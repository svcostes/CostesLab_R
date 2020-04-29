# This code will read the Thacker death and mutation data and fit them first to estimate the value
# for 0.5, 1, 2 and 4 Gy. Then, the code will read similar LET from simulations, number of DSB and 
# number of combinations. Apply the parameter to compute probability of death and mutation to
# compare to experimental data. A parameter sweep will be conducted and the minimum residual will 
# be used for parameter estimate
# S.V. Costes, August 2019, NASA - modified by E. Pariset in April 2020

require('nlme')
require('lattice')
require('ggplot2')
source("read_C_simulation_Ianik.R")
require(raster)  

# Read Thacker data from 1977 manuscript on human and hamster cells.
Thack_data =  read.csv('/Users/epariset/Desktop/Sylvain/Thacker-1977_fit.csv',header = TRUE)
# Ianik calibration dose points are more numerous for X-ray (8 instead of 4)
# Select cell line
cell_line = 'HF19' # Either HF19 (human) or V79 (hamster)
Thack_data = Thack_data[Thack_data$Cell==cell_line,]
Thack_data = Thack_data[order(Thack_data$Z,Thack_data$E..MeV.n.),] # order data so first dataset is for Z=0
num_E = dim(Thack_data)[1] # number of conditions tested for mutation and survival
num_LET = num_E
LET = 200 # Which LET is to be used to minimize residual of all fits. List of LETs with mutation data, for V79: 1, 20, 50, 90, 110, 200, for HF19: 1, 20, 28, 50, 70, 90, 110, 160, 200
optim_index = which(Thack_data$LET..keV.um.==LET) # fit is optimized only for this LET

# Dose array to create simulated curves
D_ar_lowLET = c('0.5','1','2','4')
D_ar_highLET = c('0.5','1','2') # For LETs >= 50, we remove the point at 4Gy
D_ar_X = c('0.25','0.5','1','2','3','4','5','6')
dose_ar_lowLET = c(0.5,1,2,4)
dose_ar_highLET = c(0.5,1,2)
dose_ar_X = c(0.25,0.5,1,2,3,4,5,6)
# survival parameters used for calibration 
gam1 = seq(0.,0.2,by=0.01) # dsb term
gam2 = seq(0.,0.2,by=0.01) # ncomb term
# mutation parameters used for calibration
lg1 = 10^-seq(2,9,by=0.1) # dsb term
lg2 = 10^-seq(3,9,by=0.1) # ncomb term

# Compute ln of probability death for Thacker data 1977
dose_mat_lowLET = rbind(dose_ar_lowLET,dose_ar_lowLET^2) # matrix is 2 x 4 (rows are d and d^2, columns are 4 doses)
dose_mat_highLET = rbind(dose_ar_highLET,dose_ar_highLET^2) # matrix is 2 x 3 (rows are d and d^2, columns are 3 doses)
dose_mat_X = rbind(dose_ar_X,dose_ar_X^2) # matrix is 2 x 8 (rows are d and d^2, columns are 8 doses)
param_mat = as.matrix(Thack_data[,c('surv_alpha','surv_beta')]) # matrix is num_LET x 2 (rows are num_E Carbon energies + 1 photon, col are al and bet)
pdeath_lowLET = (-param_mat[seq(2,num_LET),]%*% dose_mat_lowLET) # matrix is num_LET x 4 (rows are num_E Carbon energies, Columns are doses)
pdeath_highLET = (-param_mat[seq(2,num_LET),]%*% dose_mat_highLET) # matrix is num_LET x 3 (rows are num_E Carbon energies, Columns are doses)
pdeath_X = (-param_mat[1,]%*% dose_mat_X) # matrix is num_X x 8 (row is 1 photon, Columns are doses)


# compute ln of p_death for simulation for an array of possible parameters (gam1,gam2)
dom_str = c('0.4','0.6','0.8','1','1.25','1.5','1.75','2','2.25','2.5') #Range of domain sizes (10 domain sizes tested)
num_dom = length(dom_str)
pdeath_ar_lowLET = array(NaN,dim=c(4,num_dom,length(gam1),length(gam2))) # matrix storing the ln of the average pdeath for 4 doses and num_dom domain sizes, and num_E carbon energy
pdeath_ar_highLET = array(NaN,dim=c(3,num_dom,length(gam1),length(gam2))) # matrix storing the ln of the average pdeath for 3 doses and num_dom domain sizes, and num_E carbon energy
pdeath_ar_X = array(NaN,dim=c(8,num_dom,length(gam1),length(gam2))) # matrix storing the ln of the average pdeath for 8 doses and num_dom domain sizes, and 1 photon energy
residual_ar = array(NaN,dim=c(length(gam1),length(gam2))) # matrix storing residual for energy conditions (residual is on ln, to avoid biasing fit to low dose)
min_ar = array(NaN,dim=c(num_dom,3)) # matrix storing for each domain size (num_dom) and LET the coordinate of minimum and third value being the minimum
best_index = array(NaN,dim=c(num_dom,3)) # Store best index of coefficients gam1,gam2 for a given domain simulation, col1: index 1, col2: index 2, col3: residuals
best_index = data.frame(best_index) # initialize
colnames(best_index) = c("dim1","dim2","residual")

for (i_dom in 1:num_dom){ # index of domain size being tested to store result in min_ar
  dom_size = dom_str[i_dom]
  i_E = optim_index
    if (Thack_data[i_E,'Z']==0){ # 8 simulated doses for X-ray
      n_dose = 8
    } else if (Thack_data[i_E,'LET..keV.um.'] >= 50){ # we exclude the point at 4Gy for high LETs
      n_dose = 3
    } else {
      n_dose = 4
    }
    dsb_ar = matrix(ncol = n_dose, nrow = 1000) # matrix storing 3 (or 4 or 8) dose conditions and 1000 simulations for dsb
    ncomb_ar = matrix(ncol = n_dose, nrow = 1000) # matrix storing 3 (or 4 or 8) dose conditions and 1000 simulations for ncomb
    # Read all simulation for one complete dose response for a given energy of Carbon (or X-ray)
    for (i_D in 1:n_dose) {
      if (Thack_data[i_E,'Z']==0) { # Photon simulation
        file_path = paste('/Users/epariset/Desktop/Sylvain/data/Photon/RIFSize/Size ',dom_size,' microns/RIFs_Photon_0.1_MeV_',D_ar_X[i_D],'_Gy.dat',sep="")
      }
      if (Thack_data[i_E,'Z']==2 & Thack_data[i_E,'LET..keV.um.'] < 50) { # Helium simulation and low LET
        file_path = paste('/Users/epariset/Desktop/Sylvain/data/RIFSize_He/Size ',dom_size,' microns/RIFs_He_',Thack_data[i_E,4],'_MeV_',D_ar_lowLET[i_D],'_Gy.dat',sep="")
      }
      if (Thack_data[i_E,'Z']==2 & Thack_data[i_E,'LET..keV.um.'] >= 50) { # Helium simulation and high LET
        file_path = paste('/Users/epariset/Desktop/Sylvain/data/RIFSize_He/Size ',dom_size,' microns/RIFs_He_',Thack_data[i_E,4],'_MeV_',D_ar_highLET[i_D],'_Gy.dat',sep="")
      }
      if (Thack_data[i_E,'Z']==5) { # Boron simulation
        file_path = paste('/Users/epariset/Desktop/Sylvain/data/RIFSize_B/Size ',dom_size,' microns/RIFs_B_',Thack_data[i_E,4],'_MeV_',D_ar_highLET[i_D],'_Gy.dat',sep="")
      }      
      tmp=read_C_simulation_Ianik(file_path)
      dsb_ar[,i_D] = tmp$dsb
      ncomb_ar[,i_D] = tmp$ncomb
    }
    for (i_g1 in 1:length(gam1)) {
      for (i_g2 in 1:length(gam2)) {
        pdeath = exp(-dsb_ar * gam1[i_g1] - ncomb_ar * gam2[i_g2]) # computes survival probability for each of the 1000 cells, and each dose
        #pdeath_mean = log(apply(pdeath,2,mean,na.rm=TRUE))
        pdeath_mean = (apply(pdeath,2,mean,na.rm=TRUE)) # computes mean among all 1000 cells for each dose
        if (Thack_data[i_E,'Z']==0) { # Photon simulation
          pdeath_ar_X[,i_dom,i_g1,i_g2] = pdeath_mean
          residual_ar[i_g1,i_g2] = sum(abs(log(pdeath_mean)-pdeath_X),na.rm=TRUE)
        }
        else if (Thack_data[i_E,'LET..keV.um.'] >= 50){ # High LET simulation
          pdeath_ar_highLET[,i_dom,i_g1,i_g2] = pdeath_mean
          residual_ar[i_g1,i_g2] = sum(abs(log(pdeath_mean)-pdeath_highLET[i_E-1,]),na.rm=TRUE)
        }
        else { # Low LET simulation
          pdeath_ar_lowLET[,i_dom,i_g1,i_g2] = pdeath_mean
          residual_ar[i_g1,i_g2] = sum(abs(log(pdeath_mean)-pdeath_lowLET[i_E-1,]),na.rm=TRUE)
        }
        
      }
    }
  min_ar[i_dom,1:2] = which(residual_ar[,] == min(residual_ar[,]), arr.ind = TRUE)
  min_ar[i_dom,3] = min(residual_ar[,])
  tmp_index = which(residual_ar<max(min_ar[i_dom,3])*1.5,arr.ind = TRUE)
  if (nrow(tmp_index)==1){
    unique_pair=tmp_index
  }else{
    unique_pair=unique(tmp_index)
  }
  
  if (is.null(dim(unique_pair))){
    bin_index = data.frame(unique_pair[1],unique_pair[2]) # keeps track for each domain simulations which parameters fits data
    bin_index$avg_residual = NA # initialize
    best_pair = which(tmp_index[,1]==unique_pair[1] & tmp_index[,2]==unique_pair[2],arr.ind=TRUE)
    bin_index[1,'avg_residual'] = mean(residual_ar[tmp_index[best_pair[1],1],tmp_index[best_pair[1],2]])
    best_index[i_dom,] = bin_index[which(bin_index$avg_residual == min(bin_index$avg_residual))[1],]
  } else{
    bin_index = data.frame(unique_pair) # keeps track for each domain simulations which parameters fits data
    bin_index$avg_residual = NA # initialize
    for (i_index in 1:dim(unique_pair)[1]) {
      best_pair = which(tmp_index[,1]==unique_pair[i_index,1] & tmp_index[,2]==unique_pair[i_index,2],arr.ind=TRUE) # gives indices of tmp_index that have a couple (gamma1,gamma2) which gives a small enough residual
      bin_index[i_index,'avg_residual'] = mean(residual_ar[tmp_index[best_pair[1],1],tmp_index[best_pair[1],2]]) # averages value of residuals for all energy values chosen in optim_index for each couple (gamma1,gamma2)
    }
    best_index[i_dom,] = bin_index[which(bin_index$avg_residual == min(bin_index$avg_residual))[1],] # take the first minimum in case more than 1 solution
  }
}


# Plot fits for all possible domains with best gam1 and gam2

if(optim_index == 1){
  i_g1 = best_index[1]
  i_g2 = best_index[2]
  psurv_pred_X = data.frame(Psurv=factor(),Domain=factor(),Dose=factor())
  for (i_domain in 1:num_dom){
    psurv_pred_X_temp = data.frame(pdeath_ar_X[,i_domain,i_g1[i_domain,],i_g2[i_domain,]])
    psurv_pred_X_temp$domain = dom_str[i_domain]
    psurv_pred_X_temp$dose = dose_ar_X
    psurv_pred_X <- rbind(psurv_pred_X,psurv_pred_X_temp)
  }
  colnames(psurv_pred_X)=c("Predictions","Domain_Size","Dose")
  data_X <- data.frame(t(exp(pdeath_X)))
  data_X$dose = dose_ar_X
  colnames(data_X)=c("Data","Dose")
  plot_X <- merge(psurv_pred_X,data_X,by = "Dose")
  
  p<- ggplot(plot_X, aes(x=Dose, y=log10(Data))) + 
    geom_point() +
    geom_line(aes(x=Dose, y=log10(Predictions), color=factor(Domain_Size))) +
    xlab("Dose (Gy)") + ylab("log(P(survival))")
  print(p)
}


if(optim_index != 1 & LET < 50){
  LET_index = optim_index-1
  i_g1 = best_index[1]
  i_g2 = best_index[2]
  psurv_pred = data.frame(Psurv=factor(),Domain=factor(),Dose=factor())
  for (i_domain in 1:num_dom){
    psurv_pred_temp = data.frame(pdeath_ar_lowLET[,i_domain,i_g1[i_domain,],i_g2[i_domain,]])
    psurv_pred_temp$domain = dom_str[i_domain]
    psurv_pred_temp$dose = dose_ar_lowLET
    psurv_pred <- rbind(psurv_pred,psurv_pred_temp)
  }
  colnames(psurv_pred)=c("Predictions","Domain_Size","Dose")  
  data_LET <- data.frame(exp(pdeath_lowLET[optim_index-1,]))
  data_LET$dose = dose_ar_lowLET
  colnames(data_LET)=c("Data","Dose")
  plot_LET <- merge(psurv_pred,data_LET,by = "Dose")
  
  p<- ggplot(plot_LET, aes(x=Dose, y=log10(Data))) + 
    geom_point() +
    geom_line(aes(x=Dose, y=log10(Predictions), color=factor(Domain_Size))) +
    xlab("Dose (Gy)") + ylab("log(P(survival))")
  print(p)
}


if(optim_index != 1 & LET >= 50){
  LET_index = optim_index-1
  i_g1 = best_index[1]
  i_g2 = best_index[2]
  psurv_pred = data.frame(Psurv=factor(),Domain=factor(),Dose=factor())
  for (i_domain in 1:num_dom){
    psurv_pred_temp = data.frame(pdeath_ar_highLET[,i_domain,i_g1[i_domain,],i_g2[i_domain,]])
    psurv_pred_temp$domain = dom_str[i_domain]
    psurv_pred_temp$dose = dose_ar_highLET
    psurv_pred <- rbind(psurv_pred,psurv_pred_temp)
  }
  colnames(psurv_pred)=c("Predictions","Domain_Size","Dose")  
  data_LET <- data.frame(exp(pdeath_highLET[optim_index-1,]))
  data_LET$dose = dose_ar_highLET
  colnames(data_LET)=c("Data","Dose")
  plot_LET <- merge(psurv_pred,data_LET,by = "Dose")
  
  p<- ggplot(plot_LET, aes(x=Dose, y=log10(Data))) + 
    geom_point() +
    geom_line(aes(x=Dose, y=log10(Predictions), color=factor(Domain_Size))) +
    xlab("Dose (Gy)") + ylab("log(P(survival))")
  print(p)
}




##################################################################################################
# Repeat best fit for mutation with optimized gamma for cell death for all possible domain sizes #
##################################################################################################

# Unfortunately, mutation was incomplete and number of LET is limited compared to survival. Need to reduce Thack_data and all corresponding values
Thack_data2 = Thack_data[!is.na(Thack_data$mut_alpha),]
num_E = dim(Thack_data2)[1] # number of conditions tested for mutation and survival
num_LET = num_E
optim_index = which(Thack_data2$LET..keV.um.==LET) # fit is optimized only for these data. Re-run everything for different 3 groups


# Compute probability mutation for Thacker data 1977
dose_mat_lowLET = rbind(dose_ar_lowLET,dose_ar_lowLET^2) # matrix is 2 x 4 (rows are d and d^2, columns are 4 doses)
dose_mat_highLET = rbind(dose_ar_highLET,dose_ar_highLET^2) # matrix is 2 x 3 (rows are d and d^2, columns are 3 doses)
dose_mat_X = rbind(dose_ar_X,dose_ar_X^2) # matrix is 2 x 8 (rows are d and d^2, columns are 8 doses)
param_mat = as.matrix(Thack_data2[,c('mut_alpha','mut_beta')]) # matrix is num_LET x 2 (rows are num_E Carbon energies + 1 photon, col are al and bet)
pmut_lowLET = (param_mat[seq(2,num_LET),]%*% dose_mat_lowLET) # matrix is num_LET x 4 (rows are num_E Carbon energies, Columns are doses)
pmut_highLET = (param_mat[seq(2,num_LET),]%*% dose_mat_highLET) # matrix is num_LET x 3 (rows are num_E Carbon energies, Columns are doses)
pmut_X = (param_mat[1,]%*% dose_mat_X) # matrix is num_LET x 8 (row is 1 photon, Columns are doses)
pmut_ar_lowLET = array(NaN,dim=c(4,num_dom,length(lg1),length(lg2))) # matrix storing the ln of the average pmut for 4 doses and num_dom domain sizes
pmut_ar_highLET = array(NaN,dim=c(3,num_dom,length(lg1),length(lg2))) # matrix storing the ln of the average pmut for 3 doses and num_dom domain sizes
pmut_ar_X = array(NaN,dim=c(8,num_dom,length(lg1),length(lg2))) # matrix storing the ln of the average pmut for 8 doses and num_dom domain sizes
residual_mut_ar = array(NaN,dim=c(length(lg1),length(lg2))) # matrix storing residual for energy conditions (residual is on ln, to avoid biasing fit to low dose)
min_mut_ar = array(NaN,dim=c(num_dom,3)) # matrix storing for each domain size (num_dom) and LET the coordinate of minimum and third value being the minimum
best_index_mut = array(NaN,dim=c(num_dom,3)) # Store best index for a given domain simulation, col1: index 1, col2: index 2, col3: # of good fits
best_index_mut = data.frame(best_index_mut) # initialize
colnames(best_index_mut) = c("dim1","dim2","residual")


for (i_dom in 1:num_dom){ # index of domain size being tested to store result in min_ar
  dom_size = dom_str[i_dom]
  i_g1_surv = best_index[i_dom,1]
  i_g2_surv = best_index[i_dom,2]
  f_gam1 = gam1[i_g1_surv]
  f_gam2 = gam2[i_g2_surv]
  i_E = optim_index
    if (Thack_data[i_E,'Z']==0){ # 8 simulated doses for X-ray
      n_dose = 8
    } else if (Thack_data2[i_E,'LET..keV.um.'] >= 50){ # we exclude the point at 4Gy for high LETs
      n_dose = 3
    } else {
      n_dose = 4
    }
    dsb_ar = matrix(ncol = n_dose, nrow = 1000) # matrix storing 4 dose conditions and 1000 simulations for dsb
    ncomb_ar = matrix(ncol = n_dose, nrow = 1000) # matrix storing 4 dose conditions and 1000 simulations for dsb
    # Read all simulation for one complete dose response for a given energy of Carbon
    for (i_D in 1:n_dose) {
      if (Thack_data[i_E,'Z']==0) { # Photon simulation
        file_path = paste('/Users/epariset/Desktop/Sylvain/data/Photon/RIFSize/Size ',dom_size,' microns/RIFs_Photon_0.1_MeV_',D_ar_X[i_D],'_Gy.dat',sep="")
      }
      if (Thack_data[i_E,'Z']==2 & Thack_data2[i_E,'LET..keV.um.'] < 50) { # Helium simulation and low LET
        file_path = paste('/Users/epariset/Desktop/Sylvain/data/RIFSize_He/Size ',dom_size,' microns/RIFs_He_',Thack_data[i_E,4],'_MeV_',D_ar_lowLET[i_D],'_Gy.dat',sep="")
      }
      if (Thack_data[i_E,'Z']==2 & Thack_data2[i_E,'LET..keV.um.'] >= 50) { # Helium simulation and high LET
        file_path = paste('/Users/epariset/Desktop/Sylvain/data/RIFSize_He/Size ',dom_size,' microns/RIFs_He_',Thack_data[i_E,4],'_MeV_',D_ar_highLET[i_D],'_Gy.dat',sep="")
      }
      if (Thack_data[i_E,'Z']==5) { # Boron simulation
        file_path = paste('/Users/epariset/Desktop/Sylvain/data/RIFSize_B/Size ',dom_size,' microns/RIFs_B_',Thack_data[i_E,4],'_MeV_',D_ar_highLET[i_D],'_Gy.dat',sep="")
      }      
      tmp=read_C_simulation_Ianik(file_path)
      dsb_ar[,i_D] = tmp$dsb
      ncomb_ar[,i_D] = tmp$ncomb
    }
    for (i_g1 in 1:length(lg1)) {
      for (i_g2 in 1:length(lg2)) {
        pdeath = exp(-dsb_ar * f_gam1 - ncomb_ar * f_gam2)
        pmut = (dsb_ar * lg1[i_g1] + ncomb_ar * lg2[i_g2]) * pdeath # these are the ln of the probability of survival and mutation. Adding since ln instead of multiplying
        pmut_mean = apply(pmut,2,mean,na.rm=TRUE)/apply(pdeath,2,mean,na.rm=TRUE) # computes mean among all 1000 cells
        if (Thack_data2[i_E,'Z']==0) { # Photon simulation
          pmut_ar_X[,i_dom,i_g1,i_g2] = pmut_mean
          residual_mut_ar[i_g1,i_g2] = sum(abs(pmut_mean-pmut_X),na.rm=TRUE)
        }
        else if (Thack_data2[i_E,'LET..keV.um.'] >= 50){ # High LET simulation
          pmut_ar_highLET[,i_dom,i_g1,i_g2] = pmut_mean
          residual_mut_ar[i_g1,i_g2] = sum(abs(pmut_mean-pmut_highLET[i_E-1,]),na.rm=TRUE)
        }
        else{ # Low LET simulation
          pmut_ar_lowLET[,i_dom,i_g1,i_g2] = pmut_mean
          residual_mut_ar[i_g1,i_g2] = sum(abs(pmut_mean-pmut_lowLET[i_E-1,]),na.rm=TRUE)          
        }
      }
    }
  min_mut_ar[i_dom,1:2] = which(residual_mut_ar[,] == min(residual_mut_ar[,]), arr.ind = TRUE)
  min_mut_ar[i_dom,3] = min(residual_mut_ar[,])
  tmp_index = which(residual_mut_ar<max(min_mut_ar[i_dom,3])*1.5,arr.ind = TRUE) # make sure all fit will include every LET conditions (use max)
  if (nrow(tmp_index)==1){
    unique_pair=tmp_index
  }else{
    unique_pair=unique(tmp_index)
  }

  if (is.null(dim(unique_pair))){
    bin_index = data.frame(unique_pair[1],unique_pair[2]) # keeps track for each domain simulations which parameters fits data
    bin_index$avg_residual = NA # initialize
    best_pair = which(tmp_index[,1]==unique_pair[1] & tmp_index[,2]==unique_pair[2],arr.ind=TRUE)
    bin_index[1,'avg_residual'] = mean(residual_mut_ar[tmp_index[best_pair[1],1],tmp_index[best_pair[1],2]])
    best_index_mut[i_dom,] = bin_index[which(bin_index$avg_residual == min(bin_index$avg_residual))[1],]
  } else{
    bin_index = data.frame(unique_pair) # keeps track for each domain simulations which parameters fits data
    bin_index$avg_residual = NA # initialize
    for (i_index in 1:dim(unique_pair)[1]) {
      best_pair = which(tmp_index[,1]==unique_pair[i_index,1] & tmp_index[,2]==unique_pair[i_index,2],arr.ind=TRUE) # gives indices of tmp_index that have a couple (gamma1,gamma2) which gives a small enough residual
      bin_index[i_index,'avg_residual'] = mean(residual_mut_ar[tmp_index[best_pair[1],1],tmp_index[best_pair[1],2]]) # averages value of residuals for all energy values chosen in optim_index for each couple (gamma1,gamma2)
    }
    best_index_mut[i_dom,] = bin_index[which(bin_index$avg_residual == min(bin_index$avg_residual))[1],] # take the first minimum in case more than 1 solution
  }
}



# Plot fits for all possible domains with best gam1 and gam2

if(optim_index == 1){
  i_g1 = best_index_mut[1]
  i_g2 = best_index_mut[2]
  pmut_pred_X = data.frame(Pmut=factor(),Domain=factor(),Dose=factor())
  for (i_domain in 1:num_dom){
    pmut_pred_X_temp = data.frame(pmut_ar_X[,i_domain,i_g1[i_domain,],i_g2[i_domain,]])
    pmut_pred_X_temp$domain = dom_str[i_domain]
    pmut_pred_X_temp$dose = dose_ar_X
    pmut_pred_X <- rbind(pmut_pred_X,pmut_pred_X_temp)
  }
  colnames(pmut_pred_X)=c("Predictions","Domain_Size","Dose")
  data_X <- data.frame(t(pmut_X))
  data_X$dose = dose_ar_X
  colnames(data_X)=c("Data","Dose")
  plot_X <- merge(pmut_pred_X,data_X,by = "Dose")
  
  p<- ggplot(plot_X, aes(x=Dose, y=Data)) + 
    geom_point() +
    geom_line(aes(x=Dose, y=Predictions, color=factor(Domain_Size))) +
    xlab("Dose (Gy)") + ylab("P(mutation)")
  print(p)
}


if(optim_index != 1 & LET < 50){
  LET_index = optim_index-1
  i_g1 = best_index_mut[1]
  i_g2 = best_index_mut[2]
  pmut_pred = data.frame(Pmut=factor(),Domain=factor(),Dose=factor())
  for (i_domain in 1:num_dom){
    pmut_pred_temp = data.frame(pmut_ar_lowLET[,i_domain,i_g1[i_domain,],i_g2[i_domain,]])
    pmut_pred_temp$domain = dom_str[i_domain]
    pmut_pred_temp$dose = dose_ar_lowLET
    pmut_pred <- rbind(pmut_pred,pmut_pred_temp)
  }
  colnames(pmut_pred)=c("Predictions","Domain_Size","Dose")  
  data_LET <- data.frame(pmut_lowLET[optim_index-1,])
  data_LET$dose = dose_ar_lowLET
  colnames(data_LET)=c("Data","Dose")
  plot_LET <- merge(pmut_pred,data_LET,by = "Dose")
  
  p<- ggplot(plot_LET, aes(x=Dose, y=Data)) + 
    geom_point() +
    geom_line(aes(x=Dose, y=Predictions, color=factor(Domain_Size))) +
    xlab("Dose (Gy)") + ylab("P(mutation)")
  print(p)
}


if(optim_index != 1 & LET >= 50){
  LET_index = optim_index-1
  i_g1 = best_index_mut[1]
  i_g2 = best_index_mut[2]
  pmut_pred = data.frame(Pmut=factor(),Domain=factor(),Dose=factor())
  for (i_domain in 1:num_dom){
    pmut_pred_temp = data.frame(pmut_ar_highLET[,i_domain,i_g1[i_domain,],i_g2[i_domain,]])
    pmut_pred_temp$domain = dom_str[i_domain]
    pmut_pred_temp$dose = dose_ar_highLET
    pmut_pred <- rbind(pmut_pred,pmut_pred_temp)
  }
  colnames(pmut_pred)=c("Predictions","Domain_Size","Dose")  
  data_LET <- data.frame(pmut_highLET[optim_index-1,])
  data_LET$dose = dose_ar_highLET
  colnames(data_LET)=c("Data","Dose")
  plot_LET <- merge(pmut_pred,data_LET,by = "Dose")
  
  p<- ggplot(plot_LET, aes(x=Dose, y=Data)) + 
    geom_point() +
    geom_line(aes(x=Dose, y=Predictions, color=factor(Domain_Size))) +
    xlab("Dose (Gy)") + ylab("P(mutation)")
  print(p)
}
