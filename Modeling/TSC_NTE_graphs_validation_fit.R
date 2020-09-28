# This code compares experimental data from Thacker 1977 and our model with the optimized parameters
# Eloise Pariset - September 2020 - During COVID-19....

setwd(getwd())
require('lattice')
require('ggplot2')
require('ggpubr')

setwd(getwd())
source("fitted_death.R")
source("fitted_mutation.R")
source("LET_Z_E_Beth.R")
source("read_C_simulation_Ianik.R")

# Read Thacker data from 1977 manuscript on human and hamster cells.
Thack_data =  read.csv('Thacker-1977_fit.csv',header = TRUE)

setwd(getwd())
D_ar_lowLET = c('0.5','1','2','4')
D_ar_highLET = c('0.5','1','2') # For LETs >= 50, we remove the point at 4Gy
D_ar_X = c('0.25','0.5','1','2','3','4','5','6')
dose_ar_lowLET = c(0.5,1,2,4)
dose_ar_highLET = c(0.5,1,2)
dose_ar_X = c(0.25,0.5,1,2,3,4,5,6)

cell_line = 'V79' # Either HF19 (human) or V79 (hamster)
simp = 0 # choose 1 if you want the simplified model with bet_surv = 0 and bet_mut = 0 like Thacker for HF19
if (cell_line == 'V79') {
  dom_size = '2.25'
  Z_ar = c(0,2,2,2,2,2,5,5)
  E_ar = c(0.1,8.99,5.91,2.81,1.78,1.23,10.48,4.83)
  Zst_ar=c("Photon","He","He","He","He","He","B","B")
  num_cond = 8
} 
if (cell_line == 'HF19') {
  dom_size = '2.25'
  Z_ar = c(0,2,2,2,2,2,5,5,5)
  E_ar = c(0.1,8.99,5.91,2.81,1.78,1.23,10.48,6.50,4.83)
  Zst_ar=c("Photon","He","He","He","He","He","B","B","B")
  num_cond = 9
}

# Compute  fit for survival and mutation
cnt_ar = 1
surv_fits = data.frame(t(c(0,0,0,0,0,0,0))) # initialize
simu_df = data.frame(t(c(0,0,0,0,0,0,0,0,0,0,0,0,0))) # initialize
colnames(surv_fits) = c("Z","E_MeV","LET","al_surv","bet_surv","al_mut","bet_mut")
colnames(simu_df) = c("Z","E_MeV","LET","Dose","surv_avg","surv_se","mut_avg","mut_se","surv_fit","mut_fit","RIF","dsb","ncomb")
tmp_df = simu_df


for (i_cond in 1:num_cond) {
    LET = LET_Z_E_Beth(E_ar[i_cond],Z_ar[i_cond])/1000 # computes LET array in keV/um
    if (Thack_data[i_cond,'Z']==0){ # 8 simulated doses for X-ray
      n_dose = 8
      D_ar = D_ar_X
      dose_ar = dose_ar_X
      surv_avg = c(1,1,1,1,1,1,1,1) # initiatlize survival array
      mut_avg = c(1,1,1,1,1,1,1,1) # initialize mutation array
    } else if (Thack_data[i_cond,'LET..keV.um.'] >= 50){ # we exclude the point at 4Gy for high LETs
      n_dose = 3
      D_ar = D_ar_highLET
      dose_ar = dose_ar_highLET
      surv_avg = c(1,1,1) # initiatlize survival array
      mut_avg = c(1,1,1) # initialize mutation array
    } else {
      n_dose = 4
      D_ar = D_ar_lowLET
      dose_ar = dose_ar_lowLET
      surv_avg = c(1,1,1,1) # initiatlize survival array
      mut_avg = c(1,1,1,1) # initialize mutation array
    }
    gamma_x = fitted_death(cell_line,LET) # Parameter for death
    delta_x = fitted_mutation(cell_line,LET) # parameter for mutation
    
      for (i_D in 1:n_dose) {
          if (Thack_data[i_cond,'Z']==0) { # Photon simulation
            file_path = paste('data/Photon/RIFSize/Size ',dom_size,' microns/RIFs_Photon_0.1_MeV_',D_ar_X[i_D],'_Gy.dat',sep="")
          }
          if (Thack_data[i_cond,'Z']==2 & Thack_data[i_cond,'LET..keV.um.'] < 50) { # Helium simulation and low LET
            file_path = paste('data/RIFSize_He/Size ',dom_size,' microns/RIFs_He_',Thack_data[i_cond,4],'_MeV_',D_ar_lowLET[i_D],'_Gy.dat',sep="")
          }
          if (Thack_data[i_cond,'Z']==2 & Thack_data[i_cond,'LET..keV.um.'] >= 50) { # Helium simulation and high LET
            file_path = paste('data/RIFSize_He/Size ',dom_size,' microns/RIFs_He_',Thack_data[i_cond,4],'_MeV_',D_ar_highLET[i_D],'_Gy.dat',sep="")
          }
          if (Thack_data[i_cond,'Z']==5) { # Boron simulation
            file_path = paste('data/RIFSize_B/Size ',dom_size,' microns/RIFs_B_',Thack_data[i_cond,4],'_MeV_',D_ar_highLET[i_D],'_Gy.dat',sep="")
          }   
        
        tmp=read_C_simulation_Ianik(file_path)
        surv_avg[i_D] = mean(exp(-tmp$dsb * gamma_x[1] - tmp$ncomb * gamma_x[2]))
        surv_se = sd(exp(-tmp$dsb * gamma_x[1] - tmp$ncomb * gamma_x[2]))
        mut_avg[i_D] = mean((tmp$dsb * delta_x[1] + tmp$ncomb * delta_x[2])*exp(-tmp$dsb * gamma_x[1] - tmp$ncomb * gamma_x[2]))/surv_avg[i_D]
        mut_se = sd((tmp$dsb * delta_x[1] + tmp$ncomb * delta_x[2])*exp(-tmp$dsb * gamma_x[1] - tmp$ncomb * gamma_x[2]))/surv_avg[i_D]
        tmp_df$Dose = D_ar[i_D]
        tmp_df$Z = Z_ar[i_cond]
        tmp_df$E_MeV = E_ar[i_cond]
        tmp_df$LET = LET
        tmp_df$surv_avg = surv_avg[i_D]
        tmp_df$surv_se = surv_se
        tmp_df$mut_avg = mut_avg[i_D]
        tmp_df$mut_se = mut_se
        tmp_df$RIF = mean(tmp$rif)
        tmp_df$dsb = mean(tmp$dsb)
        tmp_df$ncomb = mean(tmp$ncomb)
        simu_df = rbind(simu_df,tmp_df)
      }
      # Fit death
      if (cell_line == 'HF19' & simp == 1){
        fits=nls(log(surv_avg) ~ -al_surv*dose_ar, start =c(al_surv=0.4))
        tmp_fits = coef(fits)
        al_surv = tmp_fits[1]
      }
      else{
        fits=nls(log(surv_avg) ~ -al_surv*dose_ar-bet_surv*dose_ar^2, start =c(al_surv=0.4,bet_surv=0.1))
        tmp_fits = coef(fits)
        if (tmp_fits[2] < 0) {# Cannot allow negative terms for beta, refit with only alpha
          fits=nls(log(surv_avg) ~ -al_surv*dose_ar, start =c(al_surv=0.4))
          tmp_fits = coef(fits)
          bet_surv = 0
        } else { 
          bet_surv = tmp_fits[2]
        }      
        al_surv = tmp_fits[1]

      
      
      # Fit mutation
      if (cell_line == 'HF19' & simp == 1){
        fits=nls(mut_avg ~ al_mut*dose_ar,  start =c(al_mut=2e-6))
        tmp_fits = coef(fits)
        al_mut = tmp_fits[1]
      }
      else{
        fits=nls(mut_avg ~ (al_mut*dose_ar+bet_mut*dose_ar^2), start =c(al_mut=2e-6,bet_mut=1e-6))
        tmp_fits = coef(fits)
        if (tmp_fits[2] < 0) {# Cannot allow negative terms for gamma2, refit with  only gamma1
          fits=nls(mut_avg ~ al_mut*dose_ar,  start =c(al_mut=2e-6))
          tmp_fits = coef(fits)
          bet_mut = 0
        } else { 
          bet_mut = tmp_fits[2]
        }
        al_mut = tmp_fits[1]
        if (al_mut<0) {
          fits=nls(mut_avg ~ bet_mut*dose_ar^2,  start =c(bet_mut=2e-6))
          tmp_fits = coef(fits)
          bet_mut = tmp_fits[1]
          al_mut = 0
        }
      }
      
      
      # Save results
      surv_fits[cnt_ar,]$al_surv = al_surv
      surv_fits[cnt_ar,]$bet_surv = bet_surv
      surv_fits[cnt_ar,]$Z = Z_ar[i_cond]
      surv_fits[cnt_ar,]$E_MeV = E_ar[i_cond]      
      surv_fits[cnt_ar,]$al_mut = al_mut
      surv_fits[cnt_ar,]$bet_mut = bet_mut
      surv_fits[cnt_ar,]$LET = LET
      cnt_ar = cnt_ar + 1
      
      # put fit in simu_df for each dose
      keep_ind = simu_df$Z == Z_ar[i_cond] & simu_df$E_MeV == E_ar[i_cond]
      simu_df[keep_ind,]$surv_fit = exp(-al_surv*dose_ar - bet_surv*dose_ar^2)
      simu_df[keep_ind,]$mut_fit = al_mut*dose_ar + bet_mut*dose_ar^2
    }
  }

Thacker_val = data.frame(t(c(0,0,0,0,0,0)))
colnames(Thacker_val) = c("LET","Z","E_MeV","Dose","surv","mut")
keep_ind = Thack_data$Cell == cell_line
Thacker = Thack_data[keep_ind,]
cnt = 0
for (i_cond in 1:num_cond) {
  for (i_D in 1:4){
    Thacker_val[cnt+i_D,]$LET = Thacker$LET[i_cond]
    Thacker_val[cnt+i_D,]$Z = Thacker$Z[i_cond]
    Thacker_val[cnt+i_D,]$E_MeV = Thacker$E..MeV.n.[i_cond]
    Thacker_val[cnt+i_D,]$Dose = i_D
    Thacker_val[cnt+i_D,]$surv = exp(-Thacker$surv_alpha[i_cond]*i_D - Thacker$surv_beta[i_cond]*i_D^2)
    Thacker_val[cnt+i_D,]$mut = Thacker$mut_alpha[i_cond]*i_D + Thacker$mut_beta[i_cond]*i_D^2
  }
  cnt = cnt+4
}

simu_val = data.frame(t(c(0,0,0,0,0,0)))
colnames(simu_val) = c("LET","Z","E_MeV","Dose","surv","mut")
cnt = 0
for (i_cond in 1:num_cond) {
  for (i_D in 1:4){
    simu_val[cnt+i_D,]$LET = Thacker$LET[i_cond]
    simu_val[cnt+i_D,]$Z = surv_fits$Z[i_cond]
    simu_val[cnt+i_D,]$E_MeV = surv_fits$E_MeV[i_cond]
    simu_val[cnt+i_D,]$Dose = i_D
    simu_val[cnt+i_D,]$surv = exp(-surv_fits$al_surv[i_cond]*i_D - surv_fits$bet_surv[i_cond]*i_D^2)
    simu_val[cnt+i_D,]$mut = surv_fits$al_mut[i_cond]*i_D + surv_fits$bet_mut[i_cond]*i_D^2
  }
  cnt = cnt+4
}

# Plot survival Thacker data (red) and predictions (black)
p<- ggplot(simu_val, aes(x=Dose, y=log10(surv))) + 
  geom_point() + geom_line() +
  geom_point(data=Thacker_val, aes(x=Dose, y=log10(surv)),color="red") + geom_line(data=Thacker_val, aes(x=Dose, y=log10(surv)),color="red") +
  scale_y_continuous(breaks = c(0,-2,-4),label=c(1,10^-2,10^-4)) + 
  facet_grid(cols=vars(LET))
print(p)

# Plot mutation Thacker data (red) and predictions (black)
p<- ggplot(simu_val, aes(x=Dose, y=log10(mut))) + 
  geom_point() + geom_line() +
  geom_point(data=Thacker_val, aes(x=Dose, y=log10(mut)),color="red") + geom_line(data=Thacker_val, aes(x=Dose, y=log10(mut)),color="red") +
  facet_grid(cols=vars(LET))
print(p)
