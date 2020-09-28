# This code uses the fitted parameters determined from survival and mutation data of Thacker 1977 to extend RBE predictions in the full range of LET (0.2 to 2,400 keV/um) in parallel and perpendicular orientations.
# Eloise Pariset - September 2020 - During COVID-19....

require('IDPmisc')
require('lattice')
require('ggplot2')
require('ggpubr')

setwd(getwd())
source("LET_Z_E_Beth.R")
source("read_C_simulation_Ianik.R")
source("fitted_death.R")
source("fitted_mutation.R")

Z_ar = c(1,6,8,10,14,18,22,26) # 8 elements
Zst_ar=c("H","C","O","Ne","Si","Ar","Ti","Fe")
E_ar = c(10,50,100,200,400,800,1000,1600) # 8 energies
D_ar = c(0.5, 1, 2, 4)
D_str = c('0.5','1','2','4') # 4 doses
cell_line = 'V79' # Either HF19 (human) or V79 (hamster)
simp = 0 # choose 1 if you want the simplified model with bet = 0 like Thacker for HF19
if (cell_line == 'V79') {
  dom_size_ar = c('2.25','2.25')
  # X-ray survival data from Thacker (V79)
  Thacker_surv = data.frame(array(c(c(0,0.25,0.5,1,2,4),c(1,0.964,0.926,0.847,0.681,0.377)),dim=c(6,2)))
  colnames(Thacker_surv) = c("dose","surv")
  Thack_al = 1.43E-1 
  Thack_bet = 2.6E-2
  Thack_RBE_al_surv = data.frame(array(c(c(20,28,50,70,90,110,200,470), c(3.4,3.8,6.5,6.7,9,5.2,7.9,6.2)), dim=c(8,2)))
  colnames(Thack_RBE_al_surv) = c("LET","RBE_al")
  # X-ray mutation data from Thacker (V79)
  Thacker_mut = data.frame(array(c(c(0,0.5,1,2,3,4,5),c(0,1.97,4.36,10.4,18.2,27.8,39.0)*1.0e-6),dim=c(7,2)))
  colnames(Thacker_mut) = c("dose","mut")
  Thack_RBE_al_mut = data.frame(array(c(c(20,50,90,110,200,470), c(2.9,14,17,18,18,15)), dim=c(6,2)))
  colnames(Thack_RBE_al_mut) = c("LET","RBE_al")
  # mutation RBE data in V79 from Kiefer, 2001 for different LETs
  Kiefer_RBE_al_mut = data.frame(array(c(c(1, 110, 256, 157, 28, 25, 14, 12, 10.8, 754, 276, 238, 46, 18, 452, 335, 294, 91, 42, 1088, 1238, 407, 218, 180, 10800, 1325), c(1, 23.7, 7.7, 16, 8.15, 5.0, 5.6, 5.7, 4.0, 4.1, 19.7, 10.8, 7.2, 2.6, 12.8, 12.2, 7.1, 8.1, 7.1, 1.7, 1.9, 3.7, 6.9, 9.4, 0.2, 0.2)), dim=c(26,2)))
  colnames(Kiefer_RBE_al_mut) = c("LET","RBE_al")
  Thack_ga1 = 0.35E-5
  Thack_ga2 = 0.86E-6
} 
if (cell_line == 'HF19') {
  dom_size_ar = c('2.25','2.25')
  # X-ray survival data from Thacker (HF19)
  Thacker_surv = data.frame(array(c(c(0,0.25,0.5,1,2,4),c(1,0.820,0.673,0.453,0.205,0.0421)),dim=c(6,2)))
  colnames(Thacker_surv) = c("dose","surv")
  Thack_al = 7.92E-1
  Thack_bet = 0
  Thack_RBE_al_surv = data.frame(array(c(c(20,28,50,70,90,110,160,200,470), c(1.4,1.6,2.5,3.2,4,3.7,3.5,2.5,1.7)), dim=c(9,2)))
  colnames(Thack_RBE_al_surv) = c("LET","RBE_al")
  # X-ray mutation data from Thacker (HF19)
  Thacker_mut = data.frame(array(c(c(0,0.5,1,2,3,4,5),c(0,15.5,31.0,62.0,93.0,124,155)*1.0e-6),dim=c(7,2)))
  colnames(Thacker_mut) = c("dose","mut")
  Thack_RBE_al_mut = data.frame(array(c(c(20,28,50,70,90,110,160,200,470), c(1.3,1.6,2.9,4.3,7.1,6.5,5.7,6.7,2.8)), dim=c(9,2)))
  colnames(Thack_RBE_al_mut) = c("LET","RBE_al")
  Thack_ga1 = 3.1E-5
  Thack_ga2 = 0
}

require(nlme)
cnt_ar = 1
surv_fits = data.frame(t(c(0,0,0,0,0,0,0,0))) # initialize
simu_df = data.frame(t(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0))) # initialize
colnames(surv_fits) = c("Z","E_MeV","LET","al","bet","gam1","gam2","orientation")
colnames(simu_df) = c("Z","E_MeV","LET","Dose","surv_avg","surv_se","mut_avg","mut_se","surv_fit","mut_fit","RIF","dsb","ncomb","orientation")
tmp_df = simu_df

# Compute death and mutation for each simulated condition
surv_avg = c(1,1,1,1) # initiatlize survival array
mut_avg = c(1,1,1,1) # initialize mutation array

for (i_z in 1:8) {
  for (i_E in 1:8) {
    LET = LET_Z_E_Beth(E_ar[i_E],Z_ar[i_z])/1000 # computes LET array in keV/um
    if (LET<28) {dom_size = dom_size_ar[1]}
    else {dom_size = dom_size_ar[2]}
    gamma_x = fitted_death(cell_line,LET) # Parameter for death
    delta_x = fitted_mutation(cell_line,LET) # parameter for mutation
      for (i_o in c("Perpendicular","Parallel")) {
        for (i_D in 1:4) {
          file_path = paste(i_o,'/RIFResults_',dom_size,'um/RIFs_',Zst_ar[i_z],"_",E_ar[i_E],'_MeV_',D_str[i_D],'_Gy_',dom_size,'um.dat',sep="")
          tmp=read_C_simulation_Ianik(file_path)
          surv_avg[i_D] = mean(exp(-tmp$dsb * gamma_x[1] - tmp$ncomb * gamma_x[2]))
          surv_se = sd(exp(-tmp$dsb * gamma_x[1] - tmp$ncomb * gamma_x[2]))
          mut_avg[i_D] = mean((tmp$dsb * delta_x[1] + tmp$ncomb * delta_x[2])*exp(-tmp$dsb * gamma_x[1] - tmp$ncomb * gamma_x[2]))/surv_avg[i_D]
          mut_se = sd((tmp$dsb * delta_x[1] + tmp$ncomb * delta_x[2])*exp(-tmp$dsb * gamma_x[1] - tmp$ncomb * gamma_x[2]))/surv_avg[i_D]
          tmp_df$Dose = D_ar[i_D]
          tmp_df$Z = Z_ar[i_z]
          tmp_df$E_MeV = E_ar[i_E]
          tmp_df$LET = LET
          tmp_df$surv_avg = surv_avg[i_D]
          tmp_df$surv_se = surv_se
          tmp_df$mut_avg = mut_avg[i_D]
          tmp_df$mut_se = mut_se
          tmp_df$RIF = mean(tmp$rif)
          tmp_df$dsb = mean(tmp$dsb)
          tmp_df$ncomb = mean(tmp$ncomb)
          tmp_df$orientation = i_o
          simu_df = rbind(simu_df,tmp_df)
        }
        # Fit death
        if (cell_line == 'HF19' & simp == 1){
          fits=nls(log(surv_avg) ~ -al*D_ar, start =c(al=0.4))
          tmp_fits = coef(fits)
          al = tmp_fits[1]
          bet = 0
        }
        else{
          fits=nls(log(surv_avg) ~ -al*D_ar-bet*D_ar^2, start =c(al=0.4,bet=0.1))
          tmp_fits = coef(fits)
          if (tmp_fits[2] < 0) {# Cannot allow negative terms for beta, refit with only alpha
            fits=nls(log(surv_avg) ~ -al*D_ar, start =c(al=0.4))
            tmp_fits = coef(fits)
            bet = 0
          } else { 
            bet = tmp_fits[2]
          }      
          al = tmp_fits[1]
        }
        
        
        # Fit mutation
        if (cell_line == 'HF19' & simp == 1){
          fits=nls(mut_avg ~ ga1*D_ar,  start =c(ga1=2e-6))
          tmp_fits = coef(fits)
          ga1 = tmp_fits[1]
          ga2 = 0
        }
        else{
          fits=nls(mut_avg ~ (ga1*D_ar+ga2*D_ar^2), start =c(ga1=2e-6,ga2=1e-6))
          tmp_fits = coef(fits)
          if (tmp_fits[2] < 0) {# Cannot allow negative terms for gamma2, refit with  only gamma1
            fits=nls(mut_avg ~ ga1*D_ar,  start =c(ga1=2e-6))
            tmp_fits = coef(fits)
            ga2 = 0
          } else { 
            ga2 = tmp_fits[2]
          }
          ga1 = tmp_fits[1]
          if (ga1<0) {
            fits=nls(mut_avg ~ ga2*D_ar^2,  start =c(ga2=2e-6))
            tmp_fits = coef(fits)
            ga2 = tmp_fits[1]
            ga1 = 0
            if(ga2<0){ga2=0}
          }
        }
        
        
        # Save results
        surv_fits[cnt_ar,]$al = al
        surv_fits[cnt_ar,]$bet = bet
        surv_fits[cnt_ar,]$Z = Z_ar[i_z]
        surv_fits[cnt_ar,]$E_MeV = E_ar[i_E]      
        surv_fits[cnt_ar,]$orientation = i_o
        surv_fits[cnt_ar,]$gam1 = ga1
        surv_fits[cnt_ar,]$gam2 = ga2
        surv_fits[cnt_ar,]$LET = LET
        cnt_ar = cnt_ar + 1
        
        # put fit in simu_df for each dose
        keep_ind = simu_df$Z == Z_ar[i_z] & simu_df$E_MeV == E_ar[i_E] & simu_df$orientation == i_o
        simu_df[keep_ind,]$surv_fit = exp(-al*D_ar - bet*D_ar^2)
        simu_df[keep_ind,]$mut_fit = ga1*D_ar + ga2*D_ar^2
      }
    }
}
keep_ind = simu_df$Z>0 #remove initialization line set to 0
simu_df= simu_df[keep_ind,]


# Plot survival predictions in parallel and perpendicular and experimental Thacker X-rays data
p<- ggplot(simu_df, aes(x=Dose, y=log10(surv_avg), color=orientation)) + 
  geom_point() +
  geom_line(aes(x=Dose, y=log10(surv_fit), color=orientation)) +
  scale_y_continuous(breaks = c(0,-1,-2),label=c(1,0.1,0.01)) +
  geom_line(data=Thacker_surv,aes(x=dose, y=log10(surv), color="Xray-Thacker"))+
  facet_grid(Z~E_MeV) +
  scale_color_manual(name="Legend",
  values = c("Perpendicular"= "#619CFF",
             "Parallel"="#F8766D",
             "Xray-Thacker"="black")) +
  labs(x ="Dose (Gy)", y = "Survival fraction") + theme(legend.position = "none")
print(p)

# Plot mutation predictions in parallel and perpendicular and experimental Thacker X-rays data
p<- ggplot(simu_df, aes(x=Dose, y=(mut_avg), color=orientation)) + 
  geom_point() + xlim(0,4) +
  geom_line(aes(x=Dose, y=(mut_fit), color=orientation)) +
  geom_line(data=Thacker_mut,aes(x=dose, y=mut, color="Xray-Thacker"))+
  scale_y_log10() +
  facet_grid(Z~E_MeV) +
  scale_color_manual(name="Legend",
                     values = c("Perpendicular"= "#619CFF",
                                "Parallel"="#F8766D",
                                "Xray-Thacker"="black")) +
  labs(x ="Dose (Gy)", y = "Mutation frequency") + theme(legend.position = "none")
print(p)

# Compute survival RBE
surv_fits$LET = LET_Z_E_Beth(surv_fits$E_MeV,surv_fits$Z)/1000
surv_fits$RBE_al = surv_fits$al/Thack_al
surv_fits$RBE_bet = surv_fits$bet/Thack_bet

# Plot RBE alpha survival with Thacker experimental data
p<- ggplot(Thack_RBE_al_surv, aes(x=log10(as.numeric(LET)),y=as.numeric(RBE_al))) +
  scale_x_continuous(breaks = c(0,1,2,3),label=c(1,10,100,1000)) + 
  geom_point(shape=15) + geom_point(data=surv_fits, aes(x=log10(LET), y=RBE_al, color=orientation), show.legend = FALSE) +
  labs(x ="LET (keV/um)", y = "RBE (alpha)") + theme_bw() +
  scale_color_manual(name="Legend",
                     values = c("Perpendicular"= "#619CFF",
                                "Parallel"="#F8766D")) + theme(legend.position = "none")
print(p)

# Plot RBE beta survival
p<- ggplot(surv_fits,aes(x=log10(LET), y=RBE_bet, color=orientation)) +
  scale_x_continuous(breaks = c(0,1,2,3),label=c(1,10,100,1000)) +
  geom_point(show.legend = FALSE) +
  labs(x ="LET (keV/um)", y = "RBE (beta)") + theme_bw() +
  scale_color_manual(name="Legend",
                     values = c("Perpendicular"= "#619CFF",
                                "Parallel"="#F8766D")) + theme(legend.position = "none")
print(p)


# Compute mutation RBE
surv_fits$RBE_mut1 = surv_fits$gam1/Thack_ga1
surv_fits$RBE_mut2 = surv_fits$gam2/Thack_ga2

# Plot RBE alpha mutation with Thacker experimental data
p<- ggplot(Thack_RBE_al_mut, aes(x=log10(as.numeric(LET)),y=as.numeric(RBE_al))) +
  scale_x_continuous(breaks = c(0,1,2,3),label=c(1,10,100,1000)) + 
  geom_point(shape=15) + geom_point(data=surv_fits, aes(x=log10(LET), y=RBE_mut1, color=orientation), show.legend = FALSE) +
  labs(x ="LET (keV/um)", y = "RBE (alpha)") + theme_bw() +
  scale_color_manual(name="Legend",
                     values = c("Perpendicular"= "#619CFF",
                                "Parallel"="#F8766D")) + theme(legend.position = "none")
print(p)

# Plot RBE alpha mutation with Thacker and Kiefer experimental data (for VF79)
p<- ggplot(Thack_RBE_al_mut, aes(x=log10(as.numeric(LET)),y=as.numeric(RBE_al))) +
  scale_x_continuous(breaks = c(0,1,2,3),label=c(1,10,100,1000)) + 
  geom_point(shape=15) + geom_point(data=surv_fits, aes(x=log10(LET), y=RBE_mut1, color=orientation), show.legend = FALSE) +
  labs(x ="LET (keV/um)", y = "RBE (alpha)") + theme_bw() +
  scale_color_manual(name="Legend",
                     values = c("Perpendicular"= "#619CFF",
                                "Parallel"="#F8766D")) + theme(legend.position = "none") + geom_point(data = Kiefer_RBE_al_mut, aes(x=log10(LET),y=RBE_al), color="black", shape=2)
print(p)

# Plot RBE beta mutation
p<- ggplot(surv_fits,aes(x=log10(LET), y=RBE_mut2, color=orientation)) +
  scale_x_continuous(breaks = c(0,1,2,3),label=c(1,10,100,1000)) +
  geom_point(show.legend = FALSE) +
  labs(x ="LET (keV/um)", y = "RBE (beta)") + theme_bw() +
  scale_color_manual(name="Legend",
                     values = c("Perpendicular"= "#619CFF",
                                "Parallel"="#F8766D")) + theme(legend.position = "none")
print(p)

# Plot RIF/Gy and Ncomb vs LET - both for perpendicular and Parallel
keep_ind = simu_df$orientation == "Perpendicular"
p<- ggplot(simu_df[keep_ind,], aes(x=log10(LET), y=log2(RIF/Dose), color=factor(E_MeV), shape=factor(Dose))) +
  geom_point() + scale_x_continuous(breaks = c(0,1,2,3),label=c(1,10,100,1000)) + scale_y_discrete(limits=c("2","4","8","16","32")) +
  labs(title="0.4um Domain Size - Perpendicular",
       x ="LET (keV/um)", y = "RIF/Gy")
print(p)

p<- ggplot(simu_df[keep_ind,], aes(x=log10(LET), y=log2(dsb/Dose), color=factor(E_MeV), shape=factor(Dose))) +
  geom_point() + scale_x_continuous(breaks = c(0,1,2,3),label=c(1,10,100,1000)) + scale_y_discrete(limits=c("2","4","8","16","32")) + 
  labs(title="0.4um Domain Size - Perpendicular",
       x ="LET (keV/um)", y = "DSB/Gy")
print(p)

p<- ggplot(simu_df[keep_ind,], aes(x=log10(LET), y=ncomb, color=factor(E_MeV), shape=factor(Dose))) +
  geom_point() + scale_x_continuous(breaks = c(0,1,2,3),label=c(1,10,100,1000)) + ylim(0,1200) +
  labs(title="0.4um Domain Size - Perpendicular",
       x ="LET (keV/um)", y = "Ncomb")
print(p)

keep_ind = simu_df$orientation == "Parallel"
p<- ggplot(simu_df[keep_ind,], aes(x=log10(LET), y=log2(RIF/Dose), color=factor(E_MeV), shape=factor(Dose))) +
  geom_point() + scale_x_continuous(breaks = c(0,1,2,3),label=c(1,10,100,1000)) + scale_y_discrete(limits=c("2","4","8","16","32")) +
  labs(title="0.4um Domain Size - Parallel",
       x ="LET (keV/um)", y = "RIF/Gy")
print(p)

p<- ggplot(simu_df[keep_ind,], aes(x=log10(LET), y=log2(dsb/Dose), color=factor(E_MeV), shape=factor(Dose))) +
  geom_point() + scale_x_continuous(breaks = c(0,1,2,3),label=c(1,10,100,1000)) + scale_y_discrete(limits=c("2","4","8","16","32")) +
  labs(title="0.4um Domain Size - Parallel",
       x ="LET (keV/um)", y = "DSB/Gy")
print(p)

p<- ggplot(simu_df[keep_ind,], aes(x=log10(LET), y=ncomb, color=factor(E_MeV), shape=factor(Dose))) +
  geom_point() + scale_x_continuous(breaks = c(0,1,2,3),label=c(1,10,100,1000)) + ylim(0,1200) +
  labs(title="0.4um Domain Size - Parallel",
       x ="LET (keV/um)", y = "Ncomb")
print(p)

# Plot impact of orientaion on DSB and RIF

p<- ggplot(simu_df, aes(x=Dose, y=dsb, color=factor(E_MeV))) +
  geom_point() +
  facet_grid(Z~orientation) +
  geom_smooth(method = "nls", se = FALSE,
              formula = y~(al*x),
              method.args = list(start=c(al=35)))
print(p)

p<- ggplot(simu_df, aes(x=Dose, y=RIF, color=factor(E_MeV))) +
  geom_point() +
  facet_grid(Z~orientation) +
  geom_smooth(method = "nls", se = FALSE,
              formula = y~(al*x+bet*x^2),
              method.args = list(start=c(al=35,bet=0)))
print(p)

# Plot all individual simulations for LET = 110 (histograms)

Zst_ar = "Ti"
E_ar = 800
for (dom_size in c(0.4, 2.25)) {
  for (i_o in c("Perpendicular","Parallel")) {
    par(mfrow=c(3,4))
    par(cex=0.4)
    par(mar = c(4, 4, 0, 0), oma = c(1, 1, 1, 1))
    for (i_D in 1:4) {
      file_path = paste(i_o,'/RIFResults_',dom_size,'um/RIFs_',Zst_ar,"_",E_ar,'_MeV_',D_str[i_D],'_Gy_',dom_size,'um.dat',sep="")
      tmp=read_C_simulation_Ianik(file_path)
      hist(tmp$dsb/as.numeric(D_str[i_D]), binwidth=5, xlab="DSB/Gy", main=NULL, xlim=c(0,150)) 
    }
    for (i_D in 1:4) {
      file_path = paste(i_o,'/RIFResults_',dom_size,'um/RIFs_',Zst_ar,"_",E_ar,'_MeV_',D_str[i_D],'_Gy_',dom_size,'um.dat',sep="")
      tmp=read_C_simulation_Ianik(file_path)
      hist(tmp$rif/as.numeric(D_str[i_D]), binwidth=2, xlab="RIF/Gy", main=NULL, xlim=c(0,100)) 
    }
    for (i_D in 1:4) {
      file_path = paste(i_o,'/RIFResults_',dom_size,'um/RIFs_',Zst_ar,"_",E_ar,'_MeV_',D_str[i_D],'_Gy_',dom_size,'um.dat',sep="")
      tmp=read_C_simulation_Ianik(file_path)
      hist(tmp$ncomb, binwidth=100, xlab="Ncomb", main=NULL, xlim=c(0,2000)) 
    }
  }
}