# Version 2:  using the actual Dose instead of the estimated dose by DSB
# Version 3: it reads the number of DSB and combination directly from output of Ianik. We don't
# need to run Matlab anymore to read compiled output
# V2: June 2019
# V3: July 2020 - During COVID-19....

setwd(getwd())
require('lattice')
require('ggplot2')
require('ggpubr')
source("LET_Z_E_Beth.R")
source("read_C_simulation_Ianik.R")
source("fitted_death.R")
source("fitted_mutation.R")

setwd(getwd())
Z_ar = c(1,6,8,10,14,18,22,26) # 8 elements
Zst_ar=c("H","C","O","Ne","Si","Ar","Ti","Fe")
E_ar = c(10,50,100,200,400,800,1000,1600) # 8 energies
D_ar = c(0.5, 1, 2, 4)
D_str = c('0.5','1','2','4') # 4 doses
cell_line = 'HF19' # Either HF19 (human) or V79 (hamster)
if (cell_line == 'V79') {
  dom_size_ar = c('0.8','2.25')
} 
if (cell_line == 'HF19') {
  dom_size_ar = c('0.4','2.25')
}


# X-ray survival (Bedford, human breast MFC10A cells) and mutation (Kiefer, Chinese hamster V79 cells) data
Bedford = data.frame(array(c(c(0,0.25,0.5,1,2,4),c(1,0.94,0.87,0.7,0.41,0.088)),dim=c(6,2)))
colnames(Bedford) = c("dose","surv")
Kiefer = data.frame(array(c(c(0,0.5,1,2,3,4,5),c(0,1.36,3.25,8.4,15.5,24.4,35.3)*1.0e-6),dim=c(7,2)))
colnames(Kiefer) = c("dose","mut")

# Check parameters for death and mutation
cnt_ar = 1
test_ar = array(NaN,dim=c(64,8))
cnt=1
for (i_z in 1:8) {
  for (i_E in 1:8) {
    LET = LET_Z_E_Beth(E_ar[i_E],Z_ar[i_z])/1000 # LET array in keV/um
    if (LET<28) {dom_size = dom_size_ar[1]}
    else {dom_size = dom_size_ar[2]}
    gamma_x = fitted_death(cell_line,LET) # Parameter for death
    delta_x = fitted_mutation(cell_line,LET) # parameter for mutation
    test_ar[cnt,]=c(E_ar[i_E],Z_ar[i_z],LET,gamma_x[1],gamma_x[2],delta_x[1],delta_x[2],as.double(dom_size))
    cnt = cnt + 1
  }
}

# Compute  fit for survival and mutation
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
        file_path = paste('/',i_o,'/RIFResults_',dom_size,'um/RIFs_',Zst_ar[i_z],"_",E_ar[i_E],'_MeV_',D_str[i_D],'_Gy_',dom_size,'um.dat',sep="")
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
      
      # Fit mutation
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

# Plot Survival with error bar in log10 scale and fit
p<- ggplot(simu_df, aes(x=Dose, y=log10(surv_avg), color=orientation)) + 
  geom_point() +
  geom_errorbar(aes(ymin=log10(surv_avg-surv_se), ymax=log10(surv_avg+surv_se)), width=.8) +
  geom_line(aes(x=Dose, y=log10(surv_fit), color=orientation)) +
  scale_y_continuous(breaks = c(0,-1,-2),label=c(1,0.1,0.01)) +
  geom_line(data=Bedford,aes(x=dose, y=log10(surv), color="Xray-Bedford"))+
  facet_grid(Z~E_MeV) +
  scale_color_manual(name="Legend",
                     values = c("Perpendicular"= "#619CFF",
                                "Parallel"="#F8766D",
                                "Xray-Bedford"="black"))
print(p)

# Plot mutation with error bar in log10 scale and fit
p<- ggplot(simu_df, aes(x=Dose, y=(mut_avg), color=orientation)) + 
  geom_point() +
  geom_errorbar(aes(ymin=(mut_avg-mut_se), ymax=(mut_avg+mut_se)), width=.8) +
  geom_line(aes(x=Dose, y=(mut_fit), color=orientation)) +
  geom_line(data=Kiefer,aes(x=dose, y=mut, color="Xray-Kiefer"))+
  scale_y_log10() +
  facet_grid(Z~E_MeV) +
  scale_color_manual(name="Legend",
                     values = c("Perpendicular"= "#619CFF",
                                "Parallel"="#F8766D",
                                "Xray-Kiefer"="black"))
print(p)

# Plot RBE based on Bedford alpha and beta
fits=nls(log(surv) ~ -al*dose-bet*dose^2, data=Bedford, start =c(al=0.4,bet=0.1))
tmp_fits = coef(fits)
Bed_al = tmp_fits[1]
Bed_bet = tmp_fits[2]
surv_fits$LET = LET_Z_E_Beth(surv_fits$E_MeV,surv_fits$Z)/1000
surv_fits$RBE_al = surv_fits$al/Bed_al
surv_fits$RBE_bet = surv_fits$bet/Bed_bet

p<- ggplot(surv_fits, aes(x=log10(LET), y=RBE_al, color=orientation)) + 
  scale_x_continuous(breaks = c(0,1,2,3),label=c(1,10,100,1000)) +
  geom_point() +
  labs(title="RBE for Survival (alpha)",
       x ="LET (keV/um)", y = "RBE")
print(p)

p<- ggplot(surv_fits,aes(x=log10(LET), y=RBE_bet, color=orientation)) +
  scale_x_continuous(breaks = c(0,1,2,3),label=c(1,10,100,1000)) +
  geom_point() +
  labs(title="RBE for Survival (beta)",
       x ="LET (keV/um)", y = "RBE")
print(p)


# Plot mutation RBE based on Kieffer gamma1 and gamma2
fits=nls(mut ~ (ga1*dose+ga2*dose^2), data=Kiefer, start =c(ga1=1e-6,ga2=1e-8))
tmp_fits = coef(fits)
Kief_ga1 = tmp_fits[1]
Kief_ga2 = tmp_fits[2]
surv_fits$RBE_mut1 = surv_fits$gam1/Kief_ga1
surv_fits$RBE_mut2 = surv_fits$gam2/Kief_ga2

p<- ggplot(surv_fits, aes(x=log10(LET), y=RBE_mut1, color=orientation)) + 
  scale_x_continuous(breaks = c(0,1,2,3),label=c(1,10,100,1000)) +
  geom_point() +
  labs(title="RBE for mut1",
       x ="LET (keV/um)", y = "RBE")

print(p)

p<- ggplot(surv_fits,aes(x=log10(LET), y=RBE_mut2, color=orientation)) +
  scale_x_continuous(breaks = c(0,1,2,3),label=c(1,10,100,1000)) +
  geom_point() +
  labs(title="RBE for mut2",
       x ="LET (keV/um)", y = "RBE")

print(p)

# Plot RIF/Gy and ncomb vs LET - both for perpendicular and Parallel
keep_ind = simu_df$orientation == "Perpendicular"
p<- ggplot(simu_df[keep_ind,], aes(x=log10(LET), y=log2(RIF/Dose), color=factor(E_MeV), shape=factor(Dose))) +
  geom_point() +
  ylim(1.5,5.5)
print(p)

p<- ggplot(simu_df[keep_ind,], aes(x=log10(LET), y=log2(dsb/Dose), color=factor(E_MeV), shape=factor(Dose))) +
  geom_point() +
  ylim(1.5,5.5)
print(p)

p<- ggplot(simu_df[keep_ind,], aes(x=log10(LET), y=ncomb, color=factor(E_MeV), shape=factor(Dose))) +
  geom_point() +
  ylim(0,1200)
print(p)

keep_ind = simu_df$orientation == "Parallel"
p<- ggplot(simu_df[keep_ind,], aes(x=log10(LET), y=log2(RIF/Dose), color=factor(E_MeV), shape=factor(Dose))) +
  geom_point() +
  ylim(1.5,5.5)
print(p)

p<- ggplot(simu_df[keep_ind,], aes(x=log10(LET), y=log2(dsb/Dose), color=factor(E_MeV), shape=factor(Dose))) +
  geom_point() +
  ylim(1.5,5.5)
print(p)

p<- ggplot(simu_df[keep_ind,], aes(x=log10(LET), y=ncomb, color=factor(E_MeV), shape=factor(Dose))) +
  geom_point() +
  ylim(0,1200)
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


# Verify no hit for last orientation

p<- ggplot(simu_df[keep_ind,], aes(x=LET, y=nohit, color=factor(Dose))) +
  geom_point() +
  geom_point(aes(x=LET, y=dpois0, color=factor(Dose), shape=factor(Dose)))
print(p)

p<- ggplot(simu_df[keep_ind,], aes(x=LET, y=nohit, color=factor(Dose))) +
  geom_point()+
  ylim(0,1)
print(p)
