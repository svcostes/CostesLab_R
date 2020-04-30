# Same as version 1, but using the actual Dose instead of the estimated dose by DSB
# June 2019
require('lattice')
require('ggplot2')
require('ggpubr')
require('tidyr')
source("LET_Z_E_Beth.R")
source("read_simulation.R")


Z_ar = c(1,6,8,10,14,18,22,26) # 8 tested Z values
Zst_ar=c("H","C","O","Ne","Si","Ar","Ti","Fe") # corresponding ions
E_ar = c(10,50,100,200,400,800,1000,1600) # 8 tested energies

mat=crossing(Z_ar,E_ar)
Z_ar_temp = mat$Z_ar
E_ar_temp = mat$E_ar
LET_ar_temp = LET_Z_E_Beth(E_ar_temp,Z_ar_temp)/1000
mat$LET = LET_ar_temp

LET_ar = LET_Z_E_Beth(E_ar,Z_ar)/1000 # LET array in keV/um
D_ar = c(0.5, 1, 2, 4) # 4 tested doses

# read data for parallel configuration using read_simulation.R
simu_df = read_simulation("RBE_parallel_wksp.mat") 
simu_df$orientation = "Parallel"

# read data for perpendicular configuration using read_simulation.R
simu_df2 = read_simulation("RBE_perpendicular_wksp.mat")
simu_df2$orientation = "Perpendicular"

# combine parallel and perpendicular data in simu_df
simu_df = rbind(simu_df,simu_df2)

# X-ray survival data (Bedford)
Bedford = data.frame(array(c(c(0,0.25,0.5,1,2,4),c(1,0.94,0.87,0.7,0.41,0.088)),dim=c(6,2)))
colnames(Bedford) = c("dose","surv")

# X-ray mutation data (Kiefer)
Kiefer = data.frame(array(c(c(0,0.5,1,2,3,4,5),c(0,1.36,3.25,8.4,15.5,24.4,35.3)*1.0e-6),dim=c(7,2)))
colnames(Kiefer) = c("dose","mut")


# Compute survival fit (alpha,beta) for survival
require(nlme)
surv_fits = data.frame(t(c(0,0,0,0,0))) # initialize
colnames(surv_fits) = c("Z","E_MeV","al","bet","orientation")
cnt_ar = 1
for (i_z in 1:8) {
  for (i_E in 1:8) {
    for (i_o in c("Parallel","Perpendicular")) {
      keep_ind = simu_df$Zi==i_z & simu_df$E_MeV == E_ar[i_E] & simu_df$orientation == i_o
      fits=nls(surv_avg ~ exp(-al*Dose-bet*Dose^2), data=simu_df[keep_ind,], start =c(al=0.4,bet=0.1))
      tmp_fits = coef(fits)
      if (tmp_fits[2] < 0) {# Cannot alow negative terms for beta, refit with  only alpha
        fits=nls(surv_avg ~ exp(-al*Dose), data=simu_df[keep_ind,], start =c(al=0.4))
        tmp_fits = coef(fits)
        bet = 0
      } else { 
        bet = tmp_fits[2]
      }
      al = tmp_fits[1]
      surv_fits[cnt_ar,]$al = al
      surv_fits[cnt_ar,]$bet = bet
      surv_fits[cnt_ar,]$Z = Z_ar[i_z]
      surv_fits[cnt_ar,]$E_MeV = E_ar[i_E]      
      surv_fits[cnt_ar,]$orientation = i_o
      # compute fitted survival
      dose_fit = simu_df$Dose[keep_ind]
      simu_df$surv_fit[keep_ind] = exp(-al*dose_fit-bet*dose_fit^2)
      cnt_ar = cnt_ar + 1
    }
  }
}

# Compute mutation fit (gam1,gam2) for survival
require(nlme)
surv_fits$gam1 = NaN
surv_fits$gam2 = NaN
cnt_ar = 1
for (i_z in 1:8) {
  for (i_E in 1:8) {
    for (i_o in c("Parallel","Perpendicular")) {
      keep_ind = simu_df$Zi==i_z & simu_df$E_MeV == E_ar[i_E] & simu_df$orientation == i_o
      fits=nls(mut_avg ~ (ga1*Dose+ga2*Dose^2), data=simu_df[keep_ind,], start =c(ga1=2e-6,ga2=1e-6))
      tmp_fits = coef(fits)
      if (tmp_fits[2] < 0) {# Cannot alow negative terms for gamma2, refit with  only gamma1
        fits=nls(mut_avg ~ (ga1*Dose), data=simu_df[keep_ind,], start =c(ga1=2e-6))
        tmp_fits = coef(fits)
        ga2 = 0
      } else { 
        ga2 = tmp_fits[2]
      }
      ga1 = tmp_fits[1]
      surv_fits[cnt_ar,]$gam1 = ga1
      surv_fits[cnt_ar,]$gam2 = ga2
      # compute fitted mutation
      dose_fit = simu_df$Dose[keep_ind]
      simu_df$mut_fit[keep_ind] = (ga1*dose_fit+ga2*dose_fit^2)
      cnt_ar = cnt_ar + 1
    }
  }
}

# Plot Survival with error bar in log10 scale and fit
p<- ggplot(simu_df, aes(x=Dose, y=log10(surv_avg), color=orientation)) + 
  geom_point() +
  geom_errorbar(aes(ymin=log10(surv_avg-surv_se), ymax=log10(surv_avg+surv_se)), width=.8) +
  geom_line(aes(x=Dose, y=log10(surv_fit), color=orientation)) +
  scale_y_continuous(breaks = c(0,-1,-2),label=c(1,0.1,0.01)) +
  geom_line(data=Bedford,aes(x=dose, y=log10(surv), color="Xray-Bedford"))+
  facet_grid(z~E_MeV) +
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
  facet_grid(z~E_MeV) +
  scale_color_manual(name="Legend",
                     values = c("Perpendicular"= "#619CFF",
                                "Parallel"="#F8766D",
                                "Xray-Kiefer"="black"))
print(p)
# Plot RBE based on Beford alpha and beta
fits=nls(surv ~ exp(-al*dose-bet*dose^2), data=Bedford, start =c(al=0.4,bet=0.1))
tmp_fits = coef(fits)
Bed_al = tmp_fits[1]
Bed_bet = tmp_fits[2]
surv_fits$LET = LET_Z_E_Beth(surv_fits$E_MeV,surv_fits$Z)/1000
surv_fits$RBE_al = surv_fits$al/Bed_al
surv_fits$RBE_bet = surv_fits$bet/Bed_bet

p<- ggplot(surv_fits, aes(x=LET, y=RBE_al, color=orientation)) + 
  geom_point() 
print(p)

p<- ggplot(surv_fits,aes(x=LET, y=RBE_bet, color=orientation)) +
  geom_point() 
print(p)

# Plot mutation RBE based on Kieffer gamma1 and gamma2
fits=nls(mut ~ (ga1*dose+ga2*dose^2), data=Kiefer, start =c(ga1=1e-6,ga2=1e-8))
tmp_fits = coef(fits)
Kief_ga1 = tmp_fits[1]
Kief_ga2 = tmp_fits[2]
surv_fits$RBE_mut1 = surv_fits$gam1/Kief_ga1
surv_fits$RBE_mut2 = surv_fits$gam2/Kief_ga2

p<- ggplot(surv_fits, aes(x=LET, y=RBE_mut1, color=orientation)) + 
  geom_point() 
print(p)

p<- ggplot(surv_fits,aes(x=LET, y=RBE_mut2, color=orientation)) +
  geom_point()
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
  facet_grid(z~orientation) +
  geom_smooth(method = "nls", se = FALSE,
              formula = y~(al*x),
              method.args = list(start=c(al=35)))
print(p)

p<- ggplot(simu_df, aes(x=Dose, y=RIF, color=factor(E_MeV))) +
  geom_point() +
  facet_grid(z~orientation) +
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
