read_simulation <- function(file_path) {
  # This funtion reads the Matlab .mat file in current directory file_path and extract both
  # Survival averages and Mutation averages
  #
  # Author: Sylvain Costes, Jan 2019, During Government Shutdown - NASA
  
  require("R.matlab")
  require("arrayhelpers")
  source("LET_Z_E_Beth.R")
  
  Z_ar = c(1,6,8,10,14,18,22,26) # 8 tested Z values
  Zst_ar=c("H","C","O","Ne","Si","Ar","Ti","Fe") # Corresponding ions
  E_ar = c(10,50,100,200,400,800,1000,1600) # 8 tested energies
  D_ar = c(0.5, 1, 2, 4) # 4 tested doses

  wksp = readMat(file_path)
  if (wksp$parallel.flag == 1){
    Xsection = pi*17*6/0.4
  } else {
    Xsection = pi*17^2/0.4
  }

  doseGy=wksp$doseGy
  survival=wksp$survival.D
  mutation=wksp$mutation
  dsb = wksp$dsb.count
  dose_avg=apply(doseGy,c(1,2,4),mean, na.rm = TRUE)
  surv_avg = apply(survival,c(1,2,4),mean,na.rm=TRUE)
  surv_std = apply(survival,c(1,2,4),sd,na.rm=TRUE)/sqrt(1000-apply(survival,c(1,2,4),function(x) sum(is.na(x))))
  surv_N = 1000-apply(doseGy,c(1,2,4),function(x) sum(is.na(x)))
  no_hit_N = apply(doseGy,c(1,2,4),function(x) sum(x<0.1 & !is.na(x)))
  mut_avg = apply(mutation,c(1,2,4),mean,na.rm=TRUE)
  mut_std = apply(mutation,c(1,2,4),sd,na.rm=TRUE)/sqrt(1000-apply(mutation,c(1,2,4),function(x) sum(is.na(x))))
  simu_df =array2df(dose_avg) 
  colnames(simu_df)=c("Avg_dose","Zi","E_MeV","Dose")
  simu_df$Z = Zst_ar[simu_df$Zi]
  simu_df$z = Z_ar[simu_df$Zi]
  simu_df$E_MeV = E_ar[simu_df$E_MeV]
  simu_df$Dose = D_ar[simu_df$Dose]
  simu_df$LET = LET_Z_E_Beth(simu_df$E_MeV,simu_df$z)/1000
  simu_df$Xsection = Xsection
  simu_df$fluence = simu_df$Avg_dose*Xsection/(simu_df$LET*1.6)
  simu_df$dpois0 = dpois(0,simu_df$fluence)
  tmp = array2df(surv_avg)
  simu_df$surv_avg = tmp$surv_avg
  tmp = array2df(surv_std)
  simu_df$surv_se = tmp$surv_std
  tmp = array2df(surv_N)
  simu_df$N = tmp$surv_N
  tmp = array2df(no_hit_N/surv_N)
  simu_df$nohit = tmp$`no_hit_N/surv_N`
  tmp = array2df(mut_avg)
  simu_df$mut_avg = tmp$mut_avg
  tmp = array2df(mut_std)
  simu_df$mut_se = tmp$mut_std
  tmp = array2df(dsb[,,,1])
  simu_df$dsb = tmp$`dsb[, , , 1]`
  tmp = array2df(dsb[,,,4])
  simu_df$RIF = tmp$`dsb[, , , 4]`
  tmp = array2df(dsb[,,,3])
  simu_df$ncomb = tmp$`dsb[, , , 3]`
  return(simu_df)
}
