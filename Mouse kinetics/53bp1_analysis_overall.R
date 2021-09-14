# This function reads the output average file from Exogen on foci analysis and generate graphs
# S. Costes, LBNL, October 2016
#

require('lattice')
require('ggplot2')
require('polynom')
require('corrplot')
require('ggpubr')
require('plyr')

theme_set(theme_grey(base_size = 18)) # See fontsize default to 18

list_group = c(seq(206,273),seq(278,301))
let_list = c('Si 350 MeV/n','Ar 350 MeV/n','Fe 600 MeV/n','X-ray')
let_val = c(63,104,170,1) # Set let of X-ray to 1
low_dose = c(0.11,0.18,0.3,0.1)
high_dose = c(0.3,0.5,0.82,1) # Anything above or equal is considered high
xray_dose = c(0.1,1,4)
num_sample = length(list_group)
num_let = length(let_list)

# Read foci data ----------------------------------------------------------

# Read animal ID conversion into gender and species - Remove all white space and put all in caps, to reduce mismatch issues due to typos
animal_ID <- read.table('C:/Users/sylvain/Documents/Sylvain/Grant work/NASA/Runs/NSRL16C/average_well/IDvsSpecies.csv',sep=",",header=TRUE)
animal_ID$Strain = toupper(animal_ID$Strain) # Turn all ID and strain to upper case in case of typo
animal_ID$Exp_ID = gsub(" ","",toupper(animal_ID$Exp_ID),fixed = TRUE) # Make all identifier upper case and w/o space
strain_list = unique(animal_ID$Strain)
num_strain = length(strain_list)

# Read first file in initial array foci, then merge new data into it.
foci <- read.table(paste('C:/Users/sylvain/Documents/Sylvain/Grant work/NASA/Runs/NSRL16C/average_well/average_well_P',list_group[1],'.txt',sep=""),sep="\t",header=TRUE)
foci$UserID = gsub(" ","",toupper(foci$UserID),fixed = TRUE) # Make all identifier upper case and w/o space
foci$duplicate = 1
foci <- merge(foci, animal_ID, by.x='UserID', by.y='Exp_ID')

# Loop over all the other files. If a duplicate file is found, mark it with duplicate string Dup2. IF not mark it with Dup1
flag_dup = 0 # flag used to keep track of when we get the first duplicate data. Once it happens set it to 1 and create new array with duplicates on same row
# Duplicate on same row will be stored in foci_dup
# Read all the other files and concatenate time, dose and let points
for(i_sample in 2:num_sample) {
  filename = paste('C:/Users/sylvain/Documents/Sylvain/Grant work/NASA/Runs/NSRL16C/average_well/average_well_P',list_group[i_sample],'.txt',sep="")
  tmp <- read.table(filename,sep="\t",header=TRUE)
  tmp$UserID = gsub(" ","",toupper(tmp$UserID), fixed = TRUE) # Make all identifier upper case and w/o space
  duplicate_ind = as.character(foci$Radiation.Type)==as.character(tmp$Radiation.Type[1]) & 
    foci$Radiation == tmp$Radiation[1] &
    as.character(tmp$Radiation.Time[1]) == as.character(foci$Radiation.Time) # See if this new treatment was entered before
  if (sum(duplicate_ind) >0) { # found a duplicate experiment
    tmp2 = merge(foci[duplicate_ind,],tmp,by.x='UserID',by.y='UserID')
    if (flag_dup>0) {
      foci_dup = rbind(foci_dup,tmp2)
    } else {
      foci_dup = tmp2
      flag_dup = 1
    }
    tmp$duplicate = 2
  } else {
    tmp$duplicate = 1
  }
  tmp <- merge(tmp, animal_ID, by.x='UserID', by.y='Exp_ID')
  foci = rbind(foci, tmp)
}

# Convert hour factor in actual number. hour is the time post IR in hour
foci$hour = as.numeric(strsplit(as.character(foci$Radiation.Time),"hr"))
# Make sur dose is considered numeric
foci$Radiation = as.numeric(foci$Radiation)


# Label, high and low dose level for each experiment
foci$dose_level = 'unknown'  # initialize all sample to nothing
foci$fluence = 0
foci$let = 1  # set let to 1 kev/um. ONly X-ray will have this value after being done.

for (i_let in 1:num_let) {
  keep_ind = foci$Radiation.Type == let_list[i_let]
  foci[keep_ind,'let'] = let_val[i_let]
  keep_ind = foci$Radiation == 0 & foci$Radiation.Type == let_list[i_let]
  foci[keep_ind,'dose_level'] = 'Control'
  keep_ind = foci$Radiation <= low_dose[i_let] & foci$Radiation >0 & foci$Radiation.Type == let_list[i_let]
  foci[keep_ind,'dose_level'] = 'Low Dose'
  foci[keep_ind,'fluence'] = 1.1;
  keep_ind = foci$Radiation >= high_dose[i_let] & foci$Radiation.Type == let_list[i_let]
  foci[keep_ind,'dose_level']= 'High Dose'
  foci[keep_ind,'fluence'] = 3;
}
keep_ind = foci$Radiation == 0.41 # this case is only true for Si, it was a mistake
foci[keep_ind,'fluence'] = 4.1;

# Make sure order of dose level goes from low to high
foci$dose_level = factor(foci$dose_level,c('Control','Low Dose','High Dose'))

# Set fluence equals to dose for X-ray
keep_ind = foci$Radiation.Type == let_list[4] & foci$Radiation > 0
foci[keep_ind,'fluence'] = foci[keep_ind,'Radiation']


# Compute bgd, slope, square term for the dose responses of all - "Bgd","FociPerGy" --------
require(nlme)
for (i_let in 1:num_let) {
  keep_ind = (foci$Radiation.Type==let_list[i_let] & foci$num_nuc>20) # Select data for given LET
  # Determine time points for this radiation type
  time_list = sort(unique(foci[keep_ind,]$Radiation.Time)) # Mazke sure time are sorted by shortest to longest
  num_time = length(time_list)
  for (i_time in 1:num_time) {
    for (i_dup in 1:2) {
      keep_ind = (foci$Radiation.Type==let_list[i_let] & foci$num_nuc>20 & foci$Radiation.Time == time_list[i_time] & foci$duplicate==i_dup) # Select data for given LET and time
      if (i_let < num_let) { # First 3 data sets are high-LET, dose dependence is linear as foci reflects # of tracks
        fits <- lmList(avg_nfoci ~ poly(Radiation, 1, raw = TRUE) | Strain, data=foci[keep_ind,])
        if (i_let == 1 & i_time == 1 & i_dup == 1) {
          tmp_fits = cbind(coef(fits),summary(fits)$r.squared)
          tmp_fits$duplicate = i_dup
          list_fits = tmp_fits
          coef_fits = cbind(tmp_fits,time_list[i_time],let_list[i_let],rownames(tmp_fits)) # collecting only the slope data in row format will allow easier plotting of data, adding time and radiation type
          colnames(list_fits) = c(paste(let_list[i_let],time_list[i_time],"Bgd"),paste(let_list[i_let],time_list[i_time],"FociPerGy"),paste(let_list[i_let],time_list[i_time],"rsq"),paste(let_list[i_let],time_list[i_time],"duplicate"))
        } else{
          tmp_fits = cbind(coef(fits),summary(fits)$r.squared)
          tmp_fits$duplicate = i_dup
          coef_fits = rbind(coef_fits,cbind(tmp_fits,time_list[i_time],let_list[i_let],rownames(tmp_fits)))
          colnames(tmp_fits) = c(paste(let_list[i_let],time_list[i_time],"Bgd"),paste(let_list[i_let],time_list[i_time],"FociPerGy"),paste(let_list[i_let],time_list[i_time],"rsq"),paste(let_list[i_let],time_list[i_time],"duplicate"))
          list_fits = cbind(list_fits,tmp_fits[match(rownames(list_fits),rownames(tmp_fits)),])
        }
      } else { # X-ray case (i.e. i_let is 4)
        if (i_time == 1) { # add asymptoptic fit for 4 hour time point
          # fits <- nlsList(avg_nfoci ~ SSasymp(Radiation, Asym, R0, lrc) | Strain, data=foci[keep_ind,seq(1,31)])
          fits <- nlsList(avg_nfoci ~ bgd+Vmax*(1-exp(-Radiation/tau)) | Strain, data=foci[keep_ind,seq(1,31)], start =c(bgd=1,tau=0.5,Vmax=5))
          fits_para = summary(fits)$parameters # format is Strain number x 4 (value, std err, t-val, p) x 3 (bgd, tau, vmax)
          # to get vmax fitted value, it will be seq(1,15)+15x8
          tmp_fits = cbind(fits_para[seq(1:num_strain)],fits_para[3*num_strain + seq(1:num_strain)], # bgd and Pval Bgd
                           fits_para[8*num_strain + seq(1:num_strain)]/fits_para[4*num_strain + seq(1:num_strain)],fits_para[7*num_strain + seq(1:num_strain)], # Vmax/tau and tau_pval (Vmax/tau is ~slope)
                           fits_para[8*num_strain + seq(1:num_strain)],fits_para[11*num_strain + seq(1:num_strain)]) # vmax and vmax_pval
          #asym_fits = fits # keep track of the asymptotic fit in Xray at 4 hours as a separate file
          rownames(tmp_fits) = rownames(coef(fits))
          tmp_fits = cbind(tmp_fits,i_dup)
          colnames(tmp_fits) = c(paste(let_list[i_let],time_list[i_time],"Bgd"),paste(let_list[i_let],time_list[i_time],"Pval_Bgd"),paste(let_list[i_let],time_list[i_time],"FociPerGy"),paste(let_list[i_let],time_list[i_time],"Pval_FociPerGy"),paste(let_list[i_let],time_list[i_time],"FociMax"),paste(let_list[i_let],time_list[i_time],"Pval_FociMax"),paste(let_list[i_let],time_list[i_time],"duplicate"))
          # Store Bgd, and slope: Fmax/tau and Rsq estimated as (1-P)^2 into coef_fits
          tmp_coef_fits = cbind(tmp_fits[,1],tmp_fits[,3],(1-tmp_fits[,2])^2,i_dup,as.character(time_list[i_time]),let_list[i_let],rownames(tmp_fits))
          colnames(tmp_coef_fits) = colnames(coef_fits)
          coef_fits = rbind(coef_fits,tmp_coef_fits)
          list_fits = cbind(list_fits,tmp_fits[match(rownames(list_fits),rownames(tmp_fits)),])
        } else {
          # for later time point, do normal linear fit
          fits <- lmList(avg_nfoci ~ poly(Radiation, 1, raw = TRUE) | Strain, data=foci[keep_ind,seq(1,31)])
          tmp_fits = cbind(coef(fits),summary(fits)$r.squared)
          tmp_fits$duplicate = i_dup
          coef_fits = rbind(coef_fits,cbind(tmp_fits,time_list[i_time],let_list[i_let],rownames(tmp_fits)))
          colnames(tmp_fits) = c(paste(let_list[i_let],time_list[i_time],"Bgd"),paste(let_list[i_let],time_list[i_time],"FociPerGy"),paste(let_list[i_let],time_list[i_time],"rsq"),paste(let_list[i_let],time_list[i_time],"duplicate"))
          list_fits = cbind(list_fits,tmp_fits[match(rownames(list_fits),rownames(tmp_fits)),])
        }
      }
    }
  }
}
colnames(coef_fits) = c("Bgd","FociPerGy","rsq","duplicate","Time","let","Strain")
coef_fits$hour = as.numeric(strsplit(as.character(coef_fits$Time),"hr"))
coef_fits$FociPerGy = as.numeric(coef_fits$FociPerGy)
coef_fits$Bgd = as.numeric(coef_fits$Bgd)
coef_fits$rsq = as.numeric(coef_fits$rsq)

coef_fits$FociPerGy[coef_fits$FociPerGy<0]=0 # Make sure there are no negative FociPerGy
coef_fits$Bgd[coef_fits$Bgd<0]=0 # Make sure there are no negative background

# Compute simple background: i.e. mean of 0 Gy, not intercept
keep_ind = foci$Radiation == 0
avg_nfoci <- ddply(foci[keep_ind,], c("Strain", "Radiation.Type", "hour","duplicate"), summarise,
                              mean = mean(avg_nfoci, na.rm=TRUE))
colnames(avg_nfoci)= c("Strain","Radiation.Type","hour","duplicate","normal_Bgd")
                           
# add background and fitted coefficients to foci data.frame
foci <- merge(foci, coef_fits, by.x=c('Strain','hour','Radiation.Type','duplicate'), 
              by.y=c('Strain','hour','let','duplicate'))

foci <- merge(foci, avg_nfoci, by=c('Strain','hour','Radiation.Type','duplicate'))


foci$nfoci_bgdsub = foci$avg_nfoci - foci$normal_Bgd
# set negative values to 0
foci$nfoci_bgdsub[foci$nfoci_bgdsub<0] = 0
foci$norm_nfoci = foci$nfoci_bgdsub/foci$fluence/foci$avg_nuc_area*5


# Reorder list_fits to have FociMax at the very end, with duplicate next to each other
tmp_list = list_fits[,!grepl("FociMax",colnames(list_fits))]
list_fits = cbind(tmp_list,list_fits[,grepl("FociMax",colnames(list_fits))])

# Compute mean norm_foci for various time points, not segregating by dose or gender
# norm_foci for LET represents the number of RIF/um with some scaling factors
require('plyr')

keep_ind = !is.infinite(foci$norm_nfoci) & !is.na(foci$norm_nfoci)  & foci$let>100
foci_test = foci[keep_ind,c("Strain", "Radiation.Type", "hour","duplicate","nfoci_bgdsub","norm_nfoci","avg_nfoci","let","Radiation")]
mean_foci_strain_let <- ddply(foci_test, c("Strain", "Radiation.Type", "hour","let"), summarise,
                              N    = sum(!is.na(norm_nfoci)),
                              mean = mean(norm_nfoci, na.rm=TRUE),
                              sd   = sd(norm_nfoci, na.rm=TRUE),
                              se   = sd / sqrt(N)
)

# Compute average number foci/cell after bgd subtraction for X-ray
keep_ind = !is.na(foci$nfoci_bgdsub) & foci$Radiation>0 & foci$Radiation.Type==let_list[4]
foci_test = foci[keep_ind,c("Strain", "Radiation.Type", "hour","duplicate","nfoci_bgdsub","norm_nfoci","avg_nfoci","let","Radiation")]
mean_foci_strain_xray <- ddply(foci_test, c("Strain", "Radiation.Type", "hour","Radiation"), summarise,
                               N    = sum(!is.na(nfoci_bgdsub)),
                               mean = mean(nfoci_bgdsub, na.rm=TRUE),
                               mean_raw = mean(avg_nfoci, na.rm=TRUE),
                               sd   = sd(nfoci_bgdsub, na.rm=TRUE),
                               se   = sd / sqrt(N)
)



# Fit average values ------------------------------------------------------
# Saturation fits, only done at 4 hour time point
keep_ind = mean_foci_strain_let$hour == 4
# compute slope of saturation between two RIF/um data point for LET
fits <- lmList(mean ~ poly(let, 1, raw = TRUE) | Strain, data=mean_foci_strain_let[keep_ind,])
saturation_fit = coef(fits)
# Compute ratio of slope between low and high let
keep_ind = mean_foci_strain_let$hour == 4 & mean_foci_strain_let$let==let_val[2]
saturation_fit$ratio = saturation_fit[,2]/mean_foci_strain_let[keep_ind,"mean"]*let_val[2]
# compute slope of saturation between 1 and 4 Gy for Xray
keep_ind = mean_foci_strain_xray$hour == 4 & mean_foci_strain_xray$Radiation > 0.1 # keep only two highest doses for X-ray
fits <- lmList(mean ~ poly(Radiation, 1, raw = TRUE) | Strain, data=mean_foci_strain_xray[keep_ind,])
saturation_fit = cbind(saturation_fit,coef(fits))
# Compute ratio of slope between low and high dose of X-ray
keep_ind = mean_foci_strain_xray$hour == 4 & mean_foci_strain_xray$Radiation <=1
fits <- lmList(mean ~ poly(Radiation, 1, raw = TRUE) | Strain, data=mean_foci_strain_xray[keep_ind,])
saturation_fit$ratio2 = saturation_fit[,5]/coef(fits)[,2]
colnames(saturation_fit) = c("RIF_LET_int","RIF_LET_slope","RIF_LET_ratio","RIF_xray_int","RIF_xray_slope","RIF_xray_ratio")
# See what parameter correlates most
res = cor(as.matrix(saturation_fit))
p.mat <- cor.mtest(as.matrix(saturation_fit))$p
corrplot(res, type = "upper", order = "original", number.font = 5, 
         tl.col = "black", tl.srt = 45, p.mat=p.mat, insig = "p-value", sig.level=-1)

# Fit time dependence
# for high-LET. Running it with RIFmax as the 0 hr intercept (i.e. Maximum of foci/um), average
# RIF/um max was 2.1 and 1.31 for 170 and 104 keV/um. This suggest linear dependence at 0 hr
# Very close to perfect. But technically at time 0, no cluster takes place and thus all RIF
# would be detected. Thus, In second iteration, using constant instead of RIFmax, equal to
# 2.1 and 1.31 for each LET respectively. LET/80 is RIFMax
keep_ind = mean_foci_strain_let$let == let_val[2]
#fits <- nlsList(mean ~ RIFmax*exp(-hour/tau) | Strain, data=mean_foci_strain_let[keep_ind,], start =c(tau=10,RIFmax=1.5))
fits <- nlsList(mean ~ let_val[2]/80*exp(-hour/tau) | Strain, data=mean_foci_strain_let[keep_ind,], start =c(tau=10))
time_fit = coef(fits)
keep_ind = mean_foci_strain_let$let == let_val[3]
#fits <- nlsList(mean ~ RIFmax*exp(-hour/tau) | Strain, data=mean_foci_strain_let[keep_ind,], start =c(tau=10,RIFmax=1.5))
fits <- nlsList(mean ~ let_val[3]/80*exp(-hour/tau) | Strain, data=mean_foci_strain_let[keep_ind,], start =c(tau=10))
time_fit =cbind(time_fit,coef(fits))
# for X-ray - Same for X-ray, replacing RIFmax by expected # of DSB: i.e. 35*dose (in Gy)
# However because of saturation at 4 hour for both 1 and 4 Gy, kinetic is complicated. 
# Ignoring 1 and 4 Gy point at 4 hours. Keeping dose dependence
for (i_dose in 1:3){
  if (i_dose == 3) {
    keep_ind = mean_foci_strain_xray$Radiation == xray_dose[i_dose] & mean_foci_strain_xray$hour>4
  }
  else {
    keep_ind = mean_foci_strain_xray$Radiation == xray_dose[i_dose]
  }
#  fits <- nlsList(mean ~ RIFmax*exp(-hour/tau) | Strain, data=mean_foci_strain_xray[keep_ind,], start =c(tau=20,RIFmax=xray_dose[i_dose]*10))
  fits <- nlsList(mean ~ xray_dose[i_dose]*35*exp(-hour/tau) | Strain, data=mean_foci_strain_xray[keep_ind,], start =c(tau=20))
  time_fit =cbind(time_fit,coef(fits))
}
#colnames(time_fit) = c("Ar_tau","Ar_RIFmax","Fe_tau","Fe_RIFmax","X0.1Gy_tau","X0.1Gy_RIFmax","X1Gy_tau","X1Gy_RIFmax","X4Gy_tau","X4Gy_RIFmax")
colnames(time_fit) = c("Ar_tau","Fe_tau","X0.1Gy_tau","X1Gy_tau","X4Gy_tau")

# Get persistent damage at 24 and 48 hour
keep_ind = mean_foci_strain_let$hour >=24 & mean_foci_strain_let$let==let_val[2] # RIF/um for high LET
persistent_RIF = aggregate(mean_foci_strain_let[keep_ind, "mean"], 
                           list(mean_foci_strain_let$Strain[keep_ind]),
                           mean)
keep_ind = mean_foci_strain_let$hour >=24 & mean_foci_strain_let$let==let_val[3]
persistent_RIF = cbind(persistent_RIF, aggregate(mean_foci_strain_let[keep_ind, "mean"], 
                           list(mean_foci_strain_let$Strain[keep_ind]),
                           mean))

keep_ind = foci$hour >=24 & foci$let==let_val[2] # RIF/cell for high LET
persistent_RIF = cbind(persistent_RIF, aggregate(foci[keep_ind, "nfoci_bgdsub"],list(foci$Strain[keep_ind]),mean))

keep_ind = foci$hour >=24 & foci$let==let_val[3] # RIF/cell for high LET
persistent_RIF = cbind(persistent_RIF, aggregate(foci[keep_ind, "nfoci_bgdsub"],list(foci$Strain[keep_ind]),mean))

for (i_xray in c(0.1, 1, 4)) {
keep_ind = mean_foci_strain_xray$hour >=24 & mean_foci_strain_xray$Radiation == i_xray
persistent_RIF = cbind(persistent_RIF, aggregate(mean_foci_strain_xray[keep_ind, "mean"], 
                           list(mean_foci_strain_let$Strain[keep_ind]),
                           mean))
}
persistent_RIF = data.frame(persistent_RIF[,c(1,2,4,6,8,10,12,14)])
colnames(persistent_RIF) = c("Strain","ArPum_per","FePum_per","Ar_per","Fe_per","X0.1Gy_per","X1Gy_per","X4Gy_per")
persistent_RIF2 = persistent_RIF[,seq(2,8)]
rownames(persistent_RIF2) = persistent_RIF$Strain


# See what time fit parameter correlates most 
# Obsolete. Only using tau - commenting out
# tau_fit = time_fit[,c(1,3,5,7,9)]
# RIFmax = time_fit[,c(2,4,6,8,10)]
# res = cor(as.matrix(tau_fit))
# p.mat <- cor.mtest(as.matrix(tau_fit))$p
# corrplot(res, type = "upper", order = "original", number.font = 5, 
#          tl.col = "black", tl.srt = 45, p.mat=p.mat, insig = "p-value", sig.level=-1)
# res = cor(as.matrix(RIFmax))
# p.mat <- cor.mtest(as.matrix(RIFmax))$p
# corrplot(res, type = "upper", order = "original", number.font = 5, 
#          tl.col = "black", tl.srt = 45, p.mat=p.mat, insig = "p-value", sig.level=-1)

time_per = cbind(time_fit,persistent_RIF2)
res = cor(as.matrix(time_per))
p.mat <- cor.mtest(as.matrix(time_per))$p
corrplot(res, type = "upper", order = "original", number.font = 5, 
         tl.col = "black", tl.srt = 45, p.mat=p.mat, insig = "p-value", sig.level=-1)


# Compute background subtracted RIF for X-ray for all dose and RIF/um for all LET as a function of time post-IR

 # Save summary files
write.csv(time_per,"Kinetic_fits.csv")
write.csv(coef_fits,"Coef_fits.csv")
write.csv(list_fits,"list_fits.csv")
write.csv(mean_foci_strain_let,'Norm_foci_summary_particle.csv')
write.csv(mean_foci_strain_xray,'Norm_foci_summary_xray.csv')

# Plotting ----------------------------------------------------------------

# Plotting average foci number as a function of time for each strain - NEED TO ADD BGD SUB

for (i_strain in 1:num_strain) {
  keep_ind = foci$Strain == strain_list[i_strain]
  imgname = paste('FociVsTime_', strain_list[i_strain],'.tif',sep="")
  #tiff(imgname,bg="transparent",width = 1500, height = 1000, units = "px", pointsize = 14, res= 300)
  print(qplot(hour, avg_nfoci, data = foci[keep_ind,], ylim=c(0,7), 
              geom=c("jitter"), color=dose_level,linetype=factor(duplicate), shape = factor(duplicate), facets=.~Radiation.Type, size=I(2.5),
              xlab="Time post exposure (hr)",ylab="Foci/cell", main=strain_list[i_strain])
        + stat_smooth(method = 'lm', formula = 'y~log(x)', size=0.3,se=TRUE) )
  #dev.off()
}


# Plot RIF vs Dose for various LET (need to plot average only per strain on separate graph with X-ray included)
keep_ind = (foci$Radiation.Type == let_list[2] | foci$Radiation.Type == let_list[3] |
              foci$Radiation.Type == let_list[4]) &
              foci$hour<48 & foci$Radiation<=1
p<- ggplot(foci[keep_ind,], aes(x=Radiation, y=nfoci_bgdsub, shape=Strain, color=Radiation.Type)) + 
  geom_point() +
  ylim(0,4.5) +
  geom_smooth(method = "nls", se = FALSE, # Asymptotic model for early time points
              formula = y~slope*x,
              method.args = list(start=c(slope=5)))
 facet(p, facet.by = c("Radiation.Time"))

# Plot LET dependence for RIF
p<- ggplot(mean_foci_strain_let, aes(x=let, y=mean, color=Strain)) + 
  geom_point()+
  xlim(0,180) +
  ylim(0,1.5) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) +
  geom_smooth(method = "nls", se = FALSE, # Asymptotic model for early time points
              formula = y~Vmax*(1-exp(-x/tau)),
              method.args = list(start=c(tau=150,Vmax=1.5)))
facet(p, facet.by = c("hour","Strain"))

# Plot time dependence for high LET
p<- ggplot(mean_foci_strain_let, aes(x=hour, y=mean, color=Radiation.Type, shape=Radiation.Type)) + 
  geom_point()+
  xlim(0,48) +
  ylim(0,1.5) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) +
  geom_smooth(method = "nls", se = FALSE, # Asymptotic model for early time points
              formula = y~RIFmax*(exp(-x/tau))+Unrep,
             method.args = list(start=c(tau=3,RIFmax=1.5,Unrep=0.25)))
facet(p, facet.by = c("Strain"))

# Plot dose dependence for RIF - X-ray
p<- ggplot(mean_foci_strain_xray, aes(x=Radiation, y=mean, color=Strain)) + 
  geom_point()+
  xlim(0,4) +
  ylim(0,10) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) +
  geom_smooth(method = "nls", se = FALSE, # Asymptotic model for early time points
              formula = y~Vmax*(1-exp(-x/tau)),
              method.args = list(start=c(tau=10,Vmax=4)))
facet(p, facet.by = c("hour","Strain"))

tiff("Xray_dose_Avg.tif",bg="transparent",width = 2000, height = 1500, units = "px", pointsize = 14, res= 300)
p<- ggplot(mean_foci_strain_xray, aes(x=Radiation, y=mean, color=factor(Strain))) + 
  geom_point()+
  xlim(0,4) +
  ylim(0,7) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) +
facet_grid(.~hour)
print(p)
dev.off()

# Plot time dependence of X-ray
p<- ggplot(mean_foci_strain_xray, aes(x=hour, y=mean, color=factor(Radiation), shape=factor(Radiation))) + 
  geom_point()+
  xlim(0,48) +
  ylim(0,7) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) +
  geom_smooth(method = "nls", se = FALSE, # Asymptotic model for early time points
              formula = y~RIFmax*(exp(-x/tau))+Unrep,
              method.args = list(start=c(tau=10,RIFmax=5,Unrep=0.1)))
facet(p, facet.by = c("Strain"))


# Plotting RIF/um for each strain as function of time or LET for ions
# First, separating gender and doses
# first all together
# Plot time dependence of X-ray
tiff("LET_Avg.tif",bg="transparent",width = 2000, height = 2000, units = "px", pointsize = 12, res= 300)
p<- ggplot(mean_foci_strain_let, aes(x=let, y=mean, color=factor(Strain))) + 
  geom_point()+
  ylim(0,2) +
  xlim(0,180) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) +
  facet_grid(.~hour)
print(p)
dev.off()

# Second, separated
for (i_strain in 1:num_strain) {
  keep_ind = foci$Strain == strain_list[i_strain] & foci$Radiation>0 & foci$Radiation.Type == let_list[c(2,3)]
  imgname = paste('NormFociVsTime_', strain_list[i_strain],'.tif',sep="")
  tiff(imgname,bg="transparent",width = 2000, height = 2000, units = "px", pointsize = 12, res= 300)
  print(qplot(hour, norm_nfoci, data = foci[keep_ind,], ylim=c(0,2), 
              geom=c("point"), color=dose_level, linetype=Gender, shape=Gender, facets=.~Radiation.Type, size=I(5),
              xlab="Time post exposure (hr)",ylab="RIF/cell/Fluence/volume", main=strain_list[i_strain])
        + stat_smooth(method = 'lm', formula = 'y~log(x)', size=1, se=FALSE) )
  imgname = paste('NormFociVsLET_', strain_list[i_strain],'.tif',sep="")
    tiff(imgname,bg="transparent",width = 2000, height = 2000, units = "px", pointsize = 12, res= 300)
  print(qplot(let, norm_nfoci, data = foci[keep_ind,], ylim=c(0,2), xlim=c(0,180), 
              geom=c("point"), color=dose_level, linetype=Gender, shape=Gender, facets=.~hour, size=I(5),
              xlab="LET (keV/um)",ylab="RIF/cell/Fluence/volume", main=strain_list[i_strain])
        + stat_smooth(method = 'lm', formula = 'y~log(x)', size=1, se=FALSE) )
  
  dev.off()
}

# Plotting RIF/cell for each strain as function of time or dose for Xray
# First, separating gender and doses
for (i_strain in 1:num_strain) {
  keep_ind = foci$Strain == strain_list[i_strain] & foci$Radiation>0 & foci$Radiation.Type == let_list[c(4)]
  imgname = paste('nfoci_bgdsubVsTime_', strain_list[i_strain],'.tif',sep="")
    tiff(imgname,bg="transparent",width = 2500, height = 2500, units = "px", pointsize = 8, res= 300)
  print(qplot(hour, nfoci_bgdsub, data = foci[keep_ind,], ylim=c(0,8), 
              geom=c("point"), color=factor(Radiation), linetype=Gender, shape=Gender, facets=.~Radiation, size=I(5),
              xlab="Time post exposure (hr)",ylab="RIF/cell-Bgd", main=strain_list[i_strain])
        + stat_smooth(method = 'lm', formula = 'y~log(x)', size=1, se=FALSE) )

  dev.off()
  imgname = paste('nfoci_bgdsubVsDose_', strain_list[i_strain],'.tif',sep="")
    tiff(imgname,bg="transparent",width = 2000, height = 2000, units = "px", pointsize = 12, res= 300)
  print(qplot(Radiation, nfoci_bgdsub, data = foci[keep_ind,], ylim=c(0,10), xlim=c(0,4), 
              geom=c("point"), color=factor(hour), linetype=Gender, shape=Gender, facets=.~hour, size=I(5),
              xlab="Dose (Gy)",ylab="RIF/cell-Bgd", main=strain_list[i_strain])
        + stat_smooth(method = 'lm', formula = 'y~(x)', size=1, se=FALSE) )
  
    dev.off()
}

# Check duplicate accuracy
qplot(avg_nfoci.x, avg_nfoci.y, data=foci_dup, size=I(2),
      facets= Radiation.Time.x~Radiation.Type.x, color=factor(Radiation.x),
      xlab="Duplicate 1", ylab="Duplicate 2", xlim=c(0,8), ylim=c(0,8)) +
  geom_smooth(method = "lm", se = TRUE)



# Correlation graph for Bgd data, FociPerGy data and FociMax vs 4h FociPerGy LET (Omit Si and Pval)
#list_FociPerGy = list_fits[,grepl(" FociPerGy",colnames(list_fits)) | grepl(" FociMax",colnames(list_fits))]
list_FociBgd = list_fits[,grepl(" Bgd",colnames(list_fits))]
list_FociBgd = list_FociBgd[,!grepl("Si 350",colnames(list_FociBgd))]
res = cor(as.matrix(list_FociBgd))
corrplot(res, type = "upper", order = "original", number.font = 5, 
         tl.col = "black", tl.srt = 45)
persistent_RIF2$Bgd = apply(list_FociBgd[,!grepl("Fe 600 MeV/n 4 hr Bgd.1",colnames(list_FociBgd))], 1, mean,na.rm = TRUE)
# see persistence RIF correlation
res = cor(as.matrix(persistent_RIF2))
p.mat <- cor.mtest(as.matrix(persistent_RIF2))$p
corrplot(res, type = "upper", order = "original", number.font = 5, 
         tl.col = "black", tl.srt = 45, p.mat=p.mat, insig = "p-value", sig.level=-1)

list_FociBgd$std = apply(list_FociBgd[,!grepl("Fe 600 MeV/n 4 hr Bgd.1",colnames(list_FociBgd))], 1, sd,na.rm = TRUE)/sqrt(12)

list_FociPerGy = list_fits[,grepl(" FociPerGy",colnames(list_fits))]
list_FociPerGy = list_FociPerGy[,!grepl("Si 350",colnames(list_FociPerGy))]
list_FociPerGy$Bgd = apply(list_FociBgd[,!grepl("Fe 600 MeV/n 4 hr Bgd.1",colnames(list_FociBgd))], 1, mean,na.rm = TRUE)
res = cor(as.matrix(list_FociPerGy))
corrplot(res, type = "upper", order = "original", number.font = 5, 
         tl.col = "black", tl.srt = 45)

list_Saturation = saturation_fit[,c(3,6)] 
colnames(list_Saturation)=c("HZE_ratio","Xray_ratio")
#list_Saturation$Fe_4hr = apply(list_FociPerGy[,grepl("Fe 600 MeV/n 4 hr FociPer",colnames(list_FociPerGy))], 1, mean,na.rm = TRUE)
list_Saturation$Fe_24hr = apply(list_FociPerGy[,grepl("Fe 600 MeV/n 24 hr FociPer",colnames(list_FociPerGy))], 1, mean,na.rm = TRUE)
#list_Saturation$Ar_4hr = apply(list_FociPerGy[,grepl("Ar 350 MeV/n 4 hr FociPer",colnames(list_FociPerGy))], 1, mean,na.rm = TRUE)
list_Saturation$Ar_24hr = apply(list_FociPerGy[,grepl("Ar 350 MeV/n 24 hr FociPer",colnames(list_FociPerGy))], 1, mean,na.rm = TRUE)
#list_Saturation$Xray_4hr = apply(list_FociPerGy[,grepl("X-ray 4 hr FociPer",colnames(list_FociPerGy))], 1, mean,na.rm = TRUE)
list_Saturation$Xray_24hr = apply(list_FociPerGy[,grepl("X-ray 24 hr FociPer",colnames(list_FociPerGy))], 1, mean,na.rm = TRUE)
res = cor(as.matrix(list_Saturation))
p.mat <- cor.mtest(as.matrix(list_Saturation))$p
corrplot(res, type = "upper", order = "original", number.font = 5, 
         tl.col = "black", tl.srt = 45, p.mat=p.mat, insig = "p-value", sig.level=-1)
write.csv(saturation_fit,"saturation_fit.csv")
write.csv(list_saturation,"list_saturation.csv")


# PCA analysis - Created all features into one csv file. Read it and see how Strains segregate
all_features = read.csv("All_Phenotypes.csv")
rownames(all_features) = all_features[,1]
all_features = all_features[,seq(2,25)]
all_features.pca <- prcomp(all_features, center = TRUE,scale. = TRUE)
summary(all_features.pca)
require('ggfortify')
set.seed(1)
autoplot(kmeans(all_features, 4),data=all_features, label=TRUE)
         # loadings = TRUE, loadings.colour = 'blue',
         # loadings.label = TRUE, loadings.label.size = 3)

main_features = read.csv("Main_phenotypes.csv")
rownames(main_features) = main_features[,1]
main_features = main_features[,seq(2,7)]
set.seed(1)
autoplot(kmeans(main_features, 4),data=main_features, label=TRUE,
loadings = TRUE, loadings.colour = 'blue',
loadings.label = TRUE, loadings.label.size = 3)
