# This function reads the output average file from Exogen on foci analysis and generate graphs
# S. Costes, LBNL, October 2016
# Same version as the original one but with output focused on fitting individual animals
# instead of grouping by strain. This is required for GWAS
# Removed duplicate separation because it would group duplicate 1 together for various dose or time points, making it an unjustified 
# assocation. Nov 2018
#

require('lattice')
require('ggplot2')
require('polynom')
require('corrplot')
require('ggpubr')
require('plyr')

theme_set(theme_grey(base_size = 12)) # See fontsize default to 18

# Removing 0 Gy Xray control 294 and 301 (24 and 4 hr respectively)
# number of foci look outliers: either too high or too low. Keeping the other control for each Time-pt
list_group = c(seq(206,273),seq(278,293), seq(295,300)) 

let_list = c('Si 350 MeV/n','Ar 350 MeV/n','Fe 600 MeV/n','X-ray')
let_list_sub = c('Si','Ar','Fe','Xray')
let_val = c(63,104,170,1) # Set let of X-ray to 1
low_dose = c(0.11,0.18,0.3,0.1)
high_dose = c(0.3,0.5,0.82,1) # Anything above or equal is considered high
xray_dose = c(0.1,1,4)
num_sample = length(list_group)
num_let = length(let_list)

# Read foci data ----------------------------------------------------------

# Read animal ID conversion into gender and species - Remove all white space and put all in caps, to reduce mismatch issues due to typos
animal_ID <- read.table('D:/Sylvain_Backup/Grant work/NASA/Runs/NSRL16C/average_well/IDvsSpecies.csv',sep=",",header=TRUE)
animal_ID$Strain = toupper(animal_ID$Strain) # Turn all ID and strain to upper case in case of typo
animal_ID$Exp_ID = gsub(" ","",toupper(animal_ID$Exp_ID),fixed = TRUE) # Make all identifier upper case and w/o space
strain_list = unique(animal_ID$Strain)
num_strain = length(strain_list)
UserID_list = unique(animal_ID$Exp_ID)
num_UserID = length(UserID_list)

# Read first file in initial array foci, then merge new data into it.
foci <- read.table(paste('../NSRL16C/average_well/average_well_P',list_group[1],'.txt',sep=""),sep="\t",header=TRUE)
foci$UserID = gsub(" ","",toupper(foci$UserID),fixed = TRUE) # Make all identifier upper case and w/o space
foci$duplicate = 1
foci <- merge(foci, animal_ID, by.x='UserID', by.y='Exp_ID')

# Loop over all the other files. If a duplicate file is found, mark it with duplicate string Dup2. IF not mark it with Dup1
flag_dup = 0 # flag used to keep track of when we get the first duplicate data. Once it happens set it to 1 and create new array with duplicates on same row
# Duplicate on same row will be stored in foci_dup
# Read all the other files and concatenate time, dose and let points
for(i_sample in 2:num_sample) {
  filename = paste('../NSRL16C/average_well/average_well_P',list_group[i_sample],'.txt',sep="")
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

# 
# # Check duplicate accuracy
# keep_ind = foci_dup$Radiation.Type.x != "Si 350 MeV/n"
# qplot(avg_nfoci.x, avg_nfoci.y, data=foci_dup[keep_ind,], size=I(2),
#       facets= Radiation.Time.x~Radiation.Type.x, color=factor(plate.x),
#       xlab="Duplicate 1", ylab="Duplicate 2", xlim=c(0,8), ylim=c(0,8)) +
#   geom_smooth(method = "lm", se = TRUE)
# 
# # check intensity bias
# qplot(plate,avg_fitc,data=foci)
# qplot(avg_fitc,avg_nfoci,data=foci,facets=Radiation~hour, color=duplicate)


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

# Create summary experimental conditions by plate#
condition_summary=unique(foci[,c(4,5,6,7)])
for (dose_level in c('Control','Low Dose','High Dose')) {
keep_ind = foci$Radiation.Type != "Si 350 MeV/n" & foci$dose_level == dose_level
p<- ggplot(foci[keep_ind,], aes(x=plate, y=avg_nfoci, fill=Radiation.Type )) + 
  geom_bar(stat = "summary", fun.y = "mean" , width=.8, position = "dodge")  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(hour~.) +
  ggtitle(dose_level) 
print(p)
}

# Compute bgd, slope, square term for the dose responses of all - "Bgd","FociPerGy" --------
require(nlme)
for (i_let in 1:num_let) {
  keep_ind = (foci$Radiation.Type==let_list[i_let] & foci$num_nuc>20) # Select data for given LET
  # Determine time points for this radiation type
  time_list = sort(unique(foci[keep_ind,]$Radiation.Time)) # Make sure time are sorted by shortest to longest
  num_time = length(time_list)
  for (i_time in 1:num_time) {
     
      keep_ind = (foci$Radiation.Type==let_list[i_let] & foci$num_nuc>20 & foci$Radiation.Time == time_list[i_time] ) # Select data for given LET and time
      if (i_let < num_let) { # First 3 data sets are high-LET, dose dependence is linear as foci reflects # of tracks
        fits <- lmList(avg_nfoci ~ poly(Radiation, 1, raw = TRUE) | UserID, data=foci[keep_ind,])
        avg_backup <- ddply(foci[keep_ind,], "UserID", summarise, mean = mean(avg_nfoci, na.rm=TRUE))
        tmp_fits = coef(fits)
        
        # identify negative slope or NA value. For those, slope is 0 and intercept is average
        keep_ind = tmp_fits[,2]<0 | is.na(tmp_fits[,1]) | is.na(tmp_fits[,2])
        tmp_fits[keep_ind,1] = avg_backup$mean[keep_ind]
        tmp_fits[keep_ind,2] = 0
        if (i_let == 1 & i_time == 1) {
          list_fits = tmp_fits
          coef_fits = cbind(tmp_fits,time_list[i_time],let_list[i_let],rownames(tmp_fits)) # collecting only the slope data in row format will allow easier plotting of data, adding time and radiation type
          colnames(list_fits) = c(paste(let_list[i_let],time_list[i_time],"Bgd"),paste(let_list[i_let],time_list[i_time],"FociPerGy"))
        } else{
          coef_fits = rbind(coef_fits,cbind(tmp_fits,time_list[i_time],let_list[i_let],rownames(tmp_fits)))
          colnames(tmp_fits) = c(paste(let_list[i_let],time_list[i_time],"Bgd"),paste(let_list[i_let],time_list[i_time],"FociPerGy"))
          list_fits2 = merge(list_fits,tmp_fits,by="row.names",all=TRUE)
          row.names(list_fits2) = list_fits2$Row.names
          list_fits = list_fits2[,seq(2,length(list_fits2))]
        }
      } else { # X-ray case (i.e. i_let is 4)
        if (i_time == 1) { # exclude 4 Gy for 4 hour time point
          keep_ind = (foci$Radiation.Type==let_list[i_let] & foci$num_nuc>20 & foci$Radiation<=1 & foci$Radiation.Time == time_list[i_time] ) # Select data for given LET and time
          } 
        # for later time point, do normal linear fit
        fits <- lmList(avg_nfoci ~ poly(Radiation, 1, raw = TRUE) | UserID, data=foci[keep_ind,seq(1,31)])
        avg_backup <- ddply(foci[keep_ind,], "UserID", summarise, mean = mean(avg_nfoci, na.rm=TRUE))
        tmp_fits = coef(fits)
        
        # identify negative slope or NA value. For those, slope is 0 and intercept is average
        keep_ind = tmp_fits[,2]<0 | is.na(tmp_fits[,1]) | is.na(tmp_fits[,2])
        tmp_fits[keep_ind,1] = avg_backup$mean[keep_ind]
        tmp_fits[keep_ind,2] = 0
        
        coef_fits = rbind(coef_fits,cbind(tmp_fits,time_list[i_time],let_list[i_let],rownames(tmp_fits)))
        colnames(tmp_fits) = c(paste(let_list[i_let],time_list[i_time],"Bgd"),paste(let_list[i_let],time_list[i_time],"FociPerGy"))
        list_fits2 = merge(list_fits,tmp_fits,by="row.names",all=TRUE)
        row.names(list_fits2) = list_fits2$Row.names
        list_fits = list_fits2[,seq(2,length(list_fits2))]
      }
    
  }
}
list_fits = merge(list_fits,animal_ID,by.x="row.names",by.y=c('Exp_ID'))

colnames(coef_fits) = c("Bgd","FociPerGy","hour","Radiation.Type","UserID")
coef_fits$Bgd = as.numeric(coef_fits$Bgd)
coef_fits$FociPerGy = as.numeric(coef_fits$FociPerGy)
coef_fits$hour = as.numeric(strsplit(as.character(coef_fits$hour),"hr"))

# Compute simple background: i.e. mean of 0 Gy, not intercept - Note background is not duplicate specific
keep_ind = foci$Radiation == 0
avg_nfoci <- ddply(foci[keep_ind,], c("UserID", "let", "Radiation.Type", "hour","Strain"), summarise,
                              mean = mean(avg_nfoci, na.rm=TRUE))
colnames(avg_nfoci)= c("UserID","let","Radiation.Type","hour","Strain","normal_Bgd")
coef_fits = merge(coef_fits,avg_nfoci,by=c('UserID','hour','Radiation.Type'), all=TRUE)


coef_fits = merge(coef_fits,animal_ID,by.x=c('UserID',"Strain"),by.y=c('Exp_ID',"Strain"))

write.csv(coef_fits,'phenotypes_by_UserID_v2.csv')

# add background and fitted coefficients to foci data.frame
foci <- merge(foci, coef_fits[,c(1,3,seq(5,8))], by.x=c('UserID','hour','let'), 
              by.y=c('UserID','hour','let'))


foci$nfoci_bgdsub = foci$avg_nfoci - foci$Bgd
foci$norm_nfoci = foci$nfoci_bgdsub/foci$fluence/foci$avg_nuc_area*5

# # Reorder list_fits to have FociMax at the very end, with duplicate next to each other
# tmp_list = list_fits[,!grepl("FociMax",colnames(list_fits))]
# list_fits = cbind(tmp_list,list_fits[,grepl("FociMax",colnames(list_fits))])

# Compute mean norm_foci for various time points, not segregating by dose or gender
# norm_foci for LET represents the number of RIF/um with some scaling factors - Exclude outliers (no_out)
# Also, do no include any point where doses lead to less than background. Remove these, make no sense.

require('plyr')


keep_ind = !is.infinite(foci$norm_nfoci) & !is.na(foci$norm_nfoci)  & foci$let>100 & foci$nfoci_bgdsub>0
mean_foci_UserID_let <- ddply(foci[keep_ind,], c("UserID", "let", "hour"), summarise,
                              N    = sum(!is.nan(norm_nfoci)),
                              mean = mean(norm_nfoci, na.rm=TRUE),
                              sd   = sd(norm_nfoci, na.rm=TRUE),
                              se   = sd / sqrt(N),
                              N_no_out = sum(!is.nan(norm_nfoci[norm_nfoci<mean+sd & norm_nfoci>mean-sd])),
                              mean_no_out = mean(norm_nfoci[norm_nfoci<mean+sd & norm_nfoci>mean-sd],na.rm=TRUE)
                              )

# Compute average number foci/cell after bgd subtraction for X-ray. Cannot exclude outliers, most groups have only 2 values
keep_ind = !is.na(foci$nfoci_bgdsub) & foci$Radiation>0 & foci$Radiation.Type==let_list[4] & foci$nfoci_bgdsub>0
mean_foci_UserID_xray <- ddply(foci[keep_ind,], c("UserID", "let", "hour","Radiation"), summarise,
                               N    = sum(!is.na(nfoci_bgdsub)),
                               mean = mean(nfoci_bgdsub, na.rm=TRUE),
                               mean_raw = mean(avg_nfoci, na.rm=TRUE),
                               sd   = sd(nfoci_bgdsub, na.rm=TRUE),
                               se   = sd / sqrt(N)
)
mean_foci_UserID_xray$meanPGy = mean_foci_UserID_xray$mean/mean_foci_UserID_xray$Radiation

# Note data are sorted by UserID alphabetically. This is critical for the next section or everything will be
# scrambled.
mean_foci_UserID_let=merge(mean_foci_UserID_let,animal_ID[,c('Exp_ID','Strain','Gender')],by.x='UserID',by.y='Exp_ID')
mean_foci_UserID_xray=merge(mean_foci_UserID_xray,animal_ID[,c('Exp_ID','Strain','Gender')],by.x='UserID',by.y='Exp_ID')
write.csv(mean_foci_UserID_let,'mean_foci_UserID_let.csv')
write.csv(mean_foci_UserID_xray,'mean_foci_UserID_xray.csv')

# Fit average values ------------------------------------------------------
# Saturation fits, only done at 4 hour time point


# Compute ratio of slope of RIF/um vs LET of slope2 over slope1: 1-slope2/slope1, if fully saturated
# slope 2 is 0 and parameter is max and equal to 1. If no saturation, slope2=slope1 and parameter is 0
# IF parameter is negative, this model is wrong and probably data are linear. Set it to 0.
keep_ind1 = mean_foci_UserID_let$hour == 4 & mean_foci_UserID_let$let==let_val[2]
keep_ind2 = mean_foci_UserID_let$hour == 4 & mean_foci_UserID_let$let==let_val[3]
slope1 = mean_foci_UserID_let$mean_no_out[keep_ind1]/let_val[2]
slope2 = (mean_foci_UserID_let$mean_no_out[keep_ind2]-mean_foci_UserID_let$mean_no_out[keep_ind1])/(let_val[3]-let_val[2])
# initialize dataframe. This will work ONLY if we have corresponding values for both LEt for each animal
saturation_fit = mean_foci_UserID_let[keep_ind1,c(1,3)] 
saturation_fit$let_slope=slope1
saturation_fit$let_saturation=1-slope2/slope1
keep_ind = saturation_fit$let_saturation <0 # search for no saturated data. Then replace by average slope and set sat to 0
saturation_fit$let_slope[keep_ind] = 0.5*(slope1[keep_ind]+slope2[keep_ind])
saturation_fit$let_saturation[keep_ind] = 0

# compute slope and saturation for Xray. Xray data are always saturated, no exception
keep_ind = mean_foci_UserID_xray$hour == 4 & mean_foci_UserID_xray$Radiation == 0.1 # 0.1 Gy point
keep_ind1 = mean_foci_UserID_xray$hour == 4 & mean_foci_UserID_xray$Radiation == 1 # 1 Gy point
keep_ind2 = mean_foci_UserID_xray$hour == 4 & mean_foci_UserID_xray$Radiation == 4 # 4 Gy point
slope1 = (mean_foci_UserID_xray$mean[keep_ind]/0.1 + mean_foci_UserID_xray$mean[keep_ind1])/2
slope2 = (mean_foci_UserID_xray$mean[keep_ind2]-mean_foci_UserID_xray$mean[keep_ind1])/3
tmp_fit = mean_foci_UserID_xray[keep_ind1,c(1,3)] 
tmp_fit$xray_slope = slope1
tmp_fit$xray_saturation=1-slope2/slope1
#if slope2 is negative, full saturation and slope 1 is adjusted
keep_ind = slope2<0
tmp_fit$xray_slope[keep_ind] = slope1[keep_ind]+0.75*slope2[keep_ind]
tmp_fit$xray_saturation[keep_ind] = 1

saturation_fit = merge(saturation_fit,tmp_fit[,c(1,3,4)],by = 'UserID',all=TRUE)
saturation_fit=merge(saturation_fit,animal_ID[,c('Exp_ID','Strain','Gender')],by.x='UserID',by.y='Exp_ID')
write.csv(saturation_fit,'saturation_fit_by_UserID.csv')

# summarize by strain
mean_foci_strain_let <- ddply(mean_foci_UserID_let, c("Strain", "let", "hour"), summarise,
                              N    = sum(!is.nan(mean_no_out)),
                              mean = mean(mean_no_out, na.rm=TRUE),
                              sd   = sd(mean_no_out, na.rm=TRUE),
                              se   = sd / sqrt(N)
)
mean_foci_strain_xray <- ddply(mean_foci_UserID_xray, c("Strain", "hour","Radiation"), summarise,
                              N    = sum(!is.nan(mean)),
                              sd   = sd(mean, na.rm=TRUE),
                              se   = sd / sqrt(N),
                              mean = mean(mean, na.rm=TRUE)
)
mean_saturation_fit <- ddply(saturation_fit, c("Strain"), summarise,
                              N_let    = sum(!is.na(let_slope)),
                             let_s_se = sd(let_slope, na.rm=TRUE)/sqrt(N_let),
                             let_sat_se = sd(let_saturation, na.rm=TRUE)/sqrt(N_let),
                             let_slope = mean(let_slope, na.rm=TRUE),
                             let_saturation = mean(let_saturation, na.rm=TRUE),
                             N_xray    = sum(!is.na(xray_slope)),
                             xray_s_se = sd(xray_slope, na.rm=TRUE)/sqrt(N_xray),
                             xray_sat_se = sd(xray_saturation, na.rm=TRUE)/sqrt(N_xray),
                             xray_slope = mean(xray_slope, na.rm=TRUE),
                             xray_saturation = mean(xray_saturation, na.rm=TRUE)
)
mean_coef_fits <- ddply(coef_fits, c("Strain","Radiation.Type", "hour"), summarise,
                        N    = sum(!is.na(FociPerGy)),
                        FociPerGy_se = sd(FociPerGy, na.rm=TRUE)/sqrt(N),
                        Bgd_se = sd(Bgd, na.rm=TRUE)/sqrt(N),
                        normal_Bgd_se =  sd(normal_Bgd, na.rm=TRUE)/sqrt(N),
                        FociPerGy = mean(FociPerGy, na.rm=TRUE),
                        Bgd = mean(Bgd, na.rm=TRUE),
                        normal_Bgd =  mean(normal_Bgd, na.rm=TRUE)
)

# create array table for correlation for coef_fits
list_hour = c(4,8,24,48) # list of hours
for (i_let in 2:num_let) {
  for (i_time in 1:4) {
    keep_ind = mean_coef_fits$Radiation.Type==let_list[i_let] & mean_coef_fits$hour == list_hour[i_time] 
    tmp_FociPerGy = unique(mean_coef_fits[keep_ind,c(1,8)])
    colnames(tmp_FociPerGy) = c("Strain",paste(let_list_sub[i_let],'_',list_hour[i_time],'Hr'))
    tmp_Bgd = unique(mean_coef_fits[keep_ind,c(1,9)])
    colnames(tmp_Bgd) = c("Strain",paste(let_list_sub[i_let],'_',list_hour[i_time],'Hr'))
    tmp_normal_Bgd = unique(mean_coef_fits[keep_ind,c(1,10)])
    colnames(tmp_normal_Bgd) = c("Strain",paste(let_list_sub[i_let],'_',list_hour[i_time],'Hr'))
    if (i_let == 2 & i_time == 1) {
      list_FociBgd = tmp_Bgd
      list_FociPerGy = tmp_FociPerGy
      list_normal_Bgd = tmp_normal_Bgd
    } else {
      list_FociBgd = merge(list_FociBgd,tmp_Bgd,by='Strain',all=TRUE)
      list_FociPerGy =  merge(list_FociPerGy, tmp_FociPerGy,by='Strain',all=TRUE)
      list_normal_Bgd =  merge(list_normal_Bgd, tmp_normal_Bgd,by='Strain',all=TRUE)
    }
  }
}
rownames(list_FociBgd) = list_FociBgd$Strain
list_FociBgd =subset(list_FociBgd, select=-c(1,11)) # Get rid of strain column and 8hr Xray since blank
rownames(list_FociPerGy) = list_FociPerGy$Strain
list_FociPerGy =subset(list_FociPerGy, select=-c(1,11))
rownames(list_normal_Bgd) = list_normal_Bgd$Strain
list_normal_Bgd =subset(list_normal_Bgd, select=-c(1,11))
## Figure 7 - Start

# Plot saturation value per strain (bargraph)
tmp_set  = mean_saturation_fit[,c(1,6,4)]
ord_ind = order(tmp_set$let_saturation)
tmp_set = tmp_set[ord_ind,]
tmp_set$Strain = factor(tmp_set$Strain, level=tmp_set$Strain)
tmp_set2 = mean_saturation_fit[,c(1,11,9)]
tmp_set$Radiation = "HZE particles"
tmp_set2$Radiation = "X-rays"
colnames(tmp_set) = c('Strain',"Saturation","Se","Radiation")
colnames(tmp_set2) = c('Strain',"Saturation","Se","Radiation")
tmp_set = rbind(tmp_set,tmp_set2)
p<- ggplot(tmp_set, aes(x=Strain, y=Saturation, fill=Radiation)) + 
  geom_bar(stat="identity", width=.8, position = "dodge")  +
  geom_errorbar(aes(ymin=Saturation-Se, ymax= Saturation+Se), width=.8, position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(p)
# See what parameter correlates most
colnames(mean_saturation_fit) = c( "Strain","N_let","let_s_se","let_sat_se","LET_slope1","LET_saturation",
                                   "N_xray","xray_s_se","xray_sat_se","Xrays_slope1","Xrays_saturation")
res = cor(as.matrix(mean_saturation_fit[,c(5,6,10,11)]))
p.mat <- cor.mtest(as.matrix(mean_saturation_fit[,c(5,6,10,11)]))$p
corrplot(res, type = "upper", order = "original", number.font = 5, 
         tl.col = "black", tl.srt = 45, p.mat=p.mat, insig = "p-value", sig.level=-1)
# Plot Figure 7 - End

# Plot RIF/Gy versus Strains - Figure 5C
# first sort strain, from more foci to less at 24 hours for Fe
keep_ind = mean_coef_fits$Radiation.Type == 'Fe 600 MeV/n' & mean_coef_fits$hour == 24
tmp_set = mean_coef_fits[keep_ind,]
ord_ind = order(tmp_set$FociPerGy)
tmp_set = tmp_set[ord_ind,]
mean_coef_fits$Strain_sorted = factor(mean_coef_fits$Strain, levels = tmp_set$Strain)
keep_ind = mean_coef_fits$Radiation.Type != 'Si 350 MeV/n' & (mean_coef_fits$hour == 4 | mean_coef_fits$hour == 24)
scale_factor = 8 # using this factor to make sure the minimu value is above 1 to avoid bar looking negatives for values between 0 and 1
p<- ggplot(mean_coef_fits[keep_ind,], aes(x=Strain_sorted, y=scale_factor*FociPerGy, fill=Radiation.Type)) + 
  geom_bar(stat="identity", width=.8, position = "dodge")  +
  geom_errorbar(aes(ymin=scale_factor*(FociPerGy-FociPerGy_se), ymax=scale_factor*(FociPerGy+FociPerGy_se)), width=.8, position = "dodge") +
  scale_y_continuous(trans='log2',limits=scale_factor*c(0.125,8),breaks = scale_factor*c(0.25,0.5,1,2,4,8), labels = c("0.25","0.5","1","2","4","8")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(hour~.)
print(p)

# Plot Bgd versus Strains - Figure 5B
# first sort strain, from more foci to less overall
mean_coef_fits_bgd <- ddply(coef_fits, c("Strain"), summarise,
                            normal_Bgd =  mean(normal_Bgd, na.rm=TRUE)
)
ord_ind = order(mean_coef_fits_bgd$normal_Bgd)
mean_coef_fits_bgd = mean_coef_fits_bgd[ord_ind,]
mean_coef_fits$Strain_sorted = factor(mean_coef_fits$Strain, levels = mean_coef_fits_bgd$Strain)
keep_ind = mean_coef_fits$Radiation.Type != 'Si 350 MeV/n' & (mean_coef_fits$hour == 4 | mean_coef_fits$hour == 24)
scale_factor = 8 # using this factor to make sure the minimu value is above 1 to avoid bar looking negatives for values between 0 and 1
p<- ggplot(mean_coef_fits[keep_ind,], aes(x=Strain_sorted, y=scale_factor*normal_Bgd, fill=Radiation.Type)) + 
  geom_bar(stat="identity", width=.8, position = "dodge")  +
  geom_errorbar(aes(ymin=scale_factor*(normal_Bgd-normal_Bgd_se), ymax=scale_factor*(normal_Bgd+normal_Bgd_se)), width=.8, position = "dodge") +
  scale_y_continuous(trans='log2',limits=scale_factor*c(0.125,8),breaks = scale_factor*c(0.25,0.5,1,2,4,8), labels = c("0.25","0.5","1","2","4","8")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(hour~.)
print(p)

# Plot Bgd versus Strains - Figure 5B(alternative - y-intercept)
# first sort strain, from more foci to less overall
mean_coef_fits_bgd <- ddply(coef_fits, c("Strain"), summarise,
                            normal_Bgd =  mean(normal_Bgd, na.rm=TRUE)
)
ord_ind = order(mean_coef_fits_bgd$normal_Bgd)
mean_coef_fits_bgd = mean_coef_fits_bgd[ord_ind,]
mean_coef_fits$Strain_sorted = factor(mean_coef_fits$Strain, levels = mean_coef_fits_bgd$Strain)
keep_ind = mean_coef_fits$Radiation.Type != 'Si 350 MeV/n' & (mean_coef_fits$hour == 4 | mean_coef_fits$hour == 24)
scale_factor = 8 # using this factor to make sure the minimu value is above 1 to avoid bar looking negatives for values between 0 and 1
p<- ggplot(mean_coef_fits[keep_ind,], aes(x=Strain_sorted, y=scale_factor*Bgd, fill=Radiation.Type)) + 
  geom_bar(stat="identity", width=.8, position = "dodge")  +
  geom_errorbar(aes(ymin=scale_factor*(Bgd-Bgd_se), ymax=scale_factor*(Bgd+Bgd_se)), width=.8, position = "dodge") +
  scale_y_continuous(trans='log2',limits=scale_factor*c(0.125,8),breaks = scale_factor*c(0.25,0.5,1,2,4,8), labels = c("0.25","0.5","1","2","4","8")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(hour~.)
print(p)
# *******************
  
# Fit time dependence
# for high-LET. Running it with RIFmax as the 0 hr intercept (i.e. Maximum of foci/um), average
# RIF/um max was 1.8 and 1.35 for 170 and 104 keV/um. This shows some saturation at 0 hr intercept (more foci at lower LET respectively)
# Using these values for fitting as RIFmax
keep_ind = mean_foci_strain_let$let == let_val[2] & mean_foci_strain_let$hour < 48
#fits <- nlsList(mean ~ RIFmax*exp(-hour/tau) | Strain, data=mean_foci_strain_let[keep_ind,], start =c(tau=10,RIFmax=1.5))
fits <- nlsList(mean ~ 1.8*exp(-hour/tau) | Strain, data=mean_foci_strain_let[keep_ind,], start =c(tau=10))
time_fit = coef(fits)


keep_ind = mean_foci_strain_let$let == let_val[3] & mean_foci_strain_let$hour < 48
#fits <- nlsList(mean ~ RIFmax*exp(-hour/tau) | Strain, data=mean_foci_strain_let[keep_ind,], start =c(tau=10,RIFmax=1.5))
fits <- nlsList(mean ~ 1.345*exp(-hour/tau) | Strain, data=mean_foci_strain_let[keep_ind,], start =c(tau=10))
time_fit =cbind(time_fit,coef(fits))

# for X-ray - Same for X-ray, replacing RIFmax by expected # of DSB: i.e. 35*dose (in Gy)
# However because of saturation at 4 hour for both 1 and 4 Gy, kinetic is complicated. 
# Ignoring 1 and 4 Gy point at 4 hours. Keeping dose dependence
for (i_dose in 1:3){
  if (i_dose > 1) {
    keep_ind = mean_foci_strain_xray$Radiation == xray_dose[i_dose] & mean_foci_strain_xray$hour>4
  }
  else {
    keep_ind = mean_foci_strain_xray$Radiation == xray_dose[i_dose]
  }
#  fits <- nlsList(mean ~ RIFmax*exp(-hour/tau) | Strain, data=mean_foci_strain_xray[keep_ind,], start =c(tau=20,RIFmax=xray_dose[i_dose]*10))
  fits <- nlsList(mean ~ xray_dose[i_dose]*35*exp(-hour/tau) | Strain, data=mean_foci_strain_xray[keep_ind,c(1,2,3,7)], start =c(tau=10))
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

# Plotting average foci number as a function of time for each strain
qplot(hour, nfoci_bgdsub, data = foci, ylim=c(0,7), 
            geom=c("jitter"), color=dose_level, facets=Strain~Radiation.Type, size=I(2.5),
            xlab="Time post exposure (hr)",ylab="Foci/cell")
    #  + stat_smooth(method = 'lm', formula = 'y~log(x)', size=0.3,se=TRUE) )

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
# Figure 3B
keep_ind = (foci$Radiation.Type == let_list[2] | foci$Radiation.Type == let_list[3]) &
              foci$hour<48 & foci$Radiation<=1
p<- ggplot(foci[keep_ind,], aes(x=Radiation, y=nfoci_bgdsub, shape=Strain, color=Radiation.Type)) + 
  geom_point() +
  ylim(-1,4.5) +
  geom_smooth(method = "nls", se = FALSE, # Asymptotic model for early time points
              formula = y~slope*x,
              method.args = list(start=c(slope=5)))
 facet(p, facet.by = c("Radiation.Time"))
# Figure 3A 
 keep_ind = (foci$Radiation.Type == let_list[2] | foci$Radiation.Type == let_list[3]) &
   foci$hour<48 & foci$Radiation<=1 & foci$Strain == "C57"
 p<- ggplot(foci[keep_ind,], aes(x=Radiation, y=avg_nfoci, color=Radiation.Type)) + 
   geom_point() +
   geom_smooth(method = "nls", se = FALSE, # Asymptotic model for early time points
               formula = y~slope*x+bgd,
               method.args = list(start=c(slope=5,bgd=0.5)))
 facet(p, facet.by = c("Radiation.Time"))
 
# Plot LET dependence for RIF/um
 keep_ind = mean_foci_strain_let$hour == 4
p<- ggplot(mean_foci_strain_let[keep_ind,], aes(x=let, y=mean, color=Strain)) + 
  geom_point()+
  ylim(0,1.5) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) +
  geom_smooth(method = "lm", se = FALSE, formula = y~x, size = 1 ) +
  scale_x_continuous(limits=c(0,180),breaks = c(0,104,170)) 
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

# Second, separated - Figure 2 for Part II
theme_set(theme_grey(base_size = 12)) # See fontsize default to 18
# Removing one obvious outlier for CC019 (48 hour - low dose) - Ideally should filter norm_nfoci for outliers, per condition
keep_ind = foci$Strain == 'CC019' & foci$Radiation==0.18 & foci$Radiation.Type == let_list[2] & foci$hour>24 &foci$norm_nfoci>0.5
foci[keep_ind,]$norm_nfoci=NA
for (i_strain in 1:num_strain) {
  keep_ind = foci$Strain == strain_list[i_strain] & foci$Radiation>0 & foci$Radiation.Type == let_list[c(2,3)]
  imgname = paste('NormFociVsTime_', strain_list[i_strain],'.tif',sep="")
  tiff(imgname,bg="transparent",width = 2000, height = 2000, units = "px", pointsize = 12, res= 400)
  print(qplot(hour, norm_nfoci, data = foci[keep_ind,], ylim=c(0,2), 
              geom=c("point"), color=dose_level, fill=dose_level, linetype=Gender, shape=Gender, facets=.~Radiation.Type, size=I(3),
              xlab="Time post exposure (hr)",ylab="RIF/cell/Fluence/volume", main=strain_list[i_strain])
        + scale_shape_manual(values=c(21,22))
        + scale_color_manual(values=c("black","black"))
        + stat_smooth(method = 'lm', formula = 'y~log(x)', size=1, se=FALSE) 
  )
  dev.off()
  imgname = paste('NormFociVsLET_', strain_list[i_strain],'.tif',sep="")
  keep_ind = foci$Strain == strain_list[i_strain] & foci$Radiation>0 & foci$Radiation.Type == let_list[c(2,3)] & foci$hour<48
    tiff(imgname,bg="transparent",width = 2000, height = 2000, units = "px", pointsize = 12, res= 400)
  print(qplot(let, norm_nfoci, data = foci[keep_ind,], ylim=c(0,2), xlim=c(0,180), 
              geom=c("point"), color=dose_level, linetype=Gender, shape=Gender, fill=Gender, facets=.~hour, size=I(3),
              xlab="LET (keV/um)",ylab="RIF/cell/Fluence/volume", main=strain_list[i_strain])
        + scale_shape_manual(values=c(16, 15))
        + stat_smooth(method = 'lm', formula = 'y~lin(x)', size=1, se=FALSE) 
  )
  dev.off()
}

# Plotting RIF/cell for each strain as function of time or dose for Xray
# First, separating gender and doses
for (i_strain in 1:num_strain) {
  keep_ind = foci$Strain == strain_list[i_strain] & foci$Radiation>0 & foci$Radiation.Type == let_list[c(4)]
  imgname = paste('nfoci_bgdsubVsTime_', strain_list[i_strain],'.tif',sep="")
    tiff(imgname,bg="transparent",width = 2700, height = 2000, units = "px", pointsize = 12, res= 400)
  print(qplot(hour, nfoci_bgdsub, data = foci[keep_ind,], ylim=c(0,8), 
              geom=c("point"), fill = I("grey"), linetype=Gender, shape=Gender, facets=.~Radiation, size=Gender,
              xlab="Time post exposure (hr)",ylab="RIF/cell-Bgd", main=strain_list[i_strain])
        + scale_shape_manual(values=c(21,22))
        + scale_size_manual(values=c(4,3))
        + stat_smooth(method = 'lm', formula = 'y~log(x)', size=1, se=FALSE, color=I("black")) )

  dev.off()
  imgname = paste('nfoci_bgdsubVsDose_', strain_list[i_strain],'.tif',sep="")
    tiff(imgname,bg="transparent",width = 2700, height = 2000, units = "px", pointsize = 12, res= 400)
  print(qplot(Radiation, nfoci_bgdsub, data = foci[keep_ind,], ylim=c(0,10), xlim=c(0,4), 
              geom=c("point"), color=factor(hour), linetype=Gender, shape=Gender, facets=.~hour, size=I(3),
              xlab="Dose (Gy)",ylab="RIF/cell-Bgd", main=strain_list[i_strain])
        + scale_shape_manual(values=c(21,22))
        + scale_size_manual(values=c(4,3))
        + stat_smooth(method = 'lm', formula = 'y~(x)', size=1, se=FALSE, color=I("black")) )
  
    dev.off()
}

# End of Figure 2, Part II



# Correlation graph for Bgd data, FociPerGy data and FociMax vs 4h FociPerGy LET (Omit Si and Pval)
# Figure 4A
res = cor(as.matrix(list_FociBgd), use="complete.obs") 
p.mat <- cor.mtest(as.matrix(list_FociBgd))$p
corrplot(res, type = "upper", number.font = 5, p.mat=p.mat, insig = "p-value",
         tl.col = "black", tl.srt = 45,  order = "original", addrect = 3, sig.level = .05)

# normal BGD is not used for figures
res = cor(as.matrix(list_normal_Bgd), use="complete.obs") 
p.mat <- cor.mtest(as.matrix(list_normal_Bgd))$p
corrplot(res, type = "upper", number.font = 5, p.mat=p.mat, insig = "p-value",
         tl.col = "black", tl.srt = 45,  order = "original", addrect = 3, sig.level = .05)

# Figure 4B
res = cor(as.matrix(list_FociPerGy))
p.mat <- cor.mtest(as.matrix(list_FociPerGy))$p
corrplot(res, type = "upper", number.font = 5, p.mat=p.mat, insig = "p-value",
         tl.col = "black", tl.srt = 45,  order = "original", addrect = 3, sig.level = .05)

persistent_RIF2$Bgd = apply(list_FociBgd[,!grepl("Fe 600 MeV/n 4 hr Bgd.1",colnames(list_FociBgd))], 1, mean,na.rm = TRUE)
# see persistence RIF correlation
res = cor(as.matrix(persistent_RIF2))
p.mat <- cor.mtest(as.matrix(persistent_RIF2))$p
corrplot(res, type = "upper", order = "original", number.font = 5, 
         tl.col = "black", tl.srt = 45, p.mat=p.mat, insig = "p-value", sig.level=-1)


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
