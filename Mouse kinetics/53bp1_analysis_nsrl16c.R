# This function reads the output average file from Exogen on foci analysis and generate graphs
# S. Costes, LBNL, October 2016
#

require('lattice')
require('ggplot2')

list_group = c(seq(206,273),seq(278,292))
let_list = c('Si 350 MeV/n','Ar 350 MeV/n','Fe 600 MeV/n','X-ray')
low_dose = c(0.11,0.18,0.3,0.1)
high_dose = c(0.3,0.5,0.82,1)
num_sample = length(list_group)
num_let = length(let_list)

# Read animal ID conversion into gender and species - Remove all white space and put all in caps, to reduce mismatch issues due to typos
animal_ID <- read.table('D:/Sylvain/Grant work/NASA/Runs/NSRL16C/average_well/IDvsSpecies.csv',sep=",",header=TRUE)
animal_ID$Strain = toupper(animal_ID$Strain) # Turn all ID and strain to upper case in case of typo
animal_ID$Exp_ID = gsub(" ","",toupper(animal_ID$Exp_ID),fixed = TRUE) # Make all identifier upper case and w/o space
strain_list = unique(animal_ID$Strain)
num_strain = length(strain_list)

# Read first file in initial array foci, then merge new data into it.
foci <- read.table(paste('D:/Sylvain/Grant work/NASA/Runs/NSRL16C/average_well/average_well_P',list_group[1],'.txt',sep=""),sep="\t",header=TRUE)
foci$UserID = gsub(" ","",toupper(foci$UserID),fixed = TRUE) # Make all identifier upper case and w/o space
foci$duplicate = 1
foci <- merge(foci, animal_ID, by.x='UserID', by.y='Exp_ID')

# Loop over all the other files. If a duplicate file is found, mark it with duplicate string Dup2. IF not mark it with Dup1
flag_dup = 0 # flag used to keep track of when we get the first duplicate data. Once it happens set it to 1 and create new array with duplicates on same row
# Duplicate on same row will be stored in foci_dup
# Read all the other files and concatenate time, dose and let points
for(i_sample in 2:num_sample) {
  filename = paste('D:/Sylvain/Grant work/NASA/Runs/NSRL16C/average_well/average_well_P',list_group[i_sample],'.txt',sep="")
  tmp <- read.table(filename,sep="\t",header=TRUE)
  tmp$UserID = gsub(" ","",toupper(tmp$UserID), fixed = TRUE) # Make all identifier upper case and w/o space
  duplicate_ind =  foci$Radiation == tmp$Radiation[1] & as.character(tmp$Radiation.Time[1]) == as.character(foci$Radiation.Time)
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
foci$dose_level = 'Control'  # initialize all sample to control

for (i_let in 1:num_let) {
  keep_ind = foci$Radiation == low_dose[i_let] & foci$Radiation.Type == let_list[i_let]
  foci[keep_ind,30] = 'Low Dose'
  keep_ind = foci$Radiation >= high_dose[i_let] & foci$Radiation.Type == let_list[i_let]
  foci[keep_ind,30]= 'High Dose'
}
# Make sure order of dose level goes from low to high



foci$dose_level = factor(foci$dose_level,c('Control','Low Dose','High Dose'))

# Check duplicate accuracy
qplot(avg_nfoci.x, avg_nfoci.y, data=foci_dup, size=I(2),
      facets= Radiation.x~., 
      xlab="Duplicate 1", ylab="Duplicate 2") +
 geom_smooth(method = "lm", se = TRUE)


# Run PCA on foci table
# source("D:/Sylvain/Grant work/NASA/Runs/NSRL16C/pca_foci.R") # Load function to concatenate file
# pca_foci(foci)


# Plotting various relationship

setwd("C:/Users/Sylvain/Documents/R_images")
for (i_strain in 1:num_strain) {
  keep_ind = foci$Strain == strain_list[i_strain]
  imgname = paste('FociVsTime_', strain_list[i_strain],'.tif',sep="")
  #tiff(imgname,bg="transparent",width = 1500, height = 1000, units = "px", pointsize = 14, res= 300)
  print(qplot(hour, avg_foci_no_outl, data = foci[keep_ind,], ylim=c(0,7), 
              geom=c("jitter"), color=dose_level,facets=.~Radiation.Type, size=I(1),
              xlab="Time post exposure (hr)",ylab="Foci/cell", main=strain_list[i_strain])
        + stat_smooth(method = 'lm', formula = 'y~log(x)', size=0.3,se=TRUE) )
  #dev.off()
}


qplot(Radiation.Time, avg_foci_no_outl, data=foci, geom=c("jitter"), shape=Radiation.Type, color=dose_level, 
            facets=.~Strain, size=I(2),
            xlab="Hours", ylab="Foci/cell") 

keep_ind = foci$Radiation.Type=="X-ray" & foci$dose_level=="High Dose"
qplot(hour, avg_foci_no_outl, data = foci[keep_ind,], 
      geom = c("point", "smooth"), color=UserID,facets=.~Strain,
      xlab="Time post exposure (hr)",ylab="Foci/cell", ylim=c(0,8))

keep_ind = (foci$Radiation.Type=="Fe 600 MeV/n" & foci$num_nuc>20)
qplot(Radiation, avg_nfoci, data = foci[keep_ind,], 
      geom = c("point", "lm"),facets=Radiation.Time~Strain,
      xlab="Dose (Gy)",ylab="Foci/cell")

qplot(Radiation, num_nuc, data = foci[keep_ind,], 
      geom = c("point", "smooth"),facets=Radiation.Time~Strain,
      xlab="Dose (Gy)",ylab="Cell number")

qplot(Radiation, avg_p2a, data = foci[keep_ind,], 
      geom = c("point", "smooth"),facets=Radiation.Time~Strain,
      xlab="Dose (Gy)",ylab="P2A")

qplot(plate, avg_nuc_dapi, data=foci,geom=c('boxplot'), color=Radiation.Type)
qplot(plate, avg_fitc_bg, data=foci,geom=c('boxplot'), color=Radiation.Type)
qplot(plate, num_nuc, data=foci,geom=c('boxplot'), color=Radiation.Type)
qplot(Radiation.Time, num_nuc, data=foci,geom=c('jitter'), color=Radiation.Type,facets=dose_level~.)
qplot(Strain, avg_foci_no_outl, data=foci, geom=c('boxplot','jitter'),color=Radiation.Type,facets=dose_level~Radiation.Time)
# Sort strain based only on background foci level
keep_ind = foci$Radiation>0 & foci$hour>8
mean_foci_strain = aggregate(foci[keep_ind, 17], list(foci$Strain[keep_ind]), mean)
mean_foci_strain = mean_foci_strain[order(mean_foci_strain$x),]
mean_foci_strain$Strain_sorted = seq(1,15)
mean_foci_strain = as.data.frame(mean_foci_strain) # put overall averages into dataframe
colnames(mean_foci_strain) = c('Strain','Overall_mean_foci','Strain_sorted')
mean_foci_strain$Strain_sorted = as.factor(mean_foci_strain$Strain_sorted)
foci = merge(foci,mean_foci_strain, by.x='Strain', by.y='Strain')
keep_ind = foci$hour>8
qplot(Strain_sorted, avg_foci_no_outl, data=foci[keep_ind,], geom=c('boxplot','jitter'),color=Radiation.Type,facets=dose_level~Radiation.Time)

# Create new dataframe where userID is the row and foci level for all three radiation condition
foci_merge = merge(foci[foci$Radiation.Type=="Si 350 MeV/n",],foci[foci$Radiation.Type=="Ar 350 MeV/n",],
                   by.x=c('UserID','hour','dose_level'), 
                   by.y=c('UserID','hour','dose_level'))
foci_merge = merge(foci_merge, foci[foci$Radiation.Type=="Fe 600 MeV/n",],
                   by.x=c('UserID','hour','dose_level'), 
                   by.y=c('UserID','hour','dose_level'))
foci_merge = merge(foci_merge, foci[foci$Radiation.Type=="X-ray",],
                   by.x=c('UserID','hour','dose_level'), 
                   by.y=c('UserID','hour','dose_level'))
qplot(avg_foci_no_outl.y,avg_foci_no_outl,data=foci_merge,color=Strain_sorted,facets=dose_level~.,
      xlab="Foci/cell for Fe", ylab="Foci/cell for Ar")

keep_ind = foci$Strain=='BALBC'&foci$Radiation>0

xyplot(avg_foci_no_outl ~ hour, groups=UserID, grid=TRUE, data=foci[keep_ind,],
       type=c("p","a"))
xyplot(avg_foci_no_outl ~ hour | UserID, grid=TRUE, data=foci[keep_ind,],
       type=c("p","a"))
print(xyplot(avg_foci_no_outl ~ hour | UserID, grid=TRUE, data=foci[keep_ind,],
             xlab = "Hours", ylab = "Foci per cell",
             panel = function(x, y) {
               panel.xyplot(x, y)
               fm <- nls(y ~ a*((1-b)*exp(-c*x) + b), start = list(a = max(y), b = min(y)/max(y), c = 0.3))
               val_fit <- summary(fm)$coefficients
               print(val_fit)
               panel.linejoin(x, predict(fm), col.line = "black")
               panel.text(mean(x),max(y),
                          paste(summary(fm),"\na = ",as.character(round(val_fit[1,1],digits=2)),
                                ", b = ",as.character(round(val_fit[2,1],digits=2)),
                                ", c = ",as.character(round(val_fit[3,1],digits=2))))
             },
             as.table=T))

# Fit time dependence of control with linear model
# Create corrected value for foci # by removing fitted control time dependence for each user ID based onkeep_ind = foci$Radiation>0.11
# require(nlme)
# keep_ind = foci$Radiation==0 & foci$Radiation.Type=='Si 350 MeV/n'
# fm <- nlsList(avg_foci_no_outl ~ a*hour + b | UserID, data=foci[keep_ind,], list(a = 0, b = 0.1))
# val_fit <- summary(fm)$coefficients
# print(val_fit)
# control_fit = as.data.frame(val_fit) # put fit for control into dataframe
# control_fit = cbind(control_fit,rownames(control_fit))
# control_fit = control_fit[,c(1,5,9)]
# colnames(control_fit) = c("ctl_slope","ctl_offset","UserID")
# foci = merge(foci, control_fit, by.x='UserID', by.y='UserID')
# foci$foci_cor = foci$avg_foci_no_outl/(foci$hour*foci$ctl_slope+foci$ctl_offset)

# Create corrected value by dividing by average controls
keep_ind = foci$Radiation==0 
bgd_avg = aggregate(x = foci[keep_ind,]$avg_foci_no_outl, by = list(foci[keep_ind,]$hour, foci[keep_ind,]$Radiation.Type, 
                                                                    foci[keep_ind,]$Strain, foci[keep_ind,]$Strain_sorted), 
                    FUN = "mean")
colnames(bgd_avg) = c("hour","Type","Strain","Strain_sorted", "Ctl_avg")
qplot(Strain_sorted, Ctl_avg, data=bgd_avg, facets=hour~Type)

foci = merge(foci, bgd_avg, by.x=c('Strain','Radiation.Type','hour'), by.y=c('Strain','Type','hour'))
foci$foci_cor = foci$avg_foci_no_outl-foci$Ctl_avg
foci$foci_cor[foci$foci_cor<0]=0 # remove any negative values for corrected foci
for (i_strain in 1:num_strain) {
  keep_ind = (foci$Strain == strain_list[i_strain] & foci$Radiation>0)
  print(qplot(Radiation.Time, foci_cor, data = foci[keep_ind,], 
        geom=c("boxplot"), color=dose_level,facets=.~Radiation.Type,ylim=c(0,4),
        xlab="Time post exposure (hr)",ylab="Foci/cell", main=strain_list[i_strain]))
}
# plot only 24 and 48 hr separately, with all strain sorted on x axis (Antoine display)
time_point = 48
keep_ind = foci$hour == time_point & foci$Radiation>0
qplot(Strain_sorted.x,foci_cor, data=foci[keep_ind,],
      geom=c("boxplot"), color = Radiation.Type, facets=.~dose_level,
      xlab = "Sorted strains", ylab="Foci/cell", main=sprintf('%d hours',time_point))

# plot foci versus LET for both high and low dose for each time point (saturation effect)

# Fit with exponential repair for irradiation high dose
keep_ind = foci$Strain=='BALBC'&foci$Radiation>0.25
fm <- nlsList(avg_foci_no_outl ~ a*((1-b)*exp(-c*hour) + b) | UserID, data=foci[keep_ind,], list(a = 7, b = 0.1, c = 0.2))
val_fit <- summary(fm)$coefficients
print(val_fit)
predict(fm)

