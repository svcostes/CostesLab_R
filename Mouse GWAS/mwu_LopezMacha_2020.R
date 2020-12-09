library(gridExtra)
library(ggbio)
library(biovizBase)
library(GenomicRanges)
library(rtracklayer)
library(snow)

# mw() - A function used to parallelize the mann-whitney and filter SNP data.
# "snp" - a table of SNPs at a given position (row) for all strains (columns)
# "phenotab" - a table of phenotypes with columns: UserID	let	hour Strain	Bgd	FociPerGy
# "thresh" - threshold for minor allele frequency. SNPs fewer than "thresh" representation 
#          in the sample of strains (not that each strain only counts once
#         , regardless of the number of animals per strain).

mw <- function(snp, phenotab, thresh) {
  p<-NA
  noNHalleles<-snp[(!(grepl("N",as.matrix(snp)))&!(grepl("H",as.matrix(snp))))]	# remove the ambiguous genotypes
  if (length(noNHalleles) > 2) {													# only use biallelic sites
    noNHalleles<-noNHalleles[3:length(noNHalleles)]							
    cleanpheno<-phenotab[which(names(phenotab)%in%names(noNHalleles))]			# subset phenotypes for only the strains that are reliably genotyped
    names(cleanpheno) <- names(phenotab)[which(names(phenotab)%in%names(noNHalleles))]
    alleles<-unique(unlist(noNHalleles[3:length(noNHalleles)]))					# get the two alleles for the SNP
    a1list<-names(snp[grepl(alleles[1],as.matrix(snp))])						# make a list of strains with allele 1
    a2list<-names(snp[grepl(alleles[2],as.matrix(snp))])						# make a list of strains with allele 2
    pheno1<-cleanpheno[,which(names(cleanpheno) %in% a1list)]					# subset phenotypes by alleles
    pheno2<-cleanpheno[,which(names(cleanpheno) %in% a2list)]
    bystrain<-noNHalleles[match(unique(names(cleanpheno)),names(noNHalleles))]	# check the minor allele frequency
    maf<-min(table(factor(unlist(bystrain))))
    if (maf>=thresh & length(pheno1)>=thresh & length(pheno2)>=thresh) {		# make sure everything passes thresh (last two tests probably unnecessary!)
      mwu<-suppressWarnings(wilcox.test(as.matrix(pheno1),as.matrix(pheno2))) # apply mann-whitney U 
      p<-mwu$p.value    
    }
  }
  return(p)																		# return p-value of the mwu test
}

# We use clusters to parallelize the mwu test corresponding to the different SNP locations
cl <- makeCluster(12)
clusterExport(cl,"mw")
t<-3 					# The minor allele frequency threshold is set to 3.
clusterExport(cl,"t")

# This next section reads in the SNP data and phenotype table, subsets the phenotype data
# and then runs the mwu function on every SNP for every phenotype. Results are stored 
# in the table "res" with each phenotype as a column. Chromosome and position are columns 1-2. 
# Note that because of filtering this will be smaller (fewer rows) than the full SNP table.

# Read SNP table
setwd('Sascha_Langley_analysis')
snps<-read.table("Sylvain_merged_GENO.txt", header=TRUE) 
snps<-snps[,2:ncol(snps)]

# Read phenotype table
pheno<-read.csv("phenotypes_by_UserID.csv",header=TRUE)

# Prepare result table for all Bgd and FociPerGy phenotypes
plist<-c("Bgd","FociPerGy")
letlist<-unique(pheno$let)
timelist<-unique(pheno$hour)
res<-data.frame(snps[,1:2]) # Variable tracking pvalue for each phenotype vs snp position
res_r<-data.frame(snps[,1:2]) # variable tracking pvalue for randomized phenotype, to determine maximum random p-value
cnt = 1 # indexation to keep track of the increment of colnames

for (i in 1:length(letlist)) { # Screen all LET conditions
  for (j in 1:length(timelist)) { # Screen all time points
    for (n in 1:length(plist)) { # Both Bgd and FociPerGy phenotypes
      subpheno<-pheno[pheno$let==letlist[i] & pheno$hour==timelist[j],plist[n]] # Select phenotypes that correspond to the LET and time point parameters
      finalpheno<-as.data.frame(t(subpheno))
      colnames(finalpheno)<-pheno[pheno$let==letlist[i] & pheno$hour==timelist[j],"Strain"]
      p<-finalpheno
      clusterExport(cl,"p")
      m<-parRapply(cl, snps, function(x) mw(x,p,t)) # Apply mw function in parallel for the different SNP locations (with the selected phenotypes)
      res<-cbind(res,m) # Create the table of p values (rows are SNP locations, columns are phenotypes)
      colnames(res)[2+cnt]<-paste(plist[n],paste(letlist[i],timelist[j],sep="_"),sep="_") # Each column of res is named by LET and time point conditions
      
      # Now we run random associations (4 times) to compare resulting p values to non-randomized ones
      for (i_r in 1:4) {
        finalpheno_r <- as.data.frame(t(sample(subpheno))) # Randomize order of strains with sample
        colnames(finalpheno_r)<-pheno[pheno$let==letlist[i] & pheno$hour==timelist[j],"Strain"]
        p<-finalpheno_r
        clusterExport(cl,"p")
        m<-parRapply(cl, snps, function(x) mw(x,p,t))
        res_r<-cbind(res_r,m)
        colnames(res_r)[2+cnt]<-paste(plist[n],paste(letlist[i],timelist[j],sep="_"),sep="_") #SVC
      }
      cnt = cnt + 1
    }
  }
}

# Determine maximum p-value for each random phenotype
mm10 <- getIdeogram(genome = "mm10",cytoband=TRUE)
p_r = array(NA,ncol(res_r)-2) # Variable tracking the maximum p-value for each randomized phenotypes
for (n in 3:ncol(res_r)) {
  sub<-res_r[,c(1,2,n)]
  sub<-sub[complete.cases(sub),]
  sub<-sub[sub[,1]!="M",]
  if (nrow(sub)>0) {
    sub[,3]<--log10(sub[,3])
    chr<-unique(sub[,1])
    name<-colnames(res_r)[n]
    title<-paste(name,"_t3",sep="")
    filename<-paste(name,"_manhattan.jpg",sep="")
    gr<-c()
    for (i in 1:length(chr)) {
      temp=sub[sub[,1]==chr[i],]
      if (nrow(temp)>0) {
        ir<-IRanges(start=temp[,2],end=temp[,2])
        tgr<-GRanges(seqnames = paste("chr",chr[i],sep=""), seqlengths = seqlengths(mm10)[paste("chr",chr[i],sep="")], ranges = ir, score = temp[,3])
        gr<-c(gr,tgr)
      }
    }
    ul<-unlist(gr)
    finalgr<-suppressWarnings(c(ul[[1]],ul[[2]],ul[[3]],ul[[4]],ul[[5]],ul[[6]],ul[[7]],ul[[8]],ul[[9]],ul[[10]],ul[[11]],ul[[12]],ul[[13]],ul[[14]],ul[[15]],ul[[16]],ul[[17]],ul[[18]],ul[[19]],ul[[20]]))
    p_r[n-2] = max(finalgr$score)    
  }
}

th_man = max(p_r,na.rm=TRUE)*1.2 # set cutoff for manhattan significance with additional 10%

# Plot Manhattan graphs
mm10 <- getIdeogram(genome = "mm10",cytoband=TRUE)
maxX<-max(-log10(unlist(res[,3:ncol(res)])[which(!is.na(unlist(res[,3:ncol(res)])))]))
maxX<-ceiling(maxX)
mycol<-c("#F15859","#4E5CA5","#C76BAD","#FAADB0","#828282","#992819","#5505D3","#29B13F","#C6C23C",
         "#16BFC1","#A826A2","#034C93","#E28E1E","#3D8964","#E05E3B","#599AE2","#A8B3C1","#AA485B",
         "#8686EA","#CEA806")

for (n in 3:ncol(res)) {
  sub<-res[,c(1,2,n)]
  sub<-sub[complete.cases(sub),]
  sub<-sub[sub[,1]!="M",]
  if (nrow(sub)>0) {
    sub[,3]<--log10(sub[,3])
    chr<-unique(sub[,1])
    #	maxX<-max(sub[,3])
    name<-colnames(res)[n]
    title<-paste(name,"_t3",sep="")
    filename<-paste(name,"_manhattan.jpg",sep="")
    listname<-paste(name,"_significant_snps.csv",sep="")
    gr<-c()
    for (i in 1:length(chr)) {
      temp=sub[sub[,1]==chr[i],]
      if (nrow(temp)>0) {
        ir<-IRanges(start=temp[,2],end=temp[,2])
        tgr<-GRanges(seqnames = paste("chr",chr[i],sep=""), seqlengths = seqlengths(mm10)[paste("chr",chr[i],sep="")], ranges = ir, score = temp[,3])
        gr<-c(gr,tgr)
      }
    }
    ul<-unlist(gr)
    finalgr<-suppressWarnings(c(ul[[1]],ul[[2]],ul[[3]],ul[[4]],ul[[5]],ul[[6]],ul[[7]],ul[[8]],ul[[9]],ul[[10]],ul[[11]],ul[[12]],ul[[13]],ul[[14]],ul[[15]],ul[[16]],ul[[17]],ul[[18]],ul[[19]],ul[[20]]))
    significant_finalgr = as.data.frame(finalgr[finalgr$score>th_man])
    if (dim(significant_finalgr)[1]>0) {
      write.table(significant_finalgr, file=listname, row.names=F, col.names=T, sep=",")
      grp<-plotGrandLinear(finalgr, aes(y = score), color = mycol) + labs(title = colnames(res)[n], y = "-log10(p-val)") + scale_y_continuous(limits = c(0, maxX)) + theme_classic() + theme(legend.position="none") + ggtitle(title)
      grp <- grp + geom_hline(yintercept=th_man, linetype="dashed", 
                              color = "red", size=1)
      jpeg(filename, width = 1200, height = 800)
      print(grp)
      dev.off()
      ggsave(filename,grp@ggplot, width = 15, height = 6)
    }
  }
  
}



