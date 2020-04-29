# mw() - A function used to parallelize the mann-whitney and filter SNP data.
# "snp" - a tab delineated table of SNPs for all strains with strain ID in the header (columns)
#	FORMAT: marker	chromosome	position(b38)	StrainA 	StrainA	StrainB	StrainB	etc.
# "phenotab" - a table of phenotypes
# 	(your) FORMAT: UserID	let	hour	Strain	intercept	Bgd	FociPerGy
# "thresh" - threshold for minor allele frequency. SNPs fewer than "thresh" respesentation 
#          in the sample of strains (not sample of samples! Each strain only counts once
#          here, regardless of the number of that strain included in the sample. This is 
#          conservative and could be changed to overall minor allele frequency in the sample,
#          but I don't recommend that.).

mw <- function(snp, phenotab, thresh) {
	p<-NA
    noNHalleles<-snp[(!(grepl("N",as.matrix(snp)))&!(grepl("H",as.matrix(snp))))]	# remove the ambiguous genotypes
    if (length(noNHalleles) > 2) {													# only use biallelic sites
        noNHalleles<-noNHalleles[3:length(noNHalleles)]							
        cleanpheno<-phenotab[which(names(phenotab)%in%names(noNHalleles))]			# subset phenotypes for only those that are reliably genotyped
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
	return(p)																		# return p-value
}

library(snow)
cl <- makeCluster(12)
clusterExport(cl,"mw")
t<-3 					# Set minor allele frequency threshold. I usually use 5 or higher, 
						# but we don't have many strains here.
clusterExport(cl,"t")


# This next section reads in your SNP data and phenotype table, subsets the phenotype data
# and then runs the above function on every SNP for every phenotype. Results are stored 
# in the table "res" with each phenotype as a column. Chromosome and position are columns 1-2. 
# Note that because of filtering this will be smaller (fewer rows) than the full SNP table ("Sylvain_merged_GENO.txt").

setwd('C:/Users/scostes/Documents/Mouse GWAS/Sascha_Langley_analysis')
snps<-read.table("Sylvain_merged_GENO.txt", header=TRUE) 
snps<-snps[,2:ncol(snps)]
plist<-c("intercept","Bgd","FociPerGy")
pheno<-read.table("phenotypes_by_UserID_SL.txt",header=TRUE)
letlist<-unique(pheno$let)
timelist<-unique(pheno$hour)
res<-data.frame(snps[,1:2])
for (i in 1:length(letlist)) {
	for (j in 1:length(timelist)) {
		for (n in 1:length(plist)) {
			subpheno<-pheno[pheno$let==letlist[i] & pheno$hour==timelist[j],plist[n]]
			finalpheno<-as.data.frame(t(subpheno))
			colnames(finalpheno)<-pheno[pheno$let==letlist[i] & pheno$hour==timelist[j],"Strain"]
			p<-finalpheno
			clusterExport(cl,"p")
			m<-parRapply(cl, snps, function(x) mw(x,p,t))
			res<-cbind(res,m)
			colnames(res)[2+i]<-paste(plist[n],paste(letlist[i],timelist[j],sep="_"),sep="_")
		}
	}
}
