mycol<-c("#F15859",
"#4E5CA5",
"#C76BAD",
"#FAADB0",
"#828282",
"#992819",
"#5505D3",
"#29B13F",
"#C6C23C",
"#16BFC1",
"#A826A2",
"#034C93",
"#E28E1E",
"#3D8964",
"#E05E3B",
"#599AE2",
"#A8B3C1",
"#AA485B",
"#8686EA",
"#CEA806")

library(gridExtra)
library(ggbio)
library(biovizBase)
library(GenomicRanges)
library(rtracklayer)
mm10 <- getIdeogram(genome = "mm10",cytoband=TRUE)


res<-read.table("Sylvain_phenotypes_by_UserID_mwu_t3.txt",header=TRUE)

maxX<-max(-log10(unlist(res[,3:ncol(res)])[which(!is.na(unlist(res[,3:ncol(res)])))]))
maxX<-ceiling(maxX)

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
	filename<-paste(name,"_manhattan.pdf",sep="")
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
	grp<-plotGrandLinear(finalgr, aes(y = score), color = mycol) + labs(title = colnames(res)[n], y = "-log10(p-val)") + scale_y_continuous(limits = c(0, maxX)) + theme_classic() + theme(legend.position="none") + ggtitle(title)
	ggsave(filename,grp@ggplot, width = 15, height = 6)
	}
}


