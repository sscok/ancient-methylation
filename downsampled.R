#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
setwd("/mnt/NAS/projects/2019_epipaleomix/epipaleo/out_shotgun/epipaleo_out/anova_marchi/ratios_10")
names= c("AKT16_replicated","Asp6_replicated","BAR25_replicated","Dil16_replicated","Ess7_replicated","Herx_replicated","Klein7_replicated","LBK_replicated","LEPE48_replicated","LEPE52_replicated","Loschbour_replicated","Nea2_replicated","Nea3_replicated","prs013_replicated","Sf12_replicated","STAR1_replicated","UstIshim_replicated","VC3-2_replicated","VLASA32_replicated","VLASA7_replicated")

#names=paste0("../", names)
for (k in c(2:20)){
for(i in names){
	k14=read.delim(i, head=F, sep="\t")
	n0=sum(k14$V3==0)
	n1=sum(k14$V3==1)
	p=n1/(n0+n1)
	m=round(ifelse(p<0.02, 49*n1, n0/49))
	deam=k14$V3
	names(deam) <- seq_along(deam)
	det=c()
	sampled=c()
	det2=c()
	if (m==49*n1) {
		det=deam[deam==0]
		sampled=sample(det, m)
		det2=c(sampled,deam[deam==1])
	} else {
		det=deam[deam==1]
		sampled=sample(det, m)
		det2=c(sampled,deam[deam==0])
	}
	index=as.integer(names(det2))
	pos=k14[index,]
	sum(sampled==1)
	p1=sum(det2==1)/(sum(det2==1)+sum(det2==0))
	print(p1)
	print(dim(pos))
	write.table(pos,paste(paste("downsampled",i,sep="_"), k, sep="_"), quote=F, col.names=F, row.names=F, sep="\t")
}
}