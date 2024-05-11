# Description: This script generates the datasets with the Methylation Score ratios per individual per CpG position.
# inputs: 
#	- outputs from the epiPALEOMIX results, individual file names
# output:
#	- processed output files per individual

individuals=c("UstIshim","Motala12","Loschbour","K14","R15","R2","R7","Sf12","R15","R7"
,"irk025","irk061","kra001","cta016","VLASA32","VLASA7","prs002","prs013","prs016","prs009","LBK", "R2","R3","R9","AKT16","Asp6","BAR25","Dil16","Ess7","Herx","Klein7","LEPE48","LEPE52","Nea2",
"Nea3","STAR1","VC3-2")

for (i in individuals){
	df=read.delim(paste0("data/",i),head=F, sep="\t")
	colnames(df) <- c("chr", "pos", "deaminated","total","bedcoord")
	filtered=df[df$total>4,]	# filter the number of reads
	result=as.data.frame(cbind(filtered, ratio=filtered$deaminated/filtered$total))
	write.table(result, paste("outputs/", i, sep=""), quote=F, col.names=T, row.names=F, sep="\t")
}