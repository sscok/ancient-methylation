individuals=c("UstIshim","Motala12","Loschbour","K14","R15","R2","R7","Sf12","R15","R7"
,"irk025","irk061","kra001","cta016","VLASA32","VLASA7","prs002","prs013","prs016","prs009","LBK", "R2","R3","R9","AKT16","Asp6","BAR25","Dil16","Ess7","Herx","Klein7","LEPE48","LEPE52","Nea2",
"Nea3","STAR1","VC3-2")

for (i in individuals){
	df=read.delim(i,head=F, sep="\t")
	colnames(df) <- c("chr", "pos", "deaminated","total","gene")
	filtered=df[df$total>9,]
	result=as.data.frame(cbind(filtered, ratio=filtered$deaminated/filtered$total))
	write.table(result, paste("ratios_10/", i, sep=""), quote=F, col.names=T, row.names=F, sep="\t")
}