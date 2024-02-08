library(ggplot2)
library(scales)
library(gridExtra)
args = commandArgs(trailingOnly=TRUE)

# args[1]: set working directory
# args[2]: sorted MS dataset input name

setwd(args[1])

                                                                # read MS dataset
input=read.table("matrix_sorted_2023",head=T, sep="\t")
inds=c("chr", "start","end", "AKT16","Asp6","BAR25","cta016","Dil16","Ess7","Herx","irk025","irk061","K14","Klein7","kra001","LBK","LEPE48","LEPE52","Loschbour",
  "Motala12","Nea2","Nea3","prs002","prs009","prs013","prs016","R15","R2","R3","R7","R9","Sf12","STAR1","Ust-Ishim","VC3-2","VLASA32","VLASA7")


hg=c("Ust-Ishim", "Motala12","Loschbour", "K14","Sf12","R15","R7", "irk025","irk061","kra001","cta016","VLASA32","VLASA7")

ne=c("prs002","prs016","prs013","prs009","LBK","R2","R3","R9","AKT16","Asp6","BAR25","Dil16","Ess7","Herx","Klein7","LEPE48","LEPE52","Nea2",
     "Nea3","STAR1","VC3-2")
a = c(hg,ne)
colnames(input) <- inds

inputhg <- input[,hg]
inputne <- input[,ne]

inputhg = inputhg[,order(colMeans(inputhg, na.rm=T))]
inputne = inputne[,order(colMeans(inputne, na.rm=T))]

inputall <- cbind(input[,1:3], inputhg,inputne)
individuals=rep(colnames(inputall)[4:ncol(inputall)], each=length(inputall$chr))

values=unlist(inputall[,4:37],use.names = FALSE)
processed=data.frame(individual=individuals,MethylationScore=values+1)


processed$individual <- factor(processed$individual, levels = unique(processed$individual))
p1=ggplot(processed, aes(x=individual, y=MethylationScore)) +
          geom_violin(color="#33638DFF", fill="steelblue",width=1, size=0.35) + scale_y_continuous(trans = 'log2',breaks = trans_breaks("log2", function(x) 2^x),
                                                                                                          labels =trans_format("log2",math_format(.x))) +
  stat_summary(fun="median", geom="point", color="black", size=1.5) + stat_summary(fun="mean", geom="point",color="#440154FF", size=1.5) + theme_bw() + labs(y="Methylation Score",size=25) +
    theme(axis.text=element_text(size=15), legend.position="none", panel.background = element_blank(),axis.title.x=element_blank(), axis.text.y = element_text(vjust = 0.5, hjust=1, color = c(rep("blue",13), rep("red", 21)))) + coord_flip()


p2=ggplot(processed, aes(x=individual, y=MethylationScore)) +
          geom_violin(color="#33638DFF", fill="steelblue",width=0.75, size=0.25) + scale_y_continuous(trans = 'log2',breaks = trans_breaks("log2", function(x) 2^x),
                                                                                                          labels =trans_format("log2",math_format(.x))) +
  stat_summary(fun="median", geom="point", color="black", size=1.5) + stat_summary(fun="mean", geom="point",color="#440154FF", size=1.5) + theme_bw() + labs(y="Methylation Score",size=25) +
    theme(axis.text=element_text(size=15),legend.position="none", panel.background = element_blank(),axis.title.x=element_blank(), axis.text.y = element_text(vjust = 0.5, hjust=1, color = c(rep("blue",13), rep("red", 21)))) + coord_flip(ylim=c(2^0, 2^0.07))
pdf("vioplots_marchi.pdf",width="12",height="20")
grid.arrange(p1,p2,ncol=2)
dev.off()
