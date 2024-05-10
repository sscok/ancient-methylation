library(cowplot)
library(grid)
library(gridExtra)
library(ggplot2)

# Description: Below scripts are used to plot the Figure 2 in the article. 
# Representative genes with the most significant differential methylation signals in linear mixed model analyses,
# related to subsistence type, tissue and sex.
# inputs: 
#	    - table: with column names corresponding to the individuals, 
#           row names representing the genomic position,
#           corresponding data being the Methylation score.
# output:
#	    - a merged scatterplot for all the candidate genes with all the categories described.

#------------------- subsistence type, laboratory-of-origin as random factor ------------------------------ 
args = commandArgs(trailingOnly=TRUE)
data=read.table(paste0(args[1],"TOX2.lab"), head=F)
result <- aggregate(V3 ~ V2 + V4 + V5, data, FUN=mean, na.rm=TRUE)
resultvalue <- aggregate(V3 ~ V4, result, FUN=mean, na.rm=TRUE)
hg=c("UstIshim","Loschbour","Sf12","R7", "irk025","irk061","kra001","VLASA32","VLASA7")

ne=c("prs002","prs009","LBK","R2","R3","R9","AKT16","Asp6","BAR25","Dil16","Ess7","Herx","Klein7","LEPE48","LEPE52","Nea2",
     "Nea3","STAR1","VC3-2")
all= c(hg,ne)
df = resultvalue
df <- df[order(match(df$V4, all)), ]
st = c(rep("HG", length(hg)), rep("NF", length(ne)) )

df_all = cbind(df, st=st)
colnames(df_all) = c("Individual", "Value", "Category")

scatterplot1 <- ggplot(df_all, aes(x = Category, y = log2(Value+1), color = Category)) +
  ggtitle("TOX2") +
  geom_jitter(position=position_jitter(0.2), size=2.5) +
  labs(x = "Category", y = "Mean Methylation Score Per Individual", tag="A") +
  scale_color_manual(values = c("HG" = "blue", "NF" = "blue"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.0001))+
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 15),axis.title.y = element_text(size=10),axis.text=element_text(size=15),plot.tag = element_text(size = 23, face = "bold"),legend.position="none",panel.background = element_blank(), axis.title.x=element_blank(), axis.text.y = element_text(vjust = 0.5, hjust=1)) +
  stat_summary(fun = median, geom = "crossbar", size = 0.5, width = 0.5, color = "black") +  
  guides(color = FALSE)

#------------------- subsistence type, individual as random factor ------------------------------ 
data=read.table(paste0(args[1],"ICAM5.indf"), head=F)
result <- aggregate(V3 ~ V2 + V4 + V5, data, FUN=mean, na.rm=TRUE)
resultvalue <- aggregate(V3 ~ V4, result, FUN=mean, na.rm=TRUE)
hg=c("UstIshim","Loschbour","Sf12","R7","R15", "cta016", "irk025","irk061","kra001","VLASA32","VLASA7")

ne=c("prs002","prs009","LBK","R2","R3","R9","AKT16","Asp6","BAR25","Dil16","Ess7","Herx","Klein7","LEPE48","LEPE52","Nea2",
     "Nea3","STAR1","VC3-2")
all= c(hg,ne)
df = resultvalue
df <- df[order(match(df$V4, all)), ]
st = c(rep("HG", length(hg)), rep("NF", length(ne)) )

df_all = cbind(df, st=st)
colnames(df_all) = c("Individual", "Value", "Category")

scatterplot2 <- ggplot(df_all, aes(x = Category, y = log2(Value+1), color = Category)) +
  ggtitle("ICAM5") +
  geom_jitter(position=position_jitter(0.2), size=2.5) +
  labs(x = "Category", y = "Mean Methylation Score Per Individual", tag="B") +
  scale_color_manual(values = c("HG" = "blue", "NF" = "blue"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.0001))+
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 15), axis.title.y = element_text(size=10),axis.text=element_text(size=15),plot.tag = element_text(size = 23, face = "bold"),legend.position="none",panel.background = element_blank(), axis.title.x=element_blank(), axis.text.y = element_text(vjust = 0.5, hjust=1)) +
  stat_summary(fun = median, geom = "crossbar", size = 0.5, width = 0.5, color = "black") +  
  guides(color = FALSE) 

#------------------- tissue type, laboratory-of-origin as random factor ------------------------------ 
data=read.table(paste0(args[1],"PCDHA2.lab"), head=F)
result <- aggregate(V3 ~ V2 + V4 + V5, data, FUN=mean, na.rm=TRUE)
resultvalue <- aggregate(V3 ~ V4, result, FUN=mean, na.rm=TRUE)
tooth=c("Loschbour" , "irk025","irk061","kra001","cta016","prs002","prs016","prs013","prs009","LBK")

bone=c("UstIshim","Sf12", "R15","VLASA32","VLASA7","R2","R3","R9","AKT16","Asp6","BAR25","Dil16","Ess7","Herx","Klein7","LEPE48","LEPE52","Nea2",
     "Nea3","STAR1","VC3-2")
all= c(bone,tooth)
df = resultvalue
df <- df[order(match(df$V4, all)), ]
st = c(rep("Bone", length(bone)), rep("Tooth", length(tooth)))

df_all = cbind(df, st=st)
colnames(df_all) = c("Individual", "Value", "Category")

scatterplot3 <- ggplot(df_all, aes(x = Category, y = log2(Value+1), color = Category)) +
  ggtitle("PCDHA2") +
  geom_jitter(position=position_jitter(0.2), size=2.5) +
  labs(x = "Category", y = "", tag="A") +
  scale_color_manual(values = c("Bone" = "red", "Tooth" = "red"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.0001))+
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 15), axis.title.y = element_text(size=10),axis.text=element_text(size=15),plot.tag = element_text(size = 23, face = "bold", color = "white"),legend.position="none",panel.background = element_blank(), axis.title.x=element_blank(), axis.text.y = element_text(vjust = 0.5, hjust=1)) +
  stat_summary(fun = median, geom = "crossbar", size = 0.5, width = 0.5, color = "black") +  
  guides(color = FALSE)

#------------------- tissue type, individual as random factor  ------------------------------ 
data=read.table(paste0(args[1],"ATP1B3.indf"), head=F)
result <- aggregate(V3 ~ V2 + V4 + V5, data, FUN=mean, na.rm=TRUE)
resultvalue <- aggregate(V3 ~ V4, result, FUN=mean, na.rm=TRUE)
tooth=c("Loschbour" , "irk025","irk061","kra001","cta016","prs002","prs016","prs013","prs009","LBK")

bone=c("UstIshim","Sf12", "R15","VLASA32","VLASA7","R2","R3","R9", "R15", "AKT16","Asp6","BAR25","Dil16","Ess7","Herx","Klein7","LEPE48","LEPE52","Nea2",
     "Nea3","STAR1","VC3-2")
all= c(bone,tooth)
df = resultvalue
df <- df[order(match(df$V4, all)), ]
st = c(rep("Bone", length(bone)), rep("Tooth", length(tooth)))

df_all = cbind(df, st=st)
colnames(df_all) = c("Individual", "Value", "Category")

scatterplot4 <- ggplot(df_all, aes(x = Category, y = log2(Value+1), color = Category)) +
  ggtitle("ATP1B3") +
  geom_jitter(position=position_jitter(0.2), size=2.5) +
  labs(x = "Category", y = "", tag="B") +
  scale_color_manual(values = c("Bone" = "red", "Tooth" = "red"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.0001))+
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 15), axis.title.y = element_text(size=10),axis.text=element_text(size=15),plot.tag = element_text(size = 23, face = "bold", color = "white"),legend.position="none",panel.background = element_blank(), axis.title.x=element_blank(), axis.text.y = element_text(vjust = 0.5, hjust=1)) +
  stat_summary(fun = median, geom = "crossbar", size = 0.5, width = 0.5, color = "black") +  
  guides(color = FALSE)

#------------------- genetic sex, laboratory-of-origin as random factor ------------------------------ 
data=read.table(paste0(args[1],"RCOR1.lab"), head=F)
result <- aggregate(V3 ~ V2 + V4 + V5, data, FUN=mean, na.rm=TRUE)
resultvalue <- aggregate(V3 ~ V4, result, FUN=mean, na.rm=TRUE)
xy=c("VC3-2","Asp6","LEPE48","LEPE52","BAR25","Dil16","Ess7","UstIshim","Motala12","Loschbour" , "K14", "R9","R15","R7","irk025","irk061","kra001","prs009","VLASA32","VLASA7")

xx=c("Sf12","cta016","prs002","LBK","R2","R3","AKT16","Herx","Klein7","Nea2",
     "Nea3","STAR1")
all= c(xx,xy)
df = resultvalue
df <- df[order(match(df$V4, all)), ]
st = c(rep("XX", length(xx)), rep("XY", length(xy)))

df_all = cbind(df, st=st)
colnames(df_all) = c("Individual", "Value", "Category")

scatterplot5 <- ggplot(df_all, aes(x = Category, y = log2(Value+1), color = Category)) +
  ggtitle("RCOR1") +
  geom_jitter(position=position_jitter(0.2), size=2.5) +
  labs(x = "Category", y = "", tag="A") +
  scale_color_manual(values = c("XX" = "purple", "XY" = "purple"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.0001))+
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 15), axis.title.y = element_text(size=10),axis.text=element_text(size=15),plot.tag = element_text(size = 23, face = "bold", color = "white"),legend.position="none",panel.background = element_blank(), axis.title.x=element_blank(), axis.text.y = element_text(vjust = 0.5, hjust=1)) +
  stat_summary(fun = median, geom = "crossbar", size = 0.5, width = 0.5, color = "black") +  
  guides(color = FALSE)
#------------------- genetic sex, individual as random factor ------------------------------ 
data=read.table(paste0(args[1],"CEP135.indf"), head=F)
result <- aggregate(V3 ~ V1 + V2 + V4 + V6, data, FUN=mean, na.rm=TRUE)
resultvalue <- aggregate(V3 ~ V4, result, FUN=mean, na.rm=TRUE)
xy=c("VC3-2","Asp6","LEPE48","LEPE52","BAR25","Dil16","Ess7","UstIshim","Loschbour", "R9","R15","irk061","kra001","prs009","prs013","VLASA32","VLASA7")

xx=c("Sf12","cta016","prs002","LBK","R3","AKT16","Herx","Klein7","Nea2",
     "Nea3","STAR1")
all= c(xx,xy)
df = resultvalue
df <- df[order(match(df$V4, all)), ]
st = c(rep("XX", length(xx)), rep("XY", length(xy)))

df_all = cbind(df, st=st)
colnames(df_all) = c("Individual", "Value", "Category")

scatterplot6 <- ggplot(df_all, aes(x = Category, y = log2(Value+1), color = Category)) +
  ggtitle("CEP135") +
  geom_jitter(position=position_jitter(0.2), size=2.5) +
  labs(x = "Category", y = "", tag="B") +
  scale_color_manual(values = c("XX" = "purple", "XY" = "purple"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.0001))+
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 15), axis.title.y = element_text(size=10),axis.text=element_text(size=15),plot.tag = element_text(size = 23, face = "bold", color = "white"),legend.position="none",panel.background = element_blank(), axis.title.x=element_blank(), axis.text.y = element_text(vjust = 0.5, hjust=1)) +
  stat_summary(fun = median, geom = "crossbar", size = 0.5, width = 0.5, color = "black") +  
  guides(color = FALSE)


# merge all the scatter plots and create the Figure 2 in the article.
pdf("plots/combined_scatterplots.pdf", width = 12, height = 6)
grid.arrange(scatterplot1, scatterplot3, scatterplot5, scatterplot2, scatterplot4, scatterplot6, ncol = 3)
dev.off()