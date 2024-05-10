# methylation dataset downsampled (normalized) 20x, and mean methylation score (MS) per gene calculated, and means of 20x replicates calculated (NAs not removed)

a = read.table("reshaped_ds_all_means_v2", row.names = 1, head = T)  # input dataset which has been created by getting the means of the 20 downsampled datasets
colnames(a)
head(a)
dim(a)

# info table

b = read.table("Shotgun_inds.tsv", row.names = 1, head = T, fill=T, sep="\t")
head(b)
setdiff(colnames(a), rownames(b)) 
setdiff(rownames(b), colnames(a))
dim(b)

b <- b[colnames(a), ]
rownames(b) == colnames(a)

levels(factor(b$Article))
b$Article <- factor(b$Article, levels=c("Lazaridis2014", "Marchi2022", "Fu2014", "Gunther2018", "Quinto2019"))

table(b$lifestyle)
table(b$tissue)
table(b$genetic_sex)
table(b$Article)

table(b$Article, b$tissue)


table(b$genetic_sex, b$tissue, b$lifestyle)

### Multidimensinal Scaling Code

par(mfrow=c(1,1))
mds <- cmdscale(dist(t(a)))
plot(mds, type="n")
text(mds, labels = rownames(mds), col=as.numeric(factor(b$lifestyle)))

# better not to use Motala

noMotala <- !(rownames(b) == "Motala12" | rownames(b) == "K14")

ax <- a[,noMotala]
bx <- b[noMotala,]

mds <- cmdscale(dist(t(ax)))

par(mfrow=c(1,1))
plot(mds, type="n")
text(mds, labels = rownames(mds), col=as.numeric(factor(bx$lifestyle)))

plot(mds, type="n")
text(mds, labels = rownames(mds), col=as.numeric(factor(bx$Article)))

plot(mds, type="n", xlab="Dimension 1", ylab="Dimension 2")
points(mds, col=as.numeric(factor(bx$Article)), pch=as.numeric(factor(bx$lifestyle))+15, cex=2)
legend(-2, 6, legend = c("HG", "NF"), pch = 1:2, bty="n")
legend(4, 6, legend = levels(bx$Article), pch = 3, bty="n", col=as.numeric(factor(levels(bx$Article))))


### analysing coverage
cor.test(bx$cov, as.numeric(as.factor(bx$lifestyle)), met="s")
boxplot(bx$cov ~ bx$lifestyle)

kruskal.test(bx$cov, as.factor(bx$Article), met="s")
summary(aov(bx$cov ~ bx$Article + bx$lifestyle))

par(mfrow=c(2,1))
boxplot(bx$cov ~ bx$lifestyle + bx$Article, las=2, xlab="", ylab="Coverage", col=rep(1:8,each=2), border=rep(1:8,each=2), log="y")

par(mfrow=c(2,1))
x <- boxplot(bx$cov ~ bx$Article + bx$lifestyle, las=2, xlab="", ylab="Coverage", col=1:8, border=1:8, log="y", ylim=c(0.1,200))
text(1:16, x$stats[5,]*3, x$n)

### analysing sum of detected genes
sums <- colSums(!is.na(ax))
names(sums) <- colnames(ax)
sums <- sums[rownames(bx)]
cbind(names(sums), rownames(bx))
bx <- cbind(bx, sums)
head(bx)

kruskal.test(bx$sums, as.factor(bx$lifestyle), met="s")
kruskal.test(bx$sums, as.factor(bx$Article), met="s")
summary(aov(bx$sums ~ bx$Article + bx$lifestyle))

par(mfrow=c(2,1))
x <- boxplot(bx$sums ~ bx$Article + bx$lifestyle, las=2, xlab="", ylab="No. detected genes", col=1:8, border=1:8, log="y", ylim=c(50,50000))
text(1:16, x$stats[5,]*3, x$n)

par(mfrow=c(1,1))
hist(rowSums(!is.na(ax)))
sum(rowSums(!is.na(ax)) >= 12) # 4717   


### analysis of variance

ax2 = ax[rowSums(!is.na(ax)) >= 20,]

dim(ax2)
ax2_aovp <- apply(ax2, 1, function(y) {
  tryCatch(summary(aov(as.numeric(y) ~ factor(bx$lifestyle)))[[1]]$Pr[1], error = function(e) "NaN") } )
hist(ax2_aovp, br=100)

ax2_aovp2 <- t(apply(ax2, 1, function(y) {
  tryCatch(summary(aov(as.numeric(y) ~ factor(bx$lifestyle) + factor(bx$tissue) + factor(bx$genetic_sex)))[[1]]$Pr,error = function(e) NA) }))

par(mfrow=c(1,3))
hist(ax2_aovp2[,1], br=20, main="Subsistence Type", xlab="ANOVA p-value", col="tomato")
hist(ax2_aovp2[,3], br=20, main="Tissue Type", xlab="ANOVA p-value", col="tomato")
hist(ax2_aovp2[,2], br=20, main="Genetic sex", xlab="ANOVA p-value", col="tomato")
# same when using 25 or 30 as cutoff

### removing marchi et al

noMarchi <- !(bx$Article == "Marchietal2022")
dim(ax[,noMarchi])
ax3 = ax[rowSums(!is.na(ax[,noMarchi])) >= 12, noMarchi]
dim(ax3)
ax3_aovp <- apply(ax3, 1, function(y) {
  summary(aov(as.numeric(y) ~ factor(bx$lifestyle[noMarchi])))[[1]]$Pr[1] })
hist(ax3_aovp, br=100)

ax3_aovp2 <- t(apply(ax3, 1, function(y) {
  summary(aov(as.numeric(y) ~ factor(bx$lifestyle[noMarchi]) + factor(bx$genetic_sex[noMarchi]) + factor(bx$tissue[noMarchi]) ))[[1]]$Pr }))

par(mfrow=c(1,3))
hist(ax3_aovp2[,1], br=20, main="Subsistence Type", xlab="ANOVA p-value", col="tomato")
hist(ax3_aovp2[,3], br=20, main="Tissue Type", xlab="ANOVA p-value", col="tomato")
hist(ax3_aovp2[,2], br=20, main="Genetic sex", xlab="ANOVA p-value", col="tomato")

cbind(bx$Article, bx$lifestyle)


### calculating means btw HG vs NF

table(bx$Article, bx$lifestyle)

ax_anto <- t(apply(
  ax[, bx$Article == "AntonioGaoMootsScience2019",], 1, function(y) 
    c(
    mean(y[ bx[bx$Article == "AntonioGaoMootsScience2019",]$lifestyle == "HG" ], na.rm = T), 
    mean(y[ bx[bx$Article == "AntonioGaoMootsScience2019",]$lifestyle == "NF" ], na.rm = T) ) ))

ax_marc <- t(apply(
  ax[, bx$Article == "Marchietal2022",], 1, function(y) 
    c(
      mean(y[ bx[bx$Article == "Marchietal2022",]$lifestyle == "HG" ], na.rm = T), 
      mean(y[ bx[bx$Article == "Marchietal2022",]$lifestyle == "NF" ], na.rm = T) ) ))

ax_laza <- t(apply(
  ax[, bx$Article == "Lazaridis2014",], 1, function(y) 
    c(
      mean(y[ bx[bx$Article == "Lazaridis2014",]$lifestyle == "HG" ], na.rm = T), 
      mean(y[ bx[bx$Article == "Lazaridis2014",]$lifestyle == "NF" ], na.rm = T) ) ))

XXX = ! (bx$Article %in% c("Lazaridis2014", "Marchietal2022", "AntonioGaoMootsScience2019"))
sum(XXX)
ax_rest <- t(apply(
  ax[, XXX], 1, function(y) 
    c(
      mean(y[ bx[XXX,]$lifestyle == "HG" ], na.rm = T), 
      mean(y[ bx[XXX,]$lifestyle == "NF" ], na.rm = T) ) ))

### pairwise comparisons between HG-NF differences across 3 datasets

# remove loci not covered in both datasets
XXX = !is.na(ax_marc[,1]+ax_marc[,2]+ ax_anto[,1]+ax_anto[,2])
sum(XXX)
# correlation between HG-NF differences
cor.test(ax_marc[XXX,1]-ax_marc[XXX,2], ax_anto[XXX,1]-ax_anto[XXX,2], m="s")

XXX = !is.na(ax_laza[,1]+ax_laza[,2]+ ax_anto[,1]+ax_anto[,2])
sum(XXX)
cor.test(ax_laza[XXX,1]-ax_laza[XXX,2], ax_anto[XXX,1]-ax_anto[XXX,2], m="s")

XXX = !is.na(ax_laza[,1]+ax_laza[,2]+ ax_marc[,1]+ax_marc[,2])
sum(XXX)
cor.test(ax_laza[XXX,1]-ax_laza[XXX,2], ax_marc[XXX,1]-ax_marc[XXX,2], m="s")

ax_nomarch <- t(apply(
  ax[, noMarchi,], 1, function(y) 
    c(
      mean(y[ bx[noMarchi,]$lifestyle == "HG" ], na.rm = T), 
      mean(y[ bx[noMarchi,]$lifestyle == "NF" ], na.rm = T) ) ))

XXX = !is.na(ax_nomarch[,1]+ax_nomarch[,2]+ ax_marc[,1]+ax_marc[,2])
sum(XXX)
cor.test(ax_nomarch[XXX,1]-ax_nomarch[XXX,2], ax_marc[XXX,1]-ax_marc[XXX,2], m="s")

############

# HG-NF differences w/o marchi
diff = ax_nomarch[,2]-ax_nomarch[,1]
diff = na.omit(diff)
length(diff)
head(diff)


### Comparison w Fagny et al.

# Eastern
d <- read.table("/mnt/NAS/projects/2019_epipaleomix/epipaleo/Methylation_AGR_RHG_DMR/wAGR_vs_wRHG_DMR.txt", head = T, fill=T, sep="\t")  
# NF-HG, from Fagny et al, shared w Dilek 
colnames(d)
head(d)
dim(d)

# multiple genes per loci: take the first gene
x = sapply(d$UCSC_REFGENE_NAME, function(x) strsplit(x, ";")[[1]][1])
xx = merge(
    cbind(names(diff), diff),
    cbind(x, d$logFC),
    by.x=1, by.y=1)
head(xx)
xxx = apply(as.matrix(xx[,-1]), 2, as.numeric)
rownames(xxx) = xx[,1]
head(xxx)
xxxx = t(sapply(unique(rownames(xxx)), function(genex) {
  apply(xxx[genex,,drop=F], 2, mean)
}))
head(xxxx)
dim(xxxx)
cor.test(xxxx[,1], xxxx[,2], m="s")

d <- read.table("/mnt/NAS/projects/2019_epipaleomix/epipaleo/Methylation_AGR_RHG_DMR/eAGR_vs_eRHG_DMR.txt", head = T, fill=T, sep="\t")  
# NF-HG, from Fagny et al, shared w Dilek 

colnames(d)
head(d)
dim(d)

# multiple genes per loci: take the first gene

x = sapply(d$UCSC_REFGENE_NAME, function(x) strsplit(x, ";")[[1]][1])
xx = merge(
  cbind(names(diff), diff),
  cbind(x, d$logFC),
  by.x=1, by.y=1)
head(xx)
xxx = apply(as.matrix(xx[,-1]), 2, as.numeric)
rownames(xxx) = xx[,1]
head(xxx)
xxxx = t(sapply(unique(rownames(xxx)), function(genex) {
  apply(xxx[genex,,drop=F], 2, mean)
}))
head(xxxx)
dim(xxxx)
cor.test(xxxx[,1], xxxx[,2], m="s")



###### PARALLEL DIFFERENCES IN DIFF MET ACROSS DATASETS ?

# compare lazaridis, antonioni, marchi, and the rest

table(bx$Article, bx$lifestyle)


# HG-NF differences in lazaridis, antonio, marchi, and all others
ax_diffs <- sapply(
  list(ax_laza, ax_anto, ax_marc, ax_rest), function(x) {
  x[,1]-x[,2]
})
colnames(ax_diffs) <- c("laza", "anto", "marc", "rest")
dim(ax_diffs)
apply(!is.na(ax_diffs), 2, sum)

# antonio has many NAs, also 0s. shift all to NA and remove
ax_diffs[ax_diffs == 0] = NA
ax_diffs <- ax_diffs[!is.na(rowSums(ax_diffs)),]
dim(ax_diffs)

ax_diffs_combs <- table ( apply( sign(ax_diffs), 1, function(x) paste(x, collapse = " ") ) )
rev( sort( ax_diffs_combs ) )

table(sign(ax_diffs))

# freq of -1
ax_diffs_minus_freq <- table(sign(ax_diffs))["-1"]/length(ax_diffs)
ax_diffs_minus_freq

binom.test(ax_diffs_combs["-1 -1 -1 -1"], n, ax_diffs_minus_freq^4)

binom.test(ax_diffs_combs["1 1 1 1"], n, (1-ax_diffs_minus_freq)^4)

mx <- matrix( c(
  ax_diffs_combs["-1 -1 -1 -1"], n*ax_diffs_minus_freq^4,
  ax_diffs_combs["1 1 1 1"], n*(1-ax_diffs_minus_freq)^4), 2, 2)
colnames(mx) = c("NF>HG", "HG>NF")
rownames(mx) = c("Obs.", "Exp.")
par(mfrow=c(1,1))
barplot(mx, bes=T, legend.text = rownames(mx))

### without "rest" (i.e. only antonio, marchi and lazaridis)

ax_diffs_combs2 <- table ( apply( sign(ax_diffs[,1:3]), 1, function(x) paste(x, collapse = " ") ) )
rev( sort( ax_diffs_combs2 ) )


table(sign(ax_diffs[,1:3]))


# freq of -1
ax_diffs_minus_freq2 <- table(sign(ax_diffs[,1:3]))["-1"]/length(ax_diffs[,1:3])
ax_diffs_minus_freq2

binom.test(ax_diffs_combs2["-1 -1 -1"], n, ax_diffs_minus_freq2^4)

binom.test(ax_diffs_combs2["1 1 1"], n, (1-ax_diffs_minus_freq2)^4)

mx <- matrix( c(
  ax_diffs_combs2["-1 -1 -1"], n*ax_diffs_minus_freq2^4,
  ax_diffs_combs2["1 1 1"], n*(1-ax_diffs_minus_freq2)^4), 2, 2)
colnames(mx) = c("NF>HG", "HG>NF")
rownames(mx) = c("Obs.", "Exp.")
par(mfrow=c(1,1))
barplot(mx, bes=T, legend.text = rownames(mx))



########

###### Compare w Fagny et al 2015 datasets

### western HG-AGR

d1 <- read.table("/mnt/NAS/projects/2019_epipaleomix/epipaleo/Methylation_AGR_RHG_DMR/wAGR_vs_wRHG_DMR.txt", head = T, fill=T, sep="\t")  
# NF-HG, from Fagny et al, shared w Dilek 
# convert to HG-NF scale
d1$logFC <- d1$logFC*-1
colnames(d1)
head(d1)
dim(d1)

head(ax_diffs)

# multiple genes per loci: take the first gene
genenamesx = as.character( sapply(d1$UCSC_REFGENE_NAME, function(x) strsplit(x, ";")[[1]][1]) )

west_diff = tapply(as.numeric(d1$logFC), genenamesx, function(y) mean(y) )
head(west_diff)
length(west_diff)
summary(west_diff)

### eastern HG-AGR

d2 <- read.table("/mnt/NAS/projects/2019_epipaleomix/epipaleo/Methylation_AGR_RHG_DMR/eAGR_vs_eRHG_DMR.txt", head = T, fill=T, sep="\t")  
# NF-HG, from Fagny et al, shared w Dilek 
# convert to HG-NF scale
d2$logFC <- d2$logFC*-1
colnames(d2)
head(d2)
dim(d2)

# multiple genes per loci: take the first gene
genenamesx = as.character( sapply(as.character(d2$UCSC_REFGENE_NAME), function(x) strsplit(x, ";")[[1]][1]) )

east_diff = tapply(as.numeric(d2$logFC), genenamesx, function(y) mean(y) )
head(east_diff)
length(east_diff)
summary(east_diff)

### both west and east

xx = merge(
  cbind(rownames(ax_diffs[,1:3]), ax_diffs[,1:3]),
  cbind(names(east_diff), east_diff),
  by.x=1, by.y=1)
xx = merge(
  xx,
  cbind(names(west_diff), west_diff),
  by.x=1, by.y=1)
dim(xx)
head(xx)
xxx = apply(as.matrix(xx[,-1]), 2, as.numeric)
rownames(xxx) = xx[,1]
head(xxx)

colnames(xxx) = c("Boston", "Stanford", "Mainz", "ECAfrica", "WCAfrica")
pairs(xxx[,c(4:5,1:3)], upper.panel = panel.cor, col = "light blue", cex=0.5, pch=19, lower.panel=panel.smooth)


round( as.dist(cor(xxx, met="s")), 3)


xxxx = apply(xxx, 1, function(x) paste(sign(x), collapse=" "))
rev(sort(table(xxxx)))

table(sign(xxx))


binom.test(table(xxxx)["-1 -1 -1 -1 -1"], length(xxxx), mean(sign(xxx) == -1)^5)

binom.test(table(xxxx)["1 1 1 1 1"], length(xxxx), (1-mean(sign(xxx) == -1))^5)


# compare shared differences w permutations

ax_diffs2 = ax_diffs[,1:3]
ax_diffs2s = sign(ax_diffs2)
head(ax_diffs2s)

apply(ax_diffs2s, 2, table)

permx <- t(sapply(1:100, function(i) {
  perm_ax_diffs2s <- apply(ax_diffs2s, 2, sample)  # sample across columns
  perm_ax_diffs_combs2 <- table ( apply( sign(perm_ax_diffs2s), 1, function(x) paste(x, collapse = " ") ) ) # summarize direction of change
  }
) )
head(permx)

par(mfrow=c(1,1))
boxplot(ax_diffs_combs2["-1 -1 -1"], permx[,1], ax_diffs_combs2["1 1 1"], permx[,2], xaxt="n", col=c(1,"lightgrey"), border=c(1,"lightgrey"), ylab="No. consistent DMR genes")
axis(1, at=c(1.5, 3.5), labels = c("NF>HG", "HG>NF"))
legend(3, 250, legend = c("Observed", "Expected"), fill = c(1,"lightgrey"), bty="n", border=c(1,"lightgrey"))

### GO ANALYSIS

# define gene sets

# all -1 genes: NF>HG
ax_diffs2_neg_genes <- names(which(rowSums( ax_diffs2s == -1) == 3))
length(ax_diffs2_neg_genes)

# all +1 genes: HG>NF
ax_diffs2_pos_genes <- names(which(rowSums( ax_diffs2s == 1) == 3))
length(ax_diffs2_pos_genes)

###

require(topGO)
require(org.Hs.eg.db)

args = commandArgs(trailingOnly=TRUE)
# args[1]: set working directory
setwd(args[1])
genes=rowSums( ax_diffs2s == -1)

# all -1 genes: NF>HG
selection <- function(i){ return(i == 3)} 
allGO2genes <- annFUN.org(whichOnto="BP", feasibleGenes=NULL, mapping="org.Hs.eg.db", ID="symbol")
GOdata <- new("topGOdata",
              ontology="BP",
              allGenes=genes,
              annot=annFUN.GO2genes,
              GO2genes=allGO2genes,
              geneSel=selection, nodeSize=5)
results.fisher <- runTest(GOdata, algorithm="elim", statistic="fisher")
d=geneData(results.fisher)[[1]]
c=geneData(results.fisher)[[2]]
goEnrichment <- GenTable(GOdata, Fisher=results.fisher, orderBy="Fisher", topNodes=length(results.fisher@score))
goEnrichment$Fisher <- as.numeric(goEnrichment$Fisher)
b=goEnrichment$Annotated
a=goEnrichment$Significant
goEnrichment$odds = (a/(c-a))/((b-a)/(d-c-b+a))
head(goEnrichment)
goEnrichment[goEnrichment$Fisher<0.01,]
sum(p.adjust(goEnrichment$Fisher)<0.1)

# all 1 genes: HG>NF
selection <- function(i){ return(i == 0)} 
allGO2genes <- annFUN.org(whichOnto="BP", feasibleGenes=NULL, mapping="org.Hs.eg.db", ID="symbol")
GOdata <- new("topGOdata",
              ontology="BP",
              allGenes=genes,
              annot=annFUN.GO2genes,
              GO2genes=allGO2genes,
              geneSel=selection, nodeSize=5)
results.fisher <- runTest(GOdata, algorithm="elim", statistic="fisher")
d=geneData(results.fisher)[[1]]
c=geneData(results.fisher)[[2]]
goEnrichment <- GenTable(GOdata, Fisher=results.fisher, orderBy="Fisher", topNodes=length(results.fisher@score))
goEnrichment$Fisher <- as.numeric(goEnrichment$Fisher)
b=goEnrichment$Annotated
a=goEnrichment$Significant
goEnrichment$odds = (a/(c-a))/((b-a)/(d-c-b+a))
head(goEnrichment)
goEnrichment[goEnrichment$Fisher<0.01,]
sum(p.adjust(goEnrichment$Fisher)<0.1)