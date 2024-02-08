# methylation dataset downsampled (normalized) 20x, and mean methylation score (MS) per gene calculated, and means of 20x replicates calculated (NAs not removed). by seda cokoglu, Aug 2022

a = read.table("all_downsampled_20_means_v2", row.names = 1, head = T)  
colnames(a)
head(a)
dim(a) # 10170    34

# info table. from seda's table at https://docs.google.com/spreadsheets/d/1s3L1OhIWnl5258WMSdMJDPveJ9igka72GIYqVZ5_lqg/edit#gid=0, but some sample and article names modified to match that in the methylation table

b = read.table("Shotgun_inds.tsv", row.names = 1, head = T, fill=T, sep="\t")
head(b)
setdiff(colnames(a), rownames(b)) # 0, i.e. no difference 
setdiff(rownames(b), colnames(a)) # "irk034"
dim(b)  # 34 13

b <- b[colnames(a), ]
rownames(b) == colnames(a)

levels(factor(b$Article))
b$Article <- factor(b$Article, levels=c("Lazaridis2014", "Marchi2022", "Fu2014", "Gunther2018", "Quinto2019"))

table(b$lifestyle)
# HG NF 
# 13 21 
table(b$tissue)
# bone tooth 
# 23    11 
table(b$genetic_sex)
# XX XY 
# 12 22 
table(b$Article)

table(b$Article, b$tissue)
#                    bone tooth
# Antonio2019           5     0
# Lazaridis2014         0     3
# Marchi2022           15     0
# Fu2014                1     0
# Gunther2018           1     0
# Kilinc2021            0     4
# Seguin-Orlando2014    1     0
# Quinto2019            0     4

table(b$genetic_sex, b$tissue, b$lifestyle)
# , ,  = HG
# bone tooth
# XX    1     1
# XY    6     5
# 
# , ,  = NF
# bone tooth
# XX    8     2
# XY    8     3

### MDS

# TO REPEAT

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
# no correlation btw group and coverage
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


### anova 

# TO REPEAT

# res: only tissue effect
# HG_columns <- c("Sf12_mean","Loschbour_mean","UstIshim_mean","VLASA7_mean", "VLASA32_mean")
# NF_columns <- c("STAR1_mean","AKT16_mean", "Asp6_mean", "BAR25_mean", "Dil16_mean", "Ess7_mean","Herx_mean", "VC3.2_mean","Klein7_mean", "LBK_mean", "LEPE48_mean", "LEPE52_mean", "Nea2_mean","Nea3_mean","prs013_mean")

# HG_presence <- apply(ax[, HG_columns], 1, function(row) any(!is.na(row)))
# NF_presence <- apply(ax[, NF_columns], 1, function(row) any(!is.na(row)))

# # Combine the conditions to filter genes with at least one representative in each group
# gene_filter <- HG_presence & NF_presence

# # Apply the filter to obtain the filtered dataset
# filtered_dataset <- ax[gene_filter, c(HG_columns, NF_columns)]

# # Print the filtered dataset
# print(filtered_dataset)


ax2 = ax[rowSums(!is.na(ax)) >= 20,]

dim(ax2) # 4717   33
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

# res: tissue effect gone

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
#                     HG NF
# Antonio2019         2  3
# Lazaridis2014       1  1
# Marchi2022          2 13
# Fu2014              1  0
# Gunther2018         1  0
# Kilinc2021          4  0
# Seguin-Orlando2014  1  0
# Quinto2019          0  4

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
dim(d) # 365401      9

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
# S = 9916086852, p-value = 0.5628
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# 0.009251208

# Eastern

d <- read.table("/mnt/NAS/projects/2019_epipaleomix/epipaleo/Methylation_AGR_RHG_DMR/eAGR_vs_eRHG_DMR.txt", head = T, fill=T, sep="\t")  
# NF-HG, from Fagny et al, shared w Dilek 

colnames(d)
head(d)
dim(d) # 365401      9

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
# S = 1.0338e+10, p-value = 0.03934
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# -0.03293041 


###### PARALLEL DIFFERENCES IN DIFF MET ACROSS DATASETS ?

# compare lazaridis, antonioni, marchi, and the rest

# res: yes

table(bx$Article, bx$lifestyle)
#                     HG NF
# Antonio2019         2  3
# Lazaridis2014       1  1
# Marchi2022          2 13
# Fu2014              1  0
# Gunther2018         1  0
# Kilinc2021          4  0
# Seguin-Orlando2014  1  0
# Quinto2019          0  4

# HG-NF differences in lazaridis, antonio, marchi, and all others
ax_diffs <- sapply(
  list(ax_laza, ax_anto, ax_marc, ax_rest), function(x) {
  x[,1]-x[,2]
})
colnames(ax_diffs) <- c("laza", "anto", "marc", "rest")
dim(ax_diffs)
apply(!is.na(ax_diffs), 2, sum)
# laza anto marc rest 
# 4909 2977 4953 4575 

# antonio has many NAs, also 0s. shift all to NA and remove
ax_diffs[ax_diffs == 0] = NA
ax_diffs <- ax_diffs[!is.na(rowSums(ax_diffs)),]
dim(ax_diffs)
# 1355    4

ax_diffs_combs <- table ( apply( sign(ax_diffs), 1, function(x) paste(x, collapse = " ") ) )
rev( sort( ax_diffs_combs ) )
# -1 -1 -1 1    1 -1 1 1   -1 -1 1 1  1 -1 -1 -1   1 -1 -1 1  -1 -1 1 -1 -1 -1 -1 -1 
# 147         134         122         119         113         102         101 
# 1 -1 1 -1  -1 1 -1 -1    1 1 -1 1   1 1 -1 -1   -1 1 -1 1     1 1 1 1    -1 1 1 1 
# 80          67          64          63          61          55          55 
# -1 1 1 -1    1 1 1 -1 
# 42          30 

# more neg than positive
table(sign(ax_diffs))
# -1    1 
# 2954 2466 

# freq of -1
ax_diffs_minus_freq <- table(sign(ax_diffs))["-1"]/length(ax_diffs)
ax_diffs_minus_freq

binom.test(ax_diffs_combs["-1 -1 -1 -1"], n, ax_diffs_minus_freq^4)
# number of successes = 101, number of trials = 1362, p-value = 0.06935
# alternative hypothesis: true probability of success is not equal to 0.0882358
# 95 percent confidence interval:
#   0.06080378 0.08937777
# sample estimates:
#   probability of success 
# 0.07415565 

binom.test(ax_diffs_combs["1 1 1 1"], n, (1-ax_diffs_minus_freq)^4)
# number of successes = 55, number of trials = 1362, p-value = 0.7377
# alternative hypothesis: true probability of success is not equal to 0.0428524
# 95 percent confidence interval:
#   0.03056337 0.05224048
# sample estimates:
#   probability of success 
# 0.04038179 

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
# -1 -1 -1  1 -1 -1  -1 -1 1   1 -1 1  -1 1 -1   1 1 -1   -1 1 1    1 1 1 
# 248      232      224      214      128      127       97       85 

table(sign(ax_diffs[,1:3]))
# -1    1 
# 2350 1715 

# freq of -1
ax_diffs_minus_freq2 <- table(sign(ax_diffs[,1:3]))["-1"]/length(ax_diffs[,1:3])
ax_diffs_minus_freq2

binom.test(ax_diffs_combs2["-1 -1 -1"], n, ax_diffs_minus_freq2^4)
# number of successes = 248, number of trials = 1362, p-value = 1.757e-14
# alternative hypothesis: true probability of success is not equal to 0.1116938
# 95 percent confidence interval:
#   0.1619260 0.2036258
# sample estimates:
#   probability of success 
# 0.1820852 

binom.test(ax_diffs_combs2["1 1 1"], n, (1-ax_diffs_minus_freq2)^4)
# number of successes = 85, number of trials = 1362, p-value = 1.013e-08
# alternative hypothesis: true probability of success is not equal to 0.03168212
# 95 percent confidence interval:
#   0.05015028 0.07659322
# sample estimates:
#   probability of success 
# 0.06240822 

mx <- matrix( c(
  ax_diffs_combs2["-1 -1 -1"], n*ax_diffs_minus_freq2^4,
  ax_diffs_combs2["1 1 1"], n*(1-ax_diffs_minus_freq2)^4), 2, 2)
colnames(mx) = c("NF>HG", "HG>NF")
rownames(mx) = c("Obs.", "Exp.")
par(mfrow=c(1,1))
barplot(mx, bes=T, legend.text = rownames(mx))



########

###### Compare w Fagny

# TO REPEAT

### western HG-AGR

d1 <- read.table("/mnt/NAS/projects/2019_epipaleomix/epipaleo/Methylation_AGR_RHG_DMR/wAGR_vs_wRHG_DMR.txt", head = T, fill=T, sep="\t")  
# NF-HG, from Fagny et al, shared w Dilek 
# convert to HG-NF scale
d1$logFC <- d1$logFC*-1
colnames(d1)
head(d1)
dim(d1) # 365401      9

head(ax_diffs)

# multiple genes per loci: take the first gene
genenamesx = as.character( sapply(d1$UCSC_REFGENE_NAME, function(x) strsplit(x, ";")[[1]][1]) )

west_diff = tapply(as.numeric(d1$logFC), genenamesx, function(y) mean(y) )
head(west_diff)
length(west_diff)
# 19124
summary(west_diff)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.700982 -0.013181  0.005287  0.004527  0.023774  1.376667 

### eastern HG-AGR

d2 <- read.table("/mnt/NAS/projects/2019_epipaleomix/epipaleo/Methylation_AGR_RHG_DMR/eAGR_vs_eRHG_DMR.txt", head = T, fill=T, sep="\t")  
# NF-HG, from Fagny et al, shared w Dilek 
# convert to HG-NF scale
d2$logFC <- d2$logFC*-1
colnames(d2)
head(d2)
dim(d2) # 365401      9

# multiple genes per loci: take the first gene
genenamesx = as.character( sapply(as.character(d2$UCSC_REFGENE_NAME), function(x) strsplit(x, ";")[[1]][1]) )

east_diff = tapply(as.numeric(d2$logFC), genenamesx, function(y) mean(y) )
head(east_diff)
length(east_diff)
# 19124
summary(east_diff)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.456509 -0.020057  0.001612  0.001636  0.023229  1.171282 

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
#             laza   anto   marc east_diff
# anto      -0.006                        
# marc       0.018 -0.062                 
# east_diff -0.002  0.005 -0.051          
# west_diff -0.012  0.010 -0.035     0.406

#

xxxx = apply(xxx, 1, function(x) paste(sign(x), collapse=" "))
rev(sort(table(xxxx)))
# -1 -1 -1 1 1  -1 -1 1 -1 -1  1 -1 -1 -1 -1    -1 -1 1 1 1 -1 -1 -1 -1 -1 
# 78             65             62             62             62 
# 1 -1 1 1 1   1 -1 1 -1 -1    1 -1 -1 1 1  -1 -1 -1 -1 1   1 -1 -1 -1 1 
# 60             58             53             46             45 
# 1 -1 1 -1 1     1 1 -1 1 1    -1 1 -1 1 1   -1 -1 1 -1 1   -1 1 1 -1 -1 
# 44             39             34             34             29 
# -1 1 -1 -1 -1      1 1 1 1 1   1 1 -1 -1 -1    1 -1 1 1 -1   -1 -1 1 1 -1 
# 28             25             25             25             25 
table(sign(xxx))
# -1    1 
# 3070 2610 

binom.test(table(xxxx)["-1 -1 -1 -1 -1"], length(xxxx), mean(sign(xxx) == -1)^5)
# number of successes = 62, number of trials = 1136, p-value = 0.1784
# alternative hypothesis: true probability of success is not equal to 0.04612647
# 95 percent confidence interval:
#   0.04209682 0.06942308
# sample estimates:
#   probability of success 
# 0.05457746 

binom.test(table(xxxx)["1 1 1 1 1"], length(xxxx), (1-mean(sign(xxx) == -1))^5)
# number of successes = 25, number of trials = 1136, p-value = 0.6749
# alternative hypothesis: true probability of success is not equal to 0.02048617
# 95 percent confidence interval:
#   0.01429122 0.03231622
# sample estimates:
#   probability of success 
# 0.02200704 

####### TO REPEAT

# compare shared differences w permutations

ax_diffs2 = ax_diffs[,1:3]
ax_diffs2s = sign(ax_diffs2)
head(ax_diffs2s)

apply(ax_diffs2s, 2, table)
#    laza anto marc
# -1  697  918  735
# 1   658  437  620

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
# 248

# all +1 genes: HG>NF
ax_diffs2_pos_genes <- names(which(rowSums( ax_diffs2s == 1) == 3))
length(ax_diffs2_pos_genes)
# 85

###

require(topGO)
require(org.Hs.eg.db)

args = commandArgs(trailingOnly=TRUE)
# args[1]: set working directory
# args[2]: corrected anova table input file
# args[3]: gene set option either "g" or "gts", g for Subsistence type significant, gts for Subsistence type significant, others non-significant 
# args[4]: output file name
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
# GO.ID                                        Term Annotated Significant
# 1 GO:0016574                      histone ubiquitination         5           4
# 2 GO:0031058 positive regulation of histone modificat...        11           6
# 3 GO:1903320 regulation of protein modification by sm...        18           8
# Expected Fisher      odds
# 1     0.91 0.0045 18.381395
# 2     1.99 0.0068  5.543662
# 3     3.26 0.0085  3.711848
sum(p.adjust(goEnrichment$Fisher)<0.1)
# 0

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
# GO.ID                         Term Annotated Significant Expected Fisher
# 1 GO:0006801 superoxide metabolic process         5           3     0.32 0.0024
# odds
# 1 22.56
sum(p.adjust(goEnrichment$Fisher)<0.1)
# 0


#######

# 34 genom - 0.02 normalize et -> anova yurut -> 
# 9967 gen, BH corr 0.05, grup 66, 25 ?
# permutasyon 100: gen sayisi ?
# GO analizi: 73 categori - gelisim
# perm 100: anlamli

