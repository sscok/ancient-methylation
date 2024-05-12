#!/usr/bin/env Rscript
# Description: This script joins individual files into a matrix for all CpG positions used in the study. 
# inputs: 
#	- input_files: individual file names containing MS ratio per individual
# - reference_data: CpG positions in the form of "chr  start  end"
# output:
#	- file output "matrix_sorted_2023", merged matrix
args = commandArgs(trailingOnly=TRUE)

input_files <- c("AKT16","Asp6","BAR25","cta016","Dil16","Ess7","Herx","irk025","irk061","K14","Klein7","kra001","LBK","LEPE48","LEPE52","Loschbour",
  "Motala12","Nea2","Nea3","prs002","prs009","prs013","prs016","R15","R2","R3","R7","R9","Sf12","STAR1","UstIshim","VC3-2","VLASA32","VLASA7")

reference_data <- read.table(args[1], header = FALSE)
colnames(reference_data) <- c("X.chrom", "genomicpos","end")
df_union = reference_data


for (file in input_files) {

  input_data <- read.table(paste0("../", file), header = F)
  colnames(input_data) <- c("X.chrom", "genomicpos","deaminated", "total", "bedcoord", "ratio")
  merged_data <- merge(df_union, input_data, by = c("X.chrom", "genomicpos"), all.x = TRUE)
  ratio_col_name <- paste(file, "ratio", sep = "_")
  df_union[ratio_col_name] <- merged_data$ratio
  rm(input_data, merged_data)
}

inds=c("chr", "start","end", "AKT16","Asp6","BAR25","cta016","Dil16","Ess7","Herx","irk025","irk061","K14","Klein7","kra001","LBK","LEPE48","LEPE52","Loschbour",
  "Motala12","Nea2","Nea3","prs002","prs009","prs013","prs016","R15","R2","R3","R7","R9","Sf12","STAR1","Ust-Ishim","VC3-2","VLASA32","VLASA7")

colnames(df_union) <- inds

print(dim(df_union[!is.na(df_union$Sf12),]))
write.table(df_union, "matrix_sorted_2024", row.names=F, col.names=TRUE, quote=F, sep="\t")

