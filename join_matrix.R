#!/usr/bin/env Rscript
# Description: This script joins individual files into a matrix for all CpG positions used in the study. 
# inputs: 
#	- input_files: individual file names containing MS ratio per individual
# - reference_data: CpG positions in the form of "chr  start  end"
# output:
#	- file output "matrix_sorted_2023", merged matrix


input_files <- c("ind.AKT16","ind.Asp6","ind.BAR25","ind.cta016","ind.Dil16","ind.Ess7","ind.Herx","ind.irk025","ind.irk061","ind.K14","ind.Klein7","ind.kra001","ind.LBK","ind.LEPE48","ind.LEPE52","ind.Loschbour",
  "ind.Motala12","ind.Nea2","ind.Nea3","ind.prs002","ind.prs009","ind.prs013","ind.prs016","ind.R15","ind.R2","ind.R3","ind.R7","ind.R9","ind.Sf12","ind.STAR1","ind.UstIshim","ind.VC3-2","ind.VLASA32","ind.VLASA7")

reference_data <- read.table("/mnt/NAS/projects/2019_epipaleomix/epipaleo/matrix/cpgsite-XY-sorted.bed", header = FALSE)
colnames(reference_data) <- c("X.chrom", "genomicpos","end")
df_union = reference_data


for (file in input_files) {
  # Read the input file
  input_data <- read.table(file, header = TRUE)
  
  # Merge the data from input_data with df_union based on 'chr' and 'start' columns
  merged_data <- merge(df_union, input_data, by = c("X.chrom", "genomicpos"), all.x = TRUE)
  
  # Update the df_union with 'ratio' values from input_data if they exist
  ratio_col_name <- paste(file, "ratio", sep = "_")
  df_union[ratio_col_name] <- merged_data$ratio
  
  # Clean up intermediate data
  rm(input_data, merged_data)
}

inds=c("chr", "start","end", "AKT16","Asp6","BAR25","cta016","Dil16","Ess7","Herx","irk025","irk061","K14","Klein7","kra001","LBK","LEPE48","LEPE52","Loschbour",
  "Motala12","Nea2","Nea3","prs002","prs009","prs013","prs016","R15","R2","R3","R7","R9","Sf12","STAR1","Ust-Ishim","VC3-2","VLASA32","VLASA7")

colnames(df_union) <- inds

print(dim(df_union[!is.na(df_union$Sf12),]))
write.table(df_union, "matrix_sorted_2023", row.names=F, col.names=TRUE, quote=F, sep="\t")

