setwd("/mnt/NAS/projects/2019_epipaleomix/epipaleo/out_shotgun/epipaleo_out/anova_marchi")
individuals = c("AKT16_replicated","Asp6_replicated","BAR25_replicated","Dil16_replicated","Ess7_replicated","Herx_replicated","Klein7_replicated","LBK_replicated","LEPE48_replicated","LEPE52_replicated","Loschbour_replicated","Nea2_replicated","Nea3_replicated","prs013_replicated","Sf12_replicated","STAR1_replicated","UstIshim_replicated","VC3-2_replicated","VLASA32_replicated","VLASA7_replicated")

for (i in individuals){

	dataset = read.table(i, head=F, sep="\t")
	dataset$concatenated <- paste(dataset$V1, dataset$V2, sep = "-")
	filtered_data <- dataset[dataset$concatenated %in% names(which(table(dataset$concatenated) >= 10)), ]
	filtered_data$concatenated = NULL
	write.table(filtered_data, paste0("ratios_10/",i), col.names=F, row.names=F, sep="\t", quote=F)
}

