result_df <- data.frame()
for (i in 1:20){
	data = read.table(paste0("reshaped_ds_means_", i), head=T)
	genes = data$V6
	data$V1_mean = NULL
	data$V2_mean = NULL
	data$V6 = NULL
	if (nrow(result_df) == 0) {
    		result_df <- data  
  	} else {
    		result_df <- result_df + data
        }
}
result_df <- result_df / 20	
all = cbind(Gene=genes, result_df)
write.table(all, file="reshaped_ds_all_means_v2",quote=F, col.names=T, row.names=F, sep="\t")
	
