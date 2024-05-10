# Description: This script calculates the mean values of 20 downsampled datasets and creates and output file for general_analysis.R
# inputs: 
#	- prefix name of the 20 rehaped downsampled matrices 
# output:
#	- an output file named "reshaped_ds_all_means_v2" which include means of 20 downsampled datasets per gene
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
write.table(all, file="reshaped_ds_all_means",quote=F, col.names=T, row.names=F, sep="\t")
	
