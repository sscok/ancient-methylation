library(dplyr)
library(reshape2)



for (i in 1:20){
	data=read.table(paste0("all_replicated_ds_sorted_a_",i), head=F)
	result <- aggregate(V3 ~ V1 + V2 + V4 + V6, data, FUN=mean, na.rm=TRUE)
	reshaped_data <- dcast(result, V1 + V2 + V6 ~ V4, value.var = "V3")
	agg_data <- reshaped_data %>%
   		group_by(V6) %>%
   		summarise_all(list(mean = ~mean(., na.rm = TRUE)))
	print(dim(agg_data))
	write.table(agg_data,file = paste0("reshaped_ds_means_",i), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
}