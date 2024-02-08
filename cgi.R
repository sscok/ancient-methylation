#!/usr/bin/env Rscript
data <- read.table("CGIs", header = FALSE, sep = "\t")
data <- data[, -c(1, 2, 3, 38, 39, 40)]  

colnames(data)[ncol(data)] <- "Category"

library(ggplot2)

categories <- c("CGI", "shores3", "shores5", "shelf3", "shelf5", "OpenSea")

filtered_data <- data[data$Category %in% categories, ]

mean_data <- aggregate(. ~ Category, data = filtered_data, FUN = function(x) mean(x, na.rm = TRUE), na.action=na.pass)

df_transposed <- data.frame(
  Category = rep(mean_data$Category, each = ncol(mean_data) - 1),
  Individual = rep(1:(ncol(mean_data) - 1), times = nrow(mean_data)),
  Value = as.vector(t(mean_data[,-1]))
)
scatterplot <- ggplot(df_transposed, aes(x = factor(Category, levels = rev(categories)), y = Value, color = Category)) +
  geom_jitter(position=position_jitter(0.2)) +
  labs(x = "Category", y = "Mean Methylation Score Per Individual") +
  scale_color_manual(values = c("CGI" = "chocolate4", "shores3" = "coral2", "shores5" = "coral2", "shelf3" = "cornflowerblue", "shelf5" = "cornflowerblue", "OpenSea" = "firebrick4")) +
  stat_summary(fun="median", geom="crossbar", color="black", size=0.5, width=0.5) +
  theme_bw() + 
  theme(axis.text=element_text(size=10.5),legend.position="none",panel.background = element_blank(), axis.title.x=element_blank(), axis.text.y = element_text(vjust = 0.5, hjust=1)) +
  coord_flip() + 
  guides(color = FALSE)
  
W <- 4  

ggsave("categorical_scatterplot1.pdf", plot = scatterplot, width = W, height = 2 * W, dpi = 1000, device = "pdf")
