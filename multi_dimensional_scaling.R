# Description: This script plots the Multi-Dimensional Scaling Analysis results 
# for both subsistence type and laboratory-of-origin (Figure 3).
# inputs: 
#	- MS means file name: "reshaped_ds_all_means_v2" which include means of 20 downsampled datasets per gene 
#      - info table file name: from the table at https://docs.google.com/spreadsheets/d/1s3L1OhIWnl5258WMSdMJDPveJ9igka72GIYqVZ5_lqg/edit#gid=0
# output:
#	- MDS plot

unique_lab_colors <- c("limegreen", "royalblue", "black", "violetred", "orange")
data = read.table("reshaped_ds_all_means_v2", row.names = 1, head = T)
info = read.table("Shotgun_inds.tsv", row.names = 1, head = T, fill=T, sep="\t")
mds <- cmdscale(dist(t(data)))

pdf("MDS_Plots.pdf", width = 14, height = 7)

par(mfrow = c(1, 2))
plot(mds, type = "n", xlab = "Dimension 1", ylab = "Dimension 2")
points(mds, col = ifelse(info$lifestyle == "HG", "blue", "red"),
       pch = 17, cex = 2)
legend(6, 3, legend = c("HG", "NF"), pch = 17, col = c("blue", "red"), bty = "n")
mtext("A", side = 3, adj = -0.12, cex = 2.15, font=2, line=2)  # Add label

plot(mds, type = "n", xlab = "Dimension 1", ylab = "Dimension 2")
points(mds, col = unique_lab_colors[as.numeric(factor(info$Lab))],
       pch = 18, cex = 2)
legend(6, 1.5, legend = levels(info$Lab), pch = 18, bty = "n", col = unique_lab_colors)
mtext("B", side = 3, adj = -0.12, cex = 2.15, font=2, line=2)  # Add label
par(mfrow = c(1, 1))

dev.off()