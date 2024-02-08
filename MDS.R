unique_lab_colors <- c("limegreen", "royalblue", "black", "violetred", "orange")

pdf("MDS_Plots.pdf", width = 14, height = 7)

par(mfrow = c(1, 2))
plot(mds, type = "n", xlab = "Dimension 1", ylab = "Dimension 2")
points(mds, col = ifelse(bx$lifestyle == "HG", "blue", "red"),
       pch = 17, cex = 2)
legend(6, 3, legend = c("HG", "NF"), pch = 17, col = c("blue", "red"), bty = "n")
mtext("A", side = 3, adj = -0.12, cex = 2.15, font=2, line=2)  # Add label "A" and left-align

plot(mds, type = "n", xlab = "Dimension 1", ylab = "Dimension 2")
points(mds, col = unique_lab_colors[as.numeric(factor(bx$Lab))],
       pch = 18, cex = 2)
legend(6, 1.5, legend = levels(bx$Lab), pch = 18, bty = "n", col = unique_lab_colors)
mtext("B", side = 3, adj = -0.12, cex = 2.15, font=2, line=2)  # Add label "B" and left-align
par(mfrow = c(1, 1))

dev.off()