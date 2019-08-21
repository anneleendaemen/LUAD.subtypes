library(Biobase)
library(Rtsne)

load("../data/Celllines_RNAseq_ClassifierGenes.RData")
data <- t(log2(exprs(eset)+1))
colors <- rep("cyan3", ncol(eset))
colors[eset$Subtype %in% "PRO"] <- "mediumpurple2"
colors[eset$Subtype %in% "MES"] <- "green3"
colors[eset$Subtype %in% "unclassified"] <- "grey"

metadata <- read.table("../data/CellLines_KRASall_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
metadata$Subtype <- factor(metadata$Subtype, levels=c("MES", "PRO", "MUC", "unclassified"))
metadata <- metadata[match(rownames(data), metadata$SAM.ID), ]
pchs <- rep(17, nrow(data))
pchs[metadata$KrasMutStatus] <- 19

pca <- prcomp(data, scale. = TRUE)

pdf("../figures/FigS4_Celllines_SubtypePCA.pdf", width = 5.5, height = 6.5)
par(mar=c(4.5, 4.5, 3, 1))
plot(pca$x[,1], pca$x[,2], xlab="PC1", ylab="PC2", 
		main = paste0("Cell lines (n=", ncol(eset), ")"),
		pch = pchs, col = colors, 
		cex.main = 1.8, cex.lab = 1.8, cex.axis = 1.5, cex = 1.5)
dev.off()

tsne <- Rtsne(data, dims=2, initial_dims=10, perplexity=30, theta = 0.2,
						pca_center = TRUE, pca_scale = TRUE,
						check_duplicates = FALSE, verbose=TRUE, max_iter=2000)

pdf("../figures/FigS4_Celllines_SubtypetSNE 10.pdf", width = 5.5, height = 6.5)
par(mar=c(4.5, 4.5, 3, 1))
plot(tsne$Y, xlab="tSNE 1", ylab="tSNE 2", 
		main = paste0("Cell lines (n=", ncol(eset), ")"),
		pch = pchs, col = colors, 
		cex.main = 1.8, cex.lab = 1.8, cex.axis = 1.5, cex = 1.5)
dev.off()
