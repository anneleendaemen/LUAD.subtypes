library(NMF)

load(file= "../data/TCGAwt_RNAseq_BimodalVariableGenes.RData")
data <- exprs(eset.data) # log2-transformed

### NMF code to assign de novo subtype labels
#NMF.opt <- nmf(data, rank=3, method='brunet', nrun=200)
#cluster <- predict(NMF.opt,prob=TRUE)
#de.novo <- as.character(cluster$predict)
#sil <- silhouette(as.numeric(cluster$predict), dmatrix=(1-consensus(NMF.opt)))
#rownames(sil) <- names(cluster$predict)
#indices.unclassifiable <- which(sil[,"sil_width"]<0.5)
#de.novo[indices.unclassifiable] <- NA
#names(de.novo) <- names(cluster$predict)

estim.r <- nmfEstimateRank(data, range=3, method='brunet', nrun=200)
annCol <- pData(eset.data)
annCol[,2] <- as.character(annCol[,2])
annCol[is.na(annCol[,2]),2] <- "unclassified"
colnames(annCol)[2] <- "de novo"
nbTumors <- nrow(annCol)

annColors <- list(subtype = c("MUC" = "cyan3", "PRO" = "mediumpurple2", "MES" = "green3", "unclassified" = "gray"),
                  `de novo` = c("1" = "cyan3", "2" = "mediumpurple2", "3" = "green3", "unclassified" = "gray"))

pdf("../figures/FigS1_TCGAwt_consensusmap_200runs.pdf", width = 6, height = 5.5, onefile=FALSE)
par(mar=c(2,2,4,1))
consensusmap(estim.r, annCol = annCol, annColors = annColors, labCol=NA, labRow=NA, tracks = NA,
						main = paste0("TCGA - KRAS wildtype (n=", nbTumors, ")"))
dev.off()
