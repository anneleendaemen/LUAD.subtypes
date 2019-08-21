library(Biobase)
library(NMF)

load("../data/TCGAmut_RNAseq_NMFGenes.RData")
data <- exprs(eset.data) # log2-transformed
nbTumors <- ncol(data)

# Estimation of the factorization rank
estim.r <- nmfEstimateRank(data, range=2:10, method='brunet', nrun=50, .opt='vp60')

# Plot NMF results across number of clusters varying from 2 to 10
pdf("../figures/FigS1_TCGAmut_NMFcoefficients.pdf", width = 6, height = 5)
par(mfrow=c(1,2), mar=c(4.5,5.5,1,0.5), oma=c(0,0,2,0))
plot(2:10, estim.r$measures[,"cophenetic"], pch=15, xlab="Number of clusters",
     ylab="Cophenetic coefficient\n", las=2, main=NULL, cex.lab=1.3, cex.axis=1.3)
lines(2:10, estim.r$measures[,"cophenetic"])

plot(2:10, estim.r$measures[,"dispersion"], pch=15, xlab="Number of clusters",
     ylab="Dispersion coefficient\n", las=2, main=NULL, cex.lab=1.3, cex.axis=1.3)
lines(2:10, estim.r$measures[,"dispersion"])

mtext(paste0("TCGA - KRAS mutant (n=",nbTumors,")"), outer=TRUE, cex=1.5)
dev.off()
