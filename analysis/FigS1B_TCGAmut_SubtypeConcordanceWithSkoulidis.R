t <- matrix(c(14,0,1,3,20,0,0,2,28), nr=3, byrow=TRUE)
rownames(t) <- c("KC", "KL", "KP")
colnames(t) <- c("1", "2", "3")

pdf("../figures/FigS1_TCGAmut_NMF_SubtypeConcordanceWithSkoulidis.pdf", width=5.5, height=6.5)
par(mar=c(4.5, 4.5, 3, 1))
barplot(t, xlab="De novo subtype", ylab="Number of tumors", col=c("cyan3", "mediumpurple2", "green3"), 
			main="TCGA - KRAS mutant (n=68)", cex.axis=1.5, cex.lab=1.8, cex.names=1.5, cex.main=1.8, ylim=c(0,30))
legend("topleft", legend=rownames(t), fill=c("cyan3", "mediumpurple2", "green3"), cex=1.6, bty="n")
dev.off()

concordance <- sum(diag(t)) / sum(t)
