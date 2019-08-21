data <- read.table("../data/GEMM_KRASmut_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
data <- data[!data$Subtype %in% "unclassified", ]
data$Subtype <- factor(as.character(data$Subtype), levels=c("PRO", "MUC"))
prev.control <- table(data$Subtype)

data <- read.table("../data/GEMM_MEKi_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
data <- data[!data$Subtype %in% "unclassified", ]
data$Subtype <- factor(as.character(data$Subtype), levels=c("PRO", "MUC"))
prev.MEKi <- table(data$Subtype)

t <- rbind(prev.control, prev.MEKi)
rownames(t) <- c("untreated", "MEKi")
prev <- t(prop.table(t, 1)*100)
fisher.test(t, alternative="less")

pdf("../figures/Fig5_GEMM_MEKi_SubtypePrevalence.pdf", width=4, height=6.5)
par(mar=c(4.5, 4.5, 3, 1), xpd=TRUE)
bp <- barplot(prev, xlab="", ylab="Tumor prevalence", col=c("mediumpurple2", "cyan3"), 
			cex.axis=1.5, cex.lab=1.8, cex.names=1.5, ylim=c(0,100), border=FALSE)
lines(x=c(bp[1], bp[2]), y=rep(102,2), lwd=2)
text(x=(bp[1]+bp[2])/2, y=105, labels="*", cex=2)
dev.off()
