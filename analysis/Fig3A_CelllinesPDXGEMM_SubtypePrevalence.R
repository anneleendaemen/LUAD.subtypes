data <- read.table("../data/CellLines_KRASall_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
data$Subtype <- factor(data$Subtype, levels=c("MES", "PRO", "MUC", "unclassified"))
prev.Celllines <- table(data$Subtype)

data <- read.table("../data/PDX_KRASall_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
data$Subtype <- factor(data$Subtype, levels=c("MES", "PRO", "MUC", "unclassified"))
prev.PDX <- table(data$Subtype)

data <- read.table("../data/GEMM_KRASmut_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
data$Subtype <- factor(data$Subtype, levels=c("MES", "PRO", "MUC", "unclassified"))
prev.GEMM <- table(data$Subtype)

prev <- rbind(prev.Celllines, prev.PDX, prev.GEMM)
rownames(prev) <- c("Cell\nlines", "PDX", "GEMM")
prev <- t(prev)
prev2 <- prop.table(prev, 2)*100

pdf("../figures/Fig3_CelllinesPDXGEMM_SubtypePrevalence.pdf", width=4, height=6.5)
par(mar=c(4.5, 4.5, 3, 1), xpd=TRUE)
bp <- barplot(prev2, xlab="", ylab="Tumor prevalence", col=c("green3", "mediumpurple2", "cyan3", "gray"), 
			cex.axis=1.5, cex.lab=1.8, cex.names=1.5, ylim=c(0,100), border=FALSE, xaxt="n")
axis(1, at=bp, pos=-5, labels=colnames(prev), lwd=0, cex.axis=1.3)
dev.off()
