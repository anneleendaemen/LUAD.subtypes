data <- read.table("../data/TCGA_KRASmut_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
data$Subtype <- factor(data$Subtype, levels=c("MES", "PRO", "MUC", "unclassified"))
prev.TCGAmut <- table(data$Subtype)

data <- read.table("../data/TCGA_KRASwt_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
data$Subtype <- factor(data$Subtype, levels=c("MES", "PRO", "MUC", "unclassified"))
prev.TCGAwt <- table(data$Subtype)

data <- read.table("../data/OAK_KRASmut_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
data$Subtype <- factor(data$Subtype, levels=c("MES", "PRO", "MUC", "unclassified"))
prev.OAKmut <- table(data$Subtype)

data <- read.table("../data/OAK_KRASwt_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
data$Subtype <- factor(data$Subtype, levels=c("MES", "PRO", "MUC", "unclassified"))
prev.OAKwt <- table(data$Subtype)

prev <- rbind(prev.TCGAmut, prev.TCGAwt, prev.OAKmut, prev.OAKwt)
rownames(prev) <- c("mut", "wt", "mut", "wt")
prev <- t(prev)
prev2 <- prop.table(prev, 2)*100

pdf("../figures/Fig2_TCGAOAK_SubtypePrevalence.pdf", width=5.5, height=6.5)
par(mar=c(4.5, 4.5, 3, 1), xpd=TRUE)
bp <- barplot(prev2, xlab="", ylab="Tumor prevalence", col=c("green3", "mediumpurple2", "cyan3", "gray"), 
			cex.axis=1.5, cex.lab=1.8, cex.names=1.5, ylim=c(0,100), border=FALSE)
text(x=1.3, y=105, labels="TCGA", cex=1.8)
lines(x=c(0.2,2.4), y=rep(102,2), lwd=2)
text(x=3.7, y=105, labels="OAK", cex=1.8)
lines(x=c(2.6,4.8), y=rep(102,2), lwd=2)
text(x=-0.4, y=-5.5, labels="KRAS", cex=1.5)
dev.off()
