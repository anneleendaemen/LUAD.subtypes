data <- read.table("../data/TCGA_KRASmut_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
data$Subtype <- factor(data$Subtype, levels=c("MES", "PRO", "MUC", "unclassified"))
prev.TCGA <- table(data$Subtype)

data <- read.table("../data/OAK_KRASmut_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
data$Subtype <- factor(data$Subtype, levels=c("MES", "PRO", "MUC", "unclassified"))
prev.OAK <- table(data$Subtype)
prev.OAK.metastatic <- table(data$Subtype[data[,"Tissue origin"] %in% "METASTATIC"])

prev <- rbind(prev.TCGA, prev.OAK)
rownames(prev) <- c("TCGA", "OAK")
prev <- t(prev)
prev2 <- prop.table(prev, 2)*100

pdf("../figures/Fig1_TCGAOAKmut_SubtypePrevalence.pdf", width=4, height=6.5)
par(mar=c(4.5, 4.5, 3, 1))
barplot(prev2, xlab="", ylab="Tumor prevalence", col=c("green3", "mediumpurple2", "cyan3", "gray"), 
			cex.axis=1.5, cex.lab=1.8, cex.names=1.5, ylim=c(0,100), border=FALSE)
dev.off()

prev.PrimaryVsMet <- rbind(prev.TCGA, prev.OAK.metastatic)
rownames(prev) <- c("TCGA", "OAKmet")
fisher.test(prev.PrimaryVsMet)
