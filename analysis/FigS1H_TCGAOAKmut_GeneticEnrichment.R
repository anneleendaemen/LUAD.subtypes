data <- read.table("../data/TCGA_KRASmut_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
data <- data[!data$Subtype %in% "unclassified",]
data$Subtype <- factor(data$Subtype, levels=c("MES", "PRO", "MUC"))
prev <- table(data$Subtype, data$STK11_mutation_or_homozygous_loss)
fisher.test(rbind(prev["PRO", ], colSums(prev[c("MUC", "MES"), ])), alternative="greater")$p.value
prev.STK11 <- prop.table(prev, 1)[,"ALT"]*100
prev <- table(data$Subtype, data$CDKN2A_mutation_or_homozygous_loss)
fisher.test(rbind(prev["MUC", ], colSums(prev[c("PRO", "MES"), ])), alternative="greater")$p.value
prev.CDKN2A <- prop.table(prev, 1)[,"ALT"]*100
prev <- table(data$Subtype, data$TP53_mutation_or_homozygous_loss)
fisher.test(rbind(prev["MES", ], colSums(prev[c("MUC", "PRO"), ])), alternative="greater")$p.value
prev.TP53 <- prop.table(prev, 1)[,"ALT"]*100
prev <- table(data$Subtype, data$NKX2.1_amplification)
fisher.test(rbind(prev["MUC", ], colSums(prev[c("PRO", "MES"), ])), alternative="less")$p.value
prev.TCGAmut <- c(prev.STK11, prev.CDKN2A, prev.TP53)

data <- read.table("../data/OAK_KRASmut_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
data <- data[!data$Subtype %in% "unclassified",]
data$Subtype <- factor(data$Subtype, levels=c("MES", "PRO", "MUC"))
prev <- table(data$Subtype, data$STK11_mutation_or_homozygous_loss)
fisher.test(rbind(prev["PRO", ], colSums(prev[c("MUC", "MES"), ])), alternative="greater")$p.value
prev.STK11 <- prop.table(prev, 1)[,"ALT"]*100
prev <- table(data$Subtype, data$CDKN2A_mutation_or_homozygous_loss)
fisher.test(rbind(prev["MUC", ], colSums(prev[c("PRO", "MES"), ])), alternative="greater")$p.value
prev.CDKN2A <- prop.table(prev, 1)[,"ALT"]*100
prev <- table(data$Subtype, data$TP53_mutation_or_homozygous_loss)
fisher.test(rbind(prev["MES", ], colSums(prev[c("MUC", "PRO"), ])), alternative="greater")$p.value
prev.TP53 <- prop.table(prev, 1)[,"ALT"]*100
prev <- table(data$Subtype, data$NKX2.1_amplification)
fisher.test(rbind(prev["MUC", ], colSums(prev[c("PRO", "MES"), ])), alternative="less")$p.value
prev.OAKmut <- c(prev.STK11, prev.CDKN2A, prev.TP53)

prev <- rbind(prev.TCGAmut[1:3], prev.OAKmut[1:3],
					prev.TCGAmut[4:6], prev.OAKmut[4:6], 
					prev.TCGAmut[7:9], prev.OAKmut[7:9])
rownames(prev) <- rep(c("TCGA", "OAK"),3)
prev <- t(prev)

pdf("../figures/FigS1_TCGAOAKmut_GeneticEnrichment.pdf", width=6.5, height=6.5)
par(mar=c(4.5, 6.5, 1, 1), xpd=TRUE)
bp <- barplot(prev, xlab="", ylab="Prevalence of genetic\nevent by subtype", col=c("green3", "mediumpurple2", "cyan3"), 
			cex.axis=1.5, cex.lab=1.8, cex.names=1.25, ylim=c(0,170), border=FALSE, space=c(rep(0.1,2),0.5,0.1,0.5,0.1))

lines(x=c(bp[1]-0.5,bp[2]+0.5), y=rep(115,2), lwd=2)
text(x=(bp[1]+bp[2])/2, y=120, labels="STK11", cex=1.8)
lines(x=c(bp[3]-0.5,bp[4]+0.5), y=rep(78,2), lwd=2)
text(x=(bp[3]+bp[4])/2, y=83, labels="CDKN2A", cex=1.8)
lines(x=c(bp[5]-0.5,bp[6]+0.5), y=rep(125,2), lwd=2)
text(x=(bp[5]+bp[6])/2, y=130, labels="TP53", cex=1.8)
text(x=-1, y=-9.5, labels="Cohort", cex=1.3)
dev.off()
