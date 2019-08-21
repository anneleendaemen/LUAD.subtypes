data.mut <- read.table("../data/TCGA_KRASmut_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
data.wt <- read.table("../data/TCGA_KRASwt_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
data.mut <- data.mut[, -match("NMF", colnames(data.mut))]
data.mut <- cbind(data.mut, KRAS="MUT")
data.wt <- cbind(data.wt, KRAS="WT")
data <- rbind(data.mut, data.wt)
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
prev.TCGA <- c(prev.STK11, prev.CDKN2A, prev.TP53)

data.mut <- read.table("../data/OAK_KRASmut_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
data.wt <- read.table("../data/OAK_KRASwt_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
data.mut <- data.mut[, -match("NMF", colnames(data.mut))]
data.mut <- cbind(data.mut, KRAS="MUT")
data.wt <- cbind(data.wt, KRAS="WT")
data <- rbind(data.mut, data.wt)
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
prev.OAK <- c(prev.STK11, prev.CDKN2A, prev.TP53)

prev <- rbind(prev.TCGA[1:3], prev.OAK[1:3],
					prev.TCGA[4:6], prev.OAK[4:6],
					prev.TCGA[7:9], prev.OAK[7:9])
rownames(prev) <- rep(c("TCGA", "OAK"),3)
prev <- t(prev)

pdf("../figures/Fig2_TCGAOAK_GeneticEnrichment.pdf", width=6.5, height=6.5)
par(mar=c(4.5, 6.5, 3, 1), xpd=TRUE)
bp <- barplot(prev, xlab="", ylab="Prevalence of genetic\nevent by subtype", col=c("green3", "mediumpurple2", "cyan3"), 
			cex.axis=1.5, cex.lab=1.8, cex.names=1.25, ylim=c(0,170), border=FALSE, space=c(rep(0.1,2),0.5,0.1,0.5,0.1))

lines(x=c(bp[1]-0.5,bp[2]+0.5), y=rep(87,2), lwd=2)
text(x=(bp[1]+bp[2])/2, y=95, labels="STK11", cex=1.8)
lines(x=c(bp[3]-0.5,bp[4]+0.5), y=rep(75,2), lwd=2)
text(x=(bp[3]+bp[4])/2, y=83, labels="CDKN2A", cex=1.8)
lines(x=c(bp[5]-0.5,bp[6]+0.5), y=rep(170,2), lwd=2)
text(x=(bp[5]+bp[6])/2, y=178, labels="TP53", cex=1.8)
text(x=-1, y=-9.5, labels="Cohort", cex=1.3)
dev.off()
