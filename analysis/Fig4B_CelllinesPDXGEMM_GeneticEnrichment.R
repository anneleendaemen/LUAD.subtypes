data <- read.table("../data/CellLines_KRASall_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
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
prev.Celllines <- c(prev.STK11, prev.CDKN2A, prev.TP53)

data <- read.table("../data/PDX_KRASall_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
data <- data[!data$Subtype %in% "unclassified",]
data$Subtype <- factor(data$Subtype, levels=c("MES", "PRO", "MUC"))
data$STK11 <- as.character(data$STK11)
data$STK11[!data$STK11 %in% "WT" & !is.na(data$STK11)] <- "ALT"
prev <- table(data$Subtype, data$STK11)
fisher.test(rbind(prev["PRO", ], colSums(prev[c("MUC", "MES"), ])), alternative="greater")$p.value
prev.STK11 <- prop.table(prev, 1)[,"ALT"]*100
data$CDKN2A <- as.character(data$CDKN2A)
data$CDKN2A[!data$CDKN2A %in% "WT" & !is.na(data$CDKN2A)] <- "ALT"
prev <- table(data$Subtype, data$CDKN2A)
fisher.test(rbind(prev["MUC", ], colSums(prev[c("PRO", "MES"), ])), alternative="greater")$p.value
prev.CDKN2A <- prop.table(prev, 1)[,"ALT"]*100
data$TP53 <- as.character(data$TP53)
data$TP53[!data$TP53 %in% "WT" & !is.na(data$TP53)] <- "ALT"
prev <- table(data$Subtype, data$TP53)
fisher.test(rbind(prev["MES", ], colSums(prev[c("MUC", "PRO"), ])), alternative="greater")$p.value
prev.TP53 <- prop.table(prev, 1)[,"ALT"]*100
prev.PDX <- c(prev.STK11, prev.CDKN2A, prev.TP53)

data <- read.table("../data/GEMM_KRASmut_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
data <- data[!data$Subtype %in% "unclassified",]
data$Subtype <- factor(data$Subtype, levels=c("PRO", "MUC"))
data$STK11 <- data$Stk11_loss
data$STK11[data$Stk11_loss] <- "ALT"
data$STK11[!data$Stk11_loss] <- "WT"
prev <- table(data$Subtype, data$STK11)
fisher.test(rbind(prev["PRO", ], prev["MUC", ]))$p.value
prev.STK11 <- prop.table(prev, 1)[,"ALT"]*100
prev.STK11["MES"] <- 0
prev.GEMM <- c(prev.STK11)

prev <- rbind(prev.Celllines[1:3], prev.PDX[1:3], prev.GEMM[1:3],
					prev.Celllines[4:6], prev.PDX[4:6], 
					prev.Celllines[7:9], prev.PDX[7:9])
rownames(prev) <- c("Cell\nlines", "PDX", "GEMM", "Cell\nlines", "PDX", "Cell\nlines", "PDX")
prev <- t(prev)

pdf("../figures/Fig4_CelllinesPDXGEMM_GeneticEnrichment.pdf", width=8, height=6.5)
par(mar=c(4.5, 6.5, 3, 1), xpd=TRUE)
bp <- barplot(prev, xlab="", ylab="Prevalence of genetic\nevent by subtype", col=c("green3", "mediumpurple2", "cyan3"), 
			cex.axis=1.5, cex.lab=1.8, ylim=c(0,250), border=FALSE, space=c(rep(0.1,3),0.5, 0.1,0.5, 0.1), xaxt="n")

lines(x=c(bp[1]-0.5,bp[3]+0.5), y=rep(58,2), lwd=2)
text(x=(bp[1]+bp[3])/2, y=67, labels="STK11", cex=1.8)
lines(x=c(bp[4]-0.5,bp[5]+0.5), y=rep(130,2), lwd=2)
text(x=(bp[4]+bp[5])/2, y=139, labels="CDKN2A", cex=1.8)
lines(x=c(bp[6]-0.5,bp[7]+0.5), y=rep(222,2), lwd=2)
text(x=(bp[6]+bp[7])/2, y=231, labels="TP53", cex=1.8)
axis(1, at=bp, pos=-8, labels=colnames(prev), lwd=0, cex.axis=1.3)
#text(x=-1, y=-15, labels="Cohort", cex=1.3)
dev.off()
