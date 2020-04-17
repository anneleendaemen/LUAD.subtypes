Study <- "GEMM KP"
data <- read.table("../data/GEMM_KRASmut_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
data <- data[!data$Subtype %in% "unclassified", ]
data <- data[!data$Subtype_13of18 %in% "unclassified", ]
nbTumors <- nrow(data)
data$Subtype <- factor(data$Subtype, levels=c("MUC", "PRO"))
data$Subtype_13of18 <- factor(data$Subtype_13of18, levels=c("KC", "KL", "KP"))
counts <- as.matrix(table(data$Subtype, data$Subtype_13of18))
rownames(counts) <- c("Murine subtype 1", "Murine subtype 2")
colnames(counts) <- c("MUC", "PRO", "MES")

pdf("../figures/FigS5_GEMM_SubtypeConcordanceWithHuman.pdf", width=5.5, height=6.5)
par(mar=c(4.5, 4.5, 3, 1))
barplot(counts, xlab="Human subtypes", ylab="Number of tumors", col=c("cyan3", "mediumpurple2"), 
			main=paste0(Study, " (n=", nbTumors, ")"), cex.axis=1.5, cex.lab=1.8, cex.names=1.5, cex.main=1.8, ylim=c(0,30))
legend("topleft", legend=rownames(counts), fill=c("cyan3", "mediumpurple2"), cex=1.6, bty="n")
dev.off()

fisher.test(counts[,c("MUC", "PRO")])
