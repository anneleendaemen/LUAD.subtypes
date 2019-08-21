data1 <- cbind(read.table("../data/TCGA_KRASmut_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE), KRAS="MUT")
data2 <- cbind(read.table("../data/TCGA_KRASwt_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE), KRAS="WT")
data1 <- data1[,-which(colnames(data1) %in% "NMF")]
data <- rbind(data1, data2)

data$Subtype <- factor(data$Subtype, levels=c("MES", "PRO", "MUC", "unclassified"))
data <- data[-which(is.na(data$Expression_Subtype)), ]
data$Expression_Subtype <- as.character(data$Expression_Subtype)
data$Expression_Subtype[data$Expression_Subtype %in% "prox.-inflam"] <- "PI"
data$Expression_Subtype[data$Expression_Subtype %in% "prox.-prolif."] <- "PP"
data$Expression_Subtype <- factor(data$Expression_Subtype, levels=c("PI", "PP", "TRU"))
nbTumors <- nrow(data)

prev <- table(data$Subtype, data$Expression_Subtype)
prev2 <- prop.table(prev, 2)*100

pdf("../figures/FigS2_TCGA_SubtypeConcordance.pdf", width=4, height=6.5)
par(mar=c(4.5, 4.5, 3, 1))
barplot(prev2, xlab="", ylab="Tumor prevalence", col=c("green3", "mediumpurple2", "cyan3", "gray"), 
			cex.axis=1.5, cex.lab=1.8, cex.names=1.5, ylim=c(0,100), border=FALSE,
			main=paste0("TCGA (n=", nbTumors, ")"), cex.main=1.8)
dev.off()
