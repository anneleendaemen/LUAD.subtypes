study <- "GEMM KP"
data <- read.table("../data/GEMM_KRASmut_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
data <- data[!data$Subtype %in% "unclassified", ]
data$Subtype <- factor(data$Subtype, levels=c("MUC", "PRO"))
nbTumors <- nrow(data)

data.gastric <- data[, c(grep("Krt20", colnames(data)), grep("Muc5ac", colnames(data)),
									grep("Muc5b", colnames(data)), grep("Pdx1", colnames(data)))]
data.gastric <- scale(data.gastric, center=TRUE, scale=TRUE)
data$Mucinous_differentiation_signature <- rowMeans(data.gastric)

colors <- rep("cyan3", nbTumors)
colors[data$Subtype %in% "PRO"] <- "mediumpurple2"
pchs <- rep(19, nbTumors)

pdf("../figures/Fig3_GEMM_STK11signatureBySubtype.pdf", width=5.5, height=6.5)
par(mar=c(4.5, 4.5, 3, 1))
boxplot(data$STK11_deficiency_signature ~ data$Subtype, outline=FALSE,
				ylab="STK11 signature", xlab="Subtype", main=paste0(study," (n=",nbTumors,")"),
				cex.lab=1.8, cex.axis=1.5, cex.main=1.8, ylim=range(data$STK11_deficiency_signature))
points(jitter(as.numeric(data$Subtype)), data$STK11_deficiency_signature,
					col=colors, pch=pchs, cex=1.5)
pval <- signif(t.test(data$STK11_deficiency_signature ~ data$Subtype)$p.value, 2)
text(x=1.5, y=max(data$STK11_deficiency_signature), labels=paste0("p=",pval), cex=1.5)
dev.off()

pdf("../figures/Fig3_GEMM_MucinoussignatureBySubtype.pdf", width=5.5, height=6.5)
par(mar=c(4.5, 4.5, 3, 1))
boxplot(data$Mucinous_differentiation_signature ~ data$Subtype, outline=FALSE,
				ylab="Mucinous differentiation signature", xlab="Subtype", main=paste0(study," (n=",nbTumors,")"),
				cex.lab=1.8, cex.axis=1.5, cex.main=1.8, ylim=range(data$Mucinous_differentiation_signature))
points(jitter(as.numeric(data$Subtype)), data$Mucinous_differentiation_signature,
					col=colors, pch=pchs, cex=1.5)
pval <- signif(t.test(data$Mucinous_differentiation_signature ~ data$Subtype)$p.value, 2)
text(x=1.5, y=max(data$Mucinous_differentiation_signature), labels=paste0("p=",pval), cex=1.5)
dev.off()
