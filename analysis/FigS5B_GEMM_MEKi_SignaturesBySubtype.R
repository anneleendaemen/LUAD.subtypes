study <- "GEMM MEKi"
data <- read.table("../data/GEMM_MEKi_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
data <- data[!data$Subtype %in% "unclassified", ]
data$Subtype <- factor(data$Subtype, levels=c("MUC", "PRO"))
nbTumors <- nrow(data)

colors <- rep("cyan3", nbTumors)
colors[data$Subtype %in% "PRO"] <- "mediumpurple2"
pchs <- rep(19, nbTumors)

pdf("../figures/FigS5_GEMM_MEKi_STK11signatureBySubtype.pdf", width=5.5, height=6.5)
par(mar=c(4.5, 4.5, 3, 1))
boxplot(data$Stk11_deficiency_signature ~ data$Subtype, outline=FALSE,
				ylab="STK11 signature", xlab="Subtype", main=paste0(study," (n=",nbTumors,")"),
				cex.lab=1.8, cex.axis=1.5, cex.main=1.8, ylim=range(data$Stk11_deficiency_signature))
points(jitter(as.numeric(data$Subtype)), data$Stk11_deficiency_signature,
					col=colors, pch=pchs, cex=1.5)
pval <- signif(t.test(data$Stk11_deficiency_signature ~ data$Subtype)$p.value, 2)
text(x=1.5, y=max(data$Stk11_deficiency_signature), labels=paste0("p=",pval), cex=1.5)
dev.off()
