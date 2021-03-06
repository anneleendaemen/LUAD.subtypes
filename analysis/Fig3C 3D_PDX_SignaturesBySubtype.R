study <- "PDX models"
data <- read.table("../data/PDX_KRASall_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
data <- data[!data$Subtype %in% "unclassified", ]
data$Subtype <- factor(data$Subtype, levels=c("MUC", "PRO", "MES"))
nbTumors <- nrow(data)

colors <- rep("cyan3", nbTumors)
colors[data$Subtype %in% "PRO"] <- "mediumpurple2"
colors[data$Subtype %in% "MES"] <- "green3"
pchs <- rep(19, nbTumors)
pchs[data$KRAS %in% "WT"] <- 17
pchs[is.na(data$KRAS)] <- 18


pdf("../figures/Fig4_PDX_STK11signatureBySubtype.pdf", width=5.5, height=6.5)
par(mar=c(4.5, 4.5, 3, 1))
boxplot(data$STK11_deficiency_signature ~ data$Subtype, outline=FALSE,
				ylab="STK11 signature", xlab="Subtype", main=paste0(study," (n=",nbTumors,")"),
				cex.lab=1.8, cex.axis=1.5, cex.main=1.8, ylim=range(data$STK11_deficiency_signature))
points(jitter(as.numeric(data$Subtype)), data$STK11_deficiency_signature,
					col=colors, pch=pchs, cex=1.5)
pval <- signif(kruskal.test(data$STK11_deficiency_signature ~ data$Subtype)$p.value, 2)
text(x=3, y=max(data$STK11_deficiency_signature), labels=paste0("p=",pval), cex=1.5)
dev.off()

t.test(data$STK11_deficiency_signature[data$Subtype %in% "PRO"],
			data$STK11_deficiency_signature[data$Subtype %in% c("MUC","MES")],
			alternative = "greater")$p.value

pdf("../figures/Fig3_PDX_MucinoussignatureBySubtype.pdf", width=5.5, height=6.5)
par(mar=c(4.5, 4.5, 3, 1))
boxplot(data$Mucinous_differentiation_signature ~ data$Subtype, outline=FALSE,
				ylab="Mucinous differentiation signature", xlab="Subtype", main=paste0(study," (n=",nbTumors,")"),
				cex.lab=1.8, cex.axis=1.5, cex.main=1.8, ylim=range(data$Mucinous_differentiation_signature))
points(jitter(as.numeric(data$Subtype)), data$Mucinous_differentiation_signature,
					col=colors, pch=pchs, cex=1.5)
pval <- signif(kruskal.test(data$Mucinous_differentiation_signature ~ data$Subtype)$p.value, 2)
text(x=3, y=max(data$Mucinous_differentiation_signature), labels=paste0("p=",pval), cex=1.5)
dev.off()

t.test(data$Mucinous_differentiation_signature[data$Subtype %in% "MUC"],
			data$Mucinous_differentiation_signature[data$Subtype %in% c("PRO","MES")],
			alternative = "greater")$p.value
