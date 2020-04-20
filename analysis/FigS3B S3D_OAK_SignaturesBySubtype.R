study <- "OAK"
data1 <- cbind(read.table("../data/OAK_KRASmut_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE), KRAS="MUT")
data2 <- cbind(read.table("../data/OAK_KRASwt_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE), KRAS="WT")
data1 <- data1[,-which(colnames(data1) %in% "NMF")]
data <- rbind(data1, data2)

data <- data[!data$Subtype %in% "unclassified", ]
data$Subtype <- factor(data$Subtype, levels=c("MUC", "PRO", "MES"))
nbTumors <- nrow(data)

colors <- rep("cyan3", nbTumors)
colors[data$Subtype %in% "PRO"] <- "mediumpurple2"
colors[data$Subtype %in% "MES"] <- "green3"
pchs <- rep(19, nbTumors) # circle
pchs[data$KRAS %in% "WT"] <- 17 # triangle

pdf("../figures/FigS3_OAK_STK11signatureBySubtype.pdf", width=5.5, height=6.5)
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
t.test(data$STK11_deficiency_signature[data$Subtype %in% "PRO" & data$KRAS %in% "MUT"],
			data$STK11_deficiency_signature[data$Subtype %in% c("MUC","MES") & data$KRAS %in% "MUT"],
			alternative = "greater")$p.value
t.test(data$STK11_deficiency_signature[data$Subtype %in% "PRO" & data$KRAS %in% "WT"],
			data$STK11_deficiency_signature[data$Subtype %in% c("MUC","MES") & data$KRAS %in% "WT"],
			alternative = "greater")$p.value

pdf("../figures/FigS3_OAK_MucinoussignatureBySubtype.pdf", width=5.5, height=6.5)
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
t.test(data$Mucinous_differentiation_signature[data$Subtype %in% "MUC" & data$KRAS %in% "MUT"],
			data$Mucinous_differentiation_signature[data$Subtype %in% c("PRO","MES") & data$KRAS %in% "MUT"],
			alternative = "greater")$p.value
t.test(data$Mucinous_differentiation_signature[data$Subtype %in% "MUC" & data$KRAS %in% "WT"],
			data$Mucinous_differentiation_signature[data$Subtype %in% c("PRO","MES") & data$KRAS %in% "WT"],
			alternative = "greater")$p.value
