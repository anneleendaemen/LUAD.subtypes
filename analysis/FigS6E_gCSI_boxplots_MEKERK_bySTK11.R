data <- read.table("../data/CellLineDrugResponse_SubtypeAssociation.tsv", sep="\t", header=TRUE, check.names=FALSE)
gCSI <- read.table("../data/CellLineDrugResponse_MeanViability.tsv", sep="\t", header=TRUE, check.names=FALSE)

# Associate MEK/ERKi response with subtype
target <- c("MEK2, MEK1", "ERK")
compounds <- as.character(data$`gCGC ID`[which(data$TARGET_NAME %in% target)])
targetsID <- as.character(data$TARGET_NAME[which(data$TARGET_NAME %in% target)])
data.gCSI <- gCSI[, compounds]
rownames(data.gCSI) <- gCSI[, "Cell line"]
colnames(data.gCSI) <- paste0(compounds, " (", targetsID, ")")
nbTumors <- nrow(data.gCSI)

STK11 <- factor(as.character(gCSI$STK11), levels=c("WT", "ALT"))
pchs <- rep(19, nbTumors)
pchs[!gCSI$KRAS] <- 17

pvalues <- c() # in all 89 cell lines
for (i in colnames(data.gCSI)) pvalues <- c(pvalues, signif(t.test(data.gCSI[,i] ~ STK11, alternative="greater")$p.value,2))
names(pvalues) <- colnames(data.gCSI)

i <- "gCGC00062 (MEK2, MEK1)"
pdf("../figures/FigS6_CellLines_CobimetinibVsSTK11.pdf", width=5.5, height=6.5)
par(mar=c(4.5,4.5,3,1))
pval <- signif(t.test(data.gCSI[,i] ~ STK11)$p.value,2)
boxplot(data.gCSI[,i] ~ STK11, ylab="Mean viability, Cobimetinib", xlab="Genetic STK11 status",
			main=paste0("Cell lines (n=",nbTumors,")"),
			cex.lab=1.8, cex.axis=1.5, cex.main=1.8, outline=FALSE,
			ylim=range(data.gCSI[,i], na.rm=TRUE))
points(jitter(as.numeric(STK11)), data.gCSI[,i], col="gray", pch=pchs, cex=1.5)
text(x=2, y=max(data.gCSI[,i], na.rm=TRUE), labels=paste0("p=",pval), cex=1.5)
dev.off()

indices.rm <- which(gCSI$Subtype %in% "unclassified")
data.gCSI <- data.gCSI[-indices.rm, ]
gCSI <- gCSI[-indices.rm, ]
STK11 <- STK11[-indices.rm]

pvalues.STK11_to_subtype <- c() # in subset of 63 cell lines with subtype 1/2/3
pvalues.subtype_to_STK11 <- c()
for (i in colnames(data.gCSI)) {
	model1a <- lm(data.gCSI[,i] ~ gCSI$Subtype)
	model1b <- lm(data.gCSI[,i] ~ STK11)
	model2 <- lm(data.gCSI[,i] ~ gCSI$Subtype + STK11)
	pvalues.STK11_to_subtype <- c(pvalues.STK11_to_subtype, anova(model1a, model2)$"Pr(>F)"[2])
	pvalues.subtype_to_STK11 <- c(pvalues.subtype_to_STK11, anova(model1b, model2)$"Pr(>F)"[2])
}
names(pvalues.STK11_to_subtype) <- names(pvalues.subtype_to_STK11) <- colnames(data.gCSI)
sort(pvalues.subtype_to_STK11[grep("MEK", names(pvalues.subtype_to_STK11))])
sort(pvalues.STK11_to_subtype[grep("MEK", names(pvalues.STK11_to_subtype))])
