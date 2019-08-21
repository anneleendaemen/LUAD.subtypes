data <- read.table("../data/CellLineDrugResponse_SubtypeAssociation.tsv", sep="\t", header=TRUE, check.names=FALSE)
gCSI <- read.table("../data/CellLineDrugResponse_MeanViability.tsv", sep="\t", header=TRUE, check.names=FALSE)
gCSI <- gCSI[!gCSI$Subtype %in% "unclassified", ]

# Associate MEK/ERKi response with subtype
target <- c("MEK2, MEK1", "ERK")
compounds <- as.character(data$`gCGC ID`[which(data$TARGET_NAME %in% target)])
targetsID <- as.character(data$TARGET_NAME[which(data$TARGET_NAME %in% target)])
data.gCSI <- gCSI[, compounds]
rownames(data.gCSI) <- gCSI[, "Cell line"]
colnames(data.gCSI) <- paste0(compounds, " (", targetsID, ")")
nbTumors <- nrow(data.gCSI)

subtypes <- factor(gCSI$Subtype, levels=c("MUC", "PRO", "MES"))
colors <- rep("cyan3", nbTumors)
colors[subtypes %in% "PRO"] <- "mediumpurple2"
colors[subtypes %in% "MES"] <- "green3"
pchs <- rep(19, nbTumors)
pchs[!gCSI$KRAS] <- 17

i <- "gCGC00062 (MEK2, MEK1)"
pdf("../figures/Fig5_CellLines_CobimetinibVsSubtype.pdf", width=5.5, height=6.5)
par(mar=c(4.5,4.5,3,1))
pval <- signif(kruskal.test(data.gCSI[,i] ~ subtypes)$p.value,2)
boxplot(data.gCSI[,i] ~ subtypes, ylab="Mean viability, Cobimetinib", xlab="Subtype",
			main=paste0("Cell lines (n=",nbTumors,")"),
			cex.lab=1.8, cex.axis=1.5, cex.main=1.8, outline=FALSE,
			ylim=range(data.gCSI[,i], na.rm=TRUE))
points(jitter(as.numeric(subtypes)), data.gCSI[,i], col=colors, pch=pchs, cex=1.5)
text(x=3, y=max(data.gCSI[,i], na.rm=TRUE), labels=paste0("p=",pval), cex=1.5)
dev.off()

subtypes <- factor(gCSI$Subtype, levels=c("MUC", "PRO", "MES"))
MAPK <- gCSI$MAPK

rsq <- ANOVA.pvalues <- c()
for (i in c(1:ncol(data.gCSI))) {
	model1 <- lm(data.gCSI[,i] ~ subtypes)
	model2 <- lm(data.gCSI[,i] ~ MAPK)
	model3 <- lm(data.gCSI[,i] ~ subtypes + MAPK)
	value <- c(summary(model1)$r.squared, summary(model2)$r.squared, summary(model3)$r.squared)
	rsq <- cbind(rsq, value)
	ANOVA.pvalues <- c(ANOVA.pvalues, anova(model2, model3)$"Pr(>F)"[2])
}
colnames(rsq) <- colnames(data.gCSI)
rownames(rsq) <- c("subtype", "MAPK", "subtype + MAPK")
names(ANOVA.pvalues) <- colnames(data.gCSI)
