sets <- rev(c("HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
				"HALLMARK_INFLAMMATORY_RESPONSE", 
				"HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", 
				"HALLMARK_APOPTOSIS", "HALLMARK_COMPLEMENT", "HALLMARK_APICAL_JUNCTION", "HALLMARK_APICAL_SURFACE",
				"HALLMARK_KRAS_SIGNALING_UP",
				"HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_E2F_TARGETS",
				"HALLMARK_G2M_CHECKPOINT", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2",
				"BYERS_MESENCHYMAL_SIGNATURE", "ROKAVEC_MESENCHYMAL_SIGNATURE"))

metadata <- read.table("../data/TCGA_KRASmut_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
metadata <- metadata[! metadata$Subtype %in% "unclassified", ]
metadata$Subtype <- factor(metadata$Subtype, levels=c("MUC", "PRO", "MES"))

scores <- subtypes <- levels <- pvalues <- c()
for (set in sets) {
	scores <- c(scores, metadata[,set], 0)
	subtypes <- c(subtypes, paste0(set,metadata$Subtype), paste0(set,"4"))
	levels <- c(levels, paste0(set,"MUC"), paste0(set, "PRO"), paste0(set, "MES"), paste0(set, "4"))
	pvalues <- c(pvalues, kruskal.test(metadata[,set] ~ metadata$Subtype)$p.value)
}
names(pvalues) <- sets
FDR <- p.adjust(pvalues, method="BH")
asterisks <- rep("NS", length(sets))
asterisks[FDR<0.05] <- "*"
asterisks[FDR<0.01] <- "**"
asterisks[FDR<0.001] <- "***"

subtypes <- gsub("HALLMARK_", "", subtypes)
subtypes <- gsub("_SIGNATURE", "", subtypes)
levels <- gsub("HALLMARK_", "", levels)
levels <- gsub("_SIGNATURE", "", levels)
sets <- gsub("HALLMARK_", "", sets)
sets <- gsub("_SIGNATURE", "", sets)

subtypes <- factor(subtypes, levels=levels)
colors <- rep(c("cyan3", "mediumpurple2", "green3", "white"), times=length(sets))

pdf("../figures/FigS2_TCGAmut_HallmarkSignatureBarplots.pdf", width=6.5, height=7.5)
par(mar=c(4.5,12.5,2.5,1))
boxplot(scores ~ subtypes, horizontal=TRUE, col=colors, outline=FALSE,
				border=c(rep("black",3),"white"), 
				cex.main = 1.8, cex.lab = 1.6, cex.axis = 1.3,
				xlab="Signature scores", yaxt="n", main="TCGA - KRAS mutant")
axis(side=2, at=seq(from=2, to=4*length(sets), by=4), labels=sets, las=2, cex.axis=0.7)
text(x=1.3, y=seq(from=2, to=4*length(sets), by=4)+0.3, labels=asterisks, cex=1.2, col="darkblue")
dev.off()
