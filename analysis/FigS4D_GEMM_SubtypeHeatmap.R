library(ComplexHeatmap)

load("../data/GEMM_RNAseq.RData")
data <- log2(exprs(eset)+1)
rownames(data) <- fData(eset)$symbol
rowannot <- fData(eset)$Subtype
rowannot[rowannot %in% "MUC"] <- "Murine subtype 1"
rowannot[rowannot %in% "PRO"] <- "Murine subtype 2"
subtype <- eset$Subtype
colannot <- data.frame(Subtype = subtype)

y <- t(scale(t(data)))
t <- 1.5
y[y>t] <- t
y[y<(-t)] <- (-t)

ht.S1 <- Heatmap(y[, subtype %in% "MUC"], name="Z-score", cluster_rows=TRUE, cluster_columns=TRUE)
ht.S2 <- Heatmap(y[, subtype %in% "PRO"], name="Z-score", cluster_rows=TRUE, cluster_columns=TRUE)
ht.NA <- Heatmap(y[, subtype %in% "unclassified"], name="Z-score", cluster_rows=TRUE, cluster_columns=TRUE)
col.order <- c(which(subtype %in% "MUC")[column_order(ht.S1)],
						which(subtype %in% "PRO")[column_order(ht.S2)],
						which(subtype %in% "unclassified")[column_order(ht.NA)])

hr <- rowAnnotation(Subtype=rowannot,
								col=list(Subtype=c("Murine subtype 1"="cyan3", "Murine subtype 2"="mediumpurple2")),
								show_annotation_name=FALSE, show_legend=FALSE)
hc <- HeatmapAnnotation(df=colannot[col.order,,drop=FALSE],
										col=list(Subtype=c("MUC"="cyan3", "PRO"="mediumpurple2", "unclassified"="gray")),
										show_annotation_name=TRUE,
										annotation_name_offset = unit(2, "mm"))
ht1 <- Heatmap(y[,col.order], name="Z-score", cluster_rows=TRUE, split=rowannot, cluster_columns=FALSE,
						show_row_names=TRUE, row_names_gp=gpar(fontsize=8),
						top_annotation=hc, show_column_names=FALSE, show_column_dend=FALSE, show_row_dend=FALSE,
						heatmap_legend_param=list(legend_direction="horizontal", legend_width=unit(4, "cm"),
						title_position="lefttop"))
pdf("../figures/FigS4_GEMM_SubtypeHeatmap.pdf", width=6, height=6)
draw(ht1+hr, heatmap_legend_side="bottom")
dev.off()
