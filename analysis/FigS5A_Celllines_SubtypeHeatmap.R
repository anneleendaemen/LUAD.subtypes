library(ComplexHeatmap)
library(Biobase)

load("../data/Celllines_RNAseq_ClassifierGenes.RData")
data <- log2(exprs(eset)+1)
rowannot <- fData(eset)$Subtype
rowannot <- factor(rowannot, levels=c("MUC up", "MUC dn", "PRO up", "PRO dn", "MES up", "MES dn"))
subtype <- eset$Subtype
colannot <- data.frame(Subtype = subtype)

y <- t(scale(t(data)))
t <- 1
y[y>t] <- t
y[y<(-t)] <- (-t)

ht.S1 <- Heatmap(y[, subtype %in% "MUC"], name="Z-score", cluster_rows=TRUE, cluster_columns=TRUE)
ht.S2 <- Heatmap(y[, subtype %in% "PRO"], name="Z-score", cluster_rows=TRUE, cluster_columns=TRUE)
ht.S3 <- Heatmap(y[, subtype %in% "MES"], name="Z-score", cluster_rows=TRUE, cluster_columns=TRUE)
ht.NA <- Heatmap(y[, subtype %in% "unclassified"], name="Z-score", cluster_rows=TRUE, cluster_columns=TRUE)
col.order <- c(which(subtype %in% "MUC")[column_order(ht.S1)],
						which(subtype %in% "PRO")[column_order(ht.S2)],
						which(subtype %in% "MES")[column_order(ht.S3)],
						which(subtype %in% "unclassified")[column_order(ht.NA)])

hr <- rowAnnotation(Subtype=rowannot,
								col=list(Subtype=c("MUC up"="cyan3", "PRO up"="mediumpurple2", "MES up"="green3",
																"MUC dn"="cyan3", "PRO dn"="mediumpurple2", "MES dn"="green3")),
								show_annotation_name = c(Subtype=FALSE), show_legend = FALSE)
hc <- HeatmapAnnotation(df=colannot[col.order,,drop=FALSE],
										col=list(Subtype=c("MUC"="cyan3", "PRO"="mediumpurple2", "MES"="green3", "unclassified"="gray")),
										show_annotation_name=TRUE,
										annotation_name_offset = unit(2, "mm"),
										annotation_legend_param = list(Subtype =
                                                             list(title = "Subtype", title_gp = gpar(fontsize = 14),
                                                                  labels = c("MUC", "PRO", "MES", "unclassified"),
                                                                  labels_gp = gpar(fontsize = 13))))
ht1 <- Heatmap(y[,col.order], column_title = paste0("Cell lines (n=",ncol(eset),")"),
						column_title_gp = gpar(fontsize = 16, fontface = "bold"),
						name="Z-score", cluster_rows=FALSE, cluster_columns=FALSE, row_split=rowannot, column_split=colannot$Subtype[col.order],
						row_title=c("up", "dn", "up", "dn", "up", "dn"), row_title_rot = 0, row_gap=unit(1.5, "mm"), column_gap=unit(1.5,"mm"),
						border=TRUE, 
						show_row_names=FALSE, 
						top_annotation=hc, show_column_names=FALSE, show_column_dend=FALSE, show_row_dend=FALSE,
						heatmap_legend_param=list(legend_direction="horizontal", legend_width=unit(4, "cm"),
						title_position="lefttop"))
pdf("../figures/FigS5_Celllines_SubtypeHeatmap.pdf", width=6, height=6)
draw(ht1+hr, heatmap_legend_side="bottom")
dev.off()
