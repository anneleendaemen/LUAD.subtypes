library(ComplexHeatmap)
library(Biobase)

load("../data/POPLAR_RNAseq_ClassifierGenes.RData")
data <-  exprs(eset)
rowannot <- data.frame(Subtype=fData(eset)$Subtype)
rowannot$Subtype <- factor(rowannot$Subtype, levels=c("MUC up", "MUC dn", "PRO up", "PRO dn", "MES up", "MES dn"))
colannot <- data.frame(Subtype = eset$Subtype)
colannot2 <- data.frame(Patient = as.character(eset$Label))

y <- t(scale(t(data)))
t <- 1.5
y[y>t] <- t
y[y<(-t)] <- (-t)

ht.S1 <- Heatmap(data[which(fData(eset)$Subtype %in% "MUC up"), ], cluster_rows = T)
ht.S2 <- Heatmap(data[which(fData(eset)$Subtype %in% "MUC dn"), ], cluster_rows = T)
ht.S3 <- Heatmap(data[which(fData(eset)$Subtype %in% "PRO up"), ], cluster_rows = T)
ht.S4 <- Heatmap(data[which(fData(eset)$Subtype %in% "PRO dn"), ], cluster_rows = T)
ht.S5 <- Heatmap(data[which(fData(eset)$Subtype %in% "MES up"), ], cluster_rows = T)
ht.S6 <- Heatmap(data[which(fData(eset)$Subtype %in% "MES dn"), ], cluster_rows = T)

r_order <- c(which(fData(eset)$Subtype %in% "MUC up")[unlist(row_order(ht.S1))],
					which(fData(eset)$Subtype %in% "MUC dn")[unlist(row_order(ht.S2))], 
					which(fData(eset)$Subtype %in% "PRO up")[unlist(row_order(ht.S3))],
					which(fData(eset)$Subtype %in% "PRO dn")[unlist(row_order(ht.S4))], 
					which(fData(eset)$Subtype %in% "MES up")[unlist(row_order(ht.S5))],
					which(fData(eset)$Subtype %in% "MES dn")[unlist(row_order(ht.S6))])
rowannot <- rowannot[r_order, , drop=FALSE]

hr <- rowAnnotation(df=rowannot,
								col=list(Subtype=c("MUC up"="cyan3", "PRO up"="mediumpurple2", "MES up"="green3",
																"MUC dn"="cyan3", "PRO dn"="mediumpurple2", "MES dn"="green3")),
								show_annotation_name = FALSE, show_legend = FALSE)
hc <- HeatmapAnnotation(df=colannot,
										col=list(Subtype=c("MUC"="cyan3", "MES"="green3", "unclassified"="gray")),
										show_annotation_name=TRUE,
										annotation_name_offset = unit(2, "mm"),
										annotation_legend_param = list(Subtype =
                                                             list(title = "Subtype", title_gp = gpar(fontsize = 14),
                                                                  labels = c("MUC", "MES", "unclassified"),
                                                                  labels_gp = gpar(fontsize = 13))))
hc2 <- HeatmapAnnotation(df = colannot2,
                        col = list(Patient = c("P1 - lung" = "aquamarine2", "P2 - mediastinal LN" = "blueviolet",
                                               "P3 - LN" = "brown1", "P3 - neck" = "brown3",
                                               "P4 - ovary" = "chartreuse3", "P5 - supraclavicular LN" = "darkgoldenrod1",
                                               "P6 - lung" = "dodgerblue")),
                        annotation_legend_param = list(Patient = list(title = "Patient", title_gp = gpar(fontsize = 14), 
                                                                      labels_gp = gpar(fontsize = 10))))

ht1 <- Heatmap(y[r_order, ], column_title = paste0("POPLAR (n=",ncol(eset),")"),
						column_title_gp = gpar(fontsize = 16, fontface = "bold"),
						name="Z-score", cluster_rows=FALSE, cluster_columns=FALSE, row_split=rowannot, column_split=substr(colannot2$Patient,1,2),
						row_title=c("up", "dn", "up", "dn", "up", "dn"), row_title_rot = 0, row_gap=unit(1.5, "mm"), column_gap=unit(1.5,"mm"),
						border=TRUE,
						show_row_names=FALSE, show_column_names=FALSE, show_column_dend=FALSE, show_row_dend=FALSE,
						top_annotation=hc, bottom_annotation=hc2,
						heatmap_legend_param=list(legend_direction="horizontal", legend_width=unit(4, "cm"), title_position="lefttop"))

pdf("../figures/Fig3_POPLAR_SubtypeHeatmap.pdf", width=6, height=6)
draw(ht1+hr, heatmap_legend_side="bottom")
dev.off()
