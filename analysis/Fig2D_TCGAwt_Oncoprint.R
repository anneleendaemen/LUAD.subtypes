library(ComplexHeatmap)

data <- read.table("../data/TCGA_KRASwt_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
samples.1 <- substr(as.character(data$SAMID[data$Subtype %in% "MUC"]), 1, 12) 
samples.2 <- substr(as.character(data$SAMID[data$Subtype %in% "PRO"]), 1, 12) 
samples.3 <- substr(as.character(data$SAMID[data$Subtype %in% "MES"]), 1, 12) 
samples.uncl <- substr(as.character(data$SAMID[data$Subtype %in% "unclassified"]), 1, 12) 
samples <- c(samples.1, samples.2, samples.3, samples.uncl)
nbTumors <- length(samples)

# Genomic alterations obtained from cBioPortal
alt <- read.table(file= "../data/TCGAwt_cBioPortal_EGFR PI3K MAPK.tsv", header = T, sep="\t")
colnames(alt) <- gsub("[.]", "-", colnames(alt))
nbGenes <- length(unique(alt$track_name))

CNA <- alt[alt$track_type %in% "CNA", ]
rownames(CNA) <- CNA$track_name
MUT <- alt[alt$track_type %in% "MUTATIONS", ]
rownames(MUT) <- MUT$track_name
FUS <- alt[alt$track_type %in% "FUSION", ]
rownames(FUS) <- FUS$track_name

mat <- matrix(data = NA, nrow = nbGenes, ncol = ncol(alt)-2)
colnames(mat) <- colnames(CNA)[-c(1:2)]
rownames(mat) <- rownames(CNA)

for(gene in rownames(mat)){
  amp <- colnames(CNA[, which(CNA[gene, ] == "Amplification")])
  del <- colnames(CNA[, which(CNA[gene, ] == "Deep Deletion")])
  mut <- colnames(MUT[, which(MUT[gene, ] == "Missense Mutation (putative driver)" |
                                MUT[gene, ] == "Inframe Mutation (putative driver)" |
                                MUT[gene, ] == "Truncating mutation (putative driver)" |
                                MUT[gene, ] == "Truncating mutation (putative passenger)" |
                                MUT[gene, ] == "Missense Mutation (putative passenger)")])
  fus <- colnames(FUS[, which(FUS[gene, ] == "Fusion")])
  
  amp_mut <- intersect(amp, mut)
  del_mut <- intersect(del, mut)
  amp_fus <- intersect(amp, fus)
  del_fus <- intersect(del, fus)
  
  mat[gene, amp] <- "AMP"
  mat[gene, del] <- "HOMDEL"
  mat[gene, mut] <- "MUT"
  mat[gene, fus] <- "FUSION"  
  mat[gene, amp_mut] <- "AMP;MUT"
  mat[gene, del_mut] <- "HOMDEL;MUT"
  mat[gene, amp_fus] <- "AMP;FUSION"
  mat[gene, del_fus] <- "HOMDEL;FUSION"  
}

mat[is.na(mat)] = ""

mat.S1 <- mat[, samples.1]
mat.S2 <- mat[, samples.2]
mat.S3 <- mat[, samples.3]
mat.Suncl <- mat[, samples.uncl]
percAltered.S1 <- (ncol(mat.S1)-length(which(apply(mat.S1, 2, function(x) length(which(!x=="")))==0)))/ncol(mat.S1)
percAltered.S2 <- (ncol(mat.S2)-length(which(apply(mat.S2, 2, function(x) length(which(!x=="")))==0)))/ncol(mat.S2)
percAltered.S3 <- (ncol(mat.S3)-length(which(apply(mat.S3, 2, function(x) length(which(!x=="")))==0)))/ncol(mat.S3)
percAltered.Suncl <- (ncol(mat.Suncl)-length(which(apply(mat.Suncl, 2, function(x) length(which(!x=="")))==0)))/ncol(mat.Suncl)

mat <- mat[, samples]
percAltered <- (ncol(mat)-length(which(apply(mat, 2, function(x) length(which(!x=="")))==0)))/ncol(mat)

# Oncoprint visualization
col = c("HOMDEL" = "blue", "AMP" = "red", "MUT" = "#008000", "FUSION" = "purple")
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  HOMDEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.9, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["HOMDEL"], col = NA))
  },
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.9, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["AMP"], col = NA))
  },
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.9, "mm"), h*0.33, 
              gp = gpar(fill = col["MUT"], col = NA))
  },
  FUSION = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.9, "mm"), h*0.5, 
              gp = gpar(fill = col["FUSION"], col = NA))
  }  
)

subt_col = c(rep("MUC", length(samples.1)),
             rep("PRO", length(samples.2)),
             rep("MES", length(samples.3)),
             rep("unclassified", length(samples.uncl)))

column_title = paste0("TCGA - KRAS wildtype (n=",nbTumors,")")
heatmap_legend_param = list(title = "Alterations", at = c("HOMDEL", "AMP", "MUT", "FUSION"), 
                            labels = c("Deep deletion", "Amplification", "Mutation", "Fusion"),
                            labels_gp = gpar(fontsize=12))
top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
                                   subtype = subt_col,
                                   col = list(subtype = c("MUC" = "cyan3", "PRO" = "mediumpurple2", "MES" = "green3", "unclassified" = "gray")),
                                   annotation_legend_param = list(subtype =
                                                                    list(title = "Subtype",
                                                                         labels = c("MUC", "PRO", "MES", "unclassified"),
                                                                  				labels_gp = gpar(fontsize = 12))))

a = oncoPrint(mat,
              alter_fun = alter_fun, col = col, 
              column_title = column_title, heatmap_legend_param = heatmap_legend_param,
              top_annotation = top_annotation,
              get_type = function(x) strsplit(x, ";")[[1]])
col_order <- column_order(a)

S1_order <- col_order[col_order %in% which(subt_col %in% "MUC")]
S2_order <- col_order[col_order %in% which(subt_col %in% "PRO")]
S3_order <- col_order[col_order %in% which(subt_col %in% "MES")]
Suncl_order <- col_order[col_order %in% which(subt_col %in% "unclassified")]

column_order <- c(S1_order, S2_order, S3_order, Suncl_order)

pdf("../figures/Fig2D_TCGAwt_Oncoprint.pdf", width=9.5, height=5.5, onefile=FALSE)
oncoPrint(mat,
          alter_fun = alter_fun, col = col, 
          column_title = column_title, heatmap_legend_param = heatmap_legend_param,
          top_annotation = top_annotation,
          column_order = column_order,
          get_type = function(x) strsplit(x, ";")[[1]])
dev.off()
