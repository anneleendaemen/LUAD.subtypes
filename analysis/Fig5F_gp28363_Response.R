data <- read.table("../data/gp28363_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
data <- data[!data$Subtype %in% "unclassified", ]
data$Subtype <- factor(data$Subtype, levels=c("MUC", "PRO", "MES"))

data$BOR <- as.character(data$BOR)
data$BOR[data$BOR == "NE"] <- NA
data$BOR[data$BOR %in% "PR"] <- "Partial Response"
data$BOR[data$BOR %in% "SD"] <- "Stable Disease"
data$BOR[data$BOR %in% "PD"] <- "Progressive Disease"
data$BOR <- factor(data$BOR, levels=c("Partial Response", "Stable Disease", "Progressive Disease"))
tbl <- table(data$BOR, data$Subtype)
nbTumors <- sum(tbl)

pdf("../figures/Fig5_gp28363_ResponseBySubtype.pdf", width=5.5, height=6.5)
par(mar=c(4.5, 4.5, 3, 1))
barplot(tbl, xlab="Subtype", ylab="Number of patient lesions", col=c("cyan", "darkblue", "brown2"), 
			main=paste0("gp28363 (n=",nbTumors,")"), cex.axis=1.5, cex.lab=1.8,
			cex.names=1.5, cex.main=1.8, ylim=c(0,10))
legend("topleft", legend=rownames(tbl), fill=c("cyan", "darkblue", "brown2"), cex=1.6, bty="n")
dev.off()
