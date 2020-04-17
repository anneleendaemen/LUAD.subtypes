study <- "OAK"
data1 <- cbind(read.table("../data/OAK_KRASmut_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE), KRAS="MUT")
data2 <- cbind(read.table("../data/OAK_KRASwt_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE), KRAS="WT")
data1 <- data1[,-which(colnames(data1) %in% "NMF")]
data <- rbind(data1, data2)
data <- data[!data$Subtype %in% "unclassified", ]
data$Subtype <- factor(data$Subtype, levels=c("MUC", "PRO", "MES"))
nbTumors <- nrow(data)

data$TMB <- as.character(data$TMB)
data$TMB[data$TMB=="<=1"] <- "1"
data$TMB <- as.numeric(data$TMB)

colors <- rep("cyan3", nbTumors)
colors[data$Subtype %in% "PRO"] <- "mediumpurple2"
colors[data$Subtype %in% "MES"] <- "green3"
pchs <- rep(19, nbTumors)
pchs[data$KRAS %in% "WT"] <- 17

pdf("../figures/FigS6_OAKall_TMBBySubtype.pdf", width=5.5, height=6.5)
par(mar=c(4.5, 4.5, 3, 1))
boxplot(data$TMB ~ data$Subtype, outline=FALSE,
				ylab="Tumor mutation load", xlab="Subtype", main=paste0(study," (n=",nbTumors,")"),
				cex.lab=1.8, cex.axis=1.5, cex.main=1.8, ylim=range(data$TMB[!is.na(data$TMB)]), log="y")
points(jitter(as.numeric(data$Subtype)), data$TMB,
					col=colors, pch=pchs, cex=1.5)
pval <- signif(kruskal.test(data$TMB ~ data$Subtype)$p.value, 2)
text(x=1, y=max(data$TMB[!is.na(data$TMB)]), labels=paste0("p=",pval), cex=1.5)
dev.off()


