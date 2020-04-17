library(survival)
library(rmeta)
library(survminer)
library(ggplot2)


data.mut <- read.table("../data/OAK_KRASmut_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
data.mut <- data.mut[!data.mut$Subtype %in% "unclassified", ]
data.mut <- data.mut[, -match("NMF", colnames(data.mut))]
data.wt <- read.table("../data/OAK_KRASwt_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
data.wt <- data.wt[!data.wt$Subtype %in% "unclassified", ]
data <- rbind(data.mut, data.wt)

indices.rm <- which(data$Subtype %in% "MES")
data <- data[-indices.rm, ]
data$Subtype <- factor(data$Subtype, levels=c("MUC", "PRO"))

data$OS.CNSR <- sign(data$TTE.OS)
data$OS.CNSR[data$OS.CNSR==(-1)] <- 0
data$OS <- abs(data$TTE.OS)

# Atezolizumab arm
data <- data[data$Category.ARM %in% "MPDL3280A",]

data <- data[!is.na(data$OS), ]
nbTumors <- nrow(data)

colors <- rep("cyan3", nrow(data))
colors[data$Subtype %in% "PRO"] <- "mediumpurple2"


# Overall survival
fit<- survfit(Surv(time=OS, event=OS.CNSR) ~ Subtype, data=data)
median.1 <- signif(summary(fit)$table[1,"median"],3)
median.2 <- signif(summary(fit)$table[2,"median"],3)

# Drawing survival curves
p <- ggsurvplot(fit, 
	data = data,
   title = paste0("OAK (n=",nbTumors,", atezolizumab)"),
   font.main = c(18, "bold", "black"),
   font.x = c(18, "bold", "black"),
   font.y = c(18, "bold", "black"),
   font.tickslab = 15, 
   color = "strata",
   palette = c("cyan3", "mediumpurple2"),
   xlab = "Overall survival (months)", 
   ylab = "Probability of survival",
   conf.int = FALSE, 
   pval=TRUE, 
   risk.table = FALSE,
   surv.median.line="h",
   legend="none", 
   size=2,
   ggtheme = theme_bw())
p$plot <- p$plot + theme(plot.title = element_text(hjust = 0.5), legend.text=element_text(size=18))
p$plot <- p$plot + ggplot2::annotate("text", x=median.1-1.5, y=0.46, label=median.1, color="cyan3", size=6)
p$plot <- p$plot + ggplot2::annotate("text", x=median.2+3, y=0.5, label=median.2, color="mediumpurple2", size=6)
p$plot <- p$plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
							panel.background = element_blank(), axis.line = element_line(colour = "black"))

pdf("../figures/FigS4A_OAK_MUCPRO_SurvivalBySubtype_Atezolizumab.pdf", width=6, height=6, onefile=FALSE)
p
dev.off()

