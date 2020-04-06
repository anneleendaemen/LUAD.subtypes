library(survival)
library(rmeta)
library(survminer)
library(ggplot2)


data <- read.table("../data/OAK_KRASmut_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
data <- data[!data$Subtype %in% "unclassified", ]
data$Subtype <- factor(data$Subtype, levels=c("MUC", "PRO", "MES"))

data$OS.CNSR <- sign(data$TTE.OS)
data$OS.CNSR[data$OS.CNSR==(-1)] <- 0
data$OS <- abs(data$TTE.OS)

# Docetaxel arm
data <- data[data$Category.ARM %in% "Docetaxel",]

data <- data[!is.na(data$OS), ]
nbTumors <- nrow(data)

colors <- rep("cyan3", nrow(data))
colors[data$Subtype %in% "PRO"] <- "mediumpurple2"
colors[data$Subtype %in% "MES"] <- "green3"

fit<- survfit(Surv(time=OS, event=OS.CNSR) ~ Subtype, data=data)
median.1 <- signif(summary(fit)$table[1,"median"],3)
median.2 <- signif(summary(fit)$table[2,"median"],3)
median.3 <- signif(summary(fit)$table[3,"median"],3)

# Drawing survival curves
p <- ggsurvplot(fit, 
	data = data,
   title = paste0("OAK - KRAS mutant (n=",nbTumors,", docetaxel)"),
   font.main = c(18, "bold", "black"),
   font.x = c(18, "bold", "black"),
   font.y = c(18, "bold", "black"),
   font.tickslab = 15, 
   color = "strata",
   palette = c("cyan3", "mediumpurple2", "green3"),
   xlab = "Overall survival (months)", 
   ylab = "Probability of survival",
   conf.int = FALSE, 
   pval=TRUE, 
   risk.table = FALSE,
   surv.median.line="h",
   legend="none", 
   ggtheme = theme_bw())
p$plot <- p$plot + theme(plot.title = element_text(hjust = 0.5), legend.text=element_text(size=18))
p$plot <- p$plot + ggplot2::annotate("text", x=median.1-1.5, y=0.47, label=median.1, color="cyan3", size=6)
p$plot <- p$plot + ggplot2::annotate("text", x=median.2-2, y=0.47, label=median.2, color="mediumpurple2", size=6)
p$plot <- p$plot + ggplot2::annotate("text", x=median.3+2.15, y=0.52, label=median.3, color="green3", size=6)
p$plot <- p$plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
							panel.background = element_blank(), axis.line = element_line(colour = "black"))

pdf("../figures/FigS3_OAKmut_SurvivalBySubtype_Docetaxel.pdf", width=6, height=6, onefile=FALSE)
p
dev.off()
