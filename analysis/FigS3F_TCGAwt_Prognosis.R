library(survival)
library(rmeta)
library(survminer)
library(ggplot2)

data <- read.table("../data/TCGA_KRASwt_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
data <- data[!data$Subtype %in% "unclassified", ]
data$Subtype <- factor(data$Subtype, levels=c("MUC", "PRO", "MES"))

data$OS <- data$days_to_death_or_last_followup
data$OS.CNSR <- NA
data$OS.CNSR[data$vital_status %in% "Alive"] <- 0
data$OS.CNSR[data$vital_status %in% "Dead"] <- 1
data$OS <- data$OS * 12/365
data <- data[!is.na(data$OS), ]
nbTumors <- nrow(data)

# Censor KRAS-wildtype data at 100 months
data$OS.CNSR[data$OS>100] <- 0 # Alive but censored
data$OS[data$OS>100] <- 100

colors <- rep("cyan3", nbTumors)
colors[data$Subtype %in% "PRO"] <- "mediumpurple2"
colors[data$Subtype %in% "MES"] <- "green3"

fit<- survfit(Surv(time=OS, event=OS.CNSR) ~ Subtype, data=data)
median.1 <- signif(summary(fit)$table[1,"median"],3)
median.2 <- signif(summary(fit)$table[2,"median"],3)
median.3 <- signif(summary(fit)$table[3,"median"],3)

# Drawing survival curves
p <- ggsurvplot(fit, 
	data = data,
   title = paste0("TCGA - KRAS wildtype (n=",nbTumors,")"),
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
   xlim=c(0,100),
   break.x.by=20,
   surv.median.line="h",
   risk.table = FALSE,
   legend="none", 
#   legend.title = "Subtype",
#   legend.labs = c('1', '2', '3'),
   ggtheme = theme_bw()) 
p$plot <- p$plot + theme(plot.title = element_text(hjust = 0.5), legend.text=element_text(size=18))
p$plot <- p$plot + ggplot2::annotate("text", x=median.1-6, y=0.46, label=median.1, color="cyan3", size=6)
p$plot <- p$plot + ggplot2::annotate("text", x=median.2+7, y=0.52, label=median.2, color="mediumpurple2", size=6)
p$plot <- p$plot + ggplot2::annotate("text", x=median.3, y=0.35, label=median.3, color="green3", size=6)
p$plot <- p$plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
							panel.background = element_blank(), axis.line = element_line(colour = "black"))

pdf("../figures/FigS3_TCGAwt_SurvivalBySubtype.pdf", width=6, height=6, onefile=FALSE)
p
dev.off()
