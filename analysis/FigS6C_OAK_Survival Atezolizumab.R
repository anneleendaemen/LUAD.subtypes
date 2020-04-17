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
data$Subtype <- factor(data$Subtype, levels=c("MUC", "PRO", "MES"))

data$OS.CNSR <- sign(data$TTE.OS)
data$OS.CNSR[data$OS.CNSR==(-1)] <- 0
data$OS <- abs(data$TTE.OS)

data$PFS.CNSR <- sign(data$TTE.PFS)
data$PFS.CNSR[data$PFS.CNSR==(-1)] <- 0
data$PFS <- abs(data$TTE.PFS)

# Atezolizumab arm
data <- data[data$Category.ARM %in% "MPDL3280A",]

data <- data[!is.na(data$OS), ]
nbTumors <- nrow(data)

colors <- rep("cyan3", nrow(data))
colors[data$Subtype %in% "PRO"] <- "mediumpurple2"
colors[data$Subtype %in% "MES"] <- "green3"


# Overall survival
fit<- survfit(Surv(time=OS, event=OS.CNSR) ~ Subtype, data=data)
median.1 <- signif(summary(fit)$table[1,"median"],3)
median.2 <- signif(summary(fit)$table[2,"median"],3)
median.3 <- signif(summary(fit)$table[3,"median"],3)

# Drawing survival curves
p <- ggsurvplot(fit, 
	data = data,
   title = paste0("OAK (n=",nbTumors,", atezolizumab)"),
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
p$plot <- p$plot + ggplot2::annotate("text", x=median.2, y=0.38, label=median.2, color="mediumpurple2", size=6)
p$plot <- p$plot + ggplot2::annotate("text", x=median.3+1.5, y=0.52, label=median.3, color="green3", size=6)
p$plot <- p$plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
							panel.background = element_blank(), axis.line = element_line(colour = "black"))

pdf("../figures/FigS6_OAK_SurvivalBySubtype_Atezolizumab.pdf", width=6, height=6, onefile=FALSE)
p
dev.off()


# Progression free survival
fit<- survfit(Surv(time=PFS, event=PFS.CNSR) ~ Subtype, data=data)
median.1 <- signif(summary(fit)$table[1,"median"],3)
median.2 <- signif(summary(fit)$table[2,"median"],3)
median.3 <- signif(summary(fit)$table[3,"median"],3)

# Drawing survival curves
p <- ggsurvplot(fit, 
	data = data,
   title = paste0("OAK (n=",nbTumors,", atezolizumab)"),
   font.main = c(18, "bold", "black"),
   font.x = c(18, "bold", "black"),
   font.y = c(18, "bold", "black"),
   font.tickslab = 15, 
   color = "strata",
   palette = c("cyan3", "mediumpurple2", "green3"),
   xlab = "Progression free survival (months)", 
   ylab = "Probability of survival",
   conf.int = FALSE, 
   pval=TRUE, 
   risk.table = FALSE,
   surv.median.line="h",
   legend="none", 
   ggtheme = theme_bw())
p$plot <- p$plot + theme(plot.title = element_text(hjust = 0.5), legend.text=element_text(size=18))
p$plot <- p$plot + ggplot2::annotate("text", x=median.1-2.2, y=0.47, label=median.1, color="cyan3", size=6)
p$plot <- p$plot + ggplot2::annotate("text", x=median.2-0.8, y=0.38, label=median.2, color="mediumpurple2", size=6)
p$plot <- p$plot + ggplot2::annotate("text", x=median.3+2.5, y=0.5, label=median.3, color="green3", size=6)
p$plot <- p$plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
							panel.background = element_blank(), axis.line = element_line(colour = "black"))

pdf("../figures/FigS6_OAK_PFSBySubtype_Atezolizumab.pdf", width=6, height=6, onefile=FALSE)
p
dev.off()
