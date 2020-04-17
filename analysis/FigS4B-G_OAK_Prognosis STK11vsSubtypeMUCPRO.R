library(survival)
library(rmeta)
library(survminer)
library(ggplot2)


### ATEZOLIZUMAB - KRAS-mutant + KRAS-wildtype ###

data1 <- read.table("../data/OAK_KRASmut_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
data2 <- read.table("../data/OAK_KRASwt_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
data1 <- data1[, -grep("NMF", colnames(data1))]
data1 <- cbind(data1, KRAS="MUT")
data2 <- cbind(data2, KRAS="WT")
data <- rbind(data1, data2)
data <- data[!data$Subtype %in% "unclassified", ]
data$Subtype <- factor(data$Subtype, levels=c("MUC", "PRO", "MES"))

data$OS.CNSR <- sign(data$TTE.OS)
data$OS.CNSR[data$OS.CNSR==(-1)] <- 0
data$OS <- abs(data$TTE.OS)
data <- data[data$Category.ARM %in% "MPDL3280A",]
data <- data[!is.na(data$OS), ]

### ANOVA comparing subtype vs. STK11 status for atezolizumab, restricted to MUC+PRO subtypes (MES = STK11-wt)
data.subtype <- data[data$Subtype %in% c("MUC","PRO"), ]
data.subtype$Subtype <- factor(as.character(data.subtype$Subtype), levels=c("MUC", "PRO"))

fit2 <- coxph(Surv(time=OS, event=OS.CNSR) ~ STK11_mutation_or_homozygous_loss, data=data.subtype)
fit3 <- coxph(Surv(time=OS, event=OS.CNSR) ~ Subtype, data=data.subtype)
fit4 <- coxph(Surv(time=OS, event=OS.CNSR) ~ Subtype + STK11_mutation_or_homozygous_loss, data=data.subtype)
fit5 <- coxph(Surv(time=OS, event=OS.CNSR) ~ Subtype * STK11_mutation_or_homozygous_loss, data=data.subtype)
anova(fit2, fit4)
anova(fit3, fit4)
anova(fit4, fit5)


data.subtype <- data[data$Subtype %in% "MUC", ]
nbTumors <- nrow(data.subtype)

fit<- survfit(Surv(time=OS, event=OS.CNSR) ~ STK11_mutation_or_homozygous_loss, data=data.subtype)
median.1 <- signif(summary(fit)$table[1,"median"],3)
median.2 <- signif(summary(fit)$table[2,"median"],3)

# Drawing survival curves
p <- ggsurvplot(fit, 
	data = data.subtype,
   title = paste0("OAK - MUC (n=",nbTumors,", atezolizumab)"),
   font.main = c(18, "bold", "black"),
   font.x = c(18, "bold", "black"),
   font.y = c(18, "bold", "black"),
   font.tickslab = 15, 
   linetype = "strata",
   palette = c("cyan3", "cyan3"),
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
p$plot <- p$plot + ggplot2::annotate("text", x=median.1-2.5, y=0.55, label=median.1, color="black", size=6)
p$plot <- p$plot + ggplot2::annotate("text", x=median.2+2.5, y=0.53, label=median.2, color="black", size=6)
p$plot <- p$plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
							panel.background = element_blank(), axis.line = element_line(colour = "black"))

pdf("../figures/FigS4B_OAK_MUC_SurvivalBySTK11_Atezo.pdf", width=6, height=6, onefile=FALSE)
p
dev.off()


data.subtype <- data[data$Subtype %in% "PRO", ]
nbTumors <- nrow(data.subtype)

fit<- survfit(Surv(time=OS, event=OS.CNSR) ~ STK11_mutation_or_homozygous_loss, data=data.subtype)
median.1 <- signif(summary(fit)$table[1,"median"],3)
median.2 <- signif(summary(fit)$table[2,"median"],3)

# Drawing survival curves
p <- ggsurvplot(fit, 
	data = data.subtype,
   title = paste0("OAK - PRO (n=",nbTumors,", atezolizumab)"),
   font.main = c(18, "bold", "black"),
   font.x = c(18, "bold", "black"),
   font.y = c(18, "bold", "black"),
   font.tickslab = 15, 
   linetype = "strata",
   palette = c("mediumpurple2", "mediumpurple2"),
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
p$plot <- p$plot + ggplot2::annotate("text", x=median.1-2.5, y=0.45, label=median.1, color="black", size=6)
p$plot <- p$plot + ggplot2::annotate("text", x=median.2+1.5, y=0.55, label=median.2, color="black", size=6)
p$plot <- p$plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
							panel.background = element_blank(), axis.line = element_line(colour = "black"))

pdf("../figures/FigS4C_OAK_PRO_SurvivalBySTK11_Atezo.pdf", width=6, height=6, onefile=FALSE)
p
dev.off()



### DOCETAXEL - KRAS-wildtype ###

data <- read.table("../data/OAK_KRASwt_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
data <- data[!data$Subtype %in% "unclassified", ]
data$Subtype <- factor(data$Subtype, levels=c("MUC", "PRO", "MES"))

data$OS.CNSR <- sign(data$TTE.OS)
data$OS.CNSR[data$OS.CNSR==(-1)] <- 0
data$OS <- abs(data$TTE.OS)
data <- data[data$Category.ARM %in% "Docetaxel",]
data <- data[!is.na(data$OS), ]


### ANOVA comparing subtype vs. STK11 status for docetaxel, restricted to MUC+PRO subtypes (MES = STK11-wt)
data.subtype <- data[data$Subtype %in% c("MUC","PRO"), ]
data.subtype$Subtype <- factor(as.character(data.subtype$Subtype), levels=c("MUC", "PRO"))
fit <- coxph(Surv(time=OS, event=OS.CNSR) ~ Subtype + STK11_mutation_or_homozygous_loss, data= data.subtype)
fit2 <- coxph(Surv(time=OS, event=OS.CNSR) ~ Subtype, data=data.subtype)
fit3 <- coxph(Surv(time=OS, event=OS.CNSR) ~ STK11_mutation_or_homozygous_loss, data=data.subtype)
anova(fit2, fit)
anova(fit3, fit)


data.subtype <- data[data$Subtype %in% "MUC", ]
nbTumors <- nrow(data.subtype)

fit<- survfit(Surv(time=OS, event=OS.CNSR) ~ STK11_mutation_or_homozygous_loss, data=data.subtype)
median.1 <- signif(summary(fit)$table[1,"median"],3)
median.2 <- signif(summary(fit)$table[2,"median"],3)

# Drawing survival curves
p <- ggsurvplot(fit, 
	data = data.subtype,
   title = paste0("OAK - KRAS wt, MUC (n=",nbTumors,", docetaxel)"),
   font.main = c(18, "bold", "black"),
   font.x = c(18, "bold", "black"),
   font.y = c(18, "bold", "black"),
   font.tickslab = 15, 
   linetype = "strata",
   palette = c("cyan3", "cyan3"),
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
p$plot <- p$plot + ggplot2::annotate("text", x=median.1+2.5, y=0.45, label=median.1, color="black", size=6)
p$plot <- p$plot + ggplot2::annotate("text", x=median.2+2.5, y=0.53, label=median.2, color="black", size=6)
p$plot <- p$plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
							panel.background = element_blank(), axis.line = element_line(colour = "black"))

pdf("../figures/FigS4D_OAKwt_MUC_SurvivalBySTK11_Docetaxel.pdf", width=6, height=6, onefile=FALSE)
p
dev.off()


data.subtype <- data[data$Subtype %in% "PRO", ]
nbTumors <- nrow(data.subtype)

fit<- survfit(Surv(time=OS, event=OS.CNSR) ~ STK11_mutation_or_homozygous_loss, data=data.subtype)
median.1 <- signif(summary(fit)$table[1,"median"],3)
median.2 <- signif(summary(fit)$table[2,"median"],3)

# Drawing survival curves
p <- ggsurvplot(fit, 
	data = data.subtype,
   title = paste0("OAK - KRAS wt, PRO (n=",nbTumors,", docetaxel)"),
   font.main = c(18, "bold", "black"),
   font.x = c(18, "bold", "black"),
   font.y = c(18, "bold", "black"),
   font.tickslab = 15, 
   linetype = "strata",
   palette = c("mediumpurple2", "mediumpurple2"),
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
p$plot <- p$plot + ggplot2::annotate("text", x=median.1-2, y=0.45, label=median.1, color="black", size=6)
p$plot <- p$plot + ggplot2::annotate("text", x=median.2+2.4, y=0.55, label=median.2, color="black", size=6)
p$plot <- p$plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
							panel.background = element_blank(), axis.line = element_line(colour = "black"))

pdf("../figures/FigS4E_OAKwt_PRO_SurvivalBySTK11_Docetaxel.pdf", width=6, height=6, onefile=FALSE)
p
dev.off()




### DOCETAXEL - KRAS-mutant ###

data <- read.table("../data/OAK_KRASmut_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
data <- data[!data$Subtype %in% "unclassified", ]
data$Subtype <- factor(data$Subtype, levels=c("MUC", "PRO", "MES"))

data$OS.CNSR <- sign(data$TTE.OS)
data$OS.CNSR[data$OS.CNSR==(-1)] <- 0
data$OS <- abs(data$TTE.OS)
data <- data[data$Category.ARM %in% "Docetaxel",]
data <- data[!is.na(data$OS), ]


### ANOVA comparing subtype vs. STK11 status for docetaxel, restricted to MUC+PRO subtypes (MES = STK11-wt)
data.subtype <- data[data$Subtype %in% c("MUC","PRO"), ]
data.subtype$Subtype <- factor(as.character(data.subtype$Subtype), levels=c("MUC", "PRO"))
fit <- coxph(Surv(time=OS, event=OS.CNSR) ~ Subtype + STK11_mutation_or_homozygous_loss, data= data.subtype)
fit2 <- coxph(Surv(time=OS, event=OS.CNSR) ~ Subtype, data=data.subtype)
fit3 <- coxph(Surv(time=OS, event=OS.CNSR) ~ STK11_mutation_or_homozygous_loss, data=data.subtype)
anova(fit2, fit)
anova(fit3, fit)


data.subtype <- data[data$Subtype %in% "MUC", ]
nbTumors <- nrow(data.subtype)

fit<- survfit(Surv(time=OS, event=OS.CNSR) ~ STK11_mutation_or_homozygous_loss, data=data.subtype)
median.1 <- signif(summary(fit)$table[1,"median"],3)
median.2 <- signif(summary(fit)$table[2,"median"],3)

# Drawing survival curves
p <- ggsurvplot(fit, 
	data = data.subtype,
   title = paste0("OAK - KRAS mut, MUC (n=",nbTumors,", docetaxel)"),
   font.main = c(18, "bold", "black"),
   font.x = c(18, "bold", "black"),
   font.y = c(18, "bold", "black"),
   font.tickslab = 15, 
   linetype = "strata",
   palette = c("cyan3", "cyan3"),
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
p$plot <- p$plot + ggplot2::annotate("text", x=median.1, y=0.55, label=median.1, color="black", size=6)
p$plot <- p$plot + ggplot2::annotate("text", x=median.2-2, y=0.46, label=median.2, color="black", size=6)
p$plot <- p$plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
							panel.background = element_blank(), axis.line = element_line(colour = "black"))

pdf("../figures/FigS4G_OAKmut_MUC_SurvivalBySTK11_Docetaxel.pdf", width=6, height=6, onefile=FALSE)
p
dev.off()


data.subtype <- data[data$Subtype %in% "PRO", ]
nbTumors <- nrow(data.subtype)

fit<- survfit(Surv(time=OS, event=OS.CNSR) ~ STK11_mutation_or_homozygous_loss, data=data.subtype)
median.1 <- signif(summary(fit)$table[1,"median"],3)
median.2 <- signif(summary(fit)$table[2,"median"],3)

# Drawing survival curves
p <- ggsurvplot(fit, 
	data = data.subtype,
   title = paste0("OAK - KRAS mut, PRO (n=",nbTumors,", docetaxel)"),
   font.main = c(18, "bold", "black"),
   font.x = c(18, "bold", "black"),
   font.y = c(18, "bold", "black"),
   font.tickslab = 15, 
   linetype = "strata",
   palette = c("mediumpurple2", "mediumpurple2"),
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
p$plot <- p$plot + ggplot2::annotate("text", x=median.1-2, y=0.46, label=median.1, color="black", size=6)
p$plot <- p$plot + ggplot2::annotate("text", x=median.2+2.4, y=0.5, label=median.2, color="black", size=6)
p$plot <- p$plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
							panel.background = element_blank(), axis.line = element_line(colour = "black"))

pdf("../figures/FigS4H_OAKmut_PRO_SurvivalBySTK11_Docetaxel.pdf", width=6, height=6, onefile=FALSE)
p
dev.off()



