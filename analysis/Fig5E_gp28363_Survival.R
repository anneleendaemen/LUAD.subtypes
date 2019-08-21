library(survival)
library(rmeta)
library(survminer)
library(ggplot2)

data <- read.table("../data/gp28363_subtypes.tsv", sep="\t", header=TRUE, check.names=FALSE)
data <- data[!data$Subtype %in% "unclassified", ]
data$Subtype <- factor(data$Subtype, levels=c("MUC", "PRO", "MES"))
nbTumors <- nrow(data)

data$PFS <- data$'PFS (months)'
data$PFS.CNSR <- NA
data$PFS.CNSR[data$'Censor PFS'==0] <- 1 # 0 = progressed ~ dead
data$PFS.CNSR[data$'Censor PFS'==1] <- 0 # 1 = not progressed ~ alive

colors <- rep("cyan3", nbTumors)
colors[data$Subtype %in% "PRO"] <- "mediumpurple2"
colors[data$Subtype %in% "MES"] <- "green3"

fit<- survfit(Surv(time=PFS, event=PFS.CNSR) ~ Subtype, data=data)
median.1 <- signif(summary(fit)$table[1,"median"],3)
median.2 <- signif(summary(fit)$table[2,"median"],3)
median.3 <- signif(summary(fit)$table[3,"median"],3)

# Drawing survival curves
p <- ggsurvplot(fit, 
	data = data,
   title = paste0("gp28363 (n=",nbTumors,"), cobimetinib + atezolizumab"),
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
   xlim=c(0,30),
   break.x.by=5,
   surv.median.line="h",
   risk.table = FALSE, 
   legend="none",
   ggtheme = theme_bw()) 
p$plot <- p$plot + theme(plot.title = element_text(hjust = 0.5), legend.text=element_text(size=12))
p$plot <- p$plot + ggplot2::annotate("text", x=median.1-1.6, y=0.47, label=median.1, color="cyan3", size=5)
p$plot <- p$plot + ggplot2::annotate("text", x=median.2+0.5, y=0.55, label=median.2, color="mediumpurple2", size=5)
p$plot <- p$plot + ggplot2::annotate("text", x=median.3, y=0.47, label=median.3, color="green3", size=5)
p$plot <- p$plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
							panel.background = element_blank(), axis.line = element_line(colour = "black"))

pdf("../figures/Fig5_gp28363_PFSBySubtype.pdf", width=6, height=6, onefile=FALSE)
p
dev.off()


model1 <- coxph(Surv(time=PFS, event=PFS.CNSR) ~ Subtype, data=data)
model2 <- coxph(Surv(time=PFS, event=PFS.CNSR) ~ MAPK_signature, data=data)
model3 <- coxph(Surv(time=PFS, event=PFS.CNSR) ~ Subtype + MAPK_signature, data=data)
summary(model1)$rsq # 0.3302
summary(model2)$rsq # 0.00427
summary(model3)$rsq # 0.3343
anova(model1, model3) # p 0.7259
anova(model2, model3) # p 0.0178
