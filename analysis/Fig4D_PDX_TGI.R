study <- "PDX models"
data <- read.table("../data/PDX_TGI.tsv", sep="\t", header=TRUE, check.names=FALSE)
data <- data[!data$Subtype %in% "unclassified", ]
data$Subtype <- factor(data$Subtype, levels=c("MUC", "PRO", "MES"))
nbTumors <- nrow(data)

model1 <- lm(data$`% TGI Cobimetinib mean` ~ data$Subtype)
model2 <- lm(data$`% TGI Cobimetinib mean` ~ data$`MAPK signature`)
model3 <- lm(data$`% TGI Cobimetinib mean` ~ data$Subtype + data$`MAPK signature`)
summary(model1)$r.squared # 0.22457
summary(model2)$r.squared # 0.09049
summary(model3)$r.squared # 0.23589
