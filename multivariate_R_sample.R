library(survival)

setwd('~/code/microarray-proportional-hazards')
data <- read.csv('test_data/GSE4412-GPL96_series_matrix.csv', header=FALSE)
time <- data[6,]
time[1] <- NULL
time <- as.numeric(as.array(as.matrix(time)))

censor <- data[5,]
censor[[1]] <- NULL
censor <- as.numeric(as.array(as.matrix(censor)))

multi <- data[2,]
multi[[1]] <- NULL
multi <- as.numeric(as.array(as.matrix(multi)))

gene <- data[11,]
gene[[1]] <- NULL
gene <- as.numeric(as.array(as.matrix(gene)))

coxdata <- data.frame(gene, time, censor, multi)

result <- coxph(formula = Surv(time, censor) ~ gene + multi, data = coxdata, model=FALSE, x=FALSE, y=FALSE, singular.ok=TRUE)
summary(result)