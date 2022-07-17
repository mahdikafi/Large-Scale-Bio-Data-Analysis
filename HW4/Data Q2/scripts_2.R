library(glasso)
library(DescTools)
library(sAIC)
library(nethet)


model <- read.csv("PrecisionMatrix.csv")
model <- model[, 2:21]
samples <- read.csv("Samples.csv")
samples <- samples[2:21]
n_rows <- dim(samples)[1]
model[model != 0] <- 1

# 1)
tprs_10 <- vector()
fprs_10 <- vector()
rows_10 <- sample(1:n_rows, size = 10)
samples_10 <- samples[rows_10, ]
x_10 <- var(samples_10)
for (rho in seq(0.0001, 1, 0.01)){
  graph_10 <- glasso(x_10, rho = rho, nobs = 10)
  graph_10 <- graph_10[['wi']]
  graph_10[graph_10 != 0] <- 1
  tp <- sum(model[model == graph_10])
  fn <- sum(model[graph_10 == 0])
  tpr <- tp/(tp+fn)
#  print(tpr)
  fp <- sum(graph_10[model==0])
  tn <- sum(graph_10[model == 0] == 0)
  fpr <- fp/(fp+tn)
#  print(fpr)
  tprs_10 <- c(tprs_10, tpr)
  fprs_10 <- c(fprs_10, fpr)
}
tprs_10 <- sort(tprs_10)
fprs_10 <- sort(fprs_10)
plot(fprs_10, tprs_10, type = 'l', col='red', xlab = 'FPR', ylab = 'TPR')
auc_10 <- AUC(fprs_10, tprs_10)
print(auc_10)

tprs_100 <- vector()
fprs_100 <- vector()
rows_100 <- sample(1:n_rows, size = 100)
samples_100 <- samples[rows_100, ]
x_100 <- var(samples_100)
for (rho in seq(0.0001, 1, 0.01)){
  graph_100 <- glasso(x_100, rho = rho, nobs = 100)
  graph_100 <- graph_100[['wi']]
  graph_100[graph_100 != 0] <- 1
  tp <- sum(model[model == graph_100])
  fn <- sum(model[graph_100 == 0])
  tpr <- tp/(tp+fn)
  #  print(tpr)
  fp <- sum(graph_100[model==0])
  tn <- sum(graph_100[model == 0] == 0)
  fpr <- fp/(fp+tn)
  #  print(fpr)
  tprs_100 <- c(tprs_100, tpr)
  fprs_100 <- c(fprs_100, fpr)
}
tprs_100 <- sort(tprs_100)
fprs_100 <- sort(fprs_100)
lines(fprs_100, tprs_100, type = 'l', col='blue')
auc_100 <- AUC(fprs_100, tprs_100)
print(auc_100)

tprs_1000 <- vector()
fprs_1000 <- vector()
x_1000 <- var(samples)
for (rho in seq(0.0001, 1, 0.01)){
  graph_1000 <- glasso(x_1000, rho = rho, nobs = 1000)
  graph_1000 <- graph_1000[['wi']]
  graph_1000[graph_1000 != 0] <- 1
  tp <- sum(model[model == graph_1000])
  fn <- sum(model[graph_1000 == 0])
  tpr <- tp/(tp+fn)
  #  print(tpr)
  fp <- sum(graph_1000[model==0])
  tn <- sum(graph_1000[model == 0] == 0)
  fpr <- fp/(fp+tn)
  #  print(fpr)
  tprs_1000 <- c(tprs_1000, tpr)
  fprs_1000 <- c(fprs_1000, fpr)
}
tprs_1000 <- sort(tprs_1000)
fprs_1000 <- sort(fprs_1000)
lines(fprs_1000, tprs_1000, type = 'l', col='green')
auc_1000 <- AUC(fprs_1000, tprs_1000)
print(auc_1000)
legend("bottomright", legend=c(paste("10:", auc_10), paste("100:", auc_100), paste("1000:", auc_1000)),
       col=c("red", "blue", "green"), lty = 1:2, cex=0.8)


# 2)
samples <- read.csv("Samples.csv")
samples <- samples[2:21]
n_rows <- dim(samples)[1]
mapped_samples <- exp(samples)

tprs_10 <- vector()
fprs_10 <- vector()
rows_10 <- sample(1:n_rows, size = 10)
samples_10 <- mapped_samples[rows_10, ]
x_10 <- var(samples_10)
for (rho in seq(0.0001, 1, 0.01)){
  graph_10 <- glasso(x_10, rho = rho, nobs = 10)
  graph_10 <- graph_10[['wi']]
  graph_10[graph_10 != 0] <- 1
  tp <- sum(model[model == graph_10])
  fn <- sum(model[graph_10 == 0])
  tpr <- tp/(tp+fn)
  #  print(tpr)
  fp <- sum(graph_10[model==0])
  tn <- sum(graph_10[model == 0] == 0)
  fpr <- fp/(fp+tn)
  #  print(fpr)
  tprs_10 <- c(tprs_10, tpr)
  fprs_10 <- c(fprs_10, fpr)
}
tprs_10 <- sort(tprs_10)
fprs_10 <- sort(fprs_10)
plot(fprs_10, tprs_10, type = 'l', col='red', xlab = 'FPR', ylab = 'TPR')
auc_10 <- AUC(fprs_10, tprs_10)
print(auc_10)


tprs_100 <- vector()
fprs_100 <- vector()
rows_100 <- sample(1:n_rows, size = 100)
samples_100 <- mapped_samples[rows_100, ]
x_100 <- var(samples_100)
for (rho in seq(0.0001, 1, 0.01)){
  graph_100 <- glasso(x_100, rho = rho, nobs = 100)
  graph_100 <- graph_100[['wi']]
  graph_100[graph_100 != 0] <- 1
  tp <- sum(model[model == graph_100])
  fn <- sum(model[graph_100 == 0])
  tpr <- tp/(tp+fn)
  #  print(tpr)
  fp <- sum(graph_100[model==0])
  tn <- sum(graph_100[model == 0] == 0)
  fpr <- fp/(fp+tn)
  #  print(fpr)
  tprs_100 <- c(tprs_100, tpr)
  fprs_100 <- c(fprs_100, fpr)
}
tprs_100 <- sort(tprs_100)
fprs_100 <- sort(fprs_100)
lines(fprs_100, tprs_100, type = 'l', col='blue')
auc_100 <- AUC(fprs_100, tprs_100)
print(auc_100)


tprs_1000 <- vector()
fprs_1000 <- vector()
x_1000 <- var(mapped_samples)
for (rho in seq(0.0001, 1, 0.01)){
  graph_1000 <- glasso(x_1000, rho = rho, nobs = 1000)
  graph_1000 <- graph_1000[['wi']]
  graph_1000[graph_1000 != 0] <- 1
  tp <- sum(model[model == graph_1000])
  fn <- sum(model[graph_1000 == 0])
  tpr <- tp/(tp+fn)
  #  print(tpr)
  fp <- sum(graph_1000[model==0])
  tn <- sum(graph_1000[model == 0] == 0)
  fpr <- fp/(fp+tn)
  #  print(fpr)
  tprs_1000 <- c(tprs_1000, tpr)
  fprs_1000 <- c(fprs_1000, fpr)
}
tprs_1000 <- sort(tprs_1000)
fprs_1000 <- sort(fprs_1000)
lines(fprs_1000, tprs_1000, type = 'l', col='green')
auc_1000 <- AUC(fprs_1000, tprs_1000)
print(auc_1000)
legend("bottomright", legend=c(paste("10:", auc_10), paste("100:", auc_100), paste("1000:", auc_1000)),
       col=c("red", "blue", "green"), lty = 1:2, cex=0.8)

# 3)
# AIC
samples <- read.csv("Samples.csv")
samples_10 <- samples[, 2:11]
model <- read.csv("PrecisionMatrix.csv")
model_10 <- model[1:10, 2:11]
model_10 <- matrix(unlist(model_10), nrow = 10, ncol = 10)
x <- var(samples_10)
best_rho <- Inf
best_aic <- Inf
for (rho in seq(0.0001, 1, 0.01)){
  graph <- glasso(x, rho = rho, nobs = 1000)
  aic <- sAIC(x = model_10, beta = graph$wi, family = "ggm")$AIC
  if (aic < best_aic){
    best_aic <- aic
    best_rho <- rho
  }
}
samples <- samples[2:21]
model <- model[2:21]
x <- var(samples)
graph_aic <- glasso(x, rho = best_rho, nobs = 1000)$wi
print(graph_aic)

# BIC
samples <- read.csv("Samples.csv")
samples_10 <- samples[, 2:11]
best_rho <- screen_bic.glasso(samples_10)$rho.opt
graph_bic <- glasso(var(samples[2:21]), rho = best_rho, nobs = 1000)$wi
print(graph_bic)



