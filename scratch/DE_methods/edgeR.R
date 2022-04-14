library(codaDE)
library(edgeR)
library(caret)
library(cowplot)
library(RSQLite)

data <- readRDS("output/datasets/122f7064-bf7a-422e-ab4f-fbc0564684a0.rds")

# Estimate true percent DE
oracle_calls <- call_DA_NB(data$simulation$abundances, data$simulation$groups)
qvals <- p.adjust(oracle_calls$pval, method = "BH")
oracle_calls <- factor(qvals < 0.05)
levels(oracle_calls) <- c("not DE", "DE")
table(oracle_calls)

nb_calls <- call_DA_NB(data$simulation$observed_counts1, data$simulation$groups)
qvals <- p.adjust(nb_calls$pval, method = "BH")
nb_calls <- factor(qvals < 0.05)
levels(nb_calls) <- c("not DE", "DE")
table(nb_calls)

noTMM_calls <- call_DA_edgeR(data = data$simulation$observed_counts1,
                             groups = data$simulation$groups,
                             use_TMM = FALSE)
qvals <- p.adjust(noTMM_calls, method = "BH")
noTMM_calls <- factor(qvals < 0.05)
levels(noTMM_calls) <- c("not DE", "DE")
table(noTMM_calls)

TMM_calls <- call_DA_edgeR(data = data$simulation$observed_counts1,
                           groups = data$simulation$groups,
                           use_TMM = TRUE)
qvals <- p.adjust(TMM_calls, method = "BH")
TMM_calls <- factor(qvals < 0.05)
levels(TMM_calls) <- c("not DE", "DE")
table(TMM_calls)

caret::confusionMatrix(nb_calls, oracle_calls, positive = "DE")

caret::confusionMatrix(noTMM_calls, oracle_calls, positive = "DE")

caret::confusionMatrix(TMM_calls, oracle_calls, positive = "DE")



