library(codaDE)

p <- 500
n <- 100
data <- simulate_RNAseq(p = p, n = n, proportion_de = 0.9, size_factor_correlation = 0.0, spike_in = TRUE)

par(mfrow = c(1,2))

# visualize counts
plot(data$abundances[,data$de_genes[1]], ylab = "orig. abundances")
plot(data$observed_counts[,data$de_genes[1]], ylab = "observed counts")

# visualize spike-in (counts)
plot(data$abundances[,(p+1)], ylab = "orig. abundances")
plot(data$observed_counts[,(p+1)], ylab = "observed counts")

logratios.abundances <- get_logratios(data, call_abundances = TRUE)
logratios.observed_counts <- get_logratios(data, call_abundances = FALSE)

# visualize noise to a DE gene added by ALR representation
plot(logratios.abundances[,data$de_genes[1]], ylab = "ALR orig. abundances")
plot(logratios.observed_counts[,data$de_genes[1]], ylab = "ALR observed counts")

use_lr <- FALSE

tp <- c()
fn <- c()
fp <- c()
tn <- c()
abundance_threshold <- 1
evaluate_gene <- apply(data$observed_counts, 2, function(x) mean(x)/nrow(data$observed_counts) > abundance_threshold)
for(feature_idx in 1:p) {
  if(evaluate_gene[feature_idx]) {
    if(use_lr) {
      gene_data <- data.frame(logratios = logratios.abundances[,feature_idx], groups = data$groups)
      fit <- lm(logratios ~ groups, data = gene_data)
      pval1 <- coef(summary(fit))[2,4]
      gene_data <- data.frame(logratios = logratios.observed_counts[,feature_idx], groups = data$groups)
      fit <- lm(logratios ~ groups, data = gene_data)
      pval2 <- coef(summary(fit))[2,4]
    } else {
      pval1 <- call_DE(data, feature_idx, call_abundances = TRUE)
      pval2 <- call_DE(data, feature_idx, call_abundances = FALSE)
    }
    if(pval1 < 0.05) {
      if(pval2 < 0.05) {
        tp <- c(tp, feature_idx)
      } else {
        fn <- c(fn, feature_idx)
      }
    } else {
      if(pval2 < 0.05) {
        fp <- c(fp, feature_idx)
      } else {
        tn <- c(tn, feature_idx)
      }
    }
  }
}

genes_evaluated <- sum(evaluate_gene)
cat("TP:",length(tp),"/",genes_evaluated,"(",round(length(tp)/genes_evaluated, 2),")\n")
cat("TN:",length(tn),"/",genes_evaluated,"(",round(length(tn)/genes_evaluated, 2),")\n")

cat("FP:",length(fp),"/",genes_evaluated,"(",round(length(fp)/genes_evaluated, 2),")\n")
cat("FN:",length(fn),"/",genes_evaluated,"(",round(length(fn)/genes_evaluated, 2),")\n")

plot_rg <- function(orig_mat, observed_mat, rg) {
  plot(orig_mat[,rg], ylab = "orig. abundances")
  plot(observed_mat[,rg], ylab = "observed counts")
}

# visualize a random TP
rg <- sample(tp)[1]
if(use_lr) {
  plot_rg(logratios.abundances, logratios.observed_counts, rg)
} else {
  plot_rg(data$abundances, data$observed_counts, rg)
}

# visualize a random TN
rg <- sample(tn)[1]
if(use_lr) {
  plot_rg(logratios.abundances, logratios.observed_counts, rg)
} else {
  plot_rg(data$abundances, data$observed_counts, rg)
}

# visualize a random FN
rg <- sample(fn)[1]
if(use_lr) {
  plot_rg(logratios.abundances, logratios.observed_counts, rg)
} else {
  plot_rg(data$abundances, data$observed_counts, rg)
}

# visualize a random FP
rg <- sample(fp)[1]
if(use_lr) {
  plot_rg(logratios.abundances, logratios.observed_counts, rg)
} else {
  plot_rg(data$abundances, data$observed_counts, rg)
}

# demonstrate the effect of zero as a FLOOR on counts
# basically: vanishingly abundant things decrease in relative abundance
#            proportionally the same amount as more abundant things but
#            their likelihood of being realizes as 
counts1 <- c(1, 2, 3, 50, 500)
props1 <- counts1/sum(counts1)

counts2 <- c(1, 2, 3, 50, 1000)
props2 <- counts2/sum(counts2)

# the proportional change is the same
props1[1]/props2[1]
props1[2]/props2[2]
props1[4]/props2[4]

# draw at same depth
counts_obs1 <- round(props1*sum(counts2))
counts_obs2 <- round(props2*sum(counts2))

# expect ratios to stay close; they basically do
cat(round((counts_obs1[1] + 0.5)/(counts_obs1[4] + 0.5), 3),"--",round((counts_obs2[1] + 0.5)/(counts_obs2[4] + 0.5), 3),"\n")
cat(round((counts_obs1[2] + 0.5)/(counts_obs1[4] + 0.5), 3),"--",round((counts_obs2[2] + 0.5)/(counts_obs2[4] + 0.5), 3),"\n")
cat(round((counts_obs1[3] + 0.5)/(counts_obs1[4] + 0.5), 3),"--",round((counts_obs2[3] + 0.5)/(counts_obs2[4] + 0.5), 3),"\n")

# draw at increased depth
counts_obs1 <- round(props1*sum(counts2)*3)
counts_obs2 <- round(props2*sum(counts2)*3)

# a tiny increase
cat(round((counts_obs1[1] + 0.5)/(counts_obs1[4] + 0.5), 3),"--",round((counts_obs2[1] + 0.5)/(counts_obs2[4] + 0.5), 3),"\n")
cat(round((counts_obs1[2] + 0.5)/(counts_obs1[4] + 0.5), 3),"--",round((counts_obs2[2] + 0.5)/(counts_obs2[4] + 0.5), 3),"\n")
cat(round((counts_obs1[3] + 0.5)/(counts_obs1[4] + 0.5), 3),"--",round((counts_obs2[3] + 0.5)/(counts_obs2[4] + 0.5), 3),"\n")

# draw at decreased depth
counts_obs1 <- round(props1*sum(counts2)*0.5)
counts_obs2 <- round(props2*sum(counts2)*0.5)

# a more noticeable increase
cat(round((counts_obs1[1] + 0.5)/(counts_obs1[4] + 0.5), 3),"--",round((counts_obs2[1] + 0.5)/(counts_obs2[4] + 0.5), 3),"\n")
cat(round((counts_obs1[2] + 0.5)/(counts_obs1[4] + 0.5), 3),"--",round((counts_obs2[2] + 0.5)/(counts_obs2[4] + 0.5), 3),"\n")
cat(round((counts_obs1[3] + 0.5)/(counts_obs1[4] + 0.5), 3),"--",round((counts_obs2[3] + 0.5)/(counts_obs2[4] + 0.5), 3),"\n")

# increasing K; proportional changes are harder to "realize" for vanishingly abundant features
# thus, it always looks like the change in the reference is LARGER
# since that change is grossly a decrease, the logratio appears to increase
#
# the protective aspect of resampling at the same size is that the least amount of distortion in the ratio is possible
#
# in general resampling at a lower abundance shouldn't product FP's because things will just drop out
a <- 0.01
b <- 10
plot(1:1000, sapply(1:1000, function(i) round(a*i)/round(b*i)))

# just evaluating DE on things about 1% abundance should fix this



