source("path_fix.R")

library(RSQLite)
library(tidyverse)
library(codaDE)
library(entropy)

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))
output_dir <- file.path("output", "datasets")
pseudocount <- 0.1

# ------------------------------------------------------------------------------
#   Sync folder and DB
# ------------------------------------------------------------------------------

simulation_fn <- list.files(output_dir, pattern = ".rds")
fs_uuids <- unname(sapply(simulation_fn, function(x) substr(x, 1, 36)))

# Remove DB entries w/ missing files
ds_uuids <- dbGetQuery(conn, "SELECT DISTINCT UUID FROM datasets")$UUID
char_uuids <- dbGetQuery(conn, "SELECT DISTINCT UUID FROM characteristics")$UUID
res_uuids <- dbGetQuery(conn, "SELECT DISTINCT UUID FROM results")$UUID
db_uuids <- unique(c(ds_uuids, char_uuids, res_uuids))

# Missing locally
unmatched_set <- setdiff(db_uuids, fs_uuids)
cat("Found",length(unmatched_set),"orphaned DB entries...\n")
if(length(unmatched_set) > 0) {
  uuid_list_str <- paste0("('",paste(unmatched_set, collapse = "','"),"')")
  discard <- dbExecute(conn, paste0("DELETE FROM datasets WHERE uuid IN ",
                                    uuid_list_str,";"))
  discard <- dbExecute(conn, paste0("DELETE FROM characteristics WHERE uuid IN ",
                                    uuid_list_str,";"))
  discard <- dbExecute(conn, paste0("DELETE FROM results WHERE uuid IN ",
                                    uuid_list_str,";"))
}

# Missing in DB
unmatched_set <- setdiff(fs_uuids, db_uuids)
cat("Found",length(unmatched_set),"orphaned files...\n")
if(length(unmatched_set) > 0) {
  for(uuid in unmatched_set) {
    unmatched_fn <- file.path(output_dir, paste0(uuid, ".rds"))
    discard <- file.remove(unmatched_fn)
  }
}

# If any features are missing, we'll update them all
update_uuids <- dbGetQuery(conn,
                           paste0("SELECT DISTINCT UUID FROM characteristics ",
                                  "WHERE ",
                                  "FOLD_CHANGE IS NULL OR ",
                                  "MEAN_CORR IS NULL OR ",
                                  "MEDIAN_CORR IS NULL OR ",
                                  "BASE_SPARSITY IS NULL OR ",
                                  "DELTA_SPARSITY IS NULL OR ",
                                  "PERCENT_STABLE IS NULL OR ",
                                  "DIR_CONSENSUS IS NULL OR ",
                                  "MAX_DELTA_REL IS NULL OR ",
                                  "MEDIAN_DELTA_REL IS NULL OR ",
                                  "BASE_ENTROPY IS NULL OR ",
                                  "DELTA_ENTROPY IS NULL;"))$UUID
for(u in 1:length(update_uuids)) {
  uuid <- update_uuids[u]
  cat(paste0("Updating ", uuid, " (", u, "/", length(update_uuids), ")\n"))
  fn <- file.path(output_dir, paste0(uuid, ".rds"))
  dataset <- readRDS(fn)
  
  for(partial in c(0, 1)) {
    if(partial == 1) {
      counts <- dataset$simulation$observed_counts2
    } else {
      counts <- dataset$simulation$observed_counts1
    }
    n <- nrow(counts)/2
    p <- ncol(counts)
    counts_A <- counts[1:n,]
    counts_B <- counts[(n+1):(n*2),]
    
    # Eliminate all-zero features (for now) by spiking in ones
    if(min(colSums(counts_A)) == 0) {
      counts_A <- spike_in_ones(counts_A)
    }
    if(min(colSums(counts_B)) == 0) {
      counts_B <- spike_in_ones(counts_B)
    }

    # --------------------------------------------------------------------------
    #   Calculate fold change
    # --------------------------------------------------------------------------
    
    mean_totalA <- mean(rowSums(counts_A))
    mean_totalB <- mean(rowSums(counts_B))
    fold_diff <- max(c(mean_totalA, mean_totalB)) / min(c(mean_totalA, mean_totalB))
    
    # --------------------------------------------------------------------------
    #   Calculate observed correlation between features (base case)
    # --------------------------------------------------------------------------
    
    observed_correlation <- cor(log(counts_A + 0.5))
    observed_correlation <- observed_correlation[upper.tri(observed_correlation,
                                                           diag = FALSE)]
    mean_cor <- mean(observed_correlation)
    median_cor <- median(observed_correlation)
    
    # --------------------------------------------------------------------------
    #   Calculate the sparsity (% zeros) in each condition
    # --------------------------------------------------------------------------
    
    sparsity_A <- sum(counts_A == 0)/(nrow(counts_A)*ncol(counts_A))
    sparsity_B <- sum(counts_B == 0)/(nrow(counts_B)*ncol(counts_B))
    delta_sparsity <- abs(sparsity_B - sparsity_A)
    
    # --------------------------------------------------------------------------
    #   Calculate percent "stable" features
    # --------------------------------------------------------------------------
    
    rel_abA <- t(apply(counts_A, 1, function(x) x/sum(x)))
    rel_abB <- t(apply(counts_B, 1, function(x) x/sum(x)))
    mean_relA <- colMeans(rel_abA)
    mean_relB <- colMeans(rel_abB)
    avg_fc <- sapply(1:p, function(x) {
      (max(c(mean_relA[x], mean_relB[x])) + pseudocount) /
        (min(c(mean_relA[x], mean_relB[x])) + pseudocount)
    })
    percent_stable <- sum(avg_fc < 1.5)/p
  
    # --------------------------------------------------------------------------
    #   Get % consensus in direction of change
    # --------------------------------------------------------------------------
    
    dir <- (mean_relB + pseudocount)/(mean_relA + pseudocount)
    agreement_dir <- max(sum(dir >= 1)/p, sum(dir < 1)/p)
    
    # --------------------------------------------------------------------------
    #   Max, median delta relative abundance
    # --------------------------------------------------------------------------
    
    delta_rel <- rel_abB - rel_abA
    max_d_rel <- max(delta_rel)
    median_d_rel <- median(delta_rel)
    
    # --------------------------------------------------------------------------
    #   Baseline entropy and change in entropy
    # --------------------------------------------------------------------------
    
    cuts <- quantile(c(rel_abA[rel_abA > 0], rel_abB[rel_abB > 0]),
                     probs = seq(from = 0, to = 1, length.out = 20))
    e1 <- entropy(table(cut(rel_abA, breaks = c(0, cuts))))
    e2 <- entropy(table(cut(rel_abB, breaks = c(0, cuts))))
    d_entropy <- abs(e2 - e1)
    
    # --------------------------------------------------------------------------
    #   Variation in sequencing depth
    # --------------------------------------------------------------------------
    
    # totals <- rowSums(counts)
    # seq_var <- sd(totals)/mean(totals)
    
    # If exists, update
    results <- dbGetQuery(conn,
                          paste0("SELECT * FROM characteristics WHERE UUID='",
                                 uuid, "' AND PARTIAL = 0;"))
    
    if(nrow(results) > 0) {
      # UPDATE
      discard <- dbExecute(conn, paste0("UPDATE characteristics SET ",
                                        "FOLD_CHANGE=",fold_diff,", ",
                                        "MEAN_CORR=",mean_cor,", ",
                                        "MEDIAN_CORR=",median_cor,", ",
                                        "BASE_SPARSITY=",sparsity_A,", ",
                                        "DELTA_SPARSITY=",delta_sparsity,", ",
                                        "PERCENT_STABLE=",percent_stable,", ",
                                        "DIR_CONSENSUS=",agreement_dir,", ",
                                        "MAX_DELTA_REL=",max_d_rel,", ",
                                        "MEDIAN_DELTA_REL=",median_d_rel,", ",
                                        "BASE_ENTROPY=",e1,", ",
                                        "DELTA_ENTROPY=",d_entropy," WHERE ",
                                        "UUID='",uuid,"' AND ",
                                        "PARTIAL=",partial,";"))
    } else {
      # INSERT
      discard <- dbExecute(conn,
                           paste0("INSERT INTO characteristics(UUID, PARTIAL, ",
                                  "FOLD_CHANGE, MEAN_CORR, MEDIAN_CORR, ",
                                  "BASE_SPARSITY, DELTA_SPARSITY, ",
                                  "PERCENT_STABLE, DIR_CONSENSUS, ",
                                  "MAX_DELTA_REL, MEDIAN_DELTA_REL, ",
                                  "BASE_ENTROPY, DELTA_ENTROPY) VALUES(",
                                  "'", uuid, "', 0, ",
                                  fold_diff, ", ",
                                  mean_cor, ", ",
                                  median_cor, ", ",
                                  sparsity_A, ", ",
                                  delta_sparsity, ", ",
                                  percent_stable, ", ",
                                  agreement_dir, ", ",
                                  max_d_rel, ", ",
                                  median_d_rel, ", ",
                                  e1, ", ",
                                  d_entropy, ");"))
             
    }
  }
}
# dbGetQuery(conn, "SELECT * FROM characteristics;")

dbDisconnect(conn)
