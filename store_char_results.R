source("path_fix.R")

library(RSQLite)
library(dplyr)
library(optparse)

option_list = list(
  make_option(c("--dir"), type = "character", default = "temp",
              help = "directory of evaluation output files", metavar = "character")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

dir <- opt$dir

if(!file.exists(dir)) {
  stop("No such output directory!\n")
}

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))

file_list <- list.files(dir, pattern = "output_char")
for(file in file_list) {
  results <- read.table(file.path(dir, file), header = TRUE)
  cat(paste0("Parsing file: ", file, "\n"))
  updates <- 0
  for(i in 1:nrow(results)) {
    job <- results[i,]
  
    if(job$partial_info == 0 & job$type == "relative_abundances") {
          updates <- updates + dbExecute(conn,
                                         paste0("UPDATE datasets SET ",
                                                "MED_ABS_TOTAL=", job$med_abs, ", ",
                                                "MED_REL_TOTAL=", job$med_rel, " ",
                                                "WHERE UUID='", job$uuid, "';"))
    }

    updates <- updates + dbExecute(conn,
                                   paste0("INSERT OR REPLACE INTO characteristics(UUID, PARTIAL, ",
                                          "TYPE, ",
                                          "TOTALS_C_FC, ",
                                          "TOTALS_C_D, ",
                                          "TOTALS_C_MAX_D, ",
                                          "TOTALS_C_MED_D, ",
                                          "TOTALS_C_SD_D, ",
                                          "CORR_RA_MED, ",
                                          "CORR_RA_SD, ",
                                          "CORR_RA_SKEW, ",
                                          "CORR_LOG_MED, ",
                                          "CORR_LOG_SD, ",
                                          "CORR_LOG_SKEW, ",
                                          "CORR_CLR_MED, ",
                                          "CORR_CLR_SD, ",
                                          "CORR_CLR_SKEW, ",
                                          "COMP_C_P0_A, ",
                                          "COMP_C_P0_B, ",
                                          "COMP_C_P1_A, ",
                                          "COMP_C_P1_B, ",
                                          "COMP_C_P5_A, ",
                                          "COMP_C_P5_B, ",
                                          "COMP_RA_P01_A, ",
                                          "COMP_RA_P01_B, ",
                                          "COMP_RA_P1_A, ",
                                          "COMP_RA_P1_B, ",
                                          "COMP_RA_P5_A, ",
                                          "COMP_RA_P5_B, ",
                                          "COMP_RA_MAX_A, ",
                                          "COMP_RA_MED_A, ",
                                          "COMP_RA_SD_A, ",
                                          "COMP_RA_SKEW_A, ",
                                          "COMP_RA_MAX_B, ",
                                          "COMP_RA_MED_B, ",
                                          "COMP_RA_SD_B, ",
                                          "COMP_RA_SKEW_B, ",
                                          "COMP_C_ENT_A, ",
                                          "COMP_C_ENT_B, ",
                                          "FW_RA_MAX_D, ",
                                          "FW_RA_MED_D, ",
                                          "FW_RA_SD_D, ",
                                          "FW_RA_PPOS_D, ",
                                          "FW_RA_PNEG_D, ",
                                          "FW_RA_PFC05_D, ",
                                          "FW_RA_PFC1_D, ",
                                          "FW_RA_PFC2_D, ",
                                          "FW_LOG_MAX_D, ",
                                          "FW_LOG_MED_D, ",
                                          "FW_LOG_SD_D, ",
                                          "FW_LOG_PPOS_D, ",
                                          "FW_LOG_PNEG_D, ",
                                          "FW_LOG_PFC05_D, ",
                                          "FW_LOG_PFC1_D, ",
                                          "FW_LOG_PFC2_D, ",
                                          "FW_CLR_MAX_D, ",
                                          "FW_CLR_MED_D, ",
                                          "FW_CLR_SD_D, ",
                                          "FW_CLR_PPOS_D, ",
                                          "FW_CLR_PNEG_D, ",
                                          "FW_CLR_PFC05_D, ",
                                          "FW_CLR_PFC1_D, ",
                                          "FW_CLR_PFC2_D) ",
                                          "VALUES(",
                                          "'", job$uuid, "', ",job$partial_info,", ",
                                          "'", job$type, "', ",
                                          job$TOTALS_C_FC, ", ",
                                          job$TOTALS_C_D, ", ",
                                          job$TOTALS_C_MAX_D, ", ",
                                          job$TOTALS_C_MED_D, ", ",
                                          job$TOTALS_C_SD_D, ", ",
                                          job$CORR_RA_MED, ", ",
                                          job$CORR_RA_SD, ", ",
                                          job$CORR_RA_SKEW, ", ",
                                          job$CORR_LOG_MED, ", ",
                                          job$CORR_LOG_SD, ", ",
                                          job$CORR_LOG_SKEW, ", ",
                                          job$CORR_CLR_MED, ", ",
                                          job$CORR_CLR_SD, ", ",
                                          job$CORR_CLR_SKEW, ", ",
                                          job$COMP_C_P0_A, ", ",
                                          job$COMP_C_P0_B, ", ",
                                          job$COMP_C_P1_A, ", ",
                                          job$COMP_C_P1_B, ", ",
                                          job$COMP_C_P5_A, ", ",
                                          job$COMP_C_P5_B, ", ",
                                          job$COMP_RA_P01_A, ", ",
                                          job$COMP_RA_P01_B, ", ",
                                          job$COMP_RA_P1_A, ", ",
                                          job$COMP_RA_P1_B, ", ",
                                          job$COMP_RA_P5_A, ", ",
                                          job$COMP_RA_P5_B, ", ",
                                          job$COMP_RA_MAX_A, ", ",
                                          job$COMP_RA_MED_A, ", ",
                                          job$COMP_RA_SD_A, ", ",
                                          job$COMP_RA_SKEW_A, ", ",
                                          job$COMP_RA_MAX_B, ", ",
                                          job$COMP_RA_MED_B, ", ",
                                          job$COMP_RA_SD_B, ", ",
                                          job$COMP_RA_SKEW_B, ", ",
                                          job$COMP_C_ENT_A, ", ",
                                          job$COMP_C_ENT_B, ", ",
                                          job$FW_RA_MAX_D, ", ",
                                          job$FW_RA_MED_D, ", ",
                                          job$FW_RA_SD_D, ", ",
                                          job$FW_RA_PPOS_D, ", ",
                                          job$FW_RA_PNEG_D, ", ",
                                          job$FW_RA_PFC05_D, ", ",
                                          job$FW_RA_PFC1_D, ", ",
                                          job$FW_RA_PFC2_D, ", ",
                                          job$FW_LOG_MAX_D, ", ",
                                          job$FW_LOG_MED_D, ", ",
                                          job$FW_LOG_SD_D, ", ",
                                          job$FW_LOG_PPOS_D, ", ",
                                          job$FW_LOG_PNEG_D, ", ",
                                          job$FW_LOG_PFC05_D, ", ",
                                          job$FW_LOG_PFC1_D, ", ",
                                          job$FW_LOG_PFC2_D, ", ",
                                          job$FW_CLR_MAX_D, ", ",
                                          job$FW_CLR_MED_D, ", ",
                                          job$FW_CLR_SD_D, ", ",
                                          job$FW_CLR_PPOS_D, ", ",
                                          job$FW_CLR_PNEG_D, ", ",
                                          job$FW_CLR_PFC05_D, ", ",
                                          job$FW_CLR_PFC1_D, ", ",
                                          job$FW_CLR_PFC2_D,
                                          ");"))
  }
  cat(paste0("Succeeded on ", updates, " rows\n"))
}

dbDisconnect(conn)
