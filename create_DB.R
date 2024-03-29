source("path_fix.R")

library(codaDE)
library(RSQLite)

# Create DB (if doesn't exist?)
conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))

#discard <- dbExecute(conn, paste0("CREATE TABLE da_realized(",
#                                  "UUID VARCHAR(36), ",
#                                  "MEASURED_BY VARCHAR(64), ",
#                                  "CALLS VARCHAR(100000), ",
#                                  "PRIMARY KEY (UUID, MEASURED_BY));"))

#discard <- dbExecute(conn, paste0("CREATE TABLE datasets(",
#                                  "UUID VARCHAR(36) PRIMARY KEY, ",
#                                  "P INTEGER, ",
#                                  "CORRP INTEGER, ",
#                                  "LOG_MEAN REAL, ",
#                                  "PERTURBATION REAL, ",
#                                  "REP_NOISE REAL, ",
#                                  "FC_ABSOLUTE REAL, ",
#                                  "FC_RELATIVE REAL, ",
#                                  "FC_PARTIAL REAL, ",
#                                  "BASELINE_CALLS VARCHAR(100000), ",
#                                  "BASELINE_BETAS VARCHAR(100000), ",
#                                  "MED_ABS_TOTAL REAL, ",
#                                  "MED_REL_TOTAL REAL, ",
#                                  "PERCENT_DIFF_SIM REAL, ",
#                                  "PERCENT_DIFF_REALIZ REAL)"))

#discard <- dbExecute(conn, paste0("INSERT INTO datasets2 SELECT UUID, P, CORRP, LOG_MEAN, PERTURBATION, REP_NOISE, ",
#                                          "FC_ABSOLUTE, FC_RELATIVE, FC_PARTIAL, BASELINE_CALLS, NULL, MED_ABS_TOTAL, ",
#                                          "MED_REL_TOTAL, PERCENT_DIFF_SIM, PERCENT_DIFF_REALIZ FROM datasets;"))

#discard <- dbExecute(conn, "DROP TABLE datasets")

#discard <- dbExecute(conn, "ALTER TABLE datasets2 RENAME TO datasets")

#discard <- dbExecute(conn, "ALTER TABLE results RENAME TO results_backup")

#discard <- dbExecute(conn, paste0("CREATE TABLE results(",
#                                  "UUID VARCHAR(36), ",
#                                  "METHOD VARCHAR(64), ",
#                                  "PARTIAL_INFO INT DEFAULT 0 NOT NULL, ",
#                                  "BASELINE_TYPE VARCHAR(16), ",
#                                  "OBSERVED_TYPE VARCHAR(64), ",
#                                  "BASELINE_CALLS VARCHAR(100000), ",
#                                  "CALLS VARCHAR(100000), ",
#                                  "BASELINE_BETAS VARCHAR(100000), ",
#                                  "BETAS VARCHAR(100000), ",
#                                  "TPR REAL, ",
#                                  "FPR REAL, ",
#                                  "FDR REAL, ",
#                                  "BETA REAL DEFAULT -1 NOT NULL, ",
#                                  "PRIMARY KEY (UUID, METHOD, PARTIAL_INFO, BASELINE_TYPE, OBSERVED_TYPE, FDR, BETA));"))

#discard <- dbExecute(conn, paste0("INSERT INTO results2 SELECT UUID, METHOD, PARTIAL_INFO, BASELINE_TYPE, OBSERVED_TYPE, ",
#                                          "BASELINE_CALLS, CALLS, NULL, NULL, TPR, FPR, FDR, BETA FROM results;"))

#discard <- dbExecute(conn, "DROP TABLE results")

#discard <- dbExecute(conn, "ALTER TABLE results2 RENAME TO results")

#discard <- dbExecute(conn, paste0("CREATE TABLE characteristics(",
#                                  "UUID VARCHAR(36),",
#                                  "PARTIAL INTEGER,",
#                                  "TYPE VARCHAR(64),",
#                                  "TOTALS_C_FC REAL,",
#                                  "TOTALS_C_D REAL,",
#                                  "TOTALS_C_MAX_D REAL,",
#                                  "TOTALS_C_MED_D REAL,",
#                                  "TOTALS_C_SD_D REAL,",
#                                  "CORR_RA_MED REAL,",
#                                  "CORR_RA_SD REAL,",
#                                  "CORR_RA_SKEW REAL,",
#                                  "CORR_LOG_MED REAL,",
#                                  "CORR_LOG_SD REAL,",
#                                  "CORR_LOG_SKEW REAL,",
#                                  "CORR_CLR_MED REAL,",
#                                  "CORR_CLR_SD REAL,",
#                                  "CORR_CLR_SKEW REAL,",
#                                  "COMP_C_P0_A REAL,",
#                                  "COMP_C_P0_B REAL,",
#                                  "COMP_C_P1_A REAL,",
#                                  "COMP_C_P1_B REAL,",
#                                  "COMP_C_P5_A REAL,",
#                                  "COMP_C_P5_B REAL,",
#                                  "COMP_RA_P01_A REAL,",
#                                  "COMP_RA_P01_B REAL,",
#                                  "COMP_RA_P1_A REAL,",
#                                  "COMP_RA_P1_B REAL,",
#                                  "COMP_RA_P5_A REAL,",
#                                  "COMP_RA_P5_B REAL,",
#                                  "COMP_RA_MAX_A REAL,",
#                                  "COMP_RA_MED_A REAL,",
#                                  "COMP_RA_SD_A REAL,",
#                                  "COMP_RA_SKEW_A REAL,",
#                                  "COMP_RA_MAX_B REAL,",
#                                  "COMP_RA_MED_B REAL,",
#                                  "COMP_RA_SD_B REAL,",
#                                  "COMP_RA_SKEW_B REAL,",
#                                  "COMP_C_ENT_A REAL,",
#                                  "COMP_C_ENT_B REAL,",
#                                  "FW_RA_MAX_D REAL,",
#                                  "FW_RA_MED_D REAL,",
#                                  "FW_RA_SD_D REAL,",
#                                  "FW_RA_PPOS_D REAL,",
#                                  "FW_RA_PNEG_D REAL,",
#                                  "FW_RA_PFC05_D REAL,",
#                                  "FW_RA_PFC1_D REAL,",
#                                  "FW_RA_PFC2_D REAL,",
#                                  "FW_LOG_MAX_D REAL,",
#                                  "FW_LOG_MED_D REAL,",
#                                  "FW_LOG_SD_D REAL,",
#                                  "FW_LOG_PPOS_D REAL,",
#                                  "FW_LOG_PNEG_D REAL,",
#                                  "FW_LOG_PFC05_D REAL,",
#                                  "FW_LOG_PFC1_D REAL,",
#                                  "FW_LOG_PFC2_D REAL,",
#                                  "FW_CLR_MAX_D REAL,",
#                                  "FW_CLR_MED_D REAL,",
#                                  "FW_CLR_SD_D REAL,",
#                                  "FW_CLR_PPOS_D REAL,",
#                                  "FW_CLR_PNEG_D REAL,",
#                                  "FW_CLR_PFC05_D REAL,",
#                                  "FW_CLR_PFC1_D REAL,",
#                                  "FW_CLR_PFC2_D REAL,",
#                                  "PRIMARY KEY(UUID, PARTIAL, TYPE));"))

# discard <- dbExecute(conn, paste0("INSERT INTO characteristics_new(",
#                                  "TYPE,",
#                                  "UUID,",
#                                  "PARTIAL,",
#                                  "TOTALS_C_FC,",
#                                  "TOTALS_C_D,",
#                                  "TOTALS_C_MAX_D,",
#                                  "TOTALS_C_MED_D,",
#                                  "TOTALS_C_SD_D,",
#                                  "CORR_RA_MED,",
#                                  "CORR_RA_SD,",
#                                  "CORR_RA_SKEW,",
#                                  "CORR_LOG_MED,",
#                                  "CORR_LOG_SD,",
#                                  "CORR_LOG_SKEW,",
#                                  "CORR_CLR_MED,",
#                                  "CORR_CLR_SD,",
#                                  "CORR_CLR_SKEW,",
#                                  "COMP_C_P0_A,",
#                                  "COMP_C_P0_B,",
#                                  "COMP_C_P1_A,",
#                                  "COMP_C_P1_B,",
#                                  "COMP_C_P5_A,",
#                                  "COMP_C_P5_B,",
#                                  "COMP_RA_P01_A,",
#                                  "COMP_RA_P01_B,",
#                                  "COMP_RA_P1_A,",
#                                  "COMP_RA_P1_B,",
#                                  "COMP_RA_P5_A,",
#                                  "COMP_RA_P5_B,",
#                                  "COMP_RA_MAX_A,",
#                                  "COMP_RA_MED_A,",
#                                  "COMP_RA_SD_A,",
#                                  "COMP_RA_SKEW_A,",
#                                  "COMP_RA_MAX_B,",
#                                  "COMP_RA_MED_B,",
#                                  "COMP_RA_SD_B,",
#                                  "COMP_RA_SKEW_B,",
#                                  "COMP_C_ENT_A,",
#                                  "COMP_C_ENT_B,",
#                                  "FW_RA_MAX_D,",
#                                  "FW_RA_MED_D,",
#                                  "FW_RA_SD_D,",
#                                  "FW_RA_PPOS_D,",
#                                  "FW_RA_PNEG_D,",
#                                  "FW_RA_PFC05_D,",
#                                  "FW_RA_PFC1_D,",
#                                  "FW_RA_PFC2_D,",
#                                  "FW_LOG_MAX_D,",
#                                  "FW_LOG_MED_D,",
#                                  "FW_LOG_SD_D,",
#                                  "FW_LOG_PPOS_D,",
#                                  "FW_LOG_PNEG_D,",
#                                  "FW_LOG_PFC05_D,",
#                                  "FW_LOG_PFC1_D,",
#                                  "FW_LOG_PFC2_D,",
#                                  "FW_CLR_MAX_D,",
#                                  "FW_CLR_MED_D,",
#                                  "FW_CLR_SD_D,",
#                                  "FW_CLR_PPOS_D,",
#                                  "FW_CLR_PNEG_D,",
#                                  "FW_CLR_PFC05_D,",
#                                  "FW_CLR_PFC1_D,",
#                                  "FW_CLR_PFC2_D) ",
#                                  "SELECT 'relative_abundances',",
#                                  "UUID,",
#                                  "PARTIAL,",
#                                  "TOTALS_C_FC,",
#                                  "TOTALS_C_D,",
#                                  "TOTALS_C_MAX_D,",
#                                  "TOTALS_C_MED_D,",
#                                  "TOTALS_C_SD_D,",
#                                  "CORR_RA_MED,",
#                                  "CORR_RA_SD,",
#                                  "CORR_RA_SKEW,",
#                                  "CORR_LOG_MED,",
#                                  "CORR_LOG_SD,",
#                                  "CORR_LOG_SKEW,",
#                                  "CORR_CLR_MED,",
#                                  "CORR_CLR_SD,",
#                                  "CORR_CLR_SKEW,",
#                                  "COMP_C_P0_A,",
#                                  "COMP_C_P0_B,",
#                                  "COMP_C_P1_A,",
#                                  "COMP_C_P1_B,",
#                                  "COMP_C_P5_A,",
#                                  "COMP_C_P5_B,",
#                                  "COMP_RA_P01_A,",
#                                  "COMP_RA_P01_B,",
#                                  "COMP_RA_P1_A,",
#                                  "COMP_RA_P1_B,",
#                                  "COMP_RA_P5_A,",
#                                  "COMP_RA_P5_B,",
#                                  "COMP_RA_MAX_A,",
#                                  "COMP_RA_MED_A,",
#                                  "COMP_RA_SD_A,",
#                                  "COMP_RA_SKEW_A,",
#                                  "COMP_RA_MAX_B,",
#                                  "COMP_RA_MED_B,",
#                                  "COMP_RA_SD_B,",
#                                  "COMP_RA_SKEW_B,",
#                                  "COMP_C_ENT_A,",
#                                  "COMP_C_ENT_B,",
#                                  "FW_RA_MAX_D,",
#                                  "FW_RA_MED_D,",
#                                  "FW_RA_SD_D,",
#                                  "FW_RA_PPOS_D,",
#                                  "FW_RA_PNEG_D,",
#                                  "FW_RA_PFC05_D,",
#                                  "FW_RA_PFC1_D,",
#                                  "FW_RA_PFC2_D,",
#                                  "FW_LOG_MAX_D,",
#                                  "FW_LOG_MED_D,",
#                                  "FW_LOG_SD_D,",
#                                  "FW_LOG_PPOS_D,",
#                                  "FW_LOG_PNEG_D,",
#                                  "FW_LOG_PFC05_D,",
#                                  "FW_LOG_PFC1_D,",
#                                  "FW_LOG_PFC2_D,",
#                                  "FW_CLR_MAX_D,",
#                                  "FW_CLR_MED_D,",
#                                  "FW_CLR_SD_D,",
#                                  "FW_CLR_PPOS_D,",
#                                  "FW_CLR_PNEG_D,",
#                                  "FW_CLR_PFC05_D,",
#                                  "FW_CLR_PFC1_D,",
#                                  "FW_CLR_PFC2_D ",
#                                  "FROM characteristics"))

dbDisconnect(conn)
