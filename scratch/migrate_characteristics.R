source("path_fix.R")

library(codaDE)
library(RSQLite)

# ------------------------------------------------------------------------------
#   Open new and old connections
# ------------------------------------------------------------------------------

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))

discard <- dbExecute(conn, "ALTER TABLE characteristics RENAME TO characteristics2;")

discard <- dbExecute(conn, paste0("CREATE TABLE characteristics(",
                                      "UUID VARCHAR(36),",
                                      "PARTIAL INTEGER,",
                                      "TOTALS_C_FC REAL,",
                                      "TOTALS_C_D REAL,",
                                      "TOTALS_C_MAX_D REAL,",
                                      "TOTALS_C_MED_D REAL,",
                                      "TOTALS_C_SD_D REAL,",
                                      "CORR_RA_MED REAL,",
                                      "CORR_RA_SD REAL,",
                                      "CORR_RA_SKEW REAL,",
                                      "CORR_LOG_MED REAL,",
                                      "CORR_LOG_SD REAL,",
                                      "CORR_LOG_SKEW REAL,",
                                      "CORR_CLR_MED REAL,",
                                      "CORR_CLR_SD REAL,",
                                      "CORR_CLR_SKEW REAL,",
                                      "COMP_C_P0_A REAL,",
                                      "COMP_C_P0_B REAL,",
                                      "COMP_C_P1_A REAL,",
                                      "COMP_C_P1_B REAL,",
                                      "COMP_C_P5_A REAL,",
                                      "COMP_C_P5_B REAL,",
                                      "COMP_RA_P01_A REAL,",
                                      "COMP_RA_P01_B REAL,",
                                      "COMP_RA_P1_A REAL,",
                                      "COMP_RA_P1_B REAL,",
                                      "COMP_RA_P5_A REAL,",
                                      "COMP_RA_P5_B REAL,",
                                      "COMP_RA_MAX_A REAL,",
                                      "COMP_RA_MED_A REAL,",
                                      "COMP_RA_SD_A REAL,",
                                      "COMP_RA_SKEW_A REAL,",
                                      "COMP_RA_MAX_B REAL,",
                                      "COMP_RA_MED_B REAL,",
                                      "COMP_RA_SD_B REAL,",
                                      "COMP_RA_SKEW_B REAL,",
                                      "COMP_C_ENT_A REAL,",
                                      "COMP_C_ENT_B REAL,",
                                      "FW_RA_MAX_D REAL,",
                                      "FW_RA_MED_D REAL,",
                                      "FW_RA_SD_D REAL,",
                                      "FW_RA_PPOS_D REAL,",
                                      "FW_RA_PNEG_D REAL,",
                                      "FW_RA_PFC05_D REAL,",
                                      "FW_RA_PFC1_D REAL,",
                                      "FW_RA_PFC2_D REAL,",
                                      "FW_LOG_MAX_D REAL,",
                                      "FW_LOG_MED_D REAL,",
                                      "FW_LOG_SD_D REAL,",
                                      "FW_LOG_PPOS_D REAL,",
                                      "FW_LOG_PNEG_D REAL,",
                                      "FW_LOG_PFC05_D REAL,",
                                      "FW_LOG_PFC1_D REAL,",
                                      "FW_LOG_PFC2_D REAL,",
                                      "FW_CLR_MAX_D REAL,",
                                      "FW_CLR_MED_D REAL,",
                                      "FW_CLR_SD_D REAL,",
                                      "FW_CLR_PPOS_D REAL,",
                                      "FW_CLR_PNEG_D REAL,",
                                      "FW_CLR_PFC05_D REAL,",
                                      "FW_CLR_PFC1_D REAL,",
                                      "FW_CLR_PFC2_D REAL,",
                                      "PRIMARY KEY(UUID, PARTIAL));"))

dbDisconnect(conn)
