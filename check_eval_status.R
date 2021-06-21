source("path_fix.R")

library(RSQLite)
library(optparse)

ps <- c(100, 1000, 5000)
methods <- c("ALDEx2", "DESeq2", "MAST", "NBGLM", "scran")

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))
corrp <- 1

for(p in ps) {
  for(method in methods) {
    res <- dbGetQuery(conn, paste0("SELECT * FROM results LEFT JOIN datasets ON results.UUID=datasets.UUID WHERE P=",p," AND CORRP=",corrp," AND method='",method,"'"))
    finished <- sum(!is.na(res$RESULT))/1
    in_progress <- sum(is.na(res$RESULT) & res$FIT_IN_PROGRESS == 1)/1
    not_started <- 400/1 - finished - in_progress
    cat(paste0("P = ", p, ", METHOD = ", method, ", FINISHED: ", finished, ", IN PROGRESS: ", in_progress, ", NOT STARTED: ", not_started, "\n"))
  }
}
dbDisconnect(conn)
