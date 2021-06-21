source("path_fix.R")

library(RSQLite)

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))
# Clear FIT_IN_PROGRESS flag
discard <- dbExecute(conn, "UPDATE results SET FIT_IN_PROGRESS = 0 WHERE TRUE;")
# Clear NA results
discard <- dbExecute(conn, "DELETE FROM results WHERE RESULT IS NULL;")
#discard <- dbExecute(conn, "DELETE FROM characteristics;")
# dbGetQuery(conn, "SELECT * FROM results;")
dbDisconnect(conn)
