get_database_transcript_table <- function(database) {

    # Connect to the database
    con <- dbConnect(SQLite(), dbname=database)

    # Fetch the table
    query <- dbSendQuery(con, "SELECT * FROM transcripts")
    transcript_table <- as.data.frame(dbFetch(query, n = -1))
    dbClearResult(query)
    dbDisconnect(con)

    return(transcript_table)
}
