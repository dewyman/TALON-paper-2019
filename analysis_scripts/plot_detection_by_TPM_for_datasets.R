
main <-function() {

    set.seed(100)
    load_packages()
    opt <- parse_options()

    # Get colors
    if (opt$color_scheme == "red") {
        color_vec <- c("white", "orange", "red2")
    } else if (opt$color_scheme == "blue") {
        color_vec <- c("white", "skyblue", "navy")
    } else if (opt$color_scheme == "green") {
        color_vec <- c("white", "yellowgreen", "springgreen4")
    }

    # Get genes and transcripts expressed in the Illumina data from the Kallisto
    # abundance file
    illumina_1 <- filter_kallisto_illumina_genes(opt$illumina_kallisto_1)
    illumina_2 <- filter_kallisto_illumina_genes(opt$illumina_kallisto_2)
    colnames(illumina_1) <- c("gene", "illumina_TPM_1")
    colnames(illumina_2) <- c("gene", "illumina_TPM_2")
    illumina_gene_table <- merge(illumina_1, illumina_2, by = "gene",
                                 all.x = T, all.y = T)

    # Remove genes only detected in one Illumina rep, then average the TPMS
    illumina_gene_table <- illumina_gene_table[complete.cases(illumina_gene_table), ]
    illumina_gene_table$tpm <- (illumina_gene_table$illumina_TPM_1 +
                                illumina_gene_table$illumina_TPM_2) / 2

    
    # Add a column to the Illumina table mapping the names to TALON IDs
    gene_name_mapping <- get_name_table(opt$database, "gene") 
    
    illumina_gene_table_with_IDs <- merge(illumina_gene_table, gene_name_mapping, by = "gene", all.x = T, all.y = F)

    # Get the names of the first and second dataset that we will be working with
    data_names <- str_split(opt$datasets, ",")[[1]]
    dataset_1 <- data_names[1]
    dataset_2 <- data_names[2]
    
    # Read in whitelist file, and get gene and transcript whitelists from that
    if (is.null(opt$whitelist)) {
        whitelisted_gene_IDs = NULL
        #whitelisted_transcript_IDs = NULL
    } else {
        whitelist <- as.data.frame(read_delim(opt$whitelist, ",", escape_double = FALSE,
                                  col_names = FALSE, trim_ws = TRUE, na = "NA"))
        whitelisted_gene_IDs <- unique(whitelist[,1])
    }

    # Now get filtered genes/transcripts detected in the datasets
    d1_genes <- get_detected_genes_for_dataset(dataset_1, whitelisted_gene_IDs, opt$database)
    d2_genes <- get_detected_genes_for_dataset(dataset_2, whitelisted_gene_IDs, opt$database)


    # Combine the Illumina tables with information about which genes/transcripts are observed in Pacbio
    illumina_gene_detection <- get_detection(illumina_gene_table_with_IDs, d1_genes, d2_genes, "gene_ID")

    # Group into buckets
    illumina_gene_detection_buckets <- get_buckets(illumina_gene_detection)    
    #illumina_transcript_detection_buckets <- get_buckets(illumina_transcript_detection)

    # Plot detection by TPM
    plot_detection(illumina_gene_detection_buckets$illumina, illumina_gene_detection_buckets$interval_labels, "gene", color_vec, opt$outdir)
    print ("--------------------------")



   #plot_length_hists_by_detection(highly_expressed, opt$outdir)
}

get_buckets <- function(illumina) {
    intervals = c(1, 2, seq(5, 20, 5), 50, 100, 500)
    illumina$group = cut(illumina$tpm, c(intervals, 100000000000))   

    cIntervals = as.character(intervals)
    interval_starts_with_dash = paste0(cIntervals, c(rep("-", length(cIntervals)-1), "+"))
    full_intervals = paste0(interval_starts_with_dash, c(as.character(intervals)[2:length(intervals)], ""))

    return(list("illumina" = illumina, "interval_labels" = full_intervals))
}


get_detection <- function(illumina, d1_items, d2_items, cat_type) {

    shared <- intersect(d1_items, d2_items)
    d1_d2_union <- union(d1_items, d2_items) 
    not_shared <- d1_d2_union[(d1_d2_union %in% shared) == F]

    illumina$detection = ifelse(illumina[,cat_type] %in% not_shared, "Detected in one PacBio rep",
                                             ifelse(illumina[,cat_type] %in% shared, "Detected in both PacBio reps", "Not detected in PacBio"))
    return(illumina)
}

get_name_table <- function(database, cat_type) {

    # Connect to the database
    con <- dbConnect(SQLite(), dbname=database)
   
    # Fetch observed starts or ends
    q_text <- paste("SELECT ID,value FROM ", cat_type, "_annotations WHERE attribute = '", cat_type, "_name'", sep="")
    
    query <- dbSendQuery(con, q_text)
    name_table <- as.data.frame(dbFetch(query, n = -1))
    dbClearResult(query)
    dbDisconnect(con)
    colnames(name_table) <- c(paste(cat_type, "ID", sep="_"), cat_type)

    return(name_table)
}

plot_detection <- function(illumina, cIntervals, cat_type, color_vec, outdir) {

    # Plot the curves
    fname <- paste(outdir, "/", cat_type, "_detection_by_TPM.png", sep="")
    xlabel <- paste(capitalize(cat_type), "expression level in Illumina data (TPM)")
    ylabel <- paste("Fraction of ", cat_type, "s", sep="")

    png(filename = fname,
        width = 2100, height = 2500, units = "px",
        bg = "white",  res = 300)
    g = ggplot(illumina, aes(x = group, fill = factor(detection, levels = c("Not detected in PacBio", "Detected in one PacBio rep", "Detected in both PacBio reps")))) +
       geom_bar(position = "fill", col = "black") + 
       scale_fill_manual("",values=color_vec)  + 
       theme_bw() + 
       scale_y_continuous(breaks=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) + 
       theme(axis.text.x = element_text(angle = 30, hjust = 1, color = "black", size=20)) + 
       xlab(xlabel)  + ylab(ylabel) + theme(text= element_text(size=22)) + 
       theme(axis.text.x = element_text(color = "black", size=21),
             axis.text.y = element_text(color = "black", size=21)) + 
       scale_x_discrete(labels=cIntervals) + 
       theme(legend.position=c(0.6,0.2),
              legend.title = element_blank(),
              legend.background = element_rect(fill="white", color = "black"),
              legend.key = element_rect(fill="transparent"),
              legend.text = element_text(colour = 'black', size = 20),
              plot.margin = margin(t = 0.5, r = 1.5, l = 0.5, b = 0.5, "cm"))

    print(g)
    dev.off()
    if (cat_type == "gene") {
        print(subset(illumina, (group == "(100,500]" | group == "(500,1e+11]") &
                            detection != "Detected in both PacBio reps"))
    
    #print(unique(as.character(subset(illumina, group == "(500,1e+11]" & detection == "Not detected in PacBio")$gene)))
    #print(unique(as.character(subset(illumina, group == "(50,100]" & detection == "Not detected in PacBio")$gene)))

 
    #print(unique(as.character(subset(illumina, (group == "(500,1e+11]" | group == "(100,500]")
    #                                            & detection == "Not detected in PacBio")$gene)))
    }
}

plot_length_hists_by_detection <- function(illumina, outdir) {

    fname <- paste(outdir, "/transcript_length_by_detection_status_100plusTPM.png", sep="")
    xlabel <- "Transcript length (kilobases)"
    ylabel <- "Count"
    illumina$length <- illumina$length/1000

    png(filename = fname,
        width = 3000, height = 2000, units = "px",
        bg = "white",  res = 300)
    g = ggplot(illumina, aes(x = length, color = factor(detection, levels = c("Not detected in PacBio", "Detected in one PacBio rep", "Detected in both PacBio reps")))) +
       geom_freqpoly(binwidth = 0.1) + 
       scale_colour_manual("",values=c("red", "skyblue", "navy"))  +
       theme_bw() + 
       theme(axis.text.x = element_text(color = "black", size=17)) +
       xlab(xlabel)  + ylab(ylabel) + theme(text= element_text(size=17)) +
       theme(axis.text.y = element_text(color = "black", size=17)) + 
       coord_cartesian(xlim=c(0, 7), ylim=c(0,60)) + 
       scale_x_continuous(breaks=seq(0,7, by = 1))

    print(g)
    dev.off()
}

get_detected_genes_for_dataset <- function(dataset, whitelisted_genes, database) {
    # Connect to the database
    con <- dbConnect(SQLite(), dbname=database)

    # Fetch observed gene IDs
    query <- dbSendQuery(con, paste("SELECT gene_ID FROM abundance LEFT JOIN transcripts ON transcripts.transcript_ID = abundance.transcript_ID WHERE dataset = '", dataset, "'", sep=""))
    geneIDs <- unique(as.data.frame(dbFetch(query, n = -1))[,1])

    dbClearResult(query)
    dbDisconnect(con)

    if (!(is.null(whitelisted_genes))) {
        filtered_geneIDs <- geneIDs[geneIDs %in% whitelisted_genes]
        return(filtered_geneIDs)
    } else {
        return(geneIDs)
    }
}

get_detected_transcripts_for_dataset <- function(dataset, whitelisted_transcripts, database) {
    # Connect to the database
    con <- dbConnect(SQLite(), dbname=database)

    # Fetch observed transcript IDs
    query <- dbSendQuery(con, paste("SELECT transcript_ID FROM abundance WHERE dataset = '", dataset, "'", sep=""))
    transcriptIDs <- as.data.frame(dbFetch(query, n = -1))[,1]

    dbClearResult(query)
    dbDisconnect(con)

    # Filter transcripts if a whetelist was provided
    if (!(is.null(whitelisted_transcripts))) {
        filtered_transcriptIDs <- transcriptIDs[transcriptIDs %in% whitelisted_transcripts]
        return(filtered_transcriptIDs)
    } else {
        return(transcriptIDs)
    }
}

load_packages <- function() {
    suppressPackageStartupMessages(library("DBI"))
    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("plyr"))
    suppressPackageStartupMessages(library("Hmisc"))
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("RSQLite"))
    suppressPackageStartupMessages(library("readr"))
    suppressPackageStartupMessages(library("reshape"))
    suppressPackageStartupMessages(library("stringr"))
    suppressPackageStartupMessages(library("data.table"))

    # Load my custom functions
    # source("/data/users/freese/mortazavi_lab/bin/TALON-paper-2019/analysis_scripts/filter_kallisto_illumina_genes.R")
    # source("/data/users/freese/mortazavi_lab/bin/TALON-paper-2019/analysis_scripts/filter_kallisto_illumina_transcripts.R")
    # source("/data/users/freese/mortazavi_lab/bin/TALON-paper-2019/analysis_scripts/get_database_transcript_table.R")

    source("/dfs2/pub/dwyman/TALON-paper-2019/analysis_scripts/filter_kallisto_illumina_genes.R")
    source("/dfs2/pub/dwyman/TALON-paper-2019/analysis_scripts/filter_kallisto_illumina_transcripts.R")
    source("/dfs2/pub/dwyman/TALON-paper-2019/analysis_scripts/get_database_transcript_table.R")
    
    return
}


parse_options <- function() {

    option_list <- list(
        make_option(c("--db"), action = "store", dest = "database",
                    default = NULL, help = "TALON database file"),
        make_option(c("--whitelist"), action = "store", dest = "whitelist",
                    default = NULL, help = "File of whitelisted transcripts for the Pacbio data"),
        make_option(c("--datasets"), action = "store", dest = "datasets",
                    default = NULL, help = "Comma-delimited list of two dataset names to include in the analysis."),
        make_option(c("--ik1"), action = "store", dest = "illumina_kallisto_1",
                    default = NULL, help = "Rep1 Illumina Kallisto file."),
        make_option(c("--ik2"), action = "store", dest = "illumina_kallisto_2",
                    default = NULL, help = "Rep2 Illumina Kallisto file."),
        make_option(c("--color"), action = "store", dest = "color_scheme",
                    default = NULL, help = "blue, red, or green"),
        make_option(c("-o","--outdir"), action = "store", dest = "outdir",
                    default = NULL, help = "Output directory for plots and outfiles")
        )

    opt <- parse_args(OptionParser(option_list=option_list))
    return(opt)
}


main()
