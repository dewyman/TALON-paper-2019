
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
        color_vec <- c("white", "olivedrab3", "springgreen4")
    }

    # Now, compute mean, median, and max transcript lengths for each gene
    gene_length_table <- compute_lengths(opt$illumina_kallisto_1, opt$illumina_kallisto_2)

    # Get genes expressed in the Illumina data from the Kallisto
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

    # Now get filtered genes detected in the datasets
    d1_genes <- get_detected_genes_for_dataset(dataset_1, whitelisted_gene_IDs, opt$database)
    d2_genes <- get_detected_genes_for_dataset(dataset_2, whitelisted_gene_IDs, opt$database)


    # Combine the Illumina tables with information about which genes are observed in Pacbio
    illumina_gene_detection <- get_detection(illumina_gene_table_with_IDs, d1_genes, d2_genes, "gene_ID")

    # Group into buckets
    illumina_gene_detection_buckets <- get_buckets(illumina_gene_detection)    

    # Plot length by detection
    detection_with_lengths <- merge(illumina_gene_detection_buckets$illumina, 
                                    gene_length_table, by = "gene", all.x = T,
                                    all.y = F)
    plot_length_by_detection(detection_with_lengths, 
                             illumina_gene_detection_buckets$interval_labels, 
                             "Median", color_vec, options$outdir)
}

plot_length_by_detection  <- function(data, cIntervals, len_type, colors, outdir) {
    data <- data[data$length_type == "Median",]
    #fname <- paste(outdir, "/length_by_detection_and_TPM_", len_type, ".png", sep="")
    #xlabel <- "Gene expression level in Illumina data (TPM)"
    #ylabel <- paste("Gene length (", len_type, " length of known transcripts)", sep="") 

    #png(filename = fname,
    #    width = 3500, height = 2500, units = "px",
    #    bg = "white",  res = 300)

    quit()
    g = ggplot(data, aes(x = group, y = tpm, fill = detection)) +
            geom_violin(alpha = 0.8) +
            geom_boxplot(width=0.1, fill="white", outlier.size=-1) +
            genom_point(size = 0.5) +
            xlab(xlabel) + ylab(ylabel) +
            theme_bw(base_family = "Helvetica", base_size = 18) +
            scale_fill_manual("", values = color_vec) +
            theme_bw(base_family = "Helvetica", base_size = 18) +
            theme(axis.line.x = element_line(color="black", size = 0.5),
                  axis.line.y = element_line(color="black", size = 0.5),
                  axis.text.x = element_text(color="black", size = rel(1.5)),
                  axis.text.y = element_text(color="black", size = rel(1.5)),
                  axis.title.x = element_text(color="black", size=rel(1.25)),
                  axis.title.y = element_text(color="black", size=rel(1.25))) +
            scale_x_discrete(labels=cIntervals)  
    print(g)
    dev.off()
}

compute_lengths <- function(ik1, ik2) {

    ik1_tab <- read_ik_file(ik1)[,c("transcript", "gene", "length")]
    ik2_tab <- read_ik_file(ik1)[,c("transcript", "gene", "length")]  
 
    all_ik <- merge(ik1_tab, ik2_tab, by = "transcript", all.x = T, all.y = T)
    all_ik$length <- (all_ik$length.x + all_ik$length.y)/2
        
    # Aggregate by gene: mean length
    print(head(all_ik))
    mean_gene_lens <- aggregate(all_ik$length, by=list(all_ik$gene.x), FUN=mean)
    print(head(mean_gene_lens))
    colnames(mean_gene_lens) <- c("gene", "length")
    mean_gene_lens$length_type <- "Mean"

    # Aggregate by gene: median length
    med_gene_lens <- aggregate(all_ik$length, by=list(all_ik$gene.x), FUN=median)    
    colnames(med_gene_lens) <- c("gene", "length")
    med_gene_lens$length_type <- "Median"

    # Aggregate by gene: max length
    max_gene_lens <- aggregate(all_ik$length, by=list(all_ik$gene.x), FUN=max)
    colnames(max_gene_lens) <- c("gene", "length")
    max_gene_lens$length_type <- "Max"

    # Combine into one table
    lengths <- rbind(mean_gene_lens, rbind(med_gene_lens, max_gene_lens)) 
    return(lengths)
}

read_ik_file <- function(kallisto_file) {
    gencode.quantitation <- as.data.frame(read_delim(kallisto_file, "\t", escape_double = FALSE,
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))
    extraCols =str_split_fixed(gencode.quantitation$target_id, "\\|",9)[,c(5,6,8,1,2)]
    colnames(extraCols) = c("transcript", "gene", "class", "t_ID", "g_ID")
    gencode.quantitation = cbind(extraCols, gencode.quantitation)

    return(gencode.quantitation)

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
    source("/dfs2/pub/dwyman/TALON-paper-2019/analysis_scripts/filter_kallisto_illumina_genes.R")
    source("/dfs2/pub/dwyman/TALON-paper-2019/analysis_scripts/filter_kallisto_illumina_transcripts.R")
    source("/pub/dwyman/TALON-paper-2019/analysis_scripts/get_database_transcript_table.R")
    
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
