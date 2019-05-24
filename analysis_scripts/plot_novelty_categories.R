main <-function() {

    load_packages()
    opt <- parse_options()
    database <- opt$database
    datasets <- opt$datasets
    outdir <- opt$outdir

    # Read in whitelist file, and get gene and transcript whitelists from that
    whitelist <- read_delim(opt$whitelist, ",", escape_double = FALSE,
                                  col_names = FALSE, trim_ws = TRUE, na = "NA")

    colnames(whitelist) <- c("gene_ID", "transcript_ID", "novelty")
    
    # Fetch observed transcripts, then filter to only those in the whitelist
    con <- dbConnect(SQLite(), dbname=database)
    query <- dbSendQuery(con, "SELECT read_name, gene_ID, transcript_ID, dataset FROM observed")
    transcript_table <- as.data.frame(dbFetch(query, n = -1))

    dbClearResult(query)
    dbDisconnect(con)

    # If datasets selected, filter the results to limit them to those datasets
    if (!is.null(datasets)) {
        datasets <- unlist(strsplit(datasets, ","))
        transcript_table <- subset(transcript_table, dataset %in% datasets)
    }

    transcript_table <- merge( x = transcript_table,
                               y = whitelist,
                               by = "transcript_ID",
                               all.x = F, all.y = F)
    
    plot_distinct_novelty(transcript_table, outdir)
    plot_novelty_on_reads(transcript_table, outdir, datasets)
}

plot_distinct_novelty <- function(observed_transcripts, outdir){
    # This function plots the number of whitelist items that belong to each 
    # transcript novelty class. So it amounts to the number of unique transcripts
    # of each type that were identified in each dataset

    # Filter table to get unique instances
    distinct_transcripts <- observed_transcripts[!duplicated(observed_transcripts[,c('transcript_ID', 'novelty')]),]
    distinct_transcripts$novelty <- as.factor(distinct_transcripts$novelty)

    # Remap levels to easier names
    distinct_transcripts$novelty = revalue(distinct_transcripts$novelty, c("FSM_transcript"="Known", 
                                            "ISM-prefix_transcript"="ISM",
                                            "ISM-suffix_transcript"="ISM",
                                            "other_ISM_transcript"="ISM",
                                            "NIC_transcript"= "NIC",
                                            "NNC_transcript"= "NNC",
                                            "antisense_transcript"= "Antisense",
                                            "intergenic_transcript"= "Intergenic"))
    distinct_transcripts$novelty <- factor(distinct_transcripts$novelty, levels = c("Known", "ISM", 
                                                                                    "NIC", "NNC", 
                                                                                    "Antisense", "Intergenic"))    
    
    # Plotting
    fname <- paste(outdir, "/distinct_isoforms_by_category.png", sep="")
    xlabel <- "Isoform category"
    ylabel <- "Number of distinct isoforms"
    ymax <- 1.02*(max(count(distinct_transcripts,"novelty")$freq))

    colors <- c("Known" = "#009E73","ISM" = "#0072B2", "NIC" = "#D55E00", 
                "NNC" = "#E69F00", "Antisense" = "#000000", 
                "Intergenic" = "#CC79A7")

    png(filename = fname,
        width = 3100, height = 3500, units = "px",
        bg = "white",  res = 300)
    g = ggplot(distinct_transcripts, aes(x = novelty, width=.6,
               fill = as.factor(novelty))) + 
               geom_bar() + 
               xlab(xlabel) + ylab(ylabel) +
               scale_fill_manual("Isoform Type", values = colors) +
               theme_bw(base_family = "Helvetica", base_size = 18) +
               theme(axis.line.x = element_line(color="black", size = 0.5),
                     axis.line.y = element_line(color="black", size = 0.5),
                     axis.text.x = element_text(color="black", size = rel(1.5)),
                     axis.text.y = element_text(color="black", size = rel(2)),
                     axis.title.x = element_text(color="black", size=rel(1.75)),
                     axis.title.y = element_text(color="black", size=rel(1.75))) +
               guides(fill = FALSE) + 
               coord_cartesian(ylim = c(0, ymax)) 
               #theme(legend.position=c(0.8,0.7),
               #      legend.title = element_blank(),
               #      legend.background = element_rect(fill="white", color = "black"),
               #      legend.key = element_rect(fill="transparent"),
               #      legend.text = element_text(colour = 'black', size = rel(2)))

    print(g)
    dev.off()

}

plot_novelty_on_reads <- function(observed_transcripts, outdir, datasets){
    # This function plots the number of reads per dataset that got assigned to
    # each novelty type.

    # Remap levels to easier names
    observed_transcripts$novelty <- revalue(as.factor(observed_transcripts$novelty), 
                                          c("FSM_transcript"="Known",
                                            "ISM-prefix_transcript"="ISM",
                                            "ISM-suffix_transcript"="ISM",
                                            "other_ISM_transcript"="ISM",
                                            "NIC_transcript"= "NIC",
                                            "NNC_transcript"= "NNC",
                                            "antisense_transcript"= "Antisense",
                                            "intergenic_transcript"= "Intergenic"))
    observed_transcripts$novelty <- factor(observed_transcripts$novelty, levels = c("Known", "ISM",
                                                                                    "NIC", "NNC",
                                                                                    "Antisense", "Intergenic"))

    # Compute percentages
    freqs_by_dataset <- count(observed_transcripts, c("dataset","novelty"))
    freqs_by_dataset <- freqs_by_dataset %>% group_by(dataset) %>% 
                        mutate(percent = round(100*freq/sum(freq),1))
    

    observed_transcripts <- merge(observed_transcripts, freqs_by_dataset, 
                                  by = c("dataset","novelty"), all.x = T, all.y = F)
       
    # Plotting
    str_datasets <- paste(datasets, collapse='-')
    fname <- paste(outdir, "/", str_datasets, "_reads_by_isoform_category.png", sep="")
    # fname <- paste(outdir, "/reads_by_isoform_category.png", sep="")
    xlabel <- "Dataset"
    ylabel <- "log2(read count)"
    colors <- c("Known" = "#009E73","ISM" = "#0072B2", "NIC" = "#D55E00",
                "NNC" = "#E69F00", "Antisense" = "#000000",
                "Intergenic" = "#CC79A7")
    ymax <- 0.9*nrow(observed_transcripts)

    n_datasets <- length(unique(observed_transcripts$dataset))
    png(filename = fname,
        width = 2000*n_datasets + 500, height = 2500, units = "px",
        bg = "white",  res = 300)
    g = ggplot(observed_transcripts, aes(x = dataset,
               fill = as.factor(novelty))) + #factor(ERCC, levels = novelty)) +
               geom_bar(position="dodge") + #custom_theme() +
               xlab(xlabel) + ylab(ylabel) +
               theme_bw(base_family = "Helvetica", base_size = 18) +
               scale_fill_manual("Isoform Type", values = colors) +
               theme_bw(base_family = "Helvetica", base_size = 18) +
               theme(axis.line.x = element_line(color="black", size = 0.5),
                     axis.line.y = element_line(color="black", size = 0.5),
                     axis.text.x = element_text(color="black", size = rel(1.5)),
                     axis.text.y = element_text(color="black", size = rel(2)),
                     axis.title.x = element_text(color="black", size=rel(1.5)),
                     axis.title.y = element_text(color="black", size=rel(1.5))) +
               theme(legend.text = element_text(color="black", size = rel(1)),
                     legend.title = element_text(color="black", size=rel(1.25)),
                     legend.position=c(0.92,0.85),
                     legend.background = element_rect(fill="white", color = "black"),
                     legend.key = element_rect(fill="transparent")) +
                yscale("log2", .format = TRUE) +
                coord_cartesian(ylim = c(1, ymax)) +
                geom_text(aes(y = ..count.., 
                  label = paste0(percent, '%')), 
                  stat = 'count', 
                  position = position_dodge(.9), 
                  size = rel(7.5), vjust=-0.25) +
                expand_limits(x = 2)
                #guides(fill = FALSE)


    print(g)
    dev.off()
    quit()

    # Write a log file
    FSM_rows = subset(observed_transcripts, novelty == "KNOWN")
    ISM_rows = subset(observed_transcripts, novelty == "ISM")
    NIC_rows = subset(observed_transcripts, novelty == "NIC")
    NNC_rows = subset(observed_transcripts, novelty == "NNC")
    antisense_rows = subset(observed_transcripts, novelty == "antisense")
    intergenic_rows = subset(observed_transcripts, novelty == "intergenic")

    FSMs_per_dataset <- table(FSM_rows$dataset)
    ISMs_per_dataset <- table(ISM_rows$dataset)
    NIC_per_dataset <- table(NIC_rows$dataset)
    NNC_per_dataset <- table(NNC_rows$dataset)
    antisense_per_dataset <- table(antisense_rows$dataset)
    intergenic_per_dataset <- table(intergenic_rows$dataset)

    reads_per_dataset <- table(observed_transcripts$dataset)
    datasets <- names(FSMs_per_dataset)

    FSM_table <- data.frame("dataset" = as.character(datasets),
                            "total_reads" = as.numeric(reads_per_dataset),
                            "percent_FSM" = 100*as.numeric(FSMs_per_dataset)/as.numeric(reads_per_dataset),
                            "percent_ISM" = 100*as.numeric(ISMs_per_dataset)/as.numeric(reads_per_dataset),
                            "percent_NIC" = 100*as.numeric(NIC_per_dataset)/as.numeric(reads_per_dataset),
                            "percent_NNC" = 100*as.numeric(NNC_per_dataset)/as.numeric(reads_per_dataset),
                            "percent_antisense" = 100*as.numeric(antisense_per_dataset)/as.numeric(reads_per_dataset),
                            "percent_intergenic" = 100*as.numeric(intergenic_per_dataset)/as.numeric(reads_per_dataset),
                            stringsAsFactors = FALSE)
    print(FSM_table)
     
} 

load_packages <- function() {
    suppressPackageStartupMessages(library("DBI"))
    suppressPackageStartupMessages(library("RSQLite"))
    suppressPackageStartupMessages(library("readr"))
    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("plyr"))
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("ggpubr"))

    return
}

parse_options <- function() {

    option_list <- list(
        make_option(c("--db"), action = "store", dest = "database",
                    default = NULL, help = "TALON database"),
        make_option(c("--w"), action = "store", dest = "whitelist",
                    default = NULL, help = "whitelist csv file"),
        make_option(c("--datasets"), action = "store", dest = "datasets",
                    default = NULL, help = "Optional: Comma-separated list of datasets to include"),
        make_option(c("-o","--outdir"), action = "store", dest = "outdir",
                    default = NULL, help = "Output directory for plots and outfiles")
        )

    opt <- parse_args(OptionParser(option_list=option_list))
    return(opt)
}

main()
