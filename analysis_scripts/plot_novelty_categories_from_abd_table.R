main <-function() {

    load_packages()
    opt <- parse_options()
    datasets <- opt$datasets
    outdir <- opt$outdir

    # Read in abundance file
    transcript_table <- as.data.frame(read_delim(opt$infile, delim = "\t",
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))

    # If datasets selected, filter the results to limit them to those datasets
    if (!is.null(datasets)) {
        datasets <- unlist(strsplit(datasets, ","))
        transcript_table <- transcript_table[,c("transcript_ID",
                                                "transcript_novelty", datasets)]
    } else {
        datasets <- names(transcript_table)[-(which(names(transcript_table) %in%
                          c("gene_ID", "transcript_ID", "annot_gene_id",
                            "annot_transcript_id", "annot_gene_name",
                            "annot_transcript_name", "n_exons",
                            "length", "gene_novelty", "transcript_novelty",
                            "ISM_subtype")))]
        transcript_table <- transcript_table[, c("transcript_ID",
                                                "transcript_novelty", datasets)] 
    }
    # Reformat the data to long
    data <- melt(transcript_table, id.vars = c("transcript_ID", "transcript_novelty"))
    colnames(data) <- c("transcript_ID", "transcript_novelty", "dataset", "count")
    data <- subset(data, count > 0)
   
    plot_distinct_novelty(data, outdir, datasets)
    plot_novelty_on_reads(data, outdir, datasets)
}

plot_distinct_novelty <- function(observed_transcripts, outdir, datasets){
    # This function plots the number of whitelist items that belong to each 
    # transcript novelty class. So it amounts to the number of unique transcripts
    # of each type that were identified in all of the datasets together

    distinct_transcripts <- unique(observed_transcripts[,c("transcript_ID", "transcript_novelty")])
    distinct_transcripts$transcript_novelty <- as.factor(distinct_transcripts$transcript_novelty)

    freqs <- distinct_transcripts %>%
                            group_by(transcript_novelty) %>% 
                            count()

    freqs$percent <- round(100*freqs$n/sum(freqs$n), 1)
    print(freqs)
 
    distinct_transcripts <- merge(distinct_transcripts, freqs, 
                                  by = c("transcript_novelty"), all.x = T, all.y = F)

    # Plotting
    str_datasets <- paste(datasets, collapse='-')
    fname <- paste(outdir, "/", str_datasets, "_distinct_isoforms_by_category.png", sep="")
    xlabel <- "Category"
    ylabel <- "Number of distinct transcript models"
    ymax <- 1.02*(max(freqs$n))
    distinct_transcripts$transcript_novelty <- factor(distinct_transcripts$transcript_novelty, 
                                               levels = c("Known", "ISM", 
                                                          "NIC", "NNC", 
                                                          "Antisense", "Intergenic"))

    colors <- c("Known" = "#009E73","ISM" = "#0072B2", "NIC" = "#D55E00", 
                "NNC" = "#E69F00", "Antisense" = "#000000", 
                "Intergenic" = "#CC79A7")

    png(filename = fname,
        width = 3000, height = 3500, units = "px",
        bg = "white",  res = 300)
    g = ggplot(distinct_transcripts, aes(x = transcript_novelty, width=.6,
               fill = as.factor(transcript_novelty))) + 
               geom_bar() + 
               ylab(ylabel) +
               scale_fill_manual("Transcript Type", values = colors) +
               theme_bw(base_family = "Helvetica", base_size = 20) +
               theme(axis.line.x = element_line(color="black", size = 0.5),
                     axis.line.y = element_line(color="black", size = 0.5),
                     axis.text.x = element_text(color="black", size = rel(1.75),
                                             angle = 25, vjust = 1, hjust=1),
                     axis.text.y = element_text(color="black", size = rel(2)),
                     axis.title.x = element_blank(),
                     axis.title.y = element_text(color="black", size=rel(1.75))) +
               geom_text(aes(y = ..count.., 
                  label = paste0(percent, '%')), 
                  stat = 'count', 
                  position = position_dodge(.9), 
                  size = rel(7.5), vjust=-0.25) +
               guides(fill=F) 

    print(g)
    dev.off()

}

plot_novelty_on_reads <- function(observed_transcripts, outdir, datasets){
    # This function plots the number of reads per dataset that got assigned to
    # each novelty type.

    observed_transcripts$transcript_novelty <- factor(observed_transcripts$transcript_novelty, 
                                                      levels = c("Known", "ISM",
                                                                 "NIC", "NNC",
                                                                 "Antisense", "Intergenic"))

    # Compute percentages
    freqs_by_dataset <- observed_transcripts %>%
                        group_by(dataset, transcript_novelty) %>%
                        mutate(freq = sum(count))
    freqs_by_dataset <- unique(freqs_by_dataset[,c("transcript_novelty", "dataset", "freq")])
    freqs_by_dataset <- freqs_by_dataset %>% group_by(dataset) %>% 
                        mutate(percent = round(100*freq/sum(freq),1))
    print(freqs_by_dataset) 

       
    # Plotting
    str_datasets <- paste(datasets, collapse='-')
    fname <- paste(outdir, "/", str_datasets, "_reads_by_isoform_category.png", sep="")
    xlabel <- "Dataset"
    ylabel <- "log2(read count)"
    colors <- c("Known" = "#009E73","ISM" = "#0072B2", "NIC" = "#D55E00",
                "NNC" = "#E69F00", "Antisense" = "#000000",
                "Intergenic" = "#CC79A7")
    ymax <- 1.3*max(freqs_by_dataset$freq)

    n_datasets <- length(unique(observed_transcripts$dataset))
    png(filename = fname,
        width = 2000*n_datasets + 500, height = 2500, units = "px",
        bg = "white",  res = 300)
    g = ggplot(freqs_by_dataset, aes(x = dataset, y = freq,
               fill = as.factor(transcript_novelty))) + 
               geom_bar(stat = "identity", position="dodge") + 
               xlab("") + ylab(ylabel) +
               theme_bw(base_family = "Helvetica", base_size = 18) +
               scale_fill_manual("Isoform Type", values = colors) +
               theme_bw(base_family = "Helvetica", base_size = 18) +
               theme(axis.line.x = element_line(color="black", size = 0.5),
                     axis.line.y = element_line(color="black", size = 0.5),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.y = element_text(color="black", size = rel(2)),
                     axis.title.x = element_text(color="black", size=rel(1.5)),
                     axis.title.y = element_text(color="black", size=rel(1.5))) +
               theme(legend.text = element_text(color="black", size = rel(1)),
                     legend.title = element_text(color="black", size=rel(1.25)),
                     legend.position=c(0.85,0.85),
                     legend.background = element_rect(fill="white", color = "black"),
                     legend.key = element_rect(fill="transparent")) +
                yscale("log2", .format = TRUE) +
                coord_cartesian(ylim = c(1, ymax)) +
                geom_text(aes(y = freq, 
                  label = paste0(percent, '%')), 
                  position = position_dodge(.9), 
                  size = rel(7.5), vjust=-0.25) +
                guides(fill = FALSE)


    print(g)
    dev.off()
     
} 

load_packages <- function() {
    suppressPackageStartupMessages(library("reshape2"))
    suppressPackageStartupMessages(library("tidyverse"))
    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("ggpubr"))

    return
}

parse_options <- function() {

    option_list <- list(
        make_option(c("--f"), action = "store", dest = "infile",
                    default = NULL, help = "TALON abundance file"),
        make_option(c("--datasets"), action = "store", dest = "datasets",
                    default = NULL, help = "Optional: Comma-separated list of datasets to include"),
        make_option(c("-o","--outdir"), action = "store", dest = "outdir",
                    default = NULL, help = "Output directory for plots and outfiles")
        )

    opt <- parse_args(OptionParser(option_list=option_list))
    return(opt)
}

main()
