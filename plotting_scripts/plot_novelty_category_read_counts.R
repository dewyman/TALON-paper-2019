main <-function() {

    load_packages()
    opt <- parse_options()
    datasets <- opt$datasets
    outdir <- opt$outdir

    # Read abundance table
    abundance_table <- as.data.frame(read_delim(opt$infile, delim = "\t",
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))

    # Filter the results to limit them to selected datasets
    datasets <- unlist(strsplit(datasets, ","))
    abundance_table <- abundance_table[,c("transcript_ID", "transcript_novelty", datasets)]

    # Remove genomic transcripts
    abundance_table <- subset(abundance_table, transcript_novelty != "Genomic")

    # Remove transcripts if they do not appear in any of the datasets
    if(length(datasets) > 1) {
        abundance_table$total <- rowSums(abundance_table[,datasets])
    } else {
        abundance_table$total <- abundance_table[,datasets]
    }
    abundance_table <- subset(abundance_table, total > 0)
    abundance_table$total <- NULL
    # Plot
    plot_novelty_on_reads(abundance_table, outdir, datasets)
}

plot_novelty_on_reads <- function(observed_transcripts, outdir, datasets){
    # This function plots the number of reads per dataset that got assigned to
    # each novelty type.

    observed_transcripts$novelty <- factor(observed_transcripts$transcript_novelty,
                                             levels = c("Known", "ISM", "NIC", "NNC",
                                             "Antisense", "Intergenic"))
    observed_transcripts <- gather(observed_transcripts, dataset, reads, datasets)
    observed_transcripts <- suppressMessages(expandRows(observed_transcripts, 
                                             "reads", count.is.col = TRUE, 
                                             drop = TRUE)[,c("novelty", "dataset")])

    # Compute percentages
    freqs_by_dataset <- observed_transcripts %>% group_by(dataset,novelty) %>% count()
    print(freqs_by_dataset)
    freqs_by_dataset <- freqs_by_dataset %>% group_by(dataset) %>% 
                        mutate(percent = round(100*n/sum(n),1))
    print(freqs_by_dataset)
    observed_transcripts <- merge(observed_transcripts, freqs_by_dataset, 
                                  by = c("dataset","novelty"), all.x = T, all.y = F)
      
    # Plotting
    str_datasets <- paste(datasets, collapse='-')
    fname <- paste(outdir, "/", str_datasets, "_reads_by_isoform_category.png", sep="")
    xlabel <- "Dataset"
    ylabel <- "Read count"
    colors <- c("Known" = "#009E73","ISM" = "#0072B2", "NIC" = "#D55E00",
                "NNC" = "#E69F00", "Antisense" = "#000000",
                "Intergenic" = "#CC79A7")
    ymax <- 0.9*nrow(observed_transcripts)

    n_datasets <- length(unique(observed_transcripts$dataset))
    png(filename = fname,
        width = 2000*n_datasets + 500, height = 2500, units = "px",
        bg = "white",  res = 300)
    g = ggplot(observed_transcripts, aes(x = dataset,
               fill = as.factor(novelty))) + 
               geom_bar(position="dodge") + 
               xlab("") + ylab(ylabel) +
               theme_bw(base_family = "Helvetica", base_size = 18) +
               scale_fill_manual("Isoform Type", values = colors) +
               theme_bw(base_family = "Helvetica", base_size = 18) +
               theme(axis.line.x = element_line(color="black", size = 0.5),
                     axis.line.y = element_line(color="black", size = 0.5),
                     axis.text.x = element_blank(), #element_text(color="black", size = rel(1.5)),
                     axis.ticks.x = element_blank(),
                     axis.text.y = element_text(color="black", size = rel(2)),
                     axis.title.x = element_text(color="black", size=rel(1.5)),
                     axis.title.y = element_text(color="black", size=rel(1.5))) +
               theme(legend.text = element_text(color="black", size = rel(1)),
                     legend.title = element_text(color="black", size=rel(1.25)),
                     legend.position=c(0.85,0.85),
                     legend.background = element_rect(fill="white", color = "black"),
                     legend.key = element_rect(fill="transparent")) +
                yscale("log10", .format = TRUE) +
                coord_cartesian(ylim = c(1, ymax)) +
                geom_text(aes(y = ..count.., 
                  label = paste0(percent, '%')), 
                  stat = 'count', 
                  position = position_dodge(.9), 
                  size = rel(7.5), vjust=-0.25) +
                guides(fill = FALSE)


    print(g)
    dev.off()
     
} 

load_packages <- function() {
    suppressPackageStartupMessages(library("tidyverse"))
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("ggpubr"))
    suppressPackageStartupMessages(library("splitstackshape"))

    return
}

parse_options <- function() {

    option_list <- list(
        make_option(c("--f"), action = "store", dest = "infile",
                    default = NULL, help = "Filtered TALON abundance file"),
        make_option(c("--datasets"), action = "store", dest = "datasets",
                    default = NULL, help = "Comma-separated list of datasets to include"),
        make_option(c("-o","--outdir"), action = "store", dest = "outdir",
                    default = NULL, help = "Output directory for plots and outfiles")
        )

    opt <- parse_args(OptionParser(option_list=option_list))
    return(opt)
}


main()
