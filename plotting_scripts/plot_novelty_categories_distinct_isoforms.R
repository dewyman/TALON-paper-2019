main <-function() {

    load_packages()
    opt <- parse_options()
    database <- opt$database
    datasets <- opt$datasets
    outdir <- opt$outdir

    # Read abundance table
    abundance_table <- as.data.frame(read_delim(opt$infile, delim = "\t",
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))

    # Filter the results to limit them to selected datasets
    datasets <- unlist(strsplit(datasets, ","))
    print(datasets)
    abundance_table <- abundance_table[,c("transcript_novelty", datasets)]
    
    # Remove genomic transcripts
    abundance_table <- subset(abundance_table, transcript_novelty != "Genomic")  

    # Remove transcripts if they do not appear in any of the datasets
    if(length(datasets) > 1) {
        abundance_table$total <- rowSums(abundance_table[,datasets])
    } else {
        abundance_table$total <- abundance_table[,datasets]
    }
    abundance_table <- subset(abundance_table, total > 0) 

    # Plot
    plot_distinct_novelty(abundance_table, outdir, datasets)
}

plot_distinct_novelty <- function(distinct_transcripts, outdir, datasets){
    # This function plots the number of transcripts that belong to each 
    # transcript novelty class. So it amounts to the number of unique transcripts
    # of each type that were identified in the datasets together

    distinct_transcripts$novelty <- factor(distinct_transcripts$transcript_novelty,
                                           levels = c("Known", "ISM", "NIC", "NNC",
                                                      "Antisense", "Intergenic"))

    # Plotting
    str_datasets <- paste(datasets, collapse='-')
    fname <- paste(outdir, "/", str_datasets, "_distinct_isoforms_by_category.png", sep="")
    xlabel <- "Category"
    ylabel <- "Number of distinct transcript models"
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
               #guides(fill = FALSE) + 
               coord_cartesian(ylim = c(0, ymax)) +
               theme(legend.position=c(0.7,0.8),
                     legend.title = element_text(colour = 'black', size = rel(2)),
                     legend.background = element_rect(fill="white", color = "black"),
                     legend.key = element_rect(fill="transparent"),
                     legend.text = element_text(colour = 'black', size = rel(2)))

    print(g)
    dev.off()

}


load_packages <- function() {
    suppressPackageStartupMessages(library("readr"))
    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("plyr"))
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("ggpubr"))

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
