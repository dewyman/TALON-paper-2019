# Given a filtered abundance file, plot the exon counts in a violin plot colored by novelty type

main <-function() {

    load_packages()
    opt <- parse_options()
    outdir <- opt$outdir

    # Read abundance table
    abundance_table <- as.data.frame(read_delim(opt$infile, delim = "\t",
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))


    # Remove genomic transcripts
    abundance_table <- subset(abundance_table, transcript_novelty != "Genomic")
    #print(subset(abundance_table[abundance_table$length > 20000,]))

    # Plot
    plot_by_novelty(abundance_table, outdir)

}

plot_by_novelty <- function(data, outdir) {

    print(nrow(data))

    fname <- paste(outdir, "/transcript_exonCount_by_novelty_type.png", sep="")
    xlabel <- "Novelty type"
    ylabel <- "Number of exons in transcript model"

    png(filename = fname,
        width = 4500, height = 2500, units = "px",
        bg = "white",  res = 300)

    color_vec <- c("Known" = "#009E73","ISM" = "#0072B2", "NIC" = "#D55E00",
                "NNC" = "#E69F00", "Antisense" = "#000000", "Intergenic" = "#CC79A7",
                "prefix ISM" = "#56B4E9", "suffix ISM" = "#698bac", "other ISM" = "#003366")

    data$transcript_novelty <- factor(data$transcript_novelty, 
                                  levels = c("Known", "ISM", "NIC", "NNC", 
                                                 "Antisense", "Intergenic"))

    # Apply jitter to outliers
    data <- data %>% group_by(transcript_novelty) %>% 
    mutate(outlier.high = n_exons > quantile(n_exons, .75) + 1.50*IQR(n_exons),
         outlier.low = n_exons < quantile(n_exons, .25) - 1.50*IQR(n_exons))

    data$outlier.color <- NA
    data$outlier.color[data$outlier.high] <- "black"
    data$outlier.color[data$outlier.low] <- "black"
    

    g = ggplot(data, 
               aes(x = transcript_novelty, y = n_exons, 
                   fill = transcript_novelty)) +
            stat_boxplot(geom ='errorbar') + 
            geom_boxplot(outlier.shape = NA) +
            geom_jitter(color = data$outlier.color, width = .2, size = 0.5, alpha = 0.5) +
            xlab(xlabel) + ylab(ylabel) +
            theme_bw(base_family = "Helvetica", base_size = 18) +
            scale_fill_manual("", values = color_vec) +
            theme_bw(base_family = "Helvetica", base_size = 18) +
            theme(axis.line.x = element_line(color="black", size = 0.5),
                  axis.line.y = element_line(color="black", size = 0.5),
                  axis.text.x = element_text(color="black", size = rel(1.5),
                                             angle = 25, vjust = 1, hjust=1),
                  axis.text.y = element_text(color="black", size = rel(1.5)),
                  axis.title.x = element_text(color="black", size=rel(1.5)),
                  axis.title.y = element_text(color="black", size=rel(1.5))) +
            guides(fill=F)
    print(g)

    dev.off()
}

load_packages <- function() {
    suppressPackageStartupMessages(library("tidyverse"))
    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("gridExtra"))
    suppressPackageStartupMessages(library("cowplot"))
    return
}

parse_options <- function() {

    option_list <- list(
        make_option(c("--f"), action = "store", dest = "infile",
                    default = NULL, help = "TALON abundance output file"),
        make_option(c("-o","--outdir"), action = "store", dest = "outdir",
                    default = NULL, help = "Output directory for plots and outfiles")
        )

    opt <- parse_args(OptionParser(option_list=option_list))
    return(opt)
}

main()
