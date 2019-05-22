# Plot a histogram of MAPQ scores

main <-function() {

    load_packages()
    opt <- parse_options()

    scores <- as.data.frame(read_delim(opt$infile, delim = ",",
                                  col_names = FALSE, trim_ws = TRUE, na = "NA"))

    colnames(scores) <- "MAPQ"

    png(filename = opt$outfile,
        width = 2500, height = 2000, units = "px",
        bg = "white",  res = 300)

    g = ggplot(scores, aes(MAPQ)) + geom_histogram(binwidth=1) + theme_bw()
    print(g)

    dev.off()
}

load_packages <- function() {
    suppressPackageStartupMessages(library("tidyverse"))
    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("optparse"))
    return
}

parse_options <- function() {

    option_list <- list(
        make_option(c("--f"), action = "store", dest = "infile",
                    default = NULL, help = "File containing MAPQ values on each line"),
        make_option(c("-o","--outfile"), action = "store", dest = "outfile",
                    default = NULL, help = "Output filename (should end in '.png'")
        )

    opt <- parse_args(OptionParser(option_list=option_list))
    return(opt)
}

main()
