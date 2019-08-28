main <-function() {

    set.seed(100)
    load_packages()
    opt <- parse_options()

    # Read the dataset names and groups
    groups <- read_delim(opt$groups, ",", escape_double = FALSE,
                             col_names = FALSE, trim_ws = TRUE, na = "NA")
    colnames(groups) <- c("group", "dataset")
    datasets <- groups$dataset


    # Read in transcripts file
    transcript_IDs <- as.data.frame(read_delim(opt$transcripts, ",", escape_double = FALSE,
                             col_names = FALSE, trim_ws = TRUE, na = "NA"))[,1]

    # Read in abundance matrix, then limit it to only the desired datasets
    abundance <- as.data.frame(read_delim(opt$infile, "\t", escape_double = FALSE,
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))
    abundance <- abundance[, c("annot_transcript_id", datasets)]

    # Reshape the data so that the datasets are rows
    abundance <- melt(abundance)
    colnames(abundance) <- c("annot_transcript_id", "dataset", "count")

    # Convert counts to TPMs for all datasets individually
    abundance <- abundance %>% group_by(dataset) %>% 
                               mutate(TPM = (count)*1000000/sum(count))

    # Now subset on the transcripts of interest
    abundance <- subset(abundance, annot_transcript_id %in% transcript_IDs)

    # Merge in the dataset groups, then average TPMs by group
    abundance <- merge(abundance, groups, by = "dataset", all.x = T, all.y = F)
    abundance <- abundance %>% group_by(annot_transcript_id, group) %>%
                               mutate(mean_TPM=mean(TPM))
    abundance <- unique(abundance[,c("annot_transcript_id", "group", "mean_TPM")])

    # Now generate the plot
    abundance$annot_transcript_id <- factor(abundance$annot_transcript_id,
                                            levels = transcript_IDs)
    abundance$group <- factor(abundance$group, 
                              levels = unique(groups$group))
    plot_TPMs_by_transcript(abundance, opt$outprefix)

}

plot_TPMs_by_transcript <- function(abundance, outprefix) {

    fname <- paste(outprefix, "/transcript_expression.png", sep="")
    xlabel <- ""
    ylabel <- "log2(mean(TPM)+1)\nof replicates"
    n_transcripts <- length(unique(abundance$annot_transcript_id))
    green = "#7dff8a"
    red = "#e26980"
    blue = "#573fe6" 
    colors <- c(blue, green, red)

    png(filename = fname,
        width = 500*n_transcripts + 1000, height = 4000, units = "px",
        bg = "white",  res = 300)

    g = ggplot(abundance, aes(x = annot_transcript_id, y = log2(mean_TPM + 1),
               fill = as.factor(group), color = as.factor(group))) + 
               geom_bar(stat="identity", position = "dodge", width = 0.6) + 
               xlab(xlabel) + ylab(ylabel) +
               scale_fill_manual("Cell line", values = colors) +
               scale_color_manual("Cell line", values = colors) +
               guides(color = FALSE) +
               theme_classic(base_family = "Helvetica", base_size = 28) +
               theme(axis.title.y = element_text(color="black", angle = 90, size = rel(3), 
                                                 margin = margin(r = 10)),
                     axis.text.y = element_text(color="black", angle = 90, size = rel(3), hjust = 0.5),
                     axis.text.x = element_blank(),
                     plot.margin = margin(t = 2, l = 1, r = 0.5, b = 0.5, "cm"))
    
    print(g)
    dev.off()
}

load_packages <- function() {
    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("dplyr"))
    #suppressPackageStartupMessages(library("cowplot"))
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("readr"))
    suppressPackageStartupMessages(library("reshape"))
    suppressPackageStartupMessages(library("stringr"))
    suppressPackageStartupMessages(library("grid"))
    return
}

parse_options <- function() {

    option_list <- list(
    make_option(c("--f"), action = "store", dest = "infile",
                    default = NULL, help = "TALON abundance file (filtered)"),
    make_option(c("--groups"), action = "store", dest = "groups",
                    default = NULL, help = "Optional: File of groupings for the datasets (ie if you want to average certain datasets). Format: group_name,dataset"),
    make_option(c("--transcripts"), action = "store", dest = "transcripts",
                    default = NULL, help = "File of transcript annot IDs for plotting."),
    make_option(c("-o","--outprefix"), action = "store", dest = "outprefix",
                default = NULL, help = "Output prefix for plots and outfiles"))

    opt <- parse_args(OptionParser(option_list=option_list))
    return(opt)
}


main()
