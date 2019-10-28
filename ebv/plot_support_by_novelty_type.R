# Plot the percentage of support across different novelty categories
    
main <-function() {

    load_packages()
    opt <- parse_options()
    ymax = opt$ymax

    # Read RNA-PET summary file
    support_data <- as.data.frame(read_delim(opt$infile, delim = ",",
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))

    # Rename column for plotting convernience
    colnames(support_data)[2] <- "support"
    
    # Read in novelty metadata file
    novelty <- as.data.frame(read_delim(opt$novelty, delim = "\t",
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))

    # Merge RNA-PET support with the novelty labels
    support_data <- merge(support_data, novelty, by = "transcript_ID", all.x = T, all.y = F)

    # Reformat the novelty information
    support_data <- support_data %>% mutate(novelty = "")
    support_data[support_data$known == 1, "novelty"] <- "Known"
    if (opt$splitISM == F) {
        support_data[support_data$ISM == 1, "novelty"] <- "ISM"
    } else {
        support_data[support_data$ISM == 1, "novelty"] <- "other ISM"
        support_data[support_data$prefix_ISM == 1, "novelty"] <- "prefix ISM"
        support_data[support_data$suffix_ISM == 1, "novelty"] <- "suffix ISM"
        support_data[support_data$prefix_ISM == 1 & support_data$suffix_ISM == 1, "novelty"] <- "other ISM"
    }
    support_data[support_data$NIC == 1, "novelty"] <- "NIC"
    support_data[support_data$NNC == 1, "novelty"] <- "NNC"
    support_data[support_data$antisense == 1, "novelty"] <- "Antisense"
    support_data[support_data$genomic == 1, "novelty"] <- "Genomic"
    support_data[support_data$intergenic == 1, "novelty"] <- "Intergenic"

    if (opt$splitISM == F) {
        support_data$novelty <- factor(support_data$novelty, 
                                  levels = rev(c("Known", "ISM", "NIC", "NNC", 
                                                 "Antisense", "Intergenic", "Genomic")))
    } else {
        support_data$novelty <- factor(support_data$novelty,
                                  levels = rev(c("Known", "prefix ISM", "suffix ISM", 
                                                 "other ISM", "NIC", "NNC",
                                                 "Antisense", "Intergenic", "Genomic")))
    }

    # Plot support
    plot_support(support_data[,c("support", "novelty")], opt$data_type, 
                 color_vec, ymax, opt$outprefix)

}


plot_support <- function(data, data_type, color, ymax, outprefix) {

    print(table(data$novelty))

    # Compute percentages
    data$support <- as.factor(data$support)
    freqs <- data %>% count(support, novelty) %>%
             group_by(novelty) %>%
             mutate(freq = n / sum(n), total = sum(n))
    freqs$novelty <- factor(freqs$novelty,
                            levels = levels(data$novelty))

    freqs$percent <- round(freqs$freq*100)
    freqs[freqs$support == "no", "percent"] <- NA
    freqs$tcolor_grp <- as.factor(ifelse(freqs$percent > 20, "white", "black"))

    colors <- c("Known" = "#009E73","ISM" = "#0072B2", "NIC" = "#D55E00", 
                "NNC" = "#E69F00", "Antisense" = "#000000", "Intergenic" = "#CC79A7",
                "prefix ISM" = "#56B4E9", "suffix ISM" = "#698bac", "other ISM" = "#003366",
                "Genomic" = "#F0E442", "Novel" = "lightblue")
    fname <- paste(outprefix, "_", data_type, "_support.png", sep="")
    xlabel <- "Transcript category"
    ylabel <- "Number of transcript models"

    png(filename = fname,
        width = 2000, height = 1500, units = "px",
        bg = "white",  res = 300)

    g = ggplot(freqs, aes(x = novelty, y = n, fill = novelty,
                          alpha = support)) +
               geom_bar(stat="identity", color = "black") +
               xlab(xlabel) + ylab(ylabel) +
               theme(legend.text = element_text(color="black", size = rel(1)),
                     legend.title = element_text(color="black", size=rel(1))) +
               theme_bw(base_family = "Helvetica", base_size = 18) +
               scale_fill_manual("", values = colors) +
               scale_alpha_manual(values=c(0,1), name = "CAGE support") +
               theme_bw(base_family = "Helvetica", base_size = 18) +
               theme(axis.line.x = element_line(color="black", size = 0.5),
                     axis.line.y = element_line(color="black", size = 0.5),
                     axis.text.x = element_text(color="black", size = rel(1.5)),
                     axis.text.y = element_text(color="black", size = rel(1.5)),
                     axis.title.x = element_text(color="black", size=rel(1.25)),
                     axis.title.y = element_blank()) +
                coord_flip(ylim=c(0,ymax)) + guides(fill=FALSE, alpha = FALSE) +
                geom_text(aes(y = ifelse(percent > 20, total + 2000, total + 2000),
                          label = paste0(percent, "%"), color = novelty),
                          position = position_dodge(0.2), size = 8) +
                scale_color_manual(values = colors) +
                guides(colour=FALSE, fill=FALSE) +
                theme(legend.position=c(0.8,0.2),
                     legend.title = element_blank(),
                     legend.background = element_rect(fill="white", color = "black"),
                     legend.key = element_rect(fill="transparent"),
                     legend.text = element_text(colour = 'black', size = 16)) +
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank())

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
                    default = NULL, help = "RNA-PET support summary file"),
        make_option(c("--t"), action = "store", dest = "data_type",
                    default = NULL, help = "Supporting data type (will be used to label plots)"),
        make_option(c("--novelty"), action = "store", dest = "novelty",
                    default = NULL, help = "File mapping transcript IDs to novelty types"),
        make_option("--ymax", action = "store", dest = "ymax",
                    default = 22000, help = "Maximum value to display on y-axis"),
        make_option(c("--splitISM"), action="store_true", dest="splitISM",
              help="Set this option to plot prefix and suffix ISMs separately", default = F),
        make_option(c("-o","--outdir"), action = "store", dest = "outprefix",
                    default = NULL, help = "Output prefix")
        )

    opt <- parse_args(OptionParser(option_list=option_list))
    opt$ymax <- as.numeric(opt$ymax)
    return(opt)
}

main()
