# Plot the percentage of RNA-PET support across different novelty categories
    
main <-function() {

    load_packages()
    opt <- parse_options()

    # Read RNA-PET summary file
    rna_pet <- as.data.frame(read_delim(opt$infile, delim = ",",
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))

    # Rename column for plotting convernience
    colnames(rna_pet)[2] <- "support"
    

    # Read in novelty metadata file
    novelty <- as.data.frame(read_delim(opt$novelty, delim = "\t",
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))

    # Merge RNA-PET support with the novelty labels
    rna_pet <- merge(rna_pet, novelty, by = "transcript_ID", all.x = T, all.y = F)

    # Reformat the novelty information
    rna_pet <- rna_pet %>% mutate(novelty = "")
    rna_pet[rna_pet$known == 1, "novelty"] <- "Known"
    rna_pet[rna_pet$ISM == 1, "novelty"] <- "ISM"
    rna_pet[rna_pet$NIC == 1, "novelty"] <- "NIC"
    rna_pet[rna_pet$NNC == 1, "novelty"] <- "NNC"
    #rna_pet[rna_pet$genomic == 1, "novelty"] <- "Genomic"
    rna_pet[rna_pet$antisense == 1, "novelty"] <- "Antisense"
    rna_pet[rna_pet$intergenic == 1, "novelty"] <- "Intergenic"
    rna_pet$novelty <- factor(rna_pet$novelty, 
                              levels = rev(c("Known", "ISM", "NIC", "NNC", 
                                         "Antisense", "Intergenic")))
 
    plot_support(rna_pet[,c("support", "novelty")], color_vec, opt$outprefix)
}

plot_support <- function(data, color, outprefix) {

    # Compute percentages
    data$support <- as.factor(data$support)
    freqs <- data %>% count(support, novelty) %>%
             group_by(novelty) %>%
             mutate(freq = n / sum(n)) %>% filter(support == "yes")
    freqs$novelty <- factor(freqs$novelty,
                           levels = rev(c("Known", "ISM", "NIC", "NNC",
                                         "Antisense", "Intergenic")))

    colors <- c("#009E73","#0072B2", "#D55E00", "#E69F00", "#000000", "#CC79A7")
    fname <- paste(outprefix, "RNA-PET_support.png", sep="_")
    xlabel <- "Transcript category"
    ylabel <- "Fraction transcripts with RNA-PET support"

    png(filename = fname,
        width = 2500, height = 2000, units = "px",
        bg = "white",  res = 300)

    g = ggplot(freqs, aes(x = novelty, y = freq, fill = novelty)) + 
               geom_bar(stat="identity") +
               xlab(xlabel) + ylab(ylabel) +
               theme(legend.text = element_text(color="black", size = rel(1)), 
                     legend.title = element_text(color="black", size=rel(1))) +
               theme_bw(base_family = "Helvetica", base_size = 18) +
               scale_fill_manual("Isoform Type", values = colors) +
               theme_bw(base_family = "Helvetica", base_size = 18) +
               theme(axis.line.x = element_line(color="black", size = 0.5),
                     axis.line.y = element_line(color="black", size = 0.5),
                     axis.text.x = element_text(color="black", size = rel(1.5)),
                     axis.text.y = element_text(color="black", size = rel(1.5)),
                     axis.title.x = element_text(color="black", size=rel(1.25)),
                     axis.title.y = element_text(color="black", size=rel(1.25))) +
                coord_flip(ylim = c(0, 1)) + guides(fill=FALSE)
                #geom_text(aes(y = ..count.., 
                #  label = paste0(percent, '%')), 
                #  stat = 'count', 
                #  position = position_dodge(.9), 
                #  size = 7, vjust=-0.25)
    print(g)
    dev.off()
} 

load_packages <- function() {
    suppressPackageStartupMessages(library("tidyverse"))
    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("optparse"))
    #suppressPackageStartupMessages(library("plyr"))
    #suppressPackageStartupMessages(library("gridExtra"))
    #suppressPackageStartupMessages(library("cowplot"))
    return
}

parse_options <- function() {

    option_list <- list(
        make_option(c("--f"), action = "store", dest = "infile",
                    default = NULL, help = "RNA-PET support summary file"),
        make_option(c("--novelty"), action = "store", dest = "novelty",
                    default = NULL, help = "File mapping transcript IDs to novelty types"),
        make_option(c("-o","--outdir"), action = "store", dest = "outprefix",
                    default = NULL, help = "Output prefix")
        )

    opt <- parse_args(OptionParser(option_list=option_list))
    return(opt)
}

main()
