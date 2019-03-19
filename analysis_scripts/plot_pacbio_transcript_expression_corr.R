main <-function() {

    load_packages()
    opt <- parse_options()

    # Read abundance table
    abundance_table <- as.data.frame(read_delim(opt$infile, delim = "\t",
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))

    # Get dataset names and check that they exist in the file
    d1 <- opt$d1
    d2 <- opt$d2
    if (d1 %in% colnames(abundance_table) == F |
        d2 %in% colnames(abundance_table) == F) {
        print("One of the provided dataset names is not in the abundance file provided. Exiting...")
        quit()
    }

    filtered <- filter_transcripts_on_options(abundance_table, opt)
    abundance_table <- filtered$abundance_table
    color_vec <- filtered$color_vec
    t_levels <- filtered$t_levels

    # Compute TPMs for each dataset
    d1_abundance <- abundance_table[,c("transcript_ID", d1)]
    d2_abundance <- abundance_table[,c("transcript_ID", d2)]
    colnames(d1_abundance) <- c("transcript_ID", "count")
    colnames(d2_abundance) <- c("transcript_ID", "count")

    d1_abundance$data1.TPM <- (d1_abundance$count)*1000000/sum(d1_abundance$count)
    d2_abundance$data2.TPM <- (d2_abundance$count)*1000000/sum(d2_abundance$count)

    merged_abundances <- merge(d1_abundance, d2_abundance, by = "transcript_ID", all.x = T, all.y = T)
    #colnames(merged_abundances) <- c("transcript_ID", "data1.TPM", "data2.TPM")
    merged_abundances[is.na(merged_abundances)] <- 0   

    # Merge in transcript novelty information
    transcript_novelty <- abundance_table[, c("transcript_ID", 
                                              "transcript_status", 
                                              "ISM_transcript",
                                              "NIC_transcript",
                                              "NNC_transcript",
                                              "antisense_transcript", 
                                              "intergenic_transcript",
                                              "genomic_transcript")]
    merged_abundances <- merge(merged_abundances, transcript_novelty, by = "transcript_ID", 
                               all.x = T, all.y = F)
    merged_abundances$novelty <- NA
    print(nrow(merged_abundances))

    merged_abundances[merged_abundances$transcript_status == "KNOWN", "novelty"] <- "Known"
    merged_abundances[merged_abundances$ISM_transcript == "ISM_transcript", "novelty"] <- "ISM"
    merged_abundances[merged_abundances$NIC_transcript == "NIC_transcript", "novelty"] <- "NIC"
    merged_abundances[merged_abundances$NNC_transcript == "NNC_transcript", "novelty"] <- "NNC"
    merged_abundances[merged_abundances$antisense_transcript == "antisense_transcript", "novelty"] <- "Antisense"
    merged_abundances[merged_abundances$intergenic_transcript == "intergenic_transcript", "novelty"] <- "Intergenic"
    merged_abundances[merged_abundances$genomic_transcript == "genomic_transcript", "novelty"] <- "Genomic"
    merged_abundances$novelty <- factor(merged_abundances$novelty, levels = t_levels)

    # Plot expression scatterplots
    expression_by_status(merged_abundances, d1, d2, opt$outdir, color_vec)
}

filter_transcripts_on_options <- function(abundance_table, opt) {

    d1 <- opt$d1
    d2 <- opt$d2

    # Get colors
    green <- "#009E73"
    blue <- "#0072B2"
    orange <- "#D55E00"
    gold <- "#E69F00"
    black <- "#000000"
    pink <- "#CC79A7"
    yellow <- "#F0E442"

    color_vec <- c(green)
    t_levels <- c("Known")

    # Restrict to transcripts that were observed in at least one of the datasets
    abundance_table <- abundance_table[abundance_table[,d1] +
                                       abundance_table[,d2] > 0, ]

    # TODO: Use the whitelist to filter transcripts if one is provided

    # Remove transcripts depending on the input options
    if (opt$ISM == F) {
        message("Removing ISM transcripts...")
        abundance_table <- subset(abundance_table,
                                  ISM_transcript == "No")
    } else {
        color_vec <- c(color_vec, blue)
        t_levels <- c(t_levels, "ISM")
    }
    if (opt$NIC == F) {
        message("Removing NIC transcripts...")
        abundance_table <- subset(abundance_table,
                                  NIC_transcript == "No")
    } else {
        color_vec <- c(color_vec, orange)
        t_levels <- c(t_levels, "NIC")
    }
    if (opt$NNC == F) {
        message("Removing NNC transcripts...")
        abundance_table <- subset(abundance_table,
                                  NNC_transcript == "No")
    } else {
        color_vec <- c(color_vec, gold)
        t_levels <- c(t_levels, "NNC")
    }
    if (opt$antisense == F) {
        message("Removing antisense transcripts...")
        abundance_table <- subset(abundance_table,
                                  antisense_transcript == "No")
    }
    else {
        color_vec <- c(color_vec, black)
        t_levels <- c(t_levels, "Antisense")
    }
    if (opt$intergenic == F) {
        message("Removing intergenic transcripts")
        abundance_table <- subset(abundance_table,
                                  intergenic_transcript == "No")
    }
    else {
        color_vec <- c(color_vec, pink)
        t_levels <- c(t_levels, "Intergenic")
    }
    if (opt$genomic == F) {
        message("Removing genomic transcripts")
        abundance_table <- subset(abundance_table,
                                  genomic_transcript == "No")
    }
    else {
        color_vec <- c(color_vec, yellow)
        t_levels <- c(t_levels, "Genomic")
    }

   
    filtered <- list(abundance_table, color_vec, t_levels)
    names(filtered) <- c("abundance_table", "color_vec", "t_levels")
    return(filtered)
}

expression_by_status <- function(merged_abundances, d1, d2, outdir, color_vec) {

    # Take log2(TPM + 1)
    merged_abundances$data1.TPM = log(merged_abundances$data1.TPM + 1, base=2)
    merged_abundances$data2.TPM = log(merged_abundances$data2.TPM + 1, base=2)
    t_levels <- levels(merged_abundances$novelty)
    
    # Plot log2(TPM + 1) for each dataset on a scatterplot. Color points according to known/novel status
    pearsonCorr = cor.test(~data1.TPM + data2.TPM, data=merged_abundances, method = "pearson", continuity = FALSE, conf.level = 0.95)$estimate
    spearmanCorr = cor.test(~data1.TPM + data2.TPM, data=merged_abundances, method = "spearman", continuity = FALSE, conf.level = 0.95)$estimate

    nov_types <- paste(t_levels, collapse = '-')
    joined_names <- paste(outdir, "/", d1, "-", d2, "_", nov_types, sep = "")
    
    fname <- paste(joined_names, "transcript", "correlationPlot.png", sep="_")
    xlabel <- paste("log2(TPM+1) of transcript in ", "Rep1", sep="")
    ylabel <- paste("log2(TPM+1) of transcript in ", "Rep2", sep="")

    png(filename = fname,
        width = 2500, height = 2500, units = "px",
        bg = "white",  res = 300)
    g = ggplot(merged_abundances, aes(x = data1.TPM, y = data2.TPM, color = novelty)) +
        geom_jitter(alpha = 0.5) + theme_bw() +
        #geom_rug(alpha = 1/2, position = "jitter") +
        xlab(xlabel)  + ylab(ylabel) + theme(text= element_text(size=24)) +
        theme(axis.text.x = element_text(color = "black", size=24),
              axis.text.y = element_text(color = "black", size=24)) +
        annotate("text", x = 5, y = 14, label = paste("Pearson r: ",
                 round(pearsonCorr, 2), "\nSpearman rho: ",
                 round(spearmanCorr, 2), sep=""),  color="black", size = 8) +
        coord_cartesian(xlim=c(0, 16), ylim=c(0, 16)) +
        scale_colour_manual("", values=color_vec) +
        theme(legend.position=c(0.8,0.2),
              legend.background = element_rect(fill="transparent"),
              legend.key = element_rect(fill="transparent"),
              legend.text = element_text(colour = 'black', size = 20))

     g = ggMarginal(g, groupColour = TRUE, groupFill = TRUE )

     print(g)
     dev.off()

}

load_packages <- function() {
    suppressPackageStartupMessages(library("tidyverse"))
    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("ggExtra"))
    return
}

parse_options <- function() {

    option_list <- list(
        make_option(c("--f"), action = "store", dest = "infile",
                    default = NULL, help = "TALON abundance output file"),
        make_option(c("--color"), action = "store", dest = "color_scheme",
                    default = NULL, help = "blue, red, or green"),
        make_option("--d1", action = "store", dest = "d1",
                    default = NULL, help = "First dataset name to use in comparison"),
        make_option("--d2", action = "store", dest = "d2",
                    default = NULL, help = "Second dataset name to use in comparison"),
        make_option("--w", action = "store", dest = "whitelist",
                    default = NULL, help = "Optional: a whitelist for filtering transcripts"),
        make_option(c("--ISM"), action="store_true", dest="ISM",
              help="Set this option to include ISM transcripts", default = F),
        make_option(c("--NIC"), action="store_true", dest="NIC",
              help="Set this option to include NIC transcripts", default = F),
        make_option(c("--NNC"), action="store_true", dest="NNC",
              help="Set this option to include ISM transcripts", default = F),
        make_option(c("--antisense"), action="store_true", dest="antisense",
              help="Set this option to include antisense transcripts", default = F),
        make_option(c("--intergenic"), action="store_true", dest="intergenic",
              help="Set this option to include intergenic transcripts", default = F),
        make_option(c("--genomic"), action="store_true", dest="genomic",
              help="Set this option to include genomic transcripts", default = F),
        make_option(c("-o","--outdir"), action = "store", dest = "outdir",
                    default = NULL, help = "Output directory for plots and outfiles")
        )

    opt <- parse_args(OptionParser(option_list=option_list))
    opt$d1 <- as.character(opt$d1)
    opt$d2 <- as.character(opt$d2)
    return(opt)
}

main()
