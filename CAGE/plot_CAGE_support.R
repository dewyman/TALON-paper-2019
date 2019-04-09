# Plot the percentage of CAGE support across different novelty categories
    
main <-function() {

    load_packages()
    opt <- parse_options()

    # Read CAGE summary file
    cage <- as.data.frame(read_delim(opt$infile, delim = ",",
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))

    # Rename column for plotting convernience
    colnames(cage)[2] <- "support"
    

    # Read in novelty metadata file
    novelty <- as.data.frame(read_delim(opt$novelty, delim = "\t",
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))

    # Read in expression level information
    abundance_table <- as.data.frame(read_delim(opt$abundance, "\t", escape_double = FALSE,
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))

    # Get dataset names and check that they exist in the file
    d1 <- opt$d1
    d2 <- opt$d2
    if (d1 %in% colnames(abundance_table) == F |
        d2 %in% colnames(abundance_table) == F) {
        print("One of the provided dataset names is not in the abundance file provided. Exiting...")
        quit()
    }

    # Read in antisense gene ID mapping file so that we can get the sense gene IDs
    antisense_2_sense <- as.data.frame(read_delim(opt$antisense, ",", escape_double = FALSE,
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))
    colnames(antisense_2_sense)[1] <- "gene_ID"

    # Replace antisense gene IDs in the abundance table with their sense counterparts
    abundance_table <- merge(abundance_table, antisense_2_sense, by = "gene_ID",
                             all.x = T, all.y = F)
    abundance_table[!is.na(abundance_table$sense_talon_ID), "gene_ID"] <- abundance_table[!is.na(abundance_table$sense_talon_ID), "sense_talon_ID"]

    # Merge CAGE support with the novelty labels
    cage <- merge(cage, novelty, by = "transcript_ID", all.x = T, all.y = F)

    # Reformat the novelty information
    cage <- cage %>% mutate(novelty = "")
    cage[cage$known == 1, "novelty"] <- "Known"
    cage[cage$ISM == 1, "novelty"] <- "ISM"
    cage[cage$NIC == 1, "novelty"] <- "NIC"
    cage[cage$NNC == 1, "novelty"] <- "NNC"
    cage[cage$antisense == 1, "novelty"] <- "Antisense"
    cage[cage$intergenic == 1, "novelty"] <- "Intergenic"
    cage$novelty <- factor(cage$novelty, 
                              levels = rev(c("Known", "ISM", "NIC", "NNC", 
                                         "Antisense", "Intergenic")))

    # Compute gene TPMs
    abundances <- compute_gene_TPMs(abundance_table, d1, d2)

    # Now merge in the TPMs
    cage <- merge(cage, abundances, by = "transcript_ID", all.x = T, all.y = F)
    cage$gene_TPM_max <- pmax(cage$gene_TPM.1, cage$gene_TPM.2)

    # Compute median TPM of gene by novelty category
    median_TPM_by_novelty <- aggregate(cage$gene_TPM_max, 
                                     by=list(cage$novelty), FUN=median)
    colnames(median_TPM_by_novelty)[2] <- "Median_gene_TPM"
    print(median_TPM_by_novelty) 

    # Compute mean TPM of gene by novelty category
    mean_TPM_by_novelty <- aggregate(cage$gene_TPM_max,
                                     by=list(cage$novelty), FUN=mean)
    colnames(mean_TPM_by_novelty)[2] <- "Mean_gene_TPM"
    print(mean_TPM_by_novelty)

    # Plot all transcripts
    plot_support(cage[,c("support", "novelty", "gene_TPM_max")], color_vec, 0, opt$outprefix)

    #print(cage[cage$gene_TPM_max > 30000,])

    # Make the same support plot, but apply a TPM cutoff
    plot_support(cage[,c("support", "novelty", "gene_TPM_max")], color_vec, 50, opt$outprefix)
    plot_support(cage[,c("support", "novelty", "gene_TPM_max")], color_vec, 100, opt$outprefix)
    plot_support(cage[,c("support", "novelty", "gene_TPM_max")], color_vec, 500, opt$outprefix) 


    # Plot gene expression level by support
    plot_expression_levels_by_support(cage[,c("support", "novelty", "gene_TPM_max")], opt$outprefix)
}

plot_expression_levels_by_support <- function(data, outprefix) {
    # Create a violin plot showing the gene expression level of transcripts
    # with and without support

    fname <- paste(outprefix, "_CAGE_support_by_gene_expression_level.png", sep="")
    xlabel <- "CAGE support for transcript"
    ylabel <- "Log2 gene expression level"
    data$novelty <- factor(data$novelty,
                           levels = c("Known", "ISM", "NIC", "NNC",
                                         "Antisense", "Intergenic"))

    # Set colors: novelty color for supported transcripts, and grey for unsupported
    data$subcat <- paste(data$novelty, data$support, sep='_')
    data$subcat <- factor(data$subcat, levels = c("Known_yes", "Known_no", 
                                                "ISM_yes", "ISM_no",
                                                "NIC_yes", "NIC_no", 
                                                "NNC_yes", "NNC_no",
                                                "Antisense_yes", "Antisense_no", 
                                                "Intergenic_yes", "Intergenic_no"))
    values = c("A" = "#E08214", "B" = "#E08214")
    color_vec <- c("Known" = "#009E73","ISM" = "#0072B2", "NIC" = "#D55E00",
                "NNC" = "#E69F00", "Antisense" = "#1A1A1A", "Intergenic" = "#CC79A7")

    png(filename = fname,
        width = 2500, height = 2000, units = "px",
        bg = "white",  res = 300)

    # Get summary stats for labels
    data %>%  group_by(novelty, support) %>% summarise(n=n(), max = max(log2(gene_TPM_max + 1))) ->Summary.data
    print(Summary.data)

    g = ggplot(data, aes(x = factor(support, levels = c("yes", "no")),
                         y = log2(gene_TPM_max + 1), fill = novelty)) +
               geom_violin(alpha = 0.8) + 
               geom_boxplot(width=0.1, fill="white", outlier.size=-1) +
               coord_cartesian(ylim = c(0, 17)) +
               xlab(xlabel) + ylab(ylabel) +
               theme_bw(base_family = "Helvetica", base_size = 18) +
               scale_fill_manual("", values = color_vec) +
               theme_bw(base_family = "Helvetica", base_size = 18) +
               theme(axis.line.x = element_line(color="black", size = 0.5),
                     axis.line.y = element_line(color="black", size = 0.5),
                     axis.text.x = element_text(color="black", size = rel(1.5)),
                     axis.text.y = element_text(color="black", size = rel(1.5)),
                     axis.title.x = element_text(color="black", size=rel(1.25)),
                     axis.title.y = element_text(color="black", size=rel(1.25))) +
                guides(fill=FALSE) + facet_wrap(~ novelty) +
                geom_text(data=Summary.data ,aes(x = support, y = max + 1, 
                          label = paste0("n=", n)), color="black", size = 5) 
    print(g)    

    dev.off()
}


compute_gene_TPMs <- function(abundance_table, d1, d2) {

    # Remove genomic transcripts
    abundance_table <- subset(abundance_table, genomic_transcript == "No")

    # Restrict to genes that were observed in at least one of the datasets
    abundance_table <- abundance_table[abundance_table[,d1] +
                                       abundance_table[,d2] > 0, ]

    # Compute gene TPMs for each dataset
    d1_abundance <- abundance_table[,c("gene_ID", d1)]
    d2_abundance <- abundance_table[,c("gene_ID", d2)]
    colnames(d1_abundance) <- c("gene_ID", "count")
    colnames(d2_abundance) <- c("gene_ID", "count")

    d1_abundance$TPM <- (d1_abundance$count)*1000000/sum(d1_abundance$count)
    d2_abundance$TPM <- (d2_abundance$count)*1000000/sum(d2_abundance$count)

    d1_abundance_agg <- aggregate(d1_abundance$TPM, by=list(d1_abundance[, "gene_ID"]), FUN=sum)
    colnames(d1_abundance_agg) <- c( "gene_ID", "TPM")
    d2_abundance_agg <- aggregate(d2_abundance$TPM, by=list(d2_abundance[, "gene_ID"]), FUN=sum)
    colnames(d2_abundance_agg) <- c("gene_ID", "TPM")      
   
    gene_abundances <- merge(d1_abundance_agg, d2_abundance_agg, by = "gene_ID", all.x = T, all.y = T)
    colnames(gene_abundances) <- c("gene_ID", "gene_TPM.1", "gene_TPM.2")
    gene_abundances[is.na(gene_abundances)] <- 0 
    
    # Create a table with the following:
    #     Transcript name (labeled as ID
    #     Gene TPM in first dataset
    #     Gene TPM in second dataset
    transcripts_with_gene_TPMs <- merge(abundance_table, gene_abundances,
                                        by = "gene_ID", all.x = T, all.y = F)
    transcripts_with_gene_TPMs <- transcripts_with_gene_TPMs[,c("annot_transcript_id", 
                                                                "gene_TPM.1", 
                                                                "gene_TPM.2")]
    colnames(transcripts_with_gene_TPMs)[1] <- c("transcript_ID")
    return(transcripts_with_gene_TPMs)
}

plot_support <- function(data, color, min_TPM, outprefix) {

    # Filter data
    data <- subset(data, gene_TPM_max >= min_TPM)

    print(summary(subset(data, novelty == "Antisense")$gene_TPM_max))
    print(table(data$novelty))

    # Compute percentages
    data$support <- as.factor(data$support)
    freqs <- data %>% count(support, novelty) %>%
             group_by(novelty) %>%
             mutate(freq = n / sum(n)) %>% filter(support == "yes")
    freqs$novelty <- factor(freqs$novelty,
                           levels = rev(c("Known", "ISM", "NIC", "NNC",
                                         "Antisense", "Intergenic")))

    colors <- c("Known" = "#009E73","ISM" = "#0072B2", "NIC" = "#D55E00", 
                "NNC" = "#E69F00", "Antisense" = "#000000", "Intergenic" = "#CC79A7")
    fname <- paste(outprefix, "_CAGE_support_minTPM-", min_TPM, ".png", sep="")
    xlabel <- "Transcript category"
    ylabel <- "Fraction transcripts with CAGE support"

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
                coord_flip(ylim = c(0, 1)) + guides(fill=FALSE) +
                geom_text(aes(y = freq - 0.12, label = paste0(round(freq*100), '%')), 
                  position = position_dodge(0.9),
                  color = "white", size = 10)
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
                    default = NULL, help = "CAGE support summary file"),
        make_option(c("--novelty"), action = "store", dest = "novelty",
                    default = NULL, help = "File mapping transcript IDs to novelty types"),
        make_option(c("--abundance"), action = "store", dest = "abundance",
                    default = NULL, help = "Abundance file (unfiltered because it is used for gene TPMs)"),
        make_option("--d1", action = "store", dest = "d1",
                    default = NULL, help = "First dataset name to use in comparison"),
        make_option("--d2", action = "store", dest = "d2",
                    default = NULL, help = "Second dataset name to use in comparison"),
        make_option("--as", action = "store", dest = "antisense",
                    default = NULL, help = "File mapping antisense TALON IDs to the sense TALON gene IDs"),
        make_option(c("-o","--outdir"), action = "store", dest = "outprefix",
                    default = NULL, help = "Output prefix")
        )

    opt <- parse_args(OptionParser(option_list=option_list))
    opt$d1 <- as.character(opt$d1)
    opt$d2 <- as.character(opt$d2)
    return(opt)
}

main()
