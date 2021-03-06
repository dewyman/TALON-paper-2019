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

    # Merge RNA-PET support with the novelty labels
    rna_pet <- merge(rna_pet, novelty, by = "transcript_ID", all.x = T, all.y = F)

    # Reformat the novelty information
    rna_pet <- rna_pet %>% mutate(novelty = "")
    rna_pet[rna_pet$known == 1, "novelty"] <- "Known"
    if (opt$splitISM == F) {
        rna_pet[rna_pet$ISM == 1, "novelty"] <- "ISM"
    } else {
        rna_pet[rna_pet$ISM == 1, "novelty"] <- "other ISM"
        rna_pet[rna_pet$prefix_ISM == 1, "novelty"] <- "prefix ISM"
        rna_pet[rna_pet$suffix_ISM == 1, "novelty"] <- "suffix ISM"
        rna_pet[rna_pet$prefix_ISM == 1 & rna_pet$suffix_ISM == 1, "novelty"] <- "other ISM"
    }
    rna_pet[rna_pet$NIC == 1, "novelty"] <- "NIC"
    rna_pet[rna_pet$NNC == 1, "novelty"] <- "NNC"
    rna_pet[rna_pet$antisense == 1, "novelty"] <- "Antisense"
    rna_pet[rna_pet$intergenic == 1, "novelty"] <- "Intergenic"

    if (opt$splitISM == F) {
        rna_pet$novelty <- factor(rna_pet$novelty, 
                                  levels = rev(c("Known", "ISM", "NIC", "NNC", 
                                                 "Antisense", "Intergenic")))
    } else {
        rna_pet$novelty <- factor(rna_pet$novelty,
                                  levels = rev(c("Known", "prefix ISM", "suffix ISM", 
                                                 "other ISM", "NIC", "NNC",
                                                 "Antisense", "Intergenic")))
    }

    # Compute gene TPMs
    abundances <- compute_gene_TPMs(abundance_table, d1, d2)

    # Now merge in the TPMs
    rna_pet <- merge(rna_pet, abundances, by = "transcript_ID", all.x = T, all.y = F)
    rna_pet$gene_TPM_max <- pmax(rna_pet$gene_TPM.1, rna_pet$gene_TPM.2)

    # Compute median TPM of gene by novelty category
    median_TPM_by_novelty <- aggregate(rna_pet$gene_TPM_max, 
                                     by=list(rna_pet$novelty), FUN=median)
    colnames(median_TPM_by_novelty)[2] <- "Median_gene_TPM"

    # Compute mean TPM of gene by novelty category
    mean_TPM_by_novelty <- aggregate(rna_pet$gene_TPM_max,
                                     by=list(rna_pet$novelty), FUN=mean)
    colnames(mean_TPM_by_novelty)[2] <- "Mean_gene_TPM"

    # Plot all transcripts
    plot_support(rna_pet[,c("support", "novelty", "gene_TPM_max")], color_vec, 0, opt$outprefix)

    # Make the same support plot, but apply a TPM cutoff
    plot_support(rna_pet[,c("support", "novelty", "gene_TPM_max")], color_vec, 50, opt$outprefix)
    plot_support(rna_pet[,c("support", "novelty", "gene_TPM_max")], color_vec, 100, opt$outprefix)
    plot_support(rna_pet[,c("support", "novelty", "gene_TPM_max")], color_vec, 500, opt$outprefix) 


    # Plot gene expression level by support
    plot_expression_levels_by_support(rna_pet[,c("support", "novelty", "gene_TPM_max")], opt$outprefix)
}

plot_expression_levels_by_support <- function(data, outprefix) {
    # Create a violin plot showing the gene expression level of transcripts
    # with and without support

    fname <- paste(outprefix, "_RNA-PET_support_by_gene_expression_level.png", sep="")
    xlabel <- "RNA-PET support for transcript"
    ylabel <- "Log2 gene expression level"

    color_vec <- c("Known" = "#009E73","ISM" = "#0072B2", "NIC" = "#D55E00",
                "NNC" = "#E69F00", "Antisense" = "#000000", "Intergenic" = "#CC79A7",
                "prefix ISM" = "#56B4E9", "suffix ISM" = "#698bac", "other ISM" = "#003366")

    png(filename = fname,
        width = 2500, height = 2000, units = "px",
        bg = "white",  res = 300)

    # Get summary stats for labels
    #print(head(data))
    #quit()
    Summary.data <- data %>%  group_by(novelty, support) %>% 
                    summarise(n=n(), max = max(log2(gene_TPM_max + 1))) %>%
                    ungroup()
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
    abundance_table <- subset(abundance_table, transcript_novelty != "Genomic")

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
             mutate(freq = n / sum(n), total = sum(n))
    freqs$novelty <- factor(freqs$novelty,
                            levels = levels(data$novelty))

    freqs$percent <- round(freqs$freq*100)
    freqs[freqs$support == "no", "percent"] <- NA
    freqs$tcolor_grp <- as.factor(ifelse(freqs$percent > 20, "white", "black"))

    colors <- c("Known" = "#009E73","ISM" = "#0072B2", "NIC" = "#D55E00", 
                "NNC" = "#E69F00", "Antisense" = "#000000", "Intergenic" = "#CC79A7",
                "prefix ISM" = "#56B4E9", "suffix ISM" = "#698bac", "other ISM" = "#003366")
    fname <- paste(outprefix, "_RNA-PET_support_minTPM-", min_TPM, ".png", sep="")
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
                coord_flip(ylim=c(0,22000)) + guides(fill=FALSE, alpha = FALSE) +
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
        make_option(c("--splitISM"), action="store_true", dest="splitISM",
              help="Set this option to plot prefix and suffix ISMs separately", default = F),
        make_option(c("-o","--outdir"), action = "store", dest = "outprefix",
                    default = NULL, help = "Output prefix")
        )

    opt <- parse_args(OptionParser(option_list=option_list))
    opt$d1 <- as.character(opt$d1)
    opt$d2 <- as.character(opt$d2)
    return(opt)
}

main()
