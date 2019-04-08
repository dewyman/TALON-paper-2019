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

    # Merge RNA-PET support with the novelty labels
    rna_pet <- merge(rna_pet, novelty, by = "transcript_ID", all.x = T, all.y = F)

    # Reformat the novelty information
    rna_pet <- rna_pet %>% mutate(novelty = "")
    rna_pet[rna_pet$known == 1, "novelty"] <- "Known"
    rna_pet[rna_pet$ISM == 1, "novelty"] <- "ISM"
    rna_pet[rna_pet$NIC == 1, "novelty"] <- "NIC"
    rna_pet[rna_pet$NNC == 1, "novelty"] <- "NNC"
    rna_pet[rna_pet$antisense == 1, "novelty"] <- "Antisense"
    rna_pet[rna_pet$intergenic == 1, "novelty"] <- "Intergenic"
    rna_pet$novelty <- factor(rna_pet$novelty, 
                              levels = rev(c("Known", "ISM", "NIC", "NNC", 
                                         "Antisense", "Intergenic")))

    # Compute gene TPMs
    abundances <- compute_gene_TPMs(abundance_table, d1, d2)

    # Now merge in the TPMs
    rna_pet <- merge(rna_pet, abundances, by = "transcript_ID", all.x = T, all.y = F)
    rna_pet$gene_TPM_max <- pmax(rna_pet$gene_TPM.1, rna_pet$gene_TPM.2)

    # Compute mean TPM of gene by novelty category
    median_TPM_by_novelty <- aggregate(rna_pet$gene_TPM_max, 
                                     by=list(rna_pet$novelty), FUN=median)
    colnames(median_TPM_by_novelty)[2] <- "Median_gene_TPM"
    print(median_TPM_by_novelty) 

    # Plot all transcripts
    plot_support(rna_pet[,c("support", "novelty", "gene_TPM_max")], color_vec, 0, opt$outprefix)

    # Make the same support plot, but apply a TPM cutoff
    plot_support(rna_pet[,c("support", "novelty", "gene_TPM_max")], color_vec, 100, opt$outprefix)
   
    # Make the same support plot, but apply a TPM cutoff
    plot_support(rna_pet[,c("support", "novelty", "gene_TPM_max")], color_vec, 500, opt$outprefix) 
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

    print(table(data$novelty))

    # Compute percentages
    data$support <- as.factor(data$support)
    freqs <- data %>% count(support, novelty) %>%
             group_by(novelty) %>%
             mutate(freq = n / sum(n)) %>% filter(support == "yes")
    freqs$novelty <- factor(freqs$novelty,
                           levels = rev(c("Known", "ISM", "NIC", "NNC",
                                         "Antisense", "Intergenic")))

    print(freqs)

    colors <- rev(c("#009E73","#0072B2", "#D55E00", "#E69F00", "#000000", "#CC79A7"))
    fname <- paste(outprefix, "_RNA-PET_support_minTPM-", min_TPM, ".png", sep="")
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
                    default = NULL, help = "RNA-PET support summary file"),
        make_option(c("--novelty"), action = "store", dest = "novelty",
                    default = NULL, help = "File mapping transcript IDs to novelty types"),
        make_option(c("--abundance"), action = "store", dest = "abundance",
                    default = NULL, help = "Abundance file (unfiltered because it is used for gene TPMs)"),
        make_option("--d1", action = "store", dest = "d1",
                    default = NULL, help = "First dataset name to use in comparison"),
        make_option("--d2", action = "store", dest = "d2",
                    default = NULL, help = "Second dataset name to use in comparison"),
        make_option(c("-o","--outdir"), action = "store", dest = "outprefix",
                    default = NULL, help = "Output prefix")
        )

    opt <- parse_args(OptionParser(option_list=option_list))
    opt$d1 <- as.character(opt$d1)
    opt$d2 <- as.character(opt$d2)
    return(opt)
}

main()
