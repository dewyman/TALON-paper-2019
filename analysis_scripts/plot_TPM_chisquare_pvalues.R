main <-function() {

    set.seed(100)
    load_packages()
    opt <- parse_options()

    # Get colors
    if (opt$color_scheme == "red") {
        fill_color <- "red2"
    } else if (opt$color_scheme == "blue") {
        fill_color <- "navy"
    } else if (opt$color_scheme == "green") {
        fill_color <- "springgreen4"
    }

    # Get the names of the first and second dataset that we will be working with
    data_names <- str_split(opt$datasets, ",")[[1]]
    dataset1 <- data_names[1]
    dataset2 <- data_names[2]

    # Get genes expressed in the Illumina data from the Kallisto
    # abundance file
    illumina_gene_table <- filter_kallisto_illumina_genes(opt$illumina_kallisto)   
    colnames(illumina_gene_table) <- c("annot_gene_name", "illumina_TPM")


    # Read PacBio abundance file
    pb_abundance <- as.data.frame(read_delim(opt$infile, "\t", escape_double = FALSE,
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))
    
    # Keep known genes only
    pb_abundance <- subset(pb_abundance, gene_status == "KNOWN")

    # Cut out unnecessary cols
    pb_abundance <- pb_abundance[, c("annot_gene_name", dataset1, dataset2)]

    # Aggregate PacBio by gene name to get gene counts
    pb_gene_abundance <- ddply(pb_abundance, c("annot_gene_name"), function(x) colSums(x[c(dataset1, dataset2)]))
    pb_gene_abundance$both_pacbio <- pb_gene_abundance[,dataset1] + pb_gene_abundance[,dataset2]

    # Merge PacBio with Illumina on annot_gene_name
    merged_illumina_pacbio <- merge(illumina_gene_table, pb_gene_abundance, by = "annot_gene_name",
                                    all.x = T, all.y = T)
    merged_illumina_pacbio[is.na(merged_illumina_pacbio)] <- 0

    IP <- merged_illumina_pacbio[, c("both_pacbio", "illumina_TPM")]
    IP$illumina_TPM <- round(IP$illumina_TPM)

    # Perform quantile normalization
    gene_counts <- as.matrix(IP[, c("both_pacbio", "illumina_TPM")])
    gene_counts <- as.data.frame(normalize.quantiles(gene_counts))
    final_table <- data.frame( norm_both_pacbio = round(gene_counts[,1]),
                                  norm_illumina_TPM = round(gene_counts[,2]))
    total_pacbio <- sum(final_table$norm_both_pacbio)
    total_illumina <- sum(final_table$norm_illumina_TPM)

    final_table[, c("pvals", "expected")] <- t(apply(final_table, 1, run_chisquare_test, total_pacbio, total_illumina))
    final_table$gene_name <- merged_illumina_pacbio$annot_gene_name
    final_table$pvals <- p.adjust(final_table$pvals, method = "bonferroni")

    volcano_plot(final_table, fill_color, opt$outdir)
    data <- plot_MA_observed_expected(final_table, fill_color, opt$outdir)

    printable <- subset(data, status == "Bonf. p-value <= 0.01")
    printable <- printable[,c("gene_name", "norm_illumina_TPM", "expected", "observed", "pvals", "A", "M")]
    colnames(printable) <- c("gene_name", "norm_illumina_TPM", "expected_TPM", "norm_observed_pacbio_TPM", "corrected_p-value", "A", "M")
    write.table(printable, paste(opt$outdir, "/MA_plot_gene_table.tsv", sep=""), row.names=F, col.names=T, quote=F, sep="\t")


    #plot_gene_longest_for_groups(data, opt$illumina_kallisto, opt$outdir)   
    #plot_gene_shortest_for_groups(data, opt$illumina_kallisto, opt$outdir)
    #plot_gene_avg_for_groups(data, opt$illumina_kallisto, opt$outdir)
}

volcano_plot <- function(data, fillcolor, outdir) {

    data$observed <- data$norm_both_pacbio
    data$pvals <- data$pvals + 10^-24

    data$log2FoldChange <- log2(data$observed + 1) - log2(data$expected + 1)
    data$fold_change = 2^data$log2FoldChange
    data$status <- as.factor(ifelse(abs(data$log2FoldChange) >= 1 & data$pvals <= 0.01, 
                             "Bonf. p-value <= 0.01", "Bonf. p-value > 0.01"))

    n_sig <- length(data$status[data$status == "Bonf. p-value <= 0.01"])
    n_no_sig <- length(data$status[data$status == "Bonf. p-value > 0.01"])

    fname <- paste(outdir, "/pacbio_illumina_gene_volcano_plot.png", sep="")
    xlabel <- "log2-fold change"
    ylabel <- "-log10 adjusted p-value"

    png(filename = fname,
        width = 2500, height = 2500, units = "px",
        bg = "white",  res = 300)

    g <- ggplot(data, aes(x=log2FoldChange, y=-log2(pvals), color = status)) +
         geom_point(alpha = 0.4, size = 2) +
         xlab(xlabel) + ylab(ylabel) + theme_bw() +
         coord_cartesian(xlim = c(-10,10)) +
         scale_color_manual(values = c("orange", fillcolor),
                                  labels = c(paste0("Significant (n = ", n_sig, ")"),
                                             paste0("Not significant (n = ", n_no_sig, ")"))) +
         theme(axis.text.x = element_text(color="black", size=20),
                     axis.text.y = element_text(color="black", size=20),
                     axis.title.x = element_text(color="black", size=16),
                     axis.title.y = element_text(color="black", size=16)) +
               guides(colour = guide_legend(override.aes = list(alpha=1,size=2.5))) +
               theme(legend.position=c(0.75,0.8),
                     legend.title = element_blank(),
                     legend.background = element_rect(fill="white", color = "black"),
                     legend.key = element_rect(fill="transparent"),
                     legend.text = element_text(colour = 'black', size = 16))        

    print(g)
    dev.off()

}


plot_MA_observed_expected <- function(data, fillcolor, outdir) {

    data$observed <- data$norm_both_pacbio

    data$M <- log2(data$observed + 1) - log2(data$expected + 1) # Same thing as the fold change
    data$A <- 0.5*(log2(data$observed + 1) + log2(data$expected + 1))

    data$status <- as.factor(ifelse(abs(data$M) >= 1 & data$pvals <= 0.01, "Bonf. p-value <= 0.01", "Bonf. p-value > 0.01"))
    n_sig <- length(data$status[data$status == "Bonf. p-value <= 0.01"])
    n_no_sig <- length(data$status[data$status == "Bonf. p-value > 0.01"])

    print(nrow(data))
    print(nrow(subset(data, status == "Bonf. p-value <= 0.01")))

    fname <- paste(outdir, "/gene_obs_expected_MA_plot.png", sep="")   
    xlabel <- "0.5*(log2(observed*expected PacBio counts))"
    ylabel <- "log2(ratio of observed to expected PacBio counts)"

    png(filename = fname,
     width = 2500, height = 2500, units = "px",
    bg = "white",  res = 300)

    p = ggplot(data, aes(x = A, y = M, color = status)) +
               geom_jitter(alpha = 0.4, size = 2) +
               xlab(xlabel) + ylab(ylabel) + theme_bw() +
               scale_color_manual(values = c("orange", fillcolor),
                                  labels = c(paste0("Significant (n = ", n_sig, ")"),
                                             paste0("Not significant (n = ", n_no_sig, ")"))) +
                                  #labels = c("Bonf. p-value <= 0.01 \nor log2 fold change > 1", "Bonf. p-value > 0.01")) +
               #coord_cartesian(xlim = c(0,1)) +
               theme(axis.text.x = element_text(color="black", size=20),
                     axis.text.y = element_text(color="black", size=20),
                     axis.title.x = element_text(color="black", size=16),
                     axis.title.y = element_text(color="black", size=16)) + 
               guides(colour = guide_legend(override.aes = list(size=2.5))) +
               theme(legend.position=c(0.75,0.1),
                     legend.title = element_blank(),
                     legend.background = element_rect(fill="white", color = "black"),
                     legend.key = element_rect(fill="transparent"),
                     legend.text = element_text(colour = 'black', size = 16))

    print(p)
    dev.off()

    return(data)
}

plot_gene_longest_for_groups <- function(data, kallisto_file, outdir) {
    # For the blues and the yellows separately, plot length distibution of genes.
    # Determine length by taking longest transcript of the gene 

    # Get table of gene lengths
    gencode.quantitation <- process_kallisto_table(kallisto_file)

    # Get longest transcript per gene
    gene_length <- aggregate(gencode.quantitation$length, by=list(gencode.quantitation$gene), FUN=max)
    colnames(gene_length) <- c("gene_name", "length")
    gene_length$length <- gene_length$length/1000

    # Merge length with data
    m <- merge(data, gene_length, by = "gene_name", all.x = T, all.y = F)
    fname <- paste(outdir, "/gene_length_by_MA_group_longest.png", sep="")
    xlabel <- "Gene length in kb (determined by longest known transcript)"
    plot_gene_length_dist(m, fname, xlabel)
}

plot_gene_shortest_for_groups <- function(data, kallisto_file, outdir) {
    # For the blues and the yellows separately, plot length distibution of genes.
    # Determine length by taking shortest transcript of the gene

    # Get table of gene lengths
    gencode.quantitation <- process_kallisto_table(kallisto_file)

    # Get longest transcript per gene
    gene_length <- aggregate(gencode.quantitation$length, by=list(gencode.quantitation$gene), FUN=min)
    colnames(gene_length) <- c("gene_name", "length")
    gene_length$length <- gene_length$length/1000

    # Merge length with data
    m <- merge(data, gene_length, by = "gene_name", all.x = T, all.y = F)
    fname <- paste(outdir, "/gene_length_by_MA_group_shortest.png", sep="")
    xlabel <- "Gene length in kb (determined by shortest known transcript)"
    plot_gene_length_dist(m, fname, xlabel)
}

plot_gene_avg_for_groups <- function(data, kallisto_file, outdir) {
    # For the blues and the yellows separately, plot length distibution of genes.
    # Determine length by taking avg transcript length

    # Get table of gene lengths
    gencode.quantitation <- process_kallisto_table(kallisto_file)

    # Get longest transcript per gene
    gene_length <- aggregate(gencode.quantitation$length, by=list(gencode.quantitation$gene), FUN=mean)
    colnames(gene_length) <- c("gene_name", "length")
    gene_length$length <- gene_length$length/1000

    # Merge length with data
    m <- merge(data, gene_length, by = "gene_name", all.x = T, all.y = F)
    fname <- paste(outdir, "/gene_length_by_MA_group_avg.png", sep="")
    xlabel <- "Gene length in kb (average of known transcripts)"
    plot_gene_length_dist(m, fname, xlabel)
}

plot_gene_length_dist <- function(data, fname, xlabel) {

    png(filename = fname,
     width = 2500, height = 2000, units = "px",
    bg = "white",  res = 300)

    p = ggplot(data, aes(length, color = status, fill = status)) +
               geom_density(alpha = 0.4, position = "identity") +
               xlab(xlabel) + theme_bw() +
               scale_color_manual(values = c("orange", "navy")) +
               scale_fill_manual(values = c("orange", "navy")) +
               #coord_cartesian(xlim = c(0,1)) +
               theme(axis.text.x = element_text(color="black", size=20),
                     axis.text.y = element_text(color="black", size=20),
                     axis.title.x = element_text(color="black", size=16),
                     axis.title.y = element_text(color="black", size=16))

    print(p)
    dev.off()

}


process_kallisto_table <- function(kallisto_file) {
    # No filtering
    # Get table of gene lengths
    gencode.quantitation <- as.data.frame(read_delim(kallisto_file, "\t", escape_double = FALSE,
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))

    # Split GENCODE transcript multi-id by '|'
    extraCols =str_split_fixed(gencode.quantitation$target_id, "\\|",9)[,c(5,6,8)]
    colnames(extraCols) = c("transcript", "gene", "class")
    gencode.quantitation = cbind(extraCols, gencode.quantitation)

    return(gencode.quantitation)
}


plot_pvalue_distribution <- function(data, fillcolor, outdir) {
    fname <- paste(outdir, "/gene_chisquare_pvalue_density.png", sep="")
    xlabel <- "P-Value"
    ylabel <- "Gene count"


#    data <- subset(data, illumina_TPM > 100 & both_pacbio == 0)
#    print(data)

    png(filename = fname,
     width = 2000, height = 2500, units = "px",
    bg = "white",  res = 300)

    top = 5
    p = ggplot(data, aes(pvals)) +
               geom_density(alpha = 0.2, fill = fillcolor, color = fillcolor) +
               xlab(xlabel) + ylab(ylabel) + theme_bw() +
               coord_cartesian(xlim = c(0,1)) +
               theme(axis.text.x = element_text(color="black", size=22),
                     axis.text.y = element_text(color="black", size=22),
                     axis.title.x = element_text(color="black", size=22),
                     axis.title.y = element_text(color="black", size=22))

    print(p)
    dev.off()
}

run_chisquare_test <- function(reads_vector, total_pacbio, total_illumina) {
    reads_pacbio <- reads_vector[1]
    reads_illumina <- reads_vector[2]

    M <- as.table(rbind(c(reads_pacbio, total_pacbio - reads_pacbio), 
                        c(reads_illumina, total_illumina - reads_illumina)))
    dimnames(M) <- list(platform = c("PacBio", "Illumina"),
                    gene = c("query_gene", "not_query_gene")) 
    Xsq <- chisq.test(M)
    return(c(Xsq$p.value, Xsq$expected[1,1]))
}

load_packages <- function() {
    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("plyr"))
    suppressPackageStartupMessages(library("Hmisc"))
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("readr"))
    suppressPackageStartupMessages(library("reshape"))
    suppressPackageStartupMessages(library("stringr"))
    suppressPackageStartupMessages(library("data.table"))
    suppressPackageStartupMessages(library("preprocessCore"))

    # Load my custom functions
    source("/pub/dwyman/TALON-paper-2019/analysis_scripts/filter_kallisto_illumina_genes.R")
    source("/pub/dwyman/TALON-paper-2019/analysis_scripts/filter_kallisto_illumina_transcripts.R")

    return
}

parse_options <- function() {

    option_list <- list(
    make_option(c("--f"), action = "store", dest = "infile",
                    default = NULL, help = "TALON abundance file (filtered)"),
    make_option(c("--datasets"), action = "store", dest = "datasets",
                    default = NULL, help = "Comma-delimited list of two dataset names to include in the analysis."),
    make_option(c("--ik"), action = "store", dest = "illumina_kallisto",
                    default = NULL, help = "Illumina Kallisto file."),
    make_option(c("--color"), action = "store", dest = "color_scheme",
                    default = NULL, help = "blue, red, or green"),
    make_option(c("-o","--outdir"), action = "store", dest = "outdir",
                default = NULL, help = "Output directory for plots and outfiles"))

    opt <- parse_args(OptionParser(option_list=option_list))
    return(opt)
}


main()
