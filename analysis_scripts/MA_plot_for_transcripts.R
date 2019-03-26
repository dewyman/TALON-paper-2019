# Generate an MA plot that compares the expected isoform-wise quantification
# base on short reads to the actual values measured in PacBio.
# Since we are looking at GENCODE-annotated transcripts only, no filtering
# of the PacBio data is needed.

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

    # Get the transcripts expressed in the Illumina data from Kallisto
    illumina_table <- filter_kallisto_illumina_transcripts(opt$illumina_kallisto)
    illumina_table <- illumina_table[,c(1,ncol(illumina_table))]
    colnames(illumina_table) <- c("annot_transcript_name", "illumina_TPM")


    # Read PacBio abundance file
    pb_abundance <- as.data.frame(read_delim(opt$infile, "\t", escape_double = FALSE,
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))
    pb_abundance <- subset(pb_abundance, transcript_status == "KNOWN")
    pb_abundance <- pb_abundance[, c("annot_transcript_name", dataset1, dataset2)] 

    # Sum together the PacBio abundances
    pb_abundance$both_pacbio <- pb_abundance[,dataset1] + pb_abundance[,dataset2]
    total_pacbio_reads <- sum(pb_abundance[,dataset1]) + sum(pb_abundance[,dataset2])

    # Merge PacBio with Illumina on annot_transcript_name
    # TODO: should we restrict to Illumina-only transcripts?
    merged_illumina_pacbio <- merge(illumina_table, pb_abundance, 
                                    by = "annot_transcript_name",
                                    all.x = T, all.y = F)
    merged_illumina_pacbio[is.na(merged_illumina_pacbio)] <- 0
    print(nrow(merged_illumina_pacbio))

    final_table <- merged_illumina_pacbio[, c("both_pacbio", "illumina_TPM")]

    # Compute the p-values
    final_table$illumina_TPM <- round(final_table$illumina_TPM)
    total_illumina <- sum(final_table$illumina_TPM) 
    final_table[, c("pvals", "expected")] <- t(apply(final_table, 1, 
                                                     run_chisquare_test, 
                                                     total_pacbio_reads, 
                                                     total_illumina))
    final_table$transcript_name <- merged_illumina_pacbio$annot_transcript_name
    final_table$pvals <- p.adjust(final_table$pvals, method = "bonferroni")

    data <- plot_MA_observed_expected(final_table, fill_color, opt$outdir)

    printable <- subset(data, status == "Bonf. p-value <= 0.01")
    printable <- printable[,c("transcript_name", "illumina_TPM", "expected", 
                              "observed", "pvals", "A", "M")]
    colnames(printable) <- c("transcript_name", "illumina_TPM", "expected_TPM", 
                             "observed_pacbio_TPM", "corrected_p-value", "A", "M")
    write.table(printable, paste(opt$outdir, "/MA_plot_gene_table.tsv", sep=""), 
                row.names=F, col.names=T, quote=F, sep="\t")

}

plot_MA_observed_expected <- function(data, fillcolor, outdir) {

    # Perform quantile normalization
    counts <- as.matrix(data[, c("both_pacbio", "expected")])
    counts <- as.data.frame(normalize.quantiles(counts))
    data$observed <- counts[,1] + 1
    data$expected <- counts[,2] + 1

    data$M <- log(data$observed, base=2) - log(data$expected, base=2)
    data$A <- 0.5*(log(data$observed, base=2) + log(data$expected, base=2))

    data$fold_change <- (data$observed - data$expected) / data$expected
    data$status <- as.factor(ifelse(abs(data$M) >= 1 & data$pvals <= 0.01, "Bonf. p-value <= 0.01", "Bonf. p-value > 0.01"))

    print(nrow(data))
    print(nrow(subset(data, status == "Bonf. p-value <= 0.01")))

    fname <- paste(outdir, "/transcript_obs_expected_MA_plot.png", sep="")
    xlabel <- "0.5*(log2(observed*expected PacBio counts))"
    ylabel <- "log2(ratio of observed to expected PacBio counts)"

    png(filename = fname,
     width = 2500, height = 2500, units = "px",
    bg = "white",  res = 300)

    p = ggplot(data, aes(x = A, y = M, color = status)) +
               geom_jitter(alpha = 0.4, size = 2.5) +
               xlab(xlabel) + ylab(ylabel) + theme_bw() +
               scale_color_manual(values = c("orange", fillcolor),
                                  labels = c("Significant", "Not significant")) +
                                  #labels = c("Bonf. p-value <= 0.01 \nor log2 fold change > 1", "Bonf. p-value > 0.01")) +
               guides(colour = guide_legend(override.aes = list(size=2.5))) +
               theme(axis.text.x = element_text(color="black", size=20),
                     axis.text.y = element_text(color="black", size=20),
                     axis.title.x = element_text(color="black", size=16),
                     axis.title.y = element_text(color="black", size=16)) +
               theme(legend.position=c(0.8,0.2),
                     legend.title = element_blank(),
                     legend.background = element_rect(fill="white", color = "black"),
                     legend.key = element_rect(fill="transparent"),
                     legend.text = element_text(colour = 'black', size = 16))

    print(p)
    dev.off()

    return(data)
}

run_chisquare_test <- function(reads_vector, total_pacbio, total_illumina) {
    reads_pacbio <- reads_vector[1]
    reads_illumina <- reads_vector[2]

    M <- as.table(rbind(c(reads_pacbio, total_pacbio - reads_pacbio),
                        c(reads_illumina, total_illumina - reads_illumina)))
    dimnames(M) <- list(platform = c("PacBio", "Illumina"),
                    transcript = c("query_transcript", "not_query_transcript"))
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
    #source("/pub/dwyman/TALON-paper-2019/analysis_scripts/filter_kallisto_illumina_genes.R")
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
