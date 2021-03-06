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
        fill_color <- "#009E73"
    }

    # Get the names of the first and second dataset that we will be working with
    data_names <- str_split(opt$datasets, ",")[[1]]
    dataset1 <- data_names[1]
    dataset2 <- data_names[2]

    # datatype label (probably ONT or PacBio)
    dtype <- opt$dtype

    # Get transcripts expressed in the Illumina data from the Kallisto
    # abundance files
    illumina_1 <- filter_kallisto_illumina_transcripts(opt$illumina_kallisto_1)[,c("t_ID","tpm", "length")]
    illumina_2 <- filter_kallisto_illumina_transcripts(opt$illumina_kallisto_2)[,c("t_ID","tpm", "length")]
    colnames(illumina_1)[1:2] <- c("annot_transcript_id", "illumina_TPM_1")
    colnames(illumina_2)[1:2] <- c("annot_transcript_id", "illumina_TPM_2")
    illumina_transcript_table <- merge(illumina_1, illumina_2, by = "annot_transcript_id",
                                 all.x = T, all.y = T)

    illumina_transcript_table[is.na(illumina_transcript_table)] <- 0
    illumina_transcript_table$length <- pmax(illumina_transcript_table$length.x, 
                                             illumina_transcript_table$length.y)
    illumina_transcript_table$length.x <- NULL
    illumina_transcript_table$length.y <- NULL
    illumina_transcript_table$illumina_TPM_1 <- round(illumina_transcript_table$illumina_TPM_1)
    illumina_transcript_table$illumina_TPM_2 <- round(illumina_transcript_table$illumina_TPM_2)


    # Read PacBio abundance file
    abundance <- as.data.frame(read_delim(opt$infile, "\t", escape_double = FALSE,
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))

    # Keep known transcripts only
    pb_transcript_abundance <- subset(abundance, transcript_novelty == "Known")

    # Cut out unnecessary cols
    pb_transcript_abundance <- pb_transcript_abundance[, c("annot_transcript_id", dataset1, dataset2)]

    # Merge PacBio with Illumina on annot_transcript_id
    merged_illumina_pacbio <- merge(illumina_transcript_table, pb_transcript_abundance, by = "annot_transcript_id",
                                    all.x = T, all.y = T)
    merged_illumina_pacbio[is.na(merged_illumina_pacbio)] <- 0
    merged_illumina_pacbio <- merged_illumina_pacbio[, c("annot_transcript_id", 
                                                         "illumina_TPM_1", 
                                                         "illumina_TPM_2",  
                                                         dataset1, dataset2)]

    # Now, format table for edgeR by setting the rownames to the transcript names
    rownames(merged_illumina_pacbio) <- merged_illumina_pacbio$annot_transcript_id
    merged_illumina_pacbio$annot_transcript_id <- NULL

    # Remove rows that only contain zeros
    merged_illumina_pacbio <- merged_illumina_pacbio[rowSums(merged_illumina_pacbio) > 0, ]

    # edgeR basics: 
    group <- factor(c("1","1","2","2")) # Indicate which group each col belongs to
    y <- DGEList(counts=merged_illumina_pacbio, group = group) # Create a DGEList object
    y <- calcNormFactors(y) # Normalize counts in the object
    design <- model.matrix(~group)
    y <- estimateDisp(y,design)

    # Pairwise testing approach for DE Genes. "classic" edgeR
    et <- exactTest(y, pair=c("1","2"))

    # Extract exact test table for plotting
    illumina_PB_et <- et$table
    illumina_PB_et$annot_transcript_id <- rownames(illumina_PB_et)
    illumina_PB_et <- merge(illumina_PB_et, 
                            abundance[, c("annot_transcript_id", "annot_transcript_name")], 
                            by = "annot_transcript_id", all.x = T, all.y = F)
    

    # Adjust p-values
    illumina_PB_et$adj_pval <- p.adjust(illumina_PB_et$PValue, method = "bonferroni")

    # MA plot
    ma_plot(illumina_PB_et, fill_color, opt$outdir, opt$dtype, opt$xmax, opt$ymax)

    # Merge the EdgeR table with the other information
    illumina_PB_et <- cbind(illumina_PB_et, merged_illumina_pacbio)
    illumina_PB_et <- illumina_PB_et[order(illumina_PB_et$adj_pval),]
    write.table(illumina_PB_et, paste(opt$outdir, "/edgeR_", dtype, "_illumina_transcripts.tsv", sep=""),
                row.names=F, col.names=T, quote=F)
}


ma_plot <- function(data, fillcolor, outdir, dtype, xmax, ymax) {

    data$status <- as.factor(ifelse(abs(data$logFC) > 1 & data$adj_pval <= 0.01,
                             "Bonf. p-value <= 0.01", "Bonf. p-value > 0.01"))

    n_sig <- length(data$status[data$status == "Bonf. p-value <= 0.01"])
    n_no_sig <- length(data$status[data$status == "Bonf. p-value > 0.01"])

    fname <- paste(outdir, "/edgeR_", dtype, "_illumina_transcript_MA_plot.png", sep="")
    xlabel <- "log(Counts per million)"
    ylabel <- paste0(dtype, " to Illumina log2-fold change")

    png(filename = fname,
        width = 3000, height = 3000, units = "px",
        bg = "white",  res = 300)

    g <- ggplot(data, aes(x=logCPM, y=logFC, color = status)) +
         geom_point(alpha = 0.4, size = 2) +
         xlab(xlabel) + ylab(ylabel) + 
         coord_cartesian(xlim=c(0,xmax), ylim = c(-1*ymax,ymax)) +
         scale_color_manual(values = c("orange", fillcolor),
                                  labels = c(paste0("Sig. (n = ", n_sig, ")"),
                                             paste0("Not sig. (n = ", n_no_sig, ")"))) +
         theme_bw(base_family = "Helvetica", base_size = 18) +
         theme(axis.text.x = element_text(color="black", size=rel(2.5)),
               axis.text.y = element_text(color="black", size=rel(2.5)),
               axis.title.x = element_text(color="black", size=rel(2)),
               axis.title.y = element_text(color="black", size=rel(2))) +
         guides(colour = guide_legend(override.aes = list(alpha=1,size=2.5))) +
         theme(legend.position=c(0.28,0.1),
               legend.title = element_blank(),
               legend.background = element_rect(fill="white", color = "black"),
               legend.key = element_rect(fill="transparent"),
               legend.text = element_text(colour = 'black', size = rel(1.5)))

    print(g)
    dev.off()
}

filter_kallisto_illumina_transcripts <- function(kallisto_file) {
    # This function takes a Kallisto abundance file and filters the transcripts
    # based on criteria designed to make the gene set comparable to what can
    # be detected using PacBio

    gencode.quantitation <- as.data.frame(read_delim(kallisto_file, "\t", escape_double = FALSE,
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))

    # Split GENCODE transcript multi-id by '|'
    extraCols <- str_split_fixed(gencode.quantitation$target_id, "\\|",9)[,c(5,6,8,1,2)]
    colnames(extraCols) <- c("transcript", "gene", "class", "t_ID", "g_ID")
    gencode.quantitation <- cbind(extraCols, gencode.quantitation)

    # Remove transcripts that are < 300 bp in length because PacBio chucks anything that size
    # Also require TPM > 1
    filtered_transcripts <- subset(gencode.quantitation, length >= 300 & tpm > 1)


    mitochondrial_blacklist <- c("MT-TF", "MT-RNR1", "MT-TV", "MT-RNR2", "MT-TL1", "MT-ND1", "MT-TI", "MT-TQ", "MT-TM", "MT-ND2", "MT-TW", "MT-TA", "MT-TN", "MT-TC", "MT-TY", "MT-CO1", "MT-TS1", "MT-TD", "MT-CO2", "MT-TK", "MT-ATP8", "MT-ATP6", "MT-CO3", "MT-TG", "MT-ND3", "MT-TR", "MT-ND4L", "MT-ND4", "MT-TH", "MT-TS2", "MT-TL2", "MT-ND5", "MT-ND6", "MT-CYB", "MTATP6P1")
    final_filtered_transcripts <- subset(filtered_transcripts, !(gene %in% mitochondrial_blacklist))

    # Normalize back to 1 million transcripts
    total_transcripts <- sum(final_filtered_transcripts$tpm)
    TPM_scaling <- 1000000/total_transcripts
    final_filtered_transcripts$tpm <- final_filtered_transcripts$tpm*TPM_scaling

    return(final_filtered_transcripts)
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
    suppressPackageStartupMessages(library("edgeR"))

    return
}

parse_options <- function() {

    option_list <- list(
    make_option(c("--f"), action = "store", dest = "infile",
                    default = NULL, help = "TALON abundance file (not filtered)"),
    make_option(c("--datasets"), action = "store", dest = "datasets",
                    default = NULL, help = "Comma-delimited list of two dataset names to include in the analysis."),
    make_option(c("--ik1"), action = "store", dest = "illumina_kallisto_1",
                    default = NULL, help = "Rep1 Illumina Kallisto file."),
    make_option(c("--ik2"), action = "store", dest = "illumina_kallisto_2",
                    default = NULL, help = "Rep2 Illumina Kallisto file."),
    make_option(c("--color"), action = "store", dest = "color_scheme",
                    default = NULL, help = "blue, red, or green"),
    make_option(c("-o","--outdir"), action = "store", dest = "outdir",
                default = NULL, help = "Output directory for plots and outfiles"),
    make_option(c("--dtype"), action = "store", dest = "dtype",
                default = "PacBio", help = "Datatype label to display on plot"),
    make_option(c("--xmax"), action = "store", dest = "xmax",
                default = 16, help = "Max x-value for plot (default = 16)"),
    make_option(c("--ymax"), action = "store", dest = "ymax",
                default = 20, help = "Max y-value for plot (default = 20)"))

    opt <- parse_args(OptionParser(option_list=option_list))
    return(opt)
}


main()
