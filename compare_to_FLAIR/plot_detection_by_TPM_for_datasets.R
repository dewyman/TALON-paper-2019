
main <-function() {

    set.seed(100)
    load_packages()
    opt <- parse_options()

    # Get colors
    if (opt$color_scheme == "red") {
        color_vec <- c("white", "orange", "red2")
    } else if (opt$color_scheme == "blue") {
        color_vec <- c("white", "skyblue", "navy")
    } else if (opt$color_scheme == "green") {
        color_vec <- c("white", "yellowgreen", "springgreen4")
    }

    # Get genes and transcripts expressed in the Illumina data from the Kallisto
    # abundance file
    illumina_1 <- filter_kallisto_illumina_genes(opt$illumina_kallisto_1)
    illumina_2 <- filter_kallisto_illumina_genes(opt$illumina_kallisto_2)
    colnames(illumina_1) <- c("gene", "illumina_TPM_1")
    colnames(illumina_2) <- c("gene", "illumina_TPM_2")
    illumina_gene_table <- merge(illumina_1, illumina_2, by = "gene",
                                 all.x = T, all.y = T)

    # Remove genes only detected in one Illumina rep, then average the TPMS
    illumina_gene_table <- illumina_gene_table[complete.cases(illumina_gene_table), ]
    illumina_gene_table$tpm <- (illumina_gene_table$illumina_TPM_1 +
                                illumina_gene_table$illumina_TPM_2) / 2

    
    # Get the names of the first and second dataset that we will be working with
    data_names <- str_split(opt$datasets, ",")[[1]]
    dataset_1 <- data_names[1]
    dataset_2 <- data_names[2]

    # Read in abundance file
    abundance_table <- as.data.frame(read_delim(opt$infile, delim = "\t",
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))

    # Get genes in each dataset
    d1_genes <- get_detected_genes_for_dataset(abundance_table, dataset_1)
    d2_genes <- get_detected_genes_for_dataset(abundance_table, dataset_2)
    

    # Combine the Illumina tables with information about which genes/transcripts are observed in Pacbio
    illumina_gene_detection <- get_detection(illumina_gene_table, d1_genes, d2_genes, "gene")

    #write.table(illumina_gene_detection[,c("gene","detection")], paste(opt$outdir, "/", "gene_detection.csv", sep=""),
    #            col.names=F, row.names=F, quote=F, sep=",")
    

    # Group into buckets
    illumina_gene_detection_buckets <- get_buckets(illumina_gene_detection)    

    # Plot detection by TPM
    plot_detection(illumina_gene_detection_buckets$illumina, illumina_gene_detection_buckets$interval_labels, "gene", color_vec, opt$outdir)
    print ("--------------------------")

   #plot_length_hists_by_detection(highly_expressed, opt$outdir)
}

get_buckets <- function(illumina) {
    intervals = c(1, 2, seq(5, 20, 5), 50, 100, 500)
    illumina$group = cut(illumina$tpm, c(intervals, 100000000000))   

    cIntervals = as.character(intervals)
    interval_starts_with_dash = paste0(cIntervals, c(rep("-", length(cIntervals)-1), "+"))
    full_intervals = paste0(interval_starts_with_dash, c(as.character(intervals)[2:length(intervals)], ""))

    return(list("illumina" = illumina, "interval_labels" = full_intervals))
}


get_detection <- function(illumina, d1_items, d2_items, cat_type) {

    shared <- intersect(d1_items, d2_items)
    d1_d2_union <- union(d1_items, d2_items) 
    not_shared <- d1_d2_union[(d1_d2_union %in% shared) == F]

    illumina$detection = ifelse(illumina[,cat_type] %in% not_shared, "Detected in one PacBio rep",
                                             ifelse(illumina[,cat_type] %in% shared, "Detected in both PacBio reps", "Not detected in PacBio"))
    return(illumina)
}

plot_detection <- function(illumina, cIntervals, cat_type, color_vec, outdir) {

    # Plot the curves
    fname <- paste(outdir, "/", cat_type, "_detection_by_TPM.png", sep="")
    xlabel <- paste(capitalize(cat_type), "expression level in Illumina data (TPM)")
    ylabel <- paste("Fraction of ", cat_type, "s", sep="")

    png(filename = fname,
        width = 2100, height = 2500, units = "px",
        bg = "white",  res = 300)
    g = ggplot(illumina, aes(x = group, fill = factor(detection, levels = c("Not detected in PacBio", "Detected in one PacBio rep", "Detected in both PacBio reps")))) +
       geom_bar(position = "fill", col = "black") + 
       scale_fill_manual("",values=color_vec)  + 
       theme_bw() + 
       scale_y_continuous(breaks=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) + 
       theme(axis.text.x = element_text(angle = 30, hjust = 1, color = "black", size=20)) + 
       xlab(xlabel)  + ylab(ylabel) + theme(text= element_text(size=22)) + 
       theme(axis.text.x = element_text(color = "black", size=21),
             axis.text.y = element_text(color = "black", size=21)) + 
       scale_x_discrete(labels=cIntervals) + 
       theme(legend.position=c(0.6,0.2),
              legend.title = element_blank(),
              legend.background = element_rect(fill="white", color = "black"),
              legend.key = element_rect(fill="transparent"),
              legend.text = element_text(colour = 'black', size = 20),
              plot.margin = margin(t = 0.5, r = 1.5, l = 0.5, b = 0.5, "cm"))

    print(g)
    dev.off()

    #print(illumina %>% group_by("detection") %>% count())

}

plot_length_hists_by_detection <- function(illumina, outdir) {

    fname <- paste(outdir, "/transcript_length_by_detection_status_100plusTPM.png", sep="")
    xlabel <- "Transcript length (kilobases)"
    ylabel <- "Count"
    illumina$length <- illumina$length/1000

    png(filename = fname,
        width = 3000, height = 2000, units = "px",
        bg = "white",  res = 300)
    g = ggplot(illumina, aes(x = length, color = factor(detection, levels = c("Not detected in PacBio", "Detected in one PacBio rep", "Detected in both PacBio reps")))) +
       geom_freqpoly(binwidth = 0.1) + 
       scale_colour_manual("",values=c("red", "skyblue", "navy"))  +
       theme_bw() + 
       theme(axis.text.x = element_text(color = "black", size=17)) +
       xlab(xlabel)  + ylab(ylabel) + theme(text= element_text(size=17)) +
       theme(axis.text.y = element_text(color = "black", size=17)) + 
       coord_cartesian(xlim=c(0, 7), ylim=c(0,60)) + 
       scale_x_continuous(breaks=seq(0,7, by = 1))

    print(g)
    dev.off()
}

get_detected_genes_for_dataset <- function(abundance, dataset) {
    data <- abundance[,c("annot_gene_id", dataset)]
    observed_genes <- unique(data[data[,dataset] > 0, "annot_gene_id"])

    return(observed_genes)
}


filter_kallisto_illumina_genes <- function(kallisto_file) {
    # This function takes a Kallisto abundance file and filters the genes
    # based on criteria designed to make the gene set comparable to what can
    # be detected using PacBio

    gencode.quantitation <- as.data.frame(read_delim(kallisto_file, "\t", escape_double = FALSE,
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))

    # Split GENCODE transcript multi-id by '|'
    extraCols =str_split_fixed(gencode.quantitation$target_id, "\\|",9)[,c(5,6,8,1,2)]
    colnames(extraCols) = c("transcript", "gene", "class", "t_ID", "g_ID")
    gencode.quantitation = cbind(extraCols, gencode.quantitation)

    # Remove transcripts that are < 300 bp in length because PacBio chucks anything that size
    gencode_quant_min300 <- subset(gencode.quantitation, length >= 300)

    # Remove genes that are on the mitochondrial blacklist
    mitochondrial_blacklist <- c("MT-TF", "MT-RNR1", "MT-TV", "MT-RNR2", "MT-TL1", 
                                 "MT-ND1", "MT-TI", "MT-TQ", "MT-TM", "MT-ND2", 
                                 "MT-TW", "MT-TA", "MT-TN", "MT-TC", "MT-TY", 
                                 "MT-CO1", "MT-TS1", "MT-TD", "MT-CO2", "MT-TK", 
                                 "MT-ATP8", "MT-ATP6", "MT-CO3", "MT-TG", "MT-ND3", 
                                 "MT-TR", "MT-ND4L", "MT-ND4", "MT-TH", "MT-TS2", 
                                 "MT-TL2", "MT-ND5", "MT-ND6", "MT-CYB","MTATP6P1")
    gencode_quant_min300_noMT <- subset(gencode_quant_min300, !(gene %in% mitochondrial_blacklist))

    # Aggregate by gene
    gene_gencode_quant_min300_noMT <- aggregate(gencode_quant_min300_noMT$tpm, by=list(gencode_quant_min300_noMT$g_ID), FUN=sum)
    colnames(gene_gencode_quant_min300_noMT) <- c("gene", "tpm")

    # Constraints: > 300 bp, TPM > 1
    final_filtered_genes <- gene_gencode_quant_min300_noMT[gene_gencode_quant_min300_noMT$tpm > 1,]

    # Normalize back to 1 million transcripts
    total_transcripts <- sum(final_filtered_genes$tpm)
    TPM_scaling <- 1000000/total_transcripts
    final_filtered_genes$tpm <- final_filtered_genes$tpm*TPM_scaling

    return(final_filtered_genes)
}

load_packages <- function() {
    suppressPackageStartupMessages(library("DBI"))
    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("plyr"))
    suppressPackageStartupMessages(library("Hmisc"))
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("RSQLite"))
    suppressPackageStartupMessages(library("tidyverse"))
    suppressPackageStartupMessages(library("reshape"))
    suppressPackageStartupMessages(library("stringr"))
    suppressPackageStartupMessages(library("data.table"))

    return
}


parse_options <- function() {

    option_list <- list(
        make_option(c("--f"), action = "store", dest = "infile",
                    default = NULL, help = "Unfiltered TALON abundance file"),
        make_option(c("--datasets"), action = "store", dest = "datasets",
                    default = NULL, help = "Comma-delimited list of two dataset names to include in the analysis."),
        make_option(c("--ik1"), action = "store", dest = "illumina_kallisto_1",
                    default = NULL, help = "Rep1 Illumina Kallisto file."),
        make_option(c("--ik2"), action = "store", dest = "illumina_kallisto_2",
                    default = NULL, help = "Rep2 Illumina Kallisto file."),
        make_option(c("--color"), action = "store", dest = "color_scheme",
                    default = NULL, help = "blue, red, or green"),
        make_option(c("-o","--outdir"), action = "store", dest = "outdir",
                    default = NULL, help = "Output directory for plots and outfiles")
        )

    opt <- parse_args(OptionParser(option_list=option_list))
    return(opt)
}


main()
