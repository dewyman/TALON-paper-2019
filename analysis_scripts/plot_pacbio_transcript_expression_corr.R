main <-function() {

    load_packages()
    opt <- parse_options()

    # Get colors
    if (opt$color_scheme == "red") {
        color_vec <- c("red2", "orange")
    } else if (opt$color_scheme == "blue") {
        color_vec <- c("navy", "orange")
    } else if (opt$color_scheme == "green") {
        color_vec <- c("springgreen4", "orange")
    }

    # Read abundance table
    abundance_table <- as.data.frame(read_delim(opt$infile, delim = "\t",
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))

    # Get dataset names and check that they exist in the file
    first_dataset <- opt$d1
    second_dataset <- opt$d2
    if (first_dataset %in% colnames(abundance_table) == F |
        second_dataset %in% colnames(abundance_table) == F) {
        print("One of the provided dataset names is not in the abundance file provided. Exiting...")
        quit()
    }

    # Read in whitelist file, and get gene and transcript whitelists from that
    whitelist <- as.data.frame(read_delim(opt$whitelist, ",", escape_double = FALSE,
                                  col_names = FALSE, trim_ws = TRUE, na = "NA"))
    whitelisted_gene_IDs <- unique(whitelist[,1])
    whitelisted_transcript_IDs <- whitelist[,2]
    
    # Plot expression scatterplots
    expression_by_status(abundance_table, first_dataset, second_dataset, whitelisted_gene_IDs, "gene", opt$outdir, color_vec)
    expression_by_status(abundance_table, first_dataset, second_dataset, whitelisted_transcript_IDs, "transcript", opt$outdir, color_vec)
}

expression_by_status <- function(abundance_table, d1, d2, whitelist, cat_type, outdir, color_vec) {

    cat_colname <- paste(cat_type, "ID", sep = "_")
    annot_colname <- paste(cat_type, "status", sep = "_")
    abundance_table <- abundance_table[abundance_table[,cat_colname] %in% whitelist,]
    annotation <- unique(abundance_table[,c(cat_colname, annot_colname)])

    d1_abundance <- abundance_table[,c(cat_colname,d1)]
    d2_abundance <- abundance_table[,c(cat_colname,d2)]
    colnames(d1_abundance) <- c(cat_colname, "count")
    colnames(d2_abundance) <- c(cat_colname, "count")

    # Convert counts to TPMs
    d1_abundance$TPM <- (d1_abundance$count)*1000000/sum(d1_abundance$count)
    d2_abundance$TPM <- (d2_abundance$count)*1000000/sum(d2_abundance$count)

    # Aggregate table by ID. This will only matter for genes, where we are summing TPMs
    if (cat_type == "gene") {
        d1_abundance_agg <- aggregate(d1_abundance$TPM, by=list(d1_abundance[,cat_colname]), FUN=sum)
        colnames(d1_abundance_agg) <- c(cat_colname, "TPM")
        d2_abundance_agg <- aggregate(d2_abundance$TPM, by=list(d2_abundance[,cat_colname]), FUN=sum)
        colnames(d2_abundance_agg) <- c(cat_colname, "TPM")
    } else {
        d1_abundance_agg = d1_abundance[,c(cat_colname, "TPM")]
        d2_abundance_agg = d2_abundance[,c(cat_colname, "TPM")]
    }

    # Merge dataset 1 and dataset 2 back together
    merged_abundances <- merge(d1_abundance_agg, d2_abundance_agg, by = cat_colname, all.x = T, all.y = T)
    colnames(merged_abundances) <- c(cat_colname, "data1.TPM", "data2.TPM")
    merged_abundances[is.na(merged_abundances)] <- 0

    # Take log2(TPM + 1)
    merged_abundances$data1.TPM = log(merged_abundances$data1.TPM + 1, base=2)
    merged_abundances$data2.TPM = log(merged_abundances$data2.TPM + 1, base=2)

    # Merge back in the known/novel info since that goes away in the aggregate
    annotated_merged_abundance <- merge(merged_abundances, annotation, by = cat_colname, all.x = T, all.Y = F)
    print(head(annotated_merged_abundance))
    
    # Break down known/novel categories further
    status_col = paste(cat_type, "status", sep="_")
    source_col = paste(cat_type, "status", "source", sep="_")
    d1_label = paste("Novel (", d1, ")", sep="")
    d2_label = paste("Novel (", d2, ")", sep="")
    annotated_merged_abundance$novelty = ifelse(annotated_merged_abundance[,annot_colname] == "KNOWN" , "Known",
                                             ifelse(annotated_merged_abundance[,annot_colname] == "NOVEL", "Novel", "Other"))

    # Plot log2(TPM + 1) for each dataset on a scatterplot. Color points according to known/novel status
    pearsonCorr = cor.test(~data1.TPM + data2.TPM, data=annotated_merged_abundance, method = "pearson", continuity = FALSE, conf.level = 0.95)$estimate
    spearmanCorr = cor.test(~data1.TPM + data2.TPM, data=annotated_merged_abundance, method = "spearman", continuity = FALSE, conf.level = 0.95)$estimate

    joined_names <- paste(outdir, "/", d1, "-", d2, sep = "")
    fname <- paste(joined_names, cat_type, "correlationPlot.png", sep="_")
    xlabel <- paste("log2(TPM+1) of ", cat_type, " in ", "Rep1", sep="")
    ylabel <- paste("log2(TPM+1) of ", cat_type, " in ", "Rep2", sep="")

    png(filename = fname,
        width = 2500, height = 2000, units = "px",
        bg = "white",  res = 300)
    g = ggplot(annotated_merged_abundance, aes(x = data1.TPM, y = data2.TPM, color = novelty)) +
        geom_jitter(alpha = 0.5) + theme_bw() +
        xlab(xlabel)  + ylab(ylabel) + theme(text= element_text(size=24)) +
        theme(axis.text.x = element_text(color = "black", size=24),
              axis.text.y = element_text(color = "black", size=24)) +
        annotate("text", x = 5, y = 14, label = paste("Pearson r: ",
                 round(pearsonCorr, 3), "\nSpearman rho: ",
                 round(spearmanCorr, 3), sep=""),  color="black", size = 10) +
                 coord_cartesian(xlim=c(0, 16), ylim=c(0, 16)) +
                 scale_colour_manual(values=color_vec)
     print(g)
     dev.off()

     #print(subset(annotated_merged_abundance,gene_ID == 20534))
     #print(head(annotated_merged_abundance$data2.TPM))

}

load_packages <- function() {
    suppressPackageStartupMessages(library("readr"))
    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("plyr"))
    suppressPackageStartupMessages(library("optparse"))

    return
}

parse_options <- function() {

    option_list <- list(
        make_option(c("--f"), action = "store", dest = "infile",
                    default = NULL, help = "TALON abundance output file"),
        make_option(c("--w"), action = "store", dest = "whitelist",
                    default = NULL, help = "whitelist csv file"),
        make_option(c("--color"), action = "store", dest = "color_scheme",
                    default = NULL, help = "blue, red, or green"),
        make_option("--d1", action = "store", dest = "d1",
                    default = NULL, help = "First dataset name to use in comparison"),
        make_option("--d2", action = "store", dest = "d2",
                    default = NULL, help = "Second dataset name to use in comparison"),
        make_option(c("-o","--outdir"), action = "store", dest = "outdir",
                    default = NULL, help = "Output directory for plots and outfiles")
        )

    opt <- parse_args(OptionParser(option_list=option_list))
    opt$d1 <- as.character(opt$d1)
    opt$d2 <- as.character(opt$d2)
    return(opt)
}

main()
