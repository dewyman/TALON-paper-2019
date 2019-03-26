main <-function() {

    load_packages()
    opt <- parse_options()

    # Get colors
    if (opt$color_scheme == "red") {
        color_vec <- c("red2", "orange")
    } else if (opt$color_scheme == "blue") {
        color_vec <- c("navy", "orange")
    } else if (opt$color_scheme == "green") {
        color_vec <- c("springgreen4", "black", "orange")
    }

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

    # Restrict to genes that were observed in at least one of the datasets
    abundance_table <- abundance_table[abundance_table[,d1] +
                                       abundance_table[,d2] > 0, ]

    # Remove transcripts depending on the input options
    if (opt$antisense == F) {
        message("Removing antisense genes")
        color_vec <- c(color_vec[1], color_vec[3])
        abundance_table <- subset(abundance_table, 
                                  antisense_gene != "antisense_gene")
    }
    if (opt$intergenic == F) {
        message("Removing intergenic genes")
        color_vec <- color_vec[-length(color_vec)]
        abundance_table <- subset(abundance_table, 
                                  intergenic_gene != "intergenic_gene")
    } 

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
   
    merged_abundances <- merge(d1_abundance_agg, d2_abundance_agg, by = "gene_ID", all.x = T, all.y = T)
    colnames(merged_abundances) <- c("gene_ID", "data1.TPM", "data2.TPM")
    merged_abundances[is.na(merged_abundances)] <- 0   

    # Merge in gene novelty information
    gene_novelty <- unique(abundance_table[, c("gene_ID", "gene_status", "antisense_gene", 
                                        "intergenic_gene")])
    merged_abundances <- merge(merged_abundances, gene_novelty, by = "gene_ID", 
                               all.x = T, all.y = F)
    merged_abundances$novelty <- NA
    merged_abundances[merged_abundances$gene_status == "KNOWN", "novelty"] <- "Known"
    merged_abundances[merged_abundances$antisense_gene == "antisense_gene", "novelty"] <- "Antisense"
    merged_abundances[merged_abundances$intergenic_gene == "intergenic_gene", "novelty"] <- "Intergenic"
    merged_abundances$novelty <- factor(merged_abundances$novelty, levels = c("Known", "Antisense", "Intergenic"))

    # Plot expression scatterplots
    expression_by_status(merged_abundances, d1, d2, opt, opt$outdir, color_vec, opt$celltype)
}

expression_by_status <- function(merged_abundances, d1, d2, options, outdir, color_vec, celltype) {

    # Take log2(TPM + 1)
    merged_abundances$data1.TPM = log(merged_abundances$data1.TPM + 1, base=2)
    merged_abundances$data2.TPM = log(merged_abundances$data2.TPM + 1, base=2)
    
    # Plot log2(TPM + 1) for each dataset on a scatterplot. Color points according to known/novel status
    pearsonCorr = cor.test(~data1.TPM + data2.TPM, data=merged_abundances, method = "pearson", continuity = FALSE, conf.level = 0.95)$estimate
    spearmanCorr = cor.test(~data1.TPM + data2.TPM, data=merged_abundances, method = "spearman", continuity = FALSE, conf.level = 0.95)$estimate

    joined_names <- paste(outdir, "/", d1, "-", d2, sep = "")
    if (options$antisense == T) {
        joined_names <- paste(joined_names, "withAntisense", sep="_")
    }
    if (options$intergenic == T) {
        joined_names <- paste(joined_names, "withIntergenic", sep="_")
    }
    fname <- paste(joined_names, "gene", "correlationPlot.png", sep="_")
    xlabel <- paste("log2(TPM+1) in ", celltype, " Rep1", sep="")
    ylabel <- paste("log2(TPM+1) in ", celltype, " Rep2", sep="")

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

     g = ggMarginal(g, groupColour = TRUE, groupFill = TRUE ) + scale_y_continuous()

     print(g)
     dev.off()

}

load_packages <- function() {
    suppressPackageStartupMessages(library("tidyverse"))
    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("ggExtra"))
    suppressPackageStartupMessages(library("scales"))
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
        make_option("--celltype", action = "store", dest = "celltype",
                    default = NULL, help = "Celltype to use in plot labels"),
        make_option(c("--antisense"), action="store_true", dest="antisense",
              help="Set this option to include antisense genes", default = F),
        make_option(c("--intergenic"), action="store_true", dest="intergenic",
              help="Set this option to include intergenic genes", default = F),
        make_option(c("-o","--outdir"), action = "store", dest = "outdir",
                    default = NULL, help = "Output directory for plots and outfiles")
        )

    opt <- parse_args(OptionParser(option_list=option_list))
    opt$d1 <- as.character(opt$d1)
    opt$d2 <- as.character(opt$d2)
    return(opt)
}

main()
