
main <-function() {

    set.seed(100)
    load_packages()
    opt <- parse_options()

    # Colorways
    if (opt$color_scheme == "blue") {
        gene_color <- "navy"
        transcript_color <- "skyblue"
    } else if (opt$color_scheme == "red") {
        gene_color <- "red"
        transcript_color <- "orange"
    } else if (opt$color_scheme == "green") {
        gene_color <- "springgreen4"
        transcript_color <- "olivedrab3"
    }

    # Process the read count file. There should be a dataset name column,
    # a platform column, and a read count column
    read_counts_by_dataset <- as.data.frame(read_delim(opt$read_ct_file, ",", 
                                         escape_double = FALSE, 
                                         trim_ws = TRUE, col_names = TRUE))       

    # Dataset names
    datasets <- read_counts_by_dataset$dataset

    # Read the abundance file and remove genomic transcripts
    abundance_table <- as.data.frame(read_delim(opt$infile, delim = "\t",
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))
    abundance_table <- subset(abundance_table, transcript_novelty != "Genomic")
    abundance_table <- subset(abundance_table, gene_novelty == "Known")

    # Restrict abundance file to gene info and only the datasets we're interested in
    abundance_table <- abundance_table[, c("gene_ID", "annot_gene_name", datasets)]

    # Aggregate genes and remove rows where the gene is not detected at all
    pb_gene_abundance <- ddply(abundance_table, c("annot_gene_name"), 
                         function(x) colSums(x[datasets]))
    pb_gene_abundance <- pb_gene_abundance[apply(pb_gene_abundance[datasets],1,function(z) any(z!=0)),]

    # Get Illumina genes in both reps, and create a unified list
    illumina_1 <- as.character(filter_kallisto_illumina_genes(opt$illumina_kallisto_1)$gene)#[,1]
    illumina_2 <- as.character(filter_kallisto_illumina_genes(opt$illumina_kallisto_2)$gene)#[,1]
    illumina_genes <- intersect(illumina_1,illumina_2)
    n_illumina_genes <- length(illumina_genes)

    # Remove rows from abundance file that are not in the Illumina file
    pb_gene_abundance <- subset(pb_gene_abundance, annot_gene_name %in% illumina_genes)

    # Now, find the genes detected in each dataset
    genes_detected_by_dataset <- lapply(X=read_counts_by_dataset$dataset, FUN=get_detected_genes_for_dataset, pb_gene_abundance)
    names(genes_detected_by_dataset) <- read_counts_by_dataset$dataset
    
    # Get all possible dataset permutations, with attached gene/transcript lists
    dataset_permutations <- permn(read_counts_by_dataset$dataset)

    # There might be a lot of them, so pick 10 at random to plot
    dataset_permutations <- sample(dataset_permutations, 10, replace=FALSE)
   
    gene_curves <- lapply(X = dataset_permutations, FUN = get_curve_points, "genes", genes_detected_by_dataset, read_counts_by_dataset)

    all_gene_curves <- do.call("rbind", gene_curves)

    # Plot the curves
    plot_discovery_curve(all_gene_curves, "genes", gene_color, n_illumina_genes, 15000, opt$outdir) 

}

plot_discovery_curve <- function(all_curves, cat_type, color, illumina_line, ymax, outdir) {

    # Plot the curves
    fname <- paste(outdir, "/discovery_curves_", cat_type, "_knownOnly.png", sep="")
    xlabel <- "Number of raw reads (millions)"
    ylabel <- paste("Number of ", cat_type, "  detected so far", sep="")
    all_curves$ID <- as.factor(all_curves$ID)
    all_curves$read_count <- round(all_curves$read_count/1000000, 2)
    if (cat_type == "genes") {
        Illumina_label <- "Number of genes detected in Illumina"
    } else if (cat_type == "transcripts") { 
        Illumina_label <- "Number of transcripts* detected in Illumina (TPM > 1)"
    }

    all_curves$percent_detected <- all_curves$items/illumina_line
    print(all_curves)

    png(filename = fname,
        width = 2500, height = 2000, units = "px",
        bg = "white",  res = 300)
    g = ggplot(all_curves, aes(x=read_count, y=items, colour=all_curves$ID)) + geom_point() + geom_line() +
               scale_colour_manual(values = c( rep(color, length(levels(all_curves$ID))))) +
               guides(colour=FALSE) + theme_bw() +
               theme(axis.text.x = element_text(color = "black", size=20),
                   axis.text.y = element_text(color = "black", size=20),
                   axis.title.x = element_text(color="black", size=20),
                   axis.title.y = element_text(color="black", size=20)) +
               xlab(xlabel) + ylab(ylabel) +
               coord_cartesian(ylim = c(0, ymax)) +
               scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7,8)) +
               geom_hline(yintercept=illumina_line, linetype="dashed", color = "black") +
               annotate("text", label = Illumina_label, x = 3, y = illumina_line + 700, color = "black", size = 7) 

    print(g)
    dev.off()
    print(paste("Total ", cat_type, ": ", max(all_curves$items), sep=""))

}


load_packages <- function() {
    suppressPackageStartupMessages(library("DBI"))
    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("plyr"))
    suppressPackageStartupMessages(library("stringr"))
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("RSQLite"))
    suppressPackageStartupMessages(library("VennDiagram"))
    suppressPackageStartupMessages(library("readr"))
    suppressPackageStartupMessages(library("ggrepel"))
    suppressPackageStartupMessages(library("reshape"))
    suppressPackageStartupMessages(library("combinat"))
    

    # Load my custom functions
    source("/pub/dwyman/TALON-paper-2019/analysis_scripts/filter_kallisto_illumina_genes.R")
    return
}


parse_options <- function() {

    option_list <- list(
        make_option(c("--f"), action = "store", dest = "infile",
                    default = NULL, help = "TALON abundance output file"),
        make_option(c("--color"), action = "store", dest = "color_scheme",
                    default = "blue", help = "Colorway to use for plot (blue, red, or green)"),
        make_option(c("--rc"), action = "store", dest = "read_ct_file",
                    default = NULL, help = "CSV file with name of dataset and number of raw reads"),
        make_option(c("--ik1"), action = "store", dest = "illumina_kallisto_1",
                    default = NULL, help = "Rep1 Illumina Kallisto file."),
        make_option(c("--ik2"), action = "store", dest = "illumina_kallisto_2",
                    default = NULL, help = "Rep2 Illumina Kallisto file."),
        make_option(c("-o","--outdir"), action = "store", dest = "outdir",
                    default = NULL, help = "Output directory for plots and outfiles")
        )

    opt <- parse_args(OptionParser(option_list=option_list))
    return(opt)
}

get_curve_points <- function(dataset_list, cat_type, items_by_dataset, read_counts_by_dataset) {

    n_datasets <- length(dataset_list)

    cat_types <- rep(cat_type, n_datasets + 1)
    perm_id_vector <- rep(paste( cat_type, "_", paste(dataset_list,collapse="-"), sep=""), n_datasets + 1)
    cumulative_read_counts <- c(0, cumsum(read_counts_by_dataset[match(dataset_list, read_counts_by_dataset$dataset),"read_count"]))
 
    # Get the number of genes/transcripts seen with each new dataset in the list 
    ordered_items <- items_by_dataset[dataset_list]

    curr_items_seen = ordered_items[[1]]
    num_items_seen = c(0, length(curr_items_seen))

    for (i in c(2:length(ordered_items))) {
        curr = ordered_items[[i]]
        prev_items_seen = curr_items_seen
        curr_items_seen = union(prev_items_seen, curr)
        num_items_seen = c(num_items_seen, length(curr_items_seen))
    }
    result <- data.frame("ID" = as.factor(perm_id_vector), "read_count" = as.numeric(cumulative_read_counts), "items" = as.numeric(num_items_seen), "cat" = as.factor(cat_types))
    return(result)
}

get_n_genes_for_dataset <- function(dataset, pb_gene_abundance) {
    # Get genes and return count
    genes <- get_detected_genes_for_dataset(dataset, whitelisted_genes, database)
    return(length(genes))
}

get_detected_genes_for_dataset <- function(dataset, pb_gene_abundance) {

    genes_detected <- pb_gene_abundance[ pb_gene_abundance[,dataset] > 0, "annot_gene_name"]
    return(genes_detected)
}




main()
