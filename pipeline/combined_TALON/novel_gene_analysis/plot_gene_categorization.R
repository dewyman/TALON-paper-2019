main <-function() {

    set.seed(100)
    load_packages()
    opt <- parse_options()

    # Read data
    data <- as.data.frame(read_delim(opt$infile, "\t",
                                         escape_double = FALSE,
                                         trim_ws = TRUE, col_names = TRUE))
    
    data$tot_length <- data$end - data$start
    data$sense_coverage <- data$known_gene_overlap/data$tot_length
    data$antisense_coverage <- data$known_gene_as_overlap/data$tot_length

    # Assign each gene to the category that is the best fit
    data$Category <- ifelse(data$chromosome == "chrEBV", "On EBV chromosome", 
                     ifelse(data$sense_coverage > data$antisense_coverage, "Overlaps known (same strand)", 
                     ifelse(data$antisense_coverage > 0 & data$n_exons == "monoexonic", "Overlaps known (antisense), is monoexonic",
                     ifelse(data$antisense_coverage > 0 & data$n_exons == "multiexonic", "Overlaps known (antisense), is multiexonic",
                     ifelse(abs(data$closest_gene_dist) < 1000 & data$closest_gene_dist < 0,
                            "< 1000 bp downstream of known (same strand)",
                     ifelse(abs(data$closest_gene_dist) < 1000 & data$closest_gene_dist > 0,
                            "< 1000 bp upstream of known (same strand)", 
                     ifelse(abs(data$closest_gene_as_dist) < 1000 & data$closest_gene_as_dist < 0,
                            "< 1000 bp downstream of known (antisense)",
                     ifelse(abs(data$closest_gene_as_dist) < 1000 & data$closest_gene_as_dist > 0,
                            "< 1000 bp upstream of known (antisense)",            
                     ifelse(data$alu_overlap > 150 | data$line1_overlap > 1000, "Repeat", "Intergenic")))))))))

    # Group similar categories
    data[data$Category == "< 1000 bp downstream of known (same strand)" |
         data$Category == "< 1000 bp downstream of known (antisense)", "Category"] <- "< 1000 bp downstream of known gene (any strand)"
    data[data$Category == "< 1000 bp upstream of known (same strand)" |
         data$Category == "< 1000 bp upstream of known (antisense)", "Category"] <- "< 1000 bp upstream of known gene (any strand)"

    plot_summary_pie_chart(data, opt$outdir)
    write.table(data, file = paste(opt$outdir, "/novel_gene_assignments.tsv", sep=""), 
                sep ="\t", quote=F, col.names=T, row.names=F)
}

plot_summary_pie_chart <- function(data, outdir) {

    fname <- paste(outdir, "/novel_genes_by_category.png", sep="")
    #print(subset(data, Category == "Intergenic")[,1:4])
    data <- as.data.frame(table(data$Category))
    colnames(data) <- c("Category", "Count")

    data %>%
    arrange(desc(Count)) %>%
    mutate(prop = percent(Count / sum(Count))) -> data 

    png(filename = fname,
        width = 3000, height = 2500, units = "px",
        bg = "white",  res = 300)

    g = ggplot(data, aes(x = "", y = Count, fill = fct_inorder(Category))) +
       geom_bar(width = 1, stat = "identity") + theme_void() +
       coord_polar("y", start = 0) +
       geom_label_repel(aes(label = prop), size=8, show.legend = F, segment.color = 'transparent') +
       guides(fill = guide_legend(title = "Category")) +
       scale_fill_brewer(palette = "Paired", type = "div") +
       theme(legend.text=element_text(color = "black", size=16)) 

    print(g)
    dev.off()

    # Plot a second time without the labels
    fname_2 <- paste(outdir, "/novel_genes_by_category_nolabels.png", sep="")
    png(filename = fname_2,
        width = 3000, height = 2500, units = "px",
        bg = "white",  res = 300)

    g = ggplot(data, aes(x = "", y = Count, fill = fct_inorder(Category))) +
       geom_bar(width = 1, stat = "identity") + theme_void() +
       coord_polar("y", start = 0) +
       guides(fill = guide_legend(title = "Category")) +
       scale_fill_brewer(palette = "Paired", type = "div") +
       theme(legend.text=element_text(color = "black", size=16))

    print(g)
    dev.off()

}

load_packages <- function() {
    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("readr"))
    suppressPackageStartupMessages(library("reshape"))
    suppressPackageStartupMessages(library("ggrepel"))
    suppressPackageStartupMessages(library("forcats"))
    suppressPackageStartupMessages(library("scales"))
    suppressPackageStartupMessages(library("RColorBrewer"))
    return
}

parse_options <- function() {

    option_list <- list(
        make_option(c("--f"), action = "store", dest = "infile",
                    default = NULL, help = "Tab-delimited file of information about the novel genes"),
        make_option(c("-o","--outdir"), action = "store", dest = "outdir",
                    default = NULL, help = "Output directory for plots and outfiles")
        )

    opt <- parse_args(OptionParser(option_list=option_list))
    return(opt)
}


main()
