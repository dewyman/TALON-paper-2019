main <-function() {

    load_packages()
    opt <- parse_options()

    # Read abundance table
    abundance_table <- as.data.frame(read_delim(opt$infile, delim = "\t",
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))

    # Get dataset names and check that they exist in the file
    d1 <- opt$d1
    d2 <- opt$d2
    d1_type <- opt$d1_type
    d2_type <- opt$d2_type
    if (d1 %in% colnames(abundance_table) == F |
        d2 %in% colnames(abundance_table) == F) {
        print("One of the provided dataset names is not in the abundance file provided. Exiting...")
        quit()
    }

    filtered <- filter_transcripts_on_options(abundance_table, opt)
    abundance_table <- filtered$abundance_table
    color_vec <- filtered$color_vec
    t_levels <- filtered$t_levels

    # Compute TPMs for each dataset
    d1_abundance <- abundance_table[,c("transcript_ID", d1)]
    d2_abundance <- abundance_table[,c("transcript_ID", d2)]
    colnames(d1_abundance) <- c("transcript_ID", "count")
    colnames(d2_abundance) <- c("transcript_ID", "count")

    d1_abundance$data1.TPM <- (d1_abundance$count)*1000000/sum(d1_abundance$count)
    d2_abundance$data2.TPM <- (d2_abundance$count)*1000000/sum(d2_abundance$count)

    merged_abundances <- merge(d1_abundance, d2_abundance, by = "transcript_ID", all.x = T, all.y = T)
    #colnames(merged_abundances) <- c("transcript_ID", "data1.TPM", "data2.TPM")
    merged_abundances[is.na(merged_abundances)] <- 0   

    # Merge in transcript novelty information
    transcript_novelty <- abundance_table[, c("transcript_ID", 
                                              "transcript_novelty")] 
    merged_abundances <- merge(merged_abundances, transcript_novelty, by = "transcript_ID", 
                               all.x = T, all.y = F)
    merged_abundances$novelty <- merged_abundances$transcript_novelty
    merged_abundances$transcript_novelty <- NULL

    merged_abundances$novelty <- factor(merged_abundances$novelty, levels = t_levels)

    # Plot expression scatterplots
    expression_by_status(merged_abundances, d1, d2, opt$outdir, color_vec, opt$celltype, opt$lsr, opt$corr_labs, opt$regression_line, d1_type, d2_type)
}

filter_transcripts_on_options <- function(abundance_table, opt) {

    d1 <- opt$d1
    d2 <- opt$d2

    # Get colors
    green <- "#009E73"
    blue <- "#0072B2"
    orange <- "#D55E00"
    gold <- "#E69F00"
    black <- "#000000"
    pink <- "#CC79A7"
    yellow <- "#F0E442"

    color_vec <- c(green)
    t_levels <- c("Known")

    # Restrict to transcripts that were observed in at least one of the datasets
    abundance_table <- abundance_table[abundance_table[,d1] +
                                       abundance_table[,d2] > 0, ]

    # TODO: Use the whitelist to filter transcripts if one is provided

    # Remove transcripts depending on the input options
    if (opt$ISM == F) {
        message("Removing ISM transcripts...")
        abundance_table <- subset(abundance_table,
                                  transcript_novelty != "ISM")
    } else {
        color_vec <- c(color_vec, blue)
        t_levels <- c(t_levels, "ISM")
    }
    if (opt$NIC == F) {
        message("Removing NIC transcripts...")
        abundance_table <- subset(abundance_table,
                                  transcript_novelty != "NIC")
    } else {
        color_vec <- c(color_vec, orange)
        t_levels <- c(t_levels, "NIC")
    }
    if (opt$NNC == F) {
        message("Removing NNC transcripts...")
        abundance_table <- subset(abundance_table,
                                  transcript_novelty != "NNC")
    } else {
        color_vec <- c(color_vec, gold)
        t_levels <- c(t_levels, "NNC")
    }
    if (opt$antisense == F) {
        message("Removing antisense transcripts...")
        abundance_table <- subset(abundance_table,
                                  transcript_novelty != "Antisense")
    }
    else {
        color_vec <- c(color_vec, black)
        t_levels <- c(t_levels, "Antisense")
    }
    if (opt$intergenic == F) {
        message("Removing intergenic transcripts")
        abundance_table <- subset(abundance_table,
                                  transcript_novelty != "Intergenic")
    }
    else {
        color_vec <- c(color_vec, pink)
        t_levels <- c(t_levels, "Intergenic")
    }
    if (opt$genomic == F) {
        message("Removing genomic transcripts")
        abundance_table <- subset(abundance_table,
                                  transcript_novelty != "Genomic")
    }
    else {
        color_vec <- c(color_vec, yellow)
        t_levels <- c(t_levels, "Genomic")
    }

   
    filtered <- list(abundance_table, color_vec, t_levels)
    names(filtered) <- c("abundance_table", "color_vec", "t_levels")
    return(filtered)
}

expression_by_status <- function(merged_abundances, d1, d2, outdir, color_vec, celltype, lsr, corr_labs, regression_line, d1_type, d2_type) {

    # Take log2(TPM + 0.1)
    merged_abundances$data1.log_TPM = log(merged_abundances$data1.TPM + 0.1, base=10)
    merged_abundances$data2.log_TPM = log(merged_abundances$data2.TPM + 0.1, base=10)
    t_levels <- levels(merged_abundances$novelty)
    
    # Plot log2(TPM + 0.1) for each dataset on a scatterplot. Color points according to known/novel status
    pearsonCorr = cor.test(~data1.TPM + data2.TPM, 
                           data=merged_abundances, method = "pearson", 
                           continuity = FALSE, conf.level = 0.95)$estimate
    spearmanCorr = cor.test(~data1.TPM + data2.TPM, 
                            data=merged_abundances, method = "spearman", 
                            continuity = FALSE, 
                            exact = FALSE, conf.level = 0.95)$estimate

    # Least-Square Regression Line
    mod<-lm(data2.TPM~data1.TPM, data=merged_abundances)
 
    nov_types <- paste(t_levels, collapse = '-')
    joined_names <- paste(outdir, "/", d1, "-", d2, "_", nov_types, sep = "")
    
    fname <- paste(joined_names, "transcript", "correlationPlot.png", sep="_")
    corr_fname <- paste(joined_names, "transcript", "correlations.txt", sep="_")

    xlabel <- paste("TPM+0.1 in ", celltype, " ", d1_type, sep="")
    ylabel <- paste("TPM+0.1 in ", celltype, " ", d2_type, sep="")
    corr_label <- paste("Pearson r: ",
                            round(pearsonCorr, 2), "\nSpearman rho: ",
                            round(spearmanCorr, 2), "\nLSR slope: ",
                            round(mod$coefficients[2], 2), sep="")
    if (lsr == T) {
        plot_label <- paste("Pearson r: ",
                            round(pearsonCorr, 2), "\nSpearman rho: ",
                            round(spearmanCorr, 2), "\nLSR slope: ",
                            round(mod$coefficients[2], 2), sep="")
    } else {
        plot_label <- paste("Pearson r: ",
                    round(pearsonCorr, 2), "\nSpearman rho: ",
                    round(spearmanCorr, 2), "\nLSR slope: ", sep="")
    }

    # write correlation numbers to outfile
    write(corr_label, corr_fname)

    png(filename = fname,
        width = 2700, height = 2500, units = "px",
        bg = "white",  res = 300)
    scatterplot <- ggplot(merged_abundances, aes(x = data1.TPM, y = data2.TPM, color = novelty)) +
        geom_jitter(alpha = 0.5) + theme_bw() +
        xlab(xlabel)  + ylab(ylabel) + theme(text= element_text(size=24)) +
        theme(axis.text.x = element_text(color = "black", size=24),
              axis.text.y = element_text(color = "black", size=24)) +
        scale_x_continuous(trans=log10_trans(), limits=c(1,32768))+
        scale_y_continuous(trans=log10_trans(), limits=c(1,32768))+
        scale_colour_manual("Transcript status", values=color_vec) +
        theme(legend.position=c(0.73,0.2),
              legend.title = element_text(colour = 'black', size = 21),
              legend.background = element_rect(fill="white", color = "black"),
              legend.key = element_rect(fill="transparent"),
              legend.text = element_text(colour = 'black', size = 20))+
        guides(colour = guide_legend(override.aes = list(alpha=1, size=3)))

    # add regression line
    if (regression_line) {
        scatterplot <- scatterplot+
              geom_abline(slope = mod$coefficients[2], intercept = mod$coefficients[1],
                   color = "gray", lwd=1, lty=2)
    }
    if (corr_labs) {
        scatterplot <- scatterplot+
              annotate("text", x = 5, y = 14, label = plot_label,
                   color="black", size = 10) 
    }


     # Find max density y value across both datasets
     vars <- unique(merged_abundances$novelty)
     xd_max <- max(sapply(vars, compute_max_density_for_var, merged_abundances, "data2.log_TPM"))
     yd_max <- max(sapply(vars, compute_max_density_for_var, merged_abundances, "data1.log_TPM"))
     plot_max <- round(max(c(xd_max, yd_max))*1.001, 2)

     # density x lims 
     density_xmin = log(0.1, base=10)
     density_xmax = log(32768, base=10)

     # Marginal density plot of x (top panel)
     xdensity <- ggplot(merged_abundances, aes(data1.log_TPM, fill=novelty, color=novelty)) + 
                        geom_density(alpha=.5) + 
                        scale_fill_manual(values = color_vec) + 
                        scale_color_manual(values = color_vec) +
                        scale_x_continuous(breaks = seq(density_xmin, density_xmax, by = density_xmax), expand = c(0,0)) +
                        scale_y_continuous(breaks = seq(0, plot_max, by = plot_max )) +
                        theme(legend.position = "none",
                              axis.title.x=element_blank(),
                              axis.text.x=element_blank(),
                              axis.ticks.x=element_blank(),
                              axis.text.y=element_text(color = "black", size=14),
                              axis.title.y=element_text(color = "black", size=20),
                              plot.margin = margin(0.75, 0, 0, 1.65, "cm"))

    # Marginal density plot of y (right panel)
    ydensity <- ggplot(merged_abundances, aes(data2.log_TPM, fill=novelty, color=novelty)) + 
                       geom_density(alpha=.5) + 
                       scale_fill_manual(values = color_vec) +
                       scale_color_manual(values = color_vec) +
                       scale_x_continuous(breaks = seq(density_xmin, density_xmax, by = density_xmax), expand = c(0,0)) +
                       scale_y_continuous(breaks = seq(0, plot_max, by = plot_max )) +
                       theme(legend.position = "none",
                             axis.title.y=element_blank(),
                             axis.text.y=element_blank(),
                             axis.ticks.y=element_blank(),
                             axis.text.x=element_text(color = "black", size=14),
                             axis.title.x=element_text(color = "black", size=20),
                             plot.margin = margin(0, 0.75, 0.3, 0, "cm")) +
                             coord_flip(ylim = c(0, plot_max))
                       

    # Blank placeholder plot
    blankPlot <- ggplot()+geom_blank(aes(1,1))+
                 theme(
                      plot.background = element_blank(), 
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), 
                      panel.border = element_blank(),
                      panel.background = element_blank(),
                      axis.title.x = element_blank(),
                      axis.title.y = element_blank(),
                      axis.text.x = element_blank(), 
                      axis.text.y = element_blank(),
                      axis.ticks = element_blank(),
                      axis.line = element_blank())

     g = grid.arrange(xdensity, blankPlot, scatterplot, ydensity, 
                      ncol=2, nrow=2, widths=c(4, 0.9), heights=c(0.9, 4))
     
     print(g)
     dev.off()

}

compute_max_density_for_var <- function(var, data, dataset) {
    # Compute the max density y value for a particular gene type
    data <- subset(data, data$novelty == var)[,dataset]
    return(max(density(data)$y))
}

load_packages <- function() {
    suppressPackageStartupMessages(library("tidyverse"))
    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("gridExtra"))
    suppressPackageStartupMessages(library("cowplot"))
    suppressPackageStartupMessages(library("scales"))


    return
}

parse_options <- function() {

    option_list <- list(
        make_option(c("--f"), action = "store", dest = "infile",
                    default = NULL, help = "FILTERED TALON abundance output file"),
        make_option(c("--color"), action = "store", dest = "color_scheme",
                    default = NULL, help = "blue, red, or green"),
       make_option("--d1", action = "store", dest = "d1",
                    default = NULL, help = "First dataset name to use in comparison"),
        make_option("--d1_type", action="store", dest="d1_type",
                    help="datatype of dataset 1 ie 'Rep1 PacBio' or 'Rep2 ONT'"),
        make_option("--d2", action = "store", dest = "d2",
                    default = NULL, help = "Second dataset name to use in comparison"),
        make_option("--d2_type", action="store", dest="d2_type",
                    help="datatype of dataset 2 ie 'Rep1 PacBio' or 'Rep2 ONT'"),
        make_option("--celltype", action = "store", dest = "celltype",
                    default = NULL, help = "Celltype to use in plot labels"),
        make_option(c("--lsr"), action="store_true", dest="lsr",
              help="Include this option if you want the LSR label on", default = F),
        make_option("--w", action = "store", dest = "whitelist",
                    default = NULL, help = "Optional: a whitelist for filtering transcripts"),
        make_option(c("--ISM"), action="store_true", dest="ISM",
              help="Set this option to include ISM transcripts", default = F),
        make_option(c("--NIC"), action="store_true", dest="NIC",
              help="Set this option to include NIC transcripts", default = F),
        make_option(c("--NNC"), action="store_true", dest="NNC",
              help="Set this option to include ISM transcripts", default = F),
        make_option(c("--antisense"), action="store_true", dest="antisense",
              help="Set this option to include antisense transcripts", default = F),
        make_option(c("--intergenic"), action="store_true", dest="intergenic",
              help="Set this option to include intergenic transcripts", default = F),
        make_option(c("--genomic"), action="store_true", dest="genomic",
              help="Set this option to include genomic transcripts", default = F),
        make_option(c("-o","--outdir"), action = "store", dest = "outdir",
                    default = NULL, help = "Output directory for plots and outfiles"),
        make_option("--correlations", action="store_true", dest = "corr_labs",
              help="Add correlation labels to plot", default = F),
        make_option("--regression_line", action="store_true", dest = "regression_line",
              help="Add regression line to plot", default = F)
        )

    opt <- parse_args(OptionParser(option_list=option_list))
    opt$d1 <- as.character(opt$d1)
    opt$d2 <- as.character(opt$d2)
    return(opt)
}

main()
