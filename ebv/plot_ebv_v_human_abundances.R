suppressMessages(library(tidyverse))
suppressMessages(library(optparse))
suppressMessages(library(stringr))

my_theme <- theme(axis.line = element_line(colour = "black"),
                    panel.background = element_blank(),
                    legend.text = element_text(color="black", size = rel(1)),
                    legend.title = element_text(color="black", size=rel(1.25)),
                    legend.position=c(0.2,0.8),
                    legend.background = element_rect(fill="white", color = "black"),
                    legend.key = element_rect(fill="transparent"),
                    plot.title = element_text(hjust = 0.5),
                    axis.line.x = element_line(color="black", size = 0.5),
                    axis.line.y = element_line(color="black", size = 0.5),
                    axis.text.x = element_text(color="black", size = rel(1.5)),
                    # axis.text.x = element_blank(),
                    axis.text.y = element_text(color="black", size = rel(2)),
                    axis.title.x = element_blank(),
                    axis.title.y = element_text(color="black", size=rel(1.5)))

# get csv files to read in data from
option_list <- list(
    make_option(c('--gene_csv'), action='store', dest='gene_csv', 
                help='CSV with gene TPMs'),
    make_option(c('--transcript_csv'), action='store', dest='transcript_csv',
                help='CSV with transcript TPMs'),
    make_option(c('--datasets'), action='store', dest='datasets'))
opt <- parse_args(OptionParser(option_list=option_list))
gene_csv = opt$gene_csv
transcript_csv = opt$transcript_csv
datasets = unlist(strsplit(opt$datasets, ","))
odir = dirname(gene_csv)

# print(gene_csv)
# print(transcript_csv)
# print(datasets)
# print(datasets[1])
# print(datasets[2])
# print(odir)

# read in infiles as dataframes
g_df <- suppressMessages(read_csv(gene_csv))
t_df <- suppressMessages(read_csv(transcript_csv))

print(distinct(g_df, gene_novelty))

# color palette 
colors <- c("Known" = "#009E73","ISM" = "#0072B2", "NIC" = "#D55E00",
               "NNC" = "#E69F00", "Antisense" = "#000000",
               "Intergenic" = "#CC79A7")

# generate figures for each dataset 
for(d in datasets){
  
  # get TPM column names 
  gt_tpm = paste('log_TPM_',d, sep="")
  print(gt_tpm)

  # get height column names
  gt_height = paste(d,'height',sep="_")
  print(gt_height)

  # get transcript/gene count column names 
  gt_count = paste('n',d,sep="_")
  print(gt_count)

  # get number of human and ebv transcripts from dataset
  t_ebv = tally(filter(t_df, UQ(as.name(d)) != 0 & ebv == 'EBV')) # ebv transcripts
  t_human = tally(filter(t_df, UQ(as.name(d)) != 0 & ebv == 'Human')) # human transcripts
  g_ebv = tally(filter(g_df, UQ(as.name(d)) != 0 & ebv == 'EBV')) # ebv genes
  g_human = tally(filter(g_df, UQ(as.name(d)) != 0 & ebv == 'Human')) # human genes

  # create output file names
  t_ofile = paste(odir, '/', d, '_transcripts_ebv_human.png', sep='')
  g_ofile = paste(odir, '/', d, '_genes_ebv_human.png', sep='')

  # refactor novelty levels to make more sense
  t_df$transcript_novelty <- factor(t_df$transcript_novelty,
                               levels = c("Known", "ISM", 
                                          "NIC", "NNC", 
                                          "Antisense", "Intergenic"))    
  g_df$gene_novelty <- factor(g_df$gene_novelty,
                          levels = c("Known", "Antisense", "Intergenic"))


  # transcript plotting
  png(filename = t_ofile,
      width = 2500, height = 2500, units = "px",
      bg = "white",  res = 300)
  p <- ggplot(data=t_df, aes(x=ebv,
        y=UQ(as.name(gt_tpm))))+
     geom_violin(fill='gray90')+
     theme_bw(base_family = "Helvetica", base_size = 18)+
     my_theme+
     geom_jitter(data=t_df, aes(color=transcript_novelty,
        size=dot_size, alpha=alpha),
        shape=16,
        position=position_jitter(0.2))+
        # alpha=0.5)+
     scale_alpha_continuous(range=c(0.2,0.5), guide=FALSE)+
     scale_size_continuous(range=c(2,3), guide=FALSE)+
     guides(colour = guide_legend(override.aes = list(alpha = 1,size=3)))+
     geom_text(data=t_df, aes(x=ebv,
        y=UQ(as.name(gt_height)),
        label=paste('n=',UQ(as.name(gt_count)),sep='')),
        size = rel(7.5))+
     scale_color_manual("Isoform Type", values=colors)+
     ylab('log2(TPM+1)')
  print(p)
  dev.off()

  # gene plotting
  png(filename = g_ofile,
      width = 2500, height = 2500, units = "px",
      bg = "white",  res = 300)
  p <- ggplot(data=g_df, aes(x=ebv,
        y=UQ(as.name(gt_tpm))))+
     geom_violin(fill='gray90')+
     theme_bw(base_family = "Helvetica", base_size = 18)+
     my_theme+
     geom_jitter(data=g_df, aes(color=gene_novelty,
        size=dot_size, alpha=alpha),
        shape=16,
        position=position_jitter(0.2))+
        # alpha=0.5)+
     scale_alpha_continuous(range=c(0.2,0.5), guide=FALSE)+
     scale_size_continuous(range=c(2,3), guide=FALSE)+
     guides(colour = guide_legend(override.aes = list(alpha = 1, size=3)))+
     geom_text(data=g_df, aes(x=ebv,
        y=UQ(as.name(gt_height)),
        label=paste('n=',UQ(as.name(gt_count)),sep='')),
        size = rel(7.5))+
     scale_color_manual("Gene Type",values=colors)+
     ylab('log2(TPM+1)')
  print(p)
  dev.off()
}



