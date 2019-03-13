
filter_kallisto_illumina_genes <- function(kallisto_file) {
    # This function takes a Kallisto abundance file and filters the genes
    # based on criteria designed to make the gene set comparable to what can
    # be detected using PacBio

    gencode.quantitation <- as.data.frame(read_delim(kallisto_file, "\t", escape_double = FALSE,
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))

    # Split GENCODE transcript multi-id by '|'
    extraCols =str_split_fixed(gencode.quantitation$target_id, "\\|",9)[,c(5,6,8)]
    colnames(extraCols) = c("transcript", "gene", "class")
    gencode.quantitation = cbind(extraCols, gencode.quantitation)

    # Remove transcripts that are < 300 bp in length because PacBio chucks anything that size, and
    # keep only transcripts that have polyA tails
    filter_set = c("protein_coding", "lincRNA", "processed_transcript", "macro_lncRNA")
    gencode_quant_min300 = subset(gencode.quantitation, eff_length >= 300 & class %in% filter_set)

    # Aggregate by gene
    gene_gencode_quant_min300 <- aggregate(gencode_quant_min300$tpm, by=list(gencode_quant_min300$gene), FUN=sum)
    colnames(gene_gencode_quant_min300) <- c("gene", "tpm")

    # Constraints: > 300 bp, polyA tail, TPM > 1
    filtered_genes <- gene_gencode_quant_min300[gene_gencode_quant_min300$tpm > 1,]

    # Remove genes that are on the mitochondrial blacklist
    mitochondrial_blacklist <- c("MT-TF", "MT-RNR1", "MT-TV", "MT-RNR2", "MT-TL1", "MT-ND1", "MT-TI", "MT-TQ", "MT-TM", "MT-ND2", "MT-TW", "MT-TA", "MT-TN", "MT-TC", "MT-TY", "MT-CO1", "MT-TS1", "MT-TD", "MT-CO2", "MT-TK", "MT-ATP8", "MT-ATP6", "MT-CO3", "MT-TG", "MT-ND3", "MT-TR", "MT-ND4L", "MT-ND4", "MT-TH", "MT-TS2", "MT-TL2", "MT-ND5", "MT-ND6", "MT-CYB")
    final_filtered_genes <- subset(filtered_genes, !(gene %in% mitochondrial_blacklist)) 

    return(final_filtered_genes)
}
