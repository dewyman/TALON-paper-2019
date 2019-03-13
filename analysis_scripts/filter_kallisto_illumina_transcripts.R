library(tidyverse)

filter_kallisto_illumina_transcripts <- function(kallisto_file) {
    # This function takes a Kallisto abundance file and filters the transcripts
    # based on criteria designed to make the gene set comparable to what can
    # be detected using PacBio

    gencode.quantitation <- as.data.frame(read_delim(kallisto_file, "\t", escape_double = FALSE,
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))

    # Split GENCODE transcript multi-id by '|'
    extraCols <- str_split_fixed(gencode.quantitation$target_id, "\\|",9)[,c(5,6,8)]
    colnames(extraCols) <- c("transcript", "gene", "class")
    gencode.quantitation <- cbind(extraCols, gencode.quantitation)

    # Remove transcripts that are < 300 bp in length because PacBio chucks anything that size, and
    # keep only transcripts that have polyA tails. Also require TPM > 1
    filter_set <- c("protein_coding", "lincRNA", "processed_transcript", "macro_lncRNA")
    filtered_transcripts <- subset(gencode.quantitation, length >= 300 & class %in% filter_set & tpm > 1)


    mitochondrial_blacklist <- c("MT-TF", "MT-RNR1", "MT-TV", "MT-RNR2", "MT-TL1", "MT-ND1", "MT-TI", "MT-TQ", "MT-TM", "MT-ND2", "MT-TW", "MT-TA", "MT-TN", "MT-TC", "MT-TY", "MT-CO1", "MT-TS1", "MT-TD", "MT-CO2", "MT-TK", "MT-ATP8", "MT-ATP6", "MT-CO3", "MT-TG", "MT-ND3", "MT-TR", "MT-ND4L", "MT-ND4", "MT-TH", "MT-TS2", "MT-TL2", "MT-ND5", "MT-ND6", "MT-CYB")
    final_filtered_transcripts <- subset(filtered_transcripts, !(gene %in% mitochondrial_blacklist))

    # Normalize back to 1 million transcripts
    total_transcripts <- sum(final_filtered_transcripts$tpm)
    TPM_scaling <- 1000000/total_transcripts
    final_filtered_transcripts$tpm <- final_filtered_transcripts$tpm*TPM_scaling

    return(final_filtered_transcripts)
}
