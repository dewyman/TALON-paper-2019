
library(IsoformSwitchAnalyzeR)
library(readr)

abundancef <- "/pub/dwyman/TALON-paper-2019/pipeline/combined_TALON/tier1_talon_abundance_filtered.tsv"
designf <- "/pub/dwyman/TALON-paper-2019/pipeline/combined_TALON/isr_design_matrix.csv"
gtf <- "/pub/dwyman/TALON-paper-2019/pipeline/combined_TALON/tier1_filtered_talon.gtf"
fasta <- "/pub/dwyman/TALON-paper-2019/pipeline/combined_TALON/tier1_filtered_talon.fasta"

# Read the TALON abundance file, then reformat as a Kallisto-style file
abundance_table <- as.data.frame(read_delim(abundancef, delim = "\t",
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))



# Read the design matrix
myDesign <- as.data.frame(read_delim(designf, delim = ",",
                                  col_names = TRUE, trim_ws = TRUE, na = "NA"))



# Create switchAnalyzeRlist
aSwitchList <- importRdata(
    isoformCountMatrix   = abundance_table$counts,
    isoformRepExpression = salmonQuant$abundance,
    designMatrix         = myDesign,
    isoformExonAnnoation = system.file(gtf, package="IsoformSwitchAnalyzeR"),
    isoformNtFasta       = system.file(fasta, package="IsoformSwitchAnalyzeR"),
    showProgress = TRUE
)
