PLOTPATH=~/mortazavi_lab/bin/TALON-paper-2019/analysis_scripts/
IPATH=~/mortazavi_lab/bin/TALON-paper-2019/Illumina/
PBPATH=~/mortazavi_lab/data/tier1_filtered/
ONTPATH=~/mortazavi_lab/data/ont_tier1/

# test
# set -e
# PLOTPATH=~/mortazavi_lab/bin/TALON-paper-2019/analysis_scripts/
# IPATH=~/mortazavi_lab/bin/TALON-paper-2019/Illumina/
# PBPATH=~/mortazavi_lab/data/test_TALON_plot_tsv_gtf_driver/
# ONTPATH=~/mortazavi_lab/data/test_TALON_plot_tsv_gtf_driver/

# GM12878 pacbio plots
Rscript ${PLOTPATH}plot_novelty_categories.R --db ${PBPATH}full_gencode_v29_2019-03-12.db --w ${PBPATH}GM12878/filtered_transcripts_GM12878.csv --datasets D8,D9 --o ${PBPATH}GM12878/
Rscript ${PLOTPATH}plot_novelty_categories.R --db ${PBPATH}full_gencode_v29_2019-03-12.db --w ${PBPATH}GM12878/filtered_transcripts_GM12878.csv --datasets D8 --o ${PBPATH}GM12878/
Rscript ${PLOTPATH}plot_novelty_categories.R --db ${PBPATH}full_gencode_v29_2019-03-12.db --w ${PBPATH}GM12878/filtered_transcripts_GM12878.csv --datasets D9 --o ${PBPATH}GM12878/
Rscript ${PLOTPATH}plot_pacbio_gene_expression_corr.R --f ${PBPATH}full_gencode_v29_pb_talon_abundance.tsv --color blue --d1 D8 --d2 D9 --celltype GM12878 -o ${PBPATH}GM12878/with_corr_labels/ --correlations --regression_line
Rscript ${PLOTPATH}plot_pacbio_gene_expression_corr.R --f ${PBPATH}full_gencode_v29_pb_talon_abundance.tsv --color blue --d1 D8 --d2 D9 --celltype GM12878 -o ${PBPATH}GM12878/with_corr_labels/ --antisense --intergenic --correlations --regression_line
Rscript ${PLOTPATH}plot_pacbio_transcript_expression_corr.R --f ${PBPATH}GM12878/GM12878_talon_abundance_filtered.tsv --d1 D8 --d2 D9 --celltype GM12878 -o ${PBPATH}GM12878/with_corr_labels/ --correlations --regression_line
Rscript ${PLOTPATH}plot_pacbio_transcript_expression_corr.R --f ${PBPATH}GM12878/GM12878_talon_abundance_filtered.tsv --d1 D8 --d2 D9 --celltype GM12878 -o ${PBPATH}GM12878/with_corr_labels/ --ISM --NIC --NNC --antisense --intergenic --correlations --regression_line
Rscript ${PLOTPATH}plot_detection_by_TPM_for_datasets.R --db ${PBPATH}full_gencode_v29_2019-03-12.db --datasets D8,D9 --ik1 ${IPATH}GM12878/Kallisto/Rep1/abundance.tsv --ik2 ${IPATH}GM12878/Kallisto/Rep2/abundance.tsv --color blue -o ${PBPATH}illumina/GM12878/
Rscript ${PLOTPATH}pacbio_v_illumina_edgeR.R --f ${PBPATH}full_gencode_v29_pb_talon_abundance.tsv --datasets D8,D9 --ik1 ${IPATH}GM12878/Kallisto/Rep1/abundance.tsv --ik2 ${IPATH}GM12878/Kallisto/Rep2/abundance.tsv --color blue -o ${PBPATH}illumina/GM12878/
Rscript ${PLOTPATH}pacbio_v_illumina_edgeR_transcripts.R --f ${PBPATH}GM12878/GM12878_talon_abundance_filtered.tsv --datasets D8,D9 --ik1 ${IPATH}GM12878/Kallisto/Rep1/abundance.tsv --ik2 ${IPATH}GM12878/Kallisto/Rep2/abundance.tsv --color green -o ${PBPATH}illumina/GM12878/

# HepG2 pacbio plots
Rscript ${PLOTPATH}plot_novelty_categories.R --db ${ONTPATH}full_gencode_v29_2019-05-24.db --w ${PBPATH}HepG2/filtered_transcripts_HepG2.csv --datasets D4,D5 --o ${PBPATH}HepG2/
Rscript ${PLOTPATH}plot_novelty_categories.R --db ${ONTPATH}full_gencode_v29_2019-05-24.db --w ${PBPATH}HepG2/filtered_transcripts_HepG2.csv --datasets D4 --o ${PBPATH}HepG2/
Rscript ${PLOTPATH}plot_novelty_categories.R --db ${ONTPATH}full_gencode_v29_2019-05-24.db --w ${PBPATH}HepG2/filtered_transcripts_HepG2.csv --datasets D5 --o ${PBPATH}HepG2/
Rscript ${PLOTPATH}plot_pacbio_gene_expression_corr.R --f ${PBPATH}full_gencode_v29_pb_talon_abundance.tsv --color green --d1 D4 --d2 D5 --celltype HepG2 -o ${PBPATH}HepG2/with_corr_labels/ --correlations --regression_line
Rscript ${PLOTPATH}plot_pacbio_gene_expression_corr.R --f ${PBPATH}full_gencode_v29_pb_talon_abundance.tsv --color green --d1 D4 --d2 D5 --celltype HepG2 -o ${PBPATH}HepG2/with_corr_labels/ --antisense --intergenic --correlations --regression_line
Rscript ${PLOTPATH}plot_pacbio_transcript_expression_corr.R --f ${PBPATH}HepG2/HepG2_talon_abundance_filtered.tsv --d1 D4 --d2 D5 --celltype HepG2 -o ${PBPATH}HepG2/with_corr_labels/ --correlations --regression_line
Rscript ${PLOTPATH}plot_pacbio_transcript_expression_corr.R --f ${PBPATH}HepG2/HepG2_talon_abundance_filtered.tsv --d1 D4 --d2 D5 --celltype HepG2 -o ${PBPATH}HepG2/with_corr_labels/ --ISM --NIC --NNC --antisense --intergenic --correlations --regression_line
Rscript ${PLOTPATH}plot_detection_by_TPM_for_datasets.R --db ${ONTPATH}full_gencode_v29_2019-05-24.db --datasets D4,D5 --ik1 ${IPATH}HepG2/Kallisto/Rep1/abundance.tsv --ik2 ${IPATH}HepG2/Kallisto/Rep2/abundance.tsv --color green -o ${PBPATH}illumina/HepG2/
Rscript ${PLOTPATH}pacbio_v_illumina_edgeR.R --f ${PBPATH}full_gencode_v29_pb_talon_abundance.tsv --datasets D4,D5 --ik1 ${IPATH}illumina/HepG2/Kallisto/Rep1/abundance.tsv --ik2 ${IPATH}illumina/HepG2/Kallisto/Rep2/abundance.tsv --color green -o ${PBPATH}illumina/HepG2/
Rscript ${PLOTPATH}pacbio_v_illumina_edgeR_transcripts.R --f ${PBPATH}HepG2/HepG2_talon_abundance_filtered.tsv --datasets D4,D5 --ik1 ${IPATH}illumina/HepG2/Kallisto/Rep1/abundance.tsv --ik2 ${IPATH}illumina/HepG2/Kallisto/Rep2/abundance.tsv --color green -o ${PBPATH}illumina/HepG2/

# K562 pacbio plots
Rscript ${PLOTPATH}plot_novelty_categories.R --db ${ONTPATH}full_gencode_v29_2019-05-24.db --w ${PBPATH}K562/filtered_transcripts_K562.csv --datasets D10,D11 --o ${PBPATH}K562/
Rscript ${PLOTPATH}plot_novelty_categories.R --db ${ONTPATH}full_gencode_v29_2019-05-24.db --w ${PBPATH}K562/filtered_transcripts_K562.csv --datasets D10 --o ${PBPATH}K562/
Rscript ${PLOTPATH}plot_novelty_categories.R --db ${ONTPATH}full_gencode_v29_2019-05-24.db --w ${PBPATH}K562/filtered_transcripts_K562.csv --datasets D11 --o ${PBPATH}K562/
Rscript ${PLOTPATH}plot_pacbio_gene_expression_corr.R --f ${PBPATH}full_gencode_v29_pb_talon_abundance.tsv --color red --d1 D10 --d2 D11 --celltype K562 -o ${PBPATH}K562/with_corr_labels/ --correlations --regression_line
Rscript ${PLOTPATH}plot_pacbio_gene_expression_corr.R --f ${PBPATH}full_gencode_v29_pb_talon_abundance.tsv --color red --d1 D10 --d2 D11 --celltype K562 -o ${PBPATH}K562/with_corr_labels/ --antisense --intergenic --correlations --regression_line
Rscript ${PLOTPATH}plot_pacbio_transcript_expression_corr.R --f ${PBPATH}K562/K562_talon_abundance_filtered.tsv --d1 D10 --d2 D11 --celltype K562 -o ${PBPATH}K562/with_corr_labels/ --correlations --regression_line
Rscript ${PLOTPATH}plot_pacbio_transcript_expression_corr.R --f ${PBPATH}K562/K562_talon_abundance_filtered.tsv --d1 D10 --d2 D11 --celltype K562 -o ${PBPATH}K562/with_corr_labels/ --ISM --NIC --NNC --antisense --intergenic --correlations --regression_line
Rscript ${PLOTPATH}plot_detection_by_TPM_for_datasets.R --db ${ONTPATH}full_gencode_v29_2019-05-24.db --datasets D10,D11 --ik1 ${IPATH}K562/Kallisto/Rep1/abundance.tsv --ik2 ${IPATH}K562/Kallisto/Rep2/abundance.tsv --color red -o ${PBPATH}illumina/K562/
Rscript ${PLOTPATH}pacbio_v_illumina_edgeR.R --f ${PBPATH}full_gencode_v29_pb_talon_abundance.tsv --datasets D10,D11 --ik1 ${IPATH}illumina/K562/Kallisto/Rep1/abundance.tsv --ik2 ${IPATH}illumina/K562/Kallisto/Rep2/abundance.tsv --color red -o ${PBPATH}illumina/K562/
Rscript ${PLOTPATH}pacbio_v_illumina_edgeR_transcripts.R --f ${PBPATH}K562/K562_talon_abundance_filtered.tsv --datasets D10,D11 --ik1 ${IPATH}illumina/K562/Kallisto/Rep1/abundance.tsv --ik2 ${IPATH}illumina/K562/Kallisto/Rep2/abundance.tsv --color green -o ${PBPATH}illumina/K562/

