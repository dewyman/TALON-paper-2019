ONTPATH=~/mortazavi_lab/data/ont_tier1/
GM_ONT=~/mortazavi_lab/data/ont_tier1/GM12878/
H_ONT=~/mortazavi_lab/data/ont_tier1/HepG2/
K_ONT=~/mortazavi_lab/data/ont_tier1/K562/
PBPATH=~/mortazavi_lab/data/tier1_filtered/
GM_PB=~/mortazavi_lab/data/tier1_filtered/GM12878/
H_PB=~/mortazavi_lab/data/tier1_filtered/HepG2/
K_PB=~/mortazavi_lab/data/tier1_filtered/K562/

# test 
# set -e
# ONTPATH=~/mortazavi_lab/data/test_TALON_plot_tsv_gtf_driver/
# GM_ONT=~/mortazavi_lab/data/test_TALON_plot_tsv_gtf_driver/GM12878/
# H_ONT=~/mortazavi_lab/data/test_TALON_plot_tsv_gtf_driver/HepG2/
# K_ONT=~/mortazavi_lab/data/test_TALON_plot_tsv_gtf_driver/K562/
# PBPATH=~/mortazavi_lab/data/test_TALON_plot_tsv_gtf_driver/
# GM_PB=~/mortazavi_lab/data/test_TALON_plot_tsv_gtf_driver/GM12878/
# H_PB=~/mortazavi_lab/data/test_TALON_plot_tsv_gtf_driver/HepG2/
# K_PB=~/mortazavi_lab/data/test_TALON_plot_tsv_gtf_driver/K562/

script=~/mortazavi_lab/bin/TALON-paper-2019/pipeline/table_figure_scripts/rename_abundance_file_datasets.py

# Table S15 
# unfiltered pacbio abundance file
python $script --f ${PBPATH}full_gencode_v29_pb_talon_abundance.tsv

# Table S28 
# filtered pacbio abundance file
python $script --f ${PBPATH}full_gencode_v29_pb_talon_abundance_filtered.tsv

# Table S3 
# GM12878 unfiltered pacbio abundance file
python $script --f ${GM_PB}GM12878_talon_abundance.tsv

# Table S4
# GM12878 filtered pacbio abundance file
python $script --f ${GM_PB}GM12878_talon_abundance_filtered.tsv

# Table S6 
# HepG2 unfiltered pacbio abundance file
python $script --f ${H_PB}HepG2_talon_abundance.tsv

# Table S7
# HepG2 filtered pacbio abundance file
python $script --f ${H_PB}HepG2_talon_abundance_filtered.tsv

# Table S9 
# K562 unfiltered pacbio abundance file
python $script --f ${K_PB}K562_talon_abundance.tsv

# Table S10
# K562 filtered pacbio abundance file
python $script --f ${K_PB}K562_talon_abundance_filtered.tsv

# Table S30
# unfiltered ONT abundance file
python $script --f ${ONTPATH}full_gencode_v29_ont_talon_abundance.tsv

# Table S31
# filtered ONT abundance file
python $script --f ${ONTPATH}full_gencode_v29_ont_talon_abundance_filtered.tsv

# Table S18
# GM12878 unfiltered ONT abundance file
python $script --f ${GM_ONT}GM12878_ont_talon_abundance.tsv

# Table S19
# GM12878 filtered ONT abundance file
python $script --f ${GM_ONT}GM12878_ont_talon_abundance_filtered.tsv

# Table S21
# HepG2 unfiltered ONT abundance file
python $script --f ${H_ONT}HepG2_ont_talon_abundance.tsv

# Table S22
# HepG2 filtered ONT abundance file
python $script --f ${H_ONT}HepG2_ont_talon_abundance_filtered.tsv

# Table S24
# K562 unfiltered ONT abundance file
python $script --f ${K_ONT}K562_ont_talon_abundance.tsv

# Table S25
# K562 filtered ONT abundance file
python $script --f ${K_ONT}K562_ont_talon_abundance_filtered.tsv

# Table S27
# unfiltered ONT + PacBio abundance file
python $script --f ${ONTPATH}full_gencode_v29_pb_ont_talon_abundance.tsv

# Table S28 
# filtered ONT + PacBio abundance file
python $script --f ${ONTPATH}full_gencode_v29_pb_ont_talon_abundance_filtered.tsv


