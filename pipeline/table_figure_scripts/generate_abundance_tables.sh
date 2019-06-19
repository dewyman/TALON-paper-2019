# paths
PBPATH=~/mortazavi_lab/data/tier1_filtered/
ONTPATH=~/mortazavi_lab/data/ont_tier1/
TALONPATH=~/mortazavi_lab/bin/TALON/post-TALON_tools/
GMPATH=~/mortazavi_lab/data/tier1_filtered/GM12878/
HPATH=~/mortazavi_lab/data/tier1_filtered/HepG2/
KPATH=~/mortazavi_lab/data/tier1_filtered/K562/
APATH=~/mortazavi_lab/bin/TALON-paper-2019/analysis_scripts/

# # test
# set -e 
# PBPATH=~/mortazavi_lab/data/test_TALON_plot_tsv_gtf_driver/
# ONTPATH=~/mortazavi_lab/data/test_TALON_plot_tsv_gtf_driver/
# TALONPATH=~/mortazavi_lab/bin/TALON/post-TALON_tools/
# GMPATH=~/mortazavi_lab/data/test_TALON_plot_tsv_gtf_driver/GM12878/
# HPATH=~/mortazavi_lab/data/test_TALON_plot_tsv_gtf_driver/HepG2/
# KPATH=~/mortazavi_lab/data/test_TALON_plot_tsv_gtf_driver/K562/
# APATH=~/mortazavi_lab/bin/TALON-paper-2019/analysis_scripts/

# PacBio tables

# Table S15 
# unfiltered pacbio abundance file
python ${TALONPATH}create_abundance_file_from_database.py --db ${ONTPATH}full_gencode_v29_2019-05-24.db --a gencode_v29 --o ${PBPATH}full_gencode_v29_pb --b hg38
printf "D4\nD5\nD8\nD9\nD10\nD11" > ${PBPATH}datasets_pb
python ${APATH}get_dataset_specific_abundance.py --infile ${PBPATH}full_gencode_v29_pb_talon_abundance.tsv --d ${PBPATH}datasets_pb

# Table S28 
# filtered pacbio abundance file
printf "D4,D5,D8,D9,D10,D11" > ${PBPATH}pairings_all
python ${TALONPATH}create_abundance_file_from_database.py --db ${ONTPATH}full_gencode_v29_2019-05-24.db --a gencode_v29 --o ${PBPATH}full_gencode_v29_pb --pairings ${PBPATH}pairings_all --filter -b hg38
python ${APATH}get_dataset_specific_abundance.py --infile ${PBPATH}full_gencode_v29_pb_talon_abundance_filtered.tsv --d ${PBPATH}datasets_pb

# Table S3 
# GM12878 unfiltered pacbio abundance file
printf "D8\nD9" > ${GMPATH}datasets_pb_GM12878
python ${APATH}get_dataset_specific_abundance.py --infile ${PBPATH}full_gencode_v29_pb_talon_abundance.tsv --outfile ${GMPATH}GM12878_talon_abundance.tsv --d ${GMPATH}datasets_pb_GM12878

# Table S4
# GM12878 filtered pacbio abundance file
printf "D8,D9" > ${GMPATH}pairings_GM12878
python ${TALONPATH}create_abundance_file_from_database.py --db ${ONTPATH}full_gencode_v29_2019-05-24.db --a gencode_v29 --o ${GMPATH}GM12878 --pairings ${GMPATH}pairings_GM12878 --filter -b hg38
python ${APATH}get_dataset_specific_abundance.py --infile ${GMPATH}GM12878_talon_abundance_filtered.tsv --d ${GMPATH}datasets_pb_GM12878

# Table S6 
# HepG2 unfiltered pacbio abundance file
printf "D4\nD5\n" > ${HPATH}datasets_pb_HepG2
python ${APATH}get_dataset_specific_abundance.py --infile ${PBPATH}full_gencode_v29_pb_talon_abundance.tsv --outfile ${HPATH}HepG2_talon_abundance.tsv --d ${HPATH}datasets_pb_HepG2

# Table S7
# HepG2 filtered pacbio abundance file
printf "D4,D5" > ${HPATH}pairings_HepG2
python ${TALONPATH}create_abundance_file_from_database.py --db ${ONTPATH}full_gencode_v29_2019-05-24.db --a gencode_v29 --o ${HPATH}HepG2 --pairings ${HPATH}pairings_HepG2 --filter -b hg38
python ${APATH}get_dataset_specific_abundance.py --infile ${HPATH}HepG2_talon_abundance_filtered.tsv --d ${HPATH}datasets_pb_HepG2

# Table S9 
# K562 unfiltered pacbio abundance file
printf "D10\nD11" > ${KPATH}datasets_pb_K562
python ${APATH}get_dataset_specific_abundance.py --infile ${PBPATH}full_gencode_v29_pb_talon_abundance.tsv --outfile ${KPATH}K562_talon_abundance.tsv --d ${KPATH}datasets_pb_K562

# Table S10
# K562 filtered pacbio abundance file
printf "D10,D11" > ${KPATH}pairings_K562
python ${TALONPATH}create_abundance_file_from_database.py --db ${ONTPATH}full_gencode_v29_2019-05-24.db --a gencode_v29 --o ${KPATH}K562 --pairings ${KPATH}pairings_K562 --filter -b hg38
python ${APATH}get_dataset_specific_abundance.py --infile ${KPATH}K562_talon_abundance_filtered.tsv --d ${KPATH}datasets_pb_K562

# ONT tables

GMPATH=~/mortazavi_lab/data/ont_tier1/GM12878/
HPATH=~/mortazavi_lab/data/ont_tier1/HepG2/
KPATH=~/mortazavi_lab/data/ont_tier1/K562/

# Table S30
# unfiltered ONT abundance file
printf "ONT21\nONT24\nONT25\nONT32\nONT33\nONT34\nONT18\nONT31" > ${ONTPATH}datasets_ont
python ${TALONPATH}create_abundance_file_from_database.py --db ${ONTPATH}full_gencode_v29_2019-05-24.db --a gencode_v29 --o ${ONTPATH}full_gencode_v29_ont -b hg38
python ${APATH}get_dataset_specific_abundance.py --infile ${ONTPATH}full_gencode_v29_ont_talon_abundance.tsv --d ${ONTPATH}datasets_ont

# Table S31
# filtered ONT abundance file
printf "ONT21,ONT24,ONT25,ONT32,ONT33,ONT34,ONT18,ONT31" > ${ONTPATH}pairings_ont
python ${TALONPATH}create_abundance_file_from_database.py --db ${ONTPATH}full_gencode_v29_2019-05-24.db --a gencode_v29 --o ${ONTPATH}full_gencode_v29_ont --pairings ${ONTPATH}pairings_ont --filter -b hg38
python ${APATH}get_dataset_specific_abundance.py --infile ${ONTPATH}full_gencode_v29_ont_talon_abundance_filtered.tsv --d ${ONTPATH}datasets_ont

# Table S18
# GM12878 unfiltered ONT abundance file
printf "ONT21\nONT24\nONT25" > ${GMPATH}datasets_ont_GM12878
python ${APATH}get_dataset_specific_abundance.py --infile ${ONTPATH}full_gencode_v29_ont_talon_abundance.tsv --outfile ${GMPATH}GM12878_ont_talon_abundance.tsv --d ${GMPATH}datasets_ont_GM12878

# Table S19
# GM12878 filtered ONT abundance file
printf "ONT21,ONT24,ONT25" > ${GMPATH}pairings_ont_GM12878
python ${TALONPATH}create_abundance_file_from_database.py --db ${ONTPATH}full_gencode_v29_2019-05-24.db --a gencode_v29 --o ${GMPATH}GM12878_ont --pairings ${GMPATH}pairings_ont_GM12878 --filter -b hg38
python ${APATH}get_dataset_specific_abundance.py --infile ${GMPATH}GM12878_ont_talon_abundance_filtered.tsv --d ${GMPATH}datasets_ont_GM12878

# Table S21
# HepG2 unfiltered ONT abundance file
printf "ONT32\nONT33\nONT34" > ${HPATH}datasets_ont_HepG2
python ${APATH}get_dataset_specific_abundance.py --infile ${ONTPATH}full_gencode_v29_ont_talon_abundance.tsv --outfile ${HPATH}HepG2_ont_talon_abundance.tsv --d ${HPATH}datasets_ont_HepG2

# Table S22
# HepG2 filtered ONT abundance file
printf "ONT32,ONT33,ONT34" > ${HPATH}pairings_ont_HepG2
python ${TALONPATH}create_abundance_file_from_database.py --db ${ONTPATH}full_gencode_v29_2019-05-24.db --a gencode_v29 --o ${HPATH}HepG2_ont --pairings ${HPATH}pairings_ont_HepG2 --filter -b hg38
python ${APATH}get_dataset_specific_abundance.py --infile ${HPATH}HepG2_ont_talon_abundance_filtered.tsv --d ${HPATH}datasets_ont_HepG2

# Table S24
# K562 unfiltered ONT abundance file
printf "ONT18\nONT31" > ${KPATH}datasets_ont_K562
python ${APATH}get_dataset_specific_abundance.py --infile ${ONTPATH}full_gencode_v29_ont_talon_abundance.tsv --outfile ${KPATH}K562_ont_talon_abundance.tsv --d ${KPATH}datasets_ont_K562

# Table S25
# K562 filtered ONT abundance file
printf "ONT18,ONT31" > ${KPATH}pairings_ont_K562
python ${TALONPATH}create_abundance_file_from_database.py --db ${ONTPATH}full_gencode_v29_2019-05-24.db --a gencode_v29 --o ${KPATH}K562_ont --pairings ${KPATH}pairings_ont_K562 --filter -b hg38
python ${APATH}get_dataset_specific_abundance.py --infile ${KPATH}K562_ont_talon_abundance_filtered.tsv --d ${KPATH}datasets_ont_K562

# ONT + PacBio tables

# Table S27
# unfiltered ONT + PacBio abundance file
python ${TALONPATH}create_abundance_file_from_database.py --db ${ONTPATH}full_gencode_v29_2019-05-24.db -a gencode_v29 --o ${ONTPATH}full_gencode_v29_pb_ont -b hg38

# Table S28 
# filtered ONT + PacBio abundance file
printf "ONT21,ONT24,ONT25,D8,D9,ONT32,ONT33,ONT34,D4,D5,ONT18,ONT31,D10,D11" > ${ONTPATH}pairings_pb_ont
python ${TALONPATH}create_abundance_file_from_database.py --db ${ONTPATH}full_gencode_v29_2019-05-24.db -a gencode_v29 --filter --pairings ${ONTPATH}pairings_pb_ont --o ${ONTPATH}full_gencode_v29_pb_ont -b hg38

# these are for pacbio + ONT comparisons
# pb+ont GM12878 biorep and platform rep filtered abundance file
printf "D8\nD9\nONT21\nONT24\nONT25" > ${GMPATH}datasets_ont_pb_GM12878
python ${TALONPATH}create_abundance_file_from_database.py --db ${ONTPATH}full_gencode_v29_2019-05-24.db -a gencode_v29 --o ${GMPATH}GM12878_pb_ont --filter --pairings ${GMPATH}pairings_ont_pb_GM12878 -b hg38
python ${APATH}get_dataset_specific_abundance.py --infile ${GMPATH}GM12878_pb_ont_talon_abundance_filtered.tsv --d ${GMPATH}datasets_ont_pb_GM12878

# pb+ont HepG2 biorep and platform rep filtered abundance file
printf "D4\nD5\nONT32\nONT33\nONT34" > ${HPATH}datasets_ont_pb_HepG2
python ${TALONPATH}create_abundance_file_from_database.py --db ${ONTPATH}full_gencode_v29_2019-05-24.db -a gencode_v29 --o ${HPATH}HepG2_pb_ont --filter --pairings ${HPATH}pairings_ont_pb_HepG2 -b hg38
python ${APATH}get_dataset_specific_abundance.py --infile ${HPATH}HepG2_pb_ont_talon_abundance_filtered.tsv --d ${HPATH}datasets_ont_pb_HepG2

# pb+ont K562 biorep and platform rep filtered abundance file
printf "D10\nD11\nONT18\nONT31" > ${KPATH}datasets_ont_pb_K562
python ${TALONPATH}create_abundance_file_from_database.py --db ${ONTPATH}full_gencode_v29_2019-05-24.db -a gencode_v29 --o ${KPATH}K562_pb_ont --filter --pairings ${KPATH}pairings_ont_pb_K562 -b hg38
python ${APATH}get_dataset_specific_abundance.py --infile ${KPATH}K562_pb_ont_talon_abundance_filtered.tsv --d ${KPATH}datasets_ont_pb_K562

