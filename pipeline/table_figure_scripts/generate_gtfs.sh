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
# HPATH=~/mortazavi_lab/test_TALON_plot_tsv_gtf_driver/HepG2/
# KPATH=~/mortazavi_lab/data/test_TALON_plot_tsv_gtf_driver/K562/
# APATH=~/mortazavi_lab/bin/TALON-paper-2019/analysis_scripts/


# PacBio GTFs

# Table S14
# filtered pacbio GTF
python ${TALONPATH}filter_talon_transcripts.py --db ${ONTPATH}full_gencode_v29_2019-05-24.db -a gencode_v29 --pairings ${PBPATH}pairings_all --o ${PBPATH}whitelist_pb
python ${TALONPATH}create_GTF_from_database.py --db ${ONTPATH}full_gencode_v29_2019-05-24.db -a gencode_v29 --whitelist ${PBPATH}whitelist_pb --o ${PBPATH}full_gencode_v29_pb_tracks -b hg38 --datasets ${PBPATH}datasets_pb

# Table S2
# GM12878 filtered pacbio GTF
python ${TALONPATH}filter_talon_transcripts.py --db ${ONTPATH}full_gencode_v29_2019-05-24.db --a gencode_v29 --pairings ${GMPATH}pairings_GM12878 --o ${GMPATH}whitelist_GM12878
python ${TALONPATH}create_GTF_from_database.py -b hg38 --db ${ONTPATH}full_gencode_v29_2019-05-24.db -a gencode_v29 --observed --datasets ${GMPATH}datasets_GM12878 --o ${GMPATH}GM12878 --whitelist ${GMPATH}whitelist_GM12878

# Table S5
# HepG2 filtered pacbio GTF
python ${TALONPATH}filter_talon_transcripts.py --db ${ONTPATH}full_gencode_v29_2019-05-24.db --a gencode_v29 --pairings ${HPATH}pairings_HepG2 --o ${HPATH}whitelist_HepG2 
python ${TALONPATH}create_GTF_from_database.py -b hg38 --db ${ONTPATH}full_gencode_v29_2019-05-24.db -a gencode_v29 --observed --datasets ${HPATH}datasets_HepG2 --o ${HPATH}HepG2 --whitelist ${HPATH}whitelist_HepG2

# Table S8
# K562 filtered pacbio GTF
python ${TALONPATH}filter_talon_transcripts.py --db ${ONTPATH}full_gencode_v29_2019-05-24.db --a gencode_v29 --pairings ${KPATH}pairings_K562 --o ${KPATH}whitelist_K562
python ${TALONPATH}create_GTF_from_database.py -b hg38 --db ${ONTPATH}full_gencode_v29_2019-05-24.db -a gencode_v29 --observed --datasets ${KPATH}datasets_K562 --o ${KPATH}K562 --whitelist ${KPATH}whitelist_K562

# ONT GTFs
GMPATH=~/mortazavi_lab/data/ont_tier1/GM12878/
HPATH=~/mortazavi_lab/data/ont_tier1/HepG2/
KPATH=~/mortazavi_lab/data/ont_tier1/K562/

# Table S29
# filtered ONT GTF
python ${TALONPATH}filter_talon_transcripts.py --db ${ONTPATH}full_gencode_v29_2019-05-24.db -a gencode_v29 --pairings ${ONTPATH}pairings_all_ont --o ${ONTPATH}whitelist_ont
python ${TALONPATH}create_GTF_from_database.py --db ${ONTPATH}full_gencode_v29_2019-05-24.db -a gencode_v29 --whitelist ${ONTPATH}whitelist_ont --observed --o ${ONTPATH}full_gencode_v29_ont_tracks -b hg38

# Table S17
# GM12878 filtered ONT GTF
python ${TALONPATH}filter_talon_transcripts.py --db ${ONTPATH}full_gencode_v29_2019-05-24.db --a gencode_v29 --pairings ${GMPATH}pairings_ont_GM12878 --o ${GMPATH}ont_whitelist_GM12878
python ${TALONPATH}create_GTF_from_database.py -b hg38 --db ${ONTPATH}full_gencode_v29_2019-05-24.db -a gencode_v29 --observed --datasets ${GMPATH}datasets_GM12878 --o ${GMPATH}GM12878_ont --whitelist ${GMPATH}ont_whitelist_GM12878

# Table S20
# HepG2 filtered ONT GTF
python ${TALONPATH}filter_talon_transcripts.py --db ${ONTPATH}full_gencode_v29_2019-05-24.db --a gencode_v29 --pairings ${HPATH}pairings_ont_HepG2 --o ${HPATH}ont_whitelist_HepG2
python ${TALONPATH}create_GTF_from_database.py -b hg38 --db ${ONTPATH}full_gencode_v29_2019-05-24.db -a gencode_v29 --observed --datasets ${HPATH}datasets_HepG2 --o ${HPATH}HepG2_ont --whitelist ${HPATH}ont_whitelist_HepG2

# Table S23
# K562 filtered ONT GTF
python ${TALONPATH}filter_talon_transcripts.py --db ${ONTPATH}full_gencode_v29_2019-05-24.db --a gencode_v29 --pairings ${KPATH}pairings_ont_K562 --o ${KPATH}ont_whitelist_K562
python ${TALONPATH}create_GTF_from_database.py -b hg38 --db ${ONTPATH}full_gencode_v29_2019-05-24.db -a gencode_v29 --observed --datasets ${KPATH}datasets_K562 --o ${KPATH}K562_ont --whitelist ${KPATH}ont_whitelist_K562

# PacBio + ONT GTFs

# Table S26
# filtered ONT + PacBio GTF
python ${TALONPATH}filter_talon_transcripts.py --db ${ONTPATH}full_gencode_v29_2019-05-24.db -a gencode_v29 --pairings ${ONTPATH}pairings_all_pb_ont --o ${ONTPATH}whitelist_pb_ont
python ${TALONPATH}create_GTF_from_database.py --db ${ONTPATH}full_gencode_v29_2019-05-24.db -a gencode_v29 --whitelist ${ONTPATH}whitelist_pb_ont --observed --o ${ONTPATH}full_gencode_v29_pb_ont_tracks -b hg38
