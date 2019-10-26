## Technology (PacBio vs. ONT vs. Illumina) splice junction comparison

1. Get the tables from the supplemental tables file that we'll be using, and set other paths that we'll be using (to TranscriptClean and to the hg38 reference genome).
```bash
# download the supplementary tables and change this path!
sup_tables=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_10-19/human_TALON/analysis/supplementary_tables/

gm_pb_gtf=${sup_tables}S2_GM12878_talon_observedOnly.gtf
gm_ont_gtf=${sup_tables}S18_GM12878_ont_talon_observedOnly.gtf

TCPATH=~/mortazavi_lab/bin/TranscriptClean/accessory_scripts/
REFPATH=~/mortazavi_lab/ref/hg38/
```

2. Extract splice junctionss from GM12878 PacBio and ONT gtfs using TranscriptClean
```bash
python ${TCPATH}get_SJs_from_gtf.py \
	--f ${gm_pb_gtf} \
	--g ${REFPATH}hg38.fa \
	--o pb_talon_GM12878_sjs.tab

python ${TCPATH}get_SJs_from_gtf.py \
	--f ${gm_ont_gtf} \
	--g ${REFPATH}hg38.fa \
	--o ont_talon_GM12878_sjs.tab
```