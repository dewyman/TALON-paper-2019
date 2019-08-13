Scripts for comparing long read and short read splice junctions.

Input:
- High-confidence SJ file from STAR (obtained by mapping short reads to the reference genome)
- SJ file from long reads (achieved by running TranscriptClean utility on TALON-generated GTF)

Output: 
- File listing the number of SJs unique o short reads, unique to long reads, and shared by both. 

1. Generate GTFs from TALON db for all ONT or all PacBio reads (tables S14 and S29 in the supplement)
```
TALONPATH=~/mortazavi_lab/bin/TALON/post-TALON_tools/
DBPATH=~/mortazavi_lab/data/ont_tier1/

printf "D10,D11,D8,D9,D4,D5" > pb_pairings_all
printf "D4\nD5\nD8\nD9\nD10\nD11" > datasets_pb

python ${TALONPATH}filter_talon_transcripts.py \
	--db ${DBPATH}full_gencode_v29_2019-05-24.db \
	-a gencode_v29 \
	--pairings pb_pairings_all \
	--o whitelist_pb

python ${TALONPATH}create_GTF_from_database.py \
	--db ${DBPATH}full_gencode_v29_2019-05-24.db \
	-a gencode_v29 \
	--whitelist whitelist_pb \
	--o full_gencode_v29_pb \
	-b hg38 \
	--datasets datasets_pb

printf "ONT21,ONT24,ONT25,ONT32,ONT33,ONT34,ONT18,ONT31" > ont_pairings_all
printf "ONT21\nONT24\nONT25\nONT32\nONT33\nONT34\nONT18\nONT31" > datasets_ont

python ${TALONPATH}filter_talon_transcripts.py \
        --db ${DBPATH}full_gencode_v29_2019-05-24.db \
        -a gencode_v29 \
        --pairings ont_pairings_all \
        --o whitelist_ont

python ${TALONPATH}create_GTF_from_database.py \
        --db ${DBPATH}full_gencode_v29_2019-05-24.db \
        -a gencode_v29 \
        --whitelist whitelist_ont \
        --o full_gencode_v29_ont \
        -b hg38 \
        --datasets datasets_ont

```

2. Extract splice junctions from gtfs using TranscriptClean
```
TCPATH=~/mortazavi_lab/bin/TranscriptClean/accessory_scripts/
REFPATH=~/mortazavi_lab/ref/hg38/

python ${TCPATH}get_SJs_from_gtf.py \
	--f full_gencode_v29_pb_talon.gtf \
	--g ${REFPATH}hg38.fa \
	--o pb_talon_sjs.bed

python ${TCPATH}get_SJs_from_gtf.py \
	--f full_gencode_v29_ont_talon.gtf \
	--g ${REFPATH}hg38.fa \
	--o ont_talon_sjs.bed
```

3. Now, let's get the splice junctions present in the Illumina data by mapping with STAR. 
```
qsub run_STAR_illumina.sh 
```
