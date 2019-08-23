## Illumina/Gencode splice junction support for whole transcripts

We want to determine what the percentage of PacBio TALON transcripts with all of their splice junctions supported by short-read data is. We decided to subset on novelty to provide some more insight. We'll show this step-by-step example for GM12878 first. 

1. Get the splice junctions from the Illumina data using STAR 
```
qsub ../run_STAR_illumina_GM12878.sh
```

2. Get splice junctions from the Gencode v29 annotation
```
conda activate python2.7
TCPATH=~/mortazavi_lab/bin/TranscriptClean/accessory_scripts/
ANNPATH=~/mortazavi_lab/ref/gencode.v29/
REFPATH=~/mortazavi_lab/ref/hg38/


python ${TCPATH}get_SJs_from_gtf.py \
    --f ${ANNPATH}gencode.v29.annotation.gtf \
	--g ${REFPATH}hg38.fa \
	--o gencode_v29_sjs.tab
```

2. We'll use the GM12878 PacBio GTF (table S2) to compare to the Illumina data for this example, and we'll plot Illumina support for the splice junctions seen in PacBio based on novelty type. 
```
python get_isoform_sj_support.py \
	-gtf ../pb_gtfs/GM12878_talon_observedOnly.gtf \
	-ref_sj_1 ../GM12878_alignedSJ.out.tab \
	-ref_sj_2 ../gencode_v29_sjs.tab \
	-sample pb_GM12878
```

3. Finally, create a summary table illustrating how many transcript isoforms have Gencode or Illumina support for all of its splice junctions.
```
python gen_isoform_support_table.py \
	-csv pb_GM12878_isoform_sj_support.csv \
	-sample pb_GM12878
```
<!-- 
testing
```
python get_isoform_sj_support.py \
	-gtf one_transcript_monoexon.gtf \
	-ref_sj_1 ../GM12878_alignedSJ.out.tab \
	-ref_sj_2 ../gencode_v29_sjs.tab \
	-sample test
``` -->

4. We can run this for the other cell lines as well:  
```
# HepG2
python get_isoform_sj_support.py \
	-gtf ../pb_gtfs/HepG2_talon_observedOnly.gtf \
	-ref_sj_1 ../HepG2_alignedSJ.out.tab \
	-ref_sj_2 ../gencode_v29_sjs.tab \
	-sample pb_HepG2

python gen_isoform_support_table.py \
	-csv pb_HepG2_isoform_sj_support.csv \
	-sample pb_HepG2

# K562
python get_isoform_sj_support.py \
	-gtf ../pb_gtfs/K562_talon_observedOnly.gtf \
	-ref_sj_1 ../K562_alignedSJ.out.tab \
	-ref_sj_2 ../gencode_v29_sjs.tab \
	-sample pb_K562

python gen_isoform_support_table.py \
	-csv pb_K562_isoform_sj_support.csv \
	-sample pb_K562
```

5. And we can run the same analysis on our Nanopore data:
```
# GM12878
python get_isoform_sj_support.py \
	-gtf ../ont_gtfs/GM12878_ont_talon_observedOnly.gtf \
	-ref_sj_1 ../GM12878_alignedSJ.out.tab \
	-ref_sj_2 ../gencode_v29_sjs.tab \
	-sample ont_GM12878

python gen_isoform_support_table.py \
	-csv ont_GM12878_isoform_sj_support.csv \
	-sample ont_GM12878

# HepG2
python get_isoform_sj_support.py \
	-gtf ../ont_gtfs/HepG2_ont_talon_observedOnly.gtf \
	-ref_sj_1 ../HepG2_alignedSJ.out.tab \
	-ref_sj_2 ../gencode_v29_sjs.tab \
	-sample ont_HepG2

python gen_isoform_support_table.py \
	-csv ont_HepG2_isoform_sj_support.csv \
	-sample ont_HepG2

# K562
python get_isoform_sj_support.py \
	-gtf ../ont_gtfs/K562_ont_talon_observedOnly.gtf \
	-ref_sj_1 ../K562_alignedSJ.out.tab \
	-ref_sj_2 ../gencode_v29_sjs.tab \
	-sample ont_K562

python gen_isoform_support_table.py \
	-csv ont_K562_isoform_sj_support.csv \
	-sample ont_K562
```

6. We also wanted to see how this works in cortex/hippocampus so we'll do the same analysis there
```
# cortex
python get_isoform_sj_support.py \
	-gtf ../mouse_brain/cortex_talon_observedOnly.gtf  \
	-ref_sj_1 ../cortex_alignedSJ.out.tab \
	-ref_sj_2 ../gencode_vM21_sjs.tab \
	-sample cortex

python gen_isoform_support_table.py \
	-csv cortex_isoform_sj_support.csv \
	-sample cortex

# hippocampus
python get_isoform_sj_support.py \
	-gtf ../mouse_brain/hippocampus_talon_observedOnly.gtf  \
	-ref_sj_1 ../hippocampus_alignedSJ.out.tab \
	-ref_sj_2 ../gencode_vM21_sjs.tab \
	-sample hippocampus

python gen_isoform_support_table.py \
	-csv hippocampus_isoform_sj_support.csv \
	-sample hippocampus
```
