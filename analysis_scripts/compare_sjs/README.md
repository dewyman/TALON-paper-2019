## Technology (PacBio vs. ONT vs. Illumina) splice junction comparison

1. We will use GM12878 to compare splice junctions across sequencing platforms. We need the GM12878 PacBio and ONT GTFs from the supplement, tables S2 and S17. These will be stored in pb_gtfs/ and ont_gtfs/ respectively.

2. Extract splice junctions from GM12878 PacBio and ONT gtfs using TranscriptClean
```
conda activate python2.7
TCPATH=~/mortazavi_lab/bin/TranscriptClean/accessory_scripts/
REFPATH=~/mortazavi_lab/ref/hg38/

python ${TCPATH}get_SJs_from_gtf.py \
	--f pb_gtfs/GM12878_talon_observedOnly.gtf \
	--g ${REFPATH}hg38.fa \
	--o pb_talon_GM12878_sjs.tab

python ${TCPATH}get_SJs_from_gtf.py \
	--f ont_gtfs/GM12878_ont_talon_observedOnly.gtf \
	--g ${REFPATH}hg38.fa \
	--o ont_talon_GM12878_sjs.tab
```

3. Now, let's get the splice junctions present in the Illumina data by mapping with STAR. 
```
qsub run_STAR_illumina_GM12878.sh
```

4. Create a venn diagram demonstrating which splice junctions are present in which dataset.
```
conda activate base
mkdir figures
python compare_sjs_venn.py \
	-pb pb_talon_GM12878_sjs.tab \
	-ont ont_talon_GM12878_sjs.tab \
	-illumina GM12878_alignedSJ.out.tab \
	-sample GM12878
```

5. We wanted to see how TranscriptClean alters SJ support by other technologies, so we pulled the SJs out of the pre-TranscriptClean sam files
```

```

5. We also ran the above analysis for the other 2 cell lines:
 ```
# HepG2 (requires tables S5 and S20)
conda activate python2.7
TCPATH=~/mortazavi_lab/bin/TranscriptClean/accessory_scripts/
REFPATH=~/mortazavi_lab/ref/hg38/

python ${TCPATH}get_SJs_from_gtf.py \
	--f pb_gtfs/HepG2_talon_observedOnly.gtf \
	--g ${REFPATH}hg38.fa \
	--o pb_talon_HepG2_sjs.tab

python ${TCPATH}get_SJs_from_gtf.py \
	--f ont_gtfs/HepG2_ont_talon_observedOnly.gtf \
	--g ${REFPATH}hg38.fa \
	--o ont_talon_HepG2_sjs.tab

qsub run_STAR_illumina_HepG2.sh

conda activate base
python compare_sjs_venn.py \
	-pb pb_talon_HepG2_sjs.tab \
	-ont ont_talon_HepG2_sjs.tab \
	-illumina HepG2_alignedSJ.out.tab \
	-sample HepG2

# K562 (requires tables S8 and S23)
conda activate python2.7
TCPATH=~/mortazavi_lab/bin/TranscriptClean/accessory_scripts/
REFPATH=~/mortazavi_lab/ref/hg38/

python ${TCPATH}get_SJs_from_gtf.py \
	--f pb_gtfs/K562_talon_observedOnly.gtf \
	--g ${REFPATH}hg38.fa \
	--o pb_talon_K562_sjs.tab

python ${TCPATH}get_SJs_from_gtf.py \
	--f ont_gtfs/K562_ont_talon_observedOnly.gtf
	--g ${REFPATH}hg38.fa \
	--o ont_talon_K562_sjs.tab

qsub run_STAR_illumina_K562.sh

conda activate base
python compare_sjs_venn.py \
	-pb pb_talon_K562_sjs.tab \
	-ont ont_talon_K562_sjs.tab \
	-illumina K562_alignedSJ.out.tab \
	-sample K562 
 ```

<img align="center" width="500" src="figures/GM12878_venn.png">
<img align="center" width="500" src="figures/HepG2_venn.png">
<img align="center" width="500" src="figures/K562_venn.png">

## Splice junction novelty types

We wanted to investigate the support for and percent makeup for splice junctions of different novelty types. In this case we define splice junction novelty as the following: 

* **Known junction:** A splice junction with the whole splice junction (splice donor and acceptor) seen in the Gencode annotation
* **NIC junction:** A splice junction where both the splice donor and acceptor individually are in the Gencode annotation, but are never seen together
* **NNC junction:** A splice junction where at least one of the splice donor or acceptor is not present in the Gencode annotation

1. Going off of this, we first obtain the splice junction file from Gencode v29 using the TranscriptClean utility.
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

2. We then use the splice junctions we extracted from the PacBio GM12878 gtf (Table S2) to label each splice junction with its novelty type
```
conda activate base
python label_sj_novelty.py \
	-sj pb_talon_GM12878_sjs.tab \
	-ref_sj gencode_v29_sjs.tab 
```

3. Now we can look at the proportions of known, nnc, and nic splice junctions present in our data in the form of a bar chart
```
python plot_sj_novelty_counts.py \
	-sj pb_talon_GM12878_sjs_novelty.tab \
	-sample "PacBio GM12878"
```

<img align="left" width="500" src="figures/PacBio_GM12878_sj_novelty.png">

4. Get splice junction novelty types for our Illumina splice junction file so we can look for support for novel splice junctions in our short-read data
```
python label_sj_novelty.py \
	-sj GM12878_alignedSJ.out.tab \
	-ref_sj gencode_v29_sjs.tab 
```

5. Let's see what percentage of novel splice junctions are supported by Illumina data. 
```
python plot_sj_novelty_counts.py \
	-sj pb_talon_GM12878_sjs_novelty.tab \
	-sample "PacBio GM12878" \
	--extra_support GM12878_alignedSJ.out_novelty.tab \
	--support_name Illumina
```

<img align="left" width="500" src="figures/PacBio_GM12878_sj_novelty_Illumina_support.png">

<!-- 
6. Look at NNC makeup with finer resolution, are both the donor and acceptor completely novel or is has one or the other been seen before?
```
# coming soon
``` -->

6. Also perform the above analysis for HepG2 and K562
```
# HepG2
python label_sj_novelty.py \
	-sj pb_talon_HepG2_sjs.tab \
	-ref_sj gencode_v29_sjs.tab 

python plot_sj_novelty_counts.py \
	-sj pb_talon_HepG2_sjs_novelty.tab \
	-sample "PacBio HepG2"

python label_sj_novelty.py \
	-sj HepG2_alignedSJ.out.tab \
	-ref_sj gencode_v29_sjs.tab 

python plot_sj_novelty_counts.py \
	-sj pb_talon_HepG2_sjs_novelty.tab \
	-sample "PacBio HepG2" \
	--extra_support HepG2_alignedSJ.out_novelty.tab \
	--support_name Illumina

# K562
python label_sj_novelty.py \
	-sj pb_talon_K562_sjs.tab \
	-ref_sj gencode_v29_sjs.tab 

python plot_sj_novelty_counts.py \
	-sj pb_talon_K562_sjs_novelty.tab \
	-sample "PacBio K562"

python label_sj_novelty.py \
	-sj K562_alignedSJ.out.tab \
	-ref_sj gencode_v29_sjs.tab 

python plot_sj_novelty_counts.py \
	-sj pb_talon_K562_sjs_novelty.tab \
	-sample "PacBio K562" \
	--extra_support K562_alignedSJ.out_novelty.tab \
	--support_name Illumina
```

<img align="left" width="500" src="figures/PacBio_HepG2_sj_novelty.png">
<img align="left" width="500" src="figures/PacBio_HepG2_sj_novelty_Illumina_support.png">

<img align="left" width="500" src="figures/PacBio_K562_sj_novelty.png">
<img align="left" width="500" src="figures/PacBio_K562_sj_novelty_Illumina_support.png">


7. Also perform this analysis for the ONT data TODO
```
# GM12878
python label_sj_novelty.py \
	-sj ont_talon_GM12878_sjs.tab \
	-ref_sj gencode_v29_sjs.tab 

python plot_sj_novelty_counts.py \
	-sj ont_talon_GM12878_sjs_novelty.tab \
	-sample "ONT GM12878"

python plot_sj_novelty_counts.py \
	-sj ont_talon_GM12878_sjs_novelty.tab \
	-sample "ONT GM12878" \
	--extra_support GM12878_alignedSJ.out_novelty.tab \
	--support_name Illumina

# HepG2
python label_sj_novelty.py \
	-sj ont_talon_HepG2_sjs.tab \
	-ref_sj gencode_v29_sjs.tab 

python plot_sj_novelty_counts.py \
	-sj ont_talon_HepG2_sjs_novelty.tab \
	-sample "ONT HepG2"

python plot_sj_novelty_counts.py \
	-sj ont_talon_HepG2_sjs_novelty.tab \
	-sample "ONT HepG2" \
	--extra_support HepG2_alignedSJ.out_novelty.tab \
	--support_name Illumina

# K562
python label_sj_novelty.py \
	-sj ont_talon_K562_sjs.tab \
	-ref_sj gencode_v29_sjs.tab 

python plot_sj_novelty_counts.py \
	-sj ont_talon_K562_sjs_novelty.tab \
	-sample "ONT K562"

python plot_sj_novelty_counts.py \
	-sj ont_talon_K562_sjs_novelty.tab \
	-sample "ONT K562" \
	--extra_support K562_alignedSJ.out_novelty.tab \
	--support_name Illumina


```

<img align="left" width="500" src="figures/PacBio_HepG2_sj_novelty.png">
<img align="left" width="500" src="figures/ONT_HepG2_sj_novelty_Illumina_support.png">

<img align="left" width="500" src="figures/PacBio_HepG2_sj_novelty.png">
<img align="left" width="500" src="figures/ONT_HepG2_sj_novelty_Illumina_support.png">

<img align="left" width="500" src="figures/PacBio_K562_sj_novelty.png">
<img align="left" width="500" src="figures/ONT_K562_sj_novelty_Illumina_support.png">
