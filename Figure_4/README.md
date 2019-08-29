# Figure 4: 5' and 3' completeness by novelty category

Files/paths used to generate the panels of this figure:
```bash
PLOTPATH=../plotting_scripts
OUTPLOTS=plots
mkdir -p ${OUTPLOTS}

PB_GTF=S2_GM12878_talon_observedOnly.gtf
ONT_GTF=S18_GM12878_ont_talon_observedOnly.gtf
CAGE=../CAGE/data/FANTOM5/hg38_CAGE.bed
```
GTF files are available as supplementary tables of the TALON paper.  
To obtain FANTOM5 CAGE data, please see instructions at https://github.com/dewyman/TALON-paper-2019/blob/master/CAGE
To obtain ENCODE RNA-PET data, please see instructions at https://github.com/dewyman/TALON-paper-2019/tree/master/RNA-PET

Software versions:  
* Python 3.7.2
* Bedtools v2.27.1
* R v3.5.1

## Panel A: Percentage of TALON transcript models with CAGE support for their 5' end by novelty category (GM12878 PacBio)
```bash
source activate mypython3.7.2
OUT=CAGE/PacBio_GM12878
mkdir -p ${OUT}
python ../CAGE/run_CAGE_analysis.py \
        --gtf ${PB_GTF} \
        --cage ${CAGE} \
        --maxdist 100 \
        --o ${OUT}/GM12878

Rscript ${PLOTPATH}/plot_support_by_novelty_type.R \
    --f ${OUT}/GM12878_CAGE_results.csv \
    --t CAGE \
    --novelty ${OUT}/transcript_beds/GM12878_novelty.csv \
    --splitISM \
    --ymax 22000 \ 
    -o ${OUTPLOTS}/GM12878_PacBio
```


## Panel B: Percentage of TALON transcript models with a poly(A) motif identified at their 3' end (GM12878 PacBio)

## Panel C: Percentage of TALON transcript models with RNA-PET support for their 5'-3' end pair (GM12878 PacBio)

## Panel D: Percentage of TALON transcript models with CAGE support for their 5' end by novelty category (GM12878 ONT)

## Panel E: Percentage of TALON transcript models with a poly(A) motif identified at their 3' end (GM12878 ONT)

## Panel F: Percentage of TALON transcript models with RNA-PET support for their 5'-3' end pair (GM12878 ONT)

