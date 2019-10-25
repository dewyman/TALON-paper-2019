# Analysis of noncanonical splice junctions in GM12878

The goal is to look at the splice junctions in GM12878 long reads before TranscriptClean and examine:
1) Their expression level compared to canonical splice junctions  
2) Their reproducibility or lack thereof across the PacBio, ONT, and Illumina platforms.

## Data processing

1) Obtain PacBio GM12878 splice junctions in pre-TranscriptClean reads (data is divided into 4 SMRT cells total from 2 libraries)
```
cat ../../../../pipeline/D8/Minimap2/3_C01/mapped_FLNC_noScaff.sam \
    ../../../../pipeline/D8/Minimap2/4_D01/mapped_FLNC_noScaff.sam \
    ../../../../pipeline/D9/Minimap2/1_A01/mapped_FLNC_noScaff.sam \
    ../../../../pipeline/D9/Minimap2/2_B01/mapped_FLNC_noScaff.sam \
    | grep -v "^@" > PacBio_GM12878_pre-TC.sam

time python ../../extract_SJs_from_sam.py \
    --sam PacBio_GM12878_pre-TC.sam \
    --genome ../../../../refs/hg38/hg38.fa \
    --o PacBio_GM12878_pre-TC

```

2) Obtain ONT GM12878 splice junctions in pre-TranscriptClean reads
```
cat ../../../../pipeline/ONT24/Minimap2/data/ONT24_mapped.sam \
    ../../../../pipeline/ONT25/Minimap2/data/ONT25_mapped.sam \
    | grep -v "^@" > ONT_GM12878_pre-TC.sam


time python ../../extract_SJs_from_sam.py \
    --sam ONT_GM12878_pre-TC.sam \
    --genome ../../../../refs/hg38/hg38.fa \
    --o ONT_GM12878_pre-TC
```

3) Use the splice junctions we extracted for Illumina using STAR:
```
/share/crsp/lab/seyedam/share/TALON_paper_data/illumina_sjs/GM12878_alignedSJ.out.tab
```

## Analysis

### Expression level


### Reproducibility

1) Extract only noncanonical SJs from each file:
```
awk '{if($5 == 0) print $0}' PacBio_GM12878_pre-TC_SJs.txt > ncsj_PacBio_GM12878_pre-TC_SJs.txt
awk '{if($5 == 0) print $0}' ONT_GM12878_pre-TC_SJs.txt > ncsj_ONT_GM12878_pre-TC_SJs.txt
awk '{if($5 == 0) print $0}' /share/crsp/lab/seyedam/share/TALON_paper_data/illumina_sjs/GM12878_alignedSJ.out.tab > ncsj_Illumina_GM12878_pre-TC_SJs.txt
```

2) Run comparison across platforms
```
mkdir figures
python ../../compare_sjs_venn.py \
        -pb ncsj_PacBio_GM12878_pre-TC_SJs.txt \
	-ont ncsj_ONT_GM12878_pre-TC_SJs.txt \
	-illumina ncsj_Illumina_GM12878_pre-TC_SJs.txt \
	-sample "GM12878 Noncanonical" \
        --log
```



