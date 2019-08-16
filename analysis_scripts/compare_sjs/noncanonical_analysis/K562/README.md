# Analysis of noncanonical splice junctions in K562

The goal is to look at the splice junctions in K562 long reads before TranscriptClean and examine:
1) Their expression level compared to canonical splice junctions  
2) Their reproducibility or lack thereof across the PacBio, ONT, and Illumina platforms.

## Data processing

1) Obtain PacBio K562 splice junctions in pre-TranscriptClean reads (data is divided into 4 SMRT cells total from 2 libraries)
```
cat ../../../../pipeline/D10/Minimap2/1_A01/mapped_FLNC_noScaff.sam \
    ../../../../pipeline/D10/Minimap2/2_B01/mapped_FLNC_noScaff.sam \
    ../../../../pipeline/D10/Minimap2/3_C01/mapped_FLNC_noScaff.sam \
    ../../../../pipeline/D10/Minimap2/4_D01/mapped_FLNC_noScaff.sam \
    ../../../../pipeline/D11/Minimap2/3_C01/mapped_FLNC_noScaff.sam \
    ../../../../pipeline/D11/Minimap2/4_D01/mapped_FLNC_noScaff.sam \
    | grep -v "^@" > PacBio_K562_pre-TC.sam

time python ../../extract_SJs_from_sam.py \
    --sam PacBio_K562_pre-TC.sam \
    --genome ../../../../refs/hg38/hg38.fa \
    --o PacBio_K562_pre-TC

```

2) Obtain ONT K562 splice junctions in pre-TranscriptClean reads
```
cat ../../../../pipeline/ONT31/Minimap2/data/ONT31_mapped.sam \
    ../../../../pipeline/ONT18/Minimap2/data/ONT18_mapped.sam \
    | grep -v "^@" > ONT_K562_pre-TC.sam


time python ../../extract_SJs_from_sam.py \
    --sam ONT_K562_pre-TC.sam \
    --genome ../../../../refs/hg38/hg38.fa \
    --o ONT_K562_pre-TC
```

3) Use the splice junctions we extracted for Illumina using STAR:
```
../../../../Illumina/K562/K562_alignedSJ.out.tab
```

## Analysis

### Expression level


### Reproducibility

1) Extract only noncanonical SJs from each file:
```
awk '{if($5 == 0) print $0}' PacBio_K562_pre-TC_SJs.txt > ncsj_PacBio_K562_pre-TC_SJs.txt
awk '{if($5 == 0) print $0}' ONT_K562_pre-TC_SJs.txt > ncsj_ONT_K562_pre-TC_SJs.txt
awk '{if($5 == 0) print $0}' ../../../../Illumina/K562/K562_alignedSJ.out.tab > ncsj_Illumina_K562_pre-TC_SJs.txt
```

2) Run comparison across platforms
```
mkdir figures
python ../../compare_sjs_venn.py \
        -pb ncsj_PacBio_K562_pre-TC_SJs.txt \
	-ont ncsj_ONT_K562_pre-TC_SJs.txt \
	-illumina ncsj_Illumina_K562_pre-TC_SJs.txt \
	-sample "K562 Noncanonical" \
        --log
```



