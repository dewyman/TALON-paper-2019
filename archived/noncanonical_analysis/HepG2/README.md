# Analysis of noncanonical splice junctions in HepG2

The goal is to look at the splice junctions in HepG2 long reads before TranscriptClean and examine:
1) Their expression level compared to canonical splice junctions  
2) Their reproducibility or lack thereof across the PacBio, ONT, and Illumina platforms.

## Data processing

1) Obtain PacBio HepG2 splice junctions in pre-TranscriptClean reads (data is divided into 4 SMRT cells total from 2 libraries)
```
cat ../../../../pipeline/D4/Minimap2/1_A01/mapped_FLNC_noScaff.sam \
    <(grep -v "^@" ../../../../pipeline/D4/Minimap2/2_B01/mapped_FLNC_noScaff.sam) \
    <(grep -v "^@" ../../../../pipeline/D5/Minimap2/3_C01/mapped_FLNC_noScaff.sam) \
    <(grep -v "^@" ../../../../pipeline/D5/Minimap2/4_D01/mapped_FLNC_noScaff.sam) \
    <(grep -v "^@" ../../../../pipeline/D5/Minimap2/3_C01_2/mapped_FLNC_noScaff.sam) \
    <(grep -v "^@" ../../../../pipeline/D5/Minimap2/4_D01_2/mapped_FLNC_noScaff.sam) \
     > PacBio_HepG2_pre-TC.sam

time python ../../extract_SJs_from_sam.py \
    --sam PacBio_HepG2_pre-TC.sam \
    --genome ../../../../refs/hg38/hg38.fa \
    --o PacBio_HepG2_pre-TC

```

2) Obtain ONT HepG2 splice junctions in pre-TranscriptClean reads
```
cat ../../../../pipeline/ONT32/Minimap2/data/ONT32_mapped.sam \
    ../../../../pipeline/ONT33/Minimap2/data/ONT33_mapped.sam \
    | grep -v "^@" > ONT_HepG2_pre-TC.sam


time python ../../extract_SJs_from_sam.py \
    --sam ONT_HepG2_pre-TC.sam \
    --genome ../../../../refs/hg38/hg38.fa \
    --o ONT_HepG2_pre-TC
```

3) Use the splice junctions we extracted for Illumina using STAR:
```
/share/crsp/lab/seyedam/share/TALON_paper_data/illumina_sjs/HepG2_alignedSJ.out.tab
```

## Analysis

### Expression level


### Reproducibility

1) Extract only noncanonical SJs from each file:
```
awk '{if($5 == 0) print $0}' PacBio_HepG2_pre-TC_SJs.txt > ncsj_PacBio_HepG2_pre-TC_SJs.txt
awk '{if($5 == 0) print $0}' ONT_HepG2_pre-TC_SJs.txt > ncsj_ONT_HepG2_pre-TC_SJs.txt
awk '{if($5 == 0) print $0}' /share/crsp/lab/seyedam/share/TALON_paper_data/illumina_sjs/HepG2_alignedSJ.out.tab > ncsj_Illumina_HepG2_pre-TC_SJs.txt
```

2) Run comparison across platforms
```
mkdir figures
python ../../compare_sjs_venn.py \
        -pb ncsj_PacBio_HepG2_pre-TC_SJs.txt \
	-ont ncsj_ONT_HepG2_pre-TC_SJs.txt \
	-illumina ncsj_Illumina_HepG2_pre-TC_SJs.txt \
	-sample "HepG2 Noncanonical" \
        --log
```



