# Analysis of noncanonical splice junctions in HepG2

The goal is to look at the splice junctions in HepG2 long reads before TranscriptClean and examine:
1) Their reproducibility or lack thereof across the PacBio, ONT, and Illumina platforms.

## Data processing

1) Obtain PacBio HepG2 splice junctions (rep1 and 2) in pre-TranscriptClean reads. Note: this is slow when you have a lot of reads.
```
data_dir=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_10-19/data
cat ${data_dir}/PacBio_HepG2_1/Minimap2/all_mapped_FLNC_noScaff.sam \
    ${data_dir}/PacBio_HepG2_2/Minimap2/all_mapped_FLNC_noScaff.sam \
    | grep -v "^@" > PacBio_HepG2_pre-TC.sam

time python ../../extract_SJs_from_sam.py \
    --sam PacBio_HepG2_pre-TC.sam \
    --genome ../../../refs/hg38/hg38.fa \
    --o PacBio_HepG2_pre-TC

```

2) Obtain ONT HepG2 splice junctions in pre-TranscriptClean reads
```
data_dir=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_10-19/data
cat ${data_dir}/ONT_HepG2_1/Minimap2/sorted_all_mapped_noScaff.sam \
    ${data_dir}/ONT_HepG2_3/Minimap2/sorted_all_mapped_noScaff.sam \
    | grep -v "^@" > ONT_HepG2_pre-TC.sam


time python ../../extract_SJs_from_sam.py \
    --sam ONT_HepG2_pre-TC.sam \
    --genome ../../../refs/hg38/hg38.fa \
    --o ONT_HepG2_pre-TC
```

3) Use the splice junctions we extracted for Illumina using STAR:
```
/share/crsp/lab/seyedam/share/TALON_paper_data/illumina_sjs/HepG2_alignedSJ.out.tab
```

## Analysis

### Reproducibility

1) Extract only noncanonical SJs from each file:
```
awk '{if($5 == 0) print $0}' PacBio_HepG2_pre-TC_SJs.txt > ncsj_PacBio_HepG2_pre-TC_SJs.txt
awk '{if($5 == 0) print $0}' ONT_HepG2_pre-TC_SJs.txt > ncsj_ONT_HepG2_pre-TC_SJs.txt
awk '{if($5 == 0) print $0}' /share/crsp/lab/seyedam/share/TALON_paper_data/illumina_sjs/HepG2_alignedSJ.out.tab > ncsj_Illumina_HepG2_pre-TC_SJs.txt
```

2) Run comparison across platforms
```
source activate mypython3.7.2
mkdir -p figures
python ../../compare_sjs_venn.py \
        -pb ncsj_PacBio_HepG2_pre-TC_SJs.txt \
	-ont ncsj_ONT_HepG2_pre-TC_SJs.txt \
	-illumina ncsj_Illumina_HepG2_pre-TC_SJs.txt \
	-sample "HepG2 Noncanonical" \
        --log
source deactivate
```
<img align="center" width="500" src="figures/HepG2_Noncanonical_venn.png">
