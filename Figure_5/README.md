# Figure 5: PacBio transcriptomes of 6-month adult male mouse cortex and hippocampus

Files/paths used to generate the panels of this figure:
```bash
PLOTPATH=../plotting_scripts/
APATH=../analysis_scripts/
# download the supplementary tables and change this path!
sup_tables=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_10-19/human_TALON/analysis/supplementary_tables/

gtf=${sup_tables}S34_mouse_brain_talon_observedOnly.gtf
```
Abundance and GTF files are available as supplementary tables of the TALON paper. 

Software versions:
* R v3.5.1

## Panel A: Isoform diversity for one representative dataset of cortex 
```bash
# TODO
```
<!-- <img align="center" width="400" src="TODO"> -->

 Pearson and Spearman correlations are recorded in ONT_GM12878_1-ONT_GM12878_2_gene_correlations.txt.

## Panel B: Isoform diversity for one representative dataset of hippocampus
```bash
# TODO
```
<!-- <img align="center" width="400" src="TODO"> -->

## Panel C: #TODO of genes with higher novelty read counts (NIC+NNC) than known, of which #TODO are only higher in cortex and #TODO higher in hippocampus
```bash
# TODO
```

<!-- <img align="center" width="400" src="TODO"> -->

## Panel D: GO semantic map for cortex genes
```bash
# TODO
```
<!-- <img align="center" width="400" src="TODO"> -->

## Panel E: GO semantic map for hippocampus genes
```bash
# TODO 
```
<!-- <img align="center" width="400" src="TODO"> -->

## Panel F: Example of Mef2a isoform expression in cortex and hippocampus
```bash
# create config file for gtf creation
# replace url with the url to your public-facing directory
url=http://crick.bio.uci.edu/freese/TALON_gtf/S34_mouse_brain_talon_observedOnly_tracks 
printf "${gtf},n+,0,none,$url" > mouse_brain_track_config
python ${APATH}gen_novelty_tracks_gtf.py \
    --c mouse_brain_track_config
mv ${sup_tables}S34_mouse_brain_talon_observedOnly_tracks/ .
```
Then copy your tracks directory to your public-facing directory and load the trackline in the genome browser, and screenshot.
<!-- <img align="center" width="400" src="TODO"> -->