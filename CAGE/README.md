# CAGE
To externally validate our long read transcript starts, we compare them to 5' end cap sites in the same cell lines derived from the CAGE assay. We start with data from ENCODE.

## ENCODE 
### Download the ENCODE data for GM12878, K562, and HepG2
```
./download_ENCODE_CAGE.sh
```
These are already mapped to hg38, so no liftover is necessary. 

### Run the pipeline
```
source activate mypython3.7.2
python run_CAGE_analysis.py --h
```

## FANTOM5
For FANTOM, we have one file of robust peaks measured across a variety of human samples rather than individual files per cell line. The peaks are mapped to hg19, so we need to use the UCSC genome browser liftover tool to convert them to hg38. I already downloaded this tool and the necessary chain file when working on the RNA-PET data.
### Download the CAGE data
```
./download_FANTOM5_CAGE.sh
```
### LiftOver to hg38
```
./liftover_FANTOM.sh
```
