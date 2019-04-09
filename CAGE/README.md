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
