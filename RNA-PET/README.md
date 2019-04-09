# RNA-PET
To externally validate our long read transcript starts and ends, we compare them to start-end pairs in the same cell lines derived from the RNA-PET assay. We chose the polyA-selected cytosol data, with clone-free library prep where possible (due to the longer read lengths). 

## Download the data from ENCODE
```
./download_RNA-PET.sh
# HepG2: clone
# GM12878: clone-free
# K562: clone-free
```

## Download Liftover
```
# Program:
Linux: wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
Mac OSX: wget http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/liftOver

# Chain file:
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
```

## Lift over each BED file that we downloaded: 
```
./liftover_RNA-PET.sh
```

## After this, run script run_RNA-PET_analysis.py to run the comparison

## Finally, run plot_RNA-PET_support.R to visualize the results 
