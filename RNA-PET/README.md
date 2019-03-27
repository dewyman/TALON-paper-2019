# RNA-PET
To externally validate our long read transcript starts and ends, we compare them to start-end pairs in the same cell lines derived from the RNA-PET assay. We chose the polyA-selected cytosol data, with clone-free library prep where possible (due to the longer read lengths). 

## Download the data from ENCODE
```
./download_RNA-PET.sh
# HepG2: clone
wget https://www.encodeproject.org/files/ENCFF001TIR/@@download/ENCFF001TIR.bed.gz

### GM12878: clone-free
wget https://www.encodeproject.org/files/ENCFF001TIL/@@download/ENCFF001TIL.bed.gz

# K562: clone-free
https://www.encodeproject.org/files/ENCFF001TJA/@@download/ENCFF001TJA.bed.gz
```

## Obtain the Liftover tool and chain file from the UCSC genome browser:
```

```

## Lift over each BED file that we downloaded
```

```
