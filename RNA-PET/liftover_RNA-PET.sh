# Liftover the RNA-PET data from hg19 to hg38
set -e

for f in data/*.bed; do
    name=$(basename $f .bed)
    ./liftOver data/${name}.bed hg19ToHg38.over.chain data/${name}_hg38.bed data/${name}_unlifted.bed
    rm data/${name}_unlifted.bed
done

