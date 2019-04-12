# Liftover the RNA-PET data from hg19 to hg38
module load jje/kent
set -e

for f in data/FANTOM5/*.bed; do
    name=$(basename $f .bed)
    liftOver data/FANTOM5/${name}.bed ../RNA-PET/hg19ToHg38.over.chain data/FANTOM5/${name}_hg38.bed data/FANTOM5/${name}_unlifted.bed
    rm data/FANTOM5/${name}_unlifted.bed
done

