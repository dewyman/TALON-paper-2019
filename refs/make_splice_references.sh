# Time consuming- run in screen

set -e
mkdir -p TranscriptClean
cd TranscriptClean

# Paths specific to Dana's system. Be sure to use Python 2.7
TC_dir=~/TranscriptClean-1.0.7
module load dwyman/anaconda/3
source activate runTC

# For human
python $TC_dir/accessory_scripts/get_SJs_from_gtf.py \
      --f ../gencode.v29.annotation.gtf \
      --g ../hg38/hg38.fa \
      --o gencode_v29_SJs.tsv

# For mouse
python $TC_dir/accessory_scripts/get_SJs_from_gtf.py \
      --f ../gencode.vM20.annotation.gtf \
      --g ../mm10/mm10.fa \
      --o gencode_vM20_SJs.tsv
