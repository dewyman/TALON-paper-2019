
# Copy over all data
WORKDIR=/pub/dwyman/TALON-paper-2019/pipeline/D5
Minimap2=${WORKDIR}/Minimap2
mkdir -p Minimap2
cd $Minimap2

mkdir -p $Minimap2/3_C01
mkdir -p $Minimap2/4_D01
mkdir -p $Minimap2/3_C01_2
mkdir -p $Minimap2/4_D01_2

cp /share/samdata/dwyman/PB70/Minimap2/3_C01/mapped_FLNC.sam $Minimap2/3_C01/. 
cp /share/samdata/dwyman/PB70/Minimap2/4_D01/mapped_FLNC.sam $Minimap2/4_D01/.
cp /share/samdata/dwyman/PB70_run2/Minimap2/3_C01/mapped_FLNC.sam $Minimap2/3_C01_2/.
cp /share/samdata/dwyman/PB70_run2/Minimap2/4_D01/mapped_FLNC.sam $Minimap2/4_D01_2/.

# Run TC
/share/samdata/dwyman/sequel_scripts/./make-TC-scripts.sh $WORKDIR

# Filter reads post-TC to keep reads that are canonical or annotated noncanonical
#/share/samdata/dwyman/sequel_scripts/./filter-TC-output.sh $WORKDIR

# Run TALON
