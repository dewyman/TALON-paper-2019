
# Copy over all data
WORKDIR=/pub/dwyman/TALON-paper-2019/pipeline/D4
Minimap2=${WORKDIR}/Minimap2
mkdir -p Minimap2
cd $Minimap2

mkdir -p $Minimap2/1_A01
mkdir -p $Minimap2/2_B01
cp /share/samdata/dwyman/PB69/Minimap2/1_A01/mapped_FLNC.sam $Minimap2/1_A01/.
cp /share/samdata/dwyman/PB69/Minimap2/2_B01/mapped_FLNC.sam $Minimap2/2_B01/.

# Run TC
/share/samdata/dwyman/sequel_scripts/./make-TC-scripts.sh $WORKDIR

# Filter reads post-TC to keep reads that are canonical or annotated noncanonical
/share/samdata/dwyman/sequel_scripts/./filter-TC-output.sh $WORKDIR

# Run TALON
