
# Copy over all data
WORKDIR=/pub/dwyman/TALON-paper-2019/pipeline/D6
Minimap2=${WORKDIR}/Minimap2
mkdir -p Minimap2
cd $Minimap2

mkdir -p $Minimap2/1_A01
mkdir -p $Minimap2/2_B01
mkdir -p $Minimap2/1_A01_2
mkdir -p $Minimap2/2_B01_2
#cp /share/samdata/dwyman/PB63/Minimap2/1_A01/mapped_FLNC.sam $Minimap2/1_A01/.
#cp /share/samdata/dwyman/PB63/Minimap2/2_B01/mapped_FLNC.sam $Minimap2/2_B01/.
cp /share/samdata/dwyman/PB63_run2/Minimap2/1_A01/mapped_FLNC.sam $Minimap2/1_A01_2/.
cp /share/samdata/dwyman/PB63_run2/Minimap2/2_B01/mapped_FLNC.sam $Minimap2/2_B01_2/.

# Run TC
#/share/samdata/dwyman/sequel_scripts/./make-TC-scripts.sh $WORKDIR

# Filter reads post-TC to keep reads that are canonical or annotated noncanonical
#/share/samdata/dwyman/sequel_scripts/./filter-TC-output.sh $WORKDIR

# Run TALON
#/share/samdata/dwyman/sequel_scripts/./make-TALON-human-script.sh $WORKDIR GM12878 PacBio-Sequel
