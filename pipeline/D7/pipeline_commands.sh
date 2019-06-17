
# Copy over all data
WORKDIR=/pub/dwyman/TALON-paper-2019/pipeline/D7
Minimap2=${WORKDIR}/Minimap2
mkdir -p Minimap2
cd $Minimap2

mkdir -p $Minimap2/3_C01
mkdir -p $Minimap2/4_D01
mkdir -p $Minimap2/1_A01
mkdir -p $Minimap2/2_B01
#cp /share/samdata/dwyman/PB64/Minimap2/3_C01/mapped_FLNC.sam $Minimap2/3_C01/.
#cp /share/samdata/dwyman/PB64/Minimap2/4_D01/mapped_FLNC.sam $Minimap2/4_D01/.
#cp /share/samdata/dwyman/PB64_run2/Minimap2/1_A01/mapped_FLNC.sam $Minimap2/1_A01/.
#cp /share/samdata/dwyman/PB64_run2/Minimap2/2_B01/mapped_FLNC.sam $Minimap2/2_B01/.

# Run TC
#/share/samdata/dwyman/sequel_scripts/./make-TC-scripts.sh $WORKDIR

# Filter reads post-TC to keep reads that are canonical or annotated noncanonical
/share/samdata/dwyman/sequel_scripts/./filter-TC-output.sh $WORKDIR

# Run TALON
#/share/samdata/dwyman/sequel_scripts/./make-TALON-human-script.sh $WORKDIR GM12878 PacBio-Sequel
