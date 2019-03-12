# Reference files
* human and mouse GENCODE annotations (GTF) 
    * Human: GENCODE v29 comprehensive, primary chromosomes only
    * Mouse: GENCODE vM20 comprehensive, primary chromosomes only
* human and mouse reference genomes  
    * Human: GRCh38 (hg38). 
    * Mouse: mm10
* STAR index for human and mouse genomes  
    * Used for mapping short reads. 
    * Downloaded from ENCODE
* Splice junctions derived from the GTF files 
    * Use with TranscriptClean
* Set of common variants from dbSNP
    * Use with TranscriptClean for human data (we don't use SNPs for mouse)
* Initialized TALON databases for human and mouse
    * These are intended to be maintained in their original state- if you want to use one for a TALON run, be sure to make a copy and run on that, rather than the original copy.
