# TALON paper repository
This repository contains custom scripts and documentation associated with the analyses in our BioRxiv manuscript concerning the application of our TALON pipeline to long read transcriptomes from PacBio and direct-RNA Oxford Nanopore. 

You can find the preprint here: https://www.biorxiv.org/content/10.1101/672931v1

To download the TALON program, please visit https://github.com/dewyman/TALON.  
To download the ENCODE DCC deployment of the TALON pipeline, please visit https://github.com/ENCODE-DCC/long-read-rna-pipeline. 

## plotting_scripts
Final versions of data visualization scripts used in the paper. 

## Figure_2
Describes exactly how the panels of Figure 2 in the paper were generated

## Figure_4
Describes exactly how the panels of Figure 4 in the paper were generated

## refs
Contains instructions for downloading reference genomes, GENCODE annotations, and variant files. Also includes scripts we used to prepare input files for TranscriptClean and TALON initializations.

## Illumina
Scripts for downloading short-read Illumina RNA-seq data from the ENCODE consortium and running Kallisto on them in order to perform short read quantification.

## analysis_scripts
Mostly visualization scripts for plotting things like gene/transcript expression correlations etc. We are in the process of deprecating this directory and moving final script versions to plotting_scripts instead.

## pipeline
Scripts and accessory files used to analyze and generate figures in our data.
