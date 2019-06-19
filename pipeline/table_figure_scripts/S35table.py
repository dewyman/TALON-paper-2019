"""
   Generates a table with gene read counts based on novelty category  """


import sys

from TALONClass import Transcript, Gene, talonResults, writeOutfile

# check that there is enough command line parameters
# sys.argv[0] is the name of the python script
# sys.argv[1] is the TALON abundance file
# sys.argv[2] is a prefix for the several output files
# sys.argv[3] is the comma-separated list of datasets to analyze 
if len(sys.argv) < 4:
	sys.exit('python %s infile outprefix datasetlist' % sys.argv[0])

talonfile = sys.argv[1]
outprefix = sys.argv[2]
datasets = sys.argv[3].split(',')


# get the TALON results as one object
myData = talonResults(talonfile, datasets)
geneDict = myData.getGenes()


# we build the filename using the prefix
outfilename4 = outprefix + '.noveltyGeneReadCount.tab'
# when we open a file with the 'w' flag, we are deleting that old file and writing over it
outfile4 = open(outfilename4,'w')


#outfile4
# we build the header line of column names in the output file
# first the standard gene info, then the novelty type columns
headerLine = 'TALONGID\tknownGID\tgeneName\tKnown\tNIC\tNNC'
outfile4.write(headerLine + '\n')

# for each gene in geneInfoDict, get the basic gene info,
# then calculate how many models of each novelty type we have

outlineList=[] 

for currentGID in geneDict:
	thisGene = geneDict[currentGID]
	outline = '%s\t%s\t%s' % (currentGID,thisGene.geneID,thisGene.geneAnnot)        
        # 
	KNOWN = 0
	NIC = 0
	NNC = 0

	for aTranscript in thisGene.getTranscripts('Known'):
		for dataset in datasets:
			KNOWN += aTranscript.getCounts(dataset)	
	for aTranscript in thisGene.getTranscripts('NIC'):
		for dataset in datasets:
			NIC += aTranscript.getCounts(dataset)	
	for aTranscript in thisGene.getTranscripts('NNC'):
		for dataset in datasets:
			NNC += aTranscript.getCounts(dataset)	

	if KNOWN < NIC + NNC:
		outline += '\t%d\t%d\t%d\n' % (KNOWN,NIC,NNC)
		Novelsum = NNC + NIC
		Novelratio = KNOWN / (NNC + NIC)
		outlineList.append((Novelratio,-1*Novelsum,outline))         
		
#outlineList.sort(reverse=True)

outlineList.sort()
for (ratio,count,line) in outlineList:
	# write each line
        outfile4.write(line)
# as a good citizen, we close the output file
outfile4.close()
