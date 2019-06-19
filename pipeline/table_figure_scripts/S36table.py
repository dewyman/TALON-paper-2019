#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Generate table with genes with higher novelty counts for cortex and hippocampus, also generates a plot to show how many of them are unique or shared between tissues

"""

import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from TALONClass import talonResults, writeOutfile

if len(sys.argv) < 2:
	sys.exit('python %s TALONabundancefile' % sys.argv[0])

talonfile = sys.argv[1]
outfile = 'noveltyReadCountVenn.tab'

cortexDatasets = ['Cortex_rep1','Cortex_rep2']
hippocampusDatasets = ['Hippocampus_rep1','Hippocampus_rep2']

cxResults = talonResults(talonfile, cortexDatasets)
hcResults = talonResults(talonfile, hippocampusDatasets)

cxGeneDict = cxResults.getGenes()
hcGeneDict = hcResults.getGenes()

cxOutGeneSet = set()
hcOutGeneSet = set()

for currentGID in cxGeneDict:
    thisGene = cxGeneDict[currentGID]
    KNOWN = 0
    NIC = 0
    NNC = 0

    for aTranscript in thisGene.getTranscripts('Known'):
        for dataset in cortexDatasets:
            KNOWN += aTranscript.getCounts(dataset)    
    for aTranscript in thisGene.getTranscripts('NIC'):
        for dataset in cortexDatasets:
            NIC += aTranscript.getCounts(dataset)    
    for aTranscript in thisGene.getTranscripts('NNC'):
        for dataset in cortexDatasets:
            NNC += aTranscript.getCounts(dataset)    

    if KNOWN < NIC + NNC:
        cxOutGeneSet.add(currentGID)       

for currentGID in hcGeneDict:
    thisGene = hcGeneDict[currentGID]
    KNOWN = 0
    NIC = 0
    NNC = 0

    for aTranscript in thisGene.getTranscripts('Known'):
        for dataset in hippocampusDatasets:
            KNOWN += aTranscript.getCounts(dataset)    
    for aTranscript in thisGene.getTranscripts('NIC'):
        for dataset in hippocampusDatasets:
            NIC += aTranscript.getCounts(dataset)    
    for aTranscript in thisGene.getTranscripts('NNC'):
        for dataset in hippocampusDatasets:
            NNC += aTranscript.getCounts(dataset)    

    if KNOWN < NIC + NNC:
        hcOutGeneSet.add(currentGID)   

vfig = venn2([cxOutGeneSet, hcOutGeneSet], set_labels = ('Cortex', 'Hippocampus'))
vfig.get_patch_by_id('10').set_color('pink')
vfig.get_patch_by_id('01').set_color('paleturquoise')
plt.title('Genes with higher NNC + NIC count than Known')
fig1 = plt.gcf()
fig1.savefig('cxhcNoveltyVenn.png')
plt.show()


#unionGeneSet = set.union(cxOutGeneSet,hcOutGeneSet)
commonGeneSet = set.intersection(cxOutGeneSet, hcOutGeneSet)
cxOnlyGeneSet = cxOutGeneSet.difference(commonGeneSet)
hcOnlyGeneSet = hcOutGeneSet.difference(commonGeneSet)

headerLine = 'group\tTALONGID\tknownGID\tgeneName\n'
outlineList = []

for (status, geneSet) in [('common',commonGeneSet), ('cortex', cxOnlyGeneSet), ('hippocampus', hcOnlyGeneSet)]:
    for currentGID in geneSet:
        thisGene = cxGeneDict[currentGID]
        outline = '%s\t%s\t%s\t%s\n' % (status, currentGID,thisGene.geneID,thisGene.geneAnnot)
        outlineList.append(outline)

writeOutfile(outfile, headerLine, outlineList)





