import sys

# this is the "class" code for creating Transcript objects. 
# all objects belong to a class, and this is the definition of it
# there will be 1 object for each Talon Transcript
# each column name is added to the object as an attribute
# we could have written this in a separate file, then imported it using import
class Transcript:
    # this is a list of column names shared between between all objects
    fieldNames = []
    
    # this is the function that is called when we create an object
    def __init__(self, inline):
        #we read the line and split it into a list of fields. .strip() removes the newline
        fields = inline.strip().split()
        self.transcript_ID = ''
        # we go through the columns in order using a numeric index
        for index in range(len(self.fieldNames)):
            attribute = self.fieldNames[index]
            value = fields[index]
            # this is the command that adds the attribute to the object for each field
            setattr(self, attribute, value)
            # This will change with the new TALON rewrite
            if attribute == 'transcript_novelty':
                self.transcriptType = value
    def __len__(self):
        return int(self.length)
    def getCounts(self,dataset):
        return int(getattr(self,dataset))
    def __str__(self):
        return "instance of Gene for " + self.transcript_ID + "\n" + str(self.__dict__)
    def __repr__(self):
        return "Talon transcript " + str(self.transcript_ID)  

class Gene:
    def __init__(self,gID,gAnnot,gTalonID,gStatus):
        self.geneID = gID
        self.status = gStatus
        self.geneAnnot = gAnnot
        self.geneTalonID = gTalonID
        self.transcriptDict = {}
    def setTypes(self,typeList):
        for aType in typeList:
            self.transcriptDict[aType]=[]
    def addTranscript(self,aTranscript):
        transcriptType    = aTranscript.transcriptType
        self.transcriptDict[transcriptType].append(aTranscript)
    def getTranscripts(self,aType=''):
        if aType in self.transcriptDict:
            return self.transcriptDict[aType]
        results=[]
        for oneType in self.transcriptDict:
            results += self.transcriptDict[oneType]
        return results
    def __str__(self):
        return "instance of Gene for " + self.geneID + "\n" + str(self.__dict__) 
        
class talonResults:
    # These are the different type of transcripts in the order that they will be recorded in the output files
    transcriptTypes = ['Known','NIC','NNC','ISM','Antisense','Intergenic','Genomic']
    def __init__(self, abundanceFile, datasetList):
        self.infilename = abundanceFile        
        self.geneDict = {}
        self.nameDict={}
        infile = open(abundanceFile)
        # the first line of the TALON abundance file has all of the column names
        firstLine = infile.readline()
        # this is where we define the column names for all of the Transcript objects
        Transcript.fieldNames = firstLine.strip().split()
        # We are going to read each line of the input file, create a transcript object, and record it in the dictionaries
        for line in infile:
            # create an object for this line
            thisTranscript = Transcript(line)
            # this is the object's TALON gene_ID
            currentGID = thisTranscript.gene_ID
            # initialize geneDict for the gene using a new Gene object if necessary
            if currentGID not in self.geneDict:
                self.geneDict[currentGID] = Gene( thisTranscript.annot_gene_id, thisTranscript.annot_gene_name,currentGID, thisTranscript.gene_novelty)
                #self.nameDict[currentGID]=thisTranscript.geneAnnot
                self.nameDict[thisTranscript.annot_gene_name.lower()] = currentGID
                self.geneDict[currentGID].setTypes(talonResults.transcriptTypes)
            # we now add the transcriptObject in the geneDict, which will track also its novelty type
            self.geneDict[currentGID].addTranscript(thisTranscript)
        infile.close()
    def __str__(self):
        return "instance of talonResults for " + self.infilename + " with " + str(len(self.geneDict)) + " genes" 
    def getGenes(self):
        return self.geneDict
    def getTranscriptTypes(self):
        return self.transcriptTypes
    def getGene(self,geneAnnot):
        thisGID = self.nameDict[geneAnnot.lower()]
        return self.geneDict[thisGID]

# just a helper function for writing a file given a name, a header line, and a list of lines
def  writeOutfile(outfilename, header, lineList):
    outfile = open(outfilename,'w')
    outfile.write(header)
    outfile.writelines(lineList)
    outfile.close()
