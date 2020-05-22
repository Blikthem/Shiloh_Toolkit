#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 09:31:06 2018

@author: claypooldj
"""

####Load dependencies
import csv
import time

"""
###------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Input variables
###------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
inpGTF="/Volumes/Untitled/Genomes/fixingLongestMeta.gtf"
outDir="/Volumes/Untitled/Genomes/"
outName="testFixes"
byID=False
"""

###------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Describe objects
###------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
class GTFRow:
    """
    GTFRow Object: Data stroage object representing a row from the anotation (GTF) file.

    Properties:
        geneName(String) - The name of the gene this row maps to
        classification(String) - The type (5'/CDS/3' etc) of this row.
        start(int) - The starting index of this row
        stop(int) - The end index of this row
        chromosome(String) - The chromosome this row belongs to
        orientation(string) - The orientation of the row, either '+' or '-'
        transcriptID(String) - The identifier for the transcript associated with the row.
    
    Methods:
        init - Initializer
        str - Conversion to string form for printing
        getGeneName - Return gene name.
        swtichToPos 
    """
    def __init__(self,chromosome,geneName,start,stop,orientation,transcriptID,classification):
        self.chromosome=chromosome
        self.classification=classification
        self.geneName=geneName
        self.start=start
        self.stop=stop
        self.orientation=orientation
        self.transcriptID=transcriptID

    def __str__(self):
        toPrint = "{GTFRow Object} Start: "+str(self.start)+"\tStop: "+str(self.stop)+"\tGene: "+str(self.geneName)+"\tClassification: "+str(self.classification)+"\tOrientation: "+str(self.orientation)+"\tChromeosome: "+str(self.chromosome)
        return(toPrint)

    def getGeneName(self):
        newName=self.geneName
        return(newName)


class transcript:
    """
    transcript Object: A transcript (splice varient) option from the gtf file.

    Properties:
        geneName(String) - The name of the gene this row maps to
        classification(String) - The type (5'/CDS/3' etc) of this row.
        start(int) - The starting index of this row
        stop(int) - The end index of this row
        chromosome(String) - The chromosome this row belongs to
        orientation(string) - The orientation of the row, either '+' or '-'
        transcriptID(String) - The identifier for the transcript associated with the row.
    
    Methods:
        init - Initializer
        str - Conversion to string form for printing
        getGeneName - Return gene name.
        switchToPos 
    """
    def __init__(self,start,stop):
        #Int start value
        self.start=start
        #Int end value
        self.stop=stop
        #List of GTFRow objects
        self.tRows = []
        #List of merged elements
        self.merged=[]
        #Sequence as a string
        self.seqString=""
        
    def __str__(self):
        toPrint="{transcript object} Start: "+str(self.start)+"\tStop: "+str(self.stop)+"\tNumber of row objects: "+str(len(self.tRows))
        return(toPrint)
    
    def addRow(self,row):
        self.tRows.append(row)

    def createTransSequence(self):
        """
        Re-express the transcript as a transcript sequence (list of read elements).
        Saves - to object property.
            List[sequenceElement]
        """
        
        #Create the initial list which has many objects of four different types
        initialList =[]
        
        counter=0
        #For each row
        for row in self.tRows:
            counter=counter+1
            added=False
            curID = row.classification
            
            #If it is a start codon
            if curID =="start_codon":
                #Make start codon object
                if (row.orientation=="+"):
                    nse = sequenceElement(row.start,row.stop,"Start_Codon +")
                if (row.orientation=="-"):
                    nse = sequenceElement(row.start,row.stop,"Start_Codon -")
                initialList.append(nse)
                added=True
                              
                    
            #If it is a stop codon
            if curID =="stop_codon":
                #Make stop codon object
                if (row.orientation=="+"):
                    nse = sequenceElement(row.start,row.stop,"Stop_Codon +")
                if (row.orientation=="-"):
                    nse = sequenceElement(row.start,row.stop,"Stop_Codon -")
                initialList.append(nse)
                added=True           
                
            #If it is a CDS
            if curID =="CDS":
                #If it is after a start
                nse = sequenceElement(row.start,row.stop,"CDS")
                initialList.append(nse)
                added=True                  
                            
            #If it is an exon
            if curID =="exon":
                #If the next item in the list is a CDS we don't want to a dd a utr
                if counter != len(self.tRows):
                    nxtEle = self.tRows[counter]
                    if (nxtEle.classification=="CDS"):    
                        added=True
            
            #Otherwise it is a noncoding region
            if (added==False):
                #Add as a UTR
                nse =sequenceElement(row.start,row.stop,"UTR")
                initialList.append(nse)

          #Create a new list that is merged
        mergedLst = mergeUTR(initialList)
        self.merged=mergedLst              
        
        
class sequenceElement:
    """
    sequenceElemnt Object: Data containner for an element within a transcript sequence.
    
    Properties:
        start(String) - Starting location.
        stop(String) - Stopping location.
        ID(String) - Identifying string for this element.
    """
    def __init__(self,start,stop,ID):
        self.start=start
        self.stop=stop
        self.ID=ID
        
    def __str__(self):
        toPrint = "{Sequence element object} ID "+str(self.ID)+"\tStart: "+str(self.start)+"\tStop: "+str(self.stop)
        return(toPrint)
        
###------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Gene object
class gtfGene:
    def __init__(self,name):
        self.name = name
        self.lstOfRows = []
        self.lstOfTranscripts = []
        #final stranscript option
        self.metaTranscript = transcript(0,1)
        #Are there any transcript objects in this object? Useful for checking if we need to generate our own.
        self.hasTranscriptLabels=False
        
    #Print statement
    def __str__(self):
        toPrint = "{gtfGene Object} Gene: "+self.name+"\tNumber of gtfRow objects included: "+str(len(self.lstOfRows))+"\tNumber of transcript objects included: "+str(len(self.lstOfTranscripts))
        return(toPrint)
            
    #Method allowing addition of GTFRow objects  
    def addRow(self,toAdd):
        self.lstOfRows.append(toAdd)

    def populateTranscriptLst(self):
        """
        This function populates the list of transcript objects based on the contents of the gene object.
        """
        #For each gtf row object in this gene
        checker=0
        for row in self.lstOfRows:
            #Get the ID of that item
            curID = row.classification
            #If this is "transcript"
            if (curID == "transcript"):
                #Create a new transcript object and add it to the list
                newTrans = transcript(row.start,row.stop)
                #Add the current gtfRow object to this list
                self.lstOfTranscripts.append(newTrans)
                #If this is NOT a new transcript
                checker=1
            if (curID!="transcript"and checker==1):
                #Add it to the current gtftranscript object
                    #The current
                toAddIndex = len(self.lstOfTranscripts)-1
                curTrans = self.lstOfTranscripts[toAddIndex]
                curTrans.addRow(row)
            
            
    def createTranscriptRows(self):
        """
        This function is for gtf files which do not contain explicit rows detailing the transcript dimensions (start/stop) for each gene.
        """
        if self.hasTranscriptLabels==True:
            return()
            
        #Loop over the rows
        curTrans=self.lstOfRows[0].transcriptID
        curStart=self.lstOfRows[0].start
        curStop=self.lstOfRows[0].stop
        startVal=0
        for i in range(1,len(self.lstOfRows)):
            cRow=self.lstOfRows[i]
            #If it is the same as the previous row
            if cRow.transcriptID==curTrans:
                if cRow.start<curStart:
                    curStart=cRow.start
                if cRow.stop>curStop:
                    curStop=cRow.stop
            if cRow.transcriptID!=curTrans or i==len(self.lstOfRows)-1:
                #Make a new row object for the row
                tRow=GTFRow(cRow.chromosome,cRow.geneName,curStart,curStop,cRow.orientation,cRow.transcriptID,"transcript")
                #Finally add this row at the start of its relevent 
                self.lstOfRows.insert(startVal,tRow)
                #And reset counters
                curTrans=cRow.transcriptID
                curStart=cRow.start
                curStop=cRow.stop
                startVal=i
    
    def getLongestTrans(self):
        """
        Gets the longest transcript object associated with a gene object
        Output:
            transcript - The longest transcript stored within this object.
        """
        lstOfLengs=[]
        lstOfCDS=[]
        if (len(self.lstOfTranscripts)!=0):
            for trans in self.lstOfTranscripts:
                #Get the sum of the exon lengths for each transcript
                counter=0
                containsCDS=False
                for ro in trans.tRows:
                    if ro.classification=="exon":
                        tAdd = int(ro.stop)-int(ro.start)
                        counter=counter+tAdd
                    if (ro.classification=="CDS"):
                        containsCDS=True
                lstOfLengs.append(tAdd)
                lstOfCDS.append(containsCDS)
        
            #Check if there is a CDS option
            isCDS=False
            for ele in lstOfCDS:
                if ele == True:
                    isCDS=True
        
            ##If there is no CDS
            if isCDS==False:
                #Select the max
                theMax = max(lstOfLengs)
                aIndex = lstOfLengs.index(theMax)
               # print(self.lstOfTranscripts[aIndex])
                return(self.lstOfTranscripts[aIndex])
            
            #If there IS a cds we need to remove all of those transcripts which do not contain a cds
            newLstOfLengs=[]
            if isCDS==True:
                for i in range(0,len(lstOfCDS)):
                    ch = lstOfCDS[i]
                    if ch==True:
                        newLstOfLengs.append(lstOfLengs[i])
        
        
            #Select the longest transcript
            theMax = max(newLstOfLengs)
            aIndex = newLstOfLengs.index(theMax)
            return(self.lstOfTranscripts[aIndex])
          
        #Otherwise return an empty transcripts    
        else:
              return(transcript(0,1))
            
    def setRowsAsLongestTrans(self):
        """
        Sets all of the rows of the gene to instead be equal to the rows in the longest transcript
        """
        theTrans = self.getLongestTrans()
        self.lstOfRows = theTrans.tRows
    
    def fixGeneNames(self):
        """
        Loops through the rows of a given gene and fixes the gene names of them.
        """
        self.name=withinQuotes(self.name)
        for ro in self.lstOfRows:
            ro.geneName=withinQuotes(ro.geneName)
        
    def getStartStart(self):
        """
        Returns the start nucleotide position of the start codon of the gene
        Outputs:
            integer
        """
        for ro in self.lstOfRows:
            if (ro.classification=="start_codon"):
                return(int(ro.start))
    
    #Returns the start nucleotide of the stop codon of the gene
    def getStopStart(self):
        for ro in self.lstOfRows:
            if (ro.classification=="stop_codon"):
                return(int(ro.start))
                
    def sortRowsByStart(self):
        """
        Sorts the row objects of this gene by their startign value
        """
        self.lstOfRows = sortByStart(self.lstOfRows)

    def firstFUTR(self,testRo):
        """
        Checks if a test ro is the first FUTR.
        Input:
            testRo (AnotRow)
        Returns:
            bool
        """
        #If the orientation is positive
        if "+" in testRo.orientation:
            for ro in self.lstOfRows:
                if ro.classification=="FUTR":
                    if int(ro.start)<int(testRo.start):
                        return(False)
            return(True)
            
        #If the orientation is negative
        if "-" in testRo.orientation:
            for ro in self.lstOfRows:
                if ro.classification=="FUTR":
                    if int(ro.start)>int(testRo.start):
                        return(False)
            return(True) 
            
    def lastTUTR(self,testRo):
        """
        Checks if a test ro is the first FUTR.
        Input:
            testRo (AnotRow)
        Returns:
            bool
        """
        #If the orientation is positive
        if "+" in testRo.orientation:
            for ro in self.lstOfRows:
                if ro.classification=="TUTR":
                    if int(ro.start)>int(testRo.start):
                        return(False)
            return(True)
            
        #If the orientation is negative
        if "-" in testRo.orientation:
            for ro in self.lstOfRows:
                if ro.classification=="TUTR":
                    if int(ro.start)<int(testRo.start):
                        return(False)
            return(True)    

    def setUTR(self):
        """
        Changes the names of the rows in the three prime utr (TPU) and five prime utr (FPU) regions of the sequence
        """
        #Get start and stop
        theStartStart = self.getStartStart()
        theStopStart = self.getStopStart()
        if (type(theStartStart)!=int) or type(theStopStart)!=int:
            return(False)
        
        exRo = self.lstOfRows[0]
        for ro in self.lstOfRows:
            #IF we have a row identified as UTR without the context (up or down stream)
            if (ro.classification=="UTR"):
                #If positive
                if "+" in exRo.orientation:
                    #Check if before the start
                    if (int(ro.start)<theStartStart):
                        ro.classification="FUTR"
                    if (int(ro.start)>=theStopStart):
                        ro.classification="TUTR"
                #If negative
                elif ("-" in exRo.orientation):
                    if (int(ro.start)<=theStopStart):
                        ro.classification="TUTR"
                    if (int(ro.start)>theStartStart):
                        ro.classification="FUTR"
           
        #Now select only the first FUTR and last TUTR, removing the others
        toKeep=[]
        for ro in self.lstOfRows:
            if ro.classification=="TUTR":
                if self.lastTUTR(ro)==True:
                    if ro not in toKeep:
                        toKeep.append(ro)
                else:
                        ro.classification="toRem"
            if ro.classification=="FUTR":
                if self.firstFUTR(ro)==True:
                    if ro not in toKeep:
                        toKeep.append(ro)
                else:
                        ro.classification="toRem"
            else:
                if ro not in toKeep:
                    toKeep.append(ro)
        #self.lstOfRows=toKeep


    def setUTR2(self):
        """
        Changes the names of the rows in the three prime utr (TPU) and five prime utr (FPU) regions of the sequence.
        """
        if len(self.lstOfRows)==0:
            return()
        exRo = self.lstOfRows[0]
        
        if "+" in exRo.orientation:
            for ro in self.lstOfRows:
                if ro.classification=="UTR":
                    ro.classification="FUTR"
                    break
                if ro.classification=="intron":
                    break
            for ro in reversed(self.lstOfRows):
                if ro.classification=="UTR":
                    ro.classification="TUTR"
                    break
                if ro.classification=="intron":
                    break               
        
        if "-" in exRo.orientation:
            for ro in reversed(self.lstOfRows):
                if ro.classification=="UTR":
                    ro.classification="FUTR"
                    break
                if ro.classification=="intron":
                    break                
            for ro in self.lstOfRows:
                if ro.classification=="UTR":
                    ro.classification="TUTR"
                    break
                if ro.classification=="intron":
                    break                   
            
            
    def removeStart(self):
        """
        Converts the start codons into coding sequences.
        """
        i=-1
        for ro in self.lstOfRows:
            i=i+1
            ro = self.lstOfRows[i]
            if (i!=0):
                prev= self.lstOfRows[i-1]
            else:
                prev=ro
            if (i!=len(self.lstOfRows)-1):
                nxt = self.lstOfRows[i+1]
            else:
                nxt=ro
                
            if (ro.classification=="start_codon"):
                ro.classification="CDS"               
                #Remove if need be
                if prev!=ro:
                    if prev.classification=="CDS":
                          self.lstOfRows.remove(ro)
                elif nxt!=ro:
                    if nxt.classification=="CDS":
                          self.lstOfRows.remove(ro)

    def removeStop(self): 
        """
        Converts the stop codons into CDS.
        """
        i=-1
        for ro in self.lstOfRows:
            i=i+1
            ro = self.lstOfRows[i]
            if (i!=0):
                prev= self.lstOfRows[i-1]
            else:
                prev=ro
            if (i!=len(self.lstOfRows)-1):
                nxt = self.lstOfRows[i+1]
            else:
                nxt=ro
                
            if (ro.classification=="stop_codon"):
                ro.classification="CDS"               
                #Remove if need be
                if prev!=ro:
                    if prev.classification=="CDS":
                          self.lstOfRows.remove(ro)
                elif nxt!=ro:
                    if nxt.classification=="CDS":
                          self.lstOfRows.remove(ro)
                    
    def addIntrons(self):
        """
        Adds the intronic sequences where needed.
        """
        lstOfIntrons = []
        counter=1
        for ro in self.lstOfRows:            
            if ro.classification=="exon":
                #Check if there is another exon downstream
                nxtRo = self.getNextExon(counter)
                if (nxtRo!=False):
                    #CREATE THE INTRON AND ADD IT TO THE LIST TO ADD
                    theStart = int(ro.stop)+1
                    theStop = int(nxtRo.start)-1
                    nIntron = GTFRow(ro.chromosome,ro.geneName,theStart,theStop,ro.orientation,ro.transcriptID,"intron")

                    lstOfIntrons.append(nIntron)
                    
            counter=counter+1
        
        #Add the introns to the gtfGene
        for intron in lstOfIntrons:
            self.addRow(intron)
        #Sort the gene
        self.sortRowsByStart
    
    def getNextExon(self,index):
        """
        Gets the next exon in a transcript after a given index or it returns false is none occurs after the index.
        Input:
            index (int) - The index where the scanning begins.
        Output:
            bool or sequenceElement
        """
        if index>len(self.lstOfRows):
            return(False)
        #Loop over each element after the index
        for i in range(index,len(self.lstOfRows)):
            ro = self.lstOfRows[i]
            if ro.classification=="intron":
                return(False)
            if ro.classification=="exon":
                return(ro)
        
        return(False)

    def mergeAllOverlaps(self):
        """
        Merge overlapping elements.
        """
        #For each element
        for ro in self.lstOfRows:
            #Compare it with each other element
            for ro2 in self.lstOfRows:
                if shouldMerge(ro,ro2)==True and ro!=ro2:
                    #Create the new row
                    nRow = mergeRowElements(ro,ro2)
                    #Add this row
                    self.lstOfRows.append(nRow)
                    #Remove the old rows
                    self.lstOfRows.remove(ro)
                    self.lstOfRows.remove(ro2)
                    return(True)
    
    def mergeAdjacent(self):
        """
        Merge adjacent regions of the same classification.
        """
        for ro in self.lstOfRows:
            for ro2 in self.lstOfRows:
                if int(ro.stop)+1==int(ro2.start) and ro.classification==ro2.classification:
                    nRow = GTFRow(ro.chromosome,ro.geneName,ro.start,ro2.stop,ro.orientation,ro.transcriptID,ro.classification)
                    #Add this row
                    self.lstOfRows.append(nRow)
                    #Remove the old rows
                    self.lstOfRows.remove(ro)
                    self.lstOfRows.remove(ro2)
                    return(True)

    def removeDegenerateExons(self):
        """
        Removes all those exons that have the same coordinates as UTR or CDS.
        """
        i=0
        for ro in self.lstOfRows:
            if ro.start==ro.stop:
                self.lstOfRows.remove(ro)
            if ro.classification=="exon":
                for ro2 in self.lstOfRows:
                    if ro!=ro2:
                        if ro2.classification=="UTR" or ro2.classification=="CDS" or ro2.classification=="FUTR" or ro2.classification=="TUTR":
                            if ro.start==ro2.start:
                                if ro.stop==ro2.stop:
                                    self.lstOfRows.remove(ro)
                                    break
                if i<=(len(self.lstOfRows)-2) and i!=0 and ro in self.lstOfRows:
                    allowed=["UTR","CDS","FUTR","TUTR"]
                    prevRo=self.lstOfRows[i-1]
                    nextRo=self.lstOfRows[i+1]
                    if prevRo.classification in allowed and nextRo.classification in allowed:
                        if int(prevRo.start)==int(ro.start) and int(nextRo.stop)==int(ro.stop):
                            self.lstOfRows.remove(ro)
            i=i+1
        
    def finalCleanUp(self,introns):
        """
        Input:
            introns(bool) - Populate fields for implied introns? Only necessary for the intron/exon splice junction itermediate file.
        """
        self.sortRowsByStart()
        #self.setUTR()
        #self.fixGeneNames()
        self.removeStart()
        self.removeStop()
        for w in range(0,len(self.lstOfRows)):
            self.mergeAllOverlaps()
            self.mergeAdjacent()
            self.sortRowsByStart()
            if introns==True:
                self.addIntrons()
        #self.removeDegenerateExons()
        self.setUTR2()
###------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Functions
###------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------               
def shouldMerge(row1,row2):
    """
    Checks whether two rows overlap and thus should be merged into a single row in the output annotation file.
    Inputs:
        row1 and row2(gtfRow)
    Returns:
        bool
    """
    #Checks the type
    if (row1.classification!=row2.classification):
        return(False)
    #Checks in the indices
    if (int(row1.start)>=int(row2.start) and int(row1.stop)>=int(row2.stop)) and int(row1.start)<=int(row2.stop):
        return(True)
    if (int(row1.start)<=int(row2.start) and int(row1.stop)<=int(row2.stop)) and int(row1.stop)>=int(row2.start):
        return(True)
    if (int(row1.start)<=int(row2.start) and int(row1.stop)>=int(row2.stop)):
        return(True)
    if (int(row1.start)>=int(row2.start) and int(row1.stop)<=int(row2.stop)):
        return(True)
    return(False)


def mergeRowElements(row1,row2):
    """
    Merges two row objects into a new row object which is returned.
    Inputs:
        row1 and row2(gtfRow)
    Returns:
        GTFRow
    """
    start1 = int(row1.start)
    start2 = int(row2.start)
    stop1 = int(row1.stop)
    stop2 = int(row2.stop)
    newStart=0
    newStop=1
    
      #Chose the smaller of the two start values
    if start1<=start2:
        newStart=start1
    else:
        newStart=start2
        
    #Chose the smaller of the two start values
    if (stop1>=stop2):
        newStop=stop1
    else:
        newStop=stop2

    mergedRow = GTFRow(row1.chromosome,row1.geneName,newStart,newStop,row1.orientation,row1.transcriptID,row1.classification)
    return(mergedRow)


def mergeSequenceElements(sEle1,sEle2):
    """
    Merges two sequence elements into a new sequence element object which is returned.
    Inputs:
        sEle1 and sEle2 (sequenceElement)
    Returns:
        sequenceElement
    """
    start1 = int(sEle1.start)
    start2 = int(sEle2.start)
    stop1 = int(sEle1.stop)
    stop2 = int(sEle2.stop)
    newStart=0
    newStop=1
    
    #Chose the smaller of the two start values
    if start1<=start2:
        newStart=start1
    else:
        newStart=start2
        
    #Chose the smaller of the two start values
    if (stop1>=stop2):
        newStop=stop1
    else:
        newStop=stop2

    #Get the ID
    newID = sEle1.ID
    
    #Return the merged object
    return(sequenceElement(newStart,newStop,newID))
    

def checkForConseqUTR(lstofSeqElem):
    """
    Takes in a list of sequence elements and returns whether there are any UTR next to each other in said list.
    Input:
        lstofSeqElem (List[sequenceElement]) - Sequence elements to check for adjacent UTR region.
    Output:
        bool
    """
    for i in range(0,len(lstofSeqElem)):
        ele1=lstofSeqElem[i]
        if ele1.ID=="UTR":
            if i<len(lstofSeqElem)-1:
                ele2 = lstofSeqElem[i+1]
                if ele2.ID=="UTR":
                    return True
    return False
    

def mergeFirst(lstofSeqElem):
    """
    Takes a list of sequence elements and merges the first copies of UTR into a single .
    Input: 
        lstofSeqElem (List[sequenceElement]) - Original list.
    Output:
        List[sequenceElement] - List with the first overlapping UTR regions merged.
    """
    w =0
    for ele1 in lstofSeqElem:
        iden1 = ele1.ID
        if (w <len(lstofSeqElem)-1):
            ele2 = lstofSeqElem[w+1]
            iden2 = ele2.ID
            #If the ID's match
            if (iden1=="UTR"and iden2=="UTR"):
                #Create a merged on at the second index
                lstofSeqElem[w+1]=mergeSequenceElements(ele1,ele2)
                #Delete the first element
                lstofSeqElem.pop(w)
                #Return the new list
                return(lstofSeqElem)
        w=w+1 


def mergeUTR(lstofSeqElem):
   """
   Takes in a list of sequence elements and merges all the UTR into single UTR.
   Input: 
        lstofSeqElem (List[sequenceElement]) - Original list.
    Output:
        List[sequenceElement] - List in which all of the UTR regions have been merged.
   """
   mL=lstofSeqElem
   conseq = checkForConseqUTR(mL)
   while conseq==True:
       #Merge one
       mL=mergeFirst(mL)
       conseq=checkForConseqUTR(mL)
   return(mL)
                     
    
def withinQuotes(myStr):
    """
    Returns the first part of a string which is within which is within ""
    Input:
        myStr (Str) - Input string containing " and "
    Output:
        newString (Str) - The substring of the originial string representing the first conetents of "". 
    """
    #Set indices to use
    startIndex = 0
    stopIndex = 0
    for w in range(len(myStr)):
        curChar = myStr[w]
        #If this is "
        if (curChar=="\""):
            if(startIndex==0):
                startIndex=w
            else:
                stopIndex=w
    newString = myStr[startIndex+1:stopIndex]
    return(newString)        


def getgtf (gtfFile,useIDInsteadOfName):
    """
    This function reads in a gtf file and re-represents it as a list of lists, where each individual list corresponds to a row in the gtf file (containing only the key information).
    Pulls lines from a tsv after skipping commented files. It then takes the relevant entries. Because GTF files store some fields in string lists (long strings with multiple entries not divided by a tab character), some creative type conversion is necessary to pull specific fields.
    Input:
        gtfFile(Str) - The name (including path if outside the working directory) of the gtf file to read.
        useIDInsteadOfName (bool) - Whether you want to use the gene ID field instead of the gene name field. Again, this should match what is going on in your binding cluster annotation process.
    Outputs:
        mstrList(List[List[Str]]) - The list of lists that represents the rows of the inputed gtf file.
    """
    mstrList=[]
    #Loop over each row, skipping those that are comment rows (rows that start with #)
    with open(gtfFile) as tsv:
        for line in csv.reader(tsv, dialect="excel-tab"): #You can also use delimiter="\t" rather than giving a dialect.
            if len(line)==0:
                continue
            toCheck=line[0][0]
            if toCheck=='#':
                print("Skipping commented line: "+str(line))
            if toCheck!="#":
                #Create a list of values that includes the seperated mess in the final entry
                indivLine=[]
                for indivEnt in line:
                    if ';' not in indivEnt:
                        indivLine.append(indivEnt)
                    #For those entires that are lists of values stored inconveniently as strings
                    elif ';' in indivEnt:
                        #Convert ;'s to ,'s
                        toAdd = indivEnt.replace(';', ',')
                        toAdd = toAdd.split(",")
                        
                        #Then pull the actual values from each string (not the header)
                        #Add two empty spacers to be replaced with transcript id and gene name
                        indivLine.append("NoTranscriptIDGiven")
                        if useIDInsteadOfName==False:
                            indivLine.append("NoGeneNameGiven")
                        else:
                            indivLine.append("NoGeneIDGiven")
                        for val in toAdd:
                            valList=val.split(" ")
                            #Remove the awkward empty spaces
                            if '' in valList:
                                valList.remove('')
                            if len(valList)>1:
                                if valList[0] =="transcript_id":
                                    toUse=valList[1]
                                    toUse=toUse.replace('"','')
                                    indivLine[8]=toUse
                                if valList[0]=="gene_name" and useIDInsteadOfName==False:
                                    toUse=valList[1]
                                    toUse=toUse.replace('"','')
                                    indivLine[9]=toUse                    
                                if valList[0]=="gene_id" and useIDInsteadOfName==True:
                                    toUse=valList[1]
                                    toUse=toUse.replace('"','')
                                    indivLine[9]=toUse  
                if len(indivLine)!=10:
                    print("Error: individual row list the wrong length (not 10). The list is:",indivLine)
                #print(indivLine)
                mstrList.append(indivLine)
                
    return(mstrList) 
    
    
def gtfObjectLst (lstofLsts):
    """
    This function takes a list of lists that represent the gtf file rows and returns a list of gtfRow objects for calculations downstream.
    Input:
        lstofLsts(List[List[Str]]) - The list of lists that represents the rows of the inputed gtf file.
    Returns:
        lstofRowObjects(List[GTFRow]) - A list of gtf row objects that represents the inputed gtf file.
    """
    #Initialize the list for the gtfRow objects
    lstofRowObjects = []
    #For each item in the individual list
    for ele in lstofLsts:
        #Get the required properties
        #Modulate the classification parameter to ensure the same term is used no matter the slight differences between inputed gtf files
        if "five_prime_utr" in ele[2]:
                classif="FUTR"
        if "three_prime_utr" in ele[2]:
                classif="TUTR"
        if ("five_prime_utr" not in ele[2]) and ("three_prime_utr" not in ele[2]):
            classif = ele[2]
        chromo = ele[0]
        geneN = ele[9]
        start = ele[3]
        stop = ele[4]
        orien = ele[6]
        tID = ele[8]
        currentgtfRow = GTFRow(chromo,geneN,start,stop,orien,tID,classif)
        lstofRowObjects.append(currentgtfRow)

    return(lstofRowObjects)
    
    
def createGeneLsts (gtfLst):
    """
    This function takes in a list of row objects and creates a list of gene objects (which have lists of gtf row objects).
    Inputs:
        gtfLst(List[GTFRow]) - A list of gtf rows objects to be converted and sorted by gene.
    Output:
        masterLst(List[gtfGene]) - A re-representation of the gtf file in terms of gene objects.
    """
    #Sort by gene name
    sorted_gtfLst=gtfLst
    sorted_gtfLst.sort(key=lambda x: x.geneName, reverse=True)
    
    counter=-1
    #Masterlist of gene objects 
    masterLst = []
    #Create a previous name section to allow for sorting by gene
    prevGeneName=""
    #Loop over each gtfRow object in the list of them
    for curObj in sorted_gtfLst:
        #Get the gene name of that item
        curGeneName = curObj.getGeneName()
        #If this is not the same as the previous we need a new gene
        if (curGeneName != prevGeneName):
            #Create a new gene object and add i
            newGene = gtfGene(curGeneName)
            #If this object is a transcript object, mark it as such so that we know we don't need to generate our own later.
            if curObj.classification=="transcript":
                newGene.hasTranscriptLabels=True
            #Add the current gtfRow object to this list
            newGene.addRow(curObj)
            #Add this gtfGene object to the master list
            masterLst.append(newGene)
            counter=counter+1
        #If this is NOT a new gene
        if (curGeneName==prevGeneName):
            #Add it to the current gtfGene object
                #The current
            curGene = masterLst[counter]
            #If this object is a transcript object, mark it as such so that we know we don't need to generate our own later.
            if curObj.classification=="transcript":
                curGene.hasTranscriptLabels=True
                #Add
            curGene.addRow(curObj)
        prevGeneName=curGeneName   
    return(masterLst)
        
    
def sortByStart(lstOfSeq):
    """
    Sorts a list of sequence element objects based on their starting values.
    Input:
        lstOfSeq ([sequenceElement])
    Returns:
        List[sequenceElement] - Sorted list based on the starting value.
    """
    #Get the item
    newLst=[]
    tupLst=[]
    newTupLst=[]
    #End length
    #Create the tuple list to sort (first = value second =index)
    for i in range(0,len(lstOfSeq)):
        ele= lstOfSeq[i]
        nTup = (int(ele.start),i)
        tupLst.append(nTup)         
    
    #Sort this list of tuples
    while (len(tupLst)>0):
        smallest = 10000000000000000000000000000
        index=0
        for j in range(0,len(tupLst)):
            curTup = tupLst[j]
            #Get current value
            curVal = curTup[0]
            #Compare to the smallest value
            if (curVal <=smallest):
                smallest=curVal
                index=j
        #Select the smallest to add to new
        newTupLst.append(tupLst[index])
        #Remove this element
        del tupLst[index]
    
    #Now that we have the tuple version simply make the list of sequence elements version
    for w in range(0,len(newTupLst)):
        #Add the corresponding element to the new list
        theTup = newTupLst[w]
        theInd = theTup[1]
        newLst.append(lstOfSeq[theInd])
    
    return(newLst)
  
    
def convertRowIntoOutput(gtfRow):
    """
    Converta a gtf row object into a lst to be outputed as a row in the final.
    Input:
        gtfRow (gtfRow) - The gtfRow object to be converted.
    Returns:
        List[str]
    """
    toRet=[gtfRow.geneName,gtfRow.classification,gtfRow.start,gtfRow.stop,gtfRow.chromosome,gtfRow.orientation]
    return(toRet)

###------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Run function
###------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def writeGTFGeneList (gtfFile,outDir,outName,geneID):
    """
    This function runs all calculations to create a list of genes objects which represents the full gtf file, which is written to an output csv file for use in the metegene generating programs.
    Input:
        gtfFile (str) - File name of the gtf file to convert into a longest transcript intermediate metagene file (include path to file if the file is not in the working directory)
        outDir (str) - Directory into which the intermediate file is to be written.
        outName (str) - Prefix used to distinguish and name the gtf intermediate file created. Standard suffix used for identification downstream.
        geneID (bool) - Should the intervals be listed by gene name or by the gene ID (often gtf regions map to a geneID but not a gene symbol).
        
    Output:
        asGenes (List[gtfGene]) - A list of gtf Genes that corresponds to the contents of the GTF file. 
    """
    start = time.time()
    #Create empty data lists for export later
    myData = []
    myDataIE=[]
    headers = ["Gene","Type","Start","Stop","Chromosome","Orientation"]
    myData.append(headers)
    myDataIE.append(headers)
    
    #Get the gtf file
    print("Reading gtf file to a list of lists...")
    gtf = getgtf(gtfFile,geneID)
    print("...done.")
    
    #Convert it to a list of gtfrow objects
    print("Convering the list of lists to a list of gtf row objects...")
    gtfasObjects = gtfObjectLst(gtf)
    print("...done.")
    
    #Seperate these gtfrow objects into genes
    print("Seperating the gtfrow objects into genes...")
    asGenes = createGeneLsts(gtfasObjects)
    print("...done.")
    
        #Populate the transcripts for these genes
    print("Populating the transcripts for each gene...")
    for gene in asGenes:
        gene.createTranscriptRows()
        gene.populateTranscriptLst()
        #Select the longest transcript to be the one defining the rows
        gene.setRowsAsLongestTrans()
        
    print("...done.")         
    
        #Create a duplicate for the I/E intermediate file
    IEDup=asGenes
    
        #Populate the values for the final csv - I/E
    print("Populating the values for the final Intron/Exon csv...")
    for gene in IEDup:
        gene.finalCleanUp(True)
        gRows = gene.lstOfRows
        for row in gRows:
            if row.classification!="toRem":
                #Convert the row into a lst to add to the myData lst of lsts
                nList = convertRowIntoOutput(row)
                #Only add it if it is NOT a transcript row
                if nList[1]!="transcript" and nList[1]!="gene":
                    myDataIE.append(nList)
    
    #Establish the name that will be used to output the results
    theName=outName+".longestTranscript.csv"
    finalOut=outDir+"/"+theName
    myFile = open(finalOut,'w')
    with myFile:
        writer = csv.writer(myFile)
        writer.writerows(myDataIE)
    
    print("...writing complete.")
   
    end = time.time()
    print("Time to run (s): "+str(end - start))

    return(asGenes)

###----------------------------------------------------------------------------
#Run
###----------------------------------------------------------------------------
#geneList = writeGTFGeneList(inpGTF,outDir,outName,byID)