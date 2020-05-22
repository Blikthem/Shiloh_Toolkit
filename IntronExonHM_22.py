#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 12:35:12 2018

@author: claypooldj
"""

####Load packages 
import csv
import os
import random
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import math
import pandas as pd
import seaborn as sns

"""
##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Inputs
##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ 
gtfFileName="/home/claypooldj/genomes/gencode_hg38.longestTranscript.csv"    #The intermediate gtf longest splice varient csv file (generated from a seperate script) to use for annotation. {Type - String}
csvDir="/home/claypooldj/myPythonScripts/testBedHM/hInputs/"    #Tlhe directory containing the cluster csv files to analyze. All csv files in this directory will be tested (so make sure there aren't other csv files in the directory or it will crash!) {Type - String}
outputDir="/home/claypooldj/myPythonScripts/testBedHM/output/"     #The directory where you want to generate all output files. {Type - String}
gtu=[]    #The list of genes which you wish to consider. Just put an empty list '[]' if you want to consider all genes. {Type - List[String] or empty List}
clustThreshold=1    #The minimum number of clusters a gene must have to be considered {Type - Int}
boundF="URC"    #The continuous field in the clusters file to use to apply filtering bounds (both upper and lower) {Type - String}
lBound=0    #The percentage of clusters to filter based on the lowest values of the bounding continuous variable (boundF) {Type - Num}
uBound=0    #The percentage of clusters to filter based on the largest values of the bounding continuous variable (boundF) {Type - Num}
wValuesByRC=False    #Whether or not the program should weight the impact of each indivudal gene on the overall metagene average linearly with the number of clusters that align to that gene {Type - Bool}
#randStates=[7211995,541995,3131994,111,222,333,444,555,888,999]   #The random states to use in the binning algorithm, which randomly deals with rounding error and then averages the results. Can be any length. {Type- List[Int]}
randStates=[999]
dpi=250    #The resolution of the output graphs in dpi {Type - Int}
imgFrmt="pdf"    #Which format to save the output graphs to. Two options, either "pdf" which saves as pdf otherwise it saves to png format. {Type - String}
wtGE="/Volumes/Untitled/Output/explorador/RNA_Seq_Background/moreWTExpressions/geneExpressionMatrix_TPM_withGeneNames.csv"    #A csv file containing the wild type gene expression levels from RNA-Seq(s) of the cell line. Takes average value. If this string is empty, we will not weigh values by WT gene expression.
wtGE=""
mainChromosomes=True #Whether we should only consider the main chromsomes (the autosomal and X/Y) in the analysis and not the strange chromosome constructs used for.
theBED_GTF="/home/claypooldj/genomes/gencode.v30.annotation.gtf"
theBEDAnotScript="/home/claypooldj/clip_metagene/bedannotator.sh"
##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""

##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Define objects and their various methods
##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ 
class anotRow:
    """
    Anotation Row Object: This object represents a row from the anotation (GTF) file.

    Properties:
        gene(String) - The name of the gene this row maps to
        ty(String) - The type (5'/CDS/3' etc) of this row
        start(int) - The starting index of this row
        stop(int) - The end index of this row
        chrom(String) - The chromosome this row belongs to
        orientation(string) - The orientation of the row, either '+' or '-'
    
    Methods:
        init - Initializer
        str - Conversion to string form for printing
        coversIndex - Checks if an index is within the anotation row object
    """
    def __init__(self, gene, ty, start, stop, chrom,orientation):
        """
        Initalizing method.cl
        """
        self.gene=gene
        self.ty=ty
        self.start=start
        self.stop=stop
        self.chrom=chrom
        self.orientation=orientation
    
    def __str__(self):
        """
        Returns a string that represents the object
        """
        toPrint="{Annotaion File Row Object} Start: "+str(self.start)+"\tStop: "+str(self.stop)+"\tGene: "+self.gene+"\tChromosome: "+self.chrom+"\tType: "+self.ty
        return(toPrint)
        
    def coversIndex(self,index):
        """
        Checks to see if a given index is within the boudns of this anotation row object
        Inputs:
            index(int) - The index to check
        Returns:
            boolean - Whether the input index is within the bounds of the object
        """
        if (self.stop>=index and index>=self.start):
            return(True)
        return(False)

class anotGene:
    """
    Anotation Gene Object: This object represents an anotation gene (a collection of row objects from the annotation gtf file that represent a single gene at the longest transcript representation).

    Properties:
        gene(String) - The name of the gene these rows map to.
        chrom(String) - The chromosome these rows map to.
        start(int) - The smallest starting index for all annotation rows in this gene.
        stop(int) - The largest ending index for all annotation rows in this gene.
        lstOfAnotRows(list[anotRows]) - A list of all of the anotRow objects in this anotationGene object
        exonIndexStarts (list[int]) - A list of all the starting positions of exons in the longest transcript splice vairent of this gene.
        exonIndexStops (list[int]) - A list of all the stoping positions of exons in the longest transcript splice vairent of this gene.
        intronIndexStarts (list[int]) - A list of all the starting positions of introns in the longest transcript splice vairent of this gene.
        intronIndexStops (list[int]) - A list of all the stoping positions of introns in the longest transcript splice vairent of this gene.
    
    Methods:
        init - Initializer
        str - Conversion to string form for printing.
        getStart - Returns the smallest(first) starting index of all the rows in this anotationGene.
        addAR - Adds a row to this anotationGene object.
        getStartandStopIndices - Populates the start and stop lists for this object.
        getType - Returns the classification (CDS/TUTR/etc) associated with a given nucleotide index in this anotationGene.
        sideOfAdjacentIntron - Check if a given exon in the object is bordered by introns and where.
        getExonAt - Returns the exon from this gene located at an inputed position.
        exonPositionIntronAdjacent - Checks if a position within an exon is within a half adjacent to an intron.
        getType - Indicates the annotation category at a given index within the anotation gene.
        
    """
    def __init__(self,gene,chrom):
        self.gene=gene
        self.chrom=chrom
        self.start=0
        self.stop=1
        self.lstOfAnotRows=[]
        self.exonIndexStarts=[]
        self.exonIndexStops=[]
        self.intronIndexStarts=[]
        self.intronIndexStops=[]
        
        
    def __str__(self):
        toPrint="{Annotaion Gene Object} Gene: "+self.gene+"\tChromosome: "+self.chrom+"\tNumber of Annotation Row Objects: "+str(len(self.lstOfAnotRows))
        return(toPrint)
    
    def repairGaps(self):
        """
        Searching for gaps between exons in the metagene and, if any are found, fills them with introns.
        """
        nList=[]
        #Sort
        sortedList=self.lstOfAnotRows
        sortedList.sort(key=lambda x: int(x.start))

        #Loop over each exon...
        for i in range(0,len(sortedList)):
            ro=sortedList[i]
            if ro.ty=="exon":
                #Check to make sure that it is immediately followed by an intron (that is, check to be sure there is no gap)
                if i+1<len(sortedList):
                    theNext=sortedList[i+1]
                    if int(theNext.start)!=(ro.stop)+1:
                        #If it is going into an intron simply expand the intron
                        if theNext.ty=="intron":
                            theNext.start=int(ro.stop)+1
                        
                        #Otherwise create a new row element
                        if theNext.ty!="intron":
                            #Create a new intron
                            nIntron=anotRow(ro.gene,"intron",int(ro.stop)+1,int(theNext.start)-1,ro.chrom,ro.orientation)
                            if nIntron not in sortedList:
                                nList.append(nIntron)
      
        for ro in sortedList:
            nList.append(ro)
        
        self.lstOfAnotRows=nList
        
     
        
    def getStart(self):
        """
        Finds the smallest starting index of the rows in the anotGene.
        Returns:
            int - The first nucleotide position in the anotGene
        """ 
        possibleVals=[]
        for row in self.lstOfAnotRows:
            possibleVals.append(row.start)
        if (len(possibleVals)==0):
            return(False)
        theMin = min(possibleVals)
        self.start=theMin
        return(theMin)


    def sideOfAdjacentIntron(self,exon):
        """
        Checks if an exon is bordered by introns on the Left (L), Right (R), or Both (B) sides or Neither (N)
        Inputs:
            exon (anotRow) - The exon to check for bordering introns.
        Returns:
            Str - Is this exon bordered by introns: on both sides("B"), one the left side ("L"), or the right side ("R"), or on neither side ("N")
        """
        
        intronBefore=False
        intronAfter=False
        exonIndex=0
        #Get the position of the exon within the list of anot rows
        for i in range(0,len(self.lstOfAnotRows)):
            ro = self.lstOfAnotRows[i] 
            #If it is the same as the exon
            if ro.start==exon.start and ro.stop==exon.stop and ro.ty=="exon":
                exonIndex=i
        #Check if it has any introns before it 
        if (exonIndex>0):
            j=exonIndex-1
            while j!=0:
                prevRo=self.lstOfAnotRows[j]  
                if (prevRo.ty=="exon"):
                    break
                if (prevRo.ty=="intron"):
                    intronBefore=True
                j=j-1
                
        #Check if it has any introns after i
        if (exonIndex!=len(self.lstOfAnotRows)):
            for k in range(exonIndex+1,len(self.lstOfAnotRows)):
                nxtRo = self.lstOfAnotRows[k]
                if (nxtRo.ty=="exon"):
                    break
                if (nxtRo.ty=="intron"):
                    intronAfter=True
                    
        
        #Now figure out what to return
        if intronBefore==True and intronAfter==True:
            return("B")
        
        if intronBefore==True:
            return("L")
            
        if intronAfter==True:
            return("R")
            
        return("N")
        
        
    def getExonAt(self,position):
        """
        Returns the exon at a particular position.
        Input:
            position (Int) - The position to check for an exon.
        Returns:
            anotRow - The anotation row that is the exon at this position.
        """
        for ro in self.lstOfAnotRows:
            #If it is exon
            if ro.ty=="exon":
                if (ro.stop>=position and ro.start<=position):
                    return(ro)
        return(False)            

    def exonPositionIntronAdjacent(self,position):
        """
        Checks if a position within an exon is within a half adjacent to an intron.
        Input:
            position (int) - The index to check if it is adjacent to an intron.
        Returns:
            bool (is this position in a half adjacent to an exon?)
        """
        #Create the bounds that represent the first and second half of the exon
        exon = self.getExonAt(position)
        #Get the length
        totLength = exon.stop-exon.start
            #Divide it by two and round down to an even number
        byTwo = int(totLength/2)
        endOfFirstHalf=exon.start+byTwo
        startOfSecondHalf=endOfFirstHalf+1
        
        #If the exon is not adjacent to anything then we know this index is no good
        if(self.sideOfAdjacentIntron(exon)=="N"):
            return(False)
        #If it is flanekd on both sides by introns then we know any index in the exon will be good
        if(self.sideOfAdjacentIntron(exon)=="B"):
            return(True)
            
        #If an intronic region flnaks on the left then the index will have to be in the first half of the exon
        if (self.sideOfAdjacentIntron(exon)=="L"):
            if position<=endOfFirstHalf and position>=exon.start:
                return(True)
            else:
                return(False)
        
        #If an intronic region flanks on the right then the index will have to be in the second half of the exon
        if (self.sideOfAdjacentIntron(exon)=="R"):
            if position>=startOfSecondHalf and position<=exon.stop:
                return(True)
            else:
                return(False)

    def addAR(self,toAdd):
        """
        Add an anotRow object to this anotGene
        Input:
            toAdd (anotRow) - The row object to add.
        """ 
        self.lstOfAnotRows.append(toAdd)
    
    def getStartandStopIndices(self):
        """
        Gets the indices for all of the individual groups with respect to the length of the whole hitlits and sets them to the object.
        """
        gStart=int(self.getStart())


        if (len(self.intronIndexStops)>0 or len(self.intronIndexStarts)>0 or len(self.exonIndexStarts)>0 or len(self.exonIndexStops)>0):
            self.intronIndexStops=[]
            self.intronIndexStarts=[]
            self.exonIndexStarts=[]
            self.exonIndexStops=[]
        firstExon=True
        for ro in self.lstOfAnotRows:
                if ro.ty=="intron":
                    #Get the starting value
                    fullVal = int(ro.start)
                    #Get the value in context of the overall
                    nVal = fullVal-gStart
                    if nVal not in self.intronIndexStarts:
                        self.intronIndexStarts.append(nVal)
                    
                    #Get the stopping value
                    fullVal2=int(ro.stop)
                    sVal=fullVal2-gStart
                    if sVal not in self.intronIndexStops:
                        self.intronIndexStops.append(sVal)
                    
                if ro.ty=="exon":
                    #If this is the first exon we just need to add one Start and stop
                    if firstExon==True:
                        #Get the start 
                        fullVal = int(ro.start)
                        #Get the original start in context of the overall
                        nVal = fullVal-gStart                      
                    
                        #Get the stop
                        fullVal2=int(ro.stop)
                        sVal=fullVal2-gStart
                        
                        #The first start value - halfway in between (Take the difference and divide by two, rounding down)
                        theDif = sVal-nVal

                        ##RANDOMLY CHOOSE BETWEEN ROUNDING UP AND ROUNDING DOWN
                        ##NEW
                        rdru=random.randint(0,1)
                            #Round down
                        if rdru==0:
                            nStart = nVal + int(theDif/2)
                            #Round up
                        elif rdru==1:
                            nStart = nVal + math.ceil(theDif/2)
                        
                        if nStart not in self.exonIndexStarts:
                            self.exonIndexStarts.append(nStart)
                        
                        #The first stop value - the end of the exon
                        if sVal not in self.exonIndexStops:
                            self.exonIndexStops.append(sVal)
                        
                        #Reset the first exon
                        firstExon=False
                    
                    #Otherwise we need to add two starts and stops
                    elif firstExon==False:
                        #Get the start
                        fullVal = int(ro.start)
                        #Get the value in context of the overall
                        nVal = fullVal-gStart
                    
                        #Get the stop
                        fullVal2=int(ro.stop)
                        sVal=fullVal2-gStart
                        
                        #Add start 1
                        if nVal not in self.exonIndexStarts:
                            self.exonIndexStarts.append(nVal)
                        
                        #Add stop 1
                        theDif = sVal-nVal
                        stop1 = nVal + int(theDif/2)
                        if stop1 not in self.exonIndexStops:
                            self.exonIndexStops.append(stop1)
                        
                        #Add start 2
                        start2=stop1+1
                        if start2 not in self.exonIndexStarts:
                            self.exonIndexStarts.append(start2)
                        
                        #Add stop 2
                        if sVal not in self.exonIndexStops:
                            self.exonIndexStops.append(sVal)
        #Remove the last start and stops
        self.exonIndexStarts.pop()
        self.exonIndexStops.pop()
        
    def getType(self,index):
        """
        Returns a character indicating the classification of the gene at a particular index
        Input:
            index (int) - The index the user wants to know the classification of.
        Output
            str  - Indicates the annotation classification at the inputed index.
            Options - "exon", "intron", "FUTR" (5' UTR), "TUTR" (3' UTR), "CDS", and "NC" (no category)
        """
        possibleOptions=[]
        #Go through each row of the anotGene
        for r in range(0,len(self.lstOfAnotRows)):
            ro = self.lstOfAnotRows[r]
            #If it is contained
            if ro.coversIndex(index)==True:
                possibleOptions.append(ro.ty)
                
        #Now pic which one to return
        if "intron" in possibleOptions:
            return("intron")
        if "exon" in possibleOptions:
            return("exon")
        if "FUTR" in possibleOptions:
            return("FUTR")
        if "TUTR" in possibleOptions:
            return("TUTR")
        if "CDS" in possibleOptions:
            return("CDS")
        return("NC")
 
        
class chromesome:
    """
    Chromosome Object: This object represents a chromosome and stores anotation elements connected to a chromosome ID.

    Properties:
        ID(String) - The name of the chromosome container
        lstOfAnotElements(list) - A list of annotation elements found in this chromosome. Intentionally nonspecific to fascilitate general use
    
    Methods:
        init - Initializer
        str - Conversion to string form for printing
    """
    def __init__(self,ID):
        self.ID=ID
        self.lstOfAnotElements=[]
        
    def __str__(self):
            toPrint="{Chromosome Object} ID: "+self.ID
            return(toPrint)
            
            
class WTGeneExpression:
    """
    WTGeneExpression Object: This object represents the gene expression profile of the wild type cell type of interest

    Properties:
        fileName (Str) - The file's name with directory if not contained in the working directory.
        DF(Pandas DataFrame) - Pandas data frame contianing the gene expression information. The gene name
        geDict (Dictionary) - A dictionary which links the gene name (str key) to the average gene 
    
    Methods:
        init - Initializer
        str - Conversion to string form for printing
    """
    def __init__(self,fileName):
        self.fileName=fileName
        self.DF=pd.read_csv(fileName)
        self.geDict=self.populateDict()
        
    def __str__(self):
            toPrint="{WTGeneExpression Object} ID: "+self.ID
            return(toPrint)

    def populateDict(self):
        """
        Uses the DataFrame stored in this object to generate a dictionary that maps gene symbol to average gene expression.
        Returns:
            geneDict (Dictionary) - Link between Str (gene symbol) and float (average gene expression in the wild type).
        """
        #Initialize the dictionary to be used to save the values
        toRet={}
        #Loop over the rows
        for index, row in self.DF.iterrows():
            theKey=str(row["gene_symbol"])
            theValue=float(row["Average"])
            toRet[theKey]=theValue
        return(toRet)

class anotFile:
    """
    anotFile Object: The object which contains and manipulates all of the genes in an annotation gtf file (loaded from the custom-made gtf longest transcript intermediate CSV file).

    Properties:
        lstOfChromosomeNames (List[str]) - The list of the unique chromosomes to be considered.
        lstOfChromosomesRows(List[chromosome]) - The list of chromosome objects which contain the anotation rows making up the gene annotation file (sorted by chrmosome).
        lstOfChromosomesGenes(List[chromosome]) - The list of chromosome objects which contain the anotation genes making up the gene annotation file (sorted by chrmosome).
    
    Methods:
        init - Initializer
        str - Conversion to string form for printing.
        getNumRows - Determines how many anotRow objects are stored in the sorted chromosome objects of this anotFile object.
        getNumGenes - Determines how many anotGene objects are stored in the sorted chromosome objects of this anotFile object.
        addAR - Adds an anotRow to the correct chromosome object in this object's list of chromosome rows.
        addAG - Adds an anotGene to the correct chromosome object in this object's list of chromosome genes.
        populateCh - Populate the lst of chromesome gene objects.
    """
    def __init__(self,chromosomeNames):
        self.lstOfChromosomeNames=chromosomeNames
        self.lstOfChromosomesRows=[]
        self.lstOfChromosomesGenes=[]
        for indivN in chromosomeNames:
            self.lstOfChromosomesRows.append(chromesome(indivN))
            self.lstOfChromosomesGenes.append(chromesome(indivN))
    
    def __str__(self):
        toPrint="{Annotaion File Object} Number of Annotation Rows: "+str(self.getNumRows())+"\tNumber of Annotation Genes:  "+str(self.getNumGenes())
        return(toPrint)
        
    def getNumRows(self):
        """
        Determines how many anotRow objects are stored in the sorted chromosome objects of this anotFile object.
        Returns:
            int
        """
        toReturn=0
        for ch in self.lstOfChromosomesRows:
            toReturn=toReturn+len(ch.lstOfAnotElements)
        return(toReturn)
            
    def getNumGenes(self):
        """
        Adds an anotGene to the correct chromosome object in this object's list of chromosome genes.
        Returns:
            int
        """
        toReturn=0
        for ch in self.lstOfChromosomesGenes:
            toReturn=toReturn+len(ch.lstOfAnotElements)
        return(toReturn)    
        
    def addAR(self,rowObj):
        """
        Adds an anotRow to the correct chromosome object in this object's list of chromosome rows.
        Input:
            rowObj(AnotRow) - The row to add to this object.
        """
        #Get the correct chromosome
        curID = rowObj.chrom
        for chrom in self.lstOfChromosomesRows:
            if (chrom.ID==curID):
                chrom.lstOfAnotElements.append(rowObj)
                return()
            
    def addAG(self,anotGeneObj):
        """
        Adds an anotGene to the correct chromosome object in this object's list of chromosome genes.
        Input:
            anotGeneObj (anotGene) - Anotation gene object to store.
        """
        curID = anotGeneObj.chrom
        for chrom in self.lstOfAnotGenes:
            if(chrom.ID==curID):
                chrom.lstOfAnotElements.append(anotGeneObj)
                return()
            
    def populateCh(self): 
        """
        Populate the lst of chromesome gene objects stored in this object.
        """
        #for each row chromosome
        for i in range(0,len(self.lstOfChromosomesRows)):
            rChr = self.lstOfChromosomesRows[i]
            gChr= self.lstOfChromosomesGenes[i]
            #For each row in the chromosome
            for indivRow in rChr.lstOfAnotElements:
                counter=0
                #Check if this is in the corresponding chromosome for the gene list
                for indivGene in gChr.lstOfAnotElements:
                    if (indivGene.gene==indivRow.gene):
                        indivGene.addAR(indivRow)
                        counter=1
                        break
                if counter==0:
                    #Make a new anotGene
                    nAnotGene=anotGene(indivRow.gene,indivRow.chrom)
                    nAnotGene.addAR(indivRow)
                    gChr.lstOfAnotElements.append(nAnotGene)
        
        
class cluster:
    """
    Cluster Object: This object represents a binding cluster.

    Properties:
        Pulled from csv rows as strings - chrom (chromosome), strand (strand oriented to), start (starting index in cluster), end (ending index of cluster), gene (gene that paralyzer aligned this cluster to if any), CS (conversion specificity), T2C (T to C fraction), URC (Unique read count for this cluster), RC (total read count)
    
    Methods:
        init - Initializer
        str - Conversion to string form for printing
        addAR - Add an anot
    """
    def __init__(self,chrom,strand,start,end,gene,CS,T2C,URC,RC):
        self.chrom=chrom
        self.strand=strand
        self.start=start
        self.end=end
        self.gene=gene
        self.CS=CS
        self.T2C=CS
        self.URC=URC
        self.RC=RC
        
    def __str__(self):
            toPrint="{Cluster Object} Gene: "+self.gene+"\tStart: "+str(self.start)+"\tEnd: "+str(self.end)+"\tChromosome: "+str(self.chrom)
            return(toPrint)                
        
class clusterGene:
    """
    ClusterGene Object: A group of clusters of the same gene

    Properties:
        gene(string) - The name of the gene around which all contained clusters are aligned
        chrom(string) - The name of the chromosome this gene falls on
        start(int) - The start location of this clustergene (ie the first start value among contained clusters)
        stop(int) - The stop location of this clustergene (ie the last stop value among contained clusters)
        lstOfClusters (list[cluster]) - A list of all of the cluster objects aligned with this gene
        hitListsIntron (List[List[Int]]) - A list of lists, where each inner list represents the hit list for an intronin this gene. That is to say, the inner list is the same length as the intron and says whether each index overlaps with a binding cluster (1) or does not (0).
        hitListsExon (List[List[Int]]) - A list of lists, where each inner list represents the hit list for an exon in this gene. That is to say, the inner list is the same length as the exon and says whether each index overlaps with a binding cluster (1) or does not (0).
        hitSums(List[int]) - The sum of hits in all intron/exon splice junctions of this gene.

    Methods:
        init - Initializer
        str - Conversion to string form for printing
        addAC - Add a cluster object to this clusterGene
        getTotalNumReads - Return the total number of reads in clusters aligned to this gene.
        doesGeneHaveNClust - Check if at least N clusters fall within this gene.
        populateNumClust - Populates the number of clusters that overlap the intronic and exonic regions of the cluster gene.
        containsIndex - Checks if a cluster in this clusterGene contains a given index.
        removeGene - Flag this gene so that it is no longer considered for downstream analyses.
        getMatchingMeta - Finds the metaGene in an anotFile object which has the same name as this clusterGene.
        populateHLTotals - Determines the total number of hits in each EIE junction, which is used to calculate the percentage of the total each overlapped nucleotide represents.
        populateHitLists - Populates the hit lists for this clusterGene object (a list which states whether each nucleotide in a region is overlapped with a cluster or not).
        checkSumTo100 - Confirms that the percentage lists in this cluster gene indeed sum to 100%.
        populateIntronBinsOfPercents - Re-express cluster overlap in the intronic region as percentages in a fixed number of bins.
        populateExonBinsOfPercents - Re-express cluster overlap in the exonic region as percentages in a fixed number of bins.
        populateGeneDistributions - Averages the percentage lists for each EIE in this gene.

    """
    def __init__(self,gene,chrom):
        self.gene=gene
        self.chrom=chrom
        self.start=0
        self.stop=1
        self.lstOfClusters=[]
        self.hitListsIntron=[]
        self.hitListsExon=[]
        self.hitSums=[]
        self.lstOfExonBinLists=[]
        self.lstOfIntronBinLists=[]
        self.percentDE=[]
        self.percentDI=[]
        self.coverageLength=1
        self.nucPercent=1
        self.firstHalfExClustCount=0
        self.secondHalfExClustCount=0
        self.intronClustCount=0
        self.goodCount=0
        self.badCount=0
        self.removed=False
        self.junctionSums=[]
        
    def __str__(self):
        """
        Scripted type conversion to string.
        Returns: 
            String representation of this object.
        """
        toPrint="{ClusterGene Object} Gene: "+self.gene+"\tChromosome: "+self.chrom+"\tNumber of Clusters: "+str(len(self.lstOfClusters))
        return(toPrint)
           
    def addAC(self,toAdd):
        """
        Add a cluster object to this clusterGene
        Inputs:
            toAdd(cluster) - Cluster object to add to the ClusterGene
        """
        self.lstOfClusters.append(toAdd)
    
    def getTotalNumReads(self):
        """
        Calculates the total number of reads aligned to all clusters aligned to this gene.
        Returns:
            int (number of reads)
        """
        toRet=0
        for indivClust in self.lstOfClusters:
            toRet=toRet+indivClust.RC
        return(toRet)    
    
    def doesGeneHaveNClust(self,nClust):
        """
        Checks if this gene has at least certain number of clusters. Used for filtering processes.
        Input:
            nClust(int) - The minimum number of clusters that must be aligned to this gene in order for this method to return true.
        Returns:
            bool (are there at least this number of clusters aligned to this particular gene?)
        """
        numC = len(self.lstOfClusters)
        if numC>=nClust:
            return(True)
        return(False)

    def populateNumClust(self,anotFile): 
        """
        Get the number of clusters that map to the first half of the exonic region of this gene and set this value to the internal count properties of the object.
        Input:
            anotFile(anotFile) - Annotation file used to map the clusters to a genome.
        """
        firstHalfT=False
        secondHalfT=False
        intronT=False
        firstHalf=0
        secondHalf=0 
        intronC=0
        matchedMG = self.getMatchingMeta(anotFile)
        if matchedMG==False or self.gene=="":
            return()
        #Loop over each cluster in this gene
        for clust in self.lstOfClusters:
            #Establish the half value (anything greater than this is part of the 3'Exon half and anything equal or less is part of the 5'exon half)
            #Find the start of the region in question (the start of the exonic cluster that contains the cluster in question)
            regStart=0
            regStop=1
            for anotRow in matchedMG.lstOfAnotRows:
                if anotRow.coversIndex(int(clust.start))==True or anotRow.coversIndex(int(clust.end))==True:
                    if anotRow.ty=="exon":
                        regStart=int(anotRow.start)
                        regStop=int(anotRow.stop)
            
            #Find the end of the region in question
            halfValue = round((regStop-regStart)/2)+regStart
            parseStart=int(clust.start)
            parseStop=int(clust.end)
            
            if clust.strand=="-":
                parseStart=regStop-int(clust.end)
                parseStop=(int(clust.end)-int(clust.start))+parseStart
                
            #Go through the indices of the cluster
            for i in range(parseStart,parseStop):
                #Check if it is exonic
                if matchedMG.getType(i)=="exon":
                    #Check if it is in the first half or the second half
                    if i<=halfValue:
                        firstHalfT=True
                    if i>halfValue:
                        secondHalfT=True
                
                if matchedMG.getType(i)=="intron":
                    intronT=True
                
            if firstHalfT==True:
                firstHalfT=False
                firstHalf=firstHalf+1
            if secondHalfT==True:
                secondHalfT=False
                secondHalf=secondHalf+1
            if intronT==True:
                intronT=False
                intronC=intronC+1
            
        self.firstHalfExClustCount=firstHalf
        self.secondHalfExClustCount=secondHalf
        self.intronClustCount=intronC
                    
    def containsIndex(self,index):
        """
        Checks if a nucleotide index is within one of the clusters of the clusterGene object
        Inputs:
            index(int) - Index to check
        Returns:
            bool - Is that index within one of the clusters in this object?
        """ 
        #For each cluster
        for clust in self.lstOfClusters:
            #Check if this index is within the bounds
            if (index>=int(clust.start) and int(clust.end)>=index):
                return(True)
        return(False)
        
    def removeGene(self):
        """
        Tag a gene for removal from consideration - ie, it will no longer be considered by the program going forward.
        """
        self.removed=True
    
    def getMatchingMeta(self,anotFile):
        """
        This function finds the matching gene from a anotFile that coresponds with this clusterGene (returned)
        Inputs:
            anotFile(anotFile) - An anotFile object containing the genes to check
        Returns:
            anotGene - The anotGene in the annotation file which matches this clusterGene object.
        """
        #Get the chromosome bin to search
        chromBin = self.chrom
        for curChrom in anotFile.lstOfChromosomesGenes:
            if curChrom.ID==chromBin:
                #Search this bin for the corresponding gene
                for mg in curChrom.lstOfAnotElements:
                    if (mg.gene==self.gene):
                        return(mg)
        return(False)
                       
    def populateHLTotals(self):
        """
        Populate the list of values which represent the total number of hits in a given intron / exon pairing. This value will be used to normalize the values of all E/I/E junctions within a given gene so that they can be averaged to create gene-level average distributions (which are then averaged for the metagene)
        """
        for i in range(0,len(self.hitListsExon)):
            su1 = sum(self.hitListsExon[i])
            su2 = sum(self.hitListsIntron[i])
            theSum = su1+su2
            self.hitSums.append(theSum)
    
    def populateHitLists(self,anotFile):
        """
        Populate the hit lists of the clusterGene (each list corresponds to a EIE hit list, either intronic or exonic). In otherwords, pass stepwise over each nucleotide in this gene and see if it is overlapped by a binding cluster (asign a value of 1) or not overlapped (assign a value of 0).
        This is the foundation for all downstream analysis. Populates the hit list properties of the gene object.
        Input:
            anotFile(anotFile) - An anotFile object containing the genes to check
        """
        #Get the matching metagene
        matchedMG = self.getMatchingMeta(anotFile)
        if matchedMG==False or self.gene=="":
            return()
        matchedMG.getStartandStopIndices()
        
        if matchedMG==False and self.gene!="" :
            print("Cant find: ",self.gene)

        theCatagories=["exon","intron"]
        for gp in theCatagories:
            if (matchedMG!=False):
                #Establish the start and stop of the specific region
                if gp=="intron":
                    endVal = len(matchedMG.intronIndexStarts)
                if gp=="exon":
                    endVal = len(matchedMG.exonIndexStarts)

                if (len(matchedMG.intronIndexStarts)==0 or len(matchedMG.exonIndexStarts)==0):
                    break

                #!!!!!!!!!
                if (2*len(matchedMG.intronIndexStarts)!=len(matchedMG.exonIndexStarts)):
                    print("Repairing for matching error")
                    print("Matched Gene: ",matchedMG.gene)
                    matchedMG.repairGaps()
                    matchedMG.getStartandStopIndices()
                
                if (2*len(matchedMG.intronIndexStarts)!=len(matchedMG.exonIndexStarts)):
                    print("-----------GTF Matching error-----------")
                    print("Matched Gene: ",matchedMG.gene)
                    print(matchedMG.intronIndexStarts)
                    print(matchedMG.exonIndexStarts)
                    break
                
                secondExon=False
                
                #POPULATE EACH INDIVIDUAL EIE - taking two start and stops for the exons and and one for the introns       
                for i in range(0,endVal):
                    hitListIntron=[]
                    hitListExon=[]
                #----------------------------                
                    self.start=matchedMG.getStart()
                
                    if (gp=="exon"): 
                        theStart=matchedMG.exonIndexStarts[i] + self.start
                        theStop=matchedMG.exonIndexStops[i] + self.start+1
                
                    if (gp=="intron"):
                        theStart=matchedMG.intronIndexStarts[i] + self.start
                        theStop=matchedMG.intronIndexStops[i] + self.start+1
                    
                    for i in range(int(theStart),theStop):
                        #We need to check if there is a cluster there
                        putAZero=True
                        if self.containsIndex(i)==True:
                            #Now that we know there is a cluster, we need to check to see if it is identified as a specific coding region
                            if (matchedMG.getType(i)==gp):
                                #Also check that the exon is adjacent to an intron at this point
                                    #ALSO CHECKS TO ENSURE THAT IT IS ADJACENT TO AN INTRON
                                if(gp=="exon" and matchedMG.exonPositionIntronAdjacent(i)==True):
                                    hitListExon.append(1)
                                    putAZero=False
                                if(gp=="intron"):
                                    hitListIntron.append(1)
                                    putAZero=False                                                
                    
                        if(putAZero==True):
                            if(gp=="exon"):
                                hitListExon.append(0)
                            if(gp=="intron"):
                                hitListIntron.append(0)
                                
                    #If this is an intron we are going to add it as a seperate hit list to the list of hit lists
                    if gp=="intron":
                        self.hitListsIntron.append(hitListIntron)
                        hitListIntron=[]
                                                     
                        
                    #If this is an exon we add it or update the counter
                    if gp=="exon":
                        if secondExon==True:
                            self.hitListsExon.append(hitListExon)
                            hitListExon=[]
                            secondExon=False
                            
                            
                        elif secondExon==False:
                            secondExon=True

        if (len(self.hitListsIntron)!=len(self.hitListsExon)):
            self.hitListIntron=[]
            self.hitListExon=[]

    def checkSumTo100(self):
        """
        Checks whether every exon/intron bin list combination (each indivdual exon/intron/exon junction) adds up to either 0 (not considered) or 100 (considered). Reports outliers by printing a warning message to console.
        """
        for i in range(0,len(self.lstOfExonBinLists)):
            #Get the sum for the intron
            theIntronSum=sum(self.lstOfIntronBinLists[i])
            #Get the sum for the exon
            theExonSum=sum(self.lstOfExonBinLists[i])
            total=theIntronSum+theExonSum
            
            if (round(total)!=100 and len(self.lstOfIntronBinLists[i])>1 and len(self.lstOfExonBinLists[i])>1):
                print("NOT REACHING 100: ",self.gene)
                intSum=0
                exSum=0
                for i in range(0,len(self.hitListsExon)):
                    intSum=intSum+sum(self.hitListsIntron[i])
                    exSum=exSum+sum(self.hitListsExon[i])

    def populateIntronBinsOfPercents(self,anotFile,numberOfBins,randomStatesList):
        """
        This function uses the hit lists to populate the list of lists that contain the percent distribution for each intron hit list seperated into a given number of bins.
        This is a critical normalization step that takes the discussion from the nucelotide level (which does not allow ready comparison between genes because of the vast differences in gene length) to a percentage discussion where all lists are the same length and can be readily compared.
        Results are set to the percentile properties of this object.
        Input:
            anotFile(anotFile) - An anotFile object containing the genes to check
            numberOfBins (int) - The number of bins into which the hit list should be spread. To choose this number, consider the average difference in length between introns and exons. Selecting a bin number that matches this difference will create a smoother curve across junctions and make comparison easier.
            randomStatesList (List[int]) - The list of random states to use in the random step where rounding must be taken into account when building the bins. The longer this list, the more times the random step will be performed before averaging. 
        """
        matchedGene = self.getMatchingMeta(anotFile)
        if matchedGene==False or self.gene=="":
            return()

        if (matchedGene==False):
            return(False)
        matchedGene.getStartandStopIndices()
        
        for h in range(0,len(self.hitListsIntron)):
            
            if self.hitSums[h]==0:
                self.lstOfIntronBinLists.append([0])
                continue
            
            #Apply loop to stabilize the heuristic algorithm and ensure an end result that is indicative of the biology.
            indivPDList=[]
            for indivSeed in randomStatesList:      
                #Determine the percent value of each indivual hit 
                nucPerc = 100.0/self.hitSums[h]
                
                #Get the hit list
                theHitList=self.hitListsIntron[h]
                
                theBinList=[]
                
                #Establish the lengths to parse of the hit list 
                totLength = len(theHitList)
                  
                if totLength==0:
                    self.lstOfIntronBinLists.append([0])
                
            
                if totLength<numberOfBins:
                    #Binsize
                    lBin=numberOfBins//totLength
                    #Figure out how much is left over
                    leftOver = numberOfBins-lBin*totLength
                    randIndices = getRandLstBounded(leftOver,totLength,indivSeed)
                           
                    #Now we need to parse over the indices of the cluster            
                    addOne=False
                    indexInHitsOld=0
                    indexInHits=0
                    while indexInHits<totLength:
                        #If this index is in the list [its +1 from the bin size]
                        if indexInHitsOld in randIndices:
                            addOne=True
                
                        val = theHitList[indexInHits]
                
                        #Add the values
                        for i in range(0,lBin):    
                            theBinList.append(val*nucPerc)                     
                    
                        if (addOne==True):
                            theBinList.append(val*nucPerc)
            
                        indexInHits=indexInHits+1
                        indexInHitsOld=indexInHitsOld+1
                        addOne=False
                    
                    #Now we need to adjust for the added amount
                        #What it should equal
                        shouldEqual=1
                        shouldEqual=sum(theHitList)*nucPerc
                        
                        #What it currently equals
                        curEquals=1
                        curEquals=sum(theBinList)
                                
                        #Find the adjusting factor
                        if curEquals!=0:
                            adjFact = shouldEqual/curEquals
                            for i in range(0,len(theBinList)):
                                newEle = theBinList[i]*adjFact
                                theBinList[i]=newEle
                              
                       
                if totLength>=numberOfBins:
                    #Determine the lengths
                    #Length of each bin before random addition
                    bLength= totLength//numberOfBins
                    #Remainder
                    rLen = totLength-(numberOfBins*bLength)
                    #Create the random numbers needed
                    randIndices = getRandLstBounded(rLen,numberOfBins,indivSeed)
            
                    #Now we need to parse over the indices of the cluster            
                    nBins=0
                    indexInHits=0
                    while nBins!=numberOfBins:
                        #If this index is in the list [its +1 from the bin size]
                        if nBins in randIndices:
                            upperCount=1+bLength

                        else:
                            upperCount=bLength
            
                        #Get the total val which represents the number of hits in this region (how many nucleotides of each type fall in the region)
                        totVal=0

                        for w in range(indexInHits,indexInHits+upperCount):
                            totVal = totVal + theHitList[w]  
                    
                        #Get the percent value (multiply by the percent represented by each in idividual nucleotide)
                        theBinList.append(totVal*nucPerc)
                                                            
                        nBins=nBins+1
                        indexInHits=indexInHits+upperCount
                
                
                #If we are dealing with a reverse transcript we will need to inverse the percentile terms
                tester = self.lstOfClusters[0]
                orient2=tester.strand
                if orient2=="-":
                    #Flip
                    theBinList.reverse()              
                
                #Add thisindividual calculation to
                indivPDList.append(theBinList)
            #Average the loop values for the list
            DFI = pd.DataFrame(indivPDList)
            
            self.lstOfIntronBinLists.append(list(DFI.mean(axis = 0)))

    def populateExonBinsOfPercents(self,anotFile,numberOfBins,randomStatesList):
        """
        This function uses the hit lists to populate the list of lists that contain the percent distribution for each exon hit list seperated into a given number of bins.
        This is a critical normalization step that takes the discussion from the nucelotide level (which does not allow ready comparison between genes because of the vast differences in gene length) to a percentage discussion where all lists are the same length and can be readily compared.
        Results are set to the percentile properties of this object.
        Input:
            anotFile(anotFile) - An anotFile object containing the genes to check
            numberOfBins (int) - The number of bins into which the hit list should be spread. To choose this number, consider the average difference in length between introns and exons. Selecting a bin number that matches this difference will create a smoother curve across junctions and make comparison easier.
            randomStatesList (List[int]) - The list of random states to use in the random step where rounding must be taken into account when building the bins. The longer this list, the more times the random step will be performed before averaging. 
        """
        matchedGene = self.getMatchingMeta(anotFile)
        if matchedGene==False or self.gene=="":
            return()

        if (matchedGene==False):
            return(False)
        matchedGene.getStartandStopIndices()
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for h in range(0,len(self.hitListsExon)):
            if self.hitSums[h]==0:
                self.lstOfExonBinLists.append([0])
                continue
            
            indivPDList=[]
            for indivSeed in randomStatesList:
                #Determine the percent value of each indivual hit 
                nucPerc = 100.0/self.hitSums[h]
                
                #Get the hit list
                theHitList=self.hitListsExon[h]
                theBinList=[]
            
                #Establish the lengths to parse of the hit list 
                totLength = len(theHitList)
                
                
                if totLength<numberOfBins:
                    #Binsize
                    #OLD
                    lBin=numberOfBins//totLength
                    #NEW
                    #Figure out how much is left over
                    leftOver = numberOfBins-lBin*totLength
                    randIndices = getRandLstBounded(leftOver,totLength,indivSeed)
                    
                    #Now we need to parse over the indices of the cluster            
                    addOne=False
                    indexInHitsOld=0
                    indexInHits=0
                    while indexInHits<totLength:
                        #If this index is in the list [its +1 from the bin size]
                        if indexInHitsOld in randIndices:
                            addOne=True
                
                        val = theHitList[indexInHits]
                
                        #Add the values
                        #NEW
                        for i in range(0,lBin):
                            theBinList.append(val*nucPerc)                     
                    
                        if (addOne==True):
                            theBinList.append(val*nucPerc)
            
                        addOne=False
                    
                    #Now we need to adjust for the added amount
                        #What it should equal
                        shouldEqual=sum(theHitList)*nucPerc
                        
                        #What it currently equals
                        curEquals=sum(theBinList)
                                                                        
                        #Find the adjusting factor
                        if curEquals!=0:
                            adjFact = shouldEqual/curEquals
                            for i in range(0,len(theBinList)):
                                newEle = theBinList[i]*adjFact
                                theBinList[i]=newEle
                                
                                
                        indexInHits=indexInHits+1
                        indexInHitsOld=indexInHitsOld+1           
                
                if totLength>=numberOfBins:
                    #Determine the lengths
                    rawLength=float(totLength)/float(numberOfBins)
                    #Length of each bin before random addition
                    bLength= totLength//numberOfBins
                    #Remainder
                    rLen = round((rawLength-bLength)*numberOfBins)
                    #Create the random numbers needed
                    randIndices = getRandLstBounded(rLen,numberOfBins,indivSeed)

                    #Now we need to parse over the indices of the cluster            
                    nBins=0
                    indexInHits=0
                    while nBins!=numberOfBins:
                        #If this index is in the list [its +1 from the bin size]
                        if nBins in randIndices:
                            upperCount=1+bLength

                        else:
                            upperCount=bLength
            
                        #Get the total val which represents the number of hits in this region (how many nucleotides of each type fall in the region)
                        totVal=0

                        for w in range(indexInHits,indexInHits+upperCount):
                            totVal = totVal + theHitList[w]  
                            #file2.write("Value at hit list index: "+str(theHitList[w])+"\n")
                    
                        #Get the percent value (multiply by the percent represented by each in idividual nucleotide)
                        theBinList.append(totVal*nucPerc)
                                                            
                        nBins=nBins+1
                        indexInHits=indexInHits+upperCount
                
                #If we are dealing with a reverse transcript we will need to inverse the percentile terms
                tester = self.lstOfClusters[0]
                orient2=tester.strand
                if orient2=="-":
                    #Flip
                    theBinList.reverse()
                
                indivPDList.append(theBinList)
                
            #Average the loop values for the list
            DFE = pd.DataFrame(indivPDList)

            self.lstOfExonBinLists.append(DFE.mean(axis = 0))

    def populateGeneDistributions(self,NIntronBins,NExonBins):
        """
        Averages the percentage lists for each EIE in this gene in order to get the average intron/exon distribution of clusters across this gene.
        It ONLY considers those binned lists which are greater than 1 (smaller values indicate errors that were set to 1 earlier in the process)
        Sets results to percent properties of the gene object.
        Input:
            NIntronBins (Int) - The number of intronic bins into which binding density should be allocated.
            NExonBins (Int) - The number of intronic bins into which binding density should be allocated.
        """
        #If the hit list sum for this gene is 0 (there is no overlap) simply set the distributions to [0] and exit
        if sum(self.hitSums)==0:
            self.percentDE=[0]
            self.percentDI=[0]
            return            
        
        intronicBinsUsed=[]
        exonicBinsUsed=[]
        
            #INTRONS
        #Loop over each bin 
        for i in range(0,NIntronBins):
            IntronEntries = []
            #For each EIE
            for binList in self.lstOfIntronBinLists:
                #If this list is not of the prescribed length, ditch it
                if len(binList)!=NIntronBins:
                    continue
                else:
                    if self.gene=="STRIP1" and binList not in intronicBinsUsed:
                        intronicBinsUsed.append(binList)
            #Otherwise add this to the entry list
                IntronEntries.append(binList[i])
                
            #Average the values from all of the genes at this bin locationF
            theVal = Average(IntronEntries)

            #Add this average to the gene object's distribution list
            self.percentDI.append(theVal)            
        
                    #Exons
        #Loop over each bin 
        for i in range(0,NExonBins):
            ExonEntries = []
            #For each EIE
            for binList in self.lstOfExonBinLists:
                #If this list is not of the prescribed length, ditch it
                if len(binList)!=NExonBins:
                    continue
                else:
                    if self.gene=="STRIP1" and list(binList) not in exonicBinsUsed:
                        exonicBinsUsed.append(list(binList))
            #Otherwise add this to the entry list
                ExonEntries.append(binList[i])
            
            #Average the values from all of the genes at this bin location
            theVal = Average(ExonEntries)
            #Add this average to the gene object's distribution list
            self.percentDE.append(theVal)
        
        if sum(self.percentDI)==0:
            self.percentDI=[0]
 
        if sum(self.percentDE)==0:
            self.percentDE=[0]
                  
            
class parclip:
    """
    Parclip Object: This object represents the output from an indivual parclip clusters.csv file.

    Properties:
        filename(String) - The name of the csv file.
        lstOfChromosomesClusters(List[Chromosomes]) - A list of the clusters in a parclip, sorted by chromosome.
        lstOfChromosomesClustersGenes (List[Chromosomes]) - A list of the cluster genes in a parclip, sorted by chromosome.
        iDistribution (List[float]) - The average distribution of binding density across the intronic region of all genes in the parclip.
        eDistribution (List[float]) - The average distribution of binding density across the exonic regions of all genes in the parclip.
        eFHDistribution (List[float]) - The average distribution of binding density across the first (5') half of all exonic regions of all genes in the parclip.
        eSHDistribution (List[float]) - The average distribution of binding density across the second (3') half of all exonic regions of all genes in the parclip.
        allCluster (List[cluster]) - A simple list of all of the clusters in this parclip object (unsorted).
    
    Methods:
        init - Initializer
        str - Conversion to string form for printing.
        addClust - Add a cluster to the parclip object (used when reading in the object).
        getNumClusters - Determines how many clusters are present in a parclip object.
        getNumGenes - Determines how many  genes are present in a parclip object.
        removeCluster - Given a cluster, this method goes through the lstOfChromsomesClusters property and removes the inputed cluster.
        applyBounds - Impose bounds on the clusters included - throwing out extremes (based on inputed property) of a specified % at the minimum and maximum end of the spectrum.
        getNumExonGenesH1 - Determines the number of genes that have cluster overlap within the 5' half of exons in this PAR-CLIP.
        getNumExonGenesH2 - Determines the number of genes that have cluster overlap within the 3' half of exons in this PAR-CLIP.
        getNumExonClustersH1 - Determines the number of clusters that fall within the 5' half of the exons in genes of this PAR-CLIP.
        getNumExonClustersH2 - Determines the number of clusters that fall within the 3' half of the exons in genes of this PAR-CLIP.
        getNumIntronGenes - Determines the number of genes that have cluster overlap within intronic regions in this PAR-CLIP.
        populateCh - Populate the chromosome lists of this object with the clusters it holds (a sorted storage step).
        populateAllHits - Populates the hit lists (and percentages of cluster distribution) for all genes in this object.
        populateExonHalves - populates the exon half lists for use later in graphing.
        
    """
    def __init__(self,filename,lstOfUniqueChromosomes):
        self.filename=filename
        theChromsomeList=[]
        for chrN in lstOfUniqueChromosomes: 
            theChromsomeList.append(chromesome(chrN))

        theChromsomeList2=[]
        for chrN in lstOfUniqueChromosomes: 
            theChromsomeList2.append(chromesome(chrN))
        
        self.lstOfChromosomesClusters=theChromsomeList
        self.lstOfChromosomesClustersGenes=theChromsomeList2
        
        self.iDistribution=[]
        self.eDistribution=[]
        self.eFHDistribution=[]
        self.eSHDistribution=[]                  
                
    def __str__(self):
        toPrint="{Parclip Object} Name: "+self.filename+"\tNumber of Clusters: "+str(self.getNumClusters())+"\tNumber of Cluster Genes: "+str(self.getNumGenes())
        return(toPrint)

    def removeClust(self, cluster):
        """
        Removes a given cluster from the lstOfChromosomesClusters value.
        Inputs:
            cluster(cluster) - The cluster to remove from the parclip object.
        """
        #Get the appropriate chromosome of this cluster object
        toRemChrom=cluster.chrom
        for chrom in self.lstOfChromosomesClusters:
            if (chrom.ID==toRemChrom):
                for clustInChrom in chrom.lstOfAnotElements:
                    if clustInChrom==cluster:
                        chrom.lstOfAnotElements.remove(cluster)
                        return() 

    def addClust (self,cluster):
        """
        Adds a given cluster to the parclip object.
        Inputs:
            cluster(cluster) - The cluster to add to the parclip object.
        """
        #Get the correct chromosome
        curID = cluster.chrom
        for chrom in self.lstOfChromosomesClusters:
            if (chrom.ID==curID):
                chrom.lstOfAnotElements.append(cluster)
                return()

    def getNumClusters(self):
        """
        Returns the number of clusters within this parclip object.
        Returns:
            int - The total count of clusters contained in this parclip object.
        """
        toReturn=0
        for ch in self.lstOfChromosomesClusters:
            toReturn=toReturn+len(ch.lstOfAnotElements)
        return(toReturn) 
        
    def getNumGenes(self):
        """
        Returns the number of genes that are overlapped by clusters within this parclip object.
        Returns:
            int - The total count of genes that are overlapped by clusters in this parclip object.
        """
        toReturn=0
        for ch in self.lstOfChromosomesClustersGenes:
            toReturn=toReturn+len(ch.lstOfAnotElements)
        return(toReturn)    
        
    def getNumExonGenesH1(self):
        """
        Determines the number of genes in this PAR-CLIP that have cluster overlap in their 5' half of thier exons.
        Returns:
            int - The number of genes.
        """
        toReturn=0
        for ch in self.lstOfChromosomesClustersGenes:
            for gene in ch.lstOfAnotElements:
                if gene.removed==True:
                    continue
                #Get the first half of the exon
                endInd=round(len(gene.percentDE)/2)
                fh=gene.percentDE[0:int(endInd)]
                if sum(fh) >0:
                    toReturn=toReturn+1
        return(toReturn)  
        
    def getNumExonGenesH2(self):
        """
        Determines the number of genes in this PAR-CLIP that have cluster overlap in their 3' half of thier exons.
        Returns:
            int - The number of genes.
        """
        toReturn=0
        for ch in self.lstOfChromosomesClustersGenes:
            for gene in ch.lstOfAnotElements:
                if gene.removed==True:
                    continue
                #Get the first half of the exon
                startInd=round(len(gene.percentDE)/2)
                sh=gene.percentDE[int(startInd):len(gene.percentDE)]
                if sum(sh) >0:
                    toReturn=toReturn+1
        return(toReturn)  
    
    def getNumExonClustersH1(self): 
        """
        Determines the number of clusters that fall within the 5' half of the exons in genes of this PAR-CLIP.
        Returns:
            int - The number of clusters.
        """
        toRet=0
        for ch in self.lstOfChromosomesClustersGenes:
            for gene in ch.lstOfAnotElements:
                if gene.removed==True:
                    continue
                toRet = toRet+ gene.firstHalfExClustCount
        return(toRet)
        
    def getNumExonClustersH2(self): 
        """
        Determines the number of clusters that fall within the 3' half of the exons in genes of this PAR-CLIP.
        Returns:
            int - The number of clusters.
        """
        toRet=0
        for ch in self.lstOfChromosomesClustersGenes:
            for gene in ch.lstOfAnotElements:
                if gene.removed==True:
                    continue
                toRet = toRet+ gene.secondHalfExClustCount
        return(toRet)
        
    def getNumIntronGenes(self):
        """
        Determines the number of genes in this PAR-CLIP with cluster density in their intronic regions.
        Returns:
            int - The number of clusters.
        """
        toReturn=0
        for ch in self.lstOfChromosomesClustersGenes:
            for gene in ch.lstOfAnotElements:
                if len(gene.percentDI)==0:
                    continue
                if sum(gene.percentDI) >0:
                    toReturn=toReturn+1
        return(toReturn)  
    
    def getNumIntronClusters(self): 
        """
        Determines the number of clusters that fall within intronic regions of genes in this PAR-CLIP.
        Returns:
            int - The number of clusters.
        """
        toRet=0
        for ch in self.lstOfChromosomesClustersGenes:
            for gene in ch.lstOfAnotElements:
                if gene.removed==True:
                    continue
                toRet = toRet+ gene.intronClustCount
        return(toRet)

    def populateCh(self):
        """
        Populate the chromosome lists with the clusters and clustergenes. You can think of this as a sorting/storage step designed to hasten downstream comparisons.
        """
        #for each cluster
        for i in range(0,len(self.lstOfChromosomesClusters)):
            cChr = self.lstOfChromosomesClusters[i]
            gChr= self.lstOfChromosomesClustersGenes[i]
                
            #For each row in the chromosome
            for indivClust in cChr.lstOfAnotElements:
                counter=0
                #Check if this is in the corresponding chromosome for the gene list
                for indivGene in gChr.lstOfAnotElements:
                    if (indivGene.gene==indivClust.gene):
                        indivGene.addAC(indivClust)
                        counter=1
                        break
                if counter==0:
                    #Make a new anotGene
                    nClustGene=clusterGene(indivClust.gene,indivClust.chrom)
                    nClustGene.addAC(indivClust)
                    gChr.lstOfAnotElements.append(nClustGene)   
        
    def populateAllHits(self,anotF,NIntronBins,NExonBins,randomStatesList):
        """
        Populate the hit lists for the intronic and exonic regions of all genes in this PAR-CLIP. In other words, this part of the program determines the nucleotide-level overlap with binding clusters before re-expressing these lists are binned for comparison later. This binning sometiems require dealing with uneven numbers. When this occurs, random bins are selected to be 1 larger than the average. This is done multiple times and the results averaged.
        Input:
            anotF (anotFile) - Annotation file that contains the information needed to construct the metagene.
            NIntronBins (Int) - THe number of bins to distribute the intronic hit lists into.
            NExonBins (Int) - THe number of bins to distribute the exonic hit lists into.
            randomStatesList (List[int]) - The list of numbers to use as random seeds in random number generation. The more random states provided, the more times the binning algorithm will execute before the averaging step.
        """
        #For each chromosome
        for chrom in self.lstOfChromosomesClustersGenes:
            #Get each gene
            for cGene in chrom.lstOfAnotElements:
                #Populate              
                cGene.populateHitLists(anotF)
                #indivHL = cGene.checkHitList()
                cGene.populateHLTotals()
                cGene.populateIntronBinsOfPercents(anotF,NIntronBins,randomStatesList)
                cGene.populateExonBinsOfPercents(anotF,NExonBins,randomStatesList)
                cGene.populateNumClust(anotF)
                cGene.checkSumTo100()
                cGene.populateGeneDistributions(NIntronBins,NExonBins)
        
    def getParclipPercentageDistributions(self,NIntronBins,NExonBins,clusterThreshold,weighValuesByRC,wtGE):
        """
        Populate the average distribution of all of the percentage hits in this parclip object. In essence, just average the distribution lists of all of the individual genes in the PAR-CLIP object.
        Input:
            NIntronBins (Int) - The number of intronic bins the hits were divided into.
            NExonBins (Int) - The number of exonic bins the hits were divided into.
            clusterThreshold (Int) - The minimum number of clusters a gene must have tobe used.
            weighValuesByRC (bool) - Whether or not we want to weight the genes with more aligned reads to be more impactful than those with fewer.
            wtGE (WTGeneExpression or False) - Whether to weight the impact of each gene on the metagene by the expression level of that gene in the wild type cell line.
        """    
        #Remove all of those lists which do not sum to 100 with their partner and ensure that those with hits in only the exonic or intronic regions are still counted.
        for chrom in self.lstOfChromosomesClustersGenes:
            #Get each gene
            for cGene in chrom.lstOfAnotElements:

                    #IF the gene has fewer than the specified number of clusters set its lists equal to zero
                if cGene.doesGeneHaveNClust(clusterThreshold)==False:
                    cGene.percentDI=[]
                    cGene.percentDE=[]
                    cGene.removeGene()
                
                    #IF the gene has only introns or only exons fix it so that it is also considered
                if len(cGene.percentDI)!=NIntronBins and len(cGene.percentDE)==NExonBins:
                    cGene.percentDI=[0]*NIntronBins
                    
                if len(cGene.percentDI)==NIntronBins and len(cGene.percentDE)!=NExonBins:
                    cGene.percentDE=[0]*NExonBins
                        
                if (round(sum(cGene.percentDI)+sum(cGene.percentDE))) != 100:
                    if (round(sum(cGene.percentDI)+sum(cGene.percentDE))) !=0 and len(cGene.percentDI)>0:
                        print("Not reaching 100 or 0: ",cGene.gene)
                        cGene.percentDI=[0]
                        cGene.percentDE=[0]
                
            #INTRONS
        #Loop over each bin 
        for i in range(0,NIntronBins):
            IntronEntries = []
            #For each gene
            for chrom in self.lstOfChromosomesClustersGenes:
                #Get each gene
                for cGene in chrom.lstOfAnotElements:
                    if round((sum(cGene.percentDI)+sum(cGene.percentDE))) != 100:
                        continue
                    #If that gene has an intron distribution list equal to the number of bins...                            
                    if len(cGene.percentDI)==NIntronBins:
                        #Add the value to the options
                        #If we are weighing by RC
                        if weighValuesByRC==True and wtGE==False:
                            for k in range(0,cGene.getTotalNumReads()):
                                IntronEntries.append(cGene.percentDI[i])
                        
                        if weighValuesByRC==False and wtGE==False:
                            IntronEntries.append(cGene.percentDI[i])
                        
                            #If we are weighing by BOTH wild type gene expression AND read count
                        if wtGE!=False:
                            #Get the value for the TPM expression
                            if cGene.gene in wtGE.geDict:
                                #Save the zero expression ones
                                if wtGE.geDict[cGene.gene]!=0:
                                    countUp=int(round(wtGE.geDict[cGene.gene]))
                                    if countUp!=0:
                                        #Get the factor to multiply each by (#CLR/TPM)
                                        multFact=cGene.getTotalNumReads()/countUp
                                        IntronEntries.append(cGene.percentDI[i]*multFact)
                                    
            #Get the averageclear
            avgVal = Average(IntronEntries)
            self.iDistribution.append(avgVal)
       
                    #ExONS
        #Loop over each bin 
        for i in range(0,NExonBins):
            ExonEntries = []
            #For each gene
            for chrom in self.lstOfChromosomesClustersGenes:
                #Get each gene
                for cGene in chrom.lstOfAnotElements:
                    if round((sum(cGene.percentDI)+sum(cGene.percentDE))) != 100:
                        continue
                    #If that gene has an intron distribution list equal to the number of bins
                    if len(cGene.percentDE)==NExonBins:
                        #Add the value to the options
                        if weighValuesByRC==True and wtGE==False:
                            for k in range(0,cGene.getTotalNumReads()):
                                ExonEntries.append(cGene.percentDE[i])
                        if weighValuesByRC==False and wtGE==False:
                            ExonEntries.append(cGene.percentDE[i])
                        
                        #If we are weighing by both wild type gene expression AND read count
                        if wtGE!=False:
                            #Get the value for the TPM expression
                            if cGene.gene in wtGE.geDict:
                                #Save the zero expression ones
                                if wtGE.geDict[cGene.gene]!=0:
                                    countUp=int(round(wtGE.geDict[cGene.gene]))
                                    if countUp!=0:
                                        #Get the factor to multiply each by (#CLR/TPM)
                                        multFact=cGene.getTotalNumReads()/countUp
                                        ExonEntries.append(cGene.percentDE[i]*multFact)    
                        
            #Get the average
            avgVal = Average(ExonEntries)
            self.eDistribution.append(avgVal)
            
        print("---------------- Final Intron Distribution: ",self.iDistribution)
        print("---------------- Final IExon Distribution: ",self.eDistribution)
        print("Length of intron final distribution: ",len(self.iDistribution))
        print("Length of exon final distribution: ",len(self.eDistribution))
        print("AVG of Exon",Average(self.eDistribution))
        print("AVG of Intron",Average(self.iDistribution))
        print("Total sum: ",sum(self.eDistribution)+sum(self.iDistribution))

    def populateExonHalves(self,NExonBins):
        """
        Populates the two half lists of the exonic distribution in order to allow them to be split for graphing later.
        Input:
            NExonBins (int) - The number of bins the exonic hit lists were compressed or stretched into.
        """
        if len(self.eDistribution)!=0:
        #Strip the first value
            self.eDistribution.pop(0)
            self.eDistribution.pop(0)
            #Strip the last value
            self.eDistribution.pop()
            self.eDistribution.pop()
        
        
        #Get the value of half of the number of bins
        hnb=(NExonBins-4)//2
        #Initiate the lists
        firstHalf = []
        secondHalf =[]
        
        for i in range(0,hnb):
            toAdd=self.eDistribution[i]
            firstHalf.append(toAdd)
            
        
        for i in range(hnb,(NExonBins-4)):
            toAdd=self.eDistribution[i]
            secondHalf.append(toAdd)
    
        self.eFHDistribution=firstHalf
        self.eSHDistribution=secondHalf
        
        print("-----5' Exon Half: ",self.eFHDistribution)
        print("-----3' Exon Half: ",self.eSHDistribution)
    
    def applyBounds(self,lowerBound,upperBound,field):
        """
        This method removes clusters from the PAR-CLIP object that are at either % end of the sorted spectrum of all clusters for an inputed property (read count, cross linked reads, etc). Rounds down through Int() typecasting at the final stage.
        Inputs:
            lowerBound (Int) - The percentage of clusters on the lower end of the value spectrum to remove.
            uppBound (Int) - The percentage of clusters on the upper end of the value spectrum to remove.
            field (Str) - The cluster.csv field to use as the sorting parameter.
        """
        #Populate the allCluster property of this parclip object
        allClusters=[]
        for indivChrom in self.lstOfChromosomesClusters:
            for clust in indivChrom.lstOfAnotElements:
                allClusters.append(clust)              
        
        if (lowerBound+upperBound)==0 or field=="":
            print("No bounds and/or filter category entered to filter by. All clusters used.")
            return()
        
        #Readout for filtering results
        cCount=len(allClusters)
            
        print("Number of clusters before: ",cCount)
        #Sort all of the clusters by the field indicated
        acceptableFields=["start","end","CS","URC","RC"]
        if field not in acceptableFields:
            print("Unable to apply filtering operation because the inputed field is invalid. It must be one of:")
            print(acceptableFields)
            return()

        if field=="start":
            print("Sorting by start location for filter...")
            sortedList=sorted(allClusters, key=lambda x: int(x.start))

        if field=="end":
            print("Sorting by end location for filter.")
            sortedList=sorted(allClusters, key=lambda x: int(x.end))

        if field=="CS":
            print("Sorting by CS (conversion specificity) for filter...")
            sortedList=sorted(allClusters, key=lambda x: int(x.CS))
            
        if field=="URC":
            print("Sorting by URC (unique read count) for filter...")
            sortedList=sorted(allClusters, key=lambda x: int(x.URC))
            
        if field=="RC":
            print("Sorting by RC (read count) for filter...")
            sortedList=sorted(allClusters, key=lambda x: int(x.RC))
        
        #Figure out how many clusters constitue a percent (round down).
        clustersInAPercent=len(sortedList)/100.0
        
        #Initialize an empty list to add the clusters to remove into
        clustersToRem=[]
        
        #Now we need to figure out how many clusters to remove for the bounds.
            #Remove for the lowest bound
                #How many clusters constitute the lower bound?
        lowerBoundCount=int(clustersInAPercent*lowerBound)
            #Loop to remove
        for i in range(0,lowerBoundCount):
            clustersToRem.append(sortedList[i])

            #Remove for the upper bound
                #How many clusters constitue the upper bound?
        upperBoundCount=int(clustersInAPercent*upperBound)
        for j in range(len(sortedList)-1,len(sortedList)-upperBoundCount-1,-1):
            clustersToRem.append(sortedList[j])
        
        #Finally, remove these clusters from the chromosome properties (The stoarge mechanism used for calculations).
        for indivClust in clustersToRem:
            self.removeClust(indivClust)
        
                #Readout for filtering results
        cCount2=0
        for indivChr in self.lstOfChromosomesClusters:
            cCount2=cCount2+len(indivChr.lstOfAnotElements)
            
        print("Number of clusters after filtering by bounds: ",cCount2)


##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Functions
##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ 
def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

def isint(value):
  try:
    int(value)
    return True
  except ValueError:
    return False

def dictionaryToCluster(cDict):
    """
    Takes a dictionary and uses it to generate a cluster object with the fields that are availbile.
    Input:
        Dictionary of cluster values.
    Returns:
        cluster
    """
    #Mandatory fields
    #Chromosome
    theChr=cDict.get("Chr")
    #Strand
    theStrand=cDict.get("Strand")
    #Start
    theStart=cDict.get("Start")
    #End
    theEnd=cDict.get("End")
    
    #Necessary field - gene name - generate if missing
    theGene=""
    if "GeneName" in cDict.keys():
        theGene=cDict.get("GeneName")
    if "GeneName" not in cDict.keys(): 
        print("Gene not found. Fetching gene...!")
        
    #Optional fields - populate with empty fields if missing
    if "ConversionSpecificity" in cDict.keys():
        if isfloat(cDict.get("ConversionSpecificity"))==True:
            theCS=float(cDict.get("ConversionSpecificity"))
        else:
            theCS=None
    else:
        theCS=None
    
    if "T2Cfraction" in cDict.keys():
        if isfloat(cDict.get("T2Cfraction"))==True:                
            T2Cf=float(cDict.get("T2Cfraction"))
        else:
            T2Cf=None
    else:
        T2Cf=None
    
    if "UniqueReads" in cDict.keys():
        if isint(cDict.get("UniqueReads"))==True:
            UR=int(cDict.get("UniqueReads"))
        else:
            UR=1
    else:
        UR=None

    if "ReadCount" in cDict.keys():
        if isint(cDict.get("ReadCount")):
            RC=int(cDict.get("ReadCount"))
        else:
            RC=1
    else:
        RC=None
    
    nClust=cluster(theChr,theStrand,theStart,theEnd,theGene,theCS,T2Cf,UR,RC)
    return(nClust)
  
def createLstOfParclips (myDir,genesToUse,lowerBound,upperBound,boundFilter,gtfObj,bedGTF,anoterScript,outDir):
    """
    Creates a list of parclip objects based on an input directory of cluster.csv files.
    Input:
        myDir(String) - File path to the directory containing all of the cluster.csv files to be inputs for the program or the name of a single clsuter.csv file to analyze.
        genesToUse(List[String]) - A list of the names of genes to be used in subsequent calclulations and visualizations. If all are desired enter an empty list.
        sort(Boolean) - Should the sort step, which keeps only those genes which are in the top quartile of total read counts, be performed?
        gtfObj (gtfObj) - GTF object that will be used to analyze the parclip clusters.
        bedGTF (str) - GTF file to use to annotate bed files.
        outDir (str) - Directory to output files into
    Returns: 
        List[parclip] - A list of parclip objects representing the csv files in the directory.
    """ 
    #Set the wd to be myDir
    if ".csv" not in myDir and ".bed" not in myDir: #If a directory
        os.chdir(myDir)
    listOfParclipObjects = []
    
    #The files to be analyzed - either a single file or a directory
    filesToAnalyze=[]
    if ".csv" in myDir or ".bed" in myDir:
        filesToAnalyze.append(myDir)
    else:
        filesToAnalyze=os.listdir(myDir)
    
    #For each file in the list of files
    for i in range(0,len(filesToAnalyze)):
            
        #Read in said file
            #Decide if it should be read
        filename= filesToAnalyze[i]
        if ".csv" not in filename and ".bed" not in filename:
            continue
        
        #Determine the type of the file
        fileType=""
        if ".bed" in filename:
            fileType="bed"
        if ".csv" in filename:
            fileType="csv"
            
        if ".csv" not in myDir and ".bed" not in myDir:
            toReadForCSV = myDir+"/"+filename
        if ".csv" in myDir or ".bed" in myDir:
            toReadForCSV=filename
        
            #Read the file in
        currentList = getFile(toReadForCSV,gtfObj,bedGTF,anoterScript,outDir)
        
            #Create a parclipobject with the name
        title = filename
        nTitle=""
        if fileType=="csv":
            nTitle = title.replace('.csv', '')
        if fileType=="bed":
            nTitle = title.replace('.bed', '')

        #Establish the chromosomes names
            #Identify the unique chromosome names needed
        uniqueChromosomes=[]
        for rDict in currentList:
            curChr = rDict.get("Chr")
            if curChr not in uniqueChromosomes:
                if curChr in gtfObj.lstOfChromosomeNames:
                    uniqueChromosomes.append(str(curChr))
        
        newPC = parclip(nTitle,uniqueChromosomes)
        #Fill the parclip with the cluster rows
        #Initialize the list ot hold the clusters
        lstOfClustersOrig = []
        #Parse over each dictionary in the list
        for cDict in currentList:
            #Take the cluster if it is a gene to be used or if there are no gene restrictions
            if (cDict.get("GeneName") in genesToUse) or (len(genesToUse))==0:
                #Create the corresponding cluster object
                nClust = dictionaryToCluster(cDict)            
                #Add it to the list of possible lusters
                lstOfClustersOrig.append(nClust) 
        
        for i in range(0,len(lstOfClustersOrig)):
            clusterToAdd = lstOfClustersOrig[i]
            newPC.addClust(clusterToAdd)
        
        #Blank the genes list
        for crClustGen in newPC.lstOfChromosomesClustersGenes:
            crClustGen.lstOfAnotElements=[]  
                
        #Apply filter(s) if you have a csv file (otherwise you don't have the property needed to apply bounds)
        if fileType=="csv":
            newPC.applyBounds(lowerBound,upperBound,boundFilter)
        
        newPC.populateCh()
        
        listOfParclipObjects.append(newPC)
    
        
    return(listOfParclipObjects)
    
def getFile (myFile,anotObj,bedGTF,anoterScript,outDir):
    """
    Takes in a directory and file name and creates a list of dictionaries, each reprenting a row. Dictionary values correspond to columns and list enteries correspond to rows
    Input:
        myFile(String) - Name of the csv file.
        anotObj (anotObj) - Annotation object to be used if bed file does not have gene names.
        bedGTF (Str) - GTF file to be used to annotate beds.
        outDir(Str) - Output directory to deposit files 
    Returns:
        List[Dictionary] - A list of dictionaries, each which represents a row from the csv to read.
    """
    #If the file is a csv cluster file.
    if ".csv" in myFile:
        #Upload CSV using the CSV package and an instantiated CSV object
        readCSV = csv.DictReader(open(myFile))
        #Create the master list to be populated
        master_list = []
    
        #poulate this list by grabbing ordered dictionaries
        for line in readCSV:
            #If this is PARALUS notation we need to make a newline
            if "chromosome" in line:
                nLine={}
                nLine["Chr"]=line["chromosome"]
                nLine["Start"]=line["start"]
                nLine["End"]=line["end"]
                nLine["GeneName"]=line["gene_name"]
                nLine["Strand"]=line["strand"]
                nLine["ClusterSequence"]=line["sequence"]
                nLine["Aligned to"]=line["annotation"]
                nLine["UniqueReads"]=line["crosslinked_reads"]
                nLine["T2Cfraction"]=line["fraction"]
                nLine["ReadCount"]=line["total_reads"]
                master_list.append(nLine)
            else:
                master_list.append(line)

        return(master_list)
    
    #Otherwise, if the file is a bed file    
    elif ".bed" in myFile:
        #See if the gene name is in the input file
        testBed= pd.read_csv(open(myFile),sep='\t')
        if "GeneName" in testBed.columns:
            readCSV = csv.DictReader(open(myFile),delimiter='\t')
            #Create the master list to be populated
            master_list = []
    
            #poulate this list by grabbing ordered dictionaries
            for line in readCSV:
                master_list.append(line)
                
            return(master_list)
        
        #Otherwise we will have to annotate
        else:
            #Run the annotation program
            print("Running annotation program on bed file...")
            #Build the command
            anotCom=anoterScript+" -i "+myFile+" -G "+bedGTF
            
            os.system(anotCom)

            #New file
            nFile=myFile.replace(".bed", ".annotated.bed", 1)
            nFile2=myFile.replace(".bed", ".annotated_withHeader.bed", 1)
            nFile3=myFile.replace(".bed", ".annotated_withHeader_spaceCorrected.bed", 1)
            
            f = open(nFile2, "w")
            writer = csv.DictWriter(f, fieldnames=["Chr", "Start","End","GeneName","AlignedTo","Strand"],delimiter='\t')
            writer.writeheader()
            f.close()
            
            #Add
            open(nFile2, "a").writelines([l for l in open(nFile).readlines()])
            
            #Replace eall of the spaces with tabs
                #input file
            fin = open(nFile2, "rt")
            #output file to write the result to
            fout = open(nFile3, "wt")
            #for each line in the input file
            for line in fin:
                #Read replace the string and write to output file
                fout.write(line.replace(' ', '\t',5))
            #close input and output files
            fin.close()
            fout.close()
            
            readCSV = csv.DictReader(open(nFile3),delimiter='\t')
            #poulate this list by grabbing ordered dictionaries
            master_list = []
            for line in readCSV:
                #Adjust for the problem in the annotating script where the "aligned to" section gets swapped with the "GeneName" section somehow
                alignValues=["CDS","UTR","gene","transcript"]
                for av in alignValues:
                    if av in line["GeneName"]:
                        nGeneName=line["Strand"]
                        nAlignedTo=line["GeneName"]
                        nStrand=line["AlignedTo"]
                        
                        line["GeneName"]=nGeneName
                        line["AlignedTo"]=nAlignedTo
                        line["Strand"]=nStrand
                        
                master_list.append(line)
            
            #move the files to the output director
            moveCommand1="rm "+nFile+" "
            os.system(moveCommand1)
            
            moveCommand="rm "+nFile2+" "
            os.system(moveCommand)
            
            moveCommand="mv "+nFile3+" "+outDir
            os.system(moveCommand)
            
            return(master_list)

            
def getCSV (myFile):
    """
    Takes in a directory and file name and creates a list of dictionaries, each reprenting a row. Dictionary values correspond to columns and list enteries correspond to rows
    Input:
        myFile(String) - Name of the csv file.
    Returns:
        List[Dictionary] - A list of dictionaries, each which represents a row from the csv to read.
    """   
    #Upload CSV using the CSV package and an instantiated CSV object
    readCSV = csv.DictReader(open(myFile))
    #Create the master list to be populated
    master_list = []
    
    #poulate this list by grabbing ordered dictionaries
    for line in readCSV:
        master_list.append(line)
    return(master_list)


def readInGTFFile (gtfFile,onlyMain):
    """
    Read in the gtf file csv file, returning an anotatedFile object
    Input:
        gtfFile(String) - Name of the csv file.
        onlyMain (bool) - Whether only the main chromosomes should be considered.
    Returns:
        anotFile - An annotation file object that contains the information from the annotation csv file.
    """ 
    #Import the data
    mData = getCSV(gtfFile)
    
    #If need be, construct the list of the main chromsomes
    mainChromosomes=[]
    if onlyMain==True:
        #Establish whether the chromosome names contain "chr"
        containC=False
        uniqueChromosomes=[]
        for rDict in mData:
            curChr = rDict.get("Chromosome")
            if "chr" in curChr:
                containC=True
                break
            
        for rDict in mData:
            curChr = rDict.get("Chromosome")
            if containC==True:
                testChr=curChr[3:len(curChr)]
            else:
                testChr=curChr
            if len(testChr)<=2:
                if curChr not in mainChromosomes and "M" not in curChr:
                    mainChromosomes.append(curChr)
    
    #Identify the unique chromosome names needed
    uniqueChromosomes=[]
    for rDict in mData:
        curChr = rDict.get("Chromosome")
        
        if curChr not in uniqueChromosomes:
            if onlyMain==True:
                if curChr in mainChromosomes:
                    uniqueChromosomes.append(str(curChr))
            else: 
                uniqueChromosomes.append(str(curChr))
    
    print("Unique Chromsomes for the gtf file: ", uniqueChromosomes)
        #Instantiate object
    myAFileObject = anotFile(uniqueChromosomes)
    
    #Parse over each dictionary in the list
    for rDict in mData:
        #Adjust the chromosome field if needed
        curChr = rDict.get("Chromosome")
        if "chr" not in curChr:
            curChr="chr"+curChr
        #Create the corresponding anotation row object
        nAnotR = anotRow(rDict.get("Gene"),rDict.get("Type"),int(rDict.get("Start")),int(rDict.get("Stop")),curChr,rDict.get("Orientation"))
        #Add this to the anotation file object
        myAFileObject.addAR(nAnotR)
        
    myAFileObject.populateCh()

    return(myAFileObject)

def getRandLst(num,seed):
    """
    This function takes in a number and generates a list of that many (semi) unique random numbers between 0 and 99.
    Input:
        num(int) - The number of random numbers to generate.
        seed(num) - Random seed to set.
    Returns:
        List[int] - A list of unique random numbers between 0 and 99 of length num.
    """
    random.seed(seed)
    lstOfRandNums = []
    #Generate initial numbers
    while len(lstOfRandNums)!=num:
        nNum = random.randint(0,99)
        if nNum not in lstOfRandNums:
            lstOfRandNums.append(nNum)
    return(lstOfRandNums)    


def getRandLstBounded(num,bound,seed):
    """
    This function takes in a number and generates a list of that many (semi) unique random numbers between 0 and an inputed bound.
    Input:
        num(int) - The number of random numbers to generate.
        bound(int) - The upper bound of possible value generation.
        seed (num) - Random seed to set.
    Returns:
        List[int] - A list of unique random numbers between 0 and the inputed upper bound of length num.
    """ 
    random.seed(seed)
    lstOfRandNums = []
    #Generate initial numbers
    while len(lstOfRandNums)!=num:
        nNum = random.randint(0,bound-1)
        if nNum not in lstOfRandNums:
            lstOfRandNums.append(nNum)
    return(lstOfRandNums)  
    
    
def createExportTable(lstOfPCs,myDir,genesToUse):
    """
    Takes in a list of parclip objects and returns a table which contains the distribution matrices for each of these. 
    This strange approach was adopted before the author learned to import sophisticated packages to expedite this process.
    Inputs:
        lstOfPCS (List[parclip]) - The list of parclip objects whose properties will be pulled and then exported in the table.
        myDir (Str) - The directory into which the table will be written.
        genesToUse (List[Str]) - A list of gene names to be considered (only these genes will have their properties identified and outputed). If empty, all genes are consodered.
    Returns:
        List[List[]] - A list of lists representing the spreadsheet.
    """ 
    toRetLst=[]
    
    genesStrToAdd=""
    for item in genesToUse:
        genesStrToAdd=genesStrToAdd+item+"_"
    
    #First create the header for percentiles
    perc = []
    perc.append("Percentile")
    for x in range(1,101):
        perc.append(x)    
    toRetLst.append(perc)
    #Now add the 5', cds, 3' sections for each group
        #For each parclip
    for pc in lstOfPCs:
        #Get the row labeling strings to add
        labelEH1 = pc.filename+": "+"5' Exon Half"
        labelI = pc.filename+": "+"Intron "
        labelEH2 = pc.filename+": "+"3' Exon Half "
            #Get the EH1'
        DEH1 = pc.eFHDistribution
            #Add the label
        DEH1.insert(0,labelEH1)
        toRetLst.append(DEH1)
            #Get the Intron
        DI = pc.iDistribution
            #Add the label
        DI.insert(0,labelI)
        toRetLst.append(DI)
            #Get the second half of the exon
        DEH2 = pc.eSHDistribution
            #Add the label
        DEH2.insert(0,labelEH2)
        toRetLst.append(DEH2)
    
      #Export as csv--------
    csvfile = myDir
    
    if len(genesToUse)==0:
        csvfile=csvfile+"/IntronExonTable.csv"
    elif len(genesToUse)!=0:
        csvfile=csvfile+"/"+genesStrToAdd+"Intron_Exon_Table.csv"
    #Assuming res is a list of lists
    with open(csvfile, "w") as output:
        writer = csv.writer(output, lineterminator='\n')
        writer.writerows(toRetLst)          
    return(toRetLst)
    

def beforeColon(myString):
    """
    This function returns a subtring of all characters in a string before the first colon in that string
    Input:
        num(int) - String with a colon from which the substring will be pulled.
    Returns:
        List[int] - A list of unique random numbers between 0 and the inputed upper bound of length num.
    """ 
    finalIndex=0
    for ch in myString:
        finalIndex=finalIndex+1
        if ch==":":
            break
    #Get the new string
    toRet=myString[0:finalIndex-1]
    return(toRet)


def indicesToCheck(exportTable):
    """
    This function returns a list of indices to check in the sub plot function.
    Input:
        exportTable(List[List[int]]) - Export table, which is a list of lists, created via the createExportTable function.
    Returns:
        List[int] - A list of indices in the table that should be checked in creating the distribution plots.
    """ 
    theLength=len(exportTable)
    theLength=theLength-1
    numEntries=int(theLength/3)
    toRet=[]
    
    #Populate with 3's
    for i in range(0,numEntries):
        toAdd = 1+i*3
        toRet.append(toAdd)       
    return(toRet)

def areAnotRowsOverlapping(row1,row2):
    """
    Checks if two anotation rows are overlapping
    Input:
        anotRow1
        anotRow2
    Returns:
        bool
    """ 
    #Type match    
    for i in range(0,1):
        if i==0:
            ro1=row1
            ro2=row2
        elif i==1:
            ro1=row2
            ro2=row1
        
        #Run the comparisons
        if ro2.start >=ro1.start:
            if ro2.start <=ro1.stop:
                return(True)
        
        if ro2.stop<=ro1.stop:
            if ro2.stop>=ro1.start:
                return(True)
        
        if ro2.start<=ro1.start and ro2.stop>=ro1.stop:
            return(True)
    
    return(False)
    
    
def getHorizontalRange(maxVal):
    """
    This function takes in the maximum domain value and returns a list of values corresponding to the y values of the horizontal lines of the distribution plots (indicating breaks across start/stop codon)
    Input:
        maxVal(int) - Maximum x value to be used in the distribution plots. 
    Returns:
        List[int] - A list of the y values of the horizontal lines to seperate the 5'UTR, CDS, and 3'UTR regions.
    """     
    #Add two
    nMax = maxVal+maxVal/4
    #Divide this into 5
    inc=nMax/5
    toRet=[]
    for w in range(0,5):
        toRet.append(w*inc)
    return(toRet)
    

def createSubPlots(exportTable,outputDir,NIntronBins,NExonBins,genesToUse,lstOfPCs,inpDPI,imgFormat):
    """
    Graphs the binding distribution (across the intronic/exonic regions) of each parclip objects in a line. These individual plots are stacked together and written to a single file.
    Input:
        exportTable(List[List]) - Representation of all parclip's binding distributions as created by the getExportTable function.
        outputDir(str) - Directory in which to creat the sub plots.
        NIntronBins(Int) - The number of bins the intronic regions were divided into.
        NExonBins(Int) - The number of bins the exonic regions were divided into.
        genesToUse(List[str]) - List of genes to consider (identified by name). If empty, all genes considered.
        lstOfPCs(List[parclip]) - The list of parclip's whose binding properties are being visualized.
        inpDPI(Int) - The resolution of the graph outputed.
        imgFormat(Str) - The format to encode the image as. If "pdf" is used, then the file will be exported as a pdf. Otherwise, it will be exported as a png.
    """
    #Grab the domain row
    #First create the header for percentiles
    genesStrToAdd=""
    for item in genesToUse:
        genesStrToAdd=genesStrToAdd+item+"_"
        
    theDomain = []
    for x in range(1,NIntronBins+NExonBins-3):
        theDomain.append(x)    

    fig=plt.figure()
    indic=0
    for w in indicesToCheck(exportTable):
        #Determine the number of intron genes in the PC
        numInt=lstOfPCs[indic].getNumIntronGenes()
        numIntDisp="Genes: "+str(numInt)
        numExo1 = lstOfPCs[indic].getNumExonGenesH1()
        numExo1Disp="Genes: "+str(numExo1)
        numExo2 = lstOfPCs[indic].getNumExonGenesH2()
        numExo2Disp="Genes: "+str(numExo2)
        numExo1Clust=lstOfPCs[indic].getNumExonClustersH1()
        numExo1ClustDisp="Clusters: "+str(numExo1Clust)
        numExo2Clust=lstOfPCs[indic].getNumExonClustersH2()
        numExo2ClustDisp="Clusters: "+str(numExo2Clust)
        numIntClust=lstOfPCs[indic].getNumIntronClusters()
        numIntClustDisp="Clusters: " +str(numIntClust)
        
        
        #Grab the range, which is actually a collection of three lists
        #Get the first list
        lst1 = exportTable[w]
        #Get the title
        sectionTitle= beforeColon(lst1[0])
        #Be sure the name of the y axis does not contain the file path
        safeTitle=""
        for c in reversed(sectionTitle):
            if c=="/":
                break
            else:
                safeTitle=safeTitle+c
        sectionTitle=safeTitle[::-1]                        
        
        del lst1[0]
        #Get the second list
        lst2 = exportTable[w+1]
        del lst2[0]
        #Get the third list
        lst3 = exportTable[w+2]
        del lst3[0]
        
        #Get the range
        theRange=lst1+lst2+lst3
        #Get the max value
        mVal = max(theRange)*1.2
        
        val=.1+.4*(indic)
        #Add the axes
        nAx = fig.add_axes([.1,val,.8,.4], ylim=(0, mVal))
                        #Add HL's
        line1Val = NExonBins//2+0.5
        line2Val = line1Val + NIntronBins
            #Get the vertical values
        vertVals = getHorizontalRange(mVal)
        nAx.plot([line1Val,line1Val,line1Val,line1Val,line1Val],vertVals,linestyle='dashed',color='lightskyblue')
        nAx.plot([line2Val,line2Val,line2Val,line2Val,line2Val],vertVals,linestyle='dashed',color='lightskyblue')
        
        #Add the number of genes used in the and exonic count
        nAx.text(line1Val+.5*NIntronBins, mVal*.926, numIntDisp,fontsize=7,color='m')
        nAx.text(0, mVal*.94, numExo1Disp,fontsize=7,color='m')
        nAx.text(line2Val, mVal*.94, numExo2Disp,fontsize=7,color='m')
        nAx.text(0,mVal*.88,numExo1ClustDisp,fontsize=7,color='m')
        nAx.text(line2Val,mVal*.88,numExo2ClustDisp,fontsize=7,color='m')
        nAx.text(line1Val+.5*NIntronBins,mVal*.88,numIntClustDisp,fontsize=7,color='m') 
        
        #Add the plot
        nAx.plot(theDomain, theRange, color='g',linewidth=0.5)
        nAx.set_ylabel(sectionTitle, color='g',fontsize=13)
        nAx.set_title("Distribution of Intron / Exon Sites",y=1.1,fontsize=15)
        nAx.axes.get_xaxis().set_visible(False)
                
        indic=indic+1
        if (indic==len(indicesToCheck(exportTable))):
            nAx.axes.text(-7, mVal*1.05, "5' Exon Half", fontsize=10,color='b')
            nAx.axes.text(line1Val+.5*NIntronBins, mVal*1.04, "Intron", fontsize=10,color='b')
            nAx.axes.text(line2Val+.15*NExonBins//2, mVal*1.04, "3' Exon Half", fontsize=10,color='b')
    
    os.chdir(outputDir)
    
    if len(genesToUse)==0:
        if imgFormat=="pdf":
            toSave='DistributionPlot_Intron_Exon.pdf'
        if imgFormat!="pdf":
            toSave='DistributionPlot_Intron_Exon.png'
    elif len(genesToUse)!=0:
        if imgFormat=="pdf":
            toSave=genesStrToAdd+"DistributionPlot_Intron_Exon.pdf"
        if imgFormat!="pdf":
            toSave=genesStrToAdd+"DistributionPlot_Intron_Exon.png"
        
    fig.savefig(toSave, bbox_inches='tight',dpi=inpDPI)


def Average(lst):
    """
    Python program to get average of a list 
    Input:
        lst(List[Number]) - A list of numbers.
    Returns:
        float - The average of the numbers in input list.
    """    
    if len(lst)!=0:
        return sum(lst) / len(lst) 
    else:
        return(0)


def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Arguments:
        data       : A 2D numpy array of shape (N,M)
        row_labels : A list or array of length N with the labels
                     for the rows
        col_labels : A list or array of length M with the labels
                     for the columns
    Optional arguments:
        ax         : A matplotlib.axes.Axes instance to which the heatmap
                     is plotted. If not provided, use current axes or
                     create a new one.
        cbar_kw    : A dictionary with arguments to
                     :meth:`matplotlib.Figure.colorbar`.
        cbarlabel  : The label for the colorbar
    All other arguments are directly passed on to the imshow call.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)
    
     # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, fraction=0.026, pad=0.01, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom",fontsize="13")


    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    ax.set_aspect(10) # X scale matches Y scale
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels,fontsize=8) #Was 15
    ax.set_yticklabels(row_labels,fontsize=16)


    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-90, ha="center",
             rotation_mode="default")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3,axis='y')
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.set_title("Exon/Intron Distribution Heat Map",fontsize=25)
    
    return im, cbar


def scaleList(inpList):
    """
    Scales a list of numbers to all be relative to the largest number in the list (all fractions) 
    Input:
        List[Number] - A list of numbers.
    Returns:
        List[Number] - A scaled list of numbers relative to the largest.
    """   
    toRetLst = []
    #Get the largest value in the list
    largestValue = max(inpList)
    #Find out what you need to multiply this by in order to get 1
    coeff = 1/largestValue
    #Now Multiply each value in the list by this scaling value
    for i in range(0,len(inpList)):
        nVal = inpList[i]*coeff
        #Set this to be the value
        toRetLst.append(nVal)
    return(toRetLst)


def s_heatmap(pdDF,outDir,xAxisLab,genesToUse,inpDPI,imgFrmt):
    """
    This function generates and saves a clustered heatmap using seaborn's cluster functionality and dendrogram visualization techinique.
    Inputs:
        pdDF (Pandas DataFrame) - The dataframe of the metagene results to visualize. Row indices are the name of the metagene.
        outDir (Str) - Output directory in which the heatmap is saved.
        xAxisLab (List[Str]) - A list of the values that constitute the x axis labels of the heatmap.
        genesToUse(List[str]) - List of genes to consider (identified by name). If empty, all genes considered.
        inpDPI(Int) - The resolution of the graph outputed.
        imgFormat(Str) - The format to encode the image as. If "pdf" is used, then the file will be exported as a pdf. Otherwise, it will be exported as a png.        
    """    
    genesStrToAdd=""
    for item in genesToUse:
        genesStrToAdd=genesStrToAdd+item+"_"
        
    sns.set(font_scale=1.6)
    #Determine the height and width to use
    height=len(pdDF.index.values)*0.6
    
    #Intitalize the cluster heat map
    g = sns.clustermap(pdDF, col_cluster=False,cmap="coolwarm",xticklabels=xAxisLab,yticklabels=pdDF.index.values,figsize=(60, height),metric="euclidean",method="centroid")
    #Move the color bar
    #g.cax.set_position((.1,.1,.1,.1))
    
    #Remove the labels on the axes
    ax = g.ax_heatmap
    ax.set_xlabel("")
    ax.set_ylabel("")
    
    cbar = ax.collections[0].colorbar
    cbar.set_ticks([0,1])
    
    #Add seperating white lines between ther rows.
    for i in range(0,len(pdDF.index.values)):
        insAt=i+1
        ax.axhline(insAt, 0, insAt, linewidth=3, c='w')
    
    #Rotate the x-axis tick marks
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)

    #Set the title
    plt.title("Intron/Exon' Cluster Distribution", fontsize = 50, loc='center')


    #Save the plot to the output directory
    
    if len(genesToUse)==0:
        if imgFrmt=="pdf":
            toSave='DistributionHeatMap_Intron_Exon.pdf'
        elif imgFrmt=="png":
            toSave='DistributionHeatMap_Intron_Exon.png'
    elif len(genesToUse)!=0:
        if imgFrmt=="pdf":
            toSave=genesStrToAdd+"DistributionHeatMap_Intron_Exon.pdf"
        elif imgFrmt=="png":
            toSave='DistributionHeatMap_Intron_Exon.png'
    g.savefig(toSave, bbox_inches='tight',dpi=inpDPI)
        

def createHeatMapSingleton(lstOfParClips,outputDir,numIntronBins,numExonBins,genesToUse,inpDPI,imgFrmt):
    """
    Creates and saves the heat map when no clustering is possible (when there is only one input file).
    Input:
        lstOfPCs(List[parclip]) - The list of parclip's whose binding properties are being visualized.
        outputDir (Str) - Output directory in which the heatmap is saved.
        NIntronBins(Int) - The number of bins the intronic regions were divided into.
        NExonBins(Int) - The number of bins the exonic regions were divided into.
        genesToUse(List[str]) - List of genes to consider (identified by name). If empty, all genes considered.
        inpDPI(Int) - The resolution of the graph outputed.
        imgFormat(Str) - The format to encode the image as. If "pdf" is used, then the file will be exported as a pdf. Otherwise, it will be exported as a png.     
        """
    genesStrToAdd=""
    for item in genesToUse:
        genesStrToAdd=genesStrToAdd+item+"_"
        
    #Establish the x axis labels
    xAxisLabels=[]
    for x in range(0,numIntronBins+numExonBins):
        if x==numExonBins//2:
            xAxisLabels.append("-Intron Start")
        elif x==numExonBins//2+numIntronBins:
            xAxisLabels.append("-Intron Stop")
        elif x!=numExonBins//2 and x!=numExonBins//2+numIntronBins:
            xAxisLabels.append("")
            
    #Establish the y catagory labels
    yAxisLabels=[]
    for pc in lstOfParClips:
        toAdd=""
        for c in reversed(pc.filename):
            if c=="/":
                break
            toAdd=toAdd+c
        yAxisLabels.append(toAdd[::-1])
    
    #Create the matrix of values
    valueMatrix=[]
    for pc in lstOfParClips:
        #Create the row to add to the value matrix
        ro=pc.eFHDistribution+pc.iDistribution+pc.eSHDistribution
        
               #Normalize this row
        nRo = scaleList(ro)
            
        valueMatrix.append(nRo)
        
    #Convert this matrix to a numpy array
    valueArray = np.array(valueMatrix)

    
    #Creates the plot
    fig, ax = plt.subplots()
    im, cbar = heatmap(valueArray, yAxisLabels, xAxisLabels, ax=ax,cmap="coolwarm", cbarlabel="Coverage")
    fig.tight_layout()
    
    #Determine height
    height=len(lstOfParClips)*0.3
    
    os.chdir(outputDir)
    fig.set_size_inches(60, height)
    
    if len(genesToUse)==0:
        if imgFrmt=="pdf":
            toSave='DistributionHeatMap_Intron_Exon.pdf'
        elif imgFrmt=="png":
            toSave='DistributionHeatMap_Intron_Exon.png'
    elif len(genesToUse)!=0:
        if imgFrmt=="pdf":
            toSave=genesStrToAdd+"DistributionHeatMap_Intron_Exon.pdf"
        elif imgFrmt=="png":
            toSave='DistributionHeatMap_Intron_Exon.png'
        
    fig.savefig(toSave, bbox_inches='tight',dpi=inpDPI)


def createHeatMap(lstOfParClips,numIntronBins,numExonBins,outputDir,genesToUse,inpDPI,imgFrmt):
    """
    Creates a heat map of metagene distribution data for all parclips considered where the intensities are scaled to be relative to the greatest intensity value in each distribution.
    Input:
        lstOfParClips(List[parclip]) - A list of parclip objects to consider.
        numIntronBins(Int) - The number of bins the intronic regions were divided into.
        numExonBins(Int) - The number of bins the exonic regions were divided into.
        outputDir (Str) - Output directory in which the heatmap is saved.
        genesToUse(List[str]) - List of genes to consider (identified by name). If empty, all genes considered.
        inpDPI(Int) - The resolution of the graph outputed.
        imgFormat(Str) - The format to encode the image as. If "pdf" is used, then the file will be exported as a pdf. Otherwise, it will be exported as a png. 
    """ 
            #Establish the x axis labels
    xAxisLab=[]
    for x in range(0,numIntronBins+numExonBins):
        if x==numExonBins//2:
            xAxisLab.append("-Intron Start")
        elif x==numExonBins//2+numIntronBins:
            xAxisLab.append("-Intron Stop")
        elif x!=numExonBins//2 and x!=numExonBins//2+numIntronBins:
            xAxisLab.append("")
    
    genesStrToAdd=""
    for item in genesToUse:
        genesStrToAdd=genesStrToAdd+item+"_"
        

    #Establish the y catagory labels
    yAxisLabels=[]
    for pc in lstOfParClips:
        yAxisLabels.append(pc.filename)
    
    #Create the matrix of values
    valueMatrix=[]
    for pc in lstOfParClips:
        #Create the row to add to the value matrix
        ro=pc.eFHDistribution+pc.iDistribution+pc.eSHDistribution
        
        #Normalize this row
        nRo = scaleList(ro)
            
        valueMatrix.append(nRo)
        
    #Convert this matrix to a pandas dataframe
    valueDF = pd.DataFrame(valueMatrix,index=yAxisLabels)

    s_heatmap(valueDF,outputDir,xAxisLab,genesToUse,inpDPI,imgFrmt)

##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------         
def run(gtfFileName,csvDir,outputDir,NIntronBins,NExonBins,genesToUse,clusterThreshold,lowerBound,upperBound,boundFilter,weighValuesByRC,randomStatesLst,dpi,imgFormat,wildTypeGE,mainChromosomesOnly,bedAnotGTF,anoterScript):
    """
    Calculate and visualize the distribution of binding clusters across the intronic and exonic regions of any number of parclip experiments.
    Input:
        gtfFileName(Str) - The name of the gtf intermediate file (generated with a seperate python program) to use as a reference for annotation purposes.
        csvDir(Str) - The directory containing all of the cluster csv files (of the correct format) to search.
        outputDir (Str) - Output directory in which the heatmap is saved.      
        NIntronBins(Int) - The number of bins the intronic regions were divided into.
        NExonBins(Int) - The number of bins the exonic regions were divided into.       
        genesToUse(List[str]) - List of genes to consider (identified by name). If empty, all genes considered.
        clusterThreshold (Int) - The minimum number of clusters that must align to a gene in order for it to be considered in the metagen calculations.
        lowerBound (num) - The lower percentage of clusters to throw out based on the boundFilter property.
        uppBound (num) - The upper percentage of clusters to throw out based on the boundFilter property.
        boundFilter(Str) - The property used to sort the clusters in a data set before the bounds are applied to remove clusters based on this property.
            Options: "start" (start coordinate of cluster),"end" (end coordinate of cluster),"CS" (Conversion Specificity of cluster),"URC" (unique read count of the cluster),"RC" (read count of the cluster)
        weighValuesByRC (Bool) - Apply a weight to each cluster percentage distribution based on the number of reads aligned to the cluster?
        randomStatesList (List[int]) - The list of numbers to use as random seeds in random number generation. The more random states provided, the more times the binning algorithm will execute before the averaging step.
        dpi(int) - The dpi to save images to.
        imgFormat (str) - "pdf" or "png". The format to save images into.
        wildTypeGE (str) - Name of a wild type gene expression file which contain "Average" and "Gene Symbol" columns.
        mainChromosomesOnly (bool) - Whether the program should only look at the "main" chromsomes (autosomal plus x/y) for analysis.
        bedAnotGTF (Str) - GTF File to be used for the bed annotation if need be.
        anoterScript (str) - The bed anotation script to use on bed files.
    """ 
    if (len(wildTypeGE)>0):
        print("Loading wild type gene expression file...")
        wtGE_File=WTGeneExpression(wildTypeGE)
        print("...loaded.")
        
    if len(wildTypeGE)==0:
        wtGE_File=False

    print("Loading gtf longest transcripts intermediate file...")
    myGTFObj = readInGTFFile(gtfFileName,mainChromosomesOnly)   
    print("...loaded.")

    print("Loading the list of parclip cluster .csv files from: "+csvDir+"...")
    lstOfParclips = createLstOfParclips(csvDir,genesToUse,lowerBound,upperBound,boundFilter,myGTFObj,bedAnotGTF,anoterScript,outputDir)
    print("...loaded.")

    print("Populating hits and overlap values for all genes in all parclip cluster files...")
    for exPC in lstOfParclips:
        exPC.populateAllHits(myGTFObj,NIntronBins,NExonBins,randomStatesLst)
        exPC.getParclipPercentageDistributions(NIntronBins,NExonBins,clusterThreshold,weighValuesByRC,wtGE_File)
        exPC.populateExonHalves(NExonBins)
    print("...populated.")

    print("Creating the export table...")
    exportTable = createExportTable(lstOfParclips,outputDir,genesToUse) 
    print("...export table created")
    ###EIlengths = outputAvgLengths(lstOfParclips,outputDir)

    print("Creating the binding distribution curves for all parclips in the directory...")
    createSubPlots(exportTable,outputDir,NIntronBins,NExonBins,genesToUse,lstOfParclips,dpi,imgFormat)
    print("...binding distribution curves generated and merged.")
    
    print("Generating the heat map...")
    if len(lstOfParclips)>1:
        createHeatMap(lstOfParclips,NIntronBins,NExonBins,outputDir,genesToUse,dpi,imgFormat)
    else:
        createHeatMapSingleton(lstOfParclips,outputDir,NIntronBins,NExonBins,genesToUse,dpi,imgFormat)
    print("... heat map generated.")

##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Run
##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ 
#run(gtfFileName,csvDir,outputDir,1250,100,gtu,clustThreshold,lBound,uBound,boundF,wValuesByRC,randStates,dpi,imgFrmt,wtGE,mainChromosomes,theBED_GTF,theBEDAnotScript)