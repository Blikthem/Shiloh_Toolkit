#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 15:40:48 2018

@author: claypooldj
"""
##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Load packages 
##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
import csv
import os
import random
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
geneListMaster=[]

"""
##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Enter parameters
##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ 
gtfFileName="/home/claypooldj/genomes/gencode_hg38.longestTranscript.csv" #The intermediate gtf longest splice varient csv file (generated from a seperate script) to use for annotation. {Type - String}
#gtfFileName="/Volumes/Untitled/Genomes/gencode_hg38.longestTranscriptBEDTEST.csv"
csvDir="/home/claypooldj/myPythonScripts/testBedHM/hInputs/" #The directory containing the cluster csv files to analyze. All csv files in this directory will be tested (so make sure there aren't other csv files in the directory or it will crash!) {Type - String}. Alternatively, a single csv file can be specified here and the analysis will be run only on that file.
outputDir="/home/claypooldj/myPythonScripts/testBedHM/output/"    #The directory where you want to generate all output files. {Type - String}
boundF="URC"    #The continuous field in the clusters file to use to apply filtering bounds (both upper and lower) {Type - String}
lBound=0    #The percentage of clusters to filter based on the lowest values of the bounding continuous variable (boundF) {Type - Num}
uBound=0    #The percentage of clusters to filter based on the largest values of the bounding continuous variable (boundF) {Type - Num}
gtu=[]    #The list of genes which you wish to consider. Just put an empty list '[]' if you want to consider all genes. {Type - List[String] or empty List}
wValuesByRC=False    #Whether or not the program should weight the impact of each indivudal gene on the overall metagene average linearly with the number of clusters that align to that gene {Type - Bool}
dpi=100    #The resolution of the output graphs in dpi {Type - Int}
imgFormat="pdf"    #Which format to save the output graphs to. Two options, either "pdf" which saves as pdf otherwise it saves to png format. {Type - String}
randStates=[7211995,541995,3131994,111,222,333,444,555,888,999]    #The random states to use in the binning algorithm, which randomly deals with rounding error and then averages the results. Can be any length. {Type- List[Int]}
randStates=[541995]
wtGE="/Volumes/Untitled/Output/explorador/RNA_Seq_Background/moreWTExpressions/geneExpressionMatrix_TPM_withGeneNames.csv"    #A csv file containing the wild type gene expression levels from RNA-Seq(s) of the cell line. Takes average value. If this string is empty, we will not weigh values by WT gene expression.
wtGE=""
mainChromosomes=True #Whether we should only consider the main chromsomes (the autosomal and X/Y) in the analysis and not the strange chromosome constructs used for.
theBED_GTF="/home/claypooldj/genomes/gencode.v30.annotation.gtf"
theBEDAnotScript="/home/claypooldj/clip_metagene/bedannotator.sh"
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
        Initalizing method.
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
        Checks to see if a given index is within the bounds of this annotation row object
        Inputs:
            index(int) - The index to check
        Returns:
            boolean - Whether the input index is within the bounds of the object
        """
        if (self.stop>=index and index>=self.start):
            return(True)
        return(False)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

class anotGene:
    """
    Anotation Gene Object: This object represents an anotation gene (a collection of row objects from the annotation gtf file that represent a single gene at the longest transcript representation).

    Properties:
        gene(String) - The name of the gene these rows map to.
        chrom(String) - The chromosome these rows map to.
        start(int) - The smallest starting index for all annotation rows in this gene.
        stop(int) - The largest ending index for all annotation rows in this gene.
        lstOfAnotRows(list[anotRows]) - A list of all of the anotRow objects in this anotationGene object (all rows aligning to this gene).
        fIndexStarts(list[int]) - The starting index values of all of the 5'UTR rows in this anotationGene object.
        tIndexStarts(list[int]) - The starting index values of all of the 3'UTR rows in this anotationGene object.
        CDSIndexStarts(list[int]) - The starting index values of all of the CDS rows in this anotationGene object.
        fIndexStops(list[int]) - The stop index values of all of the 5'UTR rows in this anotationGene object.
        tIndexStops(list[int]) - The stop index values of all of the 3'UTR rows in this anotationGene object.
        CDSIndexStops(list[int]) - The stop index values of all of the CDS UTR rows in this anotationGene object.
        getDictOverlapPercentage(dictionary) - Determines the nucleotide overlap of a dictionary to this gene as a percentage of the total length.
    
    Methods:
        init - Initializer
        str - Conversion to string form for printing.
        getStart - Returns the smallest(first) starting index of all the rows in this anotationGene.
        addAR - Adds a row to this anotationGene object.
        getStartandStopIndices - Populates the start and stop lists for this object.
        getType - Returns the classification (CDS/TUTR/etc) associated with a given nucleotide index in this anotationGene.
    """
    def __init__(self,gene,chrom):
        self.gene=gene
        self.chrom=chrom
        self.start=0
        self.stop=1
        self.lstOfAnotRows=[]
        self.fIndexStarts=[]
        self.fIndexStops=[]
        self.tIndexStarts=[]
        self.tIndexStops=[]
        self.CDSIndexStarts=[]
        self.CDSIndexStops=[]
        
    def __str__(self):
        toPrint="{Annotaion Gene Object} Gene: "+self.gene+"\tChromosome: "+self.chrom+"\tNumber of Annotation Row Objects: "+str(len(self.lstOfAnotRows))
        return(toPrint)
       
    def getStart(self):
        """
        Finds the smallest starting index of the rows in the anotGene.
        Returns:
            int - The first nucleotide position in the anotGene
        """ 
        possibleVals=[]
        for row in self.lstOfAnotRows:
            if (row.ty=="CDS" or row.ty=="FUTR" or row.ty=="TUTR"):
                possibleVals.append(row.start)
        if (len(possibleVals)==0):
            return(False)
        theMin = min(possibleVals)
        self.start=theMin
        return(theMin)
        
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
        for ro in self.lstOfAnotRows:
            if int(ro.stop)!=0: 
                ##!!!!!! Changed to not make the determining clause
                if int(ro.stop)>int(self.stop):
                    self.stop=int(ro.stop)
                if ro.ty=="FUTR":
                    #Get the starting value
                    fullVal = int(ro.start)
                    #Get the value in context of the overall
                    nVal = fullVal-gStart
                    
                    #Get the stopping value
                    fullVal2=int(ro.stop)
                    sVal=fullVal2-gStart
                    
                    #If neither are in the list already: 
                    if nVal not in self.fIndexStarts and sVal not in self.fIndexStops:
                        self.fIndexStarts.append(nVal)
                        self.fIndexStops.append(sVal)
                        
                    #If only one of them would be added...
                    if nVal not in self.fIndexStarts and sVal in self.fIndexStops:
                        #Find the index to swap 
                        for i in range(0,len(self.fIndexStops)):
                            curV=self.fIndexStops[i]
                            if curV==sVal:
                                #Now we need to resolve which start value to grab (the smaller)
                                if nVal<self.fIndexStarts[i]:
                                    self.fIndexStarts[i]=nVal
                    
                    
                    if nVal in self.fIndexStarts and sVal not in self.fIndexStops:
                        #Find the index to swap 
                        for i in range(0,len(self.fIndexStarts)):
                            curV=self.fIndexStarts[i]
                            if curV==nVal:
                                #Now we need to resolve which start value to grab (the larger)
                                if sVal>self.fIndexStops[i]:
                                    self.fIndexStops[i]=nVal
                    
                if ro.ty=="TUTR":
                    #Get the starting value
                    fullVal = int(ro.start)
                    #Get the value in context of the overall
                    nVal = fullVal-gStart
                    
                    #Get the stopping value
                    fullVal2=int(ro.stop)
                    sVal=fullVal2-gStart
                    
                    #If neither are in the list already: 
                    if nVal not in self.tIndexStarts and sVal not in self.tIndexStops:
                        self.tIndexStarts.append(nVal)
                        self.tIndexStops.append(sVal)
                        
                    #If only one of them would be added...
                    if nVal not in self.tIndexStarts and sVal in self.tIndexStops:
                        #Find the index to swap 
                        for i in range(0,len(self.tIndexStops)):
                            curV=self.tIndexStops[i]
                            if curV==sVal:
                                #Now we need to resolve which start value to grab (the smaller)
                                if nVal<self.tIndexStarts[i]:
                                    self.tIndexStarts[i]=nVal
                    
                    
                    if nVal in self.tIndexStarts and sVal not in self.tIndexStops:
                        #Find the index to swap 
                        for i in range(0,len(self.tIndexStarts)):
                            curV=self.tIndexStarts[i]
                            if curV==nVal:
                                #Now we need to resolve which start value to grab (the larger)
                                if sVal>self.tIndexStops[i]:
                                    self.tIndexStops[i]=nVal
                
                if ro.ty=="CDS":
                    #Same thing but I am not adding the values
                    #Get the starting value
                    fullVal = int(ro.start)
                    #Get the value in context of the overall
                    nVal = fullVal-gStart
                    
                    #Get the stopping value
                    fullVal2=int(ro.stop)
                    sVal=fullVal2-gStart

                    #If neither are in the list already: 
                    if nVal not in self.CDSIndexStarts and sVal not in self.CDSIndexStops:
                        self.CDSIndexStarts.append(nVal)
                        self.CDSIndexStops.append(sVal)
                        
                    #If only one of them would be added...
                    if nVal not in self.CDSIndexStarts and sVal in self.CDSIndexStops:
                        #Find the index to swap 
                        for i in range(0,len(self.CDSIndexStops)):
                            curV=self.CDSIndexStops[i]
                            if curV==sVal:
                                #Now we need to resolve which start value to grab (the smaller)
                                if nVal<self.CDSIndexStarts[i]:
                                    self.CDSIndexStarts[i]=nVal
                    
                    
                    if nVal in self.CDSIndexStarts and sVal not in self.CDSIndexStops:
                        #Find the index to swap 
                        for i in range(0,len(self.CDSIndexStarts)):
                            curV=self.CDSIndexStarts[i]
                            if curV==nVal:
                                #Now we need to resolve which start value to grab (the larger)
                                if sVal>self.CDSIndexStops[i]:
                                    self.CDSIndexStops[i]=nVal
            

    def getType(self,index):
        """
        This function will return the character indicating the classification of the gene at a particular index
        Inputs:
            index(int) - The index to check
        Returns:
            string - What type of sequence this index alignes to. NC means that it is not annotated.
        """
        possibleOptions=[]
        #Go through each row of the anotGene
        for r in range(0,len(self.lstOfAnotRows)):
            ro = self.lstOfAnotRows[r]
            #If it is contained
            if ro.coversIndex(index)==True:
                possibleOptions.append(ro.ty)
                
                
        #Now pic which one to return
        if "FUTR" in possibleOptions:
            return("FUTR")
        if "TUTR" in possibleOptions:
            return("TUTR")
        if "CDS" in possibleOptions:
            return("CDS")
        if "intron" in possibleOptions:
            return("intron")
        if "exon" in possibleOptions:
            return("exon")
        #print(possibleOptions)
        return("NC")

    def getDictOverlapPercentage(self,dictionary):
        """
        Determines the fraciton of this gene that is overlaped by a dictionary representation of a cluster.
        Input:
            dictionary (dict) - Dictionary representation of a cluster.
        Return:
            float - The fraction of overlap.
        """
        return(getNumericalOverlap(self.start,self.stop,dictionary.get("Start"),dictionary.get("End")))
        
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~         
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~         
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~       
class cluster:
    """
    Cluster Object: This object represents a cluster identified by PARALYZER - a row from a clusters.csv output file.

    Properties:
        Pulled from csv rows as strings - chrom (chromosome), strand (strand oriented to), start (starting index in cluster), end (ending index of cluster), gene (gene that paralyzer aligned this cluster to if any), CS (conversion specificity), T2C (T to C fraction), URC (Unique read count for this cluster), RC (total read count)
    
    Methods:
        init - Initializer
        str - Conversion to string form for printing
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
        
        
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
class clusterGene:
    """
    ClusterGene Object: A group of clusters of the same gene

    Properties:
        gene(string) - The name of the gene around which all contained clusters are aligned
        chrom(string) - The name of the chromosome this gene falls on
        start(int) - The start location of this clustergene (ie the first start value among contained clusters)
        stop(int) - The stop location of this clustergene (ie the last stop value among contained clusters)
        lstOfClusters (list[cluster]) - A list of all of the cluster objects aligned with this gene
        hitListCDS(list[int]) - A list of integers as long as the gene's CDS. Each number represents a nucleotide. A 1 means that there is a cluster overlapping the gene at this location that aligns to the CDS of the gene (based on the GTF file used). A 0 means there is no cluster that aligns to the CDS of the gene at that index.
        hitList3(list[int]) - A list of integers as long as the gene's 3'UTR. Each number represents a nucleotide. A 1 means that there is a cluster overlapping the gene at this location that aligns to the 3'UTR of the gene (based on the GTF file used). A 0 means there is no cluster that aligns to the 3'UTR of the gene at that index.
        hitList5(list[int]) - A list of integers as long as the gene's 5'UTR. Each number represents a nucleotide. A 1 means that there is a cluster overlapping the gene at this location that aligns to the 5'UTR of the gene (based on the GTF file used). A 0 means there is no cluster that aligns to the 5'UTR of the gene at that index.
        percentDCDS(list[float])) - The percent distribution of cluster hits aligning to the CDS across the gene's CDS region. Each number in the list represents the number of hits found in a percent bin (1%) of the total gene length. 
        percentD3(list[float])) - The percent distribution of cluster hits aligning to the 3'UTR across the gene's 3'UTR region. Each number in the list represents the number of hits found in a percent bin (1%) of the total gene length. 
        percentD5(list[float])) - The percent distribution of cluster hits aligning to the 5'UTR across the gene's 5'UTR region. Each number in the list represents the number of hits found in a percent bin (1%) of the total gene length. 
    
    Methods:
        init - Initializer
        str - Conversion to string form for printing
        addAC - Add a cluster object to this clusterGene
        getCoverageLength2 - Returns the total number of cluster nucleotides that overlap with this gene.
        containsIndex - Checks if a cluster in this clusterGene contains a given index.
        getMatchingMeta - Finds the metaGene in an anotFile object which has the same name as this clusterGene 
        getNucleotidePercent - Finds what percentage of the total gene cluster overlap a single nucleotide represents.
        populateHitLists - Populates the hit lists for this clusterGene object (a list which states whether each nucleotide in a region is overlapped with a cluster or not).
        populatePercentiles - Populates the percentile lists for this clusterGene object (the distribution of total cluster overlap, expressed as percentages, broken down by region).
        populateNumClust - Populates the number of clusters that overlap the 5', CDS, and 3' regions of the cluster gene.
        getTotalNumReads - Return the total number of reads in clusters aligned to this gene.
    """
    def __init__(self,gene,chrom):
        self.gene=gene
        self.chrom=chrom
        self.start=0
        self.stop=1
        self.lstOfClusters=[]
        self.hitListCDS=[]
        self.hitList3=[]
        self.hitList5=[]
        self.percentDCDS=[]
        self.percentD3=[]
        self.percentD5=[]
        self.nucPercent=1
        self.FClusterCount=0
        self.CDSClusterCount=0
        self.TClusterCount=0
        
    def __str__(self):
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
        toRet=0
        for indivClust in self.lstOfClusters:
            toRet=toRet+indivClust.RC
        return(toRet)
             
    def getCoverageLength2(self,anotFile):
        """
        Determines the total number of nucleotides which are covered/overlapped by clusters in the CDS, TUTR, AND FUTR regions
        Inputs:
            anotFile(cluster) - Cluster object to add to the ClusterGene
        Returns:
            int - The total number of overlapping nucleotides
        """  
        matchedMG = self.getMatchingMeta(anotFile)
        fcdc=0
        if (matchedMG!=False):                
            #Loop over each cluster ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            for clust in self.lstOfClusters:
                #Get the start of the cluster
                theStart = int(clust.start)
                #Get the end of the cluster
                theStop = int(clust.end)+1   
                    #Now we need to loop over every nucleotide in this range to populate hits
                for i in range(theStart,theStop):
                       #Now that we know there is a cluster, we need to check to see if it is identified as a specific coding region
                       if (matchedMG.getType(i)=="CDS" or matchedMG.getType(i)=="TUTR" or matchedMG.getType(i)=="FUTR"):
                           fcdc=fcdc+1 
        return(fcdc)
            
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
        
            
    def getNucleotidePercent(self,anotFile):
        """
        Gets the percent of the total gene cluster overlap that each nucleotide in that gene represents.
        Inputs:
            anotFile(anotFile) - An anotFile object containing the genes to check
        Returns:
            int - The percent of the total overlap each nucleotide in the gene represents.
            or False if the overlap length is found to be 0
        """  
        if (self.getCoverageLength2(anotFile)==0):
            return(False)
        nucPerc = 100.0/self.getCoverageLength2(anotFile)
        self.nucPercent=nucPerc
        return(nucPerc)
          
    def populateHitLists(self,anotFile):
        """
        Populate the CDS/3'UTR/5'UTR hit lists of the clusterGene by comparing clusters to a reference gtffile object
        Inputs:
            anotFile(anotFile) - An anotFile object containing the genes to check
        """
        #Get the matching metagene
        matchedMG = self.getMatchingMeta(anotFile)
        if matchedMG==False or self.gene=="":
            return()
        matchedMG.getStartandStopIndices()
        if matchedMG==False and self.gene!="" :
            print("Cant find: ",self.gene)
        
        theCatagories=["FUTR","CDS","TUTR"]
        for gp in theCatagories:
            if (matchedMG!=False):
                #Establish the start and stop of the specific region
                if gp=="FUTR":
                    endVal = len(matchedMG.fIndexStarts)
                    if len(matchedMG.fIndexStarts)!=len(matchedMG.fIndexStops):
                        print(matchedMG.gene," overlap error")
                        return()
                        
                if gp=="CDS":
                    endVal = len(matchedMG.CDSIndexStarts)
                    
                    if len(matchedMG.CDSIndexStarts)!=len(matchedMG.CDSIndexStops):
                        print(matchedMG.gene," overlap error")
                        return()
                    
                if gp=="TUTR":
                    endVal = len(matchedMG.tIndexStarts)
                    if len(matchedMG.tIndexStarts)!=len(matchedMG.tIndexStops):
                        print(matchedMG.gene," overlap error")
                        return()

                for i in range(0,endVal):
                
                    self.start=matchedMG.getStart()
                    if (gp=="FUTR"): 
                        theStart=matchedMG.fIndexStarts[i] + self.start
                        theStop=matchedMG.fIndexStops[i] + self.start+1
                
                    if (gp=="CDS"):
                        theStart=matchedMG.CDSIndexStarts[i] + self.start
                        theStop=matchedMG.CDSIndexStops[i] + self.start+1
                
                    if (gp=="TUTR"):
                        theStart=matchedMG.tIndexStarts[i] + self.start
                        theStop=matchedMG.tIndexStops[i] + self.start+1

                    for i in range(theStart,theStop):
                        #We need to check if there is a cluster there
                        putAZero=True
                        if self.containsIndex(i)==True:
                            #Now that we know there is a cluster, we need to check to see if it is identified as a specific coding region
                            if (matchedMG.getType(i)==gp):
                                if(gp=="FUTR"):
                                    self.hitList5.append(1)
                                    putAZero=False
                                if(gp=="CDS"):
                                    self.hitListCDS.append(1)
                                    putAZero=False
                                if(gp=="TUTR"):
                                    self.hitList3.append(1)  
                                    putAZero=False                                                     
                    
                        if(putAZero==True):
                            if(gp=="FUTR"):
                                self.hitList5.append(0)
                            if(gp=="CDS"):
                                self.hitListCDS.append(0)
                            if(gp=="TUTR"):
                                self.hitList3.append(0)                  

        
    def populatePercentiles(self,anotFile,randomStatesList):
        """
        Populates the CDS/3'UTR/5'UTR percentile lists of the clusterGene by looking at the respective hit lists 
        Inputs:
            anotFile(anotFile) - An anotFile object containing the genes to check
            randomStatesList(List[num]) - A list of random seeds. The length represents the number of binning processes to undergo.
        """  
        nucPerc = self.getNucleotidePercent(anotFile)
        
        check=0
        if (nucPerc==False):
            return(False)
        matchedGene=self.getMatchingMeta(anotFile)
        if (matchedGene==False):
            return(False)
        matchedGene.getStartandStopIndices()

        #Perform the binning step within a loop
        #Store the individual (each loop) hit lists in another list
        indivPDLsts5=[]
        indivPDLstsCDS=[]
        indivPDLsts3=[]
        for indivSeed in randomStatesList:
            percentD5=[]
            percentDCDS=[]
            percentD3=[]
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            theOpts=["FUTR","TUTR","CDS"]
            for grp in theOpts:
                
                #Establish the lengths to parse of the hit list
                if (grp=="FUTR"):
                    totLength = len(self.hitList5)
                    check=check+sum(self.hitList5)
            
                if (grp=="CDS"):
                    totLength = len(self.hitListCDS)
                    check=check+sum(self.hitListCDS)
                
                if (grp=="TUTR"):
                    totLength = len(self.hitList3)
                    check=check+sum(self.hitList3)
        
                if totLength==0:
                   # print("Failure at "+grp)
                    return("Failure")
                
            
                if totLength<100:
                    #Binsize
                    lBin=100//totLength
                    #Figure out how much is left over
                    leftOver = 100-lBin*totLength
                    
                    randIndices = getRandLstBounded(leftOver,totLength-1,indivSeed)
            
                    #Now we need to parse over the indices of the cluster            
                    addOne=False
                    indexInHitsOld=0
                    indexInHits=0
                    while indexInHits<totLength:
                        #If this index is in the list [its +1 from the bin size]
                        if indexInHitsOld in randIndices:
                            addOne=True
                    
                        #Get the values
                        if (grp=="FUTR"):
                            val = self.hitList5[indexInHits]
                            #Get the values
                
                        if (grp=="CDS"):
                            val = self.hitListCDS[indexInHits]
                        
                        if (grp=="TUTR"):
                            val = self.hitList3[indexInHits]
                
                        #Add the values
                        for i in range(0,lBin):    
                            if (grp=="TUTR"):
                                percentD3.append(val*nucPerc)
                            if (grp=="CDS"):
                                percentDCDS.append(val*nucPerc)
                            if (grp=="FUTR"):
                                percentD5.append(val*nucPerc)                        
                    
                        if (addOne==True):
                            if (grp=="TUTR"):
                                percentD3.append(val*nucPerc)
                            if (grp=="FUTR"):
                                percentD5.append(val*nucPerc)
                            if (grp=="CDS"):
                                percentDCDS.append(val*nucPerc)
            
                        indexInHits=indexInHits+1
                        indexInHitsOld=indexInHitsOld+1
                        addOne=False
                    
                    #Now we need to adjust for the added amount
                        #What it should equal
                        shouldEqual=1
                        if (grp=="FUTR"):
                            shouldEqual=sum(self.hitList5)*nucPerc
                
                        if (grp=="CDS"):
                            shouldEqual=sum(self.hitListCDS)*nucPerc
                        
                        if (grp=="TUTR"):
                            shouldEqual=sum(self.hitList3)*nucPerc
                        
                        #What it currently equals
                        curEquals=1
                        if (grp=="FUTR"):
                            curEquals=sum(percentD5)
                            #Get the values
                
                        if (grp=="CDS"):
                            curEquals=sum(percentDCDS)
                        
                        if (grp=="TUTR"):
                            curEquals=sum(percentD3)
                        
                        #Find the adjusting factor
                        if curEquals!=0:
                            adjFact = shouldEqual/curEquals
                            if (grp=="FUTR"):
                                for i in range(0,len(percentD5)):
                                    newEle = percentD5[i]*adjFact
                                    percentD5[i]=newEle
                                
                            if (grp=="CDS"):
                                for i in range(0,len(percentDCDS)):
                                    newEle = percentDCDS[i]*adjFact
                                    percentDCDS[i]=newEle
                        
                            if (grp=="TUTR"):
                                for i in range(0,len(percentD3)):
                                    newEle = percentD3[i]*adjFact
                                    percentD3[i]=newEle    
        
                wvals =[]
                if totLength>=100:
                    #Determine the lengths
                    rawLength=totLength/100.0
                    #Length of each bin before random addition
                    bLength= totLength//100
                    #Remainder
                    rLen = round((rawLength-bLength)*100)
                    #Create the random numbers needed
                    randIndices = getRandLst(rLen,indivSeed)
            
                    #Now we need to parse over the indices of the cluster            
                    nBins=0
                    indexInHits=0
                    while nBins!=100:
                        #If this index is in the list [its +1 from the bin size]
                        if nBins in randIndices:
                            upperCount=1+bLength

                        else:
                            upperCount=bLength
            
                        #Get the total val which represents the number of hits in this region (how many nucleotides of each type fall in the region)
                        totVal=0

                        for w in range(indexInHits,indexInHits+upperCount):
                            if (grp=="TUTR"):
                                totVal = totVal + self.hitList3[w]
                                wvals.append(w)
                            if (grp=="CDS"):
                                totVal = totVal + self.hitListCDS[w]  
                            if (grp=="FUTR"):
                                totVal = totVal + self.hitList5[w]   
                        

                        #Get the percent value (multiply by the percent represented by each in dividual nucleotide)
                        if (grp=="TUTR"):
                            percentD3.append(totVal*nucPerc)
                        if (grp=="CDS"):
                            percentDCDS.append(totVal*nucPerc)
                        if (grp=="FUTR"):
                            percentD5.append(totVal*nucPerc)
                                                            
                        nBins=nBins+1
                        indexInHits=indexInHits+upperCount
                     
            #If we are dealing with a reverse transcript we will need to inverse the percentile terms
            tester = self.lstOfClusters[0]
            orient2=tester.strand
            if orient2=="-":
                #Flip
                percentD3.reverse()
                percentD5.reverse()
                percentDCDS.reverse()
            
            #Append the individual lists into the data frame
            indivPDLsts5.append(percentD5)
            
            indivPDLstsCDS.append(percentDCDS)
            indivPDLsts3.append(percentD3)
        
        #Average the percentage lists
        DF5 = pd.DataFrame(indivPDLsts5)
        self.percentD5=list(DF5.mean(axis = 0))
        DFCDS = pd.DataFrame(indivPDLstsCDS)
        self.percentDCDS=list(DFCDS.mean(axis = 0))
        DF3 = pd.DataFrame(indivPDLsts3)
        self.percentD3=list(DF3.mean(axis = 0))

         #Check to make sure the coverage is summing to 100
        tester = self.lstOfClusters[0]
        orient2=tester.strand
        testRow=matchedGene.lstOfAnotRows[0]
        checker=[]
        for ele in matchedGene.lstOfAnotRows:
           checker.append(ele.ty)
        theTotal = sum(self.percentDCDS)+sum(self.percentD5)+sum(self.percentD3)
        
        if (round(theTotal)!=100):
            if round(theTotal)!=0 and len(self.percentDCDS)==100 and testRow.ty!="exon":  
                print("-------------------------------------------------------------------------")
                print(theTotal," Not reaching 100, ",self.gene," Matching gene: ",matchedGene.gene)
        
                self.hitList3=[]
                self.hitList5=[]
                self.hitListCDS=[]

    def populateNumClust(self,anotFile): 
        """
        Populate the number of clusters aligning to the 5'UTR, CDSn and 3'UTR regions of the gene.
        Inputs:
            anotFile(anotFile) - An anotFile object containing the genes to check
        """ 
        bool_FUTR=False
        bool_CDS=False
        bool_TUTR=False
        FUTR_count=0
        CDS_count=0 
        TUTR_count=0
        matchedMG = self.getMatchingMeta(anotFile)
        if matchedMG==False or self.gene=="":
            return()
        #Loop over each cluster in this gene
        for clust in self.lstOfClusters:
            #Go through the indices of the cluster
            for i in range(int(clust.start),int(clust.end)):
                #Check if it is 3'UTR
                if matchedMG.getType(i)=="FUTR":
                    bool_FUTR=True
                    
                if matchedMG.getType(i)=="CDS":
                    bool_CDS=True
                
                if matchedMG.getType(i)=="TUTR":
                    bool_TUTR=True
                
            if bool_FUTR==True:
                bool_FUTR=False
                FUTR_count=FUTR_count+1
            if bool_CDS==True:
                bool_CDS=False
                CDS_count=CDS_count+1
            if bool_TUTR==True:
                bool_TUTR=False
                TUTR_count=TUTR_count+1
            
        self.FClusterCount=FUTR_count
        self.CDSClusterCount=CDS_count
        self.TClusterCount=TUTR_count
        
        

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
class parclip:
    """
    Parclip Object: This object represents the output from an indivual parclip clusters.csv file.

    Properties:
        filename(String) - The name of the csv file.
        lstOfChromosomesClusters(List[Chromosomes]) - A list of the clusters in a parclip, sorted by chromosome.
        lstOfChromosomesClustersGenes - A list of the cluster genes in a parclip, sorted by chromosome.
        CDSDistribution(List[int]) - Average distribution of cluster coverage in the CDS of all of the genes in this parclip (expressed as percentages of total cluster overlap).
        fDistribution(List[int]) - Average distribution of cluster coverage in the 5'UTR of all of the genes in this parclip (expressed as percentages of total cluster overlap).
        tDistribution(List[int]) - Average distribution of cluster coverage in the 3'UTR of all of the genes in this parclip (expressed as percentages of total cluster overlap).
        totalList(List[int]) - Average distribution of cluster coverage in the 5'UTR, CDS, 3'UTR regions (the distribution lists appended in that order) for this parclip (expressed as a percentage of total cluster overlap).
        allClusters(List[Cluster]) - A list of all the clusters (in all of the chromosomes) associated with this parclip.
        
    
    Methods:
        init - Initializer
        str - Conversion to string form for printing.
        addClust - Add a cluster to the parclip object (used when reading in the object).
        getNumClusters - Determines how many clusters are present in a parclip object.
        getNumGenes - Determines how many  genes are present in a parclip object.
        populateCh - Populate the chromosome lists within parclip object.
        populateAllHits - Populate all of the hit and percentage lists for constituient clustergenes of the parclip object.
        getParclipPercentageDistributions - Populate the average distributions for the parclip object as a whole.
        getTotalList - Returns a merged list of the 5' UTR, CDS, and 3' UTR percentage distributions for parclip object.
        getNum5UTRGenes - Returns the number of genes in this PARCLIP with clusters overlapping their 5'UTR region.
        getNumCDSGenes - Returns the number of genes in this PARCLIP with clusters overlapping their CDS region.
        getNum3UTRGenes - Returns the number of genes in this PARCLIP with clusters overlapping their 3'UTR region.
        getNum5UTRClusters - Returns the number of clusters in the 5' region in this PARCLIP.
        getNumCDSClusters - Returns the number of clusters in the CDS region in this PARCLIP.
        getNum3UTRClusters - Returns the number of clusters in the 3' region in this PARCLIP.
        removeCluster - Given a cluster, this method goes through the lstOfChromsomesClusters property and removes the inputed cluster.
        applyBounds - Impose bounds on the clusters included - throwing out extremes (based on inputed property) of a specified % at the minimum and maximum end of the spectrum
    """
    def __init__(self,filename,lstOfUniqueChromosomes):
        self.filename=filename
        
        theChromsomeList=[]
        for chrN in lstOfUniqueChromosomes: 
            theChromsomeList.append(chromesome(chrN))

        theChromsomeList2=[]
        for chrN in lstOfUniqueChromosomes: 
            theChromsomeList2.append(chromesome(chrN))
            
        #self.lstOfChromosomesClusters=[chromesome("chr1"),chromesome("chr2"),chromesome("chr3"),chromesome("chr4"),chromesome("chr5"),chromesome("chr6"),chromesome("chr7"),chromesome("chr8"),chromesome("chr9"),chromesome("chr10"),chromesome("chr11"),chromesome("chr12"),chromesome("chr13"),chromesome("chr14"),chromesome("chr15"),chromesome("chr16"),chromesome("chr17"),chromesome("chr18"),chromesome("chr19"),chromesome("chr20"),chromesome("chr21"),chromesome("chr22"),chromesome("chrX"),chromesome("chrY")]
        self.lstOfChromosomesClusters=theChromsomeList
        self.lstOfChromosomesClustersGenes=theChromsomeList2
        #self.lstOfChromosomesClustersGenes=[chromesome("chr1"),chromesome("chr2"),chromesome("chr3"),chromesome("chr4"),chromesome("chr5"),chromesome("chr6"),chromesome("chr7"),chromesome("chr8"),chromesome("chr9"),chromesome("chr10"),chromesome("chr11"),chromesome("chr12"),chromesome("chr13"),chromesome("chr14"),chromesome("chr15"),chromesome("chr16"),chromesome("chr17"),chromesome("chr18"),chromesome("chr19"),chromesome("chr20"),chromesome("chr21"),chromesome("chr22"),chromesome("chrX"),chromesome("chrY")]
        self.CDSDistribution = []
        self.fDistribution = []
        self.tDistribution = []
        self.totalList=[]
        self.allCluster=[]
              
    def __str__(self):
        toPrint="{Parclip Object} Name: "+self.filename+"\tNumber of Clusters: "+str(self.getNumClusters())+"\tNumber of Cluster Genes: "+str(self.getNumGenes())
        return(toPrint)         
    
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
                        
    def writePC(self,outputDir):
        """
        Writes the clusters in a parclip to a csv file for later use. Writes other properties to a text file.
        Input:
            outputDir(str) - Directory where the csv form of clusters should be outputed.
        """
        #As CSV
        toRet=[["Chr","Start","End","Gene","Strand"]]
        #Create a list of clusters
        for chrom in self.lstOfChromosomesClustersGenes:
            for indG in chrom:
                for indClust in indG.lstOfClusters:
                    toAdd=[self.chrom,self.start,self.end,self.gene,self.strand]
                    toRet.append(toAdd)
        df=pd.DataFrame(toRet)
        toSave=outputDir+"/"+self.filename+"_internal_clusters.csv"
        df.to_csv(toSave,header=False,columns=False)
        
        #Write the other properties ot a text file
        toOpen=outputDir+"/"+self.filename+"_internal_properties.txt"
        f1=open(toOpen,"w+")
        f1.write("-------------PAR-CLIP OBJECT PROPERTIES----------\n")
        f1.write("File: "+self.filename+"\n")
        f1.write("Number of genes: "+str(self.getNumGenes)+"\n")
        f1.write("Number of clusters: "+str(self.getNumClusters)+"\n")
        f1.close()
        print("PC Written.")
    
    def getNumClusters(self):
        """
        Returns the number of clusters within this parclip object.
        Returns:
            int - The total count of clusters contained in this parclip object.
        """
        toReturn=0
        for ch in self.lstOfChromosomesClustersGenes:
            for gene in ch.lstOfAnotElements:
                toReturn=toReturn+len(gene.lstOfClusters)
        return(toReturn) 
          
    def getNumGenes(self):
        """
        Returns the number of genes that are overlapped by clusters within this parclip object.
        Returns:
            int - The total count of genes that are overlapped by clusters in this parclip object.
        """
        toReturn=0
        for ch in self.lstOfChromosomesClustersGenes:
            for gene in ch.lstOfAnotElements:
                if (len(gene.lstOfClusters)>0):
                    toReturn=toReturn+1
        return(toReturn)    
                       
    def populateCh(self):
        """
        Populate the chromosome lists with the clusters and clustergenes. This is effectively a sorting step to improve speed by reducing required comparisons downstream.
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
    
      
    def populateAllHits(self,anotF,randomStatesList):
        """
        Populate all of the hit and percentage for every cluster gene in this parclip object.
        Inputs:
            anotFile(anotFile) - An anotFile object containing the genes to check
            randomStatesList(List[num]) - A list of random seeds to use for the binning process. The length of the list corresponds to the number of binning loops that will be performed.
        """
        #For each chromosome
        #/////
        counter=0
        for chrom in self.lstOfChromosomesClustersGenes:
            #Get each gene
            for cGene in chrom.lstOfAnotElements:
                #Populate              
                cGene.populateHitLists(anotF)
                cGene.populateNumClust(anotF)
                cGene.populatePercentiles(anotF,randomStatesList)
                counter=counter+1
      
    def getParclipPercentageDistributions(self,weightByRC,wtGE,outDir):
        """
        Populate the average distribution of all of the percentage hits in this parclip object - essentially the step where we combine the many individual genes for a parclip by averaging.
        Input:
            weightByRC (Bool) - Whether to weight the impact of each gene on the metagene by the number of genes reads aligned to the gene.
            wtGE (WTGeneExpression or False) - Whether to weight the impact of each gene on the metagene by the expression level of that gene in the wild type cell line.
            outDir (The directory where the output log for those genes with zero expression should be kept).
        """        
        #Get the name of the output log for those clusters not found in the wildtype
        toSaveZeroes=outDir+"/genesNotInWildType.csv"
        if wtGE!=False:
            f= open(toSaveZeroes,"a+")
            f.write("Gene, ReadCount\n")
        #For each percentile bin
        for i in range(0,100):
            lstOfTValues=[]
            lstOfFValues=[]
            lstOfCValues=[]
            #For each gene cluster 
                #For each chromosomes
            for chrom in self.lstOfChromosomesClustersGenes:
                #Get each gene
                for gc in chrom.lstOfAnotElements:
                        if len(gc.gene)==0:
                            continue
                        #Save the genes that are not foudn in the wild type HEK293 expression levels.
                        if wtGE!=False:
                            if gc.gene not in wtGE.geDict or wtGE.geDict[gc.gene]==0:
                                #Save
                                f.write(gc.gene+", "+str(gc.getTotalNumReads())+"\n")
                            
                        if (len(gc.percentD3)==100 and len(gc.percentD5)==100 and len(gc.percentDCDS)==100) :
                            #If we are NOT weighing by RC or WT Expression levels
                            if weightByRC==False and wtGE==False:
                                #   Get the associated 3', cds, and 5' values
                                lstOfTValues.append(gc.percentD3[i])
                                lstOfFValues.append(gc.percentD5[i])
                                lstOfCValues.append(gc.percentDCDS[i])
                            
                            #If we ARE weighing by RC but not by wild type gene expression
                            if weightByRC==True and wtGE==False:
                                for k in range(0,gc.getTotalNumReads()):
                                    lstOfTValues.append(gc.percentD3[i])
                                    lstOfFValues.append(gc.percentD5[i])
                                    lstOfCValues.append(gc.percentDCDS[i])
                            
                            #If we are weighing by BOTH wild type gene expression AND read count
                            if wtGE!=False:
                                #Get the value for the TPM expression
                                if gc.gene in wtGE.geDict:
                                    #Save the zero expression ones
                                    if wtGE.geDict[gc.gene]!=0:
                                        countUp=int(round(wtGE.geDict[gc.gene]))
                                        if countUp!=0:
                                            #Get the factor to multiply each by (#CLR/TPM)
                                            multFact=gc.getTotalNumReads()/countUp
                                            #Apply the list that number of times
                                            lstOfTValues.append(gc.percentD3[i]*multFact)
                                            lstOfFValues.append(gc.percentD5[i]*multFact)
                                            lstOfCValues.append(gc.percentDCDS[i]*multFact)
                        
            #Get the average
            toAdd3 = Average(lstOfTValues)
            toAddC = Average(lstOfCValues)
            toAdd5 = Average(lstOfFValues)  

            #Add to the saved lists
            self.fDistribution.append(toAdd5)   
            self.CDSDistribution.append(toAddC)   
            self.tDistribution.append(toAdd3)  
                 
        #Compress for FUTR and CDS~~~~~~~~~~~~~~~~~~~~
            #----------FUTR
        lstOf11Bins=[]
        #Loop overthe bins
        bProgress=9.0909
        nxtBinCarryOver=0
        binTotal=0
        for pBin in self.fDistribution:
            #The fraction of a 100 bin to get
            toGet=0
            #Get an amount of the bin set by bProgress
            if bProgress>1:
                toGet=1
            if bProgress<=1:
                toGet=bProgress
            
            #Now add this amount to the binTotal
            binTotal=binTotal+pBin*toGet+nxtBinCarryOver
            
            nxtBinCarryOver=0
            
            #Now adjust the progress amount
            bProgress=bProgress-toGet
            
            #If bProgress gets to zero
            if bProgress==0:
                #If that is all of the 
                #Add the bin total
                lstOf11Bins.append(binTotal)
                #Reset the bProgress (for both bin total and for the progress through an individual )
                binTotal=0
                bProgress=9.0909 - (1-toGet)
                
                #Get the amount to carry over
                nxtBinCarryOver=(1-toGet)*pBin
        
        self.fDistribution=lstOf11Bins
        
            #-----------CDS
        lstOf70Bins=[]
        #Loop overthe bins
        bProgress=1.4285
        nxtBinCarryOver=0
            #The amount a 70 bin is equal to   
        binTotal=0
        for pBin in self.CDSDistribution:   
            #The fraction of a 100 bin to get
            toGet=0
            #Get an amount of the bin set by bProgress
            if bProgress>1:
                toGet=1
            if bProgress<=1:
                toGet=bProgress
            
            #Now add this amount to the binTotal
            binTotal=binTotal+pBin*toGet+nxtBinCarryOver
            
            nxtBinCarryOver=0
            
            #Now adjust the progress amount
            bProgress=bProgress-toGet
            
            #If bProgress gets to zero
            if bProgress==0:
                #If that is all of the 
                #Add the bin total
                lstOf70Bins.append(binTotal)
                #Reset the bProgress for both the counting ibn and for the total amount stored.
                binTotal=0
                bProgress=1.4285 - (1-toGet)
                
                #Get the amount to carry over
                nxtBinCarryOver=(1-toGet)*pBin
        
        self.CDSDistribution=lstOf70Bins
            #~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        print("---------------- Final FUTR distributions: ",self.fDistribution)
        print("---------------- Final CDS distributions: ",self.CDSDistribution)
        print("---------------- Final TUTR distributions: ",self.tDistribution)
        
        print("AVG of FUTR",Average(self.fDistribution))
        print("AVG of CDS",Average(self.CDSDistribution))
        print("AVG of TUTR",Average(self.tDistribution))
    
         
    def getTotalList(self):
        """
        Returns a single list which contains the percentage distribution of the 5'UTR, CDS, and 3'UTR appended together in that order.
        Returns:
            List[int] - The merged list of percentage distributions.
        """
        for ele in self.fDistribution:
            self.totalList.append(ele)
            
        for ele in self.CDSDistribution:
            self.totalList.append(ele)
            
        for ele in self.tDistribution:
            self.totalList.append(ele)
        
        return(self.totalList)     

    def getNum5UTRGenes(self):
        """
        Returns the number of genes with clusters in the 5' region in this PARCLIP.
        Returns:
            int - The number of genes in this object with clusters overlapping their 5'region
        """
        toRet=0
        for indivChrom in self.lstOfChromosomesClustersGenes:
            for cGene in indivChrom.lstOfAnotElements:
                #Check if it has density in the 3' region
                check=sum(cGene.hitList5)
                #If it does
                if check>0:
                    #Count
                    toRet=toRet+1
        return(toRet)
        
    def getNumCDSGenes(self):
        """
        Returns the number of genes with clusters in the CDS region in this PARCLIP.
        Returns:
            int - The number of genes in this object with clusters overlapping their CDS region
        """
        toRet=0
        for indivChrom in self.lstOfChromosomesClustersGenes:
            for cGene in indivChrom.lstOfAnotElements:
                #Check if it has density in the 3' region
                check=sum(cGene.hitListCDS)
                #If it does
                if check>0:
                    #Count
                    toRet=toRet+1
        return(toRet) 

    def getNum3UTRGenes(self):
        """
        Returns the number of genes with clusters in the 3' region in this PARCLIP.
        Returns:
            int - The number of genes in this object with clusters overlapping their 3'region
        """
        toRet=0
        for indivChrom in self.lstOfChromosomesClustersGenes:
            for cGene in indivChrom.lstOfAnotElements:
                #Check if it has density in the 3' region
                check=sum(cGene.hitList3)
                #If it does
                if check>0:
                    #Count
                    toRet=toRet+1
        return(toRet) 
        
    def getNum5UTRClusters(self):
        """
        Returns the number of cluster overlapping the 5' region in this PARCLIP.
        Returns:
            int - The number of clusters in this object overlapping their 5'region
        """
        toRet=0
        for ch in self.lstOfChromosomesClustersGenes:
            for gene in ch.lstOfAnotElements:
                toRet = toRet+ gene.FClusterCount
        return(toRet)

    def getNumCDSClusters(self): 
        """
        Returns the number of cluster overlapping the CDS region in this PARCLIP.
        Returns:
            int - The number of clusters in this object overlapping their CDSregion
        """
        toRet=0
        for ch in self.lstOfChromosomesClustersGenes:
            for gene in ch.lstOfAnotElements:
                toRet = toRet+ gene.CDSClusterCount
        return(toRet)
        
    def getNum3UTRClusters(self):
        """
        Returns the number of cluster overlapping the 3' region in this PARCLIP.
        Returns:
            int - The number of clusters in this object overlapping their 3'region
        """
        toRet=0
        for ch in self.lstOfChromosomesClustersGenes:
            for gene in ch.lstOfAnotElements:
                toRet = toRet+ gene.TClusterCount
        return(toRet)
        

    def applyBounds(self,lowerBound,upperBound,field):
        """
        This method goes through and removes clusters from the PAR-CLIP object that are at either % end of the sorted spectrum of all clusters for an inputed property (read count, cross linked reads, etc). Rounds down through Int() typecasting at the final stage.
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
        self.allClusters=allClusters        
        
        
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
            sortedList=sorted(self.allClusters, key=lambda x: int(x.start))

        if field=="end":
            print("Sorting by end location for filter.")
            sortedList=sorted(self.allClusters, key=lambda x: int(x.end))

        if field=="CS":
            print("Sorting by CS (conversion specificity) for filter...")
            sortedList=sorted(self.allClusters, key=lambda x: int(x.CS))
            
        if field=="URC":
            print("Sorting by URC (unique read count) for filter...")
            sortedList=sorted(self.allClusters, key=lambda x: int(x.URC))
            
        if field=="RC":
            print("Sorting by RC (read count) for filter...")
            sortedList=sorted(self.allClusters, key=lambda x: int(x.RC))
        
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
        
        
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
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
        dictSearch - Takes a dictionary representation of a cluster and searches for the gene (if any) it aligns to within this annotation object.
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
        Determine the total number of annotation rows stored in the chromosome objects of the object.
        Inputs:
            int - the total number of rows.
        """
        toReturn=0
        for ch in self.lstOfChromosomesRows:
            toReturn=toReturn+len(ch.lstOfAnotElements)
        return(toReturn)
                          
    def getNumGenes(self):
        """
        Determine the total number of annotation genes stored in the chromosome objects of the object.
        Inputs:
            int - the total number of genes
        """
        toReturn=0
        for ch in self.lstOfChromosomesGenes:
            toReturn=toReturn+len(ch.lstOfAnotElements)
        return(toReturn)    
                
    def addAR(self,rowObj):
        """
        Add an annotation row to the appropriate chromosome object in the anotation file object.
        """
        #Get the correct chromosome
        curID = rowObj.chrom
        for chrom in self.lstOfChromosomesRows:
            if (chrom.ID==curID):
                chrom.lstOfAnotElements.append(rowObj)
                return()
                                 
    def addAG(self,anotGeneObj):
        """
        Add an annotation gene to the appropriate chromosome object in the anotation file object.
        """ 
        curID = anotGeneObj.chrom
        for chrom in self.lstOfAnotGenes:
            if(chrom.ID==curID):
                chrom.lstOfAnotElements.append(anotGeneObj)
                return()
        
    def populateCh(self): 
        """
        This function goes through all of the annotation row objects associated with this object and adds them to their respective anotGene objects associated with this object (or makes new genes if need be) - all stored in their appropriate chromosome objects.
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
                        if (int(indivRow.start)<int(indivGene.start)):
                            indivGene.start=indivRow.start
                        if (int(indivRow.stop)>int(indivRow.stop)):
                            indivGene.stop=indivRow.stop
                        counter=1
                        break
                if counter==0:
                    #Make a new anotGene
                    nAnotGene=anotGene(indivRow.gene,indivRow.chrom)
                    nAnotGene.addAR(indivRow)
                    nAnotGene.start=indivRow.start
                    nAnotGene.stop=indivRow.stop
                    gChr.lstOfAnotElements.append(nAnotGene)
                    
    def dictSearch(self, dictionary):
        """
        Searches for the gene that a given cluster (in dictionary form) is aligning to.
        Input:
            dictionary (dict) - A dictionary representation of a cluster with the chr, start, and stop values.
        Output:
            gene (str) - The gene in the annotation file that this dictionary aligns to
        """
        #Create an options dictionary (gene name key and percentage overlap value)
        geneOptions={}
        
        #Match the chromosome
        dChrom=dictionary.get("Chr")
        for indChr in self.lstOfChromosomesGenes:
            if indChr.ID==dChrom:
                #Loop over the genes
                for anotG in indChr.lstOfAnotElements:
                    anotG.getStartandStopIndices()
                    #Get the overlaps if there is any overlap
                    toCheck=False
                    if int(anotG.start) <= int(dictionary.get("Start")) and int(dictionary.get("Start"))<=int(anotG.stop):
                        toCheck=True
                    if int(anotG.start) <= int(dictionary.get("End")) and int(dictionary.get("End"))<=int(anotG.stop):
                        toCheck=True
                    if int(anotG.start)<=int(dictionary.get("Start")) and int(anotG.stop)<=int(dictionary.get("End")):
                        toCheck=True
                    if toCheck==True:
                        theOverlap=anotG.getDictOverlapPercentage(dictionary)
                        if theOverlap>0.0:
                            geneOptions[anotG.gene]=theOverlap
        
        if len(geneOptions)==0:
            return("")                
        #Select the largest overlap
        theLargestVal=0
        theLargestGene=""
        for key in geneOptions:
            value=geneOptions.get(key)
            if value  >= theLargestVal:
                theLargestGene=key
                theLargestVal=value
        if theLargestGene not in geneListMaster:
            geneListMaster.append(theLargestGene)
        return(theLargestGene)
        
    def listAllGenes(self):
        """
        Generate a list of all the genes in this annotation object by name.
        Returns:
            List[Str]
        """
        toRet=[]
        for indChr in self.lstOfChromosomesGenes:
            for anotG in indChr.lstOfAnotElements:
                if anotG.gene not in toRet:
                    toRet.append(anotG.gene)
        return(toRet)
                            
##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Functions
##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ 
def getNumericalOverlap(refStart,refStop,testStart,testStop):
    """
    Get the overlap of two number ranges as a percentage of the test range.
    Input:
        refStart (int) - Start of the reference range.
        refStop (int) - End of the reference range.
        testStart (int) - Start of the test range.
        testStop (int) - End of the test range.
    """
    refStart=int(refStart)
    refStop=int(refStop)
    testStart=int(testStart)
    testStop=int(testStop)
    theCount=0
    anOverlap=False
    #Check if there is an overlap at all
    if refStart <= testStart and testStart<=refStop:
        anOverlap=True
        
    if refStop>=testStop and testStop>= refStart:
        anOverlap=True
        
    if testStart<=refStart and testStop>=refStop:
        anOverlap=True
    
    if anOverlap==False:
        return(0)
        
    #Test indices
    for i in range(testStart,testStop+1):
        if i in range(refStart,refStop+1):
            theCount=theCount+1
    percentage=float(theCount)/(refStop-refStart)*100.0
    return(percentage)
    
    
def getBedHeaderNames(bedFile):
    """
    Creates a list of the header names from a given bed file.
    Input:
        bedFile (str) - Bed file to be analyzed.
    Returns:
        Tuple[Str] - The column names.
    """    
    toRet=[]
    #Get an example row
    readCSV = csv.reader(open(bedFile),delimiter='\t')
    testRow=[]
    for row in readCSV:
        testRow=row
        break
    
    #Populate the optional fields
    optionalFields=["Chr","Start","End","name","score","Strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts"]
    
    for i in range(0,len(testRow)):
        toRet.append(optionalFields[i])
    return(tuple(toRet))


def fetchGene(gtfObject,inpDict):
    """
    Determine the gene within a gtf file that a given dictionary representation of a cluster aligns to.
    Input:
        gtfObject (gtfObject) - Annotation object to search.
        inpDict (Dict) - Dictionary representation of a cluster.
    """    
    return(gtfObject.dictSearch(inpDict))

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
    
    print("Unique Chromsomes: ", uniqueChromosomes)
    
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
    print("Finished Populating")
    return(myAFileObject)


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
        nNum = random.randint(0,bound)
        if nNum not in lstOfRandNums:
            lstOfRandNums.append(nNum)
    return(lstOfRandNums)     


def createExportTable(lstOfPCs,myDir,genesToUse):
    """
    This function takes in a list of parclip objects and returns a table which contains the distribution matrices for each of these
    Input:
        num(int) - The number of random numbers to generate.
        bound(int) - The upper bound of possible value generation.
    Returns:
        List[int] - A list of unique random numbers between 0 and the inputed upper bound of length num.
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
        label5 = pc.filename+": "+"5' UTR "
        labelCDS = pc.filename+": "+"CDS "
        label3 = pc.filename+": "+"3' UTR "
            #Get the 5'
        DFP = pc.fDistribution
            #Add the label
        DFP.insert(0,label5)
        toRetLst.append(DFP)
            #Get the CDS
        DCDSP = pc.CDSDistribution
            #Add the label
        DCDSP.insert(0,labelCDS)
        toRetLst.append(DCDSP)
            #Get the 5'B
        DTP = pc.tDistribution
            #Add the label
        DTP.insert(0,label3)
        toRetLst.append(DTP)
    
      #Export as csv--------
    csvfile = myDir
    #Add name to the end of directory
    if len(genesToUse)==0:
        csvfile=csvfile+"/HeatMapTableAdjustedBins.csv"
    elif len(genesToUse)!=0:
        csvfile=csvfile+"/"+genesStrToAdd+"FUTR_CDS_TUTR_Table.csv"

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
    

def createSubPlots(exportTable,outputDir,genesToUse,lstOfPCs,inpDPI,outFrmt):
    """
    This function takes in an export table of distribution values (output from the createExportTable function) and a series of subplots representing the data saved in a single figure to the specified directory.
    Input:
        exportTable(List[List[]]]) - A list of lists, with each interior list describing the distribution of clsuters from a signle parclip experiment. Output from createExportTable.
        outputDir(String) - The output directory the figure should be saved in.
        genesToUse(List[String]) - A list of the gene names for the genes to be considered. Input an empty list if all genes are to be considered.
        lstOfPCs(List[Parclip]) - A list of PARCLIP objects which are described in the exportTable.
    """ 
    #Grab the domain row
    #First create the header for percentiles
    theDomain = []
    
    genesStrToAdd=""
    for item in genesToUse:
        genesStrToAdd=genesStrToAdd+item+"_"
    
    for x in range(1,182):
        theDomain.append(x)   
        

    fig=plt.figure()
    indic=0
    for w in indicesToCheck(exportTable):        
        #Gene count lables
        num5UTRG=lstOfPCs[indic].getNum5UTRGenes()
        num5UTRGDisp="Genes: "+str(num5UTRG)
        numCDSG = lstOfPCs[indic].getNumCDSGenes()
        numCDSGDisp="Genes: "+str(numCDSG)
        num3UTRG = lstOfPCs[indic].getNum3UTRGenes()
        num3UTRGDisp="Genes: "+str(num3UTRG)
        
        #Cluster count labels
        num5UTRC=lstOfPCs[indic].getNum5UTRClusters()
        num5UTRCDisp="Clusters: "+str(num5UTRC)
        numCDSC=lstOfPCs[indic].getNumCDSClusters()
        numCDSCDisp="Clusters: "+str(numCDSC)
        num3UTRC=lstOfPCs[indic].getNum3UTRClusters()
        num3UTRCDisp="Clusters: "+str(num3UTRC)        
        
        
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
            #Get the vertical values
        vertVals = getHorizontalRange(mVal)
        nAx.plot([11.5,11.5,11.5,11.5,11.5],vertVals,linestyle='dashed',color='b')
        nAx.plot([81.5,81.5,81.5,81.5,81.5],vertVals,linestyle='dashed',color='b')
        
        #Add labels for gene and cluster counts
        nAx.text(0, mVal*.94, num5UTRGDisp,fontsize=7,color='m')
        nAx.text(36, mVal*.94, numCDSGDisp,fontsize=7,color='m')
        nAx.text(126, mVal*.94, num3UTRGDisp,fontsize=7,color='m') 
        nAx.text(0, mVal*.88, num5UTRCDisp,fontsize=7,color='m')
        nAx.text(36, mVal*.88, numCDSCDisp,fontsize=7,color='m')
        nAx.text(126, mVal*.88, num3UTRCDisp,fontsize=7,color='m')
        
        #Add the plot
        nAx.plot(theDomain, theRange, color='g')
        nAx.set_ylabel(sectionTitle, color='g',fontsize=13)
        nAx.set_title("Distribution of 5'UTR, CDS, and 3'UTR Sites",y=1.1,fontsize=13)
        nAx.axes.get_xaxis().set_visible(False)
                
        indic=indic+1
        if (indic==len(indicesToCheck(exportTable))):
            nAx.axes.text(0, mVal*1.04, "5'", fontsize=12,color='b')
            nAx.axes.text(36, mVal*1.04, "CDS", fontsize=12,color='b')
            nAx.axes.text(126, mVal*1.04, "3'", fontsize=12,color='b')
    
    if len(genesToUse)==0:
        if outFrmt=="png":
            toSave='DistributionPlot_5_CDS_3.png'
        elif outFrmt=="pdf":
            toSave='DistributionPlot_5_CDS_3.pdf'
    elif len(genesToUse)!=0:
        if outFrmt=="png":
            toSave=genesStrToAdd+"DistributionPlot_5_CDS_3.png"
        elif outFrmt=="pdf":
            toSave=genesStrToAdd+"DistributionPlot_5_CDS_3.pdf"
    
    os.chdir(outputDir)
    #old dpi = 300
    fig.savefig(toSave, bbox_inches='tight',dpi=inpDPI)
    

def outputAvgLengths(lstOfPCs,myDir):
    """
    This function outputs a CSV file which shows the average TUTR / CDS / FUTR gene lengths for all parclips within a parclip list object
    Input:
        lstOfPCs(List[parclip]) - The list of parclip objects to be analyzed.
        myDir(String) - Output directory to write the csv file to.
    Returns:
        List[List[]] - A list of lists containing the average gene lengths of each group for each gene in the parclips (the list saved as a csv).
    """  
    #Create the list to export
    toWrite=[]
    
    #For each parclip
    for pc in lstOfPCs:
        lengthsList3=[]
        lengthsList5=[]
        lengthsListCDS=[]
        for chrom in pc.lstOfChromosomesClustersGenes:
                #Get each gene
                for gc in chrom.lstOfAnotElements:
                    if len(gc.hitList3)!=0:
                        toAdd1 = gc.hitList3
                        lengthsList3.append(len(toAdd1))
                    if len(gc.hitList5)!=0:
                        toAdd1 = gc.hitList5
                        lengthsList5.append(len(toAdd1))
                    if len(gc.hitListCDS)!=0:
                        toAdd1 = gc.hitListCDS
                        lengthsListCDS.append(len(toAdd1))
                        
                        
        #Get the average lengths
        FUTRAvg = Average(lengthsList5)
        TUTRAvg = Average(lengthsList3)
        CDSAvg = Average(lengthsListCDS)
        
        #Get the multiplicative factor
        MF = 100/TUTRAvg
        
        #Multiply each of these and round to get the approximate bin count
        abcF=round(FUTRAvg*MF)
        abcT=round(TUTRAvg*MF)
        abcCDS=round(CDSAvg*MF)
        
        #Now create the seven lists needed before exporting 
        toAdd1 = [pc.filename]
        toAdd2 = ["Avg 5'UTR Length",FUTRAvg]
        toAdd3 = ["Avg CDS' Length",CDSAvg]
        toAdd4 = ["Avg 3'UTR Length",TUTRAvg]
        toAdd5 = ["# 5'UTR Bins",abcF]
        toAdd6= ["# CDS Bins",abcCDS]
        toAdd7=["# 3' UTR Bins",abcT]
                
        toWrite.append(toAdd1)
        toWrite.append(toAdd2)
        toWrite.append(toAdd3)
        toWrite.append(toAdd4)
        toWrite.append(toAdd5)
        toWrite.append(toAdd6)
        toWrite.append(toAdd7)
            
        #Export as csv--------
    csvfile = myDir
    #Add name to the end of directory
    csvfile=csvfile+"/BinLengthsIntermediate.csv"

    #Assuming res is a list of lists
    with open(csvfile, "w") as output:
        writer = csv.writer(output, lineterminator='\n')
        writer.writerows(toWrite)  
    
    return(toWrite)


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


def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Function to create a heatmap. Adapted from documentation
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
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom",fontsize="15")


    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    ax.set_aspect(2.5) # X scale matches Y scale
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels,fontsize=15)
    ax.set_yticklabels(row_labels,fontsize=18)


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
    ax.set_title("5'/CDS/3' Distribution Heat Map")
    

    return im, cbar


def createHeatMap(lstOfParClips,outputDir,genesToUse,inpDPI,imgFrmt):
    """
    Creates a heat map of metagene distribution data for all parclips considered where the intensities are scaled to be relative to the greatest intensity value in each distribution.
    Input:
        lstOfParClips(List[parclip]) - A list of parclip objects to consider.
        outputDir(String) - Directory where the heatmap is to be saved.
        genesToUse(List[String]) - List of the gene names to be used. If all are desired it is an empty list.
    """ 
    genesStrToAdd=""
    for item in genesToUse:
        genesStrToAdd=genesStrToAdd+item+"_"
        
    #Establish the x axis labels
    xAxisLabels=[]
    for x in range(0,181):
        if x==11:
            xAxisLabels.append("CDS Start")
        elif x==81:
            xAxisLabels.append("CDS Stop")
        elif x!=11 and x!=81:
            xAxisLabels.append("")
            
    #Establish the y catagory labels
    yAxisLabels=[]
    for pc in lstOfParClips:
        yAxisLabels.append(pc.filename)
    
    #Create the matrix of values
    valueMatrix=[]
    for pc in lstOfParClips:
        #Create the row to add to the value matrix
        ro=pc.fDistribution+pc.CDSDistribution+pc.tDistribution
        
        #Normalize this row
        nRo = scaleList(ro)
            
        valueMatrix.append(nRo)
        
    #Convert this matrix to a pandas dataframe
    valueDF = pd.DataFrame(valueMatrix,index=yAxisLabels)

    s_heatmap(valueDF,outputDir,xAxisLabels,genesToUse,inpDPI,imgFrmt)
    
    
def createHeatMapSingleton(lstOfParClips,outputDir,genesToUse,inpDPI,imgFrmt):
    """
    Creates a heat map of metagene distribution data for all parclips considered where the intensities are scaled to be relative to the greatest intensity value in each distribution.
    Input:
        lstOfParClips(List[parclip]) - A list of parclip objects to consider.
        outputDir(String) - Directory where the heatmap is to be saved.
        genesToUse(List[String]) - List of the gene names to be used. If all are desired it is an empty list.
    """ 
    genesStrToAdd=""
    for item in genesToUse:
        genesStrToAdd=genesStrToAdd+item+"_"
        
    #Establish the x axis labels
    xAxisLabels=[]
    for x in range(0,181):
        if x==11:
            xAxisLabels.append("CDS Start")
        elif x==81:
            xAxisLabels.append("CDS Stop")
        elif x!=11 and x!=81:
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
        ro=pc.fDistribution+pc.CDSDistribution+pc.tDistribution
        
        #Normalize this row
        nRo = scaleList(ro)
        valueMatrix.append(nRo)
        
    #Convert this matrix to a numpy array
    valueArray = np.array(valueMatrix)

    #Creates the plot
    fig, ax = plt.subplots()
    im, cbar = heatmap(valueArray, yAxisLabels, xAxisLabels, ax=ax,cmap="coolwarm", cbarlabel="Coverage")
    fig.tight_layout()
    
    os.chdir(outputDir)
    
    #Establish the height of the image - for each parclip there is .6 height
    height = .6
    fig.set_size_inches(60, height)
    
    
    if len(genesToUse)==0:
        if imgFrmt=="pdf":
            toSave='DistributionHeatMap_5_CDS_3.pdf'
        elif imgFrmt=="png":
            toSave='DistributionHeatMap_5_CDS_3.png'
    elif len(genesToUse)!=0:
        if imgFrmt=="pdf":
            toSave=genesStrToAdd+"DistributionHeatMap_5_CDS_3.pdf"
        elif imgFrmt=="png":
            toSave='DistributionHeatMap_5_CDS_3.png'
    
    fig.savefig(toSave, bbox_inches='tight',dpi=inpDPI)
    
    
def s_heatmap(pdDF,outDir,xAxisLab,genesToUse,inpDPI,imgFrmt):
    """
    This function generates and saves a clustered heatmap using seaborn's cluster functionality and dendrogram visualization techinique.
    Inputs:
        pdDF (Pandas DataFrame) - The dataframe of the metagene results to visualize. Row indices are the name of the metagene.
        outDir (Str) - Output directory in which the heatmap is saved.
        xAxisLab (List[Str]) - A list of the values that constitute the x axis labels of the heatmap.
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
    plt.title("5'/CDS/3' Cluster Distribution", fontsize = 50, loc='center')


    #Save the plot to the output directory
    if len(genesToUse)==0:
        if imgFrmt=="pdf":
            toSave='DistributionHeatMap_5_CDS_3.pdf'
        elif imgFrmt=="png":
            toSave='DistributionHeatMap_5_CDS_3.png'
    elif len(genesToUse)!=0:
        if imgFrmt=="pdf":
            toSave=genesStrToAdd+"DistributionHeatMap_5_CDS_3.pdf"
        elif imgFrmt=="png":
            toSave='DistributionHeatMap_5_CDS_3.png'
    g.savefig(toSave, bbox_inches='tight',dpi=inpDPI)
    
    
##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------         
def run(gtfFileName,csvDir,outputDir,genesToUse,lowerBound,upperBound,boundFilter,weighValuesByRC,dpi,imgFormat,randomStatesLst,wildTypeGE,mainChromosomesOnly,bedAnotGTF,anoterScript):
    """
    The primary running function for the program - it takes in the parclip data and runs all of the metegene calculations to produce distribution values, distribution plots, and a heatmap plot (all saved)
    Input:
        gtfFileDir(String) - The path to the directory which has the modified GTF csv file.
        gtfFileName(String) - The name of the modified GTF csv file to use.
        csvDir(String) - The path to the directory containing the cluster.csv files to be analyzed (if running in directory mode) or a csv file (if you are running in single csv mode).
        outputDir(String) - The path to the directory in which all files will be saved.
        genesToUse(List[String]) - A list of the genes to consider in the calculations. To consider all genes enter an empty list.
        sort(Boolean) - Do you want to include the filtering step which removes the 75% of clusters sorted to have the lowest read counts.
        lowerBound (num) - The lower percentage of clusters to throw out based on the boundFilter property.
        uppBound (num) - The upper percentage of clusters to throw out based on the boundFilter property.
        boundFilter(Str) - The property used to sort the clusters in a data set before the bounds are applied to remove clusters based on this property.
            Options: "start" (start coordinate of cluster),"end" (end coordinate of cluster),"CS" (Conversion Specificity of cluster),"URC" (unique read count of the cluster),"RC" (read count of the cluster)
        weighValuesByRC (Bool) - Apply a weight to each cluster percentage distribution based on the number of reads aligned to the cluster?
        randomStatesList (List[int]) - The list of numbers to use as random seeds in random number generation. The more random states provided, the more times the binning algorithm will execute before the averaging step.
        dpi(int) - The dpi to save images to.
        imgFormat (str) - "pdf" or "png". The format to save images into.
        wildTypeGE (str) - Name of a wild type gene expression file which contain "Average" and "Gene Symbol" columns.
        mainChromosomesOnly (bool) - Whether the program should only look at the "main" chromsomes (autosomal plus X/Y) for analysis.
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
        #Generate hits
        exPC.populateAllHits(myGTFObj,randomStatesLst)   
        #Generate percentages                
        exPC.getParclipPercentageDistributions(weighValuesByRC,wtGE_File,outputDir)
    print("...populated.")
        
    print("Creating the export table...")
    exportTable = createExportTable(lstOfParclips,outputDir,genesToUse) 
    print("...export table created.")

    print("Creating the binding distribution curves for all parclips in the directory...")
    createSubPlots(exportTable,outputDir,genesToUse,lstOfParclips,dpi,imgFormat)
    print("...binding distribution curves generated and merged.")
    #binningTable = outputAvgLengths(lstOfParclips,outputDir) 

    print("Generating the heat map...")
    #If it is a heatmap of more than one row then it can cluster
    if len(lstOfParclips)>1:
        createHeatMap(lstOfParclips,outputDir,genesToUse,dpi,imgFormat)
        print("Heat map created",genesToUse)       
    #Otherwise a single heatmap row should be not clustered.
    elif len(lstOfParclips)==1:
        createHeatMapSingleton(lstOfParclips,outputDir,genesToUse,dpi,imgFormat)
    print("... heat map generated.")
    
    print("The genes pulled: ")
    print(geneListMaster)
##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Run
##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
#run(gtfFileName,csvDir,outputDir,gtu,lBound,uBound,boundF,wValuesByRC,dpi,imgFormat,randStates,wtGE,mainChromosomes,theBED_GTF,theBEDAnotScript)