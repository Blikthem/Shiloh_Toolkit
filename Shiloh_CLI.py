#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 15:09:03 2019

@author: claypooldj
"""
##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Load dependencies
##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
import argparse
import os
import sys

    #Scripts within the same directory as this script
import MakeLongestGTF9 as MLG
import IntronExonHM_22 as HM_IE
import HeatMapper28 as HM_FCDST
import IntronExonHM_21_newWindow_5 as HM_IE2

##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Create the parser and take input
##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
my_parser = argparse.ArgumentParser(prog='metagene',description='Perform a metagene analysis on binding cluster files.')

# Add the required arguments (inputs)
my_parser.add_argument('cInput',metavar='clusterInput',type=str,help='The input cluster file(s). Either a .csv file or a directory (in this case the analysis will be run on all csv files in that directory).')
my_parser.add_argument('outDir',metavar='outputDirectory',type=str,help='The directory in which all output files will be generated.')
my_parser.add_argument('anot',metavar='anotation',type=str,help="The annotation file matching that used in cluster formation. This program runs on intermediate files (hg38 and mm10 included in download) generated from a gtf file deposited in the source directory. Pass that intermediate file or a gtf file. If a gtf file is given, the program will automatically generate the intermediate file from that gtf file, place it in the source directory, and then run the program (Note: This is quite slow so don't pass a gtf if you already have the intermediate file).")

# Add the optional arguments (inputs)
my_parser.add_argument('-gtu',metavar="genesToUse",type=list,help="A list of genes ex [gene1,gene2,gene3], the program will then ONLY consider overlap with those genes in all analysis. Default=Analyze all genes")
my_parser.add_argument('-bf',metavar="boundFilter",type=str,help="The name of the continuous variable to use when filtering clusters. Must be 'start', 'end', 'CS', 'URC', or 'RC'.  Default=ReadCount.")
my_parser.add_argument('-lb',metavar="lowerBound",type=int,help="The lower bound for the filtering step - the percentage of clusters with the lowest values of the boundFilter variable to remove. Default=0")
my_parser.add_argument('-ub',metavar="upperBound",type=int,help="The upper bound for the filtering step - the percentage of clusters with the highest values of the boundFilter variable to remove. Default=0")
my_parser.add_argument('-thresh',metavar="clusterThreshold",type=int,help="The minimum number of clusters that must be aligned to a gene in order for it to be considered in the analysis. Default=1") 
my_parser.add_argument('-rs',metavar="randomStates",type=list,help="A list of integers to use as the random seeds during bin analysis. Increasing this list size will increase the number of times the binning step is done before the results are averaged. Default=[7211995,541995,3131994,111,222,333,444,555,888,999]")
my_parser.add_argument('-dpi',metavar="dpi",type=int,help="Dots per inch (image resolution) of all graphics created. Default=250")
my_parser.add_argument('-at',metavar="analysisType",type=str,help="Whether to only run one of the metagene analysis programs. Enter IE to run only the intron/exon program or FCDST to only run the 5 prime / CDS /3 prime analyais. Default=Run both")
my_parser.add_argument('-wtexp',metavar="wildTypeExpressionFile",type=str,help="A file containing gene expression levels for the wild type cell line. If included, the program will normalize gene expression levels to the rounded average wild type gene expression.")
my_parser.add_argument('-bedGTF',metavar="bedAnnotation",type=str,help="The GTF to be used to annotate bed files without Gene fields. Note: A time-consuming indexing step occurs the first time a new GTF is used. Default=Gencode hg38.")
my_parser.add_argument('--png',help="Save images as PNG's instead of PDF's.")
my_parser.add_argument('--wbrc',help="Weight the impact of every normalized gene cluster distribution on the overall metagene by the number of reads aligning to that gene.")
my_parser.add_argument('--unusualChromosomes',help="Consider clusters that align to chromosomes other than the 'main' autosomal and sex chromosomes.")
my_parser.add_argument('--exonCentered',help="Whether to also run an exon-centric version of the intron/exon heat mapper. This feature is not as rigorously tested and also takes a lot of time and resources.")
# Execute the parse_args() method
args = my_parser.parse_args()

#--------------------------------------------------------------------------------------------------

#Check to make sure that the output directory is indeed a directory.
if not os.path.isdir(args.outDir):
    print('The output directory you specified ('+args.outDir+') does not exist.')
    sys.exit()

#Check
print("--------------------RUNNING CLIP METAGENE--------------------")
print("Positional Inputs:")
print("Cluster input: ",args.cInput)
print("Output directory: ",args.outDir)
print("Anotation File: ",args.anot)
print("-------------------------------------------------------------")

##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Establish the variables needed to run from the positional values
##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
csvDir=args.cInput    #The directory containing the cluster csv files to analyze. All csv files in this directory will be tested (so make sure there aren't other csv files in the directory or it will crash!) {Type - String}

        #The directory where you want to generate all output files. {Type - String}
outputDir=args.outDir+"/metageneOutput"

if not os.path.exists(outputDir):
    os.mkdir(outputDir)    

##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Establish the annotation file containing longest transcript information. If a gtf is given, make the requisite gtf intermediate. Otherwise, continue.
##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if '.gtf' in args.anot:
    #Run the annotation creating program on the gtf file
    print("----------GENERATING 5'/CDS/3' METAGENE----------")
    
    anotBase=os.path.basename(args.anot)
    anotBase=anotBase[0:len(anotBase)-4]
    
    #Create the intermediate file
    MLG.writeGTFGeneList(args.anot,outputDir,anotBase,False)
    
    gtfFileName=outputDir+"/"+anotBase+".longestTranscript.csv" #The intermediate gtf longest splice varient csv file (generated from a seperate script) to use for annotation. {Type - String}
    print("-------------------------------------------------")
    
else:
    gtfFileName=args.anot  

##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Establish the variables needed to run (and populate default values for these variables if needed) from the optional inputs.
##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ 
if args.bf==None:
    boundF="URC"
else:
    boundF=args.bf

if args.lb==None:
    lBound=0    
else:
    lBound=args.lb

if args.ub==None:
    uBound=0
else:
    uBound=args.ub
    
if args.gtu==None:
    gtu=[]    
else:
    gtu=args.gtu

if args.wbrc:
    wValuesByRC=True 
else:
    wValuesByRC=False
    
if args.dpi==None:
    dpi=100
else:
    dpi=args.dpi

if args.png:
    imgFormat="png"
else:
    imgFormat="pdf"    

if args.rs==None:
    randStates=[7211995,541995,3131994,111,222,333,444,555,888,999] 
else:
    randStates=args.rs
    
if args.thresh==None:
    clustThreshold=1    
else:
    clustThreshold=args.thresh

if args.unusualChromosomes:
    mainChromosomes=False
else:
    mainChromosomes=True

if args.wtexp==None:
    wtGE=""
else:
    wtGE=args.wtGE

if args.exonCentered:
    toRunNew=False
else:
    toRunNew=True

#Bed annotation information
if args.bedGTF==None:
    inp_bedGTF=os.getcwd()
    inp_bedGTF=inp_bedGTF+"/gencode.v30.annotation.gtf"

#Bed annotator script information
bedAnotScript=os.getcwd()
bedAnotScript=bedAnotScript+"/bedannotator.sh"

##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Run 5' UTR / CDS / 3' UTR Analysis
##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------     
if args.at!="IE":
    HM_FCDST.run(gtfFileName,csvDir,outputDir,gtu,lBound,uBound,boundF,wValuesByRC,dpi,imgFormat,randStates,wtGE,mainChromosomes,inp_bedGTF,bedAnotScript)

##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Run Intron / Exon Analysis
##------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ 
if args.at!="FCDST":
    HM_IE.run(gtfFileName,csvDir,outputDir,1250,100,gtu,clustThreshold,lBound,uBound,boundF,wValuesByRC,randStates,dpi,imgFormat,wtGE,mainChromosomes,inp_bedGTF,bedAnotScript)
if toRunNew==True:    
    HM_IE2.run(gtfFileName,csvDir,outputDir,1250,100,gtu,clustThreshold,lBound,uBound,boundF,wValuesByRC,randStates,dpi,imgFormat,wtGE,mainChromosomes,13,False,inp_bedGTF,bedAnotScript)