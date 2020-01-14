# clip-metagene
Suite of tools that use cluster files to visualize the binding of CLIP data across genomic features of a specified species.

## Longest Transcript Annotation File
The metagene programs do not directly use standard (GTF, RefFlat, etc.) annotation files but CSV files containing information about only the longest transcript (the transcript with the most exonic nucleotides) splice variant of each gene in an annotation. Included in this repository are such CSV files generated from the Gencode human ([gencode.v30.annotation.gtf](https://www.gencodegenes.org/human/release_30.html)) and mouse ([gencode.vM23.annotation.gtf](https://www.gencodegenes.org/mouse/release_M23.html)) chromosome assemblies. If users wish to run a metagene analysis using a different annotation source (likely whichever was used in generating their cluster files), they can provide a gtf file as the annotation input and the program will automatically generate the intermediate csv annotation file from the provided gtf file, save it to the output directory, and use it to run the metagene analyses. *Note: This step is very slow and thus a gtf file should only be used as an input when a CSV annotation file has not yet been generated.*
### Algorithm Explained
The algorithm which generates the longest transcript annotation file from a gtf file operates as follows:
1.	The GTF file is read into the program and each row is sorted into a gene by name (commented header lines are ignored). Note: Gene ID can also be used for sorting by running the writeGTFGeneList function with geneID=True.
2.	The length (number of exonic sequence nucleotides) of each possible splice variant per gene is calculated and the longest is selected to represent that gene. If no coding sequence exists, then the transcript with the most nucleotides of any classification is selected.
3.	Based on their location in the transcript, UTR regions are redefined as FUTR, TUTR, or introns.
4.	Intron rows are generated between exons within each transcript.
5.	Export to CSV.

## 2. 5' UTR / CDS / 3' UTR Metagene
