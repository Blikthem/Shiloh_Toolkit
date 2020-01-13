# clip-metagene
Suite of tools that use cluster files to visualize the binding of CLIP data across genomic features of a specified species.

## Longest Transcript Annotation
The metagene programs do not run directly from standard (GTF, RefFlat, etc.) annotation files but use CSV files containing information about only the longest transcript (the transcript with the most coding nucleotides) splice variant for each gene in an annotation. Included in this repository are examples of this type of CSV file generated from the Gencode human and mouse gtf assemblies (gencode.v30.annotation.gtf and gencode.vM23.annotation.gtf). If users wish to run a metagene analysis using a different annotation source (whichever was used in generating their cluster file), they can provide a gtf file as the annotation input and the program will automatically generate the intermediate csv annotation file (deposited in the specified output directory) from the provided gtf file and use it to run the program. Note: This step is very slow and thus should not be repeated when a CSV file is available.
###Algorithm Explained

## 5' UTR / CDS / 3' UTR Metagene
