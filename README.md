# clip-metagene
Suite of tools that use cluster files to visualize the binding of CLIP data across genomic features of a specified species.

## Longest Transcript Annotation
The metagene programs do not directly use standard (GTF, RefFlat, etc.) annotation files but CSV files containing information about only the longest transcript (the transcript with the most coding nucleotides) splice variant for each gene in an annotation. Included in this repository are such CSV files generated from the Gencode human ([gencode.v30.annotation.gtf](https://www.gencodegenes.org/human/release_30.html)) and mouse ([gencode.vM23.annotation.gtf](https://www.gencodegenes.org/mouse/release_M23.html)) chromosome assemblies. If users wish to run a metagene analysis using a different annotation source (likely whichever was used in generating their cluster files), they can provide a gtf file as the annotation input and the program will automatically generate the intermediate csv annotation file from the provided gtf file, deposit it in the output directory, and use it to run the metagene analyses. *Note: This step is very slow and thus a gtf file should only be used as an input when a CSV annotation file is unavailable.*
### Algorithm Explained

## 2. 5' UTR / CDS / 3' UTR Metagene
