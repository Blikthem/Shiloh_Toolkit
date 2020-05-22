# clip-metagene
Suite of tools that use cluster files to visualize the binding of CLIP data across genomic features of a specified species.

Perform a metagene analysis on binding cluster files.

positional arguments:
  clusterInput          The input cluster file(s). Either a .csv file or a
                        directory (in this case the analysis will be run on
                        all csv files in that directory).
  outputDirectory       The directory in which all output files will be
                        generated.
  anotation             The annotation file matching that used in cluster
                        formation. This program runs on intermediate files
                        (hg38 and mm10 included in download) generated from a
                        gtf file deposited in the source directory. Pass that
                        intermediate file or a gtf file. If a gtf file is
                        given, the program will automatically generate the
                        intermediate file from that gtf file, place it in the
                        source directory, and then run the program (Note: This
                        is quite slow so don't pass a gtf if you already have
                        the intermediate file).

optional arguments:
  -h, --help            show this help message and exit
  -gtu genesToUse       A list of genes ex [gene1,gene2,gene3], the program
                        will then ONLY consider overlap with those genes in
                        all analysis. Default=Analyze all genes
  -bf boundFilter       The name of the continuous variable to use when
                        filtering clusters. Must be 'start', 'end', 'CS',
                        'URC', or 'RC'. Default=ReadCount.
  -lb lowerBound        The lower bound for the filtering step - the
                        percentage of clusters with the lowest values of the
                        boundFilter variable to remove. Default=0
  -ub upperBound        The upper bound for the filtering step - the
                        percentage of clusters with the highest values of the
                        boundFilter variable to remove. Default=0
  -thresh clusterThreshold
                        The minimum number of clusters that must be aligned to
                        a gene in order for it to be considered in the
                        analysis. Default=1
  -rs randomStates      A list of integers to use as the random seeds during
                        bin analysis. Increasing this list size will increase
                        the number of times the binning step is done before
                        the results are averaged. Default=[7211995,541995,3131
                        994,111,222,333,444,555,888,999]
  -dpi dpi              Dots per inch (image resolution) of all graphics
                        created. Default=250
  -at analysisType      Whether to only run one of the metagene analysis
                        programs. Enter IE to run only the intron/exon program
                        or FCDST to only run the 5 prime / CDS /3 prime
                        analyais. Default=Run both
  -wtexp wildTypeExpressionFile
                        A file containing gene expression levels for the wild
                        type cell line. If included, the program will
                        normalize gene expression levels to the rounded
                        average wild type gene expression.
  -bedGTF bedAnnotation
                        The GTF to be used to annotate bed files without Gene
                        fields. Note: A time-consuming indexing step occurs
                        the first time a new GTF is used. Default=Gencode
                        hg38.
  --png PNG             Save images as PNG's instead of PDF's.
  --wbrc WBRC           Weight the impact of every normalized gene cluster
                        distribution on the overall metagene by the number of
                        reads aligning to that gene.
  --unusualChromosomes UNUSUALCHROMOSOMES
                        Consider clusters that align to chromosomes other than
                        the 'main' autosomal and sex chromosomes.
  --exonCentered EXONCENTERED
                        Whether to also run an exon-centric version of the
                        intron/exon heat mapper. This feature is not as
                        rigorously tested and also takes a lot of time and
                        resources.

## Longest Transcript Annotation File
The metagene programs do not directly use standard (GTF, RefFlat, etc.) annotation files but CSV files containing information about only the longest transcript (the transcript with the most exonic nucleotides) splice variant of each gene in an annotation. Included in this repository are such CSV files generated from the Gencode human ([gencode.v30.annotation.gtf](https://www.gencodegenes.org/human/release_30.html)) and mouse ([gencode.vM23.annotation.gtf](https://www.gencodegenes.org/mouse/release_M23.html)) chromosome assemblies. If users wish to run a metagene analysis using a different annotation source (likely whichever was used in generating their cluster files), they can provide a gtf file as the annotation input and the program will automatically generate the intermediate csv annotation file from the provided gtf file, save it to the output directory, and use it to run the metagene analyses. *Note: This step is very slow and thus a gtf file should only be used as an input when a CSV annotation file has not yet been generated.*
### Algorithm
The algorithm which generates the longest transcript annotation file from a gtf file operates as follows:
1.	The GTF file is read into the program and each row is sorted into a gene by name (commented header lines are ignored). Note: Gene ID can also be used for sorting by running the writeGTFGeneList function with geneID=True.
2.	The length (number of exonic sequence nucleotides) of each possible splice variant per gene is calculated and the longest is selected to represent that gene. If no coding sequence exists, then the transcript with the most nucleotides of any classification is selected.
3.	Based on their location in the transcript, UTR regions are redefined as FUTR, TUTR, or introns.
4.	Intron rows are generated between exons within each transcript.
5.	Export to CSV.

## 5' UTR / CDS / 3' UTR Metagene
### Algorithm
1.	The annotation csv file is read into the program.
a.	An anotFile object is instantiated to store annotation information.
b.	The program scans each row of the intermediate annotation csv file.
i.	The data in the row is saved to an anotRow object.
ii.	The row object is added to the anotFile object.

2.	The annotation information is sorted by gene.
- The rows from the intermediate annotation file are sorted into anotGene objects, which contain all the annotation information from the GTF for a given gene.

3.	The cluster files are read into the program.
a.	For each .csv file in the input directory…
Note: Be sure that all the csv files in the input directory have the correct headings and fields of a cluster file.
i.	A new parclip object is created and named after the file (minus extension).
ii.	Each row in the inputted csv file is saved to a cluster object.
iii.	Each cluster object is added to the parclip object.
iv.	If the user opts to filter clusters by a continuous variable, that filtering occurs here (all the clusters in the parclip are sorted by the trait and then the clusters falling in the selected extreme ranges are removed).
v.	The remaining clusters are sorted into the parclip’s chromosome objects for orderly retrieval downstream.

4.	Hit lists (lists specifying which nucleotides in a genetic feature are overlapped by at least one cluster and which are not) are generated.
a.	For each clusterGene
i.	The coordinates of each cluster in the gene are intersected with the anotFile object to figure out which nucleotides in the gene are overlapped by a cluster. The result is a list of 1’s (overlapped nucleotides) and 0’s (not overlapped nucleotides) for the gene’s 5’ UT, coding, and 3’ UT, regions.
ii.	The total number of overlapped nucleotides in the gene is saved for use downstream.

5.	The hit lists for each region are binned into a new list of fixed length, displaying the percentage of the region (5’ UT, coding, or 3’ UT) that falls within this fraction of the genetic feature. Allows for direct comparisons (and averages) between genes.
a.	For each clusterGene
1.	Divide the length of the 5’UT, coding, and 3’ UT regions of the gene by the specified number of bins for each region in order to divide a region into equal length bins. When a remainder (r) exists in a region, the surplus in length is dealt with by increasing the length of r random bins by one. This is done a number of times equal to the length of the random seeds listinputted (default 10) and the bin lists generated by all random states (below) are averaged.
2.	For each random seed…
a.	In each bin the number of hits (nucleotides in this region overlapped by a cluster) is multiplied by the percentage of the total number of hits in this intron-exon junction. This value is added to the bin list.
3.	Average the bin lists from each random seed.

6.	Averaging all gene distribution lists to get the average binding profile for the PAR-CLIP experiment.
a.	The average 5’ UTR, coding region, and 3’ UTR cluster distributions for all of the genes overlapped by at least one cluster are averaged to get a metagene-level depiction of the binding properties of the RNA binding protein studied.

7.	Spreadsheet created.
a.	The average 5’ UTR, coding region, and 3’ UTR distribution of clusters are saved to a csv file.

8.	Distribution graphs are generated.
a.	For each parclip in the input directory…
i.	The average distribution of binding sites from the PAR-CLIP across the 5’ UT, coding, and 3’ UT regions is displayed as a curve on a plot. The x-axis represents the length of the gene features divided into bins and the y-axis represents the percent of the total cluster overlap that falls within each bin. As such, the sum of the bins always equals 100.

9.	The heatmap is generated.
a.	Compress all of the individual distribution plots into a single row of the heat map. Blue color indicates little cluster density in a bin while red color indicates greater cluster density.
b.	Employ agglomerative hierarchical clustering to group the rows of the heatmap.

## Intron / Exon UTR Metagene
### Algorithm

1.	The annotation csv file is read into the program.
a.	An anotFile object is instantiated to store annotation information.
b.	The program scans each row of the intermediate annotation csv file.
i.	The data in the row is saved to an anotRow object.
ii.	The row object is added to the anotFile object.

2.	The annotation information is sorted by gene.
a.	The rows from the intermediate annotation file are sorted into anotGene objects, which contain all the annotation information from the GTF for a given gene.

3.	The cluster files are read into the program.
a.	For each .csv file in the input directory…
Note: Be sure that all the csv files in the input directory have the correct headings and fields of a cluster file.
i.	A new parclip object is created and named after the file (minus extension).
ii.	Each row in the inputted csv file is saved to a cluster object.
iii.	Each cluster object is added to the parclip object.
iv.	If the user opts to filter clusters by a continuous variable, that filtering occurs here (all the clusters in the parclip are sorted by the trait and then the clusters falling in the selected extreme ranges are removed).
v.	The remaining clusters are sorted into the parclip’s chromosome objects for orderly retrieval downstream.

4.	Hit lists (lists specifying which nucleotides in a genetic feature are overlapped by at least one cluster and which are not) are generated.
a.	For each clusterGene
i.	The coordinates of each cluster in the gene are intersected with the anotFile object to figure out which nucleotides in the gene are overlapped by a cluster. The result is a list of 1’s (overlapped nucleotides) and 0’s (not overlapped nucleotides) for each intron and exon in the gene.
ii.	The total number of overlapped nucleotides in each intron-exon junction (half of exon1 + full intron + half of exon 2) is saved for use downstream.

5.	The hit lists for each region are binned into a new list of fixed length, displaying the percentage of the total intron-exon junction overlap that falls within this fraction of the genetic feature. Allows for direct comparisons (and averages) between genes.
a.	For each clusterGene
i.	For each intron-exon junction (half of exon1 + full intron + half of exon 2).
1.	Divide the length of the intronic and exonic regions of the gene by the number of bins desired for each region in order to divide a region into equal length bins. When a remainder (r) exists, the surplus in length is dealt with by increasing the length of r random bins by one. This is done a number of times equal to the length of random seeds list inputted (default 10) and the bin lists generated by all random states (below) are averaged.
a.	For each random seed…
i.	In each bin the number of hits (nucleotides in this region overlapped by a cluster) is multiplied by the percentage of the total number of hits in this intron-exon junction. This value is added to the bin list.
b.	Average the bin lists from each random seed.

6.	The gene-level cluster distribution lists are calculated by averaging individual intron-exon junction bin lists.
a.	For each clusterGene
i.	Average the intronic bin lists for all intron-exon junction to get the average distribution of clusters across the intronic regions of this gene.
ii.	Average the exonic bin lists for all intron-exon junction to get the average distribution of clusters across the exonic regions of this gene.

7.	Averaging all gene-level distribution lists to get the average binding profile for the PAR-CLIP experiment.
a.	The average intronic and exonic cluster distributions for all of the genes overlapped by at least one cluster are averaged to get a metagene-level depiction of the binding properties of the RNA binding protein studied.

8.	Spreadsheet created.
a.	The average intron/exon distribution of clusters are saved to a csv file.

9.	Distribution graphs are generated.
a.	For each parclip in the input directory…
i.	The average distribution of binding sites from the PAR-CLIP across the intronic/exonic junctions displayed as a curve on a plot. The x-axis represents the length of the gene features divided into bins and the y-axis represents the percent of the total cluster overlap that falls within each bin. As such, the sum of the bins always equals 100.

10.	The heatmap is generated.
a.	Compress all of the individual distribution plots into a single row of the heat map. Blue color indicates little cluster density in a bin while red color indicates greater cluster density.
b.	Employ agglomerative hierarchical clustering to group the rows of the heatmap.

