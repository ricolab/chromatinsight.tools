# chromatinsight.tools
R tools for Chromatinsight.

A set of methods in R used in tandem with the output of Chromatinsight https://github.com/ricolab/Chromatinsight

## Install

from GitHub, using devtools:\
```devtools::install_github("ricolab/chromatinsight.tools")```\
If devtools is not installed, then\
```install.packages("devtools")```\
You might need first:\
```sudo apt-get install libgit2-dev```

## Basic usage

You can test some of the features using the following simple example (the files in the input folder are the same as those in the output folder of Chromatinsight's unit test). The output file should be identical to the file in ```./utest/output_expected/dimorphismMono_output.txt``` . Note that this code assumes connection to the server of Ensembl).

```
library(chromatinsight.tools)
library(biomaRt)

myInputFolder = "./utest/input"
myOutputFolder = "./utest/output"

gene_types = as.data.frame(read.table(header = TRUE, text = "
gene_type
3prime_overlapping_ncrna
antisense
lincRNA
miRNA
misc_RNA
processed_transcript
protein_coding
pseudogene
rRNA
sense_intronic
sense_overlapping
snoRNA
snRNA
"))

ensembl75 = useMart(host='feb2014.archive.ensembl.org',
						biomart='ENSEMBL_MART_ENSEMBL',
						dataset='hsapiens_gene_ensembl')
						
allBioMartChr = getBM(mart = ensembl75, attributes = c("ensembl_gene_id", "wikigene_name", "chromosome_name", "start_position", "end_position", "gene_biotype"))

# Loading both Chromatinsight outputs (observed and random)
dimorphismMono = getRegionData(myInputFolder, prefix = "mono", suffix = "real", totRandomStates = 11)
dimorphismMonoRnd = getRegionData(myInputFolder, prefix = "mono", suffix = "rnd", totRandomStates = 11)

# Calculating FDR
dimorphismMonoFdr = addFDR(dimorphismMono, dimorphismMonoRnd)

# Adding information from BioMart
dimorphismMonoMore = fillTableWithGeneTypes(dimorphismMonoFdr, allBioMartChr, gene_types)

# Saving output
setwd(myOutputFolder)
write.table(dimorphismMonoMore, "dimorphismMono_output.txt", sep = "\t", row.names = FALSE, quote = FALSE)
```

## Main methods

**getRegionData**\
Loads the output text files of Chromatinsight (currently Python only), using as output the median of the number of trials porformed (in the totRandomStates variable).

**addFDR**\
Generates a new table with an extra column displaying the FDR of the disparity.\
The input is two output tables from getRegionData:\
tableObserved: the table with the data including the result from the Chromatinsight analysis.\
tableRnd: the table with the same data, but randomising the category labels (in testPrediction, randomize = True).\
Example:\
```newTable = addFDR(tableObserved, tableRnd)```

**pileup**\
You give it a grouping file and key to find a set files of ChromHMM binaries and it adds up the "ones" and "zeroes" that are at the same position for the set of files, dividing by the number of files. So it gives the ratio of samples having a "one" at a specific position ie, the ratio of samples being H3K27ac or H3K4me1 at a ChromHMM bin. If the parameter "key" is omitted, then all the files will be used.\
Example:\
```dataMales = pileup(grouping = "grouping.txt", key = "male", direc = "/data", chrom = "chr10")```\
Will merge data from the files classified as "male" in grouping.txt (within the /data folder), which in the example are:\
```monocytes_S2901_mal_DATA_binary.txt```\
```monocytes_S3454_mal_DATA_binary.txt```\
```monocytes_S5715_mal_DATA_binary.txt```\
etc

**drawpile**\
You give it the frequency of a histone mark at every ChromHMM bin, ie, the output of two pileup functions (see pileup), and draws the graph comparing them.\
For H3K27ac the first group is red, the second is blue.\
For H3K4me1 the first group is pink, the second is forest green.\
To use UCSC Genome Browser coordinate format (such as chrX:10,000-25,000) see drawpilegb.

**drawpilegb**
Calls drawpile with UCSC Genome Browser coordinate format, such as chrX:10,000-25,000 (See drawpile and densitypile for more details).

**densitypile**\
You give it the output of two pileup function and it shows the graph comparing them.\
For H3K27ac the first group is red, the second is blue, when they overlap it's purple.\
For H3K4me1 the first group is pink, the second is forest green, when they overlap it's grey.\
To use UCSC Genome Browser coordinate format (such as chrX:10,000-25,000) see drawpilegb.

**fillTableWithGeneTypes**\
Merges the information from Biomart and the data from Chromatinsight, retrieved using getRegionData.
 
**compareDistributions**\
Displays the violin plots of the distributions in two analyses retrieved using getRegionData. Especially useful to compare random vs observed data.
