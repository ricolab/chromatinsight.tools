# chromatinsight.tools
R tools for Chromatinsight.

A set of methods in R used in tandem with the output of Chromatinsight https://github.com/ricolab/Chromatinsight

## Install:
from GitHub, using devtools:\
```devtools::github_install("chromatinsight.tools")```\
If devtools is not installed, then\
```install.packages("devtools")```\
You might need first:\
```sudo apt-get install libgit2-dev```

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
You give it the filename prefix to find a set files of ChromHMM binaries and it adds up the "ones" and "zeroes" that are at the same position for the set of files, dividing by the number of files. So it gives the ratio of samples having a "one" at a specific position, ie, the ratio of samples being H3K27ac or H3K4me1 at a ChromHMM bin.\
Example:\
```pileup(prefix = "mono*S*mal", direc = "/data/", chrom = "chr10")```\
Will merge data from files (within /data folder)\
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
