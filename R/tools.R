#################################
### Chromatinsight tools v1.3 ###
#################################
#
# A set of methods for R
# used in tandem with Chromatinsight
#
# Author: Marco Trevisan-Herraz, PhD
# Computational Epigenomics Laboratory
# Newcastle University, UK
# 2019-2021
#
# requires ggplot2, biomaRt
#
# Install:
# from GitHub, using devtools:
# devtools::github_install("chromatinsight.tools")
# If devtools is not installed, then
# install.packages("devtools")
# You might need first:
# sudo apt-get install libgit2-dev

# methods

#----------------------------------------------------------------------

#' Calls drawpile with UCSC Genome Browser coordinate system.
#' (such as chrX:10,000-25,000)
#' (See drawpile for more details).
#' @export
drawpilegb = function(datafem, datamal, coord = "", plusminus = 1000, sex = "both", histmod = "ac", window = 200, useDensity = TRUE, opacity = 0.5, highlight = "", label = "") {

	# note that datafem and datamal can only be part of a specific chromosome
	# so the chromosome is ignored, make sure you have got the appropriate chromosome!
	myCoord = splitCoords(coord)
	print(myCoord)
    
    chrom = myCoord$chrom
	myStart = floor(myCoord$start / window) - ceiling(plusminus / window)
	myEnd = ceiling(myCoord$end / window) + ceiling(plusminus / window)
	myPlus = myEnd - myStart
	print(paste0("Start: ", myStart * window))
	print(paste0("End: ", myEnd * window))
	print(paste0("Span: ", myPlus * window))
    
    highlightStart = 0
    highlightPlus = 0
    
    if (nchar(highlight) > 0) {
        myHighlight = splitCoords(highlight)
        highlightStart = floor(myHighlight$start / window)
        highlightPlus = ceiling(myHighlight$end / window) - highlightStart
        }
	if (useDensity) densitypile(datafem, datamal, start = myStart, plus = myPlus, sex = sex, histmod = histmod, chrom = chrom, opacity = opacity, highstart = highlightStart, highplus = highlightPlus, label = label)
    else drawpile(datafem, datamal, start = myStart, plus = myPlus, sex = sex, histmod = histmod, chrom = chrom)
}

#----------------------------------------------------------------------

#' To split the genomic coordinates
#' @examples
#' splitCoords(coords = "chrX:10,000-25,000")
#' is the same as
#' list(chrom = "chrX", start = "10000", end = "25000")
splitCoords = function(coords = "") {
    chrom = strsplit(gsub(",", "", coords), ":")[[1]][1]
    chromCoords = strsplit(gsub(",", "", coords), ":")[[1]][2]
    start = as.numeric(strsplit(chromCoords, "-")[[1]][1])
    end = as.numeric(strsplit(chromCoords, "-")[[1]][2])
    coordSplit = list(chrom = chrom, start = start, end = end)
    
    return(coordSplit)
    }

#----------------------------------------------------------------------

#' You give it the filename prefix to find a set files of ChromHMM binaries
#' and it adds up the "ones" and "zeroes" that are at the same position for the
#' set of files, dividing by the number of files.
#' So it gives the ratio of samples having a "one" at a specific position
#' ie, the ratio of samples being H3K27ac or H3K4me1 at a ChromHMM bin.
#' @export
#' @examples
#' pileup(prefix = "mono*S*mal", direc = "/data/", chrom = "chr10")
#' Will merge data from files (within /data folder)
#' monocytes_S2901_mal_DATA_binary.txt
#' monocytes_S3454_mal_DATA_binary.txt
#' monocytes_S5715_mal_DATA_binary.txt
#' etc
pileup = function(prefix = "", direc = "", chrom = ""){
	
	myChromosome = chrom
	if (nchar(myChromosome) < 4) myChromosome = paste0("chr", myChromosome)
	
	myDirec = direc
	myFilePath = paste0(myDirec, "*", prefix, "*", myChromosome, "_binary.txt")
	
	myFiles = Sys.glob(myFilePath)
	print(paste0("Checking ", length(myFiles), " files."))
	
	myList = c()
	total = 0
	for (myFile in myFiles){
		myTable = read.table(myFile, sep = "\t", header = TRUE, skip = 1)
		print(paste0(nrow(myTable), ", ", length(myTable), ": ", myFile))
		if(max(myTable) <= 1 & length(myTable) == 2){
			total = total + 1
			if (length(myList) == 0){
					myList = myTable
				} else {
					myList = myList + myTable
				}
		} else {
			print("Discarded.")
			}
		}
	
	resultList = myList / total
	
	return(resultList)
	
	}
    
#----------------------------------------------------------------------

#' You give it the frequency of a histone mark at every ChromHMM bin,
#' ie, the output of two pileup functions (see pileup).
#' and draws the graph comparing them.
#' One group is blue, the other is red, when they overlap it's purple.
#' To use UCSC Genome Browser coordinate system (such as chrX:10,000-25,000)
#' see drawpilegb.
#' @export
drawpile = function(datafem, datamal, start = 0, plus = 1000, sex = "both", histmod = "ac", chrom = ""){
	
	lineWidth = 1
	binSize = 200
	
	theEnd = nrow(datafem)
	end = start + plus
	if (end < start) end = start + 1000
	if (end > theEnd) end = theEnd
	myTitle = paste0(chrom, ":", format(start * binSize, scientific = FALSE, big.mark = ","), "-", format(end * binSize, scientific = FALSE, big.mark = ","))
	cat(myTitle)
	
	datafempart = datafem[start:end,]
	datamalpart = datamal[start:end,]
	
	myPlot <- ggplot2::ggplot(aes(x = as.numeric(row.names(datafempart)), y = 0), data = datafempart) + gglot2::geom_line() + ggplot2::ylim(0, 1) + ggplot2::scale_colour_manual(values = c("males, H3K27ac" = "blue", "females, H3K27ac" = "red",  "females, H3K4me1" = "purple", "males, H3K2me1" = "forest green"))
	myPlot <- myPlot + ggplot2::labs(title = myTitle) + ggplot2::xlab("bin in ChromHMM (1 bin = 200b)") + ggplot2::ylab("probability of histone modification")
	
	if (histmod == "ac" | histmod == "both"){
		if (sex == "fem" | sex == "both") myPlot <- myPlot + ggplot2::geom_line(aes(y = H3K27ac, color = 'females, H3K27ac'), data = datafempart, size = lineWidth)
		if (sex == "mal" | sex == "both") myPlot <- myPlot + ggplot2::geom_line(aes(y = H3K27ac, color = 'males, H3K27ac'), data = datamalpart, size = lineWidth)
		}
	
	if (histmod == "me1" | histmod == "both"){
		if (sex == "fem" | sex == "both") myPlot <- myPlot + ggplot2::geom_line(aes(y = H3K4me1, color = 'females, H3K4me1'), data = datafempart, size = lineWidth)
		if (sex == "mal" | sex == "both") myPlot <- myPlot + ggplot2::geom_line(aes(y = H3K4me1, color = 'males, H3K2me1'), data = datamalpart, size = lineWidth)
		}
		
	myPlot
	
	}
    
#----------------------------------------------------------------------

#' You give it the output of two pileup function
#' and it shows the graph comparing them.
#' @export
densitypile = function(datafem, datamal,
                    start = 0, plus = 1000,
                    highstart = 0, highplus = 0,
                    sex = "both",
                    histmod = "ac",
                    chrom = "",
                    opacity = 0.5,
                    label = ""){
	
	lineWidth = 0.25
	binSize = 200
	
	theEnd = nrow(datafem)
	end = start + plus
	if (end < start) end = start + 1000
	if (end > theEnd) end = theEnd
	myTitle = paste0(chrom, ":", format(start * binSize, scientific = FALSE, big.mark = ","), "-", format(end * binSize, scientific = FALSE, big.mark = ","))
	print(myTitle)
	
	datafempart = datafem[start:end,]
	datamalpart = datamal[start:end,]
	
	myPlot <- ggplot2::ggplot(aes(x = as.numeric(row.names(datafempart)) * binSize, y = 0), data = datafempart) + ggplot2::geom_line() + ggplot2::ylim(0, 1) + ggplot2::scale_fill_manual(values = c("males, H3K27ac" = "#0000ff", "females, H3K27ac" = "#ff0000",  "females, H3K4me1" = "purple", "males, H3K2me1" = "forest green"))
	myPlot <- myPlot + ggplot2::labs(title = myTitle) + ggplot2::xlab("bins in ChromHMM (1 bin = 200b)") + ylab("probability of histone modification")
    myPlot <- myPlot + ggplot2::theme(panel.background = element_rect(fill="white", color = "grey50", size=2), panel.grid.major = element_line(color = "grey",size=(0.2)))
    
    if (highstart > 0 & highplus > 0) {
        myPlot <- myPlot + ggplot2::geom_rect(aes(xmin = highstart * binSize, xmax = (highstart + highplus) * binSize, ymin = 0, ymax = 1), fill = "#ffff00", alpha = 0.005)
        }
	
	if (histmod == "ac" | histmod == "both"){
		if (sex == "fem" | sex == "both") myPlot <- myPlot + ggplot2::geom_density(aes(y = H3K27ac, fill = 'females, H3K27ac'), data = datafempart, size = lineWidth, stat = 'identity', alpha = opacity)
		if (sex == "mal" | sex == "both") myPlot <- myPlot + ggplot2::geom_density(aes(y = H3K27ac, fill = 'males, H3K27ac'), data = datamalpart, size = lineWidth, stat = 'identity', alpha = opacity)
		}
	
	if (histmod == "me1" | histmod == "both"){
		if (sex == "fem" | sex == "both") myPlot <- myPlot + ggplot2::geom_density(aes(y = H3K4me1, color = 'females, H3K4me1'), data = datafempart, size = lineWidth)
		if (sex == "mal" | sex == "both") myPlot <- myPlot + ggplot2::geom_density(aes(y = H3K4me1, color = 'males, H3K2me1'), data = datamalpart, size = lineWidth)
		}
    
    if (nchar(label) > 0) {
        myPlot <- myPlot + ggplot2::annotate("text", x = start * binSize, y = 0.9, label = label, family = "Ubuntu", size = 10, hjust = 0)
        }
    
	myPlot
	
	}

#----------------------------------------------------------------------

#' Loads the output text files of Chromatinsight (currently Python only),
#' Using as output the median of the number of trials porformed (in the
#' totRandomStates variable).
#' @export
getRegionData = function(myPath = "", histmod = "ac", prefix = "prefix", suffix = "", totRandomStates = 11) {
	
	globalOutput = vector()
	if (suffix != "") suffix = paste0("_", suffix)
	myFiles = list.files(myPath, pattern = paste0(prefix, "_", histmod, "_chr[X0-9]{1,2}", suffix, ".txt"), full.names = TRUE)

	for (i in 1:length(myFiles)) {
		myFile = as.data.frame(read.csv(myFiles[i], header = FALSE, sep = "\t", skip = 1))
		chrom = vector()
		init = vector()
		end = vector()
		region = vector()
		comments = vector()
		medians = vector()
		totGenes = vector()
		for (i in 1:nrow(myFile)) {
			firstElement = strsplit(as.character(myFile[1][i,]), "_|-")
			chrom = rbind(chrom, firstElement[[1]][1])
			init = rbind(init, as.numeric(firstElement[[1]][2]))
			end = rbind(end, as.numeric(firstElement[[1]][3]))
			region = rbind(region, firstElement[[1]][4])
			comments = rbind(comments, "-")
			totGenes = rbind(totGenes, 0)
			medians = rbind(medians, median(as.numeric(myFile[2:totRandomStates + 1][i,])))
		}
		outputData = data.frame(chrom = chrom, init = init, end = end, region = region, comments = comments, tot_genes = totGenes, disparity = medians)
		
		globalOutput = rbind(globalOutput, outputData)
	}
	
	return(globalOutput)
}

#----------------------------------------------------------------------

#' Merges the information from Biomart and the data from Chromatinsight,
#' retrieved using getRegionData.
#' @export
fillTableWithGeneTypes = function(data, BioMartChr, SDGeneTypes) {
geneTypeColumns = data.frame()
data$comments = as.character(data$comments)
for (i in 1:nrow(data)) {
	row = data[i,]
	mySubset = subset(BioMartChr, BioMartChr$end_position > row$init & BioMartChr$start_position < row$end & row$chrom == paste0("chr", BioMartChr$chromosome_name))
	myGeneTypes = mySubset[,6]
	myComments = mySubset[,2]
	
	commentConcat = ""
	totGenes = length(unique(mySubset[,1]))
	for (ensemblGene in unique(mySubset[,1])) {
		if(nchar(ensemblGene) > 0) commentConcat = paste(commentConcat, ensemblGene, sep = ", ")
	}
	
	if (nchar(commentConcat) > 2) commentConcat = substr(commentConcat, 3, nchar(commentConcat))
	geneTypeColumns = rbind(geneTypeColumns, countGeneTypes(myGeneTypes, SDGeneTypes))
	data[i,5] = commentConcat
	data[i,6] = totGenes
	}
	
	newData = cbind(data, geneTypeColumns)
	return(newData)
}

#----------------------------------------------------------------------

#' Counts the gene types for each of the
#' gene types dowloaded from Biomart.
countGeneTypes = function(myGeneTypes, SDGeneTypes) {
	result = data.frame()
	for (i in 1:nrow(SDGeneTypes)) {
		thisGeneTypeCount = sum(myGeneTypes == SDGeneTypes[i,])
		if(length(result) == 0) result = data.frame(thisGeneTypeCount)
		else result = cbind(result, thisGeneTypeCount)
	}
	names(result) = t(SDGeneTypes)
	return(result)
}

#----------------------------------------------------------------------

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