##################################
### Chromatinsight tools v1.21 ###
##################################
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

#' Generates a new table with an extra column displaying the FDR of the
#' disparity.
#' The input is two output tables from getRegionData:
#' tableObserved: the table with the data including the result from the
#' Chromatinsight analysis.
#' tableRnd: the table with the same data, but randomising the category
#' labels (in testPrediction, randomize = True).
#' @export
#' @examples
#' newTable = addFDR(tableObserved, tableRnd)
addBayesFactor = function(regionDataObserved, hMark, dataRndFolder, dataRndPrefix, dataRndSuffix, totRandomStates = 11) {
	dataRnd = getAllRandomStates(myPath = dataRndFolder, hMark, prefix = dataRndPrefix, suffix = dataRndSuffix, totRandomStates = totRandomStates)
	resultTable = regionDataObserved
	resultTable$BayesFactor = NA
	
	
	for (i in 1:nrow(regionDataObserved)) {
		if (i %% 100 == 0 & verbose) print(paste0(i, "/", nrow(regionDataObserved)))
		thisDisparity = regionDataObserved$disparity[i]
		totObserved = nrow(subset(regionDataObserved, disparity >= thisDisparity))
		totRnd = nrow(subset(regionDataRnd, disparity >= thisDisparity))
		# The following is the same as FP / (FP + TP)
		# Because totRND only shows false positives, but totObs includes both
		thisRnd = totRnd / totObserved
		resultTable$FDR[i] = thisRnd
		}
	
	return(resultTable)
	
	}

#----------------------------------------------------------------------

#' Generates a new table with an extra column displaying the FDR of the
#' disparity.
#' The input is two output tables from getRegionData:
#' tableObserved: the table with the data including the result from the
#' Chromatinsight analysis.
#' tableRnd: the table with the same data, but randomising the category
#' labels (in testPrediction, randomize = True).
#' @export
#' @examples
#' newTable = addFDR(tableObserved, tableRnd)
addFDR = function(regionDataObserved, regionDataRnd, verbose = TRUE) {
	resultTable = regionDataObserved
	resultTable$FDR = NA
	for (i in 1:nrow(regionDataObserved)) {
		if (i %% 100 == 0 & verbose) print(paste0(i, "/", nrow(regionDataObserved)))
		thisDisparity = regionDataObserved$disparity[i]
		totObserved = nrow(subset(regionDataObserved, disparity >= thisDisparity))
		totRnd = nrow(subset(regionDataRnd, disparity >= thisDisparity))
		# The following is the same as FP / (FP + TP)
		# Because totRND only shows false positives, but totObs includes both
		thisRnd = totRnd / totObserved
		resultTable$FDR[i] = thisRnd
		}
	
	return(resultTable)
	
	}

#----------------------------------------------------------------------

#' Calls drawpile with UCSC Genome Browser coordinate format.
#' (such as chrX:10,000-25,000)
#' (See drawpile and densitypile for more details).
#' @export
drawpilegb = function(data1,
						data2,
						coord = "",
						plusminus = 1000,
						channel = "both",
						histmod = "ac",
						window = 200,
						useDensity = TRUE,
						opacity = 0.5,
						highlight = "",
						label = "",
						filename = "",
						filedim = c(1800, 1000),
						label1 = "first",
						label2 = "second"
						) {

	# note that data1 and data2 can only be part of a specific chromosome
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
	if (useDensity) densitypile(data1, data2, start = myStart, plus = myPlus, channel = channel, histmod = histmod, chrom = chrom, opacity = opacity, highstart = highlightStart, highplus = highlightPlus, label = label, filename = filename, filedim = filedim, label1 = label1, label2 = label2)
    else drawpile(data1, data2, start = myStart, plus = myPlus, channel = channel, histmod = histmod, chrom = chrom, filename = filename, filedim = filedim, label1 = label1, label2 = label2)
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

#' You give it a grouping file and key to find a set files of ChromHMM binaries
#' and it adds up the "ones" and "zeroes" that are at the same position for the
#' set of files, dividing by the number of files.
#' So it gives the ratio of samples having a "one" at a specific position
#' ie, the ratio of samples being H3K27ac or H3K4me1 at a ChromHMM bin.
#' If the parameter "key" is omitted, then all the files will be used.
#' @export
#' @example
#' pileup(grouping = "grouping_file.txt", key = "mal", chrom = "chr10")
pileup = function(grouping = "", key = "", direc = "", chrom = ""){
	
	myChromosome = chrom
	if (nchar(myChromosome) < 4) myChromosome = paste0("chr", myChromosome)
	
	myDirec = direc
	if (nchar(myDirec) > 0) {
		if (substr(myDirec, nchar(myDirec), nchar(myDirec)) != "/") {
			myDirec = paste0(myDirec, "/")
			}
		}
	groupingFile = paste0(myDirec, grouping)
	groupingTable = as.data.frame(read.table(groupingFile, sep = "\t", header = FALSE))
	colnames(groupingTable) = c("filename", "tablekey")
	
	myFiles = groupingTable[,"filename"]
	if (key != "") myFiles = subset(groupingTable, tablekey == key)[,"filename"]
	
	print(paste0("Checking ", length(myFiles), " files."))
	
	myList = c()
	total = 0
	for (myFile in myFiles){
		myFileChr = gsub("\\*", chrom, myFile)
		myTable = read.table(myFileChr, sep = "\t", header = TRUE, skip = 1)
		print(paste0(nrow(myTable), ", ", length(myTable), ": ", myFileChr))
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

#' Deprecated version of pileup that was using as a reference the file name
#' rather than the grouping file.
#' You give it the filename prefix to find a set files of ChromHMM binaries
#' and it adds up the "ones" and "zeroes" that are at the same position for the
#' set of files, dividing by the number of files.
#' So it gives the ratio of samples having a "one" at a specific position
#' ie, the ratio of samples being H3K27ac or H3K4me1 at a ChromHMM bin.
#' @export
#' @examples
#' pileup2(prefix = "mono*S*mal", direc = "/data/", chrom = "chr10")
#' Will merge data from files (within /data folder)
#' monocytes_S2901_mal_DATA_binary.txt
#' monocytes_S3454_mal_DATA_binary.txt
#' monocytes_S5715_mal_DATA_binary.txt
#' etc
pileup2 = function(prefix = "", direc = "", chrom = ""){
	
	myChromosome = chrom
	if (nchar(myChromosome) < 4) myChromosome = paste0("chr", myChromosome)
	
	myDirec = direc
	if (nchar(myDirec) > 0) {
		if (substr(myDirec, nchar(myDirec), nchar(myDirec)) != "/") {
			myDirec = paste0(myDirec, "/")
			}
		}
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
#' Include a filename parameter if you want to save the output into a file.
#' For H3K27ac the first group is red, the second is blue,
#' For H3K4me1 the first group is pink, the second is forest green,
#' To use UCSC Genome Browser coordinate format (such as chrX:10,000-25,000)
#' see drawpilegb.
#' @export
drawpile = function(data1,
					data2,
					start = 0,
					plus = 1000,
					channel = "both",
					histmod = "ac",
					chrom = "",
					filename = "",
					filedim = c(1800, 1000),
					label1 = "first",
					label2 = "second"
					){
	
	lineWidth = 1
	binSize = 200
	
	theEnd = nrow(data1)
	end = start + plus
	if (end < start) end = start + 1000
	if (end > theEnd) end = theEnd
	myTitle = paste0(chrom, ":", format(start * binSize, scientific = FALSE, big.mark = ","), "-", format(end * binSize, scientific = FALSE, big.mark = ","))
	cat(myTitle)
	
	data1part = data1[start:end,]
	data2part = data2[start:end,]
	
	theseColors = c("red", "blue", "purple", "forest green")
	color1ac = paste0(label1, ", H3K27ac")
	color2ac = paste0(label2, ", H3K27ac")
	color1me = paste0(label1, ", H3K4me1")
	color2me = paste0(label2, ", H3K4me1")
	names(theseColors) = c(color1ac, color2ac, color1me, color2me)
	
	myPlot <- ggplot2::ggplot(ggplot2::aes(x = as.numeric(row.names(data1part)), y = 0), data = data1part) + ggplot2::geom_line() + ggplot2::ylim(0, 1) + ggplot2::scale_colour_manual(values = theseColors)
	myPlot <- myPlot + ggplot2::labs(title = myTitle) + ggplot2::xlab("bin in ChromHMM (1 bin = 200b)") + ggplot2::ylab("probability of histone modification")
	
	if (histmod == "ac" | histmod == "both"){
		if (channel == 1 | channel == "both") myPlot <- myPlot + ggplot2::geom_line(ggplot2::aes(y = H3K27ac, color = color1ac), data = data1part, size = lineWidth)
		if (channel == 2 | channel == "both") myPlot <- myPlot + ggplot2::geom_line(ggplot2::aes(y = H3K27ac, color = color2ac), data = data2part, size = lineWidth)
		}
	
	if (histmod == "me1" | histmod == "both"){
		if (channel == 1 | channel == "both") myPlot <- myPlot + ggplot2::geom_line(ggplot2::aes(y = H3K4me1, color = color1me), data = data1part, size = lineWidth)
		if (channel == 2 | channel == "both") myPlot <- myPlot + ggplot2::geom_line(ggplot2::aes(y = H3K4me1, color = color2me), data = data2part, size = lineWidth)
		}
		
	if (filename != ""){
		ggplot2::ggsave(filename, plot = myPlot, width = filedim[1] / 150, height = filedim[2] / 150, units = "in", dpi = 150)
		}
	else print(myPlot)
	
	}
    
#----------------------------------------------------------------------

#' You give it the output of two pileup function
#' and it shows the graph comparing them.
#' Include a filename parameter if you want to save the output into a file.
#' For H3K27ac the first group is red, the second is blue,
#' when they overlap it's purple.
#' For H3K4me1 the first group is pink, the second is forest green,
#' when they overlap it's grey.
#' To use UCSC Genome Browser coordinate format (such as chrX:10,000-25,000)
#' see drawpilegb.
#' @export
densitypile = function(data1, data2,
                    start = 0, plus = 1000,
                    highstart = 0, highplus = 0,
                    channel = "both", # can be 1 for data1, 2 for data2, or "both"
                    histmod = "ac",
                    chrom = "",
                    opacity = 0.5,
                    label = "",
					filename = "",
					filedim = c(1800,1000),
					label1 = "first",
					label2 = "second"
					){
	
	lineWidth = 0.25
	binSize = 200
	
	theEnd = nrow(data2)
	end = start + plus
	if (end < start) end = start + 1000
	if (end > theEnd) end = theEnd
	myTitle = paste0(chrom, ":", format(start * binSize, scientific = FALSE, big.mark = ","), "-", format(end * binSize, scientific = FALSE, big.mark = ","))
	print(myTitle)
	
	data1part = data1[start:end,]
	data2part = data2[start:end,]
	
	theseColors = c("red", "blue", "purple", "forest green")
	color1ac = paste0(label1, ", H3K27ac")
	color2ac = paste0(label2, ", H3K27ac")
	color1me = paste0(label1, ", H3K4me1")
	color2me = paste0(label2, ", H3K4me1")
	names(theseColors) = c(color1ac, color2ac, color1me, color2me)
	
	myPlot <- ggplot2::ggplot(ggplot2::aes(x = as.numeric(row.names(data1part)) * binSize, y = 0), data = data1part) + ggplot2::geom_line() + ggplot2::ylim(0, 1) + ggplot2::scale_fill_manual(values = theseColors)
	myPlot <- myPlot + ggplot2::labs(title = myTitle) + ggplot2::xlab("bins in ChromHMM (1 bin = 200b)") + ggplot2::ylab("probability of histone modification")
    myPlot <- myPlot + ggplot2::theme(panel.background = ggplot2::element_rect(fill="white", color = "grey50", size=2), panel.grid.major = ggplot2::element_line(color = "grey",size=(0.2)))
    
    if (highstart > 0 & highplus > 0) {
        myPlot <- myPlot + ggplot2::geom_rect(ggplot2::aes(xmin = highstart * binSize, xmax = (highstart + highplus) * binSize, ymin = 0, ymax = 1), fill = "#ffff00", alpha = 0.005)
        }
	
	if (histmod == "ac" | histmod == "H3K27ac" | histmod == "both"){
		if (channel == 1 | channel == "both") myPlot <- myPlot + ggplot2::geom_density(ggplot2::aes(y = H3K27ac, fill = color1ac), data = data1part, size = lineWidth, stat = 'identity', alpha = opacity)
		if (channel == 2 | channel == "both") myPlot <- myPlot + ggplot2::geom_density(ggplot2::aes(y = H3K27ac, fill = color2ac), data = data2part, size = lineWidth, stat = 'identity', alpha = opacity)
		}
	
	if (histmod == "me1" | histmod == "me" | histmod == "H3K4me1" | histmod == "both"){
		if (channel == 1 | channel == "both") myPlot <- myPlot + ggplot2::geom_density(ggplot2::aes(y = H3K4me1, fill = color1me), data = data1part, size = lineWidth, stat = 'identity', alpha = opacity)
		if (channel == 2 | channel == "both") myPlot <- myPlot + ggplot2::geom_density(ggplot2::aes(y = H3K4me1, fill = color2me), data = data2part, size = lineWidth, stat = 'identity', alpha = opacity)
		}
    
    if (nchar(label) > 0) {
        myPlot <- myPlot + ggplot2::annotate("text", x = start * binSize, y = 0.9, label = label, family = "Ubuntu", size = 10, hjust = 0)
        }
    
	if (filename != ""){
		ggplot2::ggsave(filename, plot = myPlot, width = filedim[1] / 150, height = filedim[2] / 150, units = "in", dpi = 150)
		}
	else print(myPlot)
	
	}

#----------------------------------------------------------------------

#'
#' @export
getAllRandomStates = function(myPath = "", histmod = "ac", prefix = "prefix", suffix = "", totRandomStates = 11) {
	
	globalOutput = vector()
	if (suffix != "") suffix = paste0("_", suffix)
	myFiles = list.files(myPath, pattern = paste0(prefix, "_", histmod, "_chr[X0-9]{1,2}", suffix, ".txt"), full.names = TRUE)

	allRandomStates = c()
	for (i in 1:length(myFiles)) {
		myData = as.data.frame(read.csv(myFiles[i], header = FALSE, sep = "\t", skip = 1))
		randomStates = unlist(myData)
		allRandomStates = c(allRandomStates, randomStates)
		}
	
	return(allRandomStates)
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
#' To get data from BioMart, an option can be the following:
#' library(BioMart)
#' # This applies a filter to get only data from chrX
#' ensembl75 = useMart(host='feb2014.archive.ensembl.org',
#' 						biomart='ENSEMBL_MART_ENSEMBL',
#' 						dataset='hsapiens_gene_ensembl')
#' chrX = getBM(mart = ensembl75,
#' 			attributes = c("ensembl_gene_id",
#' 						"chromosome_name",
#' 						"start_position",
#' 						"end_position",
#' 						"gene_biotype"),
#' 			filters = "chromosome_name",
#' 			values = "X")
#' # This has no filters, so gets all the data
#' allBioMartChr = getBM(mart = ensembl75,
#' 					attributes = c("ensembl_gene_id",
#' 					"wikigene_name",
#' 					"chromosome_name",
#' 					"start_position","end_position",
#' 					"gene_biotype"))
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

#' Displays the violin plots of the distributions in two analyses
#' retrieved using getRegionData. Especially useful to compare random vs
#' observed data.
#' @export
compareDistributions = function(regionData1,
								regionData2,
								label1 = "first",
								label2 = "second",
								yLabel = "data",
								graphTitle = "comparison",
								filename = "",
								filedim = c(1200, 1200)) {

	myFirst = as.data.frame(c(first = regionData1$disparity))
	myFirst$id = label1
	rownames(myFirst) = NULL
	colnames(myFirst) = c("data", "id")

	mySecond = as.data.frame(c(second = regionData2$disparity))
	mySecond$id = label2
	rownames(mySecond) = NULL
	colnames(mySecond) = c("data", "id")

	mydf = rbind(myFirst, mySecond)

	vPlot_better = ggplot2::ggplot(mydf, ggplot2::aes(x = id, y = data, fill = id), ties = min) + ggplot2::geom_violin() + ggplot2::stat_summary(fun = median, geom = "point", shape = 20, size = 7, color = "red", fill = "red") + ggplot2::theme(axis.text = ggplot2::element_text(size = 14, face = "bold"), axis.title = ggplot2::element_text(size = 16, face = "bold"), plot.title = ggplot2::element_text(size = 20, face = "bold"), legend.position = "none") + ggplot2::ylab(yLabel) + ggplot2::coord_cartesian(ylim = c(0.3, 1)) + ggplot2::ggtitle(graphTitle)
	
	if (nchar(filename) > 0) ggplot2::ggsave(filename = filename, plot = vPlot_better, width = filedim[1] / 150, height = filedim[2] / 150, units = "in", dpi = 150)
	else print(vPlot_better)
	
	}

#----------------------------------------------------------------------

#' Adds two columns with data from Chromatinsight (loaded using getRegionData)
#' to a table containing connected regions from Hi-C experiments.
#' @export
#' @examples
#' dimorphismMono = getRegionData(folderPath, prefix = "mono", totRandomStates = 5)
#' HiCDataIncludingDimorphism = addDimorphism(data = HiCData, dimorphismData = dimorphismMono)
addDimorphism = function(data, dimorphismData, suffix = "", overwrite = FALSE, counterstep = 0, sourceColumn = "disparity") {
	
	endElements = c("bait", "oe")
	totCols = ncol(data)
	
	baitDimorphismCol = paste0("baitDimorphism", suffix)
	oeDimorphismCol = paste0("oeDimorphism", suffix)	
	
	if ((baitDimorphismCol %in% colnames(data) | oeDimorphismCol %in% colnames(data)) & !overwrite) {
		print(paste0("error: columns with suffix '", suffix, "' already exists, to overwrite, please use overwrite = TRUE'"))
		return(data)
		}
		
	data[, baitDimorphismCol] = -1
	data[, oeDimorphismCol] = -1
	
	counter = 0
	totRows = nrow(data)
	for(i in 1:nrow(data)) {
		counter = counter + 1
		if (counter %% 1000 == 0 & counterstep > 0) print(paste0("Calculating row #", counter, "/", totRows, " for ", suffix, "..."))
		for (thisPart in endElements) {
			thisDimorphismCol = paste0(thisPart, sourceColumn, suffix)
			thisChr = paste0(thisPart, "Chr")
			thisStart = paste0(thisPart, "Start")
			thisEnd = paste0(thisPart, "End")
			myChrom = paste0("chr", data[i, thisChr])
			myStart = data[i, thisStart]
			myEnd = data[i, thisEnd]
			
			segmentsIncludingElement = subset(dimorphismData,
									chrom == myChrom & (
									(init <= myStart & end >= myStart) | # it starts within the segment
									(init <= myEnd & end >= myEnd) | # or it finishes within the segment
									(init > myStart & end < myEnd) # or the HiCfragment contains the segment
									)
								)
			thisDimorphism = max(segmentsIncludingElement[, sourceColumn])
			data[i, thisDimorphismCol] = thisDimorphism
			
			#if (segmentList != 1 & segmentList != 2) {
			#	print(paste0("part: ", thisPart, ", chrom: ", myChrom, ", start = ", myStart, ", end = ", myEnd, " (", segmentList, ")"))
			#	}
			
			}
		}
	
	return(data)
	
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

#' Deprecated version of the function getHistModLevels that was using the file
#' name rather than the grouping file.
#' Adds five columns to the region table, including additional ChIP-seq data
#' for all windows and samples corresponding to a genomic region:
#' 1) Ratio of windows that have the histone modification.
#' 2) Average for sample group A.
#' 3) Average for sample group B.
#' 4) Relative ratio for sample group A.
#' 5) Relative ratio for sample group B.
#' @export
#' @examples
#' getHistModLevels2(data,
#'                   direc = "",
#'                   histMod = "ac",
#'                   prefixA = "mono*S*fem*NCMLS*",
#'                   nameA = "fem",
#'                   prefixB = "mono*S*mal*NCMLS*",
#'                   nameB = "mal",
#'                   suffix = "real",
#'                   windowSize = 200,
#'                   verbose = TRUE)
getHistModLevels2 = function(data,
							direc = "",
							histMod = "ac",
							prefixA = "",
							nameA = "1st",
							prefixB = "",
							nameB = "2nd",
							suffix = "real",
							windowSize = 200,
							verbose = TRUE) {
	
	if (histMod == "ac") histMod = "H3K27ac"
	if (histMod == "me" | histMod == "me1") histMod = "H3K4me1"
	if (histMod != "me" & histMod != "me1" & histMod != "ac" & histMod != "H3K27ac" & histMod != "H3K4me1") {
		print("Unknown histone modification, aborting.")
		exit()
		}
	
	myOutput = data
	myOutput[with(myOutput, order(chrom, init)),] # just in case, but they should be sorted already
	rownames(myOutput) = NULL
	
	col_modLevel = paste0("modLevel_", suffix)
	col_avgA = paste0(nameA, "_avg_", suffix)
	col_avgB = paste0(nameB, "_avg_", suffix)
	col_relA = paste0(nameA, "_rel_", suffix)
	col_relB = paste0(nameB, "_rel_", suffix)
	
	myOutput[col_modLevel] = NA
	myOutput[col_avgA] = NA
	myOutput[col_avgB] = NA
	myOutput[col_relA] = NA
	myOutput[col_relB] = NA
	currentChr = ""
	pileA = data.frame()
	pileB = data.frame()
	
	for (i in 1:nrow(myOutput)) {
		thisChr = as.character(myOutput[i, "chrom"])
		if (thisChr != currentChr) {
			currentChr = thisChr
			print(paste0("Working on ", currentChr))
			pileA = pileup2(prefix = prefixA, direc = direc, chrom = currentChr)
			pileB = pileup2(prefix = prefixB, direc = direc, chrom = currentChr)
			}
		thisStart = as.numeric(myOutput[i, "init"])
		thisEnd = as.numeric(myOutput[i, "end"])
		
		windowStart = floor(thisStart / windowSize)
		windowEnd = ceiling(thisEnd / windowSize)
		
		# A couple of sanity checks
		if (windowStart < 1) {
			windowStart = 1
			if (verbose) print("Warning: first window was 0, readjusted to 1.")
			}
		if (windowEnd > nrow(pileA)) {
			windowEndOld = windowEnd
			windowEnd = nrow(pileA)
			if (verbose) print(paste0("Warning: last window was ", windowEndOld, ", readjusted to ", windowEnd, "(the maximum for ", thisChr, ")."))
			}
		if (windowEnd < windowStart) {
			windowEndOld = windowEnd
			windowEnd = windowStart
			if (verbose) print(paste0("Warning: last window (", windewEndOld, ") was lower than first window (", windowEnd, "), readjusted to ", windowEnd, "."))
			}
		
		avgA = sum(pileA[windowStart:windowEnd, histMod]) / (windowEnd - windowStart + 1)
		avgB = sum(pileB[windowStart:windowEnd, histMod]) / (windowEnd - windowStart + 1)
		
		modLevel = (avgA + avgB) / 2
		
		relA = 0
		relB = 0
		
		if (avgA >= 0 & avgB >= 0 & !(avgA == 0 & avgB == 0)) {
			relA = avgA / (avgA + avgB)
			relB = avgB / (avgA + avgB)
			}
		
		myOutput[i, col_modLevel] = modLevel
		myOutput[i, col_avgA] = avgA
		myOutput[i, col_avgB] = avgB
		myOutput[i, col_relA] = relA
		myOutput[i, col_relB] = relB
		
		}
	
	return(myOutput)	
	}
	
#----------------------------------------------------------------------

#' Calculates the heatmap of a set of ChromHMM output binary files.
#' Parameters:
#'   grouping: the grouping file that indicates the files to be considered.
#'   key: the key in the grouping file to filter the files to be considered
#'          (omit to use all the files).
#'   coord: the coordinates of the region to extract, such as
#'          chrX:73,040,484-73,072,558 (this parameter is mandatory)
#'   histmod: the name of the column from which extract the data, such as the
#'          histone mark. Default is "H3K27ac".
#'   windowsize: the number of bases in each ChromHMM window. Default is 200.
#'   direc: directory or folder where the data are (omit to use the working
#'          directory).
#'   filename: the filename to save a 1800x1000 px image with the heatmap. If
#'           omitted, it will just show it in the plot area.
#' @export
#' @examples
#' getheatmap2(grouping = "grouping_file.txt",
#'            key = "fem",
#'            coord = "chrX:130,822,807-130,964,890",
#'            histmod = "H3K4me1",
#'            filename = "C:/data/Firre_reduced_me.png")
getheatmap = function(grouping, key = "", coord = "", histmod = "H3K27ac", windowsize = 200, filename = "", direc = "") {
	
	if (coord == "") {
		print("Genomic coordinates, such as chrX:73,040,484-73,072,558, are mandatory")
		exit()
		}
	
	myMatrix = getmatrix(grouping = grouping, key = key, histmod = histmod, coord = coord, windowsize = windowsize, direc = direc)
	
	if (filename != "") {
		dev.new()
		png(filename, width = 1800, height = 1000)
		heatmap(data.matrix(myMatrix), Colv = NA, scale = "none", col = c("#fff0c0", "#a04040"))
		graphics.off()
		}
	else heatmap(data.matrix(myMatrix), Colv = NA, scale = "none", col = c("#fff0c0", "#a04040"))
	
	}
	
#----------------------------------------------------------------------

#' Extracts the matrix of histone modifications to be used to
#' represent the heatmap.
getmatrix = function(grouping, key = "", coord = "", start = 1, end = -1, histmod = "H3K27ac", windowsize = 200, chromosome = "", direc = ""){
	
	if (tolower(histmod) == "ac") histmod = "H3K27ac"
	if (tolower(histmod) == "me" | tolower(histmod) == "me1") histmod = "H3K4me1"
	
	# note that coord overrides start, end and chromosome
	if (coord != "") {
		myCoord = splitCoords(coord)
		start = floor(myCoord$start / windowsize)
		end = ceiling(myCoord$end / windowsize)
		chromosome = myCoord$chrom
		}
	
	myDirec = direc
	if (nchar(myDirec) > 0) {
		if (substr(myDirec, nchar(myDirec), nchar(myDirec)) != "/") {
			myDirec = paste0(myDirec, "/")
			}
		}
	
	groupingFile = paste0(myDirec, grouping)
	groupingTable = as.data.frame(read.table(groupingFile, sep = "\t", header = FALSE))
	colnames(groupingTable) = c("filename", "tablekey")
	
	myFiles = groupingTable[,"filename"]
	myKeys = groupingTable[,"tablekey"]
	if (key != "") {
		myFiles = subset(groupingTable, tablekey == key)[,"filename"]
		myKeys = subset(groupingTable, tablekey == key)[,"tablekey"]
		}
	
	print(paste0("Checking ", length(myFiles), " files."))
	
	myMatrix = data.frame()
	newKeys = data.frame()
	
	total = 0
	for (i in 1:length(myFiles)){
		myFile = myFiles[i]
		myFileChr = gsub("\\*", chromosome, myFile)
		print(paste0("Loading: ", myFileChr))
		myTable = read.table(myFileChr, sep = "\t", header = TRUE, skip = 1)
		if (end == -1) thisEnd = nrow(myTable)
		else thisEnd = end
		myRow = t(myTable[start:thisEnd, histmod])
		rownames(myRow) = paste0(myKeys[i], ":", myFileChr)
		columnNames = format((start:thisEnd) * windowsize, big.mark = ",")
		if (chromosome != "") columnNames = paste0(chromosome, ":", columnNames)
		colnames(myRow) = columnNames
		
		goodNumbers = ((min(myRow) == 0 | min(myRow) == 1) & (max(myRow) == 0 | max(myRow) == 1))
		
		if (goodNumbers) {
			myMatrix = rbind(myMatrix, myRow)
			newKeys = rbind(newKeys, myKeys[i])
			}
		else print(paste0("Warning: File ", myFileChr, " excluded [min: ", min(myRow), ", max: ", max(myRow), "]."))
		}
	
	# myKeys is not being returned for now
	return(myMatrix)
	
	}
	
#----------------------------------------------------------------------

#' Deprecated version of getheatmap that uses the pattern of the filename
#' rather than the grouping file.
#' Calculates the heatmap of a set of ChromHMM output binary files.
#' Parameters:
#'   pattern: to select the files from the filename. Asterisks (*) as wildcard.
#'   coord: the coordinates of the region to extract, such as
#'          chrX:73,040,484-73,072,558 (this parameter is mandatory)
#'   histmod: the name of the column from which extract the data, such as the
#'          histone mark. Default is "H3K27ac".
#'   windowsize: the number of bases in each ChromHMM window. Default is 200.
#'   direc: directory or folder where the data are (omit to use the working
#'          directory).
#'   filename: the filename to save a 1800x1000 px image with the heatmap. If
#'           omitted, it will just show it in the plot area.
#' @export
#' @examples
#' getheatmap2(pattern = "mono*",
#'            coord = "chrX:130,822,807-130,964,890",
#'            histmod = "H3K4me1",
#'            filename = "C:/data/Firre_reduced_me.png")
getheatmap2 = function(pattern = "*", coord = "", histmod = "H3K27ac", windowsize = 200, filename = "", direc = "") {
	
	if (coord == "") {
		print("Genomic coordinates, such as chrX:73,040,484-73,072,558, are mandatory")
		exit()
		}
	
	myMatrix = getmatrix2(direc = direc, pattern = pattern, histmod = histmod, coord = coord, windowsize = windowsize)
	
	if (filename != "") {
		dev.new()
		png(filename, width = 1800, height = 1000)
		heatmap(data.matrix(myMatrix), Colv = NA, scale = "none", col = c("#fff0c0", "#a04040"))
		graphics.off()
		}
	else heatmap(data.matrix(myMatrix), Colv = NA, scale = "none", col = c("#fff0c0", "#a04040"))
	
	}
	
#----------------------------------------------------------------------

#' Deprecated version of getmatrix that uses the pattern of the filename
#' rather than the grouping file.
#' Extracts the matrix of histone modifications to be used to
#' represent the heatmap.
getmatrix2 = function(direc = "", pattern = "*", coord = "", start = 1, end = -1, histmod = "H3K27ac", windowsize = 200, chromosome = ""){
	
	if (tolower(histmod) == "ac") histmod = "H3K27ac"
	if (tolower(histmod) == "me" | tolower(histmod) == "me1") histmod = "H3K4me1"
	
	# note that coord overrides start, end and chromosome
	if (coord != "") {
		myCoord = splitCoords(coord)
		start = floor(myCoord$start / windowsize)
		end = ceiling(myCoord$end / windowsize)
		chromosome = myCoord$chrom
		}
	
	myDirec = direc
	if (nchar(myDirec) > 0) {
		if (substr(myDirec, nchar(myDirec), nchar(myDirec)) != "/") {
			myDirec = paste0(myDirec, "/")
			}
		}
	myFilePath = paste0(myDirec, pattern) # the pattern of the *_binary.txt files
	print(myFilePath)
	
	myFiles = Sys.glob(myFilePath)
	print(paste0("Checking ", length(myFiles), " files."))
	
	myMatrix = data.frame()
	
	total = 0
	for (myFile in myFiles){
		print(paste0("Loading: ", myFile))
		myTable = read.table(myFile, sep = "\t", header = TRUE, skip = 1)
		if (end == -1) thisEnd = nrow(myTable)
		else thisEnd = end
		myRow = t(myTable[start:thisEnd, histmod])
		
		rownames(myRow) = myFile
		columnNames = format((start:thisEnd) * windowsize, big.mark = ",")
		if (chromosome != "") columnNames = paste0(chromosome, ":", columnNames)
		colnames(myRow) = columnNames
		
		goodNumbers = ((min(myRow) == 0 | min(myRow) == 1) & (max(myRow) == 0 | max(myRow) == 1))
		
		if (goodNumbers) myMatrix = rbind(myMatrix, myRow)
		else print(paste0("Warning: File ", myFile, " excluded [min: ", min(myRow), ", max: ", max(myRow), "]."))
		}
		
	return(myMatrix)
	
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