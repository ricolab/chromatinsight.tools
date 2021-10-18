##################################
### Chromatinsight tools v1.13 ###
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
drawpilegb = function(datafem,
						datamal,
						coord = "",
						plusminus = 1000,
						sex = "both",
						histmod = "ac",
						window = 200,
						useDensity = TRUE,
						opacity = 0.5,
						highlight = "",
						label = "",
						filename = "",
						filedim = c(1800, 1000)
						) {

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
	if (useDensity) densitypile(datafem, datamal, start = myStart, plus = myPlus, sex = sex, histmod = histmod, chrom = chrom, opacity = opacity, highstart = highlightStart, highplus = highlightPlus, label = label, filename = filename, filedim = filedim)
    else drawpile(datafem, datamal, start = myStart, plus = myPlus, sex = sex, histmod = histmod, chrom = chrom, filename = filename, filedim = filedim)
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
#' Include a filename parametre if you want to save the output into a file.
#' For H3K27ac the first group is red, the second is blue,
#' For H3K4me1 the first group is pink, the second is forest green,
#' To use UCSC Genome Browser coordinate format (such as chrX:10,000-25,000)
#' see drawpilegb.
#' @export
drawpile = function(datafem,
					datamal,
					start = 0,
					plus = 1000,
					sex = "both",
					histmod = "ac",
					chrom = "",
					filename = "",
					filedim = c(1800, 1000)
					){
	
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
	
	myPlot <- ggplot2::ggplot(ggplot2::aes(x = as.numeric(row.names(datafempart)), y = 0), data = datafempart) + ggplot2::geom_line() + ggplot2::ylim(0, 1) + ggplot2::scale_colour_manual(values = c("males, H3K27ac" = "blue", "females, H3K27ac" = "red",  "females, H3K4me1" = "purple", "males, H3K2me1" = "forest green"))
	myPlot <- myPlot + ggplot2::labs(title = myTitle) + ggplot2::xlab("bin in ChromHMM (1 bin = 200b)") + ggplot2::ylab("probability of histone modification")
	
	if (histmod == "ac" | histmod == "both"){
		if (sex == "fem" | sex == "both") myPlot <- myPlot + ggplot2::geom_line(ggplot2::aes(y = H3K27ac, color = 'females, H3K27ac'), data = datafempart, size = lineWidth)
		if (sex == "mal" | sex == "both") myPlot <- myPlot + ggplot2::geom_line(ggplot2::aes(y = H3K27ac, color = 'males, H3K27ac'), data = datamalpart, size = lineWidth)
		}
	
	if (histmod == "me1" | histmod == "both"){
		if (sex == "fem" | sex == "both") myPlot <- myPlot + ggplot2::geom_line(ggplot2::aes(y = H3K4me1, color = 'females, H3K4me1'), data = datafempart, size = lineWidth)
		if (sex == "mal" | sex == "both") myPlot <- myPlot + ggplot2::geom_line(ggplot2::aes(y = H3K4me1, color = 'males, H3K2me1'), data = datamalpart, size = lineWidth)
		}
		
	if (filename != ""){
		ggplot2::ggsave(filename, plot = myPlot, width = filedim[1] / 300, height = filedim[2] / 300, units = "in", dpi = 300)
		}
	else print(myPlot)
	
	}
    
#----------------------------------------------------------------------

#' You give it the output of two pileup function
#' and it shows the graph comparing them.
#' Include a filename parametre if you want to save the output into a file.
#' For H3K27ac the first group is red, the second is blue,
#' when they overlap it's purple.
#' For H3K4me1 the first group is pink, the second is forest green,
#' when they overlap it's grey.
#' To use UCSC Genome Browser coordinate format (such as chrX:10,000-25,000)
#' see drawpilegb.
#' @export
densitypile = function(datafem, datamal,
                    start = 0, plus = 1000,
                    highstart = 0, highplus = 0,
                    sex = "both",
                    histmod = "ac",
                    chrom = "",
                    opacity = 0.5,
                    label = "",
					filename = "",
					filedim = c(1800,1000)
					){
	
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
	
	myPlot <- ggplot2::ggplot(ggplot2::aes(x = as.numeric(row.names(datafempart)) * binSize, y = 0), data = datafempart) + ggplot2::geom_line() + ggplot2::ylim(0, 1) + ggplot2::scale_fill_manual(values = c("males, H3K27ac" = "#0000ff", "females, H3K27ac" = "#ff0000",  "females, H3K4me1" = "purple", "males, H3K2me1" = "forest green"))
	myPlot <- myPlot + ggplot2::labs(title = myTitle) + ggplot2::xlab("bins in ChromHMM (1 bin = 200b)") + ggplot2::ylab("probability of histone modification")
    myPlot <- myPlot + ggplot2::theme(panel.background = ggplot2::element_rect(fill="white", color = "grey50", size=2), panel.grid.major = ggplot2::element_line(color = "grey",size=(0.2)))
    
    if (highstart > 0 & highplus > 0) {
        myPlot <- myPlot + ggplot2::geom_rect(ggplot2::aes(xmin = highstart * binSize, xmax = (highstart + highplus) * binSize, ymin = 0, ymax = 1), fill = "#ffff00", alpha = 0.005)
        }
	
	if (histmod == "ac" | histmod == "H3K27ac" | histmod == "both"){
		if (sex == "fem" | sex == "both") myPlot <- myPlot + ggplot2::geom_density(ggplot2::aes(y = H3K27ac, fill = 'females, H3K27ac'), data = datafempart, size = lineWidth, stat = 'identity', alpha = opacity)
		if (sex == "mal" | sex == "both") myPlot <- myPlot + ggplot2::geom_density(ggplot2::aes(y = H3K27ac, fill = 'males, H3K27ac'), data = datamalpart, size = lineWidth, stat = 'identity', alpha = opacity)
		}
	
	if (histmod == "me1" | histmod == "me" | histmod == "H3K4me1" | histmod == "both"){
		if (sex == "fem" | sex == "both") myPlot <- myPlot + ggplot2::geom_density(ggplot2::aes(y = H3K4me1, fill = 'females, H3K4me1'), data = datafempart, size = lineWidth, stat = 'identity', alpha = opacity)
		if (sex == "mal" | sex == "both") myPlot <- myPlot + ggplot2::geom_density(ggplot2::aes(y = H3K4me1, fill = 'males, H3K2me1'), data = datamalpart, size = lineWidth, stat = 'identity', alpha = opacity)
		}
    
    if (nchar(label) > 0) {
        myPlot <- myPlot + ggplot2::annotate("text", x = start * binSize, y = 0.9, label = label, family = "Ubuntu", size = 10, hjust = 0)
        }
    
	if (filename != ""){
		ggplot2::ggsave(filename, plot = myPlot, width = filedim[1] / 300, height = filedim[2] / 300, units = "in", dpi = 300)
		}
	else print(myPlot)
	
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

#' Displays the violin plots of the distributions in two analyses
#' retrieved using getRegionData. Especially useful to compare random vs
#' observed data.
#' @export
compareDistributions = function(regionData1, regionData2, label1 = "first", label2 = "second", yLabel = "data", graphTitle = "comparison") {

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
	
	print(vPlot_better)
	
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

#' Adds five columns to the region table, including additional ChIP-seq data
#' for all windows and samples corresponding to a genomic region:
#' 1) Ratio of windows that have the histone modification.
#' 2) Average for sample group A.
#' 3) Average for sample group B.
#' 4) Relative ratio for sample group A.
#' 5) Relative ratio for sample group B.
#' @export
#' @examples
#' getHistModLevels(data,
#'                   direc = "",
#'                   histMod = "ac",
#'                   prefixA = "mono*S*fem*NCMLS*",
#'                   nameA = "fem",
#'                   prefixB = "mono*S*mal*NCMLS*",
#'                   nameB = "mal",
#'                   suffix = "real",
#'                   windowSize = 200,
#'                   verbose = TRUE)
getHistModLevels = function(data,
							direc = "",
							histMod = "ac",
							prefixA = "",
							nameA = "fem",
							prefixB = "",
							nameB = "mal",
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
			pileA = pileup(prefix = prefixA, direc = direc, chrom = currentChr)
			pileB = pileup(prefix = prefixB, direc = direc, chrom = currentChr)
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