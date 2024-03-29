% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools.R
\name{drawpile}
\alias{drawpile}
\title{You give it the frequency of a histone mark at every ChromHMM bin,
ie, the output of two pileup functions (see pileup).
and draws the graph comparing them.
Include a filename parameter if you want to save the output into a file.
For H3K27ac the first group is red, the second is blue,
For H3K4me1 the first group is pink, the second is forest green,
To use UCSC Genome Browser coordinate format (such as chrX:10,000-25,000)
see drawpilegb.}
\usage{
drawpile(
  data1,
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
)
}
\description{
You give it the frequency of a histone mark at every ChromHMM bin,
ie, the output of two pileup functions (see pileup).
and draws the graph comparing them.
Include a filename parameter if you want to save the output into a file.
For H3K27ac the first group is red, the second is blue,
For H3K4me1 the first group is pink, the second is forest green,
To use UCSC Genome Browser coordinate format (such as chrX:10,000-25,000)
see drawpilegb.
}
