% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools.R
\name{drawpile}
\alias{drawpile}
\title{You give it the frequency of a histone mark at every ChromHMM bin,
ie, the output of two pileup functions (see pileup).
and draws the graph comparing them.
One group is blue, the other is red, when they overlap it's purple.
To use UCSC Genome Browser coordinate system (such as chrX:10,000-25,000)
see drawpilegb.}
\usage{
drawpile(
  datafem,
  datamal,
  start = 0,
  plus = 1000,
  sex = "both",
  histmod = "ac",
  chrom = ""
)
}
\description{
You give it the frequency of a histone mark at every ChromHMM bin,
ie, the output of two pileup functions (see pileup).
and draws the graph comparing them.
One group is blue, the other is red, when they overlap it's purple.
To use UCSC Genome Browser coordinate system (such as chrX:10,000-25,000)
see drawpilegb.
}