% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools.R
\name{drawpilegb}
\alias{drawpilegb}
\title{Calls drawpile with UCSC Genome Browser coordinate format.
(such as chrX:10,000-25,000)
(See drawpile and densitypile for more details).}
\usage{
drawpilegb(
  datafem,
  datamal,
  coord = "",
  plusminus = 1000,
  sex = "both",
  histmod = "ac",
  window = 200,
  useDensity = TRUE,
  opacity = 0.5,
  highlight = "",
  label = ""
)
}
\description{
Calls drawpile with UCSC Genome Browser coordinate format.
(such as chrX:10,000-25,000)
(See drawpile and densitypile for more details).
}
