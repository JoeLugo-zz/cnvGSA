\name{readParamsRFile}
\alias{readParamsRFile}
\title{Read cnvGSA parameters from a file (R format)}
\description{
Reads cnvGSA parameters from a text file that has them encoded using R syntax.
}
\details{
To make it easier to integrate the association test into a larger bioinformatics pipeline, it is
convenient to read in the parameters from an external source such as a text file. One such
implementation is to record each parameter on its own line using R syntax:

\code{
grandtotals_mode <- "all"
sample_classes <- c("case", "ctrl")
fdr_iter <- 1000
extended_report <- 200
boxplot_PDFs <- FALSE
limits_type <- "DEL"
}

The package provides \code{readParamsRFile()} to parse such a file (essentially just source()ing it
and then handling the few possibilites around the cnvData\$filters parameters).
}
\usage{
readParamsRFile(filename)
}
\arguments{
  \item{filename}{Name of the file.}
}
\value{
Test parameters in a list structure that can be assigned to the \code{params} slot of a \code{CnvGSAInput} object.
}
\author{Robert Ziman \email{rziman@gmail.com}}
\examples{
paramFile <- system.file("scripts", "params_example.R", package="cnvGSA")
params <- readParamsRFile(paramFile)
}
