\name{cnvGSAexportBurdenStats}
\alias{cnvGSAexportBurdenStats}
\title{
Export burden statistics
}
\description{
Exports burden statistics into a tab-delimited text file.
}
\usage{
cnvGSAexportBurdenStats( output, filenamePrefix )
}
\arguments{
	\item{output}{Object of class \link{CnvGSAOutput}}
	\item{filenamePrefix}{Prefix to use in the names of the output files (see below)}
}
\value{
Given filenamePrefix="example", outputs "example_burdenGs_coverage.txt" and "example_burdenGs_pairs.txt".
}
\author{Robert Ziman \email{rziman@gmail.com}}
\examples{
library( "cnvGSAdata" )
data( "cnvGSA_output_example" )
cnvGSAexportBurdenStats( output, "example" )
}
