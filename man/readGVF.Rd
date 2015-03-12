\name{readGVF}
\alias{readGVF}
\title{
Read CNVs from a .gvf file
}
\description{
Reads CNV data from a .gvf file.
}
\usage{
readGVF(filename)
}
\arguments{
  \item{filename}{The name of a Genome Variation Format file.}
}
\value{
	A data frame containing the CNVs. Each row contains the following columns for one CNV:

	\describe{
		\item{\code{SampleID}}{ID assigned to the subject's DNA sample in which the CNV was
	found}
		\item{\code{Chr}}{Chromosome on which the CNV is located}
		\item{\code{Coord_i}}{Start position of the CNV on the chromosome}
		\item{\code{Coord_f}}{End position of the CNV on the chromosome}
		\item{\code{Type}}{CNV type (e.g. "DEL" or "DUP")}
		\item{\code{Genes}}{(empty string)}
		\item{\code{CnvID}}{ID assigned to the CNV}
	}

	N.B. that the \code{Genes} column is *not* assigned. Normally it should contain genes affected by the CNV -- stored in a delimited format inside a character string (e.g. "54777;255352;84435" for semicolon-delimited EntrezGene identifiers). To assign this column, use the getCnvGenes() function. The data frame can then be assigned to the \code{cnv} element of the \code{cnvData} list structure that comprises one of the slots in a \code{CnvGSAInput} object.
}
\references{
http://www.sequenceontology.org/resources/gvf.html
}
\author{Robert Ziman \email{rziman@gmail.com}}
\examples{
library("cnvGSAdata")
cnvFile <- system.file("extdata", "cnv.gvf", package="cnvGSAdata")
cnv <- readGVF(cnvFile)
}
