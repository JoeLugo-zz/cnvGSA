\name{readGMT}
\alias{readGMT}
\title{
Read gene-sets from a .gmt file
}
\description{
Reads gene-sets from a .gmt file.
}
\usage{
readGMT(filename)
}
\arguments{
  \item{filename}{The name of a file in the Gene Matrix Transposed format.}
}
\value{
A list having the following two elements:

gs2gene - A list of character vectors where each vector contains the genes for a
particular gene-set; the gene-set names ("GO:0030850" etc.) are stored as the names
of the list elements. Since gene-sets can hold different numbers of genes, the vectors
will typically have different lengths.

gs2name - A single character vector mapping each gene-set name to its description.
The descriptions are stored as the vector elements and the gene-set names are stored
as the names of the vector elements.

This list structure should be assigned to the gsData slot of a CnvGSAInput object.
}
\references{
GSEA wiki: Data formats
http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29
}
\author{Robert Ziman \email{rziman@gmail.com}}
\examples{
library("cnvGSAdata")
gsDataFile <- system.file("extdata", "gsData.gmt", package="cnvGSAdata")
gsData <- readGMT(gsDataFile)
}
