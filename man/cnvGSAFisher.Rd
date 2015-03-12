\name{cnvGSAFisher}
\alias{cnvGSAFisher}
\title{
CNV/gene-set association using Fisher Exact Test
}
\description{
Performs the CNV/gene-set association using the Fisher Exact Test.
}
\details{
Whenever a subject has at least one gene in a gene-set affected by a rare variant, a perturbation score of 1 is assigned to the (subject, gene-set) pair; for each gene-set, subject counts are tested vs. control counts using the Fisher Exact Test (FET). Significant gene-sets will have a significantly high count in cases compared to controls. Statistical reports on burden are also generated.
}
\usage{
cnvGSAFisher( input )
}
\arguments{
	\item{input}{Object of class CnvGSAInput.}
}
\value{
Outputs an object of class CnvGSAOutput containing detailed enrichment results as well as gene-centric and global burden statistics.
}
\author{Robert Ziman \email{rziman@gmail.com}}
\examples{
library("cnvGSAdata")
data("cnvGSA_input_example")
output <- cnvGSAFisher(input)
}
