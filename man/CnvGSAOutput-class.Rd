\name{CnvGSAOutput-class}
\Rdversion{1.1}
\docType{class}

\alias{CnvGSAOutput}
\alias{CnvGSAOutput-class}
%\alias{cnvData}
\alias{cnvData,CnvGSAOutput-method}
\alias{burdenSample}
\alias{burdenSample,CnvGSAOutput-method}
\alias{burdenGs}
\alias{burdenGs,CnvGSAOutput-method}
%\alias{geneData}
\alias{geneData,CnvGSAOutput-method}
\alias{enrRes}
\alias{enrRes,CnvGSAOutput-method}
\title{Class \code{"CnvGSAOutput"}}

\description{
Container class for the output structures produced by the main function (i.e. cnvGSA.Fisher()).
}
\section{Slots}{
  \describe{
    \item{\code{cnvData}:}{Object of class \code{"list"} containing original and filtered CNV data}
    \item{\code{burdenSample}:}{Object of class \code{"list"} containing burden analysis results for objects}
    \item{\code{burdenGs}:}{Object of class \code{"list"} containing burden analysis results for gene-sets}
    \item{\code{geneData}:}{Object of class \code{"list"} containing gene-centric statistics}
    \item{\code{enrRes}:}{Object of class \code{"list"} containing the gene-set enrichment results}
  }
}
\section{Methods}{
  \describe{
    \item{burdenGs}{\code{signature(obj = "CnvGSAOutput")}: Gets \code{burdenGs}.}
    \item{burdenSample}{\code{signature(obj = "CnvGSAOutput")}: Gets \code{burdenSample}.}
    \item{cnvData}{\code{signature(obj = "CnvGSAOutput")}: Gets \code{cnvData}.}
    \item{enrRes}{\code{signature(obj = "CnvGSAOutput")}: Gets \code{enrRes}.}
    \item{geneData}{\code{signature(obj = "CnvGSAOutput")}: Gets \code{geneData}.}
    \item{summary}{\code{signature(object = "CnvGSAOutput")}: Prints out several summary statistics.}
  }
}
\author{Robert Ziman \email{rziman@gmail.com}}
\examples{
## See vignette for full discussion of output elements
library("cnvGSAdata")
data("cnvGSA_output_example")
slotNames("output")
}
\keyword{classes}
