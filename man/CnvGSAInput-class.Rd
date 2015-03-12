\name{CnvGSAInput-class}
\Rdversion{1.1}
\docType{class}

\alias{CnvGSAInput-class}
\alias{CnvGSAInput}
\alias{cnvData}
\alias{cnvData,CnvGSAInput-method}
\alias{cnvData<-}
\alias{cnvData<-,CnvGSAInput-method}
\alias{gsData}
\alias{gsData,CnvGSAInput-method}
\alias{gsData<-}
\alias{gsData<-,CnvGSAInput-method}
\alias{geneData}
\alias{geneData,CnvGSAInput-method}
\alias{geneData<-}
\alias{geneData<-,CnvGSAInput-method}
\alias{params}
\alias{params,CnvGSAInput-method}
\alias{params<-}
\alias{params<-,CnvGSAInput-method}

\title{Class \code{"CnvGSAInput"}}
\description{
Container class for the input structures required by the main function (i.e. cnvGSA.Fisher()).
}
\section{Slots}{
  \describe{
    \item{\code{cnvData}:}{Object of class \code{"list"} containing CNV data}
    \item{\code{gsData}:}{Object of class \code{"list"} containing gene-set data}
    \item{\code{geneData}:}{Object of class \code{"list"} containing gene annotations (symbols and descriptive names)}
    \item{\code{params}:}{Object of class \code{"list"} containing the test parameters}
  }
}
\section{Constructor}{
      \code{CnvGSAInput( cnvData, gsData, geneData, params )}:
      Creates a CnvGSAInput object.
      \describe{
        \item{\code{cnvData}}{Structure containing CNV data along with sample-to-class information and filters}
        \item{\code{gsData}}{Structure containing gene-set data}
        \item{\code{geneData}}{Structure containing gene annotations (symbols and descriptive names)}
        \item{\code{params}}{Structure containing main test parameters}
      }
      See the vignette for complete details on each of these structures as well as a full example of how to load them.
}

\section{Methods}{
  \describe{
    \item{cnvData}{\code{signature(obj = "CnvGSAInput")}: Gets \code{cnvData}.}
    \item{cnvData<-}{\code{signature(obj = "CnvGSAInput")}: Sets \code{cnvData}.}
    \item{geneData}{\code{signature(obj = "CnvGSAInput")}: Gets \code{geneData}.}
    \item{geneData<-}{\code{signature(obj = "CnvGSAInput")}: Sets \code{geneData}.}
    \item{gsData}{\code{signature(obj = "CnvGSAInput")}: Gets \code{gsData}.}
    \item{gsData<-}{\code{signature(obj = "CnvGSAInput")}: Sets \code{gsData}.}
    \item{params}{\code{signature(obj = "CnvGSAInput")}: Gets \code{params}. }
    \item{params<-}{\code{signature(obj = "CnvGSAInput")}: Sets \code{params}. }
  }
}
\author{Robert Ziman \email{rziman@gmail.com}}
\examples{
## See vignette for full details and worked example
}
\keyword{classes}
