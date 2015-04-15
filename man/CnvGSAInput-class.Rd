\name{CnvGSAInput-class}
\Rdversion{1.1}
\docType{class}

\alias{CnvGSAInput-class}
\alias{CnvGSAInput}
\alias{config.ls}
\alias{config.ls,CnvGSAInput-method}
\alias{config.ls<-}
\alias{config.ls<-,CnvGSAInput-method}
\alias{params.ls}
\alias{params.ls,CnvGSAInput-method}
\alias{params.ls<-}
\alias{params.ls<-,CnvGSAInput-method}
\alias{cnvData.ls}
\alias{cnvData.ls,CnvGSAInput-method}
\alias{cnvData.ls<-}
\alias{cnvData.ls<-,CnvGSAInput-method}
\alias{phData.ls}
\alias{phData.ls,CnvGSAInput-method}
\alias{phData.ls<-}
\alias{phData.ls<-,CnvGSAInput-method}
\alias{gsData.ls}
\alias{gsData.ls,CnvGSAInput-method}
\alias{gsData.ls<-}
\alias{gsData.ls<-,CnvGSAInput-method}
\alias{geneID.ls}
\alias{geneID.ls,CnvGSAInput-method}
\alias{geneID.ls<-}
\alias{geneID.ls<-,CnvGSAInput-method}

\title{Class \code{"CnvGSAInput"}}
\description{
Container class for the input structures required by the main function (i.e. cnvGSALogRegTest()).
}
\section{Slots}{
  \describe{
    \item{\code{config.ls}:}{Object of class \code{"list"} containing the file names and paths.}
    \item{\code{params.ls}:}{Object of class \code{"list"} containing the test parameters.}
    \item{\code{cnvData.ls}:}{Object of class \code{"list"} containing CNV data.}
    \item{\code{phData.ls}:}{Object of class \code{"list"} containing phenotype and covariate data.}
    \item{\code{gsData.ls}:}{Object of class \code{"list"} containing gene-set data.}
    \item{\code{geneID.ls}:}{Object of class \code{"list"} containing Gene ID data.}
  }
}
\section{Constructor}{
      \code{CnvGSAInput( config.ls, params.ls, cnvData.ls, phData.ls, gsData.ls, geneID.ls)}:
      Creates a CnvGSAInput object.
      \describe{
        \item{\code{config.ls}}{Structure containing the file names and paths.}
        \item{\code{params.ls}}{Structure containing main test parameters.}
        \item{\code{cnvData.ls}}{Structure containing CNV data along with sample-to-class information and filters.}
        \item{\code{phData.ls}}{Structure containing phenotype and covariate data.}
        \item{\code{gsData.ls}}{Structure containing gene-set data.}
        \item{\code{geneID.ls}}{Structure containing Gene ID data.}
      }
      See the vignette for complete details on each of these structures as well as a full example of how to load them.
}

\section{Methods}{
  \describe{
    \item{config.ls}{\code{signature(obj = "CnvGSAInput")}: Gets \code{config.ls}. }
    \item{config.ls<-}{\code{signature(obj = "CnvGSAInput")}: Sets \code{config.ls}. }
    \item{params.ls}{\code{signature(obj = "CnvGSAInput")}: Gets \code{params.ls}. }
    \item{params.ls<-}{\code{signature(obj = "CnvGSAInput")}: Sets \code{params.ls}. }
    \item{cnvData.ls}{\code{signature(obj = "CnvGSAInput")}: Gets \code{cnvData.ls}.}
    \item{cnvData.ls<-}{\code{signature(obj = "CnvGSAInput")}: Sets \code{cnvData.ls}.}
    \item{phData.ls}{\code{signature(obj = "CnvGSAInput")}: Gets \code{phData.ls}.}
    \item{phData.ls<-}{\code{signature(obj = "CnvGSAInput")}: Sets \code{phData.ls}.}
    \item{gsData.ls}{\code{signature(obj = "CnvGSAInput")}: Gets \code{gsData.ls}.}
    \item{gsData.ls<-}{\code{signature(obj = "CnvGSAInput")}: Sets \code{gsData.ls}.}
    \item{geneID.ls}{\code{signature(obj = "CnvGSAInput")}: Gets \code{geneID.ls}.}
    \item{geneID.ls<-}{\code{signature(obj = "CnvGSAInput")}: Sets \code{geneID.ls}.}
  }
}
\author{Joseph Lugo \email{joseph.r.lugo@gmail.com}}
\examples{
## See vignette for full details and worked example
}
\keyword{classes}
