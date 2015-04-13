\name{CnvGSAOutput-class}
\Rdversion{1.1}
\docType{class}

\alias{CnvGSAOutput}
\alias{CnvGSAOutput-class}
%\alias{phData.ls}
\alias{phData.ls,CnvGSAOutput-method}
\alias{res.ls}
\alias{res.ls,CnvGSAOutput-method}
\alias{gsTables.ls}
\alias{gsTables.ls,CnvGSAOutput-method}
%\alias{gsData.ls}
\alias{gsData.ls,CnvGSAOutput-method}
\alias{config.df}
\alias{config.df,CnvGSAOutput-method}

\title{Class \code{"CnvGSAOutput"}}
\description{
Container class for the output structures produced by the main function (i.e. cnvGSAlogRegTest()).
}
\section{Slots}{
  \describe{
    \item{\code{phData.ls}:}{Object of class \code{"list"} containing phenotype and covariate data}
    \item{\code{res.ls}:}{Object of class \code{"list"} containing burden analysis results for gene-sets}
    \item{\code{gsTables.ls}:}{Object of class \code{"list"} containing the gene-set tables}
    \item{\code{gsData.ls}:}{Object of class \code{"list"} containing gene-set data}
    \item{\code{config.df}:}{Object of class \code{"list"} containing config data frame}
  }
}
\section{Methods}{
  \describe{
    \item{gsTables.ls}{\code{signature(obj = "CnvGSAOutput")}: Gets \code{gsTables.ls}.}
    \item{res.ls}{\code{signature(obj = "CnvGSAOutput")}: Gets \code{res.ls}.}
    \item{phData.ls}{\code{signature(obj = "CnvGSAOutput")}: Gets \code{phData.ls}.}
    \item{config.df}{\code{signature(obj = "CnvGSAOutput")}: Gets \code{config.df}.}
    \item{gsData.ls}{\code{signature(obj = "CnvGSAOutput")}: Gets \code{gsData.ls}.}
  }
}
\author{Joseph Lugo \email{joseph.r.lugo@gmail.com}}
\examples{
## See vignette for full details and worked example
}
\keyword{classes}
