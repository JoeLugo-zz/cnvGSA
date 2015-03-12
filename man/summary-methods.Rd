\name{summary-methods}
\docType{methods}
\alias{summary-methods}
\alias{summary,CnvGSAOutput-method}
\title{cnvGSA methods for function \code{summary}}
\description{
cnvGSA methods for function \code{summary}
}
\section{Methods}{
	\describe{
		\item{\code{signature(object = "CnvGSAOutput")}}{
			Outputs several summary statistics from the CnvGSAOutput object.
		}
	}
}
\examples{
library("cnvGSAdata")
data("cnvGSA_output_example")
summary(output)
}
