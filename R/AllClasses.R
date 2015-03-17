## --------------
## CnvGSAInput S4 
## --------------

setClass( "CnvGSAInput",
	representation(
		config.ls = "list",
		params.ls = "list",
		cnvData.ls = "list",
		phData.ls = "list",
		gsData.ls = "list"
	)
)

# CnvGSAInput constructor
CnvGSAInput <- function(
					config.ls = list(),
					params.ls = list(),
					cnvData.ls = list(),
					phData.ls = list(),
					gsData.ls = list()
				)
{
    new( "CnvGSAInput", config.ls = config.ls, params.ls = params.ls, cnvData.ls = cnvData.ls, phData.ls = phData.ls, gsData.ls = gsData.ls )
}

# CnvGSAInput accessors
setGeneric( "config.ls", function(obj) standardGeneric("config.ls") )
setGeneric( "config.ls<-", function(obj, value) standardGeneric("config.ls<-") )
setMethod( "config.ls", "CnvGSAInput", function(obj){ obj@config.ls } )
setReplaceMethod( "config.ls", "CnvGSAInput", function(obj, value){ obj@config.ls <- value } )

setGeneric( "params.ls", function(obj) standardGeneric("params.ls") )
setGeneric( "params.ls<-", function(obj, value) standardGeneric("params.ls<-") )
setMethod( "params.ls", "CnvGSAInput", function(obj){ obj@params.ls } )
setReplaceMethod( "params.ls", "CnvGSAInput", function(obj, value){ obj@params.ls <- value } )

setGeneric( "cnvData.ls", function(obj) standardGeneric("cnvData.ls") )
setGeneric( "cnvData.ls<-", function(obj, value) standardGeneric("cnvData.ls<-") )
setMethod( "cnvData.ls", "CnvGSAInput", function(obj){ obj@cnvData.ls } )
setReplaceMethod( "cnvData.ls", "CnvGSAInput", function(obj, value){ obj@cnvData.ls <- value } )

setGeneric( "phData.ls", function(obj) standardGeneric("phData.ls") )
setGeneric( "phData.ls<-", function(obj, value) standardGeneric("phData.ls<-") )
setMethod( "phData.ls", "CnvGSAInput", function(obj){ obj@phData.ls } )
setReplaceMethod( "phData.ls", "CnvGSAInput", function(obj, value){ obj@phData.ls <- value } )

setGeneric( "gsData.ls", function(obj) standardGeneric("gsData.ls") )
setGeneric( "gsData.ls<-", function(obj, value) standardGeneric("gsData.ls<-") )
setMethod( "gsData.ls", "CnvGSAInput", function(obj){ obj@gsData.ls } )
setReplaceMethod( "gsData.ls", "CnvGSAInput", function(obj, value){ obj@gsData.ls <- value } )

## ---------------
## CnvGSAOutput S4
## ---------------

setClass( "CnvGSAOutput",
	representation(
		res.ls = "list",
		gsTables.ls = "list",
		gsData.ls = "list",
		phData.ls = "list",
		config.df = "list"
	)
)

# CnvGSAOutput constructor
CnvGSAOutput <- function(
					res.ls = list(),
					gsTables.ls = list(),
					gsData.ls = list(),
					phData.ls = list(),
					config.df = list()
				)
{
    new( "CnvGSAOutput", res.ls = res.ls, gsTables.ls = gsTables.ls, gsData.ls = gsData.ls, phData.ls = phData.ls, config.df = config.df)
}

# CnvGSAOutput accessors
setGeneric( "res.ls", function(obj) standardGeneric("res.ls") )	## Already defined under CnvGSAInput accessors
setMethod( "res.ls", "CnvGSAOutput", function(obj){ obj@res.ls } )

setGeneric( "gsTables.ls", function(obj) standardGeneric("gsTables.ls") )
setMethod( "gsTables.ls", "CnvGSAOutput", function(obj){ obj@gsTables.ls } )

setGeneric( "gsData.ls", function(obj) standardGeneric("gsData.ls") )
setMethod( "gsData.ls", "CnvGSAOutput", function(obj){ obj@gsData.ls } )

setGeneric( "phData.ls", function(obj) standardGeneric("phData.ls") )
setMethod( "phData.ls", "CnvGSAOutput", function(obj){ obj@phData.ls } )

setGeneric( "config.df", function(obj) standardGeneric("config.df") )
setMethod( "config.df", "CnvGSAOutput", function(obj){ obj@config.df } )