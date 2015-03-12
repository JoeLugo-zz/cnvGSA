
##
## CnvGSAInput
##

setClass( "CnvGSAInput",
	representation(
		cnvData = "list",
		gsData = "list",
		geneData = "list",
		params = "list"
	)
)

# CnvGSAInput constructor
CnvGSAInput <- function(
					cnvData = list(),
					gsData = list(),
					geneData = list(),
					params = list()
				)
{
    new( "CnvGSAInput", cnvData = cnvData, gsData = gsData, geneData = geneData, params = params )
}

# CnvGSAInput accessors

setGeneric( "cnvData", function(obj) standardGeneric("cnvData") )
setGeneric( "cnvData<-", function(obj, value) standardGeneric("cnvData<-") )
setMethod( "cnvData", "CnvGSAInput", function(obj){ obj@cnvData } )
setReplaceMethod( "cnvData", "CnvGSAInput",	function(obj, value){ obj@cnvData <- value } )

setGeneric( "gsData", function(obj) standardGeneric("gsData") )
setGeneric( "gsData<-", function(obj, value) standardGeneric("gsData<-") )
setMethod( "gsData", "CnvGSAInput", function(obj){ obj@gsData } )
setReplaceMethod( "gsData", "CnvGSAInput",	function(obj, value){ obj@gsData <- value } )

setGeneric( "geneData", function(obj) standardGeneric("geneData") )
setGeneric( "geneData<-", function(obj, value) standardGeneric("geneData<-") )
setMethod( "geneData", "CnvGSAInput", function(obj){ obj@geneData } )
setReplaceMethod( "geneData", "CnvGSAInput",	function(obj, value){ obj@geneData <- value } )

setGeneric( "params", function(obj) standardGeneric("params") )
setGeneric( "params<-", function(obj, value) standardGeneric("params<-") )
setMethod( "params", "CnvGSAInput", function(obj){ obj@params } )
setReplaceMethod( "params", "CnvGSAInput",	function(obj, value){ obj@params <- value } )


##
## CnvGSAOutput
##

setClass( "CnvGSAOutput",
	representation(
		cnvData = "list",
		burdenSample = "list",
		burdenGs = "list",
		geneData = "list",
		enrRes = "list"
	)
)

# CnvGSAOutput accessors

#setGeneric( "cnvData", function(obj) standardGeneric("cnvData") )	## Already defined under CnvGSAInput accessors
setMethod( "cnvData", "CnvGSAOutput", function(obj){ obj@cnvData } )

setGeneric( "burdenSample", function(obj) standardGeneric("burdenSample") )
setMethod( "burdenSample", "CnvGSAOutput", function(obj){ obj@burdenSample } )

setGeneric( "burdenGs", function(obj) standardGeneric("burdenGs") )
setMethod( "burdenGs", "CnvGSAOutput", function(obj){ obj@burdenGs } )

#setGeneric( "geneData", function(obj) standardGeneric("geneData") )	## Already defined under CnvGSAInput accessors
setMethod( "geneData", "CnvGSAOutput", function(obj){ obj@geneData } )

setGeneric( "enrRes", function(obj) standardGeneric("enrRes") )
setMethod( "enrRes", "CnvGSAOutput", function(obj){ obj@enrRes } )
