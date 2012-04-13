
cnvGSAFisher <- function( input )
{
	output <- F.Pipeline(
		cnvData.ls  = input@cnvData,
		gsData.ls   = input@gsData,
		geneData.ls = input@geneData,
		params.ls	= input@params
	)
	return( output )
}

cnvGSAexportBurdenStats <- function( output, filenamePrefix )
{
	if( missing(output) ) {
		stop( "Missing 'output' argument" )
	}
	if( missing(filenamePrefix) ) {
		stop( "Missing 'filenamePrefix' argument" )
	}
	
	write.table(
		output@burdenGs$coverage,
		paste( sep="", filenamePrefix, "_burdenGs_coverage.txt" ),
		sep="\t"
	)
	
	write.table(
		output@burdenGs$pairs,
		paste( sep="", filenamePrefix, "_burdenGs_pairs.txt" ),
		sep="\t"
	)
}
