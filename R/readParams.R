
readParamsRFile <- function( filename )
{
	## Input file (see sampleparams.R) should contain something along the lines of:
	##
	##     # Main test parameters
	##     grandtotals_mode <- "all"
	##     sample_classes   <- c( "case", "ctrl" )
	##     fdr_iter         <- 2
	##     extended_report  <- 200
	##     boxplot_PDFs     <- FALSE
	##
	##     # cnvData$filters parameters
	##     limits_type      <- "DEL"
	##     rem_genes        <- c( "9696", "3106" )		# CROCC and HLA-B
	##
	## 'params', the output object, can then be built in a striaghtforwardly
	## after sourcing the file.

	## The following assignments are a workaround for "no visible binding for global variable" notes
	## in the 'R CMD check' output
	grandtotals_mode <- NULL
	sample_classes <- NULL
	fdr_iter <- NULL
	extended_report <- NULL
	boxplot_PDFs <- NULL
	limits_type <- NULL
	Type <- NULL
	Max_length <- NULL
	Max_gcount <- NULL
	rem_genes <- NULL

	source( filename, local=TRUE )

	params <- list(
		grandtotals_mode = grandtotals_mode,
		sample_classes = sample_classes,
		fdr_iter = fdr_iter,
		extended_report = extended_report,
		boxplot_PDFs = boxplot_PDFs
	)
	
	## Now deal with the parameters related to $filters...
	
	filters <- list()

	## limits_type
	if( ! is.null(limits_type) )	# Before "no visible binding for global variable" fix above, was: if( is.element( "limits_type", ls() ) )
	{
		filters$limits_type <- limits_type
	}
	
	## limits_size
	if( ! is.null(Type) )	# Before "no visible binding for global variable" fix above, was: if( is.element( "Type", ls() ) )
	{
		## First check that the other two parts are there
		if( is.null(Max_length) )	# Before "no visible binding for global variable" fix above, was: if( ! is.element( "Max_length", ls() ) )
		{
			stop( "Missing 'Max_length' parameter" )
		}
		else if( is.null(Max_gcount) )	# Before "no visible binding for global variable" fix above, was: else if( ! is.element( "Max_gcount", ls() ) )
		{
			stop( "Missing 'Max_gcount' parameter" )
		}
		
		filters$limits_size <- list(
			Type = Type,
			Max_length = Max_length,
			Max_gcount = Max_gcount
		)
	}
	
	## rem_genes
	if( ! is.null(rem_genes) )	# Before "no visible binding for global variable" fix above, was: if( is.element( "rem_genes", ls() ) )
	{
		filters$rem_genes <- rem_genes
	}
	
	params$filters <- filters

	## Finally, return the output.
	return( params )
}
