
getCnvGenes <- function( cnv, genemap, delim )
{
	## Input checks
	if( missing(cnv) ) {
		stop( "Missing 'cnv' argument" )
	}
	if( missing(genemap) ) {
		stop( "Missing 'genemap' argument" )
	}
	if( missing(delim) ) {
		stop( "Missing 'delim' argument" )
	}

	.getCnvGenes <- function( cnv.Chr, cnv.Coord_i, cnv.Coord_f, genemap )
	{
		Chr <- Coord_f <- Coord_i <- GeneID <- NULL	## workaround for "no visible binding for global variable" note in 'R CMD check' output (due to next lines of code)
		cnvGenes <- subset(
			genemap,
			(
				(Chr == cnv.Chr) &
				((cnv.Coord_i < Coord_f) & (cnv.Coord_f > Coord_i))
			),
			select = GeneID,
			drop = TRUE
		)
		return( paste( unique(cnvGenes), collapse=delim ) )
	}

	genes <- with( cnv,
		mapply( .getCnvGenes,
			Chr,
			Coord_i,
			Coord_f,
			MoreArgs = list( genemap = genemap ),
			SIMPLIFY = FALSE
		)
	)
	genes <- unname(unlist(genes))

	return( genes )
}
