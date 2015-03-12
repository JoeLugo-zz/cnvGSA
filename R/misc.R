
writeSampleGvf.SingleRecord <- function( cnv )
{
	ret <- paste( sep="\t",
		paste( "Chr", cnv$Chr, sep="" ),
		"coreCNV_AsdCt.RData",
		ifelse( cnv$Type == "DUP", "copy_number_gain", "copy_number_loss" ),
		cnv$Coord_i,
		cnv$Coord_f,
		".",
		".",
		".",
		paste( sep=";",
			paste( sep="", "ID=", cnv$CnvID ),
			paste( sep="", "Name=", cnv$CnvID, "_name" ),
			"var_origin=NA",
			paste( sep="", "Start_range=", cnv$Coord_i, ",." ),
			paste( sep="", "End_range=.,", cnv$Coord_f ),
			paste( sep="", "samples=", cnv$SampleID )
		)
	)
	return( ret )
}

writeSampleGvf <- function( cnvData, outputFile )
{
	cat( file=outputFile, sep="\n",
		writeSampleGvf.SingleRecord( cnvData$cnv )
	)
}

writeSampleGmt.SingleRow <- function( gs.name, gs.description, gs.genes )
{
	ret <- paste( sep="\t",
		toupper( gsub(" ", "_", gs.name) ),
		gs.description,
		paste( gs.genes, collapse="\t" )
	)
	return(ret)
}

writeSampleGmt <- function( gsData, outputFile )
{
	gs.genes <- unname(gsData$gs2gene)			## e.g. chr [1:42] "2736" "5176" "9241" "8626" ...
	gs.description <- names(gsData$gs2gene)		## e.g. "GO:0030850"
	gs.name <- unname(gsData$gs2name)			## e.g. "prostate gland development"

	cat( file=outputFile, "" )
	for( i in 1:length(gs.name) )
	{
		cat( file=outputFile, append=TRUE, sep="\n",
			writeSampleGmt.SingleRow(
				gs.name[i],
				gs.description[i],
				gs.genes[[i]]	## N.B. double brackets since $gs2gene is itself a list.
			)		
		)
	}
}

sort.chromosomeLabels <- function( chrs )
{
	chrs[ chrs == "X" ] <- "23"
	chrs[ chrs == "Y" ] <- "24"
	chrs <- as.character( sort(as.integer(chrs)) )
	chrs[ chrs == "23" ] <- "X"
	chrs[ chrs == "24" ] <- "Y"

	return( chrs )
}

