
setMethod( "summary",
	signature = "CnvGSAOutput",
	definition = function( object, verbose=FALSE )
	{
		cat( "\n" )
		cat( "CNV GSA Output Summary\n" )
		cat( "----------------------\n" )

		s2class <- object@cnvData$s2class
		sampleClasses <- unique(s2class$Class)
		m <- merge( object@cnvData$cnv, object@cnvData$s2class )

		## Number of sample classes (e.g. 2) and their label in the order of comparison (e.g. case, ctrl)
		cat( sep="", "Number of sample classes: ", length(sampleClasses), "\n" )
		cat( sep="", "Sample classes: ", paste(sampleClasses, collapse=" "), "\n" )

		## Number of unique samples, in total and by class, found in s2class dataframe
		## and found in cnv dataframe object (also these in total and by class)
		cat( sep="", "Number of unique samples, total: ", length(s2class$SampleID), "\n" )
		for( sc in sampleClasses ) {
			cat( sep="", "Number of unique samples, class '", sc, "' (all): ",
				length( unique( s2class$SampleID[ s2class$Class == sc ] ) ), "\n" )
			cat( sep="", "Number of unique samples, class '", sc, "' (with CNVs): ",
				length( unique( m$SampleID[ m$Class == sc ] ) ), "\n" )
		}

		## CNV type labels, and number of CNVs per type
		cnvTypes <- unique( object@cnvData$cnv$Type )
		cat( sep="", "CNV types: ", paste( collapse=" ", cnvTypes ), "\n" )
		for( type in cnvTypes ) {
			cat( sep="", "Number of CNVs of type '", type, "': ",
				length( m$Type[ m$Type == type ] ), "\n" )
		}
		
		## Chromosome label list; add min and max CNV coordinate per chromosome if "verbose details" = T
		chromosomes <- sort.chromosomeLabels( unique( m$Chr ) )
		cat( sep="", "Chromosomes: ", paste( collapse=", ", chromosomes ), "\n" )
		if( verbose == TRUE )
		{
			cat( "Min/max CNV coordinates per chromosome:\n" )
			cat( "\n" )
		
			minCoords <- data.frame(
				Coord_i = object@cnvData$full$Coord_i,
				Chr = object@cnvData$full$Chr,
				stringsAsFactors = FALSE
			)
			maxCoords <- data.frame(
				Coord_f = object@cnvData$full$Coord_i,
				Chr = object@cnvData$full$Chr,
				stringsAsFactors = FALSE
			)

			minCoords <- unlist(lapply( unstack(minCoords), min ))
			maxCoords <- unlist(lapply( unstack(maxCoords), max ))
			mm <- data.frame(
				Chr = names(minCoords),
				Min = unname(minCoords),
				Max = unname(maxCoords),
				stringsAsFactors=FALSE
			)

			mm$Chr[ mm$Chr == "X" ] <- "23"
			mm$Chr[ mm$Chr == "Y" ] <- "24"
			mm <- mm[ with(mm, order(as.integer(Chr))) , ]
				## cf. http://stackoverflow.com/questions/1296646/how-to-sort-a-dataframe-by-columns-in-r
			mm$Chr[ mm$Chr == "23" ] <- "X"
			mm$Chr[ mm$Chr == "24" ] <- "Y"

			print( mm, row.names=FALSE )
			
			cat( "\n" )
		}
			
		## Number of (unique) genes: total; by sample class; by CNV types; by sample class and CNV type

		## Clearest way to show this would be a table, e.g.:
		##
		##		Unique genes:
		##	
		##				case	ctrl	tot
		##		DEL		.		.		.
		##		DUP		.		.		.
		##		tot		.		.		.
		
		classes        <- object@cnvData$uni$class
		types          <- object@cnvData$uni$type
		selected_type  <- object@cnvData$filters$limits_type
		totals         <- object@geneData$totals

		n_types <- length(types)
		n_classes <- length(classes)
		ug <- as.data.frame(
			matrix( rep(0, n_types*n_classes), ncol = n_classes )
		)

		rownames(ug) <- types
		ug[ selected_type, ] <- totals$cnv

		ug <- cbind( ug, apply(ug, 1, sum) )
		ug <- rbind( ug, apply(ug, 2, sum) )
		names(ug) <- c( classes, "tot" )
		rownames(ug) <- c( types, "tot" )

		cat( "Unique genes:\n" )
		cat( "\n" )
		print(ug)

		## End of output
		cat( "\n" )
	}
)
