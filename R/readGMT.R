
##
## Sample data records
##
## - cf.
## 		http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29
##  	http://www.broadinstitute.org/gsea/msigdb/collections.jsp
##
## *** N.B. columns are actually tab-delimited ***
##
## PROSTATE_GLAND_DEVELOPMENT                      GO:0030850        2736    5176    9241    8626    ...
## REGULATION_OF_EPITHELIAL_CELL_DIFFERENTIATION   GO:0030856        595     54206   8626    4435    ...
## EPITHELIAL_CELL_DIFFERENTIATION                 GO:0030855        56033   2302    3713    353142  ...
## [...]
##
## 1                                               2                 3...
## gs.name                                         gs.description    gs.genes
##


readGMT <- function( filename )
{
	## Input checks
	if( missing(filename) ) {
		stop( "Missing 'filename' argument" )
	}

	## Read in the records
	gs.name <- character()
	gs.description <- character()
	gs.genes <- list()
	con <- file( filename, "r" )
	recnum <- 0
	repeat
	{
		hrec <- readLines( con, 1 )

		if( length(hrec) == 0 ) {
			break
		}
		recnum <- recnum + 1
		
		fields <- unlist( strsplit( hrec, "\t" ) )
		
		gs.name[recnum] <- fields[1]						## e.g. "prostate gland development"
		gs.description[recnum] <- fields[2]					## e.g. "GO:0030850"
		gs.genes[[recnum]] <- fields[ 3:length(fields) ]	## e.g. chr [1:42] "2736" "5176" "9241" "8626" ...
	}
	close(con)
	
	## Build the output
	names( gs.genes ) <- gs.description
	names( gs.name ) <- gs.description
	gsData <- list( gs2gene=gs.genes, gs2name=gs.name )

	return( gsData )
}
