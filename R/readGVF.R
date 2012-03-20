
##
## Sample data records
##
## - cf. http://www.sequenceontology.org/wiki/index.php/GVF_Examples (27-Jan-2012)):
##          ### check this (the description doesn't seem to match the samples in DGV)
##
## Main record (*** N.B. columns are actually tab-delimited ***):
##
## Chr17    DGVa    copy_number_gain        31450358    31501499    .   .   .   ID=essv104508;Name=essv104508;var_origin=Not tested;Start_range=31450358,.;End_range=.,31501499;samples=NA18964
## Chr17    DGVa    copy_number_variation   31450400    31496600    .   .   .   ID=esv35108;Name=esv35108;var_origin=Not tested;Start_range=31450400,.;End_range=.,31496600
## Chr10    DGVa    copy_number_loss        46363383    47212100    .   .   .   ID=essv103490;Name=essv103490;var_origin=Not tested;Start_range=46363383,.;End_range=.,47212100;samples=NA18861
## [...]
## 
## 1        2       3                       4           5           6   7   8   9
## Chr      Src_DB  Type                    Coord_i     Coord_f     ??? ??? ??? (various attributes)
##
##
## Field 9 of main record:
##
## ID=essv104508;Name=essv104508;var_origin=Not tested;Start_range=31450358,.;End_range=.,31501499;samples=NA18964
## ID=esv35108;Name=esv35108;var_origin=Not tested;Start_range=31450400,.;End_range=.,31496600
## :
## 1             2               3                     4                      5                    6
## ???           ???             ???                   Coord_i                Coord_f              SampleID


readGVF <- function( filename )
{
	## Input checks
	if( missing(filename) ) {
		stop( "Missing 'filename' argument" )
	}

	## Initialize vectors in the output data frame and read in the records
	SampleID <- character()
	Chr <- character()
	Coord_i <- integer()
	Coord_f <- integer()
	Type <- character()
	#Genes	## defined at the end
	#CnvID	## defined at the end

	con <- file( filename, "r" )
	repeat
	{
		hrec <- readLines( con, 1 )

		if( length(hrec) == 0 ) {
			break
		}
		if( substr( hrec, 1, 1 ) == "#" ) {
			next
		}

		fields <- unlist( strsplit( hrec, "\t" ) )
		
		## SampleID
		field9_6 <- unlist( strsplit( fields[9], ";" ) )[ 6 ]
		SampleID <- c( SampleID, substr( field9_6, 9, nchar(field9_6) ) )
		
		## Chr[omosome]
		Chr <- c( Chr, substr( fields[1], 4, nchar(fields[1]) ) )

		## Coord_i and Coord_f
		Coord_i <- c( Coord_i, as.integer(fields[4]) )
		Coord_f <- c( Coord_f, as.integer(fields[5]) )

		## Type
		val <- ""
		if( fields[3] == "copy_number_gain" ) {
			val <- "DUP"
		}
		else if( fields[3] == "copy_number_loss" ) {
			val <- "DEL"
		}
		else if( fields[3] == "copy_number_variation" ) {
			val <- "VAR"
		}
		Type <- c( Type, val )
	}
	close( con )

	## Genes: needs to be assigned by checking the CNVs against a gene coordinate map
	## (cf. getCnvGenes()), so just create a vector of empty strings for now.
	Genes <- rep( "", length(SampleID) )
	
	## CnvID
	CnvID <- paste( "CNV_", 1:length(SampleID), sep="" )

	## Build the output
	cnv <- data.frame( stringsAsFactors = FALSE,
		SampleID,
		Chr,
		Coord_i,
		Coord_f,
		Type,
		Genes,
		CnvID
	)

	return( cnv )
}
