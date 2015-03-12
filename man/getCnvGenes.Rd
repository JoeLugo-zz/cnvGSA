\name{getCnvGenes}
\alias{getCnvGenes}
\title{Determine genes hit by CNVs}
\description{
Takes as input a data frame of CNVs and a data frame of gene coordinates and returns the genes hit by each CNV.
}
\usage{
getCnvGenes( cnv, genemap, delim )
}
\arguments{
	\item{cnv}{
		Data frame containing the CNVs. It should contain at least the following columns:
		- Chr: Chromosome on which the CNV is found (e.g. ``1'', ``12'', ``X'')
		- Coord_i: Initial coordinate (aka start position) of the CNV on the chromosome
		- Coord_f: Final coordinate (aka end position) of the CNV on the chromosome
	}
	\item{genemap}{
		Data frame containing genes along with the chromosomes and coordinates for each. The columns should be similar to the \code{cnv} data frame:
		- Chr: Chromosome on which the gene is found (e.g. ``1'', ``12'', ``X'')
		- Coord_i: Initial coordinate (aka start position) of the gene on the chromosome
		- Coord_f: Final coordinate (aka end position) of the gene on the chromosome
		- GeneID: Gene ID
	}
	\item{delim}{
		Character to be used to separate genes in the output for each CNV.
	}
}
\value{
A vector in which each element contains a delimited string of the genes that fall within the range of the corresponding CNV in the input.
}
\author{Robert Ziman \email{rziman@gmail.com}}
\examples{
library("cnvGSAdata")

## Read in the example gene map
genemapFile <- system.file(
	"extdata",
	"merge_00k_flank_hg18_refGene_jun_2011_exon.gff",
	package = "cnvGSAdata"
)
fields <- read.table(
	genemapFile,
	sep = "\t",
	comment.char = "",
	quote = "\"",
	header = FALSE,
	stringsAsFactors = FALSE
)
genemap <- data.frame(
	Chr = fields[,1],
	Coord_i = fields[,4],
	Coord_f = fields[,5],
	GeneID = fields[,11],
	stringsAsFactors = FALSE
)
genemap$Chr <- sub(genemap$Chr, pattern = "chr", replacement = "")

## Read in a few CNVs
cnvFile <- system.file("extdata", "cnv.gvf", package="cnvGSAdata")
cnv <- readGVF(cnvFile)
cnv <- cnv[1:5,]

## Get CNV genes
delim <- ";"
genes <- getCnvGenes(cnv, genemap, delim)
}
