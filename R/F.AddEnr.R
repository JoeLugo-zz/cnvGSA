# -----------------------
# Additional data/tests
# -----------------------

# adds additional data to enrichment results
# includes additional tests
F.AddEnr <- function (enrRes.ls, cnvData.ls, gsData.ls, geneData.ls, addEnrPar.ls, assTestPar.ls)
	{
	addEnrPar.ls <- f.complete_enrpar (addEnrPar.ls)
		
	## enr.df <- subset (enrRes.ls$basic, subset = FET_fdr <= addEnrPar.ls$fdr_thr)
	enr.df <- enrRes.ls$basic[1: addEnrPar.ls$sel_n, ]
		
	# add support data
	enr.df <- f.add_support(
		enr.df,
		cnvData.ls,
		gsData.ls,
		geneData.ls,
		assTestPar.ls$test_classes[1]
	)
	enr.df <- f.add_support(
		enr.df,
		cnvData.ls,
		gsData.ls,
		geneData.ls,
		assTestPar.ls$test_classes[2]
	)

	# add samples (by class)
	enr.df <- f.add_attrib (enr.df, cnvData.ls, attrib.name = "SampleID")
	enr.df <- f.add_attrib (enr.df, cnvData.ls, attrib.name = "CnvID")
	
	# Remove top associated gene and recompute association
	enr.df <- f.rem_top (enr.df, cnvData.ls, gsData.ls, geneData.ls, assTestPar.ls, enrRes.ls)

	enrRes.ls$extended <- enr.df

	# Add tables showing subset of cnvData$full for for each gene set
	enrRes.ls$gstables <- f.add_gstables( enr.df, cnvData.ls, geneData.ls, cnvData.ls$gsep )

	return (enrRes.ls)
	}

f.complete_enrpar <- function (addEnrPar.ls)
	{
	if (is.null (addEnrPar.ls$sel_n)) addEnrPar.ls$sel_n <- 75 
	if (is.null (addEnrPar.ls$iter_comp)) addEnrPar.ls$iter_comp <- 2000
	return (addEnrPar.ls)
	}

f.add_support <- function (enr.df, cnvData.ls, gsData.ls, geneData.ls, class.ch )
	{
	sel.gs <- enr.df$GsID

#	gs_supp_id.ls <- lapply (gsData.ls$gs2gene[sel.gs], intersect, y = geneData.ls$support)
	support_class <- paste( sep="", "support_", class.ch )
	gs_supp_id.ls <- lapply( gsData.ls$gs2gene[sel.gs], intersect, y = geneData.ls[[support_class]] )

	f.conv <- function (in.genes, id2attr.chv)
		{return (id2attr.chv[in.genes])}
	gs_supp_sy.ls <- lapply (gs_supp_id.ls, f.conv, id2attr.chv = geneData.ls$ann$gene2sy)

	gs_u.genes <- unique (unlist (gsData.ls$gs2gene))
#	exp.n <- length (intersect (geneData.ls$support, gs_u.genes)) / length (gs_u.genes)
	exp.n <- length (intersect (geneData.ls[[support_class]], gs_u.genes)) / length (gs_u.genes)
	
#	enr.df$Support_size <- sapply (gs_supp_id.ls, length)
#	enr.df$Support_ratio <- (enr.df$Support_size / enr.df$GsSize) / exp.n
#	enr.df$Support_geneid <- sapply (gs_supp_id.ls, paste, collapse = cnvData.ls$gsep)
#	enr.df$Support_symbol <- sapply (gs_supp_sy.ls, paste, collapse = cnvData.ls$gsep)
	support_size_class <- paste( sep="", "Support_size_", class.ch )
	support_ratio_class <- paste( sep="", "Support_ratio_", class.ch )
	support_geneid_class <- paste( sep="", "Support_geneid_", class.ch )
	support_symbol_class <- paste( sep="", "Support_symbol_", class.ch )
	enr.df[[ support_size_class ]] <- sapply (gs_supp_id.ls, length)
	enr.df[[ support_ratio_class ]] <- (enr.df[[support_size_class]] / enr.df$GsSize) / exp.n
	enr.df[[ support_geneid_class ]] <- sapply (gs_supp_id.ls, paste, collapse = cnvData.ls$gsep)
	enr.df[[ support_symbol_class ]] <- sapply (gs_supp_sy.ls, paste, collapse = cnvData.ls$gsep)
		
	return (enr.df)
	}

f.add_attrib <- function (enr.df, cnvData.ls, attrib.name)
	{
	classes.chv <- cnvData.ls$uni$class

	f.subsetByclass <- function (class.ch, cnv.df)
		{
		Class <- NULL	## workaround for "no visible binding for global variable" note in 'R CMD check' output (due to next line of code)
		cnv_sub.df <- subset (cnv.df, subset = Class == class.ch)[, c (attrib.name, "GsID")]
		colnames (cnv_sub.df)[1] <- "Attrib"
		cnv_sub.df <- cnv_sub.df[! duplicated (cnv_sub.df), ]
		return (cnv_sub.df)
		}
	byclass.ls <- lapply (classes.chv, f.subsetByclass, cnvData.ls$full)
	names (byclass.ls) <- classes.chv

	f.getAttrib <- function (cnv.df)	
		{
		out.ls <- aggregate (Attrib ~ GsID, cnv.df, paste, collapse = ";")
		out.attrib <- out.ls$Attrib
		names (out.attrib) <- out.ls$GsID
		return (out.attrib)
		}
	byGs.ls <- lapply (byclass.ls, f.getAttrib)
	
	all.gsid <- enr.df$GsID
	c.n <- ncol (enr.df)
	l.n <- length (classes.chv)
	for (i in 1: l.n)
		{
		enr.df <- cbind (enr.df, byGs.ls[[i]][all.gsid])
		}
	colnames (enr.df)[(c.n + 1): (c.n + l.n)] <- paste (classes.chv, attrib.name, sep = "_")
	
	return (enr.df)
	}

f.rem_top <- function (enr.df, cnvData.ls, gsData.ls, geneData.ls, assTestPar.ls, enrRes.ls)
	{
	sel.gs <- enr.df$GsID
	gs2eg.ls <- lapply (gsData.ls$gs2gene[sel.gs], setdiff, y = cnvData.ls$rem_genes)
	rem.genes <- cnvData.ls$filters$rem_genes
	
	# reduce the data passed to the unit to enhance performance
	Gcount <- SampleID <- Class <- GsID <- Genes <- CnvID <- NULL	## workaround for "no visible binding for global variable" note in 'R CMD check' output (due to next line of code)
	cnv.df <- subset (cnvData.ls$full, subset = Gcount > 0, select = c (SampleID, Class, GsID, Genes, CnvID))
	
	top.genes <- sapply (gs2eg.ls, f.rem_top_topgene_unit, geneData.ls, rem.genes)

	pvalue.nv <- mapply (f.rem_top_FET_unit, sel.gs, top.genes,
						MoreArgs = list (
							cnv.df = cnv.df, 
							rem.genes = rem.genes,
							totals.nv = enrRes.ls$totals, 
							classes.chv = assTestPar.ls$test_classes)
						)
						
	enr.df$FETpv_remTop <- pvalue.nv
	enr.df$FETfdr_remTop <- sapply (pvalue.nv, f.rem_top_fdr_unit, enr.df)
	enr.df$Topgene <- geneData.ls$ann$gene2sy[top.genes]
	
	return (enr.df)
	}

f.rem_top_topgene_unit <- function (gs.genes, geneData.ls, rem.genes)
	{
	GeneID <- NULL	## workaround for "no visible binding for global variable" note in 'R CMD check' output (due to next line of code)
	top.gene  <- as.character (subset (geneData.ls$gcounts, GeneID %in% gs.genes & ! GeneID %in% rem.genes, select = GeneID, drop = T)[1])
	return (top.gene)
	}

f.rem_top_FET_unit <- function (gs.id, top.gene, cnv.df, rem.genes, totals.nv, classes.chv)
	{
	Genes <- CnvID <- SampleID <- Class <- GsID <- NULL	## workaround for "no visible binding for global variable" note in 'R CMD check' output (due to verbatim column names in sample() calls in the rest of the code here)
	
	# CNVs that carry the gene to be removed
	cnv_topgene.id <- as.character (subset (cnv.df, subset = Genes %in% top.gene, select = CnvID, drop = T))
	
	cnv.df <- subset (cnv.df, subset = ! CnvID %in% cnv_topgene.id, select = c (SampleID, Class, GsID))
	cnv.df <- cnv.df[! duplicated (cnv.df), ]
	
	c1.n <- length (unique (subset (cnv.df, subset = Class == classes.chv[1] & GsID == gs.id, select = SampleID, drop = T)))
	c2.n <- length (unique (subset (cnv.df, subset = Class == classes.chv[2] & GsID == gs.id, select = SampleID, drop = T)))

	contingency.mx <- matrix (data = c (c1.n, c2.n, totals.nv[1] - c1.n, totals.nv[2] - c2.n), ncol = 2, nrow = 2, byrow = T)

	pvalue.n <- fisher.test (contingency.mx, alternative = "greater")$p.value
	names (pvalue.n) <- top.gene
	
	return (pvalue.n)
	}

f.rem_top_fdr_unit <- function (pvalue.n, enr.df)
	{
	fdr.ix <- which.min (abs (pvalue.n - enr.df$FET_pv))
	fdr.n <- enr.df$FET_fdr[fdr.ix]

	return (fdr.n)
	}

f.add_gstables <- function( enr.df, cnvData, geneData, gsep )
{
	## Convert geneData$ann$gene2sy to a data frame so that it can be used with merge()
	gene2sy.df <- data.frame(
		Genes = names(geneData$ann$gene2sy),
		Symbols = unname(geneData$ann$gene2sy)
	)

	## Function that returns the table for a single gene-set
	f.gstable_unit <- function( gsid, cnvData, gene2sy.df, gsep )
	{
		GsID <- NULL	## workaround for "no visible binding for global variable" note in 'R CMD check' output (due to next line of code)
		gstable <- subset( cnvData$full, GsID == gsid )

		## Add column for gene symbols
		gstable <- merge( gstable, gene2sy.df )
		
		## Remove levels from these two cols; they are not needed after the merge()
		gstable$Genes <- as.character( gstable$Genes )
		gstable$Symbols <- as.character( gstable$Symbols )

		## Sort by gene symbol so that the output -- which will have the genes and their symbols
		## in unstacked form -- will have them in alphabetical order. (Code here is admittedly
		## cryptic; see man page for the order() function -- in particular, the examples)
		temp1 <- gstable[ , c("Symbols", "CnvID", "Genes")]
		temp2 <- temp1[ do.call(order, temp1), ]

		## Unstack the genes and symbols
		## N.B. if all CnvIDs have only one gene (i.e. nothing ends up getting unstacked),
		## unstack() returns a *data.frame* instead of a list-of-lists; hence the slight
		## complication with the if() blocks below. (Genes and their symbols should be
		## matched -- cf. the merge() operation above -- so a single if() statement should
		## be sufficient...)
		temp3_genes <- unstack( temp2, Genes ~ CnvID )
		temp3_symbols <- unstack( temp2, Symbols ~ CnvID )
		if( class(temp3_genes) == "data.frame" ) {
			temp4_genes <- as.list( as.character( temp3_genes[,1] ) )
			names(temp4_genes) <- rownames( temp3_genes )
			temp4_symbols <- as.list( as.character( temp3_symbols[,1] ) )
			names(temp4_symbols) <- rownames( temp3_symbols )
		}
		else { # should already be in the list-of-lists form as described above
			temp4_genes <- temp3_genes
			temp4_symbols <- temp3_symbols
		}
		
		## Collapse the list-of-lists into a simple character vector where each element
		## has the genes/symbols separated by gsep
		genes <- unlist(lapply( temp4_genes, paste, collapse=gsep ))
		symbols <- unlist(lapply( temp4_symbols, paste, collapse=gsep ))
		genes.df <- data.frame( CnvID = names(genes), Genes = unname(genes) )
		symbols.df <- data.frame( CnvID = names(symbols), Symbols = unname(symbols) )
		
		## Add these columns back into the table
		genes_and_symbols <- merge( genes.df, symbols.df )
		gstable$Genes <- NULL
		gstable$Symbols <- NULL
		gstable <- merge( unique(gstable), genes_and_symbols )

		## Sort output by $Class (again, code is cryptic -- see man page for the order() function)
		gstable <- gstable[ , c(11,1:10,12:13) ]	# Move $Class column to the front
		gstable <- gstable[ do.call(order, gstable), ]
		gstable <- gstable[ , c(2:11,1,12:13) ]		# Move $Class column back to where it was

		return( gstable )
	}

	## Apply the function to all the gene-sets and return the result
	gstables <- lapply( enr.df$GsID, f.gstable_unit, cnvData, gene2sy.df, gsep )
	names(gstables) <- enr.df$GsID
	return( gstables )
}
