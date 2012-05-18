# -----------------------
# Additional data/tests
# -----------------------

# adds additional data to enrichment results
# includes additional tests
F.AddEnr <- function (enrRes.ls, cnvData.ls, gsData.ls, geneData.ls, burdenSample.ls, assTestPar.ls, addEnrPar.ls)
	{
	addEnrPar.ls <- f.complete_enrpar (addEnrPar.ls)
	
	# Initialize enr.df using only the top 200 (or whatever addEnrPar.ls$sel_n is set to) gene-sets
	enr.df <- enrRes.ls$basic[ (1: min (nrow (enrRes.ls$basic), addEnrPar.ls$sel_n)), ]

	# Logistic regression
	if (!is.null (addEnrPar.ls$do_logistic)) {
		if (addEnrPar.ls$do_logistic == "extended") {
			cat ("Logistic regression (for only those gene-sets in the extended report)...")
			enr.df <- f.add_lrm (enr.df, cnvData.ls, burdenSample.ls, addEnrPar.ls)
			cat ("done\n")
		}
	}

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
#	support_ratio_class <- paste( sep="", "Support_ratio_", class.ch )		## 2012-04-25 RZ: Removed as per vignette corrections
	support_geneid_class <- paste( sep="", "Support_geneid_", class.ch )
	support_symbol_class <- paste( sep="", "Support_symbol_", class.ch )
	enr.df[[ support_size_class ]] <- sapply (gs_supp_id.ls, length)
#	enr.df[[ support_ratio_class ]] <- (enr.df[[support_size_class]] / enr.df$GsSize) / exp.n	## 2012-04-25 RZ: Removed as per vignette corrections
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

# logistic regression model
f.add_lrm <- function (enr.df, cnvData.ls, burdenSample.ls, addEnrPar.ls)
	{

	# 1. generate the following columns (LogLenTot and GsGene_N_Tot obtained from F.BurdenSample)
	SampleID <- c (
		burdenSample.ls$stat.ls[[1]]$LogLenMean$SampleID,	# Using [[1]] and [[2]] here since the stat.ls list elements corresponding to the case and ctrl classes
		burdenSample.ls$stat.ls[[2]]$LogLenMean$SampleID	# can have names other than "case" and "ctrl"; cf. F.BurdenSample
	)
	Class <- c (
		rep (1, length (burdenSample.ls$stat.ls[[1]]$LogLenMean$SampleID)),	# brglm requires y-values to be between 0 and 1, hence 1=case, 0=ctrl
		rep (0, length (burdenSample.ls$stat.ls[[2]]$LogLenMean$SampleID))
	)
	LogLenTot <- c (
		burdenSample.ls$stat.ls[[1]]$LogLenTot$Stat,
		burdenSample.ls$stat.ls[[2]]$LogLenTot$Stat
	)
	GsGene_N_Tot <- c (
		burdenSample.ls$stat.ls[[1]]$GsGene_N_Tot$Stat,
		burdenSample.ls$stat.ls[[2]]$GsGene_N_Tot$Stat
	)
	data.df <- data.frame (SampleID, Class, LogLenTot, GsGene_N_Tot)
	
	# Ensure that the order of data in the above columns corresponds exactly to that of cnvData.ls$tab$gen
	# (necessary because f.lrm_unit puts columns from cnvData.ls$tab$gen side-by-side with data.df)
	SampleID.df <- data.frame( SampleID=rownames(cnvData.ls$tab$gen) )
	data.df <- merge( SampleID.df, data.df )

	# Perform the logistic regression on each gene-set
	stats.mx <- t (sapply (
		1:length(enr.df$GsID),
		f.lrm_unit,
		enr.df$GsID,
		cnvData.ls,
		data.df
	))
	
	# Add Benjamini-Hochberg corrected p-values
	stats.mx[,"LM1_Dmod_bhPv"] <- p.adjust (stats.mx[,"LM1_Dmod_Pv"], method = "BH")
	stats.mx[,"LM2_Dmod_bhPv"] <- p.adjust (stats.mx[,"LM2_Dmod_Pv"], method = "BH")
		
	return (cbind (enr.df, stats.mx))
	}

f.lrm_unit <- function (gs.n, GsID.chv, cnvData.ls, data.df)
	{

	### TODO: Delete Case_MeanGsN, Ctrl_MeanGsN, Oddr_MeanGsN, Case%_GsN>0, trl%_GsN>0, 
	### Not0Oddr columns and code after confirming that they are no longer necessary

	# Progress indicator
	if (gs.n %% 25 == 0) {cat (gs.n); cat ("; ")}

	stats.names <- c (

		"LM1_GsN_Sign",
		"LM1_GsN_Pv",
		"LM1_Convg",
#		"LM1_Case_MeanGsN",
#		"LM1_Ctrl_MeanGsN",
#		"LM1_Oddr_MeanGsN",
#		"LM1_Case%_GsN>0",
#		"LM1_Ctrl%_GsN>0",
#		"LM1_Not0Oddr",
		"LM1_Dmod_Pv",
		"LM1_Dmod_bhPv",	# Placeholder for Benjamini-Hochberg corrected p-value
		"LM1_PredOddr_full",
		"LM1_PredOddr_base",

		"LM2_GsN_Sign",
		"LM2_GsN_Pv",
		"LM2_Convg",
#		"LM2_Case_MeanGsN",
#		"LM2_Ctrl_MeanGsN",
#		"LM2_Oddr_MeanGsN",
#		"LM2_Case%_GsN>0",
#		"LM2_Ctrl%_GsN>0",
#		"LM2_Not0Oddr",
		"LM2_Dmod_Pv",
		"LM2_Dmod_bhPv",	# Placeholder for Benjamini-Hochberg corrected p-value
		"LM2_PredOddr_full",
		"LM2_PredOddr_base"

	)

	stats.nv <- as.numeric (rep (NA, length (stats.names)))
	names (stats.nv) <- stats.names

	# 2-1. for each gene-set, generate gene-set specific data.frame 'data_gsi.df'
	#      by extracting corresponding line from cnvData.ls$tab$gen for that gene-set.
	#
	# data_gsi.df <- data.frame (SampleID, Class, LogLenTot, GsComp_N, Gs_N)
	# - The SampleIDs should be unique
	# - GsComp_N = GSGene_N_Tot - Gs_N

	Gs_N <- cnvData.ls$tab$gen[ , GsID.chv[gs.n] ]
	GsComp_N <- data.df$GsGene_N_Tot - Gs_N
	data_gsi.df <- data.frame (
		SampleID     = data.df$SampleID,
		Class        = data.df$Class,
		LogLenTot    = data.df$LogLenTot,
		GsGene_N_Tot = data.df$GsGene_N_Tot,
		GsComp_N,
		Gs_N,
		stringsAsFactors=FALSE	### TODO: stringsAsFactors=FALSE: necessary or not?
	)

	# -- logistic model 1 --

	# 2-2. fit logistic models

	full1.brglm <- brglm (Class ~ LogLenTot + GsComp_N + Gs_N, data = data_gsi.df, family = binomial (link = "logit"))
	base1.brglm <- brglm (Class ~ LogLenTot + GsComp_N, data = data_gsi.df, family = binomial (link = "logit"))

	# 2-3. generate/extract statistics
	
	full1.brgsm <- summary (full1.brglm)
	base1.brgsm <- summary (base1.brglm)
		
	gsn1_info.nv <- f.lrm_unit_getSmInfo (full1.brgsm, "Gs_N")

	# logistic model summaries
	stats.nv["LM1_GsN_Sign"] <- gsn1_info.nv["sign"]
	stats.nv["LM1_GsN_Pv"]   <- gsn1_info.nv["pvalue"]
	stats.nv["LM1_Convg"]    <- full1.brglm$converged

#	gsn_case.nv <- subset (data_gsi.df, subset = Class == 1, select = Gs_N, drop = T)
#	gsn_ctrl.nv <- subset (data_gsi.df, subset = Class == 0, select = Gs_N, drop = T)
#	scase.n <- length (gsn_case.nv)
#	sctrl.n <- length (gsn_ctrl.nv)
#	stats.nv["LM1_Case_MeanGsN"] <- mean (gsn_case.nv)
#	stats.nv["LM1_Ctrl_MeanGsN"] <- mean (gsn_ctrl.nv)
#	stats.nv["LM1_Oddr_MeanGsN"]  <- stats.nv["LM1_Case_MeanGsN"] / stats.nv["LM1_Ctrl_MeanGsN"]
#	stats.nv["LM1_Case%_GsN>0"] <- sum (gsn_case.nv > 0) / scase.n
#	stats.nv["LM1_Ctrl%_GsN>0"] <- sum (gsn_case.nv > 0) / sctrl.n
#	stats.nv["LM1_Not0Oddr"]    <- stats.nv["LM1_Case%_GsN>0"] / stats.nv["LM1_Ctrl%_GsN>0"]

	# comparison of logistic models
	GsN1.anova <- anova (full1.brglm, base1.brglm, test = "Chisq")
	dmodpv1.n <- GsN1.anova[["P(>|Chi|)"]][2]
	stats.nv["LM1_Dmod_Pv"] <- ifelse (is.na (dmodpv1.n), 1, dmodpv1.n)

	# prediction odds ratio (ratio between correct and incorrect predictions)
	lm1full_p1r1.n <- sum (full1.brglm$linear.predictors[full1.brglm$y == 1] > 0)
	lm1full_p0r1.n <- sum (full1.brglm$linear.predictors[full1.brglm$y == 1] < 0)
	lm1full_p0r0.n <- sum (full1.brglm$linear.predictors[full1.brglm$y == 0] < 0)
	lm1full_p1r0.n <- sum (full1.brglm$linear.predictors[full1.brglm$y == 0] > 0)
	stats.nv["LM1_PredOddr_full"] <- (lm1full_p1r1.n + lm1full_p0r0.n) / (lm1full_p1r0.n + lm1full_p0r1.n)
	lm1base_p1r1.n <- sum (base1.brglm$linear.predictors[base1.brglm$y == 1] > 0)
	lm1base_p0r1.n <- sum (base1.brglm$linear.predictors[base1.brglm$y == 1] < 0)
	lm1base_p0r0.n <- sum (base1.brglm$linear.predictors[base1.brglm$y == 0] < 0)
	lm1base_p1r0.n <- sum (base1.brglm$linear.predictors[base1.brglm$y == 0] > 0)
	stats.nv["LM1_PredOddr_base"] <- (lm1base_p1r1.n + lm1base_p0r0.n) / (lm1base_p1r0.n + lm1base_p0r1.n)

	# -- logistic model 2 --

	full2.brglm <- brglm (Class ~ LogLenTot + GsGene_N_Tot + Gs_N, data = data_gsi.df, family = binomial (link = "logit"))
	base2.brglm <- brglm (Class ~ LogLenTot + GsGene_N_Tot,        data = data_gsi.df, family = binomial (link = "logit"))

	full2.brgsm <- summary (full2.brglm)
	base2.brgsm <- summary (base2.brglm)
		
	gsn2_info.nv <- f.lrm_unit_getSmInfo (full2.brgsm, "Gs_N")
	stats.nv["LM2_GsN_Sign"] <- gsn2_info.nv["sign"]
	stats.nv["LM2_GsN_Pv"]   <- gsn2_info.nv["pvalue"]
	stats.nv["LM2_Convg"]    <- full2.brglm$converged

	GsN2.anova <- anova (full2.brglm, base2.brglm, test = "Chisq")
	dmodpv2.n <- GsN2.anova[["P(>|Chi|)"]][2]
	stats.nv["LM2_Dmod_Pv"] <- ifelse (is.na (dmodpv2.n), 1, dmodpv2.n)

	lm2full_p1r1.n <- sum (full2.brglm$linear.predictors[full2.brglm$y == 1] > 0)
	lm2full_p0r1.n <- sum (full2.brglm$linear.predictors[full2.brglm$y == 1] < 0)
	lm2full_p0r0.n <- sum (full2.brglm$linear.predictors[full2.brglm$y == 0] < 0)
	lm2full_p1r0.n <- sum (full2.brglm$linear.predictors[full2.brglm$y == 0] > 0)
	stats.nv["LM2_PredOddr_full"] <- (lm2full_p1r1.n + lm2full_p0r0.n) / (lm2full_p1r0.n + lm2full_p0r1.n)
	lm2base_p1r1.n <- sum (base2.brglm$linear.predictors[base2.brglm$y == 1] > 0)
	lm2base_p0r1.n <- sum (base2.brglm$linear.predictors[base2.brglm$y == 1] < 0)
	lm2base_p0r0.n <- sum (base2.brglm$linear.predictors[base2.brglm$y == 0] < 0)
	lm2base_p1r0.n <- sum (base2.brglm$linear.predictors[base2.brglm$y == 0] > 0)
	stats.nv["LM2_PredOddr_base"] <- (lm2base_p1r1.n + lm2base_p0r0.n) / (lm2base_p1r0.n + lm2base_p0r1.n)

	return (stats.nv)
	}

# this function just handles the exception 
# of missing x coefficient for fit problems
f.lrm_unit_getSmInfo <- function (main.sm, x.name)
	{
	output.nv <- numeric (2)
	names (output.nv) <- c ("sign", "pvalue")

	coef.mx <- main.sm$coefficients

	x.bn <- x.name %in% rownames (coef.mx)
	output.nv["sign"] <- ifelse (x.bn, sign (coef.mx[x.name, "Estimate"]), 0)
	output.nv["pvalue"] <- ifelse (x.bn, coef.mx[x.name, "Pr(>|z|)"], 1)
	
	return (output.nv)
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
