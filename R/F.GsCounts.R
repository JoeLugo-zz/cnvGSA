# -----------------------
# GS COUNTS
# -----------------------
		
F.GsCounts <- function (cnvData.ls, gsData.ls)
	{
	# Check that for the samples mapping to gs there are enough classes and samples
	f.check_cnv_gs_preTest (cnvData.ls)  

	# Table where every cell has the number of perturbed genes for the (gs; sample) pair
	cnvData.ls$tab$gen <- f.make_table_gen (cnv.df = cnvData.ls$full)

	# Table where every cell has the number of distinct perturbing CNVs for the (gs; sample) pair
	cnvData.ls$tab$cnv <- f.make_table_cnv (cnv.df = cnvData.ls$full)
	
	# Table where every cell has the number of perturbing CNVs from distinct chromosomes for the (gs; sample) pair
	cnvData.ls$tab$chr <- f.make_table_chr (cnv.df = cnvData.ls$full)
	
	# Table where every cell has the binary perturbation score
	cnvData.ls$tab$bin <- f.make_table_bin (cnv.df = cnvData.ls$full) 
	
	# return (cnvData.ls)
	
	burdGsStat.ls <- f.burden_gs_summaries (cnvData.ls, gsData.ls)
	return (list (cnvData = cnvData.ls, burdGsStat = burdGsStat.ls))
	}

f.check_cnv_gs_preTest <- function (cnvData.ls)
	{
	if (is.null (cnvData.ls$full))
	{stop ("'$full' is missing from 'cnvData.ls': check if the pre-processing function has been correctly executed")}
	GsID <- NULL	# workaround for "no visible binding or global variable" note in 'R CMD check' output due to next line of code
	cnv2gs.df <- subset (cnvData.ls$full, subset = ! is.na (GsID))
	if (length (setdiff (cnvData.ls$uni$class, cnv2gs.df$Class)))
		{stop ("when restricting to CNVs mapped to gene-sets, one of the classes is not hit by CNV any more")}
	if (length (unique (cnv2gs.df$SampleID)) < 2)
		{stop ("when restricting to CNVs mapped to gene-sets, less than two samples are hit by CNV")}
	if (length (unique (cnv2gs.df$SampleID)) < 10)
		{message ("Warning: when restricting to CNVs mapped to gene-sets, less than 10 samples are hit by CNV")}
	return ()
	}	

# Comments for all 'f.make_...'
# s2gs.tab: perturbation score for each sample x gs pair, computed according to different criteria
# will be included in the table:
# - all samples with at least a cnv (whether genic or not)
# - all gene-sets hit by at least one cnv
# will therefore be excluded from the table:
# - samples without cnv
# - gene-sets without cnv
# * $full $SampleID and $GsID column are refactored to change the level universes

f.make_table_gen <- function (cnv.df)
	{
	cnv.df$SampleID <- factor (cnv.df$SampleID)
	cnv.df$GsID <- factor (cnv.df$GsID)

	s2gs.tab <- table (cnv.df[, c ("SampleID", "GsID")])
	return (s2gs.tab)
	}	

f.make_table_cnv <- function (cnv.df)
	{
	cnv.df$SampleID <- factor (cnv.df$SampleID)
	cnv.df$GsID <- factor (cnv.df$GsID)

	# By (i) subsetting to SampleID, GsID and CnvID and (ii) pruning redundancies
	# multiple counts will be generated only when multiple CNVs hit the same gene-set in the same sample
	SampleID <- CnvID <- GsID <- NULL	# workaround for "no visible binding for global variable" note in 'R CMD check' output due to next line of code
	cnv.df <- subset (cnv.df, select = c (SampleID, CnvID, GsID))
	cnv.df <- cnv.df[! duplicated (cnv.df), ]

	s2gs.tab <- table (cnv.df[, c ("SampleID", "GsID")])
	return (s2gs.tab)	
	}

f.make_table_chr <- function (cnv.df)
	{
	cnv.df$SampleID <- factor (cnv.df$SampleID)
	cnv.df$GsID <- factor (cnv.df$GsID)

	# By (i) subsetting to SampleID, GsID and Chr and (ii) pruning redundancies
	# multiple counts will be generated only when multiple CNVs from distinct chromosomes hit the same gene-set in the same sample
	SampleID <- Chr <- GsID <- NULL	# workaround for "no visible binding for global variable" note in 'R CMD check' output due to next line of code
	cnv.df <- subset (cnv.df, select = c (SampleID, Chr, GsID))
	cnv.df <- cnv.df[! duplicated (cnv.df), ]

	s2gs.tab <- table (cnv.df[, c ("SampleID", "GsID")])
	return (s2gs.tab)	
	}

f.make_table_bin <- function (cnv.df)
	{
	cnv.df$SampleID <- factor (cnv.df$SampleID)
	cnv.df$GsID <- factor (cnv.df$GsID)

	# By (i) subsetting to SampleID, GsID and (ii) pruning redundancies
	# multiple counts will not be generated
	SampleID <- GsID <- NULL	# workaround for "no visible binding for global variable" note in 'R CMD check' output due to next line of code
	cnv.df <- subset (cnv.df, select = c (SampleID, GsID))
	cnv.df <- cnv.df[! duplicated (cnv.df), ]

	s2gs.tab <- table (cnv.df[, c ("SampleID", "GsID")])
	return (s2gs.tab)
	}

# This is an equivalent formulation
# we preferred the formulation based on 'cnv.df' input
# only for sake of elegance and simmetry with the other function
# and also because it's independent of other ways of generating the count table
# * x.tab = s2gs_gen.tab
#f.make_table_bin <- function (x.tab)	
#	{
#	bin.tab <- x.tab
#	bin.tab[bin.tab > 1] <- 1
#	return (bin.tab)
#	}

f.burden_gs_summaries <- function (cnvData.ls, gsData.ls)
	{

	# Sample subsets to compute the statistics for
	# (all, cases, controls)
	sampleset.ls <- list ()
	sampleset.ls[[1]] <- cnvData.ls$s2class$SampleID
	Class <- SampleID <- NULL	# workaround for "no visible binding for global variable" note in 'R CMD check' output due to next line of code
	sampleset.ls[[2]] <- subset (cnvData.ls$s2class, subset = Class == cnvData.ls$uni$class[1], select = SampleID, drop = T)
	sampleset.ls[[3]] <- subset (cnvData.ls$s2class, subset = Class == cnvData.ls$uni$class[2], select = SampleID, drop = T)
	names (sampleset.ls) <- c ("All", cnvData.ls$uni$class)
			
	# Coverage summaries
	# (number of samples and gene-sets hit by CNV)
	coverage_byS.mx  <- sapply (sampleset.ls, f.burden_gs_coverage_bysample_unit, cnvData.ls)
	coverage_byGS.mx <- sapply (sampleset.ls, f.burden_gs_coverage_bygset_unit,   cnvData.ls, gsData.ls) 
	coverage.mx <- rbind (coverage_byS.mx, coverage_byGS.mx)

	# Count table summaries 
	# (perturbation summaries by gene-set x sample pairs)
	gs.n <- length (gsData.ls$gs2gene)
	pairs.mx <- sapply (sampleset.ls, f.burden_gs_pairs_unit, cnvData.ls, gs.n)
		
	burdGsStat.ls <- list (coverage = coverage.mx, pairs = pairs.mx)
	
	return (burdGsStat.ls)
	}

f.burden_gs_coverage_bysample_unit <- function (set.samples, cnvData.ls)
	{
	coverage.nv <- numeric ()	
	coverage.nv[1] <- length (set.samples)
	coverage.nv[2] <- length (intersect (set.samples, cnvData.ls$cnv$SampleID))
	coverage.nv[3] <- length (intersect (set.samples, cnvData.ls$full$SampleID))
	coverage.nv[4] <- coverage.nv[3] / coverage.nv[1]
	Gcount <- SampleID <- NULL	# workaround for "no visible binding for global variable" note in 'R CMD check' output due to next line of code of code
	coverage.nv[5] <- length (intersect (set.samples, subset (cnvData.ls$full, subset = Gcount > 0, select = SampleID, drop = T)))
	coverage.nv[6] <- coverage.nv[5] / coverage.nv[1]
	GsID <- SampleID <- NULL	# workaround for "no visible binding for global variable" note in 'R CMD check' output due to next line of code
	coverage.nv[7] <- length (intersect (set.samples, subset (cnvData.ls$full, subset = ! is.na (GsID), select = SampleID, drop = T)))
	coverage.nv[8] <- coverage.nv[7] / coverage.nv[1]

	names (coverage.nv) <- c (
								"Sample N in the study, no filters",
								"Sample N with at least one cnv, no filters", 
								"Sample N with at least one cnv", 
								"Sample % with at least one cnv (on tot)", 
								"Sample N with at least one genic cnv", 
								"Sample % with at least one genic cnv (on tot)", 
								"Sample N with at least one perturbed gene-set",
								"Sample % with at least one perturbed gene-set (on tot)"
								)	

	## coverage.chv <- coverage.nv
	## count.ix <- c (1: 3, seq (from = 5, to = 7, by = 2))
	## prcnt.ix <- seq (from = 4, to = 8, by = 2)

	## coverage.chv[count.ix] <- formatC (coverage.nv[count.ix],       format = "d")
	## coverage.chv[prcnt.ix] <- formatC (coverage.nv[prcnt.ix] * 100, format = "f", digits = 2)

	## return (coverage.chv)

	return (coverage.nv)
	}
	
f.burden_gs_coverage_bygset_unit <- function (set.samples, cnvData.ls, gsData.ls)
	{
	gs_tot.n <- length (gsData.ls$gs2gene)	

	coverage.nv <- numeric ()
	SampleID <- GsID <- NULL	# workaround for "no visible binding for global variable" note in 'R CMD check' output due to next line of code
	coverage.nv[1] <- length (unique (subset (cnvData.ls$full, subset = SampleID %in% set.samples, select = GsID, drop = T)))
	coverage.nv[2] <- coverage.nv[1] / gs_tot.n	
	
	names (coverage.nv) <- c (
							"Gene-set N with at least one sample",
							"Gene-set % with at least one sample"
							)
							
	## coverage.chv <- coverage.nv
	## coverage.chv[1] <- formatC (coverage.nv[1],       format = "d")
	## coverage.chv[2] <- formatC (coverage.nv[2] * 100, format = "f", digits = 2)

	## return (coverage.chv)
	
	return (coverage.nv)
	}
	
f.burden_gs_pairs_unit <- function (set.samples, cnvData.ls, gs.n)
	{
	pairs.nv <- numeric ()
	pairs.nv[1] <- sum (cnvData.ls$tab$bin[rownames (cnvData.ls$tab$bin) %in% set.samples])
	# not all samples and all gene-sets are in the tables
	pairs.nv[2] <- pairs.nv[1] / (gs.n * length (set.samples))
	pairs.nv[3] <- sum (cnvData.ls$tab$gen[rownames (cnvData.ls$tab$gen) %in% set.samples] >= 2)
	pairs.nv[4] <- pairs.nv[3] / pairs.nv[1]
	pairs.nv[5] <- sum (cnvData.ls$tab$gen[rownames (cnvData.ls$tab$gen) %in% set.samples] >= 3)
	pairs.nv[6] <- pairs.nv[5] / pairs.nv[1]
	pairs.nv[7] <- sum (cnvData.ls$tab$cnv[rownames (cnvData.ls$tab$cnv) %in% set.samples] >= 2)
	pairs.nv[8] <- pairs.nv[7] / pairs.nv[1]
	pairs.nv[9] <- sum (cnvData.ls$tab$chr[rownames (cnvData.ls$tab$chr) %in% set.samples] >= 2)
	pairs.nv[10] <- pairs.nv[9] / pairs.nv[1]
	
	names (pairs.nv) <- c (
							"N of sample-gs pair, >= 1 CNV-perturbed gene",
							"% of sample-gs pair, >= 1 CNV-perturbed gene (on all pairs)",
							"N of sample-gs pair, >= 2 CNV-perturbed gene",
							"% of sample-gs pair, >= 2 CNV-perturbed gene (on positive pairs)",
							"N of sample-gs pair, >= 3 CNV-perturbed gene",
							"% of sample-gs pair, >= 3 CNV-perturbed gene (on positive pairs)",
							"N of sample-gs pair, >= 2 CNV",
							"% of sample-gs pair, >= 2 CNV (on positive pairs)",
							"N of sample-gs pair, >= 2 CNV on distinct chr.",
							"% of sample-gs pair, >= 2 CNV on distinct chr. (on positive pairs)"
							)
	
	## pairs.chv <- pairs.nv
	## count.ix <- seq (from = 1, to = 9, by = 2)
	## prcnt.ix <- seq (from = 2, to = 10, by = 2)
	## pairs.chv[count.ix] <- formatC (pairs.nv[count.ix],       format = "d")
	## pairs.chv[prcnt.ix] <- formatC (pairs.nv[prcnt.ix] * 100, format = "f", digits = 2)
	
	## return (pairs.chv)

	return (pairs.nv)
	}
