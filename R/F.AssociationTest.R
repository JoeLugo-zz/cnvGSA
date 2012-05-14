# -----------------------
# ASSOCIATION TEST + FDR
# -----------------------

F.AssociationTest <- function (cnvData.ls, gsData.ls, burdenSample.ls, assTestPar.ls, addEnrPar.ls)
	{
	# Different methods to generate totals and subset the count matrix
	# - 'all': 
	#   . totals: all samples in study, 
	#   . table: all samples with at least a cnv (don't subset further)
	# - 'cnv': 
	#   . totals: all samples with at least a cnv, 
	#   . table: all samples with at least a cnv (don't subset further)
	# - 'cnvGen': 
	#   . totals: all samples with at least a genic cnv, 
	#   . table: all samples with at least a genic cnv (requires subsetting)
	if (is.null (assTestPar.ls$iter)) assTestPar.ls$iter <- 2000

	methods.chv <- c ("all", "cnv", "cnvGen")
	
	f.make_test_input.ls <- list ()
	f.make_test_input.ls[[1]] <- f.make_test_input_all
	f.make_test_input.ls[[2]] <- f.make_test_input_cnv
	f.make_test_input.ls[[3]] <- f.make_test_input_cnvGen
	names (f.make_test_input.ls) <- methods.chv
	
	f.check_testpar (cnvData.ls, assTestPar.ls, methods.chv)
	
	type.ch <- assTestPar.ls$test_type
	input.ls <- f.make_test_input.ls[[type.ch]](cnvData.ls, gsData.ls, assTestPar.ls)
	
	stats.mx <- f.fet (input.ls)
	
	fdr.nv <- f.fdr (stats.mx[, "FET_pv"], input.ls, assTestPar.ls$iter)

	gs.id <- rownames (stats.mx)
	sizes.nv <- sapply (gsData.ls$gs2gene[gs.id], length)

	enr1.df <- data.frame (
				GsID = gs.id,
				GsName = gsData.ls$gs2name[gs.id],
				GsSize = sizes.nv,
				stringsAsFactors = F)
	enr2.df <- as.data.frame (stats.mx)
	enr3.df <- data.frame (FET_permFDR = fdr.nv)
	enr.df  <- cbind (enr1.df, enr2.df, enr3.df)
	
	# Logistic regression
	if( !is.null(addEnrPar.ls$do_logistic) ) {
		if( addEnrPar.ls$do_logistic == "full" ) {
			cat( "Logistic regression (for *all* gene-sets)..." )
			enr.df <- f.add_lrmstats( enr.df, cnvData.ls, burdenSample.ls, addEnrPar.ls )
			cat( "done\n" )
		}
	}
	
	enrRes.ls <- list (basic = enr.df, totals = input.ls$totals)
	
	return (enrRes.ls)
	}

f.check_testpar <- function (cnvData.ls, assTestPar.ls, methods.chv)
	{
	if (length (setdiff (assTestPar.ls$test_type, methods.chv)))
		{stop ("the association test method provided by 'assTestPar.ls$test_type' does not match with defined methods")}
	if (length (assTestPar.ls$test_type) > 1)
		{stop ("'assTestPar.ls$test_type' has more than one method")}
	if (length (setdiff (assTestPar.ls$test_classes, cnvData.ls$uni$class)))
		{stop ("the classes provided by 'assTestPar.ls$test_classes' do not match with classes in 'cnvData.ls'")}
	if (length (assTestPar.ls$test_classes) != 2)
		{stop ("'assTestPar.ls$test_classes' does not have two classes as required")}
	return ()
	}

# the make functions make the input ready for the test
# using different processing criteria (see 'F.AssociationTest' comments)
f.make_test_input_all <- function (cnvData.ls, gsData.ls, assTestPar.ls) 
	{
	data.ls <- list ()

	data.ls$tab <- cnvData.ls$tab$bin

	classes.chv <- assTestPar.ls$test_classes
	classes.ls  <- as.list (classes.chv)

	f.getSamples <- function (class.ch, s2class.df)
		{
		SampleID <- Class <- NULL	## workaround for "no visible binding for global variable" note in 'R CMD check' output (due to next line of code)
		return (unique (subset (s2class.df, select = SampleID, subset = Class == class.ch, drop = T)))
		}
	sample_byclass.ls <- lapply (classes.ls, f.getSamples, cnvData.ls$s2class)
	names (sample_byclass.ls) <- classes.chv
	
	data.ls$totals <- sapply (sample_byclass.ls, length)
	
	f.class_ix <- function (class.samples, count.tab)
		{return (which (rownames (count.tab) %in% class.samples))}
	data.ls$class_ix <- lapply (sample_byclass.ls, f.class_ix, data.ls$tab)
	
	return (data.ls)
	}
	
f.make_test_input_cnv <- function (cnvData.ls, gsData.ls, assTestPar.ls) 
	{
	data.ls <- list ()

	data.ls$tab <- cnvData.ls$tab$bin
	
	classes.chv <- assTestPar.ls$test_classes
	classes.ls  <- as.list (classes.chv)

	SampleID <- Class <- NULL	## workaround for "no visible binding for global variable" note in 'R CMD check' output (due to next line of code)
	s2class.df <- subset (cnvData.ls$full, select = c (SampleID, Class))
	s2class.df <- s2class.df[! duplicated (s2class.df), ]

	f.getSamples <- function (class.ch, s2class.df)
		{
		SampleID <- Class <- NULL	## workaround for "no visible binding for global variable" note in 'R CMD check' output (due to next line of code)
		return (unique (subset (s2class.df, select = SampleID, subset = Class == class.ch, drop = T)))
		}
	sample_byclass.ls <- lapply (classes.ls, f.getSamples, s2class.df)
	names (sample_byclass.ls) <- classes.chv
	
	data.ls$totals <- sapply (sample_byclass.ls, length)
	
	f.class_ix <- function (class.samples, count.tab)
		{return (which (rownames (count.tab) %in% class.samples))}
	data.ls$class_ix <- lapply (sample_byclass.ls, f.class_ix, data.ls$tab)
	
	return (data.ls)
	}

f.make_test_input_cnvGen <- function (cnvData.ls, gsData.ls, assTestPar.ls) 
	{
	data.ls <- list ()

	classes.chv <- assTestPar.ls$test_classes
	classes.ls  <- as.list (classes.chv)

	Gcount <- SampleID <- Class <- NULL	## workaround for "no visible binding for global variable" note in 'R CMD check' output (due to next line of code)
	s2class.df <- subset (cnvData.ls$full, subset = Gcount > 0, select = c (SampleID, Class))
	s2class.df <- s2class.df[! duplicated (s2class.df), ]

	data.ls$tab <- cnvData.ls$tab$bin[s2class.df$SampleID, ]

	f.getSamples <- function (class.ch, s2class.df)
		{
		SampleID <- Class <- NULL	## workaround for "no visible binding for global variable" note in 'R CMD check' output (due to next line of code)
		return (unique (subset (s2class.df, select = SampleID, subset = Class == class.ch, drop = T)))
		}
	sample_byclass.ls <- lapply (classes.ls, f.getSamples, s2class.df)
	names (sample_byclass.ls) <- classes.chv
	
	data.ls$totals <- sapply (sample_byclass.ls, length)
	
	f.class_ix <- function (class.samples, count.tab)
		{return (which (rownames (count.tab) %in% class.samples))}
	data.ls$class_ix <- lapply (sample_byclass.ls, f.class_ix, data.ls$tab)
	
	return (data.ls)
	}

# calls the test on each gene-set
f.fet <- function (input.ls)
	{
	gs_ix.ls <- as.list (1: ncol (input.ls$tab))

	stats.mx <- t (sapply (gs_ix.ls, f.fet_unit, input.ls))
	stats.mx[,12] <- p.adjust( stats.mx[,5], method = "BH" )	## Benjamini-Hochberg corrected p-value

	rownames (stats.mx) <- colnames (input.ls$tab)
	stats.mx <- stats.mx[order (stats.mx[, "FET_pv"], decreasing = F), ]
	
	return (stats.mx)
	}

# runs the test on a single gene-set
f.fet_unit <- function (gs.ix, input.ls)
	{
	c1_gsy.n <- sum (input.ls$tab[input.ls$class_ix[[1]], gs.ix])
	c1_gsn.n <- input.ls$totals[1] - c1_gsy.n 
	c2_gsy.n <- sum (input.ls$tab[input.ls$class_ix[[2]], gs.ix])
	c2_gsn.n <- input.ls$totals[2] - c2_gsy.n
	contingency.mx <- matrix (c (c1_gsy.n, c1_gsn.n, c2_gsy.n, c2_gsn.n), ncol = 2, nrow = 2, byrow = T)

	fet.test <- fisher.test (contingency.mx, alternative = "greater")
	fet.twosided <- fisher.test (contingency.mx, alternative = "two.sided")

	output.nv <- c (c1_gsy.n,
					c2_gsy.n,
					c1_gsy.n / input.ls$totals[1] * 100,
					c2_gsy.n / input.ls$totals[2] * 100,
					fet.test$p.value,
					fet.test$estimate,
					fet.test$conf.int,
					fet.twosided$estimate,
					fet.twosided$conf.int,
					0	## Placeholder for Benjamini-Hochberg corrected p-value
						## (see f.fet() for its assignment and n.b. the col numbers
						## in the p.adjust() line below)
					)

	names (output.nv) <- c (
						paste (names (input.ls$totals), "_N", sep = ""), 
						paste (names (input.ls$totals), "_%", sep = ""), 
						"FET_pv",
						"FET_OR",
						"FET_ORconfLow",
						"FET_ORconfHigh",
						"FET2s_OR",
						"FET2s_ORconfLow",
						"FET2s_ORconfHigh",
						"FET_bhFDR"
						)

	return (output.nv)
	}

# calls the fdr estimation on each iteration
f.fdr <- function (pv_real.nv, input.ls, iter.n) 
	{
	cat ("FDR Estimation...")
	
	fdr.nv <- list()
	
	if( iter.n == 0 )
	{
		fdr.nv <- rep( NA, length(pv_real.nv) )

		cat("SKIPPED\n")
	}
	else
	{
		iter.ls <- as.list (1: iter.n)
		# iterations at columns
		pv_rand.mx <- sapply (iter.ls, f.fdr_unit, input.ls)
		
		f.count <- function (pv_real.n, pv_real.nv, pv_rand.mx, iter.n)
			{
			op.n <- sum (pv_real.nv <= pv_real.n)
			fp.n <- sum (pv_rand.mx <= pv_real.n) / iter.n
			fdr.n <- fp.n / op.n 
			return (fdr.n)
			}
		fdr.nv <- sapply (as.list (pv_real.nv), f.count, pv_real.nv, pv_rand.mx, iter.n)
		
		cat ("done\n")
	}
	
	return (fdr.nv)
	}

# calls the test for a specific iteration
# corresponding to a single permutation of class indexes
f.fdr_unit <- function (i.n, input.ls) 
	{
	if (i.n %% 25 == 0) {cat (i.n); cat ("; ")}
	input_perm.ls <- input.ls
	
	c1.ix <- input.ls$class_ix[[1]]
	c2.ix <- input.ls$class_ix[[2]]	
	all.ix <- c (c1.ix, c2.ix)
	
	input_perm.ls$class_ix[[1]] <- sample (all.ix, size = length (c1.ix))
	input_perm.ls$class_ix[[2]] <- setdiff (all.ix, input_perm.ls$class_ix[[1]])

	pv_rand.nv <- f.fet (input_perm.ls)[, "FET_pv"]
	
	return (pv_rand.nv)
	}

