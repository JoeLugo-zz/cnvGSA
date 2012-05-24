# -----------------------
# BURDEN (SAMPLE, CLASS)
# -----------------------

F.BurdenSample <- function (cnvData.ls, assTestPar.ls)	
	{
	# Initialize output object
	burdSampleStat.ls <- list ()
	
	# Do checks
	f.check_cnv_postproc (cnvData.ls)

	# Two reports:
	# 1.
	# - statistics using all samples with at least one CNV
	# - proportions of samples without any CNV
	# 2.
	# - statistics using all samples with at least one genic CNV
	# - proportions of samples with at least one CNV but no genic ones

	## Which CNV map to gene-set genes?
	## cnvData.ls$full$GsGeneBn <- 0
	## cnvData.ls$full$GsGeneBn[! is.na (cnvData.ls$full$GsID)] <- 1

	# REPORT_1

	# Export filtered cnv data from $full
	Genes <- GsID <- NULL	# workaround for "no visible binding for global variable" note in 'R CMD check' output due to next line of code
	cnv.df <- subset (cnvData.ls$full, select = c (- Genes, - GsID))
	cnv.df <- cnv.df[!duplicated (cnv.df), ]

	# Separate cnv data by class
	f.subset <- function (class.ch, cnv.df)
		{
		Class <- NULL	# workaround for "no visible binding for global variable" note in 'R CMD check' output due to next line of code
		return (subset (cnv.df, Class == class.ch))
		}
	cnv.ls <- lapply (as.list (cnvData.ls$uni$class), f.subset, cnv.df)
	names (cnv.ls) <- cnvData.ls$uni$class

	# Make stats by sample, for the two classes separately
	stat.ls <- lapply (cnv.ls, f.burden_sample_stats)
	burdSampleStat.ls$stat.ls <- stat.ls	## Needed for f.add_lrm (the logistic test); will be NULLed afterward

	# Make summaries of the stats by sample, for the two classes separately
	summary.ls <- lapply (stat.ls, f.burden_sample_summaries)
	# names (summary.ls) <- names (stat.ls)
	
	# Compute t-test p-values for the stats-by-sample, one class vs the other
	# and make boxplot PDFs (if boxplot.bn = T)
	pv.mx <- f.burden_sample_test (stat.ls, boxplot.bn = assTestPar.ls$boxplot.bn)

	# Compute the proprtions of samples without CNV by class, and test their difference
	prop.ls <- f.burden_sample_prop_cnv (cnv.ls, cnvData.ls)

	burdSampleStat.ls$SamplesCNV <- list (summary = summary.ls, pvalue = pv.mx, no_cnv_proportion = prop.ls)

#
# RZ 2012-04-18: The following was commented out ahead of the original submission to Bioconductor (i.e. cnvGSA_1.0.0)
#
#	# REPORT_2
#
#	# Export filtered cnv and sample data from $full
#	cnv.df <- subset (cnvData.ls$full, subset = Gcount > 0, select = c (- Genes, - GsID))
#	cnv.df <- cnv.df[!duplicated (cnv.df), ]
#
#	# Separate cnv data by class
#	cnv.ls <- lapply (as.list (cnvData.ls$uni$class), f.subset, cnv.df)
#	names (cnv.ls) <- cnvData.ls$uni$class
#	
#	# Make stats by sample, for the two classes separately
#	# * CNV_N and genCNV_N will be equal
#	stat.ls <- lapply (cnv.ls, f.burden_sample_stats)
#
#	# Make summaries of the stats by sample, for the two classes separately
#	summary.ls <- lapply (stat.ls, f.burden_sample_summaries)
#
#	# Compute t-test p-values for the stats-by-sample, one class vs the other
#	pv.mx <- f.burden_sample_test (stat.ls, boxplot.bn = F)
#
#	# Compute the proprtions of samples with CNV but no gneic ones by class, and test their difference
#	prop.ls <- f.burden_sample_prop_geniccnv (cnv.ls, cnvData.ls)
#	
#	burdSampleStat.ls$SamplesGenicCNV <- list (summary = summary.ls, pvalue = pv.mx, no_cnv_proportion = prop.ls)
	
	# RETURN
	
	return (burdSampleStat.ls)
	}

f.check_cnv_postproc <- function (cnvData.ls)
	{
	if (is.null (cnvData.ls$full))
		{stop ("'$full' is missing from 'cnvData.ls': check if the pre-processing function has been correctly executed")}
	Class <- NULL	# workaround for "no visible binding for global variable" note in 'R CMD check' output due to next line of code
	if (length (setdiff (cnvData.ls$uni$class, cnvData.ls$full$Class)))
		{stop ("one of the classes is not hit by CNV any more")}
	Class <- Gcount <- NULL	# workaround for "no visible binding for global variable" note in 'R CMD check' output due to next line of code
	if (length (setdiff (cnvData.ls$uni$class, subset (cnvData.ls$full, select = Class, subset = Gcount > 0, drop = T))))
		{stop ("one of the classes is not hit by genic CNV any more")}
	if (length (unique (cnvData.ls$full$SampleID)) < 2)
		{stop ("less than two samples are hit by CNV")}
	return ()
	}	
	
f.burden_sample_stats <- function (cnv.df)
	#
	# N.B.: f.add_lrm (the subroutine that does the logistic regression) depends on the output of this function
	#
	{
	# Statistics by Sample
	stat.ls <- list ()
		
	f.count_not0 <- function (x.nv)
		{return (sum (x.nv > 0))}

	logmean <- function(v) { return (log10 (mean (v))) }
	stat.ls[[1]] <- aggregate (Length ~ SampleID, cnv.df, logmean)
	logsum <- function(v) { return (log10 (sum (v))) }
	stat.ls[[2]] <- aggregate (Length ~ SampleID, cnv.df, logsum)
	
	stat.ls[[3]] <- aggregate (CnvID ~ SampleID, cnv.df, length)
	stat.ls[[4]] <- aggregate (Gcount   ~ SampleID, cnv.df, f.count_not0)
	stat.ls[[5]] <- aggregate (GsGcount ~ SampleID, cnv.df, f.count_not0)

	stat.ls[[6]] <- aggregate (Gcount   ~ SampleID, cnv.df, mean)
	stat.ls[[7]] <- aggregate (GsGcount ~ SampleID, cnv.df, mean)
	stat.ls[[8]] <- aggregate (Gcount   ~ SampleID, cnv.df, sum)
	stat.ls[[9]] <- aggregate (GsGcount ~ SampleID, cnv.df, sum)

	names (stat.ls) <- c (
		"LogLenMean",
		"LogLenTot", 

		"CNV_N",
		"GenCNV_N",
		"GsGenCNV_N",

		"Gene_N_Mean",
		"GsGene_N_Mean",
		"Gene_N_Tot",
		"GsGene_N_Tot"
	)

	f.setStatName <- function (stat.df)
		{names (stat.df)[2] <- "Stat"; return (stat.df)}
	stat.ls <- lapply (stat.ls, f.setStatName)

	### TODO: The output now looks like:
	### 
	### List of 2
	###  $ case:List of 9
	###   ..$ LogLenMean   :'data.frame':       692 obs. of  2 variables:
	###   .. ..$ SampleID: chr [1:692] "1020_4" "1030_3" "1045_3" "1050_3" ...
	###   .. ..$ Stat    : num [1:692] 4.55 5.27 4.51 4.61 5.46 ...
	###   ..$ LogLenTot    :'data.frame':       692 obs. of  2 variables:
	###   .. ..$ SampleID: chr [1:692] "1020_4" "1030_3" "1045_3" "1050_3" ...
	###   .. ..$ Stat    : num [1:692] 4.55 5.57 4.51 4.61 5.94 ...
	###   ..$ CNV_N        :'data.frame':       692 obs. of  2 variables:
	###   .. ..$ SampleID: chr [1:692] "1020_4" "1030_3" "1045_3" "1050_3" ...
	###   .. ..$ Stat    : int [1:692] 1 2 1 1 3 2 2 1 2 4 ...
	###   [...]
	###  $ ctrl:List of 9
	###   ..$ LogLenMean   :'data.frame':       880 obs. of  2 variables:
	###   .. ..$ SampleID: chr [1:880] "B100121_1007854727" "B100331_1007873991" "B101..
	###   .. ..$ Stat    : num [1:880] 5.12 4.77 4.88 4.86 4.72 ...
	###   ..$ LogLenTot    :'data.frame':       880 obs. of  2 variables:
	###   .. ..$ SampleID: chr [1:880] "B100121_1007854727" "B100331_1007873991" "B101..
	###   .. ..$ Stat    : num [1:880] 5.12 5.07 5.35 5.16 5.02 ...
	###   ..$ CNV_N        :'data.frame':       880 obs. of  2 variables:
	###   .. ..$ SampleID: chr [1:880] "B100121_1007854727" "B100331_1007873991" "B101..
	###   .. ..$ Stat    : int [1:880] 1 2 3 2 2 3 3 3 1 4 ...
	###   [...]
	###
	### ...
	###
	### Rebuild it so that $case and $ctrl are data frames with SampleID as rows
	### and the stats as cols, as it should be...
	
	return (stat.ls)
	}

# Distribution summaries of statistics by Sample
f.burden_sample_summaries <- function (stat.ls)
	{
	summaries.chv <- c ("Min", "Q1", "Mean", "Median", "Q3", "Max")
	stat.mx <- matrix (ncol = length (stat.ls), nrow = length (summaries.chv))
	colnames (stat.mx) <- names (stat.ls)
	rownames (stat.mx) <- summaries.chv
	
	for (j in 1: length (stat.ls))
		{
		stat.mx[1, j] <- min  (stat.ls[[j]]$Stat)
		stat.mx[2, j] <- quantile (stat.ls[[j]]$Stat, p = 0.25)
		stat.mx[3, j] <- mean (stat.ls[[j]]$Stat)
		stat.mx[4, j] <- median (stat.ls[[j]]$Stat)
		stat.mx[5, j] <- quantile (stat.ls[[j]]$Stat, p = 0.75)
		stat.mx[6, j] <- max  (stat.ls[[j]]$Stat)
		}
		
	return (stat.mx)
	}

f.burden_sample_test <- function (stat.ls, boxplot.bn)
	{
	l.n <- length (stat.ls[[1]])
	pv.mx <- matrix (nrow = length (stat.ls), ncol = l.n)
	colnames (pv.mx) <- names (stat.ls[[1]])
	rownames (pv.mx) <- c (
						paste (names (stat.ls)[1], ">", names (stat.ls)[2]),
						paste (names (stat.ls)[1], "<", names (stat.ls)[2])
						)
	for (j in 1: l.n)
		{
		x1.n <- stat.ls[[1]][[j]]$Stat
		x2.n <- stat.ls[[2]][[j]]$Stat
		data.df <- data.frame (
						x = c (x1.n, x2.n), 
						class = as.factor (c (rep (1, length (x1.n)), rep (2, length (x2.n))))
						)
							
		pv.mx[1, j] <- t.test (x ~ class, data.df, alternative = "greater")$p.value
		pv.mx[2, j] <- t.test (x ~ class, data.df, alternative = "less")$p.value
		
		# Print boxplot to pdf
		if (boxplot.bn)
			{
			stat.name <- colnames (pv.mx)[j]
			col.chv <- c ("brown", "darkblue")
			pdf (paste ("BurdenSample_", stat.name, ".pdf", sep = ""))
			boxplot (x1.n, x2.n, 
					border = col.chv, 
					names = names (stat.ls), main = stat.name, 
					lwd = 2, cex = 1.5, cex.axis = 1.5)
			points (x = c (1, 2), y = c (mean (x1.n), mean (x2.n)), 
					type = "p", pch = "_", cex = 10, col = col.chv)
			dev.off ()
			}
		}
	
	# * Note: the t-test is used instead of the z-test just for ease of input
	#         however, the two tests are supposed to have very little discrepancy for large sample size
		
	return (pv.mx)
	}

f.burden_sample_prop_cnv <- function (cnv.ls, cnvData.ls)
	{
	sel1.n <- length (unique (cnv.ls[[1]]$SampleID))
	sel2.n <- length (unique (cnv.ls[[2]]$SampleID))
	
	SampleID <- Class <- NULL	# workaround for "no visible binding for global variable" note in 'R CMD check' output due to next lines of code
	tot1.n <- length (unique (subset (cnvData.ls$s2class, select = SampleID, subset = Class == names (cnv.ls)[1], drop = T)))
	tot2.n <- length (unique (subset (cnvData.ls$s2class, select = SampleID, subset = Class == names (cnv.ls)[2], drop = T)))
	
	prop.nv <- c (sel1.n / tot1.n, sel2.n / tot2.n)
	names (prop.nv) <- names (cnv.ls)
	
	# browser ()
	
	pv.n <- prop.test (c (sel1.n, sel2.n), c (tot1.n, tot2.n))$p.value
	
	return (list (PropEstimates = prop.nv, Pvalue = pv.n))
	}

# cnv.ls must have been restricted to samples with genic cnvs
f.burden_sample_prop_geniccnv <- function (cnv.ls, cnvData.ls)
	{
	sel1.n <- length (unique (cnv.ls[[1]]$SampleID))
	sel2.n <- length (unique (cnv.ls[[2]]$SampleID))
	
	SampleID <- Class <- NULL	# workaround for "no visible binding for global variable" note in 'R CMD check' output due to next lines of code
	tot1.n <- length (unique (subset (cnvData.ls$full, select = SampleID, subset = Class == names (cnv.ls)[1], drop = T)))
	tot2.n <- length (unique (subset (cnvData.ls$full, select = SampleID, subset = Class == names (cnv.ls)[2], drop = T)))
	
	prop.nv <- c (sel1.n / tot1.n, sel2.n / tot2.n)
	names (prop.nv) <- names (cnv.ls)

	# browser ()
		
	pv.n <- prop.test (c (sel1.n, sel2.n), c (tot1.n, tot2.n))$p.value
	
	return (list (PropEstimates = prop.nv, Pvalue = pv.n))
	}
