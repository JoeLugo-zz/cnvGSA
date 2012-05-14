
F.ProcessGenes <- function (cnvData.ls, gsData.ls, geneData.ls, assTestPar.ls)
	{
	f.check_genedata (geneData.ls)
		
	# Compute gene counts and association p-value
	gcounts_and_totals.ls <- f.test_genes (cnvData.ls, assTestPar.ls)

	gcounts.df <- gcounts_and_totals.ls$gcounts		# Gene counts
	totals.ls <- gcounts_and_totals.ls$totals		# Counts of total unique genes by class and CNV type
	
	# Compute support sets
	col.names <- paste (assTestPar.ls$test_classes, "_%", sep = "")
	sel1.ix <- which (gcounts.df[, col.names[1]] > gcounts.df[, col.names[2]])
	sel2.ix <- which (gcounts.df[, col.names[1]] < gcounts.df[, col.names[2]])
	support.genes.1 <- gcounts.df$GeneID[sel1.ix]
	support.genes.2 <- gcounts.df$GeneID[sel2.ix]
	support.genes.1 <- setdiff (support.genes.1, cnvData.ls$rem_genes)
	support.genes.2 <- setdiff (support.genes.2, cnvData.ls$rem_genes)
	
	geneData.ls$gcounts <- gcounts.df
	geneData.ls[[ paste( sep="", "support_", assTestPar.ls$test_classes[1] ) ]] <- support.genes.1
	geneData.ls[[ paste( sep="", "support_", assTestPar.ls$test_classes[2] ) ]] <- support.genes.2
	geneData.ls$totals <- totals.ls
		
	# Compute gene coverage statistics
	geneData.ls$coverage <- f.gene_coverage_stats (cnvData.ls, gsData.ls, geneData.ls)
	
	# add annotations to $gcounts
	geneData.ls <- f.add_ann (geneData.ls)

	return (geneData.ls)
	}

f.check_genedata <- function (geneData.ls)
	{
	if (is.null (geneData.ls))
		{stop ("Gene annotation data missing, add it to 'geneData.ls$ann'")}
	if (is.null (geneData.ls$ann))
		{stop ("Gene annotation data missing, add it to 'geneData.ls$ann'")}
	return ()
	}
	
f.test_genes <- function (cnvData.ls, assTestPar.ls)
	{
	Gcount <- GsID <- NULL	## workaround for "no visible binding for global variable" note in 'R CMD check' output (due to next line of code)
	cnv.df <- subset (cnvData.ls$full, subset = Gcount > 0, select = - GsID)
	cnv.df <- cnv.df[! duplicated (cnv.df), ]
	# refactor to exlude genes that were lost after filtering
	cnv.df$Genes <- factor (cnv.df$Genes)

	# sample counts at the gene level, for all genes
	Genes <- Class <- NULL	## workaround for "no visible binding for global variable" note in 'R CMD check' output (due to next line of code)
	gcount.tab <- table (subset (cnv.df, select = c (Genes, Class)))
	gcount.tab <- gcount.tab[, assTestPar.ls$test_classes]

	# FET can be computed using three methods for totals, as for gene-sets:
	# - all samples (s2class)
	# - all samples with at least a cnv of filtered type
	# - all samples with at least a genic cnv of filtered type
	totals.ls <- list ()
	totals.ls[["all"]] <- f.make_gene_totals (
								cnvData.ls$s2class, 
								assTestPar.ls$test_classes)
	SampleID <- Class <- NULL	## workaround for "no visible binding for global variable" note in 'R CMD check' output (due to next lines of code)
	totals.ls[["cnv"]] <- f.make_gene_totals (
								subset (cnvData.ls$full, select = c (SampleID, Class)), 
								assTestPar.ls$test_classes)
	Gcount <- SampleID <- Class <- NULL	## workaround for "no visible binding for global variable" note in 'R CMD check' output (due to next lines of code)
	totals.ls[["cnvGen"]] <- f.make_gene_totals (
								subset (cnvData.ls$full, subset = Gcount > 0, select = c (SampleID, Class)), 
								assTestPar.ls$test_classes)
	
	totals.nv <- totals.ls[[assTestPar.ls$test_type]]
	pvalue.nv <- apply (gcount.tab, 1, f.fet_gene_unit, totals.nv)
	
	r.n <- length (pvalue.nv)
	gcount1.df <- data.frame (GeneID = rownames (gcount.tab))
	gcount2.df <- as.data.frame (gcount.tab)
	gcount3.df <- as.data.frame (gcount.tab / c (rep (totals.nv[1], r.n), rep (totals.nv[2], r.n))) * 100
	gcount3.df <- round (gcount3.df, digits = 3)
	gcount4.df <- data.frame (Pvalue = pvalue.nv)
	gcount.df <- cbind (gcount1.df, gcount2.df, gcount3.df, gcount4.df)
	colnames (gcount.df)[2: 3] <- paste (assTestPar.ls$test_classes, "_N", sep = "")
	colnames (gcount.df)[4: 5] <- paste (assTestPar.ls$test_classes, "_%", sep = "")
	gcount.df <- gcount.df[order (gcount.df$Pvalue, decreasing = F), ]

	result <- list (gcounts = gcount.df, totals = totals.ls)
	
	return (result)
	}

f.make_gene_totals <- function (s2class.df, class.chv)
	{
	s2class.df <- s2class.df[! duplicated (s2class.df), ]
	s2class.df$SampleID <- factor (s2class.df$SampleID)
	
	stat.df <- aggregate (SampleID ~ Class, s2class.df, length)
	totals.nv <- stat.df$SampleID[match (class.chv, stat.df$Class)]
	names (totals.nv) <- class.chv

	return (totals.nv)
	}

f.fet_gene_unit <- function (counts.nv, totals.nv)
	{
	contingency.mx <- matrix (c (counts.nv, totals.nv), ncol = 2, nrow = 2, byrow = T)
	fet.test <- fisher.test (contingency.mx, alternative = "greater")
	return (fet.test$p.value)
	}

f.gene_coverage_stats <- function (cnvData.ls, gsData.ls, geneData.ls)
	{
	gs_u.genes <- unique (unlist (gsData.ls$gs2gene))
	Gcount <- Genes <- NULL	## workaround for "no visible binding for global variable" note in 'R CMD check' output (due to next lines of code)
	cnv_u.genes <- unique (unlist (strsplit (subset (cnvData.ls$cnv, subset = Gcount > 0, select = Genes, drop = T), split = cnvData.ls$gsep)))
	cnvf_u.genes <- unique (subset (cnvData.ls$full, subset = Gcount > 0, select = Genes, drop = T))
	
	cov.names <- c (
					"Genes in gene-set universe",
					"Genes hit by CNV (before filters)",
					"Genes hit by CNV (after filters)",
					"Genes hit by CNV (before filters) in gene-sets",
					"Genes hit by CNV (before filters) in gene-sets (%)",
					"Genes hit by CNV (after filters) in gene-sets",
					"Genes hit by CNV (after filters) in gene-sets (%)",
					"Genes in support",
					"Genes in support and also in gene-sets",
					"Genes in support and also in gene-sets (%)"
					)
					
	coverage.nv <- numeric ()
	
	coverage.nv[1] <- length (gs_u.genes)
	coverage.nv[2] <- length (cnv_u.genes)	
	coverage.nv[3] <- length (cnvf_u.genes)
	coverage.nv[4] <- length (intersect (gs_u.genes, cnv_u.genes))
	coverage.nv[5] <- coverage.nv[4] / coverage.nv[2]
	coverage.nv[6] <- length (intersect (gs_u.genes, cnvf_u.genes))
	coverage.nv[7] <- coverage.nv[6] / coverage.nv[3]

	coverage.nv[8] <- length (geneData.ls$support)
	coverage.nv[9] <- length (intersect (gs_u.genes, geneData.ls$support))
	coverage.nv[10] <- coverage.nv[9] / coverage.nv[8]
	
	names (coverage.nv) <- cov.names

	coverage.chv <- coverage.nv
	count.ix <- c (1: 4, 6, 8: 9)
	prcnt.ix <- c (5, 7, 10)
	coverage.chv[count.ix] <- formatC (coverage.nv[count.ix],       format = "d")
	coverage.chv[prcnt.ix] <- formatC (coverage.nv[prcnt.ix] * 100, format = "f", digits = 2)
	
	return (coverage.chv)
	}

f.add_ann <- function (geneData.ls)
	{
	gcount.df <- geneData.ls$gcounts
	ann.df <- data.frame (
				Symbol = geneData.ls$ann$gene2sy[as.character (gcount.df$GeneID)], 
				Name   = geneData.ls$ann$gene2name[as.character (gcount.df$GeneID)],
				stringsAsFactors = F)
	c.n <- ncol (gcount.df)
	geneData.ls$gcounts <- cbind (gcount.df$GeneID, ann.df, gcount.df[, 2: c.n])
	colnames (geneData.ls$gcounts)[1] <- "GeneID"

	return (geneData.ls)
	}
