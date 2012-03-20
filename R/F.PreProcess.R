# -----------------------
# PRE-PROCESS
# -----------------------

F.PreProcess <- function (cnvData.ls, gsData.ls)
	{
	# 0. COMPLETE DATA
	# Add lengths (CNV length, gene counts)
	gsgenes.chv <- unique (unlist (gsData.ls$gs2gene))
	cnvData.ls$cnv <- f.lengths_cnv.df (cnvData.ls, gsgenes.chv)
	
	# Add Universes
	cnvData.ls$uni <- f.add_universes (cnvData.ls)
	
	# 1. CHECKS
	f.check_cnv_s2class (cnvData.ls)
	f.check_gs (gsData.ls)
	
	# 2. MERGE ALL DATA
	cnv2g.ls  <- strsplit (cnvData.ls$cnv$Genes, split = cnvData.ls$gsep)
	names (cnv2g.ls) <- cnvData.ls$cnv$CnvID
	cnv2g.df <- stack (cnv2g.ls)
	names (cnv2g.df) <- c ("Genes", "CnvID")
	
	# CNVs without genes must be kept, hence 'all = T'
	# (mind that this will generate NA values)

	Genes <- NULL	## workaround for "no visible binding for global variable" note in 'R CMD check' output (due to next line of code)
	cnv.df <- subset (cnvData.ls$cnv, select = - Genes)
	m1.df <- merge (cnv.df, cnv2g.df, by = "CnvID", all = T)

	gs2g.df <- stack (gsData.ls$gs2gene)
	names (gs2g.df) <- c ("Genes", "GsID")

	# GSs without CNVs don't have to be kept, hence 'all.y = F'
	# (mind that this will avoid generating NA values in most of the columns)
	m2.df <- merge (m1.df, gs2g.df, by = "Genes", all.x = T, all.y = F)

	# Samples without CNVs don't have to be kept, hence 'all.y = F'
	m3.df <- merge (m2.df, cnvData.ls$s2class, by = "SampleID", all.x = T, all.y = F)

	# 3. FILTER	
	# Filter by CNV type and length
	cnvData.ls$filters <- f.init_filters (cnvData.ls)
	f.check_filters (cnvData.ls)
	f.check_genes (cnvData.ls, gsData.ls)

	# - cnvData.ls$full has all the data in one big data.frame
	cnvData.ls$full <- f.filter (m3.df, cnvData.ls)
	
	# Remove genes
	cnvData.ls$full <- f.rem_genes (cnvData.ls)
	
	return (cnvData.ls)
	}

f.lengths_cnv.df <- function (cnvData.ls, gsgenes.chv)
	{
	cnv.df <- cnvData.ls$cnv
	cnv.df$Length <- cnv.df$Coord_f - cnv.df$Coord_i

	genes.ls  <- strsplit (cnv.df$Genes, split = cnvData.ls$gsep)
	cnv.df$Gcount <- sapply (genes.ls, length)

	cnv.df$GsGcount <- sapply (lapply (genes.ls, intersect, y = gsgenes.chv), length)
	
	return (cnv.df)
	}

f.add_universes <- function (cnvData.ls)
	{
	uni.ls <- list ()
	
	uni.ls$class  <- unique (cnvData.ls$s2class$Class)
	uni.ls$sample <- unique (cnvData.ls$s2class$SampleID)
	uni.ls$type   <- unique (cnvData.ls$cnv$Type)
		
	return (uni.ls)
	}

f.check_cnv_s2class <- function (cnvData.ls)
	{
	if (length (cnvData.ls$uni$class) < 2)	
		{stop ("less than two classes in 'cnvData.ls$s2class'")}
	if (length (cnvData.ls$uni$class) > 2)	
		{stop ("more than two classes in 'cnvData.ls$s2class'")}
	if (sum (cnvData.ls$cnv$Length <= 0))
		{stop ("negative or zero value(s) in 'cnvData.ls$cnv$Length'; check '$Coord_i', '$Coord_f'")}
	if (length (setdiff (cnvData.ls$cnv$SampleID, cnvData.ls$s2class$SampleID)))
		{stop ("some of the samples in 'cnvData.ls$cnv' are not in 'cnvData.ls$s2class'")}
	}
	
f.check_gs <- function (gsData.ls)
	{
	if (length (setdiff (names (gsData.ls$gs2gene), names (gsData.ls$gs2name))))
		{stop ("gs.ls has gene-sets that are not in gs.chv")}
	return ()
	}
	
f.check_filters <- function (cnvData.ls)
	{
	if (length (setdiff (cnvData.ls$filters$limits_type, cnvData.ls$uni$type)))
		{message ("Warning: in 'cnvData': some values of '$filters$limits_type' do not match '$cnv$Type'")}
	if (length (setdiff (cnvData.ls$filters$limits_size$Type, cnvData.ls$cnv$Type)))
		{stop ("in 'cnvData': some values of '$filters$limits_size' do not match '$cnv$Type'")}
	return ()
	}

f.check_genes <- function (cnvData.ls, gsData.ls)
	{
	if (! sum (nchar (cnvData.ls$cnv$Genes) > 1))
		{stop ("in 'cnvData': no genes mapped to CNVs")}
	genes.ls <- strsplit (cnvData.ls$cnv$Genes, split = cnvData.ls$gsep)
	if (! sum (sapply (genes.ls, length) > 1))
		{message ("Warning: in 'cnvData': none of the CNVs has more than one mapped gene, check 'cnvData$gsep'")}
	cnv.genes <- unique (unlist (genes.ls))
	gs.genes  <- unique (unlist (gsData.ls$gs2gene))
	cnv_gs.genes <- intersect (cnv.genes, gs.genes)
	if (length (cnv.genes) < 50)
		{message ("Warning: in 'cnvData$cnv$Genes': less than 50 genes are hit by a CNV")}
	if (length (gs.genes) < 100)
		{message ("Warning: in 'gsData$gs2gene': less than 100 genes are present")}
	if (length (cnv_gs.genes) / length (cnv.genes) < 0.25)
		{message ("Warning: less than 25% of genes in 'cnvData$cnv$Genes' are mapped on 'gsData$gs2gene'")}
	rem.genes <- cnvData.ls$filters$rem_genes
	if (length (setdiff (rem.genes, cnv.genes)))
		{message ("Warning: some of the genes to be removed in 'cnvData.ls$filters$rem_genes' are missing from 'cnvData$cnv$Genes'")}
	return ()
	}
	
f.init_filters <- function (cnvData.ls)
	{
	limits.ls <- cnvData.ls$filters

	if (is.null (cnvData.ls$filters$limits_type))
		{cnvData.ls$filters$limits_type <- cnvData.ls$uni$type}

	if (is.null (cnvData.ls$filters$limits_size))
		{
		types.n <- length (cnvData.ls$uni$type)
		limsize.df <- data.frame (
					Type = cnvData.ls$uni$type, 
					Max_length = rep (Inf, types.n), 
					Max_gcount = rep (Inf, types.n)
					)	

		# Set max limits to 0 for types that need to be removed
		nsel.ix <- which (! limsize.df$Type %in% cnvData.ls$filters$limits_type)
		limsize.df$Max_length[nsel.ix] <- 0
		limsize.df$Max_gcount[nsel.ix] <- 0
		
		limits.ls$limits_size <- limsize.df
		}
	
	if (is.null (cnvData.ls$filters$rem_genes))
		{
		limits.ls$rem_genes <- character (0)
		}
	
	return (limits.ls)
	}
	
f.filter <- function (full.df, cnvData.ls)
	{
	limsize.df <- cnvData.ls$filters$limits_size
	
	sel.ix <- match (full.df$Type, limsize.df$Type)
	full.df$Max_length <- limsize.df$Max_length[sel.ix]
	full.df$Max_gcount <- limsize.df$Max_gcount[sel.ix]
	Length <- Max_length <- Gcount <- Max_gcount <- NULL	## workaround for "no visible binding for global variable" note in 'R CMD check' output (due to the next lines of code)
	f.df <- subset (
				full.df, 
				Length < Max_length & Gcount < Max_gcount
				)
	f.df <- subset (f.df, select = c (- Max_length, - Max_gcount))
	
	return (f.df)
	}	

f.rem_genes <- function (cnvData.ls)	
	{
	Genes <- CnvID <- NULL	## workaround for "no visible binding for global variable" note in 'R CMD check' output (due to next lines of code)
	rem.cnvid <- subset (cnvData.ls$full, 
					subset = Genes %in% cnvData.ls$filters$rem_genes, 
					select = CnvID, drop = T)
	f.df <- subset (cnvData.ls$full, ! CnvID %in% rem.cnvid)
		
	return (f.df)
	}
	
