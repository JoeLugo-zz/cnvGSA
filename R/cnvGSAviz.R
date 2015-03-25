# VIZ

#' Creating the plots from the CnvGSAOutput data
#'
#' @param cnvGSA.in A CnvGSAInput S4 object
#' @param cnvGSA.out A CnvGSAOutput S4 object
#' @return Creates the plots to better understand the output

# setwd ("/Users/josephlugo/Documents/R/PGC2_test/Results/ScriptP02_MakeLogregTable_v06_20141208h18m25_41k_size_kloci_gwas")
# load  ("wsRes_ScriptP02_MakeLogregTable_v06_20141208h18m25_41k_size_kloci_gwas.RData")

# * strange behavior of ARC / NMDAR: including known loci, ARC significant for losses and NMDAR for gains, when removing known loci, both significant for losses (but ARC more) -- need to investigate at gene level
# * synaptic sets: GO synaptic components, GO neuron projection components / morphogsListis and Darnell FMR1 target (from mouse) rank on top for losses, with as well as without known loci
# * synaptic pathway panel: glutamatargic first and significant, with and without known loci; dopaminergic borderline significant only when including known loci
# * neurophenotype: human autosomal dominant / X-linked, mouse abnormal behavior / neurological phenotype rank on top for losses, with as well as without known loci
# * mouse phenotype panel: brain and beahvior-specific pehnotype only significant after total count correction, for losses, with as well as without known loci -- need to investigate gains if any trend
# * brain expression: the set including very high and medium/high gsList without a strong pre-natal or post-natal bias is the most significant, the other are not (but need to keep in mind what the other represent)
# * brain expression: inverse correlation between loss prevalence and expression tier
# * constrain indexes show the expected correlation, with the top sets significant for losses (for GI, also wehn removing known loci); interestingly, GI index, not based on CNV data and predicting nonsynonymous intolerance, correlates better than predicted haploins.

# set groups:
# 1: synaptic and neurofunction

# frame-A-a: 100%: significance U correction
# frame-B-a: 45%: compare significance for different correction methods
# frame-C-l: 60%: U correction (rescaled), coefficient
# frame-C-r: 90%: support
# configPath = "/Users/josephlugo/Documents/R/PGC2_test/R_Works/";configFile = "vizConfig.txt";

# vizData <- f.vizPreProcess(cnvGSA.in,cnvGSA.out)
# setwd("/Users/josephlugo/Documents/R/PGC2_test/Restructure_Results/8.5/Viz/")
# load("vizInput2015031014h01m25.RData")
# f.makeViz(cnvGSA.in,cnvGSA.out)

f.makeViz <- function(cnvGSA.in,cnvGSA.out)
{
	t <- Sys.time()
	timestamp <- strftime(t,"%Y%m%d%Hh%Mm%S")

	vizData <- list(cnvGSA.in@gsData.ls$gs_all.ls,cnvGSA.out@res.ls,cnvGSA.in@cnvData.ls$cnv.df,cnvGSA.in@phData.ls$ph_TYPE.df)
	names(vizData) <- list("gs_all.ls","res.ls","cnv.df","ph_TYPE.df")

	setwd (cnvGSA.in@config.ls$outputPath)
	save(vizData,file=paste("vizInput_",timestamp,".RData",sep=""))

	config.df <- cnvGSA.in@config.ls$config.df

	Kl            <- config.df[config.df$param == "Kl","value"]
	FDRThreshold  <- as.numeric(config.df[config.df$param == "FDRThreshold","value"])
	gsList        <- config.df[config.df$param == "gsList","value"]
	cnvType       <- config.df[config.df$param == "cnvType","value"]
	plotHeight    <- as.numeric(config.df[config.df$param == "plotHeight","value"])
	outputPathViz <- config.df[config.df$param == "outputPathViz","value"]
	labelSize     <- as.numeric(config.df[config.df$param == "labelSize","value"])
	gsList        <- config.df[config.df$param == "gsList","value"]

	if (Kl == ""){Kl <- "ALL"}
	if (is.na(FDRThreshold)){FDRThreshold <- 0.1}
	if (is.na(plotHeight)){plotHeight <- 13}
	if (is.na(labelSize)){labelSize <- 0.7}

	if (gsList != ""){
			gsList <- unlist(strsplit(readLines(gsList),","))
		}
		else{
			stop("No genes to plot. Check gsList.")
		}

	gs_all.ls <- vizData$gs_all.ls
	cnv.df    <- vizData$cnv.df
	res.ls    <- vizData$res.ls

	if (length(gsList) == 0){
		stop("No gene-sets in gsList.")
	}

	gs_len.nv <- sapply (gs_all.ls, length)

	setwd (outputPathViz)
	cat(paste("Changing directory to ",outputPathViz,sep=""))
	cat("\n")

	# 1.
	# NEUROFUNCTION + SYNAPTIC
	z_set.gsid   <- gsList
	z_olp.n      <- length (z_set.gsid)
	z_olp.mx     <- matrix (data = NA, ncol = z_olp.n, nrow = z_olp.n, dimnames = list (z_set.gsid, z_set.gsid))
	# finding how many common gsList there are between the gene sets
	for (i in 1: z_olp.n)
		{
		for (j in 1: z_olp.n)
			{
			zi.gid <- gs_all.ls[[z_set.gsid[i]]]; zj.gid <- gs_all.ls[[z_set.gsid[j]]]
			z_olp.mx[i, j] <- length (intersect (zi.gid, zj.gid)) / length (zi.gid) * 100
			}
		}
	rm (i, j)

	z_olp.df <- as.data.frame(z_olp.mx)
	z_olp.df <- cbind(GeneSets = row.names(z_olp.df),z_olp.df)
	row.names(z_olp.df) <- NULL
	write.table (z_olp.df, col.names = T, row.names = F, sep = "\t", quote = F, file = paste("GsOverlap_",timestamp,".txt",sep=""))

	# cnvType: NEUROFUNCTION + SYNAPTIC
	z_set.gsid   <- gsList
	z_set.labels <- paste (z_set.gsid, gs_len.nv[z_set.gsid], sep = ": ")

	if (Kl == "ALL"){
		z_set1.col   <- rep ("gray30", length (z_set.gsid))
		z_set2.col   <- rep ("gray30", length (z_set.gsid))
		resObjectKly <- get(paste("covAll_chipAll_",cnvType,"_KLy.df",sep=""),res.ls)
		resObjectKln <- get(paste("covAll_chipAll_",cnvType,"_KLn.df",sep=""),res.ls)
		z_set1.df    <- resObjectKly[match (z_set.gsid, resObjectKly$GsID), ]
		z_set2.df    <- resObjectKln[match (z_set.gsid, resObjectKln$GsID), ]
		height_U.mx <- matrix (data = c (z_set1.df$Pvalue_U_dev_s, z_set2.df$Pvalue_U_dev_s), nrow = 2, byrow = T)
		height_Uc.mx <- matrix (data = c (z_set1.df$Coeff_U, z_set2.df$Coeff_U), nrow = 2, byrow = T)
		z_set1.col[which (z_set1.df$FDR_BH_U <= FDRThreshold)] <- "brown"
		z_set2.col[which (z_set2.df$FDR_BH_U <= FDRThreshold)] <- "brown" 
		border.vc <- c (z_set1.col, z_set2.col)
		nc.vc <- z_set1.col
	} else if (Kl == "YES"){
		z_set1.col   <- rep ("gray30", length (z_set.gsid))
		resObjectKly <- get(paste("covAll_chipAll_",cnvType,"_KLy.df",sep=""),res.ls)
		z_set1.df    <- resObjectKly[match (z_set.gsid, resObjectKly$GsID), ]
		height_U.mx <- matrix (data = c (z_set1.df$Pvalue_U_dev_s), nrow = 1, byrow = T)
		height_Uc.mx <- matrix (data = c (z_set1.df$Coeff_U), nrow = 1, byrow = T)
		z_set1.col[which (z_set1.df$FDR_BH_U <= FDRThreshold)] <- "brown"
		border.vc <- c (z_set1.col)
		nc.vc <- z_set1.col
	} else if (Kl == "NO"){
		z_set2.col   <- rep ("gray30", length (z_set.gsid))
		resObjectKln <- get(paste("covAll_chipAll_",cnvType,"_KLn.df",sep=""),res.ls)
		z_set2.df    <- resObjectKln[match (z_set.gsid, resObjectKln$GsID), ]
		height_U.mx <- matrix (data = c (z_set2.df$Pvalue_U_dev_s), nrow = 1, byrow = T)
		height_Uc.mx <- matrix (data = c (z_set2.df$Coeff_U), nrow = 1, byrow = T)
		z_set2.col[which (z_set2.df$FDR_BH_U <= FDRThreshold)] <- "brown" 
		border.vc <- c (z_set2.col)
		nc.vc <- z_set2.col
	}

	pdf(paste(cnvType,"_EffectSize_U_",timestamp,".pdf",sep=""))
	par (mar = c (plotHeight, 4, 4, 2), mgp=c(3,1,0), lwd = 2)
	barplot (
		main = paste(cnvType,": Effect Size: U",sep=""),
		height = height_Uc.mx,
		names.arg = z_set.labels, ylim = c (min (0,height_Uc.mx), ceiling(max (height_Uc.mx))), cex.names = labelSize,
		beside = T, las = 2, col = rep (c ("gray60", "gray90"), times = length (z_set.labels)), 
		border = as.character (matrix (data = border.vc, ncol = length (nc.vc), byrow = T)),
		ylab = "Coeff")
	dev.off()

	z_set.gsid   <- gsList
	z_set.labels <- paste (z_set.gsid, gs_len.nv[z_set.gsid], sep = ": ")
	z_col.names  <- c ("CASE_g1n", "CTRL_g1n", "CASE_g2n", "CTRL_g2n", "CASE_g3n", "CTRL_g3n")
	pdf(paste(cnvType,"_Support_",timestamp,".pdf",sep=""))
	par (mar = c (plotHeight, 4, 4, 2), mgp=c(3,1,0), lwd = 1)
	if (Kl == "ALL"){
		height_s.mx <- rbind (t (z_set1.df[, z_col.names]), t (z_set2.df[, z_col.names]))
	} else if (Kl == "YES"){
		height_s.mx <- t (z_set1.df[, z_col.names])
	} else if (Kl == "NO"){
		height_s.mx <- t (z_set2.df[, z_col.names])
	}
	barplot (
		main = paste(cnvType,": Support",sep=""),
		height = height_s.mx,
		names.arg = z_set.labels, ylim = c (min(0,height_s.mx), max(height_s.mx)), cex.names = labelSize, 
		beside = T, las = 2, col = rep (c ("salmon", "skyblue"), times = length (z_set.labels)), 
		border = "gray30",
		ylab = "CASE and CTRL subject %")
	dev.off()

	z_set.gsid     <- gsList
	z_set.labels   <- paste (z_set.gsid, gs_len.nv[z_set.gsid], sep = ": ")

	if (Kl == "ALL"){
		z_set1_U.col    <- rep ("gray30", length (z_set.gsid))
		z_set2_U.col    <- rep ("gray30", length (z_set.gsid))
		z_set1_TL.col    <- rep ("gray30", length (z_set.gsid))
		z_set2_TL.col    <- rep ("gray30", length (z_set.gsid))
		z_set1_CNML.col    <- rep ("gray30", length (z_set.gsid))
		z_set2_CNML.col    <- rep ("gray30", length (z_set.gsid))
		height_U.mx    <- matrix (data = c (z_set1.df$Pvalue_U_dev_s, z_set2.df$Pvalue_U_dev_s), nrow = 2, byrow = T)
		height_TL.mx   <- matrix (data = c (z_set1.df$Pvalue_TL_dev_s, z_set2.df$Pvalue_TL_dev_s), nrow = 2, byrow = T)
		height_CNML.mx <- matrix (data = c (z_set1.df$Pvalue_CNML_dev_s, z_set2.df$Pvalue_CNML_dev_s), nrow = 2, byrow = T)
		z_set1_U.col[which (z_set1.df$FDR_BH_U <= FDRThreshold)] <- "brown"
		z_set2_U.col[which (z_set2.df$FDR_BH_U <= FDRThreshold)] <- "brown" 
		z_set1_TL.col[which (z_set1.df$FDR_BH_TL <= FDRThreshold)] <- "brown"
		z_set2_TL.col[which (z_set2.df$FDR_BH_TL <= FDRThreshold)] <- "brown" 
		z_set1_CNML.col[which (z_set1.df$FDR_BH_CNML <= FDRThreshold)] <- "brown"
		z_set2_CNML.col[which (z_set2.df$FDR_BH_CNML <= FDRThreshold)] <- "brown" 
		border_U.vc <- c (z_set1_U.col, z_set2_U.col)
		border_TL.vc <- c (z_set1_TL.col, z_set2_TL.col)
		border_CNML.vc <- c (z_set1_CNML.col, z_set2_CNML.col)
		nc_U.vc <- z_set1_U.col
		nc_TL.vc <- z_set1_TL.col
		nc_CNML.vc <- z_set1_CNML.col
	} else if (Kl == "YES"){
		z_set1_U.col    <- rep ("gray30", length (z_set.gsid))
		z_set1_TL.col    <- rep ("gray30", length (z_set.gsid))
		z_set1_CNML.col    <- rep ("gray30", length (z_set.gsid))
		height_U.mx    <- matrix (data = c (z_set1.df$Pvalue_U_dev_s), nrow = 1, byrow = T)
		height_TL.mx   <- matrix (data = c (z_set1.df$Pvalue_TL_dev_s), nrow = 1, byrow = T)
		height_CNML.mx <- matrix (data = c (z_set1.df$Pvalue_CNML_dev_s), nrow = 1, byrow = T)
		z_set1_U.col[which (z_set1.df$FDR_BH_U <= FDRThreshold)] <- "brown"
		border_U.vc <- c (z_set1_U.col)
		z_set1_TL.col[which (z_set1.df$FDR_BH_TL <= FDRThreshold)] <- "brown"
		border_TL.vc <- c (z_set1_TL.col)
		z_set1_CNML.col[which (z_set1.df$FDR_BH_CNML <= FDRThreshold)] <- "brown"
		border_CNML.vc <- c (z_set1_CNML.col)
		nc_U.vc <- z_set1_U.col
		nc_TL.vc <- z_set1_TL.col
		nc_CNML.vc <- z_set1_CNML.col
	} else if (Kl == "NO"){
		z_set2_U.col    <- rep ("gray30", length (z_set.gsid))
		z_set2_TL.col    <- rep ("gray30", length (z_set.gsid))
		z_set2_CNML.col    <- rep ("gray30", length (z_set.gsid))
		height_U.mx    <- matrix (data = c (z_set2.df$Pvalue_U_dev_s), nrow = 1, byrow = T)
		height_TL.mx   <- matrix (data = c (z_set2.df$Pvalue_TL_dev_s), nrow = 1, byrow = T)
		height_CNML.mx <- matrix (data = c (z_set2.df$Pvalue_CNML_dev_s), nrow = 1, byrow = T)
		z_set2_U.col[which (z_set2.df$FDR_BH_U <= FDRThreshold)] <- "brown"
		border_U.vc <- c (z_set2_U.col)
		z_set2_TL.col[which (z_set2.df$FDR_BH_TL <= FDRThreshold)] <- "brown"
		border_TL.vc <- c (z_set2_TL.col)
		z_set2_CNML.col[which (z_set2.df$FDR_BH_CNML <= FDRThreshold)] <- "brown"
		border_CNML.vc <- c (z_set2_CNML.col)
		nc_U.vc <- z_set2_U.col
		nc_TL.vc <- z_set2_TL.col
		nc_CNML.vc <- z_set2_CNML.col
	}
	min.n          <- min (c (min (height_U.mx), min (height_TL.mx), min (height_CNML.mx)))
	max.n          <- max (c (max (height_U.mx), max (height_TL.mx), max (height_CNML.mx)))

	z_set1.col    <- rep ("gray30", length (z_set.gsid))
	z_set2.col    <- rep ("gray30", length (z_set.gsid))
	pdf(paste(cnvType,"_Significance_Compare_U_",timestamp,".pdf",sep=""))
	par (mar = c (plotHeight, 4, 4, 2), mgp=c(3,1,0), lwd = 2)
	barplot (
		main = paste(cnvType,": Compare Significance: U",sep=""),
		height = height_U.mx,
		names.arg = z_set.labels, ylim = c (min(0,min.n), ceiling(max.n)), cex.names = labelSize,
		beside = T, las = 2, col = rep (c ("gray60", "gray90"), times = length (z_set.labels)), 
		border = as.character (matrix (data = border_U.vc, ncol = length (nc_U.vc), byrow = T)),
		ylab = "-Log (Dev P-value) * sign (Coeff)")
	dev.off()

	pdf(paste(cnvType,"_Significance_Compare_TL_",timestamp,".pdf",sep=""))
	par (mar = c (plotHeight, 4, 4, 2), mgp=c(3,1,0), lwd = 2)
	barplot (
		main = paste(cnvType,": Compare Significance: TL",sep=""),
		height = height_TL.mx,
		names.arg = z_set.labels, ylim = c (min(0,min.n), ceiling(max.n)), cex.names = labelSize,
		beside = T, las = 2, col = rep (c ("gray60", "gray90"), times = length (z_set.labels)), 
		border = as.character (matrix (data = border_TL.vc, ncol = length (nc_TL.vc), byrow = T)),
		ylab = "-Log (Dev P-value) * sign (Coeff)")
	dev.off()

	pdf(paste(cnvType,"_Significance_Compare_CNML_",timestamp,".pdf",sep=""))
	par (mar = c (plotHeight, 4, 4, 2), mgp=c(3,1,0), lwd = 2)
	barplot (
		main = paste(cnvType,": Compare Significance: CNML",sep=""),
		height = height_CNML.mx,
		names.arg = z_set.labels, ylim = c (min(0,min.n), ceiling(max.n)), cex.names = labelSize,
		beside = T, las = 2, col = rep (c ("gray60", "gray90"), times = length (z_set.labels)), 
		border = as.character (matrix (data = border_CNML.vc, ncol = length (nc_CNML.vc), byrow = T)),
		ylab = "-Log (Dev P-value) * sign (Coeff)")
	dev.off()
}