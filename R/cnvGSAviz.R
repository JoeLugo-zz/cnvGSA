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
	vizData <- list(cnvGSA.in@gsData.ls$gs_all.ls,cnvGSA.out@res.ls,cnvGSA.in@cnvData.ls$cnv.df,cnvGSA.in@phData.ls$ph_TYPE.df)
	names(vizData) <- list("gs_all.ls","res.ls","cnv.df","ph_TYPE.df")

	setwd (cnvGSA.in@config.ls$outputPath)
	save(vizData,file=paste("vizInput.RData",sep=""))

	config.df <- cnvGSA.out@config.df

	FDRThreshold  <- as.numeric(config.df[config.df$param == "FDRThreshold","value"])
	gsList        <- config.df[config.df$param == "gsList","value"]
	cnvType       <- config.df[config.df$param == "cnvType","value"]
	plotHeight    <- as.numeric(config.df[config.df$param == "plotHeight","value"])
	outputPathViz <- config.df[config.df$param == "outputPathViz","value"]
	labelSize     <- as.numeric(config.df[config.df$param == "labelSize","value"])
	gsList        <- config.df[config.df$param == "gsList","value"]

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
		stop("No gsList in the gene list.")
	}

	resObjectKly <- get(paste("covAll_chipAll_",cnvType,"_KLy.df",sep=""),res.ls)
	resObjectKln <- get(paste("covAll_chipAll_",cnvType,"_KLn.df",sep=""),res.ls)

	gs_len.nv <- sapply (gs_all.ls, length)

	setwd (outputPathViz)

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
	write.table (z_olp.df, col.names = T, row.names = F, sep = "\t", quote = F, file = paste("GsOverlap_NeurofSynaptic.txt",sep=""))

	# cnvType: NEUROFUNCTION + SYNAPTIC

	z_set.gsid   <- gsList
	z_set.labels <- paste (z_set.gsid, gs_len.nv[z_set.gsid], sep = ": ")
	z_set1.df    <- resObjectKly[match (z_set.gsid, resObjectKly$GsID), ]
	z_set2.df    <- resObjectKln[match (z_set.gsid, resObjectKln$GsID), ]
	z_set1.col   <- rep ("gray30", length (z_set.gsid))
	z_set2.col   <- rep ("gray30", length (z_set.gsid))
	z_set1.col[which (z_set1.df$FDR_BH_U <= FDRThreshold)] <- "brown" # let user pick this threshold
	z_set2.col[which (z_set2.df$FDR_BH_U <= FDRThreshold)] <- "brown" 
	# leave this as default and let the user change if they want to 
	pdf(paste(cnvType,"_NeurofSyn_Significance_U.pdf",sep=""))
	par (mar = c (plotHeight, 4, 4, 2), mgp=c(3,3,0), lwd = 2) # controls the margins of the graph
	height_U.mx <- matrix (data = c (z_set1.df$Pvalue_U_dev_s, z_set2.df$Pvalue_U_dev_s), nrow = 2, byrow = T)
	barplot (
		main = paste(cnvType,": Neurof+Synaptic: Significance: U - ",sep=""),
		height = height_U.mx,
		names.arg = z_set.labels, ylim = c (min (height_U.mx), max (height_U.mx)), cex.names = labelSize,
		beside = T, las = 2, col = rep (c ("gray60", "gray90"), times = length (z_set.labels)), 
		border = as.character (matrix (data = c (z_set1.col, z_set2.col), ncol = length (z_set1.col), byrow = T)),
		ylab = "-Log (Dev P-value) * sign (Coeff)")
	dev.off()

	z_set1.col    <- rep ("gray30", length (z_set.gsid))
	z_set2.col    <- rep ("gray30", length (z_set.gsid))
	z_set1.col[which (z_set1.df$FDR_BH_U <= FDRThreshold)] <- "brown"
	z_set2.col[which (z_set2.df$FDR_BH_U <= FDRThreshold)] <- "brown" 
	pdf(paste(cnvType,"_NeurofSyn_EffectSize_U.pdf",sep=""))
	par (mar = c (plotHeight, 4, 4, 2), mgp=c(3,3,0), lwd = 2)
	height_Uc.mx <- matrix (data = c (z_set1.df$Coeff_U, z_set2.df$Coeff_U), nrow = 2, byrow = T)
	barplot (
		main = paste(cnvType,": Neurof+Synaptic: Effect Size: U",sep=""),
		height = height_Uc.mx,
		names.arg = z_set.labels, ylim = c (min (height_Uc.mx), max (height_Uc.mx)), cex.names = labelSize,
		beside = T, las = 2, col = rep (c ("gray60", "gray90"), times = length (z_set.labels)), 
		border = as.character (matrix (data = c (z_set1.col, z_set2.col), ncol = length (z_set1.col), byrow = T)),
		ylab = "Coeff")
	dev.off()

	z_set.gsid   <- gsList
	z_set.labels <- paste (z_set.gsid, gs_len.nv[z_set.gsid], sep = ": ")
	z_set1.df    <- resObjectKly[match (z_set.gsid, resObjectKly$GsID), ]
	z_set2.df    <- resObjectKln[match (z_set.gsid, resObjectKln$GsID), ]
	z_col.names  <- c ("SZ_g1n", "CT_g1n", "SZ_g2n", "CT_g2n", "SZ_g3n", "CT_g3n")
	pdf(paste(cnvType,"_NeurofSyn_Support.pdf",sep=""))
	par (mar = c (plotHeight, 4, 4, 2), mgp=c(3,3,0), lwd = 1)
	height_s.mx <- rbind (t (z_set1.df[, z_col.names]), t (z_set2.df[, z_col.names]))
	barplot (
		main = paste(cnvType,": Neurof+Synaptic: Support",sep=""),
		height = height_s.mx,
		names.arg = z_set.labels, ylim = c (0, 8), cex.names = labelSize, 
		beside = T, las = 2, col = rep (c ("salmon", "skyblue"), times = length (z_set.labels)), 
		border = "gray30",
		ylab = "SZ and CT subject %")
	dev.off()

	z_set.gsid     <- gsList
	z_set.labels   <- paste (z_set.gsid, gs_len.nv[z_set.gsid], sep = ": ")
	z_set1.df      <- resObjectKly[match (z_set.gsid, resObjectKly$GsID), ]
	z_set2.df      <- resObjectKln[match (z_set.gsid, resObjectKln$GsID), ]
	height_U.mx    <- matrix (data = c (z_set1.df$Pvalue_U_dev_s, z_set2.df$Pvalue_U_dev_s), nrow = 2, byrow = T)
	height_TL.mx   <- matrix (data = c (z_set1.df$Pvalue_TL_dev_s, z_set2.df$Pvalue_TL_dev_s), nrow = 2, byrow = T)
	height_CNML.mx <- matrix (data = c (z_set1.df$Pvalue_CNML_dev_s, z_set2.df$Pvalue_CNML_dev_s), nrow = 2, byrow = T)
	min.n          <- min (c (min (height_U.mx), min (height_TL.mx), min (height_CNML.mx)))
	max.n          <- max (c (max (height_U.mx), max (height_TL.mx), max (height_CNML.mx)))

	z_set1.col    <- rep ("gray30", length (z_set.gsid))
	z_set2.col    <- rep ("gray30", length (z_set.gsid))
	z_set1.col[which (z_set1.df$FDR_BH_U <= FDRThreshold)] <- "brown"
	z_set2.col[which (z_set2.df$FDR_BH_U <= FDRThreshold)] <- "brown" 
	pdf(paste(cnvType,"_NeurofSyn_Significance_Compare_U.pdf",sep=""))
	par (mar = c (plotHeight, 4, 4, 2), mgp=c(3,3,0), lwd = 2)
	barplot (
		main = paste(cnvType,": Neurof+Synaptic: Compare Significance: U",sep=""),
		height = height_U.mx,
		names.arg = z_set.labels, ylim = c (min.n, max.n), cex.names = labelSize,
		beside = T, las = 2, col = rep (c ("gray60", "gray90"), times = length (z_set.labels)), 
		border = as.character (matrix (data = c (z_set1.col, z_set2.col), ncol = length (z_set1.col), byrow = T)),
		ylab = "-Log (Dev P-value) * sign (Coeff)")
	dev.off()

	z_set1.col    <- rep ("gray30", length (z_set.gsid))
	z_set2.col    <- rep ("gray30", length (z_set.gsid))
	z_set1.col[which (z_set1.df$FDR_BH_TL <= FDRThreshold)] <- "brown"
	z_set2.col[which (z_set2.df$FDR_BH_TL <= FDRThreshold)] <- "brown" 
	pdf(paste(cnvType,"_NeurofSyn_Significance_Compare_TL.pdf",sep=""))
	par (mar = c (plotHeight, 4, 4, 2), mgp=c(3,3,0), lwd = 2)
	barplot (
		main = paste(cnvType,": Neurof+Synaptic: Compare Significance: TL",sep=""),
		height = height_TL.mx,
		names.arg = z_set.labels, ylim = c (min.n, max.n), cex.names = labelSize,
		beside = T, las = 2, col = rep (c ("gray60", "gray90"), times = length (z_set.labels)), 
		border = as.character (matrix (data = c (z_set1.col, z_set2.col), ncol = length (z_set1.col), byrow = T)),
		ylab = "-Log (Dev P-value) * sign (Coeff)")
	dev.off()

	z_set1.col    <- rep ("gray30", length (z_set.gsid))
	z_set2.col    <- rep ("gray30", length (z_set.gsid))
	z_set1.col[which (z_set1.df$FDR_BH_CNML <= FDRThreshold)] <- "brown"
	z_set2.col[which (z_set2.df$FDR_BH_CNML <= FDRThreshold)] <- "brown" 
	pdf(paste(cnvType,"_NeurofSyn_Significance_Compare_CNML.pdf",sep=""))
	par (mar = c (plotHeight, 4, 4, 2), mgp=c(3,3,0), lwd = 2)
	barplot (
		main = paste(cnvType,": Neurof+Synaptic: Compare Significance: CNML",sep=""),
		height = height_CNML.mx,
		names.arg = z_set.labels, ylim = c (min.n, max.n), cex.names = labelSize,
		beside = T, las = 2, col = rep (c ("gray60", "gray90"), times = length (z_set.labels)), 
		border = as.character (matrix (data = c (z_set1.col, z_set2.col), ncol = length (z_set1.col), byrow = T)),
		ylab = "-Log (Dev P-value) * sign (Coeff)")
	dev.off()

	# 2.
	# NEUROPHENOTYPE + PHENOTYPE PANEL
	z_set.gsid  <- gsList
	z_olp.n     <- length (z_set.gsid)
	z_olp.mx    <- matrix (data = NA, ncol = z_olp.n, nrow = z_olp.n, dimnames = list (z_set.gsid, z_set.gsid))
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
	write.table (z_olp.df, col.names = T, row.names = F, sep = "\t", quote = F, file = paste("GsOverlap_Phenotype.txt",sep=""))

	z_set.labels <- paste (z_set.gsid, gs_len.nv[z_set.gsid], sep = ": ")

	z_set1.df     <- resObjectKly[match (z_set.gsid, resObjectKly$GsID), ]
	z_set2.df     <- resObjectKln[match (z_set.gsid, resObjectKln$GsID), ]
	z_set1.col    <- rep ("gray30", length (z_set.gsid))
	z_set2.col    <- rep ("gray30", length (z_set.gsid))
	z_set1.col[which (z_set1.df$FDR_BH_U <= 0.1)] <- "brown"
	z_set2.col[which (z_set2.df$FDR_BH_U <= 0.1)] <- "brown" 
	pdf(paste(cnvType,"_NeurofPhen_Significance_U.pdf",sep=""))
	par (mar = c (plotHeight, 4, 4, 2), mgp=c(3,3,0), lwd = 2)
	height_U.mx <- matrix (data = c (z_set1.df$Pvalue_U_dev_s, z_set2.df$Pvalue_U_dev_s), nrow = 2, byrow = T)
	barplot (
		main = paste(cnvType,": Phenotype: Significance: U",sep=""),
		height = height_U.mx,
		names.arg = z_set.labels, ylim = c (min (height_U.mx), max (height_U.mx)), cex.names = labelSize,
		beside = T, las = 2, col = rep (c ("gray60", "gray90"), times = length (z_set.labels)), 
		border = as.character (matrix (data = c (z_set1.col, z_set2.col), ncol = length (z_set1.col), byrow = T)),
		ylab = "-Log (Dev P-value) * sign (Coeff)")
	dev.off()

	z_set1.col    <- rep ("gray30", length (z_set.gsid))
	z_set2.col    <- rep ("gray30", length (z_set.gsid))
	z_set1.col[which (z_set1.df$FDR_BH_U <= 0.1)] <- "brown"
	z_set2.col[which (z_set2.df$FDR_BH_U <= 0.1)] <- "brown" 
	pdf(paste(cnvType,"_NeurofPhen_EffectSize_U.pdf",sep=""))
	par (mar = c (plotHeight, 4, 4, 2), mgp=c(3,3,0), lwd = 2)
	height_Uc.mx <- matrix (data = c (z_set1.df$Coeff_U, z_set2.df$Coeff_U), nrow = 2, byrow = T)
	barplot (
		main = paste(cnvType,": Phenotype: Effect Size: U",sep=""),
		height = height_Uc.mx,
		names.arg = z_set.labels, ylim = c (min (height_Uc.mx), max (height_Uc.mx)), cex.names = labelSize,
		beside = T, las = 2, col = rep (c ("gray60", "gray90"), times = length (z_set.labels)), 
		border = as.character (matrix (data = c (z_set1.col, z_set2.col), ncol = length (z_set1.col), byrow = T)),
		ylab = "Coeff")
	dev.off()

	z_col.names <- c ("SZ_g1n", "CT_g1n", "SZ_g2n", "CT_g2n", "SZ_g3n", "CT_g3n")
	pdf(paste(cnvType,"_NeurofPhen_Support.pdf",sep=""))
	par (mar = c (plotHeight, 4, 4, 2), mgp=c(3,3,0), lwd = 1)
	height_s.mx <- rbind (t (z_set1.df[, z_col.names]), t (z_set2.df[, z_col.names]))
	barplot (
		main = paste(cnvType,": Phenotype: Support",sep=""),
		height = height_s.mx,
		names.arg = z_set.labels, ylim = c (min (height_s.mx), max (height_s.mx)), cex.names = labelSize, 
		beside = T, las = 2, col = rep (c ("salmon", "skyblue"), times = length (z_set.labels)), 
		border = "gray30",
		ylab = "SZ and CT subject %")
	dev.off()

	z_set.gsid     <- gsList
	z_set.labels   <- paste (z_set.gsid, gs_len.nv[z_set.gsid], sep = ": ")
	z_set1.df      <- resObjectKly[match (z_set.gsid, resObjectKly$GsID), ]
	z_set2.df      <- resObjectKln[match (z_set.gsid, resObjectKln$GsID), ]
	height_U.mx    <- matrix (data = c (z_set1.df$Pvalue_U_dev_s, z_set2.df$Pvalue_U_dev_s), nrow = 2, byrow = T)
	height_TL.mx   <- matrix (data = c (z_set1.df$Pvalue_TL_dev_s, z_set2.df$Pvalue_TL_dev_s), nrow = 2, byrow = T)
	height_CNML.mx <- matrix (data = c (z_set1.df$Pvalue_CNML_dev_s, z_set2.df$Pvalue_CNML_dev_s), nrow = 2, byrow = T)
	min.n          <- min (c (min (height_U.mx), min (height_TL.mx), min (height_CNML.mx)))
	max.n          <- max (c (max (height_U.mx), max (height_TL.mx), max (height_CNML.mx)))

	z_set1.col    <- rep ("gray30", length (z_set.gsid))
	z_set2.col    <- rep ("gray30", length (z_set.gsid))
	z_set1.col[which (z_set1.df$FDR_BH_U <= 0.1)] <- "brown"
	z_set2.col[which (z_set2.df$FDR_BH_U <= 0.1)] <- "brown" 
	pdf(paste(cnvType,"_NeurofPhen_Significance_Compare_U.pdf",sep=""))
	par (mar = c (plotHeight, 4, 4, 2), mgp=c(3,3,0), lwd = 2)
	barplot (
		main = paste(cnvType,": Phenotype: Compare Significance: U",sep=""),
		height = height_U.mx,
		names.arg = z_set.labels, ylim = c (min.n, max.n), cex.names = labelSize,
		beside = T, las = 2, col = rep (c ("gray60", "gray90"), times = length (z_set.labels)), 
		border = as.character (matrix (data = c (z_set1.col, z_set2.col), ncol = length (z_set1.col), byrow = T)),
		ylab = "-Log (Dev P-value) * sign (Coeff)")
	dev.off()

	z_set1.col    <- rep ("gray30", length (z_set.gsid))
	z_set2.col    <- rep ("gray30", length (z_set.gsid))
	z_set1.col[which (z_set1.df$FDR_BH_TL <= 0.1)] <- "brown"
	z_set2.col[which (z_set2.df$FDR_BH_TL <= 0.1)] <- "brown" 
	pdf(paste(cnvType,"_NeurofPhen_Significance_Compare_TL.pdf",sep=""))
	par (mar = c (plotHeight, 4, 4, 2), mgp=c(3,3,0), lwd = 2)
	barplot (
		main = paste(cnvType,": Phenotype: Compare Significance: TL",sep=""),
		height = height_TL.mx,
		names.arg = z_set.labels, ylim = c (min.n, max.n), cex.names = labelSize,
		beside = T, las = 2, col = rep (c ("gray60", "gray90"), times = length (z_set.labels)), 
		border = as.character (matrix (data = c (z_set1.col, z_set2.col), ncol = length (z_set1.col), byrow = T)),
		ylab = "-Log (Dev P-value) * sign (Coeff)")
	dev.off()

	z_set1.col    <- rep ("gray30", length (z_set.gsid))
	z_set2.col    <- rep ("gray30", length (z_set.gsid))
	z_set1.col[which (z_set1.df$FDR_BH_CNML <= 0.1)] <- "brown"
	z_set2.col[which (z_set2.df$FDR_BH_CNML <= 0.1)] <- "brown" 
	pdf(paste(cnvType,"_NeurofPhen_Significance_Compare_CNML.pdf",sep=""))
	par (mar = c (plotHeight, 4, 4, 2), mgp=c(3,3,0), lwd = 2)
	barplot (
		main = paste(cnvType,": Phenotype: Compare Significance: CNML",sep=""),
		height = height_CNML.mx,
		names.arg = z_set.labels, ylim = c (min.n, max.n), cex.names = labelSize,
		beside = T, las = 2, col = rep (c ("gray60", "gray90"), times = length (z_set.labels)), 
		border = as.character (matrix (data = c (z_set1.col, z_set2.col), ncol = length (z_set1.col), byrow = T)),
		ylab = "-Log (Dev P-value) * sign (Coeff)")
	dev.off()
}

# f.makeViz(configPath = "/Users/josephlugo/Documents/R/PGC2_test/R_Works/",configFile = "PGC2_config.txt")

# # ======

# # STEP-DOWN ANALYSIS:paste(cnvType,"

# covariates.chv <- c ("SEX", "CNV_metric", "CNV_platform", 
# 				     "C1", "C2", "C3", "C4", "C8") 

# f.getCoeff_sm <- function (glm.sm, var.ch) {
# 	     if (var.ch %in% rownames (glm.sm$coefficients)) {return (glm.sm$coefficients[var.ch, "Estimate"])}
# 		 else {return (NA)} }
# f.getPval_sm  <- function (glm.sm, var.ch) {
# 		if (var.ch %in% rownames (glm.sm$coefficients)) {return (glm.sm$coefficients[var.ch, "Pr(>|z|)"])}
# 		else {return (NA)} }
# f.getPval_anova  <- function (glm.anova, var.ch) {
# 		if (var.ch %in% rownames (glm.anova["Pr(>Chi)"])) {return (glm.anova["Pr(>Chi)"][var.ch, "Pr(>Chi)"])}
# 		else {return (NA)} }

# # 1. SORT BY SIGNIFICANCE

# var.chv <- NULL

# var.ch <- "GS__Neurof_GoSynaptic_paste(cnvType,""
# var.chv  <- c (var.chv, var.ch)
# glm_U_form.ch <- paste ("Condition", "~", paste (covariates.chv, collapse = " + "), "+", paste ("GS__U", "LOSS", sep = "__") , "+", paste (var.chv, collapse = " + "), sep = " ")
# x_U.glm    <- glm (as.formula (glm_U_form.ch), ph.df, family = binomial (logit))
# x_U.glm_sm <- summary (x_U.glm)
# x_U.anova  <- anova   (x_U.glm, test = "Chisq")
# f.getCoeff_sm   (x_U.glm_sm, var.ch) # 0.5208341
# f.getPval_sm    (x_U.glm_sm, var.ch) # 8.604705e-11
# f.getPval_anova (x_U.anova,  var.ch) # 2.81811e-11

# var.ch <- "GS__Neurof_UnionStringent__LOSS"
# var.chv  <- c (var.chv, var.ch)
# glm_U_form.ch <- paste ("Condition", "~", paste (covariates.chv, collapse = " + "), "+", paste ("GS__U", "LOSS", sep = "__") , "+", paste (var.chv, collapse = " + "), sep = " ")
# x_U.glm    <- glm (as.formula (glm_U_form.ch), ph.df, family = binomial (logit))
# x_U.glm_sm <- summary (x_U.glm)
# x_U.anova  <- anova   (x_U.glm, test = "Chisq")
# f.getCoeff_sm   (x_U.glm_sm, var.ch) # 0.1266792
# f.getPval_sm    (x_U.glm_sm, var.ch) # 0.09955865
# f.getPval_anova (x_U.anova,  var.ch) # 0.09764339

# var.ch <- "GS__Neurof_GoNeuronProj__LOSS"
# var.chv  <- c (var.chv, var.ch)
# glm_U_form.ch <- paste ("Condition", "~", paste (covariates.chv, collapse = " + "), "+", paste ("GS__U", "LOSS", sep = "__") , "+", paste (var.chv, collapse = " + "), sep = " ")
# x_U.glm    <- glm (as.formula (glm_U_form.ch), ph.df, family = binomial (logit))
# x_U.glm_sm <- summary (x_U.glm)
# x_U.anova  <- anova   (x_U.glm, test = "Chisq")
# f.getCoeff_sm   (x_U.glm_sm, var.ch) # -0.1416699
# f.getPval_sm    (x_U.glm_sm, var.ch) # 0.2594801
# f.getPval_anova (x_U.anova,  var.ch) # 0.2582

# var.ch <- "GS__Kirov_ARC__LOSS"
# var.chv  <- c (var.chv, var.ch)
# glm_U_form.ch <- paste ("Condition", "~", paste (covariates.chv, collapse = " + "), "+", paste ("GS__U", "LOSS", sep = "__") , "+", paste (var.chv, collapse = " + "), sep = " ")
# x_U.glm    <- glm (as.formula (glm_U_form.ch), ph.df, family = binomial (logit))
# x_U.glm_sm <- summary (x_U.glm)
# x_U.anova  <- anova   (x_U.glm, test = "Chisq")
# f.getCoeff_sm   (x_U.glm_sm, var.ch) # 0.2689707
# f.getPval_sm    (x_U.glm_sm, var.ch) # 0.1498327
# f.getPval_anova (x_U.anova,  var.ch) # 0.14704

# var.ch <- "GS__FMR1_Targets_Darnell__LOSS"
# var.chv  <- c (var.chv, var.ch)
# glm_U_form.ch <- paste ("Condition", "~", paste (covariates.chv, collapse = " + "), "+", paste ("GS__U", "LOSS", sep = "__") , "+", paste (var.chv, collapse = " + "), sep = " ")
# x_U.glm    <- glm (as.formula (glm_U_form.ch), ph.df, family = binomial (logit))
# x_U.glm_sm <- summary (x_U.glm)
# x_U.anova  <- anova   (x_U.glm, test = "Chisq")
# f.getCoeff_sm   (x_U.glm_sm, var.ch) # 0.1211829
# f.getPval_sm    (x_U.glm_sm, var.ch) # 0.1334014
# f.getPval_anova (x_U.anova,  var.ch) # 0.1310887

# # 2. SORT BY EFFECT SIZE

# var.chv <- NULL

# var.ch <- "GS__Kirov_ARC__LOSS"
# var.chv  <- c (var.chv, var.ch)
# glm_U_form.ch <- paste ("Condition", "~", paste (covariates.chv, collapse = " + "), "+", paste ("GS__U", "LOSS", sep = "__") , "+", paste (var.chv, collapse = " + "), sep = " ")
# x_U.glm    <- glm (as.formula (glm_U_form.ch), ph.df, family = binomial (logit))
# x_U.glm_sm <- summary (x_U.glm)
# x_U.anova  <- anova   (x_U.glm, test = "Chisq")
# f.getCoeff_sm   (x_U.glm_sm, var.ch) # 0.608871
# f.getPval_sm    (x_U.glm_sm, var.ch) # 0.0002645162
# f.getPval_anova (x_U.anova,  var.ch) # 0.0001808685

# var.ch <- "GS__Neurof_GoSynaptic__LOSS"
# var.chv  <- c (var.chv, var.ch)
# glm_U_form.ch <- paste ("Condition", "~", paste (covariates.chv, collapse = " + "), "+", paste ("GS__U", "LOSS", sep = "__") , "+", paste (var.chv, collapse = " + "), sep = " ")
# x_U.glm    <- glm (as.formula (glm_U_form.ch), ph.df, family = binomial (logit))
# x_U.glm_sm <- summary (x_U.glm)
# x_U.anova  <- anova   (x_U.glm, test = "Chisq")
# f.getCoeff_sm   (x_U.glm_sm, var.ch) # 0.4787264
# f.getPval_sm    (x_U.glm_sm, var.ch) # 4.372772e-08
# f.getPval_anova (x_U.anova,  var.ch) # 1.864603e-08

# var.ch <- "GS__Neurof_UnionStringent__LOSS"
# var.chv  <- c (var.chv, var.ch)
# glm_U_form.ch <- paste ("Condition", "~", paste (covariates.chv, collapse = " + "), "+", paste ("GS__U", "LOSS", sep = "__") , "+", paste (var.chv, collapse = " + "), sep = " ")
# x_U.glm    <- glm (as.formula (glm_U_form.ch), ph.df, family = binomial (logit))
# x_U.glm_sm <- summary (x_U.glm)
# x_U.anova  <- anova   (x_U.glm, test = "Chisq")
# f.getCoeff_sm   (x_U.glm_sm, var.ch) # 0.1266792
# f.getPval_sm    (x_U.glm_sm, var.ch) # 0.0958683
# f.getPval_anova (x_U.anova,  var.ch) # 0.09401131

# var.ch <- "GS__FMR1_Targets_Darnell__LOSS"
# var.chv  <- c (var.chv, var.ch)
# glm_U_form.ch <- paste ("Condition", "~", paste (covariates.chv, collapse = " + "), "+", paste ("GS__U", "LOSS", sep = "__") , "+", paste (var.chv, collapse = " + "), sep = " ")
# x_U.glm    <- glm (as.formula (glm_U_form.ch), ph.df, family = binomial (logit))
# x_U.glm_sm <- summary (x_U.glm)
# x_U.anova  <- anova   (x_U.glm, test = "Chisq")
# f.getCoeff_sm   (x_U.glm_sm, var.ch) # 0.1158742
# f.getPval_sm    (x_U.glm_sm, var.ch) # 0.1491939
# f.getPval_anova (x_U.anova,  var.ch) # 0.1469888

# var.ch <- "GS__Neurof_GoNeuronProj__LOSS"
# var.chv  <- c (var.chv, var.ch)
# glm_U_form.ch <- paste ("Condition", "~", paste (covariates.chv, collapse = " + "), "+", paste ("GS__U", "LOSS", sep = "__") , "+", paste (var.chv, collapse = " + "), sep = " ")
# x_U.glm    <- glm (as.formula (glm_U_form.ch), ph.df, family = binomial (logit))
# x_U.glm_sm <- summary (x_U.glm)
# x_U.anova  <- anova   (x_U.glm, test = "Chisq")
# f.getCoeff_sm   (x_U.glm_sm, var.ch) # -0.186064
# f.getPval_sm    (x_U.glm_sm, var.ch) # 0.1439376
# f.getPval_anova (x_U.anova,  var.ch) # 0.1430069

# # 3. FINAL AGGREGATE

# var.chv <- NULL

# var.ch <- "GS__Neurof_GoSynaptic__LOSS"
# var.chv  <- c (var.chv, var.ch)
# glm_U_form.ch <- paste ("Condition", "~", paste (covariates.chv, collapse = " + "), "+", paste ("GS__U", "LOSS", sep = "__") , "+", paste (var.chv, collapse = " + "), sep = " ")
# x_U.glm    <- glm (as.formula (glm_U_form.ch), ph.df, family = binomial (logit))
# x_U.glm_sm <- summary (x_U.glm)
# x_U.anova  <- anova   (x_U.glm, test = "Chisq")
# f.getCoeff_sm   (x_U.glm_sm, var.ch) # 0.5208341
# f.getPval_sm    (x_U.glm_sm, var.ch) # 8.604705e-11
# f.getPval_anova (x_U.anova,  var.ch) # 2.81811e-11

# var.ch <- "GS__PhMm_NeuroBehav_all__LOSS"
# var.chv  <- c (var.chv, var.ch)
# glm_U_form.ch <- paste ("Condition", "~", paste (covariates.chv, collapse = " + "), "+", paste ("GS__U", "LOSS", sep = "__") , "+", paste (var.chv, collapse = " + "), sep = " ")
# x_U.glm    <- glm (as.formula (glm_U_form.ch), ph.df, family = binomial (logit))
# x_U.glm_sm <- summary (x_U.glm)
# x_U.anova  <- anova   (x_U.glm, test = "Chisq")
# f.getCoeff_sm   (x_U.glm_sm, var.ch) # 0.1091142
# f.getPval_sm    (x_U.glm_sm, var.ch) # 0.01041729
# f.getPval_anova (x_U.anova,  var.ch) # 0.009960622

# var.ch <- "GS__BspanVHM_EqlNat__LOSS"
# var.chv  <- c (var.chv, var.ch)
# glm_U_form.ch <- paste ("Condition", "~", paste (covariates.chv, collapse = " + "), "+", paste ("GS__U", "LOSS", sep = "__") , "+", paste (var.chv, collapse = " + "), sep = " ")
# x_U.glm    <- glm (as.formula (glm_U_form.ch), ph.df, family = binomial (logit))
# x_U.glm_sm <- summary (x_U.glm)
# x_U.anova  <- anova   (x_U.glm, test = "Chisq")
# f.getCoeff_sm   (x_U.glm_sm, var.ch) # 0.06639415
# f.getPval_sm    (x_U.glm_sm, var.ch) # 0.1027003
# f.getPval_anova (x_U.anova,  var.ch) # 0.1021357

# # ======

# # WRITE GENE BURDEN WITH CORRECTION

# setwd ("/Users/daniele/Documents/Works_Research/2014/PGC2/Results/ScriptP03_Viz_20150108h09m10")

# write.table (res.ls$covAll_chipAll_KLy.df[1: 3, ], sep = "\t", quote = F, col.names = T, row.names = F, file = "UGsetBurden_KLociY__MakeLogregTable_v06_20141208h18m25__Viz_20150108h09m10.txt")

# write.table (res.ls$covAll_chipAll_KLn.df[1: 3, ], sep = "\t", quote = F, col.names = T, row.names = F, file = "UGsetBurden_KLociN__MakeLogregTable_v06_20141208h18m25__Viz_20150108h09m10.txt")
