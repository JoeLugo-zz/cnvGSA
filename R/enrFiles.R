#' Prepare the files for the enrichment maps.
#'
#' @param cnvGSA.out A CnvGSAOutput S4 object
#' @return Returns a list with the data frames of the GMT file and the generic file. 

f.enrFiles <- function(cnvGSA.out)
{
	config.df <- cnvGSA.out@config.df

	Kl              <- config.df[config.df$param == "Kl","value"]
	cnvType         <- config.df[config.df$param == "cnvType","value"]
	pVal            <- config.df[config.df$param == "pVal","value"]
	FDR             <- config.df[config.df$param == "FDR","value"]
	coeff           <- config.df[config.df$param == "coeff","value"]
	keepCoeff       <- config.df[config.df$param == "keepCoeff","value"]
	outputPathEnr   <- config.df[config.df$param == "outputPathEnr","value"]
	filtGsEnr       <- config.df[config.df$param == "filtGsEnr","value"]
	minCaseCount    <- as.numeric(config.df[config.df$param == "minCaseCount","value"])
	minControlCount <- as.numeric(config.df[config.df$param == "minControlCount","value"])
	minRatio        <- as.numeric(config.df[config.df$param == "minRatio","value"])

	geneCount.tab <- cnvGSA.out@gsData.ls$geneCount.tab
	res.ls <- cnvGSA.out@res.ls

	# MAKING GENERIC FILE
	resKL <- get(paste("covAll_chipAll_",cnvType,"_KLy.df",sep=""),res.ls)
	if (Kl == "NO"){resKL <- get(paste("covAll_chipAll_",cnvType,"_KLn.df",sep=""),res.ls)}
	resKL$Phenotype <- -1
	resKL[which(resKL[,coeff] > 0),]$Phenotype <- 1

	if (keepCoeff == "NO")
	{
		resKL <- resKL[which(resKL$Phenotype == 1),]
	}

	enrGeneric      <- subset(resKL,select = c("GsID","GsName",pVal,FDR,"Phenotype"))
	enrGeneric$GsID <- as.factor(enrGeneric$GsID); enrGeneric$GsName <- as.factor(enrGeneric$GsName);

	# MAKING GMT FILE
	gs_sel_U.df <- cnvGSA.out@gsData.ls$gs_sel_U.df
	gs.ls       <- cnvGSA.out@gsData.ls$gs.ls
	gs_info.df  <- cnvGSA.out@gsData.ls$gs_info.df

	gs_info.df <- gs_info.df[order(gs_info.df$GsKey),]
	gs.ls      <- gs.ls[order(names(gs.ls))]

	lis        <- lapply(1:length(gs.ls),function(x) gs.ls[[x]])
	gsGenes.ls <- lis
	if (filtGsEnr == "YES"){
		usrFilt.df    <- geneCount.tab[which(geneCount.tab[,"1"] >= minCaseCount),]
		usrFilt.df    <- usrFilt.df[which(usrFilt.df[,"0"] >= minControlCount),]
		usrFilt.df    <- usrFilt.df[which(usrFilt.df[,"1"]/usrFilt.df[,"0"] >= minRatio),]
		filt.vc       <- rownames(usrFilt.df) # blacklist for genes to filter out
		gsGenes.ls    <- lapply(gs.ls,setdiff,filt.vc)
	}

	names(gsGenes.ls) <- gs_info.df$GsID
	diff              <- setdiff(names(gsGenes.ls),as.character(enrGeneric$GsID))
	gsGenes.ls        <- gsGenes.ls[which(!(names(gsGenes.ls) %in% diff))]

	f.collapse <- function (input.genes){
        paste (input.genes, collapse = "\t")
    }

	f.PackGMT <- function (id2eg.ls, id2des.chv, file.name){
        genes.chv  <- unlist(lapply (id2eg.ls, f.collapse))
        output.chv <- paste(names (id2eg.ls), id2des.chv, genes.chv, sep = "\t")
         
        cat (output.chv, sep = "\n", file = file.name)
    }

    setwd(outputPathEnr)
	write.table(enrGeneric,file="enrGeneric.txt",row.names=FALSE,sep="\t",quote=FALSE)
	f.PackGMT(id2eg.ls = gsGenes.ls, id2des.chv = gs_info.df$GsName, file.name = "enr.gmt")

	# lis 	   <- lapply(gs.ls,setdiff,filt.ls)
	# gsGenes.ls <- lis
	# for (i in 1:length(gs.ls))
	# {
	# 	# filters the genes in the gene sets
	# 	currentGs <- lis[[i]][1:length(lis[[i]])]
	# 	if(filtGs == "YES"){
	# 		gsGenes.ls[[i]] <- currentGs[which(currentGs %in% filt.ls)]
	# 	}
	# 	else {
	# 		gsGenes.ls[[i]] <- currentGs
	# 	}
	# }

	# lenGS <- unlist(lapply(1:length(gsGenes.ls),function(x) length(gsGenes.ls[[x]])))
	# maxGS <- max(lenGS)
	# extra_cols.df <- as.data.frame(matrix(NA,nrow=length(gsGenes.ls),ncol=maxGS))
	# gst_info.df   <- cbind(gs_info.df,extra_cols.df)

	# for (i in 1:length(gs.ls))
	# {
	# 	# populates gst_info.df with the genes from the list
	# 	currentGs <- gsGenes.ls[[i]][1:length(gsGenes.ls[[i]])]
	# 	gsLength  <- length(currentGs)
	# 	gst_info.df[i,5:gsLength] <- currentGs
	# }

	# diff <- setdiff(as.character(gst_info.df$GsID),as.character(enrGeneric$GsID))
	# gst_info.df <- gst_info.df[which(!(gst_info.df$GsID %in% diff)),]
	# gst_info.df <- subset(gst_info.df,select= -c(GsKey,GsSize))
	
	# write.table(gst_info.df,file=paste("gst_info.gmt",sep=""),sep="\t",row.names=FALSE,quote=FALSE,col.names=FALSE)
	# enr.ls <- list(enrGeneric,gst_info.df)
	# names(enr.ls) <- list("enrGeneric.df","enrGMT.df")
	# return(enr.ls)
}

# f.enrFiles(cnvGSA.out)
# f.enrFiles(configPath = "/Users/josephlugo/Documents/R/PGC2_test/R_Works/",configFile = "PGC2_config.txt",cnvGSA.out)

# gst_info.mx <- as.matrix(gst_info.df)
# gst_info.mx[which(is.na(gst_info.df)==TRUE)] <- ""
# gst_info.df <- as.data.frame(gst_info.mx)
