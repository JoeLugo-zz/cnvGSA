#' Prepare the files for the enrichment maps.
#'
#' @param cnvGSA.in A CnvGSAInput S4 object
#' @param cnvGSA.out A CnvGSAOutput S4 object
#' @return Returns a list with the data frames of the GMT file and the generic file. 

f.enrFiles <- function(cnvGSA.in,cnvGSA.out)
{
	config.df <- cnvGSA.in@config.ls$config.df

	Kl <- config.df[config.df$param == "Kl","value"]

	if (Kl == "ALL" || Kl == "")
	{
		f.enrProcess(cnvGSA.in,cnvGSA.out,Kl = "YES")
		f.enrProcess(cnvGSA.in,cnvGSA.out,Kl = "NO")
	} else {
		f.enrProcess(cnvGSA.in,cnvGSA.out,Kl)
	}
}

f.enrProcess <- function(cnvGSA.in,cnvGSA.out,Kl)
{
	t <- Sys.time()
	timestamp <- strftime(t,"%Y%m%d%Hh%Mm%S")

	config.df <- cnvGSA.in@config.ls$config.df

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
	res.ls        <- cnvGSA.out@res.ls

	if (pVal == "")             {pVal <- "Pvalue_U_dev"}
	if (FDR == "")              {FDR <- "FDR_BH_U"}
	if (coeff == "")            {coeff <- "Coeff_U"}
	if (keepCoeff == "")        {keepCoeff <- "YES"}
	if (filtGsEnr == "")        {filtGsEnr <- "NO"}
	if (is.na(minCaseCount))    {minCaseCount <- 0}
	if (is.na(minControlCount)) {minControlCount <- 0}
	if (is.na(minRatio))        {minRatio <- 0}

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
	gsGenes.ls        <- gsGenes.ls[which(!(names(gsGenes.ls) %in% diff))] # make sure all gene-sets in gmt file are in generic file

	f.collapse <- function (input.genes){
        paste (input.genes, collapse = "\t")
    }

	f.PackGMT <- function (id2eg.ls, id2des.chv, file.name){
        genes.chv  <- unlist(lapply (id2eg.ls, f.collapse))
        output.chv <- paste(names (id2eg.ls), id2des.chv, genes.chv, sep = "\t")
         
        cat (output.chv, sep = "\n", file = file.name)
    }

    setwd(outputPathEnr)
    cat(paste("Changing directory to ",outputPathEnr,sep=""))
    cat("\n")
	write.table(enrGeneric,file=paste("enrGeneric_KL",Kl,"_",timestamp,".txt",sep=""),row.names=FALSE,sep="\t",quote=FALSE)
	f.PackGMT(id2eg.ls = gsGenes.ls, id2des.chv = gs_info.df$GsName, file.name = paste("enr_KL",Kl,"_",timestamp,".gmt",sep=""))
}

# f.enrFiles(cnvGSA.in,cnvGSA.out)