## Changes from JL
# GP20171005

# fixed as.data.frame to as.data.frame.matrix:
# as.data.frame.matrix(geneCount.tab)
# as.data.frame.matrix(sid2gs_TYPE.tab)

# assigned default values for klOlp, geneSetSizeMin, geneSetSizeMax,
# parallel, eventThreshold, fLevels, cores, CNVevents
# creates outdir if not found

# 1.0. Libraries
library(GenomicRanges)
library(doParallel)
library(foreach)
library(splitstackshape)

# 2. Reading the config file
#' Reading in the config file.
#'
#' This function is used to read in all the values from the config file and
#' # change them in the S4 objects used throughout the scripts.
#' If you would like to reload the config values, you will want to run this function.
#'
#' @param configFile The file name and full path for the config file.
#' @param cnvGSA.in A CnvGSAInput S4 object.
#' @return A CnvGSAInput object with the updated config.ls and params.ls objects.
#' @examples
#' library(cnvGSAdata)
#' data(cnvGSA_input_example)
#' ## See vignette for full details and worked example

<<<<<<< HEAD
f.readConfig <- function(configFile, cnvGSA.in) {
    # read configuration file; if not found, warn user and stop
    if(missing(configFile)) {
        stop("Missing 'configFile' arguement")
    }
    
    config.df <- read.table(
        configFile,
        header = T,
        sep = "\t",
        quote = "\"",
        strip.white = TRUE,
        stringsAsFactors = F
    )
    
    cnvFile      <-
        config.df[config.df$param == "cnvFile",      "value"]
    phFile       <-
        config.df[config.df$param == "phFile",       "value"]
    geneIDFile   <-
        config.df[config.df$param == "geneIDFile",   "value"]
    klGeneFile   <-
        config.df[config.df$param == "klGeneFile",   "value"]
    klLociFile   <-
        config.df[config.df$param == "klLociFile",   "value"]
    gsFile       <-
        config.df[config.df$param == "gsFile",       "value"]
    outputPath   <-
        config.df[config.df$param == "outputPath",   "value"]
    geneListFile <-
        config.df[config.df$param == "geneListFile", "value"]
    Kl           <-
        config.df[config.df$param == "Kl",           "value"]
    projectName  <-
        config.df[config.df$param == "projectName",  "value"]
    gsUSet       <-
        config.df[config.df$param == "gsUSet",       "value"]
    cnvType      <-
        config.df[config.df$param == "cnvType",      "value"]
    covariates   <-
        gsub(" ", "", unlist(strsplit(config.df[config.df$param == "covariates",   "value"], split = ",")), fixed = TRUE)
    klOlp        <-
        as.numeric(config.df[config.df$param == "klOlp",         "value"])
    
    if(length(klOlp) == 0 || is.na(klOlp)) {
        klOlp <- 0.5
    }
    
    corrections  <-
        gsub(" ", "", unlist(strsplit(config.df[config.df$param == "corrections", "value"], split = ",")), fixed = TRUE)
    geneSep      <-
        config.df[config.df$param == "geneSep",      "value"]
    geneSetSizeMin <-
        as.numeric(config.df[config.df$param == "geneSetSizeMin", "value"])
    
    if(length(geneSetSizeMin) == 0 ||
       is.na(geneSetSizeMin)) {
        geneSetSizeMin <-
            25
    } # fixed empty assignment in case config has empty value GP20171006
    
    geneSetSizeMax <-
        as.numeric(config.df[config.df$param == "geneSetSizeMax", "value"])
    
    if(length(geneSetSizeMax) == 0 ||
       is.na(geneSetSizeMax)) {
        geneSetSizeMax <-
            1500
    } # fixed empty assignment in case config has empty value  GP20171006
    
    filtGs       <-
        config.df[config.df$param == "filtGs",       "value"]
    covInterest  <-
        config.df[config.df$param == "covInterest",  "value"]
    parallel     <-
        config.df[config.df$param == "parallel",     "value"]
    
    cores        <-
        max(as.numeric(config.df[config.df$param == "cores",     "value"]), 1) # fixed empty assignment in case config has empty value GP20171006
    if(length(parallel) == 0 || is.na(parallel)) {
        parallel <-
            'NO'
    } # fixed empty assignment in case config has empty value GP20171006
    
    eventThreshold <-
        max(as.numeric(config.df[config.df$param == "eventThreshold", "value"]), 1) # fixed empty assignment in case config has empty value GP20171006
    fLevels      <-
        as.numeric(config.df[config.df$param == "fLevels",       "value"])
    
    if(length(fLevels) == 0 || is.na(fLevels)) {
        fLevels <- 10
    }
    
    CNVevents    <-
        max(as.numeric(config.df[config.df$param == "CNVevents", "value"]), 1) # fixed empty assignment in case config has empty value GP20171006 - not sure if it should be 1 or 0
    
    config.ls <- list(
        cnvFile,
        phFile,
        geneIDFile,
        klGeneFile,
        klLociFile,
        gsFile,
        outputPath,
        geneListFile,
        config.df
    )
    
    params.ls <- list(
        Kl,
        projectName,
        gsUSet,
        cnvType,
        covariates,
        klOlp,
        corrections,
        geneSep,
        geneSetSizeMin,
        geneSetSizeMax,
        filtGs,
        covInterest,
        eventThreshold,
        fLevels,
        cores,
        parallel,
        CNVevents
    )
    
    names(config.ls) <- list(
        "cnvFile",
        "phFile",
        "geneIDFile",
        "klGeneFile",
        "klLociFile",
        "gsFile",
        "outputPath",
        "geneListFile",
        "config.df"
    )
    
    names(params.ls) <- list(
        "Kl",
        "projectName",
        "gsUSet",
        "cnvType",
        "covariates",
        "klOlp",
        "corrections",
        "geneSep",
        "geneSetSizeMin",
        "geneSetSizeMax",
        "filtGs",
        "covInterest",
        "eventThreshold",
        "fLevels",
        "cores",
        "parallel",
        "CNVevents"
    )
    
    cnvGSA.in@config.ls <- config.ls
    cnvGSA.in@params.ls <- params.ls
    
    # creating output path if it is not present GP20171006
    cat("Checking for the presence of", outputPath, ".\n")
    if(dir.exists(outputPath)) {
        cat("Directory", outputPath, "found.\n")
    } else {
        cat("Directory", outputPath, "not found.\n")
        cat("Creating results directory.\n")
        dir.create(outputPath)
    }
    return(cnvGSA.in)
=======
f.readConfig <- function(configFile,cnvGSA.in)
{
	if(missing(configFile)){
		stop("Missing 'configFile' arguement")
	}

	config.df <- read.table (configFile, header=T, sep="\t", quote="\"", stringsAsFactors=F)

	# CONFIG.LS
	cnvFile         <- config.df[config.df$param == "cnvFile","value"]
	phFile          <- config.df[config.df$param == "phFile","value"]
	geneIDFile      <- config.df[config.df$param == "geneIDFile","value"]
	klGeneFile      <- config.df[config.df$param == "klGeneFile","value"]
	klLociFile      <- config.df[config.df$param == "klLociFile","value"]
	gsFile          <- config.df[config.df$param == "gsFile","value"]
	outputPath      <- config.df[config.df$param == "outputPath","value"]
	geneListFile    <- config.df[config.df$param == "geneListFile","value"]

	# PARAMS.LS
	Kl              <- config.df[config.df$param == "Kl","value"]
	projectName     <- config.df[config.df$param == "projectName","value"]
	gsUSet          <- config.df[config.df$param == "gsUSet","value"]
	cnvType         <- config.df[config.df$param == "cnvType","value"]
	covariates      <- gsub(" ","",unlist(strsplit(config.df[config.df$param == "covariates","value"],split = ",")),fixed=TRUE)
	klOlp           <- as.numeric(config.df[config.df$param == "klOlp","value"])
	corrections     <- gsub(" ","",unlist(strsplit(config.df[config.df$param == "corrections","value"], split = ",")),fixed=TRUE)
	geneSep         <- config.df[config.df$param == "geneSep","value"]
	# keySep          <- config.df[config.df$param == "keySep","value"]
	geneSetSizeMin  <- as.numeric(config.df[config.df$param == "geneSetSizeMin","value"])
	geneSetSizeMax  <- as.numeric(config.df[config.df$param == "geneSetSizeMax","value"])
	filtGs          <- config.df[config.df$param == "filtGs","value"] 	
	covInterest     <- config.df[config.df$param == "covInterest","value"] 	
	parallel        <- config.df[config.df$param == "parallel","value"] 	
	eventThreshold  <- as.numeric(config.df[config.df$param == "eventThreshold","value"])
	fLevels         <- as.numeric(config.df[config.df$param == "fLevels","value"])
	cores           <- as.numeric(config.df[config.df$param == "cores","value"])
	CNVevents  		<- as.numeric(config.df[config.df$param == "CNVevents","value"])

	config.ls <- list(cnvFile, phFile, geneIDFile, klGeneFile, klLociFile, gsFile, outputPath, geneListFile, config.df)
	params.ls <- list(Kl, projectName, gsUSet, cnvType, covariates, klOlp, corrections, geneSep, geneSetSizeMin, geneSetSizeMax,filtGs,covInterest, eventThreshold, fLevels,cores,parallel,CNVevents)

	names(config.ls) <- list("cnvFile","phFile","geneIDFile","klGeneFile","klLociFile","gsFile","outputPath","geneListFile","config.df")
	names(params.ls) <- list("Kl","projectName","gsUSet","cnvType","covariates","klOlp","corrections","geneSep","geneSetSizeMin","geneSetSizeMax","filtGs","covInterest","eventThreshold","fLevels","cores","parallel","CNVevents")

	#if ("" %in% config.ls || "" %in% params.ls){
	#	warning("There are empty values in the config file")
	#}

	cnvGSA.in@config.ls <- config.ls
	cnvGSA.in@params.ls <- params.ls

	return(cnvGSA.in)
>>>>>>> efcba5556174499de61f2563571b970ef3b8b388
}

# 3. Read CNV + PhenoCovar and run checks
f.readData <- function(cnvGSA.in) {
    if(missing(cnvGSA.in)) {
        stop("Missing 'cnvGSA.in' arguement")
    }
    
    cat("Reading Data")
    cat("\n")
    
    config.ls <- cnvGSA.in@config.ls
    params.ls <- cnvGSA.in@params.ls
    
    # CNV DATA
    cnv.df     <-
        read.table(
            config.ls$cnvFile,
            header = T,
            sep = "\t",
            quote = "\"",
            stringsAsFactors = F
        )
    cnv.df$CHR <- as.character(cnv.df$CHR)
    
    if(!("SID" %in% colnames(cnv.df))) {
        stop("No SID column in the CNV data frame.")
    }
    
    length(unique(cnv.df$SID)) # 34257
    
    if(params.ls$geneSep == "") {
        params.ls$geneSep <- ";"
        cnvGSA.in@params.ls$geneSep <- ";"
    }
    
    geneID.ls       <-
        strsplit(cnv.df$geneID, split = params.ls$geneSep) # list of everything in that geneID
    geneID_temp.chv <-
        setdiff(unlist(geneID.ls), NA) # everything that isnt NA in the list only has the ones with numbers
    
    geneID.ls       <- lapply(geneID.ls, setdiff, "n/a")
    geneID_temp.chv <-
        setdiff(unlist(geneID.ls), c("n/a", NA)) # everything that doesnt have "n/a" or NA
    
    cnv.df$geneID   <-
        sapply(geneID.ls, paste, collapse = params.ls$geneSep) # produces "NA" as missing value, rather than NA
    cnv.df$geneID[which(cnv.df$geneID == "NA")] <-
        NA # makes all "NA" NA
    cnv.df$geneID[which(cnv.df$geneID == "")] <- NA # makes all "" NA
    
    geneID_temp.chv <-
        setdiff(unlist(strsplit(cnv.df$geneID, split = params.ls$geneSep)), NA)
    
    keySep <- "@"
    
    cnv.df$CnvKey   <-
        with(cnv.df, paste(CHR, BP1, BP2, TYPE, sep = keySep)) # new col with all of these combined
    
    # PHENOTYPE / COVARIATE DATA
    ph.df <-
        read.table(
            config.ls$phFile,
            header = T,
            sep = "\t",
            quote = "\"",
            stringsAsFactors = F
        )
    if(!("SID" %in% colnames(ph.df))) {
        ph.df$SID <-
            with(ph.df, paste(IID, FID, sep = keySep)) # makes the SID
        ph.df     <- subset(ph.df , select = -c(IID, FID, AFF))
    }
    
    # many duplicates once you get rid of these columns
    ph.df <- ph.df[!duplicated(ph.df),]
    
    phNames.vc <- c("SID", "Condition", params.ls$covariates)
    cnvNames.vc <- colnames(cnv.df)
    
    # drop unused columns
    # only want the SIDS that are in the CNV table so we can properly analyze it
    if(length(ph.df[, which(colnames(ph.df) %in% phNames.vc &
                            !(colnames(ph.df) %in% cnvNames.vc))]) != 0) {
        cnv.df <- merge(cnv.df,
                        ph.df[, which(colnames(ph.df) %in% phNames.vc &
                                          !(colnames(ph.df) %in% colnames(subset(
                                              cnv.df, select = -c(SID)
                                          ))))],
                        by = "SID", all = F) # combines them using the SID
    }
    
    kl_gene.df <-
        read.table(
            config.ls$klGeneFile,
            sep = "",
            header = T,
            comment.char = "",
            quote = "\"",
            stringsAsFactors = F
        )
    kl_loci.df <-
        read.table(
            config.ls$klLociFile,
            sep = "",
            header = T,
            comment.char = "",
            quote = "\"",
            stringsAsFactors = F
        )
    
    kl_loci.df$locuskey <-
        with(kl_loci.df, paste("KL", CHR, BP1, BP2, paste("T:", TYPE, sep = ""), sep = keySep))
    
    # 3.2. Read GeneSets
    load(config.ls$gsFile)
    
    # if(is.na(params.ls$geneSetSizeMin)) {		params.ls$geneSetSizeMin <- 25;
    # 	cnvGSA.in@params.ls$geneSetSizeMin <- 25}
    # if(is.na(params.ls$geneSetSizeMax)) {
    # 	params.ls$geneSetSizeMax <- 1500;	cnvGSA.in@params.ls$geneSetSizeMax <- 1500}
    
    if("U" %in% names(gs_all.ls)) {
        gs.ls         <- lapply(gs_all.ls, unique)
        gs.ls         <- lapply(gs.ls, setdiff, y = NA)
        gs_lengths.nv <- sapply(gs.ls, length)
        if(params.ls$filtGs == "YES") {
            gs.ls <-
                gs.ls[which(
                    gs_lengths.nv >= params.ls$geneSetSizeMin &
                        gs_lengths.nv <= params.ls$geneSetSizeMax
                )]
        }
        gs.ls$U <- gs_all.ls$U
        cat("Already universe set in the gene-set data")
        cat("\n")
    } else {
        gs.ls         <- lapply(gs_all.ls, unique)
        gs.ls         <- lapply(gs.ls, setdiff, y = NA)
        gs_lengths.nv <- sapply(gs.ls, length)
        if(params.ls$filtGs == "YES") {
            gs.ls <-
                gs.ls[which(
                    gs_lengths.nv >= params.ls$geneSetSizeMin &
                        gs_lengths.nv <= params.ls$geneSetSizeMax
                )]
        }
        if(params.ls$gsUSet == "") {
            cat("Using all genes in cnv data as universe set")
            cat("\n")
            gs.ls$U <- geneID_temp.chv
        } else {
            cat("Using specified universe set")
            cat("\n")
        }
    }
    
    gs_info.df <-
        data.frame(
            # making new data frame and these are the columns that are included
            GsKey  = paste("GS", 1:length(gs.ls), sep = ""),
            GsID   = names(gs.ls),
            GsName = gsid2name.chv[names(gs.ls)],
            GsSize = sapply(gs.ls, length),
            row.names = paste("GS", 1:length(gs.ls), sep = ""),
            stringsAsFactors = F
        )
    names(gs.ls) <-
        paste("GS", 1:length(gs.ls), sep = "") # making new label GS1,GS2,etc...
    if(!("U" %in% gs.ls) && params.ls$gsUSet != "") {
        gs_info.df[which(grepl(params.ls$gsUSet, gs_info.df$GsID)), ]$GsID <-
            "U"
    }
    
    if(length(rownames(gs_info.df[which(gs_info.df$GsID == "U"), ])) > 1) {
        stop("There seem to be multiple universe sets")
    }
    
    gs_sel_U.ls <- gs.ls
    
    # Making sure there is at least 10% of genes in the cnv data
    gs_key_U     <- gs_info.df[which(gs_info.df$GsID == "U"), ]$GsKey
    gs_gene_list <- subset(gs_all.ls, names(gs_all.ls) != gs_key_U)
    gs_gene_list <- unlist(as.list(cbind(gs_gene_list)))
    gs_gene_list <- unique(gs_gene_list)
    if(length(geneID_temp.chv) != 0) {
        perc_gs <-
            length(Reduce(intersect, list(gs_gene_list, geneID_temp.chv))) / length(geneID_temp.chv)
    }	else if(length(geneID_temp.chv) == 0) {
        perc_gs <- 0
    }
    
    if(perc_gs < .10) {
        warning(
            "The union of all genes from the gene-sets cover less than 10% of the union of all genes from the CNV's"
        )
    }
    
    # 3.3. Compute Overlap
    f.getPercOverlap <-
        function(cnv.start,
                 cnv.end,
                 loc.start,
                 loc.end)
        {
            loc.length <-
                loc.end - loc.start + 1 # why add 1? 1-positional 1-1 leads to length 0 but thats incorrect
            
            olp.start <- max(cnv.start, loc.start)
            olp.end   <- min(cnv.end, loc.end)
            
            olp.length <- olp.end - olp.start + 1
            olp.prc    <- olp.length / loc.length
            
            return(olp.prc)
        }
    
    # 3.4. Mark known loci
    cnv.gr <-
        GRanges(
            # creates class with single start and end point on the genome
            seqnames = Rle(paste(cnv.df$CHR, cnv.df$TYPE, sep = keySep)),
            ## match by chromosome and type
            ranges   = IRanges(start = cnv.df$BP1, end = cnv.df$BP2),
            strand   = Rle(strand(rep(
                "+", nrow(cnv.df)
            ))),
            chr      = cnv.df$CHR,
            type     = cnv.df$TYPE,
            cnvkey   = cnv.df$CnvKey
        )
    
    loci.gr <- GRanges(
        seqnames = Rle(paste(
            kl_loci.df$CHR, kl_loci.df$TYPE, sep = keySep
        )),
        ## match by chromosome and type
        ranges   = IRanges(start = kl_loci.df$BP1, end = kl_loci.df$BP2),
        strand   = Rle(strand(rep(
            "+", nrow(kl_loci.df)
        ))),
        chr      = kl_loci.df$CHR,
        type     = kl_loci.df$TYPE,
        locuskey = kl_loci.df$locuskey
    )
    
    o.hits <- findOverlaps(cnv.gr, loci.gr)
    o.df   <- as.data.frame(o.hits)
    
    # ** ADD ** is there a genomicranges operation to do this?
    olp_prc.nv <- mapply(
        f.getPercOverlap,
        start(ranges(cnv.gr)) [o.df[, "queryHits"]],
        end(ranges(cnv.gr)) [o.df[, "queryHits"]],
        start(ranges(loci.gr))[o.df[, "subjectHits"]],
        end(ranges(loci.gr))[o.df[, "subjectHits"]],
        SIMPLIFY = TRUE
    )
    
    # if(is.na(params.ls$klOlp)) { # fixed in config assignment, GP20171006
    # 	params.ls$klOlp <- 0.5
    # 	cnvGSA.in@params.ls$klOlp <- 0.5}
    
    o_olp50.mx     <-
        o.df[olp_prc.nv > params.ls$klOlp,] # if overlap % > 50 then count as overlap
    
    olp50.cnvkey   <- mcols(cnv.gr)$cnvkey[o_olp50.mx[, "queryHits"]]
    
    cnv.df$OlpKL_CNV <- 0	# no overlap
    cnv.df$OlpKL_CNV[which(cnv.df$CnvKey %in% olp50.cnvkey)] <-
        1 # some overlap
    
    olp50.sid <-
        subset(cnv.df,
               subset = OlpKL_CNV == 1,
               select = SID,
               drop = T) # finding those that have overlap in cnv data frame
    
    cnv.df$OlpKL_SID <- 0 # no overlap
    cnv.df$OlpKL_SID[which(cnv.df$SID %in% olp50.sid)] <-
        1 # some overlap
    
    # 3.4.1. By gene id
    cnv2.df <- cnv.df
    cnv2.df$SubjCnvKey <-
        with(cnv2.df, paste(SID, CnvKey, sep = paste(keySep, keySep, sep = "")))
    cnv2.df <- cnv2.df[!duplicated(cnv2.df$SubjCnvKey),]
    
    cnv2gene.ls <- strsplit(cnv2.df$geneID, params.ls$geneSep)
    names(cnv2gene.ls) <- cnv2.df$SubjCnvKey
    cnv2gene.df <-
        stack(cnv2gene.ls)
    names(cnv2gene.df) <- c("geneID", "SubjCnvKey")
    cnv2gene.df <-
        merge(cnv2gene.df, cnv2.df[, c("SubjCnvKey", "CnvKey", "TYPE")], by = "SubjCnvKey", all = T)
    
    cnv2gene.df <- subset(cnv2gene.df, subset = !is.na(geneID))
    cnv2gene.df$geneID_type <-
        with(cnv2gene.df, paste(geneID, TYPE, sep = keySep))
    kl_gene.df$geneID_type  <-
        with(kl_gene.df,  paste(geneID,   TYPE, sep = keySep))
    
    # any CNV that has this gene will me marked
    kg.cnvkey <-
        subset(
            cnv2gene.df,
            subset = geneID_type %in% kl_gene.df$geneID_type,
            select = CnvKey,
            drop = T
        )
    cnv.df$OlpKL_CNV[which(cnv.df$CnvKey %in% kg.cnvkey)] <- 1
    
    # 3.4.2. Mark subjects
    klg.sid <-
        subset(cnv.df,
               subset = OlpKL_CNV == 1,
               select = SID,
               drop = T)
    
    cnv.df$OlpKL_SID <- 0
    cnv.df$OlpKL_SID[which(cnv.df$SID %in% klg.sid)] <- 1
    
    # 3.5. phenotype table
    ph.df <-
        merge(ph.df, subset(cnv.df, select = c(SID, OlpKL_SID)), by = "SID")
    
    if(nlevels(as.factor(ph.df$SEX)) >= 20) {
        warning("There are more than 20 levels in the SEX factor")
    }
    if(nlevels(as.factor(ph.df$CNV_platform)) >= 20) {
        warning("There are more than 20 levels in the CNV_platform factor")
    }
    
    ph.df <- ph.df[!duplicated(ph.df),]
    
    cnv.df$SubjCnvKey <-
        with(cnv.df, paste(SID, CnvKey, sep = paste(keySep, keySep, sep = "")))
    
    # 3.6. main table
    check_type <- params.ls$cnvType
    type.vc    <- unique(cnv.df$TYPE)
    if(params.ls$cnvType == "ALL" || params.ls$cnvType == "") {
        check_type <- type.vc
    }
    
    cnv.df             <- cnv.df[!duplicated(cnv.df$SubjCnvKey),]
    cnv2gene.ls        <-
        strsplit(cnv.df$geneID, split = params.ls$geneSep)
    names(cnv2gene.ls) <- cnv.df$SubjCnvKey
    cnv2gene.df        <-
        stack(cnv2gene.ls)
    names(cnv2gene.df) <- c("geneID", "SubjCnvKey") # sets colnames
    cnv2gene.df        <-
        merge(cnv2gene.df, cnv.df[, unlist(c(
            params.ls$covariates,
            c("CHR", "BP1", "BP2", "SubjCnvKey", "SID", "TYPE")
        ))], by = "SubjCnvKey", all = T)
    
    cnv2gene.df$geneID_TYPE  <- cnv2gene.df$geneID
    if(params.ls$cnvType != "ALL") {
        cnv2gene.df$geneID_TYPE <-
            cnv2gene.df$geneID
        cnv2gene.df$geneID_TYPE[which(cnv2gene.df$TYPE != check_type)] <-
            NA
    }
    
    # splitting up data for only loss or only gain
    cnv2gene_TYPE.df <-
        subset(cnv2gene.df, select = c(SID, geneID_TYPE, SubjCnvKey))
    
    # making list of gene sets into data frame
    gs_sel_U.df <-
        stack(gs_sel_U.ls)
    names(gs_sel_U.df) <- c("geneID", "GsKey")
    gs_sel_U.df <-
        merge(gs_sel_U.df, subset(gs_info.df, select = c(GsKey, GsID, GsName)), by =
                  "GsKey")
    
    sid2gs_TYPE.df <-
        merge(
            cnv2gene_TYPE.df,
            gs_sel_U.df,
            by.x = "geneID_TYPE",
            by.y = "geneID",
            all.x = T,
            all.y = F
        )
    
    # this is what counts the events for each gene set
    sid2gs_TYPE.tab      <- table(sid2gs_TYPE.df[, c("SID", "GsKey")])
    
    gs_colnames_TYPE.chv <-
        colnames(sid2gs_TYPE.tab)[which(apply(sid2gs_TYPE.tab > 0, 2, sum) >= params.ls$CNVevents)]
    
    geneID.df <-
        read.table(
            config.ls$geneIDFile,
            header = T,
            sep = "\t",
            quote = "\"",
            stringsAsFactors = F
        ) # "cnv_AGP_demo.txt" "PGC_41K_QC_exon.cnv.annot"
    
    # parsing out gene-ids so only one gene per row
    gene2sid.df <-
        cSplit(cnv.df,
               splitCols = "geneID",
               sep = ";",
               direction = "long")
    
    # Applying Thresholds - Gene Count Table
    geneCount.tab <- table(gene2sid.df$geneID, gene2sid.df$Condition)
    geneCount.df  <-
        as.data.frame.matrix(geneCount.tab) # fixed; GP 20171005
    colnames(geneCount.df) <- c("Controls", "Cases")
    
    geneCount.df$geneID <- rownames(geneCount.df)
    row.names(geneCount.df) <- NULL
    
    geneCount.df <-
        merge(geneCount.df,
              subset(geneID.df, select = -c(Symbol)),
              by.x = "geneID",
              by.y = "geneID")
    geneCount.df <- geneCount.df[, c(1, 4, 2, 3)]
    
    
    # gets rid of these two data frames
    rm(sid2gs_TYPE.df)
    
    gc()
    gc()
    gc()
    
    sid2gs_TYPE_tab.df     <-
        as.data.frame.matrix(sid2gs_TYPE.tab) # fixed; GP 20171005
    
    
    sid2gs_TYPE_tab.df$SID <- rownames(sid2gs_TYPE.tab)
    gc()
    gc()
    
    cat("Building Covariates")
    cat("\n")
    
    cnv.df$CnvLength_ALL  <- with(cnv.df, BP2 - BP1 + 1)
    cnv.df$CnvLength_TYPE <- with(cnv.df, BP2 - BP1 + 1)
    # all cnvlength with type 1 specifies only looking for certain type not multiplying the numbers
    # only need to run this if they want the cnv type gain or loss returns 1 or 0 if true or false
    if(params.ls$cnvType != "ALL") {
        cnv.df$CnvLength_TYPE <-
            with(cnv.df, CnvLength_ALL * as.numeric(TYPE == check_type))
    }
    # all cnvlength with type 3
    
    cnv.df$CnvCount_TYPE <-
        with(cnv.df, as.numeric(TYPE %in% c(check_type)))
    cnv.df$geneID_TYPE   <- cnv.df$geneID
    cnv.df$geneID_TYPE[which(cnv.df$TYPE != check_type)] <- NA
    
    # 4.1. COVARIATES
    cnvc_TYPE.df  <-
        aggregate(formula = CnvCount_TYPE ~ SID,
                  data = cnv.df,
                  FUN = sum)
    
    ph.df         <-
        unique(merge(ph.df, cnvc_TYPE.df, all = T, by = "SID"))
    
    tlen_TYPE.df  <-
        aggregate(formula = CnvLength_TYPE ~ SID,
                  data = cnv.df,
                  FUN = sum)
    names(tlen_TYPE.df)[2] <- "CnvTotLength_TYPE"
    
    ph.df        <-
        unique(merge(ph.df, tlen_TYPE.df, all = T, by = "SID"))
    
    mlen_TYPE.df <-
        aggregate(formula = CnvLength_TYPE ~ SID,
                  data = cnv.df,
                  FUN = mean)
    names(mlen_TYPE.df)[2] <- "CnvMeanLength_TYPE"
    
    ph.df        <-
        unique(merge(ph.df, mlen_TYPE.df, all = T, by = "SID"))
    
    ph_TYPE.df   <-
        unique(merge(ph.df, sid2gs_TYPE_tab.df, all = T, by = "SID"))
    
    if(length(ph_TYPE.df$SID) != length(unique(ph_TYPE.df$SID))) {
        warning("There are duplicated SID's in the ph_TYPE.df. May want to check this")
    }
    
    cnvData <- list(cnv.df, cnv2gene.df)
    phData  <- list(ph.df, ph_TYPE.df)
    gsData  <-
        list(gs_info.df,
             gs_sel_U.df,
             gs_colnames_TYPE.chv,
             geneCount.tab,
             gs_all.ls)
    geneID  <- list(geneID.df)
    
    cnvGSA.in@cnvData.ls         <- cnvData
    names(cnvGSA.in@cnvData.ls)  <- list("cnv.df", "cnv2gene.df")
    
    cnvGSA.in@phData.ls          <- phData
    names(cnvGSA.in@phData.ls)   <- list("ph.df", "ph_TYPE.df")
    
    cnvGSA.in@gsData.ls          <- gsData
    names(cnvGSA.in@gsData.ls)   <-
        list(
            "gs_info.df",
            "gs_sel_U.df",
            "gs_colnames_TYPE.chv",
            "geneCount.tab",
            "gs_all.ls"
        )
    
    cnvGSA.in@params.ls$check_type <- check_type
    
    cnvGSA.in@geneID.ls            <- geneID
    names(cnvGSA.in@geneID.ls)   <- list("geneID.df")
    
    return(cnvGSA.in)
}

# 5. Creating output
#' Performing the logistic regression tests on the CNV data.
#'
#' This test uses 4 different correction models and requires a case control study. It looks at odds ratios and calculates statistics for the gene-set collection.
#'
#' @param cnvGSA.in A CnvGSAInput S4 object.
#' @param cnvGSA.out A CnvGSAOutput S4 object.
#' @return A list of one or two objects depending on whether or not the user includes the known loci in the analysis. Each object in the list contains the regression results for each gene set.
#' @examples
#' library(cnvGSAdata)
#' data(cnvGSA_output_example)
#' ## See vignette for full details and worked example

cnvGSAlogRegTest <- function(cnvGSA.in, cnvGSA.out)
    # master.ls,
{
    t <- Sys.time()
    timestamp <- strftime(t, "%Y%m%d%Hh%Mm%S")
    
    cat("Running Tests")
    cat("\n")
    
    if(missing(cnvGSA.in)) {
        stop("Missing 'cnvGSA.in' argument")
    }
    if(missing(cnvGSA.out)) {
        stop("Missing 'cnvGSA.out' argument")
    }
    
    gs_info.df 			   <- cnvGSA.in@gsData.ls$gs_info.df
    gs_colnames_TYPE.chv <- cnvGSA.in@gsData.ls$gs_colnames_TYPE.chv
    ph_TYPE.df			   <- cnvGSA.in@phData.ls$ph_TYPE.df
    ph.df				      <- cnvGSA.in@phData.ls$ph.df
    phData.ls			   <- cnvGSA.in@phData.ls
    
    config.ls <- cnvGSA.in@config.ls
    params.ls <- cnvGSA.in@params.ls
    
    # 5.1. Test
    # using unlist to take apart the covariates from the config file
    if(params.ls$covariates[1] == "NONE") {
        # "SEX,CNV_metric,CNV_platform,C1,C2,C3,C4,C8"
        params.ls$covariates <- ""
    }
    
    # if(length(params.ls$eventThreshold) == 0 || is.na(params.ls$eventThreshold)) { # fixed in config, GP20171006
    # 	params.ls$eventThreshold <- 1
    # 	cnvGSA.in@params.ls$eventThreshold <- 1
    # }
    #
    # if(length(params.ls$fLevels) == 0 || is.na(params.ls$fLevels)) { # fixed in config, GP20171006
    # 	params.ls$fLevels <- 10
    # 	cnvGSA.in@params.ls$fLevels <- 10
    # }
    
    setU.gskey <-
        subset(gs_info.df, GsID == "U", select = GsKey, drop = T)
    
    # finding those rows that don't have the specific GsKey(s)
    gs_colnames_TYPE.chv <- setdiff(gs_colnames_TYPE.chv, setU.gskey)
    
    if(config.ls$geneListFile != "") {
        gs_colnames_TYPE.chv <-
            unlist(strsplit(readLines(config.ls$geneListFile), ","))
    }
    
    if(params.ls$Kl == "YES") {
        #ph_TYPE.df <- ph_TYPE.df
        kl_fn <-
            paste(
                "GsTest_",
                params.ls$projectName,
                "_",
                params.ls$cnvType,
                "_KLy_",
                timestamp,
                ".txt",
                sep = ""
            )
        dataNames <-
            list(paste("covAll_chipAll_", params.ls$cnvType, "_KLy.df", sep = ""))
        cat("Kl - YES")
        cat("\n")
    } else if(params.ls$Kl == "NO") {
        ph_TYPE.df <- subset(ph_TYPE.df, OlpKL_SID == 0)
        kl_fn2 <-
            paste(
                "GsTest_",
                params.ls$projectName,
                "_",
                params.ls$cnvType,
                "_KLn_",
                timestamp,
                ".txt",
                sep = ""
            )
        dataNames <-
            list(paste("covAll_chipAll_", params.ls$cnvType, "_KLn.df", sep = ""))
        cat("Kl - NO")
        cat("\n")
    } else if(params.ls$Kl == "ALL" || params.ls$Kl == "") {
        #ph_TYPE.df <- ph_TYPE.df
        kl_fn  <-
            paste(
                "GsTest_",
                params.ls$projectName,
                "_",
                params.ls$cnvType,
                "_KLy_",
                timestamp,
                ".txt",
                sep = ""
            )
        kl_fn2 <-
            paste(
                "GsTest_",
                params.ls$projectName,
                "_",
                params.ls$cnvType,
                "_Kln_",
                timestamp,
                ".txt",
                sep = ""
            )
        dataNames <-
            list(
                paste("covAll_chipAll_", params.ls$cnvType, "_KLy.df", sep = ""),
                paste("covAll_chipAll_", params.ls$cnvType, "_KLn.df", sep = "")
            )
        cat("KL - ALL : ")
        cat("The function will run twice for results with AND without the known loci")
        cat("\n")
    } else {
        stop("Invalid Kl value")
    }
    
    totalLenGs <- length(gs_colnames_TYPE.chv)
    
    # 5.1. Ancillary functions
    # looks at the output of the analysis
    f.getCoeff_sm <- function(glm.sm, var.ch) {
        if(var.ch %in% rownames(glm.sm$coefficients)) {
            return(glm.sm$coefficients[var.ch, "Estimate"])
        }
        else {
            return(NA)
        }
    }
    f.getPval_sm  <- function(glm.sm, var.ch) {
        if(var.ch %in% rownames(glm.sm$coefficients)) {
            return(glm.sm$coefficients[var.ch, "Pr(>|z|)"])
        }
        else {
            return(NA)
        }
    }
    f.getPval_anova  <- function(glm.anova, var.ch) {
        if(var.ch %in% rownames(glm.anova["Pr(>Chi)"])) {
            return(glm.anova["Pr(>Chi)"][var.ch, "Pr(>Chi)"])
        }
        else {
            return(NA)
        }
    }
    
    f.testGLM_wrap <-
        function(gs.colnames,
                 data.df,
                 covar.chv,
                 u.gskey,
                 covInterest,
                 eventThreshold,
                 fLevels,
                 cores) {
            # CASES
            data_sz.ix <-
                which(data.df$Condition == 1) # which rowws meet this requirement
            data_sz.df <- subset(data.df, Condition == 1) # cases
            s_sz.n <- length(data_sz.ix)
            # CONTROLS
            data_ct.ix <- which(data.df$Condition == 0)
            data_ct.df <- subset(data.df, Condition == 0) # controls
            s_ct.n <- length(data_ct.ix)
            
            data.df[, covInterest] <- as.factor(data.df[, covInterest])
            lev.ls <- levels(data.df[, covInterest])
            cat("Including", covInterest, "in test", '\n')
            
            if(length(lev.ls) > fLevels) {
                stop(paste(
                    "Number of fLevels in ",
                    covInterest,
                    " exceeds ",
                    fLevels,
                    sep = " "
                ))
            }
            
            if(length(params.ls$corrections) == 0) {
                corrections <- c("uni_gc")
            }
            
            cat("Using the following corrections:",
                params.ls$corrections,
                "\n")
            
            coreNum <- detectCores()
            
            if(is.na(cores)) {
                cores <- coreNum
            }
            
            # t returns the transpose of a matrix
            if(length(params.ls$parallel) == 0 ||
               params.ls$parallel == "NO") {
                res.mx <- t(
                    sapply(
                        gs.colnames,
                        f.testGLM_unit,
                        data.df = data.df,
                        covar.chv = covar.chv,
                        u.gskey = u.gskey,
                        sz.ix = data_sz.ix,
                        ct.ix = data_ct.ix,
                        correct.ls = params.ls$corrections,
                        covInterest = covInterest,
                        eventThreshold = eventThreshold,
                        data_sz.df = data_sz.df,
                        data_ct.df = data_ct.df,
                        s_sz.n = s_sz.n,
                        s_ct.n = s_ct.n
                    )
                )
                res.df <- as.data.frame(res.mx)
                res.df$GsKey <- gs.colnames
            }
            else {
                if(cores > coreNum) {
                    cores <- coreNum
                    cat(
                        paste(
                            "Cores specified exceeds number of cores detected. Only using ",
                            coreNum,
                            " cores.",
                            sep = ""
                        )
                    )
                    cat("\n")
                }
                registerDoParallel(cores = cores)
                cat(paste0("Using ", cores, " cores for parallelization of tests"))
                cat("\n")
                res.mx <-
                    t(as.data.frame(
                        foreach(i = 1:length(gs.colnames)) %dopar%
                            f.testGLM_unit(
                                gs.colnames[i],
                                data.df = data.df,
                                covar.chv = covar.chv,
                                u.gskey = u.gskey,
                                sz.ix = data_sz.ix,
                                ct.ix = data_ct.ix,
                                correct.ls = params.ls$corrections,
                                covInterest = covInterest,
                                eventThreshold = eventThreshold,
                                data_sz.df = data_sz.df,
                                data_ct.df = data_ct.df,
                                s_sz.n = s_sz.n,
                                s_ct.n = s_ct.n
                            )
                    ))
                res.df       <- as.data.frame(res.mx)
                res.df$GsKey <- gs.colnames
                row.names(res.df) <- NULL
            }
            
            colnames(res.df) <-
                c(
                    "Coeff",
                    "Pvalue_glm",
                    "Pvalue_dev",
                    "Pvalue_dev_s",
                    "Coeff_U",
                    "Pvalue_U_glm",
                    "Pvalue_U_dev",
                    "Pvalue_U_dev_s",
                    "Coeff_TL",
                    "Pvalue_TL_glm",
                    "Pvalue_TL_dev",
                    "Pvalue_TL_dev_s",
                    "Coeff_CNML",
                    "Pvalue_CNML_glm",
                    "Pvalue_CNML_dev",
                    "Pvalue_CNML_dev_s",
                    "CASE_g1n",
                    "CTRL_g1n",
                    "CASE_g2n",
                    "CTRL_g2n",
                    "CASE_g3n",
                    "CTRL_g3n",
                    "CASE_g4n",
                    "CTRL_g4n",
                    "CASE_g5n",
                    "CTRL_g5n",
                    "CASE_gTT",
                    "CTRL_gTT",
                    paste(c("CASE"), lev.ls, sep = "_"),
                    paste(c("CTRL"), lev.ls, sep = "_"),
                    "GsKey"
                )
            
            return(res.df)
        }
    
    # 5.2. LOGISTIC REGRESSION TEST
    f.testGLM_unit <-
        function(gs.colname,
                 data.df,
                 covar.chv,
                 u.gskey,
                 sz.ix,
                 ct.ix,
                 correct.ls,
                 covInterest,
                 eventThreshold,
                 data_sz.df,
                 data_ct.df,
                 s_sz.n,
                 s_ct.n) {
            cat(match(gs.colname, gs_colnames_TYPE.chv),
                " out of ",
                totalLenGs)
            output.nv <- numeric()
            lev.ls    <- levels(data.df[, covInterest])
            
            # no correction
            coeff <-
                NA
            pvalue_glm <- NA
            pvalue_dev <- NA
            pvalue_dev_s <- NA
            
            if("no_corr" %in% correct.ls ||
               "ALL" %in% correct.ls || length(correct.ls) == 0) {
                glm_form.ch  <-
                    paste("Condition",
                          "~",
                          paste(covar.chv, collapse = " + "),
                          "+",
                          gs.colname,
                          sep = " ")
                x.glm        <-
                    glm(as.formula(glm_form.ch), data.df, family = binomial(logit))
                x.glm_sm     <- summary(x.glm)
                x.anova      <- anova(x.glm, test = "Chisq")
                coeff        <- f.getCoeff_sm(x.glm_sm, gs.colname)
                pvalue_glm   <- f.getPval_sm(x.glm_sm, gs.colname)
                pvalue_dev   <- f.getPval_anova(x.anova, gs.colname)
                pvalue_dev_s <-
                    -log10(f.getPval_anova(x.anova, gs.colname)) * sign(f.getCoeff_sm(x.glm_sm, gs.colname))
            }
            
            # CORRECTION MODEL: universe count
            coeff_U <-
                NA
            pvalue_U_glm <- NA
            pvalue_U_dev <- NA
            pvalue_U_dev_s <- NA
            
            if("uni_gc" %in% correct.ls ||
               "ALL" %in% correct.ls || length(correct.ls) == 0) {
                # formula for the regression model
                # response variable ~ predictor variables
                glm_U_form.ch  <-
                    paste(
                        "Condition",
                        "~",
                        paste(covar.chv, collapse = " + "),
                        "+",
                        u.gskey,
                        "+",
                        gs.colname,
                        sep = " "
                    )
                x_U.glm        <-
                    glm(as.formula(glm_U_form.ch),
                        data.df,
                        family = binomial(logit))
                x_U.glm_sm     <- summary(x_U.glm)
                x_U.anova      <- anova(x_U.glm, test = "Chisq")
                coeff_U        <- f.getCoeff_sm(x_U.glm_sm, gs.colname)
                pvalue_U_glm   <- f.getPval_sm(x_U.glm_sm, gs.colname)
                pvalue_U_dev   <- f.getPval_anova(x_U.anova, gs.colname)
                pvalue_U_dev_s <-
                    -log10(f.getPval_anova(x_U.anova, gs.colname)) * sign(f.getCoeff_sm(x_U.glm_sm, gs.colname))
            }
            
            # CORRECTION MODEL: total length
            coeff_TL <-
                NA
            pvalue_TL_glm <- NA
            pvalue_TL_dev <- NA
            pvalue_TL_dev_s <- NA
            
            if("tot_l" %in% correct.ls ||
               "ALL" %in% correct.ls || length(correct.ls) == 0) {
                glm_TL_form.ch  <-
                    paste(
                        "Condition",
                        "~",
                        paste(covar.chv, collapse = " + "),
                        "+",
                        paste("CnvTotLength", "TYPE", sep = "_") ,
                        "+",
                        gs.colname,
                        sep = " "
                    )
                x_TL.glm        <-
                    glm(as.formula(glm_TL_form.ch),
                        data.df,
                        family = binomial(logit))
                x_TL.glm_sm     <- summary(x_TL.glm)
                x_TL.anova      <- anova(x_TL.glm, test = "Chisq")
                coeff_TL        <- f.getCoeff_sm(x_TL.glm_sm, gs.colname)
                pvalue_TL_glm   <- f.getPval_sm(x_TL.glm_sm, gs.colname)
                pvalue_TL_dev   <- f.getPval_anova(x_TL.anova, gs.colname)
                pvalue_TL_dev_s <-
                    -log10(f.getPval_anova(x_TL.anova, gs.colname)) * sign(f.getCoeff_sm(x_TL.glm_sm, gs.colname))
            }
            
            # CORRECTION MODEL: mean length	and number
            coeff_CNML <-
                NA
            pvalue_CNML_glm <-
                NA
            pvalue_CNML_dev <- NA
            pvalue_CNML_dev_s <- NA
            
            if("cnvn_ml" %in% correct.ls ||
               "ALL" %in% correct.ls || length(correct.ls) == 0) {
                glm_CNML_form.ch  <-
                    paste(
                        "Condition",
                        "~",
                        paste(covar.chv, collapse = " + "),
                        "+",
                        paste("CnvCount", "TYPE", sep = "_"),
                        "+",
                        paste("CnvMeanLength", "TYPE", sep = "_"),
                        "+",
                        gs.colname,
                        sep = " "
                    )
                x_CNML.glm        <-
                    glm(as.formula(glm_CNML_form.ch),
                        data.df,
                        family = binomial(logit))
                x_CNML.glm_sm     <- summary(x_CNML.glm)
                x_CNML.anova      <- anova(x_CNML.glm, test = "Chisq")
                coeff_CNML        <- f.getCoeff_sm(x_CNML.glm_sm, gs.colname)
                pvalue_CNML_glm   <- f.getPval_sm(x_CNML.glm_sm, gs.colname)
                pvalue_CNML_dev   <- f.getPval_anova(x_CNML.anova, gs.colname)
                pvalue_CNML_dev_s <-
                    -log10(f.getPval_anova(x_CNML.anova, gs.colname)) * sign(f.getCoeff_sm(x_CNML.glm_sm, gs.colname))
            }
            
            # sz = case , ct = control
            set_gn_sz.nv <- data.df[sz.ix, gs.colname]
            set_gn_ct.nv <- data.df[ct.ix, gs.colname]
            
            sz.ls <-
                unlist(lapply(lev.ls, function(x)
                    nrow(
                        subset(data_sz.df, subset = set_gn_sz.nv >= eventThreshold &
                                   data_sz.df[, covInterest] == x)
                    ) / nrow(
                        subset(data_sz.df, subset = data_sz.df[, covInterest] == x)
                    ) * 100))
            ct.ls <-
                unlist(lapply(lev.ls, function(x)
                    nrow(
                        subset(data_ct.df, subset = set_gn_ct.nv >= eventThreshold &
                                   data_ct.df[, covInterest] == x)
                    ) / nrow(
                        subset(data_ct.df, subset = data_ct.df[, covInterest] == x)
                    ) * 100))
            
            output.nv <- c(
                coeff,
                pvalue_glm,
                pvalue_dev,
                pvalue_dev_s,
                coeff_U,
                pvalue_U_glm,
                pvalue_U_dev,
                pvalue_U_dev_s,
                coeff_TL,
                pvalue_TL_glm,
                pvalue_TL_dev,
                pvalue_TL_dev_s,
                coeff_CNML,
                pvalue_CNML_glm,
                pvalue_CNML_dev,
                pvalue_CNML_dev_s,
                sum(set_gn_sz.nv >= 1) / s_sz.n * 100,
                sum(set_gn_ct.nv >= 1) / s_ct.n * 100,
                sum(set_gn_sz.nv >= 2) / s_sz.n * 100,
                sum(set_gn_ct.nv >= 2) / s_ct.n * 100,
                sum(set_gn_sz.nv >= 3) / s_sz.n * 100,
                sum(set_gn_ct.nv >= 3) / s_ct.n * 100,
                sum(set_gn_sz.nv >= 4) / s_sz.n * 100,
                sum(set_gn_ct.nv >= 4) / s_ct.n * 100,
                sum(set_gn_sz.nv >= 5) / s_sz.n * 100,
                sum(set_gn_ct.nv >= 5) / s_ct.n * 100,
                s_sz.n,
                s_ct.n,
                sz.ls[1:length(sz.ls)],
                ct.ls[1:length(ct.ls)]
            )
            
            cat("\n")
            
            rm(x.glm,
               x.glm_sm,
               x.anova,
               x_U.glm,
               x_U.glm_sm,
               x_U.anova)
            gc()
            gc()
            gc()
            
            return(output.nv)
        }
    
    res.ls <- list()
    
    {
        if(params.ls$Kl == "YES" ||
           params.ls$Kl == "ALL" || params.ls$Kl == "") {
            cat("Running test with known loci")
            cat("\n")
            res.ls$covAll_chipAll_TYPE_KLy.df <-
                f.testGLM_wrap(
                    gs.colnames = gs_colnames_TYPE.chv,
                    data.df = ph_TYPE.df,
                    covar.chv = params.ls$covariates,
                    u.gskey = setU.gskey,
                    covInterest = params.ls$covInterest,
                    eventThreshold = params.ls$eventThreshold,
                    fLevels = params.ls$fLevels,
                    cores = params.ls$cores
                )
            gc()
            gc()
            gc()
        }
        # Looking at all rows where there is no overlap according to the SID
        if(params.ls$Kl == "NO" ||
           params.ls$Kl == "ALL" || params.ls$Kl == "") {
            cat("Running test without known loci")
            cat("\n")
            res.ls$covAll_chipAll_TYPE_KLn.df <-
                f.testGLM_wrap(
                    gs.colnames = gs_colnames_TYPE.chv,
                    data.df = subset(ph_TYPE.df, OlpKL_SID == 0),
                    covar.chv = params.ls$covariates,
                    u.gskey = setU.gskey,
                    covInterest = params.ls$covInterest,
                    eventThreshold = params.ls$eventThreshold,
                    fLevels = params.ls$fLevels,
                    cores = params.ls$cores
                )
            gc()
            gc()
            gc()
        }
    }
    
    ph_TYPE.df[, params.ls$covInterest] <-
        as.factor(ph_TYPE.df[, params.ls$covInterest])
    lev.ls <- levels(ph_TYPE.df[, params.ls$covInterest])
    colOrder <-
        c(
            "GsID",
            "GsName",
            "GsSize",
            "Coeff",
            "Pvalue_glm",
            "Pvalue_dev",
            "Pvalue_dev_s",
            "FDR_BH",
            "Coeff_U",
            "Pvalue_U_glm",
            "Pvalue_U_dev",
            "Pvalue_U_dev_s",
            "FDR_BH_U",
            "Coeff_TL",
            "Pvalue_TL_glm",
            "Pvalue_TL_dev",
            "Pvalue_TL_dev_s",
            "FDR_BH_TL",
            "Coeff_CNML",
            "Pvalue_CNML_glm",
            "Pvalue_CNML_dev",
            "Pvalue_CNML_dev_s",
            "FDR_BH_CNML",
            "CASE_g1n",
            "CTRL_g1n",
            "CASE_g2n",
            "CTRL_g2n",
            "CASE_g3n",
            "CTRL_g3n",
            "CASE_g4n",
            "CTRL_g4n",
            "CASE_g5n",
            "CTRL_g5n",
            "CASE_gTT",
            "CTRL_gTT",
            paste(c("CASE"), lev.ls, sep = "_"),
            paste(c("CTRL"), lev.ls, sep = "_")
        )
    
    if(params.ls$Kl == "YES" ||
       params.ls$Kl == "ALL" || params.ls$Kl == "") {
        res.ls$covAll_chipAll_TYPE_KLy.df$FDR_BH      <-
            p.adjust(res.ls$covAll_chipAll_TYPE_KLy.df$Pvalue_dev,
                     method = "BH")
        res.ls$covAll_chipAll_TYPE_KLy.df$FDR_BH_U    <-
            p.adjust(res.ls$covAll_chipAll_TYPE_KLy.df$Pvalue_U_dev,
                     method = "BH")
        res.ls$covAll_chipAll_TYPE_KLy.df$FDR_BH_TL   <-
            p.adjust(res.ls$covAll_chipAll_TYPE_KLy.df$Pvalue_TL_dev,
                     method = "BH")
        res.ls$covAll_chipAll_TYPE_KLy.df$FDR_BH_CNML <-
            p.adjust(res.ls$covAll_chipAll_TYPE_KLy.df$Pvalue_CNML_dev,
                     method = "BH")
        res.ls$covAll_chipAll_TYPE_KLy.df <-
            merge(
                gs_info.df,
                res.ls$covAll_chipAll_TYPE_KLy.df,
                all.x = F,
                all.y = T,
                by = "GsKey"
            )
        res.ls$covAll_chipAll_TYPE_KLy.df <-
            res.ls$covAll_chipAll_TYPE_KLy.df[order(res.ls$covAll_chipAll_TYPE_KLy.df$Pvalue_U_dev_s,
                                                    decreasing = T),]
        res.ls$covAll_chipAll_TYPE_KLy.df <-
            subset(res.ls$covAll_chipAll_TYPE_KLy, select = -c(GsKey))
        res.ls$covAll_chipAll_TYPE_KLy.df <-
            res.ls$covAll_chipAll_TYPE_KLy.df[, colOrder]
    }
    
    if(params.ls$Kl == "NO" ||
       params.ls$Kl == "ALL" || params.ls$Kl == "") {
        res.ls$covAll_chipAll_TYPE_KLn.df$FDR_BH      <-
            p.adjust(res.ls$covAll_chipAll_TYPE_KLn.df$Pvalue_dev,
                     method = "BH")
        res.ls$covAll_chipAll_TYPE_KLn.df$FDR_BH_U    <-
            p.adjust(res.ls$covAll_chipAll_TYPE_KLn.df$Pvalue_U_dev,
                     method = "BH")
        res.ls$covAll_chipAll_TYPE_KLn.df$FDR_BH_TL   <-
            p.adjust(res.ls$covAll_chipAll_TYPE_KLn.df$Pvalue_TL_dev,
                     method = "BH")
        res.ls$covAll_chipAll_TYPE_KLn.df$FDR_BH_CNML <-
            p.adjust(res.ls$covAll_chipAll_TYPE_KLn.df$Pvalue_CNML_dev,
                     method = "BH")
        res.ls$covAll_chipAll_TYPE_KLn.df <-
            merge(
                gs_info.df,
                res.ls$covAll_chipAll_TYPE_KLn.df,
                all.x = F,
                all.y = T,
                by = "GsKey"
            )
        res.ls$covAll_chipAll_TYPE_KLn.df <-
            res.ls$covAll_chipAll_TYPE_KLn.df[order(res.ls$covAll_chipAll_TYPE_KLn.df$Pvalue_U_dev_s,
                                                    decreasing = T),]
        res.ls$covAll_chipAll_TYPE_KLn.df <-
            subset(res.ls$covAll_chipAll_TYPE_KLn, select = -c(GsKey))
        res.ls$covAll_chipAll_TYPE_KLn.df <-
            res.ls$covAll_chipAll_TYPE_KLn.df[, colOrder]
    }
    
    names(res.ls) <- dataNames
    
    setwd(config.ls$outputPath)
    cat("Changing directory to", config.ls$outputPath, "\n")
    
    if(params.ls$Kl == "YES" ||
       params.ls$Kl == "ALL" || params.ls$Kl == "") {
        resKLy <-
            get(paste("covAll_chipAll_", params.ls$cnvType, "_KLy.df", sep = ""),
                res.ls)
        write.table(
            resKLy,
            col.names = T,
            row.names = F,
            quote = F,
            sep = "\t",
            file = kl_fn
        )
    }
    if(params.ls$Kl == "NO" ||
       params.ls$Kl == "ALL" || params.ls$Kl == "") {
        resKLn <-
            get(paste("covAll_chipAll_", params.ls$cnvType, "_KLn.df", sep = ""),
                res.ls)
        write.table(
            resKLn,
            col.names = T,
            row.names = F,
            quote = F,
            sep = "\t",
            file = kl_fn2
        )
    }
    
    cnvGSA.out@res.ls <- res.ls
    cnvGSA.out@phData.ls <- phData.ls
    cnvGSA.out@gsData.ls <- cnvGSA.in@gsData.ls
    
    return(cnvGSA.out)
}

# 6. Creating gsTables
#' Creates the gene-set tables for each gene-set.
#'
#' Creates the gene-set tables for each gene-set.
#'
#' @param cnvGSA.in A CnvGSAInput S4 object.
#' @param cnvGSA.out A CnvGSAOutput S4 object.
#' @return A list where each object is a table corresponding to one gene-set.
#' @examples
#' library(cnvGSAdata)
#' data(cnvGSA_output_example)
#' ## See vignette for full details and worked example

cnvGSAgsTables <- function(cnvGSA.in, cnvGSA.out) {
    if (missing(cnvGSA.in)) {
        stop("Missing 'cnvGSA.in' arguement")
    }
    if (missing(cnvGSA.out)) {
        stop("Missing 'cnvGSA.out' arguement")
    }
    cat("Creating gsTables.ls")
    cat("\n")
    cnv2gene.df <- cnvGSA.in@cnvData.ls$cnv2gene.df
    gs_sel_U.df <- cnvGSA.in@gsData.ls$gs_sel_U.df
    geneID.df   <- cnvGSA.in@geneID.ls$geneID.df
    geneSep     <- cnvGSA.in@params.ls$geneSep
    params.ls   <- cnvGSA.in@params.ls
    cores       <- params.ls$cores 
    # coreNum   <- detectCores() # see below, GP20171006
    
    cnv2gene_TYPE.df <- subset(cnv2gene.df, select = c(SID, geneID_TYPE, SubjCnvKey))
    sid2gs_TYPE.df   <- merge(cnv2gene_TYPE.df, gs_sel_U.df, by.x = "geneID_TYPE", 
                              by.y = "geneID", all.x = T, all.y = F)
    
    if (params.ls$parallel == "NO") { # not using parallel
        #   original code  
        #   genes.ls <- strsplit(cnv.df$geneID,geneSep)
        #	 geneLen.chv <- 1:length(genes.ls)
        #	 geneSymbol.ls <- lapply(geneLen.chv,function(x) paste(geneID.df[which(geneID.df$geneID %in% unlist(genes.ls[x])),]$Symbol,collapse=geneSep))
        #	 cnv.df$Symbol <- geneSymbol.ls
        #
        #	 genes_TYPE.ls <- strsplit(cnv.df$geneID_TYPE,geneSep)
        #	 geneLen_TYPE.chv <- 1:length(genes_TYPE.ls)
        #	 geneSymbol_TYPE.ls <- lapply(geneLen_TYPE.chv,function(x) paste(geneID.df[which(geneID.df$geneID %in% unlist(genes_TYPE.ls[x])),]$Symbol,collapse=geneSep))
        #	 cnv.df$Symbol_TYPE <- geneSymbol_TYPE.ls
        
        # fixed code, GP20171006
        gs.df <- aggregate(geneID_TYPE ~ SID + SubjCnvKey + GsKey + GsID + GsName, data = sid2gs_TYPE.df, FUN = paste, collapse = geneSep)
        
        names(gs.df)[6] <- "geneID_GS"
        gsTables.ls <- split(gs.df, gs.df$GsID)
        
        cnv.df <- cnvGSA.in@cnvData.ls$cnv.df
        resTables.ls <- lapply(gsTables.ls, function(x) merge(cnv.df, x[, c("SubjCnvKey", "GsKey", "GsID", "GsName", "geneID_GS")],  all.y = T, by = "SubjCnvKey"))
        
    }  else { # using parallel
        
        coreNum <- detectCores() # moved from above, only makes sense if parallel is ON, GP20171006
        
        if (cores > coreNum) { # check number of cores
            cores <- coreNum
            cat(paste("Cores specified exceeds number of cores detected. Only using ", 
                      coreNum, " cores.", sep = ""))
            cat("\n")
        }
        
        registerDoParallel(cores = cores)
        cat(paste("Using ", cores, " cores for parallelization of tests", 
                  sep = ""))
        cat("\n")
        cnv.df <- cnvGSA.in@cnvData.ls$cnv.df
        genes.ls <- strsplit(cnv.df$geneID, geneSep)
        geneLen.chv <- 1:length(genes.ls)
        
        geneSymbol.ls <- foreach(i = 1:length(genes.ls)) %dopar% 
            unlist(paste(geneID.df[which(geneID.df$geneID %in% unlist(genes.ls[i])), ]$Symbol, collapse = geneSep))
        
        cnv.df$Symbol <- geneSymbol.ls
        genes_TYPE.ls <- strsplit(cnv.df$geneID_TYPE, geneSep)
        geneLen_TYPE.chv <- 1:length(genes_TYPE.ls)
        geneSymbol_TYPE.ls <- foreach(i = 1:length(genes_TYPE.ls)) %dopar% 
            unlist(paste(geneID.df[which(geneID.df$geneID %in% unlist(genes_TYPE.ls[i])), ]$Symbol, collapse = geneSep))
        cnv.df$Symbol_TYPE <- geneSymbol_TYPE.ls
    } # end using parallel
    
    # code below is the same for parallel or not
    gsTable_TYPE.df <- merge(subset(cnv.df, select = c(SID, 
                                                       CHR, BP1, BP2, TYPE, geneID, geneID_TYPE, Symbol, Symbol_TYPE)), 
                             subset(sid2gs_TYPE.df, select = -c(GsKey, geneID_TYPE, 
                                                                SubjCnvKey)), by = "SID")
    gsTable_TYPE.df        <- gsTable_TYPE.df[!duplicated(gsTable_TYPE.df), ]
    gsTable_TYPE.df$Symbol <- as.character(gsTable_TYPE.df$Symbol)
    gsTable_TYPE.df$Symbol_TYPE <- as.character(gsTable_TYPE.df$Symbol_TYPE)
    gsTable_TYPE.df$Symbol[which(gsTable_TYPE.df$Symbol == "")] <- NA
    gsTable_TYPE.df$Symbol_TYPE[which(gsTable_TYPE.df$Symbol_TYPE == 
                                          "")] <- NA
    gsTables.ls <- split(gsTable_TYPE.df, gsTable_TYPE.df$GsID)
    cnvGSA.out@gsTables.ls <- gsTables.ls
    
    return(cnvGSA.out)
}
# 7. Creating cnvGSA.in S4 object
#' Creating the input S4 object needed to run the script.
#'
#' @param configFile The file name for the config file including the full path.
#' @param cnvGSA.in A CnvGSAInput S4 object.
#' @return A cnvGSAInput S4 object.
#' @examples
#' library(cnvGSAdata)
#' data(cnvGSA_input_example)
#' ## See vignette for full details and worked example

cnvGSAIn <- function(configFile, cnvGSA.in) {
    cnvGSA.in <- f.readConfig(configFile, cnvGSA.in)
    cnvGSA.in <- f.readData(cnvGSA.in)
    return(cnvGSA.in)
}

# cnvGSA.in <- CnvGSAInput()
# cnvGSA.in <- cnvGSAIn(configFile = "/Users/josephlugo/Documents/R/PGC2_test/R_Works/PGC2_config.txt",cnvGSA.in)
# cnvGSA.in <- cnvGSAIn(configFile = "/Users/josephlugo/Documents/R/MehdiData/configFile.txt",cnvGSA.in)
# cnvGSA.out <- CnvGSAOutput()
# cnvGSA.out <- cnvGSAlogRegTest(cnvGSA.in,cnvGSA.out)