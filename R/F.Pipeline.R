
F.Pipeline <- function (cnvData.ls, gsData.ls, geneData.ls, params.ls)
{
	assTestPar.ls <- list()
	addEnrPar.ls <- list()
	assTestPar.ls$test_type    <- params.ls$grandtotals_mode
	assTestPar.ls$test_classes <- params.ls$sample_classes
	assTestPar.ls$iter         <- params.ls$fdr_iter
	assTestPar.ls$boxplot.bn   <- params.ls$boxplot_PDFs
	addEnrPar.ls$sel_n         <- params.ls$extended_report	

	cnvData.ls  <- F.PreProcess (cnvData.ls, gsData.ls)
	burdenSample.ls  <- F.BurdenSample (cnvData.ls, assTestPar.ls)
	multi.ls    <- F.GsCounts (cnvData.ls, gsData.ls)
	cnvData.ls  <- multi.ls$cnvData
	burdenGs.ls <- multi.ls$burdGsStat
	enrRes.ls   <- F.AssociationTest (cnvData.ls, gsData.ls, assTestPar.ls)
	geneData.ls <- F.ProcessGenes (cnvData.ls, gsData.ls, assTestPar.ls, geneData.ls)
	enrRes.ls   <- F.AddEnr (enrRes.ls, cnvData.ls, gsData.ls, geneData.ls, addEnrPar.ls, assTestPar.ls)

	cnvGSA.out <- new( "CnvGSAOutput", 
		cnvData = cnvData.ls,
		burdenSample = burdenSample.ls, 
		burdenGs = burdenGs.ls,
		geneData = geneData.ls,
		enrRes = enrRes.ls
	)

	return(cnvGSA.out)
}
