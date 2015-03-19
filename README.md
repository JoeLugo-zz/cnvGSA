# cnvGSA
This package is intended to facilitate gene-set association with rare CNVs in case-control studies.
cnvGSA is an R package for testing the gene-set rare variant burden in case-control studies of copy number variation (CNV).

Gene-sets need to be pre-compiled based on user-curated data or publicly available gene annotations (Gene Ontology, pathways, etc...)

"Competitive" gene-set over-representation tests are commonly used to analyze differentially gene expressed genes (e.g. Fisher's Exact Test, GSEA), but they are not suitable for rare CNV; the most appropriate choice for rare CNV is a "self-contained" burden test with global burden correction such as implemented by \texttt{cnvGSA} or other tools.

Global burden correction is important because, for many disorders implicating rare CNV (autism, schizophrenia, etc...), case subjects are enriched in large, recurrent CNV that are not observed (or observed at very low frequency) in controls, and in which only a minority of genes contribute to disease risk. In the absence of global burden correction, many gene-sets would present a biologically unspecific burden, uniquely driven by those larger and recurrent CNV. Global burden correction thus helps identifying specific pathways and functional categories implicated in disease risk by rare CNV.\\

In cnvGSA, subjects are treated as statistical sampling units. Subject-level covariates that may act as confounders can be provided by the user (sex, ethnicity, CNV genotyping platform, CNV genotyping site, array quality metrics, etc...). The gene-set burden is tested using a logistic regression approach. Two logistic regression models are fit: model A includes the subject-level covariates (described above) and a variable quantifying global CNV burden for each subject (total CNV length, or total number of CNV-overlapped genes per subject, etc...); model B includes all variables present in model A, plus the number of CNV-overlapped genes that are members of the gene-set being tested. Presence of significantly higher burden in cases compared to controls is tested by comparing the two models using a deviance test, as implemented by anova.glm
