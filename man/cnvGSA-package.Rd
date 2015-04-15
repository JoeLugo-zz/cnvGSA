\name{cnvGSA-package}
\alias{cnvGSA-package}
\alias{cnvGSA}
\docType{package}
\title{Gene-set Analysis of (Rare) Copy Number Variants}
\description{
cnvGSA is an R package for testing the rare gene-set variant burden in case-control studies of copy number variation (CNV).
}
\details{
In cnvGSA, subjects are treated as statistical sampling units. Subject-level covariates that may act as confounders can be provided by the user (e.g. sex, ethnicity, CNV genotyping platform, CNV genotyping site, array quality metrics, etc.). The gene-set burden is tested using a logistic regression approach. Two logistic regression models are fit: model A includes the subject-level covariates and a variable quantifying global CNV burden for each subject (total CNV length, or total number of CNV-overlapped genes per subject, etc.); model B includes all variables present in model A, plus the number of CNV-overlapped genes that are members of the gene-set being tested. Presence of significantly higher burden in cases compared to controls for the gene-set of interest is then tested by comparing the two models using a deviance chi-square test, as implemented by anova.glm.
}
\author{Daniele Merico \email{daniele.merico@gmail.com}, Robert Ziman \email{rziman@gmail.com}; packaged by Joseph Lugo \email{joseph.r.lugo@gmail.com}}
