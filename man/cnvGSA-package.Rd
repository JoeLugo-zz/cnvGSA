\name{cnvGSA-package}
\alias{cnvGSA-package}
\alias{cnvGSA}
\docType{package}
\title{Gene-Set Analysis of (Rare) Copy Number Variants}
\description{
cnvGSA is an R package meant to facilitate gene-set analysis of (rare) copy number variants (CNVs).
}
\details{
Known gene-sets are tested for prevalence of rare variants in case vs. control subjects. Whenever a subject has at least one gene in a gene-set affected by a rare variant, a perturbation score of 1 is assigned to the (subject, gene-set) pair; for each gene-set, subject counts are tested vs. control counts using the Fisher Exact Test (FET). Significant gene-sets will have a significantly high count in cases compared to controls. Statistical reports on burden are also generated.
}
\author{Daniele Merico \email{daniele.merico@gmail.com}; packaged by Robert Ziman \email{rziman@gmail.com}}
