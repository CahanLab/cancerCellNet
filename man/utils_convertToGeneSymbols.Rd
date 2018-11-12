% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_convertToGeneSymbols.R
\name{utils_convertToGeneSymbols}
\alias{utils_convertToGeneSymbols}
\title{Gene Symbol Conversion}
\usage{
utils_convertToGeneSymbols(expdf, typeENST = FALSE, typeENSG = FALSE,
  typeMusGene = FALSE)
}
\arguments{
\item{expdf}{the dataframe of experssion data with sample names as column names and gene names as row names}

\item{typeENSG}{TRUE if the gene identifier is in Ensembl Transcript ID}

\item{typeMusGene}{TRUE if the gene identifier is in mouse gene symbols}

\item{typeENSG}{TRUE if the gene identifier is in Ensembl Gene ID}
}
\value{
a list containing the converted experssion table as the first element and the expression table with genes that cannot be converted with our conversion dataset
}
\description{
The function can convert Ensembl Gene ID, Ensembl Transcript ID or mouse gene symbol into human gene symbol
}