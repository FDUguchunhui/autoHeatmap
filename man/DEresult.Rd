% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DESeqtools.R
\name{DEresult}
\alias{DEresult}
\title{extract log fold change information from a list of DESeqResults, return a tibble}
\usage{
DEresult(res_list, gene_list = NULL, col_name = NULL, ...)
}
\arguments{
\item{res_list}{a list of DESeqResults, must have the same length}

\item{gene_list}{a list of gene that you want to subset.
If you don't provide this, by default this function will return all the gene from DESeqResults}

\item{col_name}{a character vector that gives column name of the output}

\item{...}{other parameters can be used in res_subgroup function

reg_dir
gives which regulation direction you want to subset you gene
three options: all – up and down up – only up regulated down – only down regulated
default is all

alpha  the significant level for adjusted P-value, default is 1, i.e. no alpha requirement

reg_LFC reg gives the regulation level change in log2 fold change in absolute value
default is 0, i.e. no LFC requirement}
}
\description{
extract log fold change information from a list of DESeqResults, return a tibble
}
