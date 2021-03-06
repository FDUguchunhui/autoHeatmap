% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DESeqtools.R
\name{res_subgroup}
\alias{res_subgroup}
\title{subset the result of a differential expression analysis}
\usage{
res_subgroup(res, alpha = 0.1, reg_LFC = 1, reg_dir = "all")
}
\arguments{
\item{res}{an obj from DEseq2.result function}

\item{alpha}{the significant level for adjusted P-value}

\item{reg_LFC}{reg gives the regulation level change in log2 fold change in absolute value}

\item{reg_dir}{gives which regulation direction you want to subset you gene
three options: all -- up and down
                up  -- only up regulated
              down -- only down regulated}
}
\description{
subset the result of a differential expression analysis
}
