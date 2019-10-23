#' @title subset the result of a differential expression analysis
#'
#' @description subset the result of a differential expression analysis
#' @param res an obj from DEseq2.result function
#' @param alpha the significant level for adjusted P-value
#' @param reg_LFC reg gives the regulation level change in log2 fold change in absolute value
#' @param reg_dir gives which regulation direction you want to subset you gene
#' three options: all -- up and down
#'                 up  -- only up regulated
#'               down -- only down regulated
#' @export
res_subgroup <- function(res, alpha=0.1, reg_LFC=1, reg_dir='all'){

  res_sig_pos <- (res$padj < alpha)
  res_sig_pos[is.na(res_sig_pos)] <- F
  if(reg_dir == 'all'){
    res_LFC_pos <- (res$log2FoldChange > reg_LFC) | (res$log2FoldChange < -reg_LFC)
  }
  else if(reg_dir == 'up'){
    res_LFC_pos <- (res$log2FoldChange > reg_LFC)
  }
  else if(reg_dir == 'down'){
    res_LFC_pos <- (res$log2FoldChange < -reg_LFC)
  }

  return(res[res_LFC_pos & res_sig_pos,])

}

#'@title extract log fold change information from a list of DESeqResults, return a tibble
#'
#'@description extract log fold change information from a list of DESeqResults, return a tibble
#'@param res_list a list of DESeqResults, must have the same length
#'@param gene_list a list of gene that you want to subset.
#'If you don't provide this, by default this function will return all the gene from DESeqResults
#'@param col_name a character vector that gives column name of the output
#'@param ... other parameters can be used in res_subgroup function
#'
#'reg_dir
#'gives which regulation direction you want to subset you gene
#'three options: all – up and down up – only up regulated down – only down regulated
#'default is all
#'
#'alpha  the significant level for adjusted P-value, default is 1, i.e. no alpha requirement
#'
#'reg_LFC reg gives the regulation level change in log2 fold change in absolute value
#'default is 0, i.e. no LFC requirement
#'@export
DEresult <- function(res_list, gene_list = NULL, col_name = NULL, ...){
  # create a list that contain all result without of lfc and alpha requirement
  res_list_all <-
    lapply(
      res_list,
      res_subgroup,
      reg_dir = 'all',
      alpha = 1,
      reg_LFC = 0
    )

  # extract only the 'gene_id' and 'log2FoldChange column'
  res_list_all <-
    lapply(
      lapply(res_list_all, as_tibble, rownames = 'gene_id'),
      select,
      `gene_id`,
      `log2FoldChange`
    )

  # merge all list into one dataset that contain lfc for all the comparisons for all genes
  # this will the reference dataset to be subsetted
  # in this step the column name has duplicate, but it is okay. Ignore the warning message
  res_all <- Reduce(function(x, y){
    merge(
      x,
      y,
      by = 'gene_id',
      all = TRUE,
      no.dups = T
    )},
    x = res_list_all)

  colnames(res_all) <- c('gene_name', names(res_list_all))

  if(!is.null(col_name)){
    # rename the columns
    colnames(res_all) <- col_name
  }


  # judge whether to subset
  if (!is.null(gene_list))
    res_all <- res_all[res_all$gene_id %in% gene_list,]

  return(res_all)
}


