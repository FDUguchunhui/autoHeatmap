#'

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
