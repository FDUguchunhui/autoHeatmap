.libPaths("D:/R-3.5.1/library")
#' @import pheatmap gplots tibble tidyr readr dplyr viridis
# options
options(tibble.print_max = 10, tibble.width = Inf)

#' plot heatmap by integrate data from DESeq and enrichment analysis
#'
#' more information
#' @param fileName is the post-processed csv file from enrichment analysis
#' @return a list of tibbles of external gene names from different enrichment pathway
#' @export
get_pathway_namelist <- function(fileName){

  #read the enrichment dataset as tibble
  enrichment_csv <- file(fileName)
  string_list <- scan(enrichment_csv, what = '', sep = '\n')
  start <- which(startsWith(string_list, 'Gene/Gene Set Overlap'))
  enrichment_tibble <- readr::read_csv(file = fileName, skip = start + 1,
                                       col_types = cols(
                                         `Entrez Gene Id` = col_character())
  )
  close(enrichment_csv)
  # the enrichment.table have data structure like this
  # Entrez.Gene.Id Gene.Symbol Gene.Description                               HALLMARK_E2F_TARGETS
  # <dbl> <fct>       <fct>                                          <fct>
  #   1           7153 TOP2A       topoisomerase (DNA) II alpha 170kDa            HALLMARK_E2F_TARGETS
  # 2            991 CDC20       cell division cycle 20 homolog (S. cerevisiae) HALLMARK_E2F_TARGETS
  # 3          10733 PLK4        polo-like kinase 4                             HALLMARK_E2F_TARGETS
  # 4           6790 AURKA       aurora kinase A                                HALLMARK_E2F_TARGETS
  # 5           5347 PLK1        polo-like kinase 1                             HALLMARK_E2F_TARGETS
  # 6            983 CDK1        cyclin-dependent kinase 1                      HALLMARK_E2F_TARGETS
  # 7            332 BIRC5       baculoviral IAP repeat containing 5            HALLMARK_E2F_TARGETS
  # 8          11004 KIF2C       kinesin family member 2C                       HALLMARK_E2F_TARGETS
  # 9           9133 CCNB2       cyclin B2                                      HALLMARK_E2F_TARGETS
  # 10          24137 KIF4A       kinesin family member 4A                       HALLMARK_E2F_TARGETS

  # create a number of list corresponding to the number pathway in the enrichment analysis
  # First, get a vector of those pathway names
  pathway_names_vector <- colnames(enrichment_tibble)[-c(1:3)]
  # a list of tibble contain two column, column 1 is gene name
  #
  pathway_genes_list <- list()
  for(i in 1:length(pathway_names_vector)){
    pathway_genes_list[[i]] <- dplyr::select(enrichment_tibble, `Gene Symbol`,
                                             pathway_names_vector[i]) %>%
      dplyr::filter(.[[2]] != 'NA') %>%
      dplyr::select(`Gene Symbol`) %>%
      dplyr::pull(var = 1)
  }


  # assign the pathway name to each subset on the list
  names(pathway_genes_list) <- pathway_names_vector

  # add one additional gene names vector, which are those genes are not in any
  #of those pathway
  genes_union <- purrr::reduce(pathway_genes_list, union)
  gene_left <- filter(enrichment_tibble, !(`Gene Symbol` %in% genes_union)) %>%
    select(`Gene Symbol`)
  # add it into the list
  pathway_genes_list <- c(pathway_genes_list, left = gene_left)

  return(pathway_genes_list)
}


#' get tibbles of different subgrouped genes by pathway
#'
#' more information
#' @param count_matrix  the RNA-seq expression level count matrix, typically a
#'   normalized count matrix is used
#' @param enrichment  the modified enrichment export csv file from xlsx file
#' @param fun the fun used to get lists of genes names from enrichment analysis
#'   export
#' @return a list of tibbles of count matrix for different pathways
#' @export
#'
#'
subset_count_maxtrix <- function(count_matrix, enrichment){
  # convert count matrix as tibble
  count_matrix_tibble <- as_tibble((data.frame(gene_name = rownames(count_matrix),
                                               count_matrix)))
  pathway_genes_names_list <- get_pathway_namelist(enrichment)
  pathway_names <- names(pathway_genes_names_list)
  # get a list of lists of gene names only
  # use the lists of gene to subset the count matrix
  # initialize the subgroup
  count_matrix_subgrouped_list <- list()
  # get a list of those subgroup genes sets
  for(i in 1:length(pathway_genes_names_list)){
    count_matrix_subgrouped_list[[i]] <-
      filter(count_matrix_tibble, gene_name %in% pathway_genes_names_list[[i]])
  }
  result_list <- list('subgroups' = count_matrix_subgrouped_list, 'names' = pathway_names)
  return(result_list)
}



#' plot heatmap by integrate data from DESeq and enrichment analysis
#'
#' more information
#' @param count_matrix the RNA-seq expression level count matrix
#' @param enrichment_xlsx the post-processed Excel file from enrichment analysis
#' @param file the destination and name you want to save your heatmap, by default is the
#' @param ... other parameters used in pheatmap
#' your current work directory, the default name is heatmap.pdf
#' @return output a pdf file of heatmaps of all the pathways in the enrichment_xlsx
#' @export
#'
#'
hmplot <- function(count_matrix, enrichment, file = 'heatmap.pdf',
                   cluster_rows = TRUE,
                   cluster_cols = TRUE,
                   show_colnames = TRUE,
                   show_rownames = TRUE,
                   annotation_col = NA,
                   border_color = 'black',
                   color = inferno(10)){
  # get the a list of tibbles of subsets of genes
  genes_tibbles_list <- subset_count_maxtrix(count_matrix, enrichment)$subgroups
  # get the names of objects in list, will be used as title of heatmaps
  pathway_names <- subset_count_maxtrix(count_matrix, enrichment)$names


  # transform tibbles into count_matrices with first column as rownames
  genes_pre_matrices_list <- lapply(genes_tibbles_list, select, -1)
  genes_count_matrices_list <- lapply(genes_pre_matrices_list, data.matrix)
  for(i in 1:length(genes_tibbles_list)){
    rownames(genes_count_matrices_list[[i]]) <- pull(select(genes_tibbles_list[[i]], 1))
  }
  pdf(file)
  # plot the heatmap for subsets
  for(i in 1:(length(genes_tibbles_list))){
    #heatmap(genes_count_matrices.list[[i]], main = pathway_names[i], Rowv = F)
    #heatmap.2(genes_count_matrices.list[[i]], main = pathway_names[i], Rowv = F)
    pheatmap(genes_count_matrices_list[[i]], main = pathway_names[i],
             cluster_rows = cluster_rows,
             cluster_cols = cluster_cols,
             show_colnames = show_colnames,
             show_rownames = show_rownames,
             annotation_col = annotation_col,
             border_color = border_color,
             color = color
             )
  }

  graphics.off()

}
