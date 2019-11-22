library(tidyverse)
load(file = 'temp/normalized_matrix_no_T3.rda')
count_matrix <- normalized_matrix_no_T3

# get_pathway_namelist('data/example.xlsx')
autoHeatmap::subset_count_maxtrix(count_matrix , 'data/example.csv')
# #
autoHeatmap::hmplot(count_matrix , 'data/overlap.csv')

list <- autoHeatmap::get_pathway_namelist(fileName = 'data/overlap.csv')
reduce(list, union)


autoHeatmap::perth_table_trans(data = 'data/trans_demo.xlsx')
p <- autoHeatmap::perth_table_trans(data = 'data/trans_demo.xlsx')


