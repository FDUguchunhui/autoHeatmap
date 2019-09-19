library(tidyverse)
count_matrix <- read.csv(file = 'data/count_matrix.csv')
row.names(count_matrix) <- count_matrix[,'X']
count_matrix <- count_matrix[, -1]

# get_pathway_namelist('data/example.xlsx')
autoHeatmap::subset_count_maxtrix(count_matrix , 'data/example.csv')
# #
autoHeatmap::hmplot(count_matrix , 'data/example.csv')

list <- autoHeatmap::get_pathway_namelist(fileName = 'data/example.csv')
reduce(list, union)
