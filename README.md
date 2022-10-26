# Single-Sample-Gene-Set-Enrichment-Analysis-ssGSEA-
# List of used libraries in ssGSEA :-
library(matrixStats)
library(circlize)
library(ComplexHeatmap)

# Creating a function to perform an analysis .
ssgsea = function(X, gene_sets, alpha = 0.25, scale = T, norm = F, single = T) {
  row_names = rownames(X)
  num_genes = nrow(X)
  gene_sets = lapply(gene_sets, function(genes) {which(row_names %in% genes)})
  
  # Ranks for genes
  R = matrixStats::colRanks(X, preserveShape = T, ties.method = 'average')
  
  # Calculate enrichment score (es) for each sample (column)
  es = apply(R, 2, function(R_col) {
    gene_ranks = order(R_col, decreasing = TRUE)
    
    # Calculate es for each gene set
    es_sample = sapply(gene_sets, function(gene_set_idx) {
      # pos: match (within the gene set)
      # neg: non-match (outside the gene set)
      indicator_pos = gene_ranks %in% gene_set_idx
      indicator_neg = !indicator_pos
      
      rank_alpha  = (R_col[gene_ranks] * indicator_pos) ^ alpha
      
      step_cdf_pos = cumsum(rank_alpha)    / sum(rank_alpha)
      step_cdf_neg = cumsum(indicator_neg) / sum(indicator_neg)
      
      step_cdf_diff = step_cdf_pos - step_cdf_neg
      
      # Normalize by gene number
      if (scale) step_cdf_diff = step_cdf_diff / num_genes
      
      # Use ssGSEA or not
      if (single) {
        sum(step_cdf_diff)
      } else {
        step_cdf_diff[which.max(abs(step_cdf_diff))]
      }
    })
    unlist(es_sample)
  })
  
  if (length(gene_sets) == 1) es = matrix(es, nrow = 1)
  
  # Normalize by absolute diff between max and min
  if (norm) es = es / diff(range(es))
  
  # Prepare output
  rownames(es) = names(gene_sets)
  colnames(es) = colnames(X)
  return(es)
}

# Reading a log data created in normalization method
data = readRDS("gene_log.rds")
data[1:5,1:5]

# Converting a data into matrix form 
data = as.matrix(data)

# Provide an ideal marker's gene data
gene_set = read.csv("markers2Sep.csv")
head(gene_set)

gene_sets = as.list(as.data.frame(gene_set))
print("genes set ready")

# Assigning a created function to data 
res= ssgsea(data, gene_sets, scale = TRUE, norm = FALSE)
res1 = t(res)
head(res1)

mat=as.matrix(res1)

for(i in 1:nrow(mat)){
  vec=as.numeric(mat[i,])
  mat[i, 1:ncol(mat)]=(vec-mean(vec))/sd(vec)
}

# Plotting a heatmap 
Heatmap(t(mat), col = colorRamp2(c(-2,0,2), c("orangered", "white", "purple")))
![SSGSEA_heatmap](https://user-images.githubusercontent.com/110582335/197962860-970a5775-8c1c-48db-baa3-069870b284a7.png)

