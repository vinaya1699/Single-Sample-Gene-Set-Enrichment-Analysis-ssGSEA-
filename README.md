# Single-Sample-Gene-Set-Enrichment-Analysis-ssGSEA-
# List of used libraries in ssGSEA :-
library(matrixStats)

library(circlize)

library(ComplexHeatmap)

# Creating a gene set for analysis

data_files <- list.files("C:/Users/91973/Desktop/New folder", full.names= T)  

data_files 

gene_set= matrix(NA,ncol= length(data_files),nrow= 9)

for(i in 1:length(data_files)) {                              
  
  df= read.csv(data_files[i])
  
  gene_set[,i]= df[1:9, 1]
  
}
colnames(gene_set) <- c('M3652','M38668','p53')   #M3652 and M38668 : Genesets used for analysis

View(gene_set)

gene_sets= as.list(as.data.frame(gene_set))

# Creating a function to perform an analysis .
ssgsea = function(X, gene_sets, alpha = 0.25, scale = T, norm = F, single = T) {
    row_names = rownames(X)
    num_genes = nrow(X)
    gene_sets = lapply(gene_sets, function(genes) {which(row_names %in% genes)})
        #Ranks for genes
    R = matrixStats::colRanks(X, preserveShape = T, ties.method = 'average')

    #Calculate enrichment score (es) for each sample (column)
    es = apply(R, 2, function(R_col) {
        gene_ranks = order(R_col, decreasing = TRUE)

        #Calc es for each gene set
        es_sample = sapply(gene_sets, function(gene_set_idx) {
            #pos: match (within the gene set)
            #neg: non-match (outside the gene set)
            indicator_pos = gene_ranks %in% gene_set_idx
            indicator_neg = !indicator_pos

            rank_alpha  = (R_col[gene_ranks] * indicator_pos) ^ alpha

            step_cdf_pos = cumsum(rank_alpha)    / sum(rank_alpha)
            step_cdf_neg = cumsum(indicator_neg) / sum(indicator_neg)

            step_cdf_diff = step_cdf_pos - step_cdf_neg

            #Normalize by gene number
            if (scale) step_cdf_diff = step_cdf_diff / num_genes

            #Use ssGSEA or not
            if (single) {
                sum(step_cdf_diff)
            } else {
                step_cdf_diff[which.max(abs(step_cdf_diff))]
            }
        })
        unlist(es_sample)
    })

    if (length(gene_sets) == 1) es = matrix(es, nrow = 1)

    #Normalize by absolute diff between max and min
    if (norm) es = es / diff(range(es))

    #Prepare output
    rownames(es) = names(gene_sets)
    colnames(es) = colnames(X)
    return(es)
}

# Reading a log data created in normalization method
data1= readRDS("gene_log.rds")

data2=as.matrix(data1)

data2[1:5,1:5]
![image](https://user-images.githubusercontent.com/110582335/200167591-5108d9d5-7582-4772-b240-71058a82e8f8.png)

# Provide an ideal marker's gene data
gene_set = read.csv("markers2Sep.csv")

head(gene_set)

print("genes set ready")

# Assigning a created function to data 
res= ssgsea(data, gene_sets, scale = TRUE, norm = FALSE)

res1 = t(res)

head(res1)
![image](https://user-images.githubusercontent.com/110582335/200167697-3a415970-9729-45c1-b316-4d4c0d0c2ad5.png)

mat=as.matrix(res1)

for(i in 1:nrow(mat)){

vec=as.numeric(mat[i,])

mat[i, 1:ncol(mat)]=(vec-mean(vec))/sd(vec)

}

# Plotting a heatmap 
heatmap(t(mat))
![heatmap](https://user-images.githubusercontent.com/110582335/200167764-accd4d71-8b07-4737-9d03-7343f3e0c2e7.png)

# Interpretation :
Above Heatmap contain cancer data samples and genes sets where,colour representation is based on specific genes present in respective samples .

