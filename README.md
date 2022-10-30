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
data = readRDS("gene_log.rds")
data[1:5,1:5]

#output  :-
HT1080_FUS.DDIT3.EGFP_ruxo_p1 HT1080_FUS.DDIT3.EGFP_ruxo_p2 HT1080_FUS.DDIT3.EGFP_ruxo_p3
TSPAN6                        6.348730                      5.225328                      5.844280
TNMD                          0.000000                      0.000000                      0.000000
DPM1                          7.637741                      5.767928                      7.260560
SCYL3                         2.971615                      2.952795                      2.933835
C1orf112                      5.502201                      4.633630                      4.827602
         HT1080_FUS.DDIT3.EGFP_ruxo_p4 HT1080_FUS.DDIT3.EGFP_p1
TSPAN6                        6.655400                 5.474237
TNMD                          0.000000                 0.000000
DPM1                          7.890874                 7.585082
SCYL3                         2.746343                 3.512331
C1orf112                      5.294000                 4.694626

# Converting a data into matrix form 
data = as.matrix(data)

# Provide an ideal marker's gene data
gene_set = read.csv("markers2Sep.csv")
head(gene_set)

#output :- 
     Art1  Art2   Art3 AstroAQP1 AstroPLCG1 AstroSERPINI2 B.cell Bcell_plasma capEndo    CD4    CD8 CD8_6_1_3       DC
1  HEXIM1 ARL15   DKK2      GFAP       GPC5        ADGRV1   CD37         MZB1   CLDN5 CD40LG   CD8A      GZMA     AREG
2     ID2 MCTP1  LTBP4       TNC      HPSE2         DPP10  CD79A       JCHAIN   IFI27   CD69   CD8B      GZMK     CD74
3    IER3 MECOM    MGP  ADAMTSL3       RORA          TNIK  MS4A1         XBP1   ABCG2   CD3E   XCL1     LIME1   FCER1A
4  LRRC23 PLCB4  RAMP2      AQP1     SLC1A2       PPP2R2B RPL18A        DERL3     VWF  CXCR4 PTGER2     DUSP2   FCGR2B
5  TM4SF1 VEGFC    ELN       DST   ARHGAP24         WDR49   RPS5      HERPUD1           IL7R ZNF683     CRTAM    GRASP
6 C6orf62  ROR1 ABI3BP       TTN     CACNB2         ZNRF3   FCMR        FCRL5           IL32  KLRC4     KLRG1 HLA-DRB5

gene_sets = as.list(as.data.frame(gene_set))
print("genes set ready")

# Assigning a created function to data 
res= ssgsea(data, gene_sets, scale = TRUE, norm = FALSE)
res1 = t(res)
head(res1)

#output :- 

                                   Art1      Art2      Art3 AstroAQP1 AstroPLCG1 AstroSERPINI2    B.cell Bcell_plasma
HT1080_FUS.DDIT3.EGFP_ruxo_p1 0.3054631 0.2586226 0.1705326 0.2065896  0.1910912     0.2350131 0.2291777    0.2621320
HT1080_FUS.DDIT3.EGFP_ruxo_p2 0.3011720 0.2149711 0.1629913 0.2004283  0.1979145     0.2343365 0.2418092    0.2790913
HT1080_FUS.DDIT3.EGFP_ruxo_p3 0.3024787 0.2187553 0.1520157 0.2038492  0.1859622     0.2296506 0.2406960    0.2745030
HT1080_FUS.DDIT3.EGFP_ruxo_p4 0.3080101 0.2479575 0.1591382 0.2095333  0.1917020     0.2387244 0.2261381    0.2645014
HT1080_FUS.DDIT3.EGFP_p1      0.3136549 0.2588294 0.1601708 0.1839985  0.2120165     0.2475338 0.2395430    0.2616747
HT1080_FUS.DDIT3.EGFP_p2      0.3100692 0.2029241 0.1748144 0.1945879  0.2098940     0.2439454 0.2496574    0.2607195

mat=as.matrix(res1)

for(i in 1:nrow(mat)){
  vec=as.numeric(mat[i,])
  mat[i, 1:ncol(mat)]=(vec-mean(vec))/sd(vec)
}

# Plotting a heatmap 
Heatmap(t(mat), col = colorRamp2(c(-2,0,2), c("red", "white", "purple")))
![SSGSEA_heatmap](https://user-images.githubusercontent.com/110582335/197962860-970a5775-8c1c-48db-baa3-069870b284a7.png)

