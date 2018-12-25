Better pipeline with Seurat
================
German Novakovskiy
December 24, 2018

    ## Loading required package: ggplot2

    ## Loading required package: cowplot

    ## 
    ## Attaching package: 'cowplot'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     ggsave

    ## Loading required package: Matrix

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

Followed by [this tutorial](https://satijalab.org/seurat/pbmc3k_tutorial.html).

Setup the Seurat Object
-----------------------

``` r
#Loading the data
data.dir <-  "../../Sox17Ng_36h_German/filtered_feature_bc_matrix/"
esc.data <- Read10X(data.dir = data.dir)

#for visualization purporse
x <- row.names(esc.data)
x <- sapply(x, function(x){unlist(strsplit(x, split = "_____"))[2]})
names(x) <- NULL
x[1] <- "Neogreen"

rownames(esc.data) <- x
rm(x)

# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 7 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes
esc <- CreateSeuratObject(raw.data = esc.data, min.cells = 7, min.genes = 200, 
    project = "10X_ESC")
```

QC and selecting cells for further analysis
-------------------------------------------

``` r
# The number of genes and UMIs (nGene and nUMI) are automatically calculated
# for every object by Seurat.  For non-UMI data, nUMI represents the sum of
# the non-normalized values within a cell We calculate the percentage of
# mitochondrial genes here and store it in percent.mito using AddMetaData.
# We use object@raw.data since this represents non-transformed and
# non-log-normalized counts The % of UMI mapping to MT-genes is a common
# scRNA-seq QC metric.
mito.genes <- grep(pattern = "^MT-", x = rownames(x = esc@data), value = TRUE)
percent.mito <- Matrix::colSums(esc@raw.data[mito.genes, ])/Matrix::colSums(esc@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
esc <- AddMetaData(object = esc, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = esc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3, point.size.use = 0.005)
```

![](Exploratory_analysis_seurat_better_version_files/figure-markdown_github/unnamed-chunk-3-1.png)

``` r
# GenePlot is typically used to visualize gene-gene relationships, but can
# be used for anything calculated by the object, i.e. columns in
# object@meta.data, PC scores etc.  Since there is a rare subset of cells
# with an outlier level of high mitochondrial percentage and also low UMI
# content, we filter these as well
par(mfrow = c(1, 2))
GenePlot(object = esc, gene1 = "nUMI", gene2 = "percent.mito", pch.use = 16, cex.use = 0.5)
GenePlot(object = esc, gene1 = "nUMI", gene2 = "nGene", pch.use = 16, cex.use = 0.5)
```

![](Exploratory_analysis_seurat_better_version_files/figure-markdown_github/unnamed-chunk-4-1.png)

``` r
# We filter out cells that have unique gene counts ess than
# 200 Note that low.thresholds and high.thresholds are used to define a
# 'gate'.  -Inf and Inf should be used if you don't want a lower or upper
# threshold.
esc <- FilterCells(object = esc, subset.names = c("nGene", "percent.mito"), 
    low.thresholds = c(3000, -Inf), high.thresholds = c(Inf, 0.05))
```

Normalizing the data
--------------------

After removing unwanted cells from the dataset, the next step is to normalize the data. By default, we employ a global-scaling normalization method “LogNormalize” that normalizes the gene expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result. (Think about Scran normalization from hemberg approach).

``` r
esc <- NormalizeData(object = esc, normalization.method = "LogNormalize", 
    scale.factor = 10000)
```

Detection of variable genes across the single cells
---------------------------------------------------

Seurat calculates highly variable genes and focuses on these for downstream analysis. FindVariableGenes calculates the average expression and dispersion for each gene, places these genes into bins, and then calculates a z-score for dispersion within each bin. This helps control for the relationship between variability and average expression. This function is unchanged from (Macosko et al.), but new methods for variable gene expression identification are coming soon. We suggest that users set these parameters to mark visual outliers on the dispersion plot, but the exact parameter settings may vary based on the data type, heterogeneity in the sample, and normalization strategy. The parameters here represent typical parameter settings for UMI data that is normalized to a total of 1e4 molecules:

``` r
esc <- FindVariableGenes(object = esc, mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
```

![](Exploratory_analysis_seurat_better_version_files/figure-markdown_github/unnamed-chunk-7-1.png)

``` r
length(x = esc@var.genes)
```

    ## [1] 1548

Scaling the data and removing unwanted sources of variation
-----------------------------------------------------------

Our single cell dataset likely contains ‘uninteresting’ sources of variation. This could include not only technical noise, but batch effects, or even biological sources of variation (cell cycle stage). Regressing these signals [out](https://stats.stackexchange.com/questions/3944/explain-model-adjustment-in-plain-english) of the analysis can improve downstream dimensionality reduction and clustering.

We check the cell cycle for confounding effect. And then scale data for UMI counts, MT genes and difference in proliferating cycles (S and G2/M). We should regress out the difference between the G2M and S phase scores. This means that signals separating non-cycling cells and cycling cells will be maintained, but differences in cell cycle phase amongst proliferating cells (which are often uninteresting), will be regressed out of the data (because we have stem cells):

``` r
#cell cycle markers
cc.genes <- readLines(con = "cell_cycle_vignette_files/regev_lab_cell_cycle_genes.txt")
#cc.genes <- sapply(cc.genes, function(x){ paste("GRCh38_____", x, sep="")})
#names(cc.genes) <- NULL

# We can segregate this list into markers of G2/M phase and markers of S
# phase
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

esc <- ScaleData(object = esc, display.progress = FALSE)
#only NASP is in PC4, there are no other genes, related to cell-cycle; however a lot of are MT genes, which is expected
esc <- RunPCA(object = esc, pc.genes = esc@var.genes, pcs.print = 1:4, 
    genes.print = 10)
```

    ## [1] "PC1"
    ##  [1] "CER1"    "DKK1"    "FGF17"   "MIXL1"   "CYP26A1" "GATA6"   "LHX1"   
    ##  [8] "HAS2"    "DKK4"    "ARL4D"  
    ## [1] ""
    ##  [1] "CD24"       "AL353747.4" "PODXL"      "SEMA6A"     "TERF1"     
    ##  [6] "SOX2"       "DPPA4"      "LINC00678"  "PLA2G16"    "USP44"     
    ## [1] ""
    ## [1] ""
    ## [1] "PC2"
    ##  [1] "HMGN2" "ARL4D" "HAS2"  "CER1"  "LHX1"  "FGF17" "DKK1"  "GATA6"
    ##  [9] "APLNR" "DKK4" 
    ## [1] ""
    ##  [1] "HSPD1"    "TUBA1B"   "HSP90AA1" "ATP5F1B"  "LDHB"     "CCT5"    
    ##  [7] "CCT2"     "PDIA3"    "SSB"      "ENO1"    
    ## [1] ""
    ## [1] ""
    ## [1] "PC3"
    ##  [1] "GAL"       "LINC01356" "CNTNAP2"   "SERPINE2"  "UBE2C"    
    ##  [6] "BEX4"      "CDK1"      "AURKA"     "TOP2A"     "LINC00458"
    ## [1] ""
    ##  [1] "CD24"       "GJA1"       "AL118516.1" "PCNA"       "GINS2"     
    ##  [6] "SRSF7"      "CA2"        "FABP5"      "CTNNB1"     "LHX1"      
    ## [1] ""
    ## [1] ""
    ## [1] "PC4"
    ##  [1] "DUT"        "HIST1H1D"   "AC106864.1" "PCNA"       "GINS2"     
    ##  [6] "MCM5"       "SRSF7"      "UNG"        "FABP5"      "HIST1H1C"  
    ## [1] ""
    ##  [1] "UBE2C"   "HNRNPH1" "SOX11"   "CKS2"    "CENPF"   "AURKA"   "TOP2A"  
    ##  [8] "PLK1"    "CTNNB1"  "B3GNT7" 
    ## [1] ""
    ## [1] ""

``` r
esc <- CellCycleScoring(object = esc, s.genes = s.genes, g2m.genes = g2m.genes, 
    set.ident = TRUE)

# view cell cycle scores and phase assignments
head(x = esc@meta.data)
```

    ##                  nGene  nUMI orig.ident percent.mito      S.Score
    ## AAACCTGAGATGGGTC  4710 29613    10X_ESC   0.03438724 -0.115589141
    ## AAACCTGGTTGGAGGT  3890 21177    10X_ESC   0.04066116  0.141085671
    ## AAACGGGAGACGCTTT  4774 34914    10X_ESC   0.03242624 -0.021168221
    ## AAACGGGAGTGAAGAG  3442 17309    10X_ESC   0.04414399 -0.017462891
    ## AAACGGGCACATCCGG  4039 21412    10X_ESC   0.02699799 -0.044008609
    ## AAACGGGCACGAAGCA  4790 34842    10X_ESC   0.03019865  0.008521502
    ##                     G2M.Score Phase old.ident
    ## AAACCTGAGATGGGTC  0.057378497   G2M   10X_ESC
    ## AAACCTGGTTGGAGGT -0.089979914     S   10X_ESC
    ## AAACGGGAGACGCTTT -0.065673809    G1   10X_ESC
    ## AAACGGGAGTGAAGAG  0.004807650   G2M   10X_ESC
    ## AAACGGGCACATCCGG  0.002509739   G2M   10X_ESC
    ## AAACGGGCACGAAGCA -0.119344566     S   10X_ESC

``` r
esc <- RunPCA(object = esc, pc.genes = c(s.genes, g2m.genes), do.print = FALSE)
PCAPlot(object = esc)
```

![](Exploratory_analysis_seurat_better_version_files/figure-markdown_github/unnamed-chunk-11-1.png) For each gene, Seurat models the relationship between gene expression and the S and G2M cell cycle scores. The scaled residuals of this model represent a ‘corrected’ expression matrix, that can be used downstream for dimensional reduction.

``` r
esc@meta.data$CC.Difference <- esc@meta.data$S.Score - esc@meta.data$G2M.Score
esc <- ScaleData(object = esc, vars.to.regress = c("nUMI", "percent.mito", "CC.Difference"), display.progress = FALSE)
#esc <- ScaleData(object = esc, vars.to.regress = c("nUMI", "percent.mito"), display.progress = FALSE)

esc <- RunPCA(object = esc, pc.genes = c(s.genes, g2m.genes), do.print = FALSE)
PCAPlot(object = esc)
```

![](Exploratory_analysis_seurat_better_version_files/figure-markdown_github/unnamed-chunk-12-1.png)

Linear dimensional reduction
----------------------------

Running dimensionality reduction on highly variable genes can improve performance. However, with UMI data - particularly after regressing out technical variables, PCA returns similar (albeit slower) results when run on much larger subsets of genes, including the whole transcriptome.

``` r
esc <- RunPCA(object = esc, pc.genes = esc@var.genes, pcs.compute = 100, pcs.print = 1:4, do.print=FALSE)
```

``` r
PCAPlot(object = esc, dim.1 = 1, dim.2 = 2)
```

![](Exploratory_analysis_seurat_better_version_files/figure-markdown_github/unnamed-chunk-14-1.png)

``` r
esc <- RunPCA(object = esc, pc.genes = rownames(esc@data), pcs.compute = 100, pcs.print = 1:4, do.print=FALSE)
```

``` r
PCAPlot(object = esc, dim.1 = 1, dim.2 = 2)
```

![](Exploratory_analysis_seurat_better_version_files/figure-markdown_github/unnamed-chunk-16-1.png)

``` r
#x <- esc@meta.data
#x <- x[ , !(names(x) %in% c("Phase"))]
#esc@meta.data <- x
PCAPlot(object = esc, dim.1 = 1, dim.2 = 2, group.by='old.ident')
```

![](Exploratory_analysis_seurat_better_version_files/figure-markdown_github/unnamed-chunk-17-1.png)

PCHeatmap allows for easy exploration of the primary sources of heterogeneity in a dataset, and can be useful when trying to decide which PCs to include for further downstream analyses. Both cells and genes are ordered according to their PCA scores.

``` r
PCHeatmap(object = esc, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
```

![](Exploratory_analysis_seurat_better_version_files/figure-markdown_github/unnamed-chunk-18-1.png)

``` r
PCHeatmap(object = esc, pc.use = 1:15, cells.use = 500, do.balanced = TRUE, 
    label.columns = FALSE, use.full = FALSE)
```

![](Exploratory_analysis_seurat_better_version_files/figure-markdown_github/unnamed-chunk-19-1.png)

``` r
PCElbowPlot(object = esc)
```

![](Exploratory_analysis_seurat_better_version_files/figure-markdown_github/unnamed-chunk-20-1.png)

Clustering and t-SNE
--------------------

``` r
esc <- FindClusters(object = esc, reduction.type = "pca", dims.use = 1:40, 
    resolution = 0.6, print.output = 0, save.SNN = TRUE)

esc <- RunTSNE(object = esc, dims.use = 1:40, do.fast = TRUE)

TSNEPlot(object = esc)
```

![](Exploratory_analysis_seurat_better_version_files/figure-markdown_github/unnamed-chunk-21-1.png)

``` r
table(esc@ident) 
```

    ## 
    ##    0    1    2    3    4 
    ## 1113 1084  200   94   20

``` r
PCAPlot(object = esc, dim.1 = 1, dim.2 = 2)
```

![](Exploratory_analysis_seurat_better_version_files/figure-markdown_github/unnamed-chunk-23-1.png)

200 cells are in cluster 4, which is believed to be primitive streak

Finding of bio markerks for each cluster:

``` r
# find markers for every cluster compared to all remaining cells, report
# only the positive ones
esc.markers <- FindAllMarkers(object = esc, only.pos = TRUE, min.pct = 0.25, 
    thresh.use = 0.25)

esc.markers %>% group_by(cluster) %>% top_n(10, avg_logFC) %>% filter(cluster == 2)
```

    ## # A tibble: 10 x 7
    ## # Groups:   cluster [1]
    ##        p_val avg_logFC pct.1 pct.2 p_val_adj cluster gene   
    ##        <dbl>     <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>  
    ##  1 2.72e-282     1.16  0.845 0.046 4.93e-278 2       LHX1   
    ##  2 7.36e-218     1.18  0.92  0.109 1.33e-213 2       FGF17  
    ##  3 4.88e-157     2.25  0.94  0.232 8.84e-153 2       DKK1   
    ##  4 1.09e-154     1.05  0.92  0.196 1.97e-150 2       HAS2   
    ##  5 3.92e-148     1.65  0.765 0.118 7.10e-144 2       DKK4   
    ##  6 1.72e-121     2.27  0.975 0.446 3.11e-117 2       CER1   
    ##  7 9.10e-117     1.58  0.895 0.283 1.65e-112 2       CYP26A1
    ##  8 2.72e-109     1.29  0.995 0.726 4.92e-105 2       MIXL1  
    ##  9 3.86e- 91     1.09  0.845 0.328 6.99e- 87 2       ARL4D  
    ## 10 1.43e- 44     0.852 0.995 0.966 2.59e- 40 2       NTS

Cluster 4 is full with mesoderm and primitive streak markers.

If we check gene ontology enrichment of cluster 3 markers we see mostly negative regulation of gene expression, cell cycle and organelle organization.

Checking other key markers of primitive streak:

``` r
VlnPlot(object = esc, features.plot = c("EOMES", "CER1"))
```

![](Exploratory_analysis_seurat_better_version_files/figure-markdown_github/unnamed-chunk-25-1.png)

``` r
VlnPlot(object = esc, features.plot = c("DKK4", "PDGFRA"))
```

![](Exploratory_analysis_seurat_better_version_files/figure-markdown_github/unnamed-chunk-26-1.png)

``` r
VlnPlot(object = esc, features.plot = c("GSC", "GATA6"))
```

![](Exploratory_analysis_seurat_better_version_files/figure-markdown_github/unnamed-chunk-27-1.png)

``` r
VlnPlot(object = esc, features.plot = c("LHX1", "MIXL1"))
```

![](Exploratory_analysis_seurat_better_version_files/figure-markdown_github/unnamed-chunk-28-1.png)

Some of the markers of hESC:

``` r
VlnPlot(object = esc, features.plot = c("SOX2", "NANOG"))
```

![](Exploratory_analysis_seurat_better_version_files/figure-markdown_github/unnamed-chunk-29-1.png)

``` r
VlnPlot(object = esc, features.plot = c("POU5F1", "DPPA4"))
```

![](Exploratory_analysis_seurat_better_version_files/figure-markdown_github/unnamed-chunk-30-1.png)

``` r
VlnPlot(object = esc, features.plot = c("CDH1", "FGF2"))
```

![](Exploratory_analysis_seurat_better_version_files/figure-markdown_github/unnamed-chunk-31-1.png)

Markers of DE - we don't have DE in this culture:

``` r
#SOX17 was filtered out as undetected gene, it was not expressed
VlnPlot(object = esc, features.plot = c("FOXA2", "CXCR4"))
```

![](Exploratory_analysis_seurat_better_version_files/figure-markdown_github/unnamed-chunk-32-1.png)

``` r
#markers of primitive streak
FeaturePlot(object = esc, features.plot = c("CER1", "EOMES", "GATA6", "LHX1"), cols.use = c("grey", "red"), 
    reduction.use = "tsne")
```

![](Exploratory_analysis_seurat_better_version_files/figure-markdown_github/unnamed-chunk-33-1.png)

``` r
#markers of primitive streak
FeaturePlot(object = esc, features.plot = c("MIXL1", "DKK4", "GSC", "DKK1"), cols.use = c("grey", "red"), 
    reduction.use = "tsne")
```

![](Exploratory_analysis_seurat_better_version_files/figure-markdown_github/unnamed-chunk-34-1.png) For Neogreen

``` r
FeaturePlot(object = esc, features.plot = c("Neogreen"), cols.use = c("grey", "red"), 
    reduction.use = "tsne")
```

![](Exploratory_analysis_seurat_better_version_files/figure-markdown_github/unnamed-chunk-35-1.png)

For stem cell markers:

``` r
FeaturePlot(object = esc, features.plot = c("SOX2", "NANOG", "CDH1", "DPPA4"), cols.use = c("grey", "red"), 
    reduction.use = "tsne")
```

![](Exploratory_analysis_seurat_better_version_files/figure-markdown_github/unnamed-chunk-36-1.png)

``` r
top10 <- esc.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
# setting slim.col.label to TRUE will print just the cluster IDS instead of
# every cell name
DoHeatmap(object = esc, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
```

![](Exploratory_analysis_seurat_better_version_files/figure-markdown_github/unnamed-chunk-37-1.png)

``` r
cluster2.markers <- FindMarkers(object = esc, ident.1 = 2, min.pct = 0.25)
print(x = head(x = cluster2.markers, n = 10))
```

    ##                   p_val avg_logFC pct.1 pct.2     p_val_adj
    ## LHX1      2.724672e-282 1.1646227 0.845 0.046 4.933563e-278
    ## GATA6     9.414329e-249 0.7685805 0.850 0.062 1.704653e-244
    ## FGF17     7.358184e-218 1.1774791 0.920 0.109 1.332346e-213
    ## GATA6-AS1 1.160331e-197 0.3480916 0.485 0.011 2.101012e-193
    ## APLNR     5.805311e-189 0.6684603 0.630 0.039 1.051168e-184
    ## LHX1-DT   7.194241e-160 0.2566843 0.490 0.023 1.302661e-155
    ## DKK1      4.881593e-157 2.2544177 0.940 0.232 8.839100e-153
    ## HAS2      1.085916e-154 1.0504047 0.920 0.196 1.966267e-150
    ## LIX1      1.501319e-148 0.3939590 0.370 0.009 2.718439e-144
    ## DKK4      3.923228e-148 1.6488767 0.765 0.118 7.103790e-144

``` r
cluster3.markers <- FindMarkers(object = esc, ident.1 = 3, min.pct = 0.25)
print(x = head(x = cluster3.markers, n = 10))
```

    ##              p_val  avg_logFC pct.1 pct.2    p_val_adj
    ## GAPDH 6.626015e-44 -0.3906786 1.000 1.000 1.199773e-39
    ## PKM   3.026758e-40 -0.6881946 0.968 0.999 5.480552e-36
    ## GSTP1 3.851879e-35 -0.3511026 1.000 1.000 6.974598e-31
    ## ACTB  7.582203e-35 -0.4411506 1.000 1.000 1.372909e-30
    ## PHGDH 1.419119e-34 -0.5887536 0.915 0.995 2.569599e-30
    ## EEF2  1.068775e-33 -0.3744439 1.000 1.000 1.935231e-29
    ## ENO1  4.267650e-32 -0.4713163 1.000 1.000 7.727434e-28
    ## PRKDC 1.715239e-31 -0.5255087 0.457 0.935 3.105783e-27
    ## SNRPN 4.969095e-31 -0.4269064 1.000 1.000 8.997540e-27
    ## TKT   1.277745e-28 -0.5018311 0.926 0.995 2.313614e-24

For testing:

``` r
current.cluster.ids <- c(0, 1, 2, 3)
new.cluster.ids <- c("hESC", "hESC", "APS", "CTNNB1+ cells")
x <- esc
x@ident <- plyr::mapvalues(x = x@ident, from = current.cluster.ids, to = new.cluster.ids)
TSNEPlot(object = x, do.label = TRUE, pt.size = 0.5, label.size = 6)
```

![](Exploratory_analysis_seurat_better_version_files/figure-markdown_github/unnamed-chunk-40-1.png)

Among negative markers for APS cluster are hESC markerks: SOX2, NANOG, POU5F1, DPPA4. For CTNNB1+ cluster there is only one negative marker - CER1. Probably, CTNNB1+ cells, are intermidiate step from hESC to APS.

Additional testing
------------------

Using all 100 dimensions reveals only two clusters:

``` r
esc <- FindClusters(object = esc, reduction.type = "pca", dims.use = 1:100, 
    resolution = 0.6, print.output = 0, save.SNN = TRUE)

esc <- RunTSNE(object = esc, dims.use = 1:100, do.fast = TRUE)

TSNEPlot(object = esc)
```

![](Exploratory_analysis_seurat_better_version_files/figure-markdown_github/unnamed-chunk-41-1.png)

``` r
table(esc@ident) 
```

    ## 
    ##    0    1    2 
    ## 1199 1153  159

``` r
PCAPlot(object = esc, dim.1 = 1, dim.2 = 2)
```

![](Exploratory_analysis_seurat_better_version_files/figure-markdown_github/unnamed-chunk-43-1.png)

``` r
cluster0.markers <- FindMarkers(object = esc, ident.1 = 0, min.pct = 0.25)
print(x = head(x = cluster0.markers, n = 10))
```

    ##                 p_val  avg_logFC pct.1 pct.2     p_val_adj
    ## HSPD1   1.690253e-233 -0.4658759 1.000 1.000 3.060542e-229
    ## ATP5F1B 2.496895e-196 -0.4620275 0.999 0.999 4.521128e-192
    ## CCT2    6.321846e-186 -0.4382594 0.997 0.999 1.144697e-181
    ## HDAC2   8.974925e-186 -0.4956369 0.988 1.000 1.625090e-181
    ## ENO1    3.405592e-178 -0.3318729 1.000 1.000 6.166506e-174
    ## HSP90B1 2.030271e-174 -0.5503003 0.971 0.996 3.676213e-170
    ## TUBA1B  3.417162e-173 -0.4629375 1.000 1.000 6.187456e-169
    ## KPNA2   5.571666e-169 -0.5758295 0.952 0.990 1.008862e-164
    ## ATP5F1A 6.214712e-169 -0.3653899 0.999 1.000 1.125298e-164
    ## PDIA3   1.000018e-167 -0.4756337 0.957 0.992 1.810733e-163

``` r
cluster1.markers <- FindMarkers(object = esc, ident.1 = 1, min.pct = 0.25)
print(x = head(x = cluster1.markers, n = 10))
```

    ##                 p_val avg_logFC pct.1 pct.2     p_val_adj
    ## HSPD1   1.109860e-238 0.4575322 1.000 1.000 2.009623e-234
    ## ATP5F1B 1.636397e-204 0.4614158 0.999 0.999 2.963024e-200
    ## CCT2    3.288691e-186 0.4274616 0.999 0.998 5.954833e-182
    ## ATP5F1A 2.310694e-182 0.3737473 1.000 0.999 4.183973e-178
    ## HDAC2   2.718383e-178 0.4704625 1.000 0.990 4.922176e-174
    ## KPNA2   1.269356e-171 0.5729576 0.994 0.954 2.298423e-167
    ## CCT5    8.008010e-160 0.4023006 0.999 0.997 1.450010e-155
    ## CCT8    5.174638e-156 0.3898005 1.000 0.993 9.369717e-152
    ## TUBA1B  6.920073e-156 0.4289507 1.000 1.000 1.253018e-151
    ## NPM1    1.128485e-155 0.2680451 1.000 1.000 2.043348e-151

``` r
#markers distinguishing cluster 1 from clusters 0 
cluster1.markers <- FindMarkers(object = esc, ident.1 = 1, ident.2 = c(0), 
    min.pct = 0.25)
print(x = head(x = cluster1.markers, n = 10))
```

    ##                 p_val avg_logFC pct.1 pct.2     p_val_adj
    ## HSPD1   8.473833e-251 0.4914230 1.000 1.000 1.534357e-246
    ## ATP5F1B 3.469947e-211 0.4906356 0.999 0.999 6.283034e-207
    ## CCT2    2.443883e-198 0.4610083 0.999 0.997 4.425138e-194
    ## HDAC2   4.024502e-193 0.5159402 1.000 0.988 7.287165e-189
    ## ATP5F1A 2.216488e-186 0.3920998 1.000 0.999 4.013395e-182
    ## KPNA2   3.688376e-178 0.6103584 0.994 0.952 6.678543e-174
    ## TUBA1B  3.210953e-173 0.4772507 1.000 1.000 5.814073e-169
    ## ENO1    6.004657e-173 0.3361365 1.000 1.000 1.087263e-168
    ## PDIA3   3.851097e-170 0.4916431 0.993 0.957 6.973181e-166
    ## CCT5    1.070121e-169 0.4334976 0.999 0.997 1.937669e-165

``` r
#markers distinguishing cluster 0 from clusters 1
cluster0.markers <- FindMarkers(object = esc, ident.1 = 0, ident.2 = c(1), 
    min.pct = 0.25)
print(x = head(x = cluster0.markers, n = 10))
```

    ##                 p_val  avg_logFC pct.1 pct.2     p_val_adj
    ## HSPD1   8.473833e-251 -0.4914230 1.000 1.000 1.534357e-246
    ## ATP5F1B 3.469947e-211 -0.4906356 0.999 0.999 6.283034e-207
    ## CCT2    2.443883e-198 -0.4610083 0.997 0.999 4.425138e-194
    ## HDAC2   4.024502e-193 -0.5159402 0.988 1.000 7.287165e-189
    ## ATP5F1A 2.216488e-186 -0.3920998 0.999 1.000 4.013395e-182
    ## KPNA2   3.688376e-178 -0.6103584 0.952 0.994 6.678543e-174
    ## TUBA1B  3.210953e-173 -0.4772507 1.000 1.000 5.814073e-169
    ## ENO1    6.004657e-173 -0.3361365 1.000 1.000 1.087263e-168
    ## PDIA3   3.851097e-170 -0.4916431 0.957 0.993 6.973181e-166
    ## CCT5    1.070121e-169 -0.4334976 0.997 0.999 1.937669e-165

GO analysis
-----------

For GO analysis with GOrilla of markers of cluster 0 and 1:

``` r
x <- rownames(cluster0.markers)
y <- rownames(cluster1.markers) #they are equal

#bg list
z <- rownames(as.matrix(esc@data))

#write to file
fileConn<-file("Bg.txt")
writeLines(z, fileConn)
close(fileConn)
```

We have enriched protein folding, protein stabiliation, telomerase stability, toxin transport and others.

Continue with analysis
----------------------

``` r
library(tibble)
cluster3.markers <- FindMarkers(object = esc, ident.1 = 2, thresh.use = 0.25, 
                                 test.use = "roc", only.pos = TRUE)
cluster3.markers %>% rownames_to_column("Genes") %>% filter(power > 0.8) %>% dplyr::select("Genes")
```

    ##      Genes
    ## 1     CER1
    ## 2    MIXL1
    ## 3     DKK1
    ## 4     HAS2
    ## 5    FGF17
    ## 6     LHX1
    ## 7    ARL4D
    ## 8    GATA6
    ## 9    EOMES
    ## 10 CYP26A1
    ## 11    DKK4

Main markers.

RAZBERIS!!! - GSEA. PLUS WHAT HAPPENS AGAIN IF YOU INCLUDE RIB AND MITO GENES?

``` r
cluster2.markers <- FindMarkers(object = esc, ident.1 = 2, min.pct = 0.25)
print(x = head(x = cluster2.markers, n = 10))
```

    ##                   p_val avg_logFC pct.1 pct.2     p_val_adj
    ## LHX1      4.209939e-264 1.2861819 0.906 0.056 7.622937e-260
    ## GATA6     8.794926e-234 0.8482647 0.912 0.071 1.592497e-229
    ## GATA6-AS1 8.142031e-213 0.4077397 0.560 0.014 1.474278e-208
    ## FGF17     4.266263e-197 1.2339396 0.969 0.119 7.724922e-193
    ## APLNR     1.072462e-188 0.7498425 0.698 0.045 1.941907e-184
    ## LHX1-DT   4.733635e-178 0.3040271 0.572 0.025 8.571194e-174
    ## LIX1      1.234637e-175 0.4706836 0.447 0.010 2.235557e-171
    ## DKK4      1.786329e-148 1.7774205 0.836 0.125 3.234506e-144
    ## FAM24B    4.539038e-148 0.2967270 0.421 0.013 8.218836e-144
    ## HAS2      4.605923e-148 1.1611071 0.975 0.205 8.339944e-144

``` r
esc <- FindClusters(object = esc, reduction.type = "pca", dims.use = 1:20, 
    resolution = 0.6, print.output = 0, save.SNN = TRUE)

esc <- RunTSNE(object = esc, dims.use = 1:20, do.fast = TRUE)

TSNEPlot(object = esc)
```

![](Exploratory_analysis_seurat_better_version_files/figure-markdown_github/unnamed-chunk-51-1.png)
