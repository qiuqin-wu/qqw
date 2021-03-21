# qqw
---
  title: "scRNAseq analysis using Rmarkdown"
author: "Cankun Wang"
date: "`r format(Sys.time(), '%m/%d/%Y')`"
output: 
  html_document:
  toc: true
toc_float: true
number_sections: true
---
 ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)

cran_packages <- c(
  "Seurat","enrichR","DT"
)
cran_np <- cran_packages[!(cran_packages %in% installed.packages()[,"Package"])]

if(length(cran_np)) {
  install.packages(cran_np)
}


library(Seurat) 
library(enrichR)
library(DT)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library("org.Hs.eg.db")
library(clusterProfiler)

library(ggplot2)
library(ChIPseeker)
knitr::opts_knit$set(root.dir = "C:/Users/qiuqi/Desktop/hi-c")

# setwd()
```
# How to use this repository

GitHub for this templtate: [Click here to open](https://github.com/Wang-Cankun/scRNAseq-template)

## Step 1: Download the whole repository.

![Download repo](./tutorial/s1.png)

## Step 2: Unzip files to your local directory.

You will find these files in your directory

- count-data-scrna.csv: example gene expression matrix file.
- col-data-scrna.csv: example gene metadata file.
- scrnaseq_template_html.rmd: Rmarkdown template to generate HTML report.
- scrnaseq_template_html.html: Example HTML report.
- scrnaseq_template_pdf.rmd: Rmarkdown template to generate PDF report.
- scrnaseq_template_pdf.html: Example PDF report.
- peaks.bed: example ATAC peaks file.

## Step 3 (in Rmarkdown template): Change the working directory to your local directory.
![](./tutorial/s3.png)


# Goal

The primary goal of this study is to identify gene expression difference between conditions using a toy dataset.

Rmarkdown reference can be found in here: [Reference](https://rstudio.com/wp-content/uploads/2015/03/rmarkdown-reference.pdf)

[Mathematics in Rmarkdown](https://www.calvin.edu/~rpruim/courses/s341/S17/from-class/MathinRmd.html)

[Seurat - Guided Clustering Tutorial](https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html)

# Run Seurat strandard workflow

## Load dataset

Using raw read count as input.


```{r,echo=F,eval=T,message=F,warning=F,error=F}
## list all enrichr database
# listEnrichrDbs()

## Change to human or mouse 
dbs <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018","KEGG_2019_Human")

counts <- read.csv("count-data-scrna.csv",stringsAsFactors = F,check.names = F, row.names = 1)
meta <- read.csv("col-data-scrna.csv",stringsAsFactors = F, row.names = 1)
peak<-readPeakFile("peaks.bed")

# Overview of the dataset
counts[1:10, ]
meta[1:10,]

```

## Load dataset in Seurat

```{r,echo=T,eval=T,message=F,warning=F,error=F}
# choose whahtever column as design
obj <- CreateSeuratObject(counts = counts, project = "scRNAseq", min.cells = 3, min.features = 200)
obj <- AddMetaData(object = obj, metadata = meta)

```

## Normalizing the data

After removing unwanted cells from the dataset, the next step is to normalize the data. By default, we employ a global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result. Normalized values are stored in  obj[["RNA"]]@data.

```{r,echo=T,eval=T,message=T,warning=F,error=F}
obj <- NormalizeData(obj)
```

## Identification of highly variable features (feature selection)

We next calculate a subset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others). We and others have found that focusing on these genes in downstream analysis helps to highlight biological signal in single-cell datasets.

Our procedure in Seurat3 is described in detail here, and improves on previous versions by directly modeling the mean-variance relationship inherent in single-cell data, and is implemented in the FindVariableFeatures function. By default, we return 2,000 features per dataset. These will be used in downstream analysis, like PCA.

```{r,echo=T,eval=T,message=T,warning=F,error=F}
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)

```


## Scaling the data

Next, we apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA. The ScaleData function:
  
  * Shifts the expression of each gene, so that the mean expression across cells is 0
* Scales the expression of each gene, so that the variance across cells is 1
- This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
* The results of this are stored in obj[["RNA"]]@scale.data


```{r,echo=T,eval=T,message=T,warning=F,error=F}
all.genes <- rownames(obj)
obj <- ScaleData(obj, features = all.genes)
```


## Perform linear dimensional reduction (PCA)

Next we perform PCA on the scaled data. By default, only the previously determined variable features are used as input, but can be defined using features argument if you wish to choose a different subset.

```{r,echo=T,eval=T,message=T,warning=F,error=F}
obj <- RunPCA(obj, features = VariableFeatures(object = obj))

```


## Determine the ‘dimensionality’ of the dataset

To overcome the extensive technical noise in any single feature for scRNA-seq data, Seurat clusters cells based on their PCA scores, with each PC essentially representing a ‘metafeature’ that combines information across a correlated feature set. The top principal components therefore represent a robust compression of the dataset. However, how many componenets should we choose to include? 10? 20? 100?
  
  In Macosko et al, we implemented a resampling test inspired by the JackStraw procedure. We randomly permute a subset of the data (1% by default) and rerun PCA, constructing a ‘null distribution’ of feature scores, and repeat this procedure. We identify ‘significant’ PCs as those who have a strong enrichment of low p-value features.

```{r,echo=T,eval=T,message=T,warning=F,error=F}
ElbowPlot(obj)
```


## Cluster the cells

Seurat v3 applies a graph-based clustering approach, building upon initial strategies in (Macosko et al). Importantly, the distance metric which drives the clustering analysis (based on previously identified PCs) remains the same. However, our approach to partioning the cellular distance matrix into clusters has dramatically improved. Our approach was heavily inspired by recent manuscripts which applied graph-based clustering approaches to scRNA-seq data [SNN-Cliq, Xu and Su, Bioinformatics, 2015] and CyTOF data [PhenoGraph, Levine et al., Cell, 2015]. Briefly, these methods embed cells in a graph structure - for example a K-nearest neighbor (KNN) graph, with edges drawn between cells with similar feature expression patterns, and then attempt to partition this graph into highly interconnected ‘quasi-cliques’ or ‘communities’.

As in PhenoGraph, we first construct a KNN graph based on the euclidean distance in PCA space, and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity). This step is performed using the FindNeighbors function, and takes as input the previously defined dimensionality of the dataset (first 10 PCs).

To cluster the cells, we next apply modularity optimization techniques such as the Louvain algorithm (default) or SLM [SLM, Blondel et al., Journal of Statistical Mechanics], to iteratively group cells together, with the goal of optimizing the standard modularity function. The FindClusters function implements this procedure, and contains a resolution parameter that sets the ‘granularity’ of the downstream clustering, with increased values leading to a greater number of clusters. We find that setting this parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells. Optimal resolution often increases for larger datasets. The clusters can be found using the Idents function.

## PCA plot colored by time

```{r,echo=T,eval=T,message=T,warning=F,error=F}
obj <- FindNeighbors(obj, dims = 1:10)
obj <- FindClusters(obj, resolution = 0.5)

```

## Run non-linear dimensional reduction (UMAP/tSNE)

```{r,echo=T,eval=T,message=T,warning=F,error=F}
obj <- RunUMAP(obj, dims = 1:10)
DimPlot(obj, reduction = "umap")
```

## Plot UMAP colored by custom metadata

### UMAP plot colored by cell types
```{r,echo=T,eval=T,message=T,warning=F,error=F}
Idents(obj) <- obj$Cluster
DimPlot(obj, reduction = "umap")

```

### UMAP plot colored by time
```{r,echo=T,eval=T,message=T,warning=F,error=F}
Idents(obj) <- obj$Time
DimPlot(obj, reduction = "umap")

```

### UMAP plot colored by condition
```{r,echo=T,eval=T,message=T,warning=F,error=F}
Idents(obj) <- obj$Condition
DimPlot(obj, reduction = "umap")

```


## Finding differentially expressed features

Find cell type specific genes between control and disease group

```{r,echo=T,eval=T,message=F,warning=F,error=F}
Idents(obj) <- obj$Condition
cts_de_genes <- FindMarkers(obj, ident.1 = "Control", ident.2 = "Disease")

DT::datatable(cts_de_genes, extensions = c('FixedColumns','Buttons'),
              options = list(
                pageLength = 5,
                scrollX = TRUE,
                scrollCollapse = TRUE,
                dom = 'Bfrtip',
                buttons = c('copy', 'csv', 'excel')
              ))
```

### Violin plot for top 10 DEGs


```{r,echo=F,eval=T,message=F,warning=F,error=F, fig.width=8, fig.height=18}
Idents(obj) <- obj$Condition
VlnPlot(obj, features = rownames(cts_de_genes[which(cts_de_genes$p_val_adj < 0.05),])[1:10])

```

## Functional enrichment analysis for DEGs {.tabset}

### GO Biological Process

```{r, fig.width=10, fig.height=10,echo=T,message=FALSE,warning=F}

# This select genes of adj.p.value < 0.05, sometimes people use different threshold, like adding log2foldchange threshold.


enriched_combined <- enrichr(rownames(cts_de_genes[which(cts_de_genes$p_val_adj < 0.05 & abs(cts_de_genes$avg_logFC) > 1.5),]),dbs)

# output top 20 enriched terms
DT::datatable(head(enriched_combined$GO_Biological_Process_2018,n=20)[,c(-3,-5,-6,-7)], extensions = c('FixedColumns','Buttons'),
              options = list(
                pageLength = 5,
                scrollX = TRUE,
                scrollCollapse = TRUE,
                dom = 'Bfrtip',
                buttons = c('copy', 'csv', 'excel')
              ))
```

### GO Cellular Component

```{r, fig.width=10, fig.height=10,echo=T,message=FALSE,warning=F}

DT::datatable(head(enriched_combined$GO_Cellular_Component_2018,n=20)[,c(-3,-5,-6,-7)], extensions = c('FixedColumns','Buttons'),
              options = list(
                pageLength = 5,
                scrollX = TRUE,
                scrollCollapse = TRUE,
                dom = 'Bfrtip',
                buttons = c('copy', 'csv', 'excel')
              ))
```

### GO Molecular Function

```{r, fig.width=10, fig.height=10,echo=T,message=FALSE,warning=F}

DT::datatable(head(enriched_combined$GO_Molecular_Function_2018,n=20)[,c(-3,-5,-6,-7)], extensions = c('FixedColumns','Buttons'),
              options = list(
                pageLength = 5,
                scrollX = TRUE,
                scrollCollapse = TRUE,
                dom = 'Bfrtip',
                buttons = c('copy', 'csv', 'excel')
              ))
```

### KEGG pathway

```{r, fig.width=10, fig.height=10,echo=T,message=FALSE,warning=F}

DT::datatable(head(enriched_combined$KEGG_2019_Mouse,n=20)[,c(-3,-5,-6,-7)], extensions = c('FixedColumns','Buttons'),
              options = list(
                pageLength = 5,
                scrollX = TRUE,
                scrollCollapse = TRUE,
                dom = 'Bfrtip',
                buttons = c('copy', 'csv', 'excel')
              ))
```


# Session Infomation

```{r}
sessionInfo()
```
