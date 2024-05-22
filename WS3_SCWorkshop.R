## Workshop 3

## load
library(SingleCellExperiment)
library(Seurat)
library(readr, lib="/home/fiorini9/scratch/mjff/MF/practice_ensemblux/R")
library(tidyverse, lib="/home/fiorini9/scratch/mjff/MF/practice_ensemblux/R")
library(patchwork)
library(DoubletFinder)
library(ggplot2)
library(MetaNeighbor)
library(anndata)
library(remotes)
library(Matrix)
library(scclusteval)
library(mclust)
library(MAST, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
library(EnhancedVolcano, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")

```
Here we will run MetaNeighbor
https://github.com/maggiecrow/MetaNeighbor

Read in and properly format the data for the Lister lab. This is an example of using data in h5ad format with Seurat.

####################################
## Rhalena, I dont have enough space on my computer to to download the Lister object so I wasnt able to test this first section. This is Malrosee's code from last year. 
## If for some reason the code doesn't work we can think about just taking the trained object which is shown below.
## Looking at the code however I think it will work
####################################

```{r}
#Load as anndata, and look at the object
lister_data <- read_h5ad("data/ListerLab/Processed_data_RNA-all_full-counts-and-downsampled-CPM.h5ad")
lister_data

#Convert the basics to Seurat format, look at the object
lister_data <- CreateSeuratObject(counts = t(lister_data$X), meta.data = lister_data$obs)
```

Look at the Lister data object

```{r}
lister_data

# Look at the major cell types and clusters, also distribution of genes and UMIs
lister_data$cell_type %>% table
lister_data$major_clust %>% table
VlnPlot(lister_data, features = c("nCount_RNA", "nFeature_RNA"), pt.size = 0.001)
lister_data$nCount_RNA %>% summary
lister_data$nFeature_RNA %>% summary

```

Now we will train some classifiers using this object as the reference, but first we will do some pre-processing. 

First we wiill build the classifier using all the cells, and then we will re-run the analysis using only the non-neuronal cells. 

```{r}
Normalize, find variable features, convert to SingleCellExperiment and subset to only the variable features 
lister_data <- subset(lister_data, cell_type == "Non-Neu")
lister_data <- NormalizeData(lister_data, normalization.method = "LogNormalize", scale.factor = 10000)
lister_data <- FindVariableFeatures(lister_data, selection.method = "vst", nfeatures = 3000)
var_features <- VariableFeatures(lister_data)
Idents(lister_data) <- lister_data$major_clust
VlnPlot(lister_data, features = c("SNAP25", "GAD1", "SLC17A7"), pt.size = -1)
VlnPlot(lister_data, features = c("ALDH1L1", "MRC1"), pt.size = -1)
VlnPlot(lister_data, features = c("PLP1", "CLDN5"), pt.size = -1)
lister_data <- as.SingleCellExperiment(lister_data)
lister_data <- lister_data[var_features,]
lister_data
```

Time to train the model!

```{r}
pretrained_model_major_cluster = MetaNeighbor::trainModel(
var_genes = rownames(lister_data),
dat = lister_data,
study_id = rep("Lister", dim(lister_data)[2]),
cell_type = lister_data$major_clust)
```

Optional save trained model or load trained mode


```{r}
# training the model is computationally heavy

# Save model
saveRDS(pretrained_model_major_cluster, "data/pretainedModelLister.Rds")
```

##############
# Rhalena I tested everything from here on
##############

Load the pretrained model 
``` {r}
# load model
pretrained_model_major_cluster <- readRDS("/home/fiorini9/scratch/SCWorskshop/pretainedModelLister.Rds")

```

Remove the Lister lab object, and load our integrated Seurat object

```{r}
rm(lister_data)
```
## load the Seurat object
```{r}
## load the Seurat object
integrated_seurat_reg <- readRDS("/home/fiorini9/scratch/SCWorskshop/integrated_seurat_clusters.rds")
```

 Examine the expression data
```{r}
GetAssayData(integrated_seurat_reg, assay = "RNA", slot = "data")  %>% colSums %>% head
GetAssayData(integrated_seurat_reg, assay = "RNA", slot = "counts")  %>% colSums %>% head
```

One versus all
## Rhalena: Change [RNA_snn_res.0.2] to you actual annotation column
```{r}
aurocs = MetaNeighborUS(
trained_model = pretrained_model_major_cluster, dat = as.SingleCellExperiment(integrated_seurat_reg),
study_id = rep("integrated", dim(integrated_seurat_reg)[2]), 
cell_type = integrated_seurat_reg@meta.data$RNA_snn_res.0.2,
fast_version = TRUE
)
## print plot
pdf(paste('/home/fiorini9/scratch/SCWorskshop/temp.pdf', sep=""), width=7, height=5) 
tryCatch({plotHeatmapPretrained(aurocs, margins = c(5,5))}, error = function(error_condition) {})
dev.off()
```

One versus best
## Rhalena: Change [RNA_snn_res.0.2] to you actual annotation column
```{r}
best_hits = MetaNeighborUS(
trained_model = pretrained_model_major_cluster, dat = as.SingleCellExperiment(integrated_seurat_reg),
study_id = rep("integrated", dim(integrated_seurat_reg)[2]), 
cell_type = integrated_seurat_reg@meta.data$RNA_snn_res.0.2,
one_vs_best = TRUE,
fast_version = TRUE
)

## print plot
pdf(paste('/home/fiorini9/scratch/SCWorskshop/temp.pdf', sep=""), width=7, height=5) 
tryCatch({plotHeatmapPretrained(best_hits, margins = c(5,5))}, 
         error = function(error_condition) {})
dev.off()
```

Compare labels obtained with supervised classification of cells and unsupervised clustering, also between the two supervised clustering labels

########
## Rhalena I dont know if you added these datasets in your section or not? I saw that you did last year?
## I didnt run this because you hadn't done it yet when I wrote the code
########
```{r}
 plot1 <- PairWiseJaccardSetsHeatmap(Idents(integrated_seurat_reg), integrated_seurat@meta.data$Karolinski)
 adjustedRandIndex(Idents(integrated_seurat_reg), integrated_seurat_reg@meta.data$Karolinski) %>% print
 
 plot2 <- PairWiseJaccardSetsHeatmap(Idents(integrated_seurat_reg), integrated_seurat_reg@meta.data$Nowakowski)
 adjustedRandIndex(Idents(integrated_seurat_reg), integrated_seurat_reg@meta.data$Nowakowski) %>% print
 
 plot3 <- PairWiseJaccardSetsHeatmap(integrated_seurat_reg@meta.data$Karolinski, integrated_seurat_reg@meta.data$Nowakowski)
 adjustedRandIndex(integrated_seurat_reg@meta.data$Karolinski, integrated_seurat_reg@meta.data$Nowakowski) %>% print
```

``` {r JaccardIndexPlots, dev = c("png")}
 #pdf(file = "C:/Users/Home/Documents/GitHub/SingleCell_Workshop/JaccardIndexPlots.pdf", onefile = TRUE)
  plot1
  plot2
  plot3  
 #dev.off()
```
Assess doublets and cell cycle genes per annotated

``` {r}
table(integrated_seurat_reg@meta.data$age, integrated_seurat_reg@meta.data$DF.classifications_0.25_0.09_533)
table(integrated_seurat_reg@meta.data$age, integrated_seurat_reg@meta.data$DF.classifications_0.25_0.17_146)

#Make one variable with all the doublet classifications
integrated_seurat_reg@meta.data$DF_all <- integrated_seurat_reg@meta.data$DF.classifications_0.25_0.09_533
integrated_seurat_reg@meta.data$DF_all[integrated_seurat_reg@meta.data$age == 41] <- integrated_seurat_reg@meta.data$DF.classifications_0.25_0.17_146[integrated_seurat_reg@meta.data$age == 41]

#Make one variable with all the pANN scores
integrated_seurat_reg@meta.data$pANN_all <- integrated_seurat_reg@meta.data$pANN_0.25_0.09_533
integrated_seurat_reg@meta.data$pANN_all[integrated_seurat_reg@meta.data$age == 41] <- integrated_seurat_reg@meta.data$pANN_0.25_0.17_146[integrated_seurat_reg@meta.data$age == 41]

#Plot doublet classifications and pANN scores
table(integrated_seurat_reg@meta.data$DF_all)

DimPlot(integrated_seurat_reg, reduction = "umap.rpca", group.by = c("DF_all" ), combine = FALSE)
ggsave(paste('/home/fiorini9/scratch/SCWorskshop/temp.pdf', sep=""), height = 4, width = 12)

VlnPlot(integrated_seurat_reg, features = c("pANN_all"), pt.size = 0.001)
ggsave(paste('/home/fiorini9/scratch/SCWorskshop/temp.pdf', sep=""), height = 4, width = 12)
```

```{r}
#Plot CellCycle scores per cluster
VlnPlot(integrated_seurat_reg, features = c("S.Score", "G2M.Score"))
ggsave(paste('/home/fiorini9/scratch/SCWorskshop/temp.pdf', sep=""), height = 4, width = 12)

DimPlot(integrated_seurat_reg, reduction = "umap.rpca", group.by = c("Phase" ), combine = FALSE)
ggsave(paste('/home/fiorini9/scratch/SCWorskshop/temp.pdf', sep=""), height = 4, width = 12)
```

########   Differential gene expression with MAST    ########

Run interactive example with the group using the smallest cell type

We will begin by running differential gene expression on one of the cell types interactively. Later, we will run a loop to run DEG for all of the cell types. 

## RHalena: change "RNA_snn_res.0.2" to the column that has you final annotations

````{r}
## set default assay to RNA
DefaultAssay(integrated_seurat_reg) <- "RNA"

## Unique cell types 
unique(integrated_seurat_reg@meta.data$RNA_snn_res.0.2)

## Number of cells for each cell types
table(integrated_seurat_reg@meta.data$RNA_snn_res.0.2)

##unique subject groups: age
unique(integrated_seurat_reg@meta.data$age)

## Number of cells for each subject group
table(integrated_seurat_reg@meta.data$age)

## Number of cell types per subject group
table(integrated_seurat_reg@meta.data$RNA_snn_res.0.2, integrated_seurat_reg@meta.data$age)


## Set the idents of the cell
Idents(integrated_seurat_reg) <-  "RNA_snn_res.0.2"

## Check the idents
head(Idents(integrated_seurat_reg))

## Subset the Seurat object to only include our initial cell type of interest
celltype.sub.seu <- subset(integrated_seurat_reg, idents = 7)

## Check the idents of the subsetted Seurat object
head(Idents(celltype.sub.seu))

## Set the idents of the subsetted Seurat object to the metadata column that defines our subject groups
Idents(celltype.sub.seu) <- "age"  

## Check the new idents of the subsetted Seurat object
head(Idents(celltype.sub.seu))

## Run Differential gene expression  
DGE <- FindMarkers(celltype.sub.seu, ident.1 = 14, ident.2 = 41,  logfc.threshold = 0, test.use = "MAST")
  
## Examine the DEG results
head(DGE)

## plot Volcano plot
  #modify colour scheme
  DGE$col <- "lightgrey"
  DGE$col[DGE$avg_log2FC > 1 & DGE$p_val < 0.05] <- "indianred3"
  DGE$col[DGE$avg_log2FC < -1 & DGE$p_val < 0.05] <- "dodgerblue2"
  DGE$col[DGE$avg_log2FC > 1 & DGE$p_val_adj < 0.05] <- "red4"
  DGE$col[DGE$avg_log2FC < -1 & DGE$p_val_adj < 0.05] <- "navy"

  # fix cols
  keyvals <- DGE$col
  names(keyvals)[keyvals == 'indianred3'] <- 'Log2FC > 1; P < 0.05'
  names(keyvals)[keyvals == 'lightgrey'] <- 'NS'
  names(keyvals)[keyvals == 'dodgerblue2'] <- 'Log2FC< -1; P < 0.05'
  names(keyvals)[keyvals == 'red4'] <- 'Log2FC > 1; adj. P < 0.05'
  names(keyvals)[keyvals == 'navy'] <- 'Log2FC < -1; adj. P < 0.05'

  # volcano plot
  # volcano plot
  ## print volcano plot
      EnhancedVolcano(DGE,
                      lab = rownames(DGE),
                      x = 'avg_log2FC',
                      y = 'p_val',
                      pCutoff = 0.05,
                      FCcutoff = 1,
                      pointSize = 2.0,
                      labSize = 5, 
                      colCustom = keyvals,
                      parseLabels = TRUE,
                      drawConnectors = TRUE,
                      colAlpha = 1,
      ) +
      theme_classic() + 
      theme(plot.subtitle=element_blank(),
            plot.title=element_blank(),
            legend.position = "right",
            legend.title=element_blank())

ggsave(paste('/home/fiorini9/scratch/SCWorskshop/temp.pdf', sep=""), height = 4, width = 12)
```

Now we will run a loop to compute DGE for all of the cell types 
## RHalena: change RNA_snn_res.0.2 to your final annotation

```{r}
  ## set default assay to RNA
  DefaultAssay(integrated_seurat_reg) <- "RNA"

  ## Unique cell types 
  unique(integrated_seurat_reg@meta.data$RNA_snn_res.0.2)
  
## look into the metadata of the Seurat object
str(integrated_seurat_reg@meta.data)
unique(integrated_seurat_reg@meta.data$age)

  
  ## set up our parameters
  contrast_name <- paste0("DEG_", unique(integrated_seurat_reg@meta.data$RNA_snn_res.0.2[!is.na(integrated_seurat_reg@meta.data$RNA_snn_res.0.2)])) #### need to modify this.
  meta_data_celltype <- rep("RNA_snn_res.0.2", length(contrast_name))
  cell_type <-  unique(integrated_seurat_reg@meta.data$RNA_snn_res.0.2[!is.na(integrated_seurat_reg@meta.data$RNA_snn_res.0.2)])
  meta_data_variable <- rep("age", length(contrast_name))
  group1 <- rep(14, length(contrast_name))
  group2 <- rep(41, length(contrast_name))


  ## create a parameter data frame
  par_df <- data.frame(contrast_name,
                        meta_data_celltype,
                        cell_type,
                        meta_data_variable, 
                        group1,
                        group2)
 
list = list()

for(i in 1:nrow(par_df)){
  Idents(integrated_seurat_reg) <- par_df[i,2]    
  celltype.sub.seu <- subset(integrated_seurat_reg, idents = par_df[i,3])
  Idents(celltype.sub.seu) <- par_df[i,4]    
  
  DGE <- FindMarkers(celltype.sub.seu, ident.1 = par_df[i,5], ident.2 = par_df[i,6],  logfc.threshold = 0, test.use = "MAST")
  DGE$DEG_contrast <- par_df[i,1]
  list <- c(list, list(DGE))
}
  
combined_df <- do.call(rbind, list)
```

Now we will plot our results as a volcano plot
```{r}
#modify colour scheme
DGE <- combined_df
  DGE$col <- "lightgrey"
  DGE$col[DGE$avg_log2FC > 1 & DGE$p_val < 0.05] <- "indianred3"
  DGE$col[DGE$avg_log2FC < -1 & DGE$p_val < 0.05] <- "dodgerblue2"
  DGE$col[DGE$avg_log2FC > 1 & DGE$p_val_adj < 0.05] <- "red4"
  DGE$col[DGE$avg_log2FC < -1 & DGE$p_val_adj < 0.05] <- "navy"

  # fix cols
  keyvals <- DGE$col
  names(keyvals)[keyvals == 'indianred3'] <- 'Log2FC > 1; P < 0.05'
  names(keyvals)[keyvals == 'lightgrey'] <- 'NS'
  names(keyvals)[keyvals == 'dodgerblue2'] <- 'Log2FC< -1; P < 0.05'
  names(keyvals)[keyvals == 'red4'] <- 'Log2FC > 1; adj. P < 0.05'
  names(keyvals)[keyvals == 'navy'] <- 'Log2FC < -1; adj. P < 0.05'

  # volcano plot
  # volcano plot
  ## print volcano plot
      EnhancedVolcano(DGE,
                      lab = NA,
                      x = 'avg_log2FC',
                      y = 'p_val',
                      pCutoff = 0.05,
                      FCcutoff = 1,
                      pointSize = 2.0,
                      labSize = 5, 
                      colCustom = keyvals,
                      parseLabels = TRUE,
                      drawConnectors = TRUE,
                      colAlpha = 1,
      ) +
      theme_classic() + 
      theme(plot.subtitle=element_blank(),
            plot.title=element_blank(),
            legend.position = "right",
            legend.title=element_blank()) +
    facet_wrap(~DEG_contrast, scales = "free_y")

ggsave(paste('/home/fiorini9/scratch/SCWorskshop/temp.pdf', sep=""), height = 4, width = 12)
```


  
  

  