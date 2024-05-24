# specify Seurat 5 for this session 

# I have not specified package versions - all packages are the version on May 17, 2024

# View current library paths
.libPaths()


# Check the current library paths
old_libpaths <- .libPaths()

# Set your primary library path
second_libpath <- "/Users/rhalenathomas/Library/CustomR"
.libPaths(c(old_libpaths, second_libpath))
.libPaths(c(second_libpath,old_libpaths))
#### if you want to install everything in to your normal library 
# run 
#second_libpath <- .libPaths()

# then you won't need to change the code below 

# Install Seurat 5 to a separate library location
install.packages("Seurat", lib = second_libpath)
install.packages("SeuratObject", lib = second_libpath)
install.packages("tidyverse", lib = second_libpath)
install.packages("remotes", lib = second_libpath)
install.packages("ggplot2", lib = second_libpath)
install.packages("patchwork", lib = second_libpath)
install.packages("clustree", lib = second_libpath)
install.packages("anndata", lib = second_libpath)


remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', lib = second_libpath)  ### required some kind of password chain on my mac 
remotes::install_github("wjawaid/enrichR", lib = second_libpath)
remotes::install_github("SydneyBioX/scClassify", lib = second_libpath)


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SingleCellExperiment", lib = second_libpath, force = TRUE)


# note "Matrix" should be installed from base R but it isn't in my separate library. 

library(SingleCellExperiment, lib = second_libpath)

library(readr)
library(MetaNeighbor)
library(anndata)
library(Matrix)
library(scclusteval)
library(mclust)
library(MAST)
library(EnhancedVolcano)
library(DESeq2)








# Load Seurat 5 from the custom library location
library(Seurat, lib.loc = second_libpath)
library(tidyverse, lib.loc = second_libpath)
library(remotes, lib.loc = second_libpath)
library(ggplot2, lib.loc = second_libpath)
library(patchwork, lib = second_libpath)
library(clustree, lib.loc = second_libpath)
library(DoubletFinder, lib.loc = second_libpath)
library(enrichR, lib = second_libpath)
library(scClassify, lib = second_libpath)
library(anndata, lib = second_libpath)
library(Matrix)

# Check Seurat version
packageVersion("Seurat")
packageVersion("SeuratObject")

