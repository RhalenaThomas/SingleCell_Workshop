#Press "Run" on top right corner to install the packages

#Install packages from CRAN (https://cran.r-project.org/)
if (!require("tidyverse", quietly = TRUE))  install.packages("tidyverse")
if (!require("reticulate", quietly = TRUE))  install.packages("reticulate")
reticulate::install_miniconda(path = reticulate::miniconda_path(), update = TRUE, force = FALSE)
if (!require("Seurat", quietly = TRUE))  install.packages("Seurat")
if (!require("clustree", quietly = TRUE))  install.packages("clustree")
if (!require("devtools", quietly = TRUE))  install.packages("devtools")
if (!require("enrichR", quietly = TRUE))  install.packages("enrichR")
if (!require("Matrix", quietly = TRUE))  install.packages("Matrix")
if (!require("anndata", quietly = TRUE))  install.packages("anndata")

#Install packages from Bionconductor (https://www.bioconductor.org/)
if (!require("BiocManager", quietly = TRUE))  install.packages("BiocManager")
if (!require("scClassify", quietly = TRUE))  BiocManager::install("scClassify")
if (!require("SingleCellExperiment", quietly = TRUE))  BiocManager::install("SingleCellExperiment")
if (!require("MetaNeighbor", quietly = TRUE))  BiocManager::install("MetaNeighbor")

#Install packages from Github
if (!require("DoubletFinder", quietly = TRUE))  remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
if (!require("SeuratDisk", quietly = TRUE))  remotes::install_github("mojaveazure/seurat-disk")
