# install required packages
# R version must be  >= 4.2

# install devtools

# install bioconductor/CRAN packages
install.packages("tidyverse")
install.packages("Seurat")
install.packages("remotes")
install.packages("ggplot2")

# install packages from github
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
install_github("wjawaid/enrichR")
library(devtools)
devtools::install_github("SydneyBioX/scClassify")


# must install harmony

install.packages("harmony")
install.packages("patchwork")
install.packages("clustree")
