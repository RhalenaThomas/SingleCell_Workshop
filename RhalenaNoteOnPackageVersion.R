# specify Seurat 5 for this session 


# View current library paths
.libPaths()


# Check the current library paths
old_libpaths <- .libPaths()

# Set your primary library path
second_libpath <- "/Users/rhalenathomas/Library/CustomR"
.libPaths(c(old_libpaths, second_libpath))


#[1] "/Library/Frameworks/R.framework/Versions/4.2/Resources/library"


#My pathway

#/Users/rhalenathomas/Library/CustomR

# Install Seurat 5 to a separate library location
install.packages("Seurat", lib = "/Users/rhalenathomas/Library/CustomR")
install.packages("SeuratObject", lib = "/Users/rhalenathomas/Library/CustomR")


# Load Seurat 5 from the custom library location
library(Seurat, lib.loc = "/Users/rhalenathomas/Library/CustomR/")


# Check Seurat version
packageVersion("Seurat")
packageVersion("SeuratObject")
