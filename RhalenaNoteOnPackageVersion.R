# specify Seurat 5 for this session 


# View current library paths
.libPaths()

#[1] "/Library/Frameworks/R.framework/Versions/4.2/Resources/library"


#My pathway

#/Users/rhalenathomas/Library/CustomR

# Install Seurat 5 to a separate library location
install.packages("Seurat", lib = "/Users/rhalenathomas/Library/CustomR")

# Load Seurat 5 from the custom library location
library(Seurat, lib.loc = "/Users/rhalenathomas/Library/CustomR")


# Check Seurat version
packageVersion("Seurat")
