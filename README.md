# SingleCell_Workshop 
Brought to you by [The Neuro](https://www.mcgill.ca/neuro/) and [Single Cell Club](https://singlecellclub.openscience.mcgill.ca/)  

We present a complete workflow for single cell transcriptomics analysis - see the [2025 SingleCell Workshop Agenda](https://mcgill-my.sharepoint.com/personal/moein_yaqubi_mail_mcgill_ca/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Fmoein%5Fyaqubi%5Fmail%5Fmcgill%5Fca%2FDocuments%2FSingleCellWorkshop%20agendaJune2025%5Ffinal%5Fagenda%5FApril%5F9th%2Epdf&parent=%2Fpersonal%2Fmoein%5Fyaqubi%5Fmail%5Fmcgill%5Fca%2FDocuments&ga=1)

The workshop is presented in R notebooks which will be run in Magic Castle with  the Digital Computing Alliance of Canada

# Prepare to run analysis
- Everything will be set up for you during the workshop and we are not supporting other methods during the workshop.
- If you want to run the workbooks locally you will need to install all the required R libraries 
- See the installation script for local usage (your laptop or workstation computer) 
- If you are comfortable in base R or HPC you can also convert the workbooks and run scripts in debug/ mode make update your R with HPC_DRAC_installation.md instructions

# Single Cell Sequencing Analysis Workflow

- Section 1: Understanding data structures for scRNAseq Data
  1. Looking at the CellRanger and ParseBio Split-Pipe sequencing alignment and expression matrix outputs.
  2. Reading files into R.  
  3. Data structures and exploration in R.  

- Section 2: 
1. Loading data and generating a Seurat object.
2. scRNAseq data visualization, format and exploration in Seurat.
3. Quality Control measures
4. Filtering data

- Section 3:
1. Combining samples: merging and integration
2. Normalization
3. Feature selection
4. Dimensionality reduction: PCA, UMAP
5. Clustering: nearest neighbour network, Louvain network detection  

- Section 4: Cluster Annotation

Find cluster markers and look at reference cell type libraries
Look at expression of known cell type markers
Automated annotation

- Extra workbook
  Code for topics covered but not run in the workshop
  - Ambient RNA detection and adjustment
  - Doublet identification and removal
  - Data imputation to account for low read gene sequence drop out
  
