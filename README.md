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

- Section 1: Understanding scRNAseq Data
1. Loading data
2. Quality Control
   2.1 Cell filtering
   2.2 Doublet assessment
3. Normalization
4. Feature selection
5. Dimensionality reduction

- Section 2: Datat integration, Clustering and Annotation
6. Merging samples and batch correction
7. Clustering
8. Cluster annotation
  8.1 Find cluster markers and look at reference cell type libraries
  8.2 Look at expression of known cell type markers
  8.3 Automated annotation


- Section 3: Futher Analaysis
9.  Comparing clusters in different datasets
10. Comparing different clustering labels for the same dataset
11. Compare proportions of cell types
12. Differential Gene Expression Analysis
