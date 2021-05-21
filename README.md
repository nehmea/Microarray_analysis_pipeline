# Microarray_analysis_pipeline
An R code of an automated pipeline that download GEO dataset (raw data), then perform QC, filter low quality samples, then perform limma, and enrichment analysis on the high quality samples.

The following parameters should be specified by the user before running the pipeline (please see the first part of the code):

    *download_series_matrix = "TRUE" #should you start by downloading the series matrix? you may otherwise use your pre-prepared pdata file
    *dataset = "GSE70515" #name of desired GEO dataset to be analyzed
    *primary_variable = "disease" #the name of the variable. you may get this info from any GSM webpage related to the dataset
    *secondary_variable = NA #either specify using a character string of length 1, or NA if not available
    *tertiary_variable = NA #either specify using a character string of length 1, or NA if not available
    *write_downloaded_dataset = "TRUE" #would yu like to save the downloaded series matrix to a local file?
    *pheno_file = "pGSE157784.txt" #if download_series_matrix = "FALSE", what is the name of the phenodata file?
    *control_condition = 'Unstimulated' #the name of the control condition for your primary variable? you may get this info from a control GSM webpage related to the dataset
    *annotdb = 'hgu133plus2.db' #please specify using database name from bioconductor  

This is still a test version. Please report any error or warning message at ali.hassan.nehme@gmail.com
