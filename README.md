# "Deciphering Cancer Heterogeneity: Computational Biology Insights into Rare Cell Identification in Head and Neck Squamous Cell Carcinoma (HNSCC) "

The introduction of single-cell RNA sequencing (scRNA-seq) technology has transformed biomedical research, allowing the detailed study of individual cells. The increasing data throughput and enhanced efficiency have resulted in scRNA-seq datasets containing transcription profiles of over a million cells, accompanied by a substantial increase in available analysis tools. Despite these advancements, fully dissecting cellular heterogeneity within a large cell population remains a computational challenge. The primary goal of the project was to employ popular rare cell identification tools like giniclust2, gapclust, and fire to identify rare cells within Head and Neck Squamous Cell Carcinoma (HNSCC). This task is crucial for unraveling the intricacies within cancer populations. In addition,benchmarked these tools by fine-tuning them to ensure their optimal performance, enhancing their accuracy and effectiveness in deciphering the heterogeneity within HNSCC.

# Highlights of the project 

## Head and Neck cancer 
HNSCC occurs with an annual incidence of nearly 600,000 cases globally . HNSCC arises through either genetic alterations driven by exposure to carcinogens or through malignant transformation following high-risk HPV infection. In this study the transcriptional profile of 107738 single cells from HPV - and HPV+ HNSCC were used to decipher the heterogeneity within cancer populations.

### Step 1 :ScRNA analysis using Seurat  
The processed gene barcode matrices were collected from the GEO database (GSE139324) for HPV+ and HPV- PBMC and TIL https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139324. 





## Curation of Marker database from clustermole package
Identification of a marker database that assigns high priority to lowly expressed markers was essential first step in the project . The Clustermole package was instrumental in this process, granting  access to marker information  from a comprehensive array of eight distinct marker databases, which includes ARCHS4, CellMarker, MSigDB, PanglaoDB, SaVanT, TISSUES, and xCell. Each of these databases underwent meticulous curation and organization to prioritize markers characterized by lower expression levels by  harnessing  insights from Bulk RNA-Seq data sourced from the recount-GTEX project.

###  Marker information for each marker database in clustermole 
<img width="278" alt="image" src="https://github.com/Perezhilnagendirakumar/rare_cell_exploration/assets/97453603/e19c647b-34d2-4586-9b2f-8964f325ecc2">

###  Summary table generated based on bulk Rna seq data 
<img width="315" alt="image" src="https://github.com/Perezhilnagendirakumar/rare_cell_exploration/assets/97453603/aceef02f-428c-4c82-879c-c367d40786dc">

###  Visualization of each marker database based on Marker, Non Marker, House keeping , House keeping + Marker information 

<img width="865" alt="image" src="https://github.com/Perezhilnagendirakumar/rare_cell_exploration/assets/97453603/c4e98efd-5c31-4c55-b7fb-f4b5281839ed">
<img width="906" alt="image" src="https://github.com/Perezhilnagendirakumar/rare_cell_exploration/assets/97453603/e0585234-b7dc-46a3-b3f6-98c27b08b8d9">











