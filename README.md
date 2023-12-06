# "Deciphering Cancer Heterogeneity: Computational Biology Insights into Rare Cell Identification in Head and Neck Squamous Cell Carcinoma (HNSCC) "

The introduction of single-cell RNA sequencing (scRNA-seq) technology has transformed biomedical research, allowing the detailed study of individual cells. The increasing data throughput and enhanced efficiency have resulted in scRNA-seq datasets containing transcription profiles of over a million cells, accompanied by a substantial increase in available analysis tools. Despite these advancements, fully dissecting cellular heterogeneity within a large cell population remains a computational challenge. The primary goal of the project was to employ popular rare cell identification tools like giniclust2, gapclust, and fire to identify rare cells within Head and Neck Squamous Cell Carcinoma (HNSCC). This task is crucial for unraveling the intricacies within cancer populations. In addition,benchmarked these tools by fine-tuning them to ensure their optimal performance, enhancing their accuracy and effectiveness in deciphering the heterogeneity within HNSCC.

# Highlights of the project 

## Head and Neck cancer 
HNSCC occurs with an annual incidence of nearly 600,000 cases globally . HNSCC arises through either genetic alterations driven by exposure to carcinogens or through malignant transformation following high-risk HPV infection. In this study the transcriptional profile of 107738 single cells from HPV - and HPV+ HNSCC were used to decipher the heterogeneity within cancer populations.

### Step 1 :ScRNA analysis using Seurat  
The processed gene barcode matrices were collected from the GEO database (GSE139324) for HPV+ and HPV- PBMC and TIL https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139324. Based on the metadata information collected from GEO separate matrices were created for HPV + HNSCC (33694 genes *35k cells) and HPV - HNSCC each comprising (33694 genes * 70k cells). Seurat analysis was performed on these matrices which includes preprocessing steps like QC , Normalization , Variable selection and PCA linear dimension reduction and the two seurat object was merged together which contains almost 36694 genes * 107738 cells and further unsupervised clustering was performed . A total of 23 clusters were obtained and cells were highlighted based on the metadata information like Tissue of origin , patient information, HPV status 

![image](https://github.com/Perezhilnagendirakumar/rare_cell_exploration/assets/97453603/3f558285-7014-49bd-bcee-448339be2540)


### Step 2 : Annotation of cell type using automated and ultrafast cell type clssification tool "SC type"
Initially, Sctype was performed on the entire dataset to identify the major trends, and to identify the major immune lineages. Sctype follows marker database and cell-type identification algorithm for accurate and unsupervised cell-type annotation.Most of the clusters were classified based on the marker information from Sctype platform. 

![image](https://github.com/Perezhilnagendirakumar/rare_cell_exploration/assets/97453603/15ca7fb4-7382-45e4-8fa1-8bd72ebcf681)

### Step 3 : Exploring Rare cell type identification tools to spot the rare clusters 
Rare cells represent minor cell types in an organism. Despite low abundance, rare cell populations play an important role in determining the pathogenesis of cancer, mediating immune responses, angiogenesis in cancer and other diseases. Spotting the rare cells using unsupervised clustering approach is tedious as they include one or more parameters that can be chosen by the user to determine the resolution of the clustering. The choice of parameter often has a large effect on the outcome as it leads to underclustering or overclustering . FiRE, Giniclust2 ,Gapclust uses different strategies for optimal identification of genes associated with rare cell types compared to common ones in a complex tissue. Employed these tools on our HNSCC to spot the rare cell population. 

![image](https://github.com/Perezhilnagendirakumar/rare_cell_exploration/assets/97453603/a53c115d-1b27-4995-8d04-4f0a8b14f8a9)

### Step 4 : Subsetting the Rare Cells and  Reclustering using seurat 
Fire, Gapclust, and Giniclust identified 2500, 95, and 500 rare cells, respectively. These rare cells, pinpointed by each method, were visually highlighted in UMAP plots. Subsequently, the rare cell clusters were isolated and subjected to a re-clustering analysis to examine their distinct transcriptional states.

![image](https://github.com/Perezhilnagendirakumar/rare_cell_exploration/assets/97453603/eabd68ba-7b11-42e7-93f7-e3aa548b7ff5)
![image](https://github.com/Perezhilnagendirakumar/rare_cell_exploration/assets/97453603/4f1f5fbc-c4f5-4b3b-bcde-9835a5b47be4)

### Step 5 : Curation of Marker database from clustermole package
In our efforts to streamline the assignment of cell type labels to clusters, particularly focusing on the representation of rare cells or subcellular types, we undertook the task of identifying a marker database rich in information about lowly expressed markers. Recognizing that such markers are crucial for characterizing rare cell populations, we utilized the clustermole package, which provided access to a comprehensive array of eight distinct marker databases, including ARCHS4, CellMarker, MSigDB, PanglaoDB, SaVanT, TISSUES, and xCell. Our emphasis in this selection process was on databases that prioritize markers with lower expression levels, a feature essential for accurately identifying and characterizing rare cell types within our study, as informed by insights derived from Bulk RNA-Seq data from the recount-GTEX project.

###  Marker information for each marker database in clustermole 
<img width="278" alt="image" src="https://github.com/Perezhilnagendirakumar/rare_cell_exploration/assets/97453603/e19c647b-34d2-4586-9b2f-8964f325ecc2">

###  Summary table generated based on bulk Rna seq data 
<img width="315" alt="image" src="https://github.com/Perezhilnagendirakumar/rare_cell_exploration/assets/97453603/aceef02f-428c-4c82-879c-c367d40786dc">

###  Visualization of each marker database based on Marker, Non Marker, House keeping , House keeping + Marker information 

<img width="865" alt="image" src="https://github.com/Perezhilnagendirakumar/rare_cell_exploration/assets/97453603/c4e98efd-5c31-4c55-b7fb-f4b5281839ed">
<img width="906" alt="image" src="https://github.com/Perezhilnagendirakumar/rare_cell_exploration/assets/97453603/e0585234-b7dc-46a3-b3f6-98c27b08b8d9">

Based on the results of the density plot and interpretation on table the preferred database is ranked 

1.xCell
2.CellMarker
3.MSigDB
4.SaVanT

### Step 6: Annotating the cell types based on marker database and literature survey 
<img width="243" alt="image" src="https://github.com/Perezhilnagendirakumar/rare_cell_exploration/assets/97453603/1e1f48a6-1c23-49ea-a632-856c3916b096">

### Step 7: Optimizing Rare Cell Tools: Fine-Tuning Default Cutoffs for Accurate Rare Cell Percentage Assessment
GiniClust2 and FiRE employ a threshold-based approach for selecting informative genes contributing to rare cell populations. GiniClust2 utilizes Norm.Gini.cutoff (set to 1) and Gini.pvalue_cutoff (0.0001) to identify high Gini genes, while FiRE adopts an interquartile range (IQR) based approach (75th percentile - 25th percentile) on the rareness score to pinpoint rare cells. Instead of relying on a cutoff-based approach, genes were ordered based on these values from low to high, and the tool's performance in predicting rare cells was assessed within the range of top 100 to top 5000 genes specifically for HNSCC. PBMC 68k served as a control for performance evaluation.

![image](https://github.com/Perezhilnagendirakumar/rare_cell_exploration/assets/97453603/610163a2-f0a4-486a-85d5-81bae392929d)
![image](https://github.com/Perezhilnagendirakumar/rare_cell_exploration/assets/97453603/f29061e5-2073-4172-92ee-1e3f398adefc)

### Step 8 : Exploring Gene-PC Relationships: 
The analysis included monitoring variation at the principal component (PC) level for each gene. This was achieved through correlation analyses between individual PCs and specific gene ranges. This act as a guideline to choose appropriate pcs for downstream analysis . the results are attached in the name of correlation analysis1 to 2000_giniclust and correlation analysis1 to 2000_fire in results folder 






