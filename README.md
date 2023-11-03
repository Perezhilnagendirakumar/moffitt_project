# Optimizing the Impact of Lowly Expressed Markers on Rare Cell Type Identification: Parameter Tuning and Benchmarking of Leading Single-Cell Analysis Tools

The introduction of single-cell RNA sequencing (scRNA-seq) technology has revolutionized biomedical research by enabling the study of individual cells. With increasing data throughput and improved efficiency, scRNA-seq datasets now encompass transcription profiles of over a million cells which also lead to tremendous increase in analysis tools. Nevertheless, the absence of standardization in these tools has created ambiguity over the selection of the most suitable method for analysis. Moreover, relying on default parameters provided by these pipelines poses a significant risk, potentially leading to error-prone cell type classification and annotation.

## project overview 
The  goal of the project is to enhance the accuracy and sensitivity of cell type identification using open access scRNA-seq data by optimizing analysis tools and understanding the impact of lowly expressed markers.

# Highlights of the project 

## Curation of Marker database from clustermole package
Identification of a marker database that assigns high priority to lowly expressed markers was essential first step in the project . The Clustermole package was instrumental in this process, granting  access to marker information  from a comprehensive array of eight distinct marker databases, which includes ARCHS4, CellMarker, MSigDB, PanglaoDB, SaVanT, TISSUES, and xCell. Each of these databases underwent meticulous curation and organization to prioritize markers characterized by lower expression levels by  harnessing  insights from Bulk RNA-Seq data sourced from the recount-GTEX project.

###  Marker information for each marker database in clustermole 
<img width="278" alt="image" src="https://github.com/Perezhilnagendirakumar/rare_cell_exploration/assets/97453603/e19c647b-34d2-4586-9b2f-8964f325ecc2">

###  Summary table generated based on bulk Rna seq data 
<img width="315" alt="image" src="https://github.com/Perezhilnagendirakumar/rare_cell_exploration/assets/97453603/aceef02f-428c-4c82-879c-c367d40786dc">

###  Visualization of each marker database based on Marker, Non Marker, House keeping , House keeping + Marker information 

<img width="865" alt="image" src="https://github.com/Perezhilnagendirakumar/rare_cell_exploration/assets/97453603/c4e98efd-5c31-4c55-b7fb-f4b5281839ed">
<img width="906" alt="image" src="https://github.com/Perezhilnagendirakumar/rare_cell_exploration/assets/97453603/e0585234-b7dc-46a3-b3f6-98c27b08b8d9">











