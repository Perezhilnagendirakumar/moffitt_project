library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
##create separate seurat object for hpv+ and hpv - HNSCC
hpv_pos<- CreateSeuratObject(counts = hpvpositiveHNSCC_matrix, project = "hpv+HNSCC")
hpv_pos@meta.data$orig.ident <- factor(rep("hpv+HNSCC", nrow(hpv_pos@meta.data)))
Idents(hpv_pos) <-"hpv+HNSCC"
# load hpv - combined matrix  and creating seurat obj 
hpv_neg<- CreateSeuratObject(counts = combined_matrix, project = "hpv-HNSCC")
# merge hpv + and hpv - HNSCC together 
HNSCC <- merge(hpv_pos, y = hpv_neg, add.cell.ids = c("hpv+", "hpv-"), project = "HNSCC")
#identify the mitochondrial genes 
HNSCC[["percent.mt"]] <- PercentageFeatureSet(HNSCC, pattern = "^MT-")

## Visualize QC metrics as a violin plot
VlnPlot(HNSCC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1) + NoLegend()

##Feature Scatter is typically used to visualize feature-feature relationships, but can be used for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(HNSCC, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(HNSCC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

HNSCC <- subset(HNSCC, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 25)

### Normalize data
HNSCC<- NormalizeData(HNSCC, normalization.method = "LogNormalize", scale.factor = 10000)

### Identification of highly variable features (feature selection)
HNSCC <- FindVariableFeatures(HNSCC, selection.method = "vst", nfeatures = 2000)
length(VariableFeatures(HNSCC))

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(HNSCC), 10)
top10
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(HNSCC)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

###Scaling the data
HNSCC <- ScaleData(HNSCC, features = VariableFeatures(object = HNSCC))
### perform PCA on the scaled data
HNSCC <- RunPCA(HNSCC, features = VariableFeatures(object = HNSCC))


### Seurat provides several useful ways of visualizing both cells and features that define the PCA, including VizDimReduction, DimPlot, and DimHeatmap
# Examine and visualize PCA results a few different ways
print(HNSCC[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(HNSCC, dims = 1:10, reduction = "pca")
DimPlot(HNSCC, reduction = "pca")
DimHeatmap(HNSCC, dims = 1, cells = 1000, balanced = TRUE)
DimHeatmap(HNSCC, dims = 18:23, cells = 1000, balanced = TRUE)

## we can observe an 'elbow' around PC9-10, suggesting that the majority of true signal is captured in the first 10 PCs.
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
ElbowPlot(HNSCC)
ElbowPlot(object = HNSCC, ndims = 25)

# Determine percent of variation associated with each PC
pct <- HNSCC[["pca"]]@stdev / sum(HNSCC[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
co2
# Minimum of the two calculation
pcs <- min(co1, co2)
pcs
# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))

# Elbow plot to visualize 
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()

#cluster the cells
HNSCC <- FindNeighbors(HNSCC, dims = 1:17)
HNSCC<- FindClusters(HNSCC, resolution = 0.6)
# Look at cluster IDs of the first 5 cells
head(Idents(HNSCC), 5)
HNSCC <- RunUMAP(HNSCC, dims = 1:17) # run UMAP

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
plot1<-DimPlot(HNSCC, label = TRUE, reduction = "umap", raster=FALSE)
# adding metadata information 
metadata <- data.frame(
  SampleID = colnames(HNSCC),
  Strings<-strings,
  HPV<-HNSCC@meta.data[["orig.ident"]],
  Tissue = ifelse(strings %in% test, "Peripheral Blood", "Tumor Tissue"),
  stringsAsFactors = FALSE
)
HNSCC@meta.data$Tissue<- metadata$Tissue
HNSCC@meta.data$HPV<-metadata$HPV
# highlight based on tissue 
plot2<-DimPlot(HNSCC, reduction = "umap",group.by = "Tissue",raster=FALSE )
# highlight based on HPV status 
plot3<-DimPlot(HNSCC, reduction = "umap",group.by = "HPV",cols=c("brown", "green"),raster=FALSE )


# adding the patient information to metadata 
expression_matrix <- lapply(names(expression_matrices_GSE139324), function(name) {
  # Remove the GSM pattern from the sample name
  cleaned_name <- gsub("GSM\\d+_", "", name)
  ExpressionMatrix = expression_matrices_GSE139324[[name]]
  columns<-colnames(ExpressionMatrix)
  # Create a new list with cleaned name and corresponding expression matrix
  list(
    SampleID = cleaned_name,
    cell_names= columns
  )
})

metadata$Patient <- NA

for (i in seq_along(expression_matrix)) {
  sample_id <- expression_matrix[[i]]$SampleID
  cell_names <- expression_matrix[[i]]$cell_names
  
  # Check if the cell names are present in the metadata dataframe
  matching_indices <- match(cell_names, metadata$Strings)
  
  # Add sample ID information to the metadata dataframe in a new column 'Patient'
  metadata$Patient[matching_indices] <- sample_id
}
HNSCC@meta.data$Patient<-metadata$Patient
save(metadata,file = "metadata.RData")

# highlight based on patient information 
plot4<-DimPlot(HNSCC, reduction = "umap",group.by = "Patient",raster=FALSE )
save(HNSCC,file="HNSCC_seuratobj.RData")

plot_grid(plot1,plot2,plot3,plot4, labels = c('A', 'B','C','D'),ncol = 2,nrow=2)

cols=c("brown", "green","#FF68A1", "#00B8E7")
test<-data.frame(c(table(HNSCC_seuratobj$orig.ident),table(HNSCC_seuratobj@meta.data[["Tissue"]])))
barplot(test$freq, names.arg = rownames(test), col = cols,
        main = "Number of Cells ",width = 0.2)

### Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster 1
##cluster1.markers <- FindMarkers(CTRLPD1, ident.1 = 1, min.pct = 0.25)
cluster1.markers <- FindMarkers(HNSCC, ident.1 = 4, min.pct = 0.25)
head(cluster1.markers, n = 5)
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(HNSCC, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
HNSCC <- FindAllMarkers(HNSCC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
HNSCC.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.table(HNSCC, 'markers_allclusters.csv', sep=',', row.names = TRUE)


### Seurat has several tests for differential expression which can be set with the test.use parameter (see our DE vignette for details).
### For example, the ROC test returns the 'classification power' for any individual marker (ranging from 0 - random, to 1 - perfect).
cluster1.markers <- FindMarkers(HNSCC, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

### violin plots
VlnPlot(HNSCC, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
                            "CD8A"), pt.size = 0.1,raster=FALSE)
VlnPlot(HNSCC, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
                                   "CD8A", "NKG7", "PF4"), pt.size = 0.1, slot = "counts", log = TRUE)

# marker expression plots. Export file to pdf and set print size by 10x10 inches
features<-c( 
  "CD79A", "IGHG1", "CD14", "FCER1A", "FCGR3A", "LYZ", "CD68", "LILRA4",
  "PECAM1", "MMP2", "COL1A2", "ACTA2", "KRT17", "KRT7",
  "XCL2", "NCAM1", "CD3D", "CD3E", "CD2", "CD8A", "CD4", "FOXP3",
  "CXCL13", "LAG3", "CTLA4", "PDCD1", "TIGIT", "HAVCR2", "TCF7",
  "TPSAB1", "CD68", "CD163", "FCER1G", "IL1B", "LAMP3")
perform_feature_plot <- function(seurat_object, features, output_filename) {
  no_genes<-length(features)
  plot_list <- list()
  pdf(output_filename)
  for (i in 1:no_genes) {
    gene<-features[i]
    plot <- FeaturePlot(seurat_object, features = gene, label = TRUE,raster=FALSE)
    
    # Add the plot to the list
    plot_list[[i]] <- plot
    
    # Save the plot to the PDF file
    print(plot)
  }
  
  # Close the PDF device
  dev.off()
  
  # Return the list of plots
  return(plot_list)
}








