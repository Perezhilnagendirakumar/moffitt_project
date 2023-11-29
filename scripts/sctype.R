# cell type labelling using sc-type - automated annotation tool 

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = HNSCC[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(HNSCC@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(HNSCC@meta.data[HNSCC@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(HNSCC@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

HNSCC@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  HNSCC@meta.data$customclassif[HNSCC@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
# label the cell type for clusters 
DimPlot(HNSCC, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif',raster=FALSE)  
### Assigning cell type identity to clusters


### to reset new.cluster.ids
Idents(HNSCC)<-"customclassif"

# Rename classes based on comparision xcell and also https://www.nature.com/articles/s41467-023-37379-y
HNSCC <- RenameIdents(
  object = HNSCC,
  `Pre-B cells` = "B-cells",
  `Plasma B cells` = "Plasma cells",
  `Basophils` = "Mast cells ",
  `Myeloid Dendritic cells` = "Dendritic cells",
  `Naive CD4+ T cells`="CD4+ T cells",
  `Effector CD4+ T cells`="Effector cells",
  `Plasmacytoid Dendritic cells`="pDC",
  `Dendritic cells`="cDC",
  `Natural killer  cells`='NK cells ',
  `CD8+ NKT-like cells`="CD8+cells",
  `Classical Monocytes`="Monocytes"
)



### assign sample name to identity

library(cowplot)

p1 <- DimPlot(object = HNSCC, reduction = 'umap', label = TRUE, pt.size = 0.25, group.by = "orig.ident", repel = TRUE,raster = FALSE)
p2 <- DimPlot(object = HNSCC, reduction = 'umap', label = TRUE, pt.size = 0.5,raster = FALSE)
plot_grid(p1,p2)

save(HNSCC,file="HNSCC_seuratobj.RData")
# How many cells are in each cluster
HNSCC$CellType <- Idents(HNSCC)
table(Idents(HNSCC))
HNSCC.absolute.cells <- table(HNSCC$orig.ident)
HNSCC.absolute.cells


HNSCC.cluster.sample <- table(Idents(HNSCC), HNSCC$orig.ident)
HNSCC.cluster.sample


# What proportion of cells are in each cluster?
prop.table(table(Idents(HNSCC)))

### Subset cell types
Idents(HNSCC) <- "customclassif"
