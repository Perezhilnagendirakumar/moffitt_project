library(GEOquery)
gse<-getGEO("GSE139324",GSEMatrix = TRUE)
metadata_GSE139324<-pData(phenoData(gse[[1]]))
library(scCustomize)
data_dir <-"GSE139324_RAW"
#expression_matrices_GSE139324 <- Read10X_GEO(data_dir = data_dir)

metadata_GSE139324
load(file="~/Desktop/scrna /HNSCC/expression_matrices_GSE139324.RData")
View(expression_matrices_GSE139324)

# create HPV + and HPV - matrix by combining all cells 
hpvnegative_list<-expression_matrices_GSE139324[1:36]
hpvpositive_list<-expression_matrices_GSE139324[37:52]

# create hpv+ and Hpv - matrix separately 
hpvneg_matrix<-hpvnegative_list[[1]]
hpvpos_matrix<-hpvpositive_list[[1]]
for (i in 2:length(hpvnegative_list)){
  m2<-hpvnegative_list[[i]]
  hpvneg_matrix<-cbind2(hpvneg_matrix,m2)
}
# repeat the same function to create the hpvpos_matrix 
# file is stored in hpv-PBMC+TIL_combinedmatrix.RData and hpvpositiveHNSCC_matrix.RData
pbmchpvpositive_matrix <-hpvpositive_list[[1]]
TILpositive_matrix <- NULL

# PBMC and TIL matrices for HPV + and HPV - HNSCC
for (i in 2:length(hpvpositive_list)) {
  if (grepl("PBMC", names(hpvpositive_list)[i])) {
    m2 <-hpvpositive_list[[i]]
    pbmchpvpositive_matrix<- cbind2(pbmchpvpositive_matrix, m2)
  } else if (grepl("TIL", names(hpvnegative_list)[i])) {
    # Combine TIL matrices
    m2 <- hpvpositive_list[[i]]
    if (is.null(TILpositive_matrix)) {
      TILpositive_matrix <- m2
    } else {
      TILpositive_matrix <- cbind2(TILpositive_matrix, m2)
    }
  }
}
# repeat the same for pbmchpvnegative_matrix and TILnegative_matrix
# these files are stored in hnscc_hpv-_separatematrix.RData and hnscc_hpv+_separatematrix.RData
