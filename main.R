U24CPI3.data <- Read10X(data.dir = "C:/Bioinf/raghavi/filtered_feature_bc_matrix_CUI1/")
U24CPI3 <- CreateSeuratObject(counts = U24CPI3.data, project = "CPI1", min.cells = 3, min.features = 200 )
U24CPI3

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

ProcessSeu <- function(Seurat){
  Seurat <- NormalizeData(Seurat)
  Seurat <- FindVariableFeatures(Seurat, selection.method = "vst", nfeatures = 2000)
  Seurat <- ScaleData(Seurat, verbose = T) #could be replaced with SCTransform
  Seurat <- RunPCA(Seurat, npcs = 150)
  Seurat <- FindNeighbors(Seurat, dims = 1:150)
  Seurat <- FindClusters(Seurat, resolution = 1)
  Seurat <- RunUMAP(Seurat, dims = 1:150)
  DimPlot(object = Seurat, reduction = "umap")
  return (Seurat)
}

RDoublet <- function(tmp){
  sweep.res.list <- paramSweep_v3(tmp, PCs = 1:30, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pKopt <- as.numeric(as.character(bcmvn$pK[bcmvn$BCmetric == max(bcmvn$BCmetric)]))
  pKopt <- pKopt[order(pKopt, decreasing = TRUE) ]
  pKopt <- pKopt[1]
  homotypic.prop <- modelHomotypic(tmp$seurat_clusters) 
  nExp_poi <- round(0.1*length(colnames(tmp)))  ## Assuming 10% doublet formation rate 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  tmp <- doubletFinder_v3(tmp, PCs = 1:30, pN = 0.25, pK = pKopt, nExp = nExp_poi, reuse.pANN = FALSE)
  tmp <- doubletFinder_v3(tmp, PCs = 1:30, pN = 0.25, pK = pKopt, nExp = nExp_poi.adj, reuse.pANN = paste("pANN_0.25",pKopt,nExp_poi, sep="_"))
  return (tmp) 
}


U24CPI3[["percent.rb"]] <- PercentageFeatureSet(U24CPI3, pattern = "^RPS|^RPL|^MRPS|^MRPL", assay = 'RNA')
U24CPI3[["percent.mt"]] <- PercentageFeatureSet(U24CPI3, pattern = "^MT-")
U24CPI3 <- CellCycleScoring(U24CPI3, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

VlnPlot(U24CPI3, features = c("nFeature_RNA", "nCount_RNA","percent.mt", 'percent.rb'), ncol = 4)

U24CPI3 <- subset(U24CPI3, subset = nCount_RNA > 300 & nCount_RNA < 30000 & nFeature_RNA > 500 & nFeature_RNA < 5800 & percent.mt < 20 & percent.rb < 40)
U24CPI3 <- ScaleData(U24CPI3, verbose = T, vars.to.regress = c('percent.mt', "percent.rb","S.Score","G2M.Score")) 

U24CPI3 <- ProcessSeu(U24CPI3)

#Basic visualisation
FeaturePlot(U24CPI3, features = c('NRL'))
DimPlot(U24CPI3, reduction = 'umap', label = TRUE, repel = TRUE)

U24CPI3 <- RDoublet(U24CPI3)
U24CPI3 <- subset(U24CPI3, cells = colnames(U24CPI3 )[which(U24CPI3 [[]][12] == 'Singlet')])
U24CPI3 <- subset(U24CPI3 , cells = colnames(U24CPI3 )[which(U24CPI3 [[]][13] == 'Singlet')])

U24CPI3 <- ProcessSeu(U24CPI3)

U24CPI3$condition <- 'Uninjected'
U24CPI3$batch <- 'Batch 1'

SaveH5Seurat(U24CPI3, 'C://Bioinf/raghavi/CUI1.h5Seurat')

CUI1 <- U24CPI3

integration_list <- list(CUI1, CUI2, CUI3, CPI1, CPI2, CPI3)
features <- SelectIntegrationFeatures(object.list = integration_list)
data.anchors <- FindIntegrationAnchors(object.list = integration_list, anchor.features = features)
raghavi_integrated <- IntegrateData(anchorset = data.anchors)

ProcessInt <- function(data.integrated){
  
  data.integrated <- ScaleData(data.integrated, verbose = T)
  data.integrated <- RunPCA(data.integrated, npcs = 120, verbose = T)
  data.integrated <- FindNeighbors(data.integrated, dims = 1:200)
  data.integrated <- FindClusters(data.integrated, resolution = 1.5)
  data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:200)
}

raghavi_integrated <- ProcessInt(raghavi_integrated)

SaveH5Seurat(SB_zebrafish_integrated, 'C://Bioinf/Blackshaw_zebrafish_Hoang_science/zebrafish_Hoang_2023_EK_anno_NMDA_LD.h5Seurat')

data.combined.markers <- FindAllMarkers(SB_zebrafish_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

data.combined.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
data.combined.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(SB_zebrafish_integrated, features = top10$gene) + NoLegend()
