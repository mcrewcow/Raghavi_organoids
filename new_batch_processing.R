#load the objects
ProcessSeu <- function(Seurat){
Seurat <- NormalizeData(Seurat)
Seurat <- FindVariableFeatures(Seurat, selection.method = "vst", nfeatures = 2000)
Seurat <- ScaleData(Seurat, verbose = T, vars.to.regress = c('nCount_RNA', 'nFeature_RNA', 'percent.mt',"percent.rb","S.Score","G2M.Score"))
Seurat <- RunPCA(Seurat)
Seurat <- FindNeighbors(Seurat, dims = 1:20)
Seurat <- FindClusters(Seurat, resolution = 1.5)
Seurat <- RunUMAP(Seurat, dims = 1:20)
return (Seurat)
}

raghavi2 <- ProcessSeu(raghavi2)
raghavi_1 <- subset(raghavi2, subset = orig.ident == 'U24CPI1')
raghavi_2 <- subset(raghavi2, subset = orig.ident == 'U24CPI2')
raghavi_3 <- subset(raghavi2, subset = orig.ident == 'U24CUI1')
raghavi_4 <- subset(raghavi2, subset = orig.ident == 'U24CUI2')

integration_list <- list(raghavi_1, raghavi_2, raghavi_3,raghavi_4)
features <- SelectIntegrationFeatures(object.list = integration_list)
options(future.globals.maxSize = 64000 * 1024^2)
data.anchors <- FindIntegrationAnchors(object.list = integration_list, anchor.features = features)
data.combined <- IntegrateData(anchorset = data.anchors)
data.combined <- ProcessInt(data.combined)

data.combined <- RenameIdents(data.combined, 
                             '0' = 'Rod',
                             '1' = 'Muller glia',
                             '2' = 'Cone',
                             '3' = 'Cone',
                             '4' = 'Progenitor',
                             '5' = 'Glia',
                             '6' = 'Progenitor',
                             '7' = 'Amacrine',
                             '8' = 'RGC',
                             '9' = 'Glia',
                             '10' = 'Horizontal',
                             '11' = 'Bipolar',
                             '12' = 'Glia',
                             '13' = 'RGC',
                             '14' = 'RGC',
                             '15' = 'RPE')

DimPlot(data.combined)
data.combined$EK_anno <- data.combined@active.ident
SaveH5Seurat(data.combined, 'C://Bioinf/raghavi/human_combined_EKanno.h5Seurat')

library(SeuratDisk)
Convert('C://Bioinf/raghavi/human_combined_EKanno.h5Seurat', dest = 'h5ad')


#further analysis is in Python
