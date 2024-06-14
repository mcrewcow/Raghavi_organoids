library(Seurat)
library(SeuratDisk)
library(escape)


raghavi <- readRDS('C://Bioinf/raghavi/fetal_annotated.rds')

gene.sets1 <- getGeneSets(library = c('H',"C5"), gene.sets = c('HALLMARK_APOPTOSIS','GOBP_NECROPTOTIC_SIGNALING_PATHWAY','GOBP_NECROPTOTIC_PROCESS',
                                                           'GOBP_DNA_REPAIR', 'GOBP_PHOTORECEPTOR_CELL_DIFFERENTIATION', 'GOBP_PHOTORECEPTOR_CELL_OUTER_SEGMENT_ORGANIZATION',
                                                           'GOBP_NEUROBLAST_PROLIFERATION'), species = 'Homo sapiens')


ES <- enrichIt(obj = raghavi,
               gene.sets = gene.sets1, min.size = 0,
               groups = 1000)

head(ES)

raghavi<- AddMetaData(raghavi, ES)
ES2 <- data.frame(raghavi[[]], Idents(raghavi))
colnames(ES2)[ncol(ES2)] <- "cluster"

ridgeEnrichment(ES2, gene.set = "GOBP_NECROPTOTIC_SIGNALING_PATHWAY", group = 'condition',add.rug = FALSE) + facet_wrap(~final.id, nrow = 12)


metadata <- raghavi@meta.data

# Convert metadata to a data frame
metadata_df <- as.data.frame(metadata)

# Write the metadata to a CSV file
write.csv(metadata_df, file = "C://Bioinf/raghavi/metadata_GSEA.csv", row.names = TRUE)
