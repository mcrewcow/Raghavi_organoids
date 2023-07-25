#Read the counts
message("Reading counts...")
x <- read.csv("C://Users/Emil/10X/mgtotal_Iraw-counts.csv",header=TRUE)
rownames(x) <- x[,1]
x[,1] <- NULL
print(dim(x))
print(x[1:5,1:5])
#Read the metadata
message("Reading metadata...")
m <- read.csv("C://Users/Emil/10X/mgtotal_Imetadata.csv",header=TRUE)
rownames(m) <- m[,1]
colnames(m)[1] <- "sample"
print(dim(m))
print(head(m))
#Create and save Seurat object, could be changed to .h5Seurat
mgtotal_I <-    CreateSeuratObject(counts=t(x),meta.data=m,project="seurat",min.cells=0,min.features=0)
head(mgtotal_I)
cellchat <- createCellChat(object = mgtotal_I, group.by = "EK_anno2")
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.05, raw.use = FALSE, population.size = FALSE)
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", height = 5, width = 4)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", height = 5, width = 4)
ht1 + ht2
