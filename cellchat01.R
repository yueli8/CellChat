library(CellChat)
library(tidyverse)
library(Seurat)
library(ggalluvial)
library(anndata)
library(Seurat)
library(CellChat)
library(patchwork)
setwd("~/cellchat")
options(stringsAsFactors = FALSE)
load("data_humanSkin_CellChat.rda")
data.input = data_humanSkin$data
meta = data_humanSkin$meta
cell.use = rownames(meta)[meta$condition == "LS"] # extract the c
data.input = data.input[, cell.use]
meta = meta[cell.use, ]
#(D) Starting from an Anndata object
# read the data into R using anndata R package
#ad <- read_h5ad("scanpy_object.h5ad")# access count data matrix
#counts <- t(as.matrix(ad$X))
# normalize the count data if the normalized data is not available in the .h5ad file
#library.size <- Matrix::colSums(counts)
#data.input <- as(log1p(Matrix::t(Matrix::t(counts)/library.size)* 10000), "dgCMatrix")
# access meta data
#meta <- ad$obs
#meta$labels <- meta[["clusters"]]
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
# [1] "Create a CellChat object from a data matrix"
# Set cell identities for the new CellChat object 
# The cell groups used for CellChat analysis are  APOE+ FIB FBN1+ FIB COL11A1+ FIB Inflam. FIB cDC1 cDC2 LC Inflam. DC TC Inflam. TC CD40LG+ TC NKT 
#cellChat <- createCellChat(object = seurat.obj, group.by = "ident", assay = "RNA")
#Error: object 'seurat.obj' not found
CellChatDB <- CellChatDB.human
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")# use Secreted Signaling
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#Inference of cell-cell communication networ
cellchat <- computeCommunProb(cellchat, type = "triMean")
#triMean is used for calculating the average gene expression per cell group. 
#[1] ">>> Run CellChat on sc/snRNA-seq data <<< [2023-12-15 12:49:11.48819]"
#|=============================================================================================================| 100%
#[1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2023-12-15 12:51:34.162272]"
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, file = "cellchat_humanSkin_LS.rds")

#Visualization of cell-cell communication network 
#(A) Circle plot
pathways.show.all <- cellchat@netP$pathways
# select one pathway
pathways.show <- c("CXCL")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout =
                           "circle", color.use = NULL, sources.use = NULL, targets.use = NULL, idents.use = NULL)
#(B) Hierarchy plot
vertex.receiver = seq(1,4)
netVisual_aggregate(cellchat, signaling = pathways.show, layout =
                      "hierarchy", vertex.receiver = vertex.receiver)
#(C) Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout =
                      "chord")
par(mfrow=c(1,1))
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) #grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group =
                       group.cellType, title.name = paste0(pathways.show, " signaling network"))
#(D) Heatmap plot
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
netAnalysis_contribution(cellchat, signaling = pathways.show)
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,]
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

#(A) Bubble plot
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), signaling = c("CCL","CXCL"), remove.isolate = FALSE)
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","CXCL","FGF"))
netVisual_bubble(cellchat, sources.use = c(3,4), targets.use = c(5:8), pairLR.use = pairLR.use, remove.isolate = TRUE)
netVisual_bubble(cellchat, targets.use = c("LC","Inflam. DC","cDC2","CD40LG+ TC"), pairLR.use = pairLR.use, 
                 remove.isolate = TRUE, sort.by.target = T)
netVisual_bubble(cellchat, sources.use = c("FBN1+ FIB","APOE+ FIB","Inflam. FIB"), 
                 pairLR.use = pairLR.use, remove.isolate = TRUE,sort.by.source = T)
netVisual_bubble(cellchat, sources.use = c("FBN1+ FIB","APOE+ FIB","Inflam. FIB"), targets.use = c("LC","Inflam. DC","cDC2","CD40LG+ TC"), 
                 pairLR.use = pairLR.use, remove.isolate = TRUE, sort.by.source = T, sort.by.target = T, sort.by.source.priority = FALSE)
#(B) Chord diagram
netVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(5:11), lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = 8, legend.pos.x = 15)
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), signaling = c("CCL","CXCL"),legend.pos.x = 8)
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), slot.name = "netP", legend.pos.x = 10)
plotGeneExpression(cellchat, signaling = "CXCL", enriched.only =TRUE)

#(A)Compute and visualize the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
#(B) Visualize dominant senders (sources) and receivers (targets) in a 2D space
netAnalysis_signalingRole_scatter(cellchat, signaling = NULL)
#(C) Identify signals contributing the most to outgoing or incoming signaling of certain cell groups
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

#(A) Identify and visualize outgoing communication pattern of secreting cells
#selectK(cellchat, pattern ="outgoing")
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern ="outgoing", k = nPatterns)
netAnalysis_river(cellchat, pattern ="outgoing")
netAnalysis_dot(cellchat, pattern ="outgoing")
#(B) Identify and visualize incoming communication pattern of target cells
#selectK(cellchat, pattern = "incoming")
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
netAnalysis_river(cellchat, pattern = "incoming")
netAnalysis_dot(cellchat, pattern = "incoming")

#(A) Functional similarity
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")#install through "pip install umap-learn"
#或者下面這條命令
#cellchat <- netEmbedding(cellchat, umap.method='uwot',type = "structural")
cellchat <- netClustering(cellchat, type = "functional",do.parallel=FALSE)#加上"do.parallel=FALSE"
netVisual_embedding(cellchat, type = "functional", label.size =3.5)
netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)
#(B) Structure similarity
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")#或者下面這條命令
#cellchat <- netEmbedding(cellchat, umap.method='uwot',type = "structural")
cellchat <- netClustering(cellchat, type = "structural",do.parallel=FALSE)#加上"do.parallel=FALSE"
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)
saveRDS(cellchat, file = "cellchat_humanSkin_LS.rds")

#Part 2. Comparative analysis of cell-cell communication from pairs of scRNA-seq datasets
cellchat.LS <- readRDS("cellchat_humanSkin_LS.rds")
cellchat.NL <- readRDS("cellchat_humanSkin_NL.rds")
cellchat.NL <- updateCellChat(cellchat.NL)
cellchat.LS <- updateCellChat(cellchat.LS)
object.list <- list(NL = cellchat.NL, LS = cellchat.LS)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
save(object.list, file = "cellchat_object.list_humanSkin_NL_LS.RData")
save(cellchat, file = "cellchat_merged_humanSkin_NL_LS.RData")
gg1 <- compareInteractions(cellchat, show.legend = F, group = c (1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c (1,2), measure = "weight")
gg1 + gg2

#(A) Circle plot showing differential number of interactions or interaction strength among different cell populations across two datasets
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

#(B) Heatmap showing differential number of interactions or interaction strength among different cell populations across two datasets
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2

#(C) Circle plot showing the number of interactions or interaction strength among different cell populations across multiple datasets
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T,
                   label.edge= F, edge.weight.max = weight.max[2], edge.width.max =
                     12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

#(D) Circle plot showing the differential number of interactions or interaction strength among coarse cell types
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4))
group.cellType <- factor(group.cellType, levels = c("FIB", "DC","TC"))
object.list <- lapply(object.list, function(x) {mergeInteractions (x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T,
                   label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", names(object.list)[i]))}
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)

#(A) Identify cell populations with significant changes in sending or receiving signals
num.link <- sapply(object.list, function(x) {rowSums(x@net$count)
  + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]],
                                               title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)

#(B) Identify the signaling changes of specific cell populations
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use  = "Inflam. DC", signaling.exclude = "MIF")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "cDC1", signaling.exclude = c("MIF"))
patchwork::wrap_plots(plots = list(gg1,gg2))



source("functional.R")
cellchat <- computeNetSimilarityPairwise (cellchat, type = "functional")



cellchat <- netEmbedding(cellchat,type = "functional")
cellchat <- netClustering(cellchat, type = "functional")

