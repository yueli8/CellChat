#install.packages('devtools')
#devtools::install_github("jinworks/CellChat")
library(CellChat)
library(tidyverse)
library(ggalluvial)
#install.packages("anndata")
library(anndata)
library(anndata)
library(Seurat)
library(patchwork)
library(uwot)
library(reticulate)
library(glmGamPoi)
library(umap)
library(ComplexHeatmap)
#install.packages('NMF')
library(NMF)
library(dplyr)
library(SeuratData)
library(ggplot2)
library(svglite)
library(wordcloud)
library(wordcloud2)
library(tm)
#devtools::install_github("jokergoo/circlize")
library(circlize)
#devtools::install_github("jokergoo/ComplexHeatmap")

#source('functional.R')#千萬不要加上
setwd("~/cellchat")
rm(list=ls()) #清空所有变量
options(stringsAsFactors = FALSE)#输入数据不自动转换成因子（防止数据格式错误）

#(A) Starting from a count data matrix
load("data_humanSkin_CellChat.rda")
data.input = data_humanSkin$data
meta = data_humanSkin$meta
cell.use = rownames(meta)[meta$condition == "LS"] # extract the c
data.input = data.input[, cell.use]
meta = meta[cell.use, ]

#(B) Starting from a Seurat object
#data.input <- GetAssayData(seurat_object, assay = "RNA", slot = "data") # normalized data matrix
#labels <- Idents(seurat_object)
#meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

#(C) Starting from a SingleCellExperiment object
#data.input <- SingleCellExperiment::logcounts(object) # normalized data matrix
#meta <- as.data.frame(SingleCellExperiment::colData(object)) # extract a dataframe of the cell labels
#meta$labels <- meta[["sce.clusters"]]

#(D) Starting from an Anndata object
# read the data into R using anndata R package
#ad <- read_h5ad("scanpy_object.h5ad")# access count data matrix
#counts <- t(as.matrix(ad$X)) # normalize the count data if the normalized data is not available in the .h5ad file
#library.size <- Matrix::colSums(counts)
#data.input <- as(log1p(Matrix::t(Matrix::t(counts)/library.size)* 10000), "dgCMatrix")
#meta <- ad$obs# access meta data
#meta$labels <- meta[["clusters"]]


#2.Create a CellChat object by following option
#(A) Starting from the digital gene expression matrix and cell label information
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
# [1] "Create a CellChat object from a data matrix"
# Set cell identities for the new CellChat object 
# The cell groups used for CellChat analysis are  APOE+ FIB FBN1+ FIB COL11A1+ FIB Inflam. FIB cDC1 cDC2 LC Inflam. DC TC Inflam. TC CD40LG+ TC NKT 
#(B) Starting from a Seurat object
#cellChat <- createCellChat(object = seurat.obj, group.by = "ident", assay = "RNA")
#(C) Starting from a SingleCellExperiment object
#cellChat <- createCellChat(object = sce.obj, group.by = "sce.clusters")
#(D) Starting from an AnnData object
#sce <- zellkonverter::readH5AD(file = "adata.h5ad") # retrieve all the available assays within sce object assayNames(sce)
# add a new assay entry "logcounts" if not available
#counts <- assay(sce, "X") # make sure this is the original count data matrix
#library.size <- Matrix::colSums(counts)
#logcounts(sce) <- log1p(Matrix::t(Matrix::t(counts)/library.size) * 10000)
# extract a cell meta data
#meta <- as.data.frame(SingleCellExperiment::colData(sce)) #
#cellChat <- createCellChat(object = sce, group.by = "sce.clusters")


CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")# use Secreted Signaling
CellChatDB.use <- CellChatDB # simply use the default CellChatDB to use all CellChatDB for cell-cell communication analysis
# set the used database in the object
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#Inference of cell-cell communication network
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
#(A) Circle plot ## Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# select one pathway
pathways.show <- c("CXCL")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle", color.use = NULL, sources.use = NULL, targets.use = NULL, idents.use = NULL)
#(B) Hierarchy plot
vertex.receiver = seq(1,4)
netVisual_aggregate(cellchat, signaling = pathways.show, layout ="hierarchy", vertex.receiver = vertex.receiver)
#(C) Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout ="chord")
par(mfrow=c(1,1))
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) #grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group =group.cellType, title.name = paste0(pathways.show, " signaling network"))
#(D) Heatmap plot
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

netAnalysis_contribution(cellchat, signaling = pathways.show)
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,]
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

#(A) Bubble plot
# (1) show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
# (2) show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), signaling = c("CCL","CXCL"), remove.isolate = FALSE)
## (3) show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","CXCL","FGF"))
netVisual_bubble(cellchat, sources.use = c(3,4), targets.use = c(5:8), pairLR.use = pairLR.use, remove.isolate = TRUE)
# set the order of interacting cell pairs on x-axis
# (4) Default: first sort cell pairs based on the appearance of sources in levels(object@idents), and then based on the appearance of targets in levels(object@idents)
# (5) sort cell pairs based on the targets.use defined by users
netVisual_bubble(cellchat, targets.use = c("LC","Inflam. DC","cDC2","CD40LG+ TC"), pairLR.use = pairLR.use, 
                 remove.isolate = TRUE, sort.by.target = T)
# (6) sort cell pairs based on the sources.use defined by users
netVisual_bubble(cellchat, sources.use = c("FBN1+ FIB","APOE+ FIB","Inflam. FIB"), 
                 pairLR.use = pairLR.use, remove.isolate = TRUE,sort.by.source = T)
# (7) sort cell pairs based on the sources.use and then targets.use defined by users
netVisual_bubble(cellchat, sources.use = c("FBN1+ FIB","APOE+ FIB","Inflam. FIB"), targets.use = c("LC","Inflam. DC","cDC2","CD40LG+ TC"), pairLR.use = pairLR.use, remove.isolate = TRUE, 
                 sort.by.source = T, sort.by.target = T)
# (8) sort cell pairs based on the targets.use and then sources.use defined by users
netVisual_bubble(cellchat, sources.use = c("FBN1+ FIB","APOE+ FIB","Inflam. FIB"), targets.use = c("LC","Inflam. DC","cDC2","CD40LG+ TC"), 
                 pairLR.use = pairLR.use, remove.isolate = TRUE, sort.by.source = T, sort.by.target = T, sort.by.source.priority = FALSE)

#(B) Chord diagram
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB
netVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(5:11), lab.cex = 0.5,legend.pos.y = 30)
# show all the interactions received by Inflam.DC
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = 8, legend.pos.x = 15)
# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), signaling = c("CCL","CXCL"),legend.pos.x = 8)
# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), slot.name = "netP", legend.pos.x = 10)
plotGeneExpression(cellchat, signaling = "CXCL", enriched.only =TRUE)

#(A)Compute and visualize the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")# the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groupsnetAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
#(B) Visualize dominant senders (sources) and receivers (targets) in a 2D space
netAnalysis_signalingRole_scatter(cellchat, signaling = NULL)
#(C) Identify signals contributing the most to outgoing or incoming signaling of certain cell groups
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

#(A) Identify and visualize outgoing communication pattern of secreting cells
# infer the number of patterns.
selectK(cellchat, pattern ="outgoing")
# Both Cophenetic and Silhouette values begin to drop suddenly when the number of outgoing patterns is 3.
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern ="outgoing", k = nPatterns)
netAnalysis_river(cellchat, pattern ="outgoing")# river plot
netAnalysis_dot(cellchat, pattern ="outgoing")# dot plot

#(B) Identify and visualize incoming communication pattern of target cells
selectK(cellchat, pattern = "incoming")
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
netAnalysis_river(cellchat, pattern = "incoming")# river plot
netAnalysis_dot(cellchat, pattern = "incoming")# dot plot

#(A) Functional similarity
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")#install through "pip install umap-learn"
#或者下面三條命令都可以運行
#cellchat <- netEmbedding(cellchat, umap.method='uwot',type ="functional")
#cellchat <- netEmbedding(cellchat, umap.method='umap-learn',type ="functional")
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
netVisual_diffInteraction(cellchat, weight.scale = T)#install igraph  version1.3.5  #https://blog.csdn.net/m0_55681975/article/details/131965802
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
                     12, title.name = paste0("Number of interactions - ", names(object.list)[i]))}

#(D) Circle plot showing the differential number of interactions or interaction strength among coarse cell types
# Here, CellChat categorize the cell populations into three cell types, and then re-merge the list of CellChat objects.
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

# Similarly, CellChat can also show the differential number of interactions or interaction strength between any two cell types using circle plot.
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)

#(A) Identify cell populations with significant changes in sending or receiving signals
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]],
                                               title = names(object.list)[i], weight.MinMax = weight.MinMax)}
patchwork::wrap_plots(plots = gg)

#(B) Identify the signaling changes of specific cell populations
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use  = "Inflam. DC", signaling.exclude = "MIF")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "cDC1", signaling.exclude = c("MIF"))
patchwork::wrap_plots(plots = list(gg1,gg2))

#remotes::install_github("hafen/rminiconda")
#py <- rminiconda::find_miniconda_python("my_python")
#reticulate::use_python("/home/yueli/.local/share/r-miniconda/envs/r-reticulate/bin/python", required = TRUE)
#reticulate::py_install(packages = 'umap-learn')

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional") 
#cellchat <- netEmbedding(cellchat, umap.method="uwot",type = "functional") 
#cellchat <- netEmbedding(cellchat, umap.method="umap-learn",type = "functional")

cellchat <- netClustering(cellchat, type = "functional",do.parallel=FALSE)#
#cellchat <- netClustering(cellchat, method="umap-learn", type = "functional")#not work
#cellchat <- netClustering(cellchat, method="uwot", type = "functional")#not work
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", nCol = 2)


rankSimilarity(cellchat, type = "functional")
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

#(B) Compare outgoing (or incoming) signaling patterns associated with each cell population
i = 1
# combining all the identified signaling pathways from different datasets
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", 
                                        signaling = pathway.union, title = names(object.list)[i], width = 5, height= 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", 
                                        signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

# Compare the communication probabilities from certain cell groups to other cell groups
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),comparison= c(1, 2), angle.x = 45)
# Identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset.                 
gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in LS",
                        angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in LS",
                        angle.x = 45, remove.isolate = T)
gg1 + gg2

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "LS"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, 
                                       features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "LS",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "NL",ligand.logFC = -0.1, receptor.logFC = -0.1)
# do further deconvolution to obtain the individual signaling genes
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
# Users can also find all the significant outgoing/incoming/both signaling according to the customized DEG features and cell groups of interest
#df <- findEnrichedSignaling(object.list[[2]], features = c("CCL19", "CXCL12"), idents = c("Inflam. FIB", "COL11A1+ FIB"), pattern ="outgoing")  #It is availabe in CellChat version 2.1.1 
#Error in findEnrichedSignaling(object.list[[2]], features = c("CCL19",  : 
#                                                               could not find function "findEnrichedSignaling"pairLR.use.up = net.up[, "interaction_name", drop = F]

#(A) Bubble plot
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2), angle.x = 90, remove.isolate =  T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2), angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 + gg2

#Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 4, targets.use = c(5:11),
                     slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name =
                       paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 4, targets.use = c(5:11),
                     slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name
                     = paste0("Down-regulated signaling in ", names(object.list)[2]))

#Wordcloud plot
# visualize the enriched ligands in the first condition
computeEnrichmentScore(net.down, species = 'human')
# visualize the enriched ligands in the second condition
computeEnrichmentScore(net.up, species = 'human')

#(A) Circle plot
pathways.show <- c("CXCL")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name =
                        paste(pathways.show, names(object.list)[i]))}
#Heatmap plot
pathways.show <- c("CXCL")
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

#Visualize gene expression distribution.
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("NL", "LS")) # set factor level
plotGeneExpression(cellchat, signaling = "CXCL", split.by = "datasets", colors.ggplot = T)

save(object.list, file = "cellchat_object.list_humanSkin_NL_LS.RData")
save(cellchat, file = "cellchat_merged_humanSkin_NL_LS.RData")

#https://ndownloader.figshare.com/files/25957094
#https://ndownloader.figshare.com/files/25957634

cellchat.E13 <- readRDS("cellchat_embryonic_E13.rds")
cellchat.E13 <- updateCellChat(cellchat.E13)
cellchat.E14 <- readRDS("cellchat_embryonic_E14.rds")
cellchat.E14 <- updateCellChat(cellchat.E14)

group.new = levels(cellchat.E14@idents) # Define the cell labels to lift up
cellchat.E13 <- liftCellChat(cellchat.E13, group.new)
object.list <- list(E13 = cellchat.E13, E14 = cellchat.E14)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

# Circle plot
pathways.show <- c("WNT")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))}

save(object.list, file = "cellchat_object.list_embryonic_E13_E14.RData")
save(cellchat, file = "cellchat_merged_embryonic_E13_E14.RData")

#https://figshare.com/articles/dataset/10X_visium_data_for_spatial-informed_cell-cell_communication/23621151
load("visium_mouse_cortex_annotated.RData")

# Gene expression data
data.input = GetAssayData(visium.brain, slot = "data", assay = "SCT") # normalized data matrix
# User assigned cell labels
meta = data.frame(labels = Idents(visium.brain), row.names = names(Idents(visium.brain)))
# Spatial locations of spots from full (NOT high/low) resolution images
spatial.locs = GetTissueCoordinates(visium.brain, scale = NULL, cols = c("imagerow", "imagecol"))
# Scale factors and spot diameters of the full resolution images

#https://figshare.com/articles/dataset/spatial_imaging_data_for_the_10X_visium_brain_dataset/23709726
scale.factors = jsonlite::fromJSON('scalefactors_json.json')
scale.factors = list(spot.diameter = 65, spot = scale.factors$spot_diameter_fullres)

#Create a CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels", datatype = "spatial", coordinates = spatial.locs, scale.factors = scale.factors)
CellChatDB <- CellChatDB.mouse
# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
#CellChatDB.use <- CellChatDB # simply use the default CellChatDB to use all CellChatDB for cell-cell communication analysis
# set the used database in the object
cellchat@DB <- CellChatDB.use

#Identify over-expressed ligands or receptors.
cellchat <- subsetData(cellchat) # This step is necessary even ifusing the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

#cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, distance.use = TRUE, interaction.range = 250, scale.distance = 0.01)
#long long time
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, distance.use = TRUE,  scale.distance = 0.01)#long long time
# Filter the cell-cell communication
cellchat <- filterCommunication(cellchat, min.cells = 10)
# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
# Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, file = "cellchat_visium_mouse_cortex.Rds")

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name ="Interaction weights/strength")
pathways.show <- c("CXCL")
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
# Spatial plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout= "spatial", edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5)

#Updating the ligand-receptor interaction database CellChatDB

#https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Update-CellChat