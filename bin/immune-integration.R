#from url https://satijalab.org/seurat/v3.0/immune_alignment.html

library(tidyverse)
library(Seurat)
library(Matrix)
library(RaceID)
library(cowplot)

date <- Sys.Date()

load('/home/roman/Documents/Single cell analysis/Comparison-between-regions-Anne/RaceID/20180523-retina-vs-PV-prdata.Robj')
prdata <- prdata[,!grepl('^(RZ|Peak|Precl2)', colnames(prdata))]

hyth <- read.csv('/home/roman/Documents/Single cell analysis/20181007-Comparison-between-regions-Anne/A761.all_samples.gencode_genomic.corrected_merged.csv', stringsAsFactors = F, sep = '\t')

hyth2 <- read.csv('/home/roman/Documents/Single cell analysis/20181007-Comparison-between-regions-Anne/P215_A326.all_samples.gencode_genomic.corrected_merged.csv', stringsAsFactors = F, sep = '\t')

hyth3 <- read.csv('/home/roman/Documents/Single cell analysis/20181007-Comparison-between-regions-Anne/A558.all_samples.gencode_genomic.corrected_merged.csv', stringsAsFactors = F, sep = '\t')

hyth$GENEID <- gsub('_.*', '', hyth$GENEID)
hyth <- hyth[!duplicated(hyth$GENEID),]
hyth <- hyth[!grepl('mt-' , hyth$GENEID),]

hyth2$GENEID <- gsub('_.*', '', hyth2$GENEID)
hyth2 <- hyth2[!duplicated(hyth2$GENEID),]
hyth2<- hyth2[!grepl('mt-' , hyth2$GENEID),]

hyth3$GENEID <- gsub('_.*', '', hyth3$GENEID)
hyth3 <- hyth3[!duplicated(hyth3$GENEID),]
hyth3 <- hyth3[!grepl('mt-' , hyth3$GENEID),]

rownames(hyth) <- hyth$GENEID
rownames(hyth2) <- hyth2$GENEID
rownames(hyth3) <- hyth3$GENEID

#extract the same genes
genes <- rownames(prdata)[rownames(prdata) %in% intersect(intersect(rownames(hyth), rownames(hyth2)),rownames(hyth3))]
prdata <- prdata[genes,]
hyth <- hyth[genes,]
hyth2 <- hyth2[genes,]
hyth3 <- hyth3[genes,]

# merge datasets
eae.data <- prdata
hfd.data <- hyth[,-1]
hfd.data2 <- hyth2[,-1]
hfd.data3 <- hyth3[,-1]

# Set up control object
eae <- CreateSeuratObject(counts = eae.data, project = "IMMUNE_eae", min.cells = 5)
eae$stim <- "eae"
eae <- subset(eae, subset = nFeature_RNA > 500)
eae <- NormalizeData(eae, verbose = FALSE)
eae <- FindVariableFeatures(eae, selection.method = "vst", nfeatures = 2000)

# Set up stimulated object
hfd <- CreateSeuratObject(counts = hfd.data, project = "IMMUNE_hfd", min.cells = 5)
hfd$stim <- "hfd"
hfd <- subset(hfd, subset = nFeature_RNA > 500)
hfd <- NormalizeData(hfd, verbose = FALSE)
hfd <- FindVariableFeatures(hfd, selection.method = "vst", nfeatures = 2000)

# Set up stimulated object2
hfd2 <- CreateSeuratObject(counts = hfd.data2, project = "IMMUNE_hfd", min.cells = 5)
hfd2$stim <- "hfd"
hfd2 <- subset(hfd2, subset = nFeature_RNA > 500)
hfd2 <- NormalizeData(hfd2, verbose = FALSE)
hfd2 <- FindVariableFeatures(hfd2, selection.method = "vst", nfeatures = 2000)

# Set up stimulated object3
hfd3 <- CreateSeuratObject(counts = hfd.data3, project = "IMMUNE_hfd", min.cells = 5)
hfd3$stim <- "hfd"
hfd3 <- subset(hfd3, subset = nFeature_RNA > 500)
hfd3 <- NormalizeData(hfd3, verbose = FALSE)
hfd3 <- FindVariableFeatures(hfd3, selection.method = "vst", nfeatures = 2000)


#integrate datasets
immune.anchors <- FindIntegrationAnchors(object.list = list(eae, hfd, hfd2, hfd3), dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)

#integrated analysis
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

#plot
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

save(immune.combined, file = "data/dataset-integration-hfd-eae.Robj")


#compare fnx data
load("data/merged_prdata.Robj")

prdata_mat <- as.data.frame(as.matrix(prdata_mat))
fnx.data <- prdata_mat[,grepl("^P[1-9]", colnames(prdata_mat))]
hfd.data <- prdata_mat[,!grepl("^P[1-9]", colnames(prdata_mat))]

# Set up control object
fnx <- CreateSeuratObject(counts = fnx.data, project = "IMMUNE_fnx", min.cells = 5)
fnx$stim <- "fnx"
fnx <- subset(fnx, subset = nFeature_RNA > 500)
fnx <- NormalizeData(fnx, verbose = FALSE)
fnx <- FindVariableFeatures(fnx, selection.method = "vst", nfeatures = 2000)

# Set up stimulated object
hfd <- CreateSeuratObject(counts = hfd.data, project = "IMMUNE_hfd", min.cells = 5)
hfd$stim <- "hfd"
hfd <- subset(hfd, subset = nFeature_RNA > 500)
hfd <- NormalizeData(hfd, verbose = FALSE)
hfd <- FindVariableFeatures(hfd, selection.method = "vst", nfeatures = 2000)

#integrate datasets
immune.anchors <- FindIntegrationAnchors(object.list = list(fnx, hfd), dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)

#integrated analysis
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

#plot
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

save(immune.combined, file = "data/dataset-integration-hfd-fnx.Robj")
