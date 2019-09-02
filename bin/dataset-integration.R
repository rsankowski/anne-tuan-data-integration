#devtools::install_github(repo = "satijalab/seurat", ref = "release/3.0")
library(tidyverse)
library(Seurat)
library(Matrix)
library(RaceID)

date <- Sys.Date()

load('/home/roman/Documents/Single cell analysis/Comparison-between-regions-Anne/RaceID/20180523-retina-vs-PV-prdata.Robj')
prdata <- prdata[,!grepl('^(RZ|Peak|Precl2)', colnames(prdata))]
prdata$GENEID <- rownames(prdata)

hyth <- merge(read.csv('/home/roman/Documents/Single cell analysis/20181007-Comparison-between-regions-Anne/A761.all_samples.gencode_genomic.corrected_merged.csv', stringsAsFactors = F, sep = '\t'),
              read.csv('/home/roman/Documents/Single cell analysis/20181007-Comparison-between-regions-Anne/P215_A326.all_samples.gencode_genomic.corrected_merged.csv', stringsAsFactors = F, sep = '\t'),
              by = 'GENEID', all = T)

hyth <- merge(hyth,
              read.csv('/home/roman/Documents/Single cell analysis/20181007-Comparison-between-regions-Anne/A558.all_samples.gencode_genomic.corrected_merged.csv', stringsAsFactors = F, sep = '\t'),
              by = 'GENEID', all = T)
hyth$GENEID <- gsub('_.*', '', hyth$GENEID)


#merge them
prdata_all <- merge(prdata, hyth[-1,], by='GENEID', all = T)
rownames(prdata_all) <- prdata_all$GENEID
prdata <- na.omit(prdata_all[,-1])
cells_to_exclude <- c(read_csv("data/20181008-regions-micr-cells-to-exclude.csv")[[2]],
                      read_csv("data/20181008-regions-micr-more-cells-to-exclude.csv")[[2]],
                      read_csv("data/20181008-regions-micr-still-more-cells-to-exclude.csv")[[2]])

prdata <- prdata[, !colnames(prdata) %in% cells_to_exclude]

              
#metainfo
vars <- data.frame(row.names=names(prdata),batch=sub("(_|)\\d.+","",names(prdata)))
vars$batch <- ifelse(grepl('Onset1_7*', rownames(vars)), 'Batch1', 
                     ifelse(grepl('Precl_6*', rownames(vars)), 'Batch1', 
                            ifelse(grepl('Naive*', rownames(vars)), 'Batch1', 
                                   ifelse(grepl('Precl4d*', rownames(vars)), 'Batch1', 
                                          ifelse(grepl('A', rownames(vars)), 'Anne_Plate1', 
                                                 ifelse(grepl('B6_WT', rownames(vars)), 'Anne_Plate2', 
                                                        ifelse(grepl('B6_MG', rownames(vars)), 'Anne_Plate3', 
                                                               ifelse(grepl('^(B6_HFD|B6_CT)', rownames(vars)), 'Anne_Plate4','Batch1'))))))))

#merge metadata
              metadata <- vars
              colnames(metadata) <- "dataset"

              #extract common rows from Olah et al dataset
              micr.data <- prdata

#save datasets
              save(metadata, file = 'data/metadata.Robj')
              save(micr.data, file = 'data/merged-counts.Robj')
              
micr <- CreateSeuratObject(counts = micr.data, meta.data = metadata)
micr.list <- SplitObject(object = micr, split.by = "dataset")
for (i in 1:length(x = micr.list)) {
  micr.list[[i]] <- NormalizeData(object = micr.list[[i]], verbose = FALSE)
  micr.list[[i]] <- FindVariableFeatures(object = micr.list[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
}

reference.list <- micr.list[c(unique(vars$batch))]
micr.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
micr.integrated <- IntegrateData(anchorset = micr.anchors, dims = 1:30)

#plot
library(ggplot2)
library(cowplot)
# switch to integrated assay. The variable features of this assay are automatically set during
# IntegrateData
DefaultAssay(object = micr.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
micr.integrated <- ScaleData(object = micr.integrated, verbose = FALSE)
save(micr.integrated, file = "data/integrated-counts.Robj")
micr.integrated <- RunPCA(object = micr.integrated, npcs = 10, verbose = FALSE)
micr.integrated <- RunUMAP(object = micr.integrated, reduction = "pca", dims = 1:10)
p1 <- DimPlot(object = micr.integrated, reduction = "umap", group.by = "dataset")
p2 <- DimPlot(object = micr.integrated, reduction = "umap", group.by = "dataset", label = TRUE, 
              repel = TRUE) #+ NoLegend()
plot_grid(p1, p2)
ggsave(paste0('plots/', date, '-integration-mouse-human-dataset-C1-C9-10PCs.png'))
ggsave(paste0('plots/', date, '-integration-mouse-human-dataset-C1-C9-10PCs.pdf'))

svg(paste0('plots/', date, '-integration-mouse-human-dataset-C1-C9-10PCs.svg'))
plot_grid(p1, p2)
dev.off()

#export cell embeddings
write_csv(data.frame(ID=rownames(micr.integrated@reductions$umap@cell.embeddings),micr.integrated@reductions$umap@cell.embeddings), paste0('data/', date, '-umap-cell-embeddings.csv'))



#
micr.query <- micr.list[["human"]]
micr.anchors <- FindTransferAnchors(reference = micr.integrated, query = micr.query, 
                                        dims = 1:30)
predictions <- TransferData(anchorset = micr.anchors, refdata = micr.integrated$celltype, 
                            dims = 1:30)
micr.query <- AddMetaData(object = micr.query, metadata = predictions)
micr.query$prediction.match <- micr.query$predicted.id == micr.query$celltype
table(micr.query$prediction.match)
table(micr.query$predicted.id)

