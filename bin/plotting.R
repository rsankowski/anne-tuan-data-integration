library(tidyverse)
library(viridis)
library(readxl)
library(RaceID)
library(RColorBrewer)
library(pheatmap)

date = Sys.Date()
load('data/sc.Robj')
source("/home/roman/Documents/Single cell analysis/Advanced-plots/20181025-sankowski-et-al-functions.R")


#build data frame
if (!file.exists("data/df.Robj")) {
  df <- data.frame('Cluster' = as.numeric(sc@cpart), sc@tsne) 
df$ID <- names(sc@cpart)
rownames(df) <- names(sc@cpart)

df$Condition <- ifelse(grepl('^(P1_1|P2_7|P3_1|P3_4)', rownames(df)), 'FNc_0d', 
                       ifelse(grepl('^(P1_2|P3_2|P4_5)', rownames(df)), 'FNc_7d', 
                              ifelse(grepl('^(P1_3|P3_3|P4_8)', rownames(df)), 'FNx_7d', 
                                     ifelse(grepl('^(P1_6|P2_9|P4_9)', rownames(df)), 'FNc_30d', 
                                            ifelse(grepl('(HFD_4weeks|HFD4wks)', rownames(df)), '4wk-HFD', 
                                                   ifelse(grepl('(9wks_CT|CT_10weeks|14wks_NCD|24wks_NCD|CT_16weeks|9weeks_CT|16weeks_CT|1week_CT)', rownames(df)), 'NC_Ctrl', 
                                                                 ifelse(grepl('HFD1d', rownames(df)), '1d-HFD', 
                                                                        ifelse(grepl('(HFD_16weeks|HFD16wks)', rownames(df)), '16wks-HFD', 
                                                                               ifelse(grepl('(1week_HFD)', rownames(df)), '1w-HFD', 
                                                                                             ifelse(grepl('3days_HFD', rownames(df)), '3d-HFD', 
                                                                                                    ifelse(grepl('HFD9wks', rownames(df)), '9wk-HFD', 
                                                                                                           ifelse(grepl('^(P2_8|P2_10|P4_10)', rownames(df)), 'FNx_30d','rest'))))))))))))


df$Condition <- factor(df$Condition, levels = c('FNc_0d', 'FNc_7d',  'FNx_7d', 'FNc_30d', 'FNx_30d',   'NC_Ctrl', '1d-HFD', '3d-HFD', '1w-HFD', '4wk-HFD', '9wk-HFD', '16wks-HFD'))

df$Treatment <- ifelse(grepl("FNc", df$Condition), "FN-Ctrl",
                       ifelse(grepl("FNx", df$Condition), "FNx", 
                       ifelse(grepl("HFD", df$Condition), "HFD", "NC")))

df$Treatment <- factor(df$Treatment, levels=c("NC", "HFD", "FN-Ctrl", "FNx"))

cell_numbers <- numeric()
for (i in 1:max(df$Cluster, na.rm = T)) {
  cell_numbers[i] <- length(na.omit(df$Cluster[df$Cluster==i]))
}
names(cell_numbers) <- c(1:max(df$Cluster, na.rm = T))
retain_cl <- as.numeric(names(cell_numbers[cell_numbers > dim(sc@ndata)[2]/100]))
if (!file.exists('data/ord_clust.Robj')) {
  ord_clust <- clustheatmap(sc, final = T)
  save(ord_clust, file = 'data/ord_clust.Robj')
} else {
  load('data/ord_clust.Robj')
}
ord_clust <- ord_clust[ord_clust %in% retain_cl]
df$Cluster <- factor(df$Cluster, levels = ord_clust)

df <- df[df$Cluster %in% retain_cl,]
df$Cluster <- factor(df$Cluster, levels = ord_clust)
colnames(df)[2:3] <- c("V1", "V2")
save(df, file = "data/df.Robj")

} else {
  load("data/df.Robj")
}

#Plot tsne map
  tsne <- tsne_plot(df, FILL = df$Condition, fill_colors = c(colors_pat), point_outline = "black", point_size = 4, line_width = 0.25) 

tsne

ggsave(paste0('plots/tsne/', date, '-condition-tsne-plot.pdf'), width = 8.57, height = 5.79)  

svg(paste0('plots/tsne/', date, '-condition-tsne-plot.svg'), width = 8.57, height = 5.79)
tsne
dev.off()

#marimekko plot - patients
mosaicGG2(df, "Cluster", "Condition", c(colors_pat, colors_many), rect_col = 'black', line_width = 0.1) +
  scale_fill_brewer(palette = "Set3")
ggsave(paste0('plots/others/', date,'-compartment-marimekko-cluster-plot.pdf'))

#marimekko stat plot
mosaicGG(df, "Cluster", "Condition", rect_col = 'black', line_width = 0.1)
ggsave(paste0('plots/others/', date,'-conditions-marimekko-cluster-stat-plot.pdf'))

#plot treatment
#Plot tsne map
tsne <- tsne_plot(df, FILL = df$Treatment, fill_colors = c(colors_pat), point_outline = "black", point_size = 4, line_width = 0.25) 

tsne

ggsave(paste0('plots/tsne/', date, '-Treatment-tsne-plot.pdf'), width = 8.57, height = 5.79)  

svg(paste0('plots/tsne/', date, '-Treatment-tsne-plot.svg'), width = 8.57, height = 5.79)
tsne
dev.off()

#marimekko plot - patients
mosaicGG2(df, "Cluster", "Treatment", c(colors_pat, colors_many), rect_col = 'black', line_width = 0.1) +
  scale_fill_brewer(palette = "Set3")
ggsave(paste0('plots/others/', date,'-compartment-marimekko-cluster-plot.pdf'))

#marimekko stat plot
mosaicGG(df, "Cluster", "Treatment", rect_col = 'black', line_width = 0.1)
ggsave(paste0('plots/others/', date,'-Treatments-marimekko-cluster-stat-plot.pdf'))



#Plot cluster tsne map
tsne <- tsne_plot(df, FILL = df$Cluster, fill_colors = c(colors_pat, colors_many), point_outline = "black", point_size = 3, line_width = 0.25)

tsne

ggsave(paste0('plots/tsne/', date, '-clusters-tsne-plot.pdf'))  

svg(paste0('plots/tsne/', date, '-clusters-tsne-plot.svg'), width = 8.57, height = 5.79)
tsne
dev.off()


#plot cell signatures
signature_genes <-  read_excel('/home/roman/Documents/Single cell analysis/EAE_Final/Cluster-information-dimRed.xlsx', 'Core signature',skip = 2)

signature_genes <- data.frame("ahr"=c("Ahr", "Ahrr", "Cyp1a1", "Cyp1b1", "Entpd1", "Hif1a", "Il1b",
                                              "Il27", "Klf4", "Pparg", "Stat1", "Stat3", "Tiparp", "Vegfa" ), stringsAsFactors = F)

source('~/Documents/Single cell analysis/Advanced-plots/20190102_plot_expmap.R')

for (i in colnames(signature_genes)) {
  tryCatch({
    svg(paste0('plots/tsne/', date, as.character(i), '-gene_signature.svg'), width = 8.57, height = 5.79)
    pl <- plot_expmap(gene=c(na.omit(signature_genes[[i]])), point_size = 2.2)
    print(pl)
    dev.off()
    
    svg(paste0('plots/tsne/', date, as.character(i), '-gene_signature-logsc.svg'), width = 8.57, height = 5.79)
    pl <- plot_expmap(gene=c(na.omit(signature_genes[[i]])), point_size = 2.2, logsc = T)
    print(pl)
    dev.off()
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  #on.exit(dev.off())
}     


#plot single cell gene expression
load_data <- function(path) { 
  files <- dir(path, pattern = '\\.csv', full.names = TRUE)
  tables <- lapply(files, read.csv)
  do.call(rbind, tables)
}

up_genes <- load_data(file.path('data/Cluster specific genes/Up'))


source('~/Documents/Single cell analysis/Advanced-plots/20190102_plot_expmap.R')
plot_genes <-  as.character(unique(up_genes$GENEID))

plot_genes <- stringr::str_to_title(read_csv('data/genes-mirco.csv')[['Gene']])

for (i in plot_genes) {
  tryCatch({
    svg(paste0('plots/tsne/', i, 'logsc.svg'), width = 8.57, height = 5.79)
    pl <- plot_expmap(gene=c(i), point_size = 2.2, logsc = T)
    print(pl)
    dev.off()
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  #on.exit(dev.off())
}     

#heatmaps

load_data <- function(path) { 
  files <- dir(path, pattern = '\\.csv', full.names = TRUE)
  tables <- lapply(files, read.csv)
  do.call(rbind, tables)
}

up_genes <- load_data(file.path('data/Cluster specific genes/Up'))

up_genes <- up_genes[!grepl("(^Gm|^RP[0-9]|Rik|Hist)", up_genes$GENEID),]

gene_names <- up_genes[up_genes$padj<0.05 & up_genes$log2FoldChange >1,] %>%
  group_by(Cluster) %>%
  dplyr::arrange(Cluster, log2FoldChange) %>%
  top_n(n = 30, wt=log2FoldChange) 
genes <- unique(as.character(gene_names$GENEID))

exp <- as.data.frame(t(as.matrix(sc@ndata[genes,])+0.1))
nam <- as.data.frame(sc@cpart[sc@cpart %in% retain_cl])
nam$newid <- paste(rownames(nam),nam$`sc@cpart[sc@cpart %in% retain_cl]`,sep = "_cl")
exp.new <- merge(exp, nam, by = "row.names")
rownames(exp.new) <- exp.new$newid 
exp.new <- exp.new[,c(3:ncol(exp.new)-1)] # here you have to increase 11 to one more number per gene you will add
exp.new$`sc@cpart[sc@cpart %in% retain_cl]` <- factor(exp.new$`sc@cpart[sc@cpart %in% retain_cl]`, levels = ord_clust)
exp.new <- exp.new[order(exp.new$`sc@cpart[sc@cpart %in% retain_cl]`),]
log10exp <- log10(exp.new[,-ncol(exp.new)])
colnames(log10exp) <- gsub('_.*', '', colnames(log10exp))

#find position of column breaks for heatmap (basically where the new cluster starts)
breaks <- c(table(exp.new$`sc@cpart[sc@cpart %in% retain_cl]`))
breaks <- as.numeric(breaks)
a <- c()
for (i in 1:length(breaks)) {
  if (i==1) {
    a <- breaks[1]
  } else {a <- c(a, a[i-1]+breaks[i])}
}

#define annotation colors from url: https://stackoverflow.com/questions/33292067/pheatmap-annotation-colors-and-border
annotation_col = data.frame(ID = factor(exp.new$`sc@cpart[sc@cpart %in% retain_cl]`))
rownames(annotation_col)<-rownames(exp.new)
cols <- c(colors_pat, colors_many)[1:length(retain_cl)] #colorRampPalette(c(brewer.pal(length(breaks) ,name = 'Set3'), brewer.pal(length(breaks) - 10,name = 'Set2')))
names(cols) <- unique(annotation_col$ID)
annotation_colors <- list(ID = cols)


pheat <- pheatmap(t(log10exp),
                  show_colnames = F, 
                  #color = inferno(1000), 
                  #color = colorRampPalette(c("#FC8D59", "#FFFFBF", "#91CF60"))(length(palette.breaks) - 1),
                  #color = viridis(1000),
                  color = colorRampPalette(c('blue', 'white', 'red'))(1000),
                  cluster_cols = F, 
                  annotation_legend = T, 
                  cluster_rows = T, 
                  fontsize_row = 10, 
                  scale = 'column',
                  gaps_col=a[-length(a)],
                  annotation_col = annotation_col, 
                  annotation_colors = annotation_colors[1])


pdf('plots/heatmaps/single-cell-top-30-heatmap.pdf',height = 12, width = 12)
pheat
dev.off()


mat <- as.matrix(sc@ndata)
nd <- mat[genes,]
nd[is.na(nd)] <- 0
nd  <- t(nd[complete.cases(nd),])
clust_n <- sc@cpart[sc@cpart %in% retain_cl]
#clust_n <- clust_n[clust_n != 7 & clust_n != 8 & clust_n != 9]
mnd <- as.data.frame(cbind(clust_n,nd[rownames(nd) %in% names(clust_n),]))
mean_mnd <- aggregate(mnd[, 2:dim(mnd)[2]], list(mnd$clust_n), mean)
row.names(mean_mnd) <- paste("C",mean_mnd$Group.1,sep = "")
mean_mnd$Group.1 <- NULL
mean_mnd <- as.matrix(mean_mnd) + 0.1
gene <- as.data.frame(t(log10(mean_mnd)))
gene <- gene[complete.cases(gene),]
#row.names(gene) <- id2name(rownames(gene))
row.names(gene) <- gsub('__.*', '', row.names(gene))
pheatmap(gene, cluster_cols = T, cluster_rows=T,fontsize_row = 8, border_color = F, show_rownames = T,show_colnames = T, scale = "row")

pdf(paste0('plots/heatmaps/top-30-mean-heatmap.pdf'),height = 12, width = 6)
pheatmap(gene, cluster_cols = T, cluster_rows=T,fontsize_row = 8, border_color = F, show_rownames = T,show_colnames = T, scale = "row")
dev.off()


svg(paste0('plots/heatmaps/top-30-mean-heatmap.svg'),height = 12, width = 6)
pheatmap(gene, cluster_cols = T, cluster_rows=T,fontsize_row = 8, border_color = F, show_rownames = T,show_colnames = T, scale = "row")
dev.off()


#plot umaps from dataset integration
        load("data/dataset-integration-hfd-fnx.Robj")
        umap_embeddings <- immune.combined@reductions$umap@cell.embeddings %>%
          as.data.frame() %>%
          rownames_to_column(var="ID") 
        
        umap_embeddings$seurat_clusters <- immune.combined@meta.data$seurat_clusters
        
        df2 <- df %>% 
          left_join(umap_embeddings)
        df2$V1 <- df2$UMAP_1
        df2$V2 <- df2$UMAP_2
        #Plot tsne map
        tsne <- tsne_plot(df2, FILL = df$Treatment, fill_colors = c(colors_pat), point_outline = "black", point_size = 4, line_width = 0.25) 
        
        tsne
        
        ggsave(paste0('plots/tsne/', date, '-condition-umap-plot.pdf'), width = 8.57, height = 5.79)  
        
        svg(paste0('plots/tsne/', date, '-condition-umap-plot.svg'), width = 8.57, height = 5.79)
        tsne
        dev.off()
