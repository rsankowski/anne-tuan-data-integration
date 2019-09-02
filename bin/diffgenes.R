#Differential genes
library(tidyverse)
library(RaceID)

date = Sys.Date()
load('data/sc.Robj')

subDir <- "data/Cluster specific genes"
dir.create('data/Cluster specific genes')
dir.create('data/Cluster specific genes/Up')
dir.create('data/Cluster specific genes/Down')

for (i in unique(sc@cpart)) {
  
  tryCatch({idhwt <-names(sc@cpart[sc@cpart %in% c(i)])
  idhmut <- names(sc@cpart[sc@cpart %in% sc@cpart[sc@cpart != i]])
  diffexp <- diffexpnb(getfdata(sc,n=c(idhmut,cl)), A=idhmut, B=cl )
  diffexpgenes <- diffexp[["res"]]
  diffexpgenes <- subset(diffexpgenes, diffexpgenes$pval < 0.05)
  diffexpgenes <- subset(diffexpgenes, abs(diffexpgenes$log2FoldChange) > 0)
  diffexpgenes$GENEID <- rownames(diffexpgenes)
  diffexpgenes <- diffexpgenes %>% dplyr::arrange(padj)
  #diffexpgenes$GENEID <- gsub('__.*', '', diffexpgenes$GENEID)
  diffexpgenes$Cluster <- rep(i, nrow(diffexpgenes))
  diffgene_up <- diffexpgenes[diffexpgenes$log2FoldChange > 0, ]
  diffgene_down <- diffexpgenes[diffexpgenes$log2FoldChange < 0, ]
  write.csv(diffgene_up, file = paste0(subDir, '/Up/', as.character(date), "-diffgenes_cl", as.character(i), "_up_idhmut.csv"))
  write.csv(diffgene_down, file = paste0(subDir, '/Down/', as.character(date), "-diffgenes_cl", as.character(i), "_down_idhmut.csv"))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  tryCatch({
    png(paste0(subDir, '/Up/', as.character(date),'-MA-plot-Cl', as.character(i) ,'_idhmut.png'), res = 300, width =7, height = 7, units = 'in')
    plotdiffgenesnb(diffexp,show_names = T, pthr=.01, lthr = .0001, mthr = .0001, padj=T)
    dev.off()
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  for (j in unique(sc@cpart)) {
    tryCatch({
      idhwt <-names(sc@cpart[sc@cpart %in% i])
      idhmut <- names(sc@cpart[sc@cpart %in% j])
      diffexp <- diffexpnb(getfdata(sc,n=c(idhmut,cl)), A=idhmut, B=cl )
      diffexpgenes <- diffexp[["res"]]
      diffexpgenes <- subset(diffexpgenes, diffexpgenes$pval < 0.05)
      diffexpgenes <- subset(diffexpgenes, abs(diffexpgenes$log2FoldChange) > 0)
      diffexpgenes$GENEID <- rownames(diffexpgenes)
      diffexpgenes <- diffexpgenes %>% dplyr::arrange(padj)
      #diffexpgenes$GENEID <- gsub('__.*', '', diffexpgenes$GENEID)
      diffexpgenes$Cluster <- rep(i, nrow(diffexpgenes))
      diffgene_up <- diffexpgenes[diffexpgenes$log2FoldChange > 0, ]
      diffgene_down <- diffexpgenes[diffexpgenes$log2FoldChange < 0, ]
      write.csv(diffgene_up, file = paste0(subDir, '/Up/', as.character(date), "-diffgenes_cl", as.character(i), "_up_vs_cl", as.character(j), ".csv"))
      write.csv(diffgene_down, file = paste0(subDir, '/Down/', as.character(date), "-diffgenes_cl", as.character(i), "_down_vs_cl", as.character(j), ".csv"))
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
    tryCatch({
      png(paste0(subDir, '/Up/', as.character(date), '-MA-plot-Cl', as.character(i) , "_vs_cl", as.character(j), ".png"), res = 300, width =7, height = 7, units = 'in')
      plotdiffgenesnb(diffexp,show_names = T, pthr=.01, lthr = .0001, mthr = .0001, padj=T)
      dev.off()
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
    
  }
  
}

#export diffgenes

load_data <- function(path) { 
  files <- dir(path, pattern = '\\.csv', full.names = TRUE)
  tables <- lapply(files, read.csv)
  do.call(rbind, tables)
}

up_genes <- load_data(file.path('data/Cluster specific genes/Up'))
up_genes <- up_genes[!grepl("(^Gm|^Rp|Rik|mt-)", up_genes$GENEID),]
up_genes %>% 
  filter(up_genes$padj<0.05, up_genes$log2FoldChange > 1) %>% 
  dplyr::arrange(Cluster, desc(log2FoldChange)) %>%
  write.csv('data/up-genes.csv')


up_genes %>% 
  filter(padj<0.05, log2FoldChange > 1) %>% #, up_genes$log2FoldChange > .5
  dplyr::arrange(Cluster, desc(log2FoldChange)) %>%
  dplyr::distinct(Cluster, GENEID) %>%
  write.csv('data/unique-up-genes.csv')


#compare diagnoses

subDir <- "data/IDH-status-Cluster specific genes"
dir.create('data/IDH-status-Cluster specific genes')
dir.create('data/IDH-status-Cluster specific genes/modc-trajectory-idhwt')
dir.create('data/IDH-status-Cluster specific genes/modc-trajectory-idhmut')

#modcs
          for (i in c(18,11,15,5, 9)) {
            
            tryCatch({
              idhwt <-names(sc@cpart[sc@cpart %in% c(i) & grepl('(1|2|3)$', colnames(sc@ndata))])
            idhmut <- names(sc@cpart[sc@cpart %in% c(i) & grepl('(4|5|6)$', colnames(sc@ndata))])
            diffexp <- diffexpnb(getfdata(sc,n=c(idhmut,idhwt)), A=idhmut, B=idhwt )
            diffexpgenes <- diffexp[["res"]]
            diffexpgenes <- subset(diffexpgenes, diffexpgenes$pval < 0.05)
            diffexpgenes <- subset(diffexpgenes, abs(diffexpgenes$log2FoldChange) > 0)
            diffexpgenes$GENEID <- rownames(diffexpgenes)
            diffexpgenes <- diffexpgenes %>% dplyr::arrange(padj)
            #diffexpgenes$GENEID <- gsub('__.*', '', diffexpgenes$GENEID)
            diffexpgenes$Cluster <- rep(i, nrow(diffexpgenes))
            diffgene_up <- diffexpgenes[diffexpgenes$log2FoldChange > 0, ]
            diffgene_down <- diffexpgenes[diffexpgenes$log2FoldChange < 0, ]
            write.csv(diffgene_up, file = paste0(subDir, '/modc-trajectory-idhwt/', as.character(date), "-diffgenes_cl", as.character(i), "_up_idhwt.csv"))
            write.csv(diffgene_down, file = paste0(subDir, '/modc-trajectory-idhmut/', as.character(date), "-diffgenes_cl", as.character(i), "_up_idhmut.csv"))
            }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
            
            tryCatch({
              png(paste0(subDir, '/modc-trajectory-idhwt/', as.character(date),'-MA-plot-Cl', as.character(i) ,'_idhmut.png'), res = 300, width =7, height = 7, units = 'in')
              plotdiffgenesnb(diffexp,show_names = T, pthr=.01, lthr = .0001, mthr = .0001, padj=T)
              dev.off()
            }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
            
            
          }
          
          #export diffgenes
          
          load_data <- function(path) { 
            files <- dir(path, pattern = '\\.csv', full.names = TRUE)
            tables <- lapply(files, read.csv)
            do.call(rbind, tables)
          }
          
          up_genes <- load_data(file.path('data/IDH-status-Cluster specific genes/modc-trajectory-idhwt'))
          up_genes <- up_genes[!grepl("(^Gm|^Rp|Rik|mt-)", up_genes$GENEID),]
          up_genes %>% 
            filter(up_genes$padj<0.05) %>% 
            dplyr::arrange(Cluster, desc(log2FoldChange)) %>%
            write.csv('data/modc-trajectory-idhwt-up-genes.csv')
          
          
          up_genes %>% 
            filter(padj<0.05) %>% #, up_genes$log2FoldChange > .5
            dplyr::arrange(Cluster, desc(log2FoldChange)) %>%
            dplyr::distinct(Cluster, GENEID) %>%
            write.csv('data/modc-trajectory-idhwt-unique-up-genes.csv')
          
          up_genes <- load_data(file.path('data/IDH-status-Cluster specific genes/modc-trajectory-idhmut'))
          up_genes <- up_genes[!grepl("(^Gm|^Rp|Rik|mt-)", up_genes$GENEID),]
          up_genes %>% 
            filter(padj<0.05) %>% 
            dplyr::arrange(Cluster, desc(log2FoldChange)) %>%
            write.csv('data/modc-trajectory-idhmut-up-genes.csv')
          
          
          up_genes %>% 
            filter(padj<0.05) %>% #, up_genes$log2FoldChange > .5
            dplyr::arrange(Cluster, desc(log2FoldChange)) %>%
            dplyr::distinct(Cluster, GENEID) %>%
            write.csv('data/modc-trajectory-idhmut-unique-up-genes.csv')

#macs
            subDir <- "data/IDH-status-Cluster specific genes"
            dir.create('data/IDH-status-Cluster specific genes')
            dir.create('data/IDH-status-Cluster specific genes/mac-trajectory-idhwt')
            dir.create('data/IDH-status-Cluster specific genes/mac-trajectory-idhmut')
            
            #macs
            for (i in c(18,11,1,31)) {
              
              tryCatch({
                idhwt <-names(sc@cpart[sc@cpart %in% c(i) & grepl('(1|2|3)$', colnames(sc@ndata))])
                idhmut <- names(sc@cpart[sc@cpart %in% c(i) & grepl('(4|5|6)$', colnames(sc@ndata))])
                diffexp <- diffexpnb(getfdata(sc,n=c(idhmut,idhwt)), A=idhmut, B=idhwt )
                diffexpgenes <- diffexp[["res"]]
                diffexpgenes <- subset(diffexpgenes, diffexpgenes$pval < 0.05)
                diffexpgenes <- subset(diffexpgenes, abs(diffexpgenes$log2FoldChange) > 0)
                diffexpgenes$GENEID <- rownames(diffexpgenes)
                diffexpgenes <- diffexpgenes %>% dplyr::arrange(padj)
                #diffexpgenes$GENEID <- gsub('__.*', '', diffexpgenes$GENEID)
                diffexpgenes$Cluster <- rep(i, nrow(diffexpgenes))
                diffgene_up <- diffexpgenes[diffexpgenes$log2FoldChange > 0, ]
                diffgene_down <- diffexpgenes[diffexpgenes$log2FoldChange < 0, ]
                write.csv(diffgene_up, file = paste0(subDir, '/mac-trajectory-idhwt/', as.character(date), "-diffgenes_cl", as.character(i), "_up_idhwt.csv"))
                write.csv(diffgene_down, file = paste0(subDir, '/mac-trajectory-idhmut/', as.character(date), "-diffgenes_cl", as.character(i), "_up_idhmut.csv"))
              }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
              
              tryCatch({
                png(paste0(subDir, '/mac-trajectory-idhwt/', as.character(date),'-MA-plot-Cl', as.character(i) ,'_idhmut.png'), res = 300, width =7, height = 7, units = 'in')
                plotdiffgenesnb(diffexp,show_names = T, pthr=.01, lthr = .0001, mthr = .0001, padj=T)
                dev.off()
              }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
              
              
            }
            
            #export diffgenes
            
            load_data <- function(path) { 
              files <- dir(path, pattern = '\\.csv', full.names = TRUE)
              tables <- lapply(files, read.csv)
              do.call(rbind, tables)
            }
            
            up_genes <- load_data(file.path('data/IDH-status-Cluster specific genes/mac-trajectory-idhwt'))
            up_genes <- up_genes[!grepl("(^Gm|^Rp|Rik|mt-)", up_genes$GENEID),]
            up_genes %>% 
              filter(up_genes$padj<0.05) %>% 
              dplyr::arrange(Cluster, desc(log2FoldChange)) %>%
              write.csv('data/mac-trajectory-idhwt-up-genes.csv')
            
            
            up_genes %>% 
              filter(padj<0.05) %>% #, up_genes$log2FoldChange > .5
              dplyr::arrange(Cluster, desc(log2FoldChange)) %>%
              dplyr::distinct(Cluster, GENEID) %>%
              write.csv('data/mac-trajectory-idhwt-unique-up-genes.csv')
            
            up_genes <- load_data(file.path('data/IDH-status-Cluster specific genes/mac-trajectory-idhmut'))
            up_genes <- up_genes[!grepl("(^Gm|^Rp|Rik|mt-)", up_genes$GENEID),]
            up_genes %>% 
              filter(up_genes$padj<0.05) %>% 
              dplyr::arrange(Cluster, desc(log2FoldChange)) %>%
              write.csv('data/mac-trajectory-idhmut-up-genes.csv')
            
            
            up_genes %>% 
              filter(padj<0.05) %>% #, up_genes$log2FoldChange > .5
              dplyr::arrange(Cluster, desc(log2FoldChange)) %>%
              dplyr::distinct(Cluster, GENEID) %>%
              write.csv('data/mac-trajectory-idhmut-unique-up-genes.csv')

