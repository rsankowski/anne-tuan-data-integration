#RaceID4 
library(tidyverse)
library(viridis)
library(RaceID)
library(Seurat)
library(Matrix)

date = Sys.Date()

load("data/prdata-mesencephalon.Robj")
prdata$GENEID <- rownames(prdata)

#add hypothalamic microglia
hyth <- full_join(read.csv('/home/roman/Documents/Single cell analysis/Comparison-between-regions-Anne/ (4W 16W)all_samples.gencode_genomic.coutt_merged.csv', stringsAsFactors = F, sep = '\t'),
              read.csv('/home/roman/Documents/Single cell analysis/Comparison-between-regions-Anne/(1D1W) all_samples.gencode_genomic.coutt_merged.csv', stringsAsFactors = F, sep = '\t'))

hyth <- full_join(hyth,
              read.csv('/home/roman/Documents/Single cell analysis/Comparison-between-regions-Anne/(16W 9W 1W 3D) all_samples.gencode_genomic.coutt_merged.csv', stringsAsFactors = F, sep = '\t'))
#hyth$GENEID <- gsub('_.*', '', hyth$GENEID)


#merge them
prdata_all <- full_join(prdata, hyth)
prdata_all <- na.omit(prdata_all)
rownames(prdata_all) <- prdata_all$GENEID
prdata <-na.omit(prdata_all)
prdata <- prdata[, !grepl("16weeks_HFD|9weeks_HFD|9wks_HFD|16wks_HFD|A[1-4]", colnames(prdata))]
#prdata <- prdata[, !grepl("16weeks|9weeks|9wks|16wks|A[1-4]", colnames(prdata))]

prdata <- prdata[, !colnames(prdata) %in% read.csv('data/non-microglia-cell-ids.csv')[[2]]]
prdata_mat <- Matrix(as.matrix(prdata[,-which(colnames(prdata) == "GENEID")]), sparse = T)
rownames(prdata_mat) <- gsub("_.*", "" , prdata$GENEID)

#save merged matrix
save(prdata_mat, file="data/merged_prdata.Robj")

#RaceID
sc <- SCseq(prdata_mat)

#correct batch effects
vars <- data.frame(row.names=colnames(sc@expdata),batch=colnames(sc@expdata))
vars$batch <- ifelse(grepl('^P[1-9]', rownames(vars)), 'Batch1', 
                     ifelse(grepl('A', rownames(vars)), 'Anne_Plate1', 
                                                 ifelse(grepl('B6_WT', rownames(vars)), 'Anne_Plate2', 
                                                        ifelse(grepl('B6_MG', rownames(vars)), 'Anne_Plate3', 
                                                               ifelse(grepl('^RZ', rownames(vars)), 'Eye', 'Batch2')))))

b = list(vars[,"batch"])

# filtering of expression data
sc <- filterdata(sc, 
                 mintotal=500,
                 LBatch=b,
                 bmode="scran",
)

sc <- CCcorrect(sc, 
                dimR = T, 
                nComp = 20,
                CGenes = c('Jun',
                           'Fos',
                           'Zfp36',
                           'Atf3',
                           'Hspa1a|Hspa1b',
                           'Dusp1',
                           'Egr1',
                           'Malat1'))

sc <- compdist(sc,metric="pearson")
sc <- clustexp(sc) 

plotsaturation(sc,disp=FALSE)
plotsaturation(sc,disp=TRUE)
plotjaccard(sc)

sc <- clustexp(sc,cln=15,sat=FALSE) 
sc <- findoutliers(sc)
plotbackground(sc)
plotsensitivity(sc)
plotoutlierprobs(sc)
clustheatmap(sc)
    
    pdf(paste0('plots/heatmaps/', date,"-clustheatmap.pdf"))
    clustheatmap(sc, final = T)
    dev.off()
    
    sc <- comptsne(sc)
    sc <- compfr(sc,knn=10)

plotmap(sc)
plotmap(sc,fr=TRUE)
dev.off()


plotexpmap(sc,"Mrc1",logsc=F,fr=F)
plotexpmap(sc,"Lyve1",logsc=F,fr=F)
plotexpmap(sc,"Cd163",logsc=F,fr=F)
plotexpmap(sc,"Tmem119",logsc=F,fr=F)
plotexpmap(sc,"Cx3cr1",logsc=F,fr=F)
plotexpmap(sc,"Hexb",logsc=F,fr=F)
plotexpmap(sc,"Ptprc",logsc=F,fr=F)
plotexpmap(sc,"Cd3e",logsc=F,fr=F)
plotexpmap(sc,"Itgam",logsc=F,fr=F)
plotexpmap(sc,"Cd8a",logsc=F,fr=F)
plotexpmap(sc,"Cd4",logsc=F,fr=F)
plotexpmap(sc,"H2-Aa",logsc=F,fr=F)
plotexpmap(sc,"Zbtb46",logsc=F,fr=F)
plotexpmap(sc,"Mog",logsc=F,fr=F)
plotexpmap(sc,"Mbp",logsc=F,fr=F)
plotexpmap(sc,"Gfap",logsc=F,fr=F)
plotexpmap(sc,"Wfdc17",logsc=F,fr=F)
plotexpmap(sc,"Cd79a",logsc=F,fr=F)
plotexpmap(sc,"Cst3",logsc=F,fr=F)
plotexpmap(sc,"Nkg7",logsc=F,fr=F)
plotexpmap(sc,"Rho",logsc=F,fr=F)
plotexpmap(sc,"S100a11",logsc=F,fr=F)
plotexpmap(sc,"Fn1",logsc=F,fr=F)
plotexpmap(sc,"Gfap",logsc=F,fr=F)
plotexpmap(sc,"Cd209a",logsc=F,fr=F)
plotexpmap(sc,"Rho",logsc=F,fr=F)

plotexpmap(sc,"Mrc1",logsc=F,fr=T)
plotexpmap(sc,"Lyve1",logsc=F,fr=T)
plotexpmap(sc,"Cd163",logsc=F,fr=T)
plotexpmap(sc,"Tmem119",logsc=F,fr=T)
plotexpmap(sc,"Cx3cr1",logsc=F,fr=T)
plotexpmap(sc,"Ptprc",logsc=F,fr=T)
plotexpmap(sc,"Cd3e",logsc=F,fr=T)
plotexpmap(sc,"Itgam",logsc=F,fr=T)
plotexpmap(sc,"Cd8a",logsc=F,fr=T)
plotexpmap(sc,"H2-Aa",logsc=F,fr=T)
plotexpmap(sc,"Zbtb46",logsc=F,fr=T)
plotexpmap(sc,"Ly6c2",logsc=F,fr=T)
plotexpmap(sc,"Cd177",logsc=F,fr=T)
plotexpmap(sc,"Igkc",logsc=F,fr=T)
plotexpmap(sc,"Wfdc17",logsc=F,fr=T)
plotexpmap(sc,"Plp1",logsc=F,fr=T)
plotexpmap(sc,"Mog",logsc=F,fr=T)
#plot marker genes

#Save sc file
save(sc, file = paste0('data/sc.Robj'))

#identify myeloid cells
#write.csv(data.frame(names(sc@cpart[sc@cpart %in% c(7)])), 'data/non-microglia-cell-ids.csv')
write.csv(data.frame(names(sc@cpart[sc@cpart %in% c(5)])), 'data/more-non-microglia-cell-ids.csv')
