#Differential genes
library(tidyverse)
library(RaceID)

date = Sys.Date()
load('data/sc.Robj')

ltr <- Ltree(sc)
ltr <- compentropy(ltr)
ltr <- projcells(ltr,cthr=5,nmode=TRUE,fr=FALSE)
ltr <- projback(ltr,pdishuf=100)
ltr <- lineagegraph(ltr)
ltr <- comppvalue(ltr,pthr=0.01)

save(ltr, file = 'data/ltr.Robj')
plotgraph(ltr,scthr=0.01,showCells=FALSE,showTsne=TRUE)

svg(paste0('plots/tsne/', date, '-lineage-graph-tsne-plot.svg'), width = 8.57, height = 5.79)
plotgraph(ltr,scthr=0.01,showCells=FALSE,showTsne=TRUE)
dev.off()

pdf(paste0('plots/tsne/', date, '-lineage-graph-tsne-plot.pdf'), width = 8.57, height = 5.79)
plotgraph(ltr,scthr=0.01,showCells=FALSE,showTsne=TRUE)
dev.off()

x <- compscore(ltr,scthr=0.01)
plotdistanceratio(ltr)
plotspantree(ltr)
plotprojections(ltr)

#lineage tree for moDCs
ctrl_d7 <- c(5,4,3,11,14)
ctrl_d3 <- c(5,4,3,6)

#pseudotemporal ctrl_d7
n <- cellsfromtree(ltr,ctrl_d7)
x <- getfdata(ltr@sc)
library(FateID)
fs  <- filterset(x,n=n$f)
s1d <- getsom(fs,nb=1000,alpha=.5)
ps  <- procsom(s1d,corthr=.85,minsom=3)
y    <- ltr@sc@cpart[n$f]
fcol <- ltr@sc@fcol
plotheatmap(ps$all.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)

pdf(paste0('plots/heatmaps/', date, '-ctrl-d7-trajectory-heatmap.pdf'), width = 8.57, height = 5.79)
plotheatmap(ps$all.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
dev.off()

#export node genes
modules <- data.frame('Node' = NA, 'Genes' = NA)
for (i in 1: max(ps$nodes)) {
  gene_names <- names(ps$nodes)[ps$nodes == i]
  gene_names <- gsub('_.*', '', gene_names)
  modules2 <- data.frame('Node' = NA, 'Genes' = gene_names)
  modules2$Node <- rep(as.character(i), nrow(modules2))
  modules <- rbind(na.omit(modules), modules2)
}

write.csv(modules, paste0('data/',date,'-nodes-stemid-vector-ctrl_d7.csv'))


#pseudotemporal ordering - ctrl_d3
n <- cellsfromtree(ltr,ctrl_d3)
x <- getfdata(ltr@sc)
library(FateID)
fs  <- filterset(x,n=n$f)
s1d <- getsom(fs,nb=1000,alpha=.5)
ps  <- procsom(s1d,corthr=.85,minsom=3)
y    <- ltr@sc@cpart[n$f]
fcol <- ltr@sc@fcol
plotheatmap(ps$all.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)

pdf(paste0('plots/heatmaps/', date, '-ctrl-d3-trajectory-heatmap.pdf'), width = 8.57, height = 5.79)
plotheatmap(ps$all.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
dev.off()

#eport node genes
modules <- data.frame('Node' = NA, 'Genes' = NA)
for (i in 1: max(ps$nodes)) {
  gene_names <- names(ps$nodes)[ps$nodes == i]
  gene_names <- gsub('_.*', '', gene_names)
  modules2 <- data.frame('Node' = NA, 'Genes' = gene_names)
  modules2$Node <- rep(as.character(i), nrow(modules2))
  modules <- rbind(na.omit(modules), modules2)
}

write.csv(modules, paste0('data/',date,'-nodes-stemid-vector-ctrl_d3.csv'))


