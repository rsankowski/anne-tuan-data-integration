##### params #####
library(RaceID)

load.data  <- T
do.init <- T
mircol     <- "MIRTAR.TS.0.9"
date <- Sys.Date()


csamp <- c("Retina")

############################################ functions start ############################################

if ( load.data ){
  if (do.init){
    data   <- list()
    data.add <- list()
    cdiff <- list()
    sco   <- list()
    pca   <- list()
    entr  <- list()
    ltree <- list()
    gcl   <- list()
    lgres <- list()
    somd   <- list()
    marker <- list()
    net    <- list()
    regentr <- list()
    hyperv  <- list()
    traj    <- list()
  }
  
  gene2iso <- read.csv("/home/roman/Documents/Single cell analysis/wgEncodeGencodeBasicVM9_clean_genes2groups.tsv",sep="\t",header=FALSE)
  ercc     <- read.csv("/home/roman/Documents/Single cell analysis/ERCC_Controls_Analysis_length_ext.txt",sep="\t",header=TRUE)
  for ( n in grep("nb.cell",names(ercc)) ){ercc[,n] <- 5/6 * ercc[,n]} # original CEL-seq
  
  n1 <- 1:97
  n2 <- c(1,98:193)
  n4 <- 1:193
  
  for ( s in csamp ){  
    if (s %in% c("Retina")){
      data.add[[s]] <- list()
      for ( m in c("t","b","c") ){
        z <- list()
        if ( s == "Retina" )  { nl <- c( P1_1="D0_1", 
                                         P1_2="D7c_1", 
                                         #P1_3="D7x_1", 
                                         #P1_6="D30c_1", 
                                         P2_7="D0_2", 
                                         #P2_8="D30x_1", 
                                         #P2_9="D30c_2", 
                                         #P2_10="D30x_2", 
                                         P3_1="D0_3", 
                                         P3_2="D7c_2", 
                                         #P3_3="D7x_2", 
                                         P3_4="D0_4", 
                                         P4_5="D7c_3"#, 
                                         #P4_8="D7x_3"#, 
                                         #P4_9="D30c_3", 
                                         #P4_10="D30x_3"
        );
        fl <- list()
        for ( i in names(nl) ) fl[[i]] <- n4;  di <- "" }
        
        
        for ( sl in names(nl) ){
          if ( length(di) > 0 ){
            x <- read.csv(paste("data/counts",di,"/",sl,".cout",m,".csv",sep=""),sep="\t",header=TRUE)
          }else{
            x <- read.csv(paste("data/counts",sl,".cout",m,".csv",sep=""),sep="\t",header=TRUE)
          }
          x <- x[,fl[[sl]]]
          x <- merge(data.frame(GENEID=c(as.vector(gene2iso[,1]),as.vector(ercc[,1])),GROUP=c(as.vector(gene2iso[,2]),as.vector(ercc[,1]))),x,by="GENEID",all=TRUE)[,-1]
          names(x)[1] <- "GENEID"
          x[is.na(x[,2]),-1] <- 0
          x <- x[order(x$GENEID),]
          z[[sl]]  <- x
          names(z[[sl]])  <- c("GENEID",paste(nl[[sl]],sub("X","",names(z[[sl]])[-1]),sep="_"))
        }
        for ( i in 1:length(z) ) y <- if ( i == 1 ) z[[i]] else merge(y,z[[i]],by="GENEID")
        row.names(y) <- y$GENEID
        y <- y[,-1]
        y <- y[,apply(y,2,sum)>0]
        if ( m == "t" ){
          data[[s]] <- y
        }else{
          data.add[[s]][[m]] <- y
        }
      }
    }
  }
}


  #Create prdata file
prdata <- data[[s]][grep("^(ERCC|mt-)",row.names(data[[s]]),invert=TRUE),]
cs <- apply(prdata,2,sum)
prdata <- prdata[,cs > 500]
cs <- cs[cs > 500]
f <- t(prdata["Kcnq1ot1",])/cs < .02

#define sigcor
sigcor <- function(x,y,cthr=.4){
  if ( min(var(x),var(y)) == 0 ) return(NA)
  fit <- lm(x ~ y)
  pv <- as.data.frame(summary(fit)[4])[2,4]
  y <- as.data.frame(summary(fit)[4])[2,1]
  if ( is.na(pv) | is.na(y) ) return( NA )
  z <- sign(y)*sqrt(summary(fit)$r.square)
  if ( is.na(z) ) return(NA)
  if ( pv < .01 & abs(z) >= cthr ) return(z) else return(NA)
}


h <- rep(TRUE,nrow(prdata))
for ( g in c("Kcnq1ot1")){
  z <- apply(prdata,1,function(x,y) sigcor(x,y),y=t(prdata[g,]))
  h <- h & ( is.na(z) | z < .65 )
}
prdata <- prdata[h,]

# save prdata file
save(prdata, file =paste0('data/prdata.Robj'))


