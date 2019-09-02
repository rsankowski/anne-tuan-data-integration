
library(readr)
library(data.table)
library(tidyverse)
library(tools)

start.time <- Sys.time()

# Load necessary files
gene2iso <- read.csv("/home/roman/Documents/Single cell analysis/wgEncodeGencodeBasicVM9_clean_genes2groups.tsv",sep="\t",header=FALSE)
ercc     <- read.csv("/home/roman/Documents/Single cell analysis/ERCC_Controls_Analysis_length_ext.txt",sep="\t",header=TRUE)
geneid_ercc    <- data.frame("GENEID" = c(as.character(gene2iso$V1), as.character(ercc$Re.sort.ID)), stringsAsFactors = F )


# Specify directories
file_paths <- list("data/counts/")
names(file_paths) <- file_paths


# Search for all files in the specified directories and extract files by a given extension
files_list            <- lapply(file_paths, list.files, recursive=T)
files_by_ext          <- lapply(files_list, function(x){x[endsWith(x, suffix=".coutt.csv")]} )

# Get complete paths to all files
all_file_paths        <- unlist(lapply(seq_along(files_by_ext), function(x) {  paste(names(files_by_ext[x]), files_by_ext[[x]], sep="") } ))
names(all_file_paths) <- lapply(strsplit(all_file_paths,split="/"), function(x) { sub(".coutt.csv","",x[length(x)]) } )

# Calculate md5sums to check for duplicated
md5sums      <- lapply(all_file_paths, function(x) {md5sum(x)} )


# Check for duplicated data
if ( sum(duplicated(md5sums)) == 0 ){
  print("No duplicated data found.")
} else if ( sum(duplicated(md5sums)) >= 1 ){
  print("Warning! There is duplicated data loaded.")
  print(md5sums[duplicated(md5sums)])
}

# Check for duplicated names
file_names <- unname(unlist(lapply(strsplit(unlist(files_by_ext),split = "/"),tail,1)))
if ( sum(duplicated(file_names)) == 0 ){
  print("No duplicated file names found.")
} else if ( sum(duplicated(file_names)) >= 1 ){
  print("Warning! There are duplicated names used")
  print(file_names[duplicated(file_names)])
}

####
#### LOADING
####
# Loading data using lapply
data_list   <- lapply(all_file_paths, function(x) {fread(x, header= T)} )

# Add dataset name prefix to all columns, Merge with remaining gene names
for (d in names(data_list)) { 
  colnames(data_list[[d]]) <- c("GENEID", paste(d, "_",1:192,sep="" ))#; 
  #setkey(data_list[[d]], GENEID)#; 
  #data_list[[d]] <- merge( geneid_ercc, data_list[[d]], by="GENEID",all.x = T )  
}

# Check whether rownames order is unchanged 
# table(sapply(data_list,dim))
# table(sapply(data_list, function(x) { sum(geneid_ercc$GENEID == x[order(match(x$GENEID, geneid_ercc$GENEID)),]$GENEID) }) )

# Cbind list of data.tables and removing the GENEID column from data.tables
data_list_cbind <- reduce(data_list, full_join, by = "GENEID")
data_list_cbind <- data_list_cbind %>% full_join(geneid_ercc)
data_list_cbind[is.na(data_list_cbind)]   <- 0


# Data as tibble & replace NA values with 0
#merged_data_tib <- tibble::as_tibble(data_list_cbind)
# merged_data_tib[is.na(merged_data_tib)] <- 0


# Measure time
Sys.time() - start.time

# Remove ERCC and mitochondrial genes
#prdata           <- as.data.frame(merged_data_tib)
prdata <- data_list_cbind
rownames(prdata) <- prdata$GENEID
prdata$GENEID    <- NULL
prdata           <- prdata[grep("^(ERCC|mt|Gm|Rik)",rownames(prdata),invert=TRUE),]

save(prdata, file = "data/prdata-mesencephalon.Robj")
