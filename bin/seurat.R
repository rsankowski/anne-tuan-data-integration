library(Seurat)
library(tidyverse)

prdata <- Read10X_h5('data/Mm_CD45_IDH_SMAR_6mice_raw_feature_bc_matrix.h5', use.names = TRUE, unique.features = TRUE)
