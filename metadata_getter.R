if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")

#get metadata from GEO
library(GEOquery)
library("R.utils")
geo_id <- "GSE64810"
gse <- getGEO(geo_id,GSEMatrix=TRUE)
pdata <- pData(gse$GSE64810_series_matrix.txt.gz)

#remove characters after colon in variable values
pdata <- sapply(pdata, function(x) gsub(".*: ","",x) )
#replace character string NA with 0. These values are not recognized by is.na
pdata[pdata=="NA"]<-0
write.csv(pdata,'~/downloads/final_project/metadata.csv')
