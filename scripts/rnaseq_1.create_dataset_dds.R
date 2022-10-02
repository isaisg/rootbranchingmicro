library(ohchibi)
library(dplyr)
library(tximport)
library(DESeq2)


set.seed(130816)

#Import Salmon results 

files <- list.files("rawdata/data_rnaseq",
                    pattern = "quant.sf",recursive = T,full.names = T) %>%
  grep(pattern = "Col0",invert = T,value = T)

mnames <- files %>% 
  gsub(pattern = "rawdata/data_rnaseq/salmon_outdir_",replacement = "") %>%
  gsub(pattern = "\\/.*",replacement = "") 

names(files) <- mnames

tx2gene <- read.table("./rawdata/data_rnaseq/data_tx2gene_Arabidopsis_thaliana.TAIR10.tsv",sep = "\t") %>%
  dplyr::select(.data = .,c("V1","V3")) %>%
  dplyr::rename(.data = .,TXNAME = V1,GENEID = V3)


txi <- tximport(files = files,type = "salmon",tx2gene = tx2gene,
                countsFromAbundance = "no")

#Read metadata information
Map <- read.table(file = "./rawdata/data_rnaseq_metadata.tsv",header = T,sep = "\t")  %>%
  subset(Genotype != "POSITIVE")  %>%
  dplyr::rename(.data = .,SampleId = "Sample.ID",PlateNumber = "Plate.Number")  %>%
  dplyr::relocate(.data = .,SampleId) 

Map$Genotype <- Map$Genotype %>% 
  gsub(pattern = " ",replacement = "") %>%
  gsub(pattern = "/",replacement = "_")  %>% factor %>%
  relevel(ref = "Col0")

Map$Bacteria <- Map$Bacteria %>% factor %>%
  relevel(ref = "NB")

Map$SampleId <- Map$SampleId %>%
  gsub(pattern = " ",replacement = "_")

indices <- Map$SampleId  %>%grep(pattern = "_")
Map$SampleId[indices] <- paste0(Map$SampleId[indices],"_",Map$Position[indices])

Map$group <- paste0(Map$Genotype,"_",Map$Bacteria)

df_long <- data.frame(SampleId = mnames %>% gsub(pattern = "_L[0-9]$",replacement = ""),
           FullId = mnames)

Map <- merge(Map,df_long, by = "SampleId",all = TRUE)
rownames(Map) <- Map$FullId

#Define usable ids
Map <- match(mnames,Map$FullId) %>%
  Map[.,]

Map$group <- Map$group  %>%
  factor %>%
  relevel(ref = c("Col0_NB"))



dds <- DESeqDataSetFromTximport(txi = txi,Map,~ group)

dds <- DESeq2::collapseReplicates(object = dds, dds$SampleId)

dds <- DESeq(object = dds,betaPrior = TRUE)

#Create object to plot
mvst <- vst(object = dds,blind = F)

#Remove the batch effect from the vst matrix
mat <- assay(mvst)


Tab_z <- mat %>% t %>% scale (center = T,scale = F)


#Create matrix and project
Map_un <- Map[,1:6] %>% unique
rownames(Map_un) <- Map_un$SampleId
Map_un <- match(rownames(Tab_z),Map_un$SampleId) %>%
  Map_un[.,]

Dat <- create_dataset(Tab = Tab_z %>% t,Map = Map_un)

### Repeat with beta prior false
#Save structures
mlist_TRUE <- list(
  dds = dds,
  Dat = Dat
)


dds <- DESeqDataSetFromTximport(txi = txi,Map,~ group)

dds <- DESeq2::collapseReplicates(object = dds, dds$SampleId)

dds <- DESeq(object = dds,betaPrior = FALSE)

#Create object to plot
mvst <- vst(object = dds,blind = F)

#Remove the batch effect from the vst matrix
mat <- assay(mvst)


Tab_z <- mat %>% t %>% scale (center = T,scale = F)


#Create matrix and project
Map_un <- Map[,1:6] %>% unique
rownames(Map_un) <- Map_un$SampleId
Map_un <- match(rownames(Tab_z),Map_un$SampleId) %>%
  Map_un[.,]

Dat <- create_dataset(Tab = Tab_z %>% t,Map = Map_un)

mlist_FALSE <- list(
  dds = dds,
  Dat = Dat
)

end_list <- list(
  dds_bpFALSE = mlist_FALSE$dds,
  Dat_bpFALSE = mlist_FALSE$Dat,
  
  dds_bpTRUE = mlist_TRUE$dds,
  Dat_bpTRUE = mlist_TRUE$Dat
)



saveRDS(object = end_list,file = "./cleandata/dat_rnaseq.RDS")

rm(list=ls())
gc()
