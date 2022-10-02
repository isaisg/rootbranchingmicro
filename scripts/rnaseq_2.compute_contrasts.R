library(ohchibi)
library(dplyr)
library(DESeq2)

set.seed(130816)


Dats_rnaseq <- readRDS(file = "./cleandata/dat_rnaseq.RDS")



#### Perform contrasts ####
dds <- Dats_rnaseq$dds_bpFALSE

Tab <- dds@assays@data@listData$counts
write.table(x = Tab,
            file = "./cleandata/matrix_rnaseq_counts_rootbranching.tsv",
            append = FALSE,quote = FALSE,sep = "\t",row.names = TRUE,col.names = TRUE)

list_comparisons <- list(
  
  c("group","arf7_19_NB","Col0_NB"),
  c("group","gnom1_NB","Col0_NB"),
  c("group","lbd16_NB","Col0_NB"),
  c("group","nph4_NB","Col0_NB"),
  
  c("group","Col0_L180","Col0_NB"),
  c("group","Col0_L196","Col0_NB"),
  c("group","Col0_L264","Col0_NB"),
  c("group","Col0_L339","Col0_NB"),
  c("group","Col0_L343","Col0_NB"),
  c("group","Col0_L344","Col0_NB"),
  c("group","Col0_L359","Col0_NB"),
  c("group","Col0_L384","Col0_NB"),
  c("group","Col0_L412","Col0_NB"),
  c("group","Col0_L469","Col0_NB"),
  c("group","Col0_L76","Col0_NB"),
  c("group","Col0_RCL115","Col0_NB"),
  c("group","Col0_RMF217","Col0_NB"),
  c("group","Col0_RMF27","Col0_NB"),
  c("group","Col0_RMF331","Col0_NB"),
  c("group","Col0_RMF8","Col0_NB"),
  
  c("group","arf7_19_L344","arf7_19_NB"),
  c("group","gnom1_L196","gnom1_NB"),
  c("group","gnom1_L339","gnom1_NB"),
  c("group","gnom1_L469","gnom1_NB"),
  c("group","gnom1_RMF27","gnom1_NB"),
  c("group","lbd16_L180","lbd16_NB"),
  c("group","lbd16_L343","lbd16_NB"),
  c("group","lbd16_L359","lbd16_NB"),
  c("group","lbd16_L384","lbd16_NB"),
  c("group","lbd16_L412","lbd16_NB"),
  c("group","lbd16_L76","lbd16_NB"),
  c("group","lbd16_RMF217","lbd16_NB"),
  c("group","lbd16_RMF8","lbd16_NB"),
  c("group","nph4_L344","nph4_NB"),
  c("group","nph4_L359","nph4_NB"),
  c("group","nph4_RCL115","nph4_NB"),
  
  
  
  c("group","arf7_19_L344","Col0_L344"),
  c("group","gnom1_L196","Col0_L196"),
  c("group","gnom1_L339","Col0_L339"),
  c("group","gnom1_L469","Col0_L469"),
  c("group","gnom1_RMF27","Col0_RMF27"),
  c("group","lbd16_L180","Col0_L180"),
  c("group","lbd16_L343","Col0_L343"),
  c("group","lbd16_L359","Col0_L359"),
  c("group","lbd16_L384","Col0_L384"),
  c("group","lbd16_L412","Col0_L412"),
  c("group","lbd16_L76","Col0_L76"),
  c("group","lbd16_RMF217","Col0_RMF217"),
  c("group","lbd16_RMF8","Col0_RMF8"),
  c("group","nph4_L344","Col0_L344"),
  c("group","nph4_L359","Col0_L359"),
  c("group","nph4_RCL115","Col0_RCL115")
  
)

Res_contrasts <- NULL

for(comparison in list_comparisons){
  
  cat("Working on ",comparison,"\n")
  
  Res_contrasts <- results(object = dds,contrast = comparison) %>%
    as.data.frame  %>%
    dplyr::mutate(.data = .,gene_id = rownames(.)) %>%
    dplyr::relocate(.data = .,gene_id) %>%
    tibble::remove_rownames(.data = .) %>%
    dplyr::mutate(.data = .,
                  Contrast = comparison[2:3] %>% paste0(collapse = "_vs_")) %>%
    rbind(Res_contrasts,.)

}

Res_contrasts_FALSE <- Res_contrasts

#### True
dds <- Dats_rnaseq$dds_bpTRUE

list_comparisons <- list(
  
  c("group","arf7_19_NB","Col0_NB"),
  c("group","gnom1_NB","Col0_NB"),
  c("group","lbd16_NB","Col0_NB"),
  c("group","nph4_NB","Col0_NB"),
  
  c("group","Col0_L180","Col0_NB"),
  c("group","Col0_L196","Col0_NB"),
  c("group","Col0_L264","Col0_NB"),
  c("group","Col0_L339","Col0_NB"),
  c("group","Col0_L343","Col0_NB"),
  c("group","Col0_L344","Col0_NB"),
  c("group","Col0_L359","Col0_NB"),
  c("group","Col0_L384","Col0_NB"),
  c("group","Col0_L412","Col0_NB"),
  c("group","Col0_L469","Col0_NB"),
  c("group","Col0_L76","Col0_NB"),
  c("group","Col0_RCL115","Col0_NB"),
  c("group","Col0_RMF217","Col0_NB"),
  c("group","Col0_RMF27","Col0_NB"),
  c("group","Col0_RMF331","Col0_NB"),
  c("group","Col0_RMF8","Col0_NB"),
  
  c("group","arf7_19_L344","arf7_19_NB"),
  c("group","gnom1_L196","gnom1_NB"),
  c("group","gnom1_L339","gnom1_NB"),
  c("group","gnom1_L469","gnom1_NB"),
  c("group","gnom1_RMF27","gnom1_NB"),
  c("group","lbd16_L180","lbd16_NB"),
  c("group","lbd16_L343","lbd16_NB"),
  c("group","lbd16_L359","lbd16_NB"),
  c("group","lbd16_L384","lbd16_NB"),
  c("group","lbd16_L412","lbd16_NB"),
  c("group","lbd16_L76","lbd16_NB"),
  c("group","lbd16_RMF217","lbd16_NB"),
  c("group","lbd16_RMF8","lbd16_NB"),
  c("group","nph4_L344","nph4_NB"),
  c("group","nph4_L359","nph4_NB"),
  c("group","nph4_RCL115","nph4_NB"),
  
  
  
  c("group","arf7_19_L344","Col0_L344"),
  c("group","gnom1_L196","Col0_L196"),
  c("group","gnom1_L339","Col0_L339"),
  c("group","gnom1_L469","Col0_L469"),
  c("group","gnom1_RMF27","Col0_RMF27"),
  c("group","lbd16_L180","Col0_L180"),
  c("group","lbd16_L343","Col0_L343"),
  c("group","lbd16_L359","Col0_L359"),
  c("group","lbd16_L384","Col0_L384"),
  c("group","lbd16_L412","Col0_L412"),
  c("group","lbd16_L76","Col0_L76"),
  c("group","lbd16_RMF217","Col0_RMF217"),
  c("group","lbd16_RMF8","Col0_RMF8"),
  c("group","nph4_L344","Col0_L344"),
  c("group","nph4_L359","Col0_L359"),
  c("group","nph4_RCL115","Col0_RCL115")
  
)

Res_contrasts <- NULL

for(comparison in list_comparisons){
  
  cat("Working on ",comparison,"\n")
  
  Res_contrasts <- results(object = dds,contrast = comparison) %>%
    as.data.frame  %>%
    dplyr::mutate(.data = .,gene_id = rownames(.)) %>%
    dplyr::relocate(.data = .,gene_id) %>%
    tibble::remove_rownames(.data = .) %>%
    dplyr::mutate(.data = .,
                  Contrast = comparison[2:3] %>% paste0(collapse = "_vs_")) %>%
    rbind(Res_contrasts,.)
  
}

Res_contrasts_TRUE <- Res_contrasts

###Create list
mlist <- list(
  Res_contrasts_bpFALSE = Res_contrasts_FALSE,
  Res_contrasts_bpTRUE = Res_contrasts_TRUE
)


saveRDS(object = mlist,file = "./cleandata/res_rnaseq_contrasts.RDS")
rm(list=ls())
gc()
