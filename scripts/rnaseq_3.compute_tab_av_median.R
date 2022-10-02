library(ohchibi)
library(dplyr)
library(ggtree)
library(ComplexHeatmap)
library(org.At.tair.db)
library(clusterProfiler)


set.seed(130816)


Dat_rnaseq <- readRDS(file = "./cleandata/dat_rnaseq.RDS")

Dat <- Dat_rnaseq$Dat_bpFALSE

melted <- Dat$Tab %>% melt %>%
  dplyr::rename(.data = .,gene_id = Var1,SampleId = Var2) %>%
  merge(Dat$Map, by = "SampleId") %>%
  droplevels 



Tab_exp <- acast(data = melted,formula = gene_id~SampleId,value.var = "value") %>% t
Tab_n <- Tab_exp

for(i in 1:ncol(Tab_exp)){
  
  x <- Tab_exp[,i]
  max_val <- quantile(x,0.99)
  min_val <- quantile(x,0.01)
  x[x>max_val] <- max_val
  x[x<min_val] <- min_val
  Tab_n[,i] <- x
  
}

Tab_av <- Tab_n %>% t %>% melt %>% 
  dplyr::rename(.data = .,gene_id = Var1,SampleId = Var2) %>%
  merge(Dat$Map, by = "SampleId") %>% 
  acast(data = .,formula = gene_id~group,
        fun.aggregate = function(x)median(x,na.rm = T),value.var = "value",fill = 0)

saveRDS(object = Tab_av,file = "./cleandata/dat_tab_av_rnaseq.RDS")


rm(list=ls())
gc()
