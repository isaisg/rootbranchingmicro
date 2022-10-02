library(ohchibi)
library(xlsx)
library(dplyr)


set.seed(130816)


df <- read.xlsx(file = "./rawdata/rawdata_ionomics_col0strains_2021_12_07_Athaliana_Summary.xlsx",sheetIndex = 3)
df <- df[-1,]
df$ICP.MS.tube <- df$ICP.MS.tube %>%
  gsub(pattern = " ",replacement = "")
colnames(df)[1] <- "UiD"

df <- which(!(is.na(df$UiD))) %>%
  df[.,] %>% droplevels

Map <- df[,1:4]
Tab <- df[,5:ncol(df)] %>% as.matrix

Map$UiD <- Map$UiD %>% as.character
rownames(Tab) <- Map$UiD 
rownames(Map) <- Map$UiD

Map$Genotype <- Map$Genotype %>%
  gsub(pattern = "Col 0",replacement = "Col-0") %>%
  factor() %>%
  relevel(ref = "Col-0")

Map$Bacteria <- Map$Bacteria %>% 
  factor %>% 
  relevel(ref = "NB")


Map$GenotypeBacteria <- paste0(Map$Genotype,"|",Map$Bacteria)  %>%
  factor %>%relevel(ref = "Col-0|NB")


colnames(Tab) <- colnames(Tab) %>%
  gsub(pattern = "\\..*",replacement = "")

mode(Tab) <- "numeric"
Tab[is.na(Tab)] <- 0

#Remove outliers ###
Tab_norm <- Tab
for(i in 1:ncol(Tab)){
  
  
  max_val <- quantile(Tab[,i],0.995,na.rm = T)
  min_val <- quantile(Tab[,i],0.005,na.rm = T)
  
  x <- Tab[,i]
  
  x[x > max_val] <- max_val
  x[x < min_val] <- min_val
  
  Tab_norm[,i] <- x
  
}

Tab_norm_log10 <- log10(Tab_norm+1)

Dat_raw <- create_dataset(Tab = Tab_norm %>% t,Map = Map)
Dat_raw_log10 <- create_dataset(Tab = Tab_norm_log10 %>% t,Map = Map)

Tab_norm_z <-  Tab_norm %>% scale(center = T,scale = T) %>% as.matrix
Tab_norm_log10_z <-  Tab_norm_log10 %>% scale(center = T,scale = T) %>% as.matrix


Dat_raw <- create_dataset(Tab = Tab_norm %>% t,Map = Map)
Dar_z <- create_dataset(Tab = Tab_norm_z %>% t,Map = Map)
Dat_log10 <- create_dataset(Tab = Tab_norm_log10 %>% t,Map = Map)
Dat_log10_z <- create_dataset(Tab = Tab_norm_log10_z %>% t,Map = Map)

#Save structures ###
mlist <- list(
  
  Dat_raw = Dat_raw,
  Dar_raw_z = Dar_z,
  Dat_log10 = Dat_log10,
  Dat_log10_z = Dat_log10_z
  
  
)

saveRDS(object = mlist,file = "./cleandata/dat_rootbranching_ionome_mutants_strains.RDS")
rm(list=ls())
gc()

