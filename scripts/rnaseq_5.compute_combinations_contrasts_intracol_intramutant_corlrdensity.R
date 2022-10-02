library(ohchibi)
library(dplyr)
library(ComplexHeatmap)
library(biomaRt)

set.seed(130816)


Res <- readRDS(file = "./cleandata/res_rnaseq_contrasts.RDS") %$% Res_contrasts_bpFALSE


## Create data frame of metadata of contrasts
df_cons <- Res$Contrast %>% unique %>%
  data.frame(Contrast  = .)

df_cons$Bacteria <- df_cons$Contrast %>%
  gsub(pattern = "_vs.*",replacement = "") %>%
  gsub(pattern = "lbd16_",replacement = "") %>%
  gsub(pattern = "nph4_",replacement = "") %>%
  gsub(pattern = "gnom1_",replacement = "") %>%
  gsub(pattern = "Col0_",replacement = "") %>%
  gsub(pattern = "arf7_19_",replacement = "") 

df_cons$Type <- NA

indices <- NULL

for(i in 1:nrow(df_cons)){
  
  mbac <- df_cons[i,2]
  hit <- df_cons[i,1] %>% grep(pattern = paste0(".*",mbac,".*",mbac,".*"))
  if (length(hit)>0){
    indices <- c(indices,i)
  }
}

df_cons$Type[indices] <- "InterGenotype"

df_cons$Type[df_cons$Contrast %>% grep(pattern = ".*Col0.*Col0*",value = F)] <- "BacteriaInCol0"
df_cons$Type[df_cons$Contrast %>% grep(pattern = ".*NB.*NB.*")] <- "MutantvsCol0inNB"
df_cons$Type[which(is.na(df_cons$Type))] <- "BacteriaInMutant"

df_cons$Type <- df_cons$Type %>%
  factor(levels = c("BacteriaInCol0","BacteriaInMutant","MutantvsCol0inNB","InterGenotype"))

df_cons$Bacteria <- df_cons$Bacteria %>%  factor

#Keep solely contrasts
df_cons <- df_cons %>%
  dplyr::filter(.data = .,Type %in% c("BacteriaInCol0","BacteriaInMutant")) %>%
  droplevels


df_cons <- with(df_cons,order(Bacteria,Type)) %>%
  df_cons[.,]

### Add genotype ###
df_cons$Genotype <- NA
df_cons$Genotype[df_cons$Contrast %>% grep(pattern = "^Col0")] <- "Col0"
df_cons$Genotype[df_cons$Contrast %>% grep(pattern = "^lbd16_")] <- "lbd16"
df_cons$Genotype[df_cons$Contrast %>% grep(pattern = "^gnom1_")] <- "gnom1"
df_cons$Genotype[df_cons$Contrast %>% grep(pattern = "^nph4_")] <- "nph4"
df_cons$Genotype[df_cons$Contrast %>% grep(pattern = "^arf7_19_")] <- "arf7_19"

df_cons$BacteriaGenotype <- paste0(df_cons$Bacteria,"_",df_cons$Genotype) 


mcombinations <- df_cons$BacteriaGenotype %>%
  unique %>%
  grep(pattern = "NB",invert = T,value = T) %>%
  grep(pattern = "Col0",invert = T,value = T)

#Load results of correlation
df_cor <- readRDS(file = "./cleandata/res_cor_rnaseq_lrdensity.RDS") %$% Res_cor

df_cor$p_spearman[which(is.na(df_cor$p_spearman))] <- 1

#Use Rho
down_lr <- df_cor %>% subset(p_spearman < 0.05 & rho < 0)  %$% Gene
up_lr <- df_cor %>% subset(p_spearman < 0.05 & rho > 0)  %$% Gene
zero_lr <- df_cor %>% subset(p_spearman >= 0.05)  %$% Gene


padj_thres <- 0.05
lfc_thres <- 1

list_contrast_genes <- list()

for(mc in mcombinations){
  
  cat("Workin on",mc,"\n")
  
  df_a <- df_cons %>% subset(BacteriaGenotype == mc) %>%
    droplevels
  
  mst <- df_a$Bacteria[1] %>% as.character
  mgeno <- df_a$Genotype[1] %>% as.character
  
  df_b <- df_cons %>% subset(Bacteria == mst & Genotype == "Col0") %>%
    droplevels
  
  
  merged_c <- rbind(df_a,df_b) %>%
    dplyr::mutate(.data = .,Type = factor(Type,levels = df_cons$Type %>% levels))
  
  merged_c <- with(merged_c,order(Type)) %>%
    merged_c[.,]
  
  # Compute the contrasts based on 
  Res_sub <- Res %>% 
    dplyr::filter(.data = .,Contrast %in% (merged_c$Contrast %>% unique)) %>% 
    droplevels
  Res_sub$Direction <- 0
  Res_sub$Direction[which((Res_sub$padj < padj_thres) & (Res_sub$log2FoldChange >= lfc_thres))] <- 1
  Res_sub$Direction[which((Res_sub$padj < padj_thres) & (Res_sub$log2FoldChange <= -lfc_thres))] <- -1
  
  #Subset the universe of genes
  universe_genes <- Res_sub %>% subset(Direction == -1 | Direction == 1) %$% gene_id %>%
    as.character %>% unique
  
  Res_sub <- Res_sub %>% 
    dplyr::filter(.data = .,gene_id %in% universe_genes) %>%
    droplevels
  
  #Create a list to perform combination matrix
  list_internal <- list()
  
  for(i in 1:nrow(merged_c)){
    
    csubquery <- merged_c[i,1] %>% as.character
    
    typequery <- merged_c[i,3] %>% as.character
    
    minus <- Res_sub %>% 
      dplyr::filter(.data = .,Contrast %in% csubquery) %>%
      subset(Direction == -1) %$% gene_id
    
    zero <- Res_sub %>% 
      dplyr::filter(.data = .,Contrast %in% csubquery) %>%
      subset(Direction == 0) %$% gene_id
    
    plus <- Res_sub %>% 
      dplyr::filter(.data = .,Contrast %in% csubquery) %>%
      subset(Direction == 1) %$% gene_id
    
    list_internal[[paste0(typequery,"_-1")]] <- minus
    list_internal[[paste0(typequery,"_0")]] <- zero
    list_internal[[paste0(typequery,"_1")]] <- plus
    
  }
  
  #Append the quantitative genes
  list_internal[["CorLR_-1"]] <- down_lr
  list_internal[["CorLR_0"]] <- zero_lr
  list_internal[["CorLR_1"]] <- up_lr
  
  
  list_contrast_genes[[mc]] <- list_internal
}





#### Create combination matrices for all the dataset ####
Res_eset <- NULL
Res_size_eset <- NULL

for(mname in names(list_contrast_genes)){
  
  m = make_comb_mat(list_contrast_genes[[mname]])
  for (mc in comb_name(m)){
    
    mgenes_mc <- extract_comb(m = m,comb_name = mc) 
    
    Res_eset <- data.frame(BacteriaGenotype = mname,EnrichmentCode = mc,gene_id = mgenes_mc)   %>%
      rbind(Res_eset,.)
    
    Res_size_eset <- data.frame(BacteriaGenotype = mname,EnrichmentCode = mc, SizeEnrichmentSet = length(mgenes_mc)) %>%
      rbind(Res_size_eset ,.)
  }
  
}

Res_size_eset$EnrichmentCode <- Res_size_eset$EnrichmentCode %>%
  factor 
Res_size_eset$BacteriaGenotype <- Res_size_eset$BacteriaGenotype %>% factor


###### Load database #####
ensembl_arabidopsis <- useMart(biomart = "plants_mart",
                               dataset = "athaliana_eg_gene",
                               host = "plants.ensembl.org")


query_genes <- Res_eset$gene_id %>% as.character %>% unique

df_attr <- listAttributes(mart = ensembl_arabidopsis) 

res_bm <- getBM(attributes = c('ensembl_gene_id','tair_locus','tair_symbol','gene_biotype','description', 'go_id',"name_1006","definition_1006"),
                filters = "ensembl_gene_id",
                values = query_genes, 
                mart = ensembl_arabidopsis)

df_go_id <- aggregate(go_id~ensembl_gene_id,res_bm,FUN = function(x)paste0(x %>% unique,collapse = "|"))
df_go_name_1006 <- aggregate(name_1006~ensembl_gene_id,res_bm,FUN = function(x)paste0(x %>% unique,collapse = "|"))
df_go_definition_1006 <- aggregate(definition_1006~ensembl_gene_id,res_bm,FUN = function(x)paste0(x %>% unique,collapse = "|"))

dim(df_go_id)
dim(df_go_name_1006)
dim(df_go_definition_1006)

df_go <- merge(df_go_id,df_go_name_1006, by = "ensembl_gene_id") %>%
  merge(df_go_definition_1006,by = "ensembl_gene_id")

df_desc <- res_bm[,c('ensembl_gene_id','tair_locus','gene_biotype','description')] %>% unique 

df_tair_symbol <- res_bm[,c('ensembl_gene_id','tair_symbol')] %>% unique %>% 
  aggregate(tair_symbol~ensembl_gene_id,.,FUN = function(x)paste0(x  %>% unique,collapse = "|")) 

df_meta <- merge(df_desc,df_tair_symbol, by = "ensembl_gene_id") %>%
  merge(df_go, by = "ensembl_gene_id")  %>%
  dplyr::relocate(.data = .,ensembl_gene_id,tair_locus,tair_symbol)

##### Read gene information #####
df_gi <- read.table(file = "./rawdata/gene_families_sep_29_09_update.txt",header = T,
                    sep = "\t",comment.char = "",quote = "",fill =NA)
df_gi$Genomic_Locus_Tag <- df_gi$Genomic_Locus_Tag %>%
  gsub(pattern = " ",replacement = "") %>% toupper()

df_meta$GeneFamily <- match(df_meta$ensembl_gene_id,df_gi$Genomic_Locus_Tag) %>%
  df_gi$Gene_Family[.]

df_meta <- df_meta %>% 
  dplyr::rename(.data = .,go_name = name_1006,go_definition = definition_1006) %>%
  dplyr::relocate(.data = .,c(ensembl_gene_id,tair_locus,tair_symbol,gene_biotype,GeneFamily))


### Save structures 
mlist <- list(
  df_comb_gec = Res_eset,
  map_genes = df_meta,
  list_contrast_genes = list_contrast_genes
)

##


saveRDS(object = mlist,file = "./cleandata/res_rnaseq_combinations_contrasts_intracol_intramutant_corlrdensity.RDS")
rm(list = ls())
gc()
