library(ohchibi)
library(dplyr)
library(ComplexHeatmap)
library(biomaRt)
library(org.At.tair.db)
library(clusterProfiler)
library(ggtree)
library(xlsx)
library(DOSE)



set.seed(130816)

### Get genes ###
df_comb <- readRDS(file = "./cleandata/res_rnaseq_combinations_contrasts_intracol_intramutant_corlrdensity.RDS") %$% df_comb_gec

ids_ppp <- df_comb$EnrichmentCode %>%
  grep(pattern = "001001001")  %>%
  df_comb[.,] %$% gene_id %>% as.character %>% unique

ids_nnn <- df_comb$EnrichmentCode %>%
  grep(pattern = "100100100")  %>%
  df_comb[.,] %$% gene_id %>% as.character %>% unique



##
Res <- readRDS(file = "./cleandata/res_rnaseq_contrasts.RDS") %$% Res_contrasts_bpFALSE

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

###Get the df_cons <- Res$Contrast %>% unique %>%

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

df_cons <- df_cons %>% 
  subset(Type == "BacteriaInCol0") %>%
  droplevels

mcombinations <- df_cons$Contrast %>% as.character %>% unique

mids <- c(ids_ppp,ids_nnn) %>% unique



Res_sub <- Res %>%
  dplyr::filter(.data = .,gene_id %in% mids) %>%
  dplyr::filter(.data = .,Contrast %in% mcombinations) %>%
  droplevels



#Read regulon information
Res_regulons <- readRDS(file = "./cleandata/regulons_rnaseq.RDS")
eth <- Res_regulons$ethylene
cyto <- Res_regulons$cytokinin
aba <- readRDS("./rawdata/aba_robust_genes.RDS") %$% aba_genes_up
auxin <- readRDS(file = "./cleandata/core_auxin_genes.RDS")


#Gene set enrichment
list_hormones <- list(
  Ethylene = eth,
  Auxin = auxin,
  Cytokinin = cyto,
  ABA = aba,
  Flg22 = Res_regulons$flg22_core
)



df_gsea <- aggregate(log2FoldChange~gene_id,Res_sub,max) %>%
  dplyr::arrange(.data = .,-log2FoldChange)


df_gsea <- df_gsea %>% 
  dplyr::filter(.data = .,gene_id %in% ids_ppp) 

mvalues <- seq(0,4,0.5)
##Perform enrichment
Res_hyp_all <- NULL 
uni_genes  <- Res$gene_id %>% as.character %>% unique

for(mval in mvalues){
  
  focal_genes <- df_gsea %>%
    subset(log2FoldChange > mval) %$% gene_id %>% unique
  
  for(hormone in names(list_hormones)){
    
    mset <- list_hormones[[hormone]]
    hitInSample <- which(focal_genes %in% mset) %>% length
    hitInPop <- which(uni_genes %in% mset) %>% length
    sampleSize <- focal_genes %>% length
    failInPop <- (uni_genes  %>% length) - hitInPop
    mpval <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE);
    
    prop_mapped <- (hitInSample/sampleSize)*100
    
    Res_hyp_all <- data.frame(Direction = "+++",Hormone = hormone,
                              hitInSample = hitInSample,
                              hitInSampleProp = prop_mapped,
                              p = mpval,Cutoff = mval ) %>%
      rbind(Res_hyp_all,.)
    
    
  }
  
  

}


#Remove ABA and re order
Res_hyp_all_chosen <- Res_hyp_all %>% 
  dplyr::filter(.data = .,Hormone != "ABA") %>% droplevels %>%
  dplyr::mutate(.data = .,Hormone = Hormone %>% factor(levels = c("Ethylene","Auxin","Cytokinin","Flg22")))

Res_hyp_all_chosen$Direction <- Res_hyp_all_chosen$Direction %>%
  factor(levels = c("+++","---"))


Res_hyp_all_chosen$q <- Res_hyp_all_chosen$p %>% p.adjust(method = "fdr")
Res_hyp_all_chosen$Significance <- NA
Res_hyp_all_chosen$Significance[which(Res_hyp_all_chosen$q < 0.1)] <- "q < 0.1"

Res_hyp_all_chosen$Hormone <- Res_hyp_all_chosen$Hormone %>%
  factor(levels = c("Ethylene","Auxin","Cytokinin","Flg22"))

p1 <- ggplot(data = Res_hyp_all_chosen,aes(Cutoff,hitInSampleProp)) +
  geom_line(aes(group = Hormone)) +
  geom_point(aes(shape = Hormone,color = Significance),size = 10) +
  theme_ohchibi(size_panel_border = 0.3) +
  xlab(label = "Threshold log2 fold change") +
  ylab(label = "Proportion of genes in set") +
  scale_color_manual(values = "red",na.value = "#D9D9D9") +
  scale_shape_manual(values = 15:18) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank()
  )

#Create figure
px <- Res_hyp_all_chosen %>%
  subset(Cutoff ==3.0) %>%
  ggplot(data = .,aes(Hormone,hitInSampleProp)) +
  geom_bar(stat = "identity",aes(fill =Significance ),color = "transparent") +
  facet_grid(.~Direction,space = "free",scales = "free") +
  theme_ohchibi(size_panel_border = 0.3) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    legend.position = "top"
  ) +
  ylab(label = "Proportion of genes in set") +
  scale_fill_manual(values = "red",na.value = "black") 

### Create list of genes at cutoff 3
mval <- 3
focal_genes <- df_gsea %>%
  subset(log2FoldChange > mval) %$% gene_id %>% unique

#Retrieve information for the given genes
df_chosen <- df_gsea %>%
  dplyr::filter(.data = .,gene_id %in% focal_genes) %>%
  dplyr::arrange(.data = .,-log2FoldChange)

#Annotate genes
Map_genes <-  readRDS(file = "./cleandata/res_rnaseq_combinations_contrasts_intracol_intramutant_corlrdensity.RDS") %$% map_genes
Map_genes$IsTF <- NA
Map_genes$IsTF[Map_genes$GeneFamily %>%grep(pattern = "[Tt]ranscription [Ff]actor",ignore.case = F) ] <- "Yes"

#http://planttfdb.gao-lab.org/index.php?sp=Ath
df_tf <- read.table(file = "./rawdata/Ath_TF_list.txt",header = TRUE) 
Map_genes$IsTF[which(Map_genes$tair_locus %in% df_tf$Gene_ID)] <- "Yes"
Map_genes$FamilyTF_PTFDB <- match(Map_genes$tair_locus,df_tf$Gene_ID) %>%
  df_tf$Family[.]
Map_tf <- Map_genes %>% 
  subset(IsTF == "Yes") %>% droplevels


df_chosen$tair_symbol <- match(df_chosen$gene_id,Map_genes$ensembl_gene_id) %>%
  Map_genes$tair_symbol[.]

df_chosen$gene_biotype <- match(df_chosen$gene_id,Map_genes$ensembl_gene_id) %>%
  Map_genes$gene_biotype[.]

df_chosen$GeneFamily <- match(df_chosen$gene_id,Map_genes$ensembl_gene_id) %>%
  Map_genes$GeneFamily[.]

df_chosen$IsTF <- match(df_chosen$gene_id,Map_genes$ensembl_gene_id) %>%
  Map_genes$IsTF[.]

df_chosen$description <- match(df_chosen$gene_id,Map_genes$ensembl_gene_id) %>%
  Map_genes$description[.]

df_chosen$go_name <- match(df_chosen$gene_id,Map_genes$ensembl_gene_id) %>%
  Map_genes$go_name[.]

colnames(df_chosen)[2] <- "Avlog2FoldChange"

write.table(x = df_chosen,file = "./cleandata/table_rnaseq_3xgenes.tsv",
            append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)

### TFs 
mvalues <- seq(0,4,0.5)
##Perform enrichment
Res_hyp_all <- NULL 
uni_genes  <- Res$gene_id %>% as.character %>% unique


Map_genes <-  readRDS(file = "./cleandata/res_rnaseq_combinations_contrasts_intracol_intramutant_corlrdensity.RDS") %$% map_genes
Map_genes$IsTF <- NA
Map_genes$IsTF[Map_genes$GeneFamily %>%grep(pattern = "[Tt]ranscription [Ff]actor",ignore.case = F) ] <- "Yes"

#http://planttfdb.gao-lab.org/index.php?sp=Ath
df_tf <- read.table(file = "./rawdata/Ath_TF_list.txt",header = TRUE) 
Map_genes$IsTF[which(Map_genes$tair_locus %in% df_tf$Gene_ID)] <- "Yes"
Map_genes$FamilyTF_PTFDB <- match(Map_genes$tair_locus,df_tf$Gene_ID) %>%
  df_tf$Family[.]
Map_tf <- Map_genes %>% 
  subset(IsTF == "Yes") %>% droplevels



for(mval in mvalues){
  
  focal_genes <- df_gsea %>%
    subset(log2FoldChange > mval) %$% gene_id %>% unique
  
  for(hormone in names(list_hormones)){
    
    mset <- list_hormones[[hormone]]
    mset <- mset[which(mset %in% Map_tf$tair_locus)]
    hitInSample <- which(focal_genes %in% mset) %>% length
    hitInPop <- which(uni_genes %in% mset) %>% length
    sampleSize <- focal_genes %>% length
    failInPop <- (uni_genes  %>% length) - hitInPop
    mpval <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE);
    
    prop_mapped <- (hitInSample/sampleSize)*100
    
    Res_hyp_all <- data.frame(Direction = "+++",Hormone = hormone,
                              hitInSample = hitInSample,
                              hitInSampleProp = prop_mapped,
                              p = mpval,Cutoff = mval ) %>%
      rbind(Res_hyp_all,.)
    
    
  }
  
  
  
}



Res_hyp_all_chosen <- Res_hyp_all %>% 
  dplyr::filter(.data = .,Hormone != "ABA") %>% droplevels %>%
  dplyr::mutate(.data = .,Hormone = Hormone %>% factor(levels = c("Ethylene","Auxin","Cytokinin","Flg22")))

Res_hyp_all_chosen$Direction <- Res_hyp_all_chosen$Direction %>%
  factor(levels = c("+++","---"))


Res_hyp_all_chosen$q <- Res_hyp_all_chosen$p %>% p.adjust(method = "fdr")
Res_hyp_all_chosen$Significance <- NA
Res_hyp_all_chosen$Significance[which(Res_hyp_all_chosen$q < 0.1)] <- "q < 0.1"

Res_hyp_all_chosen$Hormone <- Res_hyp_all_chosen$Hormone %>%
  factor(levels = c("Ethylene","Auxin","Cytokinin","Flg22"))

p2 <- ggplot(data = Res_hyp_all_chosen,aes(Cutoff,hitInSampleProp)) +
  geom_line(aes(group = Hormone)) +
  geom_point(aes(shape = Hormone,color = Significance),size = 10) +
  theme_ohchibi(size_panel_border = 0.3) +
  xlab(label = "Threshold log2 fold change") +
  ylab(label = "Proportion of genes in set") +
  scale_color_manual(values = "red",na.value = "#D9D9D9") +
  scale_shape_manual(values = 15:18) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank()
  )


py <- Res_hyp_all_chosen %>%
  subset(Cutoff == 3.0) %>%
  ggplot(data = .,aes(Hormone,hitInSampleProp)) +
  geom_bar(stat = "identity",aes(fill =Significance ),color = "transparent") +
  facet_grid(.~Direction,space = "free",scales = "free") +
  theme_ohchibi(size_panel_border = 0.3) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    legend.position = "top"
  ) +
  ylab(label = "Proportion of genes in set") +
  scale_fill_manual(values = "red",na.value = "black") 

### Perform the inverted trend analysis
Res_hyp_all <- NULL 
uni_genes  <- Res$gene_id %>% as.character %>% unique

mvalues <- seq(-4,0,0.5) %>% rev

df_gsea <- aggregate(log2FoldChange~gene_id,Res_sub,min) %>%
  dplyr::arrange(.data = .,-log2FoldChange)

df_gsea <- df_gsea %>% 
  dplyr::filter(.data = .,gene_id %in% ids_nnn) 


for(mval in mvalues){
  
  focal_genes <- df_gsea %>%
    subset(log2FoldChange < mval) %$% gene_id %>% unique
  
  for(hormone in names(list_hormones)){
    
    mset <- list_hormones[[hormone]]
    hitInSample <- which(focal_genes %in% mset) %>% length
    hitInPop <- which(uni_genes %in% mset) %>% length
    sampleSize <- focal_genes %>% length
    failInPop <- (uni_genes  %>% length) - hitInPop
    mpval <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE);
    
    prop_mapped <- (hitInSample/sampleSize)*100
    
    Res_hyp_all <- data.frame(Direction = "---",Hormone = hormone,
                              hitInSample = hitInSample,
                              hitInSampleProp = prop_mapped,
                              p = mpval,Cutoff = mval ) %>%
      rbind(Res_hyp_all,.)
    
    
  }
  
  
  
}


#Remove ABA and re order
Res_hyp_all_chosen <- Res_hyp_all %>% 
  dplyr::filter(.data = .,Hormone != "ABA") %>% droplevels %>%
  dplyr::mutate(.data = .,Hormone = Hormone %>% factor(levels = c("Ethylene","Auxin","Cytokinin","Flg22")))

Res_hyp_all_chosen$Direction <- Res_hyp_all_chosen$Direction %>%
  factor(levels = c("+++","---"))


Res_hyp_all_chosen$q <- Res_hyp_all_chosen$p %>% p.adjust(method = "fdr")
Res_hyp_all_chosen$Significance <- NA
Res_hyp_all_chosen$Significance[which(Res_hyp_all_chosen$q < 0.1)] <- "q < 0.1"

Res_hyp_all_chosen$Hormone <- Res_hyp_all_chosen$Hormone %>%
  factor(levels = c("Ethylene","Auxin","Cytokinin","Flg22"))

p3 <- ggplot(data = Res_hyp_all_chosen,aes(Cutoff,hitInSampleProp)) +
  geom_line(aes(group = Hormone)) +
  geom_point(aes(shape = Hormone,color = Significance),size = 10) +
  theme_ohchibi(size_panel_border = 0.3) +
  xlab(label = "Threshold log2 fold change") +
  ylab(label = "Proportion of genes in set") +
  scale_color_manual(values = "red",na.value = "#D9D9D9") +
  scale_shape_manual(values = 15:18) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  scale_x_reverse()


### TFs inverted ####
Res_hyp_all <- NULL 
uni_genes  <- Res$gene_id %>% as.character %>% unique

mvalues <- seq(-4,0,0.5) %>% rev

df_gsea <- aggregate(log2FoldChange~gene_id,Res_sub,min) %>%
  dplyr::arrange(.data = .,-log2FoldChange)

df_gsea <- df_gsea %>% 
  dplyr::filter(.data = .,gene_id %in% ids_nnn) 


for(mval in mvalues){
  
  focal_genes <- df_gsea %>%
    subset(log2FoldChange < mval) %$% gene_id %>% unique
  
  for(hormone in names(list_hormones)){
    
    mset <- list_hormones[[hormone]]
    mset <- mset[which(mset %in% Map_tf$tair_locus)]
    hitInSample <- which(focal_genes %in% mset) %>% length
    hitInPop <- which(uni_genes %in% mset) %>% length
    sampleSize <- focal_genes %>% length
    failInPop <- (uni_genes  %>% length) - hitInPop
    mpval <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE);
    
    prop_mapped <- (hitInSample/sampleSize)*100
    
    Res_hyp_all <- data.frame(Direction = "---",Hormone = hormone,
                              hitInSample = hitInSample,
                              hitInSampleProp = prop_mapped,
                              p = mpval,Cutoff = mval ) %>%
      rbind(Res_hyp_all,.)
    
    
  }
  
  
  
}


#Remove ABA and re order
Res_hyp_all_chosen <- Res_hyp_all %>% 
  dplyr::filter(.data = .,Hormone != "ABA") %>% droplevels %>%
  dplyr::mutate(.data = .,Hormone = Hormone %>% factor(levels = c("Ethylene","Auxin","Cytokinin","Flg22")))

Res_hyp_all_chosen$Direction <- Res_hyp_all_chosen$Direction %>%
  factor(levels = c("+++","---"))


Res_hyp_all_chosen$q <- Res_hyp_all_chosen$p %>% p.adjust(method = "fdr")
Res_hyp_all_chosen$Significance <- NA
Res_hyp_all_chosen$Significance[which(Res_hyp_all_chosen$q < 0.1)] <- "q < 0.1"

Res_hyp_all_chosen$Hormone <- Res_hyp_all_chosen$Hormone %>%
  factor(levels = c("Ethylene","Auxin","Cytokinin","Flg22"))

p4 <- ggplot(data = Res_hyp_all_chosen,aes(Cutoff,hitInSampleProp)) +
  geom_line(aes(group = Hormone)) +
  geom_point(aes(shape = Hormone,color = Significance),size = 10) +
  theme_ohchibi(size_panel_border = 0.3) +
  xlab(label = "Threshold log2 fold change") +
  ylab(label = "Proportion of genes in set") +
  scale_color_manual(values = "red",na.value = "#D9D9D9") +
  scale_shape_manual(values = 15:18) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  scale_x_reverse() +
  scale_y_continuous(limits = c(0,1))

### Save graphs ###
oh.save.pdf(p = p1,outname = "rnaseq_curves_lfc_all_up.pdf",outdir = "./figures/",width = 8,height = 5)
oh.save.pdf(p = p2,outname = "rnaseq_curves_lfc_tf_up.pdf",outdir = "./figures/",width = 8,height = 5)
oh.save.pdf(p = p3,outname = "rnaseq_curves_lfc_all_down.pdf",outdir = "./figures/",width = 8,height = 5)
oh.save.pdf(p = p4,outname = "rnaseq_curves_lfc_tf_down.pdf",outdir = "./figures/",width = 8,height = 5)
oh.save.pdf(p = px,outname = "rnaseq_curves_lfc_all_up_30.pdf",outdir = "./figures/",width = 8,height = 5)
oh.save.pdf(p = py,outname = "rnaseq_curves_lfc_tf_up_30.pdf",outdir = "./figures/",width = 8,height = 5)


rm(list=ls())
gc()
