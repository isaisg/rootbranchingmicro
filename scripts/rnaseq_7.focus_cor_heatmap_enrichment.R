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

Res <- readRDS(file = "./cleandata/res_rnaseq_combinations_contrasts_intracol_intramutant_corlrdensity.RDS")

Map_genes <- Res$map_genes
Map_genes$IsTF <- NA
Map_genes$IsTF[Map_genes$GeneFamily %>%grep(pattern = "[Tt]ranscription [Ff]actor",ignore.case = F) ] <- "Yes"

#http://planttfdb.gao-lab.org/index.php?sp=Ath
df_tf <- read.table(file = "./rawdata/Ath_TF_list.txt",header = TRUE) 
Map_genes$IsTF[which(Map_genes$tair_locus %in% df_tf$Gene_ID)] <- "Yes"
Map_genes$FamilyTF_PTFDB <- match(Map_genes$tair_locus,df_tf$Gene_ID) %>%
  df_tf$Family[.]
Map_tf <- Map_genes %>% 
  subset(IsTF == "Yes") %>% droplevels


df_comb <- Res$df_comb_gec

ids_ppp <- df_comb$EnrichmentCode %>%
  grep(pattern = "001001001")  %>%
  df_comb[.,] %$% gene_id %>% as.character %>% unique

ids_nnn <- df_comb$EnrichmentCode %>%
  grep(pattern = "100100100")  %>%
  df_comb[.,] %$% gene_id %>% as.character %>% unique


#Draw heatmap
mids <- c(ids_ppp,ids_nnn) %>% unique


tf_ids_ppp <- which(ids_ppp %in% Map_tf$tair_locus) %>% ids_ppp[.] %>%
  match(Map_tf$tair_locus) %>%
  Map_tf[.,] %>% droplevels



### Read expression object to create heatmap ### 
Tab <- readRDS(file = "./cleandata/dat_tab_av_rnaseq.RDS")


###
Tab_sub <- which(rownames(Tab) %in% mids) %>% 
  Tab[.,]


res_heatmap <- chibi.heatmap(Tab = Tab_sub,
                             dist_method_rows = "pearson",dist_method_cols = "euclidean",
                             width_border_tile = 0.5,height_border_tile = 0.5,
                             size_border_tile = 0.05,
                             hclust_method_rows = "ward.D",
                             hclust_method_cols = "ward.D",
                             k_rows = 2,k_cols = 9,
                             panel_spacing = 0.1,
                             range_fill_heatmap = c(-0.5,0.5),
                             palette_border  = c("#51F5C9","#F55C7A"),
                             mtheme = theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1)),legend_position = "bottom")

melted_heat <- res_heatmap$melted

#Read regulon information
Res_regulons <- readRDS(file = "./cleandata/regulons_rnaseq.RDS")
eth <- Res_regulons$ethylene
cyto <- Res_regulons$cytokinin

### Read Lavenus 
gc()
df_lavenus <- read.xlsx(file = "./rawdata/auxin_markers/data/plcell_v27_5_1368_s1/plcell_v27_5_1368_s1/tpc132993_SupplementalDS1.xlsx",
                        sheetIndex = 2)

a <- df_lavenus %>% subset(a > 1.5) %$% AGI
b <- df_lavenus %>% subset(b > 1.5) %$% AGI
ab <- df_lavenus %>% subset(a...b > 1.5) %$% AGI
lavenus <- c(a,b,ab) %>% unique

aba <- readRDS("./rawdata/aba_robust_genes.RDS") %$% aba_genes_up


#Gene set enrichment
list_hormones <- list(
  Ethylene = eth,
  LR = lavenus,
  Cytokinin = cyto,
  ABA = aba,
  Flg22 = Res_regulons$flg22_core
)

#Create bars for hormones
melted_heat$Ethylene <- "No"
melted_heat$Ethylene[which(melted_heat$IdRows %in% list_hormones$Ethylene)] <- "Yes"

melted_heat$Auxin <- "No"
melted_heat$Auxin[which(melted_heat$IdRows %in%  list_hormones$Auxin)] <- "Yes"

melted_heat$LR <- "No"
melted_heat$LR[which(melted_heat$IdRows %in%  list_hormones$LR)] <- "Yes"


melted_heat$Cytokinin <- "No"
melted_heat$Cytokinin[which(melted_heat$IdRows %in%  list_hormones$Cytokinin)] <- "Yes"

melted_heat$ABA <- "No"
melted_heat$ABA[which(melted_heat$IdRows %in%  list_hormones$ABA)] <- "Yes"


melted_heat$Flg22 <- "No"
melted_heat$Flg22[which(melted_heat$IdRows %in%  list_hormones$Flg22)] <- "Yes"

melted_heat$Bar <- "Bar"

df_bar <- melted_heat[,c("IdRows","ClusterRows","Ethylene","Auxin","LR","Cytokinin","ABA","Flg22","Bar")] %>% unique

### Create the bars
paleta_cor <- c("white","black") 
names(paleta_cor) <- c("No","Yes")

p_eth <- ggplot(data = df_bar,aes(Bar,IdRows)) +
  geom_tile(aes(fill = Ethylene),color = "transparent") +
  facet_grid(ClusterRows~.,space = "free",scales = "free") +
  theme_ohchibi(size_panel_border = 0.5) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text.y = element_blank(),
    strip.text.x = element_text(size = 8),
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.spacing.y = unit(0.2, "lines"),
    axis.title.x = element_text(size = 8,angle = 90,vjust = 0.5,hjust = 1)
    
  ) +
  scale_fill_manual(values = paleta_cor) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  xlab("Ethylene")



p_auxin <- ggplot(data = df_bar,aes(Bar,IdRows)) +
  geom_tile(aes(fill = Auxin),color = "transparent") +
  facet_grid(ClusterRows~.,space = "free",scales = "free") +
  theme_ohchibi(size_panel_border = 0.5) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text.y = element_blank(),
    strip.text.x = element_text(size = 8),
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.spacing.y = unit(0.2, "lines"),
    axis.title.x = element_text(size = 8,angle = 90,vjust = 0.5,hjust = 1)
    
  ) +
  scale_fill_manual(values = paleta_cor) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  xlab("Auxin")

p_lr <- ggplot(data = df_bar,aes(Bar,IdRows)) +
  geom_tile(aes(fill = LR),color = "transparent") +
  facet_grid(ClusterRows~.,space = "free",scales = "free") +
  theme_ohchibi(size_panel_border = 0.5) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text.y = element_blank(),
    strip.text.x = element_text(size = 8),
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.spacing.y = unit(0.2, "lines"),
    axis.title.x = element_text(size = 8,angle = 90,vjust = 0.5,hjust = 1)
    
  ) +
  scale_fill_manual(values = paleta_cor) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  xlab("LR")

p_cyt <- ggplot(data = df_bar,aes(Bar,IdRows)) +
  geom_tile(aes(fill = Cytokinin),color = "transparent") +
  facet_grid(ClusterRows~.,space = "free",scales = "free") +
  theme_ohchibi(size_panel_border = 0.5) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text.y = element_blank(),
    strip.text.x = element_text(size = 8),
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.spacing.y = unit(0.2, "lines"),
    axis.title.x = element_text(size = 8,angle = 90,vjust = 0.5,hjust = 1)
    
  ) +
  scale_fill_manual(values = paleta_cor) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  xlab("Cytokinin")


#### Create composition
p_tree_rows <- res_heatmap$p_tree_rows
p_tree_cols <- res_heatmap$p_treee_cols
p_blank <- ggplot() + theme_void()

#Manually recreate heatmap
p_heatmap <- res_heatmap$p_heatmap +
  scale_fill_gradient2(low = "#021567",
                       mid = "white",high = "#FBDE28",midpoint = 0,limits = c(-0.5,0.5),oob = squish)



composition <- egg::ggarrange(p_blank,p_tree_cols,p_blank,p_blank,p_blank,
                              p_tree_rows,p_heatmap + theme(legend.position = "none"),p_eth,p_lr,p_cyt,
                              nrow = 2,ncol = 5,
                              heights = c(0.2,1),widths = c(0.2,1,0.01,0.01,0.01),draw = F)

dev.off()


oh.save.pdf(p = composition,outname = "rnaseq_heatmap_specific.pdf",
            outdir = "./figures/",width = 18,height = 12)
dev.off()

##### Compare clusters ########
mlist_go <- list(
  CR1 = res_heatmap$df_clust_rows  %>% subset(ClusterRows == "CR1") %$%  IdRows %>% as.character,
  CR2 = res_heatmap$df_clust_rows  %>% subset(ClusterRows == "CR2") %$%  IdRows %>% as.character
)

cg <- compareCluster(
  geneCluster=mlist_go,
  fun="enrichGO",
  keyType       = "TAIR",
  OrgDb         = org.At.tair.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.1)


p_cnet <- cnetplot(x = cg,
                   showCategory = 100, 
                   layout = "kk", 
                   circular = FALSE,
                   node_label = "category",
                   cex_label_category = 1.5,
                   cex_category = 4,
) +
  theme(legend.position = "top")



oh.save.pdf(p = p_cnet,outname = "rnaseq_cnet_specific.pdf",
            outdir = "./figures/",width = 12,height = 10)

dev.off()

### Perform hypergeometric comparisons ###
Res_hyp_all <- NULL 

for(hormone in names(list_hormones)){
  
  mset <- list_hormones[[hormone]]
  hitInSample <- which(ids_ppp %in% mset) %>% length
  hitInPop <- which(rownames(Tab) %in% mset) %>% length
  sampleSize <- ids_ppp %>% length
  failInPop <- (rownames(Tab)  %>% length) - hitInPop
  mpval <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE);
  
  prop_mapped <- (hitInSample/sampleSize)*100
  
  Res_hyp_all <- data.frame(Direction = "+++",Hormone = hormone,
                            hitInSample = hitInSample,
                            hitInSampleProp = prop_mapped,
                            p = mpval ) %>%
    rbind(Res_hyp_all,.)
  
  
}


#Inverted
for(hormone in names(list_hormones)){
  
  mset <- list_hormones[[hormone]]
  hitInSample <- which(ids_nnn %in% mset) %>% length
  hitInPop <- which(rownames(Tab) %in% mset) %>% length
  sampleSize <- ids_nnn %>% length
  failInPop <- (rownames(Tab)  %>% length) - hitInPop
  mpval <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE);
  
  prop_mapped <- (hitInSample/sampleSize)*100
  
  
  Res_hyp_all <- data.frame(Direction = "---",
                            Hormone = hormone,
                            hitInSample = hitInSample,
                            hitInSampleProp = prop_mapped,
                            p = mpval ) %>%
    rbind(Res_hyp_all,.)
  
  
}



Res_hyp_tf <- NULL
for(hormone in names(list_hormones)){
  
  mset <- list_hormones[[hormone]]
  mset <- mset[which(mset %in% Map_tf$tair_locus)]
  hitInSample <- which(ids_ppp %in% mset) %>% length
  hitInPop <- which(rownames(Tab) %in% mset) %>% length
  sampleSize <- ids_ppp %>% length
  failInPop <- (rownames(Tab)  %>% length) - hitInPop
  mpval <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE);
  
  prop_mapped <- (hitInSample/sampleSize)*100
  
  
  Res_hyp_tf <- data.frame(Direction = "+++",
                           Hormone = hormone,
                           hitInSample = hitInSample,
                           hitInSampleProp = prop_mapped,
                           p = mpval ) %>%
    rbind(Res_hyp_tf,.)
  
  
}

for(hormone in names(list_hormones)){
  
  mset <- list_hormones[[hormone]]
  mset <- mset[which(mset %in% Map_tf$tair_locus)]
  hitInSample <- which(ids_nnn %in% mset) %>% length
  hitInPop <- which(rownames(Tab) %in% mset) %>% length
  sampleSize <- ids_nnn %>% length
  failInPop <- (rownames(Tab)  %>% length) - hitInPop
  mpval <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE);
  
  prop_mapped <- (hitInSample/sampleSize)*100
  
  Res_hyp_tf <- data.frame(Direction = "---",
                           Hormone = hormone,
                           hitInSample = hitInSample,
                           hitInSampleProp = prop_mapped,
                           p = mpval ) %>%
    rbind(Res_hyp_tf,.)
  
  
}

### Perform plottting of the representation
Res_hyp_all$Direction <- Res_hyp_all$Direction %>%
  factor(levels = c("+++","---"))
Res_hyp_all$q <- Res_hyp_all$p %>% p.adjust(method = "bonferroni")

Res_hyp_all$Significance <- NA
Res_hyp_all$Significance[which(Res_hyp_all$q < 0.05)] <- "q < 0.05"

Res_hyp_tf$Direction <- Res_hyp_tf$Direction %>%
  factor(levels = c("+++","---"))
Res_hyp_tf$q <- Res_hyp_tf$p %>% p.adjust(method = "bonferroni")

Res_hyp_tf$Significance <- NA
Res_hyp_tf$Significance[which(Res_hyp_tf$q < 0.05)] <- "q < 0.05"

Res_hyp_all$Direction <- Res_hyp_all$Direction %>%
  as.character %>%
  gsub(pattern = '\\+\\+\\+',replacement = "Upregulated") %>%
  gsub(pattern = '\\-\\-\\-',replacement = "Downregulated") %>%
  factor(levels = c("Upregulated","Downregulated"))

Res_hyp_tf$Direction <- Res_hyp_tf$Direction %>%
  as.character %>%
  gsub(pattern = '\\+\\+\\+',replacement = "Upregulated") %>%
  gsub(pattern = '\\-\\-\\-',replacement = "Downregulated") %>%
  factor(levels = c("Upregulated","Downregulated"))


#Remove ABA and re order
Res_hyp_all <- Res_hyp_all %>% 
  dplyr::filter(.data = .,Hormone != "ABA") %>% droplevels %>%
  dplyr::mutate(.data = .,Hormone = Hormone %>% factor(levels = c("LR","Ethylene","Cytokinin","Flg22")))

Res_hyp_tf <- Res_hyp_tf %>% 
  dplyr::filter(.data = .,Hormone != "ABA") %>% droplevels %>%
  dplyr::mutate(.data = .,Hormone = Hormone %>% factor(levels = c("LR","Ethylene","Cytokinin","Flg22")))


p1 <- ggplot(data = Res_hyp_all,aes(Hormone,hitInSampleProp)) +
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
  scale_fill_manual(values = "red",na.value = "black") +
  scale_y_continuous(expand = c(0,0),limits = c(0,16))


p2 <- ggplot(data = Res_hyp_tf,aes(Hormone,hitInSampleProp)) +
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
  scale_fill_manual(values = "red",na.value = "black") +
  scale_y_continuous(expand = c(0,0))


#Save figures ###
oh.save.pdf(p = p1,outname = "rnaseq_hyp_all.pdf",outdir = "./figures/",width = 8,height = 5)
oh.save.pdf(p = p2,outname = "rnaseq_hyp_tf.pdf",outdir = "./figures/",width = 8,height = 5)
