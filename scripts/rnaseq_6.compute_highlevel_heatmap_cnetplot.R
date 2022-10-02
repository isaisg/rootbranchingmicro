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

df_comb <- Res$df_comb_gec

#Determine the ones in col0

ids_ppp <- df_comb$EnrichmentCode %>%
  grep(pattern = "001001001")  %>%
  df_comb[.,] %$% gene_id %>% as.character %>% unique


ids_pp0 <- df_comb$EnrichmentCode %>%
  grep(pattern = "001001010")  %>%
  df_comb[.,] %$% gene_id %>% as.character %>% unique

ids_nnn <- df_comb$EnrichmentCode %>%
  grep(pattern = "100100100")  %>%
  df_comb[.,] %$% gene_id %>% as.character %>% unique

ids_nn0 <- df_comb$EnrichmentCode %>%
  grep(pattern = "100100010")  %>%
  df_comb[.,] %$% gene_id %>% as.character %>% unique

#Draw heatmap
mids <- c(ids_ppp,ids_nnn,ids_pp0,ids_nn0) %>% unique

### Read expression object to create heatmap ### 
Tab <- readRDS(file = "./cleandata/dat_tab_av_rnaseq.RDS")

###
Tab_sub <- which(rownames(Tab) %in% mids) %>% 
  Tab[.,]


#Create df_border
mpairs <- df_comb$BacteriaGenotype %>% unique
df_border <- NULL

for(mpair in mpairs){
  
  x <- mpair %>%
    gsub(pattern = "_.*",replacement = "")
  
  y <- mpair %>% gsub(pattern = paste0(x,"_"),replacement = "")  
  idcol <- paste0(y,"_",x)
  
  df_chosen_sub <- df_comb %>% 
    subset(BacteriaGenotype == mpair) %>% droplevels
  
  if(nrow(df_chosen_sub) > 0){
    
    df_border <- data.frame(IdRows = df_chosen_sub$gene_id,
                            IdCols = idcol,
                            Significance = df_chosen_sub$EnrichmentCode) %>%
      rbind(df_border,.)
    
  }
  
}

df_border$Significance <- df_border$Significance %>%
  gsub(pattern = "001001.*",replacement = "++") %>%
  gsub(pattern = "100100.*",replacement = "--")  %>%
  factor(levels = c("++","--") %>% rev)


res_heatmap <- chibi.heatmap(Tab = Tab_sub,
                             df_border = df_border,
                             dist_method_rows = "pearson",dist_method_cols = "euclidean",
                             width_border_tile = 0.5,height_border_tile = 0.5,
                             size_border_tile = 0.05,
                             hclust_method_rows = "ward.D",
                             hclust_method_cols = "ward.D",
                             k_rows = 10,k_cols = 9,
                             panel_spacing = 0.1,
                             range_fill_heatmap = c(-0.5,0.5),
                             palette_border  = c("#51F5C9","#F55C7A"),
                             mtheme = theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1)),legend_position = "bottom")

#res_heatmap$heatmap
melted_heat <- res_heatmap$melted

melted_heat$ClusterRows <- melted_heat$ClusterRows %>% 
  as.character %>%
gsub(pattern = "CR4|CR8|CR3",replacement = "CRX") %>%
gsub(pattern = "CR5|CR10|CR2",replacement = "CRY") %>%
gsub(pattern = "CR9|CR6",replacement = "CRZ") %>%
  gsub(pattern = "CR7|CR1",replacement = "CRA")  %>%
  gsub(pattern = "CRX",replacement = "CR1")  %>%
  gsub(pattern = "CRY",replacement = "CR2")  %>%
  gsub(pattern = "CRZ",replacement = "CR3")  %>%
  gsub(pattern = "CRA",replacement = "CR4")   %>%
  factor
  

#Read correlation info
df_cor <- readRDS(file = "./cleandata/res_cor_rnaseq_lrdensity.RDS") %$% Res_cor
df_cor$p_spearman[which(is.na(df_cor$p_spearman))] <- 1

#Use Rho
down_lr <- df_cor %>% subset(p_spearman < 0.05 & rho < 0)  %$% Gene
up_lr <- df_cor %>% subset(p_spearman < 0.05 & rho > 0)  %$% Gene
zero_lr <- df_cor %>% subset(p_spearman >= 0.05)  %$% Gene


## Create map genes
Map_genes <- Res$map_genes
Map_genes$IsTF <- "No"
Map_genes$IsTF[Map_genes$GeneFamily %>% as.character %>%grep(pattern = "[Tt]ranscription [Ff]actor",ignore.case = F)] <- "Yes"

Map_genes <- Map_genes %>% 
  dplyr::relocate(.data = .,"ensembl_gene_id","tair_locus",
                  "tair_symbol","gene_biotype","GeneFamily","IsTF")  %>%
  dplyr::select(.data = .,-go_definition)


cor_pos <- up_lr
cor_neg <- down_lr
cor_null <- zero_lr



melted_heat$Cor <- "No"
melted_heat$Cor[which(melted_heat$IdRows %in% cor_pos)] <- "Positive"
melted_heat$Cor[which(melted_heat$IdRows %in% cor_neg)] <- "Negative"

melted_heat$Bar <- "Bar"

df_bar <- melted_heat[,c("IdRows","ClusterRows","Cor","Bar")] %>% unique

paleta_cor <- c("#51F5C9","white","#F55C7A") 
names(paleta_cor) <- c("Negative","No","Positive")

p_cor <- ggplot(data = df_bar,aes(Bar,IdRows)) +
  geom_tile(aes(fill = Cor),color = "transparent") +
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
  xlab("Cor LR")

#paleta_bar <- c("#E7D202FF","#5E81ACFF","#BF616AFF","#7D5329FF")
paleta_bar <- c("#8AF7FF","#EB9D67","#27A149","#A36080")
names(paleta_bar) <- c("CR1","CR2","CR3","CR4")

p_bar <- ggplot(data = df_bar,aes(Bar,IdRows)) +
  geom_tile(aes(fill = ClusterRows),color = "transparent") +
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
  scale_fill_manual(values = paleta_bar) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  xlab("Cluster")



#### Assemble panel #####

p_tree_rows <- res_heatmap$p_tree_rows
p_tree_cols <- res_heatmap$p_treee_cols
p_blank <- ggplot() + theme_void()

#Manually recreate heatmap

p_heatmap <- res_heatmap$p_heatmap +
  scale_fill_gradient2(low = "#021567",
                       mid = "white",high = "#FBDE28",midpoint = 0,limits = c(-0.5,0.5),oob = squish)


p_a <- ggplot(data = melted_heat,aes(IdCols,IdRows)) +
  geom_raster(aes(fill = value),color = "transparent") +
  geom_tile(aes(color = Border),fill = "transparent",width = 0.5,height = 0.5,size = 0.05) +
  scale_fill_gradient2(low = "#021567",
                       mid = "white",high = "#FBDE28",midpoint = 0,limits = c(-0.5,0.5),oob = squish) +
  scale_color_manual(values = c("#51F5C9","#F55C7A"),na.value = "transparent") +
  facet_grid(ClusterRows~ClusterCols,space = "free",scales = "free") +
  theme_ohchibi(size_panel_border = 0.3) +
  theme(
    strip.text.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size = 8),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8,angle = 45,vjust = 1,hjust = 1),
    axis.title.x = element_text(size = 8),
    axis.title.y= element_text(size = 8),
    strip.text.y = element_blank(),
    panel.spacing.x = unit(0.1,"lines"),
    panel.spacing.y = unit(0.1,"lines")
    
    
  ) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  xlab(label = "Treatment") +
  ylab(label = "Gene")

composition <- egg::ggarrange(p_blank,p_tree_cols,p_blank,p_blank,
                              p_tree_rows,p_a,p_bar,p_cor,
                              nrow = 2,ncol = 4,
                              heights = c(0.2,1),widths = c(0.2,1,0.01,0.01),draw = F)

dev.off()
oh.save.pdf(p = composition,outname = "rnaseq_heatmap_corlr_overall.pdf",
            outdir = "./figures/",width = 18,height = 12)


##### Compare clusters ########
mlist_go <- list(
  CR1 = df_bar  %>% subset(ClusterRows == "CR1") %$%  IdRows %>% as.character,
  CR2 = df_bar  %>% subset(ClusterRows == "CR2") %$%  IdRows %>% as.character,
  CR3 = df_bar  %>% subset(ClusterRows == "CR3") %$%  IdRows %>% as.character,
  CR4 = df_bar  %>% subset(ClusterRows == "CR4") %$%  IdRows %>% as.character
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

bp2 <- clusterProfiler::simplify(cg, cutoff = 0.8, 
                                 by = "p.adjust", select_fun = min)

p_cnet <- cnetplot(x = bp2,
         showCategory = 20, 
         layout = "kk", 
         circular = FALSE,
         node_label = "category",
         cex_label_category = 1.5,
         cex_category = 4,
         ) +
  scale_fill_manual(values = paleta_bar)  +
  theme(legend.position = "top")

oh.save.pdf(p = p_cnet,outname = "rnaseq_cnet_corlr_overall.pdf",
            outdir = "./figures/",width = 12,height = 10)

### Create R distribution stastiscic
df_cor$Type <- NA
df_cor$Type[which(df_cor$Gene %in% ids_ppp)] <- "+++"
df_cor$Type[which(df_cor$Gene %in% ids_pp0)] <- "++0"

paleta_cor_pos <- c("#F55C7A","white")


px <- df_cor %>%
  dplyr::filter(.data = .,!is.na(Type)) %>%
  ggplot(data = .,aes(rho)) +
  geom_density(aes(fill = Type,linetype = Type),alpha = 0.7,color = "black",size = 0.3) +
  theme_ohchibi(size_panel_border = 0.3,size_axis_text.x = 8,size_axis_text.y = 8,size_axis_title.y = 8,size_axis_title.x = 8,size_title_text = 10) +
  ylab(label = "Density") +
  xlab(label = "Estimate vs NB") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank() ,
    legend.position = "none"
  )  +
  scale_fill_manual(values = paleta_cor_pos) +
  scale_color_manual(values = paleta_cor_pos) +
  scale_x_continuous(limits = c(-1,1)) +
  scale_y_continuous(limits = c(0,7.5),expand = c(0,0))

df_cor$Type <- NA
df_cor$Type[which(df_cor$Gene %in% ids_nnn)] <- "---"
df_cor$Type[which(df_cor$Gene %in% ids_nn0)] <- "--0"

paleta_cor_neg <- c("#51F5C9","white")


py <- df_cor %>%
  dplyr::filter(.data = .,!is.na(Type)) %>%
  ggplot(data = .,aes(rho)) +
  geom_density(aes(fill = Type,linetype = Type),alpha = 0.7,color = "black",size = 0.3) +
  theme_ohchibi(size_panel_border = 0.3,size_axis_text.x = 8,size_axis_text.y = 8,size_axis_title.y = 8,size_axis_title.x = 8,size_title_text = 10) +
  ylab(label = "Density") +
  xlab(label = "Estimate vs NB") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank() ,
    legend.position = "none"
  )  +
  scale_fill_manual(values = paleta_cor_neg) +
  scale_color_manual(values = paleta_cor_neg) +
  scale_x_continuous(limits = c(-1,1)) +
  scale_y_continuous(limits = c(0,7.5),expand = c(0,0))



composition <- egg::ggarrange(px,py,nrow = 1)

oh.save.pdf(p = composition,outname = "rnaseq_corlr_density.pdf",
            outdir = "./figures/",width = 18,height = 12)



### Select solely usable genes
rm(list=ls())
gc()
dev.off()
