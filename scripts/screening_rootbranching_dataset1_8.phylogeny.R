library(ohchibi)
library(xlsx)
library(ggtree)

set.seed(130816)

Res_all <- readRDS(file = "./cleandata/res_cv_prop_sig_rootbranching_screening1.RDS")

#Plot counting number of modified features
df <- Res_all$df_strain_numsig_vs_nb
df$NumberFeatureSigDiff <- df$NumberFeatureSigDiff %>% factor



df_strain <- Res_all$df_strain_numsig_vs_nb 

df_iaa_acc_genome <- readRDS(file = "./cleandata/dat_rootbranching_cog_ko.RDS") %$% df_iaa_acc_genome

selected_strains <- read.xlsx(file = "./rawdata/rawdata_selaginella_screening.xlsx",sheetIndex = 1) %>%
  dplyr::rename(.data = .,Strain = Strains,
                NumberPlant = explant.nb,
                RhizophoreNumber = nb.of.Rhizophore) %>%
  dplyr::select(.data = .,c(Strain,NumberPlant,RhizophoreNumber,Bifurcation)) %>%
  dplyr::mutate(.data = .,Strain = Strain %>% gsub(pattern = "_",replacement = "") %>%
                  gsub(pattern = "MF",replacement = "RMF") %>%
                  gsub(pattern = "^CL",replacement = "RCL")) %$% Strain %>%
  grep(pattern = "NB",invert = T,value = T) %>% unique


### All 
merged <- Res_all$df_sig_vs_nb_raw 

df_est <- dcast(data = merged,formula = Strain~Feature,value.var = "estimate",fill = 0) 


#Read phylogenetic tree
tree <- readRDS(file = "./cleandata/dat_rootbranching_screening1.RDS") %$% tree

Map_strains <- readRDS(file = "./cleandata/dat_rootbranching_screening1.RDS") %$% Map_strains

paleta <- paletteer_d("ggthemes::Tableau_20",n = 19) %>% as.character
names(paleta) <- Map_strains$Order %>% table %>% names


p_tree <- ggtree(tree,ladderize = FALSE,size = 0.5) + 
  scale_y_reverse(expand = c(0.001, 0.001)) 

p_tree <- p_tree %<+% Map_strains +
   aes(color = Order) +
  geom_tippoint(aes(color = Order)) +
  scale_color_manual(values = paleta) +
  theme(
    legend.position = "none"
  )  +
  geom_tiplab(align = TRUE,size = 0)

#Determine order of tree
d <- p_tree$data %>% subset(isTip)
mtips <- with(d, label[order(y, decreasing=T)])


df_border <- merged[,c("Strain","Feature","Significance")] 
df_border$Significance <- df_border$Significance %>% gsub(pattern = "NS",replacement = NA) 


#Estimate phylogenetic signal
df_est <- df_est %>% 
  dplyr::filter(.data = .,Strain %in% mtips) %>%
  droplevels

#Melt structure
melted <- melt(data = df_est,id.vars = c("Strain")) %>%
  dplyr::rename(.data = .,Feature = variable) %>% 
  merge(df_border, by = c("Strain","Feature")) 


melted$NumberFeatureSigDiff <- match(melted$Strain,df_strain$Strain) %>%
  df_strain$NumberFeatureSigDiff[.]

melted$Sum_IAA_Bio <- match(melted$Strain,df_iaa_acc_genome$Strain) %>%
  df_iaa_acc_genome$Sum_IAA_Bio[.]

melted$Sum_IAA_Deg <- match(melted$Strain,df_iaa_acc_genome$Strain) %>%
  df_iaa_acc_genome$Sum_IAA_Deg[.]

melted$SelectedStrains <- NA
melted$SelectedStrains[which(melted$Strain %in% selected_strains )] <- "Yes"


melted <- melted %>%
  dplyr::filter(.data = .,Strain %in% mtips) %>%
  droplevels %>%
  dplyr::mutate(.data = .,Strain = Strain %>% factor(levels = mtips))


#Append information about order and clustering 
df_feat_cv_phylosig <- readRDS(file = "./cleandata/res_cv_prop_sig_rootbranching_screening1.RDS") %$% df_feat_cv_phylosig

melted <- merge(melted,df_feat_cv_phylosig, by = "Feature") 

melted$Feature <- melted$Feature %>% factor(levels = df_feat_cv_phylosig$Feature %>% levels %>% rev)
melted$Cluster <- melted$Cluster %>% factor(levels = df_feat_cv_phylosig$Cluster %>% levels)

#Read order from heatmap



p_heatmap <- ggplot(data = melted,aes(Feature,Strain)) +
  geom_raster(aes(fill = value),color = "transparent") +
  geom_tile(aes(color = Significance),fill = "transparent",width = 0.05,height = 0.05,size = 0.2) +
  #scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",limits = c(-1.5,1.5),oob = squish) +
  scale_fill_gradient2(low = "#021567",
                       mid = "white",high = "#FBDE28",midpoint = 0,limits = c(-1.5,1.5),oob = squish) +
  scale_color_manual(values = c("black"),na.value = "transparent") +
  facet_grid(.~Cluster,space = "free",scales = "free") +
  theme_ohchibi(size_panel_border = 0.3) +
  theme(
    axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
    strip.text.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +
  #scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0))


df_feature <- melted[,c("Feature","Phylosig_lambda_estimate","Phylosig_pvalue_estimate","Cluster")] %>%
  unique 
df_feature$q <- df_feature$Phylosig_pvalue_estimate %>% p.adjust(method = "fdr")
df_feature$Significance <- NA
df_feature$Significance[which(df_feature$q < 0.05)] <- "q < 0.05"

p_top <- ggplot(data = df_feature,aes(Feature,Phylosig_lambda_estimate)) +
  geom_bar(stat = "identity",aes(fill = Significance)) +
  facet_grid(.~Cluster,space = "free",scales = "free") +
  theme_ohchibi(size_panel_border = 0.3) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8)
  ) +
  #scale_x_discrete(expand = c(0,0))  +
  scale_y_continuous(breaks = seq(0,1,0.1),expand = c(0,0)) +
  ylab(label = "Phylogenetic signal") +
  scale_fill_manual(values = c("black"),na.value = "#D9D9D9")

p_blank <- ggplot() + theme_void()

#Addd bar with operon and selected strains
df_bar <- melted[,c("Strain",
                    "NumberFeatureSigDiff","Sum_IAA_Bio","Sum_IAA_Deg","SelectedStrains")] %>% unique %>%
  dplyr::mutate(.data = .,Bar = "Bar")
df_bar$Sum_IAA_Bio[df_bar$Strain %>% grep(pattern = "L404")] <- 0
df_bar$Sum_IAA_Deg[df_bar$Strain %>% grep(pattern = "L404")] <- 0

df_bar$Sum_IAA_Bio <- df_bar$Sum_IAA_Bio %>% factor
df_bar$Sum_IAA_Deg <- df_bar$Sum_IAA_Deg %>% factor

p_bio <- ggplot(data = df_bar,aes(Bar,Strain)) +
  geom_tile(aes(fill =Sum_IAA_Bio ),color = "transparent") +
  scale_x_discrete(expand = c(0,0)) +
  theme_ohchibi(size_panel_border = 0.3) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.y = element_text(angle = 0,size = 6),
    axis.title.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5,size = 8)
    
  ) +
  xlab(label = "IAA\nbiosynthesis") +
  scale_fill_manual(values = c("#F1DDBF","#84BFB9","#84BFB9","#84BFB9","#84BFB9"))

p_deg <- ggplot(data = df_bar,aes(Bar,Strain)) +
  geom_tile(aes(fill =Sum_IAA_Deg ),color = "transparent") +
  scale_x_discrete(expand = c(0,0)) +
  theme_ohchibi(size_panel_border = 0.3) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.y = element_text(angle = 0,size = 6),
    axis.title.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5,size = 8)
    
  ) +
  xlab(label = "IAA\ndegradation") +
  scale_fill_manual(values = c("#F1DDBF","#84BFB9","#84BFB9","#84BFB9","#84BFB9"))

#Plot the selected

p_chosen <- ggplot(data = df_bar,aes(Bar,Strain)) +
  geom_tile(aes(fill =SelectedStrains ),color = "transparent") +
  scale_x_discrete(expand = c(0,0)) +
  theme_ohchibi(size_panel_border = 0.3) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.y = element_text(angle = 0,size = 6),
    legend.position = "top",
    axis.title.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5,size = 8)
  ) +
  xlab(label = "Selected\nstrains") +
  scale_fill_manual(values = c("#84BFB9"),na.value = "#F1DDBF")




p_num <- ggplot(data = df_bar,aes(Bar,Strain)) +
  geom_tile(aes(fill =NumberFeatureSigDiff ),color = "transparent") +
  scale_x_discrete(expand = c(0,0)) +
  theme_ohchibi(size_panel_border = 0.3) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.y = element_text(angle = 0,size = 6),
    legend.position = "top",
    axis.title.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5,size = 8)
    
  ) +
  xlab(label = "Number of features changed") +
  scale_fill_paletteer_c("pals::kovesi.rainbow_bgyrm_35_85_c71")


p_blank <- ggplot() + theme_void()

dev.off()

#Create bar with cluster belonging
Res_heatmap <- readRDS(file = "./cleandata/res_rootbranching_screening1_heatmap_melted.RDS")
df_bar$ClusterRows <- match(df_bar$Strain,Res_heatmap$IdRows) %>%
  Res_heatmap$ClusterRows[.]

df_bar$ClusterRows <- df_bar$ClusterRows %>% factor(levels = paste0("CR",1:24))

paleta <- c(paletteer_d("ggsci::category20_d3") %>% as.character,"black","firebrick","purple","yellow")

p_cluster <- ggplot(data = df_bar,aes(Bar,Strain)) +
  geom_tile(aes(fill =ClusterRows ),color = "transparent") +
  scale_x_discrete(expand = c(0,0)) +
  theme_ohchibi(size_panel_border = 0.3) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.y = element_text(angle = 0,size = 6),
    legend.position = "top",
    axis.title.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5,size = 8)
    
  ) +
  xlab(label = "Cluster")  +
  scale_fill_manual(values = paleta)

p_cluster_print <- ggplot(data = df_bar,aes(Bar,Strain)) +
  geom_tile(aes(fill =ClusterRows ),color = "transparent") +
  scale_x_discrete(expand = c(0,0)) +
  theme_ohchibi(size_panel_border = 0.3 ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.y = element_text(angle = 0,size = 6),
    legend.position = "right",
    axis.title.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5,size = 8)
    
  ) +
  xlab(label = "Cluster")  +
  scale_fill_manual(values = paleta)

oh.save.pdf(p = p_cluster_print,outname = "rootbranching_screening1_phylogeny_clusters_legend.pdf",outdir = "./figures/",width = 12,height = 16)
dev.off()

composition <- egg::ggarrange(p_blank,p_top + theme(legend.position = "top"),p_blank,p_blank,p_blank,p_blank,
                              p_tree,p_heatmap + theme(legend.position = "none",axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank()),
                              p_cluster + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),strip.text.y = element_blank(),legend.position = "none"),
                              p_bio  + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),legend.position = "none",axis.title.y = element_blank(),strip.text.y = element_blank()),
                              p_deg  + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),legend.position = "none",axis.title.y = element_blank(),strip.text.y = element_blank()),
                              p_chosen  + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),legend.position = "none",axis.title.y = element_blank(),strip.text.y = element_blank()),
                              nrow = 2,ncol = 6,
                              heights = c(0.2,1),widths  = c(0.1,1,0.01,0.01,0.01,0.01),draw = FALSE)


oh.save.pdf(p = composition,outname = "rootbranching_screening1_heatmap_phylogeny_estimates.pdf",
            outdir = "./figures/",height = 18,width = 14)

dev.off()
rm(list=ls())
gc()

