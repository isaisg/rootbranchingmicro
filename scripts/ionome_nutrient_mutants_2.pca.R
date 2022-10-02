library(ohchibi)
library(emmeans)
library(ggrepel)
library(dplyr)
library(ggtree)


set.seed(130816)


distfun <- function(x,method) vegan::vegdist(x = x, method = "euclidean",na.rm = T)


#Read matrix
Dat_all <- readRDS(file = "./cleandata/dat_lr_ionome_nutrient_mutants.RDS")

Dat_z <- Dat_all$Dat_log10_z

#Projection
Dat_sub <- Dat_z

mpca <- oh.pca(Tab = Dat_sub$Tab %>% t,Map = Dat_sub$Map,
       retx = T,center = F,scale = F,id_var = "UiD")

chibi.pca(list_ohpca = mpca,col_val = "GenotypeBacteria",shape_val = "Medium",size = 8) +
  theme_ohchibi(legend_proportion_size = 0.01,size_panel_border = 0.3) +
  theme(
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  scale_color_paletteer_d("ggthemes::Tableau_20") +
  ggtitle(label = "Shoot ionome")

Tab_av <- Dat_sub$Tab %>% melt %>%
  dplyr::rename(.data = .,Ion = Var1,UiD = Var2) %>%
  merge(Dat_sub$Map, by = "UiD") %>%
  dplyr::mutate(.data = .,MediumGenotypeBacteria = paste0(Medium,"|",GenotypeBacteria) %>% gsub(pattern = "\\/",replacement = "")) %>%
  acast(data = .,Ion~MediumGenotypeBacteria,fun.aggregate = function(x)mean(x,na.rm = TRUE),value.var = "value")


res_heatmap <- chibi.heatmap(Tab = Tab_av,
                             dist_method_rows = "euclidean",dist_method_cols = "euclidean",
                             hclust_method_rows = "ward.D",hclust_method_cols = "ward.D",
                             range_fill_heatmap = c(-2,2),
                             k_rows = 4,k_cols = 3,axis_ticks_row = TRUE,
                           size_axis_text_col = 8,size_axis_text_row = 8)
a <- res_heatmap$heatmap

oh.save.pdf(p = a,outname = "heatmap_ionome_nutrients.pdf",outdir = "./figures/",width = 12,height = 10)
