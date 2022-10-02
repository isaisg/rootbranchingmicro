library(ohchibi)
library(xlsx)
library(emmeans)
library(multcomp)
library(ggtree)

set.seed(130816)


#Read the 
Res <- readRDS(file = "./cleandata/dat_ionome_lr_mutants_nutrients.RDS")

Dat <- Res$Dat_ionome

mpca <- oh.pca(Tab = Dat$Tab %>% t,Map = Dat$Map,
               retx = TRUE,center = FALSE,scale = FALSE,
               id_var = "SampleId")

Map_pca <- mpca$Map_pca

df_a <- Rmisc::summarySE(data = Map_pca,measurevar = "PC1",
                         groupvars = c("Genotype","Strain","Nutrient"))

df_b <- Rmisc::summarySE(data = Map_pca,measurevar = "PC2",
                         groupvars = c("Genotype","Strain","Nutrient"))


merged <- merge(df_a[,c(1,2,3,5,8)],df_b[,c(1,2,3,5,8)],
                by = c("Genotype","Strain","Nutrient"))

Map_pca$Strain <- Map_pca$Strain %>%
  factor %>% relevel(ref = "NB")
Map_pca$Genotype <- Map_pca$Genotype %>% factor
Map_pca$Nutrient <- Map_pca$Nutrient %>%
  factor

Map_pca$GenotypeStrain <- paste0(Map_pca$Genotype,"  ",Map_pca$Strain) %>% 
  factor(levels = c("Col 0  NB","Col 0  L344","Col 0  L359","Col 0  RMF27",
                    "lbd16  NB","lbd16  L359",
                    "nph4  NB","nph4  L344","nph4  L359",
                    "arf7/19  NB","arf7/19  L344",
                    "gnom1  NB","gnom1  RMF27"))

paleta <- paletteer_d("ggthemes::Tableau_10",n = 7) %>% as.character
paleta[1] <- "#FFCA64"
paleta[2] <- "#E15759FF"
paleta[3] <- "#4E79A7FF"
paleta[6] <- "#FF9DA7FF"

p <- ggplot(data = Map_pca,aes(PC1,PC2)) +
  geom_vline(xintercept = 0,linetype = "longdash",color = "#D9D9D9",size = 1)+
  geom_hline(yintercept = 0,linetype = "longdash",color = "#D9D9D9",size = 1)+
  geom_point(aes(color = GenotypeStrain,shape = Nutrient),size = 5) +
  theme_ohchibi(size_panel_border = 0.3) +
  scale_color_manual(values = paleta)+
  xlab(label = paste0("PC1 (",round(mpca$variance_explained[1],2),"%)"))+
  ylab(label = paste0("PC2 (",round(mpca$variance_explained[2],2),"%)")) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  ggtitle(label = "Shoot ionome")



#Perform heatmap 
Tab <- Dat$Tab %>%
  melt %>%
  dplyr::rename(.data = .,Ion = Var1,SampleId = Var2) %>%
  merge(Dat$Map, by = "SampleId") %>%
  acast(data =.,formula = Ion~Nutrient+Genotype+Strain,fun.aggregate = mean )

res_heatmap <- chibi.heatmap(Tab = Tab,
                             dist_method_rows = "euclidean",
                             dist_method_cols = "euclidean",
                             hclust_method_rows = "ward.D",
                             hclust_method_cols = "ward.D",
                             range_fill_heatmap = c(-1.5,1.5),
                             axis_ticks_row = TRUE,
                             size_axis_text_row = 8,
                             size_axis_text_col = 8,
                             k_rows = 5,
                             k_cols =3)

p_heatmap <- res_heatmap$p_heatmap +
  scale_fill_gradient2(low = "#021567",
                       mid = "white",high = "#FBDE28",midpoint = 0,limits = c(-1.5,1.5),oob = squish)

p_tree_rows <- res_heatmap$p_tree_rows
p_tree_cols <- res_heatmap$p_treee_cols

p_blank <- ggplot() + theme_void()

composition <- egg::ggarrange(p_blank,p_tree_cols,
                              p_tree_rows,p_heatmap + theme(legend.position = "bottom"),
                              nrow = 2,ncol = 2,
                              heights = c(0.2,1),widths  = c(0.1,1),draw = FALSE)

### Inspect the roots
df_root <- Res$df_root %>%
  subset(PR_length > 0) %>% droplevels

df_root$Genotype <- df_root$Genotype %>%
  factor(levels = c("Col0","arf7","nph4","lbd16","gnom1"))

df_root$Strain <- df_root$Strain %>% 
  factor(levels = c("NB","RMF27"))

df_root$Treatment <- paste0(df_root$Nutrient," ",df_root$Strain) %>%
  factor(levels = c("MS2 RMF27","MS2 NB","MS200 NB","MS2000 NB"))

#Perform testing within genotype
mgenotypes <- df_root$Genotype %>% levels

df_model <- NULL
for(mgeno in mgenotypes){
  
  df_temp <- df_root %>% 
    subset(Genotype == mgeno) %>% droplevels
  
  df_model <- lm(formula = PR_length~Treatment,df_temp) %>%
    emmeans(specs = "Treatment") %>%
    multcomp::cld() %>%
    as.data.frame %>%
    dplyr::mutate(.data = .,Genotype = mgeno,
                  Feature = "PrimaryRootLength") %>%
    rbind(df_model,.)
  
  df_model <- lm(formula = LR_Density~Treatment,df_temp) %>%
    emmeans(specs = "Treatment") %>%
    multcomp::cld() %>%
    as.data.frame %>%
    dplyr::mutate(.data = .,Genotype = mgeno,
                  Feature = "LRDensity") %>%
    rbind(df_model,.)
  
}

df_model %>%
  subset(Feature != "PrimaryRootLength")

#PRimary root plot
p1 <- ggplot(data = df_root,aes(Treatment,LR_Density)) +
  geom_sina(alpha = 0.3) +
  stat_summary(fun.data = mean_cl_normal,geom = "pointrange",color = "red",shape = 16,size = 1)+
  facet_grid(.~Genotype,space = "fixed",scales = "free") +
  theme_ohchibi(size_panel_border = 0.3) +
  scale_y_continuous(breaks = seq(0,6,1),limits = c(0,6),expand = c(0,0)) +
  theme(
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1) ,
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  ylab(label = "Lateral root density")


p2 <- ggplot(data = df_root,aes(Treatment,PR_length)) +
  geom_sina(alpha = 0.3) +
  stat_summary(fun.data = mean_cl_normal,geom = "pointrange",color = "red",shape = 16,size = 1)+
  facet_grid(.~Genotype,space = "fixed",scales = "free") +
  theme_ohchibi(size_panel_border = 0.3) +
  #scale_y_continuous(breaks = seq(0,6,1),limits = c(0,6),expand = c(0,0)) +
  theme(
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1) ,
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank()
    
  ) +
  ylab(label = "Primary root length") +
  scale_y_continuous(breaks = seq(0,16,2),limits = c(0,16),expand = c(0,0)) 

#Save figures
#Save plots
oh.save.pdf(p = p,outname = "fig2_ionome_mutants_nutrients_pca.pdf",outdir = "./figures/",width = 9,height = 8)

oh.save.pdf(p = composition,outname = "fig2_ionome_mutants_nutrients_heatmap.pdf",outdir = "./figures/",width = 8,height = 8)

oh.save.pdf(p = p1,outname = "fig2_ionome_mutants_nutrients_lr.pdf",outdir = "./figures/",width = 12,height = 10)

oh.save.pdf(p = p2,outname = "fig2_ionome_mutants_nutrients_pr.pdf",outdir = "./figures/",width = 12,height = 10)

rm(list=ls())
gc()
