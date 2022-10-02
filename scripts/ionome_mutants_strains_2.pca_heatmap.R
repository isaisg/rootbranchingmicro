library(ohchibi)
library(emmeans)
library(ggrepel)
library(dplyr)
library(Rmisc)


set.seed(130816)


distfun <- function(x,method) vegan::vegdist(x = x, method = "euclidean",na.rm = T)


#Read matrix
Dat_all <- readRDS(file = "./cleandata/dat_rootbranching_ionome_mutants_strains.RDS")

Dat_z <- Dat_all$Dat_log10_z

#Projection
Dat_sub <- Dat_z

mpca <- oh.pca(Tab = Dat_sub$Tab %>% t,Map = Dat_sub$Map,
               retx = T,center = F,scale = F,id_var = "UiD")

#Perform summarization
df_pc1 <- summarySE(data = mpca$Map_pca,measurevar = "PC1",groupvars = c("Bacteria","Genotype"))
df_pc2 <- summarySE(data = mpca$Map_pca,measurevar = "PC2",groupvars = c("Bacteria","Genotype"))

merged <- merge(df_pc1,df_pc2, by = c("Bacteria","Genotype")) 

p <- ggplot(data = merged,aes(PC1,PC2)) +
  geom_vline(xintercept = 0,linetype = "longdash",color = "#D9D9D9",size = 0.5) +
  geom_hline(yintercept = 0,linetype = "longdash",color = "#D9D9D9",size = 0.5) +
  geom_linerange(mapping = aes(ymin = PC2-ci.y,ymax = PC2+ci.y),alpha = 0.1)+
  geom_linerange(mapping = aes(xmin = PC1-ci.x,xmax = PC1+ci.x),alpha = 0.1)+
  geom_point(aes(color = Bacteria,shape = Genotype),size = 6) +
  theme_ohchibi(size_panel_border = 0.3) +
  scale_color_paletteer_d("ggthemes::Tableau_20") +
  ggtitle(label = "Shoot ionome") +
  scale_shape_manual(values = c(15,16,14,17,18)) +
  scale_x_continuous(limits  = c(-3,9),oob = rescale_none) +
  scale_y_continuous(limits  = c(-2.5,2.5),oob = rescale_none) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  xlab(label = paste0("PC1 (",round(mpca$variance_explained[1],2),"%)")) +
  ylab(label = paste0("PC2 (",round(mpca$variance_explained[2],2),"%)"))





### Perform complex heatmap 
melted <- Dat_z$Tab %>%
  melt %>%
  dplyr::rename(.data = .,
                Ion = Var1,UiD = Var2) %>%
  merge(Dat_z$Map, by = "UiD")

mions <- melted$Ion %>% levels

Res_contrast <- NULL
Res_em <- NULL

melted$GenotypeBacteria

for(mion in mions){
  
  m1 <- melted %>%
    subset(Ion == mion) %>%
    lm(formula = value~Genotype+Bacteria+Genotype:Bacteria,.)
  
  res_em <- emmeans(m1,pairwise~Genotype:Bacteria,adjust = "none")
  
  df_em <- res_em$emmeans %>% 
    as.data.frame
  
  df_inter <- res_em$contrasts %>%
    as.data.frame %>%
    dplyr::filter(.data = .,
                  contrast %>% grepl(pattern = ".*NB.*NB")) %>%
    dplyr::filter(.data = .,contrast %>% grepl(pattern = "Col-0")) %>%
    dplyr::mutate(.data = .,
                  Genotype = contrast %>%
                    gsub(pattern = ".* - ",replacement = "") %>%
                    gsub(pattern = "\\(",replacement = "") %>%
                    gsub(pattern = "\\)",replacement = "") %>%
                    gsub(pattern = " .*",replacement = ""),
                  Bacteria = "NB",
                  Contrast = "InterGenotype"
    )
  
  df_inter$estimate <- df_inter$estimate * -1
  
  #Pull the same genotype
  df_intra <- res_em$contrasts %>%
    as.data.frame %>%
    dplyr::filter(.data = .,contrast %>% grepl(pattern = "NB")) %>%
    dplyr::filter(.data = .,
                  contrast %>% 
                    grepl(pattern = "(.*Col-0.*Col-0|.*nph4.*nph4|.*lbd16.*lbd16|.*gnom1.*gnom1|.*arf7/19.*arf7/19)")) %>%
    dplyr::mutate(.data = .,
                  Genotype = contrast %>%
                    gsub(pattern = "\\(",replacement = "") %>%
                    gsub(pattern = "\\)",replacement = "") %>%
                    gsub(pattern = " .*",replacement = ""),
                  Bacteria = contrast %>%
                    gsub(pattern = ".* ",replacement = "") %>%
                    gsub(pattern = "\\)",replacement = ""),
                  Contrast = "IntraGenotype"
    )
  
  df_intra$estimate <- df_intra$estimate * -1
  
  
  Res_contrast <- rbind(df_inter,df_intra) %>%
    dplyr::mutate(.data = .,Ion = mion) %>%
    rbind(Res_contrast,.)
  
  
  Res_em <- df_em %>%
    dplyr::mutate(.data = .,Ion = mion) %>%
    rbind(Res_em,.)
  
}


Res_em$emmean %>% sort %>% plot
Res_em$Genotype <- Res_em$Genotype %>%
  factor(levels = c("Col-0","nph4","lbd16","arf7/19","gnom1"))

Res_em_sub <- which(!(is.na(Res_em$emmean))) %>%
  Res_em[.,] %>% droplevels

###PRocess the contrast dataset
df_contrast_inter <- Res_contrast %>% 
  subset(Contrast == "InterGenotype") %>%
  droplevels
df_contrast_inter$p.adj <- df_contrast_inter$p.value %>%
  p.adjust(method = "fdr")
df_contrast_inter$SignificanceInterGenotype <- NA
df_contrast_inter$SignificanceInterGenotype[which(df_contrast_inter$p.adj < 0.05)] <- "q < 0.05"

df_contrast_intra <- Res_contrast %>% 
  subset(Contrast != "InterGenotype") %>%
  droplevels
df_contrast_intra$p.adj <- df_contrast_intra$p.value %>%
  p.adjust(method = "fdr")
df_contrast_intra$SignificanceIntraGenotype <- NA
df_contrast_intra$SignificanceIntraGenotype[which(df_contrast_intra$p.adj < 0.05)] <- "q < 0.05"


Res_em_sub <- merge(Res_em_sub,df_contrast_inter[,c("Genotype","Bacteria","Ion","SignificanceInterGenotype")],
                    by = c("Genotype","Bacteria","Ion"),all.x = TRUE) %>%
  merge(df_contrast_intra[,c("Genotype","Bacteria","Ion","SignificanceIntraGenotype","estimate")],
        by = c("Genotype","Bacteria","Ion"),all.x = TRUE)

Res_em_sub$Label <- Res_em_sub$SignificanceIntraGenotype
Res_em_sub$Label <- Res_em_sub$Label %>%
  gsub(pattern = "q < 0.05",replacement  = "*")

##Determine the order of the heatmap
Tab <- acast(data = Res_em_sub,formula = Bacteria~Ion,fun.aggregate = mean,
             value.var = "emmean")
mclust_ion <- hclust(d = as.dist(1-cor(Tab)),method = "ward.D")
mclust_bacteria <- hclust(d = dist(Tab),method = "ward.D")

order_ions <- mclust_ion$order %>%
  mclust_ion$labels[.]
df_clust_ion <- mclust_ion %>%
  cutree(tree = .,k = 3) %>%
  data.frame(Ion = names(.),ClusterIon = paste0("ClIon",.),row.names = NULL)
df_clust_ion <- df_clust_ion[,-1]

order_bacteria <- mclust_bacteria$order %>%
  mclust_bacteria$labels[.]

Res_em_sub <- merge(Res_em_sub, df_clust_ion, by = "Ion")
Res_em_sub$Ion <- Res_em_sub$Ion %>%
  factor(levels = order_ions)
Res_em_sub$Bacteria <- Res_em_sub$Bacteria %>%
  factor(levels = order_bacteria)
Res_em_sub$Bacteria <- Res_em_sub$Bacteria %>% relevel(ref = "NB")

order_cluster_ion <- with(Res_em_sub,order(Ion)) %>%
  Res_em_sub$ClusterIon[.] %>% unique

Res_em_sub$ClusterIon <- Res_em_sub$ClusterIon %>%
  factor(levels = order_cluster_ion %>% rev)



p2 <- ggplot(data = Res_em_sub,aes(Bacteria,Ion)) +
  geom_raster(aes(fill = emmean)) +
  geom_tile(fill = "transparent",aes(color = SignificanceInterGenotype),
            width = 0.85,height = 0.85,size = 1) +
  geom_text(aes(label = Label),color = "black",size = 5)+
  facet_grid(ClusterIon~Genotype,space = "free",scales = "free") +
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-3.5,3.5),na.value = "#D9D9D9",oob = squish,name = "Abundance \nz-score") +
  scale_color_manual(values = c("black"),na.value =  "transparent",name = "Significance\nInterGenotype NB vs Col-0 NB") +
  theme_ohchibi(size_panel_border = 0.3) +
  theme(
    axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
    legend.position = "bottom",
    strip.text.y = element_blank(),
    strip.text.x = element_blank()
  ) +
  #scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0))

Res_em_sub$BacteriaGenotype <- paste0(Res_em_sub$Bacteria,"_",Res_em_sub$Genotype)

#Read the root data
df_root <- readRDS(file = "./cleandata/temp_forpresentation_mutants_lr.RDS")
df_root$Genotype <- df_root$Genotype %>% as.character %>%
  gsub(pattern = "Col0",replacement = "Col-0") %>%
  gsub(pattern = "nph4-1",replacement = "nph4") %>%
  gsub(pattern = "arf7-19",replacement = "arf7/19") %>%
  gsub(pattern = "lbd16-1",replacement = "lbd16")


df_root$BacteriaGenotype <- paste0(df_root$Strain,"_",df_root$Genotype)

df_root <- df_root %>%
  dplyr::filter(.data = .,BacteriaGenotype %in% unique(Res_em_sub$BacteriaGenotype)) %>%
  droplevels

df_root$Strain <- df_root$Strain %>% factor(levels = Res_em_sub$Bacteria %>% levels)

df_root$Genotype <- df_root$Genotype %>% factor(levels = Res_em_sub$Genotype %>% levels)

p_up <- ggplot(data = df_root,aes(Strain,value)) +
  geom_point(aes(color = Direction),alpha = 0.3) +
  geom_boxplot(aes(color = Direction),outlier.size = NA,outlier.color = NA,fill = "transparent")+
  facet_grid(.~Genotype,space = "free",scales = "free") +
  theme_ohchibi(size_panel_border = 0.3,size_axis_text.x = 8) +
  ylab(label = "Lateral root density\n(Number lateral roots by cm)") + 
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    axis.title.x = element_blank()
  )  +
  scale_fill_manual(values = c("#008C51","#E65037","#2480FF"),name = "Significance vs NB",na.value = "#B3B3B3") +
  scale_color_manual(values = c("#008C51","#E65037","#2480FF"),name = "Significance vs NB",na.value = "#B3B3B3") 

composition <- egg::ggarrange(p_up,p2,ncol = 1,heights = c(0.4,1))

oh.save.pdf(p = composition,outname = "comp_ionome_lateral_root_zscore.pdf",outdir = "./figures/",width = 16,height = 10)


p3 <- ggplot(data = Res_em_sub,aes(Bacteria,Ion)) +
  geom_raster(aes(fill = estimate)) +
  geom_tile(fill = "transparent",aes(color = SignificanceInterGenotype),
            width = 0.85,height = 0.85,size = 1) +
  geom_text(aes(label = Label),color = "black",size = 5)+
  facet_grid(ClusterIon~Genotype,space = "free",scales = "free") +
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-3.5,3.5),na.value = "#D9D9D9",oob = squish,name = "Estimate IntraGenotype\nvs NB") +
  scale_color_manual(values = c("black"),na.value =  "transparent",name = "Significance\nInterGenotype NB vs Col-0 NB") +
  theme_ohchibi(size_panel_border = 0.3) +
  theme(

    axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
    legend.position = "bottom",
    strip.text.y = element_blank(),
    strip.text.x = element_blank()
    
  ) +
  #scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0))

composition <- egg::ggarrange(p_up,p3,ncol = 1,heights = c(0.4,1))

oh.save.pdf(p = composition,outname = "comp_ionome_lateral_root_estimate.pdf",outdir = "./figures/",width = 16,height = 10)


#Save both compositions


oh.save.pdf(p = p,outname = "ionome_mutant_strains_pc.pdf",outdir = "./figures/",width = 10,height = 8)
oh.save.pdf(p = p2,outname = "ionome_mutant_strains_heatmap.pdf",outdir = "./figures/",width = 18,height = 10)

rm(list=ls())
dev.off()
gc()

