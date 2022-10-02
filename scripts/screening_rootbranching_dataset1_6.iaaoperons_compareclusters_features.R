library(ohchibi)

set.seed(130816)

### Compare distribution of IAA degraders 

melted <- readRDS(file = "./cleandata/res_rootbranching_screening1_heatmap_melted.RDS")

### Quantitative approach comparing distributions ####
mclusters <- melted$ClusterRows %>%levels

Map <- melted[,c("IdRows","IAAOperons")] %>% unique
rownames(Map) <- Map$IdRows

Res_p <- list()
Res_pval <- NULL
#paleta_operon <- c("#E6DACC","#807971")
paleta_operon <- c("#F1DDBF","#84BFB9")
names(paleta_operon) <- c("No","Yes")

df_cara <- NULL
for(mc in mclusters){
  
  Tab <- melted %>% subset(ClusterRows == mc) %>% 
    acast(data =.,formula =  IdCols~IdRows,value.var = "value")
  Map_sub <- Map
  Map_sub <- match(colnames(Tab),rownames(Map_sub)) %>%
    Map_sub[.,]
  Dat <- create_dataset(Tab = Tab,Map = Map_sub)
  
  if (Dat$Map$IAAOperons %>% as.character %>% unique %>% length == 1){
    pval <- NA

  }else{
    Tab_euc <- vegdist(Dat$Tab %>% t,method = "euclidean")
    m1 <- adonis(formula = Tab_euc ~IAAOperons,data = Dat$Map,permutations = 9999)  
    pval <- m1$aov.tab$`Pr(>F)`[1]
   
  }
  
  mpca <- oh.pca(Tab = Dat$Tab %>% t,Map = Dat$Map,center = F,retx = T,scale = F,id_var = "IdRows")
  p <- chibi.pca(list_ohpca = mpca,col_val = "IAAOperons",stroke = 0.3,size = 4,lines_zero = FALSE) +
    theme_ohchibi(size_panel_border = 0.15,size_axis_text.x = 8,size_axis_text.y = 8,size_axis_title.y = 8,size_axis_title.x = 8,size_title_text = 10) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "none",
      plot.title = element_text(size = 10)
    ) +
    scale_fill_manual(values =paleta_operon) +
    ggtitle(label = paste0(mc,"  (PERMANOVA  p = ",round(pval,3),")"))
  Res_p[[mc]] <- p
  Res_pval <- data.frame(Cluster = mc,pvalue = pval) %>%
    rbind(Res_pval,.)
  
  df_cara <- Map_sub$IAAOperons %>% table %>%
    as.data.frame %>%
    dplyr::rename(.data = .,Type = ".") %>%
    dplyr::mutate(.data = .,Cluster = mc) %>%
    rbind(df_cara,.)
  
}


#Create figure solely counting of distributions
df_bar <- melted[,c("IdRows","ClusterRows","IAAOperons")] %>% unique %>%
  dplyr::mutate(.data = .,Bar = "Bar",Dummy = 1) 
 

df_sum <- dcast(data = df_bar,formula = ClusterRows~IAAOperons,
                fun.aggregate = function(x)sum(x,na.rm = TRUE),value.var = "Dummy") %>%
  dplyr::mutate(.data = .,TotalStrains = No + Yes) %>%
  dplyr::mutate(.data = .,PropNo = No/TotalStrains,PropYes = Yes/TotalStrains)

melted_prop <- melt(data = df_sum,id.vars = "ClusterRows",measure.vars = c("PropNo","PropYes")) %>%
  dplyr::mutate(.data = .,variable = variable %>% gsub(pattern = "Prop",replacement = ""))
melted_raw <- melt(data = df_sum,id.vars = "ClusterRows",measure.vars = c("No","Yes"))

melted_united <- merge(melted_prop,melted_raw, by = c("ClusterRows","variable")) %>%
  dplyr::rename(.data = .,IAAOperons = variable,PropStrains =value.x, NumStrains = value.y)

p <- ggplot(data = melted_united,aes(IAAOperons,ClusterRows)) +
  geom_tile(aes(fill = IAAOperons),color = "black") +
  geom_text(aes(label = NumStrains))+
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0))+
  scale_fill_manual(values =paleta_operon) +
  theme_ohchibi(size_panel_border = 0.3) +
  theme(
    legend.position = "right"
  ) +
  ylab(label = "Cluster") +
  ggtitle(label = "Number of strains per cluster")



composition <- egg::ggarrange(
      Res_p$CR7,Res_p$CR9,Res_p$CR1,Res_p$CR6,
      Res_p$CR22,Res_p$CR5,Res_p$CR15,Res_p$CR17,
      Res_p$CR20,Res_p$CR24,Res_p$CR8,Res_p$CR19,
      Res_p$CR12,Res_p$CR13,Res_p$CR11,Res_p$CR10,
      Res_p$CR14,Res_p$CR21,Res_p$CR4,Res_p$CR3,
      Res_p$CR18,Res_p$CR16,Res_p$CR2,Res_p$CR23,
      nrow = 6,ncol = 4,draw = FALSE)

oh.save.pdf(p = composition,outname = "rootbranching_screening1_clusters_pca_permanova_comparison_iaaoperons.pdf",
            outdir = "./figures/",height = 18,width = 12)

dev.off()

### Test it at the variable level ####
Res_all <- readRDS(file = "./cleandata/res_cv_prop_sig_rootbranching_screening1.RDS")

chosen_features <- Res_all$df_feat_cv_phylosig %>% 
  subset(RepFeature == "Yes") %$% Feature %>% 
  as.character

df_strain <- Res_all$df_strain_numsig_vs_nb 

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
merged <- Res_all$df_sig_vs_nb_raw %>% 
  dplyr::filter(.data = .,Feature %in% chosen_features)  %>% droplevels

merged$IAAOperons <- match(merged$Strain,melted$IdRows) %>%
  melted$IAAOperons[.]

merged$SelectedStrains <- match(merged$Strain,melted$IdRows) %>%
  melted$SelectedStrains[.]

### Compare distributions  using kolgomeronov smirnof test 
merged$Feature <- merged$Feature %>% factor(levels = melted$IdCols %>% levels)

#Test the distributions
mfeatures <- merged$Feature %>% levels

Res_ks <- NULL
for(mfeature in mfeatures){
  
  merged_sub <- merged %>% subset(Feature == mfeature) %>% droplevels
  x <- merged_sub %>% subset(IAAOperons == "Yes") %$% estimate
  y <- merged_sub %>% subset(IAAOperons == "No") %$% estimate
  m1 <- ks.test(x,y,alternative = "two.sided",exact = FALSE)  
  Res_ks <- data.frame(Feature = mfeature,pvalue =m1$p.value) %>%
    rbind(Res_ks,.)
}

Res_ks$Label <- paste0("p=",round(Res_ks$pvalue ,3))
Res_ks$Feature <- Res_ks$Feature %>% factor(levels = merged$Feature %>% levels)

p <- ggplot(data = merged,aes(estimate)) +
  geom_density(aes(fill = IAAOperons,linetype = IAAOperons),alpha = 0.7,color = "black",size = 0.3) +
  facet_wrap(facets = "Feature",ncol = 5,nrow = 3) +
  geom_text(data = Res_ks,aes(x = 1.5,y = 0.3,label = Label)) +
  theme_ohchibi(size_panel_border = 0.3,size_axis_text.x = 8,size_axis_text.y = 8,size_axis_title.y = 8,size_axis_title.x = 8,size_title_text = 10) +
  ylab(label = "Density") +
  xlab(label = "Estimate vs NB") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank() ,
    legend.position = "none"
  )  +
scale_fill_manual(values = paleta_operon) +
  scale_color_manual(values = paleta_operon)

### Save figure ###
oh.save.pdf(p = p,outname = "rootbranching_screening1_features_comparison_iaaoperons.pdf",
            outdir = "./figures/",height = 8,width = 14)


rm(list=ls())
dev.off()
gc()

