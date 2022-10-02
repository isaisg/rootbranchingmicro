library(ohchibi)
library(xlsx)
set.seed(130816)

Res_all <- readRDS(file = "./cleandata/res_cv_prop_sig_rootbranching_screening1.RDS")

#Plot counting number of modified features
df <- Res_all$df_strain_numsig_vs_nb
df$NumberFeatureSigDiff <- df$NumberFeatureSigDiff %>% factor

p <- ggplot(data = df,aes(NumberFeatureSigDiff)) +
  geom_bar(color = "black",fill = "#8FC7BF",size = 0.3) +
  theme_ohchibi(size_panel_border = 0.3)+
  scale_y_continuous(breaks = seq(0,80,10),expand = c(0,0),limits = c(0,80)) +
  ylab(label = "Number of strains") +
  xlab(label = "Number of features significantly modified with respect to NB") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank()
  )

oh.save.pdf(p = p,outname = "rootbranching_screening1_numsignificant.pdf",
            outdir = "./figures/",width = 8,height = 6)


#Plot density
df <- Res_all$df_sig_vs_nb_raw %>%
  dplyr::filter(.data = .,Feature %in% c("PR_length","LR_density")) %>%
  dcast(data = .,Strain~Feature,value.var = "estimate")

df_sig <- Res_all$df_sig_vs_nb_raw %>%
  dplyr::filter(.data = .,Feature %in% c("PR_length","LR_density")) %>%
  dcast(data = .,Strain~Feature,fun.aggregate = function(x)paste0(x,collapse = " "),value.var = "Significance") %>%
  dplyr::mutate(.data = .,SignificanceBoth = paste0(PR_length," ",LR_density))

df <- merge(df,df_sig[,c("Strain","SignificanceBoth")], by = "Strain")

p1 <- ggplot(data = df,aes(PR_length,LR_density)) +
  geom_vline(xintercept = 0,linetype = "longdash",color = "#D9D9D9",size = 0.75) +
  geom_hline(yintercept = 0,linetype = "longdash",color = "#D9D9D9",size = 0.75) +
  geom_point(shape = 16,size = 2) +
  theme_ohchibi(size_panel_border = 0.3)+
  xlab(label = "Primary root length estimate with respect to NB") +
  ylab(label = "Lateral root density estimate with respect to NB") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  ggpubr::stat_cor(label.sep = "\n",size = 5) +
  scale_x_continuous(limits = c(-3.2,1.2)) +
  scale_y_continuous(limits = c(-1.75,2.1))


#Density
df <- Res_all$df_sig_vs_nb_raw %>%
  dplyr::filter(.data = .,Feature %in% c("PR_length","LR_density")) %>% droplevels

m1 <- ks.test(x = df %>% subset(Feature == "PR_length") %$% estimate,
        y = df %>% subset(Feature == "LR_density") %$% estimate,exact = FALSE)
m1$statistic

p2 <- ggplot(data = df %>% subset(Feature == "PR_length"),aes(estimate)) +
  geom_density(aes(group = Feature,fill = Feature),alpha = 0.3,color = "black")+
  theme_ohchibi(size_panel_border = 0.3)+
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.border = element_blank()
  ) +
  scale_x_continuous(limits = c(-3.2,1.2)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.7))+
  scale_fill_manual(values = c("#8568E8"))  

p3 <- ggplot(data = df %>% subset(Feature == "LR_density"),aes(estimate)) +
  geom_density(aes(group = Feature,fill = Feature),alpha = 0.3,color = "black")+
  theme_ohchibi(size_panel_border = 0.3)+
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.border = element_blank()
    
  ) +
  scale_y_continuous(limits = c(-1.75,2.1)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.7))+
  scale_fill_manual(values = c("#9C9136"))  

p_blank <- ggplot() + theme_void()
  
composition <- egg::ggarrange(p2,p_blank,
               p1,p3 + coord_flip(),nrow = 2,ncol = 2,
               widths = c(1,0.2),heights = c(0.2,1),draw = FALSE)
dev.off()

oh.save.pdf(p = composition,outname = "rootbranching_screening1_cor_primroot_lr.pdf",
            outdir = "./figures/",width = 9,height = 8)

#Graph with ks test
mtext <- paste0("D = ",round(m1$statistic,2),", p < 2.2e-16")
p4 <- ggplot(data = df,aes(estimate)) +
  geom_density(aes(group = Feature,fill = Feature),alpha = 0.3,color = "black")+
  theme_ohchibi(size_panel_border = 0.3)+
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none"
    
  ) +
  scale_y_continuous(limits = c(-1.75,2.1)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.7)) +
  scale_fill_manual(values = c("#8568E8","#9C9136") %>% rev)   +
  ylab(label = "Density") +
  xlab(label = "Estimate with respect to NB") +
  annotate(geom = "text",x = -2,y = 0.5,label =mtext,size = 5 )


oh.save.pdf(p = p4,outname = "rootbranching_screening1_density_primroot_lr.pdf",
            outdir = "./figures/",width = 9,height = 8)



#Number of features ###
chosen_features <- Res_all$df_feat_cv_phylosig %>% 
  subset(RepFeature == "Yes") %$% Feature %>% 
  as.character()

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
merged <- Res_all$df_sig_vs_nb_raw %>% 
  dplyr::filter(.data = .,Feature %in% chosen_features)  %>% droplevels

Tab <- acast(data = merged,formula = Strain~Feature,value.var = "estimate",fill = 0)




df_border <- merged[,c("Strain","Feature","Significance")] %>%
  dplyr::rename(.data = .,IdRows = Strain,IdCols = Feature) 
df_border$Significance <- df_border$Significance %>% gsub(pattern = "NS",replacement = NA)

res_heatmap <- chibi.heatmap(Tab = Tab,df_border = df_border,
                             width_border_tile = 0.05,height_border_tile = 0.05,size_border_tile = 0.2,
                             dist_method_rows = "euclidean",
                             dist_method_cols = "euclidean",
                             hclust_method_rows = "ward.D",
                             hclust_method_cols = "ward.D",range_fill_heatmap = c(-1.5,1.5),
                             k_rows = 24,k_cols = 9,
                             panel_spacing = 0.2,size_axis_text_row = 6,axis_ticks_row = FALSE,mtheme = theme(axis.ticks.y = element_blank()))

###
res_heatmap$heatmap

melted <- res_heatmap$melted
melted$NumberFeatureSigDiff <- match(melted$IdRows,df_strain$Strain) %>%
  df_strain$NumberFeatureSigDiff[.]

melted$Sum_IAA_Bio <- match(melted$IdRows,df_iaa_acc_genome$Strain) %>%
  df_iaa_acc_genome$Sum_IAA_Bio[.]

melted$Sum_IAA_Deg <- match(melted$IdRows,df_iaa_acc_genome$Strain) %>%
  df_iaa_acc_genome$Sum_IAA_Deg[.]

melted$SelectedStrains <- NA
melted$SelectedStrains[which(melted$IdRows %in% selected_strains )] <- "Yes"


df_bar <- melted[,c("IdRows","ClusterRows","NumberFeatureSigDiff","Sum_IAA_Bio","Sum_IAA_Deg","SelectedStrains")] %>% unique %>%
    dplyr::mutate(.data = .,Bar = "Bar")

df_bar$Sum_IAA_Bio <- df_bar$Sum_IAA_Bio %>% factor
df_bar$Sum_IAA_Deg <- df_bar$Sum_IAA_Deg %>% factor

p_bio <- ggplot(data = df_bar,aes(Bar,IdRows)) +
  geom_tile(aes(fill =Sum_IAA_Bio ),color = "transparent") +
  facet_grid(ClusterRows~.,space = "free",scales = "free") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_ohchibi(size_panel_border = 0.3) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.y = element_text(angle = 0,size = 6),
    axis.title.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5,size = 8),
    panel.spacing.x = unit(0.2, "lines"),
    panel.spacing.y = unit(0.2, "lines")
    
  ) +
  xlab(label = "IAA\nbiosynthesis") +
  scale_fill_manual(values = c("#F1DDBF","#84BFB9","#84BFB9","#84BFB9","#84BFB9"))



p_deg <- ggplot(data = df_bar,aes(Bar,IdRows)) +
  geom_tile(aes(fill =Sum_IAA_Deg ),color = "transparent") +
  facet_grid(ClusterRows~.,space = "free",scales = "free") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_ohchibi(size_panel_border = 0.3) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.y = element_text(angle = 0,size = 6),
    axis.title.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5,size = 8),
    panel.spacing.x = unit(0.2, "lines"),
    panel.spacing.y = unit(0.2, "lines")
    
  ) +
  xlab(label = "IAA\ndegradation") +
  scale_fill_manual(values = c("#F1DDBF","#84BFB9","#84BFB9","#84BFB9","#84BFB9"))


p_num <- ggplot(data = df_bar,aes(Bar,IdRows)) +
  geom_tile(aes(fill =NumberFeatureSigDiff ),color = "transparent") +
  facet_grid(ClusterRows~.,space = "free",scales = "free") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_ohchibi(size_panel_border = 0.3) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.y = element_text(angle = 0,size = 6),
    legend.position = "top",
    axis.title.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5,size = 8),
    panel.spacing.x = unit(0.2, "lines"),
    panel.spacing.y = unit(0.2, "lines")
    
  ) +
  xlab(label = "Number of features changed") +
  scale_fill_paletteer_c("pals::kovesi.rainbow_bgyrm_35_85_c71")

#Plot the selected

p_chosen <- ggplot(data = df_bar,aes(Bar,IdRows)) +
  geom_tile(aes(fill =SelectedStrains ),color = "transparent") +
  facet_grid(ClusterRows~.,space = "free",scales = "free") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_ohchibi(size_panel_border = 0.3) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.y = element_text(angle = 0,size = 6),
    legend.position = "top",
    axis.title.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5,size = 8),
    panel.spacing.x = unit(0.2, "lines"),
    panel.spacing.y = unit(0.2, "lines")
  ) +
  xlab(label = "Selected\nstrains") +
  scale_fill_manual(values = c("#84BFB9"),na.value = "#F1DDBF")



### Create comoposition
p_heatmap <- res_heatmap$p_heatmap +
  scale_fill_gradient2(low = "#021567",
                       mid = "white",high = "#FBDE28",midpoint = 0,limits = c(-1.5,1.5),oob = squish)

p_tree_rows <- res_heatmap$p_tree_rows
p_tree_cols <- res_heatmap$p_treee_cols

p_blank <- ggplot() + theme_void()

composition <- egg::ggarrange(p_blank,p_tree_cols,p_blank,p_blank,p_blank,p_blank,
               p_tree_rows,p_heatmap + theme(legend.position = "bottom"),
               p_num + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),strip.text.y = element_blank()),
               p_bio  + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),legend.position = "none",axis.title.y = element_blank(),strip.text.y = element_blank()),
               p_deg  + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),legend.position = "none",axis.title.y = element_blank(),strip.text.y = element_blank()),
               p_chosen  + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),legend.position = "none",axis.title.y = element_blank(),strip.text.y = element_blank()),
               nrow = 2,ncol = 6,
               heights = c(0.2,1),widths  = c(0.1,1,0.01,0.01,0.01,0.01),draw = FALSE)

oh.save.pdf(p = composition,outname = "rootbranching_screening1_heatmap_clustered_estimates.pdf",
            outdir = "./figures/",height = 18,width = 12)


dev.off()



df_bar$IAAOperons <- "Yes"
df_bar$IAAOperons[which((df_bar$Sum_IAA_Bio=="0") & (df_bar$Sum_IAA_Deg=="0"))] <- "No"
df_bar$Dummy <- 1




melted$IAAOperons <- match(melted$IdRows,df_bar$IdRows) %>%
  df_bar$IAAOperons[.]

### Save structures
saveRDS(object = melted,file = "./cleandata/res_rootbranching_screening1_heatmap_melted.RDS")
rm(list=ls())
gc()
dev.off()

