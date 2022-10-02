library(ohchibi)
library(ggtree)

set.seed(130816)


### Load dataset
melted <- readRDS(file = "./cleandata/dat_rootbranching_screening1.RDS")  %$% Dat_log10_zscore %$% Tab %>% melt %>%
  dplyr::rename(.data = .,Feature = Var1,root = Var2) %>%
  merge(readRDS(file = "./cleandata/dat_rootbranching_screening1.RDS")  %$% Dat_log10_zscore %$% Map, by = "root") %>%
  subset(Feature != "Norm_child_density") %>% droplevels

### Raw values to calculate coefficient of variation
merged_zscore <- dcast(data = melted,formula = Feature~Strain,
      fun.aggregate = function(x)mean(x,na.rm = T),value.var = "value") %>% melt %>%
  dplyr::rename(.data = .,Strain = variable)

#Matrix 
Tab <- acast(data = merged_zscore,formula = Strain~Feature,value.var = "value")
res_cor <- chibi.ggcor(Tab = Tab,hclust_method = "ward.D",cor.method = "pearson",
                       p.adjust.method = "fdr",p.adj.thres = 0.05)

df_cor <- res_cor$df_cor

df_cor$Var1 <- df_cor$Var1 %>% factor(levels = df_cor$Var1 %>% levels %>% rev)
p_cor <- ggplot(data = df_cor, aes(Var1, Var2)) + geom_raster(aes(fill = r)) + 
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37", limits = c(-1, 1)) + 
  scale_color_manual(values = c("#00000000","black")) + 
  theme_ohchibi(legend_proportion_size = NA,size_panel_border = 0.3) + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  geom_text(aes(label = rbold), parse = TRUE, 
            family = "Arial", size = 3) +
  geom_tile(aes(color = Significance), fill = "#00000000", 
            size = 0.3, width = 0.85, height = 0.85)



mclust <- res_cor$mclust 

#
Tab %>% as.data.frame() %>%
  ggplot(data = .,aes(PR_length,LR_density)) +
  geom_point()


#Load phylogenetic tree
p_tree <- ggtree(tr = as.phylo(res_cor$mclust),ladderize = TRUE,size = 0.3) +
  geom_tiplab(align = TRUE) +
  scale_x_continuous(limits = c( 0,10))

d <- p_tree$data %>% subset(isTip)
order_feature <- with(d, label[order(y, decreasing=F)])

k <- 13
cl_members <- cutree(tree = mclust, k = k)
plot(x = mclust, labels =  row.names(mclust), cex = 0.75)
rect.hclust(tree = mclust, k = k, which = 1:k, border = 1:k, cluster = cl_members)

df_cluster_feature <- cl_members %>%
  data.frame %>%
  tibble::rownames_to_column() %>%
  dplyr::rename(.data = .,Feature = rowname,Cluster = ".") %>%
  dplyr::mutate(.data = .,Cluster = paste0("Cluster",Cluster))

df_cluster_feature$Feature <- df_cluster_feature$Feature %>% factor(levels = order_feature)


order_clusters <- with(df_cluster_feature,order(Feature)) %>%
  df_cluster_feature$Cluster[.] %>% as.character %>% unique

df_cluster_feature$Cluster <- df_cluster_feature$Cluster %>% 
  factor(levels = order_clusters %>% rev)


#Read the information about the branching
df_cv <- readRDS(file = "./cleandata/dat_cv_phylosig_rootbranching_screening1.RDS") %>%
  merge(df_cluster_feature, by = "Feature") %>%
  subset(Feature != "Norm_child_density") %>% droplevels

df_cv$Feature <- df_cv$Feature %>%
  factor(levels = order_feature)


#Read the result matrix
merged <- readRDS(file = "./cleandata/res_rootbranching_screening1_vsNB.RDS") %$% Res_model_zscore %>%
  subset(Feature != "Norm_child_density") %>% droplevels


merged$Direction <- "NS"
merged$Direction[which((merged$padjFDRMW < 0.05) & (merged$estimate > 0))] <- "Higher NB, q < 0.05"
merged$Direction[which((merged$padjFDRMW < 0.05) & (merged$estimate < 0))] <- "Lower NB, q < 0.05"

merged$Significance <- merged$Direction 
merged$Significance[merged$Significance %>% grep(pattern = "q")] <- "q < 0.05"

#merged %>% subset(Feature == "LR_symmetry") %$% Significance %>% table
merged$Dummy <- 1

merged_sig <- merged


df_prop_sig <- acast(data = merged,formula = Direction~Feature,value.var = "Dummy",fun.aggregate = sum) %>%
  scale(center = F,scale = colSums(.)) %>% melt %>%
  dplyr::rename(.data = .,Direction = Var1,Feature = Var2,Prop = value)



df_raw_sig <- acast(data = merged,formula = Direction~Feature,value.var = "Dummy",fun.aggregate = sum) %>%
  melt %>%
  dplyr::rename(.data = .,Direction = Var1,Feature = Var2,Raw = value)


df_prop_sig <- merge(df_prop_sig,df_raw_sig,by = c("Direction","Feature"))

df_prop_sig$Feature <- df_prop_sig$Feature %>%
  factor(levels = order_feature)

df_prop_sig$Direction <- df_prop_sig$Direction %>%
  factor(levels = c("NS","Higher NB, q < 0.05","Lower NB, q < 0.05"))


df_prop_sig <- merge(df_prop_sig,df_cluster_feature, by = "Feature")

#Choose feature per cluster based on coefficient of variation
mclusters <- df_cv$Cluster %>% as.character %>% unique %>% sort

chosen_features <- NULL
for(mc in mclusters){
  
  chosen_features <- c(chosen_features,(df_cv %>%
                                          subset(Cluster == mc) %>%
                                          dplyr::arrange(.data = .,-CV_90) %$% Feature %>%
                                          as.character)[1])
  
}

df_cv$RepFeature <- "No"
df_cv$RepFeature[which(df_cv$Feature %in% chosen_features)] <- "Yes"


merged_chosen <- merged %>%
  dplyr::filter(.data = .,Feature %in% chosen_features) %>%
  droplevels

#### Strain perspective
df_sigdiff <- merged_chosen %>% 
  subset(Significance == "q < 0.05") %>%
  aggregate(Dummy~Strain,.,sum) %>%
  dplyr::rename(.data = .,NumberFeatureSigDiff = Dummy) %>%
  dplyr::arrange(.data = .,-NumberFeatureSigDiff) %>%
  dplyr::mutate(.data = .,PropNumberFeatureSigDiff = NumberFeatureSigDiff/length(chosen_features))

zero_st <- which(!((merged_chosen$Strain) %in% df_sigdiff$Strain)) %>% merged_chosen$Strain[.] %>% unique
df_sigdiff <- rbind(df_sigdiff,data.frame(Strain = zero_st,NumberFeatureSigDiff =0,PropNumberFeatureSigDiff =0))



#Create paneles
p_enrichment <- ggplot(data = df_prop_sig,aes(Feature,Prop)) +
  geom_bar(stat = "identity",aes(fill = Direction),alpha = 0.5) +
  geom_label(aes(label = Raw,fill = Direction)) +
  facet_grid(Cluster~.,space = "free",scales = "free") +
  theme_ohchibi(size_panel_border = 0.3) +
  theme(
    axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "top"
  ) +
  scale_y_continuous(breaks = seq(0,1,0.1),expand = c(0,0)) +
  scale_fill_manual(values = c("#B3B3B3","#E65037","#2480FF"),name = "Significance vs NB") +
  xlab(label = "Feature") +
  ylab(label = "Proportion of total strains (n = 391)") +
  coord_flip()



p_cv <- ggplot(data = df_cv,aes(Feature,CV_90)) +
  geom_bar(stat = "identity") +
  facet_grid(Cluster~.,space = "free",scales = "free") +
  theme_ohchibi(size_panel_border = 0.3) +
  coord_flip() +
  scale_y_continuous(expand = c(0,0)) +
  ylab(label = "Coefficient of variation")


#Create the heatmap of correlation
res_cor <- chibi.ggcor(Tab = Tab,hclust_method = "ward.D",cor.method = "pearson",
                       p.adjust.method = "fdr",p.adj.thres = 0.05)


df_cor <- res_cor$df_cor %>%
  dplyr::rename(.data = .,Feature = Var1,Feature2 = Var2)

df_cor <- merge(df_cor, df_cluster_feature, by = "Feature") 


df_cor$Feature <- df_cor$Feature %>% factor(levels = order_feature )
df_cor$Feature2 <- df_cor$Feature2 %>% factor(levels = order_feature %>% rev)



p_cor <- ggplot(data = df_cor, aes(Feature, Feature2)) + geom_raster(aes(fill = r)) + 
  #scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37", limits = c(-1, 1)) +
  scale_fill_paletteer_c("pals::kovesi.diverging_cwm_80_100_c22", limits = c(-1, 1)) +
  geom_tile(aes(color = Significance), fill = "#00000000", 
            size = 0.3, width = 0.85, height = 0.85) +
  scale_color_manual(values = c("#00000000", "black")) +
  facet_grid(Cluster~.,scales = "free",space = "free")+
   theme_ohchibi(legend_proportion_size = 0.5, 
                  size_panel_border = 0.3) + 
  theme(legend.position = "top", 
                           axis.title.x = element_blank(), axis.title.y = element_blank(), 
                    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  coord_flip()


## DRaw raw version
p_cor_num <- ggplot(data = df_cor, aes(Feature, Feature2)) + geom_raster(aes(fill = r)) + 
  #scale_fill_paletteer_c("pals::kovesi.diverging_bky_60_10_c30", limits = c(-1, 1)) +
  #scale_fill_gradientn(colours = c("#7EC4CF","white","#D4AFB9")) +
  #scale_fill_gradientn(colours = c("#A0CED9","white","#FFC09F")) +
  #scale_fill_gradientn(colours = c("#FCF4DD","white","#E8DFF5")) +
  #scale_fill_gradientn(colours = c("#F5E960","white","#DFB2F4")) +
  #scale_fill_gradientn(colours = c("#55D6C2","white","#DFB2F4") %>% rev) +
  #scale_fill_gradientn(colours = c("#F6BC66","white","#F55C7A") %>% rev) +
  #scale_fill_gradientn(colours = c("#F6BC66","white","#F55C7A")) +
  #scale_fill_gradientn(colours = c("#F5ED5D","white","#F55C7A")) +
  scale_fill_gradientn(colours = c("#51F5C9","white","#F55C7A"),limits = c(-1,1)) +
  
  
  #jcolors::scale_fill_jcolors_contin(palette = "pal10") +
  #scale_fill_paletteer_c("jcolors::pal12", limits = c(-1, 1),direction = -1) +
  #scale_fill_paletteer_c("scico::broc", limits = c(-1, 1),direction = 1) +
  
  #scale_fill_paletteer_c("pals::kovesi.diverging_cwm_80_100_c22", limits = c(-1, 1)) +
  #scale_fill_paletteer_c("pals::kovesi.diverging_gwv_55_95_c39", limits = c(-1, 1)) +
  geom_tile(aes(color = Significance), fill = "#00000000", 
            size = 0.3, width = 0.85, height = 0.85) +
  geom_text(aes(label = rbold),parse = TRUE,family = "Arial",size = 3) +
  scale_color_manual(values = c("#00000000", "black")) +
  facet_grid(Cluster~.,scales = "free",space = "free")+
  theme_ohchibi(legend_proportion_size = 0.5, 
                size_panel_border = 0.3) + 
  theme(legend.position = "top", 
        axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        strip.text.y = element_text(size = 9,angle = 0)) +
  coord_flip()

#Save the correlation figure with values
oh.save.pdf(p = p_cor_num,outname = "rootbranching_screening1_cor_features_numbers_pub.pdf",
            outdir = "./figures/",width = 11,height = 10)

dev.off()

###Create composition
composition <- egg::ggarrange(p_tree,
                              p_enrichment +
                                theme(strip.text.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank()),
                              p_cor +theme(strip.text.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank()),
               p_cv +
                 theme(strip.text.y = element_text(angle = 0,size = 10),axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank()),
               nrow = 1,widths = c(0.6,0.4,0.4,0.3),draw = F)



dev.off()
oh.save.pdf(p = composition,outname = "rootbranching_screening1_composition_dendrogram_cor_sig_cv.pdf",
            outdir = "./figures/",width = 24,height = 10)


##Save list
mlist <- list(
     df_feat_cv_phylosig =df_cv,
     df_feat_prop_sig_vs_nb = df_prop_sig,
     df_strain_numsig_vs_nb = df_sigdiff,
     df_sig_vs_nb_raw = merged_sig,
     df_cor = df_cor
  
     )




saveRDS(object = mlist,file = "./cleandata/res_cv_prop_sig_rootbranching_screening1.RDS")

rm(list=ls())
gc()
