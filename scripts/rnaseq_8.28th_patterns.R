library(ohchibi)
library(dplyr)
library(ComplexHeatmap)
library(biomaRt)
library(org.At.tair.db)
library(clusterProfiler)
library(ggtree)
library(xlsx)
library(DOSE)
library(ggrepel)



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


Res_regulons <- readRDS(file = "./cleandata/regulons_rnaseq.RDS")
eth <- Res_regulons$ethylene

ids_ppp <- which(ids_ppp %in% eth) %>% ids_ppp[.]


#Gene enrichment
mgo <- enrichGO(gene =ids_ppp,
                keyType       = "TAIR",
                OrgDb         = org.At.tair.db,
                ont           = "BP",
                pAdjustMethod = "fdr",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.1
)

p_cnet <- cnetplot(x = mgo,
                   showCategory = 20, 
                   layout = "kk", 
                   circular = FALSE,
                   node_label = "category",
                   cex_label_category = 1,
                   cex_category = 1,
) +
  theme(legend.position = "top")

### Save the cnet plot for the 28 genes ###
oh.save.pdf(p = p_cnet,outname = "rnaseq_cnet_28eth.pdf",
            outdir = "./figures/",width = 12,height = 10)
dev.off()

#Draw heatmap of pattern
Tab <- readRDS(file = "./cleandata/dat_tab_av_rnaseq.RDS")
Tab_sub <- which(rownames(Tab) %in% ids_ppp) %>% 
  Tab[.,]

melted <- Tab_sub %>% melt %>%
  dplyr::rename(.data = .,gene_id = Var1,GenotypeStrain = Var2)

melted$Strain <- melted$GenotypeStrain %>%
  gsub(pattern = ".*_",replacement = "") %>%
  factor %>% relevel(ref = "NB")

melted$Genotype <- melted$GenotypeStrain %>%
  gsub(pattern = "_NB.*",replacement = "") %>%
  gsub(pattern = "_RMF.*",replacement = "") %>%
  gsub(pattern = "_RCL.*",replacement = "") %>%
  gsub(pattern = "_L.*",replacement = "")  %>%
  factor %>%
  relevel(ref = "Col0")

melted_plot <- dcast(data = melted,formula = Strain+Genotype~gene_id,
      fun.aggregate = function(x)mean(x,na.rm = TRUE),value.var = "value",fill = 0) %>%
  melt(id.vars = 1:2)

melted_plot$GenotypeStrain <- paste0(melted_plot$Genotype,"|",melted_plot$Strain)


melted_plot$StrainColor <- melted_plot$Strain %>% as.character

melted_plot$StrainColor[which(!(melted_plot$Strain %in% "NB"))] <- "Other"

### Order individual plot by genotype based on the pattern of expression ###
mgenos <- melted_plot$Genotype %>% as.character %>% unique

list_plots <-  list()
for(mgeno in mgenos){
  
  order_strains <- melted_plot %>%
    subset(Genotype == mgeno) %>%
    aggregate(value~Strain,.,mean) %>%
    dplyr::arrange(.data = .,value) %$% Strain %>% as.character
  
  
  melted_plot_sub <- melted_plot %>%
    subset(Genotype == mgeno ) %>% 
    droplevels
  
  melted_plot_sub$Strain <- melted_plot_sub$Strain %>% factor(levels = order_strains) %>%
    relevel(ref = "NB")
  
  
  p_sub <- ggplot(data = melted_plot_sub,aes(Strain,value)) +
    geom_sina(shape = 16,aes(color = StrainColor)) +
    stat_summary(fun.data = mean_cl_normal,geom = "pointrange",aes(color = StrainColor)) +
    facet_wrap(facets = "Genotype",scales = "free",nrow = 1) +
    theme_ohchibi(size_panel_border = 0.3) +
    ylab(label = "Standardized expression") +
    theme(
      axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "none"
    ) +
    scale_color_manual(values = c("#5AE5ED","#FFA86E"))
  
  list_plots[[mgeno]] <- p_sub
}


composition <- egg::ggarrange(list_plots$Col0,list_plots$lbd16,list_plots$nph4,list_plots$arf7_19,list_plots$gnom1,nrow = 1)
dev.off()
oh.save.pdf(p = composition,outname = "rnaseq_boxplot_genostrain_28th.pdf",
            outdir = "./figures/",width = 26,height = 6)


### Perform the inner line plot in Col-0
df_rep <- melted_plot %>%
  subset(Genotype == "Col0") %>%
  dcast(data = .,formula = StrainColor~variable,fun.aggregate = mean,value.var = "value") %>%
  melt 

df_rep$tair_symbol <- match(df_rep$variable,Map_genes$tair_locus)  %>%
  Map_genes$tair_symbol[.]
df_rep$IsTF <- match(df_rep$variable,Map_genes$tair_locus)  %>%
  Map_genes$IsTF[.]

mindices <- df_rep$tair_symbol %>% grep(pattern = "^$")
df_rep$tair_symbol[mindices] <- df_rep$variable[mindices] %>% as.character

aggregate(value~Treatment,df_rep,mean)

df_rep$Treatment <- df_rep$StrainColor %>%
  gsub(pattern = "Other",replacement = "Average\nbacterial effect") %>%
  factor %>% relevel(ref = "NB")

p <- ggplot(df_rep,aes(Treatment,value)) +
  geom_point(aes(color = IsTF)) +
  geom_text_repel(data = df_rep ,aes(label =tair_symbol,color = IsTF ),force = TRUE,max.overlaps = 100) +
  geom_line(aes(group = variable,color = IsTF),alpha = 0.3) +
  theme_ohchibi(size_panel_border = 0.3) +
  ylab("Standardized expression in Col-0") +
  xlab(label = "Treatment")  +
  scale_color_manual(values = c("red"),na.value = "black") +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank()
  )
  
oh.save.pdf(p = p,outname = "rnaseq_pointlineplot_eth28.pdf",
            outdir = "./figures/",width = 8,height = 14)

rm(list=ls())
dev.off()
gc()
