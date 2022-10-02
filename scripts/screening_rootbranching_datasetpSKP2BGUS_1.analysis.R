library(ohchibi)
library(xlsx)
library(emmeans)
library(multcomp)

set.seed(130816)


df <- read.xlsx(file = "./rawdata/All_fig2_data/10_pSKP2BGUS data.xlsx",sheetIndex = 1)

df$Batch.nb <- paste0("Batch",df$Batch.nb) %>% 
  factor

df$Plant.number <- paste0("Plant",df$Plant.number) %>% 
  factor

df$Genotype <- df$Genotype %>%
  factor

df$Bacteria <- df$Bacteria %>%
  gsub(pattern = "_",replacement = "") %>%
  gsub(pattern = " ",replacement = "") %>%
  gsub(pattern = "MF",replacement = "RMF") %>%
  gsub(pattern = "CL",replacement = "RCL") %>%
  factor %>%
  relevel(ref = "NB")

df$Stage <- df$Stage %>%
  factor

#Rename columns
colnames(df)[4:5] <- c("PlantNumber","Batch")

df$Dummy <- 1

#COmpute relative abundance
merged <- acast(data = df,formula = Stage~Genotype+Bacteria+PlantNumber+Batch,fun.aggregate = sum,value.var = "Dummy",fill = 0) %>%
  scale(x = .,center = FALSE,scale = colSums(.))  %>%
  melt %>% 
  dplyr::rename(.data = .,Stage = Var1,UId = Var2,Prop = value)

merged <- merged$UId %>% as.character %>%
  strsplit(split = "_") %>%
  unlist %>%
  matrix(data = .,ncol = 4,byrow = TRUE) %>%
  data.frame %>%
  dplyr::rename(.data = .,Genotype = X1,Bacteria = X2,PlantNumber = X3,Batch = X4) %>%
  cbind(merged,.)

df_mean <- aggregate(Prop~Genotype+Bacteria+Stage+Batch,merged,mean)

df_mean$Stage <- df_mean$Stage %>% factor(levels = df_mean$Stage %>% levels %>% rev)
df_mean$Bacteria <- df_mean$Bacteria %>% factor %>%
  relevel(ref = "NB")

df_mean$Treatment <- paste0(df_mean$Genotype,"_",df_mean$Bacteria)

ggplot(data = df_mean,aes(Bacteria,Prop)) +
  geom_bar(stat = "identity",aes(fill = Stage))+
  facet_grid(Genotype~Batch,space = "free",scales = "free")+
  theme_ohchibi(size_panel_border = 0.3) +
  theme(
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1) 
  ) +
  scale_fill_paletteer_d("ggthemes::colorblind")


#Average values across both batches
df_av <- aggregate(Prop~Treatment+Stage,df_mean,mean)
df_av$Genotype <- match(df_av$Treatment,df_mean$Treatment) %>%
  df_mean$Genotype[.]

df_av$Bacteria <- match(df_av$Treatment,df_mean$Treatment) %>%
  df_mean$Bacteria[.]

#### Subst only ones present in previous analyse
df_av_plot <- df_av %>%
  subset(Genotype != "Col0") %>% droplevels

order_treatment <- df_av_plot %>%
  subset(Stage == "E") %>%
  dplyr::arrange(.data = .,Prop) %$% Treatment  %>%
  grep(pattern = "NB",invert = TRUE,value = TRUE) %>%
  c("SKP2B_NB",.)

df_av_plot$Treatment <- df_av_plot$Treatment %>%
  factor(levels = order_treatment)


p1 <- ggplot(data = df_av_plot,aes(Treatment,Prop)) +
  geom_bar(stat = "identity",aes(fill = Stage))+
  facet_grid(.~Genotype,space = "free",scales = "free")+
  theme_ohchibi(size_panel_border = 0.3) +
  theme(
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1) 
  ) +
  scale_fill_paletteer_d("ggthemes::colorblind") +
scale_y_continuous(breaks = seq(0,1,0.1),expand = c(0,0))


#Perform correlation with previous dataset
df_prev <- readRDS(file = "./cleandata/res_screening_rootbranching_9_pre-branching.RDS")
df_a <- df_prev %>% subset(Genotype == "Col 0")
df_b <- df_av %>% 
  subset(Genotype == "SKP2B") %>% droplevels

merged <- merge(df_a,df_b, by = c("Bacteria","Stage")) 


p2 <- ggplot(data = merged,aes(Prop.x, Prop.y)) +
  geom_point(shape = 16,size = 3) +
  geom_smooth(method = "lm",se = FALSE,color = "red")+
  ggpubr::stat_cor(label.sep = "\n")+
  theme_ohchibi(size_panel_border = 0.3)+
  xlab(label = "Primordia counting (Microscope)") +
  ylab(label = "pSKP2BGUS ")  +
  theme(panel.grid.major.x = element_blank(),
  panel.grid.major.y = element_blank())
  


#Save plots
oh.save.pdf(p = p1,outname = "fig2_stage_pSKP2BGUS_prop.pdf",
            outdir = "./figures/",width = 10,height =8)


oh.save.pdf(p = p2,outname = "fig2_stage_primordia_pSKP2BGUS_cor.pdf",
            outdir = "./figures/",width = 9,height =8)

rm(list=ls())
dev.off()
gc()



