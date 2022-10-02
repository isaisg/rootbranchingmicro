library(ohchibi)

set.seed(130816)


df <- read.csv(file = "./rawdata/rawdata_qpcr_iaa.csv")

df <- df[,colnames(df) %>%grep(pattern = "StD",invert = T)] %>%
  dplyr::select(.data = .,-ID)


df$Treatment <- df$Sample %>%
  gsub(pattern = "L_",replacement = "L") %>%
  gsub(pattern = " ",replacement = "_") %>%
  gsub(pattern = "MF",replacement = "RMF") %>%
  gsub(pattern = "CL",replacement = "RCL") %>%
  gsub(pattern = "_.*",replacement = "") %>%
  factor(levels = c("NB","IAA","L384","L404","L412","L460","L469","RMF27"))

df$Time <- df$Sample %>%
  gsub(pattern = "L_",replacement = "L") %>%
  gsub(pattern = " ",replacement = "_") %>%
  gsub(pattern = "MF",replacement = "RMF") %>%
  gsub(pattern = "CL",replacement = "RCL") %>%
  gsub(pattern = ".*_",replacement = "") %>% 
  factor(levels = c("T0","D2","D3","D4","D7"))

df$Sample <- df$Sample %>% make.unique()

#Scale
Tab <- df[,2:13]
rownames(Tab) <- df$Sample

Tab_z <- scale(x = Tab,center = FALSE,scale = TRUE)

melted <- melt(data = Tab_z) %>%
  dplyr::rename(.data = .,Sample = Var1,Gene = Var2,zscore = value) %>%
  merge(df[,c("Sample","Treatment","Time")], by = "Sample") 


p <- ggplot(data = melted,aes(Treatment,zscore))+
  geom_jitter(aes(color = Gene),width = 0.2,alpha = 0.5,shape = 16) +
  scale_color_paletteer_d("RColorBrewer::Paired") +
  stat_summary(fun = median ,geom = "point",fill = "black",shape = 23,size = 3.5) +
  facet_grid(.~Time,space = "free",scales = "free") +
  theme_ohchibi(size_panel_border = 0.3) +
  theme(
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  ylab(label = "Standardized relative expression")
#

#Save figure and go to sleep
oh.save.pdf(p = p,outname = "rnaseq_qpcr_iaa.pdf",outdir = "./figures/",width = 8,height = 6)

rm(list=ls())
dev.off()
gc()
