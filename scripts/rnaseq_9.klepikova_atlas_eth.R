library(ohchibi)
library(dplyr)



set.seed(130816)

mfiles <- list.files("./rawdata/klepikova_data/",full.names = TRUE)

Res_klepikova <- NULL
for(mfile in mfiles){
  
  mid <- mfile %>%
    gsub(pattern = ".*_",replacement = "") %>%
    gsub(pattern = "\\..*",replacement = "")
  
  Res_klepikova <- read.table(file = mfile,
             header = TRUE,sep = "\t",comment.char = "",quote = "") %>%
    dplyr::mutate(.data = .,gene_id =mid) %>%
    rbind(Res_klepikova,.)
  
}

Res_klepikova$Tissue <- Res_klepikova$Tissue %>%
  gsub(pattern = " ",replacement = "_")


#Append genes
Map_genes <- readRDS(file = "./cleandata/res_rnaseq_combinations_contrasts_intracol_intramutant_corlrdensity.RDS") %$% map_genes

Res_klepikova$tair_symbol <- match(Res_klepikova$gene_id,Map_genes$tair_locus) %>%
  Map_genes$tair_symbol[.]

mindices <- Res_klepikova$tair_symbol %>% grep(pattern = "^$")
Res_klepikova$tair_symbol[mindices] <- Res_klepikova$gene_id[mindices]

Tab <- acast(data = Res_klepikova,
             formula = Tissue~tair_symbol,
             fill = 0,value.var = "Expression.Level"
)

#Klepikova enrichment ####
res_cor <- chibi.ggcor(Tab = Tab,
            display.values = TRUE,
            display.significance = FALSE)


df_cor <- res_cor$df_cor %>%
  dplyr::rename(.data = .,Feature = Var1,Feature2 = Var2)

df_cor$Feature <- df_cor$Feature %>% factor(levels = df_cor$Feature %>% levels %>% rev)

#Plot it
p1 <- ggplot(data = df_cor, aes(Feature, Feature2)) + geom_raster(aes(fill = r)) + 
  scale_fill_gradientn(colours = c("#51F5C9","white","#F55C7A"),limits = c(-1,1)) +
  geom_tile(aes(color = Significance), fill = "#00000000", 
            size = 0.3, width = 0.85, height = 0.85) +
  #geom_text(aes(label = rbold),parse = TRUE,family = "Arial",size = 3) +
  scale_color_manual(values = c("#00000000", "black")) +
  theme_ohchibi(legend_proportion_size = 0.5, 
                size_panel_border = 0.3) + 
  theme(legend.position = "top", 
        axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        strip.text.y = element_text(size = 9,angle = 0)) 


##Get the information per gene about the most expressed genes
Res_tissue <- NULL
for(i in 1:ncol(Tab)){
  
  x <- Tab[,i] %>% sort(decreasing = TRUE) %>% head(n = 3) %>% names
  Res_tissue <- data.frame(Gene = colnames(Tab)[i],Tissues = x,Place = 1:3,Dummy = 1) %>%
    rbind(Res_tissue,.)
  
}


Res_tissue$Gene <- Res_tissue$Gene %>% factor(levels = df_cor$Feature2 %>% levels)
chosen_tissue <- Res_tissue$Tissues %>% table %>% sort(decreasing = TRUE) %>%
  data.frame %>% 
  subset(Freq > 2) %$% . %>% as.character

Res_tissue$Tissues[which(!(Res_tissue$Tissues %in% chosen_tissue))] <- "Other"

Res_tissue$Tissues <- Res_tissue$Tissues %>% 
  factor(levels = c("Germinating_seeds_1","Germinating_seeds_2","Germinating_seeds_3",
                    "Seedling_Root",
                    "Root_Apex","Root_without_apex",
                    "Internode","Other"))
paleta <- c(paletteer_d("RColorBrewer::Accent",n = 7),"#D9D9D9")

p2 <- ggplot(data = Res_tissue,aes(Gene,Place)) +
  geom_tile(aes(fill = Tissues),color = "white") +
  coord_flip() +
  ylab(label = "Top 3 Tissues") +
  theme_ohchibi(size_panel_border = 0.3) +
  scale_y_discrete(expand = c(0,0))   +
  scale_fill_manual(values = paleta) +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank()
  )
  
composition <- egg::ggarrange(p1,p2,nrow = 1,widths = c(1,0.1))

#Save the correlation figure with values
oh.save.pdf(p = composition,outname = "rnaseq_klepikova_28eth.pdf",
            outdir = "./figures/",width = 14,height = 9)

rm(list=ls())
gc()
