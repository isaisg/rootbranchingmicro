library(ohchibi)
library(dplyr)
library(ComplexHeatmap)
library(biomaRt)
library(org.At.tair.db)
library(clusterProfiler)
library(ggtree)
library(xlsx)
library(ggalluvial)
library(ggforce)

set.seed(130816)


Res <- readRDS(file = "./cleandata/res_rnaseq_combinations_contrasts_intracol_intramutant_corlrdensity.RDS")

df_comb <- Res$df_comb_gec

#Determine the ones in col
a <- df_comb$EnrichmentCode %>% grep(pattern = "^001")  %>%
  df_comb[.,] %$% gene_id %>% as.character %>% unique 

b <- df_comb$EnrichmentCode %>% grep(pattern = "^100")  %>%
  df_comb[.,] %$% gene_id %>% as.character %>% unique 

union(a,b) %>%  length

##
a <- df_comb$EnrichmentCode %>% grep(pattern = "^001001")  %>%
  df_comb[.,] %$% gene_id %>% as.character %>% unique 

b <- df_comb$EnrichmentCode %>% grep(pattern = "^100100")  %>%
  df_comb[.,] %$% gene_id %>% as.character %>% unique 

union(a,b) %>%  length

mpairs <- df_comb$BacteriaGenotype  %>% unique

df_alluvial <- NULL

for(mpair in mpairs){
  
  #mpair <- "RMF27_gnom1"
  #mpair <- "L343_lbd16"
  
  df_comb_sub <- df_comb %>% 
    subset(BacteriaGenotype == mpair) %>% droplevels
  
  
  df_freq <- df_comb_sub$EnrichmentCode %>% table %>%
    data.frame %>%
    dplyr::rename(.data = .,EnrichmentCode = .)
  
  for(ec in df_freq$EnrichmentCode){
    
    df_freq_sub  <- df_freq %>% 
      subset(EnrichmentCode == ec) %>% droplevels
    
    
    x <- gsub("(.{3})", "\\1 ", df_freq_sub$EnrichmentCode %>% as.character) %>% 
      strsplit(split = " ") %>% unlist
    
    df_alluvial <- data.frame(Pair = mpair, 
                              EnrichmentCode = df_freq_sub$EnrichmentCode %>% as.character,
                              BactInCol = x[1],
                              BactInMut = x[2],
                              CorLR = x[3],
                              Freq =df_freq_sub$Freq %>% as.numeric ) %>%
      rbind(df_alluvial,.)
    
  }
  
}

df_alluvial$BactInCol <- df_alluvial$BactInCol %>%
  gsub(pattern = "001",replacement = "+") %>%
  gsub(pattern = "010",replacement = "0") %>%
  gsub(pattern = "100",replacement = "-") %>%
  factor(levels = c("+","0","-"))

df_alluvial$BactInMut <- df_alluvial$BactInMut %>%
  gsub(pattern = "001",replacement = "+") %>%
  gsub(pattern = "010",replacement = "0") %>%
  gsub(pattern = "100",replacement = "-") %>%
  factor(levels = c("+","0","-"))


df_alluvial$CorLR <- df_alluvial$CorLR %>%
  gsub(pattern = "001",replacement = "+") %>%
  gsub(pattern = "010",replacement = "0") %>%
  gsub(pattern = "100",replacement = "-") %>%
  factor(levels = c("+","0","-"))

df_alluvial$EnrichmentCodeSymbol <- paste0(df_alluvial$BactInCol,df_alluvial$BactInMut,df_alluvial$CorLR)

df_alluvial <- na.omit(df_alluvial)

#Determine sum of counts for code
df_sum_freq <- aggregate(Freq~EnrichmentCodeSymbol,df_alluvial,sum) %>%
  dplyr::arrange(.data = .,-Freq)

#chosen_symbols <- df_sum_freq %>%
#  subset(Freq >=100) %$%  EnrichmentCodeSymbol 

chosen_symbols <- c("+++","++0","---","--0")

df_alluvial$EnrichmentCodeSymbolFill <- df_alluvial$EnrichmentCodeSymbol

df_alluvial$EnrichmentCodeSymbolFill[which(! df_alluvial$EnrichmentCodeSymbolFill %in% chosen_symbols)] <- "Other"

df_alluvial$EnrichmentCodeSymbolFill <- df_alluvial$EnrichmentCodeSymbolFill %>%
  factor(levels = c(chosen_symbols,"Other"))


paletteer_d("ochRe::lorikeet")

paleta <-c("#C03018FF","#F0A800FF","#484878FF","#609048FF","#D9D9D9")
names(paleta) <- c("+++","++0","---","--0","Other")

df_plot <- df_alluvial %>% 
  subset(BactInCol != 0) %>% droplevels 



#Test alluvial plot
p <- df_alluvial %>% 
  ggplot(data = .,aes(axis1 = BactInCol, axis2 = BactInMut, axis3 = CorLR,y = Freq)) +
  scale_x_discrete(limits = c("BactInCol", "BactInMut", "GenoInNB"), 
                   labels = c("Strain vs NB\n(Col-0)","Strain vs NB\n(Mutant)","Correlation\nLR Density ~ Expression"),
                   expand = c(.2, .05)) +
  xlab("Contrast") +
  geom_alluvium(aes(fill = EnrichmentCodeSymbolFill),alpha = 0.75,width = 1/4) +
  facet_wrap(facets = "Pair",scales = "free_y",nrow = 4,ncol = 4) +
  geom_stratum(width = 1/4) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal()  +
  theme(
    legend.position = "right"
  ) +
  scale_fill_manual(values = paleta,name = "Enrichment\npattern") +
  ylab(label = "Number of DEGs") 

##Try ggforce implementation


df_alluvial$FreqLog10 <- log10(df_alluvial$Freq)

df_alluvial %>% 
  ggplot(data = .,aes(axis1 = BactInCol, axis2 = BactInMut, axis3 = CorLR,y = FreqLog10)) +
  scale_x_discrete(limits = c("BactInCol", "BactInMut", "GenoInNB"), 
                   labels = c("Strain vs NB\n(Col-0)","Strain vs NB\n(Mutant)","Correlation\nLR Density ~ Expression"),
                   expand = c(.2, .05)) +
  xlab("Contrast") +
  geom_alluvium(aes(fill = EnrichmentCodeSymbolFill),alpha = 0.75) +
  facet_wrap(facets = "Pair",scales = "free",nrow = 4,ncol = 4) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal()  +
  theme(
    legend.position = "right"
  ) +
  scale_fill_manual(values = paleta,name = "Enrichment\npattern") +
  ylab(label = "Number of DEGs") 

  
  ##Save figure
oh.save.pdf(p = p,outname = "rnaseq_alluvial_intracol_intramutant_corlr.pdf",outdir = "./figures/",
            width = 14,height = 10)

rm(list=ls())
dev.off()
gc()
