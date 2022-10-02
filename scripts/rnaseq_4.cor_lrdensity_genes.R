library(ohchibi)
library(dplyr)


set.seed(130816)

#Read datasets
melted_tab <- readRDS(file = "./cleandata/dat_tab_av_rnaseq.RDS") %>% melt %>%
  dplyr::rename(.data = .,Gene = Var1,group = Var2)  %>%
  dplyr::filter(.data = .,group %>% grepl(pattern = "Col0")) %>% droplevels %>%
  dplyr::mutate(.data = .,Strain = group %>% gsub(pattern = "Col0_",replacement = ""))


df_pheno <- readRDS(file = "./cleandata/dat_rootbranching_screeningrnaseq.RDS")


df_av_gene <- dcast(data = melted_tab,formula = Strain~Gene,
      fun.aggregate = function(x)mean(x,na.rm = TRUE),
      value.var = "value") %>% melt %>%
  dplyr::rename(.data = .,Gene = variable,zscore = value)

df_lr <- df_pheno %>% 
  subset(Feature == "LR_density") %>% 
  aggregate(value~Strain,.,mean)
  

df_av_gene <- merge(df_av_gene,df_lr, by = "Strain") 

#Compute correlation 
mgenes <- df_av_gene$Gene %>% as.character %>% unique

Res_cor <- NULL
Res_measurements <- NULL

for(mgene in mgenes){
  
  cat("Working on",mgene,"\n")
  df_sub <- df_av_gene %>% subset(Gene == mgene) %>% droplevels
  m1 <- cor.test(df_sub$zscore,df_sub$value,method = "spearman")
  m2 <- cor.test(df_sub$zscore,df_sub$value,method = "pearson")
  
  Res_cor <- data.frame(Gene = mgene,
             rho = m1$estimate %>% as.numeric,
             p_spearman = m1$p.value %>% as.numeric,
             r = m2$estimate %>% as.numeric,
             p_pearson = m2$p.value %>% as.numeric) %>%
    rbind(Res_cor,.)
  
  rownames(df_sub) <- NULL
  df_sub <- df_sub %>%
    dplyr::rename(.data = .,GeneExpressionZ = zscore,LRDensity = value)
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  
  x <- range01(x = df_sub$GeneExpression)
  df_sub$GeneExpressionRange <- x
  df_sub <- df_sub %>% 
    dplyr::select(.data = .,c("Gene","Strain","GeneExpressionZ","GeneExpressionRange","LRDensity"))
  
  Res_measurements <- rbind(Res_measurements,df_sub)
   
}

mlist <- list(
  Res_cor = Res_cor,
  Res_measurements = Res_measurements
)

#
#Save results
saveRDS(object = mlist,file = "./cleandata/res_cor_rnaseq_lrdensity.RDS")

rm(list=ls())
dev.off()
gc()