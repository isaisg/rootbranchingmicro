library(ohchibi)

set.seed(130816)


#### Calculate coefficient of variation 
melted <- readRDS(file = "./cleandata/dat_rootbranching_screening1.RDS")  %$% Dat_log10 %$% Tab %>% melt %>%
  dplyr::rename(.data = .,Feature = Var1,root = Var2) %>%
  merge(readRDS(file = "./cleandata/dat_rootbranching_screening1.RDS")  %$% Dat_log10 %$% Map, by = "root")

### Raw values to calculate coefficient of variation
merged <- dcast(data = melted,formula = Feature~Strain,
                fun.aggregate = function(x)mean(x,na.rm = T),value.var = "value") %>% melt %>%
  dplyr::rename(.data = .,Strain = variable)


mfeatures <- merged$Feature %>% as.character %>% unique

#Calculate CV using different levels
#https://clevertap.com/blog/how-to-compare-apples-and-oranges-part-i/
Res_cv <- NULL
for(mfeature in mfeatures){
  
  x <- merged %>%
    subset(Feature == mfeature)  %>% droplevels %$% value
  
  cv_raw <- sd(x)/mean(x)
  
  top <- x %>% quantile(0.9)
  bottom <- x %>% quantile(0.1)
  x_sub <- x[which((x >bottom ) & (x < top)) ]
  cv_90 <- sd(x_sub)/mean(x_sub)
  
  top <- x %>% quantile(0.75)
  bottom <- x %>% quantile(0.25)
  x_sub <- x[which((x >bottom ) & (x < top)) ]
  cv_50 <- sd(x_sub)/mean(x_sub)
  

  Res_cv <- data.frame(Feature = mfeature,CV_Raw = cv_raw,CV_90 = cv_90,CV_50 = cv_50) %>%
    rbind(Res_cv,.)
  
}


### Calculate the phylogetic signal
tree <- readRDS(file = "./cleandata/dat_rootbranching_screening1.RDS") %$% tree

df_res <- readRDS(file = "./cleandata/res_rootbranching_screening1_vsNB.RDS") %$% Res_model_zscore %>%
  subset(Feature != "Norm_child_density") %>% droplevels %>%
  dplyr::filter(.data = .,Strain %in% tree$tip.label) %>% droplevels

Res_psig_estimate <- NULL

mfeatures <- mfeatures %>% grep(pattern = "Norm_child_density",invert = T,value = T)
for(mfeature in mfeatures){
 
  df_res_sub <- df_res %>% subset(Feature == mfeature) %>% droplevels
  df_res_sub <- match(tree$tip.label,df_res_sub$Strain) %>%
    df_res_sub[.,]
  x <- df_res_sub$estimate
  names(x) <- df_res_sub$Strain
  res_psig <- phylosig(tree = tree,x = x,method = "lambda",test = T)
  Res_psig_estimate <- data.frame(Feature = mfeature,pvalue = res_psig$P,lambda = res_psig$lambda) %>%
    rbind(Res_psig_estimate,.)
}

#Calculate lambda using raw data
Res_psig_raw <- NULL
for(mfeature in mfeatures){
  
  df_res_sub <- merged %>% subset(Feature == mfeature) %>% 
    dplyr::filter(.data = .,Strain %in% tree$tip.label) %>% droplevels
  df_res_sub <- match(tree$tip.label,df_res_sub$Strain) %>%
    df_res_sub[.,]
  x <- df_res_sub$value
  names(x) <- df_res_sub$Strain
  res_psig <- phylosig(tree = tree,x = x,method = "lambda",test = T)
  Res_psig_raw <- data.frame(Feature = mfeature,pvalue = res_psig$P,lambda = res_psig$lambda) %>%
    rbind(Res_psig_raw,.)
}




### Merge structure
merged_res <- merge(Res_cv,Res_psig_estimate,by = "Feature") %>%
  merge(Res_psig_raw, by = "Feature") %>%
  dplyr::rename(.data = .,Phylosig_pvalue_estimate = pvalue.x,Phylosig_lambda_estimate = lambda.x,
                Phylosig_pvalue_rawvalue = pvalue.y,Phylosig_lambda_rawvalue = lambda.y)

saveRDS(object = merged_res,file = "./cleandata/dat_cv_phylosig_rootbranching_screening1.RDS")
rm(list=ls())
gc()
