library(ohchibi)
library(Rmisc)
library(car)

set.seed(130816)


Dat_all <- readRDS(file = "./cleandata/dat_rootbranching_screening1.RDS")


###### Z score values ###########
Dat_z <- Dat_all$Dat_log10_zscore
Map <- Dat_z$Map
melted <- Dat_z$Tab %>% melt
colnames(melted) <- c("variable","root","value")
melted <- merge(melted,Map, by = "root")
### Test
mvariables <- melted$variable %>% as.character %>% unique
mstrains <- melted$Strain %>% as.character %>% unique %>%
  grep(pattern = "NB",value = T,invert = T)

Res_model <- NULL
for(mvar in mvariables){
  melted_var <- melted %>% subset(variable == mvar) %>% droplevels
  for(st in mstrains){
    melted_strain <- melted_var %>% subset(Strain == st)
    x <- melted_strain$value
    mbatches <- melted_strain$Batch %>% as.character %>% unique
    y <- which(melted_var$Batch %in% mbatches) %>%
      melted_var[.,] %>% subset(Strain == "NB") %$% value
    df_temp <- data.frame(Strain = st, value = x) %>%
      rbind(.,data.frame(Strain = "NB",value = y))
    df_temp$Strain <- df_temp$Strain %>% factor %>% relevel(ref = "NB")
    mres <- car::leveneTest(value~Strain,df_temp)
    if(is.nan(mres$`Pr(>F)`[1])){
      m1$p.value <- 1
      
    }else if(mres$`Pr(>F)`[1]< 0.05){
      m1 <- t.test(x,y,var.equal = FALSE)
      
    }else{
      m1 <- t.test(x,y,var.equal = TRUE)
    }
    
    res_man <- wilcox.test(x,y,alternative = "two.sided")
    mest <- mean(x) - mean(y)
    Res_model <- data.frame(Strain = st,variable = mvar,
                            pvalueTTest =m1$p.value,pvalueMW =  res_man$p.value,
                            estimate = mest) %>%
      rbind(Res_model,.) 
  }
}

colnames(Res_model)[2] <- "Feature"

#Normalize p value within variable
Res_end <- NULL
for(mvar in mvariables){
  
  Res_temp <-Res_model %>% subset(Feature ==mvar) %>% droplevels
  Res_temp <- Res_temp %>% dplyr::mutate(.data = .,padjFDRTTest = pvalueTTest %>% p.adjust(method = "fdr"),
                                         padjFDRMW = pvalueMW  %>% p.adjust(method = "fdr")) 
  Res_end <- rbind(Res_end,Res_temp)
}

Res_model <- Res_end

Res_model$padjFDRTTestG <- Res_model$pvalueTTest %>% p.adjust(method = "fdr") 
Res_model$padjFDRMWG <- Res_model$pvalueMW %>% p.adjust(method = "fdr") 

Res_model_zscore <- Res_model
rm(Res_model)



###Test using the raw values
Dat_z <- Dat_all$Dat_log10
Map <- Dat_z$Map
melted <- Dat_z$Tab %>% melt
colnames(melted) <- c("variable","root","value")
melted <- merge(melted,Map, by = "root")
### Test
mvariables <- melted$variable %>% as.character %>% unique
mstrains <- melted$Strain %>% as.character %>% unique %>%
  grep(pattern = "NB",value = T,invert = T)

Res_model <- NULL
for(mvar in mvariables){
  melted_var <- melted %>% subset(variable == mvar) %>% droplevels
  for(st in mstrains){
    melted_strain <- melted_var %>% subset(Strain == st)
    x <- melted_strain$value
    mbatches <- melted_strain$Batch %>% as.character %>% unique
    y <- which(melted_var$Batch %in% mbatches) %>%
      melted_var[.,] %>% subset(Strain == "NB") %$% value
    df_temp <- data.frame(Strain = st, value = x) %>%
      rbind(.,data.frame(Strain = "NB",value = y))
    df_temp$Strain <- df_temp$Strain %>% factor %>% relevel(ref = "NB")
    mres <- car::leveneTest(value~Strain,df_temp)
    if(is.nan(mres$`Pr(>F)`[1])){
      m1$p.value <- 1
      
    }else if(mres$`Pr(>F)`[1]< 0.05){
      m1 <- t.test(x,y,var.equal = FALSE)
      
    }else{
      m1 <- t.test(x,y,var.equal = TRUE)
    }
    
    res_man <- wilcox.test(x,y,alternative = "two.sided")
    
    mest <- mean(x) - mean(y)
    Res_model <- data.frame(Strain = st,variable = mvar,
                            pvalueTTest =m1$p.value,pvalueMW =  res_man$p.value,
                            estimate = mest) %>%
      rbind(Res_model,.)
  }
}

colnames(Res_model)[2] <- "Feature"

#Normalize p value within variable
Res_end <- NULL
for(mvar in mvariables){
  
  Res_temp <-Res_model %>% subset(Feature ==mvar) %>% droplevels
  Res_temp <- Res_temp %>% dplyr::mutate(.data = .,padjFDRTTest = pvalueTTest %>% p.adjust(method = "fdr"),
                                         padjFDRMW = pvalueMW  %>% p.adjust(method = "fdr")) 
  Res_end <- rbind(Res_end,Res_temp)
}

Res_model <- Res_end

Res_model$padjFDRTTestG <- Res_model$pvalueTTest %>% p.adjust(method = "fdr") 
Res_model$padjFDRMWG <- Res_model$pvalueMW %>% p.adjust(method = "fdr") 

Res_model_raw <- Res_model
rm(Res_model)

#Save the detaframes
mlist <- list(
  Res_model_raw = Res_model_raw,
  Res_model_zscore = Res_model_zscore
)

saveRDS(object = mlist,file = "./cleandata/res_rootbranching_screening1_vsNB.RDS")


rm(list=ls())
gc()
