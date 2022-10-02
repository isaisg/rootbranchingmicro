library(ohchibi)
library(xlsx)
library(emmeans)
library(multcomp)

set.seed(130816)


### Read Ionome data ####
df_ion <- read.xlsx(file = "./rawdata/All_fig2_data/6_Ionome data with LR phenotypes-1/2022.03.30 Mathieu_Athaliana_Summary_nutri_cleaned.xlsx",sheetIndex = 3)
colnames(df_ion) <- colnames(df_ion) %>%
  gsub(pattern = "\\.[0-9]+$",replacement = "") %>%
  gsub(pattern = "Bacteria",replacement = "Strain")
df_ion <- df_ion[-1,]

df_ion <- which(!(is.na(df_ion$Strain))) %>% 
  df_ion[.,] %>% droplevels

df_ion <- df_ion %>%
  dplyr::filter(.data = .,Strain %in% c("NB","MF27")) %>%
  droplevels

#Create Map
Map <- df_ion[,1:6] %>% 
  dplyr::mutate(.data = .,SampleId = paste0("Sample",1:nrow(df_ion))) %>%
  dplyr::relocate(.data = .,SampleId)
rownames(Map) <- Map$SampleId

Map <- Map %>%
  dplyr::rename(.data = .,Nutrient = Meduim)


Map$Nutrient <- Map$Nutrient %>%
  gsub(pattern = "\\/",replacement = "")

Map$Genotype <- Map$Genotype %>% 
  factor %>% relevel(ref = "Col 0")


Map$Strain <- Map$Strain %>%
  gsub(pattern = "MF27",replacement = "RMF27") %>%
  factor %>%
  relevel(ref = "NB")

Tab <- df_ion[,7:ncol(df_ion)] %>%
  as.matrix
mode(Tab) <- "numeric"
rownames(Tab) <- Map$SampleId

Tab[is.na(Tab)] <- 0

#Remove outliers ###
Tab_norm <- Tab
for(i in 1:ncol(Tab)){
  
  
  max_val <- quantile(Tab[,i],0.995,na.rm = T)
  min_val <- quantile(Tab[,i],0.005,na.rm = T)
  
  x <- Tab[,i]
  
  x[x > max_val] <- max_val
  x[x < min_val] <- min_val
  
  Tab_norm[,i] <- x
  
}


Tab_log <- log10(Tab_norm+1)  %>%
  scale(center = TRUE,scale = TRUE)

Dat <- create_dataset(Tab = Tab_log %>% t,Map = Map)

### Read root data ####
mpath <- "./rawdata/All_fig2_data/6_Ionome data with LR phenotypes-1/"
mfiles <- list.files(path = mpath,pattern = "csv",full.names = FALSE,recursive = TRUE)

df <- NULL
for(mfile in mfiles){
  
  
  mbatch <- "Batch1"
  
  infile <- paste0(mpath,mfile)
  df <- (read.csv(file = infile))[,1:16] %>%
    dplyr::mutate(.data = .,
                  File = mfile,Batch = mbatch) %>%
    rbind(df,.)
  
}



#Create structure of naming
df$File %>%
  gsub(pattern = ".*\\/",replacement = "") %>%
  gsub(pattern = ".*Plate_",replacement = "") %>%
  gsub(pattern = "\\.csv",replacement = "")  %>%
  gsub(pattern = "^[0-9]+_",replacement = "") %>%
  gsub(pattern = "22022203Plate86_",replacement = "") %>%
  gsub(pattern = "22022208Plate87_",replacement = "") %>%
  gsub(pattern = "22022214Plate88_",replacement = "") %>%
  gsub(pattern = " .*",replacement = "") %>% unique


df$Genotype <- NA
df$Genotype[df$File %>% grep(pattern = "Col")] <- "Col0"
df$Genotype[df$File %>% grep(pattern = "nph4")] <- "nph4"
df$Genotype[df$File %>% grep(pattern = "arf7")] <- "arf7"
df$Genotype[df$File %>% grep(pattern = "lbd16")] <- "lbd16"
df$Genotype[df$File %>% grep(pattern = "gnom1")] <- "gnom1"

df$Strain <- NA
df$Strain[df$File %>% grep(pattern = "NB")] <- "NB"
df$Strain[df$File %>% grep(pattern = "MF27")] <- "RMF27"
df$Strain[df$File %>% grep(pattern = "L344")] <- "L344"
df$Strain[df$File %>% grep(pattern = "L359")] <- "L359"

df$Nutrient <- df$File %>%
  gsub(pattern = ".*\\/",replacement = "") %>%
  gsub(pattern = ".*Plate_",replacement = "") %>%
  gsub(pattern = "\\.csv",replacement = "")  %>%
  gsub(pattern = "^[0-9]+_",replacement = "") %>%
  gsub(pattern = "22022203Plate86_",replacement = "") %>%
  gsub(pattern = "22022208Plate87_",replacement = "") %>%
  gsub(pattern = "22022214Plate88_",replacement = "")  %>%
  gsub(pattern = " .*",replacement = "") %>%
  gsub(pattern = ".*_",replacement = "") %>% factor

df <- df %>%
  dplyr::filter(.data = .,Strain %in% c("NB","RMF27")) %>%
  droplevels

#Now we need to perform algorithm to determine the number of childs within each 
##Create unique id for root
df$RootId <- paste0(df$root,"_",df$Batch,"_",df$Strain)

df$root_ontology <- df$root_ontology %>%
  gsub(pattern = " ",replacement = "")

df$root <- df$root %>%
  gsub(pattern = " ",replacement = "")

#Remove weird blank spaces in some measurements
df$root_name <- df$root_name %>%
  gsub(pattern = " ",replacement = "")

df$parent_name <- df$parent_name %>%
  gsub(pattern = " ",replacement = "")

df$parent <- df$parent %>%
  gsub(pattern = " ",replacement = "")

#Conver measurements into numeric 
df$length <- df$length %>% as.numeric
df$vector_length <- df$vector_length %>% as.numeric
df$surface <- df$surface %>% as.numeric
df$volume <- df$volume %>% as.numeric
df$direction <- df$direction %>% as.numeric
df$diameter <- df$diameter %>% as.numeric
df$insertion_position <- df$insertion_position %>% as.numeric
df$insertion_angle <- df$insertion_angle %>% as.numeric
df$n_child <- df$n_child %>% as.numeric


#Ask about the Root ontology type Root that is weird
#Check repeated ids
df_pr <- df %>% subset(root_ontology != "Lateralroot") %>% droplevels
repeated_pr <- df_pr$root %>% table %>%
  data.frame %>% subset(Freq > 1) %$% . %>% as.character

df_lr <- df %>% subset(root_ontology == "Lateralroot") %>% droplevels
repeated_lr <- df_lr$root %>% table %>%
  data.frame %>% subset(Freq > 1) %$% . %>% as.character

length(repeated_pr)
length(repeated_lr)

###Control part to know how many plants we have per strain
Map_count <- df_pr[,c("Genotype","Nutrient","Strain","root","Batch")]  %>%
  unique
Map_count$Count <- 1
Map_count$Strain %>% table %>% sort() %>% min
Map_bc <- Map_count[,c("Genotype","Nutrient","Strain","Batch","Count")] %>%unique
df_num_images <- Map_bc$Strain %>% table %>% sort(decreasing = T) %>%
  data.frame

##Counts ###
df_b <- df_pr[,c("root","length","n_child")] %>%
  dplyr::rename(.data = .,PR_number = root,
                PR_length = length,
                LR_number = n_child)
df_b$LR_Density <- df_b$LR_number/df_b$PR_length


df_b$File <- match(df_b$PR_number,df$root) %>%
  df$File[.]

#Batch
df_b$Batch <- match(df_b$PR_number,df$root) %>%
  df$Batch[.]

#Append information about  Strain and Genotype
df_b$Strain <- match(df_b$PR_number,df$root) %>%
  df$Strain[.]

df_b$Genotype <- match(df_b$PR_number,df$root) %>%
  df$Genotype[.]

df_b$Nutrient <- match(df_b$PR_number,df$root) %>%
  df$Nutrient[.]


df_csv <- df_b

rm(df_b)

#Save ionome data
mlist <- list(
  Dat_ionome = Dat,
  df_root = df_csv
)

saveRDS(object = mlist,file = "./cleandata/dat_ionome_lr_mutants_nutrients.RDS")


rm(list=ls())
gc()
