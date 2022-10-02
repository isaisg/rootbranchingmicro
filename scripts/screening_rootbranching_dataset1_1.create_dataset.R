library(ohchibi)

set.seed(130816)


### Read all the files into a structure
mpaths <- list.files(path = "./rawdata/Phenotyping_raw_data_screening_1",
                     recursive = T,pattern = "*.csv",full.names = T)
mpaths <- mpaths %>%
  grep(pattern = "data_screening_iaabiosynthesisdegradation",invert = T,value = T)

mcols <- read.csv(mpaths[1],header = T,
                  comment.char = "",quote = "") %>% colnames
df <- NULL
for(mp in mpaths){
  mflag <- "Yes"
  #Not all files have header
  mdf <- read.csv(mp,header = T,
                  comment.char = "",quote = "")
  mhit <- colnames(mdf) %>% grep(pattern = "root_name")
  if(length(mhit) == 0){
    mflag <- "No"
    mdf <- read.csv(mp,header = F,
                    comment.char = "",quote = "")
    colnames(mdf) <- mcols
  }
  #Check number of columns
  if(length (mdf) < length(mcols)){
    mdf <- read.table(file = mp,header = T,sep = ";",
                      quote = "",comment.char = "")
  }
  if(length (mdf) < length(mcols)){
    stop()
  }
  #Verify the columns have same name
  mnums <- match(mcols,colnames(mdf)) 
  morder <- mnums %>% paste0(collapse = "")
  
  mbatch <- mp %>% gsub(pattern = "./rawdata/Phenotyping_raw_data_screening_1/",replacement = "") %>%
    gsub(pattern = "/.*",replacement = "")
  mnom <- mp %>% gsub(pattern = "./rawdata/Phenotyping_raw_data_screening_1/",replacement = "") %>%
    gsub(pattern = ".*/",replacement = "") %>%
    gsub(pattern = ".csv",replacement = "")
  mdf$Path <-  mp %>% gsub(pattern = "./rawdata/Phenotyping_raw_data_screening_1/",replacement = "")
  mdf$Batch <- mbatch
  mdf$FileName <- mnom 
  mdf$Header <- mflag
  mdf$NumControl <- morder
  df <- rbind(df,mdf)
}
#Verify columns are always in the same order
df$NumControl %>% unique %>% length
df$Header %>% table

#Split the Filename intro Plates and Bacterial treatment
df$Strain <- df$FileName %>%
  gsub(pattern = "L_",replacement = "L") %>%
  gsub(pattern = "MF_",replacement = "MF") %>%
  gsub(pattern = "CL_",replacement = "CL") %>%
  gsub(pattern = " ",replacement = "") %>%
  gsub(pattern = ".*_",replacement = "") %>% 
  gsub(pattern = "!!!!",replacement = "Dob")

df$Strain %>% unique %>% length
#We hve 392 strains


df$FileNameMod <- df$FileName %>%
  gsub(pattern = "L_",replacement = "L") %>%
  gsub(pattern = "MF_",replacement = "MF") %>%
  gsub(pattern = "CL_",replacement = "CL") %>%
  gsub(pattern = " ",replacement = "") %>% 
  gsub(pattern = "!!!!",replacement = "Dob")



#Clean the space in all columns
for(i in 1:ncol(df)){
  df[,i] <- df[,i] %>%
    gsub(pattern = "\\s+",replacement = "",perl = T)
}

###Loop over each Strain and add a number id of the plates
mstrains <- df$Strain %>% unique
df_mod <- NULL
for(mst in mstrains){
  df_temp <- df %>% subset(Strain == mst) %>% droplevels
  df_temp$PlateIdRaw <- df_temp$FileNameMod %>% 
    gsub(pattern = paste0("_",mst),replacement = "") 
  df_mod <- rbind(df_mod,df_temp)
}
df_mod$PlateId <- df_mod$PlateIdRaw %>% factor %>% as.numeric %>%
  paste0("Plate",.)

#Control to know all variables are present
paste0(df_mod$PlateId,"_",df_mod$Strain) %>% unique %>% length
df_mod$FileName %>% unique %>% length

df_mod$Batch <- df_mod$Batch  %>% gsub(pattern = "_",replacement = "") %>% 
  gsub(pattern = "Batch5Warning2platesperphoto",replacement = "Batch5") %>%
  gsub(pattern = "L22",replacement = "Batch23")
df_mod$Batch %>% unique

rm(df)
df <- df_mod
rm(df_mod)


#First remove some weird measurements with solely NA values
dim(df)
df <- df$root_ontology %>% grep(pattern = "^$",invert = T) %>%
  df[.,] %>% droplevels
dim(df)

#With the revised rawdata we dont see those blank spaces anymore

###Convert to numeric the measurement columns that became character
#due to removing white space
df$length <- df$length %>% as.numeric
df$vector_length <- df$vector_length %>% as.numeric
df$surface <- df$surface %>% as.numeric
df$volume <- df$volume %>% as.numeric
df$direction <- df$direction %>% as.numeric
df$diameter <- df$diameter %>% as.numeric
df$insertion_position <- df$insertion_position %>% as.numeric
df$insertion_angle <- df$insertion_angle %>% as.numeric
df$n_child <- df$n_child %>% as.numeric
df$child_density <- df$child_density %>% as.numeric
df$insertion_first_child <- df$insertion_first_child %>% as.numeric
df$insertion_last_child <- df$insertion_last_child %>% as.numeric

###The measurement for the analysis is the individual primary roots ####
df$root_ontology %>% table

#Now we only have Primaryroot or Lateralroot

#Rename strains
df$Strain <- df$Strain %>% 
  gsub(pattern = "^MF",replacement = "RMF") %>%
  gsub(pattern = "CL",replacement = "RCL")

##Create unique id for root
df$RootId <- paste0(df$root,"_",df$Batch,"_",df$PlateId,"_",df$Strain)

#Ask about the Root ontology type Root that is weird
#Check repeated ids
df_pr <- df %>% subset(root_ontology == "Primaryroot") %>% droplevels
repeated_pr <- df_pr$root %>% table %>% 
  data.frame %>% subset(Freq > 1) %$% . %>% as.character

df_lr <- df %>% subset(root_ontology == "Lateralroot") %>% droplevels
repeated_lr <- df_lr$root %>% table %>% 
  data.frame %>% subset(Freq > 1) %$% . %>% as.character

length(repeated_pr)
length(repeated_lr)

#We dont have repeated ids anymore

###Control part to know how many plants we have per strain
Map_count <- df_pr[,c("Strain","root","Batch","PlateId")]  %>%
  unique
Map_count$Count <- 1
Map_count$Strain %>% table %>% sort() %>% min
Map_bc <- Map_count[,c("Strain","Batch","PlateId","Count")] %>%unique
df_num_images <- Map_bc$Strain %>% table %>% sort(decreasing = T) %>%
  data.frame
colnames(df_num_images) <- c("Strain","NumbeOfPlates")

#Read the list of the science strains
df_suberin <- readRDS("./cleandata/dat_suberin_mono_rootbarriers.RDS")
sts_41 <- df_suberin$SuberinScienceSum$Strain %>% 
  grep(pattern = "CDEF",invert = T,value = T)

df_num_images$ScienceStrains <- "No"
which(!(sts_41 %in% df_num_images$Strain)) %>%
  sts_41[.]
df_num_images$ScienceStrains[which(df_num_images$Strain %in% sts_41)] <- "Yes"


######################################################
######################################################
########### Part of the analysis #####################

#Determine the  of  primary roots in file
df_primaryroots <- df$root_ontology %>% grep(pattern = "Primaryroot") %>%
  df[.,] %>% droplevels

##Determine number of lateralroots per primary roots and range
#We are gonna calculate four inner statistics
#Distance between first root and last root
#Proportion of range (dividing by total root length)
#Distance between total root length and last lateral root

#First oder analysis mainly

df_tot_lr <- NULL
for(i in 1:nrow(df_primaryroots)){
  mrow <- df_primaryroots[i,]
  #Length of the primary root
  mlength <- mrow$length
  #Direction
  mdirection <- mrow$direction
  mid <- mrow$root
  df_first <- df %>% subset(parent == mid) %>% droplevels
  ids_first <- df_first$root %>% unique
  
  #Now determine the range
  tot_child_raw <- mrow$n_child
  if(tot_child_raw == 0){
    mrange_last_first <- 0
    mprop_last_first <- 0
    mrange_tip_last <- 0
    mprop_tip_last <- 0
    position_lr_diff <- 0
    interbranch <- 0
  }else if(tot_child_raw == 1){
    left_side <- which(df_first$direction < mdirection)
    right_side <- which(df_first$direction > mdirection)
    mleft <- length(left_side)/length(ids_first)
    mclock <- mleft
    
    mrange_last_first <- 0
    mprop_last_first <- 0
    mrange_tip_last <- mlength - mrow$insertion_last_child
    mprop_tip_last <- mrange_tip_last / mlength
    position_lr_diff <- 0
    interbranch <- 0
    
  }else{
    left_side <- which(df_first$direction < mdirection)
    right_side <- which(df_first$direction > mdirection)
    mleft <- length(left_side)/length(ids_first)
    mclock <- mleft
    
    mrange_last_first <- mrow$insertion_last_child - mrow$insertion_first_child
    mprop_last_first <- mrange_last_first / mlength
    mrange_tip_last <- mlength - mrow$insertion_last_child
    mprop_tip_last <- mrange_tip_last / mlength
    
    df_temp <- df %>% subset(parent == mid) %>% droplevels
    df_temp <- with(df_temp,order(insertion_position)) %>%
      df_temp[.,]
    #Position lr difference between top and down
    #0.15 is the percentage of roots to take from top and bottom 
    mnum <- ceiling(nrow(df_temp)* 0.15)
    a <- head(df_temp,mnum)
    b <- tail(df_temp,mnum)
    #Here we are comparing the differences in length
    position_lr_diff <- (mean(a$length) - mean(b$length)) 
    
    #Calculate interbranch in the meantime
    x <- df_temp$insertion_position
    
    #Normalize the interposition by total length
    #This should make it orthogonal to length
    x <- x/mlength
    
    #Loop over it calculating the intra distance
    res_id <- NULL
    for(j in 1:(length(x)-1)){
      z <- x[j+1] - x[j]
      res_id <- c(res_id,z)
    }
    interbranch <- res_id %>% mean
  }
  df_tot_lr <- data.frame(root = mid,
                          Length = mlength,
                          LR1stOrderNumRoots = length(ids_first),
                          LR1stOrderNumRootsbyCm = length(ids_first) / mlength,
                          LR1stOrderRangeLasttoFirst = mrange_last_first,
                          LR1stOrderPropRangeLasttoFirst = mprop_last_first,
                          LR1stOrderRangeTiptoLast = mrange_tip_last,
                          LR1stOrderPropRangeTiptoLast = mprop_tip_last,
                          LR1stOrderAsymmetryPattern = mclock,
                          LR1stOrderInterBranch = interbranch,
                          LR1stOrderDiffLengthTopBottom = position_lr_diff) %>%
    rbind(df_tot_lr,.)
}


#There should not be PropRange1stLRLasttoFirst > 1
df_tot_lr[which(df_tot_lr$LR1stOrderPropRangeLasttoFirst > 1),]
mindices <- which(df_tot_lr$LR1stOrderPropRangeLasttoFirst > 1)

#Change these values to set up limits for these peculiar measurements
df_tot_lr$LR1stOrderPropRangeLasttoFirst[mindices] <- 1
df_tot_lr$LR1stOrderRangeTiptoLast[mindices] <- 0
df_tot_lr$LR1stOrderPropRangeTiptoLast[mindices] <- 0

df_tot_lr$LR1stOrderRangeLasttoFirst[mindices] <- match(df_tot_lr$root[mindices],df_primaryroots$root) %>%
  df_primaryroots$length[.]


##Remove the obvious outlier
dim(df_tot_lr)
df_tot_lr <- df_tot_lr %>% subset(LR1stOrderNumRootsbyCm < 5000) %>% droplevels
dim(df_tot_lr)

#######################
#######################
#### Lateral roots ####

### First order ###
df_lateralroots <- df$root_ontology %>% grep(pattern = "Lateralroot") %>%
  df[.,] %>% droplevels

##Given that the metric is the primary root we will quantify everything in relation to it
res_lr <- NULL
for(i in 1:nrow(df_primaryroots)){
  mrow <- df_primaryroots[i,]
  mid <- mrow$root
  df_first <- df %>% subset(parent == mid) %>% droplevels
  ids_first <- df_first$root %>% unique
  tot_child_raw <- mrow$n_child
  if(tot_child_raw == 0){
    #All measurements become 0 if there is no 1st order lateral root
    mlength <- 0
    msurface <- 0
    LR1stOrderLengthAv <- 0
    LR1stOrderSurfaceAv <- 0
    LR1stOrderVolumeAv <- 0
    LR1stOrderDiameterAv <- 0
    LR1stInsertionAngleAv <- 0
    LR2ndOrderNumberRaw <- 0
    LR2ndOrderNumberNorm <- 0
    LR1stOrderRAwith2ndOrder <- 0
    LR2ndOrderLengthAv <- 0
    LR2ndOrderSurfaceAv <- 0
    LR2ndOrderVolumeAv <- 0
    LR2ndOrderDiameterAv <- 0
    
  }else if(tot_child_raw == 1){
    #Hwere we cocunt the statistic for all first order primary roots
    mlength <- df_first$length %>% mean
    msurface <- df_first$surface %>% mean
    mvolume <- df_first$volume %>% mean
    mdiameter <- df_first$diameter %>% mean
    minsertionangle <- df_first$insertion_angle %>% mean
    
    #Numer of childs is peculiar
    #We can have the average number of childs
    mnchild_raw <- df_first$n_child  %>% mean
    x <- which(df_first$n_child != 0) %>% df_first$n_child[.]
    ra_withsecond <- length(x)/length(ids_first) 
    #Average number of childs taking only into account the roots that have n_child
    mchild_norm <- 0
    if(length(x) != 0){
      mnchild_norm <- mean(x) 
      
    }
    
    raw_child_density <- df_first$child_density %>% mean
    
    notzero <- which(df_first$child_density != 0)  %>%
      df_first$child_density[.] 
    norm_child_density <- 0
    if(length(notzero) != 0){
      norm_child_density <- notzero %>% mean
    }
    
    #Proceeed with the second order
    df_second <- df %>% subset(parent %in% ids_first) %>% droplevels
    
    #Here evaluate if there are actually any second order lateral root
    if(nrow(df_second) == 0){
      mlength_2 <- 0
      msurface_2 <- 0
      mvolume_2 <- 0
      mdiameter_2 <- 0
    }else{
      mlength_2 <- df_second$length %>% mean
      msurface_2 <- df_second$surface %>% mean
      mvolume_2 <- df_second$volume %>% mean
      mdiameter_2 <- df_second$diameter %>% mean
    }
    
  }else{
    #Hwere we cocunt the statistic for all first order primary roots
    mlength <- df_first$length %>% mean
    msurface <- df_first$surface %>% mean
    mvolume <- df_first$volume %>% mean
    mdiameter <- df_first$diameter %>% mean
    minsertionangle <- df_first$insertion_angle %>% mean(na.rm = TRUE)
    
    #Numer of childs is peculiar
    #We can have the average number of childs
    mnchild_raw <- df_first$n_child  %>% mean
    x <- which(df_first$n_child != 0) %>% df_first$n_child[.]
    ra_withsecond <- length(x)/length(ids_first) 
    #Average number of childs taking only into accoutn the roots that have n_child
    mchild_norm <- 0
    if(length(x) != 0){
      mnchild_norm <- mean(x) 
      
    }
    
    raw_child_density <- df_first$child_density %>% mean
    
    notzero <- which(df_first$child_density != 0)  %>%
      df_first$child_density[.] 
    norm_child_density <- 0
    if(length(notzero) != 0){
      norm_child_density <- notzero %>% mean
    }
    
    #Proceeed with the second order
    df_second <- df %>% subset(parent %in% ids_first) %>% droplevels
    
    #Here evaluate if there are actually any second order lateral root
    if(nrow(df_second) == 0){
      mlength_2 <- 0
      msurface_2 <- 0
      mvolume_2 <- 0
      mdiameter_2 <- 0
    }else{
      mlength_2 <- df_second$length %>% mean
      msurface_2 <- df_second$surface %>% mean
      mvolume_2 <- df_second$volume %>% mean
      mdiameter_2 <- df_second$diameter %>% mean
    }
    
    
  }
  
  #Create dataframe with the results
  res_lr <- data.frame(
    root = mid,
    LR1stOrderLengthAv = mlength,
    LR1stOrderSurfaceAv = msurface,
    LR1stOrderVolumeAv = mvolume,
    LR1stOrderDiameterAv = mdiameter,
    LR1stInsertionAngleAv = minsertionangle,
    LR2ndOrderNumberRaw = mnchild_raw,
    LR2ndOrderNumberNorm = mnchild_norm,
    LR1stOrderRAwith2ndOrder = ra_withsecond,
    LR2ndOrderLengthAv = mlength_2,
    LR2ndOrderSurfaceAv = msurface_2,
    LR2ndOrderVolumeAv = mvolume_2,
    LR2ndOrderDiameterAv = mdiameter_2) %>%
    rbind(res_lr,.)
  
}


#Verification everything is a number
apply(res_lr[,-1] %>% as.matrix,MARGIN = 2,FUN = sum)

dim(df_tot_lr)
dim(res_lr)

torum <- colnames(df_tot_lr) %>% grep(pattern = "^Length")
merged <- merge(df_tot_lr[,-torum],res_lr,by = "root")


df_primaryroots <- merge(df_primaryroots,merged, by = "root")


### Remove the batch effect #####

###Normalize across the primary roots
mvariables <- c("length","surface","volume","diameter","child_density")
mvariables <- colnames(df_primaryroots) %>% grep(pattern = "LR",value = T) %>%
  c(mvariables,.)

mbatches <- df_primaryroots$Batch %>% unique

df_primaryroots_nb <- df_primaryroots%>% 
  subset(Strain == "NB") %>% droplevels

Res_pr <- data.frame(root = df_primaryroots$root)
dim(Res_pr)
Res_nfactors_pr <- NULL
for(mvar in mvariables){
  global_mean_nb <- df_primaryroots_nb[,mvar] %>% mean(na.rm = TRUE)
  df_res_interno <- NULL
  for(mbatch in mbatches){
    df_interno <- df_primaryroots_nb %>%
      subset(Batch == mbatch) %>%
      droplevels
    local_mean_nb <- df_interno[,mvar] %>% mean(na.rm = TRUE)
    nfactor <- global_mean_nb/local_mean_nb
    #Do the normalization
    df_temp <- df_primaryroots %>%
      subset(Batch == mbatch)
    if(is.infinite(nfactor)){
      cat(mvar," Infinite ",mbatch,"\n")
      df_res <- data.frame(root = df_temp$root,
                           mvar = df_temp[,mvar] + global_mean_nb
      )
    }else{
      df_res <- data.frame(root = df_temp$root,
                           mvar = df_temp[,mvar] * nfactor
      )
    }
    
    colnames(df_res)[2] <- paste0("Norm_",mvar)
    df_res_interno <- rbind(df_res_interno,df_res)
    Res_nfactors_pr <- data.frame(Variable = mvar,Batch = mbatch,NFactor = nfactor) %>%
      rbind(Res_nfactors_pr)
  }
  Res_pr <- merge(Res_pr,df_res_interno, by = "root")
}
dim(Res_pr)

#Prepare final dataframe
df <- merge(df_primaryroots,Res_pr) 


### Pepare Dat object to contain the metrics of itnerest

##Variables
mcolumns <- c("Norm_length","Norm_surface",
              "Norm_volume","Norm_diameter",
              "Norm_child_density","Norm_LR1stOrderNumRootsbyCm",
              "Norm_LR1stOrderPropRangeLasttoFirst","Norm_LR1stOrderPropRangeTiptoLast",
              "Norm_LR1stOrderAsymmetryPattern","Norm_LR1stOrderInterBranch",
              "Norm_LR1stOrderDiffLengthTopBottom","Norm_LR1stOrderLengthAv",
              "Norm_LR1stOrderSurfaceAv","Norm_LR1stOrderVolumeAv",
              "Norm_LR1stOrderDiameterAv","Norm_LR1stInsertionAngleAv",
              "Norm_LR1stOrderRAwith2ndOrder",
              "Norm_LR2ndOrderNumberNorm","Norm_LR2ndOrderLengthAv",
              "Norm_LR2ndOrderSurfaceAv","Norm_LR2ndOrderVolumeAv",
              "Norm_LR2ndOrderDiameterAv")

#PRepare structure
Tab <- df[,mcolumns] %>% as.matrix
rownames(Tab) <- df$root

Map <- df[,c("root","Strain","Batch")]
rownames(Map) <- df$root

#Appendn the GTDB information
df_gtdb <- read.table(file = "./rawdata/gtdbtk.bac120.summary.tsv",
                      header = T,sep = "\t",quote = "",comment.char = "")
df_gtdb <- df_gtdb[,c(1,2)]
df_tbind <- df_gtdb$classification %>% 
  strsplit(split = "\\;") %>% unlist %>%
  gsub(pattern = "[dpcofgs]__",replacement = "") %>%
  matrix(data = .,ncol = 7,byrow = T) %>% as.data.frame
colnames(df_tbind) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
df_gtdb <- cbind(df_gtdb,df_tbind)
colnames(df_gtdb)[1] <- "taxon_oid"

Dat <- create_dataset(Tab = Tab %>% t,Map = Map)

Tab_raw <- Dat$Tab %>% t

Tab_q <- Tab_raw
for(i in 1:ncol(Tab_raw)){
  
  x <- Tab_raw[,i]
  max_value <- quantile(x,0.95)
  min_value <- quantile(x,0.05)
  x[x>max_value] <- max_value
  x[x<min_value] <- min_value
  Tab_q[,i] <- x
}


#Perform log10 normalization
Tab_log <- apply(X = Tab_q,MARGIN = 2,FUN = function(x)log10(x+1))
Tab_log_centered <- scale(x = Tab_log,center = T,scale = T)

#Create range 0-1 based dataset
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

Tab_log_range <- apply(X = Tab_log,MARGIN = 2,FUN = normalize) 

#Create objects
Dat_raw <- create_dataset(Tab = Tab_q %>% t,Map = Map)
Dat_log10 <- create_dataset(Tab = Tab_log %>% t,Map = Map)
Dat_log10_c <- create_dataset(Tab = Tab_log_centered %>% t,Map = Map)
Dat_log10_range <- create_dataset(Tab = Tab_log_range %>% t,Map = Map)

Dat_raw_z <- create_dataset(Tab = Dat_raw$Tab %>% t %>% scale(center = T,scale = T) %>% t,Map = Map)
Tab_raw_range <- apply(X = Dat_raw$Tab %>%  t ,MARGIN = 2,FUN = normalize)
Dat_raw_range <- create_dataset(Tab = Tab_raw_range %>% t,Map = Map)


###Create the genomic dataset
tree <- read.tree(file = "./rawdata/tree_1130_isolates_plant_microbiome_fastree.rooted.newick")
df_isolates <- read.table(file = "./rawdata/373_isolates_genome.tsv",header = T,sep = "\t")
todrop <- which(!(tree$tip.label %in% df_isolates$IMG_taxon_oid)) %>%
  tree$tip.label[.]
tree <- ape::drop.tip(phy = tree,tip = todrop)
tree$tip.label <- match(tree$tip.label,df_isolates$IMG_taxon_oid) %>%
  df_isolates$FreezerId[.] %>% as.character

df_gtdb$Strain <- match(df_gtdb$taxon_oid,df_isolates$IMG_taxon_oid) %>%
  df_isolates$FreezerId[.]
Map <- df_gtdb
rownames(Map) <- Map$Strain
Map <- Map[,c(10,1,3:9)]

### Rename variables to make them consistent across all dataset
original_name <- c("Norm_length","Norm_diameter","Norm_volume","Norm_surface","Norm_LR1stOrderNumRootsbyCm",
  "Norm_LR1stInsertionAngleAv","Norm_LR1stOrderAsymmetryPattern",
  "Norm_LR1stOrderDiameterAv","Norm_LR1stOrderDiffLengthTopBottom","Norm_LR1stOrderInterBranch","Norm_LR1stOrderLengthAv",            
  "Norm_LR1stOrderPropRangeLasttoFirst","Norm_LR1stOrderPropRangeTiptoLast","Norm_LR1stOrderRAwith2ndOrder",
  "Norm_LR1stOrderSurfaceAv","Norm_LR1stOrderVolumeAv","Norm_LR2ndOrderDiameterAv",
  "Norm_LR2ndOrderLengthAv","Norm_LR2ndOrderNumberNorm","Norm_LR2ndOrderSurfaceAv","Norm_LR2ndOrderVolumeAv","Norm_child_density")

new_name <- c("PR_length","PR_diameter","PR_volume","PR_surface","LR_density","LR_angle","LR_symmetry",
"LR_diameter","LR_length_distribution","LR_interbranch_distance","LR_length",
"LR_coverage","Distance_to_the_1st_LR","Number_of_branched_LR","LR_surface","LR_volume","2ndLR_diameter",
"2ndLR_length","2ndLR_number","2ndLR_surface","2ndLR_volume","Norm_child_density")

df_rename <- data.frame(OriginalName = original_name,
            NewName = new_name)


#Rename columns
rownames(Dat_raw$Tab) <- match(rownames(Dat_raw$Tab),df_rename$OriginalName) %>%
  df_rename$NewName[.]

rownames(Dat_log10$Tab) <- match(rownames(Dat_log10$Tab),df_rename$OriginalName) %>%
  df_rename$NewName[.]

rownames(Dat_log10_c$Tab) <- match(rownames(Dat_log10_c$Tab),df_rename$OriginalName) %>%
  df_rename$NewName[.]

rownames(Dat_log10_range$Tab) <- match(rownames(Dat_log10_range$Tab),df_rename$OriginalName) %>%
  df_rename$NewName[.]


rownames(Dat_raw_z$Tab) <- match(rownames(Dat_raw_z$Tab),df_rename$OriginalName) %>%
  df_rename$NewName[.]

rownames(Dat_raw_range$Tab) <- match(rownames(Dat_raw_range$Tab),df_rename$OriginalName) %>%
  df_rename$NewName[.]


#Rename original df
mindices <- which(colnames(df) %in% df_rename$OriginalName)
torename <- match(colnames(df)[mindices],df_rename$OriginalName) %>%
  df_rename$NewName[.]
colnames(df)[mindices] <- torename


mlist <- list(df_raw = df,
              
              Dat_raw = Dat_raw,
              Dat_log10 = Dat_log10,
              Dat_log10_zscore = Dat_log10_c,
              Dat_log10_range = Dat_log10_range,
              Dat_raw_zscore = Dat_raw_z,
              Dat_raw_range = Dat_raw_range,
              
              tree = tree,
              Map_strains = Map
)



saveRDS(object = mlist,file = "./cleandata/dat_rootbranching_screening1.RDS")

rm(list=ls())
gc()

