library(ohchibi)
library(phytools)
library(RColorBrewer)
library(xlsx)



set.seed(130816)

#Load the res from the model and the metadata and tree
Dat_all <- readRDS(file = "./cleandata/dat_rootbranching_screening1.RDS")


merged_res <- Dat_all$Map_strains
tree <- Dat_all$tree
  
merged_res$Genome_Name <- paste0(merged_res$Genus, " sp.",merged_res$Strain)

#Prepare files for itol
dir.create(path = "./cleandata/itol")
write.tree(phy = tree,file = "./cleandata/itol/screening_tree_rootbranching_revisited.newick")

#LABELS
outfile <- "./cleandata/itol/itol_screening_tree_rootbranching_taxonoid2names.txt"
mline <- c("LABELS","Separator TAB","DATA")
write(x = mline,file = outfile,append = F,sep = "")
df <- data.frame(taxon_oid =merged_res$Strain,Names = merged_res$Genome_Name)
write.table(x = df,file = outfile,append = T,quote = F,sep = "\t",row.names = F,col.names = F)

#Taxonomy
mpalettes <- paste0(paletteer_d("ggthemes::Tableau_20",n = 19),collapse= "\t")
mnames <- merged_res$Order %>% table %>% names %>% paste0(.,collapse = "\t")
mshapes <- paste0(rep("2",19),collapse = "\t")
paleta <- paletteer_d("ggthemes::Tableau_20",n = 19) %>% as.character
names(paleta) <- merged_res$Order %>% table %>% names

df <- data.frame(taxon_oid = merged_res$Strain, symbol = rep(2,nrow(merged_res)),
                 size = rep(10,nrow(merged_res)),col=match(merged_res$Order,names(paleta)) %>% paleta[.] %>% as.character,
                 fill = rep(1,nrow(merged_res)),position = rep(1,nrow(merged_res)))
outfile <- "./cleandata/itol/itol_screening_tree_rootbranching_taxonoid2symbol.txt"
mline <- c("DATASET_SYMBOL","SEPARATOR TAB","DATASET_LABEL\tOrder","COLOR\tblack","MAXIMUM_SIZE\t10",
           "LEGEND_TITLE\tOrder names",
           paste0("LEGEND_SHAPES\t",mshapes,collapse = ""),
           paste0("LEGEND_COLORS\t",mpalettes,collapse = ""),
           paste0("LEGEND_LABELS\t",mnames,collapse = ""),
           "DATA")
write(x = mline,file = outfile,append = F,sep = "")

write.table(x = df,file = outfile,append = T,quote = F,sep = "\t",row.names = F,col.names = F)

#Create distribution of operons
df_iaa_acc_genome <- readRDS(file = "./cleandata/dat_rootbranching_cog_ko.RDS") %$% df_iaa_acc_genome

df_iaa_acc_genome <- data.frame(Strain = "L404",taxon_oid = "L404",
           IAA_BioIII_K00466 = 0,
           IAA_BioIV_K01721 = 0,
           IAA_BioV_K01501 = 0,
           IAA_BioVI_K04103 = 0,
           IAADegradationiac = 0,
           IAADegradatioiad = 0 ,
           ACC_K01505 = 0,
           Sum_IAA_Bio = 0,
           Sum_IAA_Deg = 0
           )  %>%
  rbind(df_iaa_acc_genome,.)



outfile <- "./cleandata/itol/itol_screening_tree_rootbranching_taxonoid2operon.txt"
mline <- c("DATASET_BINARY","SEPARATOR TAB","DATASET_LABEL\tOperons","COLOR\tpurple",
           "FIELD_SHAPES\t4\t4\t4\t4\t5\t5\t3",
           "FIELD_LABELS\tIAA_BioIII\tIAA_BioIV\tIAA_BioV\tIAA_BioVI\tIAA_Degiac\tIAA_Degiad\tSelected",
           "FIELD_COLORS\t#414141\t#414141\t#414141\t#414141\t#414141\t#414141\t#414141","HEIGHT_FACTOR\t2",
           "DATA")
write(x = mline,file = outfile,append = F,sep = "")

#df <- df_iaa_acc_genome[,c(1,3:9)]
df <- df_iaa_acc_genome[,c(1,3:8)]


#Append selected strains

selected_strains <- read.xlsx(file = "./rawdata/rawdata_selaginella_screening.xlsx",sheetIndex = 1) %>%
  dplyr::rename(.data = .,Strain = Strains,
                NumberPlant = explant.nb,
                RhizophoreNumber = nb.of.Rhizophore) %>%
  dplyr::select(.data = .,c(Strain,NumberPlant,RhizophoreNumber,Bifurcation)) %>%
  dplyr::mutate(.data = .,Strain = Strain %>% gsub(pattern = "_",replacement = "") %>%
                  gsub(pattern = "MF",replacement = "RMF") %>%
                  gsub(pattern = "^CL",replacement = "RCL")) %$% Strain %>%
  grep(pattern = "NB",invert = T,value = T) %>% unique

df$Selected <- 0
df$Selected[which(df$Strain %in% selected_strains)] <- 1



write.table(x = df,file = outfile,append = T,quote = F,sep = "\t",row.names = F,col.names = F)


#Enrich
rm(list=ls())

