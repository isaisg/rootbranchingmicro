library(ohchibi)
library(xlsx)

set.seed(130816)

mfiles <- list.files(path = "./rawdata/All_fig2_data/10_Venus",pattern = "xlsx",full.names = TRUE)

#Loop over the notebook and read
df_venus <- NULL

for(mfile in mfiles){
  
  wb <- loadWorkbook(mfile)
  msheets <- xlsx::getSheets(wb)
  
  for( msheet in names(msheets)){
    
    cat("Working on ",mfile," ",msheet,"\n")  
    if(mfile == "./rawdata/All_fig2_data/10_Venus/Data analyze day 2.xlsx" & msheet == "Project_NB_2.4_z00 "){
      
      df_venus <- (read.xlsx(file = mfile,sheetName = msheet,startRow = 4))[,1:7] %>%
        dplyr::mutate(.data = .,File = mfile,Sheet = msheet) %>% 
        rbind(df_venus,.)
      
      
    }else{
      
      df_venus <- read.xlsx(file = mfile,sheetName = msheet,startRow = 3) %>%
        dplyr::mutate(.data = .,File = mfile,Sheet = msheet) %>%
        rbind(df_venus,.)
    }
    
    gc()
    
    
  }
  
}

df_venus_raw <- df_venus

#Process the dataset and rename
df_venus <- df_venus_raw %>%
  dplyr::filter(.data = .,!Sheet %>% grepl(pattern = "LR2")) %>%
  droplevels

df_venus$Strain <- df_venus$Sheet %>%
  gsub(pattern = "Project_",replacement = "") %>%
  gsub(pattern = ".z.*",replacement = "") %>% 
  gsub(pattern = "^L_",replacement = "L") %>%
  gsub(pattern = "^CL_",replacement = "RCL") %>%
  gsub(pattern = "^MF_",replacement = "RMF") %>%
  gsub(pattern = "^MF",replacement = "RMF") %>%
  gsub(pattern = "_.*",replacement = "") %>%
  factor %>% relevel(ref = "NB")

df_venus$Batch <-   df_venus$File %>% 
  gsub(pattern = "\\.xlsx",replacement = "") %>%
  gsub(pattern = ".* ",replacement = "") %>%
  paste0("Batch",.) %>% factor()


df_venus$Plate <-   df_venus$Sheet %>%
  gsub(pattern = "Project_",replacement = "") %>%
  gsub(pattern = ".z.*",replacement = "") %>% 
  gsub(pattern = "^L_",replacement = "L") %>%
  gsub(pattern = "^CL_",replacement = "RCL") %>%
  gsub(pattern = "^MF_",replacement = "RMF") %>%
  gsub(pattern = "^MF",replacement = "RMF") %>%
  gsub(pattern = ".*_",replacement = "") %>%
  gsub(pattern = "\\..*",replacement = "") %>% 
  paste0("Plate",.)

df_venus$Plant <-   df_venus$Sheet %>%
  gsub(pattern = "Project_",replacement = "") %>%
  gsub(pattern = ".z.*",replacement = "") %>% 
  gsub(pattern = "^L_",replacement = "L") %>%
  gsub(pattern = "^CL_",replacement = "RCL") %>%
  gsub(pattern = "^MF_",replacement = "RMF") %>%
  gsub(pattern = "^MF",replacement = "RMF") %>%
  gsub(pattern = ".*_",replacement = "") %>%
  gsub(pattern = ".*\\.",replacement = "") %>% 
  paste0("Plant",.)


#Determine the measured area
df_venus$UId <- paste0(df_venus$Batch,"_",df_venus$Plate,"_",df_venus$Strain,"_",df_venus$Plant)

muids <- df_venus$UId %>% unique

df_area <- NULL
for(muid in muids){
  
  df_inner <- df_venus %>%
    subset(UId == muid)  %>% droplevels
  df_area <- data.frame(UId = muid,NumMeas = nrow(df_inner),SumArea = df_inner$area %>% sum) %>%
    rbind(df_area,.)
  
}



df_ag <- aggregate(Mean~Batch+Strain+Plant,df_venus,median) 


#Compare the distributions
m1 <- lm(formula =Mean~Strain,df_ag )
glht(model = m1,linfct = mcp(Strain="Dunnet")) %>% 
  summary

#Color based on significance and order based on agrgregate
df_ag$Color <- "NS"
df_ag$Color[df_ag$Strain %>% grep(pattern = "NB")] <- "NB"


#Perform same type of aesthetic
paleta <- c("#FFCA64","#CFB0D4","#A7A9AC","#E2F3F7")
names(paleta) <- c("NB","Up","NS","Down")

p1 <- ggplot(data = df_ag,aes(Strain,Mean)) +
  geom_jitter(alpha = 0.3,aes(color = Color)) +
  stat_summary(fun.data = mean_cl_normal,geom = "linerange",shape = 21,stroke = 0.3,size = 1,aes(color = Color))+
  stat_summary(fun = mean,geom = "point",shape = 21,stroke = 0.3,size = 3,aes(fill = Color))+
  theme_ohchibi(size_panel_border = 0.3) +
  ylab(label = "Quantified fluorescence") + 
  ggtitle(label = "No bacteria treatment") +
  theme(
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank()
  )   +
  ggtitle(label = "Venus reporter line")  +
  scale_color_manual(values = paleta,name = "Significance\nvsNB") +
  scale_fill_manual(values = paleta,name = "Significance\nvsNB") 

#Save figure
oh.save.pdf(p = p1,outname = "fig2_venus.pdf",
            outdir = "./figures/",width = 9,height =8)

rm(list=ls())
