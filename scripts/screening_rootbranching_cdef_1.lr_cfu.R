library(ohchibi)
library(xlsx)
library(emmeans)
library(multcomp)

set.seed(130816)


mpath <- "./rawdata/All_fig2_data/8_CDEF line LR density with and without bacteria/"

#Read other data in xlsx format
mfiles <- list.files(path = mpath,pattern = "xlsx",full.names = FALSE,recursive = TRUE) %>%
  grep(pattern = "Mathieu",invert = TRUE,value = TRUE)%>%
  grep(pattern = "cfu",invert = TRUE,value = TRUE) 


df <- NULL
for(mfile in mfiles){
  
  
  infile <- paste0(mpath,mfile)
  df <- (read.xlsx(file = infile,sheetIndex = 1))[,1:4] %>%
    na.omit %>%
    dplyr::mutate(.data = .,
                  File = mfile) %>%
    rbind(df,.)
  
}

df$Strain <- "NB"
df$Strain[df$File %>% grep(pattern = "MF27")] <- "RMF27"

df$Treatment <- "Col-0"
df$Treatment[df$File %>% grep(pattern = "CDEF")] <- "CDEF"

df$Genotype <- "Col-0"


df_xlsx <- df
rm(df)

df_xlsx$PR_number <- df_xlsx$PR_number %>% as.numeric
df_xlsx$PR_length <- df_xlsx$PR_length %>% as.numeric
df_xlsx$LR_number <- df_xlsx$LR_number %>% as.numeric
df_xlsx$LR_Density <- df_xlsx$LR_Density %>% as.numeric

df_xlsx$LR_Density <- df_xlsx$LR_number/df_xlsx$PR_length

df_root <- df_xlsx
rm(df_xlsx)

df_root$Group <- paste0(df_root$Strain," ",df_root$Treatment) %>%
  factor(levels = c("NB Col-0","NB CDEF",
                    "RMF27 Col-0","RMF27 CDEF"))

#Perform test
lm(formula = LR_Density~Group,df_root) %>%
  emmeans(object = .,specs = "Group") %>%
  multcomp::cld()

#PRimary root plot
p1 <- ggplot(data = df_root,aes(Group,LR_Density)) +
  geom_sina(alpha = 0.3) +
  stat_summary(fun.data = mean_cl_normal,geom = "pointrange",color = "red",shape = 16,size = 1)+
  theme_ohchibi(size_panel_border = 0.3) +
  scale_y_continuous(breaks = seq(0,6,1),limits = c(0,6),expand = c(0,0)) +
  theme(
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1) ,
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  ylab(label = "Lateral root density") +
  xlab(label = "Treatment")

### Read the CFU data
df_cfu <- read.xlsx(file = "./rawdata/All_fig2_data/8_CDEF line LR density with and without bacteria/cfu_only_cdef.xlsx",sheetIndex = 1)
df_cfu$Lvalue <- log10(df_cfu$cfu.mg)
df_cfu$Lvalue[is.infinite(df_cfu$Lvalue)] <- 0

df_cfu$Group <- paste0(df_cfu$Bacteria," ",df_cfu$Genotype)

#Remove NB fropm CFU
df_cfu <- df_cfu %>% subset(Bacteria != "NB") %>% droplevels
df_cfu$Group <- df_cfu$Group %>%
  factor(levels = c("RMF27 Col0","RMF27 CDEF"))

p2 <- ggplot(data = df_cfu,aes(Group,Lvalue)) +
  geom_jitter(alpha = 0.3) +
  stat_summary(fun = mean,geom = "point",color = "red",shape = 16,size = 3)+
  theme_ohchibi(size_panel_border = 0.3) +
  theme(
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1) ,
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  ylab(label = "log10 (CFU)") +
  xlab(label = "Treatment") +
  scale_y_continuous(breaks = seq(0,6,1),limits = c(0,6),expand = c(0,0)) 
  


composition <- egg::ggarrange(p1,p2,nrow =1,widths = c(1,0.6))

oh.save.pdf(p = composition,outname = "fig2_cdef_jitter.pdf",outdir = "./figures/",width = 14,height = 10)

rm(list=ls())
