library(ohchibi)
library(multcomp)
library(ggpubr)
library(emmeans)
library(ggrepel)

set.seed(130816)

Dat_suberin <- readRDS(file = "./cleandata/dat_screeningbarriers.RDS")

df_suberin <- Dat_suberin$Suberin

df_suberin_av <- df_suberin %>%
  subset(Localization == "Continous_expression") %>%
  aggregate(Normvalue_ra~Strains,.,mean) %>%
  dplyr::mutate(.data = .,Normvalue_ra = 1-Normvalue_ra) %>%
  dplyr::rename(.data = .,SumZones = Normvalue_ra,Strain = Strains)


#Get the lateral root density information 
Dat_all <- readRDS(file = "./cleandata/dat_rootbranching_screening1.RDS")

###### Z score values ###########
Dat <- Dat_all$Dat_raw

df_lr <- Dat$Tab %>% melt %>%
  subset(Var1 == "LR_density") %>% 
  dplyr::rename(.data = .,Feature = Var1,
                root = Var2) %>%
  merge(Dat$Map,by = "root") %>%
  aggregate(value~Strain,.,mean) %>%
  dplyr::rename(.data = .,LRDensity = value)

merged <- merge(df_lr,df_suberin_av, by = "Strain")
cor.test(merged$LRDensity,merged$SumZones,method = "pearson")

p <- merge(df_lr,df_suberin_av, by = "Strain") %>%
  ggplot(data = .,aes(SumZones,LRDensity)) +
  geom_smooth(method = "lm",se = FALSE,color = "red",size = 1)+
  ggpubr::stat_cor(method = "pearson",label.sep = "\n")+
  geom_point() +
  theme_ohchibi(size_panel_border = 0.3)  +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank()
  ) + xlab(label = "Suberin (Sum of No suberin + Discrete zones) ") +
  ylab("Lateral root density")


oh.save.pdf(p = p,outname = "fig2_cor_suberin_lr.all.pdf",outdir = "./figures/",width = 8,height = 8)
