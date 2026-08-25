library(tidyverse)
library(ggpubr)
library(ggsignif)
require("ggrepel")
library(ggbreak) 

# your folder path
setwd("")

dir.create("02_ReportedData_hRIPK3")
path="02_ReportedData_hRIPK3"

### hRIPK3
# Nucleation Score
load("nscore_df_hRIPK3.RData")

singles_stops$sig_1<-F
singles_stops[singles_stops$p.adjust<0.01, "sig_1"]<-T

singles_stops$category_1<-"WT-like"
singles_stops[singles_stops$p.adjust<0.01 & singles_stops$nscore_c >0, "category_1"]<-"NS_inc"
singles_stops[singles_stops$p.adjust<0.01 & singles_stops$nscore_c <0, "category_1"]<-"NS_dec"

singles_stops<-singles_stops[singles_stops$low_sigma == T | singles_stops$sig_1 == T, ]

singles_stops$Residue<-paste0(singles_stops$WT_AA, singles_stops$Pos)

### ddG
#PDB 7DAC
load("ddG_7DAC.RData")

# Residues classification based on side chain orientation
ddG_foldx_7DAC$orientation_7DAC<-"exposed"
ddG_foldx_7DAC[ddG_foldx_7DAC$Pos %in% c(448, 450, 452, 455, 456, 458, 459, 460, 462, 466, 468, 470), "orientation_7DAC"]<-"buried"
ddG_foldx_7DAC[ddG_foldx_7DAC$Pos %in% c(457, 461), "orientation_7DAC"]<-"glycine"

ddG_foldx_7DAC$sig_1_7DAC<-F
ddG_foldx_7DAC[ddG_foldx_7DAC$p.adjust_7DAC<0.01, "sig_1_7DAC"]<-T

ddG_foldx_7DAC$category_1_7DAC<-"WT-like"
ddG_foldx_7DAC[ddG_foldx_7DAC$p.adjust_7DAC<0.01 & ddG_foldx_7DAC$ddG_foldx_7DAC >2, "category_1_7DAC"]<-"ddG_inc"
ddG_foldx_7DAC[ddG_foldx_7DAC$p.adjust_7DAC<0.01 & ddG_foldx_7DAC$ddG_foldx_7DAC <0, "category_1_7DAC"]<-"ddG_dec"

ddG_foldx_7DAC<-ddG_foldx_7DAC[ddG_foldx_7DAC$low_ddG_foldx_SD_7DAC == T | ddG_foldx_7DAC$sig_1_7DAC == T, ]


#PDB 7DA4
load("ddG_7DA4.RData")

# Residues classification based on side chain orientation
ddG_foldx_7DA4$orientation_7DA4<-"exposed"
ddG_foldx_7DA4[ddG_foldx_7DA4$Pos %in% c(450, 452, 454, 456, 458, 459, 460, 466, 468), "orientation_7DA4"]<-"buried"
ddG_foldx_7DA4[ddG_foldx_7DA4$Pos %in% c(457, 461), "orientation_7DA4"]<-"glycine"

ddG_foldx_7DA4$sig_1_7DA4<-F
ddG_foldx_7DA4[ddG_foldx_7DA4$p.adjust_7DA4<0.01, "sig_1_7DA4"]<-T

ddG_foldx_7DA4$category_1_7DA4<-"WT-like"
ddG_foldx_7DA4[ddG_foldx_7DA4$p.adjust_7DA4<0.01 & ddG_foldx_7DA4$ddG_foldx_7DA4 >2, "category_1_7DA4"]<-"ddG_inc"
ddG_foldx_7DA4[ddG_foldx_7DA4$p.adjust_7DA4<0.01 & ddG_foldx_7DA4$ddG_foldx_7DA4 <0, "category_1_7DA4"]<-"ddG_dec"

ddG_foldx_7DA4<-ddG_foldx_7DA4[ddG_foldx_7DA4$low_ddG_foldx_SD_7DA4 == T | ddG_foldx_7DA4$sig_1_7DA4 == T, ]


#PDB 8Z94
load("ddG_8Z94.RData")

# Residues classification based on side chain orientation
ddG_foldx_8Z94$orientation_8Z94<-"exposed" 
ddG_foldx_8Z94[ddG_foldx_8Z94$Pos %in% c(450, 452, 454, 456, 458, 459, 460, 466, 468), "orientation_8Z94"]<-"buried"
ddG_foldx_8Z94[ddG_foldx_8Z94$Pos %in% c(457, 461), "orientation_8Z94"]<-"glycine"

ddG_foldx_8Z94$sig_1_8Z94<-F
ddG_foldx_8Z94[ddG_foldx_8Z94$p.adjust_8Z94<0.01, "sig_1_8Z94"]<-T

ddG_foldx_8Z94$category_1_8Z94<-"WT-like"
ddG_foldx_8Z94[ddG_foldx_8Z94$p.adjust_8Z94<0.01 & ddG_foldx_8Z94$ddG_foldx_8Z94 >2, "category_1_8Z94"]<-"ddG_inc"
ddG_foldx_8Z94[ddG_foldx_8Z94$p.adjust_8Z94<0.01 & ddG_foldx_8Z94$ddG_foldx_8Z94 <0, "category_1_8Z94"]<-"ddG_dec"

ddG_foldx_8Z94<-ddG_foldx_8Z94[ddG_foldx_8Z94$low_ddG_foldx_SD_8Z94 == T | ddG_foldx_8Z94$sig_1_8Z94 == T, ]

###
NS_ddG_RIPK3<-full_join(ddG_foldx_7DA4[c("ID", "ddG_foldx_7DA4", "ddG_foldx_SD_7DA4", "orientation_7DA4", "category_1_7DA4")], 
                        ddG_foldx_7DAC[c("ID", "ddG_foldx_7DAC", "ddG_foldx_SD_7DAC", "orientation_7DAC", "category_1_7DAC")], by="ID")

NS_ddG_RIPK3<-full_join(NS_ddG_RIPK3, ddG_foldx_8Z94[c("ID", "ddG_foldx_8Z94", "ddG_foldx_SD_8Z94", "orientation_8Z94", "category_1_8Z94")], by="ID")

NS_ddG_RIPK3<-full_join(NS_ddG_RIPK3, singles_stops[c("ID", "Pos", "Residue", "nscore_c", "sigma", "category_10", "category_1")], by="ID")

################################################################################
# 1: Correlation vs grouped mutants from the literature

#CD = cell death
CD_dec<-c("V-450-E", "V-460-E", "L-466-E", "V-460-P", "V-458-P", "G-457-D")
CD_WT<-c("N-464-D", "M-468-D")

NS_ddG_RIPK3$CD<-""
NS_ddG_RIPK3[NS_ddG_RIPK3$ID %in% CD_dec, "CD"]<-"Decreased"
NS_ddG_RIPK3[NS_ddG_RIPK3$ID %in% CD_WT, "CD"]<-"WT-like"

p_rug_1<-ggplot(NS_ddG_RIPK3[NS_ddG_RIPK3$CD != "",], aes(x=nscore_c, color=CD))+
  geom_rug(size=2)+
  scale_color_manual(values=c("red3", "grey40"))+
  labs(x="Nucleation Score", color="Cell Death")+
  theme_classic()+
  theme(legend.position = "none")
p_rug_1

ggsave(p_rug_1, file="ReportedMutants_nscore.jpg", width=4, height=4, path=path)
ggsave(p_rug_1, file="ReportedMutants_nscore.pdf", width=4, height=4, path=path)

#
p_rug_2<-ggplot(NS_ddG_RIPK3[NS_ddG_RIPK3$CD != "",], aes(x=ddG_foldx_7DA4, color=CD))+
  geom_rug(size=2)+
  scale_color_manual(values=c("red3", "grey40"))+
  labs(x="ddG 7DA4", color="Cell Death")+
  theme_classic()+
  theme(legend.position = "none")
p_rug_2

ggsave(p_rug_2, file="ReportedMutants_ddG_7DA4.jpg", width=4, height=4, path=path)
ggsave(p_rug_2, file="ReportedMutants_ddG_7DA4.pdf", width=4, height=4, path=path)

#
p_rug_3<-ggplot(NS_ddG_RIPK3[NS_ddG_RIPK3$CD != "",], aes(x=ddG_foldx_7DAC, color=CD))+
  geom_rug(size=2)+
  scale_color_manual(values=c("red3", "grey40"))+
  labs(x="ddG 7DAC", color="Cell Death")+
  theme_classic()+
  xlim(0, 20)+
  theme(legend.position = "none")
p_rug_3

ggsave(p_rug_3, file="ReportedMutants_ddG_7DAC.jpg", width=4, height=4, path=path)
ggsave(p_rug_3, file="ReportedMutants_ddG_7DAC.pdf", width=4, height=4, path=path)

#
p_rug_4<-ggplot(NS_ddG_RIPK3[NS_ddG_RIPK3$CD != "",], aes(x=ddG_foldx_8Z94, color=CD))+
  geom_rug(size=2)+
  scale_color_manual(values=c("red3", "grey40"))+
  labs(x="ddG 8Z94", color="Cell Death")+
  theme_classic()+
  theme(legend.position = "none")
p_rug_4

ggsave(p_rug_4, file="ReportedMutants_ddG_8Z94.jpg", width=4, height=4, path=path)
ggsave(p_rug_4, file="ReportedMutants_ddG_8Z94.pdf", width=4, height=4, path=path)

################################################################################
# 2: NS Correlation vs mutants from the literature
# Cell death from Li et al Cell 2012
RIPK3_mutants_HeLa_1<-data.frame("ID"=c("WT", "V-458-P", "V-460-P", "G-457-D", "N-464-D"),
                                 "cell_death"=c(24.03, 11.76, 6.14, 12.44, 18.24),
                                 "SD"=c(1.54, 2.39, 1.7, 4.26, 3.58))

singles_stops$category_1<-"WT-like"
singles_stops[singles_stops$p.adjust<0.01 & singles_stops$nscore_c >0, "category_1"]<-"NS_inc"
singles_stops[singles_stops$p.adjust<0.01 & singles_stops$nscore_c <0, "category_1"]<-"NS_dec"

singles_ns_death<-inner_join(singles_stops[c("ID", "nscore_c", "sigma", "category_1")], RIPK3_mutants_HeLa_1, by="ID")

singles_ns_death<-rbind(data.frame("ID"=c("WT"), "nscore_c"=c(0.0006341024), "sigma"=c(0.07340256), "category_1"=c("WT-like"),
                                   "cell_death"=c(24.03), "SD"=c(1.54)), singles_ns_death)

scatter_1<-ggplot(singles_ns_death, aes(x=nscore_c, y=cell_death, label=ID, color=category_1))+
  geom_errorbar(aes(ymin=cell_death-SD, ymax=cell_death+SD), 
                color="grey80", width=.05)+
  geom_errorbar(aes(xmin=nscore_c-sigma, xmax=nscore_c+sigma),
                color="grey80", width=1)+
  geom_point(size=4, color="grey40")+
  #scale_color_manual(values=c("#CD3333", "grey", "#00008B"))+
  geom_text_repel(size=4, color="grey40")+
  stat_cor(label.y = 30, label.x = -4, color="black")+ 
  labs(x="Nucleation Score", y="% Cell Death (Li et al 2012)")+
  xlim(-4, 1)+
  ylim(0, 30)+
  theme_classic()+
  theme(legend.position = "none")

scatter_1
ggsave(scatter_1, file="nscore_Lietal2012.jpg", width=4, height=4, path=path)
ggsave(scatter_1, file="nscore_Lietal2012.pdf", width=4, height=4, path=path)


# Cell-survival from Hu et al Cell Death Differentiation 2020
RIPK3_mutants_HeLa_2<-data.frame("ID"=c("WT", "N-464-D", "M-468-D"),
                                 "cell_survival"=c(12.19, 12.19, 7.53),
                                 "SD"=c(1.07, 3.22, 1.79))

singles_ns_death<-inner_join(singles_stops[c("ID", "nscore_c", "sigma", "category_1")], RIPK3_mutants_HeLa_2, by="ID")

singles_ns_death<-rbind(data.frame("ID"=c("WT"), "nscore_c"=c(0.0006341024), "sigma"=c(0.07340256), "category_1"=c("WT-like"),
                                   "cell_survival"=c(12.19), "SD"=c(1.07)), singles_ns_death)

scatter_2<-ggplot(singles_ns_death, aes(x=nscore_c, y=cell_survival, label=ID, color=category_1))+
  geom_errorbar(aes(ymin=cell_survival-SD, ymax=cell_survival+SD), 
                color="grey80", width=.05)+
  geom_errorbar(aes(xmin=nscore_c-sigma, xmax=nscore_c+sigma),
                color="grey80", width=1)+
  geom_point(size=4, color="grey40")+
  #scale_color_manual(values=c("#CD3333", "grey", "#00008B"))+
  geom_text_repel(size=4, color="grey40")+
  stat_cor(label.y = 100, label.x = -4, color="black")+ 
  xlim(-4, 1)+
  ylim(0, 25)+
  labs(x="Nucleation Score", y="% Cell Survival (Hu et al 2020)")+
  theme_classic()+
  theme(legend.position = "none")

scatter_2
ggsave(scatter_2, file="nscore_Huetal2021.jpg", width=4, height=4, path=path)
ggsave(scatter_2, file="nscore_Huetal2021.pdf", width=4, height=4, path=path)

# Wu et al. 2021
# They also demonstrated by EM and NMR that mutants L-456-D and Q-449-D change the fibril structure
mRIP3_tht_wt<-c("L-456-D")
mRIP3_tht_dec<-c("Q-449-D", "F-442-D")

# Ma et al. 2025 
# They also showed that mutant L-466-E forms a different fibril structure by EM and IF.
RIPK3_mutants_HeLa<-data.frame("ID"=c("WT", "V-450-E", "V-460-E", "L-466-E"),
                               "cell_viability"=c(51.74, 80.23, 86.63, 96.22),
                               "SD"=c(3.49, 11.05, 5.23, 6.11))

singles_ns_death<-inner_join(singles_stops[c("ID", "nscore_c", "sigma", "category_1")], RIPK3_mutants_HeLa, by="ID")

singles_ns_death<-rbind(data.frame("ID"=c("WT"), "nscore_c"=c(0.0006341024), "sigma"=c(0.07340256), "category_1"=c("WT-like"),
                                   "cell_viability"=c(51.74), "SD"=c(3.49)), singles_ns_death)

scatter_3<-ggplot(singles_ns_death, aes(x=nscore_c, y=cell_viability, label=ID, color=category_1))+
  geom_errorbar(aes(ymin=cell_viability-SD, ymax=cell_viability+SD), 
                color="grey80", width=.05)+
  geom_errorbar(aes(xmin=nscore_c-sigma, xmax=nscore_c+sigma),
                color="grey80", width=1)+
  geom_point(size=4, color="grey40")+
  #scale_color_manual(values=c("#CD3333", "grey", "#00008B"))+
  geom_text_repel(size=4, color="grey40")+
  stat_cor(label.y = 30, label.x = -4, color="black")+ 
  labs(x="Nucleation Score", y="% Cell Viability (Ma et al 2025)")+
  #xlim(-4, 1)+
  ylim(0, 110)+
  theme_classic()+
  theme(legend.position = "none")

scatter_3
ggsave(scatter_3, file="nscore_Maetal2025.jpg", width=4, height=4, path=path)
ggsave(scatter_3, file="nscore_Maetal2025.pdf", width=4, height=4, path=path)

################################################################################
# 3: ddG Correlation vs mutants from the literature

RIPK3_mutants_HeLa_1<-data.frame("ID"=c("WT", "V-458-P", "V-460-P", "G-457-D", "N-464-D"),
                                 "cell_death"=c(24.03, 11.76, 6.14, 12.44, 18.24),
                                 "SD"=c(1.54, 2.39, 1.7, 4.26, 3.58))

singles_ddG_death<-inner_join(NS_ddG_RIPK3[c("ID", "ddG_foldx_7DA4", "ddG_foldx_SD_7DA4", "ddG_foldx_7DAC", "category_1_7DA4","ddG_foldx_SD_7DAC", "category_1_7DAC", "ddG_foldx_8Z94", "ddG_foldx_SD_8Z94", "category_1_8Z94")], 
                              RIPK3_mutants_HeLa_1, by="ID")

singles_ddG_death<-rbind(data.frame("ID"=c("WT"), "ddG_foldx_7DA4"=c(0), "ddG_foldx_SD_7DA4"=c(0), "category_1_7DA4"=c("WT-like"), "ddG_foldx_7DAC"=c(0), "ddG_foldx_SD_7DAC"=c(0), "category_1_7DAC"=c("WT-like"), "ddG_foldx_8Z94"=c(0), "ddG_foldx_SD_8Z94"=c(0), "category_1_8Z94"=c("WT-like"), 
                                    "cell_death"=c(24.03), "SD"=c(1.54)), singles_ddG_death)

scatter_4<-ggplot(singles_ddG_death, aes(x=ddG_foldx_7DA4, y=cell_death, label=ID, color=category_1_7DA4))+
  geom_errorbar(aes(ymin=cell_death-SD, ymax=cell_death+SD), 
                color="grey80", width=.05)+
  geom_errorbar(aes(xmin=ddG_foldx_7DA4-ddG_foldx_SD_7DA4, xmax=ddG_foldx_7DA4+ddG_foldx_SD_7DA4),
                color="grey80", width=1)+
  geom_point(size=4, color="grey40")+
  geom_text_repel(size=4, color="grey40")+
  scale_color_manual(values=c("#b3003b", "grey"))+
  stat_cor(label.y = 30, label.x = 0, color="black")+ 
  labs(x="ddG 7DA4", y="% Cell Death (Li et al 2012)")+
  ylim(0, 30)+
  theme_classic()+
  theme(legend.position = "none")

scatter_4
ggsave(scatter_4, file="ddG_7DA4_Lietal2012.jpg", width=4, height=4, path=path)
ggsave(scatter_4, file="ddG_7DA4_Lietal2012.pdf", width=4, height=4, path=path)

#
scatter_5<-ggplot(singles_ddG_death, aes(x=ddG_foldx_7DAC, y=cell_death, label=ID, color=category_1_7DA4))+
  geom_errorbar(aes(ymin=cell_death-SD, ymax=cell_death+SD), 
                color="grey80", width=.05)+
  geom_errorbar(aes(xmin=ddG_foldx_7DAC-ddG_foldx_SD_7DAC, xmax=ddG_foldx_7DAC+ddG_foldx_SD_7DAC),
                color="grey80", width=1)+
  geom_point(size=4, color="grey40")+
  geom_text_repel(size=4, color="grey40")+
  scale_color_manual(values=c("#b3003b", "grey"))+
  stat_cor(label.y = 30, label.x = 0, color="black")+ 
  labs(x="ddG 7DAC", y="% Cell Death (Li et al 2012)")+
  ylim(0, 30)+
  theme_classic()+
  theme(legend.position = "none")

scatter_5
ggsave(scatter_5, file="ddG_7DAC_Lietal2012.jpg", width=4, height=4, path=path)
ggsave(scatter_5, file="ddG_7DAC_Lietal2012.pdf", width=4, height=4, path=path)

#
scatter_6<-ggplot(singles_ddG_death, aes(x=ddG_foldx_8Z94, y=cell_death, label=ID, color=category_1_8Z94))+
  geom_errorbar(aes(ymin=cell_death-SD, ymax=cell_death+SD), 
                color="grey80", width=.05)+
  geom_errorbar(aes(xmin=ddG_foldx_8Z94-ddG_foldx_SD_8Z94, xmax=ddG_foldx_8Z94+ddG_foldx_SD_8Z94),
                color="grey80", width=1)+
  geom_point(size=4, color="grey40")+
  geom_text_repel(size=4, color="grey40")+
  scale_color_manual(values=c("#b3003b", "grey"))+
  stat_cor(label.y = 30, label.x = 0, color="black")+ 
  labs(x="ddG 7DA4", y="% Cell Death (Li et al 2012)")+
  ylim(0, 30)+
  theme_classic()+
  theme(legend.position = "none")

scatter_6
ggsave(scatter_6, file="ddG_8Z94_Lietal2012.jpg", width=4, height=4, path=path)
ggsave(scatter_6, file="ddG_8Z94_Lietal2012.pdf", width=4, height=4, path=path)

###
# Cell-survival from Hu et al Cell Death Differentiation 2021
RIPK3_mutants_HeLa_2<-data.frame("ID"=c("WT", "N-464-D", "M-468-D"),
                                 "cell_survival"=c(12.19, 12.19, 7.53),
                                 "SD"=c(1.07, 3.22, 1.79))

singles_ddG_death<-inner_join(NS_ddG_RIPK3[c("ID", "ddG_foldx_7DA4", "ddG_foldx_SD_7DA4", "ddG_foldx_7DAC", "category_1_7DA4","ddG_foldx_SD_7DAC", "category_1_7DAC", "ddG_foldx_8Z94", "ddG_foldx_SD_8Z94", "category_1_8Z94")], 
                              RIPK3_mutants_HeLa_2, by="ID")

singles_ddG_death<-rbind(data.frame("ID"=c("WT"), "ddG_foldx_7DA4"=c(0), "ddG_foldx_SD_7DA4"=c(0), "category_1_7DA4"=c("WT-like"), "ddG_foldx_7DAC"=c(0), "ddG_foldx_SD_7DAC"=c(0), "category_1_7DAC"=c("WT-like"), "ddG_foldx_8Z94"=c(0), "ddG_foldx_SD_8Z94"=c(0), "category_1_8Z94"=c("WT-like"), 
                                    "cell_survival"=c(12.19), "SD"=c(1.07)), singles_ddG_death)

scatter_7<-ggplot(singles_ddG_death, aes(x=ddG_foldx_7DA4, y=cell_survival, label=ID, color=category_1_7DA4))+
  geom_errorbar(aes(ymin=cell_survival-SD, ymax=cell_survival+SD), 
                color="grey80", width=.05)+
  geom_errorbar(aes(xmin=ddG_foldx_7DA4-ddG_foldx_SD_7DA4, xmax=ddG_foldx_7DA4+ddG_foldx_SD_7DA4),
                color="grey80", width=1)+
  geom_point(size=4, color="grey40")+
  geom_text_repel(size=4, color="grey40")+
  scale_color_manual(values=c("#b3003b", "grey"))+
  stat_cor(color="black")+ 
  labs(x="ddG 7DA4", y="% Cell Survival (Hu et al 2021)")+
  ylim(0, 25)+
  theme_classic()+
  theme(legend.position = "none")

scatter_7
ggsave(scatter_7, file="ddG_7DA4_Huetal2021.jpg", width=4, height=4, path=path)
ggsave(scatter_7, file="ddG_7DA4_Huetal2021.pdf", width=4, height=4, path=path)

# 
scatter_8<-ggplot(singles_ddG_death, aes(x=ddG_foldx_7DAC, y=cell_survival, label=ID, color=category_1_7DAC))+
  geom_errorbar(aes(ymin=cell_survival-SD, ymax=cell_survival+SD), 
                color="grey80", width=.05)+
  geom_errorbar(aes(xmin=ddG_foldx_7DAC-ddG_foldx_SD_7DAC, xmax=ddG_foldx_7DAC+ddG_foldx_SD_7DAC),
                color="grey80", width=1)+
  geom_point(size=4, color="grey40")+
  geom_text_repel(size=4, color="grey40")+
  scale_color_manual(values=c("#b3003b", "grey"))+
  stat_cor(color="black")+ 
  labs(x="ddG 7DAC", y="% Cell Survival (Hu et al 2021)")+
  ylim(0, 25)+
  theme_classic()+
  theme(legend.position = "none")

scatter_8
ggsave(scatter_8, file="ddG_7DAC_Huetal2021.jpg", width=4, height=4, path=path)
ggsave(scatter_8, file="ddG_7DAC_Huetal2021.pdf", width=4, height=4, path=path)

# 
scatter_9<-ggplot(singles_ddG_death, aes(x=ddG_foldx_8Z94, y=cell_survival, label=ID, color=category_1_8Z94))+
  geom_errorbar(aes(ymin=cell_survival-SD, ymax=cell_survival+SD), 
                color="grey80", width=.05)+
  geom_errorbar(aes(xmin=ddG_foldx_8Z94-ddG_foldx_SD_8Z94, xmax=ddG_foldx_8Z94+ddG_foldx_SD_8Z94),
                color="grey80", width=1)+
  geom_point(size=4, color="grey40")+
  geom_text_repel(size=4, color="grey40")+
  scale_color_manual(values=c("#b3003b", "grey"))+
  stat_cor(color="black")+ 
  labs(x="ddG 8Z94", y="% Cell Survival (Hu et al 2021)")+
  ylim(0, 25)+
  theme_classic()+
  theme(legend.position = "none")

scatter_9
ggsave(scatter_9, file="ddG_8Z94_Huetal2021.jpg", width=4, height=4, path=path)
ggsave(scatter_9, file="ddG_8Z94_Huetal2021.pdf", width=4, height=4, path=path)

# Wu et al. 2021
# They also demonstrated by EM and NMR that mutants L-456-D and Q-449-D change the fibril structure
mRIP3_tht_wt<-c("L-456-D")
mRIP3_tht_dec<-c("Q-449-D", "F-442-D")

# Ma et al. 2025 
# They also showed that mutant L-466-E forms a different fibril structure by EM and IF.
RIPK3_mutants_HeLa<-data.frame("ID"=c("WT", "V-450-E", "V-460-E", "L-466-E"),
                               "cell_viability"=c(51.74, 80.23, 86.63, 96.22),
                               "SD"=c(3.49, 11.05, 5.23, 6.11))

singles_ddG_death<-inner_join(NS_ddG_RIPK3[c("ID", "ddG_foldx_7DA4", "ddG_foldx_SD_7DA4", "ddG_foldx_7DAC", "category_1_7DA4","ddG_foldx_SD_7DAC", "category_1_7DAC", "ddG_foldx_8Z94", "ddG_foldx_SD_8Z94", "category_1_8Z94")], 
                              RIPK3_mutants_HeLa, by="ID")

singles_ddG_death<-rbind(data.frame("ID"=c("WT"), "ddG_foldx_7DA4"=c(0), "ddG_foldx_SD_7DA4"=c(0), "category_1_7DA4"=c("WT-like"), "ddG_foldx_7DAC"=c(0), "ddG_foldx_SD_7DAC"=c(0), "category_1_7DAC"=c("WT-like"), "ddG_foldx_8Z94"=c(0), "ddG_foldx_SD_8Z94"=c(0), "category_1_8Z94"=c("WT-like"), 
                                    "cell_viability"=c(51.74), "SD"=c(3.49)), singles_ddG_death)

scatter_10<-ggplot(singles_ddG_death, aes(x=ddG_foldx_7DA4, y=cell_viability, label=ID, color=category_1_7DA4))+
  geom_errorbar(aes(ymin=cell_viability-SD, ymax=cell_viability+SD), 
                color="grey80", width=.05)+
  geom_errorbar(aes(xmin=ddG_foldx_7DA4-ddG_foldx_SD_7DA4, xmax=ddG_foldx_7DA4+ddG_foldx_SD_7DA4),
                color="grey80", width=1)+
  geom_point(size=4, color="grey40")+
  geom_text_repel(size=4, color="grey40")+
  scale_color_manual(values=c("#b3003b", "grey"))+
  stat_cor(label.y = 30, label.x = -4, color="black")+ 
  labs(x="ddG 7DA4", y="% Cell Viability (Ma et al 2025)")+
  ylim(0, 110)+
  theme_classic()+
  theme(legend.position = "none")

scatter_10
ggsave(scatter_10, file="ddG_7DA4_Maetal2025.jpg", width=4, height=4, path=path)
ggsave(scatter_10, file="ddG_7DA4_Maetal2025.pdf", width=4, height=4, path=path)
#

scatter_11<-ggplot(singles_ddG_death, aes(x=ddG_foldx_7DAC, y=cell_viability, label=ID, color=category_1_7DAC))+
  geom_errorbar(aes(ymin=cell_viability-SD, ymax=cell_viability+SD), 
                color="grey80", width=.05)+
  geom_errorbar(aes(xmin=ddG_foldx_7DAC-ddG_foldx_SD_7DAC, xmax=ddG_foldx_7DAC+ddG_foldx_SD_7DAC),
                color="grey80", width=1)+
  geom_point(size=4, color="grey40")+
  geom_text_repel(size=4, color="grey40")+
  scale_color_manual(values=c("#b3003b", "grey"))+
  stat_cor(label.y = 30, label.x = -4, color="black")+ 
  labs(x="ddG 7DAC", y="% Cell Viability (Ma et al 2025)")+
  ylim(0, 110)+
  theme_classic()+
  theme(legend.position = "none")

scatter_11
ggsave(scatter_11, file="ddG_7DAC_Maetal2025.jpg", width=4, height=4, path=path)
ggsave(scatter_11, file="ddG_7DAC_Maetal2025.pdf", width=4, height=4, path=path)
#

scatter_12<-ggplot(singles_ddG_death, aes(x=ddG_foldx_8Z94, y=cell_viability, label=ID, color=category_1_8Z94))+
  geom_errorbar(aes(ymin=cell_viability-SD, ymax=cell_viability+SD), 
                color="grey80", width=.05)+
  geom_errorbar(aes(xmin=ddG_foldx_8Z94-ddG_foldx_SD_8Z94, xmax=ddG_foldx_8Z94+ddG_foldx_SD_8Z94),
                color="grey80", width=1)+
  geom_point(size=4, color="grey40")+
  geom_text_repel(size=4, color="grey40")+
  scale_color_manual(values=c("#b3003b", "grey"))+
  stat_cor(label.y = 30, label.x = -4, color="black")+ 
  labs(x="ddG 8Z94", y="% Cell Viability (Ma et al 2025)")+
  ylim(0, 110)+
  theme_classic()+
  theme(legend.position = "none")

scatter_12
ggsave(scatter_12, file="ddG_8Z94_Maetal2025.jpg", width=4, height=4, path=path)
ggsave(scatter_12, file="ddG_8Z94_Maetal2025.pdf", width=4, height=4, path=path)

################################################################################
# 4: ddG vs NS with mtuants from literature

#CD = cell death
CD_dec<-c("V-450-E", "V-460-E", "L-466-E", "V-460-P", "V-458-P", "G-457-D", 
          "V-458-E", "G-461-Y", "N-464-I", "V-450-H", "L-466-E", "N-464-H")
CD_WT<-c("N-464-D", "M-468-D", "D-462-A", "L-449-A", "M-468-A", "L-466-Y", "Y-465-V", "I-452-S", "G-457-V")
CD_freq<-c("P-448-L", "V-450-I", "V-460-I", "R-447-Q")

NS_ddG_RIPK3$CD<-""
NS_ddG_RIPK3[NS_ddG_RIPK3$ID %in% CD_dec, "CD"]<-"Decreased"
NS_ddG_RIPK3[NS_ddG_RIPK3$ID %in% CD_WT, "CD"]<-"WT-like"
NS_ddG_RIPK3[NS_ddG_RIPK3$ID %in% CD_freq, "CD"]<-"AF>1e-05"

NS_ddG_RIPK3$AF<-""
NS_ddG_RIPK3[NS_ddG_RIPK3$ID %in% CD_freq, "AF"]<-"AF>1e-05"

NS_ddG_RIPK3$ID_reported<-""
NS_ddG_RIPK3[NS_ddG_RIPK3$CD != "", "ID_reported"]<-NS_ddG_RIPK3[NS_ddG_RIPK3$CD != "", "ID"]

NS_ddG_RIPK3$class<-"no_data"
NS_ddG_RIPK3[NS_ddG_RIPK3$CD ==  "Decreased", "class"]<- "Decreased"
NS_ddG_RIPK3[NS_ddG_RIPK3$CD ==  "WT-like", "class"]<- "WT-like"
NS_ddG_RIPK3[NS_ddG_RIPK3$CD ==  "AF>1e-05", "class"]<- "WT-like"

scatter_13<-ggplot(NS_ddG_RIPK3, aes(x=nscore_c, y=ddG_foldx_7DAC, label=ID_reported))+
  geom_hline(yintercept=0, size=.1)+
  geom_vline(xintercept=0, size=.1)+
  geom_point(data=NS_ddG_RIPK3[NS_ddG_RIPK3$ID_reported == "",], alpha=0.3, color="grey60")+
  geom_point(data=NS_ddG_RIPK3[NS_ddG_RIPK3$ID_reported != "",], aes(color=CD, shape=AF), size=2)+
  geom_text_repel(max.overlaps = 100)+
  scale_color_manual(values=c("grey40", "red3", "grey20"))+
  scale_shape_manual(values=c(16, 17))+
  theme_classic()
#theme(legend.position = "none")

scatter_13

ggsave(scatter_13, file="ddG_7DAC_NS_CD_all.jpg", width=5, height=4, path=path)
ggsave(scatter_13, file="ddG_7DAC_NS_CD_all.pdf", width=5, height=4, path=path)

###
# Let's measure the distance of all the dots to the origin

#NS_ddG_RIPK3<-NS_ddG_RIPK3 %>% drop_na(Pos)
# measure the length of the vector to the origin
distance_to_origin <- function(x, y){
  distance <- sqrt(x^2 + y^2)
  return(distance)
}

# iterate over the dataframe rows
# get the nscore_c (x) and ddG_foldx (y)
# calculates the distance to the origin
# store it in the column
for (i in 1:nrow(NS_ddG_RIPK3)){
  nscore_c<-NS_ddG_RIPK3[i, "nscore_c"]
  ddG_foldx_7DAC<-NS_ddG_RIPK3[i, "ddG_foldx_7DAC"]
  distance<-distance_to_origin(nscore_c, ddG_foldx_7DAC)
  NS_ddG_RIPK3[i, "distance_to_origin_7DAC"]<-distance
}

# Histogram of the distances to the origin

p_hist <- ggplot(NS_ddG_RIPK3[NS_ddG_RIPK3$class != "no_data",], 
                 aes(x=distance_to_origin_7DAC, fill=factor(class, levels=c("WT-like", "Decreased"))))+
  #geom_histogram(bins=100)+
  geom_density(alpha=0.5)+
  scale_fill_manual(values=c("grey20","red3"), guide="none")+
  labs(x="Distance to origin", y="Density")+
  theme_classic()
p_hist

ggsave(p_hist, file="Necroptosis_class_hist_7DAC_all.jpg", width = 4, height = 4, path=path)
ggsave(p_hist, file="Necroptosis_class_hist_7DAC_all.pdf", width = 4, height = 4, path=path)


# compare the distributions with Kolmogorov-Smirnov Test
ks_result<-ks.test(NS_ddG_RIPK3[NS_ddG_RIPK3$class == "WT-like", "distance_to_origin_7DAC"], 
                   NS_ddG_RIPK3[NS_ddG_RIPK3$class == "Decreased", "distance_to_origin_7DAC"])

p_ecdf<-ggplot(NS_ddG_RIPK3[NS_ddG_RIPK3$class!="no_data",], 
               aes(x=distance_to_origin_7DAC, color=factor(class, levels=c("WT-like", "Decreased")))) +
  stat_ecdf(size=1)+
  annotate("text", x=15, y=0.1, label="Kolmogorov-Smirnov Test")+
  annotate("text", x=15, y=0.02, label=paste("p-value:", round(ks_result$p.value, 2)))+
  scale_color_manual("Cell Death", values=c("grey20","red3"))+
  labs(x="Distance to origin", y="ECDF")+
  theme_classic()
p_ecdf

ggsave(p_ecdf, file="Necroptosis_class_ecdf_7DAC_all.jpg", width = 6, height = 4, path=path)
ggsave(p_ecdf, file="Necroptosis_class_ecdf_7DAC_all.pdf", width = 6, height = 4, path=path)

# Heatmap # sweet-spot
scatter_14<-ggplot(NS_ddG_RIPK3, aes(x=nscore_c, y=ddG_foldx_7DAC))+
  stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = F, bins=10) +
  scale_fill_gradientn(colours = c("white", "orange2", "red3"), name = "Density") +
  theme_pubclean()+
  theme(legend.position = "none")

scatter_14

ggsave(scatter_14, file="ddG_7DAC_NS_sweetspot.jpg", width=4, height=4, path=path)
ggsave(scatter_14, file="ddG_7DAC_NS_sweetspot.pdf", width=4, height=4, path=path)


##
scatter_15<-ggplot(NS_ddG_RIPK3, aes(x=nscore_c, y=ddG_foldx_7DA4, label=ID_reported))+
  geom_hline(yintercept=0, size=.1)+
  geom_vline(xintercept=0, size=.1)+
  geom_point(data=NS_ddG_RIPK3[NS_ddG_RIPK3$ID_reported == "",], alpha=0.3, color="grey60")+
  geom_point(data=NS_ddG_RIPK3[NS_ddG_RIPK3$ID_reported != "",], aes(color=CD, shape=AF), size=2)+
  geom_text_repel(max.overlaps = 100)+
  scale_color_manual(values=c("grey40", "red3", "grey20"))+
  scale_shape_manual(values=c(16, 17))+
  theme_classic()
#theme(legend.position = "none")

scatter_15

ggsave(scatter_15, file="ddG_7DA4_NS_CD_all.jpg", width=5, height=4, path=path)
ggsave(scatter_15, file="ddG_7DA4_NS_CD_all.pdf", width=5, height=4, path=path)

# iterate over the dataframe rows
# get the nscore_c (x) and ddG_foldx (y)
# calculates the distance to the origin
# store it in the column
for (i in 1:nrow(NS_ddG_RIPK3)){
  nscore_c<-NS_ddG_RIPK3[i, "nscore_c"]
  ddG_foldx_7DA4<-NS_ddG_RIPK3[i, "ddG_foldx_7DA4"]
  distance<-distance_to_origin(nscore_c, ddG_foldx_7DA4)
  NS_ddG_RIPK3[i, "distance_to_origin_7DA4"]<-distance
}

# Histogram of the distances to the origin

p_hist <- ggplot(NS_ddG_RIPK3[NS_ddG_RIPK3$class!="no_data",], 
                 aes(x=distance_to_origin_7DA4, fill=factor(class, levels=c("WT-like", "Decreased"))))+
  #geom_histogram(bins=100)+
  geom_density(alpha=0.5)+
  scale_fill_manual(values=c("grey20","red3"), guide="none")+
  labs(x="Distance to origin", y="Density")+
  theme_classic()
p_hist

ggsave(p_hist, file="Necroptosis_class_hist_7DA4_all.jpg", width = 4, height = 4, path=path)
ggsave(p_hist, file="Necroptosis_class_hist_7DA4_all.pdf", width = 4, height = 4, path=path)


# compare the distributions with Kolmogorov-Smirnov Test
ks_result<-ks.test(NS_ddG_RIPK3[NS_ddG_RIPK3$class == "WT-like", "distance_to_origin_7DA4"], 
                   NS_ddG_RIPK3[NS_ddG_RIPK3$class == "Decreased", "distance_to_origin_7DA4"])

p_ecdf<-ggplot(NS_ddG_RIPK3[NS_ddG_RIPK3$class!="no_data",], 
               aes(x=distance_to_origin_7DA4, color=factor(class, levels=c("WT-like", "Decreased")))) +
  stat_ecdf(size=1)+
  annotate("text", x=15, y=0.1, label="Kolmogorov-Smirnov Test")+
  annotate("text", x=15, y=0.02, label=paste("p-value:", round(ks_result$p.value, 2)))+
  scale_color_manual("Cell Death", values=c("grey20","red3"))+
  labs(x="Distance to origin", y="ECDF")+
  theme_classic()
p_ecdf

ggsave(p_ecdf, file="Necroptosis_class_ecdf_7DA4_all.jpg", width = 6, height = 4, path=path)
ggsave(p_ecdf, file="Necroptosis_class_ecdf_7DA4_all.pdf", width = 6, height = 4, path=path)


##
scatter_16<-ggplot(NS_ddG_RIPK3, aes(x=nscore_c, y=ddG_foldx_8Z94, label=ID_reported))+
  geom_hline(yintercept=0, size=.1)+
  geom_vline(xintercept=0, size=.1)+
  geom_point(data=NS_ddG_RIPK3[NS_ddG_RIPK3$ID_reported == "",], alpha=0.3, color="grey60")+
  geom_point(data=NS_ddG_RIPK3[NS_ddG_RIPK3$ID_reported != "",], aes(color=CD, shape=AF), size=2)+
  geom_text_repel(max.overlaps = 100)+
  scale_color_manual(values=c("grey40", "red3", "grey20"))+
  scale_shape_manual(values=c(16, 17))+
  theme_classic()
#theme(legend.position = "none")

scatter_16

ggsave(scatter_16, file="ddG_8Z94_NS_CD_all.jpg", width=5, height=4, path=path)
ggsave(scatter_16, file="ddG_8Z94_NS_CD_all.pdf", width=5, height=4, path=path)


# iterate over the dataframe rows
# get the nscore_c (x) and ddG_foldx (y)
# calculates the distance to the origin
# store it in the column
for (i in 1:nrow(NS_ddG_RIPK3)){
  nscore_c<-NS_ddG_RIPK3[i, "nscore_c"]
  ddG_foldx_8Z94<-NS_ddG_RIPK3[i, "ddG_foldx_8Z94"]
  distance<-distance_to_origin(nscore_c, ddG_foldx_8Z94)
  NS_ddG_RIPK3[i, "distance_to_origin_8Z94"]<-distance
}

# Histogram of the distances to the origin

p_hist <- ggplot(NS_ddG_RIPK3[NS_ddG_RIPK3$class!="no_data",], 
                 aes(x=distance_to_origin_8Z94, fill=factor(class, levels=c("WT-like", "Decreased"))))+
  #geom_histogram(bins=100)+
  geom_density(alpha=0.5)+
  scale_fill_manual(values=c("grey20","red3"), guide="none")+
  labs(x="Distance to origin", y="Density")+
  theme_classic()
p_hist

ggsave(p_hist, file="Necroptosis_class_hist_8Z94_all.jpg", width = 4, height = 4, path=path)
ggsave(p_hist, file="Necroptosis_class_hist_8Z94_all.pdf", width = 4, height = 4, path=path)


# compare the distributions with Kolmogorov-Smirnov Test
ks_result<-ks.test(NS_ddG_RIPK3[NS_ddG_RIPK3$class == "WT-like", "distance_to_origin_8Z94"], 
                   NS_ddG_RIPK3[NS_ddG_RIPK3$class == "Decreased", "distance_to_origin_8Z94"])

p_ecdf<-ggplot(NS_ddG_RIPK3[NS_ddG_RIPK3$class!="no_data",], 
               aes(x=distance_to_origin_8Z94, color=factor(class, levels=c("WT-like", "Decreased")))) +
  stat_ecdf(size=1)+
  annotate("text", x=15, y=0.1, label="Kolmogorov-Smirnov Test")+
  annotate("text", x=15, y=0.02, label=paste("p-value:", round(ks_result$p.value, 2)))+
  scale_color_manual("Cell Death", values=c("grey20","red3"))+
  labs(x="Distance to origin", y="ECDF")+
  theme_classic()
p_ecdf

ggsave(p_ecdf, file="Necroptosis_class_ecdf_8Z94_all.jpg", width = 6, height = 4, path=path)
ggsave(p_ecdf, file="Necroptosis_class_ecdf_8Z94_all.pdf", width = 6, height = 4, path=path)
