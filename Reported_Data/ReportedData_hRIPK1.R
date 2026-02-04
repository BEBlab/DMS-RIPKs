library(tidyverse)
library(ggpubr)
library(ggsignif)
require("ggrepel")
library(ggbreak) 

# your folder path
setwd("")

dir.create("02_ReportedData_hRIPK1")
path="02_ReportedData_hRIPK1"

### hRIPK1
# Nucleation Score
load("nscore_df_hRIPK1.RData")

singles_stops$sig_1<-F
singles_stops[singles_stops$p.adjust<0.01, "sig_1"]<-T

singles_stops$category_1<-"WT-like"
singles_stops[singles_stops$p.adjust<0.01 & singles_stops$nscore_c >0, "category_1"]<-"NS_inc"
singles_stops[singles_stops$p.adjust<0.01 & singles_stops$nscore_c <0, "category_1"]<-"NS_dec"

singles_stops<-singles_stops[singles_stops$low_sigma == T | singles_stops$sig_1 == T, ]

singles_stops$Residue<-paste0(singles_stops$WT_AA, singles_stops$Pos)

### ddG
#PDB 9HR6
load("ddG_9HR6.RData")

ddG_foldx_9HR6$sig_1_9HR6<-F
ddG_foldx_9HR6[ddG_foldx_9HR6$p.adjust_9HR6<0.01, "sig_1_9HR6"]<-T

ddG_foldx_9HR6$category_1_9HR6<-"WT-like"
ddG_foldx_9HR6[ddG_foldx_9HR6$p.adjust_9HR6<0.01 & ddG_foldx_9HR6$ddG_foldx_9HR6 >2, "category_1_9HR6"]<-"ddG_inc"
ddG_foldx_9HR6[ddG_foldx_9HR6$p.adjust_9HR6<0.01 & ddG_foldx_9HR6$ddG_foldx_9HR6 <0, "category_1_9HR6"]<-"ddG_dec"

ddG_foldx_9HR6<-ddG_foldx_9HR6[ddG_foldx_9HR6$low_ddG_foldx_SD_9HR6 == T | ddG_foldx_9HR6$sig_1_9HR6 == T, ]


#PDB 9HR9
load("ddG_9HR9.RData")

ddG_foldx_9HR9$sig_1_9HR9<-F
ddG_foldx_9HR9[ddG_foldx_9HR9$p.adjust_9HR9<0.01, "sig_1_9HR9"]<-T

ddG_foldx_9HR9$category_1_9HR9<-"WT-like"
ddG_foldx_9HR9[ddG_foldx_9HR9$p.adjust_9HR9<0.01 & ddG_foldx_9HR9$ddG_foldx_9HR9 >2, "category_1_9HR9"]<-"ddG_inc"
ddG_foldx_9HR9[ddG_foldx_9HR9$p.adjust_9HR9<0.01 & ddG_foldx_9HR9$ddG_foldx_9HR9 <0, "category_1_9HR9"]<-"ddG_dec"

ddG_foldx_9HR9<-ddG_foldx_9HR9[ddG_foldx_9HR9$low_ddG_foldx_SD_9HR9 == T | ddG_foldx_9HR9$sig_1_9HR9 == T, ]


#PDB 8Z93
load("ddG_8Z93.RData")

ddG_foldx_8Z93$sig_1_8Z93<-F
ddG_foldx_8Z93[is.na(ddG_foldx_8Z93$p.adjust_8Z93), "p.adjust_8Z93"]<-0
ddG_foldx_8Z93[ddG_foldx_8Z93$p.adjust_8Z93<0.01, "sig_1_8Z93"]<-T

ddG_foldx_8Z93$category_1_8Z93<-"WT-like"
ddG_foldx_8Z93[ddG_foldx_8Z93$p.adjust_8Z93<0.01 & ddG_foldx_8Z93$ddG_foldx_8Z93 >2, "category_1_8Z93"]<-"ddG_inc"
ddG_foldx_8Z93[ddG_foldx_8Z93$p.adjust_8Z93<0.01 & ddG_foldx_8Z93$ddG_foldx_8Z93 <0, "category_1_8Z93"]<-"ddG_dec"

ddG_foldx_8Z93<-ddG_foldx_8Z93[ddG_foldx_8Z93$low_ddG_foldx_SD_8Z93 == T | ddG_foldx_8Z93$sig_1_8Z93 == T, ]

###
NS_ddG_RIPK1<-full_join(ddG_foldx_9HR9[c("ID", "ddG_foldx_9HR9", "ddG_foldx_SD_9HR9", "category_1_9HR9")], 
                        ddG_foldx_9HR6[c("ID", "ddG_foldx_9HR6", "ddG_foldx_SD_9HR6", "category_1_9HR6")], by="ID")

NS_ddG_RIPK1<-full_join(NS_ddG_RIPK1, ddG_foldx_8Z93[c("ID", "ddG_foldx_8Z93", "ddG_foldx_SD_8Z93", "category_1_8Z93")], by="ID")

NS_ddG_RIPK1<-full_join(NS_ddG_RIPK1, singles_stops[c("ID", "Pos", "Residue", "nscore_c", "sigma", "category_10", "category_1")], by="ID")

################################################################################
# 1: NS Correlation vs mutants from the literature
# Cell death from Li et al Cell 2012
RIPK1_mutants_HeLa_1<-data.frame("ID"=c("WT", "I-539-P", "I-541-P", "N-545-P"),
                                 "cell_death"=c(68.38, 8.94, 21, 31.73),
                                 "SD"=c(0.89, 2.68, 2.69, 5.36))

singles_stops$category_1<-"WT-like"
singles_stops[singles_stops$p.adjust<0.01 & singles_stops$nscore_c >0, "category_1"]<-"NS_inc"
singles_stops[singles_stops$p.adjust<0.01 & singles_stops$nscore_c <0, "category_1"]<-"NS_dec"

singles_ns_death<-inner_join(singles_stops[c("ID", "nscore_c", "sigma", "category_1")], RIPK1_mutants_HeLa_1, by="ID")

singles_ns_death<-rbind(data.frame("ID"=c("WT"), "nscore_c"=c(0.115299498), "sigma"=c(0.06311079), "category_1"=c("WT-like"),
                                   "cell_death"=c(68.38), "SD"=c(0.89)), singles_ns_death)

scatter_1<-ggplot(singles_ns_death, aes(x=nscore_c, y=cell_death, label=ID, color=category_1))+
  geom_errorbar(aes(ymin=cell_death-SD, ymax=cell_death+SD), 
                color="grey80", width=.05)+
  geom_errorbar(aes(xmin=nscore_c-sigma, xmax=nscore_c+sigma),
                color="grey80", width=1)+
  geom_point(size=4, color="grey40")+
  #scale_color_manual(values=c("#CD3333", "grey", "#00008B"))+
  geom_text_repel(size=4, color="grey40")+
  stat_cor(color="black")+ 
  ylim(0, 70)+
  labs(x="Nucleation Score", y="% Cell Death (Li et al 2012)")+
  theme_classic()+
  theme(legend.position = "none")

scatter_1
ggsave(scatter_1, file="nscore_Lietal2012.jpg", width=4, height=4, path=path)
ggsave(scatter_1, file="nscore_Lietal2012.pdf", width=4, height=4, path=path)

################################################################################
# 2: ddG Correlation vs mutants from the literature

RIPK1_mutants_HeLa_1<-data.frame("ID"=c("WT", "I-539-P", "I-541-P", "N-545-P"),
                                 "cell_death"=c(68.38, 8.94, 21, 31.73),
                                 "SD"=c(0.89, 2.68, 2.69, 5.36))

singles_ddG_death<-inner_join(NS_ddG_RIPK1[c("ID", "ddG_foldx_9HR9", "ddG_foldx_SD_9HR9", "ddG_foldx_9HR6", "category_1_9HR9","ddG_foldx_SD_9HR6", "category_1_9HR6", "ddG_foldx_8Z93", "ddG_foldx_SD_8Z93", "category_1_8Z93")], 
                              RIPK1_mutants_HeLa_1, by="ID")

singles_ddG_death<-rbind(data.frame("ID"=c("WT"), "ddG_foldx_9HR9"=c(0), "ddG_foldx_SD_9HR9"=c(0), "category_1_9HR9"=c("WT-like"), "ddG_foldx_9HR6"=c(0), "ddG_foldx_SD_9HR6"=c(0), "category_1_9HR6"=c("WT-like"), "ddG_foldx_8Z93"=c(0), "ddG_foldx_SD_8Z93"=c(0), "category_1_8Z93"=c("WT-like"), 
                                    "cell_death"=c(68.38), "SD"=c(0.89)), singles_ddG_death)

scatter_4<-ggplot(singles_ddG_death, aes(x=ddG_foldx_9HR9, y=cell_death, label=ID, color=category_1_9HR9))+
  geom_errorbar(aes(ymin=cell_death-SD, ymax=cell_death+SD), 
                color="grey80", width=.05)+
  geom_errorbar(aes(xmin=ddG_foldx_9HR9-ddG_foldx_SD_9HR9, xmax=ddG_foldx_9HR9+ddG_foldx_SD_9HR9),
                color="grey80", width=1)+
  geom_point(size=4, color="grey40")+
  geom_text_repel(size=4, color="grey40")+
  scale_color_manual(values=c("#b3003b", "grey"))+
  stat_cor(color="black")+ 
  labs(x="ddG 9HR9", y="% Cell Death (Li et al 2012)")+
  ylim(0, 70)+
  theme_classic()+
  theme(legend.position = "none")

scatter_4
ggsave(scatter_4, file="ddG_9HR9_Lietal2012.jpg", width=4, height=4, path=path)
ggsave(scatter_4, file="ddG_9HR9_Lietal2012.pdf", width=4, height=4, path=path)

#
scatter_5<-ggplot(singles_ddG_death, aes(x=ddG_foldx_9HR6, y=cell_death, label=ID, color=category_1_9HR9))+
  geom_errorbar(aes(ymin=cell_death-SD, ymax=cell_death+SD), 
                color="grey80", width=.05)+
  geom_errorbar(aes(xmin=ddG_foldx_9HR6-ddG_foldx_SD_9HR6, xmax=ddG_foldx_9HR6+ddG_foldx_SD_9HR6),
                color="grey80", width=1)+
  geom_point(size=4, color="grey40")+
  geom_text_repel(size=4, color="grey40")+
  scale_color_manual(values=c("#b3003b", "grey"))+
  stat_cor(color="black")+ 
  labs(x="ddG 9HR6", y="% Cell Death (Li et al 2012)")+
  ylim(0, 70)+
  theme_classic()+
  theme(legend.position = "none")

scatter_5
ggsave(scatter_5, file="ddG_9HR6_Lietal2012.jpg", width=4, height=4, path=path)
ggsave(scatter_5, file="ddG_9HR6_Lietal2012.pdf", width=4, height=4, path=path)

#
scatter_6<-ggplot(singles_ddG_death, aes(x=ddG_foldx_8Z93, y=cell_death, label=ID, color=category_1_8Z93))+
  geom_errorbar(aes(ymin=cell_death-SD, ymax=cell_death+SD), 
                color="grey80", width=.05)+
  geom_errorbar(aes(xmin=ddG_foldx_8Z93-ddG_foldx_SD_8Z93, xmax=ddG_foldx_8Z93+ddG_foldx_SD_8Z93),
                color="grey80", width=1)+
  geom_point(size=4, color="grey40")+
  geom_text_repel(size=4, color="grey40")+
  scale_color_manual(values=c("#b3003b", "grey"))+
  stat_cor(color="black")+ 
  labs(x="ddG 8Z93", y="% Cell Death (Li et al 2012)")+
  ylim(0, 70)+
  theme_classic()+
  theme(legend.position = "none")

scatter_6
ggsave(scatter_6, file="ddG_8Z93_Lietal2012.jpg", width=4, height=4, path=path)
ggsave(scatter_6, file="ddG_8Z93_Lietal2012.pdf", width=4, height=4, path=path)
