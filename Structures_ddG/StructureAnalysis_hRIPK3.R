library(tidyverse)
library(ggpubr)
library(ggsignif)
require("ggrepel")
library(ggbreak) 
library(reshape2)

# your folder path
setwd("")

dir.create("02_StructuralComparison_hRIPK3")
path="02_StructuralComparison_hRIPK3"

### hRIPK3
# Nucleation Score
load("nscore_df_hRIPK3.RData")

singles_stops$Residue<-paste0(singles_stops$WT_AA, singles_stops$Pos)

### ddG
#PDB 7DAC
ddG_foldx_7DAC<-read.csv("ddG_foldx_7DAC_3monomers_compact.csv")
ddG_foldx_7DAC$Pos<-ddG_foldx_7DAC$Pos+446 #correct the numbering
ddG_foldx_7DAC<-ddG_foldx_7DAC[c("ID", "WT_AA", "Pos", "Mut", "ddG_foldx", "ddG_foldx_SD")] 
ddG_foldx_7DAC$ID<-paste0(ddG_foldx_7DAC$WT_AA, "-", ddG_foldx_7DAC$Pos, "-", ddG_foldx_7DAC$Mut) 
ddG_foldx_7DAC<-rename(ddG_foldx_7DAC, ddG_foldx_7DAC = ddG_foldx, ddG_foldx_SD_7DAC = ddG_foldx_SD)
ddG_foldx_7DAC$Residue<-paste0(ddG_foldx_7DAC$WT_AA, ddG_foldx_7DAC$Pos)

# Residues classification based on side chain orientation
ddG_foldx_7DAC$orientation_7DAC<-"exposed"
ddG_foldx_7DAC[ddG_foldx_7DAC$Pos %in% c(448, 450, 452, 455, 456, 458, 459, 460, 462, 466, 468, 470), "orientation_7DAC"]<-"buried"
ddG_foldx_7DAC[ddG_foldx_7DAC$Pos %in% c(457, 461), "orientation_7DAC"]<-"glycine"


#PDB 7DA4
ddG_foldx_7DA4<-read.csv("ddG_foldx_7DA4_compact.csv")
ddG_foldx_7DA4<-ddG_foldx_7DA4[c("ID", "WT_AA", "Pos", "Mut", "ddG_foldx", "ddG_foldx_SD")]
ddG_foldx_7DA4$ID<-paste0(ddG_foldx_7DA4$WT_AA, "-", ddG_foldx_7DA4$Pos, "-", ddG_foldx_7DA4$Mut)
ddG_foldx_7DA4<-rename(ddG_foldx_7DA4, ddG_foldx_7DA4 = ddG_foldx, ddG_foldx_SD_7DA4 = ddG_foldx_SD)
ddG_foldx_7DA4$Residue<-paste0(ddG_foldx_7DA4$WT_AA, ddG_foldx_7DA4$Pos)

# Residues classification based on side chain orientation
ddG_foldx_7DA4$orientation_7DA4<-"exposed"
ddG_foldx_7DA4[ddG_foldx_7DA4$Pos %in% c(450, 452, 454, 456, 458, 459, 460, 466, 468), "orientation_7DA4"]<-"buried"
ddG_foldx_7DA4[ddG_foldx_7DA4$Pos %in% c(457, 461), "orientation_7DA4"]<-"glycine"


#PDB 8Z94
ddG_foldx_8Z94<-read.csv("ddG_foldx_8Z94_3monomers_compact.csv")
ddG_foldx_8Z94<-ddG_foldx_8Z94[c("ID", "WT_AA", "Pos", "Mut", "ddG_foldx", "ddG_foldx_SD")]
ddG_foldx_8Z94$ID<-paste0(ddG_foldx_8Z94$WT_AA, "-", ddG_foldx_8Z94$Pos, "-", ddG_foldx_8Z94$Mut)
ddG_foldx_8Z94<-rename(ddG_foldx_8Z94, ddG_foldx_8Z94 = ddG_foldx, ddG_foldx_SD_8Z94 = ddG_foldx_SD)
ddG_foldx_8Z94$Residue<-paste0(ddG_foldx_8Z94$WT_AA, ddG_foldx_8Z94$Pos)

# Residues classification based on side chain orientation
ddG_foldx_8Z94$orientation_8Z94<-"exposed" 
ddG_foldx_8Z94[ddG_foldx_8Z94$Pos %in% c(450, 452, 454, 456, 458, 459, 460, 466, 468), "orientation_8Z94"]<-"buried"
ddG_foldx_8Z94[ddG_foldx_8Z94$Pos %in% c(457, 461), "orientation_8Z94"]<-"glycine"

###
NS_ddG_RIPK3<-full_join(ddG_foldx_7DA4[c("ID", "Pos", "Residue", "ddG_foldx_7DA4", "ddG_foldx_SD_7DA4", "orientation_7DA4")], 
                        ddG_foldx_7DAC[c("ID", "ddG_foldx_7DAC", "ddG_foldx_SD_7DAC", "orientation_7DAC")], by="ID")

NS_ddG_RIPK3<-full_join(NS_ddG_RIPK3, ddG_foldx_8Z94[c("ID", "ddG_foldx_8Z94", "ddG_foldx_SD_8Z94", "orientation_8Z94")], by="ID")

NS_ddG_RIPK3<-full_join(NS_ddG_RIPK3, singles_stops[c("ID", "nscore_c", "sigma")], by="ID")

# Classify the residues based on the comparison between structures
NS_ddG_RIPK3$orientation<-"exposed"
NS_ddG_RIPK3[NS_ddG_RIPK3$Pos %in% c(450, 452, 456, 458, 459, 460, 466, 468), "orientation"]<-"buried"
NS_ddG_RIPK3[NS_ddG_RIPK3$Pos %in% c(457, 461), "orientation"]<-"glycine"
NS_ddG_RIPK3[NS_ddG_RIPK3$Pos %in% c(448, 454, 455, 462), "orientation"]<-"flipped"

# drop nas
NS_ddG_RIPK3_<- NS_ddG_RIPK3 %>% drop_na(ddG_foldx_7DA4, ddG_foldx_7DAC, ddG_foldx_8Z94)

# Classify the mutations based on the ddG
NS_ddG_RIPK3_$ddG_class<-"opposite_effect"
NS_ddG_RIPK3_[NS_ddG_RIPK3_$ddG_foldx_7DA4>2 & NS_ddG_RIPK3_$ddG_foldx_7DAC>2 & NS_ddG_RIPK3_$ddG_foldx_8Z94>2, "ddG_class"]<-"destabilizing_effect"
NS_ddG_RIPK3_[NS_ddG_RIPK3_$ddG_foldx_7DA4<2 & NS_ddG_RIPK3_$ddG_foldx_7DAC<2 & NS_ddG_RIPK3_$ddG_foldx_8Z94<2, "ddG_class"]<-"stabilizing_effect"

# label the mutations with opposite effects (at least in one structure)
NS_ddG_RIPK3_$ddG_class_ID<-""
NS_ddG_RIPK3_[NS_ddG_RIPK3_$ddG_class == "opposite_effect", "ddG_class_ID"]<-NS_ddG_RIPK3_[NS_ddG_RIPK3_$ddG_class == "opposite_effect", "ID"]

# Map the structural elements
NS_ddG_RIPK3_$structure<-"Disordered-1"
NS_ddG_RIPK3_[NS_ddG_RIPK3_$Pos>=450, "structure"]<-"β-strand-1"
NS_ddG_RIPK3_[NS_ddG_RIPK3_$Pos>=453, "structure"]<-"Loop-1"
NS_ddG_RIPK3_[NS_ddG_RIPK3_$Pos>=458, "structure"]<-"β-strand-2"
NS_ddG_RIPK3_[NS_ddG_RIPK3_$Pos>=461, "structure"]<-"Loop-2"
NS_ddG_RIPK3_[NS_ddG_RIPK3_$Pos>=466, "structure"]<-"β-strand-3"
NS_ddG_RIPK3_[NS_ddG_RIPK3_$Pos>=469, "structure"]<-"Disordered-2"

NS_ddG_RIPK3_$strtype<-"Disordered"
NS_ddG_RIPK3_[NS_ddG_RIPK3_$Pos>=450, "strtype"]<-"β-strand"
NS_ddG_RIPK3_[NS_ddG_RIPK3_$Pos>=453, "strtype"]<-"Loop"
NS_ddG_RIPK3_[NS_ddG_RIPK3_$Pos>=458, "strtype"]<-"β-strand"
NS_ddG_RIPK3_[NS_ddG_RIPK3_$Pos>=461, "strtype"]<-"Loop"
NS_ddG_RIPK3_[NS_ddG_RIPK3_$Pos>=466, "strtype"]<-"β-strand"
NS_ddG_RIPK3_[NS_ddG_RIPK3_$Pos>=469, "strtype"]<-"Disordered"

NS_ddG_RIPK3_$interface<-"surface"
NS_ddG_RIPK3_[NS_ddG_RIPK3_$Pos %in% c(450, 452, 458, 460), "interface"]<-"1st-Interface"
NS_ddG_RIPK3_[NS_ddG_RIPK3_$Pos %in% c(456, 459, 466, 468), "interface"]<-"2nd-Interface"

################################################################################
# 1: compare the ddG between the different structures
scatterddG_1<-ggplot(NS_ddG_RIPK3, aes(x=ddG_foldx_7DA4, y=ddG_foldx_7DAC, label=ID))+
  geom_vline(xintercept=2, linetype="dashed", color="grey80")+
  geom_vline(xintercept=0, color="grey50")+
  geom_hline(yintercept=2, linetype="dashed", color="grey80")+
  geom_hline(yintercept=0, color="grey50")+
  geom_errorbar(aes(ymin=ddG_foldx_7DAC-ddG_foldx_SD_7DAC, ymax=ddG_foldx_7DAC+ddG_foldx_SD_7DAC), 
                color="grey80", width=.05)+
  geom_errorbar(aes(xmin=ddG_foldx_7DA4-ddG_foldx_SD_7DA4, xmax=ddG_foldx_7DA4+ddG_foldx_SD_7DA4),
                color="grey80", width=.05)+
  stat_cor()+ 
  geom_point(alpha=0.5)+
  theme_classic()+
  theme()

scatterddG_1

ggsave(scatterddG_1, path=path, file="ddG_7DA4vs7DAC.jpg", width=4, height=4)
ggsave(scatterddG_1, path=path, file="ddG_7DA4vs7DAC.pdf", width=4, height=4)

#
scatterddG_2<-ggplot(NS_ddG_RIPK3, aes(x=ddG_foldx_7DA4, y=ddG_foldx_8Z94, label=ID))+
  geom_vline(xintercept=2, linetype="dashed", color="grey80")+
  geom_vline(xintercept=0, color="grey50")+
  geom_hline(yintercept=2, linetype="dashed", color="grey80")+
  geom_hline(yintercept=0, color="grey50")+
  geom_errorbar(aes(ymin=ddG_foldx_8Z94-ddG_foldx_SD_8Z94, ymax=ddG_foldx_8Z94+ddG_foldx_SD_8Z94), 
                color="grey80", width=.05)+
  geom_errorbar(aes(xmin=ddG_foldx_7DA4-ddG_foldx_SD_7DA4, xmax=ddG_foldx_7DA4+ddG_foldx_SD_7DA4),
                color="grey80", width=.05)+
  stat_cor()+ 
  geom_point(alpha=0.5)+
  theme_classic()+
  theme()

scatterddG_2

ggsave(scatterddG_2, path=path, file="ddG_7DA4vs8Z94.jpg", width=4, height=4)
ggsave(scatterddG_2, path=path, file="ddG_7DA4vs8Z94.pdf", width=4, height=4)

#
scatterddG_3<-ggplot(NS_ddG_RIPK3, aes(x=ddG_foldx_7DAC, y=ddG_foldx_8Z94, label=ID))+
  geom_vline(xintercept=2, linetype="dashed", color="grey80")+
  geom_vline(xintercept=0, color="grey50")+
  geom_hline(yintercept=2, linetype="dashed", color="grey80")+
  geom_hline(yintercept=0, color="grey50")+
  geom_errorbar(aes(ymin=ddG_foldx_8Z94-ddG_foldx_SD_8Z94, ymax=ddG_foldx_8Z94+ddG_foldx_SD_8Z94), 
                color="grey80", width=.05)+
  geom_errorbar(aes(xmin=ddG_foldx_7DAC-ddG_foldx_SD_7DAC, xmax=ddG_foldx_7DAC+ddG_foldx_SD_7DAC),
                color="grey80", width=.05)+
  stat_cor()+ 
  geom_point(alpha=0.5)+
  theme_classic()+
  theme()

scatterddG_3

ggsave(scatterddG_3, path=path, file="ddG_7DACvs8Z94.jpg", width=4, height=4)
ggsave(scatterddG_3, path=path, file="ddG_7DACvs8Z94.pdf", width=4, height=4)

################################################################################
# 2: compare the median effect per position between structures (shape and coloring by orientation or ddG_class)
# Let's see which positions have opposing effects in each structure and also opposite orientation

# Group by position and get the median ddG
NS_ddG_RIPK3_<-NS_ddG_RIPK3_ %>% 
  group_by(Pos) %>% 
  mutate(median_ddG_7DA4_pos=median(ddG_foldx_7DA4), 
         median_ddG_7DAC_pos=median(ddG_foldx_7DAC), 
         median_ddG_8Z94_pos=median(ddG_foldx_8Z94) , 
         median_nscore_c=median(nscore_c, na.rm = TRUE))

# Remove duplicated rows
NS_ddG_RIPK3_median<-NS_ddG_RIPK3_[c("Pos", "Residue", "orientation", "median_ddG_7DA4_pos", "median_ddG_7DAC_pos", "median_ddG_8Z94_pos", "median_nscore_c")]
NS_ddG_RIPK3_median<-NS_ddG_RIPK3_median[!duplicated(NS_ddG_RIPK3_median$Pos),]

##
# Label residues with opposite ddG
NS_ddG_RIPK3_median$ddG_class<-"opposite_effect"
NS_ddG_RIPK3_median[NS_ddG_RIPK3_median$median_ddG_7DA4_pos>2 & NS_ddG_RIPK3_median$median_ddG_7DAC_pos>2, "ddG_class"]<-"destabilizing_effect"
NS_ddG_RIPK3_median[NS_ddG_RIPK3_median$median_ddG_7DA4_pos<2 & NS_ddG_RIPK3_median$median_ddG_7DAC_pos<2, "ddG_class"]<-"stabilizing_effect"

# Label them
NS_ddG_RIPK3_median$ddG_class_ID<-""
NS_ddG_RIPK3_median[NS_ddG_RIPK3_median$ddG_class == "opposite_effect", "ddG_class_ID"]<-NS_ddG_RIPK3_median[NS_ddG_RIPK3_median$ddG_class == "opposite_effect", "Residue"]

# Classify the residues based on the comparison between structures
NS_ddG_RIPK3_median$orientation<-"exposed"
NS_ddG_RIPK3_median[NS_ddG_RIPK3_median$Pos %in% c(450, 452, 456, 458, 459, 460, 466, 468), "orientation"]<-"buried"
NS_ddG_RIPK3_median[NS_ddG_RIPK3_median$Pos %in% c(457, 461), "orientation"]<-"glycine"
NS_ddG_RIPK3_median[NS_ddG_RIPK3_median$Pos %in% c(448, 454, 455, 462), "orientation"]<-"flipped"

#

NS_ddG_RIPK3_median_melted<-melt(NS_ddG_RIPK3_median[c("Residue", "median_ddG_7DA4_pos", "median_ddG_7DAC_pos", "median_ddG_8Z94_pos")], id="Residue")
NS_ddG_RIPK3_median_melted[NS_ddG_RIPK3_median_melted$value>15, "value"]<-15

min<-min(NS_ddG_RIPK3_median_melted$value)
max<-15
cols <- c(colorRampPalette(c("#4775d1", "grey95"))((-min/(-min+max)*100)-0.5), colorRampPalette("grey95")(1),
          colorRampPalette(c("grey95", "#b3003b"), bias=1)((max/(-min+max)*100)-0.5))

Residue_order<-unique(NS_ddG_RIPK3_median_melted$Residue)

heatmap_1<-ggplot(NS_ddG_RIPK3_median_melted[!NS_ddG_RIPK3_median_melted$Residue %in% c("P448", "L449"),], 
                  aes(x=factor(Residue, levels=Residue_order), y=factor(variable, levels=c("median_ddG_7DA4_pos", "median_ddG_7DAC_pos", "median_ddG_8Z94_pos"), labels=c("7DA4", "7DAC", "8Z94")), 
                      fill=value))+
  geom_tile(color='white', size=1)+
  scale_fill_gradientn(colours=cols, limits=c(min,max), na.value = "grey60", breaks=c(-1.5, 0, 5, 10, 15), labels=c(-1.5, 0, 5, 10, ">15"))+
  labs(x="", y="")+
  theme_void()+
  theme(legend.title = element_text(size=14),
        legend.text = element_text(size=14), 
        axis.title = element_text(size = 16),
        axis.text = element_text(size=14))
heatmap_1

ggsave(heatmap_1, path=path, file="heatmap_median_ddG.jpg", width=15, height=3)
ggsave(heatmap_1, path=path, file="heatmap_median_ddG.pdf", width=15, height=3)

#
scatterddG_4<-ggplot(NS_ddG_RIPK3_median, 
                     aes(x=median_ddG_7DA4_pos, y=median_ddG_7DAC_pos, label=Residue))+
  geom_vline(xintercept=2, linetype="dashed", color="grey80")+
  geom_vline(xintercept=0, color="grey50")+
  geom_hline(yintercept=2, linetype="dashed", color="grey80")+
  geom_hline(yintercept=0, color="grey50")+
  geom_point(aes(color=ddG_class, shape=factor(orientation, levels=c("buried", "exposed", "glycine", "flipped"))), 
             size=4)+
  scale_color_manual(values=c("red3", "green4", "black"))+
  scale_x_break(c(10, 20))+
  #ylim(-2, 15)+
  geom_text_repel(min.segment.length = 0.1)+
  stat_cor()+ 
  labs(shape="orientation")+
  theme_classic()+
  theme()

scatterddG_4

ggsave(scatterddG_4, path=path, file="ddG_7DA4vs7DAC_residues.jpg", width=6, height=4)
ggsave(scatterddG_4, path=path, file="ddG_7DA4vs7DAC_residues.pdf", width=6, height=4)

#
scatterddG_4_<-ggplot(NS_ddG_RIPK3_median[!NS_ddG_RIPK3_median$Residue %in% c("G457", "G461"),], 
                      aes(x=median_ddG_7DA4_pos, y=median_ddG_7DAC_pos, label=Residue))+
  geom_vline(xintercept=2, linetype="dashed", color="grey80")+
  geom_vline(xintercept=0, color="grey50")+
  geom_hline(yintercept=2, linetype="dashed", color="grey80")+
  geom_hline(yintercept=0, color="grey50")+
  geom_point(aes(color=factor(orientation, levels=c("buried", "exposed", "flipped"))), size=4)+
  scale_color_manual(values=c("grey80", "grey80", "green4"))+
  labs(color="")+
  #xlim(-2, 15)+
  #ylim(-2, 15)+
  geom_text_repel(min.segment.length = 0.1)+
  stat_cor()+ 
  labs(shape="orientation")+
  theme_classic()+
  theme()

scatterddG_4_

ggsave(scatterddG_4_, path=path, file="ddG_7DA4vs7DAC_residues_noG.jpg", width=5, height=4)
ggsave(scatterddG_4_, path=path, file="ddG_7DA4vs7DAC_residues_noG.pdf", width=5, height=4)

##
# Label residues with opposite ddG
NS_ddG_RIPK3_median$ddG_class<-"opposite_effect"
NS_ddG_RIPK3_median[NS_ddG_RIPK3_median$median_ddG_7DA4_pos>2 & NS_ddG_RIPK3_median$median_ddG_8Z94_pos>2, "ddG_class"]<-"destabilizing_effect"
NS_ddG_RIPK3_median[NS_ddG_RIPK3_median$median_ddG_7DA4_pos<2 & NS_ddG_RIPK3_median$median_ddG_8Z94_pos<2, "ddG_class"]<-"stabilizing_effect"

# Label them
NS_ddG_RIPK3_median$ddG_class_ID<-""
NS_ddG_RIPK3_median[NS_ddG_RIPK3_median$ddG_class == "opposite_effect", "ddG_class_ID"]<-NS_ddG_RIPK3_median[NS_ddG_RIPK3_median$ddG_class == "opposite_effect", "Residue"]

# Classify the residues based on the comparison between structures
NS_ddG_RIPK3_median$orientation<-"exposed"
NS_ddG_RIPK3_median[NS_ddG_RIPK3_median$Pos %in% c(450, 452, 454, 456, 458, 459, 460, 466, 468), "orientation"]<-"buried"
NS_ddG_RIPK3_median[NS_ddG_RIPK3_median$Pos %in% c(457, 461), "orientation"]<-"glycine"

scatterddG_5<-ggplot(NS_ddG_RIPK3_median, 
                     aes(x=median_ddG_7DA4_pos, y=median_ddG_8Z94_pos, label=Residue))+
  geom_vline(xintercept=2, linetype="dashed", color="grey80")+
  geom_vline(xintercept=0, color="grey50")+
  geom_hline(yintercept=2, linetype="dashed", color="grey80")+
  geom_hline(yintercept=0, color="grey50")+
  geom_point(aes(color=ddG_class, shape=factor(orientation, levels=c("buried", "exposed", "glycine", "flipped"))), size=2.5)+
  scale_color_manual(values=c("red3", "green4", "black"))+
  #scale_y_break(c(10, 30))+
  scale_x_break(c(10, 20))+
  #xlim(-2, 24)+
  #ylim(-2, 15)+
  geom_text_repel(min.segment.length = 0.1)+
  stat_cor()+ 
  labs(shape="orientation")+
  theme_classic()+
  theme()

scatterddG_5

ggsave(scatterddG_5, path=path, file="ddG_7DA4vs8Z94_residues.jpg", width=6, height=4)
ggsave(scatterddG_5, path=path, file="ddG_7DA4vs8Z94_residues.pdf", width=6, height=4)

#
#
scatterddG_5_<-ggplot(NS_ddG_RIPK3_median[!NS_ddG_RIPK3_median$Residue %in% c("G457", "G461"),], 
                      aes(x=median_ddG_7DA4_pos, y=median_ddG_8Z94_pos, label=Residue))+
  geom_vline(xintercept=2, linetype="dashed", color="grey80")+
  geom_vline(xintercept=0, color="grey50")+
  geom_hline(yintercept=2, linetype="dashed", color="grey80")+
  geom_hline(yintercept=0, color="grey50")+
  geom_point(aes(color=factor(orientation, levels=c("buried", "exposed", "flipped"))), size=4)+
  scale_color_manual(values=c("grey80", "grey80", "green4"))+
  labs(color="")+
  #xlim(-2, 15)+
  #ylim(-2, 15)+
  geom_text_repel(min.segment.length = 0.1)+
  stat_cor()+ 
  labs(shape="orientation")+
  theme_classic()+
  theme()

scatterddG_5_

ggsave(scatterddG_5_, path=path, file="ddG_7DA4vs8Z94_residues_noG.jpg", width=5, height=4)
ggsave(scatterddG_5_, path=path, file="ddG_7DA4vs8Z94_residues_noG.pdf", width=5, height=4)

##
# Label residues with opposite ddG
NS_ddG_RIPK3_median$ddG_class<-"opposite_effect"
NS_ddG_RIPK3_median[NS_ddG_RIPK3_median$median_ddG_7DAC_pos>2 & NS_ddG_RIPK3_median$median_ddG_8Z94_pos>2, "ddG_class"]<-"destabilizing_effect"
NS_ddG_RIPK3_median[NS_ddG_RIPK3_median$median_ddG_7DAC_pos<2 & NS_ddG_RIPK3_median$median_ddG_8Z94_pos<2, "ddG_class"]<-"stabilizing_effect"

# Label them
NS_ddG_RIPK3_median$ddG_class_ID<-""
NS_ddG_RIPK3_median[NS_ddG_RIPK3_median$ddG_class == "opposite_effect", "ddG_class_ID"]<-NS_ddG_RIPK3_median[NS_ddG_RIPK3_median$ddG_class == "opposite_effect", "Residue"]

# Classify the residues based on the comparison between structures
NS_ddG_RIPK3_median$orientation<-"exposed"
NS_ddG_RIPK3_median[NS_ddG_RIPK3_median$Pos %in% c(450, 452, 456, 458, 459, 460, 466, 468), "orientation"]<-"buried"
NS_ddG_RIPK3_median[NS_ddG_RIPK3_median$Pos %in% c(457, 461), "orientation"]<-"glycine"
NS_ddG_RIPK3_median[NS_ddG_RIPK3_median$Pos %in% c(448, 454, 455, 462), "orientation"]<-"flipped"

scatterddG_6<-ggplot(NS_ddG_RIPK3_median, 
                     aes(x=median_ddG_7DAC_pos, y=median_ddG_8Z94_pos, label=Residue))+
  geom_vline(xintercept=2, linetype="dashed", color="grey80")+
  geom_vline(xintercept=0, color="grey50")+
  geom_hline(yintercept=2, linetype="dashed", color="grey80")+
  geom_hline(yintercept=0, color="grey50")+
  geom_point(aes(color=ddG_class, shape=factor(orientation, levels=c("buried", "exposed", "glycine", "flipped"))), size=2.5)+
  scale_color_manual(values=c("red3", "green4", "black"))+
  scale_y_break(c(10, 25))+
  #xlim(-2, 24)+
  #ylim(-2, 15)+
  geom_text_repel(min.segment.length = 0.1)+
  stat_cor()+ 
  labs(shape="orientation")+
  theme_classic()+
  theme()

scatterddG_6

ggsave(scatterddG_6, path=path, file="ddG_7DACvs8Z94_residues.jpg", width=6, height=4)
ggsave(scatterddG_6, path=path, file="ddG_7DACvs8Z94_residues.pdf", width=6, height=4)

#
scatterddG_6_<-ggplot(NS_ddG_RIPK3_median[!NS_ddG_RIPK3_median$Residue %in% c("G457", "G461"),], 
                      aes(x=median_ddG_7DAC_pos, y=median_ddG_8Z94_pos, label=Residue))+
  geom_vline(xintercept=2, linetype="dashed", color="grey80")+
  geom_vline(xintercept=0, color="grey50")+
  geom_hline(yintercept=2, linetype="dashed", color="grey80")+
  geom_hline(yintercept=0, color="grey50")+
  geom_point(aes(color=factor(orientation, levels=c("buried", "exposed", "flipped"))), size=4)+
  scale_color_manual(values=c("grey80", "grey80", "green4"))+
  labs(color="")+
  #xlim(-2, 15)+
  #ylim(-2, 15)+
  geom_text_repel(min.segment.length = 0.1)+
  stat_cor()+ 
  labs(shape="orientation")+
  theme_classic()+
  theme()

scatterddG_6_

ggsave(scatterddG_6_, path=path, file="ddG_7DACvs8Z94_residues_noG.jpg", width=5, height=4)
ggsave(scatterddG_6_, path=path, file="ddG_7DACvs8Z94_residues_noG.pdf", width=5, height=4)


################################################################################
# 3: correlation between ddG and NS for the different structures
scatterddG_7<-ggplot(NS_ddG_RIPK3, 
                     aes(x=nscore_c, y=ddG_foldx_7DA4))+
  geom_hline(yintercept=0, color="grey50")+
  geom_vline(xintercept=0, color="grey50")+
  stat_cor()+ 
  geom_point(alpha=0.5)+
  #xlim(-9, 2)+
  theme_classic()+
  theme()

scatterddG_7

ggsave(scatterddG_7, path=path, file="ddG_7DA4_nscore.jpg", width=4, height=4)
ggsave(scatterddG_7, path=path, file="ddG_7DA4_nscore.pdf", width=4, height=4)

#
scatterddG_8<-ggplot(NS_ddG_RIPK3, 
                     aes(x=nscore_c, y=ddG_foldx_7DAC))+
  geom_hline(yintercept=0, color="grey50")+
  geom_vline(xintercept=0, color="grey50")+
  stat_cor()+ 
  geom_point(alpha=0.5)+
  #xlim(-9, 2)+
  theme_classic()+
  theme()

scatterddG_8

ggsave(scatterddG_8, path=path, file="ddG_7DAC_nscore.jpg", width=4, height=4)
ggsave(scatterddG_8, path=path, file="ddG_7DAC_nscore.pdf", width=4, height=4)

#
scatterddG_9<-ggplot(NS_ddG_RIPK3, 
                     aes(x=nscore_c, y=ddG_foldx_8Z94))+
  geom_hline(yintercept=0, color="grey50")+
  geom_vline(xintercept=0, color="grey50")+
  stat_cor()+ 
  geom_point(alpha=0.5)+
  #xlim(-9, 2)+
  theme_classic()+
  theme()

scatterddG_9

ggsave(scatterddG_9, path=path, file="ddG_8Z94_nscore.jpg", width=4, height=4)
ggsave(scatterddG_9, path=path, file="ddG_8Z94_nscore.pdf", width=4, height=4)

################################################################################
# 4: correlation between ddG and NS for the different structures facet by structural elements

NS_ddG_RIPK3_<-drop_na(NS_ddG_RIPK3_)

scatterddG_10<-ggplot(NS_ddG_RIPK3_[NS_ddG_RIPK3_$structure != "Disordered-1" & NS_ddG_RIPK3_$structure != "Disordered-2",], 
                      aes(x=nscore_c, y=ddG_foldx_7DAC))+
  geom_hline(yintercept=0, color="grey50")+
  geom_vline(xintercept=0, color="grey50")+
  stat_cor()+ 
  geom_point(alpha=0.5)+
  facet_wrap(~factor(structure, levels=c("Disordered-1", "β-strand-1", "Loop-1", "β-strand-2", "Loop-2", "β-strand-3", "Disordered-2")))+
  theme_bw()+
  theme()

scatterddG_10

ggsave(scatterddG_10, path=path, file="ddG_7DAC_nscore_strtype.jpg", width=6, height=4)
ggsave(scatterddG_10, path=path, file="ddG_7DAC_nscore_strtype.pdf", width=6, height=4)

# Correlation per structure element

corr_vector=c()
for(i in NS_ddG_RIPK3_$structure){
  correlation<-cor.test(NS_ddG_RIPK3_[NS_ddG_RIPK3_$structure==i,]$nscore_c,
                        NS_ddG_RIPK3_[NS_ddG_RIPK3_$structure==i,]$ddG_foldx_7DAC, use="complete.obs", method="pearson")
  corr<-correlation$estimate
  p_value<-correlation$p.value
  corr_vector=c(corr_vector, i, corr, p_value)
}

corr_text <- data.frame(matrix(corr_vector, ncol=3, byrow = T))
colnames(corr_text)<-c("structure", "corr","pvalue")
for(i in c(2:3)){corr_text[[i]]<-as.numeric(as.character(corr_text[[i]]))}
corr_text<-distinct(corr_text, structure, .keep_all = TRUE)

corr_text[is.na(corr_text)]<-1

corr_text$significance_pos<-""
corr_text[(corr_text$pvalue<0.05) & (corr_text$corr>0), "significance_pos"]<-"*"
corr_text[(corr_text$pvalue<0.01) & (corr_text$corr>0), "significance_pos"]<-"**"
corr_text[(corr_text$pvalue<0.001) & (corr_text$corr>0), "significance_pos"]<-"***"

corr_text$significance_neg<-""
corr_text[(corr_text$pvalue<0.05) & (corr_text$corr<0), "significance_neg"]<-"*"
corr_text[(corr_text$pvalue<0.01) & (corr_text$corr<0), "significance_neg"]<-"**"
corr_text[(corr_text$pvalue<0.001) & (corr_text$corr<0), "significance_neg"]<-"***"


p_bars_1<-ggplot(corr_text, aes(x=factor(structure, levels=c("Disordered-1", "β-strand-1", "Loop-1", "β-strand-2", 
                                                             "Loop-2", "β-strand-3", "Disordered-2")), 
                                y=corr))+
  geom_hline(yintercept = 0, color="black", linewidth=.5)+
  geom_bar(stat="identity", color="black", fill="grey", width=.5)+
  geom_text(data=corr_text, aes(label=significance_pos), vjust=0.4, size=8, colour="black")+
  geom_text(data=corr_text, aes(label=significance_neg), vjust=1, size=8, colour="black")+
  theme_void()+
  labs(x="", y="Pearson Correlation")+
  ylim(-1, 1)+
  theme_classic()

p_bars_1

ggsave(p_bars_1, file="ddG_7DAC_nscore_strtype_bars.jpg", width = 4, height = 3, path=path)
ggsave(p_bars_1, file="ddG_7DAC_nscore_strtype_bars.pdf", width = 4, height = 3, path=path)

#
scatterddG_11<-ggplot(NS_ddG_RIPK3_[NS_ddG_RIPK3_$structure != "Disordered-1" & NS_ddG_RIPK3_$structure != "Disordered-2",],
                      aes(x=nscore_c, y=ddG_foldx_7DA4))+
  geom_hline(yintercept=0, color="grey50")+
  geom_vline(xintercept=0, color="grey50")+
  stat_cor()+ 
  geom_point(alpha=0.5)+
  facet_wrap(~factor(structure, levels=c("Disordered-1", "β-strand-1", "Loop-1", "β-strand-2", "Loop-2", "β-strand-3", "Disordered-2")))+
  theme_bw()+
  theme()

scatterddG_11

ggsave(scatterddG_11, path=path, file="ddG_7DA4_nscore_strtype.jpg", width=6, height=4)
ggsave(scatterddG_11, path=path, file="ddG_7DA4_nscore_strtype.pdf", width=6, height=4)

# Correlation per structure element

corr_vector=c()
for(i in NS_ddG_RIPK3_$structure){
  correlation<-cor.test(NS_ddG_RIPK3_[NS_ddG_RIPK3_$structure==i,]$nscore_c,
                        NS_ddG_RIPK3_[NS_ddG_RIPK3_$structure==i,]$ddG_foldx_7DA4, use="complete.obs", method="pearson")
  corr<-correlation$estimate
  p_value<-correlation$p.value
  corr_vector=c(corr_vector, i, corr, p_value)
}

corr_text <- data.frame(matrix(corr_vector, ncol=3, byrow = T))
colnames(corr_text)<-c("structure", "corr","pvalue")
for(i in c(2:3)){corr_text[[i]]<-as.numeric(as.character(corr_text[[i]]))}
corr_text<-distinct(corr_text, structure, .keep_all = TRUE)

corr_text[is.na(corr_text)]<-1

corr_text$significance_pos<-""
corr_text[(corr_text$pvalue<0.05) & (corr_text$corr>0), "significance_pos"]<-"*"
corr_text[(corr_text$pvalue<0.01) & (corr_text$corr>0), "significance_pos"]<-"**"
corr_text[(corr_text$pvalue<0.001) & (corr_text$corr>0), "significance_pos"]<-"***"

corr_text$significance_neg<-""
corr_text[(corr_text$pvalue<0.05) & (corr_text$corr<0), "significance_neg"]<-"*"
corr_text[(corr_text$pvalue<0.01) & (corr_text$corr<0), "significance_neg"]<-"**"
corr_text[(corr_text$pvalue<0.001) & (corr_text$corr<0), "significance_neg"]<-"***"


p_bars_2<-ggplot(corr_text, aes(x=factor(structure, levels=c("Disordered-1", "β-strand-1", "Loop-1", "β-strand-2", 
                                                             "Loop-2", "β-strand-3", "Disordered-2")), 
                                y=corr))+
  geom_hline(yintercept = 0, color="black", linewidth=.5)+
  geom_bar(stat="identity", color="black", fill="grey", width=.5)+
  geom_text(data=corr_text, aes(label=significance_pos), vjust=0.4, size=8, colour="black")+
  geom_text(data=corr_text, aes(label=significance_neg), vjust=1, size=8, colour="black")+
  theme_void()+
  labs(x="", y="Pearson Correlation")+
  ylim(-1, 1)+
  theme_classic()

p_bars_2

ggsave(p_bars_2, file="ddG_7DA4_nscore_strtype_bars.jpg", width = 4, height = 3, path=path)
ggsave(p_bars_2, file="ddG_7DA4_nscore_strtype_bars.pdf", width = 4, height = 3, path=path)

#
scatterddG_12<-ggplot(NS_ddG_RIPK3_[NS_ddG_RIPK3_$structure != "Disordered-1" & NS_ddG_RIPK3_$structure != "Disordered-2",],
                      aes(x=nscore_c, y=ddG_foldx_8Z94))+
  geom_hline(yintercept=0, color="grey50")+
  geom_vline(xintercept=0, color="grey50")+
  stat_cor()+ 
  geom_point(alpha=0.5)+
  facet_wrap(~factor(structure, levels=c("Disordered-1", "β-strand-1", "Loop-1", "β-strand-2", "Loop-2", "β-strand-3", "Disordered-2")))+
  theme_bw()+
  theme()

scatterddG_12

ggsave(scatterddG_12, path=path, file="ddG_8Z94_nscore_strtype.jpg", width=6, height=4)
ggsave(scatterddG_12, path=path, file="ddG_8Z94_nscore_strtype.pdf", width=6, height=4)

# Correlation per structure element

corr_vector=c()
for(i in NS_ddG_RIPK3_$structure){
  correlation<-cor.test(NS_ddG_RIPK3_[NS_ddG_RIPK3_$structure==i,]$nscore_c,
                        NS_ddG_RIPK3_[NS_ddG_RIPK3_$structure==i,]$ddG_foldx_8Z94, use="complete.obs", method="pearson")
  corr<-correlation$estimate
  p_value<-correlation$p.value
  corr_vector=c(corr_vector, i, corr, p_value)
}

corr_text <- data.frame(matrix(corr_vector, ncol=3, byrow = T))
colnames(corr_text)<-c("structure", "corr","pvalue")
for(i in c(2:3)){corr_text[[i]]<-as.numeric(as.character(corr_text[[i]]))}
corr_text<-distinct(corr_text, structure, .keep_all = TRUE)

corr_text[is.na(corr_text)]<-1

corr_text$significance_pos<-""
corr_text[(corr_text$pvalue<0.05) & (corr_text$corr>0), "significance_pos"]<-"*"
corr_text[(corr_text$pvalue<0.01) & (corr_text$corr>0), "significance_pos"]<-"**"
corr_text[(corr_text$pvalue<0.001) & (corr_text$corr>0), "significance_pos"]<-"***"

corr_text$significance_neg<-""
corr_text[(corr_text$pvalue<0.05) & (corr_text$corr<0), "significance_neg"]<-"*"
corr_text[(corr_text$pvalue<0.01) & (corr_text$corr<0), "significance_neg"]<-"**"
corr_text[(corr_text$pvalue<0.001) & (corr_text$corr<0), "significance_neg"]<-"***"


p_bars_3<-ggplot(corr_text, aes(x=factor(structure, levels=c("Disordered-1", "β-strand-1", "Loop-1", "β-strand-2", 
                                                             "Loop-2", "β-strand-3", "Disordered-2")), 
                                y=corr))+
  geom_hline(yintercept = 0, color="black", linewidth=.5)+
  geom_bar(stat="identity", color="black", fill="grey", width=.5)+
  geom_text(data=corr_text, aes(label=significance_pos), vjust=0.4, size=8, colour="black")+
  geom_text(data=corr_text, aes(label=significance_neg), vjust=1, size=8, colour="black")+
  theme_void()+
  labs(x="", y="Pearson Correlation")+
  ylim(-1, 1)+
  theme_classic()

p_bars_3

ggsave(p_bars_3, file="ddG_8Z94_nscore_strtype_bars.jpg", width = 4, height = 3, path=path)
ggsave(p_bars_3, file="ddG_8Z94_nscore_strtype_bars.pdf", width = 4, height = 3, path=path)
