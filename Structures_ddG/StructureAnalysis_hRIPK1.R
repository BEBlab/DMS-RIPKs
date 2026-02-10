library(tidyverse)
library(ggpubr)
library(ggsignif)
require("ggrepel")
library(ggbreak) 
library(reshape2)

# your folder path
setwd("")

dir.create("02_StructuralComparison_hRIPK1")
path="02_StructuralComparison_hRIPK1"

### hRIPK1
# Nucleation Score
load("nscore_df_hRIPK1.RData")

### ddG
#PDB 9HR9
ddG_foldx_9HR9<-read.csv("ddG_foldx_9HR9_3monomers_compact.csv")
ddG_foldx_9HR9<-ddG_foldx_9HR9[c("ID", "WT_AA", "Pos", "Mut", "ddG_foldx", "ddG_foldx_SD")] 
ddG_foldx_9HR9$ID<-paste0(ddG_foldx_9HR9$WT_AA, "-", ddG_foldx_9HR9$Pos, "-", ddG_foldx_9HR9$Mut) 
ddG_foldx_9HR9<-rename(ddG_foldx_9HR9, ddG_foldx_9HR9 = ddG_foldx, ddG_foldx_SD_9HR9 = ddG_foldx_SD)
ddG_foldx_9HR9$Residue<-paste0(ddG_foldx_9HR9$WT_AA, ddG_foldx_9HR9$Pos)

# Residues classification based on side chain orientation
ddG_foldx_9HR9$orientation_9HR9<-"exposed"
ddG_foldx_9HR9[ddG_foldx_9HR9$Pos %in% c(529, 531, 533, 536, 539, 540, 541, 545, 547, 549), "orientation_9HR9"]<-"buried"
ddG_foldx_9HR9[ddG_foldx_9HR9$Pos %in% c(538, 542), "orientation_9HR9"]<-"glycine"


#PDB 9HR6
ddG_foldx_9HR6<-read.csv("ddG_foldx_9HR6_3monomers_compact.csv")
ddG_foldx_9HR6$Pos<-ddG_foldx_9HR6$Pos+495 #correct the numbering
ddG_foldx_9HR6<-ddG_foldx_9HR6[c("ID", "WT_AA", "Pos", "Mut", "ddG_foldx", "ddG_foldx_SD")]
ddG_foldx_9HR6$ID<-paste0(ddG_foldx_9HR6$WT_AA, "-", ddG_foldx_9HR6$Pos, "-", ddG_foldx_9HR6$Mut)
ddG_foldx_9HR6<-rename(ddG_foldx_9HR6, ddG_foldx_9HR6 = ddG_foldx, ddG_foldx_SD_9HR6 = ddG_foldx_SD)
ddG_foldx_9HR6$Residue<-paste0(ddG_foldx_9HR6$WT_AA, ddG_foldx_9HR6$Pos)

# Residues classification based on side chain orientation
ddG_foldx_9HR6$orientation_9HR6<-"exposed"
ddG_foldx_9HR6[ddG_foldx_9HR6$Pos %in% c(529, 531, 533, 536, 539, 540, 541, 545, 547, 549), "orientation_9HR6"]<-"buried"
ddG_foldx_9HR6[ddG_foldx_9HR6$Pos %in% c(538, 542), "orientation_9HR6"]<-"glycine"


#PDB 8Z93
ddG_foldx_8Z93<-read.csv("ddG_foldx_8Z93_3monomers_compact.csv")
ddG_foldx_8Z93<-ddG_foldx_8Z93[c("ID", "WT_AA", "Pos", "Mut", "ddG_foldx", "ddG_foldx_SD")]
ddG_foldx_8Z93$ID<-paste0(ddG_foldx_8Z93$WT_AA, "-", ddG_foldx_8Z93$Pos, "-", ddG_foldx_8Z93$Mut)
ddG_foldx_8Z93<-rename(ddG_foldx_8Z93, ddG_foldx_8Z93 = ddG_foldx, ddG_foldx_SD_8Z93 = ddG_foldx_SD)
ddG_foldx_8Z93$Residue<-paste0(ddG_foldx_8Z93$WT_AA, ddG_foldx_8Z93$Pos)

# Residues classification based on side chain orientation
ddG_foldx_8Z93$orientation_8Z93<-"exposed"
ddG_foldx_8Z93[ddG_foldx_8Z93$Pos %in% c(531, 533, 535, 539, 540, 541, 543, 545, 547, 549), "orientation_8Z93"]<-"buried"
ddG_foldx_8Z93[ddG_foldx_8Z93$Pos %in% c(538, 542), "orientation_8Z93"]<-"glycine"

###
NS_ddG_ripk1<-full_join(ddG_foldx_9HR6[c("ID", "Pos", "Residue", "ddG_foldx_9HR6", "ddG_foldx_SD_9HR6", "orientation_9HR6")], 
                        ddG_foldx_9HR9[c("ID", "ddG_foldx_9HR9", "ddG_foldx_SD_9HR9", "orientation_9HR9")], by="ID")

NS_ddG_ripk1<-full_join(NS_ddG_ripk1, ddG_foldx_8Z93[c("ID", "ddG_foldx_8Z93", "ddG_foldx_SD_8Z93", "orientation_8Z93")], by="ID")

NS_ddG_ripk1<-full_join(NS_ddG_ripk1, singles_stops[c("ID", "nscore_c", "sigma")], by="ID")

# Classify the residues based on the comparison between structures
NS_ddG_ripk1$orientation<-"exposed"
NS_ddG_ripk1[NS_ddG_ripk1$Pos %in% c(531, 533, 539, 540, 541, 543, 545, 547, 549), "orientation"]<-"buried"
NS_ddG_ripk1[NS_ddG_ripk1$Pos %in% c(538, 542), "orientation"]<-"glycine"
NS_ddG_ripk1[NS_ddG_ripk1$Pos %in% c(535, 536, 543), "orientation"]<-"flipped"

# drop nas
NS_ddG_ripk1_<- NS_ddG_ripk1 %>% drop_na(ddG_foldx_9HR6, ddG_foldx_9HR9, ddG_foldx_8Z93)

# Classify the mutations based on the ddG
NS_ddG_ripk1_$ddG_class<-"opposite_effect"
NS_ddG_ripk1_[NS_ddG_ripk1_$ddG_foldx_9HR6>2 & NS_ddG_ripk1_$ddG_foldx_9HR9>2 & NS_ddG_ripk1_$ddG_foldx_8Z93>2, "ddG_class"]<-"destabilizing_effect"
NS_ddG_ripk1_[NS_ddG_ripk1_$ddG_foldx_9HR6<2 & NS_ddG_ripk1_$ddG_foldx_9HR9<2 & NS_ddG_ripk1_$ddG_foldx_8Z93<2, "ddG_class"]<-"stabilizing_effect"

# label the mutations with opposite effects (at least in one structure)
NS_ddG_ripk1_$ddG_class_ID<-""
NS_ddG_ripk1_[NS_ddG_ripk1_$ddG_class == "opposite_effect", "ddG_class_ID"]<-NS_ddG_ripk1_[NS_ddG_ripk1_$ddG_class == "opposite_effect", "ID"]

# Map the structural elements
NS_ddG_ripk1_$structure<-"Disordered-1"
NS_ddG_ripk1_[NS_ddG_ripk1_$Pos>=531, "structure"]<-"β-strand-1"
NS_ddG_ripk1_[NS_ddG_ripk1_$Pos>=534, "structure"]<-"Loop-1"
NS_ddG_ripk1_[NS_ddG_ripk1_$Pos>=539, "structure"]<-"β-strand-2"
NS_ddG_ripk1_[NS_ddG_ripk1_$Pos>=542, "structure"]<-"Loop-2"
NS_ddG_ripk1_[NS_ddG_ripk1_$Pos>=547, "structure"]<-"β-strand-3"
NS_ddG_ripk1_[NS_ddG_ripk1_$Pos>=550, "structure"]<-"Disordered-2"

NS_ddG_ripk1_$strtype<-"Disordered"
NS_ddG_ripk1_[NS_ddG_ripk1_$Pos>=531, "strtype"]<-"β-strand"
NS_ddG_ripk1_[NS_ddG_ripk1_$Pos>=534, "strtype"]<-"Loop"
NS_ddG_ripk1_[NS_ddG_ripk1_$Pos>=539, "strtype"]<-"β-strand"
NS_ddG_ripk1_[NS_ddG_ripk1_$Pos>=542, "strtype"]<-"Loop"
NS_ddG_ripk1_[NS_ddG_ripk1_$Pos>=547, "strtype"]<-"β-strand"
NS_ddG_ripk1_[NS_ddG_ripk1_$Pos>=550, "strtype"]<-"Disordered"

NS_ddG_ripk1_$interface<-"surface"
NS_ddG_ripk1_[NS_ddG_ripk1_$Pos %in% c(531, 533, 539, 541), "interface"]<-"1st-Interface"
NS_ddG_ripk1_[NS_ddG_ripk1_$Pos %in% c(537, 540, 547, 549), "interface"]<-"2nd-Interface"

################################################################################
# 1: compare the ddG between the different structures
scatterddG_1<-ggplot(NS_ddG_ripk1, aes(x=ddG_foldx_9HR6, y=ddG_foldx_9HR9, label=ID))+
  geom_vline(xintercept=2, linetype="dashed", color="grey80")+
  geom_vline(xintercept=0, color="grey50")+
  geom_hline(yintercept=2, linetype="dashed", color="grey80")+
  geom_hline(yintercept=0, color="grey50")+
  geom_errorbar(aes(ymin=ddG_foldx_9HR9-ddG_foldx_SD_9HR9, ymax=ddG_foldx_9HR9+ddG_foldx_SD_9HR9), 
                color="grey80", width=.05)+
  geom_errorbar(aes(xmin=ddG_foldx_9HR6-ddG_foldx_SD_9HR6, xmax=ddG_foldx_9HR6+ddG_foldx_SD_9HR6),
                color="grey80", width=.05)+
  stat_cor()+ 
  geom_point(alpha=0.5)+
  theme_classic()+
  theme()

scatterddG_1

ggsave(scatterddG_1, path=path, file="ddG_9HR6vs9HR9.jpg", width=4, height=4)
ggsave(scatterddG_1, path=path, file="ddG_9HR6vs9HR9.pdf", width=4, height=4)

#
scatterddG_2<-ggplot(NS_ddG_ripk1, aes(x=ddG_foldx_9HR6, y=ddG_foldx_8Z93, label=ID))+
  geom_vline(xintercept=2, linetype="dashed", color="grey80")+
  geom_vline(xintercept=0, color="grey50")+
  geom_hline(yintercept=2, linetype="dashed", color="grey80")+
  geom_hline(yintercept=0, color="grey50")+
  geom_errorbar(aes(ymin=ddG_foldx_8Z93-ddG_foldx_SD_8Z93, ymax=ddG_foldx_8Z93+ddG_foldx_SD_8Z93), 
                color="grey80", width=.05)+
  geom_errorbar(aes(xmin=ddG_foldx_9HR6-ddG_foldx_SD_9HR6, xmax=ddG_foldx_9HR6+ddG_foldx_SD_9HR6),
                color="grey80", width=.05)+
  stat_cor()+ 
  geom_point(alpha=0.5)+
  theme_classic()+
  theme()

scatterddG_2

ggsave(scatterddG_2, path=path, file="ddG_9HR6vs8Z93.jpg", width=4, height=4)
ggsave(scatterddG_2, path=path, file="ddG_9HR6vs8Z93.pdf", width=4, height=4)

#
scatterddG_3<-ggplot(NS_ddG_ripk1, aes(x=ddG_foldx_9HR9, y=ddG_foldx_8Z93, label=ID))+
  geom_vline(xintercept=2, linetype="dashed", color="grey80")+
  geom_vline(xintercept=0, color="grey50")+
  geom_hline(yintercept=2, linetype="dashed", color="grey80")+
  geom_hline(yintercept=0, color="grey50")+
  geom_errorbar(aes(ymin=ddG_foldx_8Z93-ddG_foldx_SD_8Z93, ymax=ddG_foldx_8Z93+ddG_foldx_SD_8Z93), 
                color="grey80", width=.05)+
  geom_errorbar(aes(xmin=ddG_foldx_9HR9-ddG_foldx_SD_9HR9, xmax=ddG_foldx_9HR9+ddG_foldx_SD_9HR9),
                color="grey80", width=.05)+
  stat_cor()+ 
  geom_point(alpha=0.5)+
  theme_classic()+
  theme()

scatterddG_3

ggsave(scatterddG_3, path=path, file="ddG_9HR9vs8Z93.jpg", width=4, height=4)
ggsave(scatterddG_3, path=path, file="ddG_9HR9vs8Z93.pdf", width=4, height=4)

################################################################################
# 2: compare the median effect per position between structures (shape and coloring by orientation or ddG_class)
# Let's see which positions have opposing effects in each structure and also opposite orientation

# Group by position and get the median ddG
NS_ddG_ripk1_<-NS_ddG_ripk1_ %>% 
  group_by(Pos) %>% 
  mutate(median_ddG_9HR6_pos=median(ddG_foldx_9HR6), 
         median_ddG_9HR9_pos=median(ddG_foldx_9HR9), 
         median_ddG_8Z93_pos=median(ddG_foldx_8Z93) , 
         median_nscore_c=median(nscore_c, na.rm = TRUE))

# Remove duplicated rows
NS_ddG_ripk1_median<-NS_ddG_ripk1_[c("Pos", "Residue", "orientation", "median_ddG_9HR6_pos", "median_ddG_9HR9_pos", "median_ddG_8Z93_pos", "median_nscore_c")]
NS_ddG_ripk1_median<-NS_ddG_ripk1_median[!duplicated(NS_ddG_ripk1_median$Pos),]

##
# Label residues with opposite ddG
NS_ddG_ripk1_median$ddG_class<-"opposite_effect"
NS_ddG_ripk1_median[NS_ddG_ripk1_median$median_ddG_9HR6_pos>2 & NS_ddG_ripk1_median$median_ddG_9HR9_pos>2, "ddG_class"]<-"destabilizing_effect"
NS_ddG_ripk1_median[NS_ddG_ripk1_median$median_ddG_9HR6_pos<2 & NS_ddG_ripk1_median$median_ddG_9HR9_pos<2, "ddG_class"]<-"stabilizing_effect"

# Label them
NS_ddG_ripk1_median$ddG_class_ID<-""
NS_ddG_ripk1_median[NS_ddG_ripk1_median$ddG_class == "opposite_effect", "ddG_class_ID"]<-NS_ddG_ripk1_median[NS_ddG_ripk1_median$ddG_class == "opposite_effect", "Residue"]

# Classify the residues based on the comparison between structures
NS_ddG_ripk1_median$orientation<-"exposed"
NS_ddG_ripk1_median[NS_ddG_ripk1_median$Pos %in% c(531, 533, 536, 539, 540, 541, 543, 545, 547, 549), "orientation"]<-"buried"
NS_ddG_ripk1_median[NS_ddG_ripk1_median$Pos %in% c(538, 542), "orientation"]<-"glycine"

#

NS_ddG_ripk1_median_melted<-melt(NS_ddG_ripk1_median[c("Residue", "median_ddG_9HR6_pos", "median_ddG_9HR9_pos", "median_ddG_8Z93_pos")], id="Residue")
NS_ddG_ripk1_median_melted[NS_ddG_ripk1_median_melted$value>15, "value"]<-15

min<-min(NS_ddG_ripk1_median_melted$value)
max<-15
cols <- c(colorRampPalette(c("#4775d1", "grey95"))((-min/(-min+max)*100)-0.5), colorRampPalette("grey95")(1),
          colorRampPalette(c("grey95", "#b3003b"), bias=1)((max/(-min+max)*100)-0.5))

Residue_order<-unique(NS_ddG_ripk1_median_melted$Residue)

heatmap_1<-ggplot(NS_ddG_ripk1_median_melted, 
                aes(x=factor(Residue, levels=Residue_order), y=factor(variable, levels=c("median_ddG_9HR6_pos", "median_ddG_9HR9_pos", "median_ddG_8Z93_pos"), labels=c("9HR6", "9HR9", "8Z93")), 
                    fill=value))+
  geom_tile(color='white', size=1)+
  scale_fill_gradientn(colours=cols, limits=c(min,max), na.value = "grey60", breaks=c(-3, 0, 5, 10, 15), labels=c(-3, 0, 5, 10, ">15"))+
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
scatterddG_4_<-ggplot(NS_ddG_ripk1_median[!NS_ddG_ripk1_median$Residue %in% c("G538", "G542"),], 
                      aes(x=median_ddG_9HR6_pos, y=median_ddG_9HR9_pos, label=Residue))+
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

ggsave(scatterddG_4_, path=path, file="ddG_9HR6vs9HR9_residues_noG.jpg", width=5, height=4)
ggsave(scatterddG_4_, path=path, file="ddG_9HR6vs9HR9_residues_noG.pdf", width=5, height=4)

##
# Label residues with opposite ddG
NS_ddG_ripk1_median$ddG_class<-"opposite_effect"
NS_ddG_ripk1_median[NS_ddG_ripk1_median$median_ddG_9HR6_pos>2 & NS_ddG_ripk1_median$median_ddG_8Z93_pos>2, "ddG_class"]<-"destabilizing_effect"
NS_ddG_ripk1_median[NS_ddG_ripk1_median$median_ddG_9HR6_pos<2 & NS_ddG_ripk1_median$median_ddG_8Z93_pos<2, "ddG_class"]<-"stabilizing_effect"

# Label them
NS_ddG_ripk1_median$ddG_class_ID<-""
NS_ddG_ripk1_median[NS_ddG_ripk1_median$ddG_class == "opposite_effect", "ddG_class_ID"]<-NS_ddG_ripk1_median[NS_ddG_ripk1_median$ddG_class == "opposite_effect", "Residue"]

# Classify the residues based on the comparison between structures
NS_ddG_ripk1_median$orientation<-"exposed"
NS_ddG_ripk1_median[NS_ddG_ripk1_median$Pos %in% c(531, 533, 539, 540, 541, 543, 545, 547, 549), "orientation"]<-"buried"
NS_ddG_ripk1_median[NS_ddG_ripk1_median$Pos %in% c(538, 542), "orientation"]<-"glycine"
NS_ddG_ripk1_median[NS_ddG_ripk1_median$Pos %in% c(535, 536, 543), "orientation"]<-"flipped"

#
scatterddG_5_<-ggplot(NS_ddG_ripk1_median[!NS_ddG_ripk1_median$Residue %in% c("G538", "G542"),], 
                     aes(x=median_ddG_9HR6_pos, y=median_ddG_8Z93_pos, label=Residue))+
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

ggsave(scatterddG_5_, path=path, file="ddG_9HR6vs8Z93_residues_noG.jpg", width=5, height=4)
ggsave(scatterddG_5_, path=path, file="ddG_9HR6vs8Z93_residues_noG.pdf", width=5, height=4)

cor.test(NS_ddG_ripk1_median[!NS_ddG_ripk1_median$Residue %in% c("G538", "G542"),]$median_ddG_8Z93_pos,
         NS_ddG_ripk1_median[!NS_ddG_ripk1_median$Residue %in% c("G538", "G542"),]$median_ddG_9HR6_pos)

##
# Label residues with opposite ddG
NS_ddG_ripk1_median$ddG_class<-"opposite_effect"
NS_ddG_ripk1_median[NS_ddG_ripk1_median$median_ddG_9HR9_pos>2 & NS_ddG_ripk1_median$median_ddG_8Z93_pos>2, "ddG_class"]<-"destabilizing_effect"
NS_ddG_ripk1_median[NS_ddG_ripk1_median$median_ddG_9HR9_pos<2 & NS_ddG_ripk1_median$median_ddG_8Z93_pos<2, "ddG_class"]<-"stabilizing_effect"

# Label them
NS_ddG_ripk1_median$ddG_class_ID<-""
NS_ddG_ripk1_median[NS_ddG_ripk1_median$ddG_class == "opposite_effect", "ddG_class_ID"]<-NS_ddG_ripk1_median[NS_ddG_ripk1_median$ddG_class == "opposite_effect", "Residue"]

# Classify the residues based on the comparison between structures
NS_ddG_ripk1_median$orientation<-"exposed"
NS_ddG_ripk1_median[NS_ddG_ripk1_median$Pos %in% c(531, 533, 539, 540, 541, 543, 545, 547, 549), "orientation"]<-"buried"
NS_ddG_ripk1_median[NS_ddG_ripk1_median$Pos %in% c(538, 542), "orientation"]<-"glycine"
NS_ddG_ripk1_median[NS_ddG_ripk1_median$Pos %in% c(535, 536, 543), "orientation"]<-"flipped"


#
scatterddG_6_<-ggplot(NS_ddG_ripk1_median[!NS_ddG_ripk1_median$Residue %in% c("G538", "G542"),], 
                      aes(x=median_ddG_9HR9_pos, y=median_ddG_8Z93_pos, label=Residue))+
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

ggsave(scatterddG_6_, path=path, file="ddG_9HR9vs8Z93_residues_noG.jpg", width=5, height=4)
ggsave(scatterddG_6_, path=path, file="ddG_9HR9vs8Z93_residues_noG.pdf", width=5, height=4)


################################################################################
# 3: correlation between ddG and NS for the different structures
scatterddG_7<-ggplot(NS_ddG_ripk1, 
                     aes(x=nscore_c, y=ddG_foldx_9HR6))+
  geom_hline(yintercept=0, color="grey50")+
  geom_vline(xintercept=0, color="grey50")+
  stat_cor()+ 
  geom_point(alpha=0.5)+
  xlim(-9, 2)+
  theme_classic()+
  theme()

scatterddG_7

ggsave(scatterddG_7, path=path, file="ddG_9HR6_nscore.jpg", width=4, height=4)
ggsave(scatterddG_7, path=path, file="ddG_9HR6_nscore.pdf", width=4, height=4)

#
scatterddG_8<-ggplot(NS_ddG_ripk1, 
                     aes(x=nscore_c, y=ddG_foldx_9HR9))+
  geom_hline(yintercept=0, color="grey50")+
  geom_vline(xintercept=0, color="grey50")+
  stat_cor()+ 
  geom_point(alpha=0.5)+
  xlim(-9, 2)+
  theme_classic()+
  theme()

scatterddG_8

ggsave(scatterddG_8, path=path, file="ddG_9HR9_nscore.jpg", width=4, height=4)
ggsave(scatterddG_8, path=path, file="ddG_9HR9_nscore.pdf", width=4, height=4)

#
scatterddG_9<-ggplot(NS_ddG_ripk1, 
                     aes(x=nscore_c, y=ddG_foldx_8Z93))+
  geom_hline(yintercept=0, color="grey50")+
  geom_vline(xintercept=0, color="grey50")+
  stat_cor()+ 
  geom_point(alpha=0.5)+
  xlim(-9, 2)+
  theme_classic()+
  theme()

scatterddG_9

ggsave(scatterddG_9, path=path, file="ddG_8Z93_nscore.jpg", width=4, height=4)
ggsave(scatterddG_9, path=path, file="ddG_8Z93_nscore.pdf", width=4, height=4)

################################################################################
# 4: correlation between ddG and NS for the different structures facet by structural elements

NS_ddG_ripk1_<-drop_na(NS_ddG_ripk1_)

scatterddG_10<-ggplot(NS_ddG_ripk1_[NS_ddG_ripk1_$structure != "Disordered-2",], 
                      aes(x=nscore_c, y=ddG_foldx_9HR9))+
  geom_hline(yintercept=0, color="grey50")+
  geom_vline(xintercept=0, color="grey50")+
  stat_cor()+ 
  geom_point(alpha=0.5)+
  facet_wrap(~factor(structure, levels=c("Disordered-1", "β-strand-1", "Loop-1", "β-strand-2", "Loop-2", "β-strand-3", "Disordered-2")))+
  theme_bw()+
  theme()

scatterddG_10

ggsave(scatterddG_10, path=path, file="ddG_9HR9_nscore_strtype.jpg", width=6, height=4)
ggsave(scatterddG_10, path=path, file="ddG_9HR9_nscore_strtype.pdf", width=6, height=4)

# Correlation per structure element

corr_vector=c()
for(i in NS_ddG_ripk1_$structure){
  correlation<-cor.test(NS_ddG_ripk1_[NS_ddG_ripk1_$structure==i,]$nscore_c,
                        NS_ddG_ripk1_[NS_ddG_ripk1_$structure==i,]$ddG_foldx_9HR9, use="complete.obs", method="pearson")
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

ggsave(p_bars_1, file="ddG_9HR9_nscore_strtype_bars.jpg", width = 4, height = 3, path=path)
ggsave(p_bars_1, file="ddG_9HR9_nscore_strtype_bars.pdf", width = 4, height = 3, path=path)

#
scatterddG_11<-ggplot(NS_ddG_ripk1_[NS_ddG_ripk1_$structure != "Disordered-2",], 
                      aes(x=nscore_c, y=ddG_foldx_9HR6))+
  geom_hline(yintercept=0, color="grey50")+
  geom_vline(xintercept=0, color="grey50")+
  stat_cor()+ 
  geom_point(alpha=0.5)+
  facet_wrap(~factor(structure, levels=c("Disordered-1", "β-strand-1", "Loop-1", "β-strand-2", "Loop-2", "β-strand-3", "Disordered-2")))+
  theme_bw()+
  theme()

scatterddG_11

ggsave(scatterddG_11, path=path, file="ddG_9HR6_nscore_strtype.jpg", width=6, height=4)
ggsave(scatterddG_11, path=path, file="ddG_9HR6_nscore_strtype.pdf", width=6, height=4)

# Correlation per structure element

corr_vector=c()
for(i in NS_ddG_ripk1_$structure){
  correlation<-cor.test(NS_ddG_ripk1_[NS_ddG_ripk1_$structure==i,]$nscore_c,
                        NS_ddG_ripk1_[NS_ddG_ripk1_$structure==i,]$ddG_foldx_9HR6, use="complete.obs", method="pearson")
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

ggsave(p_bars_2, file="ddG_9HR6_nscore_strtype_bars.jpg", width = 4, height = 3, path=path)
ggsave(p_bars_2, file="ddG_9HR6_nscore_strtype_bars.pdf", width = 4, height = 3, path=path)

#
scatterddG_12<-ggplot(NS_ddG_ripk1_[NS_ddG_ripk1_$structure != "Disordered-2",], 
                      aes(x=nscore_c, y=ddG_foldx_8Z93))+
  geom_hline(yintercept=0, color="grey50")+
  geom_vline(xintercept=0, color="grey50")+
  stat_cor()+ 
  geom_point(alpha=0.5)+
  facet_wrap(~factor(structure, levels=c("Disordered-1", "β-strand-1", "Loop-1", "β-strand-2", "Loop-2", "β-strand-3", "Disordered-2")))+
  theme_bw()+
  theme()

scatterddG_12

ggsave(scatterddG_12, path=path, file="ddG_8Z93_nscore_strtype.jpg", width=6, height=4)
ggsave(scatterddG_12, path=path, file="ddG_8Z93_nscore_strtype.pdf", width=6, height=4)

# Correlation per structure element

corr_vector=c()
for(i in NS_ddG_ripk1_$structure){
  correlation<-cor.test(NS_ddG_ripk1_[NS_ddG_ripk1_$structure==i,]$nscore_c,
                        NS_ddG_ripk1_[NS_ddG_ripk1_$structure==i,]$ddG_foldx_8Z93, use="complete.obs", method="pearson")
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

ggsave(p_bars_3, file="ddG_8Z93_nscore_strtype_bars.jpg", width = 4, height = 3, path=path)
ggsave(p_bars_3, file="ddG_8Z93_nscore_strtype_bars.pdf", width = 4, height = 3, path=path)

