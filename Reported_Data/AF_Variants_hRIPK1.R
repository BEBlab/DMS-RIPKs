library(tidyverse)
require(reshape2)
library(ggpubr)
library(seqinr)
require("ggrepel")
library(ggbreak) 

# your folder path
setwd("")

dir.create("03_AF_Variants_hRIPK1")
path="03_AF_Variants_hRIPK1"

### hRIPK1
# Nucleation Score
load("nscore_df_hRIPK1.RData")

singles_stops$sig_1<-F
singles_stops[singles_stops$p.adjust<0.01, "sig_1"]<-T

singles_stops$category_1<-"WT-like"
singles_stops[singles_stops$p.adjust<0.01 & singles_stops$nscore_c >0, "category_1"]<-"NS_inc"
singles_stops[singles_stops$p.adjust<0.01 & singles_stops$nscore_c <0, "category_1"]<-"NS_dec"

# Let's be super stringent with this comaprisons
singles_stops<-singles_stops[singles_stops$low_sigma == T | singles_stops$sig_1 == T, ]

singles_stops$Residue<-paste0(singles_stops$WT_AA, singles_stops$Pos)

### ddG
#PDB 9HR6
load("ddG_9HR6.RData")

# Residues classification based on side chain orientation
ddG_foldx_9HR6$orientation_9HR6<-"exposed"
ddG_foldx_9HR6[ddG_foldx_9HR6$Pos %in% c(448, 450, 452, 455, 456, 458, 459, 460, 462, 466, 468, 470), "orientation_9HR6"]<-"buried"
ddG_foldx_9HR6[ddG_foldx_9HR6$Pos %in% c(457, 461), "orientation_9HR6"]<-"glycine"

ddG_foldx_9HR6$sig_1_9HR6<-F
ddG_foldx_9HR6[ddG_foldx_9HR6$p.adjust_9HR6<0.01, "sig_1_9HR6"]<-T

ddG_foldx_9HR6$category_1_9HR6<-"WT-like"
ddG_foldx_9HR6[ddG_foldx_9HR6$p.adjust_9HR6<0.01 & ddG_foldx_9HR6$ddG_foldx_9HR6 >2, "category_1_9HR6"]<-"ddG_inc"
ddG_foldx_9HR6[ddG_foldx_9HR6$p.adjust_9HR6<0.01 & ddG_foldx_9HR6$ddG_foldx_9HR6 <0, "category_1_9HR6"]<-"ddG_dec"

ddG_foldx_9HR6<-ddG_foldx_9HR6[ddG_foldx_9HR6$low_ddG_foldx_SD_9HR6 == T | ddG_foldx_9HR6$sig_1_9HR6 == T, ]


#PDB 9HR9
load("ddG_9HR9.RData")

# Residues classification based on side chain orientation
ddG_foldx_9HR9$orientation_9HR9<-"exposed"
ddG_foldx_9HR9[ddG_foldx_9HR9$Pos %in% c(450, 452, 454, 456, 458, 459, 460, 466, 468), "orientation_9HR9"]<-"buried"
ddG_foldx_9HR9[ddG_foldx_9HR9$Pos %in% c(457, 461), "orientation_9HR9"]<-"glycine"

ddG_foldx_9HR9$sig_1_9HR9<-F
ddG_foldx_9HR9[ddG_foldx_9HR9$p.adjust_9HR9<0.01, "sig_1_9HR9"]<-T

ddG_foldx_9HR9$category_1_9HR9<-"WT-like"
ddG_foldx_9HR9[ddG_foldx_9HR9$p.adjust_9HR9<0.01 & ddG_foldx_9HR9$ddG_foldx_9HR9 >2, "category_1_9HR9"]<-"ddG_inc"
ddG_foldx_9HR9[ddG_foldx_9HR9$p.adjust_9HR9<0.01 & ddG_foldx_9HR9$ddG_foldx_9HR9 <0, "category_1_9HR9"]<-"ddG_dec"

ddG_foldx_9HR9<-ddG_foldx_9HR9[ddG_foldx_9HR9$low_ddG_foldx_SD_9HR9 == T | ddG_foldx_9HR9$sig_1_9HR9 == T, ]


#PDB 8Z93
load("ddG_8Z93.RData")

# Residues classification based on side chain orientation
ddG_foldx_8Z93$orientation_8Z93<-"exposed" 
ddG_foldx_8Z93[ddG_foldx_8Z93$Pos %in% c(450, 452, 454, 456, 458, 459, 460, 466, 468), "orientation_8Z93"]<-"buried"
ddG_foldx_8Z93[ddG_foldx_8Z93$Pos %in% c(457, 461), "orientation_8Z93"]<-"glycine"

ddG_foldx_8Z93$sig_1_8Z93<-F
ddG_foldx_8Z93[ddG_foldx_8Z93$p.adjust_8Z93<0.01, "sig_1_8Z93"]<-T

ddG_foldx_8Z93$category_1_8Z93<-"WT-like"
ddG_foldx_8Z93[ddG_foldx_8Z93$p.adjust_8Z93<0.01 & ddG_foldx_8Z93$ddG_foldx_8Z93 >2, "category_1_8Z93"]<-"ddG_inc"
ddG_foldx_8Z93[ddG_foldx_8Z93$p.adjust_8Z93<0.01 & ddG_foldx_8Z93$ddG_foldx_8Z93 <0, "category_1_8Z93"]<-"ddG_dec"

ddG_foldx_8Z93<-ddG_foldx_8Z93[ddG_foldx_8Z93$low_ddG_foldx_SD_8Z93 == T | ddG_foldx_8Z93$sig_1_8Z93 == T, ]

###
NS_ddG_RIPK1<-full_join(ddG_foldx_9HR9[c("ID", "ddG_foldx_9HR9", "ddG_foldx_SD_9HR9", "orientation_9HR9", "category_1_9HR9")], 
                        ddG_foldx_9HR6[c("ID", "ddG_foldx_9HR6", "ddG_foldx_SD_9HR6", "orientation_9HR6", "category_1_9HR6")], by="ID")

NS_ddG_RIPK1<-full_join(NS_ddG_RIPK1, ddG_foldx_8Z93[c("ID", "ddG_foldx_8Z93", "ddG_foldx_SD_8Z93", "orientation_8Z93", "category_1_8Z93")], by="ID")

NS_ddG_RIPK1<-full_join(NS_ddG_RIPK1, singles_stops[c("ID", "Pos", "Residue", "nscore_c", "sigma", "low_sigma", "sig_10", "category_10", "sig_1", "category_1")], by="ID")


##############################################################################
gnomad_variants<-read.csv("gnomAD_v4.1.0_ENSG00000137275_2025_10_30_09_15_47_hRIPK1_filtered.csv")

singles_stops_gnomad<-full_join(NS_ddG_RIPK1, gnomad_variants, by="ID")

singles_stops_gnomad_1<-singles_stops_gnomad[singles_stops_gnomad$low_sigma == T | singles_stops_gnomad$sig_10 == T,]
singles_stops_gnomad_1<-singles_stops_gnomad_1[!is.na(singles_stops_gnomad_1),]
singles_stops_gnomad_2<-singles_stops_gnomad[singles_stops_gnomad$low_sigma == T | singles_stops_gnomad$sig_1 == T,]
singles_stops_gnomad_2<-singles_stops_gnomad_2[!is.na(singles_stops_gnomad_2),]

p_gnomad_ns<-ggplot(singles_stops_gnomad_2, aes(x=nscore_c, y=X.log10.AF.))+
  geom_point(aes(color=category_1), size=3)+
  geom_text_repel(aes(color=category_1, label=ID), size=3, min.segment.length = 0.2)+
  scale_color_manual(values=c("#CD3333", "#00008B", "grey"))+
  annotate("text", x = Inf, y = -Inf, vjust = -1, hjust = 1,  label = paste0("FDR 1%"))+
  labs(x="Nucleation Score", y="Allele Frequency")+
  scale_y_continuous(limits=c(4, 6.5), breaks=c(4, 5, 6), 
                     labels=c(1e-04, 1e-05, 1e-06))+
  scale_x_continuous(limits=c(-4, 3), breaks=c(-4, -2, 0, 2))+
  theme_classic()+
  theme(legend.position = "none")
p_gnomad_ns

ggsave(p_gnomad_ns, path=path, file="p_gnomad_ns_lowsigma_classified.jpg", width=4, height=4)
ggsave(p_gnomad_ns, path=path, file="p_gnomad_ns_lowsigma_classified.pdf", width=4, height=4)


#

RGC_variants<-read.csv("RGC_hRIPK1_251030_filtered.csv")

singles_stops_RGC<-full_join(NS_ddG_RIPK1, RGC_variants, by="ID")

singles_stops_RGC_1<-singles_stops_RGC[singles_stops_RGC$low_sigma == T | singles_stops_RGC$sig_10 == T,]
singles_stops_RGC_2<-singles_stops_RGC[singles_stops_RGC$low_sigma == T | singles_stops_RGC$sig_1 == T,]

p_RGC_ns<-ggplot(singles_stops_RGC_2, aes(x=nscore_c, y=X.log10.AAF.))+
  geom_point(aes(color=category_1), size=3)+
  geom_text_repel(aes(color=category_1, label=ID), size=3, min.segment.length = 0.2)+
  scale_color_manual(values=c("#CD3333", "#00008B", "grey"))+
  annotate("text", x = Inf, y = -Inf, vjust = -1, hjust = 1,  label = paste0("FDR 1%"))+
  labs(x="Nucleation Score", y="Allele Frequency")+
  scale_y_continuous(limits=c(4, 6.5), breaks=c(4, 5, 6), 
                     labels=c(1e-04, 1e-05, 1e-06))+
  scale_x_continuous(limits=c(-4, 2.5), breaks=c(-4, -2, 0, 2))+
  theme_classic()+
  theme(legend.position = "none")
p_RGC_ns

ggsave(p_RGC_ns, path=path, file="p_RGC_ns_lowsigma_classified.jpg", width=4, height=4)
ggsave(p_RGC_ns, path=path, file="p_RGC_ns_lowsigma_classified.pdf", width=4, height=4)

################################################################################
# 9HR9
p_gnomad_ns_2<-ggplot(singles_stops_gnomad_2, aes(x=ddG_foldx_9HR9, y=X.log10.AF.))+
  geom_point(aes(color=category_1_9HR9), size=3)+
  geom_text_repel(aes(color=category_1_9HR9, label=ID), size=3, min.segment.length = 0.2)+
  scale_color_manual(values=c("#4775d1", "#b3003b", "grey"))+
  annotate("text", x = Inf, y = -Inf, vjust = -1, hjust = 1,  label = paste0("FDR 1%"))+
  labs(x="ddG 9HR9", y="Allele Frequency")+
  scale_y_continuous(limits=c(4, 6.5), breaks=c(4, 5, 6), 
                     labels=c(1e-04, 1e-05, 1e-06))+
  scale_x_continuous(limits=c(-5, 40), breaks=c(-5, 0, 10, 20, 30))+
  theme_classic()+
  theme(legend.position = "none")
p_gnomad_ns_2

ggsave(p_gnomad_ns_2, path=path, file="p_gnomad_ddG_9HR9_lowsigma_classified.jpg", width=4, height=4)
ggsave(p_gnomad_ns_2, path=path, file="p_gnomad_ddG_9HR9_lowsigma_classified.pdf", width=4, height=4)

#

p_RGC_ns_2<-ggplot(singles_stops_RGC_2, aes(x=ddG_foldx_9HR9, y=X.log10.AAF.))+
  geom_point(aes(color=category_1_9HR9), size=3)+
  geom_text_repel(aes(color=category_1_9HR9, label=ID), size=3, min.segment.length = 0.2)+
  scale_color_manual(values=c("#4775d1", "#b3003b", "grey"))+
  annotate("text", x = Inf, y = -Inf, vjust = -1, hjust = 1,  label = paste0("FDR 1%"))+
  labs(x="ddG 9HR9", y="Allele Frequency")+
  scale_y_continuous(limits=c(4, 6.5), breaks=c(4, 5, 6), 
                     labels=c(1e-04, 1e-05, 1e-06))+
  scale_x_continuous(limits=c(-5, 25), breaks=c(-5, 0, 10, 20))+
  theme_classic()+
  theme(legend.position = "none")
p_RGC_ns_2

ggsave(p_RGC_ns_2, path=path, file="p_RGC_ddG_9HR9_lowsigma_classified.jpg", width=4, height=4)
ggsave(p_RGC_ns_2, path=path, file="p_RGC_ddG_9HR9_lowsigma_classified.pdf", width=4, height=4)

################################################################################
# 9HR6
p_gnomad_ns_3<-ggplot(singles_stops_gnomad_2, aes(x=ddG_foldx_9HR6, y=X.log10.AF.))+
  geom_point(aes(color=category_1_9HR6), size=3)+
  geom_text_repel(aes(color=category_1_9HR6, label=ID), size=3, min.segment.length = 0.2)+
  scale_color_manual(values=c("#4775d1", "#b3003b", "grey"))+
  annotate("text", x = Inf, y = -Inf, vjust = -1, hjust = 1,  label = paste0("FDR 1%"))+
  labs(x="ddG 9HR6", y="Allele Frequency")+
  scale_y_continuous(limits=c(4, 6.5), breaks=c(4, 5, 6), 
                     labels=c(1e-04, 1e-05, 1e-06))+
  scale_x_continuous(limits=c(-5, 22), breaks=c(-5, 0, 5, 10, 20))+
  #scale_x_break(c(11, 25))+
  theme_classic()+
  theme(legend.position = "none")
p_gnomad_ns_3

ggsave(p_gnomad_ns_3, path=path, file="p_gnomad_ddG_9HR6_lowsigma_classified.jpg", width=4, height=4)
ggsave(p_gnomad_ns_3, path=path, file="p_gnomad_ddG_9HR6_lowsigma_classified.pdf", width=4, height=4)

#

p_RGC_ns_3<-ggplot(singles_stops_RGC_2, aes(x=ddG_foldx_9HR6, y=X.log10.AAF.))+
  geom_point(aes(color=category_1_9HR6), size=3)+
  geom_text_repel(aes(color=category_1_9HR6, label=ID), size=3, min.segment.length = 0.2)+
  scale_color_manual(values=c("#4775d1", "#b3003b", "grey"))+
  annotate("text", x = Inf, y = -Inf, vjust = -1, hjust = 1,  label = paste0("FDR 1%"))+
  labs(x="ddG 9HR6", y="Allele Frequency")+
  scale_y_continuous(limits=c(4, 6.5), breaks=c(4, 5, 6), 
                     labels=c(1e-04, 1e-05, 1e-06))+
  scale_x_continuous(limits=c(-5, 30), breaks=c(-5, 0, 5, 10, 25, 30))+
  scale_x_break(c(11, 25))+
  theme_classic()+
  theme(legend.position = "none")
p_RGC_ns_3

ggsave(p_RGC_ns_3, path=path, file="p_RGC_ddG_9HR6_lowsigma_classified.jpg", width=4, height=4)
ggsave(p_RGC_ns_3, path=path, file="p_RGC_ddG_9HR6_lowsigma_classified.pdf", width=4, height=4)

################################################################################
# 8Z93
p_gnomad_ns_4<-ggplot(singles_stops_gnomad_2, aes(x=ddG_foldx_8Z93, y=X.log10.AF.))+
  geom_point(aes(color=category_1_8Z93), size=3)+
  geom_text_repel(aes(color=category_1_8Z93, label=ID), size=3, min.segment.length = 0.2)+
  scale_color_manual(values=c("#4775d1", "#b3003b", "grey"))+
  annotate("text", x = Inf, y = -Inf, vjust = -1, hjust = 1,  label = paste0("FDR 1%"))+
  labs(x="ddG 8Z93", y="Allele Frequency")+
  scale_y_continuous(limits=c(4, 6.5), breaks=c(4, 5, 6), 
                     labels=c(1e-04, 1e-05, 1e-06))+
  scale_x_continuous(limits=c(-5, 20), breaks=c(-5, 0, 5, 10, 15, 20))+
  theme_classic()+
  theme(legend.position = "none")
p_gnomad_ns_4

ggsave(p_gnomad_ns_4, path=path, file="p_gnomad_ddG_8Z93_lowsigma_classified.jpg", width=4, height=4)
ggsave(p_gnomad_ns_4, path=path, file="p_gnomad_ddG_8Z93_lowsigma_classified.pdf", width=4, height=4)

#

p_RGC_ns_4<-ggplot(singles_stops_RGC_2, aes(x=ddG_foldx_8Z93, y=X.log10.AAF.))+
  geom_point(aes(color=category_1_8Z93), size=3)+
  geom_text_repel(aes(color=category_1_8Z93, label=ID), size=3, min.segment.length = 0.2)+
  scale_color_manual(values=c("#4775d1", "#b3003b", "grey"))+
  annotate("text", x = Inf, y = -Inf, vjust = -1, hjust = 1,  label = paste0("FDR 1%"))+
  labs(x="ddG 8Z93", y="Allele Frequency")+
  scale_y_continuous(limits=c(4, 6.5), breaks=c(4, 5, 6), 
                     labels=c(1e-04, 1e-05, 1e-06))+
  scale_x_continuous(limits=c(-5, 15), breaks=c(-5, 0, 5, 10, 15))+
  theme_classic()+
  theme(legend.position = "none")
p_RGC_ns_4

ggsave(p_RGC_ns_4, path=path, file="p_RGC_ddG_8Z93_lowsigma_classified.jpg", width=4, height=4)
ggsave(p_RGC_ns_4, path=path, file="p_RGC_ddG_8Z93_lowsigma_classified.pdf", width=4, height=4)

