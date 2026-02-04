library(tidyverse)
require(reshape2)
library(ggpubr)
library(seqinr)
require("ggrepel")
library(ggbreak) 

# your folder path
setwd("")

dir.create("03_AF_Variants_hRIPK3")
path="03_AF_Variants_hRIPK3"

### hRIPK3
# Nucleation Score
load("nscore_df_hRIPK3.RData")

singles_stops$sig_1<-F
singles_stops[singles_stops$p.adjust<0.01, "sig_1"]<-T

singles_stops$category_1<-"WT-like"
singles_stops[singles_stops$p.adjust<0.01 & singles_stops$nscore_c >0, "category_1"]<-"NS_inc"
singles_stops[singles_stops$p.adjust<0.01 & singles_stops$nscore_c <0, "category_1"]<-"NS_dec"

# Let's be super stringent with this comaprisons
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

NS_ddG_RIPK3<-full_join(NS_ddG_RIPK3, singles_stops[c("ID", "Pos", "Residue", "nscore_c", "sigma", "low_sigma", "sig_10", "category_10", "sig_1", "category_1")], by="ID")


##############################################################################
gnomad_variants<-read.csv("gnomAD_v4.1.0_ENSG00000129465_2025_10_30_09_15_23_hRIPK3_filtered.csv")

singles_stops_gnomad<-full_join(NS_ddG_RIPK3, gnomad_variants, by="ID")

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
  scale_x_continuous(limits=c(-2, 1.5), breaks=c(-2, -1, 0, 1))+
  theme_classic()+
  theme(legend.position = "none")
p_gnomad_ns

ggsave(p_gnomad_ns, path=path, file="p_gnomad_ns_lowsigma_classified.jpg", width=4, height=4)
ggsave(p_gnomad_ns, path=path, file="p_gnomad_ns_lowsigma_classified.pdf", width=4, height=4)


#

RGC_variants<-read.csv("RGC_hRIPK3_251030_filtered.csv")

singles_stops_RGC<-full_join(NS_ddG_RIPK3, RGC_variants, by="ID")

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
  scale_x_continuous(limits=c(-2, 1.5), breaks=c(-2, -1, 0, 1))+
  theme_classic()+
  theme(legend.position = "none")
p_RGC_ns

ggsave(p_RGC_ns, path=path, file="p_RGC_ns_lowsigma_classified.jpg", width=4, height=4)
ggsave(p_RGC_ns, path=path, file="p_RGC_ns_lowsigma_classified.pdf", width=4, height=4)

################################################################################
# 7DA4
p_gnomad_ns_2<-ggplot(singles_stops_gnomad_2, aes(x=ddG_foldx_7DA4, y=X.log10.AF.))+
  geom_point(aes(color=category_1_7DA4), size=3)+
  geom_text_repel(aes(color=category_1_7DA4, label=ID), size=3, min.segment.length = 0.2)+
  scale_color_manual(values=c("#4775d1", "#b3003b", "grey"))+
  annotate("text", x = Inf, y = -Inf, vjust = -1, hjust = 1,  label = paste0("FDR 1%"))+
  labs(x="ddG 7DA4", y="Allele Frequency")+
  scale_y_continuous(limits=c(4, 6.5), breaks=c(4, 5, 6), 
                     labels=c(1e-04, 1e-05, 1e-06))+
  scale_x_continuous(limits=c(-5, 15), breaks=c(-5, 0, 5, 10, 15))+
  theme_classic()+
  theme(legend.position = "none")
p_gnomad_ns_2

ggsave(p_gnomad_ns_2, path=path, file="p_gnomad_ddG_7DA4_lowsigma_classified.jpg", width=4, height=4)
ggsave(p_gnomad_ns_2, path=path, file="p_gnomad_ddG_7DA4_lowsigma_classified.pdf", width=4, height=4)

#

p_RGC_ns_2<-ggplot(singles_stops_RGC_2, aes(x=ddG_foldx_7DA4, y=X.log10.AAF.))+
  geom_point(aes(color=category_1_7DA4), size=3)+
  geom_text_repel(aes(color=category_1_7DA4, label=ID), size=3, min.segment.length = 0.2)+
  scale_color_manual(values=c("#4775d1", "#b3003b", "grey"))+
  annotate("text", x = Inf, y = -Inf, vjust = -1, hjust = 1,  label = paste0("FDR 1%"))+
  labs(x="ddG 7DA4", y="Allele Frequency")+
  scale_y_continuous(limits=c(4, 6.5), breaks=c(4, 5, 6), 
                     labels=c(1e-04, 1e-05, 1e-06))+
  scale_x_continuous(limits=c(-5, 15), breaks=c(-5, 0, 5, 10, 15))+
  theme_classic()+
  theme(legend.position = "none")
p_RGC_ns_2

ggsave(p_RGC_ns_2, path=path, file="p_RGC_ddG_7DA4_lowsigma_classified.jpg", width=4, height=4)
ggsave(p_RGC_ns_2, path=path, file="p_RGC_ddG_7DA4_lowsigma_classified.pdf", width=4, height=4)

################################################################################
# 7DAC
p_gnomad_ns_3<-ggplot(singles_stops_gnomad_2, aes(x=ddG_foldx_7DAC, y=X.log10.AF.))+
  geom_point(aes(color=category_1_7DAC), size=3)+
  geom_text_repel(aes(color=category_1_7DAC, label=ID), size=3, min.segment.length = 0.2)+
  scale_color_manual(values=c("#4775d1", "#b3003b", "grey"))+
  annotate("text", x = Inf, y = -Inf, vjust = -1, hjust = 1,  label = paste0("FDR 1%"))+
  labs(x="ddG 7DAC", y="Allele Frequency")+
  scale_y_continuous(limits=c(4, 6.5), breaks=c(4, 5, 6), 
                     labels=c(1e-04, 1e-05, 1e-06))+
  scale_x_continuous(limits=c(-5, 20), breaks=c(-5, 0, 5, 10, 15, 20))+
  scale_x_break(c(11, 18))+
  theme_classic()+
  theme(legend.position = "none")
p_gnomad_ns_3

ggsave(p_gnomad_ns_3, path=path, file="p_gnomad_ddG_7DAC_lowsigma_classified.jpg", width=4, height=4)
ggsave(p_gnomad_ns_3, path=path, file="p_gnomad_ddG_7DAC_lowsigma_classified.pdf", width=4, height=4)

#

p_RGC_ns_3<-ggplot(singles_stops_RGC_2, aes(x=ddG_foldx_7DAC, y=X.log10.AAF.))+
  geom_point(aes(color=category_1_7DAC), size=3)+
  geom_text_repel(aes(color=category_1_7DAC, label=ID), size=3, min.segment.length = 0.2)+
  scale_color_manual(values=c("#4775d1", "#b3003b", "grey"))+
  annotate("text", x = Inf, y = -Inf, vjust = -1, hjust = 1,  label = paste0("FDR 1%"))+
  labs(x="ddG 7DAC", y="Allele Frequency")+
  scale_y_continuous(limits=c(4, 6.5), breaks=c(4, 5, 6), 
                     labels=c(1e-04, 1e-05, 1e-06))+
  scale_x_continuous(limits=c(-5, 20), breaks=c(-5, 0, 5, 10, 15, 20))+
  scale_x_break(c(11, 18))+
  theme_classic()+
  theme(legend.position = "none")
p_RGC_ns_3

ggsave(p_RGC_ns_3, path=path, file="p_RGC_ddG_7DAC_lowsigma_classified.jpg", width=4, height=4)
ggsave(p_RGC_ns_3, path=path, file="p_RGC_ddG_7DAC_lowsigma_classified.pdf", width=4, height=4)

################################################################################
# 8Z94
p_gnomad_ns_4<-ggplot(singles_stops_gnomad_2, aes(x=ddG_foldx_8Z94, y=X.log10.AF.))+
  geom_point(aes(color=category_1_8Z94), size=3)+
  geom_text_repel(aes(color=category_1_8Z94, label=ID), size=3, min.segment.length = 0.2)+
  scale_color_manual(values=c("#4775d1", "#b3003b", "grey"))+
  annotate("text", x = Inf, y = -Inf, vjust = -1, hjust = 1,  label = paste0("FDR 1%"))+
  labs(x="ddG 8Z94", y="Allele Frequency")+
  scale_y_continuous(limits=c(4, 6.5), breaks=c(4, 5, 6), 
                     labels=c(1e-04, 1e-05, 1e-06))+
  scale_x_continuous(limits=c(-5, 15), breaks=c(-5, 0, 5, 10, 15))+
  theme_classic()+
  theme(legend.position = "none")
p_gnomad_ns_4

ggsave(p_gnomad_ns_4, path=path, file="p_gnomad_ddG_8Z94_lowsigma_classified.jpg", width=4, height=4)
ggsave(p_gnomad_ns_4, path=path, file="p_gnomad_ddG_8Z94_lowsigma_classified.pdf", width=4, height=4)

#

p_RGC_ns_4<-ggplot(singles_stops_RGC_2, aes(x=ddG_foldx_8Z94, y=X.log10.AAF.))+
  geom_point(aes(color=category_1_8Z94), size=3)+
  geom_text_repel(aes(color=category_1_8Z94, label=ID), size=3, min.segment.length = 0.2)+
  scale_color_manual(values=c("#4775d1", "#b3003b", "grey"))+
  annotate("text", x = Inf, y = -Inf, vjust = -1, hjust = 1,  label = paste0("FDR 1%"))+
  labs(x="ddG 8Z94", y="Allele Frequency")+
  scale_y_continuous(limits=c(4, 6.5), breaks=c(4, 5, 6), 
                     labels=c(1e-04, 1e-05, 1e-06))+
  scale_x_continuous(limits=c(-5, 15), breaks=c(-5, 0, 5, 10, 15))+
  theme_classic()+
  theme(legend.position = "none")
p_RGC_ns_4

ggsave(p_RGC_ns_4, path=path, file="p_RGC_ddG_8Z94_lowsigma_classified.jpg", width=4, height=4)
ggsave(p_RGC_ns_4, path=path, file="p_RGC_ddG_8Z94_lowsigma_classified.pdf", width=4, height=4)

