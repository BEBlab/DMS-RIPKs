library(tidyverse)
library(ggpubr)
library(ggrepel)

# your folder path
setwd("")

# Create a folder for the figures
dir.create("02_Distribution_RepsCorrelation")
path="02_Distribution_RepsCorrelation"

# open processed data
load("nscore_df_mRIP3.RData")

#################################################################################
# Nucleation score distribution per variant type (Missense, Nonsense and Synonymous)

dist_singles<-singles_stops[,c("Mut", "nscore_c")]
dist_singles<-rbind(dist_singles, silent[silent$Nmut_codons==1,c("Mut", "nscore_c")])
dist_singles$type<-"Missense"
dist_singles[dist_singles$Mut=="*",]$type<-"Nonsense"
dist_singles[dist_singles$Mut=="silent",]$type<-"Synonymous"

p_hist<-ggplot(dist_singles, aes(x=factor(type, levels=c("Synonymous", "Nonsense", "Missense")), y=nscore_c))+
  geom_jitter(color="grey20", height = 0, width = 0.1, size=0.5)+
  geom_violin(fill="grey90", alpha=0.5, scale="width")+
  geom_boxplot(width = 0.02, outlier.shape=NA)+
  labs(x="", y="Nucleation Score", fill="")+
  theme_classic()

p_hist

ggsave(p_hist, file="p_hist.jpg", path=path, width = 4, height = 4)
ggsave(p_hist, file="p_hist.pdf", path=path, width = 4, height = 4)

################################################################################
# correlation between replicates

pair_plot <- ggpairs(
  singles_stops[c("nscore1_c", "nscore2_c", "nscore3_c")],
  diag = list(continuous = "blankDiag"),
  lower = list(continuous = wrap("points", color = "grey80", alpha = 0.5)))+
  theme_bw()

pair_plot

ggsave(pair_plot, path=path, file="pair_plot.pdf", width=4, height=4)
###

#1 and 2
corr<-cor.test(singles_stops$nscore1_c, singles_stops$nscore2_c, use="complete.obs")
R<-corr$estimate
p<-corr$p.value

p_corr_12<-ggplot(singles_stops, aes(x=nscore1_c, y=nscore2_c) )+
  geom_point(size=2, alpha=0.5, color="black")+
  theme_classic()+
  annotate("text", x = Inf, y = -Inf, vjust = -2, hjust = 1.5,  label = paste0("R=", round(R, 2)))+
  annotate("text", x = Inf, y = -Inf, vjust = -0.5, hjust = 1.2, label = paste0("p=",format(p, digits = 2, scientific = T)))+
  geom_smooth(method = "lm", linetype = 2, size=1, se = F, color = "grey")

p_corr_12

ggsave(p_corr_12,path=path, file="p_corr_rep1vs2.jpg", width=5, height=4)

#1 and 3
corr<-cor.test(singles_stops$nscore1_c, singles_stops$nscore3_c, use="complete.obs")
R<-corr$estimate
p<-corr$p.value

p_corr_13<-ggplot(singles_stops, aes(x=nscore1_c, y=nscore3_c) )+
  geom_point(size=2, alpha=0.5, color="black")+
  theme_classic()+
  annotate("text", x = Inf, y = -Inf, vjust = -2, hjust = 1.5,  label = paste0("R=", round(R, 2)))+
  annotate("text", x = Inf, y = -Inf, vjust = -0.5, hjust = 1.2, label = paste0("p=",format(p, digits = 2, scientific = T)))+
  geom_smooth(method = "lm", linetype = 2, size=1, se = F, color = "grey")

p_corr_13

ggsave(p_corr_13, path=path, file="p_corr_reps1vs3.jpg", width=5, height=4)

#2 and 3
corr<-cor.test(singles_stops$nscore2_c, singles_stops$nscore3_c, use="complete.obs")
R<-corr$estimate
p<-corr$p.value

p_corr_23<-ggplot(singles_stops, aes(x=nscore2_c, y=nscore3_c) )+
  geom_point(size=2, alpha=0.5, color="black")+
  theme_classic()+
  annotate("text", x = Inf, y = -Inf, vjust = -2, hjust = 1.5,  label = paste0("R=", round(R, 2)))+
  annotate("text", x = Inf, y = -Inf, vjust = -0.5, hjust = 1.2, label = paste0("p=",format(p, digits = 2, scientific = T)))+
  geom_smooth(method = "lm", linetype = 2, size=1, se = F, color = "grey")

p_corr_23

ggsave(p_corr_23,path=path, file="p_corr_reps2vs3.jpg", width=5, height=4)

#
pair_plot <- ggpairs(
  singles_stops[c("nscore1_c", "nscore2_c", "nscore3_c")],
  diag = list(continuous = "blankDiag"),
  lower = list(continuous = wrap("points", color = "grey80", alpha = 0.5)))+
  theme_bw()

pair_plot

ggsave(pair_plot, path=path, file="pair_plot.pdf", width=4, height=4)
