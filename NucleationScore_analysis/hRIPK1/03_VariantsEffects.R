library(tidyverse)
require(reshape2)
library(ggpubr)
library(seqinr)

# your folder path
setwd("")

# Create a folder for the figures
dir.create("03_VariantsEffects")
path="03_VariantsEffects"

# open processed data
load("nscore_df_hRIPK1.RData")

# position where the library starts
start_position<-522
  
# amino acid sequence mutagenised
peptide_seq<-"PPTDESIKYTIYNSTGIQIGAYNYMEIGGTSSSL"
peptide_seq<-c(strsplit(peptide_seq, '')[[1]])

peptide_seq_pos<-c()
for (n in seq_along(peptide_seq)){
  peptide_seq_pos<-c(peptide_seq_pos, paste0(peptide_seq[n], '\n', n+start_position))
}

# all the amino acids + stop codon
vectorAA <- c("G","A","V","L","M","I","F","Y","W","K","R","D","E","S","T","C","N","Q","H", "P", "*")

# color scale to be used
min<-min(singles_stops$nscore_c)
max<-max(singles_stops$nscore_c)
cols <- c(colorRampPalette(c( "brown3", "grey95"))((-min/(-min+max)*100)-0.5), colorRampPalette("grey95")(1),
          colorRampPalette(c("grey95",  "darkblue"), bias=1)((max/(-min+max)*100)-0.5))

##################################################################################

#####  heatmap (contains stops)

#add syn
positions<-(523:556)
syn.df<-data.frame(
  "aa_seq"="PPTDESIKYTIYNSTGIQIGAYNYMEIGGTSSSL",
  "WT_AA"= peptide_seq,
  "Mut"= peptide_seq,
  "Pos"= positions, 
  "sigma"=0, 
  "nscore_c"=0,
  "ID"="syn",
  "mean_count"=NA
  
)

heatmap_df<-rbind(singles_stops[,c("aa_seq","WT_AA", "Mut", "Pos", "sigma", "nscore_c", "ID", "mean_count")], syn.df)

heatmap_df$WT<-""
heatmap_df[heatmap_df$ID=="syn",]$WT<-"WT"

p_heatmap<-ggplot(heatmap_df)+
  geom_tile(aes(Pos,factor(Mut, levels=rev(vectorAA)), fill=nscore_c), color='white', size=1)+
  geom_text(aes(Pos, Mut, label=WT), size=3)+
  scale_x_continuous(breaks=seq(523,556), labels = peptide_seq_pos, expand=c(0,0))+
  theme_minimal()+
  labs(x="hRIPK1", y="Mutant amino acid", fill="Nucleation score")+
  theme(legend.title = element_text(size=14, face = "bold"),
        legend.text = element_text(size=14), 
        axis.title = element_text(size = 16),
        axis.text = element_text(size=14))+
  scale_fill_gradientn(colours=cols, limits=c(min,max), na.value = "grey60",breaks=c(-10, -5, -2, 0, 2))
p_heatmap

ggsave(p_heatmap, path=path, file="p_heatmap_nscore.jpg", width=15, height=8)

##### Violin plot per position
subsp_median_df<-as.data.frame(heatmap_df %>% group_by(Pos) %>% dplyr::summarise(median_p=median(nscore_c)))
heatmap_df<-left_join(heatmap_df, subsp_median_df)

p_violin_p<-ggplot(heatmap_df, aes(x=as.factor(Pos), y=nscore_c))+
  geom_hline(yintercept = 0, size=0.1)+
  geom_violin(scale = "width", size=0.2, aes(fill=median_p))+
  geom_boxplot(width=0.15, outlier.shape = NA, size=0.2)+
  scale_x_discrete(breaks=seq(523,556), labels = peptide_seq_pos, expand=c(0,0))+
  theme_bw()+
  labs(x="hRIPK1", y="Nucleation score", fill="Median NS")+
  theme(legend.title = element_text(size=14, face = "bold"),
        legend.text = element_text(size=14), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        plot.margin = unit(c(0,0,0,0), 'cm'))+
  scale_fill_gradientn(colours=cols, limits=c(min,max), na.value = "grey60", breaks=c(-10, -5, -2, 0, 2))

p_violin_p

ggsave(p_violin_p, path=path, file="p_violinplotp_nscore.jpg", width=15, height=8)

##### violin per change

subsm_median_df<-as.data.frame(heatmap_df %>% group_by(Mut) %>% dplyr::summarise(median_m=median(nscore_c)))
heatmap_df<-left_join(heatmap_df, subsm_median_df)

p_violin_m<-ggplot(heatmap_df, aes(x=factor(Mut, levels=rev(vectorAA)), y=nscore_c))+
  geom_hline(yintercept = 0, size=0.1)+
  geom_violin(scale = "width", size=0.2, aes(fill=median_m))+
  geom_boxplot(width=0.15, outlier.shape = NA, size=0.2)+
  theme_bw()+
  labs(x="hRIPK1", y="Nucleation score", fill="Median NS")+
  theme(legend.title = element_text(size=14, face = "bold"),
        legend.text = element_text(size=14), 
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=16),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        plot.margin = unit(c(0,0,0,0), 'cm'))+
  scale_fill_gradientn(colours=cols, limits=c(min,max), na.value = "grey60", breaks=c(-10, -5, -2, 0, 2))+ 
  coord_flip() 

p_violin_m

ggsave(p_violin_m, file="p_violinplot_mut.jpg",width=4, height=8, path=path)

#### stacked bars
# Considering stops
fdr_categories<-singles_stops[,c("aa_seq", "Pos", "WT_AA", "Mut", "ID", "nscore_c", "p.adjust", "category_10")]

fdr_categories$category<-"WT-like"
fdr_categories[(fdr_categories$p.adjust<0.1 & fdr_categories$nscore_c<0),]$category<- "NS-"
fdr_categories[(fdr_categories$p.adjust<0.1 & fdr_categories$nscore_c>0),]$category<- "NS+"


categories_p <- fdr_categories %>% group_by(Pos,category) %>% dplyr::summarise(Freq_p=n()) 

levels_cat = c("NS-", "WT-like", "NS+")
colors_cat<-c("#DF9292", "#F2F2F2", "#7979BE")

p_categories_p<-ggplot(categories_p, aes(fill=factor(category, levels=levels_cat), x=factor(Pos), y=Freq_p))+
  scale_x_discrete(breaks=seq(523,556), labels = peptide_seq_pos, expand=c(0,0))+
  theme_bw()+
  theme(legend.title=element_text(size=14, face = "bold"),
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color='black', size=0.25),
        axis.title = element_text(size=16),
        axis.text = element_text(size=14),
        plot.margin = unit(c(0,0,0,0), 'cm'))+
  geom_bar(position="stack", stat="identity", alpha=0.8, width=0.8)+
  scale_fill_manual("Nucleation (FDR=0.1)", values=colors_cat)+
  labs(x= "hRIPK1", y="Counts")

p_categories_p

ggsave(p_categories_p, path=path, file="p_categoriesp_nscore.jpg", width=15, height=8)

# changes from

categories_from <- fdr_categories %>% group_by(WT_AA, category) %>% dplyr::summarise(n=n())  %>% mutate(freq = n / sum(n) )

p_categories_from<-ggplot(categories_from, aes(fill=factor(category, levels=levels_cat), x=factor(WT_AA), y=freq))+
  theme_bw()+
  theme(legend.title=element_text(size=14, face = "bold"),
        legend.text = element_text(size=14),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color='black', size=0.25),
        axis.title = element_text(size=16),
        axis.text = element_text(size=14),
        plot.margin = unit(c(0,0,0,0), 'cm'))+
  geom_bar(position="stack", stat="identity", alpha=0.8, width=0.8)+
  scale_fill_manual("Nucleation (FDR=0.1)", values=colors_cat)+
  labs(x= "WT AA", y="Frequency")

p_categories_from

ggsave(p_categories_from, path=path, file="p_categoriesfrom_nscore.jpg", width=5, height=4)

######

categories_m <- fdr_categories %>% group_by(Mut, category) %>% dplyr::summarise(Freq_m=n()) 

p_categories_m<-ggplot(categories_m, aes(fill=factor(category, levels=levels_cat), y=factor(Mut, levels=rev(vectorAA)), x=Freq_m))+
  theme_bw()+
  theme(legend.title=element_text(size=14, face = "bold"),
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.title.y = element_blank(),
        axis.line = element_line(color='black', size=0.25),
        axis.title = element_text(size=16),
        axis.text = element_text(size=14),
        plot.margin = unit(c(0,0,0,0), 'cm'))+
  geom_bar(position="stack", stat="identity", alpha=0.8, width=0.8)+
  scale_fill_manual("Nucleation (FDR=0.1)", values=colors_cat)+
  labs(x= "Counts")

p_categories_m

ggsave(p_categories_m, path=path, file="p_categoriesm_nscore.jpg", width=15, height=8)
##################

p_hv <- ggarrange(p_heatmap + rremove("x.text"), p_violin_m + rremove("y.text"), p_categories_m + rremove("y.text"), 
                  ncol=3, align = "h",  legend="none", widths = c(1, 0.2, 0.2), common.legend=TRUE) 

p_hv

ggsave(p_hv, path=path, file="p_heatmap_violinm_categoreism.jpg", width=15, height=8)
ggsave(p_hv, path=path, file="p_heatmap_violinm_categoreism.pdf", width=15, height=8)


### Arrange plots
p_heatmap <- p_heatmap + theme(axis.title.x = element_blank())
p_hvc<-ggarrange(p_heatmap + rremove("x.text"), p_violin_p + rremove("x.text"), p_categories_p, 
                 nrow=3, align = "v",  legend="none", heights = c(1, 0.3, 0.3))

p_hvc

ggsave(p_hvc, path=path, file="p_heatmap_violin_categories_nscore.jpg", width=14, height=12)
ggsave(p_hvc, path=path, file="p_heatmap_violin_categories_nscore.pdf", width=14, height=12)

###
# RHIM in vs RHIM out

singles_stops$region<-"flanking"
singles_stops[singles_stops$Pos>=531 & singles_stops$Pos<=547, "region"]<-"RHIM"

categories_r <- singles_stops %>%
  filter(low_sigma == TRUE | sig_10 == TRUE) %>%
  group_by(region, category_10) %>%
  summarise(Freq_p = n(), .groups = "drop_last") %>%
  mutate(
    percent = 100 * Freq_p / sum(Freq_p)
  ) %>%
  ungroup()

levels_cat<-c("NS_dec", "WT-like", "NS_inc")

region_p<-ggplot(categories_r, aes(x=region, y=percent, fill=factor(category_10, levels=levels_cat)))+
  geom_bar(position="stack", stat="identity", alpha=0.8, width=0.8)+
  theme_classic()+
  scale_fill_manual("Nucleation (FDR=0.1)", values=colors_cat)+
  labs(x= "Region", y="Counts")
region_p

ggsave(region_p, path=path, file="Region_RHIM_stackedbars.jpg", width=4, height=4)
ggsave(region_p, path=path, file="Region_RHIM_stackedbars.pdf", width=4, height=4)
#

categories_r[categories_r$category_10 != "WT-like", "category_10"]<-"NS_disrupted"
categories_r_ <- categories_r %>% group_by(region, category_10) %>% dplyr::summarise(Freq_p=sum(Freq_p)) 

# Make a contingency table
tbl <- xtabs(Freq_p ~ region + category_10, data = categories_r_)

# Run chi-squared test
chisq.test(tbl)

##################################################################################
####  heatmap FDR

# work only with those variants that can be classified based on FDR and low sigma (normalized to fitness range)
singles_stops<-singles_stops[singles_stops$low_sigma == T | singles_stops$sig_10 == T,]

heatmap_fdr<-singles_stops[,c("aa_seq", "Pos", "WT_AA", "Mut", "ID", "nscore_c", "p.adjust", "category_10", "low_sigma")]

#add syn
syn.df<-data.frame(
  "aa_seq"="PPTDESIKYTIYNSTGIQIGAYNYMEIGGTSSSL",
  "WT_AA"= peptide_seq,
  "Mut"= peptide_seq,
  "Pos"= c(523:556), 
  "nscore_c"=0,
  "ID"="syn",
  "p.adjust"=NA,
  "category_10"="WT-like",
  "low_sigma"=FALSE
)

heatmap_fdr<-rbind(heatmap_fdr, syn.df)

heatmap_fdr$category<-"WT-like"

heatmap_fdr[(heatmap_fdr$p.adjust<0.25 & heatmap_fdr$nscore_c<0),]$category<- "NS- 25%"
heatmap_fdr[(heatmap_fdr$p.adjust<0.1 & heatmap_fdr$nscore_c<0),]$category<- "NS- 10%"
heatmap_fdr[(heatmap_fdr$p.adjust<0.05 & heatmap_fdr$nscore_c<0),]$category<- "NS- 5%"
heatmap_fdr[(heatmap_fdr$p.adjust<0.01 & heatmap_fdr$nscore_c<0),]$category<- "NS- 1%"

heatmap_fdr[(heatmap_fdr$p.adjust<0.25 & heatmap_fdr$nscore_c>0),]$category<- "NS+ 25%"
heatmap_fdr[(heatmap_fdr$p.adjust<0.1 & heatmap_fdr$nscore_c>0),]$category<- "NS+ 10%"
heatmap_fdr[(heatmap_fdr$p.adjust<0.05 & heatmap_fdr$nscore_c>0),]$category<- "NS+ 5%"
heatmap_fdr[(heatmap_fdr$p.adjust<0.01 & heatmap_fdr$nscore_c>0),]$category<- "NS+ 1%"

heatmap_fdr$WT<-""
heatmap_fdr[heatmap_fdr$ID=="syn",]$WT<-"WT"

#heatmap_fdr$SNVs<-""
#heatmap_fdr[(heatmap_fdr$aa_seq %in% SNVs$aa_seq) & (heatmap_fdr$ID != "syn"), "SNVs"]<-"|"

levels = c("NS- 1%", "NS- 5%", "NS- 10%", "NS- 25%", "WT-like","NS+ 25%", "NS+ 10%", "NS+ 5%", "NS+ 1%")

colors<-c( "NS- 1%"="#CD3333", 
           "NS- 5%"="#D66262", 
           "NS- 10%"="#DF9292", 
           "NS- 25%"="#E8C2C2",
           "WT-like"="#F2F2F2", 
           "NS+ 25%"="#B5B5D8",
           "NS+ 10%"="#7979BE",
           "NS+ 5%"="#3C3CA4",
           "NS+ 1%"="#00008B")

p_heatmap_fdr<-ggplot(heatmap_fdr)+
  geom_tile(aes(Pos,factor(Mut, levels=rev(vectorAA)),fill=factor(category, levels=levels)), size=1, color='white')+
  geom_text(aes(Pos, Mut, label=WT), size=3)+
  #geom_text(aes(Pos, Mut, label=SNVs), color="black", size=3, nudge_x=0.32, nudge_y=0.37, angle=45, fontface="bold")+
  scale_fill_manual("Category (FDR)", values= colors, labels=levels)+
  scale_x_continuous(breaks=seq(523,556), labels = peptide_seq_pos, expand=c(0,0))+
  theme_minimal()+
  labs(x="", y="Mutant amino acid", fill="Nucleation score")+
  theme(legend.title = element_text(size=14, face = "bold"),
        legend.text = element_text(size=14), 
        axis.title = element_text(size = 16),
        axis.text = element_text(size=14))
p_heatmap_fdr

ggsave(p_heatmap_fdr,path=path, file="p_heatmap_FDR.jpg", width=15, height=8)

#####
fdr_categories<-singles_stops[,c("Pos", "WT_AA", "Mut", "ID", "nscore_c", "p.adjust", "category_10")]

fdr_categories$category<-"WT-like"
fdr_categories[(fdr_categories$p.adjust<0.25 & fdr_categories$nscore_c<0),]$category<- "NS- 25%"
fdr_categories[(fdr_categories$p.adjust<0.1 & fdr_categories$nscore_c<0),]$category<- "NS- 10%"
fdr_categories[(fdr_categories$p.adjust<0.05 & fdr_categories$nscore_c<0),]$category<- "NS- 5%"
fdr_categories[(fdr_categories$p.adjust<0.01 & fdr_categories$nscore_c<0),]$category<- "NS- 1%"

fdr_categories[(fdr_categories$p.adjust<0.25 & fdr_categories$nscore_c>0),]$category<- "NS+ 25%"
fdr_categories[(fdr_categories$p.adjust<0.1 & fdr_categories$nscore_c>0),]$category<- "NS+ 10%"
fdr_categories[(fdr_categories$p.adjust<0.05 & fdr_categories$nscore_c>0),]$category<- "NS+ 5%"
fdr_categories[(fdr_categories$p.adjust<0.01 & fdr_categories$nscore_c>0),]$category<- "NS+ 1%"


categories_fdr <- fdr_categories %>% group_by(Pos, category) %>% dplyr::summarise(Freq=n()) 

levels_cat_fdr = c("NS- 1%", "NS- 5%", "NS- 10%", "NS- 25%", "WT-like","NS+ 25%", "NS+ 10%", "NS+ 5%", "NS+ 1%")
colors_cat_fdr<-c( "NS- 1%"="#CD3333", 
                   "NS- 5%"="#D66262", 
                   "NS- 10%"="#DF9292", 
                   "NS- 25%"="#E8C2C2",
                   "WT-like"="#F2F2F2", 
                   "NS+ 25%"="#B5B5D8",
                   "NS+ 10%"="#7979BE",
                   "NS+ 5%"="#3C3CA4",
                   "NS+ 1%"="#00008B")

p_categories_fdr<-ggplot(categories_fdr, aes(fill=factor(category, levels=levels_cat_fdr), x=factor(Pos), y=Freq))+
  scale_x_discrete(breaks=seq(523,556), labels = peptide_seq_pos, expand=c(0,0))+
  theme_bw()+
  theme(legend.title=element_text(size=14, face = "bold"),
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color="black", size=0.25),
        axis.title = element_text(size=16),
        axis.text = element_text(size=14),
        plot.margin = unit(c(0,0,0,0), 'cm'))+
  geom_bar(position="stack", stat="identity", alpha=0.8, width=0.8)+
  scale_fill_manual("Nucleation FDR(%)", values=colors_cat_fdr)+
  labs(x= "hRIPK1 peptide", y="Counts")

p_categories_fdr

ggsave(p_categories_fdr,path=path, file="p_categories_FDR.jpg", width=15, height=8)

####
p_hfdrc<-ggarrange(p_heatmap_fdr + rremove("x.text"), p_categories_fdr, 
                   nrow=2, align = "v",  legend="right", heights = c(1, 0.4),
                   common.legend = TRUE)

p_hfdrc

ggsave(p_hfdrc, path=path, file="p_heatmap_fdr_categories.jpg", width=15, height=10)
ggsave(p_hfdrc, path=path, file="p_heatmap_fdr_categories.pdf", width=15, height=10)

#####


#### stacked bars - all variants
# Considering stops
fdr_categories_all<-singles_stops[,c("Pos", "WT_AA", "Mut", "ID", "nscore_c", "p.adjust", "category_10", "low_sigma")]

fdr_categories_all<-fdr_categories_all[fdr_categories_all$low_sigma == T | fdr_categories_all$category_10 != "WT-like",]

fdr_categories_all$category<-"WT-like"
#fdr_categories_all[fdr_categories_all$low_sigma == F,]$category<- "unknown"
fdr_categories_all[(fdr_categories_all$p.adjust<0.1 & fdr_categories_all$nscore_c<0),]$category<- "NS-"
fdr_categories_all[(fdr_categories_all$p.adjust<0.1 & fdr_categories_all$nscore_c>0),]$category<- "NS+"

categories_p_all <- fdr_categories_all %>% group_by(category) %>% dplyr::summarise(Freq_p=n()) 
categories_p_all$library<-"hRIPK1"
categories_p_all$percentage<-round((categories_p_all$Freq_p/sum(categories_p_all$Freq_p))*100, 2)

levels_cat = c("NS-", "WT-like", "NS+", "unknown")
colors_cat<-c("#DF9292", "#F2F2F2", "#7979BE", "white")

p_categories_p_all<-ggplot(categories_p_all, aes(fill=factor(category, levels=levels_cat), x=library, y=percentage))+
  geom_bar(position="stack", stat="identity", alpha=0.8, width=0.8)+
  scale_fill_manual("Nucleation (FDR=0.1)", values=colors_cat)+
  labs(x= "", y="% of variants")+
  annotate("text", x="hRIPK1", y=10, label="NS+", size=5)+
  annotate("text", x="hRIPK1", y=35, label="WT-like", size=5)+
  annotate("text", x="hRIPK1", y=75, label="NS-", size=5)+
  theme_bw()+
  theme(legend.title=element_text(size=14, face = "bold"),
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color='black', size=0.25),
        axis.title = element_text(size=16),
        axis.text = element_text(size=14),
        plot.margin = unit(c(0,0,0,0), 'cm'))

p_categories_p_all

ggsave(p_categories_p_all, path=path, file="p_categories_p_all.jpg", width=4, height=4)
ggsave(p_categories_p_all, path=path, file="p_categories_p_all.pdf", width=4, height=4)





