library(tidyverse)
require(reshape2)
library(ggpubr)

# your folder path
setwd("")

# Create a folder for the figures
dir.create("09_VariantsEffectsIndels")
path="09_VariantsEffectsIndels"

# open processed data
load("nscore_df_hRIPK3.RData")

# start and end positions
start_p<-442
end_p<-475

# mutagenized sequence
peptide_seq<-'NPVTGRPLVNIYNCSGVQVGDNNYLTMQQTTALP'
peptide_seq<-c(strsplit(peptide_seq, '')[[1]])

peptide_seq_pos<-c()
for (n in seq_along(peptide_seq)){
  peptide_seq_pos<-c(peptide_seq_pos, paste0(peptide_seq[n], '\n', n+start_p-1))
}

vectorAA <- c("G","A","V","L","M","I","F","Y","W","K","R","D","E","S","T","C","N","Q","H", "P")

min<-min(insertions.df$nscore_c)
max<-max(insertions.df$nscore_c)
cols <- c(colorRampPalette(c( "brown3", "grey95"))((-min/(-min+max)*100)-0.5), colorRampPalette("grey95")(1),
          colorRampPalette(c("grey95",  "darkblue"), bias=1)((max/(-min+max)*100)-0.5))
#####
# single AA insertions

insertions<-insertions.df[,c("aa_seq", "ins_pos", "ins_aa", "nscore_c", "sigma", "p.adjust")]

heatmap_insertions<-insertions

heatmap_insertions$ID<-paste0("Ins_k1_", heatmap_insertions$ins_pos, "_", heatmap_insertions$ins_aa)


p_heatmap_inser<-ggplot(heatmap_insertions)+
  geom_tile(aes(factor(ins_pos, levels=c(start_p:(end_p-1)), labels=c(start_p:(end_p-1))),
                factor(ins_aa, levels=rev(vectorAA)), fill=nscore_c), color="white", size=1)+
  theme_minimal()+
  labs(x="hRIPK3", y="Mutant amino acid", fill="Nucleation score")+
  theme(legend.title = element_text(size=14),
        legend.text = element_text(size=14), 
        axis.title = element_text(size = 16),
        axis.text = element_text(size=14))+
  scale_fill_gradientn(colours=cols, limits=c(min,max), na.value = "grey60", breaks=c(-3, -1, 1))
p_heatmap_inser

ggsave(p_heatmap_inser, file="p_heatmap_nscore_insertions.jpg", width=15, height=8, path=path)

##### Violin plot per position
subsp_median_df<-as.data.frame(heatmap_insertions %>% group_by(ins_pos) %>% dplyr::summarise(median_p=median(nscore_c)))
heatmap_insertions<-left_join(heatmap_insertions, subsp_median_df)

p_violin_p<-ggplot(heatmap_insertions, aes(x=factor(ins_pos, levels=c(start_p:(end_p-1)), labels=c(start_p:(end_p-1))), y=nscore_c))+
  geom_hline(yintercept = 0, size=0.1)+
  geom_violin(scale = "width", size=0.2, aes(fill=median_p))+
  geom_boxplot(width=0.15, outlier.shape = NA, size=0.2)+
  scale_x_discrete(breaks=seq(442,475), labels = peptide_seq_pos, expand=c(0,0))+
  theme_bw()+
  labs(x="hRIPK3", y="Nucleation score", fill="Median NS")+
  theme(legend.title = element_text(size=14),
        legend.text = element_text(size=14), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        plot.margin = unit(c(0,0,0,0), 'cm'))+
  scale_fill_gradientn(colours=cols, limits=c(min,max), na.value = "grey60", breaks=c(-3, -1, 1))

p_violin_p

ggsave(p_violin_p, path=path, file="p_violinplotp_nscore_insertions.jpg", width=15, height=8)

##### violin per change

subsm_median_df<-as.data.frame(heatmap_insertions %>% group_by(ins_aa) %>% dplyr::summarise(median_m=median(nscore_c)))
heatmap_insertions<-left_join(heatmap_insertions, subsm_median_df)

p_violin_m<-ggplot(heatmap_insertions, aes(x=factor(ins_aa, levels=rev(vectorAA)), y=nscore_c))+
  geom_hline(yintercept = 0, size=0.1)+
  geom_violin(scale = "width", size=0.2, aes(fill=median_m))+
  geom_boxplot(width=0.15, outlier.shape = NA, size=0.2)+
  theme_bw()+
  labs(x="hRIPK3", y="Nucleation score", fill="Median NS")+
  theme(legend.title = element_text(size=14),
        legend.text = element_text(size=14), 
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=16),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        plot.margin = unit(c(0,0,0,0), 'cm'))+
  scale_fill_gradientn(colours=cols, limits=c(min,max), na.value = "grey60", breaks=c(-3, -1, 1))+ 
  coord_flip() 

p_violin_m

ggsave(p_violin_m, path=path, file="p_violinplot_inset.jpg", width=4, height=8)


#### stacked bars
# Considering stops
fdr_categories<-insertions.df[,c("ins_pos", "ins_aa", "nscore_c", "p.adjust", "category_10")]

fdr_categories$category<-"WT-like"
fdr_categories[(fdr_categories$p.adjust<0.1 & fdr_categories$nscore_c<0),]$category<- "NS-"
fdr_categories[(fdr_categories$p.adjust<0.1 & fdr_categories$nscore_c>0),]$category<- "NS+"

categories_p <- fdr_categories %>% group_by(ins_pos,category) %>% dplyr::summarise(Freq_p=n()) 

levels_cat = c("NS-", "WT-like", "NS+")
colors_cat<-c("#DF9292", "#F2F2F2", "#7979BE")

p_categories_p<-ggplot(categories_p, aes(fill=factor(category, levels=levels_cat), x=factor(ins_pos, levels=c(start_p:(end_p-1)), labels=c(start_p:(end_p-1))), y=Freq_p))+
  theme_bw()+
  theme(legend.title=element_text(size=14),
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
  labs(x= "hRIPK3", y="Counts")

p_categories_p

ggsave(p_categories_p, path=path, file="p_categories_nscore_insertions.jpg", width=15, height=8)


######

categories_m <- fdr_categories %>% group_by(ins_aa, category) %>% dplyr::summarise(Freq_m=n()) 


p_categories_m<-ggplot(categories_m, aes(fill=factor(category, levels=levels_cat), y=factor(ins_aa, levels=rev(vectorAA)), x=Freq_m))+
  theme_bw()+
  theme(legend.title=element_text(size=14),
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

ggsave(p_categories_m, path=path, file="p_categoriesm_nscore_insertions.jpg", width=15, height=8)


##################

p_hv <- ggarrange(p_heatmap_inser + rremove("x.text"), p_violin_m + rremove("y.text"), p_categories_m + rremove("y.text"), 
                  ncol=3, align = "h",  legend="none", widths = c(1, 0.2, 0.2), common.legend=TRUE) 

p_hv

ggsave(p_hv, path=path, file="RIPK3_heatmap_violinm_categoreism_ins.jpg", width=15, height=8)
ggsave(p_hv, path=path, file="RIPK3_heatmap_violinm_categoreism_ins.pdf", width=15, height=8)


### Arrange plots
p_heatmap_inser <- p_heatmap_inser + theme(axis.title.x = element_blank())
p_hvc<-ggarrange(p_heatmap_inser + rremove("x.text"), p_violin_p + rremove("x.text"), p_categories_p, 
                 nrow=3, align = "v",  legend="none", heights = c(1, 0.3, 0.3))

p_hvc

ggsave(p_hvc, path=path, file="p_heatmap_violin_categories_nscore_ins.jpg", width=14, height=12)
ggsave(p_hvc, path=path, file="p_heatmap_violin_categories_nscore_ins.pdf", width=14, height=12)

##################################################################################
#### FDR
insertions.df<-insertions.df[insertions.df$low_sigma == T | insertions.df$sig_10 == T, ]
insertions<-insertions.df[,c("ins_pos", "ins_aa", "nscore_c", "sigma", "p.adjust", "low_sigma")]

heatmap_insertions_fdr<-rbind(insertions)

heatmap_insertions_fdr$category<-"WT-like"

heatmap_insertions_fdr$p.adjust<-as.numeric(heatmap_insertions_fdr$p.adjust)

heatmap_insertions_fdr[(heatmap_insertions_fdr$category!="missing" &  heatmap_insertions_fdr$p.adjust<0.25 & heatmap_insertions_fdr$nscore_c<0),]$category<- "NS- 25%"
heatmap_insertions_fdr[(heatmap_insertions_fdr$category!="missing"  & heatmap_insertions_fdr$p.adjust<0.1 & heatmap_insertions_fdr$nscore_c<0),]$category<- "NS- 10%"
heatmap_insertions_fdr[(heatmap_insertions_fdr$category!="missing"  & heatmap_insertions_fdr$p.adjust<0.05 & heatmap_insertions_fdr$nscore_c<0),]$category<- "NS- 5%"
heatmap_insertions_fdr[(heatmap_insertions_fdr$category!="missing"  & heatmap_insertions_fdr$p.adjust<0.01 & heatmap_insertions_fdr$nscore_c<0),]$category<- "NS- 1%"

heatmap_insertions_fdr[(heatmap_insertions_fdr$category!="missing"  & heatmap_insertions_fdr$p.adjust<0.25 & heatmap_insertions_fdr$nscore_c>0),]$category<- "NS+ 25%"
heatmap_insertions_fdr[(heatmap_insertions_fdr$category!="missing"  & heatmap_insertions_fdr$p.adjust<0.1 & heatmap_insertions_fdr$nscore_c>0),]$category<- "NS+ 10%"
heatmap_insertions_fdr[(heatmap_insertions_fdr$category!="missing"  & heatmap_insertions_fdr$p.adjust<0.05 & heatmap_insertions_fdr$nscore_c>0),]$category<- "NS+ 5%"
heatmap_insertions_fdr[(heatmap_insertions_fdr$category!="missing"  & heatmap_insertions_fdr$p.adjust<0.01 & heatmap_insertions_fdr$nscore_c>0),]$category<- "NS+ 1%"

heatmap_insertions_fdr$significance = ""
heatmap_insertions_fdr[(heatmap_insertions_fdr$low_sigma == TRUE) & (heatmap_insertions_fdr$category == "WT-like"),]$significance = "|"

ramp_colors <- c(colorRampPalette(c("#DF9292", "#F2F2F2", "#7979BE"))(9))
colors<-c(ramp_colors,  "grey30" )

levels = c("NS- 1%", "NS- 5%", "NS- 10%", "NS- 25%", "WT-like","NS+ 25%", "NS+ 10%", "NS+ 5%", "NS+ 1%", "missing")

### 
p_fdr_heatmap_ins<-ggplot(heatmap_insertions_fdr)+
  geom_tile(aes(x=factor(ins_pos, levels=c(start_p:(end_p-1)), labels=c(start_p:(end_p-1))),
                y=factor(ins_aa, levels=rev(vectorAA)),fill=factor(category, levels=levels)), 
                alpha=0.8, color='white', lwd = 1)+
  scale_fill_manual("Category (FDR)", values= colors, labels=levels)+
  theme_minimal()+
  labs(x="", y="Mutant amino acid", fill="Nucleation score")+
  theme(legend.title = element_text(size=14, face = "bold"),
        legend.text = element_text(size=14), 
        axis.title = element_text(size = 16),
        axis.text = element_text(size=14))

p_fdr_heatmap_ins

ggsave(p_fdr_heatmap_ins, file="p_heatmap_fdr_insertions.jpg", width=15, height=8, path=path)

#####
categories_fdr <- heatmap_insertions_fdr %>% group_by(ins_pos, category) %>% dplyr::summarise(Freq=n()) 

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

p_categories_fdr<-ggplot(categories_fdr, aes(fill=factor(category, levels=levels_cat_fdr), x=factor(ins_pos, levels=c(start_p:(end_p-1)), labels=c(start_p:(end_p-1))), y=Freq))+
  theme_bw()+
  theme(legend.title=element_text(size=14),
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
  labs(x= "hRIPK3 peptide", y="Counts")

p_categories_fdr

ggsave(p_categories_fdr,path=path, file="p_categories_FDR.jpg", width=15, height=8)


####
p_hfdrc<-ggarrange(p_fdr_heatmap_ins + rremove("x.text"), p_categories_fdr, 
                   nrow=2, align = "v",  legend="right", heights = c(1, 0.4),
                   common.legend = TRUE)

p_hfdrc

ggsave(p_hfdrc, path=path, file="p_heatmap_fdr_categories_insertions.jpg", width=15, height=10)
ggsave(p_hfdrc, path=path, file="p_heatmap_fdr_categories_insertions.pdf", width=15, height=10)


#####################################################################################
#### single AA deletions

p_single_deletion<-ggplot(deletions.df, aes(x=factor(del_pos, levels = c(start_p:end_p)), y=nscore_c ))+
  geom_hline(yintercept = 0, size=0.05)+
  geom_errorbar(aes(ymin=nscore_c-1.96*sigma, ymax=nscore_c+1.96*sigma), width=0, size=0.1)+
  geom_point(data=deletions.df[deletions.df$category_10 %in% c("NS_inc", "NS_dec"),],
             aes(fill=nscore_c), size=4, shape=21, stroke=1.5)+
  geom_point(data=deletions.df[!deletions.df$category_10 %in% c("NS_inc", "NS_dec"),],
             aes(fill=nscore_c), size=4, shape=21, stroke=0.5)+
  scale_x_discrete(breaks=seq(start_p,end_p), labels = peptide_seq_pos)+
  theme_classic()+
  theme(legend.title = element_text(size=14),
        legend.text = element_text(size=14), 
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        plot.margin = unit(c(0,0,0,0), 'cm'))+
  scale_fill_gradientn(colours=cols, limits=c(min,max), breaks=c(-3, -1, 1))+
  labs(y="Nucleation score", x="Deleted amino acid and position", fill="Nucleation score")
  
p_single_deletion

ggsave(p_single_deletion, file="p_deletions_effect.jpg", width = 15, height = 4, path=path)
ggsave(p_single_deletion, file="p_deletions_effect.pdf", width = 15, height = 4, path=path)

###
singles_indels_df<-rbind(singles_stops[,c("aa_seq", "ID", "nscore_c", "sigma", "zscore", "p.adjust", "dataset", "low_sigma", "category_10")],
                         insertions.df[,c("aa_seq", "ID", "nscore_c", "sigma", "zscore", "p.adjust", "dataset", "low_sigma", "category_10")],
                         deletions.df[,c("aa_seq", "ID", "nscore_c", "sigma", "zscore", "p.adjust", "dataset", "low_sigma", "category_10")])

fdr_categories<-singles_indels_df

fdr_categories$category<-"WT-like"
fdr_categories[(fdr_categories$p.adjust<0.25 & fdr_categories$nscore_c<0),]$category<- "NS- 25%"
fdr_categories[(fdr_categories$p.adjust<0.1 & fdr_categories$nscore_c<0),]$category<- "NS- 10%"
fdr_categories[(fdr_categories$p.adjust<0.05 & fdr_categories$nscore_c<0),]$category<- "NS- 5%"
fdr_categories[(fdr_categories$p.adjust<0.01 & fdr_categories$nscore_c<0),]$category<- "NS- 1%"
fdr_categories[(fdr_categories$p.adjust<0.25 & fdr_categories$nscore_c>0),]$category<- "NS+ 25%"
fdr_categories[(fdr_categories$p.adjust<0.1 & fdr_categories$nscore_c>0),]$category<- "NS+ 10%"
fdr_categories[(fdr_categories$p.adjust<0.05 & fdr_categories$nscore_c>0),]$category<- "NS+ 5%"
fdr_categories[(fdr_categories$p.adjust<0.01 & fdr_categories$nscore_c>0),]$category<- "NS+ 1%"

fdr_categories<-fdr_categories[fdr_categories$low_sigma == T | fdr_categories$category != "WT-like",]

# Merge overall categories and splitted by region
categories <- as.data.frame(fdr_categories %>% group_by(dataset, category) %>% dplyr::summarise(n=n())  %>% mutate(freq = n / sum(n) ))
all_variants_categories <- as.data.frame(fdr_categories %>% group_by(category) %>% dplyr::summarise(n=n())  %>% mutate(freq = n / sum(n) ))
all_variants_categories$dataset<-"All"

categories<-rbind(categories, all_variants_categories)
categories<-categories[categories$dataset != "WT",]

colors_fdr <- c(colorRampPalette(c( "#DF9292", "#F2F2F2", "#7979BE"))(9))

levels = c("NS- 1%", "NS- 5%", "NS- 10%", "NS- 25%", "WT-like","NS+ 25%", "NS+ 10%", "NS+ 5%", "NS+ 1%")

p_categories_variants<-ggplot(categories, aes(fill=factor(category, levels=rev(levels)), 
                            y=factor(dataset, levels=c("All", "missense", "insertion", "deletion", "STOP"), labels=c("All", "Missense", "Insertion", "Deletion", "Nonsense")), 
                            x=freq)) + 
  geom_bar(position="stack", stat="identity", alpha=0.8, width=0.8)+
  theme_classic()+
  theme(legend.title=element_text(size=14),
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color='black', size=0.25),
        axis.title = element_text(size=16),
        axis.text = element_text(size=14),
        plot.margin = unit(c(0,0,0,0), 'cm'))+
  scale_fill_manual("FDR category", values=rev(colors_fdr))+
  labs( x= "Frequency",y="")
p_categories_variants

ggsave(p_categories_variants, file="p_categories_all.jpg", width = 6, height = 5, path=path)
ggsave(p_categories_variants, file="p_categories_all.pdf", width = 6, height = 5, path=path)

###
fdr_categories[fdr_categories$dataset %in% c("STOP", "missense"), "dataset"]<-"substitution"
categories <- as.data.frame(fdr_categories %>% group_by(dataset, category_10) %>% dplyr::summarise(n=n())  %>% mutate(freq = n / sum(n) ))
all_variants_categories <- as.data.frame(fdr_categories %>% group_by(category_10) %>% dplyr::summarise(n=n())  %>% mutate(freq = n / sum(n) ))
all_variants_categories$dataset<-"All"

categories<-rbind(categories, all_variants_categories)
categories<-categories[categories$dataset != "WT",]

levels_cat = c("NS_dec", "WT-like", "NS_inc")
colors_cat<-c("#DF9292", "#F2F2F2", "#7979BE")

p_categories_variants<-ggplot(categories, aes(fill=factor(category_10, levels=rev(levels_cat)), 
                                            y=factor(dataset, levels=c("All", "substitution", "insertion", "deletion"), labels=c("All", "Substitutions", "Insertions", "Deletions")), 
                                            x=freq)) + 
  geom_bar(position="stack", stat="identity", alpha=0.8, width=0.8)+
  theme_classic()+
  theme(legend.title=element_text(size=14),
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color='black', size=0.25),
        axis.title = element_text(size=16),
        axis.text = element_text(size=14),
        plot.margin = unit(c(0,0,0,0), 'cm'))+
  scale_fill_manual("FDR category", values=rev(colors_cat))+
  labs( x= "Frequency",y="")
p_categories_variants

ggsave(p_categories_variants, file="p_categories_all_cat10.jpg", width = 6, height = 5, path=path)
ggsave(p_categories_variants, file="p_categories_all_cat10.pdf", width = 6, height = 5, path=path)