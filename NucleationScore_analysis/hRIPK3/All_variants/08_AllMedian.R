library(tidyverse)
library(ggpubr)
require("ggrepel")
library(stringr)
library(shadowtext)

# your folder path
setwd("")

# Create a folder for the figures
dir.create("08_AllMedian")
path="08_AllMedian"

# Load ripk3 data

load("nscore_df_hRIPK3.RData")

start_p<-442
end_p<-475

# Pre-processing of singles
singles_ripk3<-singles_stops[singles_stops$Mut != "*",]
singles_ripk3$protein<-"ripk3"
singles_ripk3$Pos_c<-singles_ripk3$Pos-(start_p-1) # Substract the initial position. All the libraries will start at 1 for comparisons
singles_ripk3$nscore_c_ripk3<-singles_ripk3$nscore_c
singles_ripk3$wt_ripk3<-singles_ripk3$WT_AA
singles_ripk3$pos_ripk3<-singles_ripk3$Pos
singles_ripk3$mut_ripk3<-singles_ripk3$Mut

singles_ripk3$residue<-paste0(singles_ripk3$WT_AA,"\n " ,singles_ripk3$Pos)
singles_ripk3$change<-paste0(singles_ripk3$Pos, singles_ripk3$Mut)
singles_ripk3$change_type<-paste0(singles_ripk3$Pos, singles_ripk3$type_mut)

singles_ripk3$dataset<-"hRIPK3"

singles_ripk3 <- singles_ripk3 %>% group_by(residue) %>% mutate(median_nscore=median(nscore_c_ripk3))

singles_ripk3 <- singles_ripk3 %>% group_by(residue) %>% mutate(median_zscore=median(zscore))

###
# Load mRIP1 data
load("nscore_df_mRIP1.RData")

start_p_m <- 512
end_p_m <- 545

# Pre-processing of singles
singles_mRIP1<-singles_stops[singles_stops$Mut != "*",]
singles_mRIP1$protein<-"mRIP1"
singles_mRIP1$Pos_c<-singles_mRIP1$Pos-(start_p_m-1) # Substract the initial position. All the libraries will start at 1 for comparisons
singles_mRIP1$nscore_c_mRIP1<-singles_mRIP1$nscore_c
singles_mRIP1$wt_mRIP1<-singles_mRIP1$WT_AA
singles_mRIP1$pos_mRIP1<-singles_mRIP1$Pos
singles_mRIP1$mut_mRIP1<-singles_mRIP1$Mut

singles_mRIP1$residue<-paste0(singles_mRIP1$WT_AA, "\n ", singles_mRIP1$Pos)
singles_mRIP1$change<-paste0(singles_mRIP1$Pos, singles_mRIP1$Mut)

singles_mRIP1$dataset<-"mRIPK1"

singles_mRIP1 <- singles_mRIP1 %>% group_by(residue) %>% mutate(median_nscore=median(nscore_c_mRIP1))

singles_mRIP1 <- singles_mRIP1 %>% group_by(residue) %>% mutate(median_zscore=median(zscore))

###
# Load mRIP3 data
load("nscore_df_mRIP3.RData")

start_p_m <- 432
end_p_m <- 465

# Pre-processing of singles
singles_mRIP3<-singles_stops[singles_stops$Mut != "*",]
singles_mRIP3$protein<-"mRIP3"
singles_mRIP3$Pos_c<-singles_mRIP3$Pos-(start_p_m-1) # Substract the initial position. All the libraries will start at 1 for comparisons
singles_mRIP3$nscore_c_mRIP3<-singles_mRIP3$nscore_c
singles_mRIP3$wt_mRIP3<-singles_mRIP3$WT_AA
singles_mRIP3$pos_mRIP3<-singles_mRIP3$Pos
singles_mRIP3$mut_mRIP3<-singles_mRIP3$Mut

singles_mRIP3$residue<-paste0(singles_mRIP3$WT_AA, "\n ", singles_mRIP3$Pos)
singles_mRIP3$change<-paste0(singles_mRIP3$Pos, singles_mRIP3$Mut)

singles_mRIP3$dataset<-"mRIPK3"

singles_mRIP3 <- singles_mRIP3 %>% group_by(residue) %>% mutate(median_nscore=median(nscore_c_mRIP3))

singles_mRIP3 <- singles_mRIP3 %>% group_by(residue) %>% mutate(median_zscore=median(zscore))

###
# Load hripk1 data
load("nscore_df_hRIPK1.RData")

start_p_m <- 523
end_p_m <- 556

# Pre-processing of singles
singles_ripk1<-singles_stops[singles_stops$Mut != "*",]
singles_ripk1$protein<-"ripk1"
singles_ripk1$Pos_c<-singles_ripk1$Pos-(start_p_m-1) # Substract the initial position. All the libraries will start at 1 for comparisons
singles_ripk1$nscore_c_ripk1<-singles_ripk1$nscore_c
singles_ripk1$wt_ripk1<-singles_ripk1$WT_AA
singles_ripk1$pos_ripk1<-singles_ripk1$Pos
singles_ripk1$mut_ripk1<-singles_ripk1$Mut

singles_ripk1$residue<-paste0(singles_ripk1$WT_AA, "\n ", singles_ripk1$Pos)
singles_ripk1$change<-paste0(singles_ripk1$Pos, singles_ripk1$Mut)

singles_ripk1$dataset<-"hRIPK1"

singles_ripk1 <- singles_ripk1 %>% group_by(residue) %>% mutate(median_nscore=median(nscore_c_ripk1))

singles_ripk1 <- singles_ripk1 %>% group_by(residue) %>% mutate(median_zscore=median(zscore))

##
cols_interest<-c("dataset", "residue", "ID", "WT_AA", "Mut", "Pos_c", "nscore_c", "median_nscore", "median_zscore", "category_10", "low_sigma")

# make a single df with all the data
singles_ripks<-rbind(singles_ripk3[cols_interest], 
                     singles_ripk1[cols_interest],
                     singles_mRIP3[cols_interest], 
                     singles_mRIP1[cols_interest])

min<-min(singles_ripks$median_nscore)
max<-max(singles_ripks$median_nscore)
cols <- c(colorRampPalette(c( "brown3", "grey95"))((-min/(-min+max)*100)-0.5), colorRampPalette("grey95")(1),
          colorRampPalette(c("grey95",  "darkblue"), bias=1)((max/(-min+max)*100)-0.5))

singles_ripks <- singles_ripks[!duplicated(singles_ripks),]

align_p<-ggplot(singles_ripks, aes(x=Pos_c, y=factor(dataset, levels=c("mRIPK3", "mRIPK1", "hRIPK3", "hRIPK1")), 
                                        fill=median_nscore, label=WT_AA))+
  geom_tile(color='white', size=1)+
  geom_shadowtext(size=3, bg.color="white", color="grey20")+
  scale_fill_gradientn(colours=cols, limits=c(min,max), na.value = "grey60")+
  labs(x="", y="")+
  theme_void()+
  theme(legend.title = element_text(size=14),
        legend.text = element_text(size=14), 
        axis.title = element_text(size = 16),
        axis.text = element_text(size=14),
        axis.text.x = element_blank())
align_p

ggsave(align_p, path=path, file="p_heatmap_median_nscore_all.jpg", width=15, height=4)
ggsave(align_p, path=path, file="p_heatmap_median_nscore_all.pdf", width=15, height=4)
