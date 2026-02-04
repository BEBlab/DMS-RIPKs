library(tidyverse)
library(ggpubr)
require("ggrepel")
library(stringr)

# your folder path
setwd("")

# Create a folder for the figures
dir.create("01_Correlation_Twist_IDT")
path="01_Correlation_Twist_IDT"

# Let's see how the results obtained with IDT (singles) and Twist (indels + singles) libraries correlate

load("hRIPK3_indels_singles.RData")
singles_indels_df<-rename(singles_indels_df, nscore_c_twist = nscore_c)

load("nscore_df_hRIPK3.RData")
singles_stops<-rename(singles_stops, nscore_c_idt = nscore_c)

#
singles_ripk3<-inner_join(singles_stops, singles_indels_df, by="aa_seq")

correlation_all<-cor.test(singles_ripk3$nscore_c_twist, singles_ripk3$nscore_c_idt, use="complete.obs")
correlation_all

y <- singles_ripk3$nscore_c_idt
x <- singles_ripk3$nscore_c_twist

# Fit the regression # then i can predict the values of x (as y)
LR<-lm(y ~ x)

singles_ripk3$residual_abs<-as.numeric(abs(LR$residuals))

SD2<-2*sd(resid(LR))

singles_ripk3$outliers<-FALSE
singles_ripk3[singles_ripk3$residual_abs>SD2,]$outliers<-TRUE


p_corr_all<-ggplot(singles_ripk3, aes(y=nscore_c_idt, x=nscore_c_twist))+
  geom_hline(yintercept = 0, color="black", linewidth=.2)+
  geom_vline(xintercept = 0, color="black", linewidth=.2)+
  geom_point(aes(color=outliers), alpha=.6)+
  scale_color_manual(values=c("black", "red"))+
  geom_text_repel(data=singles_ripk3[singles_ripk3$outliers == T,], aes(label=ID.x), size=3, min.segment.length=0.1)+
  theme_classic()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  annotate("text", label=paste0("R=", round(correlation_all$estimate,2)), x=-Inf, y=Inf, hjust=-0.5, vjust=2, size=4, color="black")+
  annotate("text", label=paste0("p=", format(correlation_all$p.value, digits = 2, scientific = T)), x=-Inf, y=Inf, hjust=-0.2, vjust=3, size=4, color="black")+
  geom_smooth(method = "lm", linetype = 2, size=1, se = F, color = "grey")+
  labs(y="Nucleation Score IDT", x="Nucleation Score Twist")

p_corr_all

ggsave(p_corr_all, file="Corr_ripk3_twist_idt_all.jpg", width = 5, height = 4, path=path)

# There's good correlation between libraries!

# Transform the values
insertions.df$nscore_c<-predict(LR, newdata = data.frame(x = insertions.df$nscore_c))
insertions.df$nscore1_c<-predict(LR, newdata = data.frame(x = insertions.df$nscore1_c))
insertions.df$nscore2_c<-predict(LR, newdata = data.frame(x = insertions.df$nscore2_c))
insertions.df$nscore3_c<-predict(LR, newdata = data.frame(x = insertions.df$nscore3_c))
# Error propagation
insertions.df$sigma<-insertions.df$sigma*LR$coefficients[2]
insertions.df<-rename(insertions.df, sig_10 = sig_fdr, category_10 = category_fdr, category_10_sigma = category_fdr_sigma)

deletions.df$nscore_c<-predict(LR, newdata = data.frame(x = deletions.df$nscore_c))
deletions.df$nscore1_c<-predict(LR, newdata = data.frame(x = deletions.df$nscore1_c))
deletions.df$nscore2_c<-predict(LR, newdata = data.frame(x = deletions.df$nscore2_c))
deletions.df$nscore3_c<-predict(LR, newdata = data.frame(x = deletions.df$nscore3_c))
deletions.df$sigma<-deletions.df$sigma*LR$coefficients[2]
deletions.df<-rename(deletions.df, sig_10 = sig_fdr, category_10 = category_fdr, category_10_sigma = category_fdr_sigma)

singles.df$nscore_c<-predict(LR, newdata = data.frame(x = singles.df$nscore_c))
singles.df$nscore1_c<-predict(LR, newdata = data.frame(x = singles.df$nscore1_c))
singles.df$nscore2_c<-predict(LR, newdata = data.frame(x = singles.df$nscore2_c))
singles.df$nscore3_c<-predict(LR, newdata = data.frame(x = singles.df$nscore3_c))
singles.df$sigma<-singles.df$sigma*LR$coefficients[2]
singles.df$ID<-paste(singles.df$WT_AA, singles.df$Pos, singles.df$Mut, sep = "-")
singles.df<-rename(singles.df, sig_10 = sig_fdr, category_10 = category_fdr, category_10_sigma = category_fdr_sigma)
singles.df<-singles.df[c("nt_seq", "aa_seq", "WT_AA", "Pos", "Mut", "ID", "dataset", "mean_count",
                         "nscore_c", "nscore1_c", "nscore2_c", "nscore3_c", "sigma", "zscore", 
                         "p.adjust", "sig_10", "category_10", "sigma_norm_iqr", "sigma_norm_first_toWT", 
                         "low_sigma", "category_10_sigma")]
singles<-singles.df
#
stops<-singles_stops[singles_stops$STOP == T,]
stops$nscore_c<-stops$nscore_c_idt
stops$dataset<-"STOP"
stops<-stops[c("nt_seq", "aa_seq", "WT_AA", "Pos", "Mut", "ID", "dataset", "mean_count",
               "nscore_c", "nscore1_c", "nscore2_c", "nscore3_c", "sigma", "zscore", 
               "p.adjust", "sig_10", "category_10", "sigma_norm_iqr", "sigma_norm_first_toWT", 
               "low_sigma", "category_10_sigma")]
singles_stops<-rbind(singles, stops)

############
# Singles
####################################################
# FDR=0.01 correction and assignment into categories

singles_stops$zscore<-singles_stops$nscore_c/singles_stops$sigma
singles_stops$p.adjust<-p.adjust(2*pnorm(-abs(singles_stops$zscore)), method = "BH")

singles_stops$sig_10<-FALSE
singles_stops[singles_stops$p.adjust<0.1,]$sig_10<-TRUE

singles_stops$category_10<-"WT-like"
singles_stops[singles_stops$sig_10==T & singles_stops$nscore_c<0,]$category_10<-"NS_dec"
singles_stops[singles_stops$sig_10==T & singles_stops$nscore_c>0,]$category_10<-"NS_inc"

#################################################################################
# interquartile range == fitness range

summary(singles_stops$nscore_c)
iqr<-IQR(singles_stops$nscore_c)

first_to_wt<-summary(singles_stops$nscore_c)[[2]]

p_iqr<-ggplot(singles_stops, aes(x=nscore_c))+
  geom_histogram(color="black", fill="grey", bins=100)+
  theme_classic()+
  labs(x="Nucleation score", title="Fitness range")+
  
  geom_vline(xintercept=0)+
  annotate("text", label="WT fitness= 0", x=0.1, y=110, color="black")+
  
  geom_vline(xintercept=summary(singles_stops$nscore_c)[[2]], color="red")+
  annotate("text", label=paste0("1st Qu= ", round(summary(singles_stops$nscore_c)[[2]], 2)), 
           x=summary(singles_stops$nscore_c)[[2]]-0.8, y=140, color="red")+
  
  geom_vline(xintercept=summary(singles_stops$nscore_c)[[5]], color="red")+
  annotate("text", label=paste0("3rd Qu= ", round(summary(singles_stops$nscore_c)[[5]], 2)), 
           x=summary(singles_stops$nscore_c)[[5]]+0.8, y=140, color="red")+
  
  geom_vline(xintercept=summary(singles_stops$nscore_c)[[3]], color="blue")+
  annotate("text", label=paste0("median= ",round(summary(singles_stops$nscore_c)[[3]], 2) ), 
           x=summary(singles_stops$nscore_c)[[3]]+0.1, y=100, color="blue")
p_iqr

ggsave(p_iqr, file="p_iqr.jpg", path=path, width = 4, height = 6)


#sigma distribution
p_sigma_dist<-ggplot(singles_stops, aes(x=sigma))+
  geom_histogram(color="black", fill="grey", bins=100)+
  theme_classic()+
  scale_x_continuous(limits = c(0,3))
p_sigma_dist

ggsave(p_sigma_dist, file="p_sigma_dist.jpg", path=path, width = 5, height = 3)

# normalise sigmas to fitness range - if fitness range is IQR

singles_stops$sigma_norm_iqr<-""
fitness_range_iqr=abs(IQR(singles_stops$nscore_c))

# or if fitness range is 1rst to WT fitness

singles_stops$sigma_norm_first_toWT<-""
fitness_range_first_toWT=abs(summary(singles_stops$nscore_c)[[2]])

singles_stops$sigma_norm_iqr<-singles_stops$sigma / fitness_range_iqr
singles_stops$sigma_norm_first_toWT<-singles_stops$sigma / fitness_range_first_toWT

singles_stops$sigma_norm_iqr<-as.numeric(singles_stops$sigma_norm_iqr)
singles_stops$sigma_norm_first_toWT<-as.numeric(singles_stops$sigma_norm_first_toWT)

#sigma_norm_iqr

p1<-ggplot(singles_stops, aes(x=sigma_norm_iqr))+
  geom_histogram(bins=200, aes(fill=factor(category_10, levels=c("NS_inc", "NS_dec", "WT-like"))))+
  scale_fill_manual(values=c("#7979BE", "#DF9292", "#F2F2F2"))+
  facet_wrap(~category_10, ncol=1)+
  labs(title="sigma normalised to fitness range (1st Qu to 3r Qu)", fill="Category_FDR")+
  scale_x_continuous(limits = c(0,1))+
  geom_vline(xintercept = 0.3, linetype="dashed", color="grey20")+
  theme_classic()
p1

ggsave(p1, file="p_sigma_normalised.jpg", path=path, width = 6, height = 6)

## exclude those that are over 30% of the fitness range
singles_stops$low_sigma<-FALSE
singles_stops[singles_stops$sigma_norm_iqr<=0.3,]$low_sigma<-TRUE

# Categories based on sigma normalized to fitness range
singles_stops$category_10_sigma<-"unclassified"

singles_stops[singles_stops$category_10=="NS_inc"& singles_stops$low_sigma==T,]$category_10_sigma<-"NS_inc"
singles_stops[singles_stops$category_10=="NS_dec" & singles_stops$low_sigma==T,]$category_10_sigma<-"NS_dec"
singles_stops[singles_stops$category_10=="WT-like" & singles_stops$low_sigma==T,]$category_10_sigma<-"WT-like"
#

############
# Insertions
####################################################
# FDR=0.01 correction and assignment into categories

insertions.df$zscore<-insertions.df$nscore_c/insertions.df$sigma
insertions.df$p.adjust<-p.adjust(2*pnorm(-abs(insertions.df$zscore)), method = "BH")

insertions.df$sig_10<-FALSE
insertions.df[insertions.df$p.adjust<0.1,]$sig_10<-TRUE

insertions.df$category_10<-"WT-like"
insertions.df[insertions.df$sig_10==T & insertions.df$nscore_c<0,]$category_10<-"NS_dec"
insertions.df[insertions.df$sig_10==T & insertions.df$nscore_c>0,]$category_10<-"NS_inc"

#################################################################################
# interquartile range == fitness range

summary(insertions.df$nscore_c)
iqr<-IQR(insertions.df$nscore_c)

# normalise sigmas to fitness range - if fitness range is IQR

insertions.df$sigma_norm_iqr<-""
fitness_range_iqr=abs(IQR(insertions.df$nscore_c))

# or if fitness range is 1rst to WT fitness

insertions.df$sigma_norm_first_toWT<-""
fitness_range_first_toWT=abs(summary(insertions.df$nscore_c)[[2]])

insertions.df$sigma_norm_iqr<-insertions.df$sigma / fitness_range_iqr
insertions.df$sigma_norm_first_toWT<-insertions.df$sigma / fitness_range_first_toWT

insertions.df$sigma_norm_iqr<-as.numeric(insertions.df$sigma_norm_iqr)
insertions.df$sigma_norm_first_toWT<-as.numeric(insertions.df$sigma_norm_first_toWT)

## exclude those that are over 30% of the fitness range
insertions.df$low_sigma<-FALSE
insertions.df[insertions.df$sigma_norm_iqr<=0.3,]$low_sigma<-TRUE

# Categories based on sigma normalized to fitness range
insertions.df$category_10_sigma<-"unclassified"

insertions.df[insertions.df$category_10=="NS_inc"& insertions.df$low_sigma==T,]$category_10_sigma<-"NS_inc"
insertions.df[insertions.df$category_10=="NS_dec" & insertions.df$low_sigma==T,]$category_10_sigma<-"NS_dec"
insertions.df[insertions.df$category_10=="WT-like" & insertions.df$low_sigma==T,]$category_10_sigma<-"WT-like"

############
# Deletions
####################################################
# FDR=0.01 correction and assignment into categories

deletions.df$zscore<-deletions.df$nscore_c/deletions.df$sigma
deletions.df$p.adjust<-p.adjust(2*pnorm(-abs(deletions.df$zscore)), method = "BH")

deletions.df$sig_10<-FALSE
deletions.df[deletions.df$p.adjust<0.1,]$sig_10<-TRUE

deletions.df$category_10<-"WT-like"
deletions.df[deletions.df$sig_10==T & deletions.df$nscore_c<0,]$category_10<-"NS_dec"
deletions.df[deletions.df$sig_10==T & deletions.df$nscore_c>0,]$category_10<-"NS_inc"

#################################################################################
# interquartile range == fitness range

summary(deletions.df$nscore_c)
iqr<-IQR(deletions.df$nscore_c)

# normalise sigmas to fitness range - if fitness range is IQR

deletions.df$sigma_norm_iqr<-""
fitness_range_iqr=abs(IQR(deletions.df$nscore_c))

# or if fitness range is 1rst to WT fitness

deletions.df$sigma_norm_first_toWT<-""
fitness_range_first_toWT=abs(summary(deletions.df$nscore_c)[[2]])

deletions.df$sigma_norm_iqr<-deletions.df$sigma / fitness_range_iqr
deletions.df$sigma_norm_first_toWT<-deletions.df$sigma / fitness_range_first_toWT

deletions.df$sigma_norm_iqr<-as.numeric(deletions.df$sigma_norm_iqr)
deletions.df$sigma_norm_first_toWT<-as.numeric(deletions.df$sigma_norm_first_toWT)

## exclude those that are over 30% of the fitness range
deletions.df$low_sigma<-FALSE
deletions.df[deletions.df$sigma_norm_iqr<=0.3,]$low_sigma<-TRUE

# Categories based on sigma normalized to fitness range
deletions.df$category_10_sigma<-"unclassified"

deletions.df[deletions.df$category_10=="NS_inc"& deletions.df$low_sigma==T,]$category_10_sigma<-"NS_inc"
deletions.df[deletions.df$category_10=="NS_dec" & deletions.df$low_sigma==T,]$category_10_sigma<-"NS_dec"
deletions.df[deletions.df$category_10=="WT-like" & deletions.df$low_sigma==T,]$category_10_sigma<-"WT-like"

################################################################################

save(silent, singles, singles_stops, insertions.df, deletions.df, file="nscore_df_hRIPK3.RData")

################################################################################
library("writexl")
write_xlsx(singles_stops, "hRIPK3_singles.xlsx")
write_xlsx(silent, "hRIPK3_synonymous.xlsx")
write_xlsx(insertions.df, "hRIPK3_insertions.xlsx")
write_xlsx(deletions.df, "hRIPK3_deletions.xlsx")

write.csv(singles_stops, "hRIPK3_singles.csv", row.names=FALSE)
write.csv(silent, "hRIPK3_synonymous.csv", row.names=FALSE)
write.csv(insertions.df, "hRIPK3_insertions.csv", row.names=FALSE)
write.csv(deletions.df, "hRIPK3_deletions.csv", row.names=FALSE)
################################################################################