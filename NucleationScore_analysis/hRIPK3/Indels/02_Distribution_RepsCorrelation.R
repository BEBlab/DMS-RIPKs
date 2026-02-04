library(tidyverse)
library(ggpubr)

# your folder path
setwd("")

# Create a folder for the figures
dir.create("02_Distribution_RepsCorrelation")
path="02_Distribution_RepsCorrelation"

name<-"hRIPK3_indels_singles"

load(paste0(name, ".RData"))

###

# replicate correlations 
# considering only designed sequences

combinations<-c("1_2", "1_3", "2_3")

my_list<-list()

for(comb in combinations){
  
  first_rep<-unlist(strsplit(comb, "_"))[1]
  second_rep<-unlist(strsplit(comb, "_"))[2]
  
  subset<-singles_indels_df[!is.na(singles_indels_df[[paste0("nscore", first_rep, "_c")]]),]
  subset<-subset[!is.na(subset[[paste0("nscore",second_rep, "_c")]]),]
  
  print(paste0("n(",comb,")= ",length(subset$ID)))
  
  corr<-cor.test(singles_indels_df[[paste0("nscore", first_rep, "_c")]], 
                 singles_indels_df[[paste0("nscore", second_rep, "_c")]], use="complete.obs")
  R<-corr$estimate
  p<-corr$p.value
  
  p_corr<-ggplot(singles_indels_df, aes(x=.data[[paste0("nscore", first_rep, "_c")]], 
                               y=.data[[paste0("nscore", second_rep, "_c")]] ))+
    stat_binhex()+
    theme_bw()+
    labs(x=paste0("Replicate ", first_rep), y=paste0("Replicate ", second_rep))+
    theme(  panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(color='black'),
            axis.title  = element_text(size = 20),
            axis.text = element_text(size=16))+
    annotate("text", x = -4, y = 3, label = paste0("R=", round(R, 2)), size=6)+
    annotate("text", x = -4, y = 2, label = paste0("p=",format(p, digits = 2, scientific = T)), size=8)+
    scale_fill_gradient(high="grey30", low="grey90")
  p_corr
  
  my_list[[comb]]<-p_corr
  
}

p_all<-ggarrange(my_list$`1_2`, my_list$`1_3`, my_list$`2_3`, ncol=3, common.legend = T, legend = "right")
p_all

ggsave(p_all, file="reps_correlation.jpg", width = 12, height = 3.5, path=path)
ggsave(p_all, file="reps_correlation.pdf", width = 12, height = 3.5, path=path)

#################################################################################
# nucleation socre distribution (Synonymous, Missense, Insertions and Deletions)

ns_dist_df<-rbind(synonymous[,c("nt_seq", "nscore_c", "dataset")], singles_indels_df[,c("nt_seq", "nscore_c", "dataset")]) 
ns_dist_df<-ns_dist_df[ns_dist_df$dataset != "WT",]

p_hist<-ggplot(ns_dist_df, aes(x=nscore_c))+
  geom_histogram(binwidth = 0.2, position="identity", alpha=0.5, color=NA, aes(fill=factor(dataset, levels = c("missense", "insertion", "deletion","synonymous"))))+
  geom_vline(aes(xintercept=0), color="black", linetype="dashed", size=0.5)+
  scale_fill_manual(values=c("grey70", "darkorange", "darkgreen", "red"))+
  theme_bw()+
  labs(x="Nucleation score", y="Counts")+
  theme(legend.position = c(0.1,0.8),
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        legend.key.width = unit(0.4,"cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color='black'),
        axis.title  = element_text(size = 14),
        axis.text = element_text(size=14),
        strip.text = element_text(size=14),
        plot.margin = margin(0,0,0,0, unit = 'cm')
        )+
  xlim(-6,3.5)
p_hist

ggsave(p_hist, file="distribution_hist.jpg", path=path, width = 12, height = 8)
ggsave(p_hist, file="distribution_hist.pdf", path=path, width = 12, height = 8)

###

# designed sequences

n_text<-as.data.frame(table(ns_dist_df$dataset))
colnames(n_text)<-c("dataset", "n")

p_dist<-ggplot(ns_dist_df, aes(x=nscore_c))+
  geom_histogram(bins=50, fill="grey70")+
  geom_vline(xintercept = 0,linetype="dashed", size=0.2)+
  facet_wrap(~factor(dataset, levels=c("synonymous", "missense", "insertion", "deletion"), labels=c("Synonymous", "Missense", "Insertion", "Deletion")),
             ncol=4)+
  theme_bw()+
  geom_text(data=n_text[n_text$dataset!="wt",], aes(label=paste0("n=",n), x=-2.5, y=60), size=5, colour="black")+
  theme(panel.grid = element_blank(),
        strip.text = element_text(size=20),
        strip.background = element_blank(),
        panel.border = element_blank(),
        axis.title  = element_text(size = 20),
        axis.text = element_text(size=16),
        axis.line = element_line(color='black'),
        axis.ticks = element_line(size=0.1))+
  labs(x="Nucleation score", y="Counts")
p_dist

ggsave(p_dist, file="distribution_hist_wrap.jpg", path=path, width = 12, height = 4)
ggsave(p_dist, file="distribution_hist_wrap.pdf", path=path, width = 12, height = 4)

