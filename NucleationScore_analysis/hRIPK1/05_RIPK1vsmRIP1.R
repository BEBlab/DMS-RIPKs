library(tidyverse)
library(ggpubr)
require("ggrepel")
library(stringr)

# your folder path
setwd("")

# Create a folder for the figures
dir.create("05_hRIPK1VsmRIP1")
path="05_hRIPK1VsmRIP1"

# all the amino acid
all_aa<-c("G","A","V","L","M","I","F","Y","W","K","R","D","E","S","T","C","N","Q","H","P")

AA_type<-data.frame("wt_aa"= all_aa,
                    "mut_aa" = all_aa,
                    "name_AA"=c("Glycine", "Alanine","Valine","Leucine","Methionine","Isoleucine","Phenylalanine",
                                "Tyrosine","Tryptophan","Lysine","Arginine","Aspartic acid","Glutamic acid","Serine","Threonine",
                                "Cysteine","Asparagine","Glutamine","Histidine","Proline"),
                    "type"=c("glycine",rep("aliphatic",5),rep("aromatic",3),rep("positive",2),rep("negative",2),rep("polar",6),"proline"))

color_axis_x=c("N442"="black", "P443"="black", "V444"="black", "T445"="black", "G446"="black", "R447"="black",
               "P448"="black", "L449"="black", "V450"="red", "N451"="red", "I452"="red", "Y453"="red",
               "N454"="red", "C455"="red", "S456"="red", "G457"="red", "V458"="red", "Q459"="red",
               "V460"="red", "G461"="red", "D462"="red", "N463"="red", "N464"="red", "Y465"= "red",
               "L466"="red", "T467"="black", "M468"="black", "Q469"="black","Q470"="black", "T471"="black",
               "T472"="black", "A473"="black", "L468"="black", "P469"="black")

###
# Load RIPK1 data
load("nscore_df_hRIPK1.RData")

start_p<-523
end_p<-556

# Pre-processing of singles
singles_ripk1<-singles_stops[singles_stops$Mut != "*",]
singles_ripk1$protein<-"RIPK1"
singles_ripk1$Pos_c<-singles_ripk1$Pos-(start_p-1) # Substract the initial position. All the libraries will start at 1 for comparisons
singles_ripk1$nscore_c_ripk1<-singles_ripk1$nscore_c
singles_ripk1$wt_ripk1<-singles_ripk1$WT_AA
singles_ripk1$pos_ripk1<-singles_ripk1$Pos
singles_ripk1$mut_ripk1<-singles_ripk1$Mut

singles_ripk1<-left_join(singles_ripk1, AA_type[,c("wt_aa", "type")], by=c("wt_ripk1"="wt_aa"))
singles_ripk1<-rename(singles_ripk1, type_wt = type)
singles_ripk1<-left_join(singles_ripk1, AA_type[,c("wt_aa", "type")], by=c("mut_ripk1"="wt_aa"))
singles_ripk1<-rename(singles_ripk1, type_mut = type)

singles_ripk1$residue_ripk1<-paste0(singles_ripk1$WT_AA, singles_ripk1$Pos)
singles_ripk1$change<-paste0(singles_ripk1$Pos, singles_ripk1$Mut)
singles_ripk1$change_type<-paste0(singles_ripk1$Pos, singles_ripk1$type_mut)

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

singles_mRIP1<-left_join(singles_mRIP1, AA_type[,c("wt_aa", "type")], by=c("wt_mRIP1"="wt_aa"))
singles_mRIP1<-rename(singles_mRIP1, type_wt = type)
singles_mRIP1<-left_join(singles_mRIP1, AA_type[,c("wt_aa", "type")], by=c("mut_mRIP1"="wt_aa"))
singles_mRIP1<-rename(singles_mRIP1, type_mut = type)

singles_mRIP1$residue_mRIP1<-paste0(singles_mRIP1$WT_AA, singles_mRIP1$Pos)
singles_mRIP1$change<-paste0(singles_mRIP1$Pos, singles_mRIP1$Mut)


################################################################################
# Singles correlation 

singles_ripk1$change<-paste0(singles_ripk1$Pos_c, singles_ripk1$Mut)
singles_ripk1$change_pos<-paste0(singles_ripk1$WT_AA, singles_ripk1$Pos, singles_ripk1$Mut)
singles_ripk1$change_pos_c<-paste0(singles_ripk1$WT_AA, singles_ripk1$Pos_c, singles_ripk1$Mut)
singles_ripk1$change_aa<-paste0(singles_ripk1$wt_ripk1, singles_ripk1$mut_ripk1)
singles_ripk1$change_type<-paste0(singles_ripk1$Pos_c, singles_ripk1$type_mut)
singles_ripk1$change_pos_type<-paste0(singles_ripk1$type_wt, singles_ripk1$Pos_c, singles_ripk1$type_mut)

singles_mRIP1$change<-paste0(singles_mRIP1$Pos_c, singles_mRIP1$Mut)
singles_mRIP1$change_pos<-paste0(singles_mRIP1$WT_AA, singles_mRIP1$Pos, singles_mRIP1$Mut)
singles_mRIP1$change_pos_c<-paste0(singles_mRIP1$WT_AA, singles_mRIP1$Pos_c, singles_mRIP1$Mut)
singles_mRIP1$change_aa<-paste0(singles_mRIP1$wt_mRIP1, singles_mRIP1$mut_mRIP1)
singles_mRIP1$change_type<-paste0(singles_mRIP1$Pos_c, singles_mRIP1$type_mut)
singles_mRIP1$change_pos_type<-paste0(singles_mRIP1$type_wt, singles_mRIP1$Pos_c, singles_mRIP1$type_mut)

# Correlation per position - mut
singles_ripk1m1<-inner_join(singles_ripk1, singles_mRIP1, by="change")
singles_ripk1m1<-rename(singles_ripk1m1, Pos = Pos_c.x, 
                        change_ripk1 = change_pos.x, change_mRIP1 = change_pos.y,
                        change_ripk1_c = change_pos_c.x, change_mRIP1_c = change_pos_c.y)
singles_ripk1m1<-singles_ripk1m1[c("change", "change_ripk1", "change_mRIP1", "change_ripk1_c", "change_mRIP1_c", 
                                   "Pos", "residue_ripk1", "residue_mRIP1", "nscore_c_ripk1", "nscore_c_mRIP1")]
singles_ripk1m1<-arrange(singles_ripk1m1, Pos)
singles_ripk1m1$RHIM<-FALSE
singles_ripk1m1[(singles_ripk1m1$Pos>8) & (singles_ripk1m1$Pos<26),]$RHIM<-TRUE
singles_ripk1m1$shared_pos<-FALSE
singles_ripk1m1[singles_ripk1m1$change_ripk1_c == singles_ripk1m1$change_mRIP1_c,]$shared_pos<-TRUE

correlation_all<-cor.test(singles_ripk1m1$nscore_c_ripk1, singles_ripk1m1$nscore_c_mRIP1, use="complete.obs")

correlation_rhim<-cor.test(singles_ripk1m1[singles_ripk1m1$RHIM == TRUE,]$nscore_c_ripk1, 
                           singles_ripk1m1[singles_ripk1m1$RHIM == TRUE,]$nscore_c_mRIP1, use="complete.obs")

correlation_out<-cor.test(singles_ripk1m1[singles_ripk1m1$RHIM == FALSE,]$nscore_c_ripk1, 
                          singles_ripk1m1[singles_ripk1m1$RHIM == FALSE,]$nscore_c_mRIP1, use="complete.obs")

p_corr_all<-ggplot(singles_ripk1m1, aes(x=nscore_c_ripk1, y=nscore_c_mRIP1))+
  geom_hline(yintercept = 0, color="black", linewidth=.2)+
  geom_vline(xintercept = 0, color="black", linewidth=.2)+
  geom_point(aes(color=RHIM), alpha=.6)+
  scale_color_manual(values=c("grey70", "coral"))+
  theme_classic()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  annotate("text", label=paste0("R=", round(correlation_all$estimate,2)), x=-Inf, y=Inf, hjust=-0.5, vjust=2, size=4, color="black")+
  annotate("text", label=paste0("p=", format(correlation_all$p.value, digits = 2, scientific = T)), x=-Inf, y=Inf, hjust=-0.2, vjust=3, size=4, color="black")+
  annotate("text", label=paste0("R=", round(correlation_rhim$estimate,2)), x=-Inf, y=Inf,hjust=-0.5, vjust=4.5, size=4, color="coral")+
  annotate("text", label=paste0("p=", format(correlation_rhim$p.value, digits = 2, scientific = T)), x=-Inf, y=Inf,hjust=-0.2, vjust=5.5, size=4, color="coral")+
  annotate("text", label=paste0("R=", round(correlation_out$estimate,2)), x=-Inf, y=Inf,hjust=-0.5, vjust=7, size=4, color="grey60")+
  annotate("text", label=paste0("p=", format(correlation_out$p.value, digits = 2, scientific = T)), x=-Inf, y=Inf,hjust=-0.2, vjust=8, size=4, color="grey60")+
  labs(x="Nucleation Score RIPK1", y="Nucleation Score mRIP1")

p_corr_all

ggsave(p_corr_all, file="Corr_singles_nscore_pos_RIPK1_mRIP1.jpg", width = 5, height = 4, path=path)
ggsave(p_corr_all, file="Corr_singles_nscore_pos_RIPK1_mRIP1.pdf", width = 5, height = 4, path=path)

# Correlation per position
corr_vector=c()
for(i in singles_ripk1m1$Pos){
  correlation<-cor.test(singles_ripk1m1[singles_ripk1m1$Pos==i,][["nscore_c_ripk1"]],
                        singles_ripk1m1[singles_ripk1m1$Pos==i,][["nscore_c_mRIP1"]], use="complete.obs")
  corr<-correlation$estimate
  p_value<-correlation$p.value
  corr_vector=c(corr_vector, i, corr, p_value)
}

corr_text <- data.frame(matrix(corr_vector, ncol=3, byrow = T))
colnames(corr_text)<-c("Pos", "corr","pvalue")
for(i in c(2:3)){corr_text[[i]]<-as.numeric(as.character(corr_text[[i]]))}
corr_text<-distinct(corr_text, Pos, .keep_all = TRUE)

#
p_corr_pos<-ggplot(singles_ripk1m1, aes(x=nscore_c_ripk1, y=nscore_c_mRIP1))+
  geom_rect(data = singles_ripk1m1, aes(fill = RHIM), xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.05)+
  scale_fill_manual(values=c("white", "mistyrose"))+
  geom_hline(yintercept = 0, linetype="dashed", color="darkgrey")+
  geom_vline(xintercept = 0, linetype="dashed", color="darkgrey")+
  geom_point()+
  theme_bw()+
  facet_wrap(~factor(Pos, label=unique(singles_ripk1m1$Pos+start_p-1)))+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  geom_text(data=corr_text, aes(label=paste0("R=",round(corr, 2)), x=-Inf, y=Inf), hjust=-0.05, vjust=1.5, size=4, colour="red")+
  geom_text(data=corr_text, aes(label=paste0("p=", format(pvalue, digits = 2, scientific = T)), x=-Inf, y=Inf), hjust=-0.05, vjust=3, size=4, colour="red")+
  geom_text(data=singles_ripk1m1, aes(label=residue_ripk1, x=-Inf, y=-Inf, color="orange"), hjust=-0.2, vjust=-0.5)+
  geom_text(data=singles_ripk1m1, aes(label=residue_mRIP1, x=Inf, y=-Inf, color="brown"), hjust=1.1, vjust=-0.5)+
  scale_color_manual(name = "Residue", labels = c("Mouse", "Human"), values=c("brown", "orange"))+
  labs(x="Nucleation Score RIPK1", y="Nucleation Score mRIP1")
p_corr_pos

ggsave(p_corr_pos, file="Corr_singles_nscore_poseach_RIPK1_mRIP1.jpg", width = 12, height = 10, path=path)
ggsave(p_corr_pos, file="Corr_singles_nscore_poseach_RIPK1_mRIP1.pdf", width = 12, height = 10, path=path)

# Bar plot of correlation scores
corr_text[is.na(corr_text)]<-1

corr_text$significance_pos<-""
corr_text[(corr_text$pvalue<0.05) & (corr_text$corr>0), "significance_pos"]<-"*"
corr_text[(corr_text$pvalue<0.01) & (corr_text$corr>0), "significance_pos"]<-"**"

corr_text$significance_neg<-""
corr_text[(corr_text$pvalue<0.05) & (corr_text$corr<0), "significance_neg"]<-"*"
corr_text[(corr_text$pvalue<0.01) & (corr_text$corr<0), "significance_neg"]<-"**"

peptide_seq<-'PVADDLIKYTIFNSSGIQIGNHNYMDVGLNSQPP'
peptide_seq<-c(strsplit(peptide_seq, '')[[1]])

peptide_seq_pos<-c()
for (n in seq_along(peptide_seq)){
  peptide_seq_pos<-c(peptide_seq_pos, paste0(peptide_seq[n], '\n', n+start_p_m-1))
}

p_bars<-ggplot(corr_text, aes(x=factor(Pos, labels=peptide_seq_pos), y=corr))+
  geom_hline(yintercept = 0, color="black", linewidth=.5)+
  geom_bar(stat="identity", color="black", fill="grey")+
  geom_text(data=corr_text, aes(label=significance_pos), vjust=0.4, size=8, colour="black")+
  geom_text(data=corr_text, aes(label=significance_neg), vjust=1, size=8, colour="black")+
  theme_void()+
  labs(x="RIPK1", y="Pearson Correlation", title="NS vs mRIP1")+
  ylim(-1,1)+
  theme(legend.title = element_text(size=12, face = "bold"),
        legend.text = element_text(size=12), 
        axis.title.y = element_text(size = 14, angle = 90),
        axis.text.x = element_text(color=color_axis_x, size=12),
        axis.text.y = element_text(color="black", size=12),
        axis.line.y = element_line(color="black"))

p_bars

ggsave(p_bars, file="Corr_nscore_mRIP1_bars.jpg", width = 15, height = 3, path=path)
ggsave(p_bars, file="Corr_nscore_mRIP1_bars.pdf", width = 15, height = 3, path=path)

################################################################################
### Correlation per aa change

singles_ripk1_grouped<-as.data.frame(singles_ripk1 %>% group_by(change_aa) %>% dplyr::summarise(mean_ns_ripk1=mean(nscore_c_ripk1)))
singles_mRIP1_grouped<-as.data.frame(singles_mRIP1 %>% group_by(change_aa) %>% dplyr::summarise(mean_ns_mRIP1=mean(nscore_c_mRIP1)))

singles_ripk1m1<-inner_join(singles_ripk1_grouped, singles_mRIP1_grouped, by="change_aa")

correlation_all<-cor.test(singles_ripk1m1$mean_ns_ripk1, singles_ripk1m1$mean_ns_mRIP1, use="complete.obs")


p_corr_all<-ggplot(singles_ripk1m1, aes(x=mean_ns_ripk1, y=mean_ns_mRIP1))+
  geom_hline(yintercept = 0, color="black", linewidth=.2)+
  geom_vline(xintercept = 0, color="black", linewidth=.2)+
  geom_point(color="grey60", alpha=.6)+
  theme_classic()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  annotate("text", label=paste0("R=", round(correlation_all$estimate,2)), x=-Inf, y=Inf, hjust=-0.5, vjust=2, size=5, color="black")+
  annotate("text", label=paste0("p=", format(correlation_all$p.value, digits = 2, scientific = T)), x=-Inf, y=Inf, hjust=-0.2, vjust=3, size=5, color="black")+
  labs(x="Nucleation Score RIPK1", y="Nucleation Score mRIP1")
p_corr_all

ggsave(p_corr_all, file="Corr_singles_nscore_aachange_RIPK1_mRIP1.jpg", width = 5, height = 4, path=path)
ggsave(p_corr_all, file="Corr_singles_nscore_aachange_RIPK1_mRIP1.pdf", width = 5, height = 4, path=path)
