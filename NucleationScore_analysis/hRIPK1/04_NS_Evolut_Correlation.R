library(tidyverse)
require(reshape2)
library(ggpubr)
library(viridis)

# your folder path
setwd("")

# Create a folder for the figures
dir.create("04_NS_evol_corrpearson")
path="04_NS_evol_corrpearson"

# open processed data
load("nscore_df_hRIPK1.RData")

# start and end position in the library
start_p<-523
end_p<-556

# amino acid sequence mutagenised
peptide_seq<-'PPTDESIKYTIYNSTGIQIGAYNYMEIGGTSSSL'
peptide_seq<-c(strsplit(peptide_seq, '')[[1]])

peptide_seq_pos<-c()
for (n in seq_along(peptide_seq)){
  peptide_seq_pos<-c(peptide_seq_pos, paste0(peptide_seq[n], '\n', n+start_p-1))
}

# all the amino acid
all_aa<-c("G","A","V","L","M","I","F","Y","W","K","R","D","E","S","T","C","N","Q","H","P")

AA_type<-data.frame("wt_aa"= all_aa,
                    "mut_aa" = all_aa,
                    "name_AA"=c("Glycine", "Alanine","Valine","Leucine","Methionine","Isoleucine","Phenylalanine",
                                "Tyrosine","Tryptophan","Lysine","Arginine","Aspartic acid","Glutamic acid","Serine","Threonine",
                                "Cysteine","Asparagine","Glutamine","Histidine","Proline"),
                    "type"=c("glycine",rep("aliphatic",5),rep("aromatic",3),rep("positive",2),rep("negative",2),rep("polar",6),"proline"))

color_axis_x=c("N523"="black", "P443"="black", "V444"="black", "T445"="black", "G446"="black", "R447"="black",
               "P448"="black", "L449"="black", "V531"="red", "N451"="red", "I452"="red", "Y453"="red",
               "N454"="red", "C455"="red", "S456"="red", "G457"="red", "V458"="red", "Q459"="red",
               "V460"="red", "G461"="red", "D462"="red", "N463"="red", "N464"="red", "Y547"= "red",
               "L466"="red", "T467"="black", "M468"="black", "Q469"="black","Q470"="black", "T471"="black",
               "T472"="black", "A473"="black", "L468"="black", "P469"="black")

# Load GEMME data 
gemme_score<-read_tsv("RIPK1_GEMME.tsv")
gemme_score<-gemme_score %>% drop_na(GEMME_score)

# Generate new columns based on aa type
gemme_score<-left_join(gemme_score, AA_type[,c("wt_aa", "type")], by=c("WT_AA"="wt_aa"))
gemme_score<-rename(gemme_score, type_wt = type)
gemme_score<-left_join(gemme_score, AA_type[,c("wt_aa", "type")], by=c("variable"="wt_aa"))
gemme_score<-rename(gemme_score, type_mut = type)

gemme_score$change_type<-paste0(gemme_score$type_wt, gemme_score$Pos, gemme_score$type_mut)


# Generate new columns based on aa type
singles_stops<-arrange(singles_stops, Pos)

singles_stops$change<-paste0(singles_stops$WT_AA, singles_stops$Pos, singles_stops$Mut)
singles_stops<-left_join(singles_stops, AA_type[,c("wt_aa", "type")], by=c("WT_AA"="wt_aa"))
singles_stops<-rename(singles_stops, type_wt = type)
singles_stops<-left_join(singles_stops, AA_type[,c("wt_aa", "type")], by=c("Mut"="wt_aa"))
singles_stops<-rename(singles_stops, type_mut = type)

singles_stops$Residue<-paste0(singles_stops$WT_AA, singles_stops$Pos)
singles_stops$change_type<-paste0(singles_stops$type_wt, singles_stops$Pos, singles_stops$type_mut)

# join both datasets
RIPK1_gemme_nscore<-inner_join(singles_stops, gemme_score, by="change")
RIPK1_gemme_nscore<-rename(RIPK1_gemme_nscore, "WT_AA" = "WT_AA.x", "Pos" = "Pos.x")
RIPK1_gemme_nscore<-RIPK1_gemme_nscore[,c("change", "WT_AA", "Pos", "Mut", "nscore_c", "GEMME_score")]
RIPK1_gemme_nscore$Residue<-paste0(RIPK1_gemme_nscore$WT_AA, RIPK1_gemme_nscore$Pos)
RIPK1_gemme_nscore$RHIM<-FALSE
RIPK1_gemme_nscore[(RIPK1_gemme_nscore$Pos>=531) & (RIPK1_gemme_nscore$Pos<=547),]$RHIM<-TRUE

RIPK1_gemme_nscore$Pos<-as.integer(RIPK1_gemme_nscore$Pos)
RIPK1_gemme_nscore<-arrange(RIPK1_gemme_nscore, Pos)

### Correlation GEMME score vs NScore
correlation_all<-cor.test(RIPK1_gemme_nscore$nscore_c, RIPK1_gemme_nscore$GEMME_score, use="complete.obs")

correlation_rhim<-cor.test(RIPK1_gemme_nscore[RIPK1_gemme_nscore$RHIM == TRUE,]$nscore_c, 
                           RIPK1_gemme_nscore[RIPK1_gemme_nscore$RHIM == TRUE,]$GEMME_score, use="complete.obs")

correlation_out<-cor.test(RIPK1_gemme_nscore[RIPK1_gemme_nscore$RHIM == FALSE,]$nscore_c, 
                          RIPK1_gemme_nscore[RIPK1_gemme_nscore$RHIM == FALSE,]$GEMME_score, use="complete.obs")


p_corr_all<-ggplot(RIPK1_gemme_nscore, aes(x=nscore_c, y=GEMME_score))+
  geom_hline(yintercept = 0, color="black", linewidth=.2)+
  geom_vline(xintercept = 0, color="black", linewidth=.2)+
  geom_point(aes(color=RHIM), alpha=.6)+
  scale_color_manual(values=c("grey70", "coral"))+
  theme_classic()+
  theme(panel.grid = element_blank())+
  annotate("text", label=paste0("R=", round(correlation_all$estimate,2)), x=-Inf, y=Inf, hjust=-0.5, vjust=1, size=4, color="black")+
  annotate("text", label=paste0("p=", format(correlation_all$p.value, digits = 2, scientific = T)), x=-Inf, y=Inf, hjust=-0.2, vjust=2, size=4, color="black")+
  annotate("text", label=paste0("R=", round(correlation_rhim$estimate,2)), x=-Inf, y=Inf,hjust=-0.5, vjust=3.5, size=4, color="coral")+
  annotate("text", label=paste0("p=", format(correlation_rhim$p.value, digits = 2, scientific = T)), x=-Inf, y=Inf,hjust=-0.2, vjust=4.5, size=4, color="coral")+
  annotate("text", label=paste0("R=", round(correlation_out$estimate,2)), x=-Inf, y=Inf,hjust=-0.5, vjust=6, size=4, color="grey60")+
  annotate("text", label=paste0("p=", format(correlation_out$p.value, digits = 2, scientific = T)), x=-Inf, y=Inf,hjust=-0.2, vjust=7, size=4, color="grey60")+
  labs(x="Nucleation Score", y="GEMME Score")
p_corr_all

ggsave(p_corr_all, file="Corr_nscore_gemme_.jpg", width = 5, height = 4, path=path)
ggsave(p_corr_all, file="Corr_nscore_gemme_.pdf", width = 5, height = 4, path=path)

### Correlation GEMME score vs NScore per position

corr_vector=c()
for(i in RIPK1_gemme_nscore$Pos){
  correlation<-cor.test(RIPK1_gemme_nscore[RIPK1_gemme_nscore$Pos==i,]$nscore_c,
                        RIPK1_gemme_nscore[RIPK1_gemme_nscore$Pos==i,]$GEMME_score, use="complete.obs")
  corr<-correlation$estimate
  p_value<-correlation$p.value
  corr_vector=c(corr_vector, i, corr, p_value)
}

corr_text <- data.frame(matrix(corr_vector, ncol=3, byrow = T))
colnames(corr_text)<-c("Pos", "corr","pvalue")
for(i in c(2:3)){corr_text[[i]]<-as.numeric(as.character(corr_text[[i]]))}
corr_text<-distinct(corr_text, Pos, .keep_all = TRUE)

#
p_corr_pos<-ggplot(RIPK1_gemme_nscore, aes(x=nscore_c, y=GEMME_score))+
  geom_rect(data = RIPK1_gemme_nscore, aes(fill = RHIM), xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.05)+
  scale_fill_manual(values=c("white", "mistyrose"))+
  geom_hline(yintercept = 0, linetype="dashed", color="grey")+
  geom_vline(xintercept = 0, linetype="dashed", color="grey")+
  geom_point(color="grey20")+
  theme_bw()+
  facet_wrap(~factor(Pos, labels=distinct(RIPK1_gemme_nscore, Residue)[[1]]))+
  theme(panel.grid = element_blank())+
  geom_text(data=corr_text, aes(label=paste0("R=",round(corr, 2)), x=-Inf, y=Inf), hjust=-0.1, vjust=1.5, size=4, colour="red")+
  geom_text(data=corr_text, aes(label=paste0("p=", format(pvalue, digits = 2, scientific = T)), x=-Inf, y=Inf), hjust=-0.1, vjust=3, size=4, colour="red")+
  scale_color_manual(name = "Residue", labels = c("Mouse", "Human"), values=c("brown", "orange"))+
  labs(x="Nucleation Score", y="GEMME Score")
p_corr_pos

ggsave(p_corr_pos, file="Corr_nscore_gemme_pos.jpg", width = 12, height = 10, path=path)
ggsave(p_corr_pos, file="Corr_nscore_gemme_pos.pdf", width = 12, height = 10, path=path)

# make a bar plot of the correlations
corr_text[is.na(corr_text)]<-1

corr_text$significance_pos<-""
corr_text[(corr_text$pvalue<0.05) & (corr_text$corr>0), "significance_pos"]<-"*"
corr_text[(corr_text$pvalue<0.01) & (corr_text$corr>0), "significance_pos"]<-"**"

corr_text$significance_neg<-""
corr_text[(corr_text$pvalue<0.05) & (corr_text$corr<0), "significance_neg"]<-"*"
corr_text[(corr_text$pvalue<0.01) & (corr_text$corr<0), "significance_neg"]<-"**"

p_bars<-ggplot(corr_text, aes(x=factor(Pos, labels=peptide_seq_pos), y=corr))+
  geom_hline(yintercept = 0, color="black", linewidth=.5)+
  geom_bar(stat="identity", color="black", fill="grey")+
  geom_text(data=corr_text, aes(label=significance_pos), vjust=0.4, size=8, colour="black")+
  geom_text(data=corr_text, aes(label=significance_neg), vjust=1, size=8, colour="black")+
  theme_void()+
  labs(x="RIPK1", y="Pearson Correlation", title="NS vs GEMME")+
  ylim(-1,1)+
  theme(legend.title = element_text(size=12, face = "bold"),
        legend.text = element_text(size=12), 
        axis.title.y = element_text(size = 14, angle = 90),
        axis.text.x = element_text(color=color_axis_x, size=12),
        axis.text.y = element_text(color="black", size=12),
        axis.line.y = element_line(color="black"))

p_bars

ggsave(p_bars, file="Corr_nscore_gemme_bars.jpg", width = 15, height = 3, path=path)
ggsave(p_bars, file="Corr_nscore_gemme_bars.pdf", width = 15, height = 3, path=path)

# GEMME heatmap
vectorAA <- c("G","A","V","L","M","I","F","Y","W","K","R","D","E","S","T","C","N","Q","H", "P")

gemme_score$Pos<-substr(gemme_score$residue, 2, 4)
gemme_score$Pos<-as.integer(gemme_score$Pos)
gemme_score<-rename(gemme_score, Mut = variable, Residue = residue)

p_gemme_heatmap<-ggplot(gemme_score)+
  geom_tile(aes(Pos, factor(Mut, levels=rev(vectorAA)), fill=GEMME_score), color='white', size=1)+
  scale_x_continuous(breaks=seq(start_p, end_p), labels = peptide_seq_pos, expand=c(0,0))+
  scale_fill_viridis(option = "plasma", direction = -1)+
  theme_minimal()+
  labs(x="RIPK1", y="Mutant amino acid", fill="GEMME score")+
  theme(legend.title = element_text(size=12),
        legend.text = element_text(size=12), 
        axis.title = element_text(size = 14),
        axis.text.x = element_text(color=c(rep("black", 8), rep("red", 17), rep("black", 9)), size=12),
        axis.text.y = element_text(color="black", size=12))
p_gemme_heatmap

ggsave(p_gemme_heatmap, file="gemme_heatmap.jpg", width=15, height=8, path=path)

#################################################################################
### Correlation per position considering amino acid type

singles_stops_group<-as.data.frame(singles_stops %>% group_by(Pos, Residue, change_type) %>% dplyr::summarise(mean_ns=mean(nscore_c)))

singles_gemme_group<-as.data.frame(gemme_score %>% group_by(Pos, Residue, change_type) %>% dplyr::summarise(mean_gemme=mean(GEMME_score)))

singles_nscore_gemme<-inner_join(singles_stops_group, singles_gemme_group, by="change_type")
singles_nscore_gemme<-rename(singles_nscore_gemme, "Pos" = "Pos.x")
singles_nscore_gemme<-singles_nscore_gemme[c("change_type", "Pos", "mean_ns", "mean_gemme")]
singles_nscore_gemme<-arrange(singles_nscore_gemme, Pos)
singles_nscore_gemme$RHIM<-FALSE
singles_nscore_gemme[(singles_nscore_gemme$Pos>=531) & (singles_nscore_gemme$Pos<=547),]$RHIM<-TRUE


correlation_all<-cor.test(singles_nscore_gemme$mean_ns, singles_nscore_gemme$mean_gemme, use="complete.obs")

correlation_rhim<-cor.test(singles_nscore_gemme[singles_nscore_gemme$RHIM == TRUE,]$mean_ns, 
                           singles_nscore_gemme[singles_nscore_gemme$RHIM == TRUE,]$mean_gemme, use="complete.obs")

correlation_out<-cor.test(singles_nscore_gemme[singles_nscore_gemme$RHIM == FALSE,]$mean_ns, 
                          singles_nscore_gemme[singles_nscore_gemme$RHIM == FALSE,]$mean_gemme, use="complete.obs")

p_corr_all<-ggplot(singles_nscore_gemme, aes(x=mean_ns, y=mean_gemme))+
  geom_hline(yintercept = 0, color="black", linewidth=.2)+
  geom_vline(xintercept = 0, color="black", linewidth=.2)+
  geom_point(aes(color=RHIM), alpha=.6)+
  scale_shape_manual(values=c(19, 6))+
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
  labs(x="Nucleation Score", y="GEMME Score")
p_corr_all

ggsave(p_corr_all, file="Corr_singles_nscore_postype_RIPK1_gemme.jpg", width = 5, height = 4, path=path)
ggsave(p_corr_all, file="Corr_singles_nscore_postype_RIPK1_gemme.pdf", width = 5, height = 4, path=path)

# Correlation per position based on amino acid type
corr_vector=c()
for(i in singles_nscore_gemme$Pos){
  correlation<-cor.test(singles_nscore_gemme[singles_nscore_gemme$Pos==i,][["mean_ns"]],
                        singles_nscore_gemme[singles_nscore_gemme$Pos==i,][["mean_gemme"]], use="complete.obs")
  corr<-correlation$estimate
  p_value<-correlation$p.value
  corr_vector=c(corr_vector, i, corr, p_value)
}

corr_text <- data.frame(matrix(corr_vector, ncol=3, byrow = T))
colnames(corr_text)<-c("Pos", "corr","pvalue")
for(i in c(2:3)){corr_text[[i]]<-as.numeric(as.character(corr_text[[i]]))}
corr_text<-distinct(corr_text, Pos, .keep_all = TRUE)

#
p_corr_pos<-ggplot(singles_nscore_gemme, aes(x=mean_ns, y=mean_gemme))+
  geom_rect(data = singles_nscore_gemme, aes(fill = RHIM), xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.05)+
  scale_fill_manual(values=c("white", "mistyrose"))+
  geom_hline(yintercept = 0, linetype="dashed", color="darkgrey")+
  geom_vline(xintercept = 0, linetype="dashed", color="darkgrey")+
  geom_point()+
  theme_bw()+
  facet_wrap(~factor(Pos, label = unique(singles_stops$Residue)))+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  geom_text(data=corr_text, aes(label=paste0("R=",round(corr, 2)), x=-Inf, y=Inf), hjust=-0.05, vjust=1.5, size=4, colour="red")+
  geom_text(data=corr_text, aes(label=paste0("p=", format(pvalue, digits = 2, scientific = T)), x=-Inf, y=Inf), hjust=-0.05, vjust=3, size=4, colour="red")+
  scale_color_manual(name = "Residue", labels = c("Mouse", "Human"), values=c("brown", "orange"))+
  labs(x="Nucleation Score", y="GEMME Score")
p_corr_pos

ggsave(p_corr_pos, file="Corr_singles_nscore_postypeeach_RIPK1_gemme.jpg", width = 12, height = 10, path=path)
ggsave(p_corr_pos, file="Corr_singles_nscore_postypeeach_RIPK1_gemme.pdf", width = 12, height = 10, path=path)

# bar plot based on correlation scores
corr_text[is.na(corr_text)]<-1

corr_text$significance_pos<-""
corr_text[(corr_text$pvalue<0.05) & (corr_text$corr>0), "significance_pos"]<-"*"
corr_text[(corr_text$pvalue<0.01) & (corr_text$corr>0), "significance_pos"]<-"**"

corr_text$significance_neg<-""
corr_text[(corr_text$pvalue<0.05) & (corr_text$corr<0), "significance_neg"]<-"*"
corr_text[(corr_text$pvalue<0.01) & (corr_text$corr<0), "significance_neg"]<-"**"


p_bars<-ggplot(corr_text, aes(x=factor(Pos, labels=peptide_seq_pos), y=corr))+
  geom_hline(yintercept = 0, color="black", linewidth=.5)+
  geom_bar(stat="identity", color="black", fill="grey")+
  geom_text(data=corr_text, aes(label=significance_pos), vjust=0.4, size=8, colour="black")+
  geom_text(data=corr_text, aes(label=significance_neg), vjust=1, size=8, colour="black")+
  theme_void()+
  labs(x="RIPK1", y="Pearson Correlation", title="NS vs gemme - per aa type")+
  ylim(-1,1)+
  theme(legend.title = element_text(size=12, face = "bold"),
        legend.text = element_text(size=12), 
        axis.title.y = element_text(size = 14, angle = 90),
        axis.text.x = element_text(color=color_axis_x, size=12),
        axis.text.y = element_text(color="black", size=12),
        axis.line.y = element_line(color="black"))

p_bars

ggsave(p_bars, file="Corr_nscore_gemme_postype_bars.jpg", width = 15, height = 3, path=path)
ggsave(p_bars, file="Corr_nscore_gemme_postype_bars.pdf", width = 15, height = 3, path=path)


##################################################################################
# Same approach as for GEMME but for ESM1b
# Load ESM1b data

ESM1b_score<-read_csv("RIPK1_ESM1b.csv")
ESM1b_score<-rename(ESM1b_score, change = variant, ESM1b_score = score)
ESM1b_score<-left_join(ESM1b_score, AA_type[,c("wt_aa", "type")], by=c("WT_AA"="wt_aa"))
ESM1b_score<-rename(ESM1b_score, type_wt = type)
ESM1b_score<-left_join(ESM1b_score, AA_type[,c("wt_aa", "type")], by=c("Mut"="wt_aa"))
ESM1b_score<-rename(ESM1b_score, type_mut = type)

ESM1b_score$change_type<-paste0(ESM1b_score$type_wt, ESM1b_score$pos, ESM1b_score$type_mut)


RIPK1_ESM1b_nscore<-inner_join(singles_stops, ESM1b_score, by="change")
RIPK1_ESM1b_nscore<-rename(RIPK1_ESM1b_nscore, "WT_AA" = "WT_AA.x", "Mut" = "Mut.x")
RIPK1_ESM1b_nscore<-RIPK1_ESM1b_nscore[,c("change", "WT_AA", "Pos", "Mut", "nscore_c", "ESM1b_score")]
RIPK1_ESM1b_nscore$Residue<-paste0(RIPK1_ESM1b_nscore$WT_AA, RIPK1_ESM1b_nscore$Pos)
RIPK1_ESM1b_nscore$RHIM<-FALSE
RIPK1_ESM1b_nscore[(RIPK1_ESM1b_nscore$Pos>=531) & (RIPK1_ESM1b_nscore$Pos<=547),]$RHIM<-TRUE

RIPK1_ESM1b_nscore$Pos<-as.integer(RIPK1_ESM1b_nscore$Pos)
RIPK1_ESM1b_nscore<-arrange(RIPK1_ESM1b_nscore, Pos)

### Correlation ESM1b score vs NScore
correlation_all<-cor.test(RIPK1_ESM1b_nscore$nscore_c, RIPK1_ESM1b_nscore$ESM1b_score, use="complete.obs")

correlation_rhim<-cor.test(RIPK1_ESM1b_nscore[RIPK1_ESM1b_nscore$RHIM == TRUE,]$nscore_c, 
                           RIPK1_ESM1b_nscore[RIPK1_ESM1b_nscore$RHIM == TRUE,]$ESM1b_score, use="complete.obs")

correlation_out<-cor.test(RIPK1_ESM1b_nscore[RIPK1_ESM1b_nscore$RHIM == FALSE,]$nscore_c, 
                          RIPK1_ESM1b_nscore[RIPK1_ESM1b_nscore$RHIM == FALSE,]$ESM1b_score, use="complete.obs")

p_corr_all<-ggplot(RIPK1_ESM1b_nscore, aes(x=nscore_c, y=ESM1b_score))+
  geom_hline(yintercept = 0, color="black", linewidth=.2)+
  geom_vline(xintercept = 0, color="black", linewidth=.2)+
  geom_point(aes(color=RHIM), alpha=.6)+
  scale_color_manual(values=c("grey70", "coral"))+
  theme_classic()+
  theme(panel.grid = element_blank())+
  annotate("text", label=paste0("R=", round(correlation_all$estimate,2)), x=-Inf, y=Inf, hjust=-0.5, vjust=1, size=4, color="black")+
  annotate("text", label=paste0("p=", format(correlation_all$p.value, digits = 2, scientific = T)), x=-Inf, y=Inf, hjust=-0.2, vjust=2, size=4, color="black")+
  annotate("text", label=paste0("R=", round(correlation_rhim$estimate,2)), x=-Inf, y=Inf,hjust=-0.5, vjust=3.5, size=4, color="coral")+
  annotate("text", label=paste0("p=", format(correlation_rhim$p.value, digits = 2, scientific = T)), x=-Inf, y=Inf,hjust=-0.2, vjust=4.5, size=4, color="coral")+
  annotate("text", label=paste0("R=", round(correlation_out$estimate,2)), x=-Inf, y=Inf,hjust=-0.5, vjust=6, size=4, color="grey60")+
  annotate("text", label=paste0("p=", format(correlation_out$p.value, digits = 2, scientific = T)), x=-Inf, y=Inf,hjust=-0.2, vjust=7, size=4, color="grey60")+
  labs(x="Nucleation Score", y="ESM1b Score")
p_corr_all

ggsave(p_corr_all, file="Corr_nscore_ESM1b.jpg", width = 5, height = 4, path=path)
ggsave(p_corr_all, file="Corr_nscore_ESM1b.pdf", width = 5, height = 4, path=path)

### Correlation ESM1b score vs NScore per position

corr_vector=c()
for(i in RIPK1_ESM1b_nscore$Pos){
  correlation<-cor.test(RIPK1_ESM1b_nscore[RIPK1_ESM1b_nscore$Pos==i,]$nscore_c,
                        RIPK1_ESM1b_nscore[RIPK1_ESM1b_nscore$Pos==i,]$ESM1b_score, use="complete.obs")
  corr<-correlation$estimate
  p_value<-correlation$p.value
  corr_vector=c(corr_vector, i, corr, p_value)
}

corr_text <- data.frame(matrix(corr_vector, ncol=3, byrow = T))
colnames(corr_text)<-c("Pos", "corr","pvalue")
for(i in c(2:3)){corr_text[[i]]<-as.numeric(as.character(corr_text[[i]]))}
corr_text<-distinct(corr_text, Pos, .keep_all = TRUE)

#
p_corr_pos<-ggplot(RIPK1_ESM1b_nscore, aes(x=nscore_c, y=ESM1b_score))+
  geom_rect(data = RIPK1_ESM1b_nscore, aes(fill = RHIM), xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.05)+
  scale_fill_manual(values=c("white", "mistyrose"))+
  geom_hline(yintercept = 0, linetype="dashed", color="grey")+
  geom_vline(xintercept = 0, linetype="dashed", color="grey")+
  geom_point(color="grey20")+
  theme_bw()+
  facet_wrap(~factor(Pos, labels=distinct(RIPK1_ESM1b_nscore, Residue)[[1]]))+
  theme(panel.grid = element_blank())+
  geom_text(data=corr_text, aes(label=paste0("R=",round(corr, 2)), x=-Inf, y=Inf), hjust=-0.1, vjust=1.5, size=4, colour="red")+
  geom_text(data=corr_text, aes(label=paste0("p=", format(pvalue, digits = 2, scientific = T)), x=-Inf, y=Inf), hjust=-0.1, vjust=3, size=4, colour="red")+
  labs(x="Nucleation Score", y="ESM1b Score")
p_corr_pos

ggsave(p_corr_pos, file="Corr_nscore_ESM1b_pos.jpg", width = 12, height = 10, path=path)
ggsave(p_corr_pos, file="Corr_nscore_ESM1b_pos.pdf", width = 12, height = 10, path=path)

# bar plot
corr_text[is.na(corr_text)]<-1

corr_text$significance_pos<-""
corr_text[(corr_text$pvalue<0.05) & (corr_text$corr>0), "significance_pos"]<-"*"
corr_text[(corr_text$pvalue<0.01) & (corr_text$corr>0), "significance_pos"]<-"**"

corr_text$significance_neg<-""
corr_text[(corr_text$pvalue<0.05) & (corr_text$corr<0), "significance_neg"]<-"*"
corr_text[(corr_text$pvalue<0.01) & (corr_text$corr<0), "significance_neg"]<-"**"

p_bars<-ggplot(corr_text, aes(x=factor(Pos, labels=peptide_seq_pos), y=corr))+
  geom_hline(yintercept = 0, color="black", linewidth=.5)+
  geom_bar(stat="identity", color="black", fill="grey")+
  geom_text(data=corr_text, aes(label=significance_pos), vjust=0.4, size=8, colour="black")+
  geom_text(data=corr_text, aes(label=significance_neg), vjust=1, size=8, colour="black")+
  theme_void()+
  labs(x="RIPK1", y="Pearson Correlation", title= "NS vs ESM1b")+
  ylim(-1,1)+
  theme(legend.title = element_text(size=12, face = "bold"),
        legend.text = element_text(size=12), 
        axis.title.y = element_text(size = 14, angle = 90),
        axis.text.x = element_text(color=color_axis_x, size=12),
        axis.text.y = element_text(color="black", size=12),
        axis.line.y = element_line(color="black"))

p_bars

ggsave(p_bars, file="Corr_nscore_ESM1b_bars.jpg", width = 15, height = 3, path=path)
ggsave(p_bars, file="Corr_nscore_ESM1b_bars.pdf", width = 15, height = 3, path=path)

# ESM1b heatmap

ESM1b_score<-ESM1b_score[ESM1b_score$pos <= end_p & ESM1b_score$pos >= start_p,]
ESM1b_score$Residue<-substr(ESM1b_score$change, 1, 4)
ESM1b_score$Mut<-substr(ESM1b_score$change, 5, 5)
ESM1b_score$Pos<-as.integer(ESM1b_score$pos)


p_ESM1b_heatmap<-ggplot(ESM1b_score)+
  geom_tile(aes(Pos, factor(Mut, levels=rev(vectorAA)), fill=ESM1b_score), color='white', size=1)+
  scale_x_continuous(breaks=seq(start_p, end_p), labels = peptide_seq_pos, expand=c(0,0))+
  #scale_fill_gradient2(low="gold", high="purple", mid="white")+
  scale_fill_viridis(direction=-1)+
  theme_minimal()+
  labs(x="RIPK1", y="Mutant amino acid", fill="ESM1b score")+
  theme(legend.title = element_text(size=12),
        legend.text = element_text(size=12), 
        axis.title = element_text(size = 14),
        axis.text.x = element_text(color=c(rep("black", 8), rep("red", 17), rep("black", 9)), size=12),
        axis.text.y = element_text(color="black", size=12))
p_ESM1b_heatmap

ggsave(p_ESM1b_heatmap, file="ESM1b_heatmap.jpg", width=15, height=8, path=path)

#################################################################################
### Correlation per position considering amino acid type

singles_stops_group<-as.data.frame(singles_stops %>% group_by(Pos, Residue, change_type) %>% dplyr::summarise(mean_ns=mean(nscore_c)))

singles_ESM1b_group<-as.data.frame(ESM1b_score %>% group_by(pos, residue, change_type) %>% dplyr::summarise(mean_ESM1b=mean(ESM1b_score)))

singles_nscore_ESM1b<-inner_join(singles_stops_group, singles_ESM1b_group, by="change_type")
singles_nscore_ESM1b<-singles_nscore_ESM1b[c("change_type", "Pos", "mean_ns", "mean_ESM1b")]
singles_nscore_ESM1b<-arrange(singles_nscore_ESM1b, Pos)
singles_nscore_ESM1b$RHIM<-FALSE
singles_nscore_ESM1b[(singles_nscore_ESM1b$Pos>=531) & (singles_nscore_ESM1b$Pos<=547),]$RHIM<-TRUE

correlation_all<-cor.test(singles_nscore_ESM1b$mean_ns, singles_nscore_ESM1b$mean_ESM1b, use="complete.obs")

correlation_rhim<-cor.test(singles_nscore_ESM1b[singles_nscore_ESM1b$RHIM == TRUE,]$mean_ns, 
                           singles_nscore_ESM1b[singles_nscore_ESM1b$RHIM == TRUE,]$mean_ESM1b, use="complete.obs")

correlation_out<-cor.test(singles_nscore_ESM1b[singles_nscore_ESM1b$RHIM == FALSE,]$mean_ns, 
                          singles_nscore_ESM1b[singles_nscore_ESM1b$RHIM == FALSE,]$mean_ESM1b, use="complete.obs")

p_corr_all<-ggplot(singles_nscore_ESM1b, aes(x=mean_ns, y=mean_ESM1b))+
  geom_hline(yintercept = 0, color="black", linewidth=.2)+
  geom_vline(xintercept = 0, color="black", linewidth=.2)+
  geom_point(aes(color=RHIM), alpha=.6)+
  scale_shape_manual(values=c(19, 6))+
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
  labs(x="Nucleation Score", y="ESM1b Score")
p_corr_all

ggsave(p_corr_all, file="Corr_singles_nscore_postype_RIPK1_ESM1b.jpg", width = 5, height = 4, path=path)
ggsave(p_corr_all, file="Corr_singles_nscore_postype_RIPK1_ESM1b.pdf", width = 5, height = 4, path=path)

# Correlation per position considering amino acid type
corr_vector=c()
for(i in singles_nscore_ESM1b$Pos){
  correlation<-cor.test(singles_nscore_ESM1b[singles_nscore_ESM1b$Pos==i,][["mean_ns"]],
                        singles_nscore_ESM1b[singles_nscore_ESM1b$Pos==i,][["mean_ESM1b"]], use="complete.obs")
  corr<-correlation$estimate
  p_value<-correlation$p.value
  corr_vector=c(corr_vector, i, corr, p_value)
}

corr_text <- data.frame(matrix(corr_vector, ncol=3, byrow = T))
colnames(corr_text)<-c("Pos", "corr","pvalue")
for(i in c(2:3)){corr_text[[i]]<-as.numeric(as.character(corr_text[[i]]))}
corr_text<-distinct(corr_text, Pos, .keep_all = TRUE)

#
p_corr_pos<-ggplot(singles_nscore_ESM1b, aes(x=mean_ns, y=mean_ESM1b))+
  geom_rect(data = singles_nscore_ESM1b, aes(fill = RHIM), xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.05)+
  scale_fill_manual(values=c("white", "mistyrose"))+
  geom_hline(yintercept = 0, linetype="dashed", color="darkgrey")+
  geom_vline(xintercept = 0, linetype="dashed", color="darkgrey")+
  geom_point()+
  theme_bw()+
  facet_wrap(~factor(Pos, label = unique(singles_stops$Residue)))+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  geom_text(data=corr_text, aes(label=paste0("R=",round(corr, 2)), x=-Inf, y=Inf), hjust=-0.05, vjust=1.5, size=4, colour="red")+
  geom_text(data=corr_text, aes(label=paste0("p=", format(pvalue, digits = 2, scientific = T)), x=-Inf, y=Inf), hjust=-0.05, vjust=3, size=4, colour="red")+
  scale_color_manual(name = "Residue", labels = c("Mouse", "Human"), values=c("brown", "orange"))+
  labs(x="Nucleation Score", y="ESM1b Score")
p_corr_pos

ggsave(p_corr_pos, file="Corr_singles_nscore_postypeeach_RIPK1_ESM1b.jpg", width = 12, height = 10, path=path)
ggsave(p_corr_pos, file="Corr_singles_nscore_postypeeach_RIPK1_ESM1b.pdf", width = 12, height = 10, path=path)

# Bar plot
corr_text[is.na(corr_text)]<-1

corr_text$significance_pos<-""
corr_text[(corr_text$pvalue<0.05) & (corr_text$corr>0), "significance_pos"]<-"*"
corr_text[(corr_text$pvalue<0.01) & (corr_text$corr>0), "significance_pos"]<-"**"

corr_text$significance_neg<-""
corr_text[(corr_text$pvalue<0.05) & (corr_text$corr<0), "significance_neg"]<-"*"
corr_text[(corr_text$pvalue<0.01) & (corr_text$corr<0), "significance_neg"]<-"**"

p_bars<-ggplot(corr_text, aes(x=factor(Pos, labels=peptide_seq_pos), y=corr))+
  geom_hline(yintercept = 0, color="black", linewidth=.5)+
  geom_bar(stat="identity", color="black", fill="grey")+
  geom_text(data=corr_text, aes(label=significance_pos), vjust=0.4, size=8, colour="black")+
  geom_text(data=corr_text, aes(label=significance_neg), vjust=1, size=8, colour="black")+
  theme_void()+
  labs(x="RIPK1", y="Pearson Correlation", title="NS vs ESM1b - per aa type")+
  ylim(-1,1)+
  theme(legend.title = element_text(size=12, face = "bold"),
        legend.text = element_text(size=12), 
        axis.title.y = element_text(size = 14, angle = 90),
        axis.text.x = element_text(color=color_axis_x, size=12),
        axis.text.y = element_text(color="black", size=12),
        axis.line.y = element_line(color="black"))

p_bars

ggsave(p_bars, file="Corr_nscore_ESM1b_postype_bars.jpg", width = 15, height = 3, path=path)
ggsave(p_bars, file="Corr_nscore_ESM1b_postype_bars.pdf", width = 15, height = 3, path=path)

##################################################################################
# same as for GEMME and ESM1b but for Alphamissense
#Load alphamissense data

alphamissense_score<-read_csv("RIPK1_alphamissense.csv")
alphamissense_score<-rename(alphamissense_score, change = variant, alphamissense_score = score)
alphamissense_score<-left_join(alphamissense_score, AA_type[,c("wt_aa", "type")], by=c("WT_AA"="wt_aa"))
alphamissense_score<-rename(alphamissense_score, type_wt = type)
alphamissense_score<-left_join(alphamissense_score, AA_type[,c("wt_aa", "type")], by=c("Mut"="wt_aa"))
alphamissense_score<-rename(alphamissense_score, type_mut = type)

alphamissense_score$change_type<-paste0(alphamissense_score$type_wt, alphamissense_score$pos, alphamissense_score$type_mut)

RIPK1_alphamissense_nscore<-inner_join(singles_stops, alphamissense_score, by="change")
RIPK1_alphamissense_nscore<-rename(RIPK1_alphamissense_nscore, "WT_AA" = "WT_AA.x", "Mut" = "Mut.x")
RIPK1_alphamissense_nscore<-RIPK1_alphamissense_nscore[,c("change", "WT_AA", "Pos", "Mut", "nscore_c", "alphamissense_score")]
RIPK1_alphamissense_nscore$Residue<-paste0(RIPK1_alphamissense_nscore$WT_AA, RIPK1_alphamissense_nscore$Pos)
RIPK1_alphamissense_nscore$RHIM<-FALSE
RIPK1_alphamissense_nscore[(RIPK1_alphamissense_nscore$Pos>=531) & (RIPK1_alphamissense_nscore$Pos<=547),]$RHIM<-TRUE

RIPK1_alphamissense_nscore$Pos<-as.integer(RIPK1_alphamissense_nscore$Pos)
RIPK1_alphamissense_nscore<-arrange(RIPK1_alphamissense_nscore, Pos)

### Correlation alphamissense score vs NScore
correlation_all<-cor.test(RIPK1_alphamissense_nscore$nscore_c, RIPK1_alphamissense_nscore$alphamissense_score, use="complete.obs")

correlation_rhim<-cor.test(RIPK1_alphamissense_nscore[RIPK1_alphamissense_nscore$RHIM == TRUE,]$nscore_c, 
                           RIPK1_alphamissense_nscore[RIPK1_alphamissense_nscore$RHIM == TRUE,]$alphamissense_score, use="complete.obs")

correlation_out<-cor.test(RIPK1_alphamissense_nscore[RIPK1_alphamissense_nscore$RHIM == FALSE,]$nscore_c, 
                          RIPK1_alphamissense_nscore[RIPK1_alphamissense_nscore$RHIM == FALSE,]$alphamissense_score, use="complete.obs")

p_corr_all<-ggplot(RIPK1_alphamissense_nscore, aes(x=nscore_c, y=alphamissense_score))+
  geom_hline(yintercept = 0, color="black", linewidth=.2)+
  geom_vline(xintercept = 0, color="black", linewidth=.2)+
  geom_point(aes(color=RHIM), alpha=.6)+
  scale_color_manual(values=c("grey70", "coral"))+
  theme_classic()+
  theme(panel.grid = element_blank())+
  annotate("text", label=paste0("R=", round(correlation_all$estimate,2)), x=-Inf, y=Inf, hjust=-0.5, vjust=1, size=4, color="black")+
  annotate("text", label=paste0("p=", format(correlation_all$p.value, digits = 2, scientific = T)), x=-Inf, y=Inf, hjust=-0.2, vjust=2, size=4, color="black")+
  annotate("text", label=paste0("R=", round(correlation_rhim$estimate,2)), x=-Inf, y=Inf,hjust=-0.5, vjust=3.5, size=4, color="coral")+
  annotate("text", label=paste0("p=", format(correlation_rhim$p.value, digits = 2, scientific = T)), x=-Inf, y=Inf,hjust=-0.2, vjust=4.5, size=4, color="coral")+
  annotate("text", label=paste0("R=", round(correlation_out$estimate,2)), x=-Inf, y=Inf,hjust=-0.5, vjust=6, size=4, color="grey60")+
  annotate("text", label=paste0("p=", format(correlation_out$p.value, digits = 2, scientific = T)), x=-Inf, y=Inf,hjust=-0.2, vjust=7, size=4, color="grey60")+
  labs(x="Nucleation Score", y="AlphaMissense Score")
p_corr_all

ggsave(p_corr_all, file="Corr_nscore_alphamissense.jpg", width = 5, height = 4, path=path)
ggsave(p_corr_all, file="Corr_nscore_alphamissense.pdf", width = 5, height = 4, path=path)

# Correlation alphamissense score vs NScore per position
corr_vector=c()
for(i in RIPK1_alphamissense_nscore$Pos){
  correlation<-cor.test(RIPK1_alphamissense_nscore[RIPK1_alphamissense_nscore$Pos==i,]$nscore_c,
                        RIPK1_alphamissense_nscore[RIPK1_alphamissense_nscore$Pos==i,]$alphamissense_score, use="complete.obs")
  corr<-correlation$estimate
  p_value<-correlation$p.value
  corr_vector=c(corr_vector, i, corr, p_value)
}

corr_text <- data.frame(matrix(corr_vector, ncol=3, byrow = T))
colnames(corr_text)<-c("Pos", "corr","pvalue")
for(i in c(2:3)){corr_text[[i]]<-as.numeric(as.character(corr_text[[i]]))}
corr_text<-distinct(corr_text, Pos, .keep_all = TRUE)

#
p_corr_pos<-ggplot(RIPK1_alphamissense_nscore, aes(x=nscore_c, y=alphamissense_score))+
  geom_rect(data = RIPK1_alphamissense_nscore, aes(fill = RHIM), xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.05)+
  scale_fill_manual(values=c("white", "mistyrose"))+
  geom_hline(yintercept = 0, linetype="dashed", color="grey")+
  geom_vline(xintercept = 0, linetype="dashed", color="grey")+
  geom_point(color="grey20")+
  theme_bw()+
  facet_wrap(~factor(Pos, labels=distinct(RIPK1_alphamissense_nscore, Residue)[[1]]))+
  theme(panel.grid = element_blank())+
  geom_text(data=corr_text, aes(label=paste0("R=",round(corr, 2)), x=-Inf, y=Inf), hjust=-0.1, vjust=1.5, size=4, colour="red")+
  geom_text(data=corr_text, aes(label=paste0("p=", format(pvalue, digits = 2, scientific = T)), x=-Inf, y=Inf), hjust=-0.1, vjust=3, size=4, colour="red")+
  labs(x="Nucleation Score", y="AlphaMissense Score")
p_corr_pos

ggsave(p_corr_pos, file="Corr_nscore_alphamissense_pos.jpg", width = 12, height = 10, path=path)
ggsave(p_corr_pos, file="Corr_nscore_alphamissense_pos.pdf", width = 12, height = 10, path=path)

#
corr_text[is.na(corr_text)]<-1

corr_text$significance_pos<-""
corr_text[(corr_text$pvalue<0.05) & (corr_text$corr>0), "significance_pos"]<-"*"
corr_text[(corr_text$pvalue<0.01) & (corr_text$corr>0), "significance_pos"]<-"**"

corr_text$significance_neg<-""
corr_text[(corr_text$pvalue<0.05) & (corr_text$corr<0), "significance_neg"]<-"*"
corr_text[(corr_text$pvalue<0.01) & (corr_text$corr<0), "significance_neg"]<-"**"

p_bars<-ggplot(corr_text, aes(x=factor(Pos, labels=peptide_seq_pos), y=corr))+
  geom_hline(yintercept = 0, color="black", linewidth=.5)+
  geom_bar(stat="identity", color="black", fill="grey")+
  geom_text(data=corr_text, aes(label=significance_pos), vjust=0.4, size=8, colour="black")+
  geom_text(data=corr_text, aes(label=significance_neg), vjust=1, size=8, colour="black")+
  theme_void()+
  labs(x="RIPK1", y="Pearson Correlation", title= "NS vs AlphaMissense")+
  ylim(-1,1)+
  theme(legend.title = element_text(size=12, face = "bold"),
        legend.text = element_text(size=12), 
        axis.title.y = element_text(size = 14, angle = 90),
        axis.text.x = element_text(color=color_axis_x, size=12),
        axis.text.y = element_text(color="black", size=12),
        axis.line.y = element_line(color="black"))

p_bars

ggsave(p_bars, file="Corr_nscore_alphamissense_bars.jpg", width = 15, height = 3, path=path)
ggsave(p_bars, file="Corr_nscore_alphamissense_bars.pdf", width = 15, height = 3, path=path)

# alphamissense heatmap

alphamissense_score<-alphamissense_score[alphamissense_score$pos <= end_p & alphamissense_score$pos >= start_p,]
alphamissense_score$Residue<-substr(alphamissense_score$change, 1, 4)
alphamissense_score$Mut<-substr(alphamissense_score$change, 5, 5)
alphamissense_score$Pos<-as.integer(alphamissense_score$pos)

p_alphamissense_heatmap<-ggplot(alphamissense_score)+
  geom_tile(aes(Pos, factor(Mut, levels=rev(vectorAA)), fill=alphamissense_score), color='white', size=1)+
  scale_x_continuous(breaks=seq(start_p, end_p), labels = peptide_seq_pos, expand=c(0,0))+
  #scale_fill_gradient2(low="gold", high="purple", mid="white")+
  scale_fill_viridis(direction=-1, option="E")+
  theme_minimal()+
  labs(x="RIPK1", y="Mutant amino acid", fill="AlphaMissense score")+
  theme(legend.title = element_text(size=12),
        legend.text = element_text(size=12), 
        axis.title = element_text(size = 14),
        axis.text.x = element_text(color=c(rep("black", 8), rep("red", 17), rep("black", 9)), size=12),
        axis.text.y = element_text(color="black", size=12))
p_alphamissense_heatmap

ggsave(p_alphamissense_heatmap, file="alphamissense_heatmap.jpg", width=15, height=8, path=path)

#################################################################################
### Correlation per position considering amino acid type

singles_stops_group<-as.data.frame(singles_stops %>% group_by(Pos, Residue, change_type) %>% dplyr::summarise(mean_ns=mean(nscore_c)))

singles_alphamissense_group<-as.data.frame(alphamissense_score %>% group_by(pos, residue, change_type) %>% dplyr::summarise(mean_alphamissense=mean(alphamissense_score)))

singles_nscore_alphamissense<-inner_join(singles_stops_group, singles_alphamissense_group, by="change_type")
singles_nscore_alphamissense<-singles_nscore_alphamissense[c("change_type", "Pos", "mean_ns", "mean_alphamissense")]
singles_nscore_alphamissense<-arrange(singles_nscore_alphamissense, Pos)
singles_nscore_alphamissense$RHIM<-FALSE
singles_nscore_alphamissense[(singles_nscore_alphamissense$Pos>=531) & (singles_nscore_alphamissense$Pos<=547),]$RHIM<-TRUE


correlation_all<-cor.test(singles_nscore_alphamissense$mean_ns, singles_nscore_alphamissense$mean_alphamissense, use="complete.obs")

correlation_rhim<-cor.test(singles_nscore_alphamissense[singles_nscore_alphamissense$RHIM == TRUE,]$mean_ns, 
                           singles_nscore_alphamissense[singles_nscore_alphamissense$RHIM == TRUE,]$mean_alphamissense, use="complete.obs")

correlation_out<-cor.test(singles_nscore_alphamissense[singles_nscore_alphamissense$RHIM == FALSE,]$mean_ns, 
                          singles_nscore_alphamissense[singles_nscore_alphamissense$RHIM == FALSE,]$mean_alphamissense, use="complete.obs")

p_corr_all<-ggplot(singles_nscore_alphamissense, aes(x=mean_ns, y=mean_alphamissense))+
  geom_hline(yintercept = 0, color="black", linewidth=.2)+
  geom_vline(xintercept = 0, color="black", linewidth=.2)+
  geom_point(aes(color=RHIM), alpha=.6)+
  scale_shape_manual(values=c(19, 6))+
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
  labs(x="Nucleation Score", y="AlphaMissense Score")
p_corr_all

ggsave(p_corr_all, file="Corr_singles_nscore_postype_RIPK1_alphamissense.jpg", width = 5, height = 4, path=path)
ggsave(p_corr_all, file="Corr_singles_nscore_postype_RIPK1_alphamissense.pdf", width = 5, height = 4, path=path)

# correlation per position
corr_vector=c()
for(i in singles_nscore_alphamissense$Pos){
  correlation<-cor.test(singles_nscore_alphamissense[singles_nscore_alphamissense$Pos==i,][["mean_ns"]],
                        singles_nscore_alphamissense[singles_nscore_alphamissense$Pos==i,][["mean_alphamissense"]], use="complete.obs")
  corr<-correlation$estimate
  p_value<-correlation$p.value
  corr_vector=c(corr_vector, i, corr, p_value)
}

corr_text <- data.frame(matrix(corr_vector, ncol=3, byrow = T))
colnames(corr_text)<-c("Pos", "corr","pvalue")
for(i in c(2:3)){corr_text[[i]]<-as.numeric(as.character(corr_text[[i]]))}
corr_text<-distinct(corr_text, Pos, .keep_all = TRUE)

#
p_corr_pos<-ggplot(singles_nscore_alphamissense, aes(x=mean_ns, y=mean_alphamissense))+
  geom_rect(data = singles_nscore_alphamissense, aes(fill = RHIM), xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.05)+
  scale_fill_manual(values=c("white", "mistyrose"))+
  geom_hline(yintercept = 0, linetype="dashed", color="darkgrey")+
  geom_vline(xintercept = 0, linetype="dashed", color="darkgrey")+
  geom_point()+
  theme_bw()+
  facet_wrap(~factor(Pos, label = unique(singles_stops$Residue)))+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  geom_text(data=corr_text, aes(label=paste0("R=",round(corr, 2)), x=-Inf, y=Inf), hjust=-0.05, vjust=1.5, size=4, colour="red")+
  geom_text(data=corr_text, aes(label=paste0("p=", format(pvalue, digits = 2, scientific = T)), x=-Inf, y=Inf), hjust=-0.05, vjust=3, size=4, colour="red")+
  scale_color_manual(name = "Residue", labels = c("Mouse", "Human"), values=c("brown", "orange"))+
  labs(x="Nucleation Score", y="AlphaMissense Score")
p_corr_pos

ggsave(p_corr_pos, file="Corr_singles_nscore_postypeeach_RIPK1_alphamissense.jpg", width = 12, height = 10, path=path)
ggsave(p_corr_pos, file="Corr_singles_nscore_postypeeach_RIPK1_alphamissense.pdf", width = 12, height = 10, path=path)

# bar plot
corr_text[is.na(corr_text)]<-1

corr_text$significance_pos<-""
corr_text[(corr_text$pvalue<0.05) & (corr_text$corr>0), "significance_pos"]<-"*"
corr_text[(corr_text$pvalue<0.01) & (corr_text$corr>0), "significance_pos"]<-"**"

corr_text$significance_neg<-""
corr_text[(corr_text$pvalue<0.05) & (corr_text$corr<0), "significance_neg"]<-"*"
corr_text[(corr_text$pvalue<0.01) & (corr_text$corr<0), "significance_neg"]<-"**"

p_bars<-ggplot(corr_text, aes(x=factor(Pos, labels=peptide_seq_pos), y=corr))+
  geom_hline(yintercept = 0, color="black", linewidth=.5)+
  geom_bar(stat="identity", color="black", fill="grey")+
  geom_text(data=corr_text, aes(label=significance_pos), vjust=0.4, size=8, colour="black")+
  geom_text(data=corr_text, aes(label=significance_neg), vjust=1, size=8, colour="black")+
  theme_void()+
  labs(x="RIPK1", y="Pearson Correlation", title="NS vs AlphaMissense - per aa type")+
  ylim(-1,1)+
  theme(legend.title = element_text(size=12, face = "bold"),
        legend.text = element_text(size=12), 
        axis.title.y = element_text(size = 14, angle = 90),
        axis.text.x = element_text(color=color_axis_x, size=12),
        axis.text.y = element_text(color="black", size=12),
        axis.line.y = element_line(color="black"))

p_bars

ggsave(p_bars, file="Corr_nscore_alphamissense_postype_bars.jpg", width = 15, height = 3, path=path)
ggsave(p_bars, file="Corr_nscore_alphamissense_postype_bars.pdf", width = 15, height = 3, path=path)

###
# let's plot median score per position
RIPK1_VEPs<-inner_join(RIPK1_alphamissense_nscore, RIPK1_ESM1b_nscore[c("change", "ESM1b_score")], by="change")
RIPK1_VEPs<-inner_join(RIPK1_VEPs, RIPK1_gemme_nscore[c("change", "GEMME_score")], by="change")

RIPK1_VEPs<- RIPK1_VEPs %>% group_by(Residue) %>% mutate("median_nscore"=median(nscore_c),
                                                         "median_alphamissense"=median(alphamissense_score), 
                                                         "median_ESM1b"=median(ESM1b_score),
                                                         "median_GEMME"=median(GEMME_score))

RIPK1_VEPs_resmedian<-RIPK1_VEPs[c("Residue", "median_nscore", "median_alphamissense", "median_ESM1b", "median_GEMME")]
RIPK1_VEPs_resmedian<-RIPK1_VEPs_resmedian[!duplicated(RIPK1_VEPs_resmedian),]

RIPK1_VEPs_median_pivot <- 
  RIPK1_VEPs_resmedian %>%
  pivot_longer(
    cols = !Residue,
    names_to = "score_name",
    values_to = "score"
  )

## Median NS
# color scale to be used
min<-min(RIPK1_VEPs_median_pivot[RIPK1_VEPs_median_pivot$score_name == "median_nscore", "score"])
max<-max(RIPK1_VEPs_median_pivot[RIPK1_VEPs_median_pivot$score_name == "median_nscore", "score"])
cols <- c(colorRampPalette(c( "brown3", "grey95"))((-min/(-min+max)*100)-0.5), colorRampPalette("grey95")(1),
          colorRampPalette(c("grey95",  "darkblue"), bias=1)((max/(-min+max)*100)-0.5))

p_heatmap_ns<-ggplot(RIPK1_VEPs_median_pivot[RIPK1_VEPs_median_pivot$score_name == "median_nscore",])+
  geom_tile(aes(x=factor(Residue, levels=unique(RIPK1_VEPs_median_pivot$Residue)), y=score_name, fill=score), color='white', size=1)+
  scale_fill_gradientn(colours=cols, limits=c(min,max))+
  labs(fill="Nucleation Score", x="", y="")+
  theme_minimal()+
  theme(legend.title = element_text(size=12),
        legend.text = element_text(size=12), 
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(color="black", size=12))
p_heatmap_ns

## Median alphamissense
# color scale to be used
min<-min(RIPK1_VEPs_median_pivot[RIPK1_VEPs_median_pivot$score_name == "median_alphamissense", "score"])
max<-max(RIPK1_VEPs_median_pivot[RIPK1_VEPs_median_pivot$score_name == "median_alphamissense", "score"])


p_heatmap_alphamissense<-ggplot(RIPK1_VEPs_median_pivot[RIPK1_VEPs_median_pivot$score_name == "median_alphamissense",])+
  geom_tile(aes(x=factor(Residue, levels=unique(RIPK1_VEPs_median_pivot$Residue)), y=score_name, fill=score), color='white', size=1)+
  scale_fill_gradient(low = "gold3", high = "grey95", limits=c(min,max))+
  labs(fill="Alphamissense Score", x="", y="")+
  theme_minimal()+
  theme(legend.title = element_text(size=12),
        legend.text = element_text(size=12), 
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(color="black", size=12))
p_heatmap_alphamissense

## Median ESM1b
# color scale to be used
min<-min(RIPK1_VEPs_median_pivot[RIPK1_VEPs_median_pivot$score_name == "median_ESM1b", "score"])
max<-max(RIPK1_VEPs_median_pivot[RIPK1_VEPs_median_pivot$score_name == "median_ESM1b", "score"])

p_heatmap_ESM1b<-ggplot(RIPK1_VEPs_median_pivot[RIPK1_VEPs_median_pivot$score_name == "median_ESM1b",])+
  geom_tile(aes(x=factor(Residue, levels=unique(RIPK1_VEPs_median_pivot$Residue)), y=score_name, fill=score), color='white', size=1)+
  scale_fill_gradient(low = "orange1", high = "grey95", limits=c(min,max))+
  labs(fill="ESM1b Score", x="", y="")+
  theme_minimal()+
  theme(legend.title = element_text(size=12),
        legend.text = element_text(size=12), 
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(color="black", size=12))
p_heatmap_ESM1b

## Median GEMME
# color scale to be used
min<-min(RIPK1_VEPs_median_pivot[RIPK1_VEPs_median_pivot$score_name == "median_GEMME", "score"])
max<-max(RIPK1_VEPs_median_pivot[RIPK1_VEPs_median_pivot$score_name == "median_GEMME", "score"])

p_heatmap_GEMME<-ggplot(RIPK1_VEPs_median_pivot[RIPK1_VEPs_median_pivot$score_name == "median_GEMME",])+
  geom_tile(aes(x=factor(Residue, levels=unique(RIPK1_VEPs_median_pivot$Residue)), y=score_name, fill=score), color='white', size=1)+
  scale_fill_gradient(low = "brown", high = "grey95", limits=c(min,max))+
  labs(fill="GEMME Score", x="", y="")+
  theme_minimal()+
  theme(legend.title = element_text(size=12),
        legend.text = element_text(size=12), 
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(color="black", size=12))
p_heatmap_GEMME

# arrange all the median heatmaps
arrange_heatmap <- ggarrange(p_heatmap_ns + rremove("x.text"), p_heatmap_GEMME + rremove("x.text"), p_heatmap_ESM1b + rremove("x.text"), p_heatmap_alphamissense, 
                             ncol=1, nrow=4, align = "v")

arrange_heatmap

ggsave(arrange_heatmap, file="VEPs_median_heatmap.jpg", width=15, height = 7, path=path)
ggsave(arrange_heatmap, file="VEPs_median_heatmap.pdf", width=15, height = 7, path=path)
