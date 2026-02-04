library(tidyverse)

# your folder path
setwd("")

# Create a folder for the figures
dir.create("04_mRIP3vsTheOthers")
path="04_mRIP3vsTheOthers"

# Load ripk3 data
load("nscore_df_hRIPK3.RData")

# library start and end positions
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
singles_ripk3$ID_ripk3<-singles_ripk3$ID
singles_ripk3$residue_ripk3<-paste0(singles_ripk3$WT_AA, singles_ripk3$Pos)

singles_ripk3<-singles_ripk3[c("Pos_c", "residue_ripk3", "wt_ripk3", "pos_ripk3", "mut_ripk3", "ID_ripk3", "nscore_c_ripk3")]
singles_ripk3$change<-paste0(singles_ripk3$Pos_c, singles_ripk3$mut_ripk3)

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
singles_mRIP3$ID_mRIP3<-singles_mRIP3$ID
singles_mRIP3$residue_mRIP3<-paste0(singles_mRIP3$WT_AA, singles_mRIP3$Pos)

singles_mRIP3<-singles_mRIP3[c("Pos_c", "residue_mRIP3", "wt_mRIP3", "pos_mRIP3", "mut_mRIP3", "ID_mRIP3", "nscore_c_mRIP3")]
singles_mRIP3$change<-paste0(singles_mRIP3$Pos_c, singles_mRIP3$mut_mRIP3)

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
singles_ripk1$ID_ripk1<-singles_ripk1$ID
singles_ripk1$residue_ripk1<-paste0(singles_ripk1$WT_AA, singles_ripk1$Pos)

singles_ripk1<-singles_ripk1[c("Pos_c", "residue_ripk1", "wt_ripk1", "pos_ripk1", "mut_ripk1", "ID_ripk1", "nscore_c_ripk1")]
singles_ripk1$change<-paste0(singles_ripk1$Pos_c, singles_ripk1$mut_ripk1)

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
singles_mRIP1$ID_mRIP1<-singles_mRIP1$ID
singles_mRIP1$residue_mRIP1<-paste0(singles_mRIP1$WT_AA, singles_mRIP1$Pos)

singles_mRIP1<-singles_mRIP1[c("Pos_c", "residue_mRIP1", "wt_mRIP1", "pos_mRIP1", "mut_mRIP1", "ID_mRIP1", "nscore_c_mRIP1")]
singles_mRIP1$change<-paste0(singles_mRIP1$Pos_c, singles_mRIP1$mut_mRIP1)

################################################################################
singles_all<-inner_join(singles_ripk3, singles_mRIP3, by="change")
singles_all<-inner_join(singles_all, singles_ripk1, by="change")
singles_all<-inner_join(singles_all, singles_mRIP1, by="change")

# Let's see how's the correlation between equivalent position E vs G
singles_mRIP3_E447<-singles_all[singles_all$residue_mRIP3 == "E447",]

corr_vector=c()
# vs RIPK3
correlation<-cor.test(singles_mRIP3_E447$nscore_c_mRIP3,
                      singles_mRIP3_E447$nscore_c_ripk3, use="complete.obs", method="pearson")
corr<-correlation$estimate
p_value<-correlation$p.value
corr_vector=c(corr_vector, "vsRIPK3", corr, p_value)
# vs RIPK1
correlation<-cor.test(singles_mRIP3_E447$nscore_c_mRIP3,
                      singles_mRIP3_E447$nscore_c_ripk1, use="complete.obs", method="pearson")
corr<-correlation$estimate
p_value<-correlation$p.value
corr_vector=c(corr_vector, "vsRIPK1", corr, p_value)
# vs mRIP1
correlation<-cor.test(singles_mRIP3_E447$nscore_c_mRIP3,
                      singles_mRIP3_E447$nscore_c_mRIP1, use="complete.obs", method="pearson")
corr<-correlation$estimate
p_value<-correlation$p.value
corr_vector=c(corr_vector, "vsmRIP1", corr, p_value)

corr_text <- data.frame(matrix(corr_vector, ncol=3, byrow = T))
colnames(corr_text)<-c("sequence", "corr","pvalue")
for(i in c(2:3)){corr_text[[i]]<-as.numeric(as.character(corr_text[[i]]))}
corr_text<-distinct(corr_text, sequence, .keep_all = TRUE)

corr_text[is.na(corr_text)]<-1

corr_text$significance_pos<-""
corr_text[(corr_text$pvalue<0.05) & (corr_text$corr>0), "significance_pos"]<-"*"
corr_text[(corr_text$pvalue<0.01) & (corr_text$corr>0), "significance_pos"]<-"**"
corr_text[(corr_text$pvalue<0.001) & (corr_text$corr>0), "significance_pos"]<-"***"

corr_text$significance_neg<-""
corr_text[(corr_text$pvalue<0.05) & (corr_text$corr<0), "significance_neg"]<-"*"
corr_text[(corr_text$pvalue<0.01) & (corr_text$corr<0), "significance_neg"]<-"**"
corr_text[(corr_text$pvalue<0.001) & (corr_text$corr<0), "significance_neg"]<-"***"


p_bars<-ggplot(corr_text, aes(x=factor(sequence, levels=c("vsRIPK3", "vsRIPK1", "vsmRIP1"), labels=c("RIPK3", "RIPK1", "mRIP1")), 
                              y=corr))+
  geom_hline(yintercept = 0, color="black", linewidth=.5)+
  geom_bar(stat="identity", color="black", fill="grey", width=.5)+
  geom_text(data=corr_text, aes(label=significance_pos), vjust=0.4, size=8, colour="black")+
  geom_text(data=corr_text, aes(label=significance_neg), vjust=1, size=8, colour="black")+
  theme_void()+
  labs(x="", y="Pearson Correlation")+
  theme_classic()

p_bars

ggsave(p_bars, file="E447vsG_allRIPKs.jpg", width = 4, height = 3, path=path)
ggsave(p_bars, file="E447vsG_allRIPKs.pdf", width = 4, height = 3, path=path)
