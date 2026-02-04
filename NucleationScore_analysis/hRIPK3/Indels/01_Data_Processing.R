library(tidyverse)
library(ggpubr)
library(readxl)

# your folder path
setwd("")

# Create a folder for the figures
dir.create("01_Data_Processing")
path="01_Data_Processing"

# Open the dimsum output with the fitness estimates
name<-"hRIPK3_indels"
load(paste0(name, '_fitness_replicates.RData'))

# Read desgined library file
indels_library<-read_excel("RIPK3_library.xlsx")

# Read dimsum output file
variant_data_merge_df<-read_tsv(paste0(name,"_variant_data_merge.tsv"))

###
all_variants_df<-all_variants

###

all_variants_df<-left_join(all_variants_df, indels_library, by=c("aa_seq"="AA_sequence"))

#centering to synonymous
synonymous$dataset<-"synonymous"
synonymous<-rename(synonymous, nscore_c = fitness)
mean_syn_1codon<-weighted.mean(synonymous$nscore_c, synonymous$sigma^-2, na.rm = T)

all_variants_df$nscore_c<-as.numeric(all_variants_df$fitness-mean_syn_1codon)
all_variants_df$nscore1_c<-as.numeric(all_variants_df$fitness1_uncorr-mean_syn_1codon)
all_variants_df$nscore2_c<-as.numeric(all_variants_df$fitness2_uncorr-mean_syn_1codon)
all_variants_df$nscore3_c<-as.numeric(all_variants_df$fitness3_uncorr-mean_syn_1codon)

all_variants_df$WT<-replace_na(all_variants_df$WT, FALSE)
all_variants_df[all_variants_df$WT == TRUE,]$dataset<-"WT"
#all_variants_df[all_variants_df$aa_seq %in% singles[singles$STOP == TRUE,]$aa_seq, ]$dataset<-"STOP"

all_variants_df<-all_variants_df[,c("nt_seq", "aa_seq", "ID", "dataset", "mean_count", "nscore_c", "nscore1_c", "nscore2_c","nscore3_c", "sigma")]

all_variants_df$dataset <- all_variants_df$dataset %>% replace_na("not_designed")

# Keep alive designed sequences
singles_indels_df<-drop_na(all_variants_df, ID)

###

# find non-nucleating variants (they have input reads but 0 output reads- NS not calculated)
# each AA variant is only resulting from one codon change by design (non nuc nt seq results in non nuc aa seq)
n_reads<-200
non_nuc_df<-variant_data_merge_df[variant_data_merge_df$input1_e1_s0_bNA_count >= n_reads &
                                    variant_data_merge_df$input2_e2_s0_bNA_count >= n_reads &
                                    variant_data_merge_df$input3_e3_s0_bNA_count >= n_reads &
                                    (variant_data_merge_df$output1_e1_s1_b1_count == 0 &
                                       variant_data_merge_df$output2_e2_s1_b1_count == 0 &
                                       variant_data_merge_df$output3_e3_s1_b1_count == 0),]
### There are no variants like this

# test variants against WT at FDR=0.01
singles_indels_df$zscore<-singles_indels_df$nscore_c/singles_indels_df$sigma
singles_indels_df$p.adjust<-p.adjust(2*pnorm(-abs(singles_indels_df$zscore)), method = "BH")

singles_indels_df$sig_fdr<-FALSE
singles_indels_df[singles_indels_df$p.adjust<0.1,]$sig_fdr<-TRUE

singles_indels_df$category_fdr<-"WT-like"
singles_indels_df[singles_indels_df$sig_fdr==T & singles_indels_df$nscore_c<0,]$category_fdr<-"NS_dec"
singles_indels_df[singles_indels_df$sig_fdr==T & singles_indels_df$nscore_c>0,]$category_fdr<-"NS_inc"

###

# interquartile range == fitness range

summary(singles_indels_df$nscore_c)
iqr<-IQR(singles_indels_df$nscore_c)

first_to_wt<-summary(singles_indels_df$nscore_c)[[2]]

p_iqr<-ggplot(singles_indels_df, aes(x=nscore_c))+
  geom_histogram(color="black", fill="grey", bins=100)+
  theme_bw()+
  labs(x="Nucleation score", title="Fitness range")+
  
  geom_vline(xintercept=0)+
  annotate("text", label="WT fitness= 0", x=0.55, y=140, color="black")+
  
  geom_vline(xintercept=summary(singles_indels_df$nscore_c)[[2]], color="red")+
  annotate("text", label=paste0("1st Qu= ", round(summary(singles_indels_df$nscore_c)[[2]], 2)), 
           x=summary(singles_indels_df$nscore_c)[[2]]-0.8, y=150, color="red")+
  
  geom_vline(xintercept=summary(singles_indels_df$nscore_c)[[5]], color="red")+
  annotate("text", label=paste0("3rd Qu= ", round(summary(singles_indels_df$nscore_c)[[5]], 2)), 
           x=summary(singles_indels_df$nscore_c)[[5]]+0.8, y=150, color="red")+
  
  geom_vline(xintercept=summary(singles_indels_df$nscore_c)[[3]], color="blue")+
  annotate("text", label=paste0("median= ",round(summary(singles_indels_df$nscore_c)[[3]], 2) ), 
           x=summary(singles_indels_df$nscore_c)[[3]]-0.55, y=150, color="blue")
p_iqr

ggsave(p_iqr, file="p_iqr.jpg", path=path, width = 5, height = 3)

# sigma distribution
p_sigma_dist<-ggplot(singles_indels_df, aes(x=sigma))+
  geom_histogram(color="black", fill="grey", bins=100)+
  theme_bw()+
  scale_x_continuous(limits = c(0,3))
p_sigma_dist

ggsave(p_sigma_dist, file="p_sigma_dist.jpg", path=path, width = 5, height = 3)

# normalise sigmas to fitness range 
singles_indels_df$sigma_norm_iqr<-""
fitness_range_iqr=abs(IQR(singles_indels_df$nscore_c))

# or if fitness range is 1rst to WT fitness
singles_indels_df$sigma_norm_first_toWT<-""
fitness_range_first_toWT=abs(summary(singles_indels_df$nscore_c)[[2]])

singles_indels_df$sigma_norm_iqr<-singles_indels_df$sigma / fitness_range_iqr
singles_indels_df$sigma_norm_first_toWT<-singles_indels_df$sigma / fitness_range_first_toWT

singles_indels_df$sigma_norm_iqr<-as.numeric(singles_indels_df$sigma_norm_iqr)
singles_indels_df$sigma_norm_first_toWT<-as.numeric(singles_indels_df$sigma_norm_first_toWT)

#sigma_norm_iqr

# sigma_norm_iqr
p1<-ggplot(singles_indels_df, aes(x=sigma_norm_iqr))+
  geom_histogram(bins=200, aes(fill=factor(category_fdr, levels=c("NS_inc", "NS_dec", "WT-like"))))+
  scale_fill_manual(values=c("#7979BE", "#DF9292", "#F2F2F2"))+
  facet_wrap(~category_fdr, ncol=1)+
  labs(title="sigma normalised to fitness range (1st Qu to 3r Qu)", fill="Category_FDR")+
  scale_x_continuous(limits = c(0,1))+
  geom_vline(xintercept = 0.3, linetype="dashed", color="grey20")+
  theme_classic()
p1

ggsave(p1, file="p_sigma_normalised.jpg", path=path, width = 6, height = 6)

## exclude those that are over 30% of the fitness range
singles_indels_df$low_sigma<-FALSE
singles_indels_df[singles_indels_df$sigma_norm_iqr<=0.3,]$low_sigma<-TRUE

# Categories based on sigma normalized to fitness range
singles_indels_df$category_fdr_sigma<-"unclassified"

singles_indels_df[singles_indels_df$category_fdr=="NS_inc"& singles_indels_df$low_sigma==T,]$category_fdr_sigma<-"NS_inc"
singles_indels_df[singles_indels_df$category_fdr=="NS_dec" & singles_indels_df$low_sigma==T,]$category_fdr_sigma<-"NS_dec"
singles_indels_df[singles_indels_df$category_fdr=="WT-like" & singles_indels_df$low_sigma==T,]$category_fdr_sigma<-"WT-like"

###

#divide into mutation types: single AA substitutions, single AA insertions, single AA deletions, truncations and internal deletions
singles.df<-singles_indels_df[singles_indels_df$dataset=="missense",]
insertions.df<-singles_indels_df[singles_indels_df$dataset=="insertion",]
deletions.df<-singles_indels_df[singles_indels_df$dataset=="deletion",]

###########
# singles #
###########

singles.df$WT_AA<-""
singles.df$Pos<-""
singles.df$Mut<-""

for(i in 1:nrow(singles.df)){
  singles.df[i,]$WT_AA<-unlist(strsplit(singles.df[i,]$ID, "_"))[2]
  singles.df[i,]$Pos<-unlist(strsplit(singles.df[i,]$ID, "_"))[3]
  singles.df[i,]$Mut<-unlist(strsplit(singles.df[i,]$ID, "_"))[4]
}

singles.df$Pos<-as.numeric(singles.df$Pos)

##############
# insertions #
##############

# one coding sequence can come from more than one insertion- duplicate rows

insertions.df$ins_pos<-""
insertions.df$ins_aa<-""

for (i in 1:nrow(insertions.df)){
  
  if(is.na(str_extract(insertions.df[i,]$ID, ";"))){
    
    ID<-unlist(str_split(insertions.df[i,]$ID, "_"))
    
    insertions.df[i,]$ins_pos<-ID[4]
    insertions.df[i,]$ins_aa<-ID[5]
    
  }else{
    
    all_ID<-unlist(str_split(insertions.df[i,]$ID, ";"))
    num_ID<-length(all_ID)
    
    for(j in all_ID){
      
      new_row<- insertions.df[i,]
      new_row$ID<-j
      
      ID<-unlist(str_split(j, "_"))
      
      new_row$ins_pos<-ID[4]
      new_row$ins_aa<-ID[5]
      
      insertions.df<-rbind(insertions.df,new_row)
      
    }
  }
}

#remove rows with more than one ID (they are already split)
insertions.df<-insertions.df[is.na(str_extract(insertions.df$ID, ";")),]

#######################
# single AA deletions #
#######################

# one coding sequence can come from more than one single deletion- duplicate rows

deletions.df<-deletions.df

deletions.df$del_pos<-""

for (i in 1:nrow(deletions.df)){
  
  if(is.na(str_extract(deletions.df[i,]$ID, ";"))){
    
    ID<-unlist(str_split(deletions.df[i,]$ID, "_"))
    
    deletions.df[i,]$del_pos<-unlist(strsplit(ID[4],"-"))[1]
    
  }else{
    
    all_ID<-unlist(str_split(deletions.df[i,]$ID, ";"))
    
    for(j in all_ID){
      
      new_row<- deletions.df[i,]
      new_row$ID<-j
      
      ID<-unlist(str_split(j, "_"))
      
      new_row$del_pos<-unlist(strsplit(ID[4],"-"))[1]
      
      deletions.df<-rbind(deletions.df,new_row)
      
    }
  }
}

#remove rows with more than one ID (they are already split)
deletions.df<-deletions.df[is.na(str_extract(deletions.df$ID, ";")),]

###

summary_singles_indels_1<-as.data.frame(singles_indels_df %>% group_by(category_fdr) %>% dplyr::summarise(n=n())  %>% mutate(freq = n / sum(n) ))
summary_singles_indels_2<-as.data.frame(singles_indels_df %>% group_by(category_fdr_sigma) %>% dplyr::summarise(n=n())  %>% mutate(freq = n / sum(n) ))

###
# correct the residue numbering
singles.df$Pos<-singles.df$Pos+441
insertions.df$ins_pos<-as.numeric(insertions.df$ins_pos)
insertions.df$ins_pos<-insertions.df$ins_pos+441
deletions.df$del_pos<-as.numeric(deletions.df$del_pos)
deletions.df$del_pos<-deletions.df$del_pos+441

###

save(synonymous, all_variants_df, singles_indels_df, singles.df, insertions.df, deletions.df, file = paste0(name, "_singles.RData"))

save(non_nuc_df, file="non_nucleating_variants.RData")

write.table(all_variants_df, file="RIPK3_indels_processed_data.tsv", sep="\t", quote = F, row.names = F)
write.table(singles_indels_df, file="RIPK3_indels_processed_data_designed.tsv", sep="\t", quote = F, row.names = F)

