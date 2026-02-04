library(tidyverse)
library(ggpubr)

# your folder path
setwd("")

# Create a folder for the figures
dir.create("01_Centering_Correction")
path="01_Centering_Correction"

# Open the dimsum output with the fitness estimates
name<-'hRIPK1'
load(paste0(name, '_fitness_replicates.RData'))

# rename the fitness as nucleation score (nscore)
all_variants<-rename(all_variants, 
                     nscore = fitness, 
                     nscore1 = fitness1_uncorr, nscore2 = fitness2_uncorr, nscore3 = fitness3_uncorr)
                   
all_variants<-select(all_variants, aa_seq, starts_with('nscore') | starts_with('count'))

synonymous<-rename(synonymous, 
                   nscore = fitness, 
                   nscore1 = fitness1_uncorr, nscore2 = fitness2_uncorr, nscore3 = fitness3_uncorr)

singles<-select(singles, !starts_with('fitness'))
singles<-inner_join(singles, all_variants, by='aa_seq')

# correct the residue position
singles$Pos<-singles$Pos+522

# centering to the weighted mean of synonymous with 1 mut codon
silent<-synonymous
mean_syn_1codon<-weighted.mean(silent[silent$Nmut_codons==1,]$nscore, silent[silent$Nmut_codons==1,]$sigma^-2, na.rm = T)

# silents centering
silent$nscore_c<-as.numeric(silent$nscore-mean_syn_1codon)
silent$ID<-"silent"
silent$Mut<-"silent"

# singles centering
singles$nscore_c<-as.numeric(singles$nscore-mean_syn_1codon)
singles$nscore1_c<-as.numeric(singles$nscore1-mean_syn_1codon)
singles$nscore2_c<-as.numeric(singles$nscore2-mean_syn_1codon)
singles$nscore3_c<-as.numeric(singles$nscore3-mean_syn_1codon)

# Generate the variant ID
singles$ID<-paste(singles$WT_AA, singles$Pos, singles$Mut, sep = "-")

# FDR=0.01 correction and assignment into categories
singles$zscore<-singles$nscore_c/singles$sigma
singles$p.adjust<-p.adjust(2*pnorm(-abs(singles$zscore)), method = "BH")

singles$sig_10<-FALSE
singles[singles$p.adjust<0.1,]$sig_10<-TRUE

singles$category_10<-"WT-like"
singles[singles$sig_10==T & singles$nscore_c<0,]$category_10<-"NS_dec"
singles[singles$sig_10==T & singles$nscore_c>0,]$category_10<-"NS_inc"

# singles with stops
singles_stops<-singles

# singles without stops
singles<-singles[singles$Mut!="*",]

# mean input and ouput reads counts
singles_stops <- singles_stops %>%
  rowwise() %>%
  mutate(input_mean_count=mean(c(count_e1_s0, count_e2_s0, count_e3_s0)),
         output_mean_count=mean(c(count_e1_s1, count_e2_s1, count_e3_s1)))

#################################################################################
# Check how's the noise in this library
# Use the sigma normalize to the  interquartile range == fitness range

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

# sigma distribution
p_sigma_dist<-ggplot(singles_stops, aes(x=sigma))+
  geom_histogram(color="black", fill="grey", bins=100)+
  theme_classic()+
  scale_x_continuous(limits = c(0,3))
p_sigma_dist

ggsave(p_sigma_dist, file="p_sigma_dist.jpg", path=path, width = 5, height = 3)

# normalise sigmas to fitness range 
singles_stops$sigma_norm_iqr<-""
fitness_range_iqr=abs(IQR(singles_stops$nscore_c))

# or if fitness range is 1rst to WT fitness
singles_stops$sigma_norm_first_toWT<-""
fitness_range_first_toWT=abs(summary(singles_stops$nscore_c)[[2]])

singles_stops$sigma_norm_iqr<-singles_stops$sigma / fitness_range_iqr
singles_stops$sigma_norm_first_toWT<-singles_stops$sigma / fitness_range_first_toWT

singles_stops$sigma_norm_iqr<-as.numeric(singles_stops$sigma_norm_iqr)
singles_stops$sigma_norm_first_toWT<-as.numeric(singles_stops$sigma_norm_first_toWT)

# sigma_norm_iqr
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

# exclude those that are over 30% of the fitness range
singles_stops$low_sigma<-FALSE
singles_stops[singles_stops$sigma_norm_iqr<=0.3,]$low_sigma<-TRUE

# Categories based on sigma normalized to fitness range
singles_stops$category_10_sigma<-"unclassified"

singles_stops[singles_stops$category_10=="NS_inc"& singles_stops$low_sigma==T,]$category_10_sigma<-"NS_inc"
singles_stops[singles_stops$category_10=="NS_dec" & singles_stops$low_sigma==T,]$category_10_sigma<-"NS_dec"
singles_stops[singles_stops$category_10=="WT-like" & singles_stops$low_sigma==T,]$category_10_sigma<-"WT-like"

# save this R file
save(silent, singles, singles_stops, file="nscore_df_hRIPK1.RData")

################################################################################
# Check how many variants are present in the input and in the output
variant_data_merge<-read_tsv("hRIPK1_variant_data_merge.tsv")

variant_data_merge_singles<-variant_data_merge[variant_data_merge$Nham_aa == 1 & 
                                               variant_data_merge$indel == F ,]

input_output_reads<-function(n_reads){
  aa_unique_input<-unique(variant_data_merge_singles[variant_data_merge_singles$input1_e1_s0_bNA_count >= n_reads &
                                                     variant_data_merge_singles$input2_e2_s0_bNA_count >= n_reads &
                                                     variant_data_merge_singles$input3_e3_s0_bNA_count >= n_reads, "aa_seq"])
  
  aa_unique_output<-unique(variant_data_merge_singles[variant_data_merge_singles$input1_e1_s0_bNA_count >= n_reads &
                                                      variant_data_merge_singles$input2_e2_s0_bNA_count >= n_reads &
                                                      variant_data_merge_singles$input3_e3_s0_bNA_count >= n_reads &
                                                      (variant_data_merge_singles$output1_e1_s1_b1_count > 0 |
                                                       variant_data_merge_singles$output2_e2_s1_b1_count > 0 |
                                                       variant_data_merge_singles$output3_e3_s1_b1_count > 0), "aa_seq"])
  
  print("***********************************************************************")
  print(paste0("Variants in the input (>=", n_reads, "):"))
  print(nrow(aa_unique_input))
  print(paste0("Variants with NS (input >=", n_reads, "):"))
  print(nrow(aa_unique_output))
  print("***********************************************************************")

}

input_output_reads(10)
input_output_reads(100)
input_output_reads(200)
input_output_reads(300) # mininputcountall used in dimsum

################################################################################
# save in excel format
library("writexl")
write_xlsx(singles_stops, "hRIPK1_singles.xlsx")
write_xlsx(silent, "hRIPK1_synonymous.xlsx")

write.csv(singles_stops, "hRIPK1_singles.csv", row.names=FALSE)
write.csv(silent, "hRIPK1_synonymous.csv", row.names=FALSE)

################################################################################
