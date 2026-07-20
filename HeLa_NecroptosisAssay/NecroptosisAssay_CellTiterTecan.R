library(tidyverse)
library(ggpubr)
library(ggsignif)
require("ggrepel")

# your folder path:
setwd("")

file_name = "HeLa_NecroptosisAssay.txt"
dir.create("NecroptosisAssayPlots")
path = "NecroptosisAssayPlots"

tecan_results<-read_delim(file_name, locale=locale(encoding="latin1"), delim="\t")
tecan_results<-as.data.frame(t(tecan_results))
names(tecan_results)<-tecan_results[1,]
tecan_results<-tecan_results[2:nrow(tecan_results),]
tecan_results$sample_rep<-rownames(tecan_results)

tecan_results_pivot <- 
  tecan_results %>% 
  pivot_longer(
    cols = !sample_rep,
    names_to = "rows", 
    values_to = "Luminiscence"
  )

# Whole plate was treated with 20 ng/ml Doxy
# Add the stimulus used
tecan_results_pivot[tecan_results_pivot$rows %in% c("B", "C", "D"), "stimulus"]<-"DMSO"
tecan_results_pivot[tecan_results_pivot$rows %in% c("E", "F", "G"), "stimulus"]<-"TSI"

# Add a column with the sample
tecan_results_pivot[tecan_results_pivot$sample_rep %in% c("WT (1)", "WT (2)", "WT (3)", "WT (4)", "WT (5)", "WT (6)", "WT (7)", "WT (8)", "WT (9)"), "sample"]<-"WT"
tecan_results_pivot[tecan_results_pivot$sample_rep %in% c("I452S (1)", "I452S (2)", "I452S (3)"), "sample"]<-"I452S"
tecan_results_pivot[tecan_results_pivot$sample_rep %in% c("V458E (1)", "V458E (2)", "V458E (3)"), "sample"]<-"V458E"
tecan_results_pivot[tecan_results_pivot$sample_rep %in% c("G461Y (1)", "G461Y (2)", "G461Y (3)"), "sample"]<-"G461Y"
tecan_results_pivot[tecan_results_pivot$sample_rep %in% c("D462A (1)", "D462A (2)", "D462A (3)"), "sample"]<-"D462A"
tecan_results_pivot[tecan_results_pivot$sample_rep %in% c("L449A (1)", "L449A (2)", "L449A (3)"), "sample"]<-"L449A"
tecan_results_pivot[tecan_results_pivot$sample_rep %in% c("G457V (1)", "G457V (2)", "G457V (3)"), "sample"]<-"G457V"
tecan_results_pivot[tecan_results_pivot$sample_rep %in% c("M468A (1)", "M468A (2)", "M468A (3)"), "sample"]<-"M468A"
tecan_results_pivot[tecan_results_pivot$sample_rep %in% c("N464I (1)", "N464I (2)", "N464I (3)"), "sample"]<-"N464I"
tecan_results_pivot[tecan_results_pivot$sample_rep %in% c("L466Y (1)", "L466Y (2)", "L466Y (3)"), "sample"]<-"L466Y"
tecan_results_pivot[tecan_results_pivot$sample_rep %in% c("Y465V (1)", "Y465V (2)", "Y465V (3)"), "sample"]<-"Y465V"
tecan_results_pivot[tecan_results_pivot$sample_rep %in% c("V450H (1)", "V450H (2)", "V450H (3)"), "sample"]<-"V450H"
tecan_results_pivot[tecan_results_pivot$sample_rep %in% c("N464H (1)", "N464H (2)", "N464H (3)"), "sample"]<-"N464H"
tecan_results_pivot[tecan_results_pivot$sample_rep %in% c("L466E (1)", "L466E (2)", "L466E (3)"), "sample"]<-"L466E"

tecan_results_pivot[tecan_results_pivot$sample_rep %in% c("WT (1)", "I452S (1)", "V458E (1)", "G461Y (1)", "D462A (1)", "L449A (1)", "G457V (1)", "M468A (1)", "N464I (1)", "L466Y (1)", "Y465V (1)", "V450H (1)", "N464H (1)", "L466E (1)"), "rep"]<-1
tecan_results_pivot[tecan_results_pivot$sample_rep %in% c("WT (2)", "I452S (2)", "V458E (2)", "G461Y (2)", "D462A (2)", "L449A (2)", "G457V (2)", "M468A (2)", "N464I (2)", "L466Y (2)", "Y465V (2)", "V450H (2)", "N464H (2)", "L466E (2)" ), "rep"]<-2
tecan_results_pivot[tecan_results_pivot$sample_rep %in% c("WT (3)", "I452S (3)", "V458E (3)", "G461Y (3)", "D462A (3)", "L449A (3)", "G457V (3)", "M468A (3)", "N464I (3)", "L466Y (3)", "Y465V (3)", "V450H (3)", "N464H (3)", "L466E (3)" ), "rep"]<-3
tecan_results_pivot[tecan_results_pivot$sample_rep %in% c("WT (4)" ), "rep"]<-4
tecan_results_pivot[tecan_results_pivot$sample_rep %in% c("WT (5)" ), "rep"]<-5
tecan_results_pivot[tecan_results_pivot$sample_rep %in% c("WT (6)" ), "rep"]<-6
tecan_results_pivot[tecan_results_pivot$sample_rep %in% c("WT (7)" ), "rep"]<-7
tecan_results_pivot[tecan_results_pivot$sample_rep %in% c("WT (8)" ), "rep"]<-8
tecan_results_pivot[tecan_results_pivot$sample_rep %in% c("WT (9)" ), "rep"]<-9
#
tecan_results_pivot$sample_stimulus<-paste0(tecan_results_pivot$sample, "-", tecan_results_pivot$stimulus)
tecan_results_pivot$sample_rep_stimulus<-paste0(tecan_results_pivot$sample_rep, "-", tecan_results_pivot$stimulus)

#
tecan_results_pivot$Luminiscence<-as.numeric(tecan_results_pivot$Luminiscence)

### Normalize each sample to the well with no stimulus
#
wt_lux_1<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "WT (1)-DMSO",]$Luminiscence)
I452S_lux_1<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "I452S (1)-DMSO",]$Luminiscence)
V458E_lux_1<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "V458E (1)-DMSO",]$Luminiscence)
G461Y_lux_1<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "G461Y (1)-DMSO",]$Luminiscence)
D462A_lux_1<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "D462A (1)-DMSO",]$Luminiscence)
L449A_lux_1<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "L449A (1)-DMSO",]$Luminiscence)
G457V_lux_1<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "G457V (1)-DMSO",]$Luminiscence)
M468A_lux_1<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "M468A (1)-DMSO",]$Luminiscence)
N464I_lux_1<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "N464I (1)-DMSO",]$Luminiscence)
L466Y_lux_1<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "L466Y (1)-DMSO",]$Luminiscence)
Y465V_lux_1<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "Y465V (1)-DMSO",]$Luminiscence)
V450H_lux_1<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "V450H (1)-DMSO",]$Luminiscence)
N464H_lux_1<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "N464H (1)-DMSO",]$Luminiscence)
L466E_lux_1<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "L466E (1)-DMSO",]$Luminiscence)

wt_lux_2<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "WT (2)-DMSO",]$Luminiscence)
I452S_lux_2<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "I452S (2)-DMSO",]$Luminiscence)
V458E_lux_2<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "V458E (2)-DMSO",]$Luminiscence)
G461Y_lux_2<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "G461Y (2)-DMSO",]$Luminiscence)
D462A_lux_2<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "D462A (2)-DMSO",]$Luminiscence)
L449A_lux_2<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "L449A (2)-DMSO",]$Luminiscence)
G457V_lux_2<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "G457V (2)-DMSO",]$Luminiscence)
M468A_lux_2<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "M468A (2)-DMSO",]$Luminiscence)
N464I_lux_2<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "N464I (2)-DMSO",]$Luminiscence)
L466Y_lux_2<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "L466Y (2)-DMSO",]$Luminiscence)
Y465V_lux_2<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "Y465V (2)-DMSO",]$Luminiscence)
V450H_lux_2<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "V450H (2)-DMSO",]$Luminiscence)
N464H_lux_2<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "N464H (2)-DMSO",]$Luminiscence)
L466E_lux_2<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "L466E (2)-DMSO",]$Luminiscence)

wt_lux_3<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "WT (3)-DMSO",]$Luminiscence)
I452S_lux_3<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "I452S (3)-DMSO",]$Luminiscence)
V458E_lux_3<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "V458E (3)-DMSO",]$Luminiscence)
G461Y_lux_3<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "G461Y (3)-DMSO",]$Luminiscence)
D462A_lux_3<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "D462A (3)-DMSO",]$Luminiscence)
L449A_lux_3<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "L449A (3)-DMSO",]$Luminiscence)
G457V_lux_3<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "G457V (3)-DMSO",]$Luminiscence)
M468A_lux_3<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "M468A (3)-DMSO",]$Luminiscence)
N464I_lux_3<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "N464I (3)-DMSO",]$Luminiscence)
L466Y_lux_3<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "L466Y (3)-DMSO",]$Luminiscence)
Y465V_lux_3<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "Y465V (3)-DMSO",]$Luminiscence)
V450H_lux_3<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "V450H (3)-DMSO",]$Luminiscence)
N464H_lux_3<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "N464H (3)-DMSO",]$Luminiscence)
L466E_lux_3<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "L466E (3)-DMSO",]$Luminiscence)

wt_lux_4<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "WT (4)-DMSO",]$Luminiscence)
wt_lux_5<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "WT (5)-DMSO",]$Luminiscence)
wt_lux_6<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "WT (6)-DMSO",]$Luminiscence)
wt_lux_7<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "WT (7)-DMSO",]$Luminiscence)
wt_lux_8<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "WT (8)-DMSO",]$Luminiscence)
wt_lux_9<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "WT (9)-DMSO",]$Luminiscence)

#
tecan_results_pivot[tecan_results_pivot$sample_rep == "WT (1)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "WT (1)", "Luminiscence"]/wt_lux_1*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "I452S (1)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "I452S (1)", "Luminiscence"]/I452S_lux_1*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "V458E (1)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "V458E (1)", "Luminiscence"]/V458E_lux_1*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "G461Y (1)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "G461Y (1)", "Luminiscence"]/G461Y_lux_1*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "D462A (1)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "D462A (1)", "Luminiscence"]/D462A_lux_1*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "L449A (1)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "L449A (1)", "Luminiscence"]/L449A_lux_1*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "G457V (1)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "G457V (1)", "Luminiscence"]/G457V_lux_1*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "M468A (1)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "M468A (1)", "Luminiscence"]/M468A_lux_1*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "N464I (1)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "N464I (1)", "Luminiscence"]/N464I_lux_1*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "L466Y (1)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "L466Y (1)", "Luminiscence"]/L466Y_lux_1*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "Y465V (1)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "Y465V (1)", "Luminiscence"]/Y465V_lux_1*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "V450H (1)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "V450H (1)", "Luminiscence"]/V450H_lux_1*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "N464H (1)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "N464H (1)", "Luminiscence"]/N464H_lux_1*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "L466E (1)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "L466E (1)", "Luminiscence"]/L466E_lux_1*100


tecan_results_pivot[tecan_results_pivot$sample_rep == "WT (2)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "WT (2)", "Luminiscence"]/wt_lux_2*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "I452S (2)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "I452S (2)", "Luminiscence"]/I452S_lux_2*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "V458E (2)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "V458E (2)", "Luminiscence"]/V458E_lux_2*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "G461Y (2)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "G461Y (2)", "Luminiscence"]/G461Y_lux_2*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "D462A (2)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "D462A (2)", "Luminiscence"]/D462A_lux_2*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "L449A (2)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "L449A (2)", "Luminiscence"]/L449A_lux_2*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "G457V (2)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "G457V (2)", "Luminiscence"]/G457V_lux_2*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "M468A (2)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "M468A (2)", "Luminiscence"]/M468A_lux_2*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "N464I (2)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "N464I (2)", "Luminiscence"]/N464I_lux_2*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "L466Y (2)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "L466Y (2)", "Luminiscence"]/L466Y_lux_2*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "Y465V (2)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "Y465V (2)", "Luminiscence"]/Y465V_lux_2*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "V450H (2)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "V450H (2)", "Luminiscence"]/V450H_lux_2*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "N464H (2)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "N464H (2)", "Luminiscence"]/N464H_lux_2*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "L466E (2)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "L466E (2)", "Luminiscence"]/L466E_lux_2*100

tecan_results_pivot[tecan_results_pivot$sample_rep == "WT (3)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "WT (3)", "Luminiscence"]/wt_lux_3*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "I452S (3)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "I452S (3)", "Luminiscence"]/I452S_lux_3*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "V458E (3)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "V458E (3)", "Luminiscence"]/V458E_lux_3*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "G461Y (3)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "G461Y (3)", "Luminiscence"]/G461Y_lux_3*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "D462A (3)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "D462A (3)", "Luminiscence"]/D462A_lux_3*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "L449A (3)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "L449A (3)", "Luminiscence"]/L449A_lux_3*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "G457V (3)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "G457V (3)", "Luminiscence"]/G457V_lux_3*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "M468A (3)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "M468A (3)", "Luminiscence"]/M468A_lux_3*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "N464I (3)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "N464I (3)", "Luminiscence"]/N464I_lux_3*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "L466Y (3)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "L466Y (3)", "Luminiscence"]/L466Y_lux_3*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "Y465V (3)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "Y465V (3)", "Luminiscence"]/Y465V_lux_3*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "V450H (3)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "V450H (3)", "Luminiscence"]/V450H_lux_3*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "N464H (3)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "N464H (3)", "Luminiscence"]/N464H_lux_3*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "L466E (3)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "L466E (3)", "Luminiscence"]/L466E_lux_3*100

tecan_results_pivot[tecan_results_pivot$sample_rep == "WT (4)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "WT (4)", "Luminiscence"]/wt_lux_4*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "WT (5)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "WT (5)", "Luminiscence"]/wt_lux_5*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "WT (6)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "WT (6)", "Luminiscence"]/wt_lux_6*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "WT (7)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "WT (7)", "Luminiscence"]/wt_lux_7*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "WT (8)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "WT (8)", "Luminiscence"]/wt_lux_8*100
tecan_results_pivot[tecan_results_pivot$sample_rep == "WT (9)", "LuxNorm"]<-tecan_results_pivot[tecan_results_pivot$sample_rep == "WT (9)", "Luminiscence"]/wt_lux_9*100



### Plot for each replicate the result
#
p_box_1 <- ggplot(tecan_results_pivot, aes(x = stimulus, y = Luminiscence)) +
  geom_boxplot(position = position_dodge(width = 0.9), alpha=0.5) +
  geom_jitter(alpha = 0.8) +
  facet_grid(rep~factor(sample, levels=c("WT", "I452S", "V458E", "G461Y", 
                                         "D462A", "L449A", "G457V", "M468A", 
                                         "N464I", "L466Y", "Y465V",
                                         "V450H", "N464H", "L466E")))+
  labs(y="Counts/s (Viable Cells)", x="")+
  theme_bw()
p_box_1

ggsave(path=path, p_box_1, file="HeLa_NecroptosisAssay_RHIM_1.jpg", width=10, height=5)

#
p_box_2 <- ggplot(tecan_results_pivot, aes(x = stimulus, y = LuxNorm)) +
  geom_boxplot(position = position_dodge(width = 0.9), alpha=0.5) +
  geom_jitter() +
  facet_grid(rep~factor(sample, levels=c("WT", "I452S", "V458E", "G461Y", 
                                         "D462A", "L449A", "G457V", "M468A", 
                                         "N464I", "L466Y", "Y465V",
                                         "V450H", "N464H", "L466E")))+
  labs(y="Counts/s (Viable Cells)", x="")+
  theme_bw()
p_box_2

ggsave(path=path, p_box_2, file="HeLa_NecroptosisAssay_RHIM_1_LuxNorm.jpg", width=10, height=5)

### Get the mean value from each experiment and do only 1 plot
#
wt_lux_1_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "WT (1)-DMSO",]$LuxNorm)
wt_lux_1_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "WT (1)-TSI",]$LuxNorm)
I452S_lux_1_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "I452S (1)-DMSO",]$LuxNorm)
I452S_lux_1_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "I452S (1)-TSI",]$LuxNorm)
V458E_lux_1_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "V458E (1)-DMSO",]$LuxNorm)
V458E_lux_1_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "V458E (1)-TSI",]$LuxNorm)
G461Y_lux_1_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "G461Y (1)-DMSO",]$LuxNorm)
G461Y_lux_1_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "G461Y (1)-TSI",]$LuxNorm)
D462A_lux_1_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "D462A (1)-DMSO",]$LuxNorm)
D462A_lux_1_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "D462A (1)-TSI",]$LuxNorm)
L449A_lux_1_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "L449A (1)-DMSO",]$LuxNorm)
L449A_lux_1_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "L449A (1)-TSI",]$LuxNorm)
G457V_lux_1_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "G457V (1)-DMSO",]$LuxNorm)
G457V_lux_1_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "G457V (1)-TSI",]$LuxNorm)
M468A_lux_1_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "M468A (1)-DMSO",]$LuxNorm)
M468A_lux_1_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "M468A (1)-TSI",]$LuxNorm)
N464I_lux_1_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "N464I (1)-DMSO",]$LuxNorm)
N464I_lux_1_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "N464I (1)-TSI",]$LuxNorm)
L466Y_lux_1_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "L466Y (1)-DMSO",]$LuxNorm)
L466Y_lux_1_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "L466Y (1)-TSI",]$LuxNorm)
Y465V_lux_1_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "Y465V (1)-DMSO",]$LuxNorm)
Y465V_lux_1_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "Y465V (1)-TSI",]$LuxNorm)
V450H_lux_1_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "V450H (1)-DMSO",]$LuxNorm)
V450H_lux_1_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "V450H (1)-TSI",]$LuxNorm)
N464H_lux_1_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "N464H (1)-DMSO",]$LuxNorm)
N464H_lux_1_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "N464H (1)-TSI",]$LuxNorm)
L466E_lux_1_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "L466E (1)-DMSO",]$LuxNorm)
L466E_lux_1_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "L466E (1)-TSI",]$LuxNorm)

#
wt_lux_2_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "WT (2)-DMSO",]$LuxNorm)
wt_lux_2_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "WT (2)-TSI",]$LuxNorm)
I452S_lux_2_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "I452S (2)-DMSO",]$LuxNorm)
I452S_lux_2_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "I452S (2)-TSI",]$LuxNorm)
V458E_lux_2_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "V458E (2)-DMSO",]$LuxNorm)
V458E_lux_2_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "V458E (2)-TSI",]$LuxNorm)
G461Y_lux_2_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "G461Y (2)-DMSO",]$LuxNorm)
G461Y_lux_2_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "G461Y (2)-TSI",]$LuxNorm)
D462A_lux_2_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "D462A (2)-DMSO",]$LuxNorm)
D462A_lux_2_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "D462A (2)-TSI",]$LuxNorm)
L449A_lux_2_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "L449A (2)-DMSO",]$LuxNorm)
L449A_lux_2_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "L449A (2)-TSI",]$LuxNorm)
G457V_lux_2_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "G457V (2)-DMSO",]$LuxNorm)
G457V_lux_2_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "G457V (2)-TSI",]$LuxNorm)
M468A_lux_2_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "M468A (2)-DMSO",]$LuxNorm)
M468A_lux_2_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "M468A (2)-TSI",]$LuxNorm)
N464I_lux_2_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "N464I (2)-DMSO",]$LuxNorm)
N464I_lux_2_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "N464I (2)-TSI",]$LuxNorm)
L466Y_lux_2_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "L466Y (2)-DMSO",]$LuxNorm)
L466Y_lux_2_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "L466Y (2)-TSI",]$LuxNorm)
Y465V_lux_2_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "Y465V (2)-DMSO",]$LuxNorm)
Y465V_lux_2_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "Y465V (2)-TSI",]$LuxNorm)
V450H_lux_2_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "V450H (1)-DMSO",]$LuxNorm)
V450H_lux_2_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "V450H (2)-TSI",]$LuxNorm)
N464H_lux_2_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "N464H (2)-DMSO",]$LuxNorm)
N464H_lux_2_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "N464H (2)-TSI",]$LuxNorm)
L466E_lux_2_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "L466E (2)-DMSO",]$LuxNorm)
L466E_lux_2_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "L466E (2)-TSI",]$LuxNorm)
#
wt_lux_3_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "WT (3)-DMSO",]$LuxNorm)
wt_lux_3_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "WT (3)-TSI",]$LuxNorm)
I452S_lux_3_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "I452S (3)-DMSO",]$LuxNorm)
I452S_lux_3_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "I452S (3)-TSI",]$LuxNorm)
V458E_lux_3_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "V458E (3)-DMSO",]$LuxNorm)
V458E_lux_3_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "V458E (3)-TSI",]$LuxNorm)
G461Y_lux_3_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "G461Y (3)-DMSO",]$LuxNorm)
G461Y_lux_3_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "G461Y (3)-TSI",]$LuxNorm)
D462A_lux_3_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "D462A (3)-DMSO",]$LuxNorm)
D462A_lux_3_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "D462A (3)-TSI",]$LuxNorm)
L449A_lux_3_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "L449A (3)-DMSO",]$LuxNorm)
L449A_lux_3_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "L449A (3)-TSI",]$LuxNorm)
G457V_lux_3_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "G457V (3)-DMSO",]$LuxNorm)
G457V_lux_3_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "G457V (3)-TSI",]$LuxNorm)
M468A_lux_3_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "M468A (3)-DMSO",]$LuxNorm)
M468A_lux_3_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "M468A (3)-TSI",]$LuxNorm)
N464I_lux_3_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "N464I (3)-DMSO",]$LuxNorm)
N464I_lux_3_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "N464I (3)-TSI",]$LuxNorm)
L466Y_lux_3_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "L466Y (3)-DMSO",]$LuxNorm)
L466Y_lux_3_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "L466Y (3)-TSI",]$LuxNorm)
Y465V_lux_3_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "Y465V (3)-DMSO",]$LuxNorm)
Y465V_lux_3_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "Y465V (3)-TSI",]$LuxNorm)
V450H_lux_3_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "V450H (3)-DMSO",]$LuxNorm)
V450H_lux_3_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "V450H (3)-TSI",]$LuxNorm)
N464H_lux_3_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "N464H (3)-DMSO",]$LuxNorm)
N464H_lux_3_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "N464H (3)-TSI",]$LuxNorm)
L466E_lux_3_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "L466E (3)-DMSO",]$LuxNorm)
L466E_lux_3_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "L466E (3)-TSI",]$LuxNorm)

#
wt_lux_4_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "WT (4)-DMSO",]$LuxNorm)
wt_lux_4_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "WT (4)-TSI",]$LuxNorm)

wt_lux_5_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "WT (5)-DMSO",]$LuxNorm)
wt_lux_5_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "WT (5)-TSI",]$LuxNorm)

wt_lux_6_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "WT (6)-DMSO",]$LuxNorm)
wt_lux_6_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "WT (6)-TSI",]$LuxNorm)

wt_lux_7_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "WT (7)-DMSO",]$LuxNorm)
wt_lux_7_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "WT (7)-TSI",]$LuxNorm)

wt_lux_8_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "WT (8)-DMSO",]$LuxNorm)
wt_lux_8_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "WT (8)-TSI",]$LuxNorm)

wt_lux_9_dmso<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "WT (9)-DMSO",]$LuxNorm)
wt_lux_9_tsi<-mean(tecan_results_pivot[tecan_results_pivot$sample_rep_stimulus == "WT (9)-TSI",]$LuxNorm)



df_triplicates<-data.frame("sample"=rep(c("WT", "I452S", "V458E", "G461Y", "D462A", "L449A", "G457V", "M468A", "N464I", "L466Y", "Y465V", "V450H", "N464H","L466E"), 6),
                           "rep"=rep(c(rep(c("rep1"), 14), rep(c("rep2"), 14), rep(c("rep3"), 14)), 2),
                           "stimulus"=c(rep(c("DMSO"), 42), rep(c("TSI"), 42)),
                           "Cell_survival"=c(wt_lux_1_dmso, I452S_lux_1_dmso, V458E_lux_1_dmso, G461Y_lux_1_dmso, D462A_lux_1_dmso, L449A_lux_1_dmso, G457V_lux_1_dmso, M468A_lux_1_dmso, N464I_lux_1_dmso, L466Y_lux_1_dmso, Y465V_lux_1_dmso, V450H_lux_1_dmso, N464H_lux_1_dmso, L466E_lux_1_dmso,
                                             wt_lux_2_dmso, I452S_lux_2_dmso, V458E_lux_2_dmso, G461Y_lux_2_dmso, D462A_lux_2_dmso, L449A_lux_2_dmso, G457V_lux_2_dmso, M468A_lux_2_dmso, N464I_lux_2_dmso, L466Y_lux_2_dmso, Y465V_lux_2_dmso, V450H_lux_2_dmso, N464H_lux_2_dmso, L466E_lux_2_dmso,
                                             wt_lux_3_dmso, I452S_lux_3_dmso, V458E_lux_3_dmso, G461Y_lux_3_dmso, D462A_lux_3_dmso, L449A_lux_3_dmso, G457V_lux_3_dmso, M468A_lux_3_dmso, N464I_lux_3_dmso, L466Y_lux_3_dmso, Y465V_lux_3_dmso, V450H_lux_3_dmso, N464H_lux_3_dmso, L466E_lux_3_dmso,
                                             wt_lux_1_tsi, I452S_lux_1_tsi, V458E_lux_1_tsi, G461Y_lux_1_tsi, D462A_lux_1_tsi, L449A_lux_1_tsi, G457V_lux_1_tsi, M468A_lux_1_tsi, N464I_lux_1_tsi, L466Y_lux_1_tsi, Y465V_lux_1_tsi, V450H_lux_1_tsi, N464H_lux_1_tsi, L466E_lux_1_tsi,
                                             wt_lux_2_tsi, I452S_lux_2_tsi, V458E_lux_2_tsi, G461Y_lux_2_tsi, D462A_lux_2_tsi, L449A_lux_2_tsi, G457V_lux_2_tsi, M468A_lux_2_tsi, N464I_lux_2_tsi, L466Y_lux_2_tsi, Y465V_lux_2_tsi, V450H_lux_2_tsi, N464H_lux_2_tsi, L466E_lux_2_tsi,
                                             wt_lux_3_tsi, I452S_lux_3_tsi, V458E_lux_3_tsi, G461Y_lux_3_tsi, D462A_lux_3_tsi, L449A_lux_3_tsi, G457V_lux_3_tsi, M468A_lux_3_tsi, N464I_lux_3_tsi, L466Y_lux_3_tsi, Y465V_lux_3_tsi, V450H_lux_3_tsi, N464H_lux_3_tsi, L466E_lux_3_tsi
                           )) 
df_wt<-data.frame("sample"=rep(c("WT"), 12),
                  "rep"=c(rep(c("rep4"), 2), rep(c("rep5"), 2), rep(c("rep6"), 2), rep(c("rep7"), 2), rep(c("rep8"), 2), rep(c("rep9"), 2)),
                  "stimulus"=rep(c("DMSO", "TSI"), 6),
                  "Cell_survival"=c(wt_lux_4_dmso, wt_lux_4_tsi, wt_lux_5_dmso, wt_lux_5_tsi, wt_lux_6_dmso, wt_lux_6_tsi, wt_lux_7_dmso, wt_lux_7_tsi, wt_lux_8_dmso, wt_lux_8_tsi, wt_lux_9_dmso, wt_lux_9_tsi))

df_triplicates<-rbind(df_triplicates, df_wt)

# function to summarize data
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df_summary<-data_summary(df_triplicates, "Cell_survival", c("sample", "stimulus"))

# Create the bar plot with error bars using the summarized data
p_bar_1 <- ggplot(df_summary, aes(x=factor(sample, levels=c("V458E", "I452S", "WT", "L449A", "M468A", "D462A", "G457V", "G461Y", "V450H", "N464H", "L466E")), 
                                  y=Cell_survival, fill=stimulus)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5, color="black") +
  geom_errorbar(aes(ymin = Cell_survival, ymax = Cell_survival + sd), width = 0.2, position = position_dodge(0.5)) +
  scale_fill_manual(values=c("grey20", "grey80"))+
  scale_y_continuous(limits = c(0, 120), breaks=c(0, 30, 60, 90, 120))+
  labs(x="")+
  theme_classic()
p_bar_1

ggsave(path=path, p_bar_1, file="HeLa_NecroptosisAssay_RHIM_1_CellTiter_bars.jpg", width=10, height=3)
ggsave(path=path, p_bar_1, file="HeLa_NecroptosisAssay_RHIM_1_CellTiter_bars.pdf", width=10, height=3)


#########################################################################################################
# Necreoptosis vs Nucleation Score
load("nscore_df_hRIPK3.RData")

singles_stops$category_1<-"WT-like"
singles_stops[singles_stops$p.adjust<0.01 & singles_stops$nscore_c >0, "category_1"]<-"NS_inc"
singles_stops[singles_stops$p.adjust<0.01 & singles_stops$nscore_c <0, "category_1"]<-"NS_dec"


df_summary$ID<-paste0(str_sub(df_summary$sample, 1, 1), "-", str_sub(df_summary$sample, 2, 4), "-", str_sub(df_summary$sample, 5, 5))

df_lum_necorptosis<-df_summary[df_summary$stimulus == "TSI", c("ID", "Cell_survival", "sd")]

singles_ns_death<-inner_join(singles_stops[c("ID", "nscore_c", "sigma", "category_1")], df_lum_necorptosis, by="ID")

singles_ns_death<-rbind(data.frame("ID"=c("WT"), "nscore_c"=c(0.0006341024), "sigma"=c(0.07340256), "category_1"=c("WT-like"),
                                   "Cell_survival"=c(45.16202), "sd"=c(1.045397e+01)), singles_ns_death)

singles_ns_death<-singles_ns_death[!duplicated(singles_ns_death), ]

scatter_ns<-ggplot(singles_ns_death, aes(x=nscore_c, y=Cell_survival, label=ID, color=category_1))+
  geom_errorbar(aes(ymin=Cell_survival-sd, ymax=Cell_survival+sd), 
                color="grey80", width=.1)+
  geom_errorbar(aes(xmin=nscore_c-sigma, xmax=nscore_c+sigma),
                color="grey80", width=1)+
  geom_point(size=4)+
  scale_color_manual(values=c("#CD3333", "#00008B", "grey"))+
  geom_text_repel(size=4, max.overlaps = Inf)+
  stat_cor(label.y = 30, label.x = -5.7, color="black")+ 
  labs(x="Nucleation Score", y="% Cell survival")+
  theme_classic()+
  theme(legend.position = "none")

scatter_ns

ggsave(path=path, scatter_ns, file="HeLa_NecroptosisAssay_RHIM_1_nscore.jpg", width=4, height=4)
ggsave(path=path, scatter_ns, file="HeLa_NecroptosisAssay_RHIM_1_nscore.pdf", width=4, height=4)

###
min<-min(singles_stops$nscore_c)
max<-max(singles_stops$nscore_c)
cols <- c(colorRampPalette(c( "brown3", "grey95"))((-min/(-min+max)*100)-0.5), colorRampPalette("grey95")(1),
          colorRampPalette(c("grey95",  "darkblue"), bias=1)((max/(-min+max)*100)-0.5))


# Create the bar plot with error bars using the summarized data and the color scale of the nucleation score
p_bar_2 <- ggplot(singles_ns_death, aes(x=factor(ID, levels=c("V-458-E", "V-450-H" ,"I-452-S", "WT", "L-449-A", "M-468-A", 
                                                              "Y-465-V", "L-466-Y", "D-462-A", "G-457-V", "N-464-I", 
                                                              "N-464-H", "G-461-Y", "L-466-E")), 
                                        y=Cell_survival, fill=nscore_c)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5, color="black") +
  geom_errorbar(aes(ymin = Cell_survival, ymax = Cell_survival + sd), width = 0.2, position = position_dodge(0.5)) +
  #scale_fill_manual(values=c("grey20", "grey80"))+
  scale_fill_gradientn(colours=cols, limits=c(min,max))+
  labs(x="")+
  theme_classic()
p_bar_2

ggsave(path=path, p_bar_2, file="HeLa_NecroptosisAssay_RHIM_1_nscore_bar.jpg", width=8, height=3)
ggsave(path=path, p_bar_2, file="HeLa_NecroptosisAssay_RHIM_1_nscore_bar.pdf", width=8, height=3)


### Test the difference between the mutants and the wt
df_triplicates$ID_experiment<-factor(paste0(df_triplicates$sample, "-", df_triplicates$stimulus), 
                                     levels=c("WT-TSI", "I452S-TSI", "V458E-TSI", 
                                              "G461Y-TSI", "D462A-TSI", "L449A-TSI", "G457V-TSI", 
                                              "M468A-TSI", "N464I-TSI", "L466Y-TSI", 
                                              "Y465V-TSI", "V450H-TSI", "N464H-TSI", "L466E-TSI"))
# ANOVA
ANOVA=aov(Cell_survival ~ ID_experiment, data=df_triplicates[df_triplicates$stimulus == "TSI",])
summary(ANOVA)
# Perform Dunnett's post-hoc test
library(multcomp)
summary(glht(ANOVA, linfct = mcp(ID_experiment = "Dunnett")))
# Perform Tukey's post-hoc test
TukeyHSD(ANOVA)

###################################################################################
# Necreoptosis vs ddG
#PDB 7DAC
load("ddG_7DAC.RData")

ddG_foldx_7DAC$sig_1_7DAC<-F
ddG_foldx_7DAC[ddG_foldx_7DAC$p.adjust_7DAC<0.01, "sig_1_7DAC"]<-T

ddG_foldx_7DAC$category_1_7DAC<-"WT-like"
ddG_foldx_7DAC[ddG_foldx_7DAC$p.adjust_7DAC<0.01 & ddG_foldx_7DAC$ddG_foldx_7DAC >2, "category_1_7DAC"]<-"ddG_inc"
ddG_foldx_7DAC[ddG_foldx_7DAC$p.adjust_7DAC<0.01 & ddG_foldx_7DAC$ddG_foldx_7DAC <0, "category_1_7DAC"]<-"ddG_dec"

ddG_foldx_7DAC<-ddG_foldx_7DAC[ddG_foldx_7DAC$low_ddG_foldx_SD_7DAC == T | ddG_foldx_7DAC$sig_1_7DAC == T, ]

singles_ddG_7DAC_death<-inner_join(ddG_foldx_7DAC[c("ID", "ddG_foldx_7DAC", "ddG_foldx_SD_7DAC", "category_1_7DAC")], df_lum_necorptosis, by="ID")

singles_ddG_7DAC_death<-rbind(data.frame("ID"=c("WT"), "ddG_foldx_7DAC"=c(0), "ddG_foldx_SD_7DAC"=c(0), "category_1_7DAC"=c("WT-like"),
                                         "Cell_survival"=c(45.16202), "sd"=c(1.045397e+01)), singles_ddG_7DAC_death)

singles_ddG_7DAC_death<-singles_ddG_7DAC_death[!duplicated(singles_ddG_7DAC_death), ]

scatter_ddg_7DAC<-ggplot(singles_ddG_7DAC_death, aes(x=ddG_foldx_7DAC, y=Cell_survival, label=ID, color=category_1_7DAC))+
  geom_errorbar(aes(ymin=Cell_survival-sd, ymax=Cell_survival+sd), 
                color="grey80", width=.2)+
  geom_errorbar(aes(xmin=ddG_foldx_7DAC-ddG_foldx_SD_7DAC, xmax=ddG_foldx_7DAC+ddG_foldx_SD_7DAC),
                color="grey80", width=1)+
  geom_point(size=4)+
  geom_text_repel(size=4)+
  scale_color_manual(values=c("#4775d1", "#b3003b", "grey"))+
  stat_cor(label.y = 110, color="black")+ 
  labs(x="ddG FoldX (PDB:7DAC)", y="% Cell survival")+
  theme_classic()+
  theme(legend.position = "none")

scatter_ddg_7DAC

ggsave(path=path, scatter_ddg_7DAC, file="ddG_cellsurvival_7DAC.jpg", width=4, height=4)
ggsave(path=path, scatter_ddg_7DAC, file="ddG_cellsurvival_7DAC.pdf", width=4, height=4)


#PDB 7DA4
load("ddG_7DA4.RData")

ddG_foldx_7DA4$sig_1_7DA4<-F
ddG_foldx_7DA4[ddG_foldx_7DA4$p.adjust_7DA4<0.01, "sig_1_7DA4"]<-T

ddG_foldx_7DA4$category_1_7DA4<-"WT-like"
ddG_foldx_7DA4[ddG_foldx_7DA4$p.adjust_7DA4<0.01 & ddG_foldx_7DA4$ddG_foldx_7DA4 >2, "category_1_7DA4"]<-"ddG_inc"
ddG_foldx_7DA4[ddG_foldx_7DA4$p.adjust_7DA4<0.01 & ddG_foldx_7DA4$ddG_foldx_7DA4 <0, "category_1_7DA4"]<-"ddG_dec"

ddG_foldx_7DA4<-ddG_foldx_7DA4[ddG_foldx_7DA4$low_ddG_foldx_SD_7DA4 == T | ddG_foldx_7DA4$sig_1_7DA4 == T, ]

singles_ddG_7DA4_death<-inner_join(ddG_foldx_7DA4[c("ID", "ddG_foldx_7DA4", "ddG_foldx_SD_7DA4", "category_1_7DA4")], df_lum_necorptosis, by="ID")

singles_ddG_7DA4_death<-rbind(data.frame("ID"=c("WT"), "ddG_foldx_7DA4"=c(0), "ddG_foldx_SD_7DA4"=c(0), "category_1_7DA4"=c("WT-like"),
                                         "Cell_survival"=c(45.16202), "sd"=c(1.045397e+01)), singles_ddG_7DA4_death)

singles_ddG_7DA4_death<-singles_ddG_7DA4_death[!duplicated(singles_ddG_7DA4_death), ]

scatter_ddg_7DA4<-ggplot(singles_ddG_7DA4_death, aes(x=ddG_foldx_7DA4, y=Cell_survival, label=ID, color=category_1_7DA4))+
  geom_errorbar(aes(ymin=Cell_survival-sd, ymax=Cell_survival+sd), 
                color="grey80", width=.1)+
  geom_errorbar(aes(xmin=ddG_foldx_7DA4-ddG_foldx_SD_7DA4, xmax=ddG_foldx_7DA4+ddG_foldx_SD_7DA4),
                color="grey80", width=1)+
  geom_point(size=4)+
  geom_text_repel(size=4)+
  scale_color_manual(values=c("#b3003b", "grey"))+
  stat_cor(label.y = 110, color="black")+ 
  labs(x="ddG FoldX (PDB:7DA4)", y="% Cell survival")+
  theme_classic()+
  theme(legend.position = "none")

scatter_ddg_7DA4

ggsave(path=path, scatter_ddg_7DA4, file="ddG_cellsurvival_7DA4.jpg", width=4, height=4)
ggsave(path=path, scatter_ddg_7DA4, file="ddG_cellsurvival_7DA4.pdf", width=4, height=4)

#PDB 8Z94
load("ddG_8Z94.RData")

ddG_foldx_8Z94$sig_1_8Z94<-F
ddG_foldx_8Z94[ddG_foldx_8Z94$p.adjust_8Z94<0.01, "sig_1_7DA4"]<-T

ddG_foldx_8Z94$category_1_8Z94<-"WT-like"
ddG_foldx_8Z94[ddG_foldx_8Z94$p.adjust_8Z94<0.01 & ddG_foldx_8Z94$ddG_foldx_8Z94 >2, "category_1_8Z94"]<-"ddG_inc"
ddG_foldx_8Z94[ddG_foldx_8Z94$p.adjust_8Z94<0.01 & ddG_foldx_8Z94$ddG_foldx_8Z94 <0, "category_1_8Z94"]<-"ddG_dec"

ddG_foldx_8Z94<-ddG_foldx_8Z94[ddG_foldx_8Z94$low_ddG_foldx_SD_8Z94 == T | ddG_foldx_8Z94$sig_1_8Z94 == T, ]

singles_ddG_8Z94_death<-inner_join(ddG_foldx_8Z94[c("ID", "ddG_foldx_8Z94", "ddG_foldx_SD_8Z94", "category_1_8Z94")], df_lum_necorptosis, by="ID")

singles_ddG_8Z94_death<-rbind(data.frame("ID"=c("WT"), "ddG_foldx_8Z94"=c(0), "ddG_foldx_SD_8Z94"=c(0), "category_1_8Z94"=c("WT-like"),
                                         "Cell_survival"=c(45.16202), "sd"=c(1.045397e+01)), singles_ddG_8Z94_death)

singles_ddG_8Z94_death<-singles_ddG_8Z94_death[!duplicated(singles_ddG_8Z94_death), ]

scatter_ddg_8Z94<-ggplot(singles_ddG_8Z94_death, aes(x=ddG_foldx_8Z94, y=Cell_survival, label=ID, color=category_1_8Z94))+
  geom_errorbar(aes(ymin=Cell_survival-sd, ymax=Cell_survival+sd), 
                color="grey80", width=.1)+
  geom_errorbar(aes(xmin=ddG_foldx_8Z94-ddG_foldx_SD_8Z94, xmax=ddG_foldx_8Z94+ddG_foldx_SD_8Z94),
                color="grey80", width=1)+
  geom_point(size=4)+
  geom_text_repel(size=4)+
  scale_color_manual(values=c("#b3003b", "grey"))+
  stat_cor(label.y = 110, color="black")+ 
  labs(x="ddG FoldX (PDB:8Z94)", y="% Cell survival")+
  theme_classic()+
  theme(legend.position = "none")

scatter_ddg_8Z94

ggsave(path=path, scatter_ddg_8Z94, file="ddG_cellsurvival_8Z94.jpg", width=4, height=4)
ggsave(path=path, scatter_ddg_8Z94, file="ddG_cellsurvival_8Z94.pdf", width=4, height=4)

###