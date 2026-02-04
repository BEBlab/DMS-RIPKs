library(tidyverse)
library(ggpubr)
library(ggsignif)
require("ggrepel")
library(ggbreak) 
library(reshape2)

# your folder path
setwd("")

dir.create("02_ReportedData_mRIP3")
path="02_ReportedData_mRIP3"

### mRIP3
# Nucleation Score
load("nscore_df_mRIP3.RData")

singles_stops$sig_1<-F
singles_stops[singles_stops$p.adjust<0.01, "sig_1"]<-T

singles_stops$category_1<-"WT-like"
singles_stops[singles_stops$p.adjust<0.01 & singles_stops$nscore_c >0, "category_1"]<-"NS_inc"
singles_stops[singles_stops$p.adjust<0.01 & singles_stops$nscore_c <0, "category_1"]<-"NS_dec"

# Let's be super stringent with this comaprisons
singles_stops["low_sigma"]<-F
singles_stops[singles_stops$sigma_norm_iqr<0.2, "low_sigma"]<-T
singles_stops<-singles_stops[singles_stops$low_sigma == T | singles_stops$sig_1 == T, ]

singles_stops$Residue<-paste0(singles_stops$WT_AA, singles_stops$Pos)

### ddG
#PDB 6JPD
load("ddG_6JPD.RData")

ddG_foldx_6JPD$sig_1_6JPD<-F
ddG_foldx_6JPD[ddG_foldx_6JPD$p.adjust_6JPD<0.01, "sig_1_6JPD"]<-T

ddG_foldx_6JPD$category_1_6JPD<-"WT-like"
ddG_foldx_6JPD[ddG_foldx_6JPD$p.adjust_6JPD<0.01 & ddG_foldx_6JPD$ddG_foldx_6JPD >2, "category_1_6JPD"]<-"ddG_inc"
ddG_foldx_6JPD[ddG_foldx_6JPD$p.adjust_6JPD<0.01 & ddG_foldx_6JPD$ddG_foldx_6JPD <0, "category_1_6JPD"]<-"ddG_dec"

ddG_foldx_6JPD<-ddG_foldx_6JPD[ddG_foldx_6JPD$low_ddG_foldx_SD_6JPD == T | ddG_foldx_6JPD$sig_1_6JPD == T, ]

###
NS_ddG_mRIP3<-full_join(ddG_foldx_6JPD[c("ID", "ddG_foldx_6JPD", "ddG_foldx_SD_6JPD", "category_1_6JPD")],
                        singles_stops[c("ID", "Pos", "Residue", "nscore_c", "sigma", "category_10", "category_1")], by="ID")


################################################################################
# 1: NS Correlation vs mutants from the literature

mRIP3_mutants_HeLa<-data.frame("ID"=c("WT", "F-442-D", "Q-449-D", "L-456-D", "F-442-A", "Q-449-A", "L-456-A"),
                               "cell_survival"=c(34.35, 100.33, 100.46, 99.3, 59.44, 64.05, 53.09),
                               "SD"=c(0.96, 3.27, 0.79, 3.43, 5.39, 6.38, 3.61))

singles_ns_death<-inner_join(singles_stops[c("ID", "nscore_c", "sigma", "category_1")], mRIP3_mutants_HeLa, by="ID")

singles_ns_death<-rbind(data.frame("ID"=c("WT"), "nscore_c"=c(0.056252746), "sigma"=c(0.01378721), "category_1"=c("WT-like"),
                                   "cell_survival"=c(34.35), "SD"=c(0.96)), singles_ns_death)

scatter_1<-ggplot(singles_ns_death, aes(x=nscore_c, y=cell_survival, label=ID, color=category_1))+
  geom_errorbar(aes(ymin=cell_survival-SD, ymax=cell_survival+SD), 
                color="grey80", width=.05)+
  geom_errorbar(aes(xmin=nscore_c-sigma, xmax=nscore_c+sigma),
                color="grey80", width=1)+
  geom_point(size=4)+
  scale_color_manual(values=c("#CD3333", "#00008B", "grey"))+
  geom_text_repel(size=4)+
  stat_cor(label.y = 40, label.x = -4, color="black")+ 
  labs(x="Nucleation Score", y="% Cell Survival (Wu et al 2021)")+
  theme_classic()+
  theme(legend.position = "none")

scatter_1

ggsave(scatter_1, path=path, file="nscore_cellsurvival_Wuetal2021.jpg", width=4, height=4)
ggsave(scatter_1, path=path, file="nscore_cellsurvival_Wuetal2021.pdf", width=4, height=4)

#################################################################################
# 3: ddG Correlation vs mutants from the literature

singles_ddG_death<-inner_join(NS_ddG_mRIP3[c("ID", "ddG_foldx_6JPD", "ddG_foldx_SD_6JPD", "category_1_6JPD")], mRIP3_mutants_HeLa, by="ID")

singles_ddG_death<-rbind(data.frame("ID"=c("WT"), "ddG_foldx_6JPD"=c(0), "ddG_foldx_SD_6JPD"=c(0), "category_1_6JPD"=c("WT-like"), 
                                    "cell_survival"=c(34.35), "SD"=c(0.96)), singles_ddG_death)

scatter_2<-ggplot(singles_ddG_death, aes(x=ddG_foldx_6JPD, y=cell_survival, label=ID, color=category_1_6JPD))+
  geom_errorbar(aes(ymin=cell_survival-SD, ymax=cell_survival+SD), 
                color="grey80", width=.05)+
  geom_errorbar(aes(xmin=ddG_foldx_6JPD-ddG_foldx_SD_6JPD, xmax=ddG_foldx_6JPD+ddG_foldx_SD_6JPD),
                color="grey80", width=1)+
  geom_point(size=4)+
  geom_text_repel(size=4)+
  scale_color_manual(values=c("#b3003b", "grey"))+
  stat_cor(label.y = 30, label.x = 0, color="black")+ 
  labs(x="ddG 6JPD", y="% Cell Surivival (Wu et al 2021)")+
  theme_classic()+
  theme(legend.position = "none")

scatter_2

ggsave(scatter_2, file="ddG_6JPD_cellsurvival_Wuetal2021.jpg", width=4, height=4, path=path)
ggsave(scatter_2, file="ddG_6JPD_cellsurvival_Wuetal2021.pdf", width=4, height=4, path=path)

#################################################################################
# 4: ddG vs NS with mutants from literature

singles_death<-full_join(NS_ddG_mRIP3[c("ID", "ddG_foldx_6JPD","nscore_c")], mRIP3_mutants_HeLa, by="ID")

singles_death$ID_reported<-""
singles_death[is.na(singles_death$cell_survival) == F, "ID_reported"]<-singles_death[is.na(singles_death$cell_survival) == F, "ID"]

scatter_3<-ggplot(singles_death, aes(x=nscore_c, y=ddG_foldx_6JPD, label=ID_reported))+
  geom_hline(yintercept=2, linetype="dashed")+
  geom_point(data=singles_death[singles_death$ID_reported == "",], alpha=0.3, color="grey60")+
  geom_point(data=singles_death[singles_death$ID_reported != "",], aes(color=cell_survival), size=2)+
  scale_color_continuous()+
  geom_text_repel(max.overlaps = 100)+
  theme_classic()
  #theme(legend.position = "none")

scatter_3

ggsave(scatter_3, file="ddG_6JPD_NS_cellsurvival_Wuetal2021.jpg", width=5, height=4, path=path)
ggsave(scatter_3, file="ddG_6JPD_NS_cellsurvival_Wuetal2021.pdf", width=5, height=4, path=path)

#################################################################################
# 5: ThT data

# Open raw data from Wu et al 2021
ThT_data<-read_tsv("Wu2021_ThT_data.txt")
# substract the blank
for (i in c("Wild type - rep1", "Wild type - rep2", "Wild type - rep3", "mRIPK3-F442D-rep1",
            "mRIPK3-F442D-rep2", "mRIPK3-F442D-rep3", "mRIPK3-Q449D-rep1", "mRIPK3-Q449D-rep2", 
            "mRIPK3-Q449D-rep3", "mRIPK3-L456D-rep1", "mRIPK3-L456D-rep2", "mRIPK3-L456D-rep3")){
  ThT_data[paste0(i, "- corrected")]<-ThT_data[i]-ThT_data["Blank"]
}

ThT_data_final<-ThT_data[c("Time(min)", "Wild type - rep1- corrected", "Wild type - rep2- corrected",
                           "Wild type - rep3- corrected", "mRIPK3-F442D-rep1- corrected",
                           "mRIPK3-F442D-rep2- corrected", "mRIPK3-F442D-rep3- corrected", 
                           "mRIPK3-Q449D-rep1- corrected", "mRIPK3-Q449D-rep2- corrected", 
                           "mRIPK3-Q449D-rep3- corrected", "mRIPK3-L456D-rep1- corrected", 
                           "mRIPK3-L456D-rep2- corrected", "mRIPK3-L456D-rep3- corrected")]

# melt the dataframe
ThT_data_melted<-melt(ThT_data_final, id="Time(min)")
# drop rows with na in value
ThT_data_melted<-ThT_data_melted %>% drop_na(value)
# rename Time column
ThT_data_melted<-ThT_data_melted %>% rename(Time = `Time(min)`)
# Tag the resp with same ID
ThT_data_melted$ID<-"WT"
ThT_data_melted[ThT_data_melted$variable %in% c("mRIPK3-F442D-rep1- corrected", "mRIPK3-F442D-rep2- corrected", "mRIPK3-F442D-rep3- corrected"), "ID"]<-"F-442-D"
ThT_data_melted[ThT_data_melted$variable %in% c("mRIPK3-Q449D-rep1- corrected", "mRIPK3-Q449D-rep2- corrected", "mRIPK3-Q449D-rep3- corrected"), "ID"]<-"Q-449-D"
ThT_data_melted[ThT_data_melted$variable %in% c("mRIPK3-L456D-rep1- corrected", "mRIPK3-L456D-rep2- corrected", "mRIPK3-L456D-rep3- corrected"), "ID"]<-"L-456-D"

# Get the mean and the sd for each ID
ThT_data_melted_grouped <- ThT_data_melted %>% group_by(Time, ID) %>% mutate(mean_value=mean(value), std_value=sd(value))

# Plot it
tht_p<-ggplot(ThT_data_melted_grouped, aes(x=Time, y=mean_value, color=ID))+
  geom_point()+
  geom_errorbar(aes(ymin=mean_value-std_value, ymax=mean_value+std_value), width=.2, alpha=0.1)+
  #scale_colour_brewer(palette="Set2")+
  scale_colour_manual(values=c("#D66262", "#7979BE", "#CD3333", "grey80"))+
  geom_smooth(se=F)+
  xlim(0, 210)+
  theme_classic()+
  labs(x="Time (min)", y="ThT Fluorescence", color="")
tht_p

ggsave(tht_p, file="ThT_mRIP3_asitis_Wu2021.jpg", width=5.5, height=4, path=path)
ggsave(tht_p, file="ThT_mRIP3_asitis_Wu2021.pdf", width=5.5, height=4, path=path)

###
# Normalized to each trace
ThT_data_melted_grouped <- ThT_data_melted_grouped %>% group_by(variable) %>% mutate(value_norm=(value/max(value))*100)
ThT_data_melted_grouped <- ThT_data_melted_grouped %>% group_by(Time, ID) %>% mutate(mean_value_norm=mean(value_norm), std_value_norm=sd(value_norm))
ThT_data_melted_grouped <- ThT_data_melted_grouped %>% group_by(variable) %>% mutate(value_trace=value_norm-min(value_norm))
ThT_data_melted_grouped <- ThT_data_melted_grouped %>% group_by(Time, ID) %>% mutate(mean_value_trace=mean(value_trace), std_value_trace=sd(value_trace))
ThT_data_melted_grouped <- ThT_data_melted_grouped %>% group_by(variable) %>% mutate(value_trace_norm=(value_trace/max(value_trace))*100)
ThT_data_melted_grouped <- ThT_data_melted_grouped %>% group_by(Time, ID) %>% mutate(mean_value_trace_norm=mean(value_trace_norm), std_value_trace_norm=sd(value_trace_norm))

scatter_tht_normalized_trace<-ggplot(ThT_data_melted_grouped, 
                                     aes(y=mean_value_trace_norm, x=Time, 
                                         color=factor(ID)))+
  geom_point()+
  geom_errorbar(aes(ymin=mean_value_trace_norm-std_value_trace_norm, ymax=mean_value_trace_norm+std_value_trace_norm), width=.2,
                alpha=0.05)+
  geom_smooth(se=F)+
  scale_y_continuous(breaks=c(0, 100, 200))+
  theme_bw()+
  labs(x="Time (min)", y="Normalized\n ThT Fluorescence", color="")
scatter_tht_normalized_trace

ggsave(scatter_tht_normalized_trace, file="ThT_mRIP3_Wu2021.jpg", width=8, height=6, path=path)
ggsave(scatter_tht_normalized_trace, file="ThT_mRIP3_Wu2021.pdf", width=8, height=6, path=path)

#
scatter_tht_normalized_trace<-ggplot(ThT_data_melted_grouped, 
                                     aes(y=value_trace_norm, x=Time, 
                                         color=factor(variable)))+
  geom_point()+
  geom_smooth(se=F)+
  facet_wrap(~ID)+
  scale_y_continuous(breaks=c(0, 100, 200))+
  theme_bw()+
  labs(x="Time (min)", y="Normalized\n ThT Fluorescence", color="")
scatter_tht_normalized_trace

ggsave(scatter_tht_normalized_trace, file="ThT_mRIP3_reps_Wu2021.jpg", width=8, height=6, path=path)
ggsave(scatter_tht_normalized_trace, file="ThT_mRIP3_reps_Wu2021.pdf", width=8, height=6, path=path)

# wt fit
loess_model_wt_1 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<200 & ThT_data_melted_grouped$variable == "Wild type - rep1- corrected",])
t1_2_wt_1<-predict(loess_model_wt_1, 50)
print(t1_2_wt_1)

loess_model_wt_2 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<200 & ThT_data_melted_grouped$variable == "Wild type - rep2- corrected",])
t1_2_wt_2<-predict(loess_model_wt_2, 50)
print(t1_2_wt_2)

loess_model_wt_3 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<200 & ThT_data_melted_grouped$variable == "Wild type - rep3- corrected",])
t1_2_wt_3<-predict(loess_model_wt_3, 50)
print(t1_2_wt_3)

mean_t1_2_wt<-mean(c(t1_2_wt_1, t1_2_wt_2, t1_2_wt_3))
sd_t1_2_wt<-sd(c(t1_2_wt_1, t1_2_wt_2, t1_2_wt_3))

# Q449D fit
loess_model_Q449D_1 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-Q449D-rep1- corrected",])
t1_2_Q449D_1<-predict(loess_model_Q449D_1, 50)
print(t1_2_Q449D_1)

loess_model_Q449D_2 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-Q449D-rep2- corrected",])
t1_2_Q449D_2<-predict(loess_model_Q449D_2, 50)
print(t1_2_Q449D_2)

loess_model_Q449D_3 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-Q449D-rep3- corrected",])
t1_2_Q449D_3<-predict(loess_model_Q449D_3, 50)
print(t1_2_Q449D_3)

mean_t1_2_Q449D<-mean(c(t1_2_Q449D_1, t1_2_Q449D_2, t1_2_Q449D_3))
sd_t1_2_Q449D<-sd(c(t1_2_Q449D_1, t1_2_Q449D_2, t1_2_Q449D_3))

# L456D fit
loess_model_L456D_1 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-L456D-rep1- corrected",])
t1_2_L456D_1<-predict(loess_model_L456D_1, 50)
print(t1_2_L456D_1)

loess_model_L456D_2 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-L456D-rep2- corrected",])
t1_2_L456D_2<-predict(loess_model_L456D_2, 50)
print(t1_2_L456D_2)

loess_model_L456D_3 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-L456D-rep3- corrected",])
t1_2_L456D_3<-predict(loess_model_L456D_3, 50)
print(t1_2_L456D_3)

mean_t1_2_L456D<-mean(c(t1_2_L456D_1, t1_2_L456D_2, t1_2_L456D_3))
sd_t1_2_L456D<-sd(c(t1_2_L456D_1, t1_2_L456D_2, t1_2_L456D_3))

# F442D fit
loess_model_F442D_1 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-F442D-rep1- corrected",])
t1_2_F442D_1<-predict(loess_model_F442D_1, 50)
print(t1_2_F442D_1)

loess_model_F442D_2 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-F442D-rep2- corrected",])
t1_2_F442D_2<-predict(loess_model_F442D_2, 50)
print(t1_2_F442D_2)

loess_model_F442D_3 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-F442D-rep3- corrected",])
t1_2_F442D_3<-predict(loess_model_F442D_3, 50)
print(t1_2_F442D_3)

mean_t1_2_F442D<-mean(c(t1_2_F442D_1, t1_2_F442D_2, t1_2_F442D_3))
sd_t1_2_F442D<-sd(c(t1_2_F442D_1, t1_2_F442D_2, t1_2_F442D_3))

#
kinetics_df<-data.frame(ID=c("WT", "Q-449-D", "L-456-D", "F-442-D"),
                        mean_t1_2=c(round(mean_t1_2_wt, 2), round(mean_t1_2_Q449D, 2), round(mean_t1_2_L456D, 2), round(mean_t1_2_F442D, 2)),
                        sd_t1_2=c(round(sd_t1_2_wt, 2), round(sd_t1_2_Q449D, 2), round(sd_t1_2_L456D, 2), round(sd_t1_2_F442D, 2))
)

kinetics_df<-left_join(kinetics_df, NS_ddG_mRIP3, by="ID")
kinetics_df<-kinetics_df[c("ID", "mean_t1_2", "sd_t1_2", "nscore_c", "sigma", "ddG_foldx_6JPD", "ddG_foldx_SD_6JPD")]

kinetics_df[kinetics_df$ID == "WT", "nscore_c"]<-0.056252746
kinetics_df[kinetics_df$ID == "WT", "sigma"]<-0.01378721
kinetics_df[kinetics_df$ID == "WT", "ddG_foldx_6JPD"]<-0
kinetics_df[kinetics_df$ID == "WT", "ddG_foldx_SD_6JPD"]<-0


kinetics_fitness_plot<-ggplot(kinetics_df, aes(x=mean_t1_2, y=nscore_c))+
  geom_point(size=3)+
  geom_errorbar(aes(xmin=mean_t1_2-sd_t1_2, xmax=mean_t1_2+sd_t1_2), width=0.2, alpha=0.5)+
  geom_errorbar(aes(ymin=nscore_c-sigma, ymax=nscore_c+sigma), width=0.2, alpha=0.5)+
  geom_text_repel(aes(label=ID))+
  geom_smooth(method = "lm", se = F, color="black", linetype = "dashed")+
  stat_cor(method="pearson")+
  labs(title="T1/2 estimation")+
  theme_classic()
kinetics_fitness_plot

ggsave(kinetics_fitness_plot, file="ThT_NS_halftime.jpg", width=4, height=4, path=path)
ggsave(kinetics_fitness_plot, file="ThT_NS_halftime.pdf", width=4, height=4, path=path)

###
# cutting the traces at 150 min
ThT_data_melted_grouped<-ThT_data_melted_grouped[ThT_data_melted_grouped$Time<150,]
ThT_data_melted_grouped <- ThT_data_melted_grouped %>% group_by(variable) %>% mutate(value_norm=(value/max(value))*100)
ThT_data_melted_grouped <- ThT_data_melted_grouped %>% group_by(Time, ID) %>% mutate(mean_value_norm=mean(value_norm), std_value_norm=sd(value_norm))
ThT_data_melted_grouped <- ThT_data_melted_grouped %>% group_by(variable) %>% mutate(value_trace=value_norm-min(value_norm))
ThT_data_melted_grouped <- ThT_data_melted_grouped %>% group_by(Time, ID) %>% mutate(mean_value_trace=mean(value_trace), std_value_trace=sd(value_trace))
ThT_data_melted_grouped <- ThT_data_melted_grouped %>% group_by(variable) %>% mutate(value_trace_norm=(value_trace/max(value_trace))*100)
ThT_data_melted_grouped <- ThT_data_melted_grouped %>% group_by(Time, ID) %>% mutate(mean_value_trace_norm=mean(value_trace_norm), std_value_trace_norm=sd(value_trace_norm))

# wt fit
loess_model_wt_1 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<200 & ThT_data_melted_grouped$variable == "Wild type - rep1- corrected",])
t1_2_wt_1<-predict(loess_model_wt_1, 50)
print(t1_2_wt_1)

loess_model_wt_2 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<200 & ThT_data_melted_grouped$variable == "Wild type - rep2- corrected",])
t1_2_wt_2<-predict(loess_model_wt_2, 50)
print(t1_2_wt_2)

loess_model_wt_3 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<200 & ThT_data_melted_grouped$variable == "Wild type - rep3- corrected",])
t1_2_wt_3<-predict(loess_model_wt_3, 50)
print(t1_2_wt_3)

mean_t1_2_wt<-mean(c(t1_2_wt_1, t1_2_wt_2, t1_2_wt_3))
sd_t1_2_wt<-sd(c(t1_2_wt_1, t1_2_wt_2, t1_2_wt_3))

# Q449D fit
loess_model_Q449D_1 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-Q449D-rep1- corrected",])
t1_2_Q449D_1<-predict(loess_model_Q449D_1, 50)
print(t1_2_Q449D_1)

loess_model_Q449D_2 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-Q449D-rep2- corrected",])
t1_2_Q449D_2<-predict(loess_model_Q449D_2, 50)
print(t1_2_Q449D_2)

loess_model_Q449D_3 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-Q449D-rep3- corrected",])
t1_2_Q449D_3<-predict(loess_model_Q449D_3, 50)
print(t1_2_Q449D_3)

mean_t1_2_Q449D<-mean(c(t1_2_Q449D_1, t1_2_Q449D_2, t1_2_Q449D_3))
sd_t1_2_Q449D<-sd(c(t1_2_Q449D_1, t1_2_Q449D_2, t1_2_Q449D_3))

# L456D fit
loess_model_L456D_1 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-L456D-rep1- corrected",])
t1_2_L456D_1<-predict(loess_model_L456D_1, 50)
print(t1_2_L456D_1)

loess_model_L456D_2 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-L456D-rep2- corrected",])
t1_2_L456D_2<-predict(loess_model_L456D_2, 50)
print(t1_2_L456D_2)

loess_model_L456D_3 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-L456D-rep3- corrected",])
t1_2_L456D_3<-predict(loess_model_L456D_3, 50)
print(t1_2_L456D_3)

mean_t1_2_L456D<-mean(c(t1_2_L456D_1, t1_2_L456D_2, t1_2_L456D_3))
sd_t1_2_L456D<-sd(c(t1_2_L456D_1, t1_2_L456D_2, t1_2_L456D_3))

# F442D fit
loess_model_F442D_1 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-F442D-rep1- corrected",])
t1_2_F442D_1<-predict(loess_model_F442D_1, 50)
print(t1_2_F442D_1)

loess_model_F442D_2 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-F442D-rep2- corrected",])
t1_2_F442D_2<-predict(loess_model_F442D_2, 50)
print(t1_2_F442D_2)

loess_model_F442D_3 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-F442D-rep3- corrected",])
t1_2_F442D_3<-predict(loess_model_F442D_3, 50)
print(t1_2_F442D_3)

mean_t1_2_F442D<-mean(c(t1_2_F442D_1, t1_2_F442D_2, t1_2_F442D_3))
sd_t1_2_F442D<-sd(c(t1_2_F442D_1, t1_2_F442D_2, t1_2_F442D_3))

#
kinetics_df<-data.frame(ID=c("WT", "Q-449-D", "L-456-D", "F-442-D"),
                        mean_t1_2=c(round(mean_t1_2_wt, 2), round(mean_t1_2_Q449D, 2), round(mean_t1_2_L456D, 2), round(mean_t1_2_F442D, 2)),
                        sd_t1_2=c(round(sd_t1_2_wt, 2), round(sd_t1_2_Q449D, 2), round(sd_t1_2_L456D, 2), round(sd_t1_2_F442D, 2))
)

kinetics_df<-left_join(kinetics_df, NS_ddG_mRIP3, by="ID")
kinetics_df<-kinetics_df[c("ID", "mean_t1_2", "sd_t1_2", "nscore_c", "sigma", "ddG_foldx_6JPD", "ddG_foldx_SD_6JPD")]

kinetics_df[kinetics_df$ID == "WT", "nscore_c"]<-0.056252746
kinetics_df[kinetics_df$ID == "WT", "sigma"]<-0.01378721
kinetics_df[kinetics_df$ID == "WT", "ddG_foldx_6JPD"]<-0
kinetics_df[kinetics_df$ID == "WT", "ddG_foldx_SD_6JPD"]<-0


kinetics_fitness_plot<-ggplot(kinetics_df, aes(x=mean_t1_2, y=nscore_c))+
  geom_point(size=3)+
  geom_errorbar(aes(xmin=mean_t1_2-sd_t1_2, xmax=mean_t1_2+sd_t1_2), width=0.2, alpha=0.5)+
  geom_errorbar(aes(ymin=nscore_c-sigma, ymax=nscore_c+sigma), width=0.2, alpha=0.5)+
  geom_text_repel(aes(label=ID))+
  geom_smooth(method = "lm", se = F, color="black", linetype = "dashed")+
  stat_cor(method="pearson")+
  labs(title="T1/2 estimation until 150 min")+
  theme_classic()
kinetics_fitness_plot

ggsave(kinetics_fitness_plot, file="ThT_NS_halftime_150.jpg", width=4, height=4, path=path)
ggsave(kinetics_fitness_plot, file="ThT_NS_halftime_150.pdf", width=4, height=4, path=path)

###
# cutting the traces at 75 min
ThT_data_melted_grouped<-ThT_data_melted_grouped[ThT_data_melted_grouped$Time<75,]
ThT_data_melted_grouped <- ThT_data_melted_grouped %>% group_by(variable) %>% mutate(value_norm=(value/max(value))*100)
ThT_data_melted_grouped <- ThT_data_melted_grouped %>% group_by(Time, ID) %>% mutate(mean_value_norm=mean(value_norm), std_value_norm=sd(value_norm))
ThT_data_melted_grouped <- ThT_data_melted_grouped %>% group_by(variable) %>% mutate(value_trace=value_norm-min(value_norm))
ThT_data_melted_grouped <- ThT_data_melted_grouped %>% group_by(Time, ID) %>% mutate(mean_value_trace=mean(value_trace), std_value_trace=sd(value_trace))
ThT_data_melted_grouped <- ThT_data_melted_grouped %>% group_by(variable) %>% mutate(value_trace_norm=(value_trace/max(value_trace))*100)
ThT_data_melted_grouped <- ThT_data_melted_grouped %>% group_by(Time, ID) %>% mutate(mean_value_trace_norm=mean(value_trace_norm), std_value_trace_norm=sd(value_trace_norm))

# wt fit
loess_model_wt_1 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<200 & ThT_data_melted_grouped$variable == "Wild type - rep1- corrected",])
t1_2_wt_1<-predict(loess_model_wt_1, 50)
print(t1_2_wt_1)

loess_model_wt_2 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<200 & ThT_data_melted_grouped$variable == "Wild type - rep2- corrected",])
t1_2_wt_2<-predict(loess_model_wt_2, 50)
print(t1_2_wt_2)

loess_model_wt_3 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<200 & ThT_data_melted_grouped$variable == "Wild type - rep3- corrected",])
t1_2_wt_3<-predict(loess_model_wt_3, 50)
print(t1_2_wt_3)

mean_t1_2_wt<-mean(c(t1_2_wt_1, t1_2_wt_2, t1_2_wt_3))
sd_t1_2_wt<-sd(c(t1_2_wt_1, t1_2_wt_2, t1_2_wt_3))

# Q449D fit
loess_model_Q449D_1 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-Q449D-rep1- corrected",])
t1_2_Q449D_1<-predict(loess_model_Q449D_1, 50)
print(t1_2_Q449D_1)

loess_model_Q449D_2 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-Q449D-rep2- corrected",])
t1_2_Q449D_2<-predict(loess_model_Q449D_2, 50)
print(t1_2_Q449D_2)

loess_model_Q449D_3 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-Q449D-rep3- corrected",])
t1_2_Q449D_3<-predict(loess_model_Q449D_3, 50)
print(t1_2_Q449D_3)

mean_t1_2_Q449D<-mean(c(t1_2_Q449D_1, t1_2_Q449D_2, t1_2_Q449D_3))
sd_t1_2_Q449D<-sd(c(t1_2_Q449D_1, t1_2_Q449D_2, t1_2_Q449D_3))

# L456D fit
loess_model_L456D_1 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-L456D-rep1- corrected",])
t1_2_L456D_1<-predict(loess_model_L456D_1, 50)
print(t1_2_L456D_1)

loess_model_L456D_2 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-L456D-rep2- corrected",])
t1_2_L456D_2<-predict(loess_model_L456D_2, 50)
print(t1_2_L456D_2)

loess_model_L456D_3 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-L456D-rep3- corrected",])
t1_2_L456D_3<-predict(loess_model_L456D_3, 50)
print(t1_2_L456D_3)

mean_t1_2_L456D<-mean(c(t1_2_L456D_1, t1_2_L456D_2, t1_2_L456D_3))
sd_t1_2_L456D<-sd(c(t1_2_L456D_1, t1_2_L456D_2, t1_2_L456D_3))

# F442D fit
loess_model_F442D_1 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-F442D-rep1- corrected",])
t1_2_F442D_1<-predict(loess_model_F442D_1, 50)
print(t1_2_F442D_1)

loess_model_F442D_2 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-F442D-rep2- corrected",])
t1_2_F442D_2<-predict(loess_model_F442D_2, 50)
print(t1_2_F442D_2)

loess_model_F442D_3 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-F442D-rep3- corrected",])
t1_2_F442D_3<-predict(loess_model_F442D_3, 50)
print(t1_2_F442D_3)

mean_t1_2_F442D<-mean(c(t1_2_F442D_1, t1_2_F442D_2, t1_2_F442D_3))
sd_t1_2_F442D<-sd(c(t1_2_F442D_1, t1_2_F442D_2, t1_2_F442D_3))

#
kinetics_df<-data.frame(ID=c("WT", "Q-449-D", "L-456-D", "F-442-D"),
                        mean_t1_2=c(round(mean_t1_2_wt, 2), round(mean_t1_2_Q449D, 2), round(mean_t1_2_L456D, 2), round(mean_t1_2_F442D, 2)),
                        sd_t1_2=c(round(sd_t1_2_wt, 2), round(sd_t1_2_Q449D, 2), round(sd_t1_2_L456D, 2), round(sd_t1_2_F442D, 2))
)

kinetics_df<-left_join(kinetics_df, NS_ddG_mRIP3, by="ID")
kinetics_df<-kinetics_df[c("ID", "mean_t1_2", "sd_t1_2", "nscore_c", "sigma", "ddG_foldx_6JPD", "ddG_foldx_SD_6JPD")]

kinetics_df[kinetics_df$ID == "WT", "nscore_c"]<-0.056252746
kinetics_df[kinetics_df$ID == "WT", "sigma"]<-0.01378721
kinetics_df[kinetics_df$ID == "WT", "ddG_foldx_6JPD"]<-0
kinetics_df[kinetics_df$ID == "WT", "ddG_foldx_SD_6JPD"]<-0


kinetics_fitness_plot<-ggplot(kinetics_df, aes(x=mean_t1_2, y=nscore_c))+
  geom_point(size=3)+
  geom_errorbar(aes(xmin=mean_t1_2-sd_t1_2, xmax=mean_t1_2+sd_t1_2), width=0.2, alpha=0.5)+
  geom_errorbar(aes(ymin=nscore_c-sigma, ymax=nscore_c+sigma), width=0.2, alpha=0.5)+
  geom_text_repel(aes(label=ID))+
  geom_smooth(method = "lm", se = F, color="black", linetype = "dashed")+
  stat_cor(method="pearson")+
  labs(title="T1/2 estimation until 75 min")+
  theme_classic()
kinetics_fitness_plot

ggsave(kinetics_fitness_plot, file="ThT_NS_halftime_75.jpg", width=4, height=4, path=path)
ggsave(kinetics_fitness_plot, file="ThT_NS_halftime_75.pdf", width=4, height=4, path=path)

###
# cutting the traces at 50 min
ThT_data_melted_grouped<-ThT_data_melted_grouped[ThT_data_melted_grouped$Time<50,]
ThT_data_melted_grouped <- ThT_data_melted_grouped %>% group_by(variable) %>% mutate(value_norm=(value/max(value))*100)
ThT_data_melted_grouped <- ThT_data_melted_grouped %>% group_by(Time, ID) %>% mutate(mean_value_norm=mean(value_norm), std_value_norm=sd(value_norm))
ThT_data_melted_grouped <- ThT_data_melted_grouped %>% group_by(variable) %>% mutate(value_trace=value_norm-min(value_norm))
ThT_data_melted_grouped <- ThT_data_melted_grouped %>% group_by(Time, ID) %>% mutate(mean_value_trace=mean(value_trace), std_value_trace=sd(value_trace))
ThT_data_melted_grouped <- ThT_data_melted_grouped %>% group_by(variable) %>% mutate(value_trace_norm=(value_trace/max(value_trace))*100)
ThT_data_melted_grouped <- ThT_data_melted_grouped %>% group_by(Time, ID) %>% mutate(mean_value_trace_norm=mean(value_trace_norm), std_value_trace_norm=sd(value_trace_norm))

# wt fit
loess_model_wt_1 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<200 & ThT_data_melted_grouped$variable == "Wild type - rep1- corrected",])
t1_2_wt_1<-predict(loess_model_wt_1, 50)
print(t1_2_wt_1)

loess_model_wt_2 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<200 & ThT_data_melted_grouped$variable == "Wild type - rep2- corrected",])
t1_2_wt_2<-predict(loess_model_wt_2, 50)
print(t1_2_wt_2)

loess_model_wt_3 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<200 & ThT_data_melted_grouped$variable == "Wild type - rep3- corrected",])
t1_2_wt_3<-predict(loess_model_wt_3, 50)
print(t1_2_wt_3)

mean_t1_2_wt<-mean(c(t1_2_wt_1, t1_2_wt_2, t1_2_wt_3))
sd_t1_2_wt<-sd(c(t1_2_wt_1, t1_2_wt_2, t1_2_wt_3))

# Q449D fit
loess_model_Q449D_1 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-Q449D-rep1- corrected",])
t1_2_Q449D_1<-predict(loess_model_Q449D_1, 50)
print(t1_2_Q449D_1)

loess_model_Q449D_2 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-Q449D-rep2- corrected",])
t1_2_Q449D_2<-predict(loess_model_Q449D_2, 50)
print(t1_2_Q449D_2)

loess_model_Q449D_3 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-Q449D-rep3- corrected",])
t1_2_Q449D_3<-predict(loess_model_Q449D_3, 50)
print(t1_2_Q449D_3)

mean_t1_2_Q449D<-mean(c(t1_2_Q449D_1, t1_2_Q449D_2, t1_2_Q449D_3))
sd_t1_2_Q449D<-sd(c(t1_2_Q449D_1, t1_2_Q449D_2, t1_2_Q449D_3))

# L456D fit
loess_model_L456D_1 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-L456D-rep1- corrected",])
t1_2_L456D_1<-predict(loess_model_L456D_1, 50)
print(t1_2_L456D_1)

loess_model_L456D_2 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-L456D-rep2- corrected",])
t1_2_L456D_2<-predict(loess_model_L456D_2, 50)
print(t1_2_L456D_2)

loess_model_L456D_3 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-L456D-rep3- corrected",])
t1_2_L456D_3<-predict(loess_model_L456D_3, 50)
print(t1_2_L456D_3)

mean_t1_2_L456D<-mean(c(t1_2_L456D_1, t1_2_L456D_2, t1_2_L456D_3))
sd_t1_2_L456D<-sd(c(t1_2_L456D_1, t1_2_L456D_2, t1_2_L456D_3))

# F442D fit
loess_model_F442D_1 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-F442D-rep1- corrected",])
t1_2_F442D_1<-predict(loess_model_F442D_1, 50)
print(t1_2_F442D_1)

loess_model_F442D_2 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-F442D-rep2- corrected",])
t1_2_F442D_2<-predict(loess_model_F442D_2, 50)
print(t1_2_F442D_2)

loess_model_F442D_3 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-F442D-rep3- corrected",])
t1_2_F442D_3<-predict(loess_model_F442D_3, 50)
print(t1_2_F442D_3)

mean_t1_2_F442D<-mean(c(t1_2_F442D_1, t1_2_F442D_2, t1_2_F442D_3))
sd_t1_2_F442D<-sd(c(t1_2_F442D_1, t1_2_F442D_2, t1_2_F442D_3))

#
kinetics_df<-data.frame(ID=c("WT", "Q-449-D", "L-456-D", "F-442-D"),
                        mean_t1_2=c(round(mean_t1_2_wt, 2), round(mean_t1_2_Q449D, 2), round(mean_t1_2_L456D, 2), round(mean_t1_2_F442D, 2)),
                        sd_t1_2=c(round(sd_t1_2_wt, 2), round(sd_t1_2_Q449D, 2), round(sd_t1_2_L456D, 2), round(sd_t1_2_F442D, 2))
)

kinetics_df<-left_join(kinetics_df, NS_ddG_mRIP3, by="ID")
kinetics_df<-kinetics_df[c("ID", "mean_t1_2", "sd_t1_2", "nscore_c", "sigma", "ddG_foldx_6JPD", "ddG_foldx_SD_6JPD")]

kinetics_df[kinetics_df$ID == "WT", "nscore_c"]<-0.056252746
kinetics_df[kinetics_df$ID == "WT", "sigma"]<-0.01378721
kinetics_df[kinetics_df$ID == "WT", "ddG_foldx_6JPD"]<-0
kinetics_df[kinetics_df$ID == "WT", "ddG_foldx_SD_6JPD"]<-0


kinetics_fitness_plot<-ggplot(kinetics_df, aes(x=mean_t1_2, y=nscore_c))+
  geom_point(size=3)+
  geom_errorbar(aes(xmin=mean_t1_2-sd_t1_2, xmax=mean_t1_2+sd_t1_2), width=0.2, alpha=0.5)+
  geom_errorbar(aes(ymin=nscore_c-sigma, ymax=nscore_c+sigma), width=0.2, alpha=0.5)+
  geom_text_repel(aes(label=ID))+
  geom_smooth(method = "lm", se = F, color="black", linetype = "dashed")+
  stat_cor(method="pearson")+
  labs(title="T1/2 estimation until 50 min")+
  theme_classic()
kinetics_fitness_plot

ggsave(kinetics_fitness_plot, file="ThT_NS_halftime_50.jpg", width=4, height=4, path=path)
ggsave(kinetics_fitness_plot, file="ThT_NS_halftime_50.pdf", width=4, height=4, path=path)

###
# cutting the traces at 25 min
ThT_data_melted_grouped<-ThT_data_melted_grouped[ThT_data_melted_grouped$Time<25,]
ThT_data_melted_grouped <- ThT_data_melted_grouped %>% group_by(variable) %>% mutate(value_norm=(value/max(value))*100)
ThT_data_melted_grouped <- ThT_data_melted_grouped %>% group_by(Time, ID) %>% mutate(mean_value_norm=mean(value_norm), std_value_norm=sd(value_norm))
ThT_data_melted_grouped <- ThT_data_melted_grouped %>% group_by(variable) %>% mutate(value_trace=value_norm-min(value_norm))
ThT_data_melted_grouped <- ThT_data_melted_grouped %>% group_by(Time, ID) %>% mutate(mean_value_trace=mean(value_trace), std_value_trace=sd(value_trace))
ThT_data_melted_grouped <- ThT_data_melted_grouped %>% group_by(variable) %>% mutate(value_trace_norm=(value_trace/max(value_trace))*100)
ThT_data_melted_grouped <- ThT_data_melted_grouped %>% group_by(Time, ID) %>% mutate(mean_value_trace_norm=mean(value_trace_norm), std_value_trace_norm=sd(value_trace_norm))

# wt fit
loess_model_wt_1 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<200 & ThT_data_melted_grouped$variable == "Wild type - rep1- corrected",])
t1_2_wt_1<-predict(loess_model_wt_1, 50)
print(t1_2_wt_1)

loess_model_wt_2 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<200 & ThT_data_melted_grouped$variable == "Wild type - rep2- corrected",])
t1_2_wt_2<-predict(loess_model_wt_2, 50)
print(t1_2_wt_2)

loess_model_wt_3 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<200 & ThT_data_melted_grouped$variable == "Wild type - rep3- corrected",])
t1_2_wt_3<-predict(loess_model_wt_3, 50)
print(t1_2_wt_3)

mean_t1_2_wt<-mean(c(t1_2_wt_1, t1_2_wt_2, t1_2_wt_3))
sd_t1_2_wt<-sd(c(t1_2_wt_1, t1_2_wt_2, t1_2_wt_3))

# Q449D fit
loess_model_Q449D_1 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-Q449D-rep1- corrected",])
t1_2_Q449D_1<-predict(loess_model_Q449D_1, 50)
print(t1_2_Q449D_1)

loess_model_Q449D_2 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-Q449D-rep2- corrected",])
t1_2_Q449D_2<-predict(loess_model_Q449D_2, 50)
print(t1_2_Q449D_2)

loess_model_Q449D_3 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-Q449D-rep3- corrected",])
t1_2_Q449D_3<-predict(loess_model_Q449D_3, 50)
print(t1_2_Q449D_3)

mean_t1_2_Q449D<-mean(c(t1_2_Q449D_1, t1_2_Q449D_2, t1_2_Q449D_3))
sd_t1_2_Q449D<-sd(c(t1_2_Q449D_1, t1_2_Q449D_2, t1_2_Q449D_3))

# L456D fit
loess_model_L456D_1 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-L456D-rep1- corrected",])
t1_2_L456D_1<-predict(loess_model_L456D_1, 50)
print(t1_2_L456D_1)

loess_model_L456D_2 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-L456D-rep2- corrected",])
t1_2_L456D_2<-predict(loess_model_L456D_2, 50)
print(t1_2_L456D_2)

loess_model_L456D_3 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-L456D-rep3- corrected",])
t1_2_L456D_3<-predict(loess_model_L456D_3, 50)
print(t1_2_L456D_3)

mean_t1_2_L456D<-mean(c(t1_2_L456D_1, t1_2_L456D_2, t1_2_L456D_3))
sd_t1_2_L456D<-sd(c(t1_2_L456D_1, t1_2_L456D_2, t1_2_L456D_3))

# F442D fit
loess_model_F442D_1 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-F442D-rep1- corrected",])
t1_2_F442D_1<-predict(loess_model_F442D_1, 50)
print(t1_2_F442D_1)

loess_model_F442D_2 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-F442D-rep2- corrected",])
t1_2_F442D_2<-predict(loess_model_F442D_2, 50)
print(t1_2_F442D_2)

loess_model_F442D_3 <- loess(Time ~ value_trace_norm, data = ThT_data_melted_grouped[ThT_data_melted_grouped$Time<220 & ThT_data_melted_grouped$variable == "mRIPK3-F442D-rep3- corrected",])
t1_2_F442D_3<-predict(loess_model_F442D_3, 50)
print(t1_2_F442D_3)

mean_t1_2_F442D<-mean(c(t1_2_F442D_1, t1_2_F442D_2, t1_2_F442D_3))
sd_t1_2_F442D<-sd(c(t1_2_F442D_1, t1_2_F442D_2, t1_2_F442D_3))

#
kinetics_df<-data.frame(ID=c("WT", "Q-449-D", "L-456-D", "F-442-D"),
                        mean_t1_2=c(round(mean_t1_2_wt, 2), round(mean_t1_2_Q449D, 2), round(mean_t1_2_L456D, 2), round(mean_t1_2_F442D, 2)),
                        sd_t1_2=c(round(sd_t1_2_wt, 2), round(sd_t1_2_Q449D, 2), round(sd_t1_2_L456D, 2), round(sd_t1_2_F442D, 2))
)

kinetics_df<-left_join(kinetics_df, NS_ddG_mRIP3, by="ID")
kinetics_df<-kinetics_df[c("ID", "mean_t1_2", "sd_t1_2", "nscore_c", "sigma", "ddG_foldx_6JPD", "ddG_foldx_SD_6JPD")]

kinetics_df[kinetics_df$ID == "WT", "nscore_c"]<-0.056252746
kinetics_df[kinetics_df$ID == "WT", "sigma"]<-0.01378721
kinetics_df[kinetics_df$ID == "WT", "ddG_foldx_6JPD"]<-0
kinetics_df[kinetics_df$ID == "WT", "ddG_foldx_SD_6JPD"]<-0


kinetics_fitness_plot<-ggplot(kinetics_df, aes(x=mean_t1_2, y=nscore_c))+
  geom_point(size=3)+
  geom_errorbar(aes(xmin=mean_t1_2-sd_t1_2, xmax=mean_t1_2+sd_t1_2), width=0.2, alpha=0.5)+
  geom_errorbar(aes(ymin=nscore_c-sigma, ymax=nscore_c+sigma), width=0.2, alpha=0.5)+
  geom_text_repel(aes(label=ID))+
  geom_smooth(method = "lm", se = F, color="black", linetype = "dashed")+
  stat_cor(method="pearson")+
  labs(title="T1/2 estimation until 25 min")+
  theme_classic()
kinetics_fitness_plot

ggsave(kinetics_fitness_plot, file="ThT_NS_halftime_25.jpg", width=4, height=4, path=path)
ggsave(kinetics_fitness_plot, file="ThT_NS_halftime_25.pdf", width=4, height=4, path=path)
