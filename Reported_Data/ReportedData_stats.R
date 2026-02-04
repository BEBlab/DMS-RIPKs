# Script for testing the data reported in the literature
library(multcomp)

# t.test without reps data
t_test_summary <- function(ID1, mean1, sd1, n1, ID2, mean2, sd2, n2) {
  se_diff <- sqrt(sd1^2/n1 + sd2^2/n2)
  t_value <- (mean1 - mean2) / se_diff
  df <- ((sd1^2/n1 + sd2^2/n2)^2) /
    ((sd1^2/n1)^2/(n1-1) + (sd2^2/n2)^2/(n2-1))
  p_value <- 2 * pt(-abs(t_value), df)
  return(data.frame(t_value, df, p_value))
}

###
# Li et a. 2012 Cell RIPK3
df <- data.frame("ID"=factor(c("WT", "V-458-P", "V-460-P", "V458P/V460P", "G-457-D", "N-464-D", "VQVG>AAAA"), levels=c("WT", "V-458-P", "V-460-P", "V458P/V460P", "G-457-D", "N-464-D", "VQVG>AAAA")),
                 "mean"=c(24.03, 11.76, 6.14, 8.86, 12.44, 18.24, 6.65),
                 "sd"=c(1.54, 2.39, 1.7, 2.22, 4.26, 3.58, 0.51),
                 "n"=c(3, 3, 3, 3, 3, 3, 3))
# t.test
for (id in df$ID){
  result_t_test<-t_test_summary("WT", df[df$ID == "WT", "mean"], df[df$ID == "WT", "sd"], df[df$ID == "WT", "n"],
                                df[df$ID == id, "ID"], df[df$ID == id, "mean"], df[df$ID == id, "sd"], df[df$ID == id, "n"])
  print(paste("WT vs", id))
  print(result_t_test)
  print("***********************************************************************")
}

##
# simulated ANOVA with Tukey posthoc
set.seed(123)
sim_data <- df %>%
  group_by(ID) %>%
  do(data.frame(value = rnorm(.$n, mean = .$mean, sd = .$sd)))

# Run ANOVA
anova_res <- aov(value ~ ID, data = sim_data)
summary(anova_res)

# Perform Dunnett's post-hoc test
summary(glht(anova_res, linfct = mcp(ID = "Dunnett")))

# Tukey post-hoc test
TukeyHSD(anova_res)



###
# Li et al. 2012 Cell RIPK1
df <- data.frame("ID"=factor(c("WT", "IQIG>AAAA", "I-539-P", "I-541-P", "N-545-P"), levels=c("WT", "IQIG>AAAA", "I-539-P", "I-541-P", "N-545-P")),
                 "mean"=c(68.83, 25.03, 8.94, 21, 31.73),
                 "sd"=c(0.89, 5.36, 2.68, 2.69, 5.36),
                 "n"=c(3, 3, 3, 3, 3))
# t.test
for (id in df$ID){
  result_t_test<-t_test_summary("WT", df[df$ID == "WT", "mean"], df[df$ID == "WT", "sd"], df[df$ID == "WT", "n"],
                                df[df$ID == id, "ID"], df[df$ID == id, "mean"], df[df$ID == id, "sd"], df[df$ID == id, "n"])
  print(paste("WT vs", id))
  print(result_t_test)
  print("***********************************************************************")
}

##
# simulated ANOVA with Tukey posthoc
set.seed(123)
sim_data <- df %>%
  group_by(ID) %>%
  do(data.frame(value = rnorm(.$n, mean = .$mean, sd = .$sd)))

# Run ANOVA
anova_res <- aov(value ~ ID, data = sim_data)
summary(anova_res)

# Perform Dunnett's post-hoc test
summary(glht(anova_res, linfct = mcp(ID = "Dunnett")))

# Tukey post-hoc test
TukeyHSD(anova_res)

###
# Hu et al. 2020 Cell death diff RIPK3
# timepoint 10 hs
df <- data.frame("ID"=factor(c("WT", "N464D", "M468D", "N464D/M468D"), levels=c("WT", "N464D", "M468D", "N464D/M468D")),
                 "mean"=c(35.48, 25.45, 29.75, 95.34),
                 "sd"=c(2.51, 1.07, 2.15, 3.23),
                 "n"=c(2, 2, 2, 2))
# t.test
for (id in df$ID){
  result_t_test<-t_test_summary("WT", df[df$ID == "WT", "mean"], df[df$ID == "WT", "sd"], df[df$ID == "WT", "n"],
                                df[df$ID == id, "ID"], df[df$ID == id, "mean"], df[df$ID == id, "sd"], df[df$ID == id, "n"])
  print(paste("WT vs", id))
  print(result_t_test)
  print("***********************************************************************")
}

##
# simulated ANOVA with Tukey posthoc
set.seed(123)
sim_data <- df %>%
  group_by(ID) %>%
  do(data.frame(value = rnorm(.$n, mean = .$mean, sd = .$sd)))

# Perform Dunnett's post-hoc test
summary(glht(anova_res, linfct = mcp(ID = "Dunnett")))

# Run ANOVA
anova_res <- aov(value ~ ID, data = sim_data)
summary(anova_res)

# Tukey post-hoc test
TukeyHSD(anova_res)

# timepoint 24 hs
df <- data.frame("ID"=factor(c("WT", "N464D", "M468D", "N464D/M468D"), levels=c("WT", "N464D", "M468D", "N464D/M468D")),
                 "mean"=c(12.19, 12.19, 7.53, 96.77),
                 "sd"=c(1.07, 3.22, 1.79, 2.15),
                 "n"=c(2, 2, 2, 2))
# t.test
for (id in df$ID){
  result_t_test<-t_test_summary("WT", df[df$ID == "WT", "mean"], df[df$ID == "WT", "sd"], df[df$ID == "WT", "n"],
                                df[df$ID == id, "ID"], df[df$ID == id, "mean"], df[df$ID == id, "sd"], df[df$ID == id, "n"])
  print(paste("WT vs", id))
  print(result_t_test)
  print("***********************************************************************")
}

##
# simulated ANOVA with Tukey posthoc
set.seed(123)
sim_data <- df %>%
  group_by(ID) %>%
  do(data.frame(value = rnorm(.$n, mean = .$mean, sd = .$sd)))

# Run ANOVA
anova_res <- aov(value ~ ID, data = sim_data)
summary(anova_res)

# Perform Dunnett's post-hoc test
summary(glht(anova_res, linfct = mcp(ID = "Dunnett")))

# Tukey post-hoc test
TukeyHSD(anova_res)

###
# Ma et al. PNAS 2025
# timepoint 2 hs
df <- data.frame("ID"=factor(c("WT", "V450E", "V460E", "L466E", "V450E/V460E"), levels=c("WT", "V450E", "V460E", "L466E", "V450E/V460E")),
                 "mean"=c(51.74, 80.23, 86.63, 96.22, 87.21),
                 "sd"=c(3.49, 11.05, 5.23, 6.11, 3.2),
                 "n"=c(4, 4, 4, 4, 4))
# t.test
for (id in df$ID){
  result_t_test<-t_test_summary("WT", df[df$ID == "WT", "mean"], df[df$ID == "WT", "sd"], df[df$ID == "WT", "n"],
                                df[df$ID == id, "ID"], df[df$ID == id, "mean"], df[df$ID == id, "sd"], df[df$ID == id, "n"])
  print(paste("WT vs", id))
  print(result_t_test)
  print("***********************************************************************")
}

##
# simulated ANOVA with Tukey posthoc
set.seed(123)
sim_data <- df %>%
  group_by(ID) %>%
  do(data.frame(value = rnorm(.$n, mean = .$mean, sd = .$sd)))

# Run ANOVA
anova_res <- aov(value ~ ID, data = sim_data)
summary(anova_res)

# Perform Dunnett's post-hoc test
summary(glht(anova_res, linfct = mcp(ID = "Dunnett")))

# Tukey post-hoc test
TukeyHSD(anova_res)


