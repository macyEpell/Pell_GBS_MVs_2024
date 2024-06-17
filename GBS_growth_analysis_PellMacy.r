###Analysis and visualization of growth and viable CFU data of Group B Streptococcus cultures across various antibiotic treatment groups
##Author: Pell, M.E.
##Date Created: 2023-09-11
##Date Modified: 2024-06-12

####Load required packages####
require(ggplot2)
require(tidyverse)
require(tidyNano)
require(tidyr)
require(rstatix)
require(viridis)
library(scales)
require(FSA)
require(asbio)

####Import and tidy data####
###This file includes 
  ## 1) all raw absorbance readings (OD600) for each biological replicate 
  ## 2) Average CFU/mL results (calculated from viable colony counts, averaged across 2-3 technical replicates) for each biological replicate 
raw_df <- read_csv("OD_CFU_raw_df.csv") %>%
  mutate_at(vars("Treatment", "Growth_timepoint(hr)", "Bio_rep"), as.factor) %>%
  mutate_at(vars("CFU/mL"), as.numeric)

##Summarize OD600 values across biological replicates for each treatment group
avg_OD600_df <- raw_df %>%
  nanolyze(name= "avg_OD600", param_var = OD600, Treatment, `Growth_timepoint(hr)`)

##remove "NAs"
tidy_CFU_df <- raw_df %>%
  drop_na()

##Summarize CFU/mL data across biological replicates for each treatment group
avg_CFU_df <- tidy_CFU_df %>%
  nanolyze(name = "avg_CFU/mL", param_var = `CFU/mL`, Treatment, `Growth_timepoint(hr)`, Treatment_timepoint)

##Join OD600 and CFU/mL summary data
all_data_summary <- left_join(avg_OD600_df, avg_CFU_df, by = c("Treatment", "Growth_timepoint(hr)"))


####Data Visualization: Figure 1####


###Figure 1A
##Plot growth curve (OD600 values over time)
growth_curve_OD600 <- all_data_summary %>%
  ggplot(aes(x = `Growth_timepoint(hr)`,
             y = `avg_OD600_mean`, group = Treatment, color = Treatment))+
  geom_errorbar(aes(ymin=avg_OD600_mean-avg_OD600_se, ymax=avg_OD600_mean+avg_OD600_se), width=.1)+
  geom_line(aes(color = as.factor(Treatment)))+
  geom_point(aes(color = as.factor(Treatment)))+
  scale_y_log10()+
  scale_color_manual(values = c("#60CEACFF", "#0B0405FF", "#395D9CFF"))+
  ylab("OD600")+
  xlab("Time (hrs)") +
  theme_bw()
growth_curve_OD600


###Figure 1B
##Plot viable CFUs over time (CFU/mL over time)
growth_curve_CFUs <- all_data_summary %>%
  ggplot(aes(x = Treatment_timepoint,
             y = `avg_CFU/mL_mean`, group = Treatment, color = Treatment))+
  geom_errorbar(aes(ymin=`avg_CFU/mL_mean`-`avg_CFU/mL_se`, ymax=`avg_CFU/mL_mean`+`avg_CFU/mL_se`), width=.1)+
  geom_line(aes(color = Treatment))+
  geom_point(aes(color = Treatment))+
  scale_y_log10()+
  scale_color_manual(values = c("#60CEACFF", "#0B0405FF", "#395D9CFF"))+
  ylab("Average CFU/mL")+
  xlab("Time post-treatment (hrs)") +
  theme_bw()
growth_curve_CFUs


####Statistical Analysis####

###Evaluate statistical differences across treatment groups at each timepoint

##Growth timepoint T0 (OD600 only)##
T0 <- raw_df %>% 
  filter(`Growth_timepoint(hr)` ==0) %>% 
  select(Treatment, `Growth_timepoint(hr)`, OD600, `CFU/mL`) %>% 
  unite("Sample", Treatment, `Growth_timepoint(hr)`, sep= "_")
T0

##Shapiro-Wilk test for normality (T0 - OD600)
shapiro.test(T0$OD600)
# W=0.89724, p-value=0.1461
# interpretation = normally distributed

##ANOVA and Tukey-HSD (T0 - OD600)
aov(`OD600`~Sample, data = T0) %>%
  tukey_hsd()
# No significant differences across treatment groups (p-values >0.05)




##Growth timepoint T1 (OD600 only)##
T1 <- raw_df %>% 
  filter(`Growth_timepoint(hr)` ==1) %>% 
  select(Treatment, `Growth_timepoint(hr)`, OD600, `CFU/mL`) %>% 
  unite("Sample", Treatment, `Growth_timepoint(hr)`, sep= "_")
T1

##Shapiro-Wilk test for normality (T1 - OD600)
shapiro.test(T1$OD600)
# W=0.90727, p-value=0.1968
# interpretation = normally distributed

##ANOVA and Tukey-HSD (T1 - OD600)
aov(`OD600`~Sample, data = T1) %>%
  tukey_hsd()
#No significant differences (p-values >0.05)





##Growth timepoint T2 (OD600 only)##
T2 <- raw_df %>% 
  filter(`Growth_timepoint(hr)` ==2) %>% 
  select(Treatment, `Growth_timepoint(hr)`, OD600, `CFU/mL`) %>% 
  unite("Sample", Treatment, `Growth_timepoint(hr)`, sep= "_")
T2


##Shapiro-Wilk test for normality (T2 - OD600)
shapiro.test(T2$OD600)
# W=0.95888, p-value=0.7677
# interpretation = normally distributed

##ANOVA and Tukey-HSD (T2 - OD600)
aov(`OD600`~Sample, data = T2) %>%
  tukey_hsd()
#No significant differences (p-values >0.05)






##Growth timepoint T3##
T3 <- raw_df %>% 
  filter(`Growth_timepoint(hr)` ==3) %>% 
  select(Treatment, `Growth_timepoint(hr)`, OD600, `CFU/mL`) %>% 
  unite("Sample", Treatment, `Growth_timepoint(hr)`, sep= "_")
T3


##Shapiro-Wilk tests for normality
#OD600#
shapiro.test(T3$OD600)
# W = 0.7484, p-value = 0.002572
# Interpretation = NOT normally distributed

#CFU/mL#
shapiro.test(T3$`CFU/mL`)
# W = 0.93823, p-value = 0.4754
# Interpretation = normally distributed

##Kruskal-Wallis and Dunn's tests (T3 - OD600)##
kruskal.test(OD600~Sample, data = T3)
# Kruskal-Wallis chi-squared = 1.2254, df = 2, p-value = 0.5419
dunn_test(OD600~Sample, data = T3)
# No significant differences (p-values >0.05)


##ANOVA and Tukey-HSD tests (T3 - CFU/mL)##
aov(`CFU/mL`~Sample, data = T3) %>%
  tukey_hsd()
# No significant differences (p-values >0.05)




##Growth timepoint T4##
T4 <- raw_df %>% 
  filter(`Growth_timepoint(hr)` ==4) %>% 
  select(Treatment, `Growth_timepoint(hr)`, OD600, `CFU/mL`) %>% 
  unite("Sample", Treatment, `Growth_timepoint(hr)`, sep= "_")
T4

##Shapiro-Wilk tests for normality
#OD600#
shapiro.test(T4$OD600)
# W = 0.83985, p-value = 0.02757
# Interpretation = NOT normally distributed

#CFU/mL#
shapiro.test(T4$`CFU/mL`)
# W = 0.96267, p-value = 0.8211
# Interpretation = normally distributed

##Kruskal-Wallis and Dunn's tests (T4 - OD600)##
kruskal.test(OD600~Sample, data = T4)
# Kruskal-Wallis chi-squared = 0.93926, df = 2, p-value = 0.6252
dunn_test(OD600~Sample, data = T4)
# No significant differences (p-values >0.05)


##ANOVA and Tukey-HSD tests (T4 - CFU/mL)##
aov(`CFU/mL`~Sample, data = T4) %>%
  tukey_hsd()
# Significant difference between Ampicillin and Control group (p=0.0281)




##Growth timepoint T5##
T5 <- raw_df %>% 
  filter(`Growth_timepoint(hr)` ==5) %>% 
  select(Treatment, `Growth_timepoint(hr)`, OD600, `CFU/mL`) %>% 
  unite("Sample", Treatment, `Growth_timepoint(hr)`, sep= "_")
T5

##Shapiro-Wilk tests for normality
#OD600#
shapiro.test(T5$OD600)
# W = 0.92955, p-value = 0.3754
# Interpretation = normally distributed

#CFU/mL#
shapiro.test(T5$`CFU/mL`)
# W = 0.79062, p-value = 0.007343
# Interpretation = NOT normally distributed


##ANOVA and Tukey-HSD tests (T5 - OD600)##
aov(OD600~Sample, data = T5) %>%
  tukey_hsd()
# Significant difference between Ampicillin and Control group (p=0.00979)

##Kruskal-Wallis and Dunn's tests (T5 - CFU/mL)##
kruskal.test(`CFU/mL`~Sample, data = T5)
# Kruskal-Wallis chi-squared = 7.7308, df = 2, p-value = 0.02095
dunn_test(`CFU/mL`~Sample, data = T5)
# Significant differences between Ampicillin-Control (p=0.008) and Erythromycin-Control (p=0.0395)



##Growth timepoint T6##
T6 <- raw_df %>% 
  filter(`Growth_timepoint(hr)` ==6) %>% 
  select(Treatment, `Growth_timepoint(hr)`, OD600, `CFU/mL`) %>% 
  unite("Sample", Treatment, `Growth_timepoint(hr)`, sep= "_")
T6

##Shapiro-Wilk tests for normality
#OD600#
shapiro.test(T6$OD600)
# W = 0.93556, p-value = 0.4427
# Interpretation = normally distributed

#CFU/mL#
shapiro.test(T6$`CFU/mL`)
# W = 0.75007, p-value = 0.002677
# Interpretation = NOT normally distributed


##ANOVA and Tukey-HSD tests (T6 - OD600)##
aov(OD600~Sample, data = T6) %>%
  tukey_hsd()
# Significant difference between Ampicillin and Control group (p=0.00259)

##Kruskal-Wallis and Dunn's tests (T6 - CFU/mL)##
kruskal.test(`CFU/mL`~Sample, data = T6)
# Kruskal-Wallis chi-squared = 7.4231, df = 2, p-value = 0.02444
dunn_test(`CFU/mL`~Sample, data = T6)
# Significant differences between Ampicillin-Control (p=0.0241) and Erythromycin-Control (p=0.0142)



##Growth timepoint T7##
T7 <- raw_df %>% 
  filter(`Growth_timepoint(hr)` ==7) %>% 
  select(Treatment, `Growth_timepoint(hr)`, OD600, `CFU/mL`) %>% 
  unite("Sample", Treatment, `Growth_timepoint(hr)`, sep= "_")
T7

##Shapiro-Wilk tests for normality
#OD600#
shapiro.test(T7$OD600)
# W = 0.92776, p-value = 0.357
# Interpretation = normally distributed

#CFU/mL#
shapiro.test(T7$`CFU/mL`)
# W = 0.75038, p-value = 0.002697
# Interpretation = NOT normally distributed


##ANOVA and Tukey-HSD tests (T7 - OD600)##
aov(OD600~Sample, data = T7) %>%
  tukey_hsd()
# Significant difference between Ampicillin-Control (p=0.000932) and Erythromycin-Control (p=0.0214)

##Kruskal-Wallis and Dunn's tests (T7 - CFU/mL)##
kruskal.test(`CFU/mL`~Sample, data = T7)
# Kruskal-Wallis chi-squared = 8.3462, df = 2, p-value = 0.0154
dunn_test(`CFU/mL`~Sample, data = T7)
# Significant differences between Ampicillin-Control (p=0.00446)




##Growth timepoint T8##
T8 <- raw_df %>% 
  filter(`Growth_timepoint(hr)` ==8) %>% 
  select(Treatment, `Growth_timepoint(hr)`, OD600, `CFU/mL`) %>% 
  unite("Sample", Treatment, `Growth_timepoint(hr)`, sep= "_")
T8

##Shapiro-Wilk tests for normality
#OD600#
shapiro.test(T8$OD600)
# W = 0.92175, p-value = 0.3007
# Interpretation = normally distributed

#CFU/mL#
shapiro.test(T8$`CFU/mL`)
# W = 0.77203, p-value = 0.004584
# Interpretation = NOT normally distributed


##ANOVA and Tukey-HSD tests (T8 - OD600)##
aov(OD600~Sample, data = T8) %>%
  tukey_hsd()
# Significant difference between Ampicillin-Control (p=0.00075) and Erythromycin-Control (p=0.0124)

##Kruskal-Wallis and Dunn's tests (T8 - CFU/mL)##
kruskal.test(`CFU/mL`~Sample, data = T8)
# Kruskal-Wallis chi-squared = 7.7308, df = 2, p-value = 0.02095
dunn_test(`CFU/mL`~Sample, data = T8)
# Significant differences between Ampicillin-Control (p=0.00811) and Erythromycin-Control (p=0.0395)


