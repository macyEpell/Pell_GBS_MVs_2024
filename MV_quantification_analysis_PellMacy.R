###Analysis and visualization of nanoparticle tracking data to quantify and compare production of Group B Streptococcus membrane vesicles
#across antibiotic treatment groups

###Author: Pell, M.E.
###Date Created: 2023.06.21
###Date Modified: 2024.06.11

####Load packages####
require(tidyverse)
# devtools::install_github("nguyens7/tidyNano", force = TRUE)
require(tidyNano)
require(ggrepel)
require(ggplot2)
require(viridis)
require(rstatix)
require(scales)

####Data Import and Tidying####
#read in the key csv file (CFU and volume information for each sample)
key_df <- read.csv("112_amp_erm_4hr/112_amp-erm_key.csv") %>% 
  mutate_at(vars(Treatment), as.factor) %>% 
  mutate_at(vars(Bio_Rep), as.factor) %>%
  mutate_at(vars(Strain), as.factor)
key_df

#Read in all raw data output from the NanoSight in .csv format
#The raw .csv output for each replicate experiment is organized into one directory called "experiment_summaries"
nanosight_master_df <- nanocombine("112_amp_erm_4hr/experiment_summaries/")
nanosight_master_df

#tidy raw data by separating the components for each sample name into columns  
tidy_df <- nanosight_master_df %>% nanotidy(sep_var = c("Strain","Treatment", "Dilution", "Bio_Rep", "Tech_Rep"))
View(tidy_df)

#join the tidy data output and the sample key 
tidy_df_full <- tidy_df  %>% 
  inner_join(key_df, by = c("Strain", "Treatment", "Bio_Rep"))
View(tidy_df_full)


####Data analysis and Normalization####

#Summarize particle counts across technical replicates
#returns an average count for each particle size across the 5 technical replicates per biological replicate 
summary.tech_reps <- tidy_df_full %>%
  nanolyze(particle_size, Treatment, Bio_Rep, Prep_Vol, Res_Vol, 
           param_var = True_count,
           name = "avg_MV_count")
View(summary.tech_reps)

##Normalize counts by the resuspension volume for each sample (recorded in sample key dataframe)
#Note: The raw NanoSight output is in "particles per mL"
normalized_counts <- summary.tech_reps %>% 
  mutate(particles_per_uL = avg_MV_count_mean/1000) %>%  #convert to uL (nanosight output is in particles/mL)
  mutate(avg_total_MVs = particles_per_uL*Res_Vol) %>%  #normalize to resuspension volume (uL) to give total MVs isolated from culture
  mutate(avg_MVs_per_mL = avg_total_MVs/Prep_Vol)       #divide total MVs by prep volume to give normalized MVs/mL of culture volume
View(normalized_counts)


##Calculate summary of counts across "bio_reps" (biological replicates) for each treatment from normalized data
summary.bio_rep.total_normalized <- normalized_counts %>%
  nanolyze(particle_size, Treatment, param_var = avg_total_MVs, name="total_MVs_bio_rep")
View(summary.bio_rep.total_normalized)


####Figure 3A####
#Plot the size distribution of average total counts across biological replicates
#Note: geom_ribbon (grey) displays standard error
PLOT_avg_bio_rep_distribution <- summary.bio_rep.total_normalized %>%
  ggplot(aes(x=particle_size, y=total_MVs_bio_rep_mean, group = Treatment, color = Treatment)) +
  geom_ribbon(aes(ymin=total_MVs_bio_rep_mean-total_MVs_bio_rep_se, ymax= total_MVs_bio_rep_mean+total_MVs_bio_rep_se), 
              fill="grey70", color = "grey") +
  geom_line(linewidth = 1) +
  theme_bw() +
  scale_x_continuous(limits = c(0,1000)) +
  scale_color_manual(values = c("#60CEACFF", "#0B0405FF", "#395D9CFF" ),
                     name = "Treatment", labels=c("Ampicillin", "Control", "Erythromycin")) +
  labs(x = "\nParticle Size (nm)",
       y = "Average Particle Count\n") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
PLOT_avg_bio_rep_distribution



####Figure 3B####

##Calculate the sum of the average total-particle-counts across all particle sizes for each biological replicate
sum_avg_total_mvs <- normalized_counts %>%
  nanocount(Treatment, Bio_Rep, 
            param_var = avg_total_MVs, 
            name = "Total_all_MVs")
View(sum_avg_total_mvs)

##Annotations of statistical significance for Boxplot (**See "Statistical Analysis" Section below for details)
MVperPrep_anno <- data.frame(x1 = c(1, 2, 1), 
                        x2 = c(2, 3, 3), 
                        y1 = c(5.1e+10, 4e+10, 5.5e+10), 
                        y2 = c(5.2e+10, 4.1e+10, 5.6e+10), 
                        xstar = c(1.5, 2.5, 2), 
                        ystar = c(5.25e+10, 4.15e+10, 5.65e+10),
                        lab = c("***", "*","**")
                        )
MVperPrep_anno

###Figure 3B
##Boxplot of total membrane vesicles (MVs) across treatment groups, with 6 biological replicates each
MVperPrep_plot <- sum_avg_total_mvs %>% 
  ggplot(aes(x=Treatment, y = Total_all_MVs_count, group = Treatment, color = Treatment))+
  geom_boxplot(fill = NA, outlier.color = "black") +
  geom_jitter(width = 0.1, alpha=0.5, size = 2) +
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_text(inherit.aes = FALSE, data = MVperPrep_anno, aes(x = xstar,  y = ystar, label = lab), size=6) +
  scale_color_manual(values = c("#60CEACFF", "#0B0405FF", "#395D9CFF"),
                     name = "Treatment", labels=c("Ampicillin", "Control", "Erythromycin")) +
  geom_segment(inherit.aes = FALSE, data = MVperPrep_anno, aes(x = x1, xend = x1, 
                                     y = y1, yend = y2),
               colour = "black") +
  geom_segment(inherit.aes = FALSE, data = MVperPrep_anno, aes(x = x2, xend = x2, 
                                     y = y1, yend = y2),
               colour = "black") +
  geom_segment(inherit.aes = FALSE, data = MVperPrep_anno, aes(x = x1, xend = x2, 
                                     y = y2, yend = y2),
               colour = "black")+
  labs(y = "Total Number of MVs\n",
       x = "Treatment Group\n")

MVperPrep_plot



####Statistical Analysis####


###Shapiro-Wilk normality test

shapiro.test(sum_avg_total_mvs$Total_all_MVs_count)
## W = 0.90567, p-value = 0.07218
## interpretation = data is normally distributed



###Anova with Tukey-HSD

aov(Total_all_MVs_count~Treatment, data = sum_avg_total_mvs) %>%
  tukey_hsd()
# term      group1     group2       null.value      estimate      conf.low     conf.high      p.adj p.adj.signif
# * <chr>     <chr>      <chr>             <dbl>         <dbl>         <dbl>         <dbl>      <dbl> <chr>       
#   1 Treatment ampicillin control               0 -31545984123. -42257691247. -20834276999. 0.00000423 ****        
#   2 Treatment ampicillin erythromycin          0 -17687853486. -28399560610.  -6976146362. 0.00176    **          
#   3 Treatment control    erythromycin          0  13858130637.   3146423513.  24569837761. 0.0113     * 
## Interpretation: Total MV production for the following comparisons are significantly different; 
##  1) Ampicillin vs. Control (p=0.00000423)
##  2) Ampicillin vs. Erythromycin (p=0.00176)
##  3) Control vs. Erythromycin (p=0.0113)
