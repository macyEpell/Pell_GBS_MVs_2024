###Analysis and visualization of proteomics data from GBS membrane vesicles across various antibiotic treatment groups.


##Date Created: 2023.07.09
##Created by: Pell, M.E.
##Date Modified: 2024-06-13


####Load required packages####
require(tidyverse)
require(tidyNano)
require(devtools)
require(factoextra)
require("viridis")
require(pairwiseAdonis)

require(hrbrthemes)
require(tm)
require(proustr)
require(VennDiagram)
require(vegan)

require(FSA)
require(dunn.test)
require(dplyr)
require(tidyr)
require(purrr)
require(stringr)
require(asbio)

require(ggbeeswarm)
require(ggsci)

require(EnhancedVolcano)
require(readr)

####Import and Tidy data####

##Import metadata with isolate and treatment information
key_df <- read_csv("Strain_Key.csv") %>%
  mutate_at(vars(Treatment, Replicate, Antibiotic), as.factor)

##Import the "Total Spectrum Counts" and respective spectrum identifications for all samples (processed and filtered in Scaffold)
Raw_df <- read_csv("proteomics_raw_2023-08-15/proteomics_raw_scaffold_output_2023-09-21.csv") %>% 
  select(`Accession Number`,`Identified Proteins (493)`,
         `AMP 1`, `AMP 2`, `AMP 3`, `AMP 4`,
         `Control 1`, `Control 2`, `Control 3`, `Control 4`,
         `ERM 1`, `ERM 2`, `ERM 3`, `ERM 4`) %>%
  gather(Sample, Spectral_Counts, `AMP 1`:`ERM 4`) %>%
  separate(Sample, into = c("Treatment", "Replicate"), sep = " ") %>% 
  mutate(Treatment = as.factor(Treatment))  %>%
  inner_join(key_df, by =  c("Treatment", "Replicate"))


##Import data key for contaminants and exclude contaminant spectra from raw data
Contaminants <-  read_csv("Contaminants.csv")

No_contaminants_Raw_df <- Raw_df %>% 
  anti_join(Contaminants) 

##Tidy Accession numbers of spectra identifications
Corrected_Accession_Raw_df <- No_contaminants_Raw_df %>% 
  select(`Accession Number`, `Identified Proteins (493)`, `Treatment`, `Replicate`, Antibiotic, `Spectral_Counts`) %>% 
  mutate(`Accession Number` = str_replace(`Accession Number`, "(?s) .*", "")) #Removes (number)


##Filter data for identified spectra that are present at least twice per treatment group
Hits_Present_Twice_Raw <- Corrected_Accession_Raw_df %>% 
  select(`Accession Number`, Treatment, Replicate, Antibiotic, Spectral_Counts) %>% 
  drop_na() %>% 
  mutate(Is_Detected = case_when(Spectral_Counts == 0 ~ "Not_Present",
                                 Spectral_Counts > 0 ~ "Present"
  )) %>% 
  filter(Is_Detected == "Present") %>% 
  group_by(`Treatment`,`Accession Number` ) %>% 
  tally() %>% 
  filter(n > 1)

##Determine the number of treatment groups that contain each accession number
Accessions_of_Present_Twice_Raw <- Hits_Present_Twice_Raw %>% 
  group_by(`Accession Number`) %>% 
  tally() # n = number of treatment groups that contain the accession number
colnames(Accessions_of_Present_Twice_Raw)[colnames(Accessions_of_Present_Twice_Raw) == 'n'] <- 'Number_of_Treatments_Occurred_In'

##Join raw data with accessions present twice to filter for spectra of interest 
#(i.e. those that are present in at least two replicates for each treatment group)
HOI_Raw <- Corrected_Accession_Raw_df %>% 
  inner_join(Accessions_of_Present_Twice_Raw)

raw.2 <- Corrected_Accession_Raw_df %>%
  inner_join(Hits_Present_Twice_Raw)

##Calculate the average spectral counts for each identified peptide across treatment groups
Averaged_by_Treatment_Raw <- raw.2 %>% 
  nanolyze(`Accession Number`, `Treatment`, param_var = Spectral_Counts, name = "Spectral_Counts") %>% 
  select(-(Spectral_Counts_sd),-(Spectral_Counts_se))



###Figure 5A: #Venn Diagram####

##Assign main dataframe to filtered spectra that were present in at least two samples
my_data <- Hits_Present_Twice_Raw

##Organize and bin data separately by treatment group: "Control" (untreated), "Ampicillin" (Amp), "Erythromycin" (Erm)
Control_1 <- my_data %>% filter(Treatment == "Control")
Control_list <- c(Control_1$`Accession Number`)
length(Control_list) #Total number of distinct peptides identified in untreated vesicles

Amp_1 <- my_data %>% filter(Treatment == "AMP")
Amp_list <- c(Amp_1$`Accession Number`)
length(Amp_list) #Total number of distinct peptides identified in ampicillin-treated vesicles

Erm_1 <- my_data %>% filter(Treatment == "ERM")
Erm_list <- c(Erm_1$`Accession Number`)
length(Erm_list) #Total number of distinct peptides identified in erythromycin-treated vesicles

###Figure 5A
##Create VennDiagram to compare quantities of distinct peptides/proteins identified in vesicles across treatment groups
venn.diagram(
  x = list(Control_list, Amp_list, Erm_list
   ),
  category.names = c("Control (330)" , "Ampicillin (276)" , "Erythromycin (319)"),
  filename = 'VennDiagram_Proteomics_2023-09-24.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 400,
  compression = "lzw",
  lwd = 1,
  col=c("#60CEACFF", "#0B0405FF", "#395D9CFF"),
  fill = c(alpha("#60CEACFF",0.3), alpha("#0B0405FF",0.3), alpha('#395D9CFF',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.27,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  cat.col = c("#60CEACFF", "#0B0405FF", "#395D9CFF"),
  rotation = 1
)

####Figure 5B:PCoA plot####

##Organize and reshape data to appropriate format for PCoA
data.les.pcoa <- HOI_Raw %>% 
  unite("Sample", Treatment, Replicate, sep = "_" ) %>% 
  select(c(`Sample`, Antibiotic, `Accession Number`, `Spectral_Counts`)) %>% 
  spread(key=`Accession Number`, value = `Spectral_Counts`)

data_cols <- c(3:419)
data.les.pcoa[ , data_cols] <- apply(data.les.pcoa[ , data_cols], 2, function(x) as.numeric(as.character(x)))

##Compute dissimilarity indices (Euclidean distance)
dis <- vegdist(data.les.pcoa[,3:419], method = "euclidean")

##Create groups based on treatment
groups <- data.les.pcoa$Antibiotic

##Perform analysis of multivariate homogeneity of group dispersions (variance)
bdisp <- betadisper(dis, groups)
##Output:
# Homogeneity of multivariate dispersions
# 
# Call: betadisper(d = dis, group = groups)
# 
# No. of Positive Eigenvalues: 11
# No. of Negative Eigenvalues: 0
# 
# Average distance to median:
#   Ampicillin Erythromycin         None 
# 45.92        51.24        64.01 
# 
# Eigenvalues for PCoA axes:
#   (Showing 8 of 11 eigenvalues)
# PCoA1 PCoA2 PCoA3 PCoA4 PCoA5 PCoA6 PCoA7 PCoA8 
# 76855 19522  9920  6570  3390  2771  2067  1793 

plot(bdisp)

##Perform principal coordinates analysis (PCoA) using dissimilarity matrix
pcoa <- cmdscale(dis)

##Plot PCoA ordination (Figure 5B)
pcoa_plot <- plot(pcoa, col = c("#60CEACFF", "#0B0405FF", "#395D9CFF")[groups], pch = c(19,19,19)[groups],
     xlab = "PCoA 1", ylab = "PCoA 2")
#Add 99% confidence ellipses to plot
pcoa_plot <- ordiellipse(pcoa, groups, display = "sites", col = c("#60CEACFF", "#0B0405FF", "#395D9CFF"),
              kind = "se", lty=c(1), conf = 0.99, alpha = 0.05, lwd = 1.5) 
#Add legend to plot
pcoa_plot<- legend(x="bottomright", legend = c("Ampicillin", "Control", "Erythromycin"),
         fill = c("#60CEACFF", "#0B0405FF", "#395D9CFF"),
         border = "black",
         cex = 0.75)

###Statistical Analysis: permanova and pairwise tests from vegan::adonis2 and pairwise.adonis package
set.seed(13)
adonis2(formula = dis ~ Antibiotic, data = data.les.pcoa, permutations = 999)
#Output:
#Significant differences in dataset (p-value=0.001)


##Pairwise statistical analysis with 999 permutations
pairwise.adonis2(dis ~ Antibiotic, data.les.pcoa, nperm = 999)
#OUTPUT:
#Amp vs Control = significantly different (p=0.038)
#Amp vs Erm = significantly different (p=0.028)
#Control vs Erm = significantly different (p=0.024)



####Analysis of differential expression of proteins (i.e. spectral abundance) across treatment groups####

#identify proteins significantly enriched in given treatment groups
HOI_stats <- HOI_Raw %>% 
  select(`Accession Number`, Antibiotic, Spectral_Counts)

HOI_stats[is.na(HOI_stats)] =  0


###Statistical analysis for differences in protein abundance across treatment groups
##Pairwise Kruskal-Wallis tests

Abx_KW <- HOI_stats %>% 
  group_by(Antibiotic) %>% 
  nest(-`Accession Number`) %>%
  mutate(test = map(data, ~pairw.kw(.$Spectral_Counts, .$Antibiotic, conf = 0.95))) %>% 
  mutate(summary = map(test, pluck, "summary")) %>% 
  mutate(summary = map(summary, mutate_if, is.factor, as.character)) %>% 
  unnest(summary)
View(Abx_KW)

#Tidy statistical output by adding identifiers for pairwise comparisons
Abx_KW_clean <- Abx_KW %>% 
  group_by(`Accession Number`) %>% 
  mutate(Identifier = str_c("ID", row_number(), sep='_'))
View(Abx_KW_clean)

##Create metadata key to identify pairwise comparisons in dataset
Identifier <- c("ID_1","ID_2","ID_3")
pairs <- c("A & E","A & C","E & C")
df <- data.frame(pairs, Identifier) %>% 
  select(Identifier, pairs)

##Join ID data frame with ST_KW_clean
Abx_post_stats <- inner_join(Abx_KW_clean, df, by = "Identifier") 

##Filter data for comparisons with a significant p-value (p<0.05)
Abx_post_stats_clean <- Abx_post_stats %>% 
  dplyr::select(`Accession Number`,Diff:`Adj. P-value`,pairs) %>% 
  filter(`Adj. P-value` <0.05)
#Add protein description information to dataframe
Stats_df_Supp <- inner_join(HOI_Raw, Abx_post_stats_clean, by = "Accession Number")  %>% 
  select(`Accession Number`, Diff, Lower, Upper, Decision, `Adj. P-value`, pairs, `Identified Proteins (493)`) %>% 
  unique()


##Summary of all pairwaise comparisons with protein description information
Stats_df_all <- inner_join(HOI_Raw, Abx_post_stats, by = "Accession Number")  %>% 
  select(`Accession Number`, Diff, Lower, Upper, Decision, `Adj. P-value`, pairs, `Identified Proteins (493)`) %>% 
  unique()

##Calculate average spectral counts for each distinct peptide by treatment group
Averaged_by_Treatment_stats <- HOI_stats %>% 
  nanolyze(`Accession Number`, `Antibiotic`, param_var = Spectral_Counts, name = "Spectral_Counts") %>% 
  select(-(Spectral_Counts_sd),-(Spectral_Counts_se), -(Spectral_Counts_N)) %>%
  pivot_wider(names_from = "Antibiotic", values_from = "Spectral_Counts_mean")

##Summary of peptides with statistically different abundance across treatment groups (pairwise) with average spectral count information by treatment group
Significant_comparisons_summary <- inner_join(Stats_df_Supp, Averaged_by_Treatment_stats, by = "Accession Number")


##Summary of all peptides with statistical output and average spectral count information by treatment group
All_comparisons_summary <- inner_join(Stats_df_all, Averaged_by_Treatment_stats, by = "Accession Number")


####Figure 6: Volcano plots####

##Format protein identity information from the "Hits Present Twice" dataframe for volcano plot
HPtwice <- Hits_Present_Twice_Raw %>%
  pivot_wider(names_from = Treatment, values_from = n)


##Filter for proteins that are shared only between Control and Ampicillin groups
A.C_shared <- HPtwice %>%
  select(-ERM) %>%
  drop_na()

##Filter for proteins that are shared only between Control and Erythromycin groups
E.C_shared <- HPtwice %>%
  select(-AMP) %>%
  drop_na()

##Filter for proteins that are shared only between Ampicillin and Erythromycin groups
A.E_shared <- HPtwice %>%
  select(-Control) %>%
  drop_na()


##Select data from pairwise Kruskal-Wallis statistical analysis summary
stats_df <- All_comparisons_summary %>%
  select("Accession Number", "pairs", "Adj. P-value", "Ampicillin", "Erythromycin", "None")
stats_df$`Adj. P-value` <- as.numeric(as.character(stats_df$`Adj. P-value`))

###Figure 6A: Ampicillin vs. Control
##calculate log2-Fold-Change of spectral counts for each identified protein across treatment groups
AvC_df <- inner_join(stats_df, A.C_shared, by="Accession Number") %>%
  select(-AMP, -Control, -Erythromycin) %>%
  filter(pairs == "A & C") %>%
  mutate(log2FoldChange = log2(Ampicillin)-log2(None))

##Construct volcano plot to visualize enriched proteins in Ampicillin vs. Control groups
AvC_VP <- EnhancedVolcano(AvC_df, lab=AvC_df$`Accession Number`, x='log2FoldChange', y='Adj. P-value', 
                          xlim=c(-7,6), ylim=c(0,3), pCutoff=0.05, FCcutoff=1.5, pointSize=2.0, labSize=3.0, 
                          title= 'Ampicillin vs. Control', cutoffLineCol = 'black', cutoffLineWidth = 0.2, hline=c(), 
                          hlineCol=c('grey0', 'grey25', 'grey50'), col=c('grey75', 'grey75', 'grey75', 'blueviolet'),
                          hlineType='longdash', hlineWidth=0.2, drawConnectors = TRUE, widthConnectors = 0.2, arrowheads = FALSE, 
                          gridlines.major=FALSE, gridlines.minor=FALSE, legendPosition = "none", max.overlaps = Inf)
AvC_VP


###Figure 6B: Erythromycin vs. Control Analysis
##calculate log2-Fold-Change of spectral counts for each identified protein across treatment groups
EvC_df <- inner_join(stats_df, E.C_shared, by="Accession Number") %>%
  select(-ERM, -Control, -Ampicillin) %>%
  filter(pairs == "E & C") %>%
  mutate(log2FoldChange = log2(Erythromycin)-log2(None))

##Construct volcano plot to visualize enriched proteins in Erytrhomycin vs. Control groups
EvC_VP <- EnhancedVolcano(EvC_df, lab=EvC_df$`Accession Number`, x='log2FoldChange', y='Adj. P-value', 
                          xlim=c(-7,6), ylim=c(0,3), pCutoff=0.05, FCcutoff=1.5, pointSize=2.0, labSize=3.0, 
                          title= 'Erythromycin vs. Control', cutoffLineCol = 'black', cutoffLineWidth = 0.2, hline=c(), 
                          hlineCol=c('grey0', 'grey25', 'grey50'), col=c('grey75', 'grey75', 'grey75', 'blueviolet'),
                          hlineType='longdash', hlineWidth=0.2, drawConnectors = TRUE, widthConnectors = 0.2, arrowheads = FALSE, 
                          gridlines.major=FALSE, gridlines.minor=FALSE, legendPosition = "none", max.overlaps = Inf)
EvC_VP


###Figure 6C: Ampicillin vs. Erythromycin Analysis
##calculate log2-Fold-Change of spectral counts for each identified protein across treatment groups
AvE_df <- inner_join(stats_df, A.E_shared, by="Accession Number") %>%
  select(-AMP, -ERM, -None) %>%
  filter(pairs == "A & E") %>%
  mutate(log2FoldChange = log2(Ampicillin)-log2(Erythromycin))

##Construct volcano plot to visualize enriched proteins in Ampicillin vs. Erythromycin groups
AvE_VP <- EnhancedVolcano(AvE_df, lab=AvE_df$`Accession Number`, x='log2FoldChange', y='Adj. P-value', 
                          xlim=c(-7,6), ylim=c(0,3), pCutoff=0.05, FCcutoff=1.5, pointSize=2.0, labSize=3.0, 
                          title= 'Ampicillin vs. Eryhtromycin', cutoffLineCol = 'black', cutoffLineWidth = 0.2, hline=c(), 
                          hlineCol=c('grey0', 'grey25', 'grey50'), col=c('grey75', 'grey75', 'grey75', 'blueviolet'),
                          hlineType='longdash', hlineWidth=0.2, drawConnectors = TRUE, widthConnectors = 0.2, arrowheads = FALSE, 
                          gridlines.major=FALSE, gridlines.minor=FALSE, legendPosition = "none", max.overlaps = Inf)

AvE_VP

####Figure 7 and S2: Plot spectral abundance for proteins of interest across treatment groups####

####Figure 7

##Figure 7A: Plot penicillin binding proteins
pbps_list <- c("Q8E1Q5", "Q8E240", "Q8DWZ3", "Q8E0G8", "Q8E1R6", "Q8E153")
pbps_df <-  filter(HOI_Raw, `Accession Number` %in% pbps_list)
pbps_df

pbps_plot <- pbps_df %>% 
  ggplot(aes(x = `Accession Number`,
             y = `Spectral_Counts`, group = Treatment, color = Treatment))+
  geom_boxplot(fill = NA, position = "dodge", outlier.colour = "black")+
  geom_point(position=position_jitterdodge(),aes(group=Treatment), alpha = 0.5)+
  facet_wrap(vars(`Accession Number`), scales = "free_x", nrow = 1)+
  theme_bw()+
  theme(strip.text.x = element_blank(),
        legend.position = "top")+
  scale_x_discrete(labels=c("Q8E1Q5"="Pbp1A", "Q8E240"="Pbp1B", "Q8DWZ3"="Pbp2A", 
                            "Q8E0G8"="Pbp2B", "Q8E1R6"="Pbp2X", "Q8E153"="FibA"))+
  ylab("Spectral Counts")+
  xlab("")+
  scale_color_manual(values = c("#60CEACFF", "#0B0405FF", "#395D9CFF"),
                     name = "", labels=c("Ampicillin", "Control", "Erythromycin"))+
  theme(text = element_text(size = 16))

pbps_plot


##Figure 7B: 50S ribosomal proteins
rib_list <- c("Q8E2B5", "Q8E2C8", "Q8E2C9", "Q8E2D1", "Q8E2B9")
rib_df <-  filter(HOI_Raw, `Accession Number` %in% rib_list)

rib_plot <- rib_df %>% 
  ggplot(aes(x = `Accession Number`,
             y = `Spectral_Counts`, group = Treatment, color = Treatment))+
  geom_boxplot(fill = NA, position = "dodge", outlier.colour = "black")+
  geom_point(position=position_jitterdodge(),aes(group=Treatment), alpha = 0.5)+
  facet_wrap(vars(`Accession Number`), scales = "free_x", nrow = 1)+
  theme_bw()+
  theme(strip.text.x = element_blank(),
        legend.position = "top")+
  scale_x_discrete(labels=c("Q8E2B5"="RplO", "Q8E2C8"="RplB", "Q8E2C9"="RplW",
                            "Q8E2D1"="RplC", "Q8E2B9"="RplF"))+
  ylab("Spectral Counts")+
  xlab("")+
  scale_color_manual(values = c("#60CEACFF", "#0B0405FF", "#395D9CFF"),
                     name = "", labels=c("Ampicillin", "Control", "Erythromycin"))+
  theme(text = element_text(size = 16))

rib_plot



###Figure S2

##Figure S2A: Virulence proteins
vir_list <- c("Q8CX01", "Q8DWP1", "Q9AFI1")
vir_df <-  filter(HOI_Raw, `Accession Number` %in% vir_list)

vir_plot <- vir_df %>% 
  ggplot(aes(x = `Accession Number`,
             y = `Spectral_Counts`, group = Treatment, color = Treatment))+
  geom_boxplot(fill = NA, position = "dodge", outlier.colour = "black")+
  geom_point(position=position_jitterdodge(),aes(group=Treatment), alpha = 0.5)+
  facet_wrap(vars(`Accession Number`), scales = "free_x", nrow = 1)+
  theme_bw()+
  theme(strip.text.x = element_blank(),
        legend.position = "top")+
  scale_x_discrete(labels=c("Q8CX01"="cAMP \nfactor", "Q8DWP1"="Serine \nprotease", "Q9AFI1" = "CpsD"))+
  ylab("Spectral Counts")+
  xlab("")+
  scale_color_manual(values = c("#60CEACFF", "#0B0405FF", "#395D9CFF"),
                     name = "", labels=c("Ampicillin", "Control", "Erythromycin"))+
  theme(text = element_text(size = 16))

vir_plot


##Figure S2B: Stress-associated proteins
stress_list <- c("P0A3J3", "Q8DZG4", "Q8DZT2", "Q8E0J3")
stress_df <-  filter(HOI_Raw, `Accession Number` %in% stress_list)

stress_plot <- stress_df %>% 
  ggplot(aes(x = `Accession Number`,
             y = `Spectral_Counts`, group = Treatment, color = Treatment))+
  geom_boxplot(fill = NA, position = "dodge", outlier.colour = "black")+
  geom_point(position=position_jitterdodge(),aes(group=Treatment), alpha = 0.5)+
  facet_wrap(vars(`Accession Number`), scales = "free_x", nrow = 1)+
  theme_bw()+
  theme(strip.text.x = element_blank(),
        legend.position = "top")+
  scale_x_discrete(labels=c("P0A3J3"="DnaK","Q8DZG4"="Gls24", 
                            "Q8DZT2"="CtsA", "Q8E0J3"= "General \nstress protein" ))+
  ylab("Spectral Counts")+
  xlab("")+
  scale_color_manual(values = c("#60CEACFF", "#0B0405FF", "#395D9CFF"),
                     name = "", labels=c("Ampicillin", "Control", "Erythromycin"))+
  theme(text = element_text(size = 16))

stress_plot

##Figure S2C: cell wall-associated proteins
wall_list <- c("Q8E2H1", "Q8DY79", "Q8E1I6")
wall_df <-  filter(HOI_Raw, `Accession Number` %in% cell_list)

wall_plot <- wall_df %>% 
  ggplot(aes(x = `Accession Number`,
             y = `Spectral_Counts`, group = Treatment, color = Treatment))+
  geom_boxplot(fill = NA, position = "dodge", outlier.colour = "black")+
  geom_point(position=position_jitterdodge(),aes(group=Treatment), alpha = 0.5)+
  facet_wrap(vars(`Accession Number`), scales = "free_x", nrow = 1)+
  theme_bw()+
  theme(strip.text.x = element_blank(),
        legend.position = "top")+
  scale_x_discrete(labels=c("Q8E2H1"="PcsB", "Q8DY79"="MltG", "Q8E1I6" = "LytR-cpsA-psr"))+
  ylab("Spectral Counts")+
  xlab("")+
  scale_color_manual(values = c("#60CEACFF", "#0B0405FF", "#395D9CFF"),
                     name = "", labels=c("Ampicillin", "Control", "Erythromycin"))+
  theme(text = element_text(size = 16))

wall_plot


##Figure S2D: cell division-associated proteins
cell_list <- c("Q8E144", "Q8E178", "Q8E1R7", "Q8E143", "Q8CX05")
cell_df <-  filter(HOI_Raw, `Accession Number` %in% cell_list)

cell_plot <- cell_df %>%
  ggplot(aes(x = `Accession Number`,
             y = `Spectral_Counts`, group = Treatment, color = Treatment))+
  geom_boxplot(fill = NA, position = "dodge", outlier.colour = "black")+
  geom_point(position=position_jitterdodge(),aes(group=Treatment), alpha = 0.5)+
  facet_wrap(vars(`Accession Number`), scales = "free_x", nrow = 1)+
  theme_bw()+
  theme(strip.text.x = element_blank(),
        legend.position = "top")+
  scale_x_discrete(labels=c("Q8E144"="FtsE", "Q8E178"="DivIVA", 
                            "Q8E1R7"="FtsL", "Q8E143"="FtsX", "Q8CX05"="FtsK"))+
  ylab("Spectral Counts")+
  xlab("")+
  scale_color_manual(values = c("#60CEACFF", "#0B0405FF", "#395D9CFF"),
                     name = "", labels=c("Ampicillin", "Control", "Erythromycin"))+
  theme(text = element_text(size = 16))

cell_plot










