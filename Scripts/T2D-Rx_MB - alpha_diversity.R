# T2D-Rx_MB_analysis - effect of T2D-Rx towards microbiome - alpha diversity analysis

# Load packages
library(here)
library(dplyr)
library(flextable)
library(tidyverse)
library(tibble)
library(stringr)
library(ggplot2)
library(mia)
library(miaViz)
library(purrr)
library(tidyr)
library(vegan)
library(pheatmap)
library(ggcorrplot)
library(ggpubr)
library(paletteer)
library(ggthemes)
library(rstatix)

# ___________________________________________________________________________ #

# Load input files

## Load counts (otu-matrix)
counts <- read.csv("otu_matrix.csv", sep = ";", row.names = 1, check.names = FALSE)

## Load taxonomy (rowData)
tax <- read.csv("taxonomy.csv", sep = ";", row.names = 1, check.names = FALSE)

## Load metadata (colData)
samples <- read.csv("metadata.csv", row.names = 1, sep = ",", check.names = FALSE)

## Order columns in counts based on rows in samples
counts <- counts[ , rownames(samples)]

## Convert input files into right format
counts <- as.matrix(counts)

## Make a TreeSummarizedExperiment
tse_taxa <- TreeSummarizedExperiment(assays =  SimpleList(counts = counts),
                                     colData = DataFrame(samples),
                                     rowData = DataFrame(tax))

## Remove features containing human DNA
tse_prelim <- tse_taxa[!rowData(tse_taxa)$Phylum == "" & 
                         !rowData(tse_taxa)$Class == "" & 
                         !rowData(tse_taxa)$Kingdom == "Unassigned", ]

## Exclude samples
tse <- tse_prelim[ , !colnames(tse_prelim) %in% c("GLP1RA-5-2", "GLP1RA-5-3", 
                                                  "GLP1RA-5-4", "GLP1RA-6-2", 
                                                  "NEG-Valida", "GLP1RA-15-2-2", 
                                                  "Elini-proov", "AL")]


# ___________________________________________________________________________ #


# Convert counts into relative abundances
tse <- transformCounts(tse, assay.type = "counts", method = "relabundance")

# Convert relative abundances into CLR-transformed values
mat <- assay(tse, "relabundance")
tse <- transformCounts(tse, assay.type = "relabundance", method = "clr", 
                       pseudocount = min(mat[mat>0]))

# Collapse into Genus level
tse_genus <- agglomerateByRank(tse, rank = "Genus")

# Remove "Genus:" label
rownames(tse_genus) <- sub("Genus:", "", rownames(tse_genus))

# Separate by medication
tse_glp <- tse_genus[ , colData(tse_genus)$Medication == "GLP-1-RA"]
tse_sglt <- tse_genus[ , colData(tse_genus)$Medication == "SGLT-2"]

# ___________________________________________________________________________ #

# Alpha diversity analysis with repeated measures ANOVA and t-tests

## GLP-1-RA

### Calculate alpha diversity
tse_glp <- estimateDiversity(tse_glp, assay.type = "relabundance", index = "shannon")
tse_glp <- estimateEvenness(tse_glp, assay.type = "relabundance", index = "pielou")
tse_glp <- estimateRichness(tse_glp, assay.type = "relabundance", index = "observed")

## Prepare data

### Extract CLR-transformed abundance values
tse_glp_clr <- as.data.frame(assay(tse_glp, "clr"))

### Extract metadata with alpha diversity indices
tse_glp_samples <- as.data.frame(colData(tse_glp))

### Perform PCA 
pca_glp <- prcomp(t(tse_glp_clr))

### Combine PCs with metadata and alpha diversity
glp_alpha_beta_raw <- data.frame(pca_glp$x[ , 1:5],
                                 tse_glp_samples)

## Prepare a for-loop for repeated measures ANOVA

### Create a list for alpha diversity indices
alpha_div <- c("shannon", "pielou", "observed")

### Create a tibble for results
glp_alpha_anova_results <- tibble(x = 1:3) %>% 
  add_column(alpha_div, .before = "x") %>% 
  dplyr::rename(p_value = x)

# Execute the for-loop

for (j in 1:length(alpha_div)){
  
  glp_temp_alpha <- glp_alpha_beta_raw %>% 
    rownames_to_column(var = "SampleID") %>% 
    select(c(7:8, 72:74)) %>% 
    pivot_longer(cols = 3:5, names_to = "alpha", values_to = "value") %>% 
    filter(alpha == alpha_div[j])
  
  glp_alpha_test <- glp_temp_alpha %>% 
    group_by(PatientID) %>% 
    filter(n() != 1) %>% 
    arrange(Timepoint, PatientID) %>% 
    ungroup()
    
    model_glp <- rstatix::anova_test(data = glp_alpha_test, dv = value, wid = PatientID, within = Timepoint)
    
    glp_alpha_anova_results$p_value[j] <- model_glp$ANOVA$p
    
}

## Prepare a for-loop for t-tests

### Create a list for timepoints
timepoints <- c("II", "III", "IV")

### Create a list for alpha diversity indices
alpha_div <- c("shannon", "pielou", "observed")

### Make a tibble for estimates
glp_alpha_ttest_estimates <- tibble(x = 1:3, y = 1:3, z = 1:3) %>% 
  add_column(alpha_div, .before = "x") %>% 
  dplyr::rename(II = x,
                III = y,
                IV = z) %>% 
  column_to_rownames(var = "alpha_div")

### Make a tibble for p-values
glp_alpha_ttest_pvalues <- tibble(x = 1:3, y = 1:3, z = 1:3) %>% 
  add_column(alpha_div, .before = "x") %>% 
  dplyr::rename(II = x,
                III = y,
                IV = z) %>% 
  column_to_rownames(var = "alpha_div")

### Execute the for-loop
for (j in 1:length(alpha_div)){
  
  temp_alpha_glp <- glp_alpha_beta_raw %>% 
    rownames_to_column(var = "SampleID") %>% 
    select(c(7:8, 72:74)) %>% 
    pivot_longer(cols = 3:5, names_to = "alpha_div", values_to = "value") %>% 
    filter(alpha_div == alpha_div[j])
  
  for (i in 1:length(timepoints)){
    alpha_test_glp <- temp_alpha_glp %>% 
      filter(Timepoint %in% c("I", timepoints[i])) %>% 
      group_by(PatientID) %>% 
      filter(n() != 1) %>% 
      arrange(Timepoint, PatientID) %>% 
      ungroup()
    
    alpha_test_glp_raw <- t.test(value ~ Timepoint, data = alpha_test_glp, paired = TRUE) 
    
    alpha_test_glp <- alpha_test_glp_raw %>% 
      broom::tidy()
    
    glp_alpha_ttest_estimates[alpha_div[j], timepoints[i]] <- alpha_test_glp$estimate
    glp_alpha_ttest_pvalues[alpha_div[j], timepoints[i]] <- alpha_test_glp$p.value
    
    glp_alpha_ttest_results <- cbind(glp_alpha_ttest_estimates, glp_alpha_ttest_pvalues)
    
  }
}

# ___________________________________________________________________________ #

## SGLT-2

### Calculate alpha diversity
tse_sglt <- estimateDiversity(tse_sglt, assay.type = "relabundance", index = "shannon")
tse_sglt <- estimateEvenness(tse_sglt, assay.type = "relabundance", index = "pielou")
tse_sglt <- estimateRichness(tse_sglt, assay.type = "relabundance", index = "observed")

## Prepare data

### Extract CLR-transformed abundance values
tse_sglt_clr <- as.data.frame(assay(tse_sglt, "clr"))

### Extract metadata with alpha diversity indices
tse_sglt_samples <- as.data.frame(colData(tse_sglt))

### Perform PCA 
pca_sglt <- prcomp(t(tse_sglt_clr))

### Combine PCs with metadata and alpha diversity
sglt_alpha_beta_raw <- data.frame(pca_sglt$x[ , 1:5],
                                  tse_sglt_samples)

## Prepare a for-loop for repeated measures ANOVA

### Create a list for timepoints
timepoints <- c("II", "III", "IV")

### Create a list for alpha diversity indices
alpha_div <- c("shannon", "pielou", "observed")

### Create a tibble for results
sglt_alpha_anova_results <- tibble(x = 1:3) %>% 
  add_column(alpha_div, .before = "x") %>% 
  dplyr::rename(p_value = x)

# Execute the for-loop

for (j in 1:length(alpha_div)){
  
  sglt_temp_alpha <- sglt_alpha_beta_raw %>% 
    rownames_to_column(var = "SampleID") %>% 
    select(c(7:8, 72:74)) %>% 
    pivot_longer(cols = 3:5, names_to = "alpha", values_to = "value") %>% 
    filter(alpha == alpha_div[j])
  
  sglt_alpha_test <- sglt_temp_alpha %>% 
    group_by(PatientID) %>% 
    filter(n() != 1) %>% 
    arrange(Timepoint, PatientID) %>% 
    ungroup()
  
  model_sglt <- rstatix::anova_test(data = sglt_alpha_test, dv = value, wid = PatientID, within = Timepoint)
  
  sglt_alpha_anova_results$p_value[j] <- model_sglt$ANOVA$p
  
}

## Prepare a for-loop for t-tests

### Create a list for timepoints
timepoints <- c("II", "III", "IV")

### Create a list for alpha diversity indices
alpha_div <- c("shannon", "pielou", "observed")

### Make a tibble for estimates
sglt_alpha_ttest_estimates <- tibble(x = 1:3, y = 1:3, z = 1:3) %>% 
  add_column(alpha_div, .before = "x") %>% 
  dplyr::rename(II = x,
                III = y,
                IV = z) %>% 
  column_to_rownames(var = "alpha_div")

### Make a tibble for p-values
sglt_alpha_ttest_pvalues <- tibble(x = 1:3, y = 1:3, z = 1:3) %>% 
  add_column(alpha_div, .before = "x") %>% 
  dplyr::rename(II = x,
                III = y,
                IV = z) %>% 
  column_to_rownames(var = "alpha_div")

### Execute the for-loop
for (j in 1:length(alpha_div)){
  
  temp_alpha_sglt <- sglt_alpha_beta_raw %>% 
    rownames_to_column(var = "SampleID") %>% 
    select(c(7:8, 72:74)) %>% 
    pivot_longer(cols = 3:5, names_to = "alpha_div", values_to = "value") %>% 
    filter(alpha_div == alpha_div[j])
  
  for (i in 1:length(timepoints)){
    alpha_test_sglt <- temp_alpha_sglt %>% 
      filter(Timepoint %in% c("I", timepoints[i])) %>% 
      group_by(PatientID) %>% 
      filter(n() != 1) %>% 
      arrange(Timepoint, PatientID) %>% 
      ungroup()
    
    alpha_test_sglt_raw <- t.test(value ~ Timepoint, data = alpha_test_sglt, paired = TRUE) 
    
    alpha_test_sglt <- alpha_test_sglt_raw %>% 
      broom::tidy()
    
    sglt_alpha_ttest_estimates[alpha_div[j], timepoints[i]] <- alpha_test_sglt$estimate
    sglt_alpha_ttest_pvalues[alpha_div[j], timepoints[i]] <- alpha_test_sglt$p.value
    
    sglt_alpha_ttest_results <- cbind(sglt_alpha_ttest_estimates, sglt_alpha_ttest_pvalues)
    
  }
}

# ADD P-VALUE CORRECTION!

# ___________________________________________________________________________ #

# Visualize the results

## GLP-1-RA

glp_alpha_plot_data <- glp_alpha_beta_raw %>% 
  rownames_to_column(var = "SampleID") %>% 
  select(c(7:9, 72:74)) %>% 
  pivot_longer(cols = 4:6, names_to = "alpha_div", values_to = "value") %>% 
  pivot_wider(names_from = Timepoint, values_from = value) %>% 
  mutate(I_II = I - II, 
         I_III = I - III,
         I_IV = I - IV) %>% 
  pivot_longer(cols = 4:7, names_to = "Timepoint", values_to = "Value") %>% 
  pivot_longer(cols = 4:6, names_to = "chg", values_to = "chg_value")


glp_alpha_shannon <- glp_alpha_plot_data %>% 
  filter(alpha_div == "shannon") %>% 
  na.omit() %>% 
  group_by(chg) %>% 
  ggplot(aes(x = chg, y = chg_value)) +
  facet_grid(Medication~alpha_div, switch = "y") +
  geom_boxplot() +
  geom_point() +
  ggpubr::stat_compare_means(comparisons = my_comparisons, method = "t.test") +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.3) +
  geom_hline(aes(yintercept = 0, colour = "red"), show.legend = F, linetype = "dashed") +
  labs(x = NULL,
       y = "Change in value")

glp_alpha_pielou <- glp_alpha_plot_data %>% 
  filter(alpha_div == "pielou") %>% 
  na.omit() %>% 
  group_by(chg) %>% 
  ggplot(aes(x = chg, y = chg_value)) +
  facet_grid(Medication~alpha_div, switch = "y") +
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.3) +
  geom_hline(aes(yintercept = 0, colour = "red"), show.legend = F, linetype = "dashed") +
  labs(x = "Change in timepoint compared to BL",
       y = NULL)  

glp_alpha_observed <- glp_alpha_plot_data %>% 
  filter(alpha_div == "observed") %>% 
  na.omit() %>% 
  group_by(chg) %>% 
  ggplot(aes(x = chg, y = chg_value)) +
  facet_grid(Medication~alpha_div, switch = "y") +
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.3) +
  geom_hline(aes(yintercept = 0, colour = "red"), show.legend = F, linetype = "dashed") +
  labs(x = NULL,
       y = NULL)

## SGLT-2

sglt_alpha_plot_data <- sglt_alpha_beta_raw %>% 
  rownames_to_column(var = "SampleID") %>% 
  select(c(7:9, 72:74)) %>% 
  pivot_longer(cols = 4:6, names_to = "alpha_div", values_to = "value") %>% 
  pivot_wider(names_from = Timepoint, values_from = value) %>% 
  mutate(I_II = I - II, 
         I_III = I - III,
         I_IV = I - IV) %>% 
  pivot_longer(cols = 4:7, names_to = "Timepoint", values_to = "Value") %>% 
  pivot_longer(cols = 4:6, names_to = "chg", values_to = "chg_value")

sglt_alpha_shannon <- sglt_alpha_plot_data %>% 
  filter(alpha_div == "shannon") %>% 
  na.omit() %>% 
  group_by(chg) %>% 
  ggplot(aes(x = chg, y = chg_value)) +
  facet_grid(Medication~alpha_div, switch = "y") +
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.3) +
  geom_hline(aes(yintercept = 0, colour = "red"), show.legend = F, linetype = "dashed") +
  labs(x = NULL,
       y = "Change in value")    

sglt_alpha_pielou <- sglt_alpha_plot_data %>% 
  filter(alpha_div == "pielou") %>% 
  na.omit() %>% 
  group_by(chg) %>% 
  ggplot(aes(x = chg, y = chg_value)) +
  facet_grid(Medication~alpha_div, switch = "y") +
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.3) +
  geom_hline(aes(yintercept = 0, colour = "red"), show.legend = F, linetype = "dashed") +
  labs(x = "Change in timepoint compared to BL",
       y = NULL)

sglt_alpha_observed <- sglt_alpha_plot_data %>% 
  filter(alpha_div == "observed") %>% 
  na.omit() %>% 
  group_by(chg) %>% 
  ggplot(aes(x = chg, y = chg_value)) +
  facet_grid(Medication~alpha_div, switch = "y") +
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.3) +
  geom_hline(aes(yintercept = 0, colour = "red"), show.legend = F, linetype = "dashed") +
  labs(x = NULL,
       y = NULL)

common_alpha <- ggpubr::ggarrange(glp_alpha_shannon, glp_alpha_pielou, glp_alpha_observed, 
                                  sglt_alpha_shannon, sglt_alpha_pielou, sglt_alpha_observed, 
                                  ncol = 3, nrow = 2, common.legend = TRUE, legend = "right")
