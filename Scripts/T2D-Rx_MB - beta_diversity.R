# T2D-Rx_MB_analysis - effect of T2D-Rx towards microbiome - beta diversity analysis

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

# Beta diversity analysis with repeated measures ANOVA and t-tests

## GLP-1-RA

### Extract CLR-transformed abundance values
tse_glp_clr <- as.data.frame(assay(tse_glp, "clr"))

### Extract metadata
tse_glp_samples <- as.data.frame(colData(tse_glp))

### Perform PCA 
pca_glp <- prcomp(t(tse_glp_clr))

### Combine PCs with metadata
glp_beta_raw <- data.frame(pca_glp$x[ , 1:5],
                           tse_glp_samples)


## Prepare a for-loop for repeated measures ANOVA

### Create a list for beta diversity indices
beta_div <- c("PC1", "PC2", "PC3", "PC4", "PC5")

### Create a tibble for results
glp_beta_anova_results <- tibble(x = 1:5) %>% 
  add_column(beta_div, .before = "x") %>% 
  dplyr::rename(p_value = x)

# Execute the for-loop

for (j in 1:length(beta_div)){
  
  glp_temp_beta <- glp_beta_raw %>% 
    rownames_to_column(var = "SampleID") %>% 
    select(c(2:8)) %>% 
    pivot_longer(cols = 1:5, names_to = "beta", values_to = "value") %>% 
    filter(beta == beta_div[j])
  
  glp_beta_test <- glp_temp_beta %>% 
    group_by(PatientID) %>% 
    filter(n() != 1) %>% 
    arrange(Timepoint, PatientID) %>% 
    ungroup()
  
  model_glp <- rstatix::anova_test(data = glp_beta_test, dv = value, wid = PatientID, within = Timepoint)
  
  glp_beta_anova_results$p_value[j] <- model_glp$ANOVA$p
  
}

### Correct p-values for multiple testing w/ Benjamini-Hochberg method
glp_anova_beta_results <- glp_beta_anova_results %>% 
  mutate(p_value_BH = p.adjust(p_value, method = "BH"))

## Prepare a for-loop for t-tests

### Create a list for timepoints
timepoints <- c("II", "III", "IV")

### Create a list for beta diversity indices
beta_div <- c("PC1", "PC2", "PC3", "PC4", "PC5")

### Make a tibble for estimates
glp_beta_ttest_estimates <- tibble(x = 1:5, y = 1:5, z = 1:5) %>% 
  add_column(beta_div, .before = "x") %>% 
  dplyr::rename(II = x,
                III = y,
                IV = z) %>% 
  column_to_rownames(var = "beta_div")

### Make a tibble for p-values
glp_beta_ttest_pvalues <- tibble(x = 1:5, y = 1:5, z = 1:5) %>% 
  add_column(beta_div, .before = "x") %>% 
  dplyr::rename(II = x,
                III = y,
                IV = z) %>% 
  column_to_rownames(var = "beta_div")

### Execute the for-loop
for (j in 1:length(beta_div)){
  
  temp_beta_glp <- glp_beta_raw %>% 
    rownames_to_column(var = "SampleID") %>% 
    select(c(2:8)) %>% 
    pivot_longer(cols = 1:5, names_to = "beta_div", values_to = "value") %>% 
    filter(beta_div == beta_div[j])
  
  for (i in 1:length(timepoints)){
    beta_test_glp <- temp_beta_glp %>% 
      filter(Timepoint %in% c("I", timepoints[i])) %>% 
      group_by(PatientID) %>% 
      filter(n() != 1) %>% 
      arrange(Timepoint, PatientID) %>% 
      ungroup()
    
    beta_test_glp_raw <- t.test(value ~ Timepoint, data = beta_test_glp, paired = TRUE) 
    
    beta_results_glp <- beta_test_glp_raw %>% 
      broom::tidy()
    
    glp_beta_ttest_estimates[beta_div[j], timepoints[i]] <- beta_results_glp$estimate
    glp_beta_ttest_pvalues[beta_div[j], timepoints[i]] <- beta_results_glp$p.value
    
    glp_beta_estimates <- glp_beta_ttest_estimates %>% 
      rownames_to_column(var = "Beta_div") %>% 
      pivot_longer(cols = 2:4, names_to = "Timepoint", values_to = "Estimate")
    
    glp_beta_pvalues <- glp_beta_ttest_pvalues %>% 
      rownames_to_column(var = "Beta_div") %>% 
      pivot_longer(cols = 2:4, names_to = "Timepoint", values_to = "p_value")
    
    glp_beta_ttest_results <- glp_beta_estimates %>% 
      mutate(p_value = glp_beta_pvalues$p_value)
    
  }
}

### Correct p-values for multiple testing w/ Benjamini-Hochberg method
glp_beta_ttest_BH <- glp_beta_ttest_results %>% 
  mutate(p_value_BH = p.adjust(p_value, method = "BH"))

# ___________________________________________________________________________ #


## SGLT-2

### Extract CLR-transformed abundance values
tse_sglt_clr <- as.data.frame(assay(tse_sglt, "clr"))

### Extract metadata
tse_sglt_samples <- as.data.frame(colData(tse_sglt))

### Perform PCA 
pca_sglt <- prcomp(t(tse_sglt_clr))

### Combine PCs with metadata
sglt_beta_raw <- data.frame(pca_sglt$x[ , 1:5],
                            tse_sglt_samples)

## Prepare a for-loop for repeated measures ANOVA

### Create a list for beta diversity indices
beta_div <- c("PC1", "PC2", "PC3", "PC4", "PC5")

### Create a tibble for results
sglt_beta_anova_results <- tibble(x = 1:5) %>% 
  add_column(beta_div, .before = "x") %>% 
  dplyr::rename(p_value = x)

# Execute the for-loop

for (j in 1:length(beta_div)){
  
  sglt_temp_beta <- sglt_beta_raw %>% 
    rownames_to_column(var = "SampleID") %>% 
    select(c(2:8)) %>% 
    pivot_longer(cols = 1:5, names_to = "beta", values_to = "value") %>% 
    filter(beta == beta_div[j])
  
  sglt_beta_test <- sglt_temp_beta %>% 
    group_by(PatientID) %>% 
    filter(n() != 1) %>% 
    arrange(Timepoint, PatientID) %>% 
    ungroup()
  
  model_sglt <- rstatix::anova_test(data = sglt_beta_test, dv = value, wid = PatientID, within = Timepoint)
  
  sglt_beta_anova_results$p_value[j] <- model_sglt$ANOVA$p
  
}

### Correct p-values for multiple testing w/ Benjamini-Hochberg method
sglt_anova_beta_results <- sglt_beta_anova_results %>% 
  mutate(p_value_BH = p.adjust(p_value, method = "BH"))

## Prepare a for-loop for t-tests - only testing PC2 and PC5 (ANOVA < 0.05)

### Create a list for timepoints
timepoints <- c("II", "III", "IV")

### Create a list for beta diversity indices
beta_div <- c("PC1", "PC2", "PC3", "PC4", "PC5")

### Make a tibble for estimates
sglt_beta_ttest_estimates <- tibble(x = 1:2, y = 1:2, z = 1:2) %>% 
  add_column(beta_div, .before = "x") %>% 
  dplyr::rename(II = x,
                III = y,
                IV = z) %>% 
  column_to_rownames(var = "beta_div")

### Make a tibble for p-values
sglt_beta_ttest_pvalues <- tibble(x = 1:2, y = 1:2, z = 1:2) %>% 
  add_column(beta_div, .before = "x") %>% 
  dplyr::rename(II = x,
                III = y,
                IV = z) %>% 
  column_to_rownames(var = "beta_div")

### Execute the for-loop
for (j in 1:length(beta_div)){
  
  temp_beta_sglt <- sglt_beta_raw %>% 
    rownames_to_column(var = "SampleID") %>% 
    select(c(2:9)) %>% 
    pivot_longer(cols = 1:5, names_to = "beta_div", values_to = "value") %>% 
    filter(beta_div == beta_div[j])
  
  for (i in 1:length(timepoints)){
    beta_test_sglt <- temp_beta_sglt %>% 
      filter(Timepoint %in% c("I", timepoints[i])) %>% 
      group_by(PatientID) %>% 
      filter(n() != 1) %>% 
      arrange(Timepoint, PatientID) %>% 
      ungroup()
    
    beta_test_sglt_raw <- t.test(value ~ Timepoint, data = beta_test_sglt, paired = TRUE) 
    
    beta_resuls_sglt <- beta_test_sglt_raw %>% 
      broom::tidy()
    
    sglt_beta_ttest_estimates[beta_div[j], timepoints[i]] <- beta_resuls_sglt$estimate
    sglt_beta_ttest_pvalues[beta_div[j], timepoints[i]] <- beta_resuls_sglt$p.value
    
    sglt_beta_estimates <- sglt_beta_ttest_estimates %>% 
      rownames_to_column(var = "Beta_div") %>% 
      pivot_longer(cols = 2:4, names_to = "Timepoint", values_to = "Estimate")
    
    sglt_beta_pvalues <- sglt_beta_ttest_pvalues %>% 
      rownames_to_column(var = "Beta_div") %>% 
      pivot_longer(cols = 2:4, names_to = "Timepoint", values_to = "p_value")
    
    sglt_beta_ttest_results <- sglt_beta_estimates %>% 
      mutate(p_value = sglt_beta_pvalues$p_value)
    
  }
}

### Correct p-values for multiple testing w/ Benjamini-Hochberg method
sglt_beta_ttest_BH <- sglt_beta_ttest_results %>% 
  mutate(p_value_BH = p.adjust(p_value, method = "BH"))

# ___________________________________________________________________________ #

# Summary results

beta_anova_results <- cbind(glp_anova_beta_results, sglt_anova_beta_results)

t_test_alpha_results <- cbind(glp_beta_ttest_BH, sglt_beta_ttest_BH)
## Cannot combine them - tested PC2 and PC5 only for SGLT-2 (ANOVA < 0.05)

# ___________________________________________________________________________ #

# Visualize the results

## GLP-1-RA

glp_beta_plot_data <- glp_alpha_beta_raw %>% 
  rownames_to_column(var = "SampleID") %>% 
  select(c(2:9)) %>% 
  pivot_longer(cols = 1:5, names_to = "beta_div", values_to = "value") %>% 
  pivot_wider(names_from = Timepoint, values_from = value) %>% 
  mutate(I_II = I - II, 
         I_III = I - III,
         I_IV = I - IV) %>% 
  pivot_longer(cols = 4:7, names_to = "Timepoint", values_to = "Value") %>% 
  pivot_longer(cols = 4:6, names_to = "chg", values_to = "chg_value")


glp_PC1 <- glp_beta_plot_data %>% 
  filter(beta_div == "PC1") %>% 
  na.omit() %>% 
  group_by(chg) %>% 
  ggplot(aes(x = chg, y = chg_value)) +
  facet_grid(Medication~beta_div, switch = "y") +
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.3) +
  geom_hline(aes(yintercept = 0, colour = "red"), show.legend = F, linetype = "dashed") +
  scale_x_discrete(labels = c("BL vs M1", "BL vs M3", "BL vs M12")) +
  labs(x = NULL,
       y = "Change in value") +
  annotate(geom = "text", label = "ns", x = 1, y = 0.75) +
  annotate(geom = "text", label = "ns", x = 2, y = 0.75) +
  annotate(geom = "text", label = "ns", x = 3, y = 0.75)

glp_PC2 <- glp_beta_plot_data %>% 
  filter(beta_div == "PC2") %>% 
  na.omit() %>% 
  group_by(chg) %>% 
  ggplot(aes(x = chg, y = chg_value)) +
  facet_grid(Medication~beta_div, switch = "y") +
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.3) +
  geom_hline(aes(yintercept = 0, colour = "red"), show.legend = F, linetype = "dashed") +
  scale_x_discrete(labels = c("BL vs M1", "BL vs M3", "BL vs M12")) +
  labs(x = NULL,
       y = "Change in value") +
  annotate(geom = "text", label = "ns", x = 1, y = 0.1) +
  annotate(geom = "text", label = "ns", x = 2, y = 0.1) +
  annotate(geom = "text", label = "ns", x = 3, y = 0.1)

glp_PC3 <- glp_beta_plot_data %>% 
  filter(beta_div == "PC3") %>% 
  na.omit() %>% 
  group_by(chg) %>% 
  ggplot(aes(x = chg, y = chg_value)) +
  facet_grid(Medication~beta_div, switch = "y") +
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.3) +
  geom_hline(aes(yintercept = 0, colour = "red"), show.legend = F, linetype = "dashed") +
  scale_x_discrete(labels = c("BL vs M1", "BL vs M3", "BL vs M12")) +
  labs(x = NULL,
       y = "Change in value") +
  annotate(geom = "text", label = "ns", x = 1, y = 30) +
  annotate(geom = "text", label = "ns", x = 2, y = 30) +
  annotate(geom = "text", label = "ns", x = 3, y = 30)

glp_PC4 <- glp_beta_plot_data %>% 
  filter(beta_div == "PC4") %>% 
  na.omit() %>% 
  group_by(chg) %>% 
  ggplot(aes(x = chg, y = chg_value)) +
  facet_grid(Medication~beta_div, switch = "y") +
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.3) +
  geom_hline(aes(yintercept = 0, colour = "red"), show.legend = F, linetype = "dashed") +
  scale_x_discrete(labels = c("BL vs M1", "BL vs M3", "BL vs M12")) +
  labs(x = NULL,
       y = "Change in value") 

glp_PC5 <- glp_beta_plot_data %>% 
  filter(beta_div == "PC5") %>% 
  na.omit() %>% 
  group_by(chg) %>% 
  ggplot(aes(x = chg, y = chg_value)) +
  facet_grid(Medication~beta_div, switch = "y") +
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.3) +
  geom_hline(aes(yintercept = 0, colour = "red"), show.legend = F, linetype = "dashed") +
  scale_x_discrete(labels = c("BL vs M1", "BL vs M3", "BL vs M12")) +
  labs(x = NULL,
       y = "Change in value")

## SGLT-2

sglt_beta_plot_data <- sglt_alpha_beta_raw %>% 
  rownames_to_column(var = "SampleID") %>% 
  select(c(2:9)) %>% 
  pivot_longer(cols = 1:5, names_to = "beta_div", values_to = "value") %>% 
  pivot_wider(names_from = Timepoint, values_from = value) %>% 
  mutate(I_II = I - II, 
         I_III = I - III,
         I_IV = I - IV) %>% 
  pivot_longer(cols = 4:7, names_to = "Timepoint", values_to = "Value") %>% 
  pivot_longer(cols = 4:6, names_to = "chg", values_to = "chg_value")


sglt_PC1 <- sglt_beta_plot_data %>% 
  filter(beta_div == "PC1") %>% 
  na.omit() %>% 
  group_by(chg) %>% 
  ggplot(aes(x = chg, y = chg_value)) +
  facet_grid(Medication~beta_div, switch = "y") +
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.3) +
  geom_hline(aes(yintercept = 0, colour = "red"), show.legend = F, linetype = "dashed") +
  scale_x_discrete(labels = c("BL vs M1", "BL vs M3", "BL vs M12")) +
  labs(x = NULL,
       y = "Change in value") +
  annotate(geom = "text", label = "ns", x = 1, y = 0.6) +
  annotate(geom = "text", label = "ns", x = 2, y = 0.6) +
  annotate(geom = "text", label = "ns", x = 3, y = 0.6)

sglt_PC2 <- sglt_beta_plot_data %>% 
  filter(beta_div == "PC2") %>% 
  na.omit() %>% 
  group_by(chg) %>% 
  ggplot(aes(x = chg, y = chg_value)) +
  facet_grid(Medication~beta_div, switch = "y") +
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.3) +
  geom_hline(aes(yintercept = 0, colour = "red"), show.legend = F, linetype = "dashed") +
  scale_x_discrete(labels = c("BL vs M1", "BL vs M3", "BL vs M12")) +
  labs(x = NULL,
       y = "Change in value") +
  annotate(geom = "text", label = "ns", x = 1, y = 0.15) +
  annotate(geom = "text", label = "ns", x = 2, y = 0.15) +
  annotate(geom = "text", label = "ns", x = 3, y = 0.15)

sglt_PC3 <- sglt_beta_plot_data %>% 
  filter(beta_div == "PC3") %>% 
  na.omit() %>% 
  group_by(chg) %>% 
  ggplot(aes(x = chg, y = chg_value)) +
  facet_grid(Medication~beta_div, switch = "y") +
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.3) +
  geom_hline(aes(yintercept = 0, colour = "red"), show.legend = F, linetype = "dashed") +
  scale_x_discrete(labels = c("BL vs M1", "BL vs M3", "BL vs M12")) +
  labs(x = NULL,
       y = "Change in value") + 
  annotate(geom = "text", label = "ns", x = 1, y = 10) +
  annotate(geom = "text", label = "ns", x = 2, y = 10) +
  annotate(geom = "text", label = "ns", x = 3, y = 10) 

sglt_PC4 <- sglt_beta_plot_data %>% 
  filter(beta_div == "PC4") %>% 
  na.omit() %>% 
  group_by(chg) %>% 
  ggplot(aes(x = chg, y = chg_value)) +
  facet_grid(Medication~beta_div, switch = "y") +
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.3) +
  geom_hline(aes(yintercept = 0, colour = "red"), show.legend = F, linetype = "dashed") +
  scale_x_discrete(labels = c("BL vs M1", "BL vs M3", "BL vs M12")) +
  labs(x = NULL,
       y = "Change in value")

sglt_PC5 <- sglt_beta_plot_data %>% 
  filter(beta_div == "PC5") %>% 
  na.omit() %>% 
  group_by(chg) %>% 
  ggplot(aes(x = chg, y = chg_value)) +
  facet_grid(Medication~beta_div, switch = "y") +
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.3) +
  geom_hline(aes(yintercept = 0, colour = "red"), show.legend = F, linetype = "dashed") +
  scale_x_discrete(labels = c("BL vs M1", "BL vs M3", "BL vs M12")) +
  labs(x = NULL,
       y = "Change in value")

common_alpha <- ggpubr::ggarrange(glp_alpha_shannon, glp_alpha_pielou, glp_alpha_observed, 
                                  sglt_alpha_shannon, sglt_alpha_pielou, sglt_alpha_observed, 
                                  ncol = 3, nrow = 2)

