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

## Remove features containing human DNA and unassigned taxa
tse_prelim <- tse_taxa[!rowData(tse_taxa)$Phylum == "" & 
                         !rowData(tse_taxa)$Genus == "" &
                         !rowData(tse_taxa)$Kingdom == "Unassigned", ]

## Exclude samples
tse <- tse_prelim[ , !colnames(tse_prelim) %in% c("GLP1RA-5-2", "GLP1RA-5-3", 
                                                  "GLP1RA-5-4", "GLP1RA-6-2", 
                                                  "NEG-Valida", "GLP1RA-15-2-2", 
                                                  "Elini-proov", "AL")]


# ___________________________________________________________________________ #

# Convert counts into relative abundances
tse <- transformAssay(tse, assay.type = "counts", method = "relabundance")

# Convert relative abundances into CLR-transformed values
mat <- assay(tse, "relabundance")
tse <- transformAssay(tse, assay.type = "relabundance", method = "clr", 
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

### Create a tibble for general effect sizes (ges)
glp_beta_anova_ges <- tibble(x = 1:5) %>% 
  add_column(beta_div, .before = "x") %>% 
  dplyr::rename(ges = x)

glp_beta_anova_pvalues <- tibble(x = 1:5) %>% 
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
  
  glp_beta_anova_ges$ges[j] <- model_glp$ANOVA$ges
  glp_beta_anova_pvalues$p_value[j] <- model_glp$ANOVA$p
  
  glp_beta_anova_results <- cbind(glp_beta_anova_ges, glp_beta_anova_pvalues)
  
}

### Correct p-values for multiple testing w/ Benjamini-Hochberg method

glp_beta_anova_results <- glp_beta_anova_results[-3]

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

### Create a tibble for general effect sizes (ges)
sglt_beta_anova_ges <- tibble(x = 1:5) %>% 
  add_column(beta_div, .before = "x") %>% 
  dplyr::rename(ges = x)

sglt_beta_anova_pvalues <- tibble(x = 1:5) %>% 
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
  
  sglt_beta_anova_ges$ges[j] <- model_sglt$ANOVA$ges
  sglt_beta_anova_pvalues$p_value[j] <- model_sglt$ANOVA$p
  
  sglt_beta_anova_results <- cbind(sglt_beta_anova_ges, sglt_beta_anova_pvalues)
  
}

### Correct p-values for multiple testing w/ Benjamini-Hochberg method

sglt_beta_anova_results <- sglt_beta_anova_results[-3]

sglt_anova_beta_results <- sglt_beta_anova_results %>% 
  mutate(p_value_BH = p.adjust(p_value, method = "BH"))

## Prepare a for-loop for t-tests - only testing PC2 and PC5 (ANOVA < 0.05)

### Create a list for timepoints
timepoints <- c("II", "III", "IV")

### Create a list for beta diversity indices
beta_div <- c("PC1", "PC2", "PC3", "PC4", "PC5")

### Make a tibble for estimates
sglt_beta_ttest_estimates <- tibble(x = 1:5, y = 1:5, z = 1:5) %>% 
  add_column(beta_div, .before = "x") %>% 
  dplyr::rename(II = x,
                III = y,
                IV = z) %>% 
  column_to_rownames(var = "beta_div")

### Make a tibble for p-values
sglt_beta_ttest_pvalues <- tibble(x = 1:5, y = 1:5, z = 1:5) %>% 
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

t_test_beta_results <- cbind(glp_beta_ttest_BH, sglt_beta_ttest_BH)

# writexl::write_xlsx(beta_anova_results, "anova_beta_results.xlsx")
# writexl::write_xlsx(t_test_beta_results, "t_test_beta_results.xlsx")

# ___________________________________________________________________________ #

# Visualize the results

## GLP-1-RA

# paletteer::paletteer_d("fishualize::Thunnus_obesus") 

#  #669ABFFF #D4D05BFF #FDF6E2FF #747798FF #534D56FF 

glp_beta_data <- glp_beta_raw %>% 
  rownames_to_column(var = "SampleID") %>% 
  select(c(2:9)) %>%
  pivot_longer(cols = 1:5, names_to = "beta_div", values_to = "value")

glp_PC1_plot <- glp_beta_data %>% 
  filter(beta_div == "PC1") %>% 
  ggplot(aes(x = Timepoint, y = value, fill = beta_div)) +
  scale_fill_manual(values = "#7AD151FF") +
  facet_grid(Medication~beta_div, switch = "y", scale = "free") +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.3) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "",
       y = "") +
  ggsignif::geom_signif(
    y_position = c(150, 180, 210), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.14", "0.04", "0.16"), tip_length = 0.02) +
  guides(fill = "none") +
  theme_classic() +
  theme(strip.text.x = element_text(size = 16), 
        strip.text.y = element_text(size = 16), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))

glp_PC2_plot <- glp_beta_data %>% 
  filter(beta_div == "PC2") %>% 
  ggplot(aes(x = Timepoint, y = value, fill = beta_div)) +
  scale_fill_manual(values = "#7AD151FF") +
  facet_grid(.~beta_div, switch = "y", scale = "free") +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "",
       y = "") +
  ggsignif::geom_signif(
    y_position = c(100, 120, 140), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.62", "0.29", "0.38"), tip_length = 0.02) +
  guides(fill = "none") +
  theme_classic() +
  theme(strip.text.x = element_text(size = 16), 
        strip.text.y = element_text(size = 16), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))

glp_PC3_plot <- glp_beta_data %>% 
  filter(beta_div == "PC3") %>% 
  ggplot(aes(x = Timepoint, y = value, fill = beta_div)) +
  scale_fill_manual(values = "#7AD151FF") +
  facet_grid(.~beta_div, switch = "y", scale = "free") +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "",
       y = "") +
  ggsignif::geom_signif(
    y_position = c(100, 120, 140), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.84", "0.87", "0.75"), tip_length = 0.02) +
  guides(fill = "none") +
  theme_classic() +
  theme(strip.text.x = element_text(size = 16), 
        strip.text.y = element_text(size = 16), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))

glp_PC4_plot <- glp_beta_data %>% 
  filter(beta_div == "PC4") %>% 
  ggplot(aes(x = Timepoint, y = value, fill = beta_div)) +
  scale_fill_manual(values = "#7AD151FF") +
  facet_grid(.~beta_div, switch = "y", scale = "free") +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "",
       y = "") +
  ggsignif::geom_signif(
    y_position = c(100, 120, 140), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.56", "0.7", "0.18"), tip_length = 0.02) +
  guides(fill = "none") +
  theme_classic() +
  theme(strip.text.x = element_text(size = 16), 
        strip.text.y = element_text(size = 16), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))

glp_PC5_plot <- glp_beta_data %>% 
  filter(beta_div == "PC5") %>% 
  ggplot(aes(x = Timepoint, y = value, fill = beta_div)) +
  scale_fill_manual(values = "#7AD151FF") +
  facet_grid(.~beta_div, switch = "y", scale = "free") +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "",
       y = "") +
  ggsignif::geom_signif(
    y_position = c(50, 60, 70), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.6", "0.64", "0.59"), tip_length = 0.02) +
  guides(fill = "none") +
  theme_classic() +
  theme(strip.text.x = element_text(size = 16), 
        strip.text.y = element_text(size = 16), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))

glp_beta_comb <- ggpubr::ggarrange(glp_PC1_plot, glp_PC2_plot, glp_PC3_plot,
                                   nrow = 1, ncol = 3)

## SGLT-2

# paletteer::paletteer_d("fishualize::Thunnus_obesus") 

#  #669ABFFF #D4D05BFF #FDF6E2FF #747798FF #534D56FF 

sglt_beta_data <- sglt_beta_raw %>% 
  rownames_to_column(var = "SampleID") %>% 
  select(c(2:9)) %>%
  pivot_longer(cols = 1:5, names_to = "beta_div", values_to = "value") %>% 
  mutate(Medication = recode(Medication, 
                             "SGLT-2" = "SGLT-2i"))

sglt_PC1_plot <- sglt_beta_data %>% 
  filter(beta_div == "PC1") %>% 
  ggplot(aes(x = Timepoint, y = value, fill = beta_div)) +
  scale_fill_manual(values = "#669ABFFF") +
  facet_grid(Medication~beta_div, switch = "y", scale = "free") +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "",
       y = "") +
  guides(fill = "none") +
  ggsignif::geom_signif(
    y_position = c(130, 160, 190), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.66", "0.12", "0.46"), tip_length = 0.02) +
  theme_classic() +
  theme(strip.text.x = element_text(size = 16), 
        strip.text.y = element_text(size = 16), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))

sglt_PC2_plot <- sglt_beta_data %>% 
  filter(beta_div == "PC2") %>% 
  ggplot(aes(x = Timepoint, y = value, fill = beta_div)) +
  scale_fill_manual(values = "#669ABFFF") +
  facet_grid(.~beta_div, switch = "y", scale = "free") +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "",
       y = "") +
  guides(fill = "none") +
  ggsignif::geom_signif(
    y_position = c(100, 120, 140), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.51", "0.47", "0.34"), tip_length = 0.02) +
  theme_classic() +
  theme(strip.text.x = element_text(size = 16), 
        strip.text.y = element_text(size = 16), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))

sglt_PC3_plot <- sglt_beta_data %>% 
  filter(beta_div == "PC3") %>% 
  ggplot(aes(x = Timepoint, y = value, fill = beta_div)) +
  scale_fill_manual(values = "#669ABFFF") +
  facet_grid(.~beta_div, switch = "y", scale = "free") +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "",
       y = "") +
  guides(fill = "none") +
  ggsignif::geom_signif(
    y_position = c(100, 120, 140), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.32", "0.44", "0.15"), tip_length = 0.02) +
  theme_classic() +
  theme(strip.text.x = element_text(size = 16), 
        strip.text.y = element_text(size = 16), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))

sglt_PC4_plot <- sglt_beta_data %>% 
  filter(beta_div == "PC4") %>% 
  ggplot(aes(x = Timepoint, y = value, fill = beta_div)) +
  scale_fill_manual(values = "#669ABFFF") +
  facet_grid(.~beta_div, switch = "y", scale = "free") +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "",
       y = "") +
  guides(fill = "none") +
  ggsignif::geom_signif(
    y_position = c(50, 70, 90), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.72", "0.87", "0.56"), tip_length = 0.02) +
  theme_classic() +
  theme(strip.text.x = element_text(size = 16), 
        strip.text.y = element_text(size = 16), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))

sglt_PC5_plot <- sglt_beta_data %>% 
  filter(beta_div == "PC5") %>% 
  ggplot(aes(x = Timepoint, y = value, fill = beta_div)) +
  scale_fill_manual(values = "#669ABFFF") +
  facet_grid(.~beta_div, switch = "y", scale = "free") +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "",
       y = "") +
  guides(fill = "none") +
  ggsignif::geom_signif(
    y_position = c(80, 100, 120), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.57", "0.009", "0.10"), tip_length = 0.02) +
  theme_classic() +
  theme(strip.text.x = element_text(size = 16), 
        strip.text.y = element_text(size = 16), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))

sglt_beta_comb <- ggpubr::ggarrange(sglt_PC1_plot, sglt_PC2_plot, sglt_PC3_plot, 
                                    nrow = 1, ncol = 3)

beta_comb <- ggpubr::ggarrange(glp_beta_comb, sglt_beta_comb, 
                                nrow = 2, ncol = 1)

ggsave("combined_beta_plot.svg", device = "svg", dpi = 300, width = 11, height = 11)
