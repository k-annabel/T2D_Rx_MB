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
library(emmeans)

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

### Create a tibble for general effect sizes (ges)
glp_alpha_anova_ges <- tibble(x = 1:3) %>% 
  add_column(alpha_div, .before = "x") %>% 
  dplyr::rename(ges = x)

glp_alpha_anova_pvalues <- tibble(x = 1:3) %>% 
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
    
    glp_alpha_anova_ges$ges[j] <- model_glp$ANOVA$ges
    glp_alpha_anova_pvalues$p_value[j] <- model_glp$ANOVA$p
    
    glp_alpha_anova_results <- cbind(glp_alpha_anova_ges, glp_alpha_anova_pvalues)
    
}

### Correct p-values for multiple testing w/ Benjamini-Hochberg method

glp_alpha_anova_results <- glp_alpha_anova_results[-3]

glp_anova_alpha_results <- glp_alpha_anova_results %>% 
  mutate(p_value_BH = p.adjust(p_value, method = "BH"))

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
    
    glp_alpha_estimates <- glp_alpha_ttest_estimates %>% 
      rownames_to_column(var = "Alpha_div") %>% 
      pivot_longer(cols = 2:4, names_to = "Timepoint", values_to = "Estimate")
    
    glp_alpha_pvalues <- glp_alpha_ttest_pvalues %>% 
      rownames_to_column(var = "Alpha_div") %>% 
      pivot_longer(cols = 2:4, names_to = "Timepoint", values_to = "p_value")
    
    glp_alpha_ttest_results <- glp_alpha_estimates %>% 
      mutate(p_value = glp_alpha_pvalues$p_value)
    
  }
}

### Correct p-values for multiple testing w/ Benjamini-Hochberg method
glp_alpha_ttest_BH <- glp_alpha_ttest_results %>% 
  mutate(p_value_BH = p.adjust(p_value, method = "BH"))

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

### Create tibbles for results
sglt_alpha_anova_ges <- tibble(x = 1:3) %>% 
  add_column(alpha_div, .before = "x") %>% 
  dplyr::rename(ges = x)

sglt_alpha_anova_pvalues <- tibble(x = 1:3) %>% 
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
  
  sglt_alpha_anova_ges$ges[j] <- model_sglt$ANOVA$ges
  sglt_alpha_anova_pvalues$p_value[j] <- model_sglt$ANOVA$p
  
  sglt_alpha_anova_results <- cbind(sglt_alpha_anova_ges, sglt_alpha_anova_pvalues)
  
}

sglt_alpha_anova_results <- sglt_alpha_anova_results[-3]

### Correct p-values for multiple testing w/ Benjamini-Hochberg method
sglt_anova_alpha_results <- sglt_alpha_anova_results %>% 
  mutate(p_value_BH = p.adjust(p_value, method = "BH"))

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
    
    
    sglt_alpha_estimates <- sglt_alpha_ttest_estimates %>% 
      rownames_to_column(var = "Alpha_div") %>% 
      pivot_longer(cols = 2:4, names_to = "Timepoint", values_to = "Estimate")
    
    sglt_alpha_pvalues <- sglt_alpha_ttest_pvalues %>% 
      rownames_to_column(var = "Alpha_div") %>% 
      pivot_longer(cols = 2:4, names_to = "Timepoint", values_to = "p_value")
    
    
    sglt_alpha_ttest_results <- sglt_alpha_estimates %>% 
      mutate(p_value = sglt_alpha_pvalues$p_value)
    
  }
}

### Correct p-values for multiple testing w/ Benjamini-Hochberg method
sglt_alpha_ttest_BH <- sglt_alpha_ttest_results %>% 
  mutate(p_value_BH = p.adjust(p_value, method = "BH"))

# ___________________________________________________________________________ #

# Summary results

anova_alpha_results <- cbind(glp_anova_alpha_results, sglt_anova_alpha_results)

t_test_alpha_results <- cbind(glp_alpha_ttest_BH, sglt_alpha_ttest_BH)

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
  pivot_longer(cols = 4:6, names_to = "chg", values_to = "chg_value") %>% 
  mutate(alpha_div = recode(alpha_div, 
                            "shannon" = "Shannon",
                            "pielou" = "Pielou",
                            "observed" = "Observed"))


glp_alpha_shannon <- glp_alpha_plot_data %>% 
  filter(alpha_div == "shannon") %>% 
  na.omit() %>% 
  group_by(chg) %>% 
  ggplot(aes(x = chg, y = chg_value)) +
  facet_grid(Medication~alpha_div, switch = "y") +
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
  
glp_alpha_pielou <- glp_alpha_plot_data %>% 
  filter(alpha_div == "Pielou") %>% 
  na.omit() %>% 
  group_by(chg) %>% 
  ggplot(aes(x = chg, y = chg_value)) +
  facet_grid(.~alpha_div, switch = "y") +
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.3) +
  geom_hline(aes(yintercept = 0, colour = "red"), show.legend = F, linetype = "dashed") +
  scale_x_discrete(labels = c("BL vs M1", "BL vs M3", "BL vs M12")) +
  labs(x = NULL,
       y = NULL)  +
  annotate(geom = "text", label = "ns", x = 1, y = 0.1) +
  annotate(geom = "text", label = "ns", x = 2, y = 0.1) +
  annotate(geom = "text", label = "ns", x = 3, y = 0.1)

glp_alpha_observed <- glp_alpha_plot_data %>% 
  filter(alpha_div == "Observed") %>% 
  na.omit() %>% 
  group_by(chg) %>% 
  ggplot(aes(x = chg, y = chg_value)) +
  facet_grid(.~alpha_div, switch = "y") +
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.3) +
  geom_hline(aes(yintercept = 0, colour = "red"), show.legend = F, linetype = "dashed") +
  scale_x_discrete(labels = c("BL vs M1", "BL vs M3", "BL vs M12")) +
  labs(x = NULL,
       y = NULL) +
  annotate(geom = "text", label = "ns", x = 1, y = 30) +
  annotate(geom = "text", label = "ns", x = 2, y = 30) +
  annotate(geom = "text", label = "ns", x = 3, y = 30)

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
  pivot_longer(cols = 4:6, names_to = "chg", values_to = "chg_value") %>% 
  mutate(alpha_div = recode(alpha_div, 
                            "shannon" = "Shannon",
                            "pielou" = "Pielou",
                            "observed" = "Observed"))


sglt_alpha_shannon <- sglt_alpha_plot_data %>% 
  filter(alpha_div == "Shannon") %>% 
  na.omit() %>% 
  group_by(chg) %>% 
  ggplot(aes(x = chg, y = chg_value)) +
  facet_grid(Medication~alpha_div, switch = "y") +
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

sglt_alpha_pielou <- sglt_alpha_plot_data %>% 
  filter(alpha_div == "Pielou") %>% 
  na.omit() %>% 
  group_by(chg) %>% 
  ggplot(aes(x = chg, y = chg_value)) +
  facet_grid(.~alpha_div, switch = "y") +
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.3) +
  geom_hline(aes(yintercept = 0, colour = "red"), show.legend = F, linetype = "dashed") +
  xlab("Change in timepoint compared to BL") +
  theme(axis.title.x = element_text(vjust = -1, size = 11), axis.title.y = element_blank()) +
  scale_x_discrete(labels = c("BL vs M1", "BL vs M3", "BL vs M12")) +
  annotate(geom = "text", label = "ns", x = 1, y = 0.15) +
  annotate(geom = "text", label = "ns", x = 2, y = 0.15) +
  annotate(geom = "text", label = "ns", x = 3, y = 0.15)

sglt_alpha_observed <- sglt_alpha_plot_data %>% 
  filter(alpha_div == "Observed") %>% 
  na.omit() %>% 
  group_by(chg) %>% 
  ggplot(aes(x = chg, y = chg_value)) +
  facet_grid(.~alpha_div, switch = "y") +
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.3) +
  geom_hline(aes(yintercept = 0, colour = "red"), show.legend = F, linetype = "dashed") +
  scale_x_discrete(labels = c("BL vs M1", "BL vs M3", "BL vs M12")) +
  labs(x = NULL,
       y = NULL) +
  annotate(geom = "text", label = "ns", x = 1, y = 10) +
  annotate(geom = "text", label = "ns", x = 2, y = 10) +
  annotate(geom = "text", label = "ns", x = 3, y = 10)

common_alpha <- ggpubr::ggarrange(glp_alpha_shannon, glp_alpha_pielou, glp_alpha_observed, 
                                  sglt_alpha_shannon, sglt_alpha_pielou, sglt_alpha_observed, 
                                  ncol = 3, nrow = 2)

# New version

## GLP-1-RA

# paletteer_d("ggprism::viridis")

  #440154FF #414487FF #2A788EFF #22A884FF #7AD151FF #FDE725FF

glp_alpha_data <- glp_alpha_beta_raw %>% 
  rownames_to_column(var = "SampleID") %>% 
  select(c(7:9, 72:74)) %>% 
  pivot_longer(cols = 4:6, names_to = "alpha_div", values_to = "value") %>% 
  mutate(alpha_div = recode(alpha_div, 
                            "shannon" = "Shannon",
                            "pielou" = "Pielou",
                            "observed" = "Observed"))

glp_alpha_observed <- glp_alpha_data %>% 
  filter(alpha_div == "Observed") %>% 
  ggplot(aes(x = Timepoint, y = value, fill = alpha_div)) +
  facet_grid(Medication~alpha_div, switch = "y", scales = "free") +
  scale_fill_manual(values = "#7AD151FF") +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.3) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "",
       y = "Value") +
  guides(fill = "none") +
  ggsignif::geom_signif(
    y_position = c(90, 95, 100), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.68", "0.37", "0.29"), tip_length = 0.02) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))


glp_alpha_shannon <- glp_alpha_data %>% 
  filter(alpha_div == "Shannon") %>% 
  ggplot(aes(x = Timepoint, y = value, fill = alpha_div)) +
  scale_fill_manual(values = "#7AD151FF") +
  facet_grid(.~alpha_div, switch = "y", scales = "free") +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.3) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  guides(fill = "none") +
  labs(x = "",
       y = "") +
  ggsignif::geom_signif(
    y_position = c(3.6, 3.8, 4.0), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.47", "0.12", "0.51"), tip_length = 0.02) +
  theme(strip.text = element_text(size = 12))

glp_alpha_pielou <- glp_alpha_data %>% 
  filter(alpha_div == "Pielou") %>% 
  ggplot(aes(x = Timepoint, y = value, fill = alpha_div)) +
  facet_grid(.~alpha_div, switch = "y", scales = "free") +
  scale_fill_manual(values = "#7AD151FF") +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.3) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "",
       y = "") +
  guides(fill = "none") +
  ggsignif::geom_signif(
    y_position = c(0.85, 0.9, 0.95), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.5", "0.07", "0.67"), tip_length = 0.02) +
  theme(strip.text = element_text(size = 12))

glp_alpha_comb <- ggpubr::ggarrange(glp_alpha_observed, glp_alpha_shannon,
                                    nrow = 1, ncol = 2)

## SGLT-2

sglt_alpha_data <- sglt_alpha_beta_raw %>% 
  rownames_to_column(var = "SampleID") %>% 
  select(c(7:9, 72:74)) %>% 
  pivot_longer(cols = 4:6, names_to = "alpha_div", values_to = "value") %>% 
  mutate(alpha_div = recode(alpha_div, 
                            "shannon" = "Shannon",
                            "pielou" = "Pielou",
                            "observed" = "Observed"))

sglt_alpha_observed <- sglt_alpha_data %>% 
  filter(alpha_div == "Observed") %>% 
  ggplot(aes(x = Timepoint, y = value, fill = alpha_div)) +
  scale_fill_manual(values = "#7AD151FF") +
  facet_grid(Medication~alpha_div, switch = "y", scales = "free") +
  geom_boxplot() +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.3) +
  geom_point(position = position_jitter(width = 0.02)) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "Timepoint",
       y = "Value") +
  guides(fill = "none") +
  ggsignif::geom_signif(
    y_position = c(102, 106, 110), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.87", "0.58", "0.93"), tip_length = 0.02) +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))

sglt_alpha_shannon <- sglt_alpha_data %>% 
  filter(alpha_div == "Shannon") %>% 
  ggplot(aes(x = Timepoint, y = value, fill = alpha_div)) +
  scale_fill_manual(values = "#7AD151FF") +
  facet_grid(.~alpha_div, switch = "y", scales = "free") +
  geom_boxplot() +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  geom_point(position = position_jitter(width = 0.02)) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "",
       y = "") +
  guides(fill = "none") +
  ggsignif::geom_signif(
    y_position = c(3.6, 3.7, 3.8), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.55", "0.68", "0.08"), tip_length = 0.02) +
  theme(strip.text = element_text(size = 12))

sglt_alpha_pielou <- sglt_alpha_data %>% 
  filter(alpha_div == "Pielou") %>% 
  ggplot(aes(x = Timepoint, y = value, fill = alpha_div)) +
  scale_fill_manual(values = "#7AD151FF") +
  facet_grid(.~alpha_div, switch = "y", scales = "free") +
  geom_boxplot() +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  geom_point(position = position_jitter(width = 0.02)) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "",
       y = "") +
  guides(fill = "none") +
  ggsignif::geom_signif(
    y_position = c(0.82, 0.84, 0.86), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.53", "0.71", "0.08"), tip_length = 0.02) +
  theme(strip.text = element_text(size = 12))

sglt_alpha_comb <- ggpubr::ggarrange(sglt_alpha_observed, sglt_alpha_shannon,
                                     nrow = 1, ncol = 2)

alpha_comb <- ggpubr::ggarrange(glp_alpha_comb, sglt_alpha_comb, 
                                nrow = 2, ncol = 1)

ggsave("combined_alpha_plot.png", device = "png", dpi = 300, width = 12, height = 11)
