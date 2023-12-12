# T2D-Rx_MB_analysis - effect of T2D-Rx towards microbiome - differential abundance analysis

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
tse <- transformAssay(tse, assay.type = "counts", method = "relabundance")

# Convert relative abundances into CLR-transformed values
mat <- assay(tse, "relabundance")
tse <- transformAssay(tse, assay.type = "relabundance", method = "clr", 
                       pseudocount = min(mat[mat>0]))

# Collapse into Genus level
tse_genus <- agglomerateByRank(tse, rank = "Genus")

# Remove "Genus:" label
rownames(tse_genus) <- sub("Genus:", "", rownames(tse_genus))

# Remove taxa with Genus observation empty
tse_genus <- tse_genus[rowData(tse_genus)$Genus != "", ]

# Separate by medication
tse_glp <- tse_genus[ , colData(tse_genus)$Medication == "GLP-1-RA"]
tse_sglt <- tse_genus[ , colData(tse_genus)$Medication == "SGLT-2"]

# ___________________________________________________________________________ #

# Differential abundance analysis with repeated measures ANOVA and t-tests

# GLP-1-RA

## Gather top taxa (present in at least 10 samples)
tse_glp <- transformAssay(tse_glp, method = "pa")
glp_top_taxa <- names(rowSums(assay(tse_glp, "pa")))[rowSums(assay(tse_glp, "pa")) > 10]

# Clean and transform relative abundance data corresponding to most prevalent genera
glp_genera_comparisons <- assay(tse_glp, "clr") %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Genus") %>% 
  mutate(Genus = rownames(assay(tse_glp, "clr"))) %>% 
  filter(Genus %in% glp_top_taxa) %>% 
  pivot_longer(cols = 2:36, 
               names_to = "SampleID", 
               values_to = "clr") %>%
  mutate(Timepoint = str_extract(SampleID, "(\\d+)$")) %>% 
  mutate(Timepoint = recode(Timepoint, 
                            "1" = "I", 
                            "2" = "II", 
                            "3" = "III", 
                            "4" = "IV")) %>% 
  mutate(PatientID = str_extract(SampleID, "GLP1RA-\\d+")) %>% 
  relocate(PatientID, .before = "SampleID") %>% 
  select(-SampleID)


# Remove duplicate genera
# glp_top_taxa2 <- glp_top_taxa[!(glp_top_taxa %in% c("uncultured", "uncultured_1"))]

# Create a tibble for results
glp_results_estimates <- tibble(x = 1:128) %>% 
  rownames_to_column(var = "Genus") %>% 
  dplyr::rename(p_value = x) %>% 
  add_column(glp_top_taxa, .before = "p_value") %>% 
  select(-Genus) %>% 
  dplyr::rename(Genus = glp_top_taxa)

glp_results_pvalue <- tibble(x = 1:128) %>% 
  rownames_to_column(var = "Genus") %>% 
  dplyr::rename(p_value = x) %>% 
  add_column(glp_top_taxa, .before = "p_value") %>% 
  select(-Genus) %>% 
  dplyr::rename(Genus = glp_top_taxa)


# Execute the for-loop for repeated measures ANOVA

for (j in 1:length(glp_top_taxa)){
  
  temp_data <- glp_genera_comparisons %>% 
    dplyr::filter(Genus == glp_top_taxa[j])
  
  df_glp <- temp_data %>% 
    group_by(PatientID) %>% 
    filter(n() != 1) %>% 
    arrange(Timepoint, PatientID) %>% 
    ungroup()
  
  #GLP-1-RA
  model_glp <- rstatix::anova_test(data = df_glp, dv = clr, wid = PatientID, within = Timepoint)
  
  glp_results_da$p_value[j] <- model_glp$ANOVA$p
  
}

### Correct p-values for multiple testing w/ Benjamini-Hochberg method
glp_anova_da_results <- glp_results_da %>% 
  filter(p_value < 0.05) %>% 
  mutate(p_value_BH = p.adjust(p_value, method = "BH"))

# ___________________________________________________________________________ #

# Perform t-tests for five significant genera in ANOVA results

### Initialize empty vector for results
timepoints <- c("II", "III", "IV")

### Pull five significant genera
glp_anova_genera <- glp_results_da %>% 
  filter(p_value < 0.05) %>% 
  pull(Genus)

### Make a tibble for estimates
glp_da_ttest_estimates <- tibble(x = 1:5, y = 1:5, z = 1:5) %>% 
  add_column(glp_anova_genera, .before = "x") %>% 
  dplyr::rename(II = x,
                III = y,
                IV = z) %>% 
  column_to_rownames(var = "glp_anova_genera")

### Make a tibble for p-values
glp_da_ttest_pvalues <- tibble(x = 1:5, y = 1:5, z = 1:5) %>% 
  add_column(glp_anova_genera, .before = "x") %>% 
  dplyr::rename(II = x,
                III = y,
                IV = z) %>% 
  column_to_rownames(var = "glp_anova_genera")


for (j in 1:length(glp_anova_genera)){
  
  temp_data <- glp_genera_comparisons %>% 
    filter(Genus == glp_anova_genera[j])
  
  for (i in 1:length(timepoints)){
    df_glp <- temp_data %>% 
      filter(Timepoint %in% c("I", timepoints[i])) %>% 
      group_by(PatientID) %>% 
      filter(n() != 1) %>% 
      arrange(Timepoint, PatientID) %>% 
      ungroup()
    
    da_test_glp_raw <- t.test(clr ~ Timepoint, data = df_glp, paired = TRUE) 
    
    da_results_glp <- da_test_glp_raw %>% 
      broom::tidy()
    
    glp_da_ttest_estimates[glp_anova_genera[j], timepoints[i]] <- da_results_glp$estimate
    glp_da_ttest_pvalues[glp_anova_genera[j], timepoints[i]] <- da_results_glp$p.value

  }
}

### Correct p-values for multiple testing w/ Benjamini-Hochberg method
glp_da_estim <- glp_da_ttest_estimates %>% 
  rownames_to_column(var = "Genus_GLP") %>% 
  pivot_longer(cols = 2:4, names_to = "timepoint_GLP", values_to = "estimate_GLP")

glp_da_pvalues <- glp_da_ttest_pvalues %>% 
  rownames_to_column(var = "Genus") %>% 
  pivot_longer(cols = 2:4, names_to = "timepoint", values_to = "p_value_GLP") %>% 
  select(-c(Genus, timepoint)) %>% 
  mutate(p_value_BH_GLP = p.adjust(p_value_GLP, method = "BH"))

glp_da_ttest_BH <- bind_cols(glp_da_estim, glp_da_pvalues)

# ___________________________________________________________________________ #

# SGLT-2

## Gather top taxa (present in at least 10 samples)
tse_sglt <- transformAssay(tse_sglt, method = "pa")
sglt_top_taxa <- names(rowSums(assay(tse_sglt, "pa")))[rowSums(assay(tse_sglt, "pa")) > 10]

# Clean and transform relative abundance data corresponding to most prevalent genera
sglt_genera_comparisons <- assay(tse_sglt, "clr") %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Genus") %>% 
  mutate(Genus = rownames(assay(tse_sglt, "clr"))) %>%
  pivot_longer(cols = 2:39,
               names_to = "SampleID", 
               values_to = "clr") %>%
  mutate(Timepoint = str_extract(SampleID, "(\\d+)$")) %>% 
  mutate(Timepoint = recode(Timepoint, 
                            "1" = "I", 
                            "2" = "II", 
                            "3" = "III", 
                            "4" = "IV")) %>% 
  mutate(PatientID = str_extract(SampleID, "GLP1RA-\\d+")) %>% 
  relocate(PatientID, .before = "SampleID") %>% 
  select(-SampleID)

# Perform repeated measures ANOVA

# Remove duplicate genera
# sglt_top_taxa2 <- sglt_top_taxa[!(sglt_top_taxa %in% c("uncultured", "uncultured_1"))]

# Create a tibble for results
sglt_results_da <- tibble(x = 1:147) %>% 
  rownames_to_column(var = "Genus") %>% 
  dplyr::rename(p_value = x) %>% 
  add_column(sglt_top_taxa, .before = "p_value") %>% 
  select(-Genus) %>% 
  dplyr::rename(Genus = sglt_top_taxa)

# Execute the for-loop for repeated measures ANOVA

for (j in 1:length(sglt_top_taxa)){
  
  temp_data <- sglt_genera_comparisons %>% 
    filter(Genus == sglt_top_taxa[j])
  
  df_sglt <- temp_data %>% 
    group_by(PatientID) %>% 
    filter(n() != 1) %>% 
    arrange(Timepoint, PatientID) %>% 
    ungroup()
  
  #SGLT-2
  model_sglt <- rstatix::anova_test(data = df_sglt, dv = clr, wid = PatientID, within = Timepoint)
  
  sglt_results_da$p_value[j] <- model_sglt$ANOVA$p
  
}

### Correct p-values for multiple testing w/ Benjamini-Hochberg method
sglt_anova_da_results <- sglt_results_da %>% 
  filter(p_value <= 0.05) %>% 
  mutate(p_value_BH = p.adjust(p_value, method = "BH"))

# Perform t-tests

# Initialize empty vector for results
timepoints <- c("II", "III", "IV")

# Pull five significant genera
sglt_anova_genera <- sglt_results_da %>% 
  filter(p_value <= 0.05) %>% 
  pull(Genus)

### Make a tibble for estimates
sglt_da_ttest_estimates <- tibble(x = 1:5, y = 1:5, z = 1:5) %>% 
  add_column(sglt_anova_genera, .before = "x") %>% 
  dplyr::rename(II = x,
                III = y,
                IV = z) %>% 
  column_to_rownames(var = "sglt_anova_genera")

### Make a tibble for p-values
sglt_da_ttest_pvalues <- tibble(x = 1:5, y = 1:5, z = 1:5) %>% 
  add_column(sglt_anova_genera, .before = "x") %>% 
  dplyr::rename(II = x,
                III = y,
                IV = z) %>% 
  column_to_rownames(var = "sglt_anova_genera")

for (j in 1:length(sglt_anova_genera)){
  
  temp_data <- sglt_genera_comparisons %>% 
    filter(Genus == sglt_anova_genera[j])
  
  for (i in 1:length(timepoints)){
    df_sglt <- temp_data %>% 
      filter(Timepoint %in% c("I", timepoints[i])) %>% 
      group_by(PatientID) %>% 
      filter(n() != 1) %>% 
      arrange(Timepoint, PatientID) %>% 
      ungroup()
    
    da_test_sglt_raw <- t.test(clr ~ Timepoint, data = df_sglt, paired = TRUE) 
    
    da_results_sglt <- da_test_sglt_raw %>% 
      broom::tidy()
    
    sglt_da_ttest_estimates[sglt_anova_genera[j], timepoints[i]] <- da_results_sglt$estimate
    sglt_da_ttest_pvalues[sglt_anova_genera[j], timepoints[i]] <- da_results_sglt$p.value
    
  }
}

### Correct p-values for multiple testing w/ Benjamini-Hochberg method
sglt_da_estim <- sglt_da_ttest_estimates %>% 
  rownames_to_column(var = "Genus_SGLT") %>% 
  pivot_longer(cols = 2:4, names_to = "timepoint_SGLT", values_to = "estimate_SGLT")


sglt_da_pvalues <- sglt_da_ttest_pvalues %>% 
  rownames_to_column(var = "Genus") %>% 
  pivot_longer(cols = 2:4, names_to = "timepoint", values_to = "p_value_SGLT") %>% 
  select(-c(Genus, timepoint)) %>% 
  mutate(p_value_BH_SGLT = p.adjust(p_value_SGLT, method = "BH"))

sglt_da_ttest_BH <- bind_cols(sglt_da_estim, sglt_da_pvalues)

# ___________________________________________________________________________ #

# Summary results

anova_da_results <- bind_cols(glp_da_ttest_BH, sglt_da_ttest_BH)

t_test_da_results <- bind_cols(glp_da_ttest_BH, sglt_da_ttest_BH)

# ___________________________________________________________________________ #

# Visualize results for significant genera after ANOVA

# paletteer::paletteer_d("calecopal::superbloom1")
# 
# #B9C7E2FF #ECAB99FF #F1C100FF #5B6530FF #9484B1FF 

# Add a medication facet label
glp_med_string <- replicate(4480, "GLP-1-RA")

# Add this string as a column to a data frame
glp_genera_comparisons <- glp_genera_comparisons %>% 
  mutate(Medication = glp_med_string)

glp_da_plot_rb <- glp_genera_comparisons %>% 
  filter(Genus == "Romboutsia") %>% 
  ggplot(aes(x = Timepoint, y = clr, fill = Genus)) +
  scale_fill_manual(values = "#B9C7E2FF") +
  guides(fill = "none") +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  facet_grid(Medication~Genus, scales = "free", switch = "y") +
  ggsignif::geom_signif(
    y_position = c(6.2, 7.0, 7.8), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.67", "0.59", "0.11"), tip_length = 0.02) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "", 
       y = "CLR") +
  theme(strip.text.x = element_text(size = 12, face = "italic"),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12))

glp_da_plot_im <- glp_genera_comparisons %>% 
  filter(Genus == "Intestinimonas") %>%
  ggplot(aes(x = Timepoint, y = clr, fill = Genus)) +
  scale_fill_manual(values = "#ECAB99FF") +
  guides(fill = "none") +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  facet_wrap(~Genus, scales = "free") +
  ggsignif::geom_signif(
    y_position = c(6, 7, 8), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.77", "0.9", "0.55"), tip_length = 0.02) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "",
       y = "") +
  theme(strip.text = element_text(size = 12, face = "italic"),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))

glp_da_plot_hf <- glp_genera_comparisons %>% 
  filter(Genus == "Haemophilus") %>% 
  ggplot(aes(x = Timepoint, y = clr, fill = Genus)) +
  scale_fill_manual(values = "#F1C100FF") +
  guides(fill = "none") +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  facet_wrap(~Genus, scales = "free") +
  ggsignif::geom_signif(
    y_position = c(35, 39, 43), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.17", "0.16", "0.34"), tip_length = 0.02) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "",
       y = "") +
  theme(strip.text = element_text(size = 12, face = "italic"),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))

glp_da_plot_mb <- glp_genera_comparisons %>% 
  filter(Genus == "Marvinbryantia") %>% 
  ggplot(aes(x = Timepoint, y = clr, fill = Genus)) +
  scale_fill_manual(values = "#5B6530FF") +
  guides(fill = "none") +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  facet_wrap(~Genus, scales = "free") +
  ggsignif::geom_signif(
    y_position = c(8, 9, 10), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.2", "0.16", "0.67"), tip_length = 0.02) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "",
       y = "") +
  theme(strip.text = element_text(size = 12, face = "italic"),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))

glp_da_plot_rt <- glp_genera_comparisons %>% 
  filter(Genus == "Rothia") %>% 
  ggplot(aes(x = Timepoint, y = clr, fill = Genus)) +
  scale_fill_manual(values = "#9484B1FF") +
  guides(fill = "none") +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  facet_wrap(~Genus, scales = "free") +
  ggsignif::geom_signif(
    y_position = c(13, 15, 17), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.03", "0.36", "0.86"), tip_length = 0.02) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "",
       y = "") +
  theme(strip.text = element_text(size = 12, face = "italic"),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))

glp_da <- ggpubr::ggarrange(glp_da_plot_rb, glp_da_plot_im, 
                            glp_da_plot_hf, 
                            glp_da_plot_mb, glp_da_plot_rt,
                            ncol = 5, nrow = 1)

# SGLT-2

# paletteer::paletteer_d("calecopal::bixby")
# #286A81FF #045CB4FF #7F6F43FF #748B75FF #B8B196FF

# Add a medication facet label
sglt_med_string <- replicate(14896, "SGLT-2")

# Add this string as a column to a data frame
sglt_genera_comparisons <- sglt_genera_comparisons %>% 
  mutate(Medication = sglt_med_string)

sglt_da_plot_as <- sglt_genera_comparisons %>% 
  filter(Genus == "Alistipes") %>%
  ggplot(aes(x = Timepoint, y = clr, fill = Genus)) +
  scale_fill_manual(values = "#045CB4FF") +
  guides(fill = "none") +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  facet_grid(Medication~Genus, scales = "free", switch = "y") +
  ggsignif::geom_signif(
    y_position = c(40, 45, 50), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.18", "0.19", "0.97"), tip_length = 0.02) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "",
       y = "CLR") +
  guides(fill = "none") +
  theme(strip.text.x = element_text(size = 12, face = "italic"),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))

sglt_da_plot_cc <- sglt_genera_comparisons %>% 
  filter(Genus == "Coprococcus 2") %>% 
  ggplot(aes(x = Timepoint, y = clr, fill = Genus)) +
  scale_fill_manual(values = "#286A81FF") +
  guides(fill = "none") +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  facet_wrap(~Genus, scales = "free") +
  ggsignif::geom_signif(
    y_position = c(16, 18, 20), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("1", "0.72", "0.04"), tip_length = 0.02) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "", 
       y = "") +
  guides(fill = "none") +
  theme(strip.text = element_text(size = 12, face = "italic"), 
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))

sglt_da_plot_im <- sglt_genera_comparisons %>% 
  filter(Genus == "Intestinimonas") %>% 
  ggplot(aes(x = Timepoint, y = clr, fill = Genus)) +
  scale_fill_manual(values = "#ECAB99FF") +
  guides(fill = "none") + 
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  facet_wrap(~Genus, scales = "free") +
  ggsignif::geom_signif(
    y_position = c(13, 15, 17), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("1", "0.24", "0.03"), tip_length = 0.02) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(y = "") +
  guides(fill = "none") +
  theme(strip.text = element_text(size = 12, face = "italic"),
        axis.text.x = element_text(size = 12), 
        axis.title.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))


sglt_da_plot_tb <- sglt_genera_comparisons %>% 
  filter(Genus == "Turicibacter") %>% 
  ggplot(aes(x = Timepoint, y = clr, fill = Genus)) +
  scale_fill_manual(values = "#748B75FF") +
  guides(fill = "none") +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  facet_wrap(~Genus, scales = "free") + 
  ggsignif::geom_signif(
    y_position = c(17, 20, 23), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.46", "0.35", "0.08"), tip_length = 0.02) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "",
       y = "") +
  guides(fill = "none") +
  theme(strip.text = element_text(size = 12, face = "italic"), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))


sglt_da_plot_dt <- sglt_genera_comparisons %>% 
  filter(Genus == "DTU089") %>% 
  ggplot(aes(x = Timepoint, y = clr, fill = Genus)) +
  scale_fill_manual(values = "#B8B196FF") +
  guides(fill = "none") +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  facet_wrap(~Genus, scales = "free") +
  ggsignif::geom_signif(
    y_position = c(8, 9, 10), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.4", "0.05", "0.21"), tip_length = 0.02) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "",
       y = "") +
  guides(fill = "none") +
  theme(strip.text = element_text(size = 12), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))

sglt_da <- ggpubr::ggarrange(sglt_da_plot_as, sglt_da_plot_cc, 
                             sglt_da_plot_im, sglt_da_plot_tb,
                             sglt_da_plot_dt,
                             ncol = 5, nrow = 1)

combined_da_plot <- ggpubr::ggarrange(glp_da, sglt_da, nrow = 2, ncol = 1)

ggsave("combined_da_plot.svg", device = "svg", width = 18, height = 11)

