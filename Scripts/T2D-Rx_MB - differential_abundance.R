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

# Collapse into Genus level
tse_genus <- agglomerateByRank(tse, rank = "Genus")

# Convert counts into relative abundances
tse_genus <- transformAssay(tse_genus, assay.type = "counts", method = "relabundance")

# Convert relative abundances into CLR-transformed values
mat <- assay(tse_genus, "relabundance")
tse_genus <- transformAssay(tse_genus, assay.type = "relabundance", method = "clr", 
                      pseudocount = min(mat[mat>0]))

# Remove "Genus:" label
rownames(tse_genus) <- sub("Genus:", "", rownames(tse_genus))

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

# Create tibbles for results
glp_da_anova_ges <- tibble(x = 1:128) %>% 
  add_column(glp_top_taxa, .before = "x") %>% 
  dplyr::rename(ges = x)

glp_da_anova_pvalues <- tibble(x = 1:128) %>% 
  add_column(glp_top_taxa, .before = "x") %>% 
  dplyr::rename(p_value = x)

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
  
  glp_da_anova_ges$ges[j] <- model_glp$ANOVA$ges
  glp_da_anova_pvalues$p_value[j] <- model_glp$ANOVA$p
  
  glp_da_anova_results_raw <- cbind(glp_da_anova_ges, glp_da_anova_pvalues)
  
}

### Correct p-values for multiple testing w/ Benjamini-Hochberg method

glp_da_anova_results <- glp_da_anova_results_raw[-3]

glp_anova_da_results_BH <- glp_da_anova_results %>% 
  mutate(p_value_BH = p.adjust(p_value, method = "BH")) %>% 
  filter(p_value <= 0.053)

# ___________________________________________________________________________ #

# Perform t-tests for five significant genera in ANOVA results

### Initialize empty vector for results
timepoints <- c("II", "III", "IV")

### Pull five significant genera
glp_anova_genera <- glp_anova_da_results_BH %>% 
  filter(p_value <= 0.053) %>% 
  pull(glp_top_taxa)

### Make a tibble for estimates
glp_da_ttest_estimates <- tibble(x = 1:7, y = 1:7, z = 1:7) %>% 
  add_column(glp_anova_genera, .before = "x") %>% 
  dplyr::rename(II = x,
                III = y,
                IV = z) %>% 
  column_to_rownames(var = "glp_anova_genera")

### Make a tibble for p-values
glp_da_ttest_pvalues <- tibble(x = 1:7, y = 1:7, z = 1:7) %>% 
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
  select(-c(Genus, timepoint))# %>% 
  #mutate(p_value_BH_GLP = p.adjust(p_value_GLP, method = "BH"))

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
sglt_da_anova_ges <- tibble(x = 1:147) %>% 
  add_column(sglt_top_taxa, .before = "x") %>% 
  dplyr::rename(ges = x)

sglt_da_anova_pvalues <- tibble(x = 1:147) %>% 
  add_column(sglt_top_taxa, .before = "x") %>% 
  dplyr::rename(p_value = x)

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
  
  sglt_da_anova_ges$ges[j] <- model_sglt$ANOVA$ges
  sglt_da_anova_pvalues$p_value[j] <- model_sglt$ANOVA$p
  
  sglt_da_anova_results_raw <- cbind(sglt_da_anova_ges, sglt_da_anova_pvalues)
  
}

### Correct p-values for multiple testing w/ Benjamini-Hochberg method

sglt_da_anova_results <- sglt_da_anova_results_raw[-3]

sglt_anova_da_results_BH <- sglt_da_anova_results %>% 
  mutate(p_value_BH = p.adjust(p_value, method = "BH"))

# Perform t-tests

# Initialize empty vector for results
timepoints <- c("II", "III", "IV")

# Pull five significant genera (there are more... (N = 10))
sglt_anova_genera <- sglt_anova_da_results_BH %>% 
  filter(p_value <= 0.051) %>% 
  pull(sglt_top_taxa)

### Make a tibble for estimates
sglt_da_ttest_estimates <- tibble(x = 1:10, y = 1:10, z = 1:10) %>% 
  add_column(sglt_anova_genera, .before = "x") %>% 
  dplyr::rename(II = x,
                III = y,
                IV = z) %>% 
  column_to_rownames(var = "sglt_anova_genera")

### Make a tibble for p-values
sglt_da_ttest_pvalues <- tibble(x = 1:10, y = 1:10, z = 1:10) %>% 
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
  select(-c(Genus, timepoint))
  #mutate(p_value_BH_SGLT = p.adjust(p_value_SGLT, method = "BH"))

sglt_da_ttest_BH <- bind_cols(sglt_da_estim, sglt_da_pvalues)

# ___________________________________________________________________________ #

# Summary results

anova_da_results <- bind_cols(glp_anova_da_results_BH, sglt_anova_da_results_BH)

t_test_da_results <- bind_cols(glp_da_ttest_BH, sglt_da_ttest_BH)

# writexl::write_xlsx(anova_da_results, "anova_da_results.xlsx")
# writexl::write_xlsx(t_test_da_results, "t_test_da_results.xlsx")

# ___________________________________________________________________________ #

# Visualize results for significant genera after ANOVA

# paletteer::paletteer_d("ltc::minou")
# 
# #00798CFF #D1495BFF #EDAE49FF #66A182FF #2E4057FF #8D96A3FF

# Add a medication facet label
glp_med_string <- replicate(4480, "GLP-1-RA")

# Add this string as a column to a data frame
glp_genera_comparisons <- glp_genera_comparisons %>% 
  mutate(Medication = glp_med_string)

glp_da_plot_am <- glp_genera_comparisons %>% 
  filter(Genus == "Actinomyces") %>%
  ggplot(aes(x = Timepoint, y = clr, fill = Genus)) +
  scale_fill_manual(values = "#00798CFF") +
  guides(fill = "none") +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  facet_wrap(~Genus, scales = "free") +
  ggsignif::geom_signif(
    y_position = c(7, 9, 11), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.71", "0.41", "0.38"), tip_length = 0.02) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "",
       y = "CLR") +
  theme(strip.text = element_text(size = 12, face = "italic"),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))

glp_da_plot_cs <- glp_genera_comparisons %>% 
  filter(Genus == "Candidatus Soleaferrea") %>%
  ggplot(aes(x = Timepoint, y = clr, fill = Genus)) +
  scale_fill_manual(values = "#D1495BFF") +
  guides(fill = "none") +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  facet_wrap(~Genus, scales = "free") +
  ggsignif::geom_signif(
    y_position = c(3, 4, 5), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.18", "0.18", "0.02"), tip_length = 0.02) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "",
       y = "") +
  theme(strip.text = element_text(size = 12, face = "italic"),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))

glp_da_plot_fb <- glp_genera_comparisons %>% 
  filter(Genus == "Fusobacterium") %>%
  ggplot(aes(x = Timepoint, y = clr, fill = Genus)) +
  scale_fill_manual(values = "#EDAE49FF") +
  guides(fill = "none") +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  facet_wrap(~Genus, scales = "free") +
  ggsignif::geom_signif(
    y_position = c(8, 10, 12), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.48", "0.25", "0.03"), tip_length = 0.02) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "",
       y = "") +
  theme(strip.text = element_text(size = 12, face = "italic"),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))

glp_da_plot_hp <- glp_genera_comparisons %>% 
  filter(Genus == "Haemophilus") %>%
  ggplot(aes(x = Timepoint, y = clr, fill = Genus)) +
  scale_fill_manual(values = "#66A182FF") +
  guides(fill = "none") +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  facet_wrap(~Genus, scales = "free") +
  ggsignif::geom_signif(
    y_position = c(7, 9, 11), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.19", "0.16", "0.08"), tip_length = 0.02) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "",
       y = "") +
  theme(strip.text = element_text(size = 12, face = "italic"),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))

glp_da_plot_pb <- glp_genera_comparisons %>% 
  filter(Genus == "Pseudobutyrivibrio") %>%
  ggplot(aes(x = Timepoint, y = clr, fill = Genus)) +
  scale_fill_manual(values = "#2E4057FF") +
  guides(fill = "none") +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  facet_wrap(~Genus, scales = "free") +
  ggsignif::geom_signif(
    y_position = c(6, 8, 10), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.30", "0.09", "0.07"), tip_length = 0.02) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "",
       y = "") +
  theme(strip.text = element_text(size = 12, face = "italic"),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))

glp_da_plot_vl <- glp_genera_comparisons %>% 
  filter(Genus == "Veillonella") %>%
  ggplot(aes(x = Timepoint, y = clr, fill = Genus)) +
  scale_fill_manual(values = "#8D96A3FF") +
  guides(fill = "none") +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  facet_wrap(~Genus, scales = "free") +
  ggsignif::geom_signif(
    y_position = c(5, 7, 9), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.33", "0.06", "0.02"), tip_length = 0.02) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "",
       y = "") +
  theme(strip.text = element_text(size = 12, face = "italic"),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))

glp_da <- ggpubr::ggarrange(glp_da_plot_am, glp_da_plot_cs, 
                            glp_da_plot_fb, glp_da_plot_hp, 
                            glp_da_plot_pb, glp_da_plot_vl, 
                            ncol = 6, nrow = 1)

# ___________________________________________________________________________ #

glp_da_plot_rb <- glp_genera_comparisons %>% 
  filter(Genus == "Romboutsia") %>% 
  ggplot(aes(x = Timepoint, y = clr, fill = Genus)) +
  scale_fill_manual(values = "#00798CFF") +
  guides(fill = "none") +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  facet_grid(Medication~Genus, scales = "free", switch = "y") +
  ggsignif::geom_signif(
    y_position = c(6.2, 7.0, 7.8), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.67", "0.57", "0.11"), tip_length = 0.02) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "", 
       y = "CLR") +
  theme(strip.text.x = element_text(size = 12, face = "italic"),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12))

glp_da_plot_sl <- glp_genera_comparisons %>% 
  filter(Genus == "Slackia") %>%
  #filter(PatientID != "GLP1RA-10") %>% 
  ggplot(aes(x = Timepoint, y = clr, fill = Genus)) +
  scale_fill_manual(values = "#D1495BFF") +
  guides(fill = "none") +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  facet_wrap(~Genus, scales = "free") +
  ggsignif::geom_signif(
    y_position = c(30, 34, 38), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.84", "0.99", "0.08"), tip_length = 0.02) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "",
       y = "") +
  theme(strip.text = element_text(size = 12, face = "italic"),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))

# SGLT-2

# paletteer::paletteer_d("Manu::Hoiho")
# #CABEE9FF #7C7189FF #FAE093FF #D04E59FF #BC8E7DFF #2F3D70FF 

# Add a medication facet label
sglt_med_string <- replicate(14896, "SGLT-2")

# Add this string as a column to a data frame
sglt_genera_comparisons <- sglt_genera_comparisons %>% 
  mutate(Medication = sglt_med_string)

sglt_da_plot_ab <- sglt_genera_comparisons %>% 
  filter(Genus == "Agathobacter") %>%
  ggplot(aes(x = Timepoint, y = clr, fill = Genus)) +
  scale_fill_manual(values = "#CABEE9FF") +
  guides(fill = "none") +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  facet_grid(Medication~Genus, scales = "free", switch = "y") +
  ggsignif::geom_signif(
    y_position = c(10, 12, 14), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.24", "0.52", "0.03"), tip_length = 0.02) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "",
       y = "CLR") +
  guides(fill = "none") +
  theme(strip.text.x = element_text(size = 12, face = "italic"),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))

sglt_da_plot_am <- sglt_genera_comparisons %>% 
  filter(Genus == "Akkermansia") %>% 
  ggplot(aes(x = Timepoint, y = clr, fill = Genus)) +
  scale_fill_manual(values = "#7C7189FF") +
  guides(fill = "none") +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  facet_wrap(~Genus, scales = "free") +
  ggsignif::geom_signif(
    y_position = c(10, 12, 14), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.96", "0.26", "0.27"), tip_length = 0.02) +
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
  scale_fill_manual(values = "#FAE093FF") +
  guides(fill = "none") +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  facet_wrap(~Genus, scales = "free") +
  ggsignif::geom_signif(
    y_position = c(3, 4, 5), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.44", "0.02", "0.16"), tip_length = 0.02) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "", 
       y = "") +
  guides(fill = "none") +
  theme(strip.text = element_text(size = 12, face = "italic"), 
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))

sglt_da_plot_ps <- sglt_genera_comparisons %>% 
  filter(Genus == "Parasutterella") %>% 
  ggplot(aes(x = Timepoint, y = clr, fill = Genus)) +
  scale_fill_manual(values = "#D04E59FF") +
  guides(fill = "none") + 
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  facet_wrap(~Genus, scales = "free") +
  ggsignif::geom_signif(
    y_position = c(8, 10, 12), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.92", "0.30", "0.08"), tip_length = 0.02) +
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
  scale_fill_manual(values = "#BC8E7DFF") +
  guides(fill = "none") +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  facet_wrap(~Genus, scales = "free") + 
  ggsignif::geom_signif(
    y_position = c(5, 7, 9), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.86", "0.71", "0.06"), tip_length = 0.02) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "",
       y = "") +
  guides(fill = "none") +
  theme(strip.text = element_text(size = 12, face = "italic"), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))


sglt_da_plot_uc <- sglt_genera_comparisons %>% 
  filter(Genus == "uncultured bacterium_8") %>% 
  ggplot(aes(x = Timepoint, y = clr, fill = Genus)) +
  scale_fill_manual(values = "#2F3D70FF") +
  guides(fill = "none") +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.02)) +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  facet_wrap(~Genus, scales = "free") +
  ggsignif::geom_signif(
    y_position = c(4, 5, 6), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.87", "0.04", "0.46"), tip_length = 0.02) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "",
       y = "") +
  guides(fill = "none") +
  theme(strip.text = element_text(size = 12), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))

sglt_da <- ggpubr::ggarrange(sglt_da_plot_ab, sglt_da_plot_am, sglt_da_plot_dt,
                             sglt_da_plot_ps, sglt_da_plot_tb,
                             sglt_da_plot_uc,
                             ncol = 6, nrow = 1)

combined_da_plot <- ggpubr::ggarrange(glp_da, sglt_da, nrow = 2, ncol = 1)

ggsave("combined_da_plot.svg", device = "svg", width = 18, height = 11)

