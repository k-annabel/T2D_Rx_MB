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

# Alpha diversity analysis with repeated measures ANOVA and t-tests

## GLP-1-RA

### Calculate alpha diversity
tse_glp <- estimateDiversity(tse_glp, assay.type = "relabundance", index = "shannon")
tse_glp <- estimateEvenness(tse_glp, assay.type = "relabundance", index = "pielou")
tse_glp <- estimateRichness(tse_glp, assay.type = "relabundance", index = "observed")

## Prepare data

### Extract metadata with alpha diversity indices
tse_glp_metadata <- as.data.frame(colData(tse_glp))

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
  
  glp_temp_alpha <- tse_glp_metadata %>% 
    rownames_to_column(var = "SampleID") %>% 
    select(c(2:4, 67:69)) %>% 
    pivot_longer(cols = 4:6, names_to = "alpha", values_to = "value") %>% 
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
  
  temp_alpha_glp <- tse_glp_metadata %>% 
    rownames_to_column(var = "SampleID") %>% 
    select(c(2:4, 67:69)) %>% 
    pivot_longer(cols = 4:6, names_to = "alpha_div", values_to = "value") %>% 
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

### Extract metadata with alpha diversity indices
tse_sglt_metadata <- as.data.frame(colData(tse_sglt))

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
  
  sglt_temp_alpha <- tse_sglt_metadata %>% 
    rownames_to_column(var = "SampleID") %>% 
    select(c(2:4, 67:69)) %>% 
    pivot_longer(cols = 4:6, names_to = "alpha", values_to = "value") %>% 
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
  
  temp_alpha_sglt <- tse_sglt_metadata %>% 
    rownames_to_column(var = "SampleID") %>% 
    select(c(2:4, 67:69)) %>% 
    pivot_longer(cols = 4:6, names_to = "alpha_div", values_to = "value") %>% 
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

# writexl::write_xlsx(anova_alpha_results, "anova_alpha_results.xlsx")
# writexl::write_xlsx(t_test_alpha_results, "t_test_alpha_results.xlsx")

# ___________________________________________________________________________ #

# Visualize the results

## GLP-1-RA

# paletteer_d("ggprism::viridis")

  #440154FF #414487FF #2A788EFF #22A884FF #7AD151FF #FDE725FF

glp_alpha_data <- tse_glp_metadata %>% 
  rownames_to_column(var = "SampleID") %>% 
  select(c(2:4, 67:69)) %>% 
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
       y = "") +
  guides(fill = "none") +
  ggsignif::geom_signif(
    y_position = c(150, 160, 170), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.62", "0.36", "0.09"), tip_length = 0.02) +
  theme_classic() +
  theme(strip.text.x = element_text(size = 16), 
        strip.text.y = element_text(size = 16), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))


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
    y_position = c(4.0, 4.2, 4.4), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.79", "0.29", "0.7"), tip_length = 0.02) +
  theme_classic() +
  theme(strip.text.x = element_text(size = 16), 
        strip.text.y = element_text(size = 16), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))

glp_alpha_pielou <- glp_alpha_data %>% 
  filter(alpha_div == "Pielou") %>% 
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
    y_position = c(0.85, 0.9, 0.95), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.93", "0.28", "0.9"), tip_length = 0.02) +
  theme_classic() +
  theme(strip.text.x = element_text(size = 16), 
        strip.text.y = element_text(size = 16), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))

glp_alpha_comb <- ggpubr::ggarrange(glp_alpha_observed, glp_alpha_shannon, glp_alpha_pielou, 
                                    nrow = 1, ncol = 3)

## SGLT-2

sglt_alpha_data <- tse_sglt_metadata %>% 
  rownames_to_column(var = "SampleID") %>% 
  select(c(2:4, 67:69)) %>% 
  pivot_longer(cols = 4:6, names_to = "alpha_div", values_to = "value") %>% 
  mutate(alpha_div = recode(alpha_div, 
                            "shannon" = "Shannon",
                            "pielou" = "Pielou",
                            "observed" = "Observed")) %>% 
  mutate(Medication = recode(Medication, 
                             "SGLT-2" = "SGLT-2i"))
  

sglt_alpha_observed <- sglt_alpha_data %>% 
  filter(alpha_div == "Observed") %>% 
  ggplot(aes(x = Timepoint, y = value, fill = alpha_div)) +
  scale_fill_manual(values = "#669ABFFF") +
  facet_grid(Medication~alpha_div, switch = "y", scales = "free") +
  geom_boxplot() +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.3) +
  geom_point(position = position_jitter(width = 0.05)) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "",
       y = "") +
  guides(fill = "none") +
  ggsignif::geom_signif(
    y_position = c(150, 160, 170), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.35", "0.05", "0.3"), tip_length = 0.02) +
  theme_classic() +
  theme(strip.text.x = element_text(size = 16), 
        strip.text.y = element_text(size = 16), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))

sglt_alpha_shannon <- sglt_alpha_data %>% 
  filter(alpha_div == "Shannon") %>% 
  ggplot(aes(x = Timepoint, y = value, fill = alpha_div)) +
  scale_fill_manual(values = "#669ABFFF") +
  facet_grid(.~alpha_div, switch = "y", scales = "free") +
  geom_boxplot() +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  geom_point(position = position_jitter(width = 0.02)) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "",
       y = "") +
  guides(fill = "none") +
  ggsignif::geom_signif(
    y_position = c(4.0, 4.2, 4.4), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.61", "0.8", "0.21"), tip_length = 0.02) +
  theme_classic() +
  theme(strip.text.x = element_text(size = 16), 
        strip.text.y = element_text(size = 16), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))

sglt_alpha_pielou <- sglt_alpha_data %>% 
  filter(alpha_div == "Pielou") %>% 
  ggplot(aes(x = Timepoint, y = value, fill = alpha_div)) +
  scale_fill_manual(values = "#669ABFFF") +
  facet_grid(.~alpha_div, switch = "y", scales = "free") +
  geom_boxplot() +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.2) +
  geom_point(position = position_jitter(width = 0.02)) +
  scale_x_discrete(labels = c("BL", "M1", "M3", "M12")) +
  labs(x = "",
       y = "") +
  guides(fill = "none") +
  ggsignif::geom_signif(
    y_position = c(0.8, 0.84, 0.88), xmin = c(1, 1, 1), xmax = c(2, 3, 4),
    annotation = c("0.47", "0.91", "0.12"), tip_length = 0.02) +
  theme_classic() +
  theme(strip.text.x = element_text(size = 16), 
        strip.text.y = element_text(size = 16), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))

sglt_alpha_comb <- ggpubr::ggarrange(sglt_alpha_observed, sglt_alpha_shannon, sglt_alpha_pielou,
                                     nrow = 1, ncol = 3)

alpha_comb <- ggpubr::ggarrange(glp_alpha_comb, sglt_alpha_comb, 
                                nrow = 2, ncol = 1)

ggsave("combined_alpha_plot.svg", dpi = 300, width = 11, height = 11)
