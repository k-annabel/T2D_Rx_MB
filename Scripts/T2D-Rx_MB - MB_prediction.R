# T2D-Rx_MB_analysis - MB prediction ability of T2D treatment endpoints

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

# GLP-1-RA

# Transform data for the correlation matrix

# Calculate alpha diversity
tse_glp <- estimateDiversity(tse_glp, assay.type = "relabundance", index = "shannon")
tse_glp <- estimateEvenness(tse_glp, assay.type = "relabundance", index = "pielou")
tse_glp <- estimateRichness(tse_glp, assay.type = "relabundance", index = "observed")

# Extract CLR-values for PCA
tse_glp_clr <- as.data.frame(assay(tse_glp, "clr"))

# Perform PCA 
pca_glp <- prcomp(t(tse_glp_clr))

# Pull out top genera based on CLR-values and prevalence (5 of 20 samples)
top_genera_glp <- mia::getTopFeatures(tse_glp,
                                      assay_name = "clr",
                                      top = 5)

# Extract CLR-values for the heatmap
heatmap_glp_clr_raw <- as.data.frame(t(assay(tse_glp, "clr"))) 

# Filter out less prevalent genera
heatmap_glp_clr <- heatmap_glp_clr_raw %>% 
  select(all_of(top_genera_glp)) %>% 
  rownames_to_column(var = "SampleID")

# Extract metadata 
heatmap_metadata_glp_raw <- as.data.frame(colData(tse_glp)) %>% 
  rownames_to_column(var = "SampleID")

# Add PCA data
heatmap_metadata_glp <- data.frame(pca_glp$x[ , 1:5], heatmap_metadata_glp_raw)

# Merge two data frames
heatmap_glp_raw <- merge(heatmap_glp_clr, heatmap_metadata_glp)

# Subset the data frame
glp_corr_data <- heatmap_glp_raw %>% 
  select(c(1:12, 13, 21, 25, 63:66, 77:79)) %>% 
  relocate(PatientID, .before = Blautia) %>% 
  relocate(Timepoint, .before = Blautia) %>% 
  relocate(c("shannon", "pielou", "observed"), .before = PC1) %>% 
  mutate(across(.cols = 4:22, .fns=as.numeric)) %>% 
  pivot_longer(cols = 4:22, names_to = "Parameter", values_to = "Value")

# Pull out BL-data
glp_corr_BL <- glp_corr_data %>% 
  filter(Timepoint == "I") %>% 
  select(-Timepoint) %>% 
  dplyr::rename(Value_BL = Value)

# Pull out other timepoints' data
glp_corr_TP <- glp_corr_data %>% 
  filter(!Timepoint %in% c("I"))

# Merge two data frames together
glp_heatmap_comb <- merge(glp_corr_BL, glp_corr_TP, 
                          by = c("PatientID", "Parameter"),
                          all = TRUE)

# Tidy the merged data frame
glp_heatmap <- glp_heatmap_comb %>% 
  arrange(PatientID, Parameter, SampleID.y) %>% 
  mutate(chg = Value - Value_BL) %>% 
  select(-Value) %>% 
  pivot_wider(names_from = Parameter, values_from = c(Value_BL, chg)) %>% 
  filter(!Timepoint %in% c("II", "IV")) %>% 
  select(-c(PatientID, SampleID.x, SampleID.y, Timepoint))
# dim(glp_heatmap) - 10 x 38

# Calculate the correlation matrix
glp_corr <- round(cor(glp_heatmap, use = "complete.obs"), 3)
#glp_corr <- round(cor(glp_heatmap), 1)
glp_corr[is.na(glp_corr)] <- 0

# Compute a matrix of correlation p-values
p.mat_glp <- cor_pmat(glp_heatmap)
p.mat_glp[is.na(p.mat_glp)] <- 0
  

# Visualize the inital correlation matrix
pheatmap(glp_corr)

# Make a submatrix for analysis (all TPs included)
glp_corr_alam <- glp_corr[c(1, 3:5, 10:15, 17:19), c(21, 25:28, 35)]
glp_pmat_alam <- p.mat_glp[c(1, 3:5, 10:15, 17:19), c(21, 25:28, 35)]
glp_labels_row <- c("Alistipes", "Dorea", "Roseburia", "Bacteroides", "Blautia")
glp_labels_col <- c("Blood glucose", "BMI", "HbA1c (%)")
glp_alam_hm <- pheatmap(glp_corr_alam, 
                        #labels_row = glp_labels_row, 
                        #labels_col = glp_labels_col, 
                        fontsize = 18, angle_col = 45)


corrplot::corrplot(glp_corr_alam,
                   method = "square")


# Option 2 with ggplot2

# Select correct columns and rows
glp_corr_plot <- glp_corr %>% 
  as.data.frame() %>% 
  select(c(21, 24, 26:28, 35)) %>% 
  slice(c(1, 3:4, 6, 10:15, 17:19)) %>% 
  rownames_to_column(var = "BL_Parameter") %>% 
  pivot_longer(cols = 2:7, names_to = "chg_Parameter", values_to = "Value")

glp_corr_plot_p <- p.mat_glp %>% 
  as.data.frame() %>% 
  column_to_rownames(var = "rowname") %>% 
  select(c(21, 24, 26:28, 35)) %>% 
  slice(c(1, 3:4, 6, 10:15, 17:19)) %>% 
  rownames_to_column(var = "BL_Parameter") %>% 
  pivot_longer(cols = 2:7, names_to = "chg_Parameter", values_to = "p_value")

glp_corr_plot <- glp_corr_plot %>% 
  mutate(p_value = glp_corr_plot_p$p_value) %>% 
  mutate(p_value_BH = p.adjust(p_value, method = "BH")) %>% 
  mutate(sig_p = ifelse(p_value_BH < 0.05, T, F)) %>% 
  mutate(p_if_sig_BH = ifelse(sig_p, p_value_BH, NA))

ggplot(glp_corr_plot, aes(x = BL_Parameter, y = chg_Parameter, fill = Value)) +
  geom_tile(colour = "black") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
  coord_fixed() +
  #geom_text(aes(label = round(p_if_sig_BH,3)), colour = "white", size = 3) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  annotate("point", x = 0.85, y = 4, shape = "*", size = 7, colour = "black") +
  annotate("point", x = 1.2, y = 4, shape = "*", size = 7, colour = "black") +
  annotate("point", x = 1, y = 5, shape = "*", size = 7, colour = "black") +
  annotate("point", x = 1, y = 6, shape = "*", size = 7, colour = "black") +
  scale_x_discrete(labels = c("Alistipes", "Blautia", "Dorea", "Lachnoclostridium",
                   "Observed", "PC1", "PC2", "PC3", "PC4", "PC5", "Pielou", "Roseburia", "Shannon")) +
  scale_y_discrete(labels = c("BMI", "HbA1c", "Lymphocytes", "Neutrophils", "NLR", "White blood cells")) +
  labs(x = "Baseline MB parameter", 
       y = NULL, 
       fill = "Pearson correlation")

ggsave("glp_heatmap.svg", device = "svg", width = 8, height = 10)
  


# SGLT-2

# Transform data for the correlation matrix

# Calculate alpha diversity
tse_sglt <- estimateDiversity(tse_sglt, assay.type = "relabundance", index = "shannon")
tse_sglt <- estimateEvenness(tse_sglt, assay.type = "relabundance", index = "pielou")
tse_sglt <- estimateRichness(tse_sglt, assay.type = "relabundance", index = "observed")

# Extract CLR-values for PCA
tse_sglt_clr <- as.data.frame(assay(tse_sglt, "clr"))

# Perform PCA 
pca_sglt <- prcomp(t(tse_sglt_clr))

# Pull out top genera based on CLR-values and prevalence (5 of 20 samples)
top_genera_sglt <- mia::getTopFeatures(tse_sglt,
                                      assay_name = "clr",
                                      top = 5)

# Extract CLR-values for the heatmap
heatmap_sglt_clr_raw <- as.data.frame(t(assay(tse_sglt, "clr"))) 

# Filter out less prevalent genera
heatmap_sglt_clr <- heatmap_sglt_clr_raw %>% 
  select(all_of(top_genera_sglt)) %>% 
  rownames_to_column(var = "SampleID")

# Convert colData into a data frame
heatmap_metadata_sglt_raw <- as.data.frame(colData(tse_sglt)) %>% 
  rownames_to_column(var = "SampleID")

# Add PCA data
heatmap_metadata_sglt <- data.frame(pca_sglt$x[ , 1:5], heatmap_metadata_sglt_raw)

# Merge two data frames
heatmap_sglt_raw <- merge(heatmap_sglt_clr, heatmap_metadata_sglt)

# Subset the data frame
sglt_corr_data <- heatmap_sglt_raw %>% 
  relocate(PatientID, .before = Blautia) %>% 
  relocate(Timepoint, .before = Blautia) %>% 
  select(c(1:13, 21, 25, 39, 77:79)) %>% 
  relocate(c("shannon", "pielou", "observed"), .before = PC1) %>% 
  mutate(across(.cols = 4:19, .fns=as.numeric)) %>% 
  pivot_longer(cols = 4:19, names_to = "Parameter", values_to = "Value")

# Pull out BL-data
sglt_corr_BL <- sglt_corr_data %>% 
  filter(Timepoint == "I") %>% 
  select(-Timepoint) %>% 
  dplyr::rename(Value_BL = Value)

# Pull out other timepoints' data
sglt_corr_TP <- sglt_corr_data %>% 
  filter(!Timepoint %in% c("I"))

# Merge two data frames together
sglt_heatmap_comb <- merge(sglt_corr_BL, sglt_corr_TP, 
                          by = c("PatientID", "Parameter"),
                          all = TRUE)

# Tidy the merged data frame
sglt_heatmap <- sglt_heatmap_comb %>% 
  arrange(PatientID, Parameter, SampleID.y) %>% 
  mutate(chg = Value - Value_BL) %>% 
  select(-Value) %>% 
  pivot_wider(names_from = Parameter, values_from = c(Value_BL, chg)) %>% 
  filter(!Timepoint %in% c("II", "IV")) %>% 
  select(-c(PatientID, SampleID.x, SampleID.y, Timepoint))
# dim(sglt_heatmap) - 10 x 32

#### Calculate the correlation matrix
sglt_corr <- round(cor(sglt_heatmap, use = "complete.obs"), 3)
#glp_corr <- round(cor(glp_heatmap), 1)
sglt_corr[is.na(sglt_corr)] <- 0

#### Compute a matrix of correlation p-values
p.mat_sglt <- cor_pmat(sglt_heatmap)
p.mat_sglt[is.na(p.mat_sglt)] <- 0

#### Visualize the inital correlation matrix
pheatmap(sglt_corr)

# Select correct columns and rows
sglt_corr_plot <- sglt_corr %>% 
  as.data.frame() %>% 
  select(c(17, 20:21)) %>% 
  slice(c(2:3, 6:16)) %>% 
  rownames_to_column(var = "BL_Parameter") %>% 
  pivot_longer(cols = 2:4, names_to = "chg_Parameter", values_to = "Value")

sglt_corr_plot_p <- p.mat_sglt %>% 
  as.data.frame() %>% 
  column_to_rownames(var = "rowname") %>% 
  select(c(17, 20:21)) %>% 
  slice(c(2:3, 6:16))%>% 
  rownames_to_column(var = "BL_Parameter") %>% 
  pivot_longer(cols = 2:4, names_to = "chg_Parameter", values_to = "p_value")

sglt_corr_plot <- sglt_corr_plot %>% 
  mutate(p_value = sglt_corr_plot_p$p_value) %>% 
  mutate(p_value_BH = p.adjust(p_value, method = "BH")) %>% 
  mutate(sig_p = ifelse(p_value_BH < 0.05, T, F)) %>% 
  mutate(p_if_sig_BH = ifelse(sig_p, p_value_BH, NA))

ggplot(sglt_corr_plot, aes(x = BL_Parameter, y = chg_Parameter, fill = Value)) +
  geom_tile(colour = "black") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
  coord_fixed() +
  geom_text(aes(label = round(p_if_sig_BH,3)), colour = "white", size = 3) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  scale_x_discrete(labels = c("Blautia", "Christensenellaceae R-7 group", "Observed",
                              "PC1", "PC2", "PC3", "PC4", "PC5", "Pielou", "Roseburia", 
                              "Ruminococcaceae UCG-002", "Ruminococcaceae UCG-005", "Shannon")) +
  scale_y_discrete(labels = c("BMI", "GFR", "HbA1c")) +
  labs(x = "Baseline MB parameter", 
       y = NULL, 
       fill = "Pearson correlation")

ggsave("sglt_heatmap.svg", device = "svg", width = 11, height = 8)
