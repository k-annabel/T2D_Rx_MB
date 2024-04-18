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
  relocate(PatientID, .before = Bacteroides) %>% 
  relocate(Timepoint, .before = Bacteroides) %>% 
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

# Option 2 with ggplot2

# Select correct columns and rows
glp_corr_plot <- glp_corr %>% 
  as.data.frame() %>% 
  dplyr::select(c(21, 26:29, 35)) %>% 
  dplyr::slice(c(1, 3:6, 11:15, 17:19)) %>% 
  rownames_to_column(var = "BL_Parameter") %>% 
  pivot_longer(cols = 2:7, names_to = "chg_Parameter", values_to = "Value")

glp_corr_plot_p <- p.mat_glp %>% 
  as.data.frame() %>% 
  column_to_rownames(var = "rowname") %>% 
  dplyr::select(c(21, 26:29, 35)) %>% 
  dplyr::slice(c(1, 3:6, 11:15, 17:19)) %>% 
  rownames_to_column(var = "BL_Parameter") %>% 
  pivot_longer(cols = 2:7, names_to = "chg_Parameter", values_to = "p_value")

glp_corr_plot <- glp_corr_plot %>% 
  mutate(p_value = glp_corr_plot_p$p_value) %>% 
  mutate(p_value_BH = p.adjust(p_value, method = "BH")) %>% 
  mutate(sig_p = ifelse(p_value_BH < 0.05, T, F)) %>% 
  mutate(p_if_sig_BH = ifelse(sig_p, p_value_BH, NA))

glp_heatmap <- ggplot(glp_corr_plot, aes(x = BL_Parameter, y = chg_Parameter, fill = Value)) +
  geom_tile(colour = "black") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
  coord_fixed() +
#  geom_text(aes(label = round(p_if_sig_BH,3)), colour = "white", size = 3) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  annotate("point", x = 0.85, y = 4, shape = "*", size = 7, colour = "black") +
  annotate("point", x = 1.2, y = 4, shape = "*", size = 7, colour = "black") +
  annotate("point", x = 1, y = 5, shape = "*", size = 7, colour = "black") +
#  annotate("point", x = 1, y = 6, shape = "*", size = 7, colour = "black") +
  scale_x_discrete(labels = c("Alistipes", "Bacteroides", "Blautia", "Collinsella", "Faecalibacterium",
                   "Observed", "PC1", "PC2", "PC3", "PC4", "PC5", "Pielou", "Shannon")) +
  scale_y_discrete(labels = c("BMI", "HbA1c", "Lymphocytes", "Neutrophils", "NLR", "White blood cells")) +
  labs(x = "Baseline MB parameter", 
       y = "Change in the parameter value", 
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
  relocate(PatientID, .before = Bacteroides) %>% 
  relocate(Timepoint, .before = Bacteroides) %>% 
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
  select(c(18, 23:24)) %>% 
  slice(c(1, 3:6, 9:16)) %>% 
  rownames_to_column(var = "BL_Parameter") %>% 
  pivot_longer(cols = 2:4, names_to = "chg_Parameter", values_to = "Value")

sglt_corr_plot_p <- p.mat_sglt %>% 
  as.data.frame() %>% 
  column_to_rownames(var = "rowname") %>% 
  select(c(18, 23:24)) %>% 
  slice(c(1, 3:6, 9:16))%>% 
  rownames_to_column(var = "BL_Parameter") %>% 
  pivot_longer(cols = 2:4, names_to = "chg_Parameter", values_to = "p_value")

sglt_corr_plot <- sglt_corr_plot %>% 
  mutate(p_value = sglt_corr_plot_p$p_value) %>% 
  mutate(p_value_BH = p.adjust(p_value, method = "BH")) %>% 
  mutate(sig_p = ifelse(p_value_BH < 0.05, T, F)) %>% 
  mutate(p_if_sig_BH = ifelse(sig_p, p_value_BH, NA))

sglt_heatmap <- ggplot(sglt_corr_plot, aes(x = BL_Parameter, y = chg_Parameter, fill = Value)) +
  geom_tile(colour = "black") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
  coord_fixed() +
  geom_text(aes(label = round(p_if_sig_BH,3)), colour = "white", size = 3) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  scale_x_discrete(labels = c("Alistipes", "Bacteroides", "Blautia", "Collinsella", 
                              "Faecalibacterium", "Observed", "PC1", "PC2", "PC3", 
                              "PC4", "PC5", "Pielou", "Shannon")) +
  scale_y_discrete(labels = c("BMI", "GFR", "HbA1c")) +
  labs(x = "Baseline MB parameter", 
       y = "Change in the parameter value", 
       fill = "Pearson correlation")

ggsave("sglt_heatmap.svg", device = "svg", width = 11, height = 8)

comb_prediction <- ggarrange(glp_heatmap, sglt_heatmap, ncol = 1, nrow = 2, common.legend = TRUE)

# ___________________________________________________________________________ #

# New analyses for all genera

# GLP-1-RA

# Leave in genera found differentially abundant by ANOVA
heatmap_glp_clr <- heatmap_glp_clr_raw %>% 
  select(all_of(glp_top_taxa)) %>% 
  rownames_to_column(var = "SampleID")

# Extract metadata 
heatmap_metadata_glp_raw <- as.data.frame(colData(tse_glp)) %>% 
  rownames_to_column(var = "SampleID")

# Add PCA data
heatmap_metadata_glp <- data.frame(pca_glp$x[ , 1:3], heatmap_metadata_glp_raw)

# Merge two data frames
heatmap_glp_raw_v2 <- merge(heatmap_glp_clr, heatmap_metadata_glp)

# Subset the data frame
glp_corr_data <- heatmap_glp_raw_v2 %>% 
  select(c(1:134, 142, 146, 184:187, 198:200)) %>% 
  relocate(PatientID, .before = "Prevotella 9") %>% 
  relocate(Timepoint, .before = "Prevotella 9") %>% 
  relocate(c("shannon", "pielou", "observed"), .before = PC1) %>% 
  mutate(across(.cols = 4:143, .fns=as.numeric)) %>% 
  pivot_longer(cols = 4:143, names_to = "Parameter", values_to = "Value")

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

# Calculate the correlation matrix
glp_corr <- round(cor(glp_heatmap, use = "complete.obs"), 3)
#glp_corr <- round(cor(glp_heatmap), 1)
#glp_corr[is.na(glp_corr)] <- 0

# Compute a matrix of correlation p-values
p.mat_glp <- cor_pmat(glp_heatmap)
#p.mat_glp[is.na(p.mat_glp)] <- 0

# Visualize the inital correlation matrix
pheatmap(glp_corr)

# Select correct columns and rows
glp_corr_plot <- glp_corr %>% 
  as.data.frame() %>% 
  select(c(151, 191, 208, 212, 214, 257)) %>% 
  slice(c(1:10, 12:50, 52:67, 69:71, 73, 75:116, 118:140)) %>% 
  rownames_to_column(var = "BL_Parameter") %>% 
  pivot_longer(cols = 2:7, names_to = "chg_Parameter", values_to = "Value")

glp_corr_plot_p <- p.mat_glp %>% 
  as.data.frame() %>% 
  column_to_rownames(var = "rowname") %>% 
  select(c(151, 191, 208, 212, 214, 257)) %>% 
  slice(c(1:10, 12:50, 52:67, 69:71, 73, 75:116, 118:140)) %>% 
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
  geom_text(aes(label = round(p_value,3)), colour = "black", size = 3) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

# Only stat signf association: BL_Alistipes ~ chg_Neutr

# SGLT-2

# Extract CLR-values of prevalrnt genera
heatmap_sglt_clr <- heatmap_sglt_clr_raw %>% 
  select(all_of(sglt_top_taxa)) %>% 
  rownames_to_column(var = "SampleID")

# Convert colData into a data frame
heatmap_metadata_sglt_raw <- as.data.frame(colData(tse_sglt)) %>% 
  rownames_to_column(var = "SampleID")

# Add PCA data
heatmap_metadata_sglt <- data.frame(pca_sglt$x[ , 1:3], heatmap_metadata_sglt_raw)

# Merge two data frames
heatmap_sglt_raw <- merge(heatmap_sglt_clr, heatmap_metadata_sglt)

# Subset the data frame
sglt_corr_data <- heatmap_sglt_raw %>% 
  relocate(PatientID, .before = "Prevotella 9") %>% 
  relocate(Timepoint, .before = "Prevotella 9") %>% 
  select(c(1:153, 161, 165, 179, 203:206, 217:219)) %>% 
  relocate(c("shannon", "pielou", "observed"), .before = PC1) %>% 
  mutate(across(.cols = 4:163, .fns=as.numeric)) %>% 
  pivot_longer(cols = 4:163, names_to = "Parameter", values_to = "Value")

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

#### Calculate the correlation matrix
sglt_corr <- round(cor(sglt_heatmap, use = "complete.obs"), 3)
#sglt_corr <- round(cor(sglt_heatmap), 1)
#sglt_corr[is.na(sglt_corr)] <- 0

#### Compute a matrix of correlation p-values
p.mat_sglt <- cor_pmat(sglt_heatmap)
#p.mat_sglt[is.na(p.mat_sglt)] <- 0

#### Visualize the inital correlation matrix
pheatmap(sglt_corr)

# Select correct columns and rows
sglt_corr_plot <- sglt_corr %>% 
  as.data.frame() %>% 
  select(c(171, 207, 211, 229, 236, 238, 285)) %>% 
  slice(c(1:10, 12:46, 48:50, 52:68, 70:75, 77, 79:124, 126:160)) %>% 
  rownames_to_column(var = "BL_Parameter") %>% 
  pivot_longer(cols = 2:8, names_to = "chg_Parameter", values_to = "Value")

sglt_corr_plot_p <- p.mat_sglt %>% 
  as.data.frame() %>% 
  column_to_rownames(var = "rowname") %>% 
  select(c(171, 207, 211, 229, 236, 238, 285)) %>% 
  slice(c(1:10, 12:46, 48:50, 52:68, 70:75, 77, 79:124, 126:160))%>% 
  rownames_to_column(var = "BL_Parameter") %>% 
  pivot_longer(cols = 2:8, names_to = "chg_Parameter", values_to = "p_value")

sglt_corr_plot <- sglt_corr_plot %>% 
  mutate(p_value = sglt_corr_plot_p$p_value) %>% 
  mutate(p_value_BH = p.adjust(p_value, method = "BH")) %>% 
  mutate(sig_p = ifelse(p_value_BH < 0.05, T, F)) %>% 
  mutate(p_if_sig_BH = ifelse(sig_p, p_value_BH, NA))

ggplot(sglt_corr_plot, aes(x = BL_Parameter, y = chg_Parameter, fill = Value)) +
  geom_tile(colour = "black") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
  coord_fixed() +
  geom_text(aes(label = round(p_value,3)), colour = "black", size = 1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))

# ___________________________________________________________________________ #

# Correlation plots for selected associations

glp_heatmap %>% 
  ggplot(aes(x = Value_BL_Alistipes, y = chg_Neutr)) +
  geom_point() +
  geom_smooth(method = "lm") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf, fill = "red", alpha = 0.2) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "green", alpha = 0.2) +
  geom_abline(aes(slope = 0, intercept = 0)) +
  ggpubr::stat_regline_equation(label.y = -0.5, aes(label = after_stat(eq.label))) +
  ggpubr::stat_regline_equation(label.y = -0.75, aes(label = after_stat(rr.label)))


sglt_heatmap %>% 
  ggplot(aes(x = Value_BL_Bilophila, y = chg_GFR)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggpubr::stat_regline_equation(label.y = 0, aes(label = after_stat(eq.label))) +
  ggpubr::stat_regline_equation(label.y = -0.5, aes(label = after_stat(rr.label))) +
  labs(x = "", y = "")

# ___________________________________________________________________________ #

# Top 3 for HbA1c

# GLP-1-RA

# Make data longer and remove "Value_BL_" string from genera names

glp_plot_data <- glp_heatmap %>% 
  pivot_longer(cols = 1:140, names_to = "BL_Parameter", values_to = "BL_Value") %>% 
  pivot_longer(cols = 1:140, names_to = "chg_Parameter", values_to = "chg_Value") %>% 
  mutate(Parameter_BL = gsub("Value_BL_", "", BL_Parameter)) %>% 
  mutate(Parameter_chg = gsub("chg_", "", chg_Parameter)) %>% 
  select(-c(BL_Parameter, chg_Parameter)) %>% 
  relocate(Parameter_BL, .before = BL_Value) %>% 
  relocate(Parameter_chg, .before = chg_Value)

e_s_hba1c_plot <- glp_plot_data %>% 
  filter(Parameter_BL == "Escherichia-Shigella" & Parameter_chg == "HbA1c") %>% 
  ggplot(aes(x = BL_Value, y = chg_Value)) +
  facet_grid(.~Parameter_BL, switch = "y", scales = "free") +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "",
       y = "Change in HbA1c (%, DCCT)",
       subtitle = "p = 0.00372 **, r = -0.819, p(adj) = 0.32 ns") +
  theme_classic() +
  theme(strip.text = element_text(size = 11, face = "italic"),
        axis.title.y = element_text(size = 12), 
        title = element_text(size = 14))

ls_nd3007_hba1c_plot <- glp_plot_data %>% 
  filter(Parameter_BL == "Lachnospiraceae ND3007 group" & Parameter_chg == "HbA1c") %>% 
  ggplot(aes(x = BL_Value, y = chg_Value)) +
  facet_grid(.~Parameter_BL, switch = "y", scales = "free") +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "",
       y = "", 
       subtitle = "p = 0.0106 *, r = 0.761, p(adj) = 0.62 ns") +
  theme_classic() +
  theme(strip.text = element_text(size = 12, face = "italic"),
        axis.title.y = element_text(size = 12), 
        title = element_text(size = 14))

catenibacterium_hba1c_plot <- glp_plot_data %>% 
  filter(Parameter_BL == "Catenibacterium" & Parameter_chg == "HbA1c") %>% 
  ggplot(aes(x = BL_Value, y = chg_Value)) +
  facet_grid(.~Parameter_BL, switch = "y", scales = "free") +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "", 
       y = "",
       subtitle = "p = 0.017 *, r = 0.728, p(adj) = 0.63 ns") +
  theme_classic() +
  theme(strip.text = element_text(size = 12, face = "italic"),
        axis.title.y = element_text(size = 12), 
        title = element_text(size = 14))

# SGLT-2i

# Make data longer and remove "Value_BL_" string from genera names

sglt_plot_data <- sglt_heatmap %>% 
  pivot_longer(cols = 1:160, names_to = "BL_Parameter", values_to = "BL_Value") %>% 
  pivot_longer(cols = 1:160, names_to = "chg_Parameter", values_to = "chg_Value") %>% 
  mutate(Parameter_BL = gsub("Value_BL_", "", BL_Parameter)) %>% 
  mutate(Parameter_chg = gsub("chg_", "", chg_Parameter)) %>% 
  select(-c(BL_Parameter, chg_Parameter)) %>% 
  relocate(Parameter_BL, .before = BL_Value) %>% 
  relocate(Parameter_chg, .before = chg_Value)

p_nk_hba1c_plot <- sglt_plot_data %>% 
  filter(Parameter_BL == "Prevotellaceae NK3B31 group" & Parameter_chg == "HbA1c") %>% 
  ggplot(aes(x = BL_Value, y = chg_Value)) +
  facet_grid(.~Parameter_BL, switch = "y", scales = "free") +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "",
       y = "",
       subtitle = "p = 0.005 **, r = 0.836, p(adj) = 0.76 ns") +
  theme_classic() +
  theme(strip.text = element_text(size = 12, face = "italic"),
        axis.title.y = element_text(size = 12), 
        title = element_text(size = 14))

cs_hba1c_plot <- sglt_plot_data %>% 
  filter(Parameter_BL == "Candidatus Soleaferrea" & Parameter_chg == "HbA1c") %>% 
  ggplot(aes(x = BL_Value, y = chg_Value)) +
  facet_grid(.~Parameter_BL, switch = "y", scales = "free") +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "",
       y = "",
       subtitle = "p = 0.006 **, r = -0.823, p(adj) = 0.76 ns") +
  theme_classic() +
  theme(strip.text = element_text(size = 12, face = "italic"),
        axis.title.y = element_text(size = 12), 
        title = element_text(size = 14))

uba_hba1c_plot <- sglt_plot_data %>% 
  filter(Parameter_BL == "UBA1819" & Parameter_chg == "HbA1c") %>% 
  ggplot(aes(x = BL_Value, y = chg_Value)) +
  facet_grid(.~Parameter_BL, switch = "y", scales = "free") +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "",
       y = "",
       subtitle = "p = 0.018 *, r = -0.758, p(adj) = 0.99 ns") +
  theme_classic() +
  theme(strip.text = element_text(size = 12, face = "italic"),
        axis.title.y = element_text(size = 12), 
        title = element_text(size = 14))

hba1c_comb_plot <- ggpubr::ggarrange(e_s_hba1c_plot, ls_nd3007_hba1c_plot, 
                                     catenibacterium_hba1c_plot, p_nk_hba1c_plot,
                                     cs_hba1c_plot, uba_hba1c_plot,
                                     nrow = 2, ncol = 3)

ggsave("hba1c_corr_plot.svg", dpi = 300, width = 13, height = 11)
  
# ___________________________________________________________________________ #

# Top 3 for BMI

# GLP-1-RA

tb_bmi_plot <- glp_plot_data %>% 
  filter(Parameter_BL == "Turicibacter" & Parameter_chg == "BMI") %>% 
  ggplot(aes(x = BL_Value, y = chg_Value)) +
  facet_grid(.~Parameter_BL, switch = "y", scales = "free") +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "",
       y = "Change in BMI",
       subtitle = "p = 0.00738 **, r = -0.783, p(adj) = 0.54 ns") +
  theme_classic() +
  theme(strip.text = element_text(size = 12, face = "italic"),
        axis.title.y = element_text(size = 12), 
        title = element_text(size = 14))

al_bmi_plot <- glp_plot_data %>% 
  filter(Parameter_BL == "Allisonella" & Parameter_chg == "BMI") %>% 
  ggplot(aes(x = BL_Value, y = chg_Value)) +
  facet_grid(.~Parameter_BL, switch = "y", scales = "free") +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "",
       y = "",
       subtitle = "p = 0.01670 *, r = 0.729, p(adj) = 0.62 ns") +
  theme_classic() +
  theme(strip.text = element_text(size = 12, face = "italic"),
        axis.title.y = element_text(size = 12), 
        title = element_text(size = 14))

de_bmi_plot <- glp_plot_data %>% 
  filter(Parameter_BL == "Defluviitaleaceae UCG-011" & Parameter_chg == "BMI") %>% 
  ggplot(aes(x = BL_Value, y = chg_Value)) +
  facet_grid(.~Parameter_BL, switch = "y", scales = "free") +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "",
       y = "",
       subtitle = "p = 0.02350 *, r = -0.703, p(adj) = 0.63 ns") +
  theme_classic() +
  theme(strip.text = element_text(size = 12, face = "italic"),
        axis.title.y = element_text(size = 12), 
        title = element_text(size = 14))

# SGLT-2i

eu_bmi_plot <- sglt_plot_data %>% 
  filter(Parameter_BL == "[Eubacterium] ventriosum group" & Parameter_chg == "BMI") %>% 
  ggplot(aes(x = BL_Value, y = chg_Value)) +
  facet_grid(.~Parameter_BL, switch = "y", scales = "free") +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "",
       y = "",
       subtitle = "p = 0.00184 **, r = 0.878, p(adj) = 0.57 ns") +
  theme_classic() +
  theme(strip.text = element_text(size = 12, face = "italic"),
        axis.title.y = element_text(size = 12), 
        title = element_text(size = 14))

ad_bmi_plot <- sglt_plot_data %>% 
  filter(Parameter_BL == "Adlercreutzia" & Parameter_chg == "BMI") %>% 
  ggplot(aes(x = BL_Value, y = chg_Value)) +
  facet_grid(.~Parameter_BL, switch = "y", scales = "free") +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "",
       y = "",
       subtitle = "p = 0.00368 **, r = -0.850, p(adj) = 0.76 ns") +
  theme_classic() +
  theme(strip.text = element_text(size = 12, face = "italic"),
        axis.title.y = element_text(size = 12), 
        title = element_text(size = 14))

pa_bmi_plot <- sglt_plot_data %>% 
  filter(Parameter_BL == "Parabacteroides" & Parameter_chg == "BMI") %>% 
  ggplot(aes(x = BL_Value, y = chg_Value)) +
  facet_grid(.~Parameter_BL, switch = "y", scales = "free") +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "",
       y = "",
       subtitle = "p = 0.01190 *, r = 0.787, p(adj) = 0.96 ns") +
  theme_classic() +
  theme(strip.text = element_text(size = 12, face = "italic"),
        axis.title.y = element_text(size = 12), 
        title = element_text(size = 14))

bmi_comb_plot <- ggpubr::ggarrange(tb_bmi_plot, al_bmi_plot, de_bmi_plot, 
                                   eu_bmi_plot, ad_bmi_plot, pa_bmi_plot,
                                   nrow = 2, ncol = 3)

ggsave("bmi_corr_plot.svg", dpi = 300, width = 13, height = 11)

# ___________________________________________________________________________ #

# Additional visualizations

# Pull out BL-data
glp_corr_BL <- glp_corr_data %>% 
  filter(Timepoint == "I") %>% 
  dplyr::rename(Parameter_BL = Parameter, 
                Value_BL = Value)

# Pull out other timepoints' data
glp_corr_TP3 <- glp_corr_data %>% 
  filter(Timepoint %in% c("III")) %>% 
  dplyr::rename(Parameter_TP3 = Parameter, 
                Value_TP3 = Value)

# Merge two data frames together
glp_plot_comb <- merge(glp_corr_BL, glp_corr_TP3, by = "PatientID", all = TRUE) %>% 
  select(-c(SampleID.x, SampleID.y)) %>% 
  pivot_wider(names_from = Parameter_BL, 
              values_from = Value_BL, names_prefix = "BL_") %>% 
  pivot_wider(names_from = Parameter_TP3, values_from = Value_TP3, 
              names_prefix = "TP3_")
  

glp_plot_comb %>% 
  ggplot(aes(x = BL_Neutr, y = TP3_Neutr)) +
  geom_point(aes(colour = BL_Alistipes, size = 5)) +
  scale_colour_gradient2(low="#22c1c3", mid = "white", midpoint = 5, high="#fdbb2d") +
#  geom_smooth(method = "lm") +
  geom_abline(aes(slope = , intercept = 0)) + #slope: y = "a"x+b, intercept: y = ax + "b"
  ggpubr::stat_regline_equation(label.y = -0.25, aes(label = after_stat(eq.label))) +
  ggpubr::stat_regline_equation(label.y = -0.5, aes(label = after_stat(rr.label))) +
  guides(size = "none")
