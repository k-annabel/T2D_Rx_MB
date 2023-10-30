## STANDARD SCRIPT FOR MICROBIOME STATISTICAL ANALYSIS

### Last modified: Sys.Date()

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

# Make a TreeSummarizedExperiment
tse_taxa <- TreeSummarizedExperiment(assays =  SimpleList(counts = counts),
                                     colData = DataFrame(samples),
                                     rowData = DataFrame(tax))

# Remove features containing human DNA
tse_prelim <- tse_taxa[!rowData(tse_taxa)$Phylum == "" & 
                  !rowData(tse_taxa)$Class == "" & 
                  !rowData(tse_taxa)$Kingdom == "Unassigned", ]

# Exclude samples
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

# Filter genera by prevalence
tse_genus <- tse_genus[rowSums(assay(tse_genus, "clr")) > 0, ]
  
# Separate by medication
tse_glp <- tse_genus[ , colData(tse_genus)$Medication == "GLP-1-RA"]
tse_sglt <- tse_genus[ , colData(tse_genus)$Medication == "SGLT-2"]

# ___________________________________________________________________________ #

# Make taxaplots

## GLP-1-RA

### Get top 20 genera
top_taxa_list <- getTopFeatures(tse_genus,
                                top = 19,
                                assay.type = "relabundance")

### Remove "Genus:" string
top_taxa_list <- stringr::str_remove_all(top_taxa_list, "Genus:")

### Rename less abundant genera as "Other"
genus_renamed <- lapply(rowData(tse_glp)$Genus,
                        function(x){if (x %in% top_taxa_list) {x} else {"Other"}})

### Make changes in rowData
rowData(tse_glp)$Genus <- as.character(genus_renamed)

# __________________________________ V1 ______________________________________ #

### Melt data
glp_melt_raw <- meltAssay(tse_glp, 
                      assay.type = "relabundance", 
                      add_row_data = TRUE, 
                      add_col_data = TRUE) 

### Tidy and transform
glp_melt <- glp_melt_raw %>% 
  # filter(!str_detect(FeatureID, "Family:")) %>% 
  select(!c(4:8, 10, 14:75)) %>% 
  mutate(FeatureID = case_when(FeatureID %in% top_taxa_list ~ FeatureID,
                           TRUE ~ "Other")) %>%
  mutate(Timepoint2 = case_match(Timepoint,
                   "I" ~ "Baseline",
                   "II" ~ "1st month",
                   "III" ~ "3rd month",
                   "IV" ~ "12th month"))

### Pull all genera
glp_taxaplot_genera <- glp_melt %>% 
  arrange(desc(relabundance)) %>% 
  pull(FeatureID)

### Leave in only unique
glp_taxaplot_genera <- unique(glp_taxaplot_genera)

### Move "Other" to 1st position
glp_taxaplot_genera <- c("Other", glp_taxaplot_genera)

### Remove duplicate "Other"
glp_taxaplot_genera <- glp_taxaplot_genera[-10]

### Reorder genera as "Other" as first factor level
glp_melt$FeatureID <- fct_relevel(glp_melt$FeatureID, glp_taxaplot_genera)

### Rename PatientIDs for the plot
glp_taxaplot_patientids <- c("Patient 1", "Patient 2", "Patient 3", "Patient 4",
                             "Patient 5", "Patient 6", "Patient 7", "Patient 8",
                             "Patient 9", "Patient 10")

### Make a plot
glp_taxaplot <- glp_melt %>% 
  ggplot(aes(x = PatientID, y = relabundance, fill = FeatureID)) +
  geom_col(position = "fill") +
  # geom_bar(aes(fill = FeatureID), stat = "identity", position = "fill", colour = NA) +
  facet_grid(~factor(Timepoint2, levels = c("Baseline", "1st month", "3rd month", "12th month")), scales = "free_x") +
  scale_fill_paletteer_d("ggthemes::Tableau_20") +
  scale_y_continuous(expand = c(0,0), labels = scales::percent) +
  theme(legend.text = element_text(face = "italic")) +
  theme(axis.text.x = element_text(angle = 33, vjust = 0.75, hjust = 0.75)) +
  coord_cartesian(ylim=c(0,1)) +
  scale_x_discrete(labels = glp_taxaplot_patientids) + 
  labs(x = "",
       y = "Relative abundance (%)", 
       fill = "Genus") +
  guides(fill=guide_legend(ncol=1)) +
  ggtitle("Semaglutide (GLP1-RA)")

# __________________________________ V2 ______________________________________ #

### Get initial data for the plot
top_genera_plot_glp <- miaViz::plotAbundance(tse_glp, 
                                             rank = "Genus",
                                             assay_name = "relabundance", 
                                             order_rank_by = "abund", 
                                             one_facet = FALSE)

### Convert into data frame
top_genera_glp_data <- top_genera_plot_glp$data %>% 
  as.data.frame() %>% 
  mutate(Timepoint = case_when(
    grepl("-1$", X) ~ "I",
    grepl("-2$", X) ~ "II",
    grepl("-3$", X) ~ "III",
    grepl("-4$", X) ~ "IV",
    TRUE ~ NA
  )) %>% 
  dplyr::rename(Genus = colour_by,
                SampleID = X, 
                rel_abund = Y) %>% 
  mutate(Timepoint2 = case_match(Timepoint,
                                 "I" ~ "Baseline",
                                 "II" ~ "1st month",
                                 "III" ~ "3rd month",
                                 "IV" ~ "12th month"))

### Make a plot
glp_genera_plot <- top_genera_glp_data %>% 
  ggplot(aes(x = SampleID, y = rel_abund, fill = Genus)) +
  geom_col() +
  facet_grid(~factor(Timepoint2, levels = c("Baseline", "1st month", "3rd month", "12th month")), scales = "free_x") +
  scale_fill_paletteer_d("ggthemes::Tableau_20") +
  scale_y_continuous(expand = c(0,0), labels = scales::percent) +
  theme(legend.text = element_text(face = "italic")) +
  theme(axis.text.x = element_text(angle = 33, vjust = 0.75, hjust = 0.75)) +
  coord_cartesian(ylim=c(0,1)) +
  scale_x_discrete(labels = glp_taxaplot_patientids) + 
  labs(x = "",
       y = "Relative abundance (%)", 
       fill = "Genus") +
  guides(fill=guide_legend(ncol=1)) +
  ggtitle("Semaglutide (GLP1-RA)")

# ___________________________________________________________________________ #

## SGLT-2

### Get top 20 genera
top_taxa_list <- getTopFeatures(tse_genus,
                                top = 19,
                                assay.type = "relabundance")

### Remove "Genus:" string
top_taxa_list <- stringr::str_remove_all(top_taxa_list, "Genus:")

### Leave in names for top 20 genera
genus_renamed_sglt <- lapply(rowData(tse_sglt)$Genus,
                        function(x){if (x %in% top_taxa_list) {x} else {"Other"}})

### Make changes in rowData
rowData(tse_sglt)$Genus <- as.character(genus_renamed_sglt)


# __________________________________ V1 ______________________________________ #

### Melt data
sglt_melt_raw <- meltAssay(tse_sglt, 
                          assay.type = "relabundance", 
                          add_row_data = TRUE, 
                          add_col_data = TRUE) 

### Tidy and transform
sglt_melt <- sglt_melt_raw %>% 
  # filter(!str_detect(FeatureID, "Family:")) %>% 
  select(!c(4:8, 10, 14:75)) %>% 
  mutate(FeatureID = case_when(FeatureID %in% top_taxa_list ~ FeatureID, 
                           TRUE ~ "Other")) %>% 
  mutate(Timepoint2 = case_match(Timepoint,
                                 "I" ~ "Baseline",
                                 "II" ~ "1st month",
                                 "III" ~ "3rd month",
                                 "IV" ~ "12th month"))

### Pull all genera
sglt_taxaplot_genera <- sglt_melt %>% 
  arrange(desc(relabundance)) %>% 
  pull(FeatureID)

### Leave in only unique
sglt_taxaplot_genera <- unique(sglt_taxaplot_genera)

### Move "Other" to 1st position
sglt_taxaplot_genera <- c("Other", sglt_taxaplot_genera)

### Remove duplicate "Other"
sglt_taxaplot_genera <- sglt_taxaplot_genera[-13]

### Reorder genera as "Other" as first factor level
sglt_melt$FeatureID <- fct_relevel(sglt_melt$FeatureID, sglt_taxaplot_genera)

### Rename PatientIDs for the plot
sglt_taxaplot_patientids <- c("Patient 11", "Patient 12", "Patient 13", "Patient 14",
                             "Patient 15", "Patient 16", "Patient 17", "Patient 18",
                             "Patient 19", "Patient 20")

### Make a plot
sglt_taxaplot <- sglt_melt %>% 
  ggplot(aes(x = PatientID, y = relabundance, fill = FeatureID)) +
  geom_col(position = "fill") +
#  geom_bar(aes(fill = FeatureID), stat = "identity", position = "fill", colour = NA) +
  facet_grid(~factor(Timepoint2, levels = c("Baseline", "1st month", "3rd month", "12th month")), scales = "free_x") +
  scale_fill_manual(values = tableau_20_sglt) +
  scale_y_continuous(expand = c(0,0), labels = scales::percent) +
  theme(legend.text = element_text(face = "italic")) +
  theme(axis.text.x = element_text(angle = 33, vjust = 0.75, hjust = 0.75)) +
  coord_cartesian(ylim=c(0,1)) +
  scale_x_discrete(labels = sglt_taxaplot_patientids) + 
  labs(x = "",
       y = "Relative abundance (%)", 
       fill = "Genus") +
  guides(fill=guide_legend(ncol=1)) +
  ggtitle("Empagliflozin (SGLT-2i)")

tableau_20_sglt <- c("#4E79A7", "#A0CBE8", "#F28E2B", "#8CD17D", "#FFBE7D", 
                     "#B6992D", "#F1CE63", "#D37295", "#86BCB6", "#499894", 
                     "#FF9D9A", "#79706E", "#59A14F", "#FABFD2", "#D4A6C8", 
                     "#BAB0AC", "#9D7660", "#B07AA1", "#D37295", "#D7B5A6")

# __________________________________ V2 ______________________________________ #

### Get initial data for the plot
top_genera_plot_sglt <- miaViz::plotAbundance(tse_sglt, 
                                             rank = "Genus",
                                             assay_name = "relabundance", 
                                             order_rank_by = "abund", 
                                             one_facet = FALSE)

### Convert into data frame
top_genera_sglt_data <- top_genera_plot_sglt$data %>% 
  as.data.frame() %>% 
  mutate(Timepoint = case_when(
    grepl("-1$", X) ~ "I",
    grepl("-2$", X) ~ "II",
    grepl("-3$", X) ~ "III",
    grepl("-4$", X) ~ "IV",
    TRUE ~ NA
  )) %>% 
  dplyr::rename(Genus = colour_by,
                SampleID = X, 
                rel_abund = Y) %>% 
  mutate(Timepoint2 = case_match(Timepoint,
                                 "I" ~ "Baseline",
                                 "II" ~ "1st month",
                                 "III" ~ "3rd month",
                                 "IV" ~ "12th month"))

### Make a plot
sglt_genera_plot <- top_genera_sglt_data %>% 
  ggplot(aes(x = SampleID, y = rel_abund, fill = Genus)) +
  geom_col() +
  facet_grid(~factor(Timepoint2, levels = c("Baseline", "1st month", "3rd month", "12th month")), scales = "free_x") +
  scale_fill_manual(values = tableau_20_sglt) +
  scale_y_continuous(expand = c(0,0), labels = scales::percent) +
  theme(legend.text = element_text(face = "italic")) +
  theme(axis.text.x = element_text(angle = 33, vjust = 0.75, hjust = 0.75)) +
  coord_cartesian(ylim=c(0,1)) +
  scale_x_discrete(labels = sglt_taxaplot_patientids) + 
  labs(x = "",
       y = "Relative abundance (%)", 
       fill = "Genus") +
  guides(fill=guide_legend(ncol=1)) +
  ggtitle("Empagliflozin (SGLT-2i)")

# ___________________________________________________________________________ #

### Combine two plots
glp_sglt_taxaplot <- ggpubr::ggarrange(glp_genera_plot, sglt_genera_plot, 
                                       nrow = 2, common.legend = TRUE, 
                                       legend = "right")

### Export
ggplot2::ggsave("glp_sglt_taxaplot.svg", width = 12, height = 7, dpi = 600)

# ___________________________________________________________________________ #

# I. Effect of T2D medications towards MB

# Repeated measures ANOVA for a general comparison between timepoints

# https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r/ (rstatix)
# https://www.statology.org/repeated-measures-anova-in-r/

## GLP-1-RA

### Calculate alpha diversity
tse_glp <- estimateDiversity(tse_glp, assay.type = "relabundance", index = "shannon")
tse_glp <- estimateEvenness(tse_glp, assay.type = "relabundance", index = "pielou")
tse_glp <- estimateRichness(tse_glp, assay.type = "relabundance", index = "observed")

### Prepare data

#### Extract CLR-transformed abundance values for top genera per sample
tse_glp_clr <- as.data.frame(assay(tse_glp, "clr"))

#### Extract metadata with alpha diversity indices
tse_glp_samples <- as.data.frame(colData(tse_glp))

#### Perform PCA 
pca_glp <- prcomp(t(tse_glp_clr))

#### Combine PCs with metadata and alpha diversity
glp_alpha_beta_raw <- data.frame(pca_glp$x[ , 1:5],
                                 tse_glp_samples)


#### Perform ANOVA
model_glp <- aov(observed~factor(Timepoint), data = glp_alpha_beta_raw)
summary(model_glp)

# ___________________________________________________________________________ #

## SGLT-2

### Calculate alpha diversity
tse_sglt <- estimateDiversity(tse_sglt, assay.type = "relabundance", index = "shannon")
tse_sglt <- estimateEvenness(tse_sglt, assay.type = "relabundance", index = "pielou")
tse_sglt <- estimateRichness(tse_sglt, assay.type = "relabundance", index = "observed")

### Prepare data

#### Extract CLR-transformed abundance values for top genera per sample
tse_sglt_clr <- as.data.frame(assay(tse_sglt, "clr"))

#### Extract metadata with alpha diversity indices
tse_sglt_samples <- as.data.frame(colData(tse_sglt))

#### Perform PCA 
pca_sglt <- prcomp(t(tse_sglt_clr))

#### Combine PCs with metadata and alpha diversity
sglt_alpha_beta_raw <- data.frame(pca_sglt$x[ , 1:5],
                                 tse_sglt_samples)


#### Perform ANOVA
model_sglt <- aov(observed~factor(Timepoint), data = sglt_alpha_beta_raw)
summary(model_sglt)

# ___________________________________________________________________________ #

## Create a list for timepoints
timepoints <- c("II", "III", "IV")

## Create a list for alpha diversity indices
alpha_div <- c("Shannon", "Pielou", "Observed")

## Create an empty vector for results
alpha_diversity_results <- matrix(nrow = length(timepoints), 
                                  ncol = length(alpha_div), 
                                  dimnames = list(timepoints, alpha_div))

## Execute the for-loop

for (i in length(timepoints)){
    tse_glp_test <- tse_glp[ ,colData(tse_glp)$Timepoint == timepoints[i]]
  
}

# ___________________________________________________________________________ #

# Differential abundance analysis

## Gather top taxa
top_taxa <- names(rowSums(assay(tse_genus, "clr")) > 0)
top_taxa <- stringr::str_remove_all(top_taxa, "Genus:")
top_taxa <- stringr::str_remove_all(top_taxa, "Family:")

## Convert rowData into data frame
tse_glp_rowData <- as.data.frame(rowData(tse_glp))

## Combine tables 
diff_abund_data <- tse_glp_metadata %>% 
  rownames_to_column(var = "SampleID") %>% 
  select(1:4)

## Convert data
tse_glp_melt_prelim <- meltAssay(tse_glp, 
                          assay.type = "clr", 
                          add_row_data = TRUE, 
                          add_col_data = TRUE)

tse_glp_melt <- tse_glp_melt_prelim %>% 
  select(2:3, 9, 11:12) %>% 
  group_by(PatientID) %>% 
  filter(n() != 1) %>% 
  arrange(Timepoint, PatientID) %>% 
  ungroup()

tse_glp_melt <- tse_glp_melt_prelim %>% 
  select(1:3, 9, 11:12) %>% 
  group_by(PatientID) %>% 
  filter(all(c("I", "II", "III", "IV") %in% Timepoint)) %>% 
  ungroup()
  
## Create an empty matrix for results
diff_abund_results <- matrix(nrow = length(top_taxa), 
                             ncol = length(timepoints),
                             dimnames = list(top_taxa, timepoints))

for (i in length(top_taxa)){
  
  for(j in length(timepoints)){
  
  diff_abund_test_BL <- tse_glp_melt %>% 
    filter(Genus == top_taxa[i]) %>% 
    filter(Timepoint == "I")
    
  diff_abund_test <- tse_glp_melt %>% 
    filter(Genus == top_taxa[i]) %>% 
    filter(Timepoint == timepoints[j])
    
    t.test_taxa <- t.test(diff_abund_test_BL$clr, 
                          diff_abund_test$clr, 
                          paired = TRUE)$p.value
    
    diff_abund_results <- diff_abund_results[t.test_taxa, timepoints[j]]
  }
}

# ___________________________________________________________________________ #

# II. Effect of MB towards the efficacy of T2D medications

## GLP-1-RA

### Transform data for the correlation matrix

#### Pull out top genera based on CLR-values and prevalence (5 of 20 samples)
top_genera_glp <- mia::getTopTaxa(tse_glp,
                                  assay_name = "clr",
                                  detection = 0,
                                  prevalence = 5/20)
# There's only five genera now? - Bacteroides, Blautia, Alistipes, Roseburia, Dorea

#### Extract CLR-values for the heatmap
heatmap_glp_clr_raw <- as.data.frame(t(assay(tse_glp, "clr"))) 

#### Filter out less prevalent genera
heatmap_glp_clr <- heatmap_glp_clr_raw %>% 
  select(all_of(top_genera_glp)) %>% 
  rownames_to_column(var = "SampleID")

#### Extract metadata 
heatmap_metadata_glp_raw <- as.data.frame(colData(tse_glp)) %>% 
  rownames_to_column(var = "SampleID")

#### Add PCA data
heatmap_metadata_glp <- data.frame(pca_glp$x[ , 1:5], heatmap_metadata_glp_raw)

#### Merge two data frames
heatmap_glp_raw <- merge(heatmap_glp_clr, heatmap_metadata_glp)

#### Subset the data frame
glp_corr_data <- heatmap_glp_raw %>% 
  select(!c(14:18, 67:76)) %>% 
  relocate(PatientID, .before = Bacteroides) %>% 
  relocate(Timepoint, .before = Bacteroides) %>% 
  relocate(c("shannon", "pielou", "observed"), .before = Weight) %>% 
  mutate(across(.cols = 4:64, .fns=as.numeric)) %>% 
  pivot_longer(cols = 4:64, names_to = "Parameter", values_to = "Value")


#### Pull out BL-data
glp_corr_BL <- glp_corr_data %>% 
  filter(Timepoint == "I") %>% 
  select(-Timepoint) %>% 
  rename("Value" = "Value_BL")

#### Pull out other timepoints' data
glp_corr_TP <- glp_corr_data %>% 
  filter(!Timepoint %in% c("I"))

#### Merge two data frames together
glp_heatmap_comb <- merge(glp_corr_BL, glp_corr_TP, 
                          by = c("PatientID", "Parameter"),
                          all = TRUE)

#### Tidy the merged data frame
glp_heatmap <- glp_heatmap_comb %>% 
  arrange(PatientID, Parameter, SampleID.y) %>% 
  mutate(chg = Value - Value_BL) %>% 
  select(-Value) %>% 
  pivot_wider(names_from = Parameter, values_from = c(Value_BL, chg)) %>% 
  filter(!Timepoint %in% c("II", "IV")) %>% 
  select(-c(PatientID, SampleID.x, SampleID.y, Timepoint))
# dim(glp_heatmap) - 17 x 122

#### Calculate the correlation matrix
glp_corr <- round(cor(glp_heatmap, use = "complete.obs"), 3)
#glp_corr <- round(cor(glp_heatmap), 1)
glp_corr[is.na(glp_corr)] <- 0

#### Compute a matrix of correlation p-values
p.mat_glp <- cor_pmat(glp_heatmap)
p.mat_glp[is.na(p.mat_glp)] <- 0

#### Visualize the inital correlation matrix
pheatmap(glp_corr)

### Make a submatrix for analysis (all TPs included)
glp_corr_alam <- glp_corr[c(5, 7:8, 12, 41), c(67, 78, 80)]
glp_pmat_alam <- p.mat_glp[c(5, 7:8, 12, 41), c(67, 78, 80)]
glp_labels_row <- c("Alistipes", "Dorea", "Roseburia", "Bacteroides", "Blautia")
glp_labels_col <- c("Blood glucose", "BMI", "HbA1c (%)")
glp_alam_hm <- pheatmap(glp_corr_alam, 
                        labels_row = glp_labels_row, 
                        labels_col = glp_labels_col, 
                        fontsize = 18, angle_col = 45, filename = "glp_alam_hm.png")
                        
                        
                        
#### Tidy the merged data frame
glp_heatmap_HbA1c <- glp_heatmap_comb %>% 
  arrange(PatientID, Parameter, SampleID.y) %>% 
  mutate(chg = Value - Value_BL) %>% 
  select(-Value) %>% 
  filter(Parameter %in% c("Alistipes", "Bacteroides", "Blautia", "Dorea", 
                          "HbA1c", "PC1", "PC2", "PC3", "PC4", "PC5", 
                          "Roseburia", "observed", "pielou", "shannon")) %>% 
  pivot_wider(names_from = Parameter, values_from = c(Value_BL, chg)) %>% 
  pivot_wider(names_from = Timepoint, values_from = chg_HbA1c) %>% 
  select(-c(PatientID, SampleID.x, SampleID.y))

#### Calculate the correlation matrix
glp_corr_HbA1c <- round(cor(glp_heatmap_HbA1c, use = "complete.obs"), 3)
glp_corr_HbA1c <- round(cor(glp_heatmap_HbA1c), 3)

alam_glp_corr_HbA1c <- glp_corr_HbA1c[c(28:30), c(1:14)]

heatmap(alam_glp_corr_HbA1c)


glp_corr_HbA1c <- glp_corr[c(5, 7:8, 12, 33:37, 41, 59:61), c(80)]
corrplot::corrplot(glp_corr_HbA1c)


## SGLT-2

### Transform data for the correlation matrix

#### Convert colData into a data frame
heatmap_metadata_sglt_raw <- as.data.frame(colData(tse_sglt)) %>% 
  rownames_to_column(var = "SampleID")

#### Subset the data frame
sglt_corr_data <- sglt_coldata_df %>% 
  rownames_to_column(var = "SampleID") %>% 
  select(!c(4:8, 56:66))

#### Pull out BL-metadata
sglt_corr_BL <- sglt_corr_data %>% 
  filter(Timepoint == "I") %>% 
  select(-Timepoint) %>% 
  pivot_longer(cols = 3:52, names_to = "Parameter", values_to = "Value", 
               values_transform = as.numeric) %>% 
  dplyr::rename(BL_Value = Value)

#### Pull out other timepoints' data
sglt_corr_TP <- sglt_corr_data %>% 
  filter(Timepoint != "I") %>% 
  pivot_longer(cols = 4:53, names_to = "Parameter", values_to = "Value", 
               values_transform = as.numeric)

#### Merge two data frames together
sglt_corr_comb <- merge(sglt_corr_BL, sglt_corr_TP, by = c("PatientID", "Parameter")) %>% 
  mutate(Timepoint = factor(Timepoint)) %>% 
  mutate(Timepoint = fct_relevel(Timepoint, c("II", "III", "IV"))) %>% 
  arrange(Timepoint) %>% 
  mutate(chg_Value = Value - BL_Value) %>% 
  select(!c(SampleID.x, SampleID.y, Value)) %>% 
  pivot_wider(names_from = "Parameter", values_from = c(BL_Value, chg_Value)) %>% 
  select(!c(PatientID, Timepoint))
# dim(sglt_corr_comb) - 28 x 100

#### Calculate the correlation matrix
sglt_corr <- round(cor(sglt_corr_comb, use = "complete.obs"), 3)
#glp_corr <- round(cor(glp_heatmap), 1)
sglt_corr[is.na(sglt_corr)] <- 0

#### Compute a matrix of correlation p-values
p.mat_sglt <- cor_pmat(sglt_corr_comb)
p.mat_sglt[is.na(p.mat_sglt)] <- 0

#### Visualize the inital correlation matrix
pheatmap(sglt_corr)

# ___________________________________________________________________________ #

## HbA1c and BMI plots

### GLP-1-RA

comp <- list(c("I", "III"), c("I", "IV"))

hba1c_plot_glp <- heatmap_metadata_glp_raw %>% 
  filter(Timepoint != "II") %>% 
  ggplot(aes(x = Timepoint, y = HbA1c)) + 
  facet_wrap(. ~ Medication) +
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.4) +
  ggpubr::stat_compare_means(comparisons = comp) +
  labs(x = "", y = "") +
  scale_x_discrete(labels = c("Enne ravi", "3. kuu", "12. kuu")) +
  theme(
    axis.text.x = element_text(
      angle = 33,
      hjust = 1,
      vjust = 0.85, 
      size = rel(1.2)
    ))

ggsave("glp_hba1c_plot.png", width = 5)

bmi_plot_glp <- heatmap_metadata_glp_raw %>% 
  filter(Timepoint != "II") %>% 
  ggplot(aes(x = Timepoint, y = BMI)) + 
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.4) +
  ggpubr::stat_compare_means(comparisons = comp) +
  labs(x = "", y = "") +
  scale_x_discrete(labels = c("Enne ravi", "3. kuu", "12. kuu")) +
  theme(
    axis.text.x = element_text(
      angle = 33,
      hjust = 1,
      vjust = 0.85
    ))

### SGLT-2

hba1c_plot_sglt <- heatmap_metadata_sglt_raw %>% 
  filter(Timepoint != "II") %>% 
  ggplot(aes(x = Timepoint, y = HbA1c)) + 
  facet_wrap(. ~ Medication) +
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.4) +
  ggpubr::stat_compare_means(comparisons = comp) +
  labs(x = "", y = "") +
  scale_x_discrete(labels = c("Enne ravi", "1. kuu", "3. kuu")) +
  theme(
    axis.text.x = element_text(
      angle = 33,
      hjust = 1,
      vjust = 0.85,
      size = rel(1.2)
    ))

ggsave("sglt_hba1c_plot.png", width = 5)

bmi_plot_sglt <- heatmap_metadata_sglt_raw %>% 
  filter(Timepoint != "II") %>% 
  ggplot(aes(x = Timepoint, y = BMI)) + 
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.4) +
  ggpubr::stat_compare_means(comparisons = comp) +
  labs(x = "", y = "") +
  scale_x_discrete(labels = c("Enne ravi", "1. kuu", "3. kuu")) +
  theme(
    axis.text.x = element_text(
      angle = 33,
      hjust = 1,
      vjust = 0.85
    ))

### Both medications together as facet plots

metadata_tse <- as.data.frame(colData(tse))

hba1c_plot <- metadata_tse %>% 
  filter(Timepoint != "II") %>% 
  ggplot(aes(x = Timepoint, y = HbA1c)) + 
  facet_wrap(. ~ Medication) +
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.4) +
  ggpubr::stat_compare_means(comparisons = comp) +
  labs(x = "", y = "HbA1c (%)") +
  scale_x_discrete(labels = c("Baseline", "3rd month", "12th month")) +
  theme(
    axis.text.x = element_text(
      angle = 33,
      hjust = 1,
      vjust = 0.85, 
      size = rel(1.2)
    ))
  
ggsave("hba1c_plot.png", height = 5, width = 6.5)

bmi_plot <- metadata_tse %>% 
  filter(Timepoint != "II") %>% 
  ggplot(aes(x = Timepoint, y = BMI)) + 
  facet_wrap(. ~ Medication) +
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group = PatientID), color = "grey", linewidth = 0.4) +
  ggpubr::stat_compare_means(comparisons = comp) +
  labs(x = "", y = "BMI") +
  scale_x_discrete(labels = c("Baseline", "3rd month", "12th month")) +
  theme(
    axis.text.x = element_text(
      angle = 33,
      hjust = 1,
      vjust = 0.85, 
      size = rel(1.2)
    ))

ggsave("bmi_plot.png", height = 5, width = 6.5)
      

# ___________________________________________________________________________ #

# Scatterplots for presentation

## GLP-1-RA

scatter_blautia_bmi <- glp_heatmap %>% 
  ggplot(aes(x = Value_BL_Blautia, y = chg_BMI)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggpubr::stat_regline_equation(label.y = 0, aes(label = after_stat(eq.label))) +
  ggpubr::stat_regline_equation(label.y = -0.5, aes(label = after_stat(rr.label))) +
  labs(x = "", y = "")

ggsave("scatter_blautia_bmi.png")
