# =====================================================================
# PROJECT: CYP450-AF-PCA
# PURPOSE: PCA-Driven Diagnostic Pipeline for AlphaFold Structural Auditing
# AUTHOR: [Your Name] | Ohio Supercomputer Center (OSC)
# VERSION: 1.0.1 (Production/GitHub Release)
# =====================================================================

# ---------------------------------------------------------------------
# Section 1: DATA PREPARATION & CLEANING
# ---------------------------------------------------------------------
# Load dataset
df <- read.csv("alphafold_data_mining_set_1000.csv")

# Clean pocket volume data (remove '%' and cast to numeric)
df$Pocket_Volume_Change_vs_PDB_Percent <- as.numeric(gsub("%", "", df$Pocket_Volume_Change_vs_PDB_Percent))

# Filter for valid structural benchmarks (exclude rows without crystal structure comparisons)
df_clean <- df[!is.na(df$Pocket_Distortion_Status) & df$Pocket_Distortion_Status != "", ]
df_clean$Pocket_Distortion_Status <- as.factor(df_clean$Pocket_Distortion_Status)

# Define Core Features for PCA/Clustering
# Using Local vs. Global pLDDT as the primary variance drivers
feature_cols <- c("Pocket_Mean_pLDDT", "Pocket_Min_pLDDT", "Global_Mean_pLDDT")
features <- df_clean[, feature_cols]

# Normalize features (Z-score scaling) - Critical for PCA and K-Means
scaled_features <- scale(features)

# ---------------------------------------------------------------------
# Section 2: DIAGNOSTIC METRIC: SCAFFOLD-POCKET DIVERGENCE
# ---------------------------------------------------------------------
# Quantify the "Confidence Gap" between the scaffold and the active site
df_clean$Scaffold_Pocket_Divergence <- df_clean$Global_Mean_pLDDT - df_clean$Pocket_Mean_pLDDT

hist(df_clean$Scaffold_Pocket_Divergence, 
     main="Metric Distribution: Structural Dissonance", 
     xlab="Divergence (Global Mean - Pocket Mean pLDDT)", 
     col="skyblue", border="white", breaks=20)

# Save it to a FILE
dev.copy(pdf, "Scaffold_Pocket_Divergence_Distribution.pdf", width=7, height=5)
dev.off() # This will show "null device 1" again, but your file is now saved!

# ---------------------------------------------------------------------
# Section 3: HIERARCHICAL CLUSTERING (Dendrogram Analysis)
# ---------------------------------------------------------------------
# 1. Load Data and Clean
df <- read.csv("alphafold_data_mining_set_1000.csv")

# FIX 1: Strip the '%' symbol and convert text to actual numbers
df$Pocket_Volume_Change_vs_PDB_Percent <- as.numeric(gsub("%", "", df$Pocket_Volume_Change_vs_PDB_Percent))

# FIX 2: Filter out rows that do not have real crystal structure calculations
# This prevents K-Means and Hierarchical distance math from crashing!
df_clean <- df[!is.na(df$Pocket_Distortion_Status) & df$Pocket_Distortion_Status != "", ]

# Convert target label to a clean factor category for analysis later
df_clean$Pocket_Distortion_Status <- as.factor(df_clean$Pocket_Distortion_Status)

# Update feature selection to use the clean dataset rows
features <- df_clean[, c("Pocket_Mean_pLDDT", "Pocket_Min_pLDDT", "Global_Mean_pLDDT")]

# Normalize/Scale the data (Critical step taught in Week 2 for K-Means)
scaled_features <- scale(features)
```

```{r}
# =====================================================================
# Section 2: EXPLORATORY DANA ANALYSIS (EDA)
# =====================================================================
# Exploratory Data Analysis (EDA)
# Create a new column measuring the "Pocket Deficit" 
df_clean$Scaffold_Pocket_Divergence <- df_clean$Global_Mean_pLDDT - df_clean$Pocket_Mean_pLDDT

# Check distribution of divergence (Dissonance check)
# Plot to the SCREEN (RStudio window)
hist(df_clean$Scaffold_Pocket_Divergence, 
     main="Metric Distribution: Structural Dissonance", 
     xlab="Divergence (Global Mean - Pocket Mean pLDDT)", 
     col="skyblue", border="white", breaks=20)

# save it to a FILE
dev.copy(pdf, "Scaffold_Pocket_Divergence_Distribution.pdf", width=7, height=5)
dev.off() # This will show "null device 1" again, but your file is now saved!
```

```{r}
# =====================================================================
# Section 3: HIERARCHICAL CLUSTERING (Tuned for 14x8 High-Res)
# =====================================================================
# TIGHTENED MARGINS: Bottom reduced from 8 to 5.5
par(mar = c(5.5, 4, 4, 2)) 

# Get tree metrics
max_h <- max(hc_result$height)
tot_l <- length(hc_result$order)

# Plot the tree
plot(hc_result, 
     labels = df_clean$UniProt_ID, 
     main = "Hierarchical Clustering of Structural Conformations",
     sub = "", xlab = "", ylab = "Clustering Height", 
     cex = 0.5,    
     hang = -1)    

# Draw red boxes
rect.hclust(hc_result, k=3, border="red")

# Place labels (Keep original logic - it's looking great)
y_pos <- max_h * 0.93 

text(x = tot_l * 0.08, y = y_pos, labels = "STRUCTURAL\nDISSONANCE", 
     col = "red", font = 2, cex = 0.8, adj = 0.5)

text(x = tot_l * 0.35, y = y_pos, labels = "TRANSITIONAL\nSTABILITY", 
     col = "red", font = 2, cex = 0.8, adj = 0.5)

text(x = tot_l * 0.72, y = y_pos, labels = "STRUCTURAL\nCONSONANCE", 
     col = "red", font = 2, cex = 0.8, adj = 0.5)

# PULL X-AXIS UP:
mtext("Protein Entry ID (UniProt)", side = 1, line = 1.0, cex = 1)

# Save to a png file
dev.copy(png, "Figure1_Dendrogram_CYP450.png", width=2000, height=1400, res=300)
dev.off()

# Save to a pdf file
dev.copy(pdf, "Figure1_Dendrogram_CYP450.pdf", width=10, height=7)
dev.off()

# =====================================================================
# Section 4: k-MEANS CLUSTERING & VALIDATION
# =====================================================================
# 1. The Elbow Method: Finding Optimal K
wcss <- vector()
for (i in 1:10) {
  wcss[i] <- sum(kmeans(scaled_features, centers=i)$withinss)
}
plot(1:10, wcss, type="b", main="The Elbow Method for Optimal K",
     xlab="Number of Clusters (K)", ylab="Within-Cluster Sum of Squares")

# 2. Run K-means with K = 3 (Pathway A, Pathway B, and Consonance)
set.seed(42)
km_result <- kmeans(scaled_features, centers=3, nstart=25)
df_clean$Cluster_Label <- as.factor(km_result$cluster)

# 3. Cluster Summary Statistics
cluster_summary <- aggregate(
  df_clean[, c("Global_Mean_pLDDT", "Pocket_Mean_pLDDT", "Pocket_Min_pLDDT", "Scaffold_Pocket_Divergence")],
  by = list(Cluster = df_clean$Cluster_Label),
  FUN = mean
)
print("=== CLUSTER PHENOTYPE SUMMARY ===")
print(cluster_summary)

# 4. Cross-tabulation (Comparison of HC and K-Means)
print("=== CLUSTERING AGREEMENT (HC vs K-Means) ===")
print(table(hc_cluster, km_result$cluster))

# ---------------------------------------------------------------------
# BIOLOGICAL VALIDATION (Data Mining Insights)
# ---------------------------------------------------------------------
# Check if AlphaFold confidence correlates with Taxonomic Domain
print("=== TAXONOMIC DISTRIBUTION BY CLUSTER ===")
print(table(df_clean$Taxonomic_Domain, df_clean$Cluster_Label))

# Verify if statistical clusters match physical structural distortion
print("=== PHYSICAL DISTORTION VALIDATION ===")
print(table(df_clean$Pocket_Distortion_Status, df_clean$Cluster_Label))

# ---------------------------------------------------------------------
# Section 5: PRINCIPAL COMPONENT ANALYSIS (PCA) & PCA-LOADING MAP
# ---------------------------------------------------------------------
# Compute PCA
pca_result <- prcomp(scaled_features, center = FALSE, scale. = FALSE)

# Prepare PCA Plot Data
pca_data <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Cluster = df_clean$Cluster_Label,
  Status = df_clean$Pocket_Distortion_Status
)

# Axis Labels (Variance Explained)
var_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100
pc1_label <- sprintf("PC1 (%.1f%% Variance)", var_explained[1])
pc2_label <- sprintf("PC2 (%.1f%% Variance)", var_explained[2])

# Pearson Correlation check for PC2 as the "Divergence Axis"
message(sprintf("Pearson Correlation (Divergence vs PC2): %.4f", 
                cor(df_clean$Scaffold_Pocket_Divergence, pca_data$PC2)))

# Generate the PCA Biplot
# Colors: Blue = Transitional, Green = Consonance, Red = Dissonance
cluster_colors <- c("#3498db", "#2ecc71", "#e74c3c") 

plot(pca_data$PC1, pca_data$PC2,
     col = cluster_colors[pca_data$Cluster],
     pch = 19, cex = 1.2, xlim = c(-2.8, 7.5),
     main = "PCA Diagnostic Projection: Structural Consonance vs. Dissonance",
     xlab = pc1_label, ylab = pc2_label)
grid(col = "lightgray", lty = "dotted")

# PCA LOADING VECTORS (Arrows)
print("=== PCA LOADING MATRIX ===")
loadings <- pca_result$rotation
print(loadings)

arrow_scale <- 1.5
for (i in 1:nrow(loadings)) {
  arrows(x0 = 0, y0 = 0, 
         x1 = loadings[i, 1] * arrow_scale, 
         y1 = loadings[i, 2] * arrow_scale, 
         col = "black", lwd = 2.0, length = 0.08)
}

# SPACE-OPTIMIZED VECTOR LABELS
text(loadings["Global_Mean_pLDDT", 1] * arrow_scale, 
     loadings["Global_Mean_pLDDT", 2] * arrow_scale, 
     labels = "Global_Mean", col = "black", font = 2, cex = 0.7, pos = 3)

text(loadings["Pocket_Mean_pLDDT", 1] * arrow_scale - 0.2, 
     loadings["Pocket_Mean_pLDDT", 2] * arrow_scale + 0.15, 
     labels = "Pocket_Mean", col = "black", font = 2, cex = 0.7, pos = 4)

text(loadings["Pocket_Min_pLDDT", 1] * arrow_scale - 0.2, 
     loadings["Pocket_Min_pLDDT", 2] * arrow_scale - 0.15, 
     labels = "Pocket_Min", col = "black", font = 2, cex = 0.7, pos = 4)

legend("topright", 
       legend = c("Transitional Stability", "Structural Consonance", "Structural Dissonance"),
       col = cluster_colors, pch = 19, cex = 0.8, bty = "o", bg = "white")

# Save a high-resolution copy to a .png file: 
dev.copy(png, "cyp_pca_diagnostic_plot.png", width=1100, height=800, res=150)
dev.off()

# ---------------------------------------------------------------------
# Section 6: LOGISTIC REGRESSION (Risk Calibration)
# ---------------------------------------------------------------------
# Prepare numeric target (0/1) for GLM
df_clean$Pocket_Distortion_Status_Numeric <- as.numeric(as.character(df_clean$Pocket_Distortion_Status))

# Define model: Can pLDDT metrics predict structural collapse?
logit_model <- glm(Pocket_Distortion_Status_Numeric ~ Pocket_Mean_pLDDT + Pocket_Min_pLDDT + Global_Mean_pLDDT, 
                   data = df_clean, family = "binomial")

print("=== LOGISTIC REGRESSION SUMMARY ===")
print(summary(logit_model))

# Calculate Odds Ratios (Risk Multipliers)
print("=== CALCULATED ODDS RATIOS (RISK PER UNIT pLDDT) ===")
odds_ratios <- exp(coef(logit_model))
print(odds_ratios)

# Confidence Intervals for Manuscript Table
print("=== 95% CONFIDENCE INTERVALS ===")
print(exp(confint(logit_model)))
