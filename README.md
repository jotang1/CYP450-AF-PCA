# CYP450-AF-PCA
A PCA-driven diagnostic pipeline for auditing structural consonance and functional viability in AlphaFold-predicted Cytochrome P450 enzymes.

1. Structural Consonance Diagnostic Pipeline
   AlphaFold (AF2/AF3) produces high-confidence scaffolds but can 'hallucinate' pocket geometries in cofactor-dependent enzymes. This pipeline uses Principal Component Analysis (PCA) to decouple global structural confidence (PC1) from local pocket divergence (PC2), identifying 124 distinct outlier phenotypes (Loop Failures and Stiff Cages).

2. Key Features
   - HPC-Ready: Optimized for high-throughput scanning of thousands of PDB files.
   - Automated Triage: Uses KNN-clustering to categorize targets into Structural Consonance or Structural Dissonance (Pathway A/B).
   - Multi-Language: Python-based pLDDT extraction + R-based statistical visualization.

3. The Diagnostic Framework
   - Structural Consonance: High global/local agreement (n=463).
   - Pathway A (Loop Failure): Low local confidence despite stable scaffold.
   - Pathway B (Stiff Cage): High confidence but biophysically rigid/inaccessible pocket.

4. Installation & Usage
   Include a requirements.txt for Python and a list of R libraries (ggplot2, factoextra, class).
  
5. System Architecture
   Processed 993 targets in XX seconds on OSC clusters.

## cyp_pca_pipeline.R | Structural Consonance Diagnostic Script
1. Overview
   This R-based analytical pipeline is designed to perform dimensionality reduction (PCA) and unsupervised clustering (Hierarchical/K-Means) on AlphaFold- predicted protein structures. Specifically, it identifies "Structural Dissonance"—cases where global scaffold confidence (pLDDT) remains high while the functional active site (substrate recognition sites) exhibits structural collapse or physical impossibility.

2. Theoretical Framework
   The script utilizes a dual-axis Principal Component Analysis (PCA) to audit AlphaFold outputs:
   - PC1 (Global Confidence): Captures the overall folding stability of the protein scaffold.
   - PC2 (Scaffold-to-Pocket Divergence): Maps the decoupling of local active-site confidence from the global fold.
   By projecting 999 Cytochrome P450 targets into this space, we categorize models into three distinct phenotypes:
   - Structural Consonance: High-fidelity models where the scaffold and pocket are in alignment.
   - Pathway A (Loop Failure): High-confidence scaffolds with collapsed/volatile substrate recognition loops.
   - Pathway B (Stiff Cage): Overly rigid, non-productive geometries that fail to accommodate heme-thiolate biophysics.

3. Key Technical Features
   - Automated Data Sanitization: Handles raw CSV exports, cleans percentage strings, and filters missing structural benchmarks.
   - Multi-Model Clustering: Executes both Ward’s D2 Hierarchical Clustering (for structural taxonomy) and K-Means (for phenotype partitioning).
   - Publication-Quality Artifacts: Automatically generates high-resolution .png and .pdf diagnostic plots, including PCA loading vectors (Biplots) and cluster  dendrograms.
   - Statistical Risk Assessment: Runs a logistic regression model to calculate Odds Ratios for structural distortion based on pLDDT variance.

4. Prerequisites & Installation
   The script is written in base R to ensure maximum portability across HPC environments (e.g., OSC Pitzer/Owens clusters) with minimal dependency bloat.
   Required R Libraries:
   - stats (Standard)
   - graphics (Standard)
   - grDevices (Standard)

5. Usage
   (1) Ensure the dataset (e.g., alphafold_data_mining_set_XXX.csv) is in the working directory.
   (2) Run the script via RStudio or the terminal:
           Rscript cyp_pca_pipeline.R
   (3) The pipeline will generate the following outputs:
       - cyp_pca_diagnostic_plot.png: The primary PCA cluster map.
       - dendrogram_structural_conformations.png: The hierarchical structural tree.
       - logistic_regression_results.csv: A table of Odds Ratios for manuscript reporting.

6. Data Availability Note
   The current repository includes a sample_data.csv containing 20 representative targets for code verification. The full benchmarking dataset (n=999) is  currently withheld and will be released upon peer-reviewed publication of the accompanying manuscript.
