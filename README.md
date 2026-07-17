# CYP450-AF-PCA: An Automated Diagnostic Pipeline for AlphaFold Structural Auditing
## Project Overview
CYP450-AF-PCA is a high-throughput computational framework designed to audit the structural reliability of AlphaFold 2 (AF2) predictions within the Cytochrome P450 superfamily.

While AF2 provides high-confidence global scaffolds, it frequently fails to reconcile the structural plasticity of the P450 active site with the volumetric requirements of the heme prosthetic group. This pipeline identifies "Structural Dissonance"—cases where high global confidence masks localized functional failures—using a combination of PCA-driven diagnostics and supervised machine learning.

## System Architecture & Pipeline
The pipeline is optimized for High-Performance Computing (HPC) environments (e.g., Ohio Supercomputer Center) and is divided into three functional layers:
1. Data Mining & Benchmarking (Python)
   - get_ids.py: Automated retrieval of 999 UniProtKB Cytochrome P450 identifiers.
   - build_dataset.py: Core engine for fetching AlphaFold coordinates and cross-referencing against experimental PDB crystal structures. It calculates local pLDDT metrics and geometric RMSD deviations.
2. Unsupervised Auditing (R)
   - cyp_pca_pipeline.R: Performs dimensionality reduction (PCA) to decouple global confidence (PC1) from local pocket stability (PC2). Includes Hierarchical Clustering (Ward’s D2) to map the structural taxonomy of the dataset.
3. Supervised Predictive Modeling (Python/Scikit-Learn)
   - compare_cyp_models.py: Benchmarks Random Forest, SVM, and KNN classifiers to identify the most accurate predictor of structural failure.
   - cyp_ml_pipeline.py: Implements the validated Random Forest model and extracts feature importance rankings.
   - run_knn_discovery.py: The "Discovery Engine" that partitions uncharacterized targets into functional tiers (Consonance vs. Dissonance).

## Key Findings (999 Target Scan)
| Structural Phenotype | Count (n=913) | Biophysical Interpretation |
| :--- | :---: | :--- |
| **Structural Consonance** | 463 | Validated high-fidelity models; scaffold/pocket alignment. |
| **Transitional Stability** | 326 | Moderate risk; requires conformational relaxation (MD). |
| **Scaffold-Pocket Divergence** | 124 | **Outlier Phenotype**: Likely structural hallucination. |

  - Pathway A (113 targets): Localized Loop Failures (SRS regions).
  - Pathway B (11 targets): "Stiff Cage" artifacts; biophysically non-productive geometries.

## Getting Started
### Data Availability Note
  - For academic evaluation, a representative subset (n=20) is provided in sample_data.csv. 
  - The full benchmarking dataset (n=999) is currently under research embargo pending manuscript submission.
### Prerequisites
  - Python 3.8+: biopython, requests, pandas, scikit-learn, numpy, matplotlib
  - R 4.0+: stats, graphics
### Usage
  1. Initialize Dataset: python build_dataset.py
  2. Run PCA Diagnostic: Rscript cyp_pca_pipeline.R
  3. Execute ML Discovery: python run_knn_discovery.py

### Author
Joseph Tang <br>
HPC Storage Engineer, Ohio Supercomputer Center <br>
Specializing in the intersection of HPC Systems and Structural Bioinformatics
