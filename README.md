# CYP450-AF-PCA: An Automated Diagnostic Pipeline for AlphaFold Structural Auditing
## Project Overview
CYP450-AF-PCA is a high-throughput computational framework designed to audit the structural reliability of AlphaFold 2 (AF2) predictions within the Cytochrome P450 (CYP450) superfamily.

While AlphaFold 2 (AF2) provides high-confidence global scaffolds, it frequently fails to reconcile the structural plasticity of the CYP450 active site with the volumetric requirements of the heme prosthetic group. This pipeline identifies "Structural Dissonance"—cases where high global confidence masks localized functional failures—using a combination of PCA-driven diagnostics and supervised machine learning.

## System Architecture & Pipeline
The pipeline is optimized for High-Performance Computing (HPC) environments (e.g., Ohio Supercomputer Center) and is divided into three functional layers:
1. **Data Mining & Benchmarking** (Python)
   * **get_ids.py**: Automated retrieval of 999 UniProtKB Cytochrome P450 identifiers.
   * **build_dataset.py**: Core engine for fetching AlphaFold coordinates and cross-referencing against experimental PDB crystal structures. It calculates local pLDDT metrics and geometric RMSD deviations.
2. **Unsupervised Auditing** (R)
   * **cyp_pca_pipeline.R**: Performs dimensionality reduction (PCA) to decouple global confidence (PC1) from local pocket stability (PC2). Includes Hierarchical Clustering (Ward’s D2) to map the structural taxonomy of the dataset.
3. **Supervised Predictive Modeling** (Python/Scikit-Learn)
   * **compare_cyp_models.py**: Benchmarks Random Forest, SVM, and KNN classifiers to identify the most accurate predictor of structural failure.
   * **cyp_ml_pipeline.py**: Implements the validated Random Forest model and extracts feature importance rankings.
   * **run_knn_discovery.py**: The "Discovery Engine" that partitions uncharacterized targets into functional tiers (Consonance vs. Dissonance).

## Key Findings (999 Target Scan)
| Structural Phenotype | Count (n=913) | Biophysical Interpretation |
| :--- | :---: | :--- |
| **Structural Consonance** | 463 | Validated high-fidelity models; scaffold/pocket alignment. |
| **Transitional Stability** | 326 | Moderate stability; requires conformational relaxation (MD). |
| **Scaffold-Pocket Divergence** | 124 | **Outlier Phenotype**: Likely structural hallucination. |

### Structural Outlier Phenotypes
Our PCA-KNN diagnostic identified 124 targets exhibiting **Structural Dissonance**, bifurcating into two distinct pathological trajectories:

* **Pathway A: Stiff Cages (n=113)**
  * **Metric**: Pocket pLDDT > Global pLDDT (Negative PC2 Divergence).
  * **Mechanism**: **Hallucinated Rigidity.** In the absence of the heme-thiolate group, the model over-stabilizes the apo-form active site, resulting in a hyper-rigid and non-physical geometry.

* **Pathway B: Loop Failures (n=11)**
  * **Metric**: Global pLDDT > Pocket pLDDT (Positive PC2 Divergence).
  * **Mechanism**: **Structural Insecurity.** The substrate recognition sites (SRS) exhibit localized distortion and low confidence, failing to define a stable conformational state without the steric cues of the missing cofactor.

The 10:1 ratio of Stiff Cages to Loop Failures suggests that AF2 exhibits a systematic bias toward hallucinated stability in the CYP450 active site.

## Getting Started
### Data Availability Note
  - For academic evaluation, a representative subset (n=20) is provided in sample_data.csv. 
  - The full benchmarking dataset (n=999) is available upon request.
### Prerequisites
  - Python 3.8+: biopython, requests, pandas, scikit-learn, numpy, matplotlib
  - R 4.0+: stats, graphics
### Usage
  1. Initialize Dataset: python build_dataset.py
  2. Run PCA Diagnostic: Rscript cyp_pca_pipeline.R
  3. Execute ML Discovery: python run_knn_discovery.py

### Author
Joseph Tang <br>
Ohio Supercomputer Center
