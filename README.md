# CYP450-AF-PCA
A PCA-driven diagnostic pipeline for auditing structural consonance and functional viability in AlphaFold-predicted Cytochrome P450 enzymes.

1. Structural Consonance Diagnostic Pipeline
   AlphaFold (AF2/AF3) produces high-confidence scaffolds but can 'hallucinate' pocket geometries in cofactor-dependent enzymes. This pipeline uses Principal Component Analysis (PCA) to decouple global structural confidence (PC1) from local pocket divergence (PC2), identifying 124 distinct outlier phenotypes (Loop Failures and Stiff Cages).

2. Key Features
   HPC-Ready: Optimized for high-throughput scanning of thousands of PDB files.
   Automated Triage: Uses KNN-clustering to categorize targets into Structural Consonance or Structural Dissonance (Pathway A/B).
   Multi-Language: Python-based pLDDT extraction + R-based statistical visualization.

3. The Diagnostic Framework
   Structural Consonance: High global/local agreement (n=463).
   Pathway A (Loop Failure): Low local confidence despite stable scaffold.
   Pathway B (Stiff Cage): High confidence but biophysically rigid/inaccessible pocket.

4. Installation & Usage
   Include a requirements.txt for Python and a list of R libraries (ggplot2, factoextra, class).
  
5. System Architecture
   Processed 993 targets in XX seconds on OSC clusters.
