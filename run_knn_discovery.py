#!/usr/bin/env python3
import os
import pandas as pd
from sklearn.neighbors import KNeighborsClassifier

# Load Data
if os.path.exists("alphafold_data_mining_set_1000.csv"):
    data_path = "alphafold_data_mining_set_1000.csv"
    print("System: Full 1000-target dataset detected.")
elif os.path.exists("sample_data.csv"):
    data_path = "sample_data.csv"
    print("System: Sample dataset detected for verification.")
else:
    print("Error: No valid dataset (Full or Sample) found in directory.")
    exit()

df = pd.read_csv(data_path)

# Filter Calibration vs Discovery sets
df_benchmarks = df.dropna(subset=['Pocket_Distortion_Status']).copy()
df_benchmarks['Pocket_Distortion_Status'] = df_benchmarks['Pocket_Distortion_Status'].astype(int)
df_discovery = df[df['Pocket_Distortion_Status'].isna() | (df['Pocket_Distortion_Status'] == "")].copy()

features = ['Pocket_Mean_pLDDT', 'Pocket_Min_pLDDT', 'Global_Mean_pLDDT']
X_train = df_benchmarks[features]
y_train = df_benchmarks['Pocket_Distortion_Status']
X_discover = df_discovery[features]

# Train KNN Model across all 86 benchmarks
print("Calibrating KNN model across 86 benchmarks")
knn = KNeighborsClassifier(n_neighbors=3, weights='distance')
knn.fit(X_train, y_train)

# Predict on the 913 Unknowns
df_discovery['KNN_Predicted_Class'] = knn.predict(X_discover)
df_discovery['KNN_Warp_Probability'] = knn.predict_proba(X_discover)[:, 1]

df_flags = df_discovery[df_discovery['KNN_Predicted_Class'] == 1]

# Export knn discovery report
output_csv = "alphafold_knn_discovery_predictions.csv"
df_discovery.to_csv(output_csv, index=False)

# =====================================================================
# DATA PARTITIONING (Replacing manual file reading with Pandas)
# =====================================================================
# Structural Consonance: Probability is 0.0 AND Predicted Class is 0
df_consonance = df_discovery[(df_discovery['KNN_Warp_Probability'] == 0.0) & (df_discovery['KNN_Predicted_Class'] == 0)]

# Scaffold Pocket Divergence: Predicted Class is 1 (distorted)
df_divergence = df_discovery[df_discovery['KNN_Predicted_Class'] == 1]

# Pathway A (Pocket_Mean > Global_Mean): Stiff_Cage
df_stiff_cages = df_divergence[df_divergence['Pocket_Mean_pLDDT'] > df_divergence['Global_Mean_pLDDT']]

# Pathway B: The remaining 11 - loop failure
df_loop_failures = df_divergence[df_divergence['Pocket_Mean_pLDDT'] <= df_divergence['Global_Mean_pLDDT']]

# Intermediate Stability: Everything else
moderate_stability = len(df_discovery) - len(df_consonance) - len(df_divergence)

# Export the sub-lists for your records
df_consonance.to_csv("structural_consonance.csv", index=False)
df_divergence.to_csv("scaffold_pocket_divergence.csv", index=False)
df_stiff_cages.to_csv("pathway_a_stiff_cages.csv", index=False)
df_loop_failures.to_csv("pathway_b_loop_failures.csv", index=False)

# =====================================================================
# FINAL DISCOVERY RESULTS OUTPUT
# =====================================================================
print("\n" + "="*40)
print("        KNN DISCOVERY PHASE RESULTS")
print("="*40)
print(f"Total Uncharacterized Proteins Scanned: {len(df_discovery)}")
print(f"Structural Consonance Proteins:         {len(df_consonance)}")
print(f"Intermediate Stability Proteins:        {moderate_stability}")
print(f"Scaffold Pocket Divergence Proteins:    {len(df_divergence)}")
print("-" * 50)
print(f"   --> Pathway A (Stiff Cages):         {len(df_stiff_cages)}")
print(f"   --> Pathway B (Loop Failures):       {len(df_loop_failures)}")
print("="*50)

print(f"\nFull discovery predictions exported to '{output_csv}'")
print(f"Phenotype reports exported to 'pathway_a_loop_failures.csv' and 'pathway_b_stiff_cages.csv'")
