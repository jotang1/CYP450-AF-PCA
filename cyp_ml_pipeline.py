#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np

# Force matplotlib to print figures to files without a screen monitor
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score

# Load Data Matrix Table
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

# Ensure data types are uniform
features = ['Pocket_Mean_pLDDT', 'Pocket_Min_pLDDT', 'Global_Mean_pLDDT']

# =====================================================================
# PHASE 1 & 2: CALIBRATION & VALIDATION POOL (The 86 Benchmark Proteins)
# =====================================================================
# Filter out rows that have reported crystal structure answers
df_benchmarks = df.dropna(subset=['Pocket_Distortion_Status']).copy()
df_benchmarks['Pocket_Distortion_Status'] = df_benchmarks['Pocket_Distortion_Status'].astype(int)

X_benchlines = df_benchmarks[features]
y_benchlines = df_benchmarks['Pocket_Distortion_Status']

# Split the 86 proteins: 80% for Phase 1 Training, 20% for Phase 2 Validation Testing
X_train, X_val, y_train, y_val = train_test_split(
    X_benchlines, y_benchlines, test_size=0.2, random_state=42, stratify=y_benchlines
)

print("="*60)
print(f"PHASE 1: TRAINING MODEL ON {len(X_train)} CALIBRATION PROTEINS")
print("="*60)
clf = RandomForestClassifier(n_estimators=100, max_depth=5, random_state=42)
clf.fit(X_train, y_train)
print("Model calibration complete.")

print("\n" + "="*60)
print(f"PHASE 2: VALIDATING PERFORMANCE ON {len(X_val)} HIDDEN BENCHMARKS")
print("="*60)
# Test the model on data it has never seen before to calculate true course metrics
y_val_pred = clf.predict(X_val)

print(f"Validation Accuracy: {accuracy_score(y_val, y_val_pred):.2%}")
print("\nConfusion Matrix:")
print(confusion_matrix(y_val, y_val_pred))
print("\nClassification Report (Precision, Recall, F1):")
print(classification_report(y_val, y_val_pred, target_names=["Stable (0)", "Warped (1)"]))

# Calculate Feature Importances from the validated model
importances = clf.feature_importances_
feat_imp = pd.DataFrame({'Metric': features, 'Importance': importances}).sort_values(by='Importance', ascending=False)
print("\nFeature Importance Rankings:")
print(feat_imp.to_string(index=False))

# =====================================================================
# PHASE 3: THE DISCOVERY PHASE (The 913 Uncharacterized Proteins)
# =====================================================================
print("\n" + "="*60)
print(f"PHASE 3: DEPLOYING ENGINE TO SCAN {df['Pocket_Distortion_Status'].isna().sum()} UNCHARACTERIZED TARGETS")
print("="*60)

# Filter for the rows that are missing crystal structures
df_discovery = df[df['Pocket_Distortion_Status'].isna() | (df['Pocket_Distortion_Status'] == "")].copy()
X_discover = df_discovery[features]

# Use the trained model to predict the hidden states
df_discovery['Predicted_Class'] = clf.predict(X_discover)
df_discovery['Collapse_Probability'] = clf.predict_proba(X_discover)[:, 1]

df_flags = df_discovery[df_discovery['Predicted_Class'] == 1]
print(f"Scan Complete! Found {len(df_flags)} proteins predicted to be secretly warped/collapsed.")

# Export discovery results to a clean CSV report
output_csv = "alphafold_discovery_predictions.csv"
df_discovery.to_csv(output_csv, index=False)
print(f"Discovery report exported to '{output_csv}'")
print("="*60)
