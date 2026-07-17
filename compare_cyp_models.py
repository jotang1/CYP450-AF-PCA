#!/usr/bin/env python3

import os
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score

# Load data:
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

# Filter out rows missing crystal results using clean column title
df_benchmarks = df.dropna(subset=['Pocket_Distortion_Status']).copy()
df_benchmarks['Pocket_Distortion_Status'] = df_benchmarks['Pocket_Distortion_Status'].astype(int)

features = ['Pocket_Mean_pLDDT', 'Pocket_Min_pLDDT', 'Global_Mean_pLDDT']
X = df_benchmarks[features]
y = df_benchmarks['Pocket_Distortion_Status']

# Stratified 80/20 Train-Test Split
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42, stratify=y
)

# Initialize 3 models: RF. SVM, and KNN
models = {
    "Random Forest": RandomForestClassifier(n_estimators=100, max_depth=5, class_weight= 'balanced', random_state=42),
    "Support Vector Machine (SVM)": SVC(kernel='linear', C=1.0, class_weight= 'balanced', random_state=42),
    "K-Nearest Neighbors (KNN)": KNeighborsClassifier(n_neighbors=5, weights= 'distance')
}

performance_records = []

# Train, predict, and extract metrics side-by-side
print("Training and benchmarking supervised classifiers across 86 metrics...")
for name, model in models.items():
    model.fit(X_train, y_train)
    predictions = model.predict(X_test)

    # Calculate performance metrics required by your grading criteria
    acc = accuracy_score(y_test, predictions)
    prec = precision_score(y_test, predictions, zero_division=0)
    rec = recall_score(y_test, predictions, zero_division=0)
    f1 = f1_score(y_test, predictions, zero_division=0)

    performance_records.append({
        "Model Name": name,
        "Accuracy": f"{acc:.2%}",
        "Precision": f"{prec:.2%}",
        "Recall": f"{rec:.2%}",
        "F1-Score": f"{f1:.2%}"
    })

# Print a comparison table to the terminal
summary_df = pd.DataFrame(performance_records)
print("\n" + "="*70)
print("             SUPERVISED MODEL COMPARISON MATRIX")
print("="*70)
print(summary_df.to_string(index=False))
print("="*70)

# Export comparison to a summary CSV for your report tables
summary_df.to_csv("supervised_model_comparison_matrix.csv", index=False)
print("\nSummary matrix exported cleanly to 'supervised_model_comparison_matrix.csv'")
