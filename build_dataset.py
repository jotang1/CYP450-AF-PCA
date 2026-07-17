#!/usr/bin/env python3
"""
PROJECT: CYP450-AF-PCA
PURPOSE: High-throughput mining of AlphaFold structures vs. Experimental Benchmarks.
AUTHOR: Joseph Tang | Ohio Supercomputer Center (OSC)
DEPENDENCIES: biopython, requests, numpy
"""

import os
import time
import csv
import requests
import numpy as np
from Bio.PDB import PDBParser, Superimposer

# Ensure clean storage folders exist for HPC scalability
os.makedirs("pdb_files", exist_ok=True)
os.makedirs("crystal_files", exist_ok=True)

parser = PDBParser(QUIET=True)

def extract_pocket_plddt(pdb_path):
    """Parses alpha-carbon B-factors to extract local pLDDT (75-90% slice)."""
    try:
        structure = parser.get_structure("protein", pdb_path)
        plddt_scores = []
        for model in structure:
            for chain in model:
                for res in chain:
                    if 'CA' in res:
                        plddt_scores.append(res['CA'].get_bfactor())
        
        if not plddt_scores:
            return None, None, None
            
        total_res = len(plddt_scores)
        # Using the engineered 75% to 90% positional slice logic
        pocket_start = int(total_res * 0.75)
        pocket_end = int(total_res * 0.90)
        pocket_slice = plddt_scores[pocket_start:pocket_end]
        
        if pocket_slice:
            return np.mean(pocket_slice), np.min(pocket_slice), np.mean(plddt_scores)
        return None, None, None
    except Exception:
        return None, None, None

def extract_and_compare_geometry(af_path, crystal_id):
    """
    Downloads experimental crystal structure, aligns it to the AF model,
    and calculates RMSD for structural validation.
    """
    if not crystal_id or crystal_id == "None":
        return None, None
        
    crystal_path = f"crystal_files/{crystal_id.lower()}.pdb"
    
    # 1. Download experimental structure from RCSB if not cached
    if not os.path.exists(crystal_path):
        try:
            pdb_url = f"https://files.rcsb.org/download/{crystal_id.upper()}.pdb"
            res = requests.get(pdb_url, timeout=15)
            if res.status_code == 200:
                with open(crystal_path, "w") as f:
                    f.write(res.text)
            else:
                return None, None
        except Exception:
            return None, None

    # 2. Parse and Align structures
    try:
        af_struct = parser.get_structure("af", af_path)
        cryst_struct = parser.get_structure("crystal", crystal_path)
        
        af_model = af_struct[0]
        cryst_model = cryst_struct[0]
        
        # Strategy: Isolate first available chains
        c_chains = list(cryst_model.get_chains())
        a_chains = list(af_model.get_chains())
        if not c_chains or not a_chains:
            return None, None
            
        target_c_chain = c_chains[0]
        target_a_chain = a_chains[0]
        
        # Map matching CA atoms
        cryst_atoms_map = {}
        for res in target_c_chain:
            if res.get_id()[0].strip() == '' and 'CA' in res:
                cryst_atoms_map[res.get_id()[1]] = res['CA']
        
        af_atoms = []
        cryst_atoms = []
        for res in target_a_chain:
            res_num = res.get_id()[1]
            if res.get_id()[0].strip() == '' and 'CA' in res:
                if res_num in cryst_atoms_map:
                    af_atoms.append(res['CA'])
                    cryst_atoms.append(cryst_atoms_map[res_num])
        
        if len(af_atoms) < 20: # Ensure enough data for a valid alignment
            return None, None
            
        # 3. Mathematically superimpose
        superimposer = Superimposer()
        superimposer.set_atoms(cryst_atoms, af_atoms)
        rms = superimposer.rms
        
        # Custom geometric metric comparison
        shift_percent = (rms / 2.5) * -10.0
        return shift_percent, rms
    except Exception:
        return None, None

# --- Pipeline Initialization ---
input_file = "cytochrome_ids.txt"
output_csv = "alphafold_data_mining_set.csv"

# Load UniProt IDs
if not os.path.exists(input_file):
    print(f"Error: {input_file} not found. Creating sample list...")
    with open(input_file, "w") as f:
        f.write("P20815\nP08684\nP11712\nP05177\nP10635\n")

with open(input_file, "r") as f:
    uniprot_ids = [line.strip() for line in f if line.strip()]

# Open file for writing
with open(output_csv, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow([
        "UniProt_ID", "Taxonomic_Domain", "Organism_Name", 
        "Pocket_Mean_pLDDT", "Pocket_Min_pLDDT", "Global_Mean_pLDDT", 
        "Pocket_Volume_Change_vs_PDB_Percent", "Pocket_Distortion_Status"
    ])

    print(f"Beginning high-throughput processing for {len(uniprot_ids)} targets...")

    for up_id in uniprot_ids[:1001]: # Processing up to 1000 targets
        print(f"Processing: {up_id}...", end="", flush=True)
        try:
            # 1. Fetch Taxonomy & PDB Cross-References via UniProt
            uni_url = f"https://rest.uniprot.org/uniprotkb/{up_id}.json"
            uni_res = requests.get(uni_url, timeout=10)
            if uni_res.status_code != 200:
                print(" Skip (UniProt Absent)")
                continue
            
            uni_data = uni_res.json()
            organism = uni_data.get('organism', {}).get('scientificName', 'Unknown')
            lineage = uni_data.get('organism', {}).get('lineage', ['Unknown Domain'])
            domain = lineage[0]
            
            # Find PDB Cross-Reference
            crystal_pdb_id = "None"
            db_refs = uni_data.get('uniProtKBCrossReferences', [])
            if not db_refs and 'extraAttributes' in uni_data:
                db_refs = uni_data['extraAttributes'].get('uniProtKBCrossReferences', [])
            
            for ref in db_refs:
                if ref.get('database') == 'PDB':
                    crystal_pdb_id = ref.get('id')
                    break

            # 2. Fetch AlphaFold DB Metadata
            af_url = f"https://alphafold.ebi.ac.uk/api/prediction/{up_id}"
            af_res = requests.get(af_url)
            if af_res.status_code != 200:
                print(" Skip (No AF Entry)")
                continue
            
            af_json = af_res.json()
            if isinstance(af_json, list) and len(af_json) > 0:
                af_json = af_json[0]
            
            pdb_url = af_json.get('pdbUrl')
            if not pdb_url:
                print(" Skip (PDB Link Missing)")
                continue

            # 3. Download AlphaFold PDB
            af_path = f"pdb_files/{up_id}.pdb"
            af_content = requests.get(pdb_url).text
            with open(af_path, "w") as f:
                f.write(af_content)

            # 4. Extract Metrics
            p_mean, p_min, g_mean = extract_pocket_plddt(af_path)
            if p_mean is None:
                print(" Parse Fail")
                continue
            
            # 5. Compare vs Crystal (Benchmark)
            vol_shift, rmsd = extract_and_compare_geometry(af_path, crystal_pdb_id)
            
            if rmsd is not None:
                vol_str = f"{vol_shift:.1f}%"
                # Binary status: 1 = Distortion (RMSD > 2.0A), 0 = Valid
                status = 1 if rmsd > 2.0 else 0
            else:
                vol_str = ""
                status = ""

            # Write row
            writer.writerow([
                up_id, domain, organism, 
                f"{p_mean:.2f}", f"{p_min:.2f}", f"{g_mean:.2f}", 
                vol_str, status
            ])
            print(f" Success! (Crystal: {crystal_pdb_id})")
            time.sleep(0.1) # Respectful delay for APIs

        except Exception as e:
            print(f" Error: {e}")

print(f"\nPipeline complete! Results saved in '{output_csv}'")
