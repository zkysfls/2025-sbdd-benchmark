import yaml
from rdkit import Chem
from rdkit.Chem import AllChem
import meeko
#import vina
import os
import csv
from pathlib import Path
import time
from tdc import Oracle

#import yaml



# models = ['dst', 'graph_ga', 'GraphAF', 'jt_vae', 'mimosa', 'moldqn', 'pasithea', 'reinvent',
#         'smiles_ga', 'smiles_lstm_hc', 'smiles_vae']

def process_txt_file(txt_path):
    """
    Process a single txt file containing SMILES strings
    Returns only non-empty lines
    """
    try:
        with open(txt_path, 'r') as f:
            data = f.readlines()
            # Filter out empty lines and strip whitespace
            return [line.strip() for line in data if line.strip()]
    except Exception as e:
        print(f"Error processing {txt_path}: {str(e)}")
        return []  # Return empty list instead of dict if error
    

def process_all_folders(root_directory):
    """
    Recursively process all folders and extract SMILES data from files
    """
   
    sa_oracle = Oracle(name = 'SA')
    qed_oracle = Oracle(name = 'QED')
    logp_oracle = Oracle(name = 'LogP')

    # Walk through all directories
    for model_dir in os.listdir(root_directory + '/Generated_molecules_archeive'):
        model_path = os.path.join(root_directory + '/Generated_molecules_archeive', model_dir)
       
        # Skip if not a directory
        if not os.path.isdir(model_path):
            continue
           
        print(f"Processing model directory: {model_dir}")
       
        # Walk through target directories
        for target_dir in os.listdir(model_path):
            target_path = os.path.join(model_path, target_dir)
           
            # Skip if not a directory
            if not os.path.isdir(target_path):
                continue

            # Check if SMILES.txt exists
            smiles_path = os.path.join(target_path, 'SMILES.txt')
            if not os.path.exists(smiles_path):
                print(f"SMILES.txt not found in {target_path} - skipping")
                continue
                
            # Ensure output directory exists
            csv_dir = '/data2/kai/SBDD/results/oracle_csv_files'
            os.makedirs(csv_dir, exist_ok=True)

            # Create CSV file with headers
            csv_filename = os.path.join(csv_dir, f"{model_dir}_{target_dir}.csv")
            with open(csv_filename, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(['Input smiles', 'SA', 'QED', 'LogP'])
                
                # Process the YAML file
                data = process_txt_file(smiles_path)
               
                # Store the data with model and target information
                for smiles in data:
                    try:
                        sa = sa_oracle(smiles)
                        qed = qed_oracle(smiles)
                        logp = logp_oracle(smiles)
                        writer.writerow([smiles, sa, qed, logp])
                    except Exception as e:
                        print(f"Error processing SMILES {smiles}: {str(e)}")
                        continue

def main():
    # Replace with your root directory path
    root_directory = '/data2/kai/SBDD/results'
   
    # Process all folders
    print("Starting to process folders...")
    process_all_folders(root_directory)

if __name__ == '__main__':
    main()                 

    # To run this script:
    # 1. Ensure you have Python 3 installed
    # 2. Install required packages: pip install meeko vina pandas
    # 3. Make sure you have the proper directory structure:
    #    - sbdd-improvement/ (this script's location)
    #    - ../vina-prepared-receptors/ (containing receptor PDBQT files and config files)
    #    - ../Pocket2Mol/ or ../3D-Generative-SBDD2/ (containing generated molecules)
    # 4. Run the script: python vina_test.py
            