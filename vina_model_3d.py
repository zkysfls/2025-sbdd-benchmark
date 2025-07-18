import yaml
from rdkit import Chem
from rdkit.Chem import AllChem
import meeko
import vina
import os
import csv
from pathlib import Path
import time

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
    
def read_config(config_file):
    """Read and parse docking configuration file"""
    with open(config_file, 'r') as f:
        config = dict(line.strip().split(' = ') for line in f if '=' in line)
    
    center = [float(config['center_x']), float(config['center_y']), float(config['center_z'])]
    box_size = [float(config['size_x']), float(config['size_y']), float(config['size_z'])]
    exhaustiveness = int(config['exhaustiveness'])
    energy_range = float(config['energy_range'])
    
    return center, box_size, exhaustiveness, energy_range

def prepare_ligand(smiles):
    """Prepare ligand for docking"""
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    
    meeko_prep = meeko.MoleculePreparation()
    meeko_prep.prepare(mol)
    lig_pdbqt = meeko_prep.write_pdbqt_string()
    
    return lig_pdbqt

def dock_molecule(receptor_pdbqt, lig_pdbqt, center, box_size, energy_range, exhaustiveness=8):
    """Perform docking using Vina"""
    v = vina.Vina(sf_name='vina', verbosity=0)
    v.set_receptor(receptor_pdbqt)
    v.set_ligand_from_string(lig_pdbqt)
    v.compute_vina_maps(center=center, box_size=box_size)
    
    # Perform docking
    v.dock(exhaustiveness=exhaustiveness, n_poses=20)
    
    poses = v.poses(n_poses=1, energy_range=energy_range)
    energies = v.energies(n_poses=1, energy_range=energy_range)
    # Extract only the first energy value from the array
    if len(energies) > 0:
        energies = energies[0][0]  # Get the first energy value from the first pose
    print(f"Found {len(poses)} poses within the energy range with energy: {energies}")

    
    return energies, poses

def process_molecules(root_directory, vina_prepared_receptors_dir):
    # Create log directory if not exists
    log_dir = '../results/logs'
    os.makedirs(log_dir, exist_ok=True)
    
    for model_dir in os.listdir(root_directory + '/Generated_molecules'):
        model_path = os.path.join(root_directory + '/Generated_molecules', model_dir)
        
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
                
            # Add file existence check FIRST
            smiles_path = os.path.join(target_path, 'SMILES.txt')
            if not os.path.exists(smiles_path):
                print(f"SMILES.txt not found in {target_path} - skipping")
                continue
            
            # Start timing only if SMILES.txt exists
            start_time = time.time()
            
            smiles_list = process_txt_file(smiles_path)
            
            receptor_dir_path = os.path.join(vina_prepared_receptors_dir, target_dir)
            if not os.path.isdir(receptor_dir_path):
                continue
            
            docking_results = []
            poses_results = []
            
            receptor_pdbqt = None
            config_file = None
            for f in os.listdir(receptor_dir_path):
                if f.endswith('.pdbqt'):
                    receptor_pdbqt = os.path.join(receptor_dir_path, f)
                elif f.endswith('.txt'):
                    config_file = os.path.join(receptor_dir_path, f) 
            
            if not (receptor_pdbqt and config_file):
                print(f"Missing required files in {receptor_dir_path}")
                continue

            print(f"Processing target: {target_dir}")
            print(f"Using receptor: {receptor_pdbqt}")
            print(f"Using config: {config_file}")

            # Read docking configuration
            center, box_size, exhaustiveness, energy_range = read_config(config_file)

            # Create individual log file path
            individual_log_path = os.path.join(log_dir, f"{model_dir}_{target_dir}.log")
            
            # Initialize log file with parameters
            with open(individual_log_path, 'w') as log_f:
                log_f.write(f"Model: {model_dir}\nTarget: {target_dir}\n")
                log_f.write(f"Receptor: {os.path.basename(receptor_pdbqt)}\n")
                log_f.write(f"Center: {center}\nBox size: {box_size}\n")
                log_f.write(f"Exhaustiveness: {exhaustiveness}\nEnergy range: {energy_range}\n")
                log_f.write(f"Total molecules: {len(smiles_list)}\n\n")

            for smiles in smiles_list:
                print(f"Processing smiles: {smiles}")
                try:
                    lig_pdbqt = prepare_ligand(smiles)
                    energies, poses = dock_molecule(receptor_pdbqt, lig_pdbqt, center, box_size, energy_range, exhaustiveness)
                    docking_results.append(energies)
                    poses_results.append(poses)
                except Exception as e:
                    with open(individual_log_path, 'a') as log_f:
                        log_f.write(f"Failed to process SMILES: {smiles} - {str(e)}\n")
            
            # Calculate elapsed time
            elapsed_time = time.time() - start_time
            
            # Update individual log with results
            with open(individual_log_path, 'a') as log_f:
                log_f.write(f"\nSuccessful dockings: {len(docking_results)}\n")
                log_f.write(f"Failed dockings: {len(smiles_list)-len(docking_results)}\n")
                log_f.write(f"Total time: {elapsed_time:.2f} seconds\n")
            
            if not os.path.exists('../results'):
                os.makedirs('../results')
            # Create separate folders for CSV and PDBQT files
            csv_dir = f"../results/csv_files"
            pdbqt_dir = f"../results/pdbqt_files"
            os.makedirs(csv_dir, exist_ok=True)
            os.makedirs(pdbqt_dir, exist_ok=True)
            
            docking_csv = f"{csv_dir}/{model_dir}_{target_dir}_docking_score.csv"
            poses_file = f"{pdbqt_dir}/{model_dir}_{target_dir}_poses.pdbqt"
            with open(docking_csv, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(['Input smiles', 'Docking score'])
                for smiles, score in zip(smiles_list, docking_results):
                    writer.writerow([smiles, score])
            with open(poses_file, 'w') as f:
                # Join all pose strings with newlines
                f.write('\n'.join(str(pose) for pose in poses_results))

def main():
    
    root_directory = '/data2/kai/SBDD/sbdd-improvement'
    vina_prepared_receptors_dir = './vina-prepared-receptors_copy'
    
    process_molecules(root_directory, vina_prepared_receptors_dir)       

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
            