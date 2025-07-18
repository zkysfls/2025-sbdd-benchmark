import yaml
from rdkit import Chem
from rdkit.Chem import AllChem
import meeko
import vina
import os
import csv
from pathlib import Path



def process_yaml_file(yaml_path):
    """
    Process a single YAML file containing SMILES strings
    """
    try:
        with open(yaml_path, 'r') as f:
            data = yaml.safe_load(f)
            return list(data.keys())
    except Exception as e:
        print(f"Error processing {yaml_path}: {str(e)}")
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
    # print(f"Found {len(poses)} poses within the energy range.")

    
    return energies[0], poses

def process_molecules(root_directory, oracle_dir, model_name):
    # Path("../results").mkdir(parents=True, exist_ok=True)

    # Process only the specified model directory
    model_path = os.path.join(root_directory, model_name)
    if not os.path.isdir(model_path):
        print(f"Model directory not found: {model_name}")
        return

    print(f"Processing model directory: {model_name}")
    
    # Create output directory for this model
    output_dir = os.path.join("..", "results", f"{model_name}_vina_pose")
    os.makedirs(output_dir, exist_ok=True)

    # Process targets in specified model directory
    for target_dir in os.listdir(model_path):
        target_path = os.path.join(model_path, target_dir)
        if not os.path.isdir(target_path):
            continue

        # Process yaml file for this target
        yaml_file = next((f for f in os.listdir(target_path) if f.endswith('.yaml')), None)
        if not yaml_file:
            print(f"No YAML file found in {target_path}")
            continue
            
        yaml_path = os.path.join(target_path, yaml_file)
        smiles_list = process_yaml_file(yaml_path)
        
        # Find corresponding oracle target
        oracle_target_path = os.path.join(oracle_dir, target_dir)
        if not os.path.exists(oracle_target_path):
            print(f"Oracle target not found: {target_dir}")
            continue
        
        # Find receptor and config
        receptor_pdbqt, config_file = None, None
        for f in os.listdir(oracle_target_path):
            if f.endswith('.pdbqt'):
                receptor_pdbqt = os.path.join(oracle_target_path, f)
            elif f.endswith('.txt'):
                config_file = os.path.join(oracle_target_path, f)
        
        if not (receptor_pdbqt and config_file):
            print(f"Missing files in {oracle_target_path}")
            continue

        # Read config once per target
        try:
            center, box_size, exhaustiveness, energy_range = read_config(config_file)
        except Exception as e:
            print(f"Error reading config {config_file}: {str(e)}")
            continue

        # Prepare output files for this target
        docking_csv = os.path.join(output_dir, f"{model_name}_{target_dir}_docking_score.csv")
        poses_file = os.path.join(output_dir, f"{model_name}_{target_dir}_poses.pdbqt")
        
        # Open files in append mode to write results incrementally
        with open(docking_csv, 'w', newline='') as csvfile, open(poses_file, 'w') as posefile:
            writer = csv.writer(csvfile)
            writer.writerow(['SMILES', 'Docking Score'])
            
            # Process molecules one by one
            for i, smiles in enumerate(smiles_list):
                try:
                    lig_pdbqt = prepare_ligand(smiles)
                    energies, poses = dock_molecule(
                        receptor_pdbqt, lig_pdbqt, 
                        center, box_size, 
                        energy_range, exhaustiveness
                    )
                    
                    # Write results immediately
                    writer.writerow([smiles, energies[0]])
                    
                    if poses:
                        posefile.write(poses + "\n")
                    posefile.flush()  # Ensure immediate write
                    
                except Exception as e:
                    print(f"Error processing {smiles}: {str(e)}")
                    writer.writerow([smiles, 'ERROR'])
                    posefile.write(f"ERROR: {smiles}\n")
                
                # Print progress every 10 molecules
                if (i + 1) % 10 == 0:
                    print(f"Processed {i+1}/{len(smiles_list)} molecules for {target_dir}")

def main():
    import sys
    if len(sys.argv) < 2:
        print("Usage: python vina_test_model_pmo.py <model_name>")
        sys.exit(1)
        
    model_name = sys.argv[1]
    root_directory = '../'
    receptor_dir = './recetpor'
    process_molecules(root_directory, receptor_dir, model_name)       

if __name__ == '__main__':
    main()                 

            
            