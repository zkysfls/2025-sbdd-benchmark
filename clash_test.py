from posebusters import PoseBusters
from pathlib import Path
from posecheck import PoseCheck
import os
import pandas as pd

def process_pose_busters(root_directory):
    pc = PoseCheck()
    clash_check_dir = os.path.join(root_directory, 'clash_check_results')
    # Create output directory if not exists - add exist_ok for safety
    os.makedirs(clash_check_dir, exist_ok=True)
    
    # Fix potential path issues using pathlib
    receptors_dir = Path(root_directory) / 'vina-prepared-receptors_copy'
    sdf_dir = Path(root_directory) / 'pdbqt_files' / 'sdf_files'

    # Check if input directories exist
    if not receptors_dir.exists():
        raise FileNotFoundError(f"Receptors directory not found: {receptors_dir}")
    if not sdf_dir.exists():
        raise FileNotFoundError(f"SDF directory not found: {sdf_dir}")

    # Walk through target directories
    for target_dir in os.listdir(receptors_dir):
        target_path = os.path.join(receptors_dir, target_dir)
        
        # Skip if not a directory
        if not os.path.isdir(target_path):
            continue

        # Generate true file path
        true_file_path = os.path.join(target_path, f"{target_dir}.pdb")
        if not os.path.exists(true_file_path):
            print(f"Warning: Protein file not found: {true_file_path}")
            continue
            
        try:
            pc.load_protein_from_pdb(true_file_path)
        except Exception as e:
            print(f"Error loading protein {true_file_path}: {str(e)}")
            continue

        # Process each model's SDF files
        for model_dir in ["Pocket2Mol", "PocketFlow", "ResGen"]:
            # Generate pred file path
            pred_file_name = f"{model_dir}_{target_dir}_poses.sdf"
            pred_path = os.path.join(sdf_dir, pred_file_name)
            if not os.path.exists(pred_path):
                print(f"Warning: SDF file not found: {pred_path}")
                continue
                
            pred_file = Path(pred_path)
            try:
                pc.load_ligands_from_sdf(pred_path)
            except Exception as e:
                print(f"Error loading ligands from {pred_path}: {str(e)}")
                continue

            # Run PoseBusters with error handling
            print(f"Processing {model_dir} for target {target_dir}")
            try:
                # Check for clashes
                clashes = pc.calculate_clashes()
                if not clashes:
                    print(f"Warning: No clash data returned for {pred_path}")
                    continue
                print(f"Number of clashes per molecule: {clashes[0]}")

                # Check for strain
                strain = pc.calculate_strain_energy()
                if not strain:
                    print(f"Warning: No strain energy data returned for {pred_path}")
                    continue
                print(f"Strain energy per molecule: {strain[0]}")

                # Create DataFrame from results
                results = list(zip(clashes, strain))  # Pair corresponding entries
                df = pd.DataFrame(results, columns=['clashes', 'strain_energy'])
                
                # Save results if DataFrame is not empty
                if not df.empty:
                    output_file = Path(clash_check_dir) / f"{model_dir}_{target_dir}_clash_check_results.csv"
                    try:
                        df.to_csv(output_file)
                        print(f"Saved results to {output_file}")
                    except Exception as e:
                        print(f"Error saving results to {output_file}: {str(e)}")
                else:
                    print(f"No results generated for {pred_file}")
            except Exception as e:
                print(f"Error processing {pred_file}: {str(e)}")

def main():
    root_directory = "/data2/kai/SBDD/results"
    try:
        process_pose_busters(root_directory)
    except Exception as e:
        print(f"Fatal error in main process: {str(e)}")

if __name__ == "__main__":
    main()
