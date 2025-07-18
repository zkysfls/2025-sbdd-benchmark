from posebusters import PoseBusters
from pathlib import Path
import os

def process_pose_busters(root_directory):
    pose_busters_dir = os.path.join(root_directory, 'pose_busters_results')
    # Create output directory if not exists - add exist_ok for safety
    os.makedirs(pose_busters_dir, exist_ok=True)
    
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
        cond_file = Path(true_file_path)
        # Process each model's SDF files
        for model_dir in ["Pocket2Mol", "PocketFlow", "ResGen"]:
            # Generate pred file path
            pred_file_name = f"{model_dir}_{target_dir}_poses.sdf"
            pred_path = os.path.join(sdf_dir, pred_file_name)
            if not os.path.exists(pred_path):
                continue
            pred_file = Path(pred_path)

            # Run PoseBusters with error handling
            print(f"Processing {model_dir} for target {target_dir}")
            try:
                buster = PoseBusters(config="dock")
                df = buster.bust([pred_file], None, cond_file)
                
                # Save results if DataFrame is not empty
                if not df.empty:
                    output_file = Path(pose_busters_dir) / f"{model_dir}_{target_dir}_pose_busters_results.csv"
                    df.to_csv(output_file)
                    print(f"Saved results to {output_file}")
                else:
                    print(f"No results generated for {pred_file}")
            except Exception as e:
                print(f"Error processing {pred_file}: {str(e)}")

def main():
    root_directory = "/data2/kai/SBDD/results"
    process_pose_busters(root_directory)

if __name__ == "__main__":
    main()
