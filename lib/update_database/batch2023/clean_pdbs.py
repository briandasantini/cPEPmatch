import subprocess
import os
import shutil

def read_pdbs(file_path):
    """Reads PDB IDs from a given file."""
    with open(file_path, 'r') as file:
        return file.read().split(',')

def process_pdbs(pdbs, archive_dir):
    """Processes each PDB ID."""
    failed_pdbs = []  # List to keep track of failed PDB IDs

    for pdb_id in pdbs:
        try:
            # Define file names
            pdb_gz_file = f"{pdb_id}.pdb.gz"
            pdb_file = f"{pdb_id}.pdb"
            pdb_renum_file = f"{pdb_id}_renum.pdb"

            # Unzip the .pdb.gz file
            subprocess.run(f"gunzip -k {pdb_gz_file}", shell=True, check=True)

            # Run pdb4amber
            subprocess.run(f"pdb4amber -i {pdb_file} -o {pdb_renum_file} --nohyd", shell=True, check=True)

            # Move original files to archive directory
            if not os.path.exists(archive_dir):
                os.makedirs(archive_dir)
            shutil.move(pdb_gz_file, os.path.join(archive_dir, pdb_gz_file))
            shutil.move(pdb_file, os.path.join(archive_dir, pdb_file))

        except subprocess.CalledProcessError:
            print(f"Failed to process {pdb_id}")
            failed_pdbs.append(pdb_id)

    # Handling failed PDB IDs
    if failed_pdbs:
        with open('failed_pdbs.txt', 'w') as f:
            for pdb_id in failed_pdbs:
                f.write(f"{pdb_id}\n")

# Path to your pdbs.txt file
pdbs_file = 'pdbs.txt'  # Replace with the actual path to your pdbs.txt file
archive_dir = 'archive'  # Directory to store original files

# Process
pdbs = read_pdbs(pdbs_file)
process_pdbs(pdbs, archive_dir)

