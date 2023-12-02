#Python script to find duplicate PDB IDs

def read_pdbs(file_path):
    with open(file_path, 'r') as file:
        return file.read().split(',')

# Read PDB IDs from both files
batch23_pdb = read_pdbs('batch2019/pdbs.txt')
batch19_pdb = read_pdbs('batch_by-hand2023/pdbs.txt')

# Combine and find duplicates
all_pdbs = batch23_pdb + batch19_pdb
duplicates = set([pdb for pdb in all_pdbs if all_pdbs.count(pdb) > 1])

print("Duplicates:", duplicates)
