def read_duplicates(file_path):
    """Reads duplicate PDB IDs from a given file."""
    with open(file_path, 'r') as file:
        # Extracting the set of duplicates from the file
        duplicates_line = file.readline().strip()
        duplicates = set(eval(duplicates_line.split(":")[1].strip()))
        return duplicates

def filter_duplicates(batch_file, duplicates, output_file):
    """Filters out duplicate PDB IDs from a batch file and writes to a new file."""
    with open(batch_file, 'r') as file:
        ids = file.read().split(',')

    # Removing duplicates
    filtered_ids = [id for id in ids if id not in duplicates]

    with open(output_file, 'w') as file:
        file.write(','.join(filtered_ids))

# Paths to your files
duplicates_file = 'duplicates.txt'  # Update with the path to your duplicates file
batch23_file = 'batch_by-hand2023/pdbs.txt'               # Path to the batch23 pdbs file
output_file = 'batch_by-hand2023/pdbs_filtered.txt'       # Path for the new file without duplicates

# Processing
duplicates = read_duplicates(duplicates_file)
filter_duplicates(batch23_file, duplicates, output_file)

