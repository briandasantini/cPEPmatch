# cPEPmatch

**Created by**: Brianda L. Santini

**Supervised by**: Prof. Dr. Martin Zacharias

**Affiliation**: Physics Department T38, Technical University of Munich, Garching, Germany.


Webserver version: [cPEPmatch](https://t38webservices.nat.tum.de/cpepmatch/), Cite: [Santini & Zacharias](https://pubmed.ncbi.nlm.nih.gov/33134275/)

## General Description

cPEPmatch is a Python tool designed to find cyclic peptides that mimic proteins and target their binding partners. The program consists of four steps:

1. **Cyclic Peptide Database Creation**:
   - Reads from a pre-constructed database of 400 PDB structures downloaded from the PDB databank.
   - Characterizes motifs by calculating the CA carbon distances along the whole structure.
   - Users can select the motif size (4-7), consecutive or not, and specific features such as cyclization type or presence of non-standard amino acids.

2. **Characterization of the Protein-Protein or Protein-Ligand Interface**:
   - The motif is characterized using the same CA distance matrix as the cyclic peptide database.
   - Users can select the interface cutoff in angstroms and an optional parameter for protein-specific residues for targeted matching.

3. **Backbone Match Motif**:
   - Matches the CA distance motifs of the protein interface and the cyclic peptide database.
   - Selects matches that fit within the given frmsd_threshold and returns a Distance RMSD value in angstroms.

4. **Superimpose and Mutate**:
   - Matches are read and superimposed to match the structure of the protein mimic.
   - Returns a new Fit RMSD value in angstroms and mutates matching side-chains to mimic the protein.
   - Main limitation: Modeller cannot process non-standard residues, so these peptides are excluded from the mutation step.
   - Outputs a match_list.txt file and all mutated or aligned structures.

5. **Results Analysis**:
   - Output table contains match number, pdb name, Dist-RMSD, matched cyclic peptide residues, matched protein residues, Fit-RMSD, and file names.
   - Notebook version allows structure visualization by changing the `match_to_view` value.

## Usage

### Required Parameters

- `pdb_name`: PDB file name containing both the protein and the target.
- `protein`: Protein chain name.
- `target`: Target chain name.
- `motif_size`: Number of CA carbons to match.
- `consecutive`: Motif type (consecutive or non-consecutive).
- `interface_cutoff`: Distance in angstroms to define contact interface residues.
- `frmsd_threshold`: Fit-RMSD value in angstroms for matching precision.
- `working_location`: Path to the protein-target pdb file location.
- `database_location`: Path to the cPEPmatch database.

### Optional Parameters

- `protein_specific_residues`: Residue numbers to match if known (e.g., from hot-spot analysis).
- `cyclization_type`: Specific cyclization type for cyclic peptides (e.g., 'head to tail').
- `exclude_non_standard`: Exclude non-standard amino acid containing peptides (default: False).

## Directory Structure
```
cPEPmatch
├── cPEPmatch.py
├── cpep_modules
│ ├── backbone_match.py
│ ├── cpep_database.py
│ ├── protein_target_characterization.py
│ ├── superimpose_mutate.py
│ └── init.py
├── database
│ ├── cyclo_pep.csv
│ └── *.pdb
├── lib
│ ├── old_databases
│ ├── test_system
│ └── update_database
├── README.md
└── requirements.txt
```

## Installation and Setup

## Prerequisites
- Anaconda (recommended over Miniconda for compatibility reasons)

Creating a Conda Environment
Create a dedicated Conda environment for cPEPmatch:
```bash
conda create -n "cpepmatch" python=3.7.11

conda activate cpepmatch
``` 
## Install Required Python Packages
Install the necessary Python packages using the requirements.txt file:
```bash
  pip install -r requirements.txt
``` 
## Install Additional Conda Modules
Some modules need to be installed specifically via Conda:
 ```bash
 conda install -c conda-forge vmd.python

conda config --add channels salilab

conda install modeller
```
However, it will instruct you to modify /home/user/miniconda3/envs/server/lib/modeller-10.3/modlib/modeller/config.py and replace XXXX with your Modeller licence key (in this case: MODELIRANJE). Execute the command again. </br>
If the installation of vmd-python doesn't work, go to https://github.com/Eigenstate/vmd-python, download the repo, and run it.

VMD-Python Installation Issues
If you encounter issues installing vmd-python, visit vmd-python GitHub for alternative installation methods.

## Making the Script Executable and Setting Up PATH

1. Update the Shebang Line
Make sure the cpepmatch environment is activated
Find the Path to the Python Interpreter:
Use the which command to find out the path to the Python interpreter in the activated environment:
```bash
which python
```
This command will output the path to the Python interpreter, something like:
/home/user/anaconda3/envs/cpepmatch/bin/python

Edit the cPEPmatch.py file to include the shebang line with the path to the Anaconda Python interpreter at the very top:
```bash
#!/home/user/anaconda3/envs/cpepmatch/bin/python
```
Replace the path with the one you obtained from the which python command.


2. Make cPEPmatch Executable
Navigate to the cPEPmatch directory and make the main script executable:
```bash
chmod +x cPEPmatch.py
```

3. Add cPEPmatch to PATH
For convenience, add cPEPmatch to your system's PATH:
```bash
echo 'export PATH="$PATH:/path/to/cPEPmatch"' >> ~/.bashrc
source ~/.bashrc
```

## Usage
After installation, you can run cPEPmatch from any directory in your terminal:
```bash
cPEPmatch.py [arguments]
```
Replace [arguments] with appropriate command-line arguments for the script.

To learn about the various command-line arguments and options that cPEPmatch supports, you can use the help option:
```bash
cPEPmatch.py -h
```
