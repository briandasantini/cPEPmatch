#!/home/user/anaconda3/envs/cpepmatch/bin/python
# # cPEPmatch
# Python program designed to find cyclic peptides that mimic proteins and target their binding partners.

#     Created by Brianda L. Santini  
#     Supervised by Prof. Dr. Martin Zacharias
#     Physics Department T38, Technical University of Munich, Garching, Germany.

import os
import sys
import argparse


#cPEPmatch modules
from cpep_modules import cpep_database
from cpep_modules import protein_target_characterization
from cpep_modules import backbone_match
from cpep_modules import superimpose_mutate


def update_progress(progress):
    bar_length = 50  
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
    if progress < 0:
        progress = 0
    if progress >= 1:
        progress = 1

    block = int(round(bar_length * progress))
    text = "\rProgress: [{0}] {1:.1f}%".format("#" * block + "-" * (bar_length - block), progress * 100)
    sys.stdout.write(text)
    sys.stdout.flush()


def main(args):

    print("\n\ncPEPmatch process initiated.\n")
    # 1. Create cPEP library
    print("\nStep 1: Updating the cPEP database.")
    database = cpep_database.create_cyclolib_database(args.database_location, args.motif_size, args.consecutive, 
                                        args.cyclization_type, args.exclude_non_standard,
                                        csv_database_file='cyclo_pep.csv')
    update_progress(0.25)
    print("\nStep 1 complete.\n")

    print("\nStep 2: Characterizing the protein interface.")
    # 2. Create and Characterize Protein to Mimic
    interface_residues = protein_target_characterization.create_interfaces(args.working_location, args.protein, args.target, args.pdb_name, 
                                        args.interface_cutoff, residues=args.protein_specific_residues)
    protein_motifs = protein_target_characterization.process_interface(args.working_location, args.pdb_name, args.motif_size, motif_type=args.consecutive)
    update_progress(0.5)
    print("\nStep 2 complete.\n")

    print("\nStep 3: Finding matches. \nNote: This step might take longer for non-consecutive motifs.")
    # 3. Find Matches
    all_match = backbone_match.find_match(args.database_location, protein_motifs, args.motif_size, args.consecutive, args.frmsd_threshold)
    all_match.sort()
    update_progress(0.75)
    print("\nStep 3 complete.\n")
    
    print("\nStep 4: Finalizing - mutating and refining your matches.")
    original_stdout = sys.stdout
    with open('md.log', 'w') as f:
        sys.stdout = f
        real_match_list = superimpose_mutate.superimpose_and_mutate(args.working_location, all_match, args.motif_size, 
                                                args.interface_cutoff, args.frmsd_threshold, args.database_location, 
                                                args.protein, args.pdb_name, args.target)
        sys.stdout = original_stdout 

    update_progress(1.0)

    print("\nYay, all steps completed successfully! \n\n\n")
    try:
        with open(args.working_location + 'match_list.txt', 'r') as file:
            for line in file:
                print(line, end='')
    except IOError:
        print(f"Check your working directory for your results.")
    
    print("\n \nYour results are now available in the following directory: {}" .format(args.working_location))
    print("\n\nThank you for using cPEPmatch!\n\n")

if __name__ == "__main__":

    script_dir = os.path.dirname(os.path.realpath(__file__))
    default_database_location = os.path.join(script_dir, 'database/')
    default_test_system_location = os.path.join(script_dir, 'lib', 'test_system/')
    
    
    parser = argparse.ArgumentParser(description="cPEPmatch: Python program designed to find cyclic peptides that mimic proteins and target their binding partners.") #The program consists of four steps: \n\n1: Cyclic Peptide Database Creation\nThe program reads from a pre-constructed database of 400 PDB structures downloaded from the PDBdatabank. It characterizes the motifs by calculating the CA carbon distances along the whole structure. Users can select the motif_size (4-7), consecutive (consecutive or not consecutive), and whether to select specific features such as cyclization type (e.g., 'head to tail') or whether they contain non-standard amino acids or not. \n\n2: Characterization of the Protein-Protein or Protein-Ligand Interface  \nThe motif is characterized using the same CA distance matrix as the cyclic peptide database, with the same given parameters (motif_size, motif_type). Users can select the interface_cutoff which defines the distance in angstroms to select the contact interface residues that will be matched. An optional parameter is protein_specific_residues, which takes a given selection of residue numbers to match, for example, if hot-spot analysis has been done prior to the matching.\n\n3: Backbone Match Motif\nThis section matches the CA distance motifs of the protein interface and the cyclic peptide database. It selects matches that fit within the frmsd_threshold given, and returns a Distance RMSD value in angstroms.\n\n4: Superimpose and Mutate\nAll matches are read and superimposed to match the structure of the protein mimic. A new Fit RMSD value in angstroms is returned that determines how well the structures' overall backbones fit. Matching side-chains, except for disulphide bonded cysteines, are mutated into the matching residues of the protein to mimic. Modeller is used for side-chain mutation and optimization. A main limitation of this step is that modeller is unable to process non-standard residues, so all non-standard cyclic peptides are excluded from the mutation step. A match_list.txt file is outputted, along with all the mutated (match_pdbname.pdb) or aligned (match_pdbname-not_mutated.pdb) structures.""")

    parser.add_argument("-n", "--pdb_name", default='system', help="PDB file name, should contain both the protein and the target")
    parser.add_argument("-p", "--protein", default='A', help="Protein to mimic chain name")
    parser.add_argument("-t", "--target", default='B', help="Target chain name")
    parser.add_argument("-wl", "--working_location", default=default_test_system_location, help="Path to where your protein-target pdb file is located. Results will be outputted here")
    parser.add_argument("-dl", "--database_location", default=default_database_location, help="Path to where the cPEPmatch database is located")
   
    parser.add_argument("-ic", "--interface_cutoff", type=int, default=6, help="Distance in angstroms to define the contact interface residues that will be matched")
    parser.add_argument("-ft", "--frmsd_threshold", type=float, default=0.5, help="Fit-RMSD value in angstroms of how precise you want the matching")
    parser.add_argument("-ms", "--motif_size", type=int, default=5, help="Number of CA carbons to match")
    parser.add_argument("-cs", "--consecutive", type=lambda x: (str(x).lower() == 'true'), default=True, help="Motif type: consecutive or non-consecutive motifs")
    parser.add_argument("-psr", "--protein_specific_residues", default='', help="Selection of given residue numbers to match, for example, if hot-spot analysis has been done prior to the matching")
   
    parser.add_argument("-ct", "--cyclization_type", default='', help="Select only cyclic peptides from the database that are cyclized in a specific type, e.g., 'head to tail'")
    parser.add_argument("-ens", "--exclude_non_standard", type=lambda x: (str(x).lower() == 'true'), default=False, help="Exclude cyclic peptides from the database containing non-standard amino acids.")

    args = parser.parse_args()
    main(args)
