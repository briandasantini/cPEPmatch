# superimpose_and_mutate.py
# Part of the cPEPmatch tool developed by Brianda L. Santini, supervised by Prof. Dr. Martin Zacharias
# Physics Department T38, Technical University of Munich, Garching, Germany.
#
# Module Description:
# This module is responsible for the final stages of the cPEPmatch process. It superimposes matched 
# structures from the cyclic peptide database onto protein targets and then performs mutations 
# on the aligned structures. This module's functions ensure that the output includes the superimposed 
# PDB files and a comprehensive list of matches, facilitating the analysis of potential peptide mimics.

import sys
import os
import re
import glob
import operator
import pickle
import time
import warnings
import shutil
import csv
import argparse
import numpy as np
import pandas as pd
from scipy import spatial
from collections import defaultdict
from itertools import product, combinations, permutations
from Bio.PDB import PDBParser, parse_pdb_header
from modeller import *
from modeller.optimizers import MolecularDynamics, ConjugateGradients
from modeller.automodel import autosched
import Bio
import vmd
from vmd import molecule, atomsel
    
# Functions for the superimposition and mutation processes

# Modeller Optimization

def optimize(atmsel, sched):
    #conjugate gradient
    for step in sched:
        step.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)
    #md
    refine(atmsel)
    cg = ConjugateGradients()
    cg.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)


# Molecular Dynamics

def refine(atmsel):
    # at T=1000, max_atom_shift for 4fs is cca 0.15 A.
    md = MolecularDynamics(cap_atom_shift=0.39, md_time_step=4.0,
                           md_return='FINAL')
    init_vel = True
    for (its, equil, temps) in ((200, 20, (150.0, 250.0, 400.0, 700.0, 1000.0)),
                                (200, 600,
                                 (1000.0, 800.0, 600.0, 500.0, 400.0, 300.0))):
        for temp in temps:
            md.optimize(atmsel, init_velocities=init_vel, temperature=temp,
                         max_iterations=its, equilibrate=equil)
            init_vel = False


def model_side_chains(match_path, modelname, match, chain, respos, restyp, mutation_count):  
    #Modeller's standard sidechain mutation (TODO:insert link)

    #log.verbose
    env = Environ(rand_seed=-49837)
    env.io.hetatm = True
    env.edat.dynamic_sphere=False
    env.edat.dynamic_lennard=True
    env.edat.contact_shell = 4.0
    env.edat.update_dynamic = 0.39
    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')
    mdl1 = Model(env, file=modelname)
    ali = Alignment(env)
    ali.append_model(mdl1, atom_files=modelname, align_codes=modelname)

    s = Selection(mdl1.chains[chain].residues[respos])
    s.mutate(residue_type=restyp)
    ali.append_model(mdl1, align_codes=modelname)
    mdl1.clear_topology()
    mdl1.generate_topology(ali[-1])

    mdl1.transfer_xyz(ali)
    mdl1.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')
    mdl2 = Model(env, file=modelname)

    mdl1.res_num_from(mdl2,ali)
    mdl1.write(file=modelname+restyp+respos+'.tmp')
    mdl1.read(file=modelname+restyp+respos+'.tmp')
    make_restraints(mdl1, ali)
    mdl1.env.edat.nonbonded_sel_atoms=1

    sched = autosched.loop.make_for_model(mdl1)
    s = Selection(mdl1.chains[chain].residues[respos])

    mdl1.restraints.unpick_all()
    mdl1.restraints.pick(s)

    s.energy()
    s.randomize_xyz(deviation=4.0)

    mdl1.env.edat.nonbonded_sel_atoms=2
    optimize(s, sched)
    mdl1.env.edat.nonbonded_sel_atoms=1
    optimize(s, sched)
    s.energy()

    #Output file
    mdl1.write(match_path + ('/mutation{}.pdb'.format(mutation_count)))
    #delete the temporary file
    os.remove(modelname+restyp+respos+'.tmp')

    shutil.copyfile((match_path + ('/mutation{}.pdb'.format(mutation_count))), 
                    (match_path + ('/mutated_match.pdb')))

    modelname = (match_path + ('/mutated_match.pdb'))
    return(modelname)
    


# Homologs and dihedral library for dihedral angle restraints

def make_restraints(mdl1, aln):
    rsr = mdl1.restraints
    rsr.clear()
    s = Selection(mdl1)
    for typ in ('stereo', 'phi-psi_binormal'):
        rsr.make(s, restraint_type=typ, aln=aln, spline_on_site=True)
    for typ in ('omega', 'chi1', 'chi2', 'chi3', 'chi4'):
        rsr.make(s, restraint_type=typ+'_dihedral', spline_range=4.0,
                 spline_dx=0.3, spline_min_points = 5, aln=aln,
                 spline_on_site=True)
      

# Superimpose and Mutate

def mutate_matched_side_chains(match_path, match, real_match_count, motif_size, ppi_chain, cyclo_chain):
    modelname = (match_path + '/cpep_aligned.pdb') 
    motif_length = len(match[2]) - 1
    mutation_count = 0
    chain = 'A'
    disulfide = 'CYX'

    for matched_residue in match[2]:
        respos = str(matched_residue)
        residue_number_to_match = match[3][mutation_count]
        restyp = ppi_chain[residue_number_to_match].get_resname()     
            
        try:
            cyclic_restyp = cyclo_chain[matched_residue].get_resname()
        except Exception:
            mutation_count = mutation_count + 1
            break
            
        #Do not mutate disulfide-bonded cysteines
    
        if cyclic_restyp != disulfide:
            try:
                #Standard Mutation
                modelname = model_side_chains(match_path, modelname, match, chain, respos, restyp, mutation_count)
                mutation_count = mutation_count + 1
            except:
                #Non-Standard Amino Acids
                mutation_count = mutation_count + 1
                pass 
        else:
            #Disulfide-bonded cysteines
            mutation_count = mutation_count + 1
    try:
        shutil.copyfile((match_path + ('/mutated_match.pdb')), 
                        ('match{}_{}.pdb'.format(real_match_count, match[0])))
        match_name = 'match{}_{}.pdb'.format(real_match_count, match[0])
    except Exception:
        shutil.copyfile((match_path + ('/cpep_aligned.pdb')), 
                        ('match{}_{}-NotMutated.pdb'.format(real_match_count, match[0])))
        match_name = 'match{}_{}-NotMutated.pdb'.format(real_match_count, match[0])

    real_match_count = real_match_count + 1 

    for fname in os.listdir(match_path):
        if fname.startswith("mutation"):
            os.remove(os.path.join(match_path, fname))
            
    return(real_match_count, modelname, match_name)



def write_match_list_header(motif_size, interface_cutoff,frmsd_threshold, working_location):
    os.chdir('{}'.format(working_location))
    with open('match_list.txt', 'w') as f:
        f.write("List of matches found by cPEPmatch.\nMotif Size: {} \nInterface Cutoff: {} \nFit-RMSD Threshold: {} \n\n" .format(motif_size, interface_cutoff,frmsd_threshold))
        f.write(" Match  PDB   Dist-RMSD      cPep Residues              PP-Interface Residues      Fit-RMSD\n\n")


# Function to convert a list of numbers to a string separated by a given separator
def int_list(lyst, sep=' ', quote=None):
    """Return the list of numbers in 'lyst' as a string separated by 'sep'."""
    f = str if quote is None else lambda i: quote+str(i)+quote
    return sep.join(map(f, lyst))


def remove_old_match_directories(working_location):
    dir_list = os.listdir(working_location)
    for dir_name in dir_list:
        if "match" in dir_name:
            dir_path = os.path.join(working_location, dir_name)
            if os.path.isdir(dir_path):
                try:
                    shutil.rmtree(dir_path)
                    print(f"Removed: {dir_path}")
                except Exception as e:
                    print(f"Error removing {dir_path}: {e}")



def eliminate_overlap(match, pdb_name, c1, working_location, match_path, target):
    no_overlap_match = []
    
    #Eliminates the matched cyclic peptides that overlap with the receptor.    
    cyclo_rlist = int_list(match[2], quote='"')
    matched_cyclo_residues = cyclo_rlist
    ppi_path = working_location + ("/{}.pdb".format(pdb_name))
    overlap = overlap_finder(match_path, ppi_path, matched_cyclo_residues, target)

    if overlap == [[],[]]:
        no_overlap_match = match

    #Write output file
    with open('match_list.txt', 'a') as f:
        if no_overlap_match != []:
            c2 = -1
            f.write("{:>6}".format(c1))
            for i in no_overlap_match:
                c2 = c2 + 1
                if c2 < 1:
                    f.write("{:>6} ".format(i))
                elif c2 < 2:
                    f.write("{:>9.4f}".format(i))
                elif c2 < 3:
                    f.write("   ")
                    for r in i:
                        f.write("{:>5}".format(r))
                elif c2 < 4:
                    f.write("   ")
                    for r in i:
                        f.write("{:>5}".format(r))
                else:
                    f.write("   ")
                    f.write("{:>9.4f}".format(i))
                
            f.write("\n")
            c1 = c1 + 1  
            
    return(no_overlap_match, c1)


def overlap_finder(match_path, complex_fname, matched_sequence, target):
    
    parser = PDBParser()
    ppi = molecule.load("pdb", complex_fname)    
    receptor_chain = atomsel("chain {}" .format(target)) #change this if the receptor has a different chain name.
    
    parser = PDBParser()
    cyclo = molecule.load("pdb", (match_path + '/cpep_aligned.pdb'))
    cyclo_chain = atomsel('chain A and not resid {}'.format(matched_sequence))
    
    overlap = atomsel.contacts(cyclo_chain, receptor_chain, 1)
    
    return(overlap)   


def superimpose(match_path, match, ppi_atoms, cyclo_atoms, cyclo_model, ppi_model, cyclo_structure, ppi_structure):
    # Superimposer:
    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ppi_atoms, cyclo_atoms)
    super_imposer.apply(cyclo_model.get_atoms())

    # Add to fit_RMSD the results:
    superimpose_rmsd = super_imposer.rms
    lst = list(match)
    lst = lst + [superimpose_rmsd]
    match = tuple(lst)

    # Save the aligned version 
    io = Bio.PDB.PDBIO()
    io.set_structure(cyclo_structure) 
    io.save(match_path + "/cpep_aligned.pdb")
    return(match)


def make_list_of_atoms_to_align(match, cyclo_model, ppi_model):
    cyclo_atoms = []
    ppi_atoms = []
    cyclo_residues_to_be_aligned = match[2]
    ppi_residues_to_be_aligned = match[3]

    for cyclo_chain in cyclo_model:
        for cyclo_res in cyclo_chain:
            if cyclo_res.get_id()[1] in cyclo_residues_to_be_aligned: # Check if residue number ( .get_id() ) is in the list
                cyclo_atoms.append(cyclo_res['CA'])


    for ppi_chain in ppi_model:
        for ppi_res in ppi_chain:
            if ppi_res.get_id()[1] in ppi_residues_to_be_aligned:
                ppi_atoms.append(ppi_res['CA'])
                
    return(cyclo_atoms, ppi_atoms)


def superimpose_and_mutate(working_location, all_match, motif_size, interface_cutoff, frmsd_threshold, database_location, protein, pdb_name, target):
    
    match_directory_list = []
    real_match_list= []
    all_mutations = dict()
    motif_end = motif_size - 1
    real_match_count = 1
    cr = 1
    c1 = 1
    
    remove_old_match_directories(working_location)           
    write_match_list_header(motif_size, interface_cutoff, frmsd_threshold, working_location)

    for match in all_match: 

        no_overlap_match = []
        pdb_parser = Bio.PDB.PDBParser(QUIET = True)
        matched_cyclopep = match[0]
        match_dir = ("match{}_{}".format(cr, matched_cyclopep))
        match_directory_list.append(match_dir)

        #Create a directory for each match
        try:
            os.mkdir("{}".format(match_dir))
            cr  =  cr + 1
        except RuntimeError as e:
            warn(e)
            pass

        cpep = ('{}{}.pdb'.format(database_location, matched_cyclopep))
        match_folder = (working_location + ('{}/cpep_match.pdb'.format(match_dir)))
        match_path = (working_location + ('{}'.format(match_dir)))

        try:
            shutil.copy(cpep, match_folder, follow_symlinks=True)
        except RuntimeError as e:
            warn(e)
            pass

        #Superimpose every match  
        cyclo_structure = pdb_parser.get_structure("reference", (working_location + ("{}/cpep_match.pdb".format(match_dir))))  
        ppi_structure = pdb_parser.get_structure("sample", (working_location + "interface.pdb"))

        cyclo_model = cyclo_structure[0]
        cyclo_chain = cyclo_model['A']
        ppi_model = ppi_structure[0]
        ppi_chain = ppi_model['{}'.format(protein)]
        
        cyclo_atoms, ppi_atoms = make_list_of_atoms_to_align(match, cyclo_model, ppi_model)
        match = superimpose(match_path, match, ppi_atoms, cyclo_atoms, cyclo_model, ppi_model, cyclo_structure, ppi_structure)

        #Find overlaping structures
        no_overlap_match, c1 = eliminate_overlap(match, pdb_name, c1, working_location, match_path, target)
        
        #Mutate residues of matches with no overlap
        if no_overlap_match != []:
            real_match_count, modelname, match_name = mutate_matched_side_chains(match_path, match, real_match_count, motif_size, ppi_chain, cyclo_chain)
            real_match_list.append(match_name)
        
        shutil.rmtree(match_dir)
    return(real_match_list)


def read_match_list(working_location, real_match_list, motif_size):
    C_list = []
    R_list = [] 
    for res in range(1, motif_size + 1):
        C_list.append("C" + str(res))
        R_list.append("R" + str(res))

    header=['Match', 'PDB', 'Dist-RMSD'] + C_list + R_list + ['Fit-RMSD']
    match_list = pd.read_csv(working_location + 'match_list.txt',  names=header, skiprows=6, sep='\s+')
    match_list['PDB Output Name']=real_match_list
    styles = [{'selector': 'th', 'props': [('text-align', 'left')]},
              {'selector': 'td', 'props': [('text-align', 'left')]}]
    match_list['Match'] = match_list['Match'].astype(str)
    display(match_list.style.hide_index().set_table_styles(styles))
    return(match_list)
