# protein_characterization.py
# Part of the cPEPmatch tool developed by Brianda L. Santini, supervised by Prof. Dr. Martin Zacharias
# Physics Department T38, Technical University of Munich, Garching, Germany.
#
# Module Description:
# This module focuses on the characterization of protein-target interactions within the cPEPmatch tool. 
# It provides essential functionalities for identifying interface residues between proteins and their 
# targets and for constructing and analyzing motifs within these interfaces. The module facilitates 
# the understanding of protein-protein or protein-ligand interactions, crucial for the subsequent steps 
# of matching cyclic peptides to these protein interfaces.


import sys
import os
import pickle
import warnings
import numpy as np
import pandas as pd
from scipy import spatial
from collections import defaultdict
from itertools import product, combinations, permutations
from Bio.PDB import PDBParser, parse_pdb_header
import Bio
import vmd
from vmd import molecule, atomsel

# Main function to create interfaces based on the given parameters
def create_interfaces(working_location, protein, target, pdb_name, interface_cutoff, residues=''): 
    os.chdir('{}' .format(working_location))
    parser = PDBParser()
    try:
        pdb_fname = ("{}.pdb" .format(pdb_name))
        molid  = molecule.load("pdb", pdb_fname) #from vmd library
        
        if residues:
            Iall = interfaces_for_specific_residues(molid, residues, protein)
        else:
            Iall = interfaces_interface_cutoff(molid, interface_cutoff, target, protein)
            
        int_fname = 'interface.pdb'
        print('{} -> {}' .format(pdb_fname, int_fname))
        molecule.write(molid, 'pdb', int_fname, selection=Iall)
        molecule.delete(molid)
        
    except RuntimeError as e:
        warn(e)
        pass
    return(Iall)

# Function to determine interface residues based on a distance cutoff
def interfaces_interface_cutoff(molid, cutoff, target, protein, minresidues=4):
    I = []
    I = atomsel("chain {} and same residue as exwithin {} of chain {}".format(protein, cutoff, target))
    if len(set(I.residue)) < minresidues: 
        I = []
    return I
# Function to determine interface residues based on specific residue IDs
def interfaces_for_specific_residues(molid, protein_resid, protein, minresidues=4):
    I = []
    I = atomsel("chain {} and resid {}".format(protein, protein_resid))
    if len(set(I.residue)) < minresidues: 
        I = []
    return I

# Function to convert a list of numbers into a string
def int_list(lyst, sep=' ', quote=None):
    """Return the list of numbers in 'lyst' as a string separated by 'sep'."""
    f = str if quote is None else lambda i: quote+str(i)+quote
    return sep.join(map(f, lyst))

# Function to load interface data from a PDB file
def load_interface(pdb_fname): 
    """Return a vmd molecule for the interface found in pdb_fname in a tuple of (vmd-mol, pdb-code, tag, list-of-chains)."""
    d,f = os.path.split(pdb_fname) # d = path/.../, f = x.pdb
    pdb,tag,chains,ext = f.split('.') # e.g. 4gko, pp, B-G, pdb
    chains = chains.split('-')
    mol = molecule.load("pdb", pdb_fname) 
    return (mol, pdb, tag, chains) 

# Function to calculate a distance matrix for consecutive motifs
def distance_matrix_consecutive(nres, R, fname):
    M = dict()
    for i,res in enumerate(R[:-nres+1]):
        resids = R[i:i+nres]
        if not np.all(np.diff(resids)==1): continue
        rlist = int_list(resids, quote='"')      
        resid_query = atomsel("name CA and resid {}" .format(rlist))
        clz = resid_query.structure
        r = np.array(resid_query.centerperresidue())
        
        # Ignore sequences that are either near the ends of a chain, or are near ends of a broken chain
        
        if len(r) != nres: continue
        assert len(r) == nres, (len(r), r, nres, res, rnbrs, resid_query, fname)
        
        dr = spatial.distance_matrix(r,r)
        m = dr[np.triu_indices(nres,k=2)]
        M[resids[0]] = dict(
            m=m,
            resid=resids,
        )
    return M
    
# Function to calculate a general distance matrix
def distance_matrix(nres, R, fname):
    M = dict()
    combination = combinations(R, nres)
    count = 1
    for i in combination:
        resids = i
        #if not np.all(np.diff(resids)==1): continue
        rlist = int_list(resids, quote='"')
        #if not rnbrs: continue
        resid_query = atomsel("name CA and resid {}" .format(rlist))
        r = np.array(resid_query.centerperresidue())
        # Ignore sequences that are either near the ends of a chain, or are near ends of a broken chain
        #if len(r) != nres: continue
        #assert len(r) == nres, (len(r), r, nres, res, rnbrs, resid_query, fname)
        dr = spatial.distance_matrix(r,r)
        m = dr[np.triu_indices(nres,k=1)]
        M[str(count)] = dict(m=m, resid=resids)
        count = count + 1
    return M

# Function to convert a list of numbers into a string
def motifs(molid, motif_size, inter_fname, consecutive):
    P = atomsel("all")
    P_res = list(sorted(set(P.resid)))
    if consecutive:
        MP = distance_matrix_consecutive(motif_size, P_res, inter_fname)
    else:
        MP = distance_matrix(motif_size, P_res, inter_fname)
    return (MP)

# Function to process the interface and return motifs
def process_interface(working_location, pdb_name, motif_size, motif_type):
    os.chdir('{}'.format(working_location, pdb_name))
    parser = PDBParser()
    inter_fname = "interface.pdb"
    molid  = molecule.load("pdb", inter_fname)
    mtfs = motifs(molid, motif_size, inter_fname, motif_type)
    return(mtfs)