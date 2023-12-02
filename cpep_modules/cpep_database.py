# cpep_database.py
# Part of the cPEPmatch tool developed by Brianda L. Santini, supervised by Prof. Dr. Martin Zacharias
# Physics Department T38, Technical University of Munich, Garching, Germany.
# 
# Module Description:
# This module is dedicated to the creation and characterization of the cyclic peptide (cPEP) database.
# It includes functions for reading and processing cyclic peptide data, calculating distance matrices,
# and characterizing peptide motifs. The module plays a critical role in the cPEPmatch tool by 
# providing the foundational data needed for matching and analysis of cyclic peptides in the context
# of protein mimicry and target binding.


import sys
import os
import operator
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

# Function to calculate distance matrix for cyclic peptides
def cyclo_distance_matrix(nres, R, fname):
    M = dict()
    combination = combinations(R, nres)
    count = 1
    for i in combination:
        resids = i
        #if not np.all(np.diff(resids)==1): continue
        rlist = int_list(resids, quote='"')
        #if not rnbrs: continue
        QA = '''(name CA and (altloc "A" or altloc "") and resid %(rlist)s)'''
        QA = ' '.join(str(QA).split())
        
        resid_query = atomsel(QA % dict(rlist=rlist))
        r = np.array(resid_query.centerperresidue())
        # Ignore sequences that are either near the ends of a chain, or are near ends of a broken chain
        if len(r) != nres: continue
        assert len(r) == nres, (len(r), r, nres, res, rnbrs, resid_query, fname)
        dr = spatial.distance_matrix(r,r)
        m = dr[np.triu_indices(nres,k=1)]
        M[str(count)] = dict(m=m, resid=resids)
        count = count + 1
    return M

# Function to convert a list of numbers to a string separated by a given separator
def int_list(lyst, sep=' ', quote=None):
    """Return the list of numbers in 'lyst' as a string separated by 'sep'."""
    f = str if quote is None else lambda i: quote+str(i)+quote
    return sep.join(map(f, lyst))


# Function to calculate distance matrix for consecutive cyclic peptides
def cyclo_distance_matrix_consecutive(nres, R, fname):
    M = dict()
    
    for i,res in enumerate(R[:-nres+1]):
        #print(res, res+nres, R[i:i+nres])
        resids = R[i:i+nres]
        if not np.all(np.diff(resids)==1): continue
        rlist = int_list(resids, quote='"')
        #if not rnbrs: continue
        QA = '''(name CA and (altloc "A" or altloc "") and resid %(rlist)s)'''
        QA = ' '.join(str(QA).split())
        
        resid_query = atomsel(QA % dict(rlist=rlist))
        r = np.array(resid_query.centerperresidue())

        # Ignore sequences that are either near the ends of a chain, or are near ends of a broken chain
        if len(r) != nres: continue
        assert len(r) == nres, (len(r), r, nres, res, rnbrs, resid_query, fname)
        
        dr = spatial.distance_matrix(r,r)
        m = dr[np.triu_indices(nres,k=2)]
        M[resids[0]] = dict(
            m=m,
            #r=r,
            resid=resids,
        )
    return M


# Function to characterize motifs in cyclic peptides
def cyclo_motifs(molid, cyclo, motif_size, consecutive):
    MC = dict()
    C = atomsel('protein and backbone and (altloc "A" or altloc "") ')
    C_res = list(sorted(set(C.resid)))
    if consecutive:
        MC = cyclo_distance_matrix_consecutive(motif_size, C_res, cyclo)
    else:
        MC = cyclo_distance_matrix(motif_size, C_res, cyclo)
    return MC

# Function to create a library database of cyclic peptides
def create_cyclolib_database(database_location, motif_size, consecutive,  
                             cyclization_type, exclude_non_standard,
                             csv_database_file='cyclo_pep.csv'):
    
    os.chdir(database_location)
    parser = PDBParser()
    cyclo_mtfs = []
    lib = []
    
    if consecutive:
        motif_type = 'consecutive'
    else:
        motif_type = 'not_consecutive'
        
    complete_database = (pd.read_csv(csv_database_file)).to_numpy()
    
    #Take only chosen peptides
    
    for peptide in complete_database:
        
        if exclude_non_standard:
            if cyclization_type:    
                if peptide[5] == 'None':
                    if peptide[2] == cyclization_type:
                        lib.append(peptide[1])
            else:
                if peptide[5] == 'None':
                    lib.append(peptide[1])
                    
        elif cyclization_type:
            if peptide[2] == cyclization_type:
                lib.append(peptide[1])
        else:
            lib.append(peptide[1])
            

    for standard_peptide in lib:
        motifs = dict()
        cyclo = ('{}.pdb' .format(standard_peptide))
        molid = molecule.load("pdb", cyclo)
        c_mtfs = cyclo_motifs(molid, cyclo, motif_size, consecutive)
        
        if c_mtfs is not None:
            motifs.update({standard_peptide: c_mtfs})
            cyclo_mtfs.append(motifs)
            
    
    database = open("database_{}-{}.pkl" .format(motif_size, motif_type), "wb")
    pickle.dump(cyclo_mtfs, database)
    database.close()
    database_name = "database_{}-{}.pkl" .format(motif_size, motif_type)
    return(database_name)