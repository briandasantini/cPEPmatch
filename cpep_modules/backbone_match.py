# backbone_match.py
# Part of the cPEPmatch tool developed by Brianda L. Santini, supervised by Prof. Dr. Martin Zacharias
# Physics Department T38, Technical University of Munich, Garching, Germany.
#
# Module Description:
# This module is designed to handle the backbone matching process in the cPEPmatch tool. 
# It focuses on matching the CA (carbon alpha) distance motifs of protein-protein interfaces (PPI) 
# with those in the cyclic peptide database. Key functionalities include reading and interpreting 
# motif data from the database, comparing these motifs with the PPI motifs, and calculating 
# the fit-RMSD (Root Mean Square Deviation) to identify potential matches.

import sys
import os
import pickle
import csv
import argparse
import numpy as np
import pandas as pd
from scipy import spatial
from collections import defaultdict
from itertools import product, combinations, permutations
from Bio.PDB import PDBParser, parse_pdb_header
import Bio
import vmd
from vmd import molecule, atomsel

# Read and return motifs from the cyclic peptide database.
def read_database_motifs(database_location, motif_size, motif_type):
    cyclo_mtfs = [ ]
    with open("{}database_{}-{}.pkl" .format(database_location, motif_size, motif_type),'rb') as database:
        cyclo_mtfs = pickle.load(database)
    return(cyclo_mtfs)
    
#Find matches between PPI motifs and cyclic peptide database motifs based on RMSD.
def find_match(database_location, mtfs, motif_size, motif_type, frmsd_threshold):
    if motif_type:
        motif_type = 'consecutive'
    else:
        motif_type = 'not_consecutive'
    cyclo_mtfs = read_database_motifs(database_location, motif_size, motif_type)
    database_length = len(cyclo_mtfs)
    
    match = dict()
    all_match = []
    matches = []
    d_ppi = []
    
    for pp in mtfs:
        m_count = len(mtfs[pp]["m"])
        d_ppi = []
        for j in range(0, m_count):
            d_ppi_single = mtfs[pp]["m"][j]
            d_ppi.append(d_ppi_single)
        for i in range(0,database_length):
            for c in cyclo_mtfs[i]:
                peptide = c
                for d in (cyclo_mtfs[i][c]):
                    d_cyclo = []
                    for k in range(0, m_count):
                        d_cyclo_single = cyclo_mtfs[i][c][d]["m"][k]
                        d_cyclo.append(d_cyclo_single)
                    sum = 0
                    for r in range(0, m_count):
                        dd = (d_ppi[r]-d_cyclo[r])**2
                        sum = sum + dd
                    rmsd = np.sqrt(sum)
                    if rmsd < frmsd_threshold:
                        rmsd = rmsd.round(decimals=4)
                        match = ( peptide , rmsd,  cyclo_mtfs[i][c][d]["resid"], mtfs[pp]["resid"] )
                        all_match.append(match) 
    return(all_match)