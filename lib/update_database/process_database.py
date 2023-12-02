#!/usr/bin/env python
# coding: utf-8

# # cPEP Database Construction


import os
import re
import glob
import numpy as np
import pandas as pd
import subprocess
import sys
from numpy.linalg import norm
import pytraj as pt
import mdtraj as md 
from itertools import groupby
from scipy.spatial.distance import squareform, pdist  
import collections
import statistics
#from statistics import multimode
from pathlib import Path


# ## Functions

# In[33]:


def abbreviate_seq(ints, presort=True):
    """
    Given `ints`, an array of int, intervalize into lists of 2-tuples.

    Example:
        Input: [1,2,3,5,6,8]
        Output: [[1, 3], [5, 6], [8, 8]]
    """
    if presort:
        ints = sorted(ints)

    ret = []
    temp = [ints[0], ints[0]]
    if len(ints) == 1:
        return [temp]

    for i in ints[1:]:
        if temp[1] + 1 != i:
            # new interval
            ret.append(temp)
            temp = [i, i]
        else:
            # update interval
            temp[1] = i
    ret.append(temp)

    return ret


def expand_tupseq(tupseq):
    """
    Expand 2-tuples into a sequence of numbers.
    Used in conjunction with `abbreviate_seq`.

    Example: [(1, 3), (5, 5), (7, 8)] -> [1, 2, 3, 5, 7, 8]
    """
    seq = []
    for tup in tupseq:
        seq.extend(range(tup[0], tup[1]+1))
    return seq


def seq_to_string(seq, offset=0, presort=True):
    """
    Abbreviate a sequence of numbers into a comma-separated list of numbers
    and hyphenized number pairs.

    Example: [1,2,3,4,6,8,9] -> '1-4,6,8-9'
    """
    if offset != 0:
        seq = [el+offset for el in seq]

    tuples = abbreviate_seq(seq, presort=presort)
    ret = []
    for tup in tuples:
        if tup[0] == tup[1]:
            ret.append(str(tup[0]))
        else:
            ret.append("{}-{}".format(tup[0], tup[1]))
    return ",".join(ret)


def create_chain_dict(top: pt.Topology, only_resid=False, key="chain"):
    """
    Create a dictionary of residue objects, grouped by their chain index

    Parameters
    ----------
    only_resid : bool
        If true, residues are only stored as residue indices (counting from 0),
        else residue objects (pytraj.core.topology_objects.Residue) are stored
    key : str ["chain", "residue"]
        If "chain", create a chain:list_of_residues dictionary.
        If "residue", create a resid:chain dictionary (automatically assumes only_resid=True)
    """
    if key not in ["chain", "residue"]:
        raise ValueError("`key` has to be either 'chain' or 'residue'")
    f_chaindict = key == "chain"

    chain_dict = {}
    for res in top.residues:
        atom = top.atom(res.first)
        chain = atom.chain

        if f_chaindict:
            reslist = chain_dict.setdefault(chain, [])
            reslist.append(res.index) if only_resid else reslist.append(res)
        else:
            chain_dict[res.index] = chain
    return chain_dict


# only select files/directories (including pathway) that contain a certain string in title

def find_matching_files(patterns, file_dir):
   matches = []
   for pattern in patterns:
      search_path = os.path.join(file_dir, '*{}*'.format(pattern))
      for match in glob.iglob(search_path):
         matches.append(match)
   return matches 


def remove_water_and_H(root, fn):
    # original pdb file
    absfn = os.path.join(root, fn)
    
    # pdb4amber output file 
    newfn = os.path.join(root, "dry_nohyd_"+ fn)

        # Check if the newfn file already exists
    if not os.path.exists(newfn):
        cmd = "pdb4amber -i " + absfn + " -o " + newfn + " --dry --nohyd"

        # call pdb4amber command to remove water and hydrogen atoms
        subprocess.run(cmd, shell=True)
    else:
        print(f"File {newfn} already exists. Skipping pdb4amber processing.")

    return absfn, newfn


def load_topology_of_pdb(newfn):

    # load pdb file with mdtraj 
    pdb_md = md.load(newfn)

    # get topology from pdb file
    initial_topology = pdb_md.topology
 
    return pdb_md, initial_topology 



def read_connects():

    # output non-peptide covalent bonds using CONECTS 

    non_pep_cov=[]
    disulfide_bonds=0
    print(fn[-8:])
    cleanf = os.path.join(root,"clean_"+fn[-8:])
    # fout=open(cleanf, "w")
    with open(newfn, "r") as orig:
        for line in orig.readlines():
            # return number for atoms that engage in non-peptide covalent bonds
            if "CONECT" in line:
                CON=line.split()[1:]
    #             print([int(x) for x in CON])
                non_pep_cov.append([int(x) for x in CON])
                #print(non_pep_cov)
                #print(len(non_pep_cov))
                #print(non_pep_cov[0][0])

    # print atoms between which a non-peptide covalent bond exists (not hydrogens)
    # CON_no is a set of lists each containing the numbers of the CONECT atoms
    # len(CON_no) is the number of atoms in each CONECT 
    # CON_atom is the atom name that corresponds to the number of the atom
    # non_pep_bond_atoms is the name of all the atoms in one CONECT (list emptied after each CONECT)
    # print(len(atom_name))
    for CON_no in non_pep_cov:
    #     print("new CON")
        non_pep_bond_atoms=[]
        non_pep_bond_residue=[]
        for i in range(len(CON_no)):
    #         print(CON_no[i])
            # get atom name and residue ID from serial atom number
            CON_atom_df=table.loc[table['serial']==CON_no[i], 'name']
            CON_residSeq_df=table.loc[table['serial']==CON_no[i], 'resSeq']
            non_pep_bond_residue.append(np.array(CON_residSeq_df))
            CON_atom=np.array(CON_atom_df)
            non_pep_bond_atoms.append(CON_atom[0])
        print(non_pep_bond_atoms)
        non_pep_bond="-".join(non_pep_bond_atoms)

        # exclude non-pep bonds within same residue 
        # count how many times the first residue turns up in all the connected atoms - if it turns up as much as len(CON_no), then all connected atoms are from the same residue

        if len(CON_no) == non_pep_bond_residue.count(non_pep_bond_residue[0]): 
            continue
                
        # exclude peptide bonds 
        
        elif set(non_pep_bond_atoms).issubset(['C', 'CA', 'N', 'O'])==True:
            continue 
        
        # determine which kind of non-pep bond it is 

        else: 
            
            # print non-peptide bond
            print(non_pep_bond + " bond")

            # print disulfide bonds 
            if "SG-SG" == non_pep_bond:
                disulfide_bonds+=1
                print("disulfide bridge = %i"%disulfide_bonds)
            
            elif "SG-SG" not in non_pep_bond:
                
                # print side-chain to backbone bond 
                # if the connected atoms contain at least one backbone atom but are not fully comprised of them
                
                if 1<=len([el in non_pep_bond_atoms for el in ['C', 'CA', 'N', 'O']])<len(CON_no): 
                    print("side-chain to backbone")
                    
                # print side-chain to side chain bond 
                # if none of the non_pep_bond_atoms are in ['C', 'CA', 'N', 'O']
            
                elif any(ele in non_pep_bond_atoms for ele in ['C', 'CA', 'N', 'O'])==False: 
                    print("side-chain to side-chain")

                    # if one residue == "not amino acid" 
                    # print staple bond 

                    # else:
                    # print side-chain to side-chain bond 
                

# fout.close()


# make initial table and remake table after changing topology 

def rewrite_table(initial_topology):

    # display topology in table 

    initial_table, bonds = initial_topology.to_dataframe()

    return initial_table


# groupby chainID and then by resSeq

def groupby_chain_res(initial_table, residSeq):

    peptide_chains = {}

    for chain_id, chain_table in initial_table.groupby("chainID"):
        peptide_chains[chain_id] = chain_table

    peptide_residues = {}

    for chain_id, chain_table in peptide_chains.items():
        chain_residuedict = {}

        dfgroupby = chain_table.groupby(residSeq)

        for resSeq, res_table in dfgroupby:
            chain_residuedict[resSeq] = res_table
        peptide_residues[chain_id] = chain_residuedict
    
    return peptide_chains, peptide_residues


# if topology is empty, append pdb name to delete_pdb 
def remove_pdb(batch):  

    # select all files that contain pdb code of input file 

    non_cyclic_peptide_files = find_matching_files([fn[3:7]], batch + "/renum/")
    # print(non_cyclic_peptide_files)    
    
    # move files that do no contain cyclic peptides into the trash 

    for non_cp_file_path in non_cyclic_peptide_files:
        
        # split path to file so that you can move file into new directory 
        split_path = non_cp_file_path.split("/")

        # move non-cylic peptide file into trash directory
        os.rename(non_cp_file_path , batch + "/trash/" + split_path[-1])

    # append empty pdb file to delete_pdb 
    delete_pdb.append(pdb[PDB_count])

    return delete_pdb
    

# delete atoms of a peptide/protein that is longer or shorter than the required length of a cyclic peptide or is a duplicate

def initial_delete_peptide(peptide_chains, peptide_residues, initial_topology):
 
    # delete duplicates

    # residue sequences of all chains 
    all_resSeq=[]

    # index atoms that need to be deleted
    delete_atoms=[]

    for m,unsplit_chain in peptide_chains.items():
        # print(unsplit_chain[['serial','resName']])
        resSequence = unsplit_chain['resName'].values.tolist()
        
        # if chain is already present 
        if resSequence in all_resSeq: 

            # get atom indexes of duplicate chains
            delete_atoms.append(unsplit_chain.index[unsplit_chain["chainID"]==m].tolist())

        # if chain is not present 
        else: 
            all_resSeq.append(resSequence)

    # delete cyclic peptides with length diverging from standard cp length 

    for k,chain in peptide_residues.items():

        # iterate through residues
        for r, residue in chain.items():

            # delete chains that contain less than 4 and more than 30 residues
            if not 4 < len(chain) < 30:

                # append atom index in each residue to delete_atoms list 
                delete_atoms.append(residue.index[residue["chainID"]==k].tolist())

    # combine all atoms to delete (from duplicate and unwanted length) into one list
    concat_delete_atoms=[j for i in delete_atoms for j in i]

    # delete atom indexes that are in the list twice 
    concat_delete_atoms= list(set(concat_delete_atoms))

    # sort list from highest to lowest number (delete from back)
    concat_delete_atoms.sort(reverse=True)

    # delete atoms with index from delete_atoms list (mdtraj method)

    for atom in concat_delete_atoms:
        initial_topology.delete_atom_by_index(atom)

    # create new topology after deleting atoms
    topology = initial_topology

    # if topology does not contain any atoms, delete pdb name from list and move all files containing pdb to a trash directory  
    if topology.n_atoms == 0:
        delete_pdb = remove_pdb(batch)  

    return topology, concat_delete_atoms


def adjust_coordinates(input_coordindates, delete_atoms_list): 

    # delete coordinates of atoms that were deleted and are thus no longer in topology 
    coordinates_out=np.delete(input_coordindates, delete_atoms_list, axis=0)

    # rewrite table 
    table = rewrite_table(initial_topology)

    # adjust table
    x, y, z = map(list, zip(*coordinates_out))
    table['x'], table['y'], table['z'] = [x, y, z]

    return table, coordinates_out 


# groupby only residue 

def groupby_res(table):

    peptide_resids = {}

    for res_id, res_table in table.groupby("resSeq"):
        peptide_resids[res_id] = res_table
    
    return peptide_resids


# In[44]:


# determine if atoms are bonded based on distances between atoms 

def calculate_distance_matrix(table): 

    # Van der Waals Radii of atoms found in amino acids
    vdw_rad = {
        "N": .155,
        "C": .170,
        "O": .152,
        "S": .180,
        "Li": .182,
        "Be": .153,
        "B": .192,
        "F": .135,
        "Na": .227,
        "Mg": .173,
        "Al": .184,
        "Si": .210,
        "P": .180,
        "S": .180,
        "Cl": .175,
        "K": .275, 
        "Ca": .231,
        "Sc": .211, 
        "Ti": .187,
        "V": .179,
        "Cr": .189, 
        "Mn": .197, 
        "Fe": .194,
        "Co": .192, 
        "Ni": .163,
        "Cu": .140, 
        "Zn": .139, 
        "Ga": .187, 
        "Ge": .211, 
        "As": .185,
        "Se": .190, 
        "Br": .183, 
        "Rb": .303,
        "Sr": .249,
        "Rh": .195,
        "Pd": .202, 
        "Ag": .172,
        "Cd": .158,
        "In": .193,
        "Sn": .217, 
        "Sb": .206,
        "Te": .206, 
        "I": .198, 
        "Cs": .343,
        "Ba": .268, 
        "Ir": .202,
        "Pt": .209,
        "Au": .166,
        "Hg": .209, 
        "Tl": .196,  
        "Pb": .202, 
        "Bi": .207, 
        "Po": .197,
        "At": .202
    }

    # determine if atoms are close enough to be bonded 

    # append a column with the Van der Waals radii of the atoms 
    table['VDW_rad'] = table.apply(lambda x: vdw_rad[x.element], axis=1)

    # calculate the distance matrix (distances between all atoms)
    dist_mat=pd.DataFrame(squareform(pdist(table.iloc[:,7:10], 'euclid')))

    # calculate the radii matrix (cutoff in order to determine if atoms are bonded)
    radii_sum_mat = pd.DataFrame([[(x + y)*0.6 for x in table['VDW_rad']] for y in table['VDW_rad']])

    # bond exists if the calculated distance is lower than the radii matrix values
    bond_mat_true = dist_mat < radii_sum_mat

    #combine distance matrix with the 'True' or 'False' bond table
    bond_mat = dist_mat.combine(bond_mat_true, np.multiply)

    # only take lower triangle of bond matrix (don't repeat the same bond twice)
    bond_mat.values[np.triu_indices_from(bond_mat, k=1)] = np.nan 
    # fill diagonals with zero 
    # distances between atom and self are zero
    np.fill_diagonal(radii_sum_mat.values, 0)
    # replace 0 with nan
    bond_mat = bond_mat.replace(0, np.nan)
    # drop all nan 
    bond_mat = bond_mat.unstack().dropna()

    # combine bonded atom indexes 
    pair_indices = np.array([list(pairkey) for pairkey in bond_mat.keys()])
    
    # if there are not any bonds present
    if len(pair_indices) == 0: 
        delete_pdb = remove_pdb(batch)

    return pair_indices 


def create_bond_table(pair_indices, table):
    
    # create table of bonds between all atoms within each residue 
    bond_df = pd.DataFrame(pair_indices, columns=['at1', 'at2'])

    # add atom name columns 
    bond_df['at1_name'] = table.name[pair_indices[:, 0]].to_numpy()
    bond_df['at2_name'] = table.name[pair_indices[:, 1]].to_numpy()

    # add atom element columns 
    bond_df['at1_element'] = table.element[pair_indices[:, 0]].to_numpy()
    bond_df['at2_element'] = table.element[pair_indices[:, 1]].to_numpy()

    # add residue ID columns    
    bond_df['resSeq1'] = table.resSeq[pair_indices[:, 0]].to_numpy()
    bond_df['resSeq2'] = table.resSeq[pair_indices[:, 1]].to_numpy()

    # add chain ID column so that you can delete atoms by chain later on (if necessary)
    bond_df['chainID_1'] = table.chainID[pair_indices[:, 0]].to_numpy()
    bond_df['chainID_2'] = table.chainID[pair_indices[:, 1]].to_numpy()

    # if bond between atoms in different chains then delete bonds 
    bond_df = bond_df[bond_df["chainID_1"]==bond_df["chainID_2"]]

    # remove chainID_2 from dataframe
    bond_df.drop("chainID_2", axis=1, inplace=True)
    
    # rename chain column 
    bond_df.rename(columns = {'chainID_1':'chainID'}, inplace = True)

    # select bonds within residue 
    bonds_within_resid = bond_df.loc[(bond_df['resSeq1']== bond_df['resSeq2'])]

    # bonds between residues 
    bonds_betwix_resid= bond_df.loc[~(bond_df['resSeq1']== bond_df['resSeq2'])]

    #  duplicate bonds between residues so that each residue contains all bonds
    duplicate_bonds = bonds_betwix_resid.copy()

    # invert at1 and at2, at1name and at2name, at2_element, at1_element, resSeq1 and resSeq2
    new_column_titles = ["at2", "at1", "at2_name", "at1_name", "at2_element", "at1_element", "resSeq2", "resSeq1", "chainID"]
    reordered_columns=duplicate_bonds.reindex(columns=new_column_titles)
    reordered_columns.columns=["at1", "at2", "at1_name", "at2_name", "at1_element", "at2_element", "resSeq1", "resSeq2", "chainID"]

    # combine dataframes within resid, betwix residues and reordered betwix residue   
    all_bonds = pd.concat([bonds_within_resid, bonds_betwix_resid, reordered_columns], ignore_index=True, sort=False)

    return all_bonds


def check_peptide_is_cyclic(peptide_residues_bonds, table):

    # check if peptide is cyclic 

    # list of non_cyclic peptide atoms if present 
    non_cyclic_peptide_atoms =[]

    # concatenate list of bonds between residues for each chain 
    concat_bonds_betwix_res = []
    
    # chain ID of cyclic peptide 
    chain_ID = None

    # chain ID of non-cyclic peptide 
    chain_ID_delete = []

    # number of cyclic peptide chains (only take one cyclic peptide)
    cyclic = []

    # iterate through chains
    for k, chain_bonds in peptide_residues_bonds.items(): 

        # list of bonds between each residue and other residues 
        bonds_betwix_res_per_chain = []
        
        # iterate through residues 
        for k2, bonds_per_residue in chain_bonds.items(): 

                # count number of residues with which residue is bonded 
                diff_residues = 0 

                # check how many residues each residue is bonded with
                
                for index, row in bonds_per_residue.iterrows():
                    if row['resSeq1']!=row['resSeq2']:
                        diff_residues+=1 
                
                bonds_betwix_res_per_chain.append(diff_residues)

        # check if all residues bind to two other residues
        if all(elem == 2 for elem in bonds_betwix_res_per_chain) and 1 not in cyclic:
            cyclic.append(1)
            
            # if cyclic, save chainID of peptide 
            chain_ID = k
        
        # if one residue has 3 bonds  

        elif 3 in bonds_betwix_res_per_chain and 1 not in cyclic: 
            cyclic.append(1)

            chain_ID = k

        # if not cyclic 

        else:
            cyclic.append(0)

            # remove sublist
            bonds_betwix_res_per_chain = []
            
            # get atom indexes of non-cyclic peptide chain atoms in respective chain from table   
            cyclic_peptide_atom_indexes = list(table[table["chainID"]==k].index.values)

            # append atom indexes from non-cyclic peptide chain to non_cyclic_peptide atoms 
            non_cyclic_peptide_atoms.extend(cyclic_peptide_atom_indexes)
            
            # append chain ID of non-cyclic peptide to chainID_delete
            chain_ID_delete.append(k)

        # append sublist to bonds_betwix_res list 
        concat_bonds_betwix_res.append(bonds_betwix_res_per_chain)

    # list of bonds between each residue and other residues only in correct chain  
    bonds_betwix_res = [j for i in concat_bonds_betwix_res for j in i]

    return chain_ID, non_cyclic_peptide_atoms, chain_ID_delete, bonds_betwix_res
            

def write_cp_final_table(non_cyclic_peptide_atoms, topology):

    # if list of non cyclic peptide atoms to be deleted is not empty delete these atoms
    
    if len(non_cyclic_peptide_atoms) != 0: 
        
        # delete chain from peptide_residues_bonds 
        for k in chain_ID_delete:
            peptide_residues_bonds.pop(k)

        # reverse order of atom indexes of non-cyclic peptides (all chains) to be deleted 
        non_cyclic_peptide_atoms.sort(reverse=True)

        # delete atoms from topology with index from non_cyclic_peptide_atoms list (mdtraj method)

        for atom in non_cyclic_peptide_atoms:
            topology.delete_atom_by_index(atom)
    
        # create new topology after deleting non-cyclic peptide atoms 
        cp_topology = topology 

        # if topology does not contain any atoms, delete pdb name from list and move all files containing pdb to a trash directory  
        if cp_topology.n_atoms == 0:
            delete_pdb = remove_pdb(batch)
            
            # final table and coordinates are empty
            final_table = None
            final_coordinates = None

        # if topology is not empty, rewrite table 

        else:

            # delete coordinates of atoms that were deleted and are thus no longer in topology and rewrite table 
            final_table, final_coordinates = adjust_coordinates(coordinates, non_cyclic_peptide_atoms)
    
    else: 
        # if table, coordinates and topology don't change, then keep the same
        final_table = table
        final_coordinates = coordinates
        cp_topology = topology
    
    return peptide_residues_bonds, final_table, final_coordinates, cp_topology

# note that you can no longer compare atom indexes in peptide_residues_bonds with atom serials in final_table (no longer the same)


# iterate through peptide residues in cyclic peptide chain and determine features of cyclic peptide
def determine_cp_features(peptide_residues_bonds, bonds_betwix_res_num):

    # number of amino acids (out of all residues in one cyclic peptide chain)
    amino_acid_num = 0

    # number of residues in one cyclic peptide chain
    residue_num = 0

    # list of standard backbone bonds (also any of these bonds could be reversed)
    std_bb_bonds = ["N-CA", "CA-C", "C-O", "CA-N", "C-CA", "O-C"]
    # list of standard peptide bond 
    std_pep_bond = ["N-C", "C_N"]
    # list of standard amino acids 
    std_amino_acids = ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'CYS', 'SEC', 'GLY', 'PRO', 'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP']
    bonds_betwix_resid_name=[]
    cp_type=[]
    Caps=[]
    non_standard_aa=[]
    other_cp_features=[]

    # iterate through all bonds for each residue
    for bpr, final_bonds_per_resid in peptide_residues_bonds[chain_ID].items(): 

        # number of residues
        residue_num += 1

        # number of atoms N binds with within residue 
        N_bonds=0
        # number of atoms CA binds with within residue 
        CA_bonds=0
        # create list of C atom indexes that bond with N or CA 
        C_index=[]

        # list all bonds within each residue (name)
        bonds_within_resid_name=[]
        # list all bonds within each residue (element)
        bonds_within_resid_element=[]

        # check if each residue contains correctly bonded amino acid backbone
        
        # iterate through each bond for each residue        
        for index, row in final_bonds_per_resid.iterrows():

            # check bonds within residue
            if row['resSeq1']==row['resSeq2']:
                
                # get atom names from bonded atom indexes 
                bonded_at1_name = row['at1_name']
                bonded_at2_name = row['at2_name']
                bonded_at1_element = row['at1_element']
                bonded_at2_element = row['at2_element']

                # append bonded atom name pairs to list  
                bond_within_resid_name = bonded_at1_name + "-" + bonded_at2_name
                bonds_within_resid_name.append(bond_within_resid_name)

                # get atom elements from bonded atom indexes
                bond_within_resid_element = bonded_at1_element + "-" + bonded_at2_element
                bonds_within_resid_element.append(bond_within_resid_element)

                # if N in bond, enumerate 
                if bonded_at1_name=="N" or bonded_at2_name=="N":

                    N_bonds+=1

                    # index C atom that is bound to N atom 
                    if "C" in bond_within_resid_element:

                        # replace comma so that you can easily index atoms
                        bond_elements = bond_within_resid_element.replace("-","")
                        
                        # find at1/at2 of atom of C element atoms that bind to N
                        C_index.append(row['at'+str(bond_elements.index("C")+1)])

                elif bonded_at1_name=="CA" or bonded_at2_name=="CA":

                    CA_bonds+=1

                    # index C atom that is bound to CA atom 
                    if "C" in bond_within_resid_element:

                        # replace comma so that you can easily index atoms
                        bond_elements = bond_within_resid_element.replace("-","")

                        # find at1/at2 of atom of C element atoms that bind to CA
                        C_index.append(row['at'+str(bond_elements.index("C")+1)])

                # give number of rows that contain C_index
                

            # check if N is bonded to two atoms and if CA is bonded to 3 atoms (not including N-CA) and the C element is not bonded to any 

            # check bonds between residues (where resSeq1 is smaller than resSeq2)
            elif row['resSeq1']<row['resSeq2']:  
            
                # get atom names from bonded atom indexes
                bond_betwix_resid_name = row['at1_name'] + "-" + row['at2_name']
                bonds_betwix_resid_name.append(bond_betwix_resid_name)
            
            # ignore duplicate bonds
            else:
                pass

        # print(bonds_within_resid_name)
        # print(bonds_within_resid_element)
        
        # if residue contains correctly bonded backbone atoms to be a standard amino acid
        if len([i for i in bonds_within_resid_name if i in std_bb_bonds]) == 3:
            amino_acid_num+=1
            # print("standard amino acid")
            # print(row['resSeq1'])

            # select resname of non-standard amino acids 
            amino_acid_name = final_table.loc[final_table["resSeq"]==bpr, "resName"].iloc[0]
            if amino_acid_name not in std_amino_acids:
                non_standard_aa.append(amino_acid_name)
                
                # iterate through C element atoms that are bonded to N or CA name atoms  
                for ind in C_index:

                    # get number of C element atom bonds
                    C_bonds = len(final_bonds_per_resid[(final_bonds_per_resid["at1"]==ind) | (final_bonds_per_resid["at2"]==ind)])
                
                    # if N name atom bonds to two other atoms (CA and C) within resid and this C atom binds to only the N  
                    if N_bonds == 2 and C_bonds == 1:
                        
                        # append methyl group attached to N to other_features list   
                        other_cp_features.append(amino_acid_name + "(methylated at N)")

                    # if CA name atom bonds to 3 other atoms (R group, Carboxyl group and methyl group - (N not included again)) within resid and the Methyl group C atom binds to only CA (would not include methylated glycine) 
                    elif CA_bonds == 3 and C_bonds ==1:

                        # append methyl group attached to CA to other_features list
                        other_cp_features.append(amino_acid_name + "(methylated at CA)")    
                
                # print("C_bonds: " + str(C_bonds))
        
        # if residue is not std. amino acid
        else:
            # if residue is only bonded to one other residue, it is a cap
            if bonds_betwix_res_num[residue_num-1] == 1:
                Caps.append(final_table.loc[final_table["resSeq"]==bpr, "resName"].iloc[0])

            # if residue bonded to two other residues, it is a staple
            elif bonds_betwix_res_num[residue_num-1] == 2:

                # if cyclic peptide stapled, make this cyclic peptide type
                cp_type.append(final_table.loc[final_table["resSeq"]==bpr, "resName"].iloc[0] + " STAPLE")

                # if staple, change chainID to 0 in final_table

                # if Staple name in line, make new chain with value = 0
                final_table.loc[final_table["resSeq"]==bpr, 'chainID'] = 0
                
    # split strings in bonds_betwix_resid_name 

    # if disulfide bridge in peptide, append to cyclic peptide type
    if "SG-SG" in bonds_betwix_resid_name:
        cp_type.append(str(bonds_betwix_resid_name.count("SG-SG")) + "x disulfide bridge")

    # if there are as many standard peptide bonds as there are residues then the cyclic peptide is head to tail
    if (bonds_betwix_resid_name.count("C-N") + bonds_betwix_resid_name.count("N-C")) == residue_num:
        cp_type.append("head to tail")

    else:
        # iterate through each bond between residues 
        for name_bonds in bonds_betwix_resid_name:
            
            # if not peptide bond
            if name_bonds != "C-N" and name_bonds != "N-C":
                # print(name_bonds)
                
                # split name bonds 
                split_name_bonds = name_bonds.split("-")
                atom1=split_name_bonds[0]
                atom2=split_name_bonds[1]
                
                # if STAPLE in cyclic peptide, do not mention how its bonds 
                if len(cp_type) != 0: 
                    pass

                # if one bond atom is C or N it is side-chain to backbone  
                elif atom1=="C" or atom1=="N" or atom2=="C" or atom2=="N":
                    cp_type.append("side-chain to backbone")
                
                # if bond atoms are both not C or N 
                else: 
                    cp_type.append("side-chain to side chain")

    # append length of cyclic peptide to a_a_length list 
    a_a_length.append(amino_acid_num)            

    # append non-standard amino acids to non_standard_aa list
    if len(non_standard_aa)== 0: 
        non_standard_aa.append("None")

    # append None to Caps if list empty
    if len(Caps)== 0: 
        Caps.append("None")

    # append number of cp features to other_cp_features
    if len(other_cp_features)== 0: 
        other_cp_features.append("None")

    return final_table, cp_type, a_a_length, non_standard_aa, Caps, other_cp_features


# renumber serial, resSeq in final_table

def renumber_final_table(final_table):

    # sort final table by chain_ID and if duplicate value present then by serial number
    cp_table1 = final_table.sort_values(by = ['chainID', 'serial'])
    
    # reset index 
    cp_table2 = cp_table1.reset_index()

    #make index into new serial column 
    cp_table2['serial'] = cp_table2.index

    # if the chainID in row is == chain_ID then change chainID to 1
    cp_table2.loc[cp_table2["chainID"]==chain_ID, 'chainID'] = 1

    # if resName changes, add 1 to the resSeq 
    cp_table2['resSeq'] = (cp_table2["resSeq"] != cp_table2["resSeq"].shift()).cumsum()

    # remove index column 
    cp_final_table = cp_table2.drop(['index'], axis=1)

    return cp_final_table


def save_pdb(cp_final_table, final_coordinates):

    # convert final_table into final_topology 
    final_topology = md.Topology.from_dataframe(cp_final_table)

    # convert topology and coordinates into pdb file
    final_pdb = md.Trajectory(final_coordinates, final_topology)
    final_pdb.save_pdb(batch + "/libraryready_coded/" + fn[0:4] + '.pdb')

    return final_pdb


def save_pdb_try(cp_final_table, final_coordinates, batch, fn):
    # Construct the file path
    file_path = os.path.join(batch, "libraryready_coded", fn[0:4] + '.pdb')

    # Convert final_table into final_topology 
    final_topology = md.Topology.from_dataframe(cp_final_table)

    # Convert topology and coordinates into a PDB file
    final_pdb = md.Trajectory(final_coordinates, final_topology)

    # Check if the file already exists
    if os.path.exists(file_path):
        print(f"File {file_path} already exists. Skipping save.")
        return final_pdb  # Return the final_pdb object even if not saving

    try:
        # Attempt to save the PDB file
        final_pdb.save_pdb(file_path)
        print(f"File {file_path} has been successfully saved.")
    except Exception as e:
        print(f"Failed to save the file {file_path}. Error: {e}")

    return final_pdb



# find secondary structure of cyclic peptide from md trajectory

def find_secondary_structure(final_pdb):
    sec_structure = md.compute_dssp(final_pdb)

    # replace secondary structure symbols with actual secondary structure name 
    replacements = {
        'H':'alpha helix',
        'G':'alpha helix',
        'I':'alpha helix',
        'B':'hairpin',
        'E':'hairpin', 
        'C':'random coil',
        'T':'random coil',
        'S':'random coil'
        }
    replacer = replacements.get

    structure_prevalence = [replacer(n, n) for n in sec_structure[0]]

    # get frequency of each secondary structure 
    secondary_structure.append(dict(zip(*np.unique(structure_prevalence, return_counts=True))))
    
    return secondary_structure


# In[53]:


# combine all strings in sub-list of cp_type, non_standard_aa, Caps and other_cp_features list 
def list_comprehension(feature_list):
    new_feature_list= [', '.join(sub_list) for sub_list in feature_list]
    return new_feature_list


def format_structure_dict(structure_dict):
    if structure_dict:
        # Format the dictionary into 'key: value' pairs
        formatted_items = [f"{key}: {value}" for key, value in structure_dict.items()]
        return ', '.join(formatted_items)
    return None


# ## The Actual Code

# working directory 
cp_directory = r'.'

# select all batches from working directory by calling function
all_batches = find_matching_files(["batch2023"], cp_directory)
all_batches

# create list for pdb names, amino acid length, secondary structure, non-standard amino acids, CAPs and other features for cyclic peptide table 

pdb=[]
a_a_length=[]
secondary_structure=[]
cp_type_new=[]
non_standard_aa_new=[]
Caps_new=[]
other_cp_features_new=[]

# create a list of pdbs that do not contain cyclic peptides and remove these from pdb list (won't be in cp_dataframe)
delete_pdb=[]

# iterate through batches 
for batch in all_batches:

    print(batch)
    
    PDB_count = 0
    # directory through which we iterate
    root = batch+""

    # iterate through pdb files
    for PDB in os.listdir(root):
        
        # only use files that start with "pdb" as input files
        if PDB[4:14] == "_renum.pdb":
            pdb.append(PDB[0:4])
        
            fn = PDB

            # run everything 

            # remove hydrogen atoms and water from pdb 
            absfn, newfn = remove_water_and_H(root, fn)
            
            # if the pdb4amber command outputs a file 
            if Path(newfn).is_file():

                # get topology from pdb file
                pdb_md, initial_topology = load_topology_of_pdb(newfn)
            
                # write topology into table
                initial_table = rewrite_table(initial_topology)
            
                # call function groupby chains and residues
                peptide_chains, peptide_residues = groupby_chain_res(initial_table, "resSeq") 

                # call function to delete peptides that are not correct length  
                topology, concat_delete_atoms = initial_delete_peptide(peptide_chains, peptide_residues, initial_topology)

                # if peptide with correct length
                if topology.n_atoms != 0: 
                
                    # adjust coordinates and rewrite table 
                    table, coordinates = adjust_coordinates(pdb_md.xyz[0], concat_delete_atoms)

                    # call function to group updated peptide by only residue (not chain)
                    peptide_resids = groupby_res(table)

                    # call distance matrix function 
                    pair_indices = calculate_distance_matrix(table)
                    
                    # if bonds exist
                    if len(pair_indices) != 0: 

                        # call function that creates table of bonded atoms in residues 
                        all_bonds = create_bond_table(pair_indices, table)

                        # call function groupby chains and residues
                        peptide_chains_bonds, peptide_residues_bonds = groupby_chain_res(all_bonds, "resSeq1")

                        # call function that finds non-cyclic peptides
                        chain_ID, non_cyclic_peptide_atoms, chain_ID_delete, bonds_betwix_res_num = check_peptide_is_cyclic(peptide_residues_bonds, table)
                        
                        # delete found non-cyclic peptides and rewrite table
                        peptide_residues_bonds, final_table, final_coordinates, cp_topology = write_cp_final_table(non_cyclic_peptide_atoms, topology)

                        # after second deletion 
                        if cp_topology.n_atoms != 0: 
                            
                            # iterate through peptide residues in cyclic peptide chain and determine features of cyclic peptide
                            final_table, cp_type, a_a_length, non_standard_aa, Caps, other_cp_features = determine_cp_features(peptide_residues_bonds, bonds_betwix_res_num)

                            # renumber serial, resSeq in final_table
                            cp_final_table = renumber_final_table(final_table)
                            
                            # convert table into trajectory and save as pdb file 
                            
                            final_pdb = save_pdb_try(cp_final_table, final_coordinates, batch, fn)

                            # call function that finds secondary structure of cyclic peptide 
                            secondary_structure = find_secondary_structure(final_pdb)
                            
                            # append lists for dataframe
                            cp_type_new.append(', '.join(cp_type))
                            non_standard_aa_new.append(', '.join(non_standard_aa))
                            Caps_new.append(', '.join(Caps))
                            other_cp_features_new.append(', '.join(other_cp_features))

            else:
                delete_pdb = remove_pdb(batch)

            PDB_count+=1           

formatted_secondary_structure = [format_structure_dict(structure) for structure in secondary_structure]

# create dataframe of cyclic peptides
cyclic_peptide_df = pd.DataFrame(list(zip(pdb, cp_type_new, a_a_length, formatted_secondary_structure, non_standard_aa_new, Caps_new, other_cp_features_new)), columns=['PDB', 'type', 'aminoacid_length', 'secondary_structure', 'non-standard_a.a.', 'CAPS', 'other features' ])

# save dafaframe as csv
cyclic_peptide_df.to_csv("cyclo_pep2023.csv")
            



