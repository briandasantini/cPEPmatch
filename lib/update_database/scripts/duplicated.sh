#!/bin/bash

# Read PDB IDs from both files into arrays
mapfile -t batch23_pdb < batch23/pdbs.txt
mapfile -t batch19_pdb < batch19/pdbs.txt

# Concatenate arrays and sort them
all_pdbs=("${batch23_pdb[@]}" "${batch19_pdb[@]}")
sorted_pdbs=($(echo "${all_pdbs[@]}" | tr ' ' '\n' | sort))

# Find duplicates
echo "Duplicates:"
uniq -d <<< "${sorted_pdbs[@]}"
