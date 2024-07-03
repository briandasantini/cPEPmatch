**Overview**\
This project contains a Jupyter Notebook designed to expand the preexisting cyclic peptide database. This program analyses and characterises new cyclic peptides structures based on various attributes such as secondary structure, length, connectivity and the presence of non-standard amino acids and caps. 
The notebook processes multiple PDB files, extracts relevant features, and compiles the data into a comprehensive CSV file.

**Requirements**\
To run this notebook, you need the following:

+ **Python**: Version 3.7 or higher.
+ **Jupyter Notebook**: An interactive environment to run the notebook.
+ **Libraries**:
  + pandas
  + numpy
  + BioPython
  + mdtraj

These libraries can be installed via pip (bash code):
pip install pandas numpy biopython

**Files**
+ Cpep_analysis.ipynb: The main Jupyter Notebook containing the code for analysing cyclic peptides.
+ PDB files: The input data files of potential cyclic peptides in PDB format located in specified directories.

**Structure of the Notebook**\
The notebook is divided into several sections, each responsible for a specific part of the analysis:

1. Imports and Setup:  Importing necessary libraries and setting up the initial configurations.
2. Defining Functions: Functions to process PDB files, extract features, and analyse cyclic peptides.
3. Processing PDB Files: Iterating over PDB files in specified directories to extract and compute features.
4. Creating the DataFrame: Compiling extracted features into a Pandas DataFrame.
5. Saving Results: Saving the DataFrame to a CSV file for further analysis.

**Key Functions**
+ find_secondary_structure(pdb_path): Analyses the secondary structure of a cyclic peptide from a PDB file.
+ list_comprehension(data_list): Combines substrings in a sublist for cleaner data representation.
+ process_pdb_files(root_dir): Main function to iterate over PDB files in the specified directory, extract features, and compile them into lists.
+ remove_pdb(batch_dir): Function to clean up and remove temporary PDB files e.g. of non cyclic peptides.

**Usage**
1. Setup Directory: Place all your PDB files in the specified directory structure.
2. Run Notebook: Open Cpep_analysis.ipynb in Jupyter Notebook and run all cells sequentially.
3. Review Output: The notebook will process the PDB files, extract features, and save the results to cyclo_pep.csv.

**Example**\
To run the notebook on a new set of PDB files:

1. Place the PDB files in the directory /home/yourusername/projects/cPEP_library/batchX/renum/. (insert whatever pathway leads to the input files)
2. Update the batch variable in the notebook to point to the correct directory.
3. Execute all cells in the notebook.
4. Check the generated cyclo_pep.csv for the results.

**Conclusion**\
This notebook provides a comprehensive tool for characterising cyclic peptides from PDB files. By following the setup and usage instructions, users can efficiently analyse their datasets and obtain detailed information on peptide properties.
