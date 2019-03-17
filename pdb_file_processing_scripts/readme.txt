We will have a set of scripts which will parse the PDB files provided by a user to get the all-atom files and C-alpha files to be used later in our calculations.
The native PDB file provided by the user will be processed to extract the chain of interest as follows.
1. Remove all ligand/Het atoms
2. Remove alternate locations for residues (if any) by retaining only the 'A' location
3. Make sure that all the residue names are in the set of the standard 20 amino acid names.
 
