# RDKit-Custom-Sanitization
This program addresses the issue of incorrect covalent bond recognition that can occur when RDKit and MDAnalysis convert molecules into their native software representations.

MDAnalysis does not have built-in functions for sanitizing molecules, checking valences, or removing small fragmetns, as RDKit does.Therefore, these erroneous covalent bonds remain in the MDAnalysis universe unless explicitly handled. RDKit, while capable of identifying atoms with incorrect valences, lacks a sophisticated mechanism for managing such instances. This can lead to crashes during the conversion of a pdb file to an RDKit molecule.

This function establishes a custom method for managing incorrect valences between and within protein residues. It uses `amino_acid_connectivity.json`, a dictionary containing the "ground truth" covalent connectivity of all 20 amino acids and the protein backbone, to screen the local bonded neighborhood of an atom RDKit flags as having incorrect valence, and subsequently removes the incorrect connections from _BOTH_ the RDKit molecule and MD Universe objects. 

Current functionality alters the MDA and RDKit "bonded representations" of the protein only, not the underlying location of atoms in the pdb file that may be causing the sanitization failures. This program cannot resolve the misplacement of atoms, move atoms into more plausible locations, or reorient/rebuild side chains. If the valence error thrown is a result of erroneous bond degree, the function does not currently support the recharacterization of bonding degrees, although this can be added in a future update.

###########################################################################################################
# Approach Philosophy:
You will notice that the function `custom_sanitization()` is recursive. Oftentimes, when RDKit sanitization fails, the pdb file(s) in questions has many underlying structural impossibilites resulting in the recognition of erroneous covalent bonds. The Protein Data Bank is a public repository, and while efforts have been made to standardize the formatting of pdb files, many inconsistencies and errors of varying nature persist throughout. Sometimes it is impossible to resolve all heavy atoms for AA side chains from experimental electron densities. Maybe you forgot to minimize the energy of your pdb structure before performing analysis. Whatever the reason, the point is that it is unlikely that a single atom in your pdb structure will have incorrect valence. This necessitated a looping procedure to remove erroneous covalent bonds from _ALL_ atoms identified by RDKit. 

The current approach uses a recursive call to `custom_sanitization()` with a depth > 25 termination criteria. In other words, if your pdb file necessitates more than 25 rounds of custom sanitization, you have some serious problems on your hands, and you might want to take another look at your file before trying to perform analysis. I realize this is a sort of 'brute force' implementation, but it works. RDKit is an open source community, and I felt the need to contribute a solution from my own analysis because I saw one was lacking. I hope others can take this code an improve upon it!

###########################################################################################################
# USAGE:
In your process of digitizing your first molecule and performing analysis on it with either RDKit of MDAnalysis, you will do something like: 
'''
 # Create MD Universe objects from pdb_file. Utilze the in-built method for determining covalent bonds
 u = mda.Universe(pdb_file, guess_bonds=True)

 # Create rdkit molecule from pdb_file while preserving explicit hydrogen atoms
 rdkit = Chem.MolFromPDBFile(pdb_file, removeHs=False, sanitize=False) 
'''

Perhaps, as I did, you want to utilize both MDA and RDKit representations for the various useful methods present in each library. Great! As I mentioned earlier, MDAnalysis does not have built-in functions for sanitizing molecules. It is only RDKit that will complain about the valence issues if they exist. This could be an oversight if you plan to just use MDAnalysis.

Anyways, you decide to proceed with sanitization:
'''
 # Manually sanitize the molecules
 try:
   Chem.SanitizeMol(rdkit, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL)
   
 # If the sanitization fails, utilize the custom sanitization process to remove erroneous covalent bonds
 except Chem.AtomValenceException as e:
   print(f"Sanitization error: {e}")
   error_message = str(e)
   match = re.search(r'# (\d+)', error_message) # Parse the error message for the atom index throwing the error
   atom_index = int(match.group(1))
   rdkit, u = custom_sanitization(pdb_ID, rdkit, u, atom_index)
'''

And you're done! `rdkit` and `u` should now be cleaned and ready for analysis. Take a look at `RDKits_SanitizationFlags_andExplanations.txt` for an explanation of RDKit's various sanitization operations. We apply ALL in this implementation. 
