import json
import os
import re
from rdkit import Chem

# `error_log` : Stores a log of sanitization errors and their corresponding pdb_ID
error_log = os.path.join(os.getcwd(), 'error_log.txt')

def custom_sanitization(pdb_ID, molecule, u, atom_index, depth=0, max_depth=25):
    '''
    Inputs:
    `pdb_ID`: unique identifier for PPI
    `molecule`: RDKit molecular representation we want to modify
    `u`: MDAnalysis Universe object we want to modify
    `atom index`: Unique atom index from RDKit where the valence error occurs. You obtain this by parsing the error message string
    `depth`: Current depth for recursion
    `max_depth`: breaking condition for recursion. If the function fails to sanitize your molecule after 25 iterations, you have some
                 serious issues with your file and custom sanitization fails. Sorry!

    Outputs:
    `molecule`: The updated RDKit molecule with erroneous covalent bonds removed
    `u`: The updated MDAnalysis Universe object with erroneous covalent bonds removed.
    '''
    # The problematic_atom has a valence greater than permitted due to erroneous bond(s)
    problematic_atom = molecule.GetAtomWithIdx(atom_index)
    print(problematic_atom.GetFormalCharge())

    # Define res_info in the RDKit ecosystem for the residue containing the problematic_atom
    res_info = problematic_atom.GetPDBResidueInfo()
    residue = res_info.GetResidueName()
    res_id = res_info.GetResidueNumber()
    chain_id =  res_info.GetChainId()
    atom_name = res_info.GetName().strip()
    tag = f'{residue}{res_id}.{chain_id}'
    
    # Find the equivalent atom in the MDA Universe object `u`
    mda_atom = u.select_atoms(f"resname {residue} and resid {res_id} and segid {chain_id} and name {atom_name}")
    
    # Loop through the bonds associated with the problematic_atom and 
    # create a mapping of {connection: atom_indices} like {'O <-> C': '768<->767', 'OG <-> O': '778<->768'}
    mapping = {}
    for bond in problematic_atom.GetBonds():
        atom1 = molecule.GetAtomWithIdx(bond.GetBeginAtomIdx())
        atom2 = molecule.GetAtomWithIdx(bond.GetEndAtomIdx())
        atom_indices = f'{atom1.GetIdx()}<->{atom2.GetIdx()}'

        res_info1 = atom1.GetPDBResidueInfo()
        res_info2 = atom2.GetPDBResidueInfo()
        connection = f'{res_info1.GetName().strip()} <-> {res_info2.GetName().strip()}'
        
        mapping[connection] = atom_indices

    # Load in the json file where valid amino_acid_connectivity information is stored
    with open('data_prep/amino_acid_connectivity.json', 'r') as f:
        bond_dict = json.load(f)

    # Obtain the valid connectivity for the residue
    valid_connectivity = bond_dict.get(residue, [])
    
        # Obtain valid backbone bonds
    backbone_bonds = bond_dict.get('Backbone', [])
    
    # Obtain the indices for all atoms in the RDKit residue
    rdkit_residue_indices = []
    for atom in molecule.GetAtoms():
        if (atom.GetPDBResidueInfo().GetResidueName() == residue and
            atom.GetPDBResidueInfo().GetResidueNumber() == res_id and
            atom.GetPDBResidueInfo().GetChainId() == chain_id):
            rdkit_residue_indices.append(atom.GetIdx())
    print(f'RDKit Indices of [{tag}]: {rdkit_residue_indices}')

    # Obtain the indices for all atoms in the MDA residue
    residue_of_atom = mda_atom.residues
    atoms_in_residue = residue_of_atom.atoms
    mda_residue_indices = atoms_in_residue.indices
    print(f'MDA Indices of [{tag}]: {mda_residue_indices}\n')
          

    # Remove erroneous bonds from both 'u' and 'molecule' 
    # either within the same residue or between neighboring residues
    # based on camparison to the `valid_connectivity` list and the indices of the bonded atoms
    for connection, indices in mapping.items():
        index_source = int(indices.split('<->')[0])
        index_target = int(indices.split('<->')[1])

        # Invalid intra-residue connection
        if (index_source in rdkit_residue_indices) and (index_target in rdkit_residue_indices):
            if (connection not in valid_connectivity):
                if connection in backbone_bonds:
                    if connection == "N <-> C" or connection == "C <-> N":
                        # Remove erroneous backbone-backbone connection
                        molecule = update_rdkit(molecule, index_source, index_target, connection, reason=f'Erroneous backbone-backbone bond')
                        u = update_mda(u, connection, mda_atom, res_id, reason=f'Erroneous backbone-backbone bond')
                        print('-------------------------------------------------------------------------------------------------------------------------')
                    else:  
                        continue # DO NOT remove the valid backbone-backbone connections
                else:
                    # Remove erroneous sidechain-backbone and sidechain-sidechain bonds 
                    molecule = update_rdkit(molecule, index_source, index_target, connection, reason=f'Connection not found in {residue}')
                    u = update_mda(u, connection, mda_atom, res_id, reason=f'Connection not found in {residue}')
                    print('-------------------------------------------------------------------------------------------------------------------------')

        # Invalid cross residue connection - "N <-> C" and "C <-> N" are the only valid residue linking bonds
        if (index_source not in rdkit_residue_indices) or (index_target not in rdkit_residue_indices):
            if (connection != "N <-> C") and (connection != "C <-> N"): 
                # Remove erroneous backbone-backbone connection
                molecule = update_rdkit(molecule, index_source, index_target, connection, reason=f'Erroneous cross residue connection')
                u = update_mda(u, connection, mda_atom, res_id, reason=f'Erroneous cross residue connection')
                print('-------------------------------------------------------------------------------------------------------------------------')
    
    # Try the sanitization process again
    try:
        Chem.SanitizeMol(molecule, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL)
        print('Successfully sanitized RDKit molecule via custom process.')
        return molecule, u
    except Chem.AtomValenceException as e:
        print(f"Sanitization error: {e}")
        error_message = str(e)
        match = re.search(r'# (\d+)', error_message)
        atom_index = int(match.group(1))
        
        if depth > max_depth:
            # If the maximum depth is reached, stop the recursion
            print("[ERROR] Maximum recursion depth reached, sanitization failed.")
            with open(error_log, "a+") as file:
                # Go to the beginning of the file in case there is existing content
                file.seek(0)
                file.write(pdb_ID + "\n" + str(e) + "\n")
        else:
            # Recursive call to clear up additional downstream sanitization issues
            return custom_sanitization(pdb_ID, molecule, u, atom_index, depth+1)


def update_rdkit(molecule, index_source, index_target, connection, reason=None):
    # Convert the molecule to a RWMol
    rw_mol = Chem.RWMol(molecule)

    # Now you can remove the bond
    rw_mol.RemoveBond(index_source, index_target)
    print(f'[RDKit] Removed bond: {connection} between atom {index_source} and atom {index_target} Reason: {reason}')

    # If you want to convert it back to a regular molecule
    molecule = rw_mol.GetMol()
    
    return molecule

def update_mda(u, connection, mda_atom, res_id, reason=None):
    '''
    The `else` clause in a `for` loop is executed when the loop has 
    exhausted iterating the list. If a `break` statement is executed 
    inside the `for` loop then the `else` block is skipped. 
    '''
    valid_connections = []
    for bond in mda_atom.bonds:
        atom1, atom2 = bond.atoms
        if atom1.name == connection.split(' <-> ')[0] and atom2.name == connection.split(' <-> ')[1]:
            if (reason == 'Erroneous backbone-backbone bond') and atom1.residue.resid == atom2.residue.resid:
                    found = True
                    u.delete_bonds([bond])
                    print(f'[MDA Universe] Removed bond: {atom1.name} <-> {atom2.name} between atom {atom1.ix} and atom {atom2.ix} Reason: {reason}')
                    break      
            else:
                found = True
                u.delete_bonds([bond])
                print(f'[MDA Universe] Removed bond: {atom1.name} <-> {atom2.name} between atom {atom1.ix} and atom {atom2.ix} Reason: {reason}')
                break
        elif atom2.name  == connection.split(' <-> ')[0] and atom1.name == connection.split(' <-> ')[1]:
            if (reason == 'Erroneous backbone-backbone bond') and atom2.residue.resid == atom1.residue.resid:
                found = True
                u.delete_bonds([bond])
                print(f'[MDA Universe] Removed bond: {atom2.name} <-> {atom1.name} between atom {atom2.ix} and atom {atom1.ix} Reason: {reason}')
                break
            else:
                found = True
                u.delete_bonds([bond])
                print(f'[MDA Universe] Removed bond: {atom2.name} <-> {atom1.name} between atom {atom2.ix} and atom {atom1.ix} Reason: {reason}')
                break
        else:
            valid_connections.append(f'{atom1.name} <-> {atom2.name}')
    else:
        print(f'[MDA Universe] {connection} not found for the problematic atom [{mda_atom[0].name}] with connections: {valid_connections}')
        return u
    
    return u
