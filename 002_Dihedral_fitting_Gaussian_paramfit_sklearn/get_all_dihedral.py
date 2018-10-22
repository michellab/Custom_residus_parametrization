import Bio.PDB
import numpy
parser = Bio.PDB.PDBParser(QUIET=True) # QUIET=True avoids comments on errors in the pdb.

structures = parser.get_structure('1prot', '1prot.pdb')
structure = structures[0] # 'structures' may contain several proteins in this case only one.

def get_close_atoms(target_atom): 

	close_atoms = ns.search(target_atom.coord, 2)
	return close_atoms



parser = Bio.PDB.PDBParser(QUIET=True) # QUIET=True avoids comments on errors in the pdb.

structures = parser.get_structure('1prot', '1prot.pdb')
structure = structures[0] # 'structures' may contain several proteins in this case only one.
atoms  = Bio.PDB.Selection.unfold_entities(structure, 'A')
ns = Bio.PDB.NeighborSearch(atoms)


neighboors={}

list_bond_heavy=[]
target_atoms = [structure['A'][1]['CA']]
done=[]

while target_atoms != [] :
	target_atom=target_atoms.pop()
	list_heavy_atoms.append(atom) 
	close_atoms =get_close_atoms(target_atom)
	close_atoms.remove(target_atom)
	if target_atom.element in ['C' , 'N', 'S' , 'O'] : 
		for atm in close_atoms : if  atm.element in  ['C' , 'N', 'S' ,'O'] :
			list_bond_heavy.append([target_atom,atm])
	neighboors.update({target_atom:close_atoms})
	for a in close_atoms:  
		if a not in neighboors.keys() and target_atoms = structure['A'][1]:
			target_atoms.append(a)
	
all_dihedrals=[]
all_dihedrals_type=[]

for bond in list_bond_heavy:
	alldihedral.append(np.cartesian(( [a for a in neighboors(bond[0]) if a!=neighboors(bond[1]) ] , [a for a in neighboors(bond[1]) if a!=neighboors(bond[0]) ] )))
	alldihedral_type.append(np.cartesian(( [a.type for a in neighboors(bond[0]) if a!=neighboors(bond[1]) ] , [a.type for a in neighboors(bond[1]) if a!=neighboors(bond[0]) ] )))



