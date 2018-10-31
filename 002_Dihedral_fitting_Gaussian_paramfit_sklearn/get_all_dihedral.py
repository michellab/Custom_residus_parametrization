import Bio.PDB
import numpy
#import MDAnalysis as MDA
import parmed
import os
#def get_close_atoms(target_atom):

	#close_atoms = ns.search(target_atom.coord, 2)
	#return close_atoms




def make_missing_parms_parmck2(inputfile):
    os.system( 'parmchk2 -i  %s  -o frcmod -f mol2 -s 2 -o amber.frcmod' %(inputfile))
    tleapinput=open('temp.in','w')
    tleapinput.writelines(['loadAmberParams amber.frcmod \n',  'source leaprc.protein.ff14SB \n'  ,  'input = loadmol2 %s \n '  %(inputfile) ,\
    'saveAmberParm  input  input.prmtop input.inpcrd\n ', 'quit'])
    tleapinput.close()
    os.system( 'tleap -f temp.in' )

def find_dihedrals(inputfile):
    make_missing_parms_parmck2(inputfile)
    all_dihedrals=[]
    dihedrals_heavy=[]
    dihedrals_heavy_index=[]
    dihe=[]
    heavy_bonds=[]
    dihedrals_heavy_centralbond_name=[]
    dihedrals_heavy_name=[]
    topol = parmed.load_file( 'input.prmtop')
    for resid in topol.residues:
             for atom in resid.atoms:

                 if atom.atomic_number in [6, 7, 8, 16]: # C, N , O and S

                    for bond in atom.bonds :
                        b = [bond.atom1.idx , bond.atom2.idx]
                        b.sort()
                        if b not in heavy_bonds and bond.atom1.atomic_number in [6, 7, 8, 16] and   bond.atom2.atomic_number in[6, 7, 8, 16] :
                            heavy_bonds.append(b)


    dihe= topol.dihedrals

    for bond in heavy_bonds:
        all_dihedrals.append([])

        dihedrals_heavy.append([])
        for i  in range (len(dihe)) :

             d = [dihe[i].atom2.idx,dihe[i].atom3.idx]
             d.sort()

             if bond ==d :

                 if dihe[i] not in  all_dihedrals[-1]:
                     all_dihedrals[-1].append(dihe[i])
                     ans = [' "%s" "%s" "%s" "%s" ' %( dihe[i].atom1.name, dihe[i].atom2.name,dihe[i].atom3.name,dihe[i].atom4.name)]
                     an = [ dihe[i].atom1.name, dihe[i].atom2.name,dihe[i].atom3.name,dihe[i].atom4.name]

                     ai = [ dihe[i].atom1.idx, dihe[i].atom2.idx,dihe[i].atom3.idx,dihe[i].atom4.idx]
                     #ai = [ ' "%s" "%s" "%s" "%s" ' %(dihe[i].atom1.idx, dihe[i].atom2.idx,dihe[i].atom3.idx,dihe[i].atom4.idx)]

                     if dihe[i].atom1.atomic_number in [6, 7, 8, 16] and  dihe[i].atom4.atomic_number in [6, 7, 8, 16] and ai not in  dihedrals_heavy_index :
                            dihedrals_heavy[-1].append([dihe[i]])
                            if dihe[i].atom2.type != 'CA' and dihe[i].atom3.type != 'CA'  :  # no aromatic or it will be a mess!
                                dihedrals_heavy_name.append( ans)
                                dihedrals_heavy_centralbond_name.append([an[1], an[2]])
                                dihedrals_heavy_index.append(ai)
    torsion_names=[]
    i=1
    all_dihedrals_type=[]
    for bond in dihedrals_heavy_centralbond_name:
        bond.sort()
        all_dihedrals_type.append([])
        for i  in range (len(dihe)) :
            d = [dihe[i].atom2.idx,dihe[i].atom3.idx]
            d.sort()
            if bond ==d :
                at = [dihe[i].atom1.type,dihe[i].atom2.type,dihe[i].atom3.type,dihe[i].atom4.type,dihe[i].multiplicity]
                if at not in all_dihedrals_type[-1] and ai not in  dihedrals_heavy_index:
                    all_dihedrals_type[-1].append(at)

        if ['CA' ,'N'] == bond:
            torsion_names.append('PHI')
        elif [ "C" , "CA" ] == bond :
            torsion_names.append('PSI')
        else :
            torsion_names.append('CHI'+str(i) )
            i+=1
    return all_dihedrals,   all_dihedrals_type, dihedrals_heavy, dihedrals_heavy_name, torsion_names, dihedrals_heavy_index


if __name__ == '__main__' :




    make_missing_parms_parmck2('Mol-sm_m1-c1.mol2')






















'''
#   traj = md.load('PSI_-70.mol2')
#    inputfile= 'capped-CLS-input-gaussian.mol2'
#    inputfile='capped-CLS-input-gaussian.pdb'
    inputfile='Mol-sm_m1-c1.mol2'
#    parser = Bio.PDB.PDBParser(QUIET=True) # QUIET=True avoids comments on errors in the pdb.
#    structures = parser.get_structure('input', inputfile)
#    structure = structures[0] # 'structures' may contain several proteins in this case only one.
#    atoms  = Bio.PDB.Selection.unfold_entities(structure, 'A')
#    ns = Bio.PDB.NeighborSearch(atoms)

    neighboors={}
    list_bond_heavy=[]
#    target_atoms = [structure['*'][1]['CA']]
    target_atoms=['CA']
    done=[]
    topology= MDA.Universe(inputfile)

    #topology= MDA.coordinates.MOL2.MOL2Reader(inputfile)
    list_heavy_atoms=[]
    while target_atoms != [] :
        target_atom=target_atoms.pop()
        list_heavy_atoms.append(target_atom)
        #close_atoms =get_close_atoms(target_atom)
        close_atoms=topology.select_atoms('sphzone 2.0 (name %s) and not name %s ' %(target_atom, target_atom) )

        if target_atom[0] in ['C' , 'N', 'S' , 'O'] :
            for atm in close_atoms :
                if  atm.element in  ['C' , 'N', 'S' ,'O'] :
                    list_bond_heavy.append([target_atom,atm])
                neighboors.update({target_atom:close_atoms})
        for a in close_atoms:
            if a not in neighboors.keys() and target_atoms == structure['A'][1]:
                target_atoms.append(a)

    all_dihedrals=[]
    all_dihedrals_type=[]

    for bond in list_bond_heavy:
    	alldihedral.append(np.cartesian(( [a for a in neighboors(bond[0]) if a!=neighboors(bond[1]) ] , [a for a in neighboors(bond[1]) if a!=neighboors(bond[0]) ] )))
    	alldihedral_type.append(np.cartesian(( [a.type for a in neighboors(bond[0]) if a!=neighboors(bond[1]) ] , [a.type for a in neighboors(bond[1]) if a!=neighboors(bond[0]) ] )))
'''
