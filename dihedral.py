from MDAnalysis import Universe
from MDAnalysis.analysis.dihedrals import Dihedral
import pandas as pd
import sys

#Note: 
# - MDAnalysis gives direct computation of phi, omega, and chi 1 angles
# - See https://userguide.mdanalysis.org/1.0.0/examples/analysis/structure/dihedrals.html
# - AtomGroup object can be created by "u.select_atoms()" 
# - See https://docs.mdanalysis.org/1.0.0/documentation_pages/analysis/dihedrals.html#MDAnalysis.analysis.dihedrals.Dihedral.angles

def chi2_selection(u):
    """
    Args
      u: MDAnalysis.Universe
    Return:
      atomgroups: AtomGroups in MDAnalysis
    """
    res_letters  = ['PHE','TYR','TRP','HIS','MET','LEU','ILE','ARG','GLU','GLN','LYS']    
    chi2_atoms   = {'PHE':['CA', 'CB', 'CG' , 'CD1'],
                    'TYR':['CA', 'CB', 'CG' , 'CD1'],
                    'TRP':['CA', 'CB', 'CG' , 'CD1'],
                    'HIS':['CA', 'CB', 'CG' , 'ND1'],
                    'MET':['CA', 'CB', 'CG' , 'SD' ],
                    'LEU':['CA', 'CB', 'CG' , 'CD1'],
                    'ILE':['CA', 'CB', 'CG1', 'CD'],
                    'ARG':['CA', 'CB', 'CG' , 'CD' ],
                    'GLU':['CA', 'CB', 'CG' , 'CD' ],
                    'GLN':['CA', 'CB', 'CG' , 'CD' ],
                    'LYS':['CA', 'CB', 'CG' , 'CD' ]}

    #Get residue ids and 3-letter
    res_sni = [] # stands for RESidue Segment, Name, and Index
    for residue in u.residues:
        resn = residue.resname
        if resn in res_letters:
            res_sni.append([resn, residue.resid, residue.segid])
    
    #Make a list of atomgroups
    atomgroups = []
    for resn, resid, chainid in res_sni:
        atom_names = ' '.join(chi2_atoms[resn])
        atomgroups.append(u.select_atoms(f'segid {chainid} and resid {resid} and name {atom_names}'))

    n_atoms_in_each = set(map(len, atomgroups))
    for a in atomgroups:
        if len(a) != 4: 
            print(f'Stop) All Atom groups must have four atoms, but the No. of atoms stored are:')
            print(a.resnames, a.names)
            exit()
    #if len(n_atoms_in_each) != 1: 
    #    print(f'Stop) All Atom groups must have four atoms, but the No. of atoms stored are: {n_atoms_in_each}')
    #    exit()
    
    return atomgroups 

def rotCH3(u):
    """
    Args
      u: MDAnalysis.Universe
    Return:
      atomgroups: AtomGroups in MDAnalysis
    
    Note: This function is for dihedral angle of N-Ca-Cb-HB1 of Ala.  
    """
    res_letters  = ['ALA', 'MET', 'LEU', 'ILE', 'THR', 'VAL']    
    chi2_atoms   = {'ALA':['N', 'CA', 'CB' , 'HB1'],
                    'MET':['CG','SD','CE','HE1'],
                    'LEU':['CB','CG','CD1','HD11'], #or ['CB','CG','CD2','HD21']
                    'ILE':['CB','CG1','CD','HD1'],  #or ['CA','CB','CG2','HG21']
                    'THR':['CA','CB','CG2','HG21'],
                    'VAL':['CA','CB','CG1','HG11']  #or ['CA','CB','CG2','HG21']
                    }

    #Get residue ids and 3-letter
    res_sni = [] # stands for RESidue Segment, Name, and Index
    for residue in u.residues:
        resn = residue.resname
        if resn in res_letters:
            res_sni.append([resn, residue.resid, residue.segid])
    
    #Make a list of atomgroups
    atomgroups = []
    for resn, resid, chainid in res_sni:
        atom_names = ' '.join(chi2_atoms[resn])
        atomgroups.append(u.select_atoms(f'segid {chainid} and resid {resid} and name {atom_names}'))

    n_atoms_in_each = set(map(len, atomgroups))
    for a in atomgroups:
        if len(a) != 4: 
            print(f'Stop) All Atom groups must have four atoms, but the No. of atoms stored are:')
            print(a.resnames, a.names)
            exit()
    
    return atomgroups     

def parser():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('-r','--ref', required=True) 
    p.add_argument('-t','--trj')
    p.add_argument('-at','--angle_type', default='X2')

    args = p.parse_args()
    return args.ref, args.trj, args.angle_type

def make_columns(atomgroups):
    selected_residues = []
    print('---Selected residues and atoms for dihedral angle calculation---')
    
    for i in atomgroups:
        segid = i.segids[0]
        if segid == 'SYSTEM': segid = 'A'

        selected_residues.append(f'{i.resnames[0]}{i.resids[0]}{segid}')
        print(i.resnames[0], i.resids[0], segid, *i.names)

    return selected_residues

def main():
    ref, traj, at  = parser()
    u = Universe(ref, traj) 

    if   at=='X2':
        atomgroups = chi2_selection(u)

    elif at=='X1':
        atomgroups = [res.chi1_selection() for res in u.residues]
        atomgroups = [ele for ele in atomgroups if ele] # this is for eliminating None element from atomgroups. "if ele" works only if ele is not empty (None).

    elif at=='CH3':
        atomgroups = rotCH3(u)

    columns = make_columns(atomgroups) #for table's column

    # Each atom group must contain 4 atoms.
    R = Dihedral(atomgroups).run() # values are stored in R.angles
    df = pd.DataFrame(R.angles, columns = columns)
    df.to_csv(f'{at}.csv')

if __name__ == '__main__':
    main()
