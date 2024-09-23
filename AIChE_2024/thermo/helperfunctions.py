from rmgpy import settings
from rmgpy.data.thermo import ThermoDatabase, ThermoData, remove_thermo_data, add_thermo_data
import rmgpy.molecule.group as gr
from rmgpy.molecule.group import Group, GroupAtom, GroupBond
from rmgpy.data.base import Entry
from rmgpy.data.base import LogicOr
from rmgpy.molecule.atomtype import ATOMTYPES

import logging
import os
import sys
import numpy as np
np.set_printoptions(threshold=sys.maxsize)
from copy import copy, deepcopy
import itertools
import matplotlib.pyplot as plt
import matplotlib as mpl
plt.style.use("seaborn-poster")
import numpy as np

from copy import deepcopy
from sklearn.linear_model import RidgeCV,LassoCV,ElasticNetCV,LinearRegression
from sklearn.linear_model import Ridge,Lasso,ElasticNet
from sklearn.metrics import mean_squared_error, mean_absolute_error


def check_data(label=None,index=None):
    """
    Helper functions for checking outlier species in the parity plot. 
    User can input either the label or the index of a species.
    This function prints out
    - GAE with estimated GAV: group additivity estimation (GAE) with not exact fit node 
                              (node data is string or None)
    - missing groups
    - GAE without estimated GAV: group additivity estimation with only the exact fit node
    - GAE with fitted GAV: GAE with newly fitted GAVs. 
                           Calculated by directly adding the new groups to the GAE without estimated GAV value.
                           Does not descend the tree again.
    - New GAE: group additivity estimation using the new database. Descending from the new tree.
    """
    if index is None:
        index = label_index_dict[label]
    
    entry = all_data["entry"][index]
    display(entry.item)
    print(entry.item.smiles)
    entry_thermo = all_data["entry thermo"][index]
    
    print("Entry")
    print(entry.short_desc)
    print(entry.label)
    print(entry_thermo)
    print('')
    print("GAE with estimated GAV")
    print(all_data["GAE with estimated GAV"][index])
    print('')
    print("missing")
    print(all_data["missing groups"][index])
    print('')
    print('GAE without estimated GAV')
    print(all_data["GAE without estimated GAV"][index])
    print('')
    print("GAE with fitted GAV")
    print(all_data["GAE with fitted GAV"][index])
    print('')
    if "New GAE" in all_data:
        print("New GAE")
        print(all_data["New GAE"][index])
    
def add_thermo_data_uncertainty(thermo_data1, thermo_data2, group_additivity=True, verbose=False):
    """
    Adapted from the add_thermo_data function from `rmgpy.data.thermo`.
    Add the thermodynamic data, including the uncertainties, from `thermo_data2` to the data `thermo_data1`,
    and return `thermo_data1`.

    If `group_additivity` is True, append comments related to group additivity estimation
    If `verbose` is False, omit the comments from a "zero entry", whose H298, S298, and Cp are all 0.
    If `verbose` is True, or thermo_data2 is not a zero entry, add thermo_data2.comment to thermo_data1.comment.
    """
    if (len(thermo_data1.Tdata.value_si) != len(thermo_data2.Tdata.value_si) or
            any([T1 != T2 for T1, T2 in zip(thermo_data1.Tdata.value_si, thermo_data2.Tdata.value_si)])):
        raise ValueError('Cannot add these ThermoData objects due to their having different temperature points.')

    for i in range(thermo_data1.Tdata.value_si.shape[0]):
        thermo_data1.Cpdata.value_si[i] += thermo_data2.Cpdata.value_si[i]
        thermo_data1.Cpdata.uncertainty_si[i] += thermo_data2.Cpdata.uncertainty_si[i]
    thermo_data1.H298.value_si += thermo_data2.H298.value_si
    thermo_data1.H298.uncertainty_si += thermo_data2.H298.uncertainty_si
    thermo_data1.S298.value_si += thermo_data2.S298.value_si
    thermo_data1.S298.uncertainty_si += thermo_data2.S298.uncertainty_si

    test_zero = sum(abs(value) for value in
                    [thermo_data2.H298.value_si, thermo_data2.S298.value_si] + thermo_data2.Cpdata.value_si.tolist())
    # Used to check if all of the entries in thermo_data2 are zero

    if group_additivity:
        if verbose or test_zero != 0:
            if thermo_data1.comment:
                thermo_data1.comment += ' + {0}'.format(thermo_data2.comment)
            else:
                thermo_data1.comment = 'Thermo group additivity estimation: ' + thermo_data2.comment

    return thermo_data1

def get_neighbors(atom, parent_node, group_atoms, n_degree_neighbor=1):
    """
    Get neighbors within n degree for the center atom recursively.
    
    Args:
        atom (Atom): the center atom
        group_atoms (dict): an empty dictionary
        n_degree_neighbor (int): If a neighbor is n bond apart from the center atom, 
                                    it's defined as n degree neighbor.
        
    
    Returns:
        group_atoms (dict): a dictionary that has Atom object as key and the corresponding 
                            GroupAtom object as value
    """
    #let's look at the specifics of the atoms on the parent node
    on_carbon_lp = False
    on_carbon_c = False
    on_fluorine_lp = False
    on_fluorine_c = False
    on_oxygen_lp = False
    on_oxygen_c = False
    on_sulfur_lp = False
    on_sulfur_c = False
    on_chlorine_lp = False
    on_chlorine_c = False
    on_bromine_lp = False
    on_bromine_c = False


    for atm in parent_node.item.atoms:
        if atm.atomtype[0].label not in ["O2s", "S2s", "O4tc", "Cs"]:
            if atm.is_carbon(): 
                if atm.lone_pairs!=[]:  # if the parent node does not specify lone pairs, then the child shouldn't have to
                    on_carbon_lp = True
                if atm.charge!=[]:  # if the parent node does not specify charge, then the child shouldn't have to
                    on_carbon_c = True
            if atm.is_fluorine(): 
                if atm.lone_pairs!=[]:  # if the parent node does not specify lone pairs, then the child shouldn't have to
                    on_fluorine_lp = True
                if atm.charge!=[]:  # if the parent node does not specify charge, then the child shouldn't have to
                    on_fluorine_c = True #but if the parent group is specific, make the child group specific
            if atm.is_oxygen(): 
                if atm.lone_pairs!=[]:  # if the parent node does not specify lone pairs, then the child shouldn't have to
                    on_oxygen_lp = True
                if atm.charge!=[]:  # if the parent node does not specify charge, then the child shouldn't have to
                    on_oxygen_c = True #but if the parent group is specific, make the child group specific
            if atm.is_sulfur(): 
                if atm.lone_pairs!=[]:  # if the parent node does not specify lone pairs, then the child shouldn't have to
                    on_sulfur_lp = True
                if atm.charge!=[]:  # if the parent node does not specify charge, then the child shouldn't have to
                    on_sulfur_c = True #but if the parent group is specific, make the child group specific
            if atm.is_chlorine(): 
                if atm.lone_pairs!=[]:  # if the parent node does not specify lone pairs, then the child shouldn't have to
                    on_chlorine_lp = True
                if atm.charge!=[]:  # if the parent node does not specify charge, then the child shouldn't have to
                    on_chlorine_c = True #but if the parent group is specific, make the child group specific
            if atm.is_bromine(): 
                if atm.lone_pairs!=[]:  # if the parent node does not specify lone pairs, then the child shouldn't have to
                    on_bromine_lp = True
                if atm.charge!=[]:  # if the parent node does not specify charge, then the child shouldn't have to
                    on_bromine_c = True #but if the parent group is specific, make the child group specific
                                
                    
                    
    #save the flag values and pass them on
    flags = [on_carbon_lp, 
            on_carbon_c,
            on_fluorine_lp,
            on_fluorine_c,
            on_oxygen_lp,
            on_oxygen_c,
            on_sulfur_lp,
            on_sulfur_c,  
            on_chlorine_lp,
            on_chlorine_c,
            on_bromine_lp,
            on_bromine_c]
    
    # print(f'In get_neighbers, flags based on parent node are: {flags}')
    if n_degree_neighbor == 1:
    
        for atm in atom.edges:
            if atm not in group_atoms:
                # print(atm.atomtype)
                if atm.is_carbon(): 
                    #let's only worry about the first two flags
                    flags_to_look_at = [flags[0], flags[1]]

                if atm.is_fluorine(): 
                    flags_to_look_at = [flags[2], flags[3]]

                if atm.is_oxygen(): 
                    flags_to_look_at = [flags[4], flags[5]]

                #hydrogens should never have lone pairs or charge in a molecule
                if not atm.is_non_hydrogen():
                    flags_to_look_at = [False, False]

                if atm.is_sulfur(): 
                    #let's only worry about the last two flags
                    flags_to_look_at = [flags[6], flags[7]]
                #editing here for chlorine: 
                if atm.is_chlorine():
                    flags_to_look_at = [flags[8], flags[9]]                    
                if atm.is_bromine():
                    flags_to_look_at = [flags[10], flags[11]]                     
                    


                if flags_to_look_at == [True, True]: #lp, charge both True on carbon
                    group_atoms[atm] = GroupAtom(atomtype=[atm.atomtype],
                                                radical_electrons=[atm.radical_electrons],
                                                lone_pairs=[atm.lone_pairs],
                                                charge=[atm.charge],
                                                label='')
                if flags_to_look_at == [True, False]: #only lp True on carbon
                    group_atoms[atm] = GroupAtom(atomtype=[atm.atomtype],
                                                radical_electrons=[atm.radical_electrons],
                                                lone_pairs=[atm.lone_pairs],
                                                label='')
                if flags_to_look_at == [False, True]: #only charge on carbon
                    group_atoms[atm] = GroupAtom(atomtype=[atm.atomtype],
                                                radical_electrons=[atm.radical_electrons],
                                                charge=[atm.charge],
                                                label='')
                if flags_to_look_at == [False, False]: #lp, charge both False
                    group_atoms[atm] = GroupAtom(atomtype=[atm.atomtype],
                                                radical_electrons=[atm.radical_electrons],
                                                label='')                    

                
    else:
        
        for atm in atom.edges:
            if atm not in group_atoms:
                if atm.is_carbon(): 
                    #let's only worry about the first two flags
                    flags_to_look_at = [flags[0], flags[1]]

                if atm.is_fluorine(): 
                    #let's only worry about the last two flags
                    flags_to_look_at = [flags[2], flags[3]]

                if atm.is_oxygen(): 
                    #let's only worry about the last two flags
                    flags_to_look_at = [flags[4], flags[5]]

                #hydrogens should never have lone pairs or charge in a molecule
                if not atm.is_non_hydrogen():
                    flags_to_look_at = [False, False]

                if atm.is_sulfur(): 
                    #let's only worry about the last two flags
                    flags_to_look_at = [flags[6], flags[7]]
                #editing here for chlorine: 
                if atm.is_chlorine():
                    flags_to_look_at = [flags[8], flags[9]] 
                if atm.is_bromine():
                    flags_to_look_at = [flags[10], flags[11]] 

                if flags_to_look_at == [True, True]: #lp, charge both True on carbon
                    group_atoms[atm] = GroupAtom(atomtype=[atm.atomtype],
                                                radical_electrons=[atm.radical_electrons],
                                                lone_pairs=[atm.lone_pairs],
                                                charge=[atm.charge],
                                                label='')
                if flags_to_look_at == [True, False]: #only lp True on carbon
                    group_atoms[atm] = GroupAtom(atomtype=[atm.atomtype],
                                                radical_electrons=[atm.radical_electrons],
                                                lone_pairs=[atm.lone_pairs],
                                                label='')
                if flags_to_look_at == [False, True]: #only charge on carbon
                    group_atoms[atm] = GroupAtom(atomtype=[atm.atomtype],
                                                radical_electrons=[atm.radical_electrons],
                                                charge=[atm.charge],
                                                label='')
                if flags_to_look_at == [False, False]: #lp, charge both False
                    group_atoms[atm] = GroupAtom(atomtype=[atm.atomtype],
                                                radical_electrons=[atm.radical_electrons],
                                                label='')   
            assert  n_degree_neighbor<10
            get_neighbors(atm, parent_node, group_atoms, n_degree_neighbor-1)
                
    return group_atoms
        
def make_bonds(atom, group, group_atoms, n_degree_neighbor=1):
    """
    Make bonds between the group atoms in the group object recursively, using the edges of the center atom
    Arg:
        atom (Atom): the center atom for the group
        group (Group): the group object where GroupAtoms are filled but GroupBonds are not added yet
        group_atoms (dict): a dictionary that has Atom object as key and the corresponding 
                            GroupAtom object as value
        n_degree_neighbor (int): If a neighbor is n bond apart from the center atom, 
                                    it's defined as n degree neighbor.
    Returns:
        group (Group): the group object with GroupBonds added
    """
    
    if n_degree_neighbor == 1:
        
        for bonded_atom, bond in atom.edges.items():
            if not group.has_bond(group_atoms[atom],group_atoms[bonded_atom]):
                group.add_bond(GroupBond(group_atoms[atom],group_atoms[bonded_atom],order=[bond.order]))
            else:
                pass
                
    else:
        
        for bonded_atom, bond in atom.edges.items():
            if not group.has_bond(group_atoms[atom],group_atoms[bonded_atom]):
                group.add_bond(GroupBond(group_atoms[atom],group_atoms[bonded_atom],order=[bond.order]))
            else:
                pass
            make_bonds(bonded_atom, group, group_atoms, n_degree_neighbor-1)
            
    return group

def make_group(atom, parent_node, n_degree_neighbor=1):
    """
    Make a group with `atom` as the center atom and `n_degree_neighbor` degree of neighbor, 
    using `get_neighbors` and `make_bonds` functions.
    
    Args:
        atom (Atom): the center atom for the group
        n_degree_neighbor (Int): If a neighbor is n bond apart from the center atom, 
                                    it's defined as n degree neighbor.
    Returns:
        group (Group): the group containing the center atom and the n degrees of neighbor
    """

    group_atoms = {}
    bonds = []
    
    # print(f'This is the parent node\n{parent_node.item.to_adjacency_list()}')
    group_atoms = get_neighbors(atom, parent_node, group_atoms, n_degree_neighbor=n_degree_neighbor) #this is a dictionary
    
    #now let's look at each of the group atoms 

    #let's look at the specifics of the atoms on the parent node
    on_carbon_lp = False
    on_carbon_c = False
    on_fluorine_lp = False
    on_fluorine_c = False
    on_oxygen_lp = False
    on_oxygen_c = False 
    on_sulfur_lp = False
    on_sulfur_c = False 
    on_chlorine_lp = False
    on_chlorine_c = False
    on_bromine_lp = False
    on_bromine_c = False
    
    
    for atm in parent_node.item.atoms:
        if atm.atomtype[0].label not in ["O2s", "S2s", "O4tc", "Cs"]:
            if atm.is_carbon(): 
                if atm.lone_pairs!=[]:  # if the parent node does not specify lone pairs, then the child shouldn't have to
                    on_carbon_lp = True
                if atm.charge!=[]:  # if the parent node does not specify charge, then the child shouldn't have to
                    on_carbon_c = True
            if atm.is_fluorine(): 
                if atm.lone_pairs!=[]:  # if the parent node does not specify lone pairs, then the child shouldn't have to
                    on_fluorine_lp = True
                if atm.charge!=[]:  # if the parent node does not specify charge, then the child shouldn't have to
                    on_fluorine_c = True #but if the parent group is specific, make the child group specific
            if atm.is_oxygen(): 
                if atm.lone_pairs!=[]:  # if the parent node does not specify lone pairs, then the child shouldn't have to
                    on_oxygen_lp = True
                if atm.charge!=[]:  # if the parent node does not specify charge, then the child shouldn't have to
                    on_oxygen_c = True 
            if atm.is_sulfur():
                if atm.lone_pairs!=[]:  # if the parent node does not specify lone pairs, then the child shouldn't have to
                    on_sulfur_lp = True
                if atm.charge!=[]:  # if the parent node does not specify charge, then the child shouldn't have to
                    on_sulfur_c = True                 
            if atm.is_chlorine():
                if atm.lone_pairs!=[]:  # if the parent node does not specify lone pairs, then the child shouldn't have to
                    on_chlorine_lp = True
                if atm.charge!=[]:  # if the parent node does not specify charge, then the child shouldn't have to
                    on_chlorine_c = True  
            if atm.is_bromine():
                if atm.lone_pairs!=[]:  # if the parent node does not specify lone pairs, then the child shouldn't have to
                    on_bromine_lp = True
                if atm.charge!=[]:  # if the parent node does not specify charge, then the child shouldn't have to
                    on_bromine_c = True                     
                    
                    
    #save the flag values and pass them on
    flags = [on_carbon_lp, 
            on_carbon_c,
            on_fluorine_lp,
            on_fluorine_c,
            on_oxygen_lp,
            on_oxygen_c,
            on_sulfur_lp,
            on_sulfur_c,
            on_chlorine_lp, 
            on_chlorine_c,
            on_bromine_lp, 
            on_bromine_c] 

    # print(f'In make_groups, flags based on parent node are: {flags}')

    if atom.atomtype.label in ["O2s", "S2s", "O4tc", "Cs"]:
        group_atoms[atom] = GroupAtom(atomtype=[atom.atomtype],
                                     radical_electrons=[atom.radical_electrons],
                                     label='*')
    else:
        # print(atom.atomtype)
        if atom.is_carbon(): 
            # print('1.')
            #let's only worry about the first two flags
            flags_to_look_at = [flags[0], flags[1]]

        if atom.is_fluorine(): 
            # print('2.')
            #let's only worry about the last two flags
            flags_to_look_at = [flags[2], flags[3]]
        
        if atom.is_oxygen():
            # print('3.') 
            #let's only worry about the last two flags
            flags_to_look_at = [flags[4], flags[5]]

        #hydrogens should never have lone pairs or charge in a molecule
        if not atom.is_non_hydrogen():
            # print('4.')
            flags_to_look_at = [False, False]

        if atom.is_sulfur():
            # print('5.')
            flags_to_look_at = [flags[6], flags[7]]
            
        if atom.is_chlorine():
            # print('6.')
            flags_to_look_at = [flags[8], flags[9]]
        if atom.is_bromine():
            # print('7.')
            flags_to_look_at = [flags[10], flags[11]]
            

        if flags_to_look_at == [True, True]: #lp, charge both True
            # print('1')
            group_atoms[atom] = GroupAtom(atomtype=[atom.atomtype],
                                        radical_electrons=[atom.radical_electrons],
                                        lone_pairs=[atom.lone_pairs],
                                        charge=[atom.charge],
                                        label='*')
        if flags_to_look_at == [True, False]: #only lp True
            # print('2')
            group_atoms[atom] = GroupAtom(atomtype=[atom.atomtype],
                                        radical_electrons=[atom.radical_electrons],
                                        lone_pairs=[atom.lone_pairs],
                                        label='*')
        if flags_to_look_at == [False, True]: #only charge
            # print('3')
            group_atoms[atom] = GroupAtom(atomtype=[atom.atomtype],
                                        radical_electrons=[atom.radical_electrons],
                                        charge=[atom.charge],
                                        label='*')
        if flags_to_look_at == [False, False]: #lp, charge both False
            # print('4')
            group_atoms[atom] = GroupAtom(atomtype=[atom.atomtype],
                                        radical_electrons=[atom.radical_electrons],
                                        label='*')

    group = Group(atoms=list(group_atoms.values()))
    
    group = make_bonds(atom, group, group_atoms, n_degree_neighbor=n_degree_neighbor)

    group.update()

    return group    

def make_neighbor_name(atom,n_degree_neighbor=1,exclude=[]):
    """
    Make the name for n degree of neighbor recursively
    Args:
        atom (Atom): the center atom for the group
        n_degree_neighbor (Int): If a neighbor is n bond apart from the center atom, 
                                    it's defined as n degree neighbor.
        exclude (list of Atom): helper list to record atoms that have been named to avoid repetition 
        
    Returns:
        neighbors (string): name for n degrees of neighbors
        exclude (list of Atom): helper list to record atoms that have been named to avoid repetition 
        
    """
    
    if n_degree_neighbor == 1:

        names = []

        for atom2 in atom.edges.keys():
            
            if atom2.atomtype.label != 'H':
                
                if atom2 not in exclude:
                    atom_neighbor = atom2.atomtype.label
                    names.append(atom_neighbor)
                    exclude.append(atom2)
                
        neighbors = ''.join(sorted(names))
        neighbors += 'H' * len(['H' for atom2 in atom.edges.keys() if atom2.atomtype.label == 'H' and atom2 not in exclude])
        
        return neighbors, exclude

        
    else:
        
        names = []
        
        if atom not in exclude:
            exclude.append(atom)
        
        for atom2 in atom.edges.keys():
            
            if atom2.atomtype.label != 'H':
                
                atom_neighbor = ''
                
                if atom2 not in exclude:
                    
                    atom_neighbor += atom2.atomtype.label
                    exclude.append(atom2)
                    
                names.append(atom_neighbor)
                    
        for (i,atom2) in enumerate(atom.edges.keys()):
                
            atom2_neighbor, exclude = make_neighbor_name(atom2, n_degree_neighbor-1,exclude=exclude)

            if atom2_neighbor:
                names.append(f"({atom2_neighbor})")
                
        neighbors = ''.join(sorted(names))
        neighbors += 'H' * len(['H' for atom2 in atom.edges.keys() if atom2.atomtype.label == 'H' and atom2 not in exclude])
        
        for atom2 in atom.edges.keys():
            if atom2 not in exclude:
                exclude.append(atom2)
        
        return neighbors, exclude
    
def make_group_name(atom,n_degree_neighbor=1):
    """
    Make the group string for a group with atom as its center atom and n degree of neighbor, 
    using `make_neighbor_name`
    Args:
        atom (Atom): the center atom for the group
        n_degree_neighbor (Int): If a neighbor is n bond apart from the center atom, 
                                    it's defined as n degree neighbor.
        
    Returns:
        group_str (string): name for the group
    """
    
    
    exclude = []
    group_str = atom.atomtype.label
    
    neighbors = make_neighbor_name(atom,n_degree_neighbor=n_degree_neighbor,exclude=exclude)[0]
    
    if atom.radical_electrons==1:
        group_str += "J"

    if neighbors:
        group_str += f'-{neighbors}'
        
    return group_str