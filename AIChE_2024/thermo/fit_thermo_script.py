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
from helperfunctions import get_neighbors, make_bonds, make_group, make_neighbor_name, make_group_name

libraries = ["C1_C2_Fluorine",
            "C1_C3_hydrofluorocarbons_NIST",  #obtained via email by Linteris on 3/10/23
             "NCSU_C2_C8_PFAS", 
             "PFCA_thermo", 
             
             #used for training by David
             "CHOF_G4",
             "CHOFCl_G4",
             "CHOFClBr_G4",
             "CHOFBr_G4",
            ]

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


special_list = [
    'O0sc-N5dc',
    'S2d-CS',
    'S2d-C2d',
    'O2d-N5dc',
    'O2d-CO',
    'O2d-S6dd',
    'O2d-N3d',
    'O2d-S4dd',
    'O2d-S4d',
    'O2d-Cdd',
    'Cb-CSCbCb',
    'Cdd-CdO2d',
]

database = ThermoDatabase()
database.load(os.path.join(settings['database.directory'],"thermo"),
             libraries = libraries,
             depository=False)



#for storing all of the data
all_data = dict()
all_data["entry"] = list()
all_data["entry thermo"] = list()
all_data["GAE without estimated GAV"] = list()
all_data["GAE with estimated GAV"] = list()
all_data["missing GAE"] = list()
all_data["missing groups"] = list()

missing_group_index_dict = dict()
missing_group_dict = dict()
missing_group_index = 0

good_estimates_with_matched_node = []
missing_group_generated = []
matched_node_with_bad_estimate = []

nodes_that_matched = {}
for library in database.libraries:
    print(f'Starting to look at library: {library}')
    entries=list(database.libraries[library].entries.items())
    
    for item, entry in entries:
        
        if entry.data is not None:
            
            if not isinstance(entry.data, ThermoData):
                try:
                    entry_thermo = entry.data.to_thermo_data()
                except:
                    continue
            else:
                entry_thermo = entry.data

            molecule = entry.item
            print("\n========================================================")
            print(f"Current species is {entry} ({molecule.get_formula()})\n")
            pass_this_molecule_through_as_missing_group=False
            
            if molecule.smiles in ["[Ar]","[He]","[Ne]"]:
                #Current group additivity tree only contains C, N, S, O, and some halogen chemistry
                #Use `primaryThermoLibrary` for noble gas thermo
                continue
            if molecule.is_radical():
                print('GAV not built for radicals. Passing this species.')
            if not molecule.is_radical():
                #Current notebook only works on non-radical species.
                #Should be able to extend for radical groups with careful adaptation
                
                estimated_thermo = database.estimate_thermo_via_group_additivity(molecule)
                
                if (estimated_thermo.H298.value_si-entry_thermo.H298.value_si)/4180 < 2:
                    print('Estimated thermo via GAV is already good.')                        
                    molecule.sort_atoms()
                    for atom in molecule.atoms:
                        if atom.is_non_hydrogen() and not atom.is_halogen() and not atom.is_nitrogen():
                            node0 = database.groups['group'].descend_tree(molecule, {'*': atom}, None) #this is the node that it matches to. 
                            data = node0.data
                            print(f"The matched node used to estimate the GAV was {node0}.")
                            if node0.long_desc!= '': 
                                print(f"Long desc attached with this node:\n {node0.long_desc}")
                            if node0.long_desc=='': 
                                print('No long desc attached to this node.')
                            #these can be retrained, but the estimates were good so it won't do much
                            if estimated_thermo.H298.value_si!=0.0:
                                good_estimates_with_matched_node.append(node0)
                            if estimated_thermo.H298.value_si==0.0:
                                pass_this_molecule_through_as_missing_group=True
                        else: 
                            print('This is a hydrogen or halogen or nitrogen')

                if (estimated_thermo.H298.value_si-entry_thermo.H298.value_si)/4180 > 2 or pass_this_molecule_through_as_missing_group==True:                    
                    print('Estimated thermo is off.')
                    
                    missing_grp = list()
                    missing = 0 #let's count the missing groups

                    real_data_thermo = ThermoData(
                            Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
                            Cpdata=([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], "J/(mol*K)"),
                            H298=(0.0, "kJ/mol"),
                            S298=(0.0, "J/(mol*K)"),
                        )

                    molecule.sort_atoms()
                    
                    for atom in molecule.atoms:
                            
                        #here we loop through each atom. If atom is hydrogen or halogen, we don't descend
                        #it down tree because we don't consider hydrogen or halogen as center atom.
                        #If the node data is None, string, or all zeros, we go in the branch of generating
                        #missing group structure. If the node data is just normal data, we add it to
                        #`real_data_thermo`.
    
                        if atom.is_non_hydrogen() and not atom.is_halogen():
                            #Hydrogen and halogen are not considered as center atom

                            node0 = database.groups['group'].descend_tree(molecule, {'*': atom}, None) #this is the best matched node
                            node = node0
                            data = node.data

                            add_to_real_data_thermo = True
                            print(f"The best matched node is: {node}")
                            
                            ###############################################################################
                                #Start identifying missing group
                            if data is None or isinstance(data,str) or data.is_all_zeros():
                                print('There is no data in this best matched group. \nThis is probably a parent node with no specific child.')
                                
                                #make group structure and group string
                                n_degree_neighbor = 1
                                #let's pass in the parent node to make_group
                                parent_node = node0
                                group = make_group(atom, parent_node, n_degree_neighbor=n_degree_neighbor)
                                group_str = make_group_name(atom, n_degree_neighbor=n_degree_neighbor)
                                group_0 = group
                                group_str_0 = group_str
                                print(f"First initial generated group is: {group_str_0}")
                                print(group_0.to_adjacency_list())
                                
                                if group_str.split("_")[-1] not in special_list:
                                    
                                    while not group.is_subgraph_isomorphic(node.item):
                                    #while not group.is_subgraph_isomorphic(node.item, generate_initial_map = True):#Su suggested this
                                        #new group has to be the child-node of the originally matched group
                                        #to have correct parent-node child-node relation
                                        print(f"newly made group \n{group.to_adjacency_list()} is not subgraph isomorphic to the parent node {node} \n{node.item.to_adjacency_list()}")
                                        print("Increasing n_degree_neighbor and trying again.")
                                        n_degree_neighbor+=1
                                        group = make_group(atom, parent_node, n_degree_neighbor=n_degree_neighbor)
                                        group_str = make_group_name(atom, n_degree_neighbor=n_degree_neighbor)
                                    print('Newly made child group passed check for being sub-iso to parent node')
                                    while any([(group.make_sample_molecule()).is_subgraph_isomorphic(child.item, generate_initial_map = True) for child in node.children]):
                                        #Child-node can't be the child of other children
                                        #to avoid ambiguous group selection
                                        n_degree_neighbor+=1
                                        group = make_group(atom, parent_node, n_degree_neighbor=n_degree_neighbor)
                                        group_str = make_group_name(atom, n_degree_neighbor=n_degree_neighbor)
                                    print('Newly made child group passed check to not being sub-iso to other children of the parent node.')
                                if group_str.split("_")[-1] not in special_list:
                                    add_to_real_data_thermo = False 

                                    group.sort_atoms()
                                                  
                                                  #parent label  #new label
                                    group_str = f'{node.label}_{group_str}'

                                    missing += 1 #increase count for missing groups for this molecule
                                    missing_grp.append(group_str) 
                                    missing_group_generated.append(group_str)
                                    if group_str not in missing_group_index_dict:
                                        missing_group_index_dict[group_str] = missing_group_index
                                        missing_group_index+=1

                                        missing_group_dict[group_str] = dict()
                                        missing_group_dict[group_str]["group"] = [group]
                                        missing_group_dict[group_str]["atom"] = [atom]
                                        missing_group_dict[group_str]["molecule"] = [molecule]
                                        missing_group_dict[group_str]["label"] = [entry.label]
                                    else:
                                        missing_group_dict[group_str]["group"].append(group)
                                        missing_group_dict[group_str]["atom"].append(atom)
                                        missing_group_dict[group_str]["molecule"].append(molecule)
                                        missing_group_dict[group_str]["label"].append(entry.label)

                ###############################################################################
                #calculate real data thermo value
                            if add_to_real_data_thermo: #this only happens if it is not a missing group.
                                while node is not None and node.data is None:
                                    node = node.parent #if there is no data in this node, go up to the parent with data. 
                                if node is None:
                                    raise DatabaseError(f'Unable to determine thermo parameters for atom {atom} in molecule {molecule}: '
                                                        f'no data for node {node0} or any of its ancestors in database {database.label}.')

                                data = node.data
                                comment = node.label
                                loop_count = 0
                                while isinstance(data, str):
                                    loop_count += 1
                                    if loop_count > 100:
                                        raise DatabaseError("Maximum iterations reached while following thermo group data pointers. A circular"
                                                            f" reference may exist. Last node was {node.label} pointing to group called {data} in "
                                                            f"database {database.label}")

                                    for entr in database.groups["group"].entries.values():
                                        if entr.label == data:
                                            data = entr.data
                                            comment = entr.label
                                            break
                                    else:
                                        raise DatabaseError(f"Node {node.label} points to a non-existing group called {data} "
                                                            f"in database {database.label}")

                                data.comment = '{0}({1})'.format(database.groups['group'].label, comment)
                                print('Adding thermo data ') 
                                matched_node_with_bad_estimate.append(node)
                                add_thermo_data(real_data_thermo, data, group_additivity=True)

                    cyclic = molecule.is_cyclic()

                    if cyclic:
                        sssr = molecule.get_smallest_set_of_smallest_rings()
                        for ring in sssr:
                            for atomPair in itertools.permutations(ring, 2):
                                try:
                                    database._add_group_thermo_data(real_data_thermo, database.groups['longDistanceInteraction_cyclic'], molecule,
                                                                {'*1': atomPair[0], '*2': atomPair[1]})
                                except KeyError:
                                    pass

                    # Do ring corrections separately because we only want to match
                    # each ring one time

                    if cyclic:
                        monorings, polyrings = molecule.get_disparate_cycles()
                        for ring in monorings:
                            # Make a temporary structure containing only the atoms in the ring
                            # NB. if any of the ring corrections depend on ligands not in the ring, they will not be found!
                            try:
                                database._add_ring_correction_thermo_data_from_tree(real_data_thermo, database.groups['ring'], molecule, ring)
                            except KeyError:
                                logging.error("Couldn't find a match in the monocyclic ring database even though "
                                              "monocyclic rings were found.")
                                logging.error(molecule)
                                logging.error(molecule.to_adjacency_list())
                                raise
                        for polyring in polyrings:
                            # Make a temporary structure containing only the atoms in the ring
                            # NB. if any of the ring corrections depend on ligands not in the ring, they will not be found!
                            try:
                                database._add_polycyclic_correction_thermo_data(real_data_thermo, molecule, polyring)
                            except KeyError:
                                logging.error("Couldn't find a match in the polycyclic ring database even though "
                                              "polycyclic rings were found.")
                                logging.error(molecule)
                                logging.error(molecule.to_adjacency_list())

    ########################################################################################################
                    if missing > 0:
            #If there are missing groups identified in this molecule, add to all_data (dict)
                        entry.short_desc = library
                        all_data["entry"].append(entry)
                        all_data["GAE without estimated GAV"].append(real_data_thermo)
                        all_data["GAE with estimated GAV"].append(estimated_thermo)
                        all_data["missing groups"].append(missing_grp)
                        all_data["entry thermo"].append(entry_thermo)

                        try:
                            #we remove the contribution from `real_data_thermo` in the old thermo estimation
                            #the rest of thermo are contributed by the missing groups
                            missing_group_thermo = remove_thermo_data(deepcopy(entry_thermo),real_data_thermo)
                            all_data["missing GAE"].append(missing_group_thermo)
                        except (ValueError,IndexError):
                            #We need Cp0 and CpInf to perform the inversion from Nasa to ThermoData
                            if entry_thermo.Cp0 is None:
                                cp_0 = molecule.calculate_cp0()
                                entry_thermo.Cp0 = (cp_0, "J/(mol*K)")
                            if entry_thermo.CpInf is None:
                                cp_inf = molecule.calculate_cpinf()
                                entry_thermo.CpInf = (cp_inf, "J/(mol*K)")

                            nasa = entry_thermo.to_nasa(Tmin=10.0, Tmax=3000.0, Tint=500.0)
                            entry_thermo = nasa.to_thermo_data()
                            missing_group_thermo = remove_thermo_data(deepcopy(entry_thermo),real_data_thermo)
                            all_data["missing GAE"].append(missing_group_thermo)
                    if missing == 0:
                        #if there are no missing groups, let's isolate these species and do something with them later
                        print('no groups missing')
        # print('completed')
spc_num = len(all_data["entry"])
grp_num = len(missing_group_index_dict.keys())
print(f"Fitting {grp_num} of new groups with {spc_num} of species")


