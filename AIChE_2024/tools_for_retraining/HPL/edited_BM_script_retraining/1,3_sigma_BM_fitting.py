import pyximport
pyximport.install(reload_support=True)
from matplotlib.widgets import Button, Slider
import rmgpy
import numpy as np
from rmgpy.molecule.molecule import *
from rmgpy.species import *
from rmgpy.chemkin import *
from rmgpy.data.rmg import RMGDatabase
from IPython.display import display
from rmgpy.data.thermo import ThermoLibrary
from rmgpy.rmg.react import react
from rmgpy.species import Species
from rmgpy.reaction import Reaction
from rmgpy.data.rmg import get_db
from rmgpy.molecule.group import Group
from rmgpy.kinetics.arrhenius import ArrheniusBM
from rmgpy.data.kinetics.family import _make_rule
from rmgpy import settings
import time
import matplotlib.pyplot as plt
import matplotlib
import importlib
import sys
#matplotlib.use('qtagg')
matplotlib.use('TkAgg')
# %matplotlib inline

settings

thermo_libs = [
'C1_C2_Fluorine', #adding Siddha's as first most trusted because this is the thermo library that Franklin used
'PFCA_thermo',
'NCSU_C2_C8_PFAS', #adding Westmoreland's thermo as the second most trusted
'primaryThermoLibrary',
'Fluorine',
'FFCM1(-)',
'halogens',
'CHOF_G4',
'CHOCl_G4',
'CHOBr_G4',
'CHOFCl_G4',
'CHOFBr_G4',
'CHOFClBr_G4',
'DFT_QCI_thermo',
'2-BTP_G4',
'thermo_DFT_CCSDTF12_BAC',
'SulfurHaynes'
]

kin_families = ['1,3_sigmatropic_rearrangement']

database = RMGDatabase()
database.load(
            path = settings['database.directory'],
            thermo_libraries = thermo_libs,
            transport_libraries = [],
            reaction_libraries = [],
            seed_mechanisms = [],#['BurkeH2O2inN2','ERC-FoundationFuelv0.9'],
            kinetics_families = kin_families,
            kinetics_depositories = ['training'],
            #frequenciesLibraries = self.statmechLibraries,
            depository = False, # Don't bother loading the depository information, as we don't use it
        )

database.kinetics.families

# 'Lactone_Formation' 
family_to_train = "1,3_sigmatropic_rearrangement"
family = database.kinetics.families[family_to_train]

family.clean_tree()

start = time.time()
family.generate_tree(thermo_database=database.thermo,
                     nprocs=1,
                     new_fraction_threshold_to_reopt_node=0.25,
                     max_batch_size=800,
                     extension_iter_max=2,
                     extension_iter_item_cap=100)

end = time.time()
print(end-start)

len(family.groups.entries)

#let's see the group entries
family.groups.entries

#check the tree
start = time.time()
family.check_tree()
end = time.time()
print(end-start)

#regularize
start = time.time()
family.regularize(thermo_database=database.thermo)
end = time.time()
print(end-start)

#get reaction matches
start = time.time()
templateRxnMap = family.get_reaction_matches(thermo_database=database.thermo,remove_degeneracy=True,
                                             get_reverse=True,exact_matches_only=False,fix_labels=True)
end = time.time()
print(end-start)

len(templateRxnMap)

family.clean_tree_rules()

#make sure its clean
family.rules.entries

temperatures = np.array([ 300., 500.,  600.,  700. , 800.,  900., 1000., 1100., 1200., 1500., 2000.])
inverse_temps = [1000/T for T in temperatures]

def get_tr_rate(tr_rxn): 
    """Get training reaction rates over the range of temperatures"""
    log_tr_rates = []
    for temp in temperatures: 
        tr_rate = tr_rxn.get_rate_coefficient(temp)
        log_tr_rates.append(np.log(tr_rate))
    return log_tr_rates


def plot_node(node_label): 
    """
    Plot all of the rxns in the node and the fit for that node. 
    """
    #beginning of plotting 
    plt.ion()  # Enable interactive mode
    plt.figure()
    fig, ax = plt.subplots()

    
    #first get the BM rule of the node
    BM_rule = family.retrieve_original_entry(node_label)[0].data
    
    #second, let's get the training reactions of that node
    training_rxns = templateRxnMap[node_label]
    
    #calculate the rates for all of training reactions across all temps 
    for tr_rxn in training_rxns: 
        log_tr_rates = get_tr_rate(tr_rxn)
        ax.plot(inverse_temps, log_tr_rates, '-', c= 'black', label = str(tr_rxn))

        
    #define initial parameters 
    init_index_of_training_rxn = 0 #initialize with any of the training reaction indices
        
    def get_BM_rate_at_this_dH(index_of_training_rxn):
        log_BM_rates = []
        dH = training_rxns[index_of_training_rxn].get_enthalpy_of_reaction(298)
        for temp in temperatures: 
            BM_rate = BM_rule.get_rate_coefficient(temp, dH)
            log_BM_rates.append(np.log(BM_rate))
        return log_BM_rates #returns ydata
    
    #plot the BM fit 
    line = ax.plot(inverse_temps, get_BM_rate_at_this_dH(init_index_of_training_rxn), c='r', lw=2, label=str(training_rxns[init_index_of_training_rxn]))
    
    
    # adjust the main plot to make room for the sliders
    fig.subplots_adjust(bottom=0.25)
    
    # Make a horizontal slider to control the dH of the BM fitting.
    axTR_Ind = fig.add_axes([0.25, 0.1, 0.65, 0.03])
    dH_slider = Slider(
        ax=axTR_Ind,
        label='dH',
        valmin=0,
        valmax=len(training_rxns),
        valinit=init_index_of_training_rxn,
        valstep=1,  # Ensures only integer steps are allowed
    )
    

    
    # The function to be called anytime a slider's value changes
    def update(val):
        line.set_ydata(get_BM_rate_at_this_dH(index_of_training_rxn.val))
        #highlight the training reaction that its getting the dH from
        log_tr_rates_to_highlight = get_tr_rate(training_rxns[index_of_training_rxn.val])
        ax.plot(inverse_temps, log_tr_rates_to_highlight, ':', c= 'r', label = str(training_rxns[index_of_training_rxn.val]))
        fig.canvas.draw_idle()
        
    # register the update function with each slider
    dH_slider.on_changed(update)
    
    
    # Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
    resetax = fig.add_axes([0.8, 0.025, 0.1, 0.04])
    button = Button(resetax, 'Reset', hovercolor='0.975')
    
    
    def reset(event):
        dH_slider.reset()
    button.on_clicked(reset)
    
    ax.set_xlabel('1000/T')
    ax.set_ylabel('ln(k)')
    plt.title(node_label)
    plt.legend()
    plt.show()
    
    
#doing one node at a time 
#testing
entries = list(family.groups.entries.values())
Tref=1000.0
fmax=1.0e5
completed_nodes = 0 


rule_keys = family.rules.entries.keys()
for entry in family.groups.entries.values():
    if entry.label not in rule_keys:
        family.rules.entries[entry.label] = []


for node_name, rxns_list in templateRxnMap.items():
    
    #print(f'STARTING NODE {node_name}')
    #make a subset dictionary
    subset_dict = {}
    subset_dict[node_name] = rxns_list
    # family.make_bm_rules_from_template_rxn_map(subset_dict)
    
    
    rxnlists = [(subset_dict[entry.label], entry.label) for entry in entries if entry.label in subset_dict.keys()]
    # rxnlists = [(templateRxnMap[entry.label], entry.label) for entry in entries if entry.label in templateRxnMap.keys()]

    inputs = np.array([(family.forward_recipe.actions, rxns, Tref, fmax, label, [r.rank for r in rxns]) for rxns, label in rxnlists])
    inds = np.arange(len(inputs))
    np.random.shuffle(inds)  # want to parallelize in random order
    inds = inds.tolist()
    revinds = [inds.index(x) for x in np.arange(len(inputs))]

    #Modifying slightly so its not in random order
    kinetics_list = []
    for ind in inds: 
        node_kinetics = _make_rule(inputs[ind])
        kinetics_list.append(node_kinetics)
    kinetics_list = np.array(kinetics_list)
    kinetics_list = kinetics_list[revinds]  # fix order

    for i, kinetics in enumerate(kinetics_list): 
        
        for entry_ in entries: 
            if entry_.label==node_name:
                entry = entry_
        new_entry = Entry(
                    index=entry.index,
                    label=entry.label,
                    item=family.forward_template,
                    data=kinetics,
                )
        family.rules.entries[entry.label].append(new_entry)

    print(f'FINISHED NODE {node_name}')
    plot_node(node_name)
    completed_nodes+=1
    
print(f'{completed_nodes} were successfully fit to ArrheniusBM')