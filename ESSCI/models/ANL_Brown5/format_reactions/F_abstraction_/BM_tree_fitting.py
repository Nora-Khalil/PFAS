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
from rmgpy import settings
import time
import matplotlib.pyplot as plt
import matplotlib

print(settings)

database = RMGDatabase()
database.load(
            path = settings['database.directory'],
            thermo_libraries = ['Klippenstein_Glarborg2016', 'BurkeH2O2', 'thermo_DFT_CCSDTF12_BAC', 
                               'DFT_QCI_thermo',
                           'primaryThermoLibrary', 'primaryNS', 'NitrogenCurran', 'NOx2018', 'FFCM1(-)',
'SulfurLibrary', 'SulfurGlarborgH2S','SABIC_aromatics'],
            transport_libraries = [],
            reaction_libraries = [],
            seed_mechanisms = [],#['BurkeH2O2inN2','ERC-FoundationFuelv0.9'],
            kinetics_families = 'F_Abstraction',
            kinetics_depositories = ['training'],
            #frequenciesLibraries = self.statmechLibraries,
            depository = False, # Don't bother loading the depository information, as we don't use it
        )

family = database.kinetics.families["F_Abstraction"]

family.clean_tree()

start = time.time()
family.generate_tree(thermo_database=database.thermo,
                     nprocs=1,
                     new_fraction_threshold_to_reopt_node=0.25,
                     max_batch_size=800,
                     extension_iter_max=2,
                     extension_iter_item_cap=300)
end = time.time()
print(f'time to generate tree {end-start}')

print(len(family.groups.entries))
print(family.groups.entries)


start = time.time()
family.check_tree()
end = time.time()
print(f'time to check the tree {end-start}')

start = time.time()
family.regularize(thermo_database=database.thermo)
end = time.time()
print(f'time to regularize the tree {end-start}')

start = time.time()
templateRxnMap = family.get_reaction_matches(thermo_database=database.thermo,remove_degeneracy=True,
                                             get_reverse=True,exact_matches_only=False,fix_labels=True)
end = time.time()
print(f'time to make template map {end-start}')

print(len(templateRxnMap))

family.clean_tree_rules()

start = time.time()
family.make_bm_rules_from_template_rxn_map(templateRxnMap)#,nprocs=6)
end = time.time()
print(f'time to make the bm rules {end-start}')

start = time.time()
family.check_tree()
end = time.time()
print(f'check tree again {end-start}')

save_path = os.path.join(settings['database.directory'], 'kinetics', 'families', family.name)
print(save_path)

family.save(save_path)