'''
Papas, Zhang, p. 1150: "The stoichiometry for the hydrofluorocarbon-air systems considered was determined by taking the combustion products to be CO2, HF, and H2O. If there was insufficient hydrogen available for formation of HF and H2O, then the formation of HF took preference over H2O formation. If there was insufficient hydrogen available for all the fluorine to form HF, then the remaining fluorine was assumed to produce CF2O in preference of carbon forming CO2. "

Assumed stoichiometric combustion for CH2F2 is: 


CH2F2 + (O2 + 3.76 N2) = CO2 + 2HF + 0 H2O + 3.76 N 2

'''

import cantera as ct
import numpy as np
import pandas as pd
import csv
import sys
import re

print("Running Cantera Version: " + str(ct.__version__ ))

To = 298
Po = ct.one_atm


path_to_mech = '/work/westgroup/nora/Code/projects/PFAS/simulations/WPI_flamespeeds/models/ANL_brown5.yaml'


gas_full = ct.Solution(path_to_mech)


#####
# use this toggle to use either the full mechanism (True) or a smaller version (False)
FullMech = False
#FullMech = True

if FullMech==False:
    
    #exclude = ['CH2CCHCH2CHCH2'] #add specific species here. (for illustration purposes, since the Max_N_oxygen below would eliminate it anyway)
    max_N_carbon = 4
    max_N_oxygen = 3
    
    all_species = ct.Species.listFromFile(path_to_mech)
    species = []
    # Filter species
    for S in all_species:
#         if S.name in exclude:
#             print( "excluding %s"%(S.name) )
#             continue #skip this species
        comp = S.composition
        if S.name in ['C6H9', 'C6H10']:
            continue
        if 'C' in comp and comp['C']> max_N_carbon:
            print( "excluding %s"%(S.name) )
            continue
        if 'O' in comp and comp['O']> max_N_oxygen:
            print( "excluding %s"%(S.name) )
            continue            
#         if 'He' in comp: # Exclude Helium
#             print( "excluding %s"%(S.name) )
#             continue
#         if 'Ar' in comp: # Exclude Argon
#             print( "excluding %s"%(S.name) )
#             continue
#         if 'Kr' in comp: # Exclude Krypton
#             print( "excluding %s"%(S.name) )
#             continue   
        #print S.name    
        species.append(S)

    species_names = {S.name for S in species}
    # Filter reactions, keeping only those that only involve the selected species
    all_reactions = gas_full.reactions()
    reactions = []

    for R in all_reactions:
        if not all(reactant in species_names for reactant in R.reactants):
            continue
        if not all(product in species_names for product in R.products):
            continue
        reactions.append(R)

    gas_small = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                       species=species, reactions=reactions)
    gas = gas_small    
else:
    gas = gas_full
    
print( "final mechanism has %d species and %d reactions"%(gas.n_species, gas.n_reactions) )
print(path_to_mech)



mole_fractions = list(np.linspace(0.4, 1.6, 30))   #mole fraction of CH2F2
results = {}

try: 

    x = mole_fractions[int(eval(sys.argv[-1]))]

    mole_frac_dict = {'CH2F2': x, 'O2':1, 'N2':3.76} 

    # phi calculated as follows: phi = (F/A)_ac / (F/A)_stoic, where (F/A)_stoic assumed to be one, and A in (F/A)_ac is 1.
    # F is x, the mole fraction of CH2F2

    phi = x 
    string = f'****************************starting phi: {phi}  **************************'

    gas.TPX = To, Po, mole_frac_dict

    width = 0.12
    flame = ct.FreeFlame(gas, width=width)
    flame.set_max_grid_points(flame.domains[flame.domain_index("flame")], 1e4)
    flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1) 
    flame.max_time_step_count = 2600
    loglevel = 1 


#         try:
#             phi_before = mole_fractions[i-1]
#             if i!=0:
#                 d = f'{phi_before}_difluoromethane.csv'
#                 if os.path.exists(d):  
#                     arr2 = ct.SolutionArray(gas)
#                     arr2.read_csv(d)
#                     flame.set_initial_guess(data=arr2)
#                     print('initial guess has been set')
#         except: 
#             print('initial guess not set for this volume fraction')
        

    flame.solve(loglevel=loglevel, auto=False)
    Su = flame.velocity[0]
    results[x] = Su
    sltn = flame.to_solution_array()
    df1 = sltn.to_pandas()
    df1.to_csv(f'{phi}.csv', index=False)            

except Exception as e: 
    print(f'********************passed phi:{phi}, error: {e}*************************************')
    pass


phis = list(results.keys())
flame_speeds = list(results.values())


print("Phi is:")
print(phi)

print("flame speed is:")
print(Su)

        