'''
Papas, Zhang, p. 1150: "The stoichiometry for the hydrofluorocarbon-air systems considered was determined by taking the combustion products to be CO2, HF, and H2O. If there was insufficient hydrogen available for formation of HF and H2O, then the formation of HF took preference over H2O formation. If there was insufficient hydrogen available for all the fluorine to form HF, then the remaining fluorine was assumed to produce CF2O in preference of carbon forming CO2. "

Assumed stoichiometric combustion for propane is: 


C3H8 + 5(O2 + 3.76 N2) = 3CO2 + 4 H2O + 3.76*5 N2



calculated phi as: 

at stoichiometric: propane: 0.2, O2: 1, N2: 3.76

therefore phi = (F/M ac) / (F/M stoich) = (x/(x+1+3.76))/(0.2/(0.2+1+3.76)) = (24.8)*(x/(x+1+3.76)) #WRONG

phi = x/.2

'''
import cantera as ct
import numpy as np
import pandas as pd
import csv
import sys
import os
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


mole_fractions_list = []
percents_of_CH2F2 = [0.2, 0.3, 0.5, 0.7, 1, 2, 3, 4, 5, 6]
for x in percents_of_CH2F2: 
    # if x<=2: 
    y = (x/100)*4.96 # this will give you y mole fraction of CH2F2 in the mixture
    fractions = {'CH2F2': y, 'C3H8': .2, 'O2': 1, 'N2': 3.76-y}
    mole_fractions_list.append(fractions)


results = {}


try: 


    mole_frac_dict = mole_fractions_list[int(eval(sys.argv[-1]))]
    print(mole_frac_dict)
    

    vf = percents_of_CH2F2[int(eval(sys.argv[-1]))]


    string = f'****************************starting vf: {vf}  **************************'

    gas.TPX = To, Po, mole_frac_dict

    width = 0.12
    flame = ct.FreeFlame(gas, width=width)
    flame.set_max_grid_points(flame.domains[flame.domain_index("flame")], 1e4)
    flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1) 
    flame.max_time_step_count = 2600
    loglevel = 1 
      
    flame.solve(loglevel=loglevel, auto=False)
    Su = flame.velocity[0]
    results[phi] = Su
    sltn = flame.to_solution_array()
    df1 = sltn.to_pandas()
    df1.to_csv(f'{vf}.csv', index=False)            
            
except Exception as e: 
    print(f'********************passed vf:{vf}, error: {e}*************************************')
    pass

print("vf is:")
print(vf)

print("flame speeds is:")
print(Su)

    
        