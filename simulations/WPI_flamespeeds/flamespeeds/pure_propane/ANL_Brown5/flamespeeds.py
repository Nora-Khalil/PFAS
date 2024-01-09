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



######  skip phis that are already calculated 

csv_files = [file for file in os.listdir('.') if re.match('([0-9]\.[0-9]+)\.csv', file)]
phis_already_calculated = [float(re.match('([0-9]\.[0-9]+)\.csv', file).group(1)) for file in csv_files]
#########

mole_fractions = list(np.linspace(0.1, 0.4, 40))   #mole fraction of propane
results = {}

for pos_in_list in range(0, len(mole_fractions)):

    x = mole_fractions[pos_in_list]

    phi=x/.2

    if phi in phis_already_calculated: 
        print(f'this phi was already calculated: {phi}')
        continue

    try: 
        mole_frac_dict = {'C3H8': x, 'O2':1, 'N2':3.76} 

        string = f'****************************starting phi: {phi}  **************************'

        gas.TPX = To, Po, mole_frac_dict

        #width = 0.08 #changing width, have to "make sure that the domain is significantly wider than the flame thickness" https://groups.google.com/g/cantera-users/c/G9eEfVkeiM0
        width = 0.12
        flame = ct.FreeFlame(gas, width=width)
        #Set maxiumum number of grid points to be very high (otherwise default is 1000)
        flame.set_max_grid_points(flame.domains[flame.domain_index("flame")], 1e4)
        #flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1)
        flame.set_refine_criteria(ratio=5, slope=0.2, curve=0.3) #using middle criteria
        #flame.set_refine_criteria(ratio=10, slope=0.8, curve=0.8) 
        flame.max_time_step_count = 2600
        loglevel = 1 


        try:
            if pos_in_list!=0:
                previous_phi = (mole_fractions[pos_in_list-1])/.2
                d = f'{previous_phi}_edited_names.csv'
                print(f'Looking for a data file that matches {d}')
                if os.path.exists(d):  
                    arr2 = ct.SolutionArray(gas)
                    arr2.read_csv(d)
                    flame.set_initial_guess(data=arr2)
                    print('initial guess has been set')
                else: 
                    print(f'Could not find data file of {d} to set initial guess')
                    print('initial guess not set for this phi')
        except: 
            print('There was an error. Initial guess not set for this phi')


        flame.solve(loglevel=loglevel, auto=False) #turning off auto=False
        Su = flame.velocity[0]
        results[phi] = Su
        print(f' at {phi} phi, speed is {Su}')
        sltn = flame.to_solution_array()
        df1 = sltn.to_pandas()
        df1.to_csv(f'{phi}.csv', index=False)  

        #make sure you also rewrite the solution array so it can be used as a guess in the future
        rewritten = []
        current_csv_file = f'{phi}.csv'
        with open(current_csv_file, 'r') as f: 
            lines = f.readlines() 
            for line in lines: 
                new_line = line.replace('Y_', '')
                rewritten.append(new_line)
        new_file_name = current_csv_file.replace('.csv', '_edited_names.csv')
        with open(new_file_name, 'w') as f: 
            f.writelines(rewritten)

    except Exception as e: 
        print(f'********************passed phi:{phi}, error: {e}*************************************')
        pass

    print("Phi is:")
    print(phi)

    print("flame speed is:")
    print(Su)
     
