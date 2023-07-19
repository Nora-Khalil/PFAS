from os import uname
import cantera as ct
import numpy as np
import pandas as pd
import os
import csv

print("Running Cantera Version: " + str(ct.__version__))

To = 298
Po = 1e5 # ct.one_atm

directory = '/work/westgroup/nora/Code/projects/PFAS/models/testing_Siddhas_chemistry/CH3F/Siddhas_halogen_pdep_fluoromethane_switched/cantera/chem_annotated.yaml'

gas = ct.Solution(directory)


i = 0.025
mole_frac_list = list(np.linspace(0.025, 0.25, 50))


results = {}

for i in range(len(mole_frac_list)): 
    try: 
        x = mole_frac_list[i]
        string = f'****************************starting new volume fraction: {x} **************************'
        print(string)
        

        norm_ox = (1-x)*.21
        mole_frac_dict = {'CH3F(1)': (x/norm_ox), 'O2(2)':((1-x)*.21)/norm_ox, 'N2':((1-x)*0.79)/norm_ox } 
        gas.TPX = To, Po, mole_frac_dict
        width = 0.08
        flame = ct.FreeFlame(gas, width=width)
        flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1) 
        flame.max_time_step_count = 900
        loglevel = 1 
        try:
            if i!=0:
                d = f'data/{mole_frac_list[i-1]}_107.csv'
                if os.path.exists(d):  
                    arr2 = ct.SolutionArray(gas)
                    arr2.read_csv(d)
                    flame.set_initial_guess(data=arr2)
                    print(' initial guess has been set')
        except: 
            print('initial guess not set for this volume fraction')
            
        
        flame.solve(loglevel=loglevel, auto=False)
        Su = flame.velocity[0]
        results[x] = Su
        sltn = flame.to_solution_array()
        pd = sltn.to_pandas()
        pd.to_csv(f'data/{x}_107.csv', index=False)
        
    except Exception as e: 
        print(f'********************passed volume fraction:{mole_frac_list[i]}, error: {e}*************************************')
        pass

vol_fracs = list(results.keys())
flame_speeds = list(results.values())


print("volume fractions are:")
print(vol_fracs)

print("flame speeds are:")
print(flame_speeds)

with open('speeds_107.csv', 'w+') as g:
    g.write(directory)
    writers = csv.writer(g)
    writers.writerow(vol_fracs)
    writers.writerow(flame_speeds)
        


