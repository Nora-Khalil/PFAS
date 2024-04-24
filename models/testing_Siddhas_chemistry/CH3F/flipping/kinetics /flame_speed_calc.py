from os import uname
import cantera as ct
import numpy as np
import pandas as pd
import os
import csv
import sys

print("Running Cantera Version: " + str(ct.__version__))

To = 298
Po = 1e5 # ct.one_atm

file_name=sys.argv[-1]

directory = f'/work/westgroup/nora/Code/projects/PFAS/models/testing_Siddhas_chemistry/CH3F/flipping/flipped_models/copies/{file_name}'

gas = ct.Solution(directory)



mole_frac_list = [0.125]


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
                d = f'data/copy_1006.cti.csv'
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
        pd.to_csv(f'data/{file_name}.csv', index=False)
        
    except Exception as e: 
        print(f'********************passed volume fraction:{mole_frac_list[i]}, error: {e}*************************************')
        pass

vol_fracs = list(results.keys())
flame_speeds = list(results.values())


print("volume fractions are:")
print(vol_fracs)

print("flame speeds are:")
print(flame_speeds)


        


