#########calculates flamespeeds for the .cti files at different equivalence ratios. Uses initial guess from the previous model #############################


#####will be used to calculate the burning velocity of ******stoichiometric******* 2,3,3,3-tetrafluoropropene-difluoromethane mix increasing mole fraction of tetra 

import cantera as ct
import numpy as np
import pandas as pd
import os
import csv 
import sys
import re


print("Running Cantera Version: " + str(ct.__version__))

file_name='220'

To = 298
Po = ct.one_atm
directory ='/work/westgroup/nora/Code/projects/PFAS/ESSCI/models/RMG/after_ESSCI_second_try/4FP_CH2F2/species_in_core/chemkin/copies/copy_chem0220.cti'


gas = ct.Solution(directory)

#### has tuples containing (tetra mole fraction, stoichiometric dictionary at that mole fraction)
# master_dictionary = [(0.0, {'CH2F2(39)': 1.0, 'C3H2F4(1)': 0.0, 'O2(3)': 1.0, 'N2': 3.76}),
#  (0.041666666666666664,
#   {'CH2F2(39)': 0.9019607843137255,
#    'C3H2F4(1)': 0.0392156862745098,
#    'O2(3)': 1.0,
#    'N2': 3.76}),
#  (0.08333333333333333,
#   {'CH2F2(39)': 0.8148148148148148,
#    'C3H2F4(1)': 0.07407407407407407,
#    'O2(3)': 1.0,
#    'N2': 3.76}),
#  (0.125,
#   {'CH2F2(39)': 0.7368421052631579,
#    'C3H2F4(1)': 0.10526315789473684,
#    'O2(3)': 1.0,
#    'N2': 3.76}),
#  (0.16666666666666666,
#   {'CH2F2(39)': 0.6666666666666667,
#    'C3H2F4(1)': 0.13333333333333333,
#    'O2(3)': 1.0,
#    'N2': 3.7599999999999993}),
#  (0.20833333333333331,
#   {'CH2F2(39)': 0.6031746031746033,
#    'C3H2F4(1)': 0.15873015873015872,
#    'O2(3)': 1.0,
#    'N2': 3.76}),
#  (0.25,
#   {'CH2F2(39)': 0.5454545454545454,
#    'C3H2F4(1)': 0.18181818181818182,
#    'O2(3)': 1.0,
#    'N2': 3.76})]
                     
                     
                     
#  (0.29166666666666663,
#   {'CH2F2(39)': 0.4927536231884058,
#    'C3H2F4(1)': 0.20289855072463767,
#    'O2(3)': 1.0,
#    'N2': 3.7599999999999993}),
#  (0.3333333333333333,
#   {'CH2F2(39)': 0.4444444444444445,
#    'C3H2F4(1)': 0.2222222222222222,
#    'O2(3)': 1.0,
#    'N2': 3.76}),
#  (0.375, {'CH2F2(39)': 0.4, 'C3H2F4(1)': 0.24, 'O2(3)': 1.0, 'N2': 3.76}),
#  (0.41666666666666663,
#   {'CH2F2(39)': 0.358974358974359,
#    'C3H2F4(1)': 0.25641025641025633,
#    'O2(3)': 1.0,
#    'N2': 3.76}),
#  (0.4583333333333333,
#   {'CH2F2(39)': 0.3209876543209877,
#    'C3H2F4(1)': 0.2716049382716049,
#    'O2(3)': 1.0,
#    'N2': 3.76}),
#  (0.5,
#   {'CH2F2(39)': 0.2857142857142857,
#    'C3H2F4(1)': 0.2857142857142857,
#    'O2(3)': 1.0,
#    'N2': 3.7600000000000002}),


# master_dictionary = [(0.5416666666666666,
#   {'CH2F2(39)': 0.2528735632183908,
#    'C3H2F4(1)': 0.29885057471264365,
#    'O2(3)': 1.0,
#    'N2': 3.76}),
#  (0.5833333333333333,
#   {'CH2F2(39)': 0.22222222222222227,
#    'C3H2F4(1)': 0.31111111111111106,
#    'O2(3)': 1.0,
#    'N2': 3.76}),
#  (0.625,
#   {'CH2F2(39)': 0.1935483870967742,
#    'C3H2F4(1)': 0.3225806451612903,
#    'O2(3)': 1.0,
#    'N2': 3.76}),
#  (0.6666666666666666,
#   {'CH2F2(39)': 0.16666666666666669,
#    'C3H2F4(1)': 0.3333333333333333,
#    'O2(3)': 1.0,
#    'N2': 3.76}),
#  (0.7083333333333333,
#   {'CH2F2(39)': 0.1414141414141414,
#    'C3H2F4(1)': 0.3434343434343433,
#    'O2(3)': 1.0,
#    'N2': 3.76}),
#  (0.75,
#   {'CH2F2(39)': 0.11764705882352941,
#    'C3H2F4(1)': 0.35294117647058826,
#    'O2(3)': 1.0,
#    'N2': 3.76})]
master_dictionary = [(0.7916666666666666,
  {'CH2F2(39)': 0.09523809523809526,
   'C3H2F4(1)': 0.3619047619047619,
   'O2(3)': 1.0,
   'N2': 3.76}),
 (0.8333333333333333,
  {'CH2F2(39)': 0.0740740740740741,
   'C3H2F4(1)': 0.37037037037037024,
   'O2(3)': 1.0,
   'N2': 3.76}),
  (0.875,
  {'CH2F2(39)': 0.05405405405405406,
   'C3H2F4(1)': 0.3783783783783784,
   'O2(3)': 1.0,
   'N2': 3.7600000000000002}),
 (0.9166666666666666,
  {'CH2F2(39)': 0.03508771929824563,
   'C3H2F4(1)': 0.38596491228070173,
   'O2(3)': 1.0,
   'N2': 3.76}),
 (0.9583333333333333,
  {'CH2F2(39)': 0.01709401709401712,
   'C3H2F4(1)': 0.39316239316239304,
   'O2(3)': 1.0,
   'N2': 3.76}),
 (1.0,
  {'CH2F2(39)': 0.0,
   'C3H2F4(1)': 0.4,
   'O2(3)': 1.0,
   'N2': 3.7599999999999993})]

results = {}

for fraction, dictionary in master_dictionary:
    try: 
        
        tol_ss = [1.0e-13, 1.0e-9]  #abs and rel tolerances for steady state problem
        tol_ts = [1.0e-13, 1.0e-9]  #abs and rel tie tolernces for time step function
        
       
        print(f'****************************starting new mole fraction for 2,3,3,3-tetrafluoropropene: {fraction}**************************')

        mole_frac_dict = dictionary

        print(mole_frac_dict)
        gas.TPX = To, Po, mole_frac_dict
        width = 0.08
        flame = ct.FreeFlame(gas, width=width)
        flame.flame.set_steady_tolerances(default=tol_ss)   #set tolerances
        flame.flame.set_transient_tolerances(default=tol_ts)
        #flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1) 
        flame.set_refine_criteria(ratio=5, slope=0.25, curve=0.27)
        flame.max_time_step_count = 1500
        loglevel = 1 
        
        ######################################################################################
#         if i!=0:
#             d = f'./data/range10pts_test_{BTPmole_list[i-1]}.csv'
#             if os.path.exists(d):  
#                 arr2 = ct.SolutionArray(gas)
#                 arr2.read_csv(d)
#                 flame.set_initial_guess(data=arr2)
#                 print(' initial guess has been set')
        #######################################################################################        
                
        #"False" stops the calculation from retrying over and over, thanks Chao 
        #flame.solve(loglevel=loglevel, auto=False)
        flame.solve(loglevel=loglevel, auto=True)
        Su = flame.velocity[0]

        results[fraction] = Su
        sltn = flame.to_solution_array()
        df1 = sltn.to_pandas()
        #edited this here!! index=False
        df1.to_csv(f'{fraction}_{file_name}.csv', index=False)
    except Exception as e: 
        print(f'********************   passed mole fraction:{fraction}, error: {e}    *************************************')
        pass


mole_fraction_of_tetra = list(results.keys())
flame_speeds = list(results.values())


print("mole fractions of 2,3,3,3-tetrafluoropropene are:")
print(mole_fraction_of_tetra)

print("flame speeds are:")
print(flame_speeds)


with open(f'flame_speeds_75_100_{file_name}_species.csv', 'w+') as g:
    g.write(directory)
    g.write('\n')
    writers = csv.writer(g)
    writers.writerow(mole_fraction_of_tetra)
    writers.writerow(flame_speeds)

        

