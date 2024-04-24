import cantera as ct
from PIL import Image
from subprocess import run
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import re 
import os
import numpy as np
import csv
import pandas as pd


###################### get the model #######################################


directory = '/work/westgroup/nora/Code/Advanced_Kinetics_Class/project/models/RMG/Fluorine_lib/chemkin/copies/copy_chem_annotated_208.cti'

file_name = directory.split('/')[-1]

# Z=.50, phi = 1
x=.4
phi = x/.4
species_con =  {'C3H2F4(1)': x, 'O2(3)':1, 'N2':1}



########### make the gas, specify the TPX

gas = ct.Solution(directory)
gas.TPX = 298, 1e5, species_con

########### create the FreeFlame object in Cantera, solve the flame
width = 0.08
f = ct.FreeFlame(gas, width=width)
f.set_refine_criteria(ratio=3, slope=0.1, curve=0.1) 
f.max_time_step_count = 900
f.solve()
    
# ########## save a flux diagram 
# diagram = ct.ReactionPathDiagram(f.gas,'H')

# diagram.title = 'Reaction path diagram following H'
# diagram.label_threshold = 0.01

# dot_file = 'rxnpath.dot'
# img_file = 'rxnpath.png'
# img_path = Path.cwd().joinpath(img_file)

# diagram.write_dot(dot_file)
# #print(diagram.get_data())

# print("Wrote graphviz input file to '{0}'.".format(Path.cwd().joinpath(dot_file)))

# run('dot {0} -Tpng -o{1} -Gdpi=200'.format(dot_file, img_file).split())
# print("Wrote graphviz output file to '{0}'.".format(img_path))


########## get sensitivities of each reaction rate on flamespeed
'''
.get_flame_speed_reaction_sensitivities(): Compute the normalized sensitivities of the laminar flame speed S with respect to the reaction rate constants k:

                    s_i =   k_i    dS_u
                           -----  -------
                            S_U    dk_i

'''
######## sort the sensitivities by magnitude #######################################

sens = f.get_flame_speed_reaction_sensitivities()

sensitivity = {}
for m in range(gas.n_reactions):
    sensitivity[m] = abs(sens[m])

sorted_sensitivity = dict(sorted(sensitivity.items(), key=lambda item: item[1], reverse=True)) #sort with highest magnitude first

######### revert the sensitivity values back to original sign  ####################

#sorted_sensitivity_list = [[k,sens[k],gas.reaction(k)] for k,v in sorted_sensitivity.items() ]

data = {
    'k_s': [k for k,v in sorted_sensitivity.items()], #this is number of reaction in gas.reactions list
    'sensitivity': [sens[k] for k,v in sorted_sensitivity.items()], #sensitivity
    'cantera equation': [gas.reaction(k).equation for k,v in sorted_sensitivity.items()],
    'cantera products': [gas.reaction(k).products for k,v in sorted_sensitivity.items()],
    'cantera reactants': [gas.reaction(k).reactants for k,v in sorted_sensitivity.items()],
#     'type': [gas.reaction(k).type for k,v in sorted_sensitivity.items()],
}

df = pd.DataFrame(data)
df.to_csv('sensitivities.csv', index=False)




