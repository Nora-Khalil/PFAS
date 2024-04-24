import numpy as np
import re
import os
import csv
import pandas as pd
import copy
import cantera as ct
import sys

#functions 
#ign delay function 

def find_id(time,data):
    m = len(data)
    diff = np.zeros(m)
    diff_l = np.zeros(m)
    diff_r = np.zeros(m)
    a = 0
    for i in range(m):
        if i>1 and i<m-1:
            diff_l[i] = (data[i] - data[i-1])/(time[i]-time[i-1])
            diff_r[i] = (data[i+1] - data[i])/(time[i+1]-time[i])
            diff[i] = (diff_l[i]+diff_r[i])/2
    a = np.max(diff) #slope,m
    b = np.argmax(diff) #location of slope, x1
    
    
    #making of the function (y-y1 = m(x-x1)) 
    y = data[b] #location of y1
    x = time[b]  #location of x1 in terms of the modified time s
    
    #added by Nora 
    print(f'slope (a): {a}, location of slope (b): {b}, y: {y}, x: {x}')
    
    t_id = -y/a+x #ignition delay location
    
    
    #returning the slope approximation function 
    
    lo = x-0.005
    hi = x+0.005
   
    x_0  = np.linspace(lo,hi)
    y_0 = a*(x_0-x)+y
    print(f'in the find_id function, found :{t_id}')
    return diff,t_id,x_0,y_0

#from sys.args 
from_command_line = str(sys.argv[-1]).split('_')
[temperature, x_value] = [float(string) for string in from_command_line]

particular_temps = [temperature]
columns_ = ['reaction_index', f'Temp {temperature}']

#percentages of CH2F2 in propane 
X_list = [x_value]

conc_names = ['CH2F2(1): 0, C3H8(2): 1, O2(3): 6, Ar:93',   #0%
              'CH2F2(1):0.11,C3H8(2):1.1,O2(3):6,Ar:92.13', #0.1%
              'CH2F2(1): 1.99, C3H8(2):0.5, O2(3): 3, Ar:94.7'] #4%

pressure_dict = {0: 1,  #I'm just assuming this is 1 atm for no CH2F2, there's actually no csv file for it
    0.1:  1.0127232051380075,
    0.5:  1.0158269373681306, 
    2: 1.012855580622881,
    4: 0.9891502432791955}

master_dict={}
data_for_csv = []

#load in the gas
directory = '/work/westgroup/nora/Code/projects/PFAS/models/C3H8_CH2F2/C3H8_CH2F2/second_try_on_rebase/chemkin/copies/copy_98_chem_annotated.cti' #initial RMG model
gas = ct.Solution(directory)

dk = 1e-2

#now calculate the initial ign delay: 

for index, x in enumerate(X_list):
    
    if x==4: 
        mole_frac_dict = conc_names[-1]
    if x==0: 
        mole_frac_dict = conc_names[0]
    if x==0.1: 
        mole_frac_dict = conc_names[1]
        
    t_id_nuig = np.zeros(len(particular_temps))
    for i in range(len(particular_temps)):
        print(f'starting {particular_temps[i]} for {x}% of CH2F2')
        
        #now let's perturb
        new_ign_delays = {} #will store the sensitivities for each reaction in this list
        perturb_rxn_indices = list(range(0,len(gas.reactions())))#perturb all reactions one at a time
        # noras_gas.set_multiplier(1)
        for rxn_index in perturb_rxn_indices: 

        
            gas.TPX = particular_temps[i],pressure_dict[x]*ct.one_atm, mole_frac_dict
            gas.set_multiplier(1) #reset all multipliers
            
            print(f'perturbing reaction {rxn_index}')            
            #perturb the rates of the chosen reactions 
            gas.set_multiplier(1+dk, rxn_index) #perturb the k for just this reaction

            #react    
            r = ct.Reactor(contents=gas)
            sim = ct.ReactorNet([r])
            states = ct.SolutionArray(gas, extra=['t'])
        
            dt_max = 5e-6
            t_end = (2000 * dt_max)*4

            while sim.time < t_end:
                sim.advance(sim.time + dt_max)
                states.append(r.thermo.state, t=sim.time*1e3)

            diff,perturbed_t_id,x_0,y_0 = find_id(states.t,states.X[:, gas.species_index('OH(7)')]/np.max(states.X[:, gas.species_index('OH(7)')]))
            data_for_csv.append([rxn_index, perturbed_t_id])
            new_ign_delays[rxn_index] = perturbed_t_id

        
    master_dict[x] = new_ign_delays

print(master_dict)
df = pd.DataFrame(data_for_csv, columns=columns_)
df.to_csv(f'perturbed_ign_delays_{temperature}_{x_value}',index=False)
