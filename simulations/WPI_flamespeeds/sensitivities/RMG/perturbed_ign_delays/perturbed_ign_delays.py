import cantera as ct
import numpy as np
import time
import pandas as pd
import os 
import sys

mechfile ='/work/westgroup/nora/Code/projects/PFAS/simulations/WPI_flamespeeds/models/RMG/with_franklins_edits/chemkin/copy_chem_annotated.cti'

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
floats_from_command_line = [float(string) for string in from_command_line]
[temperature, x_value] = floats_from_command_line
particular_temps = [temperature]
columns_ = ['reaction_index', f'Temp {temperature}']


#percentages of CH2F2 in propane 
X_list = [x_value]

#to store sensitivities and data
data_for_csv = []
noras_master_dict_perturbed={}

#perturb k
dk = 1e-2

concentrations_dict = {0.0: 'CH2F2(1):0, C3H8(2): 1, O2(3): 6, Ar:93',
                       0.1: 'CH2F2(1):0.11,C3H8(2):1.1,O2(3):6,Ar:92.13',
                       0.5: 'CH2F2(1):0.52,C3H8(2):1,O2(3):6,Ar:92.5',
                       2.0: 'CH2F2(1):2,C3H8(2):1.1,O2(3):6,Ar:90.9',
                       4.0: 'CH2F2(1):1.99,C3H8(2):0.5,O2(3):3,Ar:94.7'}

pressure_dict = {0.0: 1,  #I'm just assuming this is 1 atm for no CH2F2, there's actually no csv file for it
    0.1:  1.0127232051380075,
    0.5:  1.0158269373681306, 
    2.0: 1.012855580622881,
    4.0: 0.9891502432791955}

for x in X_list:

    mole_fract = concentrations_dict[x]
    

    for i in range(len(particular_temps)):
        print(f'starting {particular_temps[i]} for {x}% of CH2F2')
        
        #now let's perturb
        new_ign_delays = {} #will store the sensitivities for each reaction in this list
        perturb_rxn_indices = list(range(0,1095))#I have 1095 reactions
        # noras_gas.set_multiplier(1)
        for rxn_index in perturb_rxn_indices: 
            
            noras_gas=ct.Solution(mechfile)
            
            noras_gas.TPX = particular_temps[i],pressure_dict[x]*ct.one_atm, mole_fract

            
            print(f'perturbing reaction {rxn_index}')
            #noras_gas.set_multiplier(1) #reset all multipliers
            
            
            #perturb the rates of the chosen reactions 
            noras_gas.set_multiplier(1+dk, rxn_index) #perturb the k for just this reaction

            #react    

            noras_r = ct.Reactor(contents=noras_gas)
            sim_noras = ct.ReactorNet([noras_r])
            #sim_noras.verbose = True
            noras_states = ct.SolutionArray(noras_gas, extra=['t'])

            dt_max = 5e-6
            t_end = (2000 * dt_max)*4
            
            while sim_noras.time < t_end:
                sim_noras.advance(sim_noras.time + dt_max)
                noras_states.append(noras_r.thermo.state, t=sim_noras.time*1e3)


            #find new ignition delay time
            noras_diff_perturbed,noras_t_id_perturbed,noras_x_0_perturbed,noras_y_0_perturbed = find_id(noras_states.t,noras_states.X[:, noras_gas.species_index('OH(7)')]/np.max(noras_states.X[:, noras_gas.species_index('OH(7)')]))
     
            data_for_csv.append([rxn_index, noras_t_id_perturbed])
            new_ign_delays[rxn_index] = noras_t_id_perturbed
            
        noras_gas.set_multiplier(1)
        
    #save the data 
    noras_master_dict_perturbed[x] = new_ign_delays
    
print(noras_master_dict_perturbed)
df = pd.DataFrame(data_for_csv, columns=columns_)
df.to_csv(f'perturbed_ign_delays_{temperature}_{x_value}',index=False)

