import cantera as ct
import numpy as np
import time
import pandas as pd
import os 
import re
import sys
import json

print('Runnning Cantera version: ' + ct.__version__)

cti_index = int(eval(sys.argv[-1]))
full_path_to_ctis = '/work/westgroup/nora/Code/projects/PFAS/models/C3H8_CH2F2/C3H8_CH2F2/with_franklins_edited_rates/chemkin/copies/'


cti_files = [file for file in os.listdir(full_path_to_ctis) if 'cti' in file] 

cti_file = cti_files[cti_index] 
print(cti_file)

mechfile = full_path_to_ctis + cti_file

gas=ct.Solution(mechfile)


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


#doping data, eduardos experiments 

no_CH2F2_temps_inv, no_CH2F2_taus, no_CH2F2_label, no_CH2F2_color, no_CH2F2_shape  = [[0.63269666, 0.56358202, 0.61326454, 0.62338108, 0.6588242,  0.74256908,
 0.77973157, 0.74591045, 0.65343648, 0.6813308], 
 [0.140397, 0.03592654, 0.08382547, 0.1015549,  0.21059732, 1.74050565,
 3.16531584, 1.67252058, 0.19990897, 0.36081361], '0% CH$_2$F$_2$', 'k', 's']
                                                                      
                    
pt1_CH2F2_temps_inv, pt1_CH2F2_taus, pt1_CH2F2_label, pt1_CH2F2_color,  pt1_CH2F2_shape =[[0.71625797, 0.56876523, 0.58520137, 0.6216291, 0.70274495, 0.68888998,
 0.74213088], 
 [1.21382961, 0.0471075,  0.07372893, 0.12092486, 0.67064511, 0.55730096,
 2.02355779], '0.1% CH$_2$F$_2$', 'm', 's']
                                                          
    
pt5_CH2F2_temps_inv, pt5_CH2F2_taus, pt5_CH2F2_label, pt5_CH2F2_color,  pt5_CH2F2_shape  = [[0.64051452, 0.74256407, 0.75806182, 0.59921964, 0.55737473, 0.69714701,
 0.71456021, 0.68253075], 
 [0.19242997, 2.08706923, 2.24174715, 0.08639935, 0.04194953, 0.63637152,
 1.06022888, 0.4540505 ], '0.5% CH$_2$F$_2$', 'c', 'v'] 

_2_CH2F2_temps_inv, _2_CH2F2_taus, _2_CH2F2_label, _2_CH2F2_color,  _2_CH2F2_shape = [[0.6014515, 0.72489571, 0.76260487, 0.63965302, 0.64586785, 0.66499352,
 0.66447004, 0.69194404, 0.55187601], 
 [0.12278134, 2.15618223, 3.8948634,  0.28953096, 0.36193926, 0.57132118,
 0.62294087, 1.16312355, 0.051226],  '2% CH$_2$F$_2$', 'b', '^']

_4_CH2F2_temps_inv, _4_CH2F2_taus, _4_CH2F2_label, _4_CH2F2_color,  _4_CH2F2_shape = [[0.59342286, 0.56586408, 0.61717534, 0.66810691, 0.69578168, 0.71392225,
 0.52988636], 
 [0.2396125, 0.13106904, 0.40859232, 1.46789449, 2.3113845,  2.8922329,
 0.06264836], '4% CH$_2$F$_2$', 'g', 'o']

#this gets specific temperatures
temps_for_0 = sorted([1000/T_inv for T_inv in no_CH2F2_temps_inv])
temps_for_0_1 = sorted([1000/T_inv for T_inv in pt1_CH2F2_temps_inv])
temps_for_0_5 = sorted([1000/T_inv for T_inv in pt5_CH2F2_temps_inv])
temps_for_2 = sorted([1000/T_inv for T_inv in _2_CH2F2_temps_inv])
temps_for_4 = sorted([1000/T_inv for T_inv in _4_CH2F2_temps_inv])

#doped propane (my model) 

#specific Ts

specific_temps = {0: temps_for_0,
0.1: temps_for_0_1,
0.5: temps_for_0_5,
2: temps_for_2,
4: temps_for_4}

X_list = [0,0.1,0.5,2,4]
concentrations = ['0','0.1', '0.5', '2', '4']
conc_names = ['CH2F2(1):0, C3H8(2): 1, O2(3): 6, Ar:93','CH2F2(1):0.11,C3H8(2):1.1,O2(3):6,Ar:92.13','CH2F2(1):0.52,C3H8(2):1,O2(3):6,Ar:92.5','CH2F2(1):2,C3H8(2):1.1,O2(3):6,Ar:90.9','CH2F2(1):1.99,C3H8(2):0.5,O2(3):3,Ar:94.7']

pressure_dict = {0: 1,  #I'm just assuming this is 1 atm for no CH2F2, there's actually no csv file for it
    0.1:  1.0127232051380075,
    0.5:  1.0158269373681306, 
    2: 1.012855580622881,
    4: 0.9891502432791955}

master_dict={}

for index, x in enumerate(X_list): 
    
    temps_we_want = specific_temps[x]
    print(temps_we_want)

    mole_frac_dict = conc_names[index]
    print(conc_names)
    
    t_id_nuig = np.zeros(len(temps_we_want))
    for i in range(len(temps_we_want)):
        print(f'starting {temps_we_want[i]} of {concentrations[index]}% CH2F2')
                        
        gas.TPX = temps_we_want[i],pressure_dict[x]*ct.one_atm, mole_frac_dict
        r = ct.Reactor(contents=gas)


        sim = ct.ReactorNet([r])
        sim.verbose = True

        states = ct.SolutionArray(gas, extra=['t'])
        
                
        dt_max = 5e-6
        t_end = (2000 * dt_max)*4

        while sim.time < t_end:
            sim.advance(sim.time + dt_max)
            states.append(r.thermo.state, t=sim.time*1e3)


        diff,t_id,x_0,y_0 = find_id(states.t,states.X[:, gas.species_index('OH(7)')]/np.max(states.X[:, gas.species_index('OH(7)')]))
        t_id_nuig[i] = t_id #id in ms
        

        #plt.figure()
        #plt.plot(states.t,states.X[:, gas.species_index('OH(7)')]/np.max(states.X[:, gas.species_index('OH(7)')]))
#     #save the data   
    master_dict[x] = t_id_nuig

    print(master_dict)
# with open(f'{cti_file}_ign_delays', 'w') as f:
#     json.dump(master_dict, f)