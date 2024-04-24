'''
Papas, Zhang, p. 1150: "The stoichiometry for the hydrofluorocarbon-air systems considered was determined by taking the combustion products to be CO2, HF, and H2O. If there was insufficient hydrogen available for formation of HF and H2O, then the formation of HF took preference over H2O formation. If there was insufficient hydrogen available for all the fluorine to form HF, then the remaining fluorine was assumed to produce CF2O in preference of carbon forming CO2. "

Stoichiometric equation for 2,3,3,3-tetrafluoropropene is provided by https://www.sciencedirect.com/science/article/pii/S0022113917301586?via%3Dihub  (Babushok, Linteris): 

Assumed stoichiometric combustion for 2,3,3,3-tetrafluoropropene is: 

CH2CFCF3 + 2.5 (O2 + 3.76 N2) = 2 CO2 + CF20 + 2 HF + 2.5*3.76 N2


calculated phi as: 

at stoichiometric: CH2CFCF3: 0.4, O2: 1, N2: 3.76

therefore phi = (F/M ac) / (F/M stoich) = (x/(x+1+3.76))/(0.4/(0.4+1+3.76)) = (12.9)*(x/(x+1+3.76)) #is this wrong?? why am I doing molar fraction? 

'''

#39% O2 instead of 21%



import cantera as ct
import numpy as np
import pandas as pd
import csv

print("Running Cantera Version: " + str(ct.__version__ ))

To = 298
Po = 1e5 # ct.one_atm
gas = ct.Solution('./chemkin/copies/copy_chem0192.cti')

mole_fractions = list(np.linspace(0.2, 1.0, 25))   #mole fraction of tetra


results = {}

for i in  range(len(mole_fractions)):

    try: 

        x = mole_fractions[i]
        
        mole_frac_dict = {'C3H2F4(1)': x, 'O2(3)':1, 'N2':1.56} # 1/(1+1.56) = .39
        
        #doc string shows how phi equation was found 
#         phi = (12.9)*(x/(x+1+3.76))
        phi = x/.4
        string = f'****************************starting phi: {phi}  **************************'

        gas.TPX = To, Po, mole_frac_dict
        
        width = 0.08
        flame = ct.FreeFlame(gas, width=width)
        flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1) 
        flame.max_time_step_count = 900
        loglevel = 1 
        
        
        try:
            if i!=0:
                d = f'data_39/192/{mole_fractions[i-1]}.csv'
                if os.path.exists(d):  
                    arr2 = ct.SolutionArray(gas)
                    arr2.read_csv(d)
                    flame.set_initial_guess(data=arr2)
                    print('initial guess has been set')
        except: 
            print('initial guess not set for this volume fraction')
        

        flame.solve(loglevel=loglevel, auto=False)
        Su = flame.velocity[0]
        results[phi] = Su
        print(f'at {phi} phi, velocity is {Su}')
        sltn = flame.to_solution_array()
        df1 = sltn.to_pandas()
        df1.to_csv(f'data_39/192/{phi}.csv', index=False)            
            
    except Exception as e: 
        print(f'********************passed phi:{phi}, error: {e}*************************************')
        pass


phis = list(results.keys())
flame_speeds = list(results.values())


print("Phis are:")
print(phis)

print("flame speeds are:")
print(flame_speeds)

    
with open('flamespeeds_192_species.csv', 'w+') as g:
    writers = csv.writer(g)
    writers.writerow(phis)
    writers.writerow(flame_speeds)
        