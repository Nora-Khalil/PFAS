
import cantera as ct
import numpy as np
import pandas as pd
import os 
import re
from subprocess import getoutput
import sys
import csv


print("Running Cantera Version: " + str(ct.__version__))



files = ['chem.inp',
 'chem0006.inp',
 'chem0007.inp',
 'chem0008.inp',
 'chem0009.inp',
 'chem0010.inp',
 'chem0011.inp',
 'chem0012.inp',
 'chem0013.inp',
 'chem0014.inp',
 'chem0015.inp',
 'chem0016.inp',
 'chem0017.inp',
 'chem0018.inp',
 'chem0019.inp',
 'chem0020.inp',
 'chem0021.inp',
 'chem0022.inp',
 'chem0023.inp',
 'chem0024.inp',
 'chem0025.inp',
 'chem0026.inp',
 'chem0027.inp',
 'chem0028.inp',
 'chem0029.inp',
 'chem0030.inp',
 'chem0031.inp',
 'chem0032.inp',
 'chem0033.inp',
 'chem0034.inp',
 'chem0035.inp',
 'chem0036.inp',
 'chem0037.inp',
 'chem0038.inp',
 'chem0039.inp',
 'chem0040.inp',
 'chem0041.inp',
 'chem0042.inp',
 'chem0043.inp',
 'chem0044.inp',
 'chem0045.inp',
 'chem0046.inp',
 'chem0047.inp',
 'chem0048.inp',
 'chem0049.inp',
 'chem0050.inp',
 'chem0051.inp',
 'chem0052.inp',
 'chem0053.inp',
 'chem0054.inp',
 'chem0055.inp',
 'chem0056.inp',
 'chem0057.inp',
 'chem0058.inp',
 'chem0059.inp',
 'chem0060.inp',
 'chem0061.inp',
 'chem0062.inp',
 'chem0063.inp',
 'chem_annotated.inp']



def convert(file_name):
    ############### copies chemkin files to dups folder #############################
                    


    #os.system('source activate ct_env') #if on local, this is cantera 2.6 beta
    os.system('source activate cantera_env') #if on discovery, will be cantera 2.5


    #copy folders so i dont screw up the original, and change into the new directory with copies
    os.makedirs('copies', exist_ok=True)

    #copy chem.inp file into dups folder, and will then convert this copy into .cti
    command = f'scp {file_name} copies/copy_184_{file_name}'
    os.system(command)
    os.system('scp tran.dat copies/tran.dat')

    #now look in the dups folder
    os.chdir('./copies')









    ############################# converts the dup_chem.inp files to .cti files #############################

    x = 0
    while x == 0: 
    #this is a string of the output when I try to convert this file to a cti file. Will probably produce an error
        output = getoutput( f'ck2cti --input=copy_184_{file_name} --transport=tran.dat') 
        #if command passed without an error
        if re.search('PASSED',output):
            print('**************************command passed, converting to cti***************************')
            x += 1 
        #if command generated the Duplicate error
        else:
            print('*******************************needs some work****************************************')
            if re.search('Encountered\sunmarked\sduplicate\sreaction',output):
                match = re.search('See\slines\s([0-9]+)\sand\s([0-9]+)\sof\sthe\sinput\sfile',output)
                #capture the line numbers with the duplicate reactions
                line_numbers = [int(match.group(1)), int(match.group(2))]
                print(f'Unmarked duplicates on lines {match.group(1)} and {match.group(2)}')
                print('Editing chemkin file to allow conversion to .cti')
                #write the lines of the chemkin input file to a list so that I can insert the "DUPLICATE" statement
                with open(f'./copy_184_{file_name}','r') as f:
                    data = f.readlines()
                    print(data[line_numbers[0]-1], data[line_numbers[1]-1])

                #start editing the .inp file below

                #'adjustments' will make sure that, even when I add an element in 'data', my index will still be correct
                adjustments = [0,1]
                for i,adjust in zip(line_numbers,adjustments): 
                    start = i+adjust-1
                    count = 0 
                    while count == 0: 
                        #if you don't see a blank line after the duplicated reaction line, keep going until you do
                        if not re.search('^\n', data[start]): 
                            print('no match')
                            print(start)
                            print(data[start])
                            start += 1
                        #when we get to the blank line after the reaction block, insert "DUPLICATE" and stop the loop for this line number
                        else: 
                            print('there is a match')
                            data.insert(start,'DUPLICATE')
                            count = 1 
                #now overwrite the input file with the change 
                with open(f'copy_184_{file_name}','w+') as f: 
                    for l in data: 
                        f.write(l)
                        x==0

            #if the command generated an error that is not the Duplicate error
            else:
                #if code ever gets to here, just cry
                print('There is another error, see Output')
                print(output)
                x += 1
    os.chdir('./..')

for file_name in files: 
    convert(file_name)

            
