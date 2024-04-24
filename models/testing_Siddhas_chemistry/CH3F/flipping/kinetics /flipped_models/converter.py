
import cantera as ct
import numpy as np
import pandas as pd
import os 
import re
from subprocess import getoutput
import sys
import csv


print("Running Cantera Version: " + str(ct.__version__))


files = ['1006.inp',
 '1007.inp',
 '1008.inp',
 '1009.inp',
 '1010.inp',
 '1011.inp',
 '1012.inp',
 '1013.inp',
 '1015.inp',
 '1024.inp',
 '1026.inp',
 '1032.inp',
 '1036.inp',
 '1041.inp',
 '1044.inp',
 '1047.inp',
 '1049.inp',
 '1050.inp',
 '1067.inp',
 '1076.inp',
 '1081.inp',
 '1092.inp',
 '1101.inp',
 '1103.inp',
 '1105.inp',
 '1108.inp',
 '1110.inp',
 '1112.inp',
 '1116.inp',
 '1120.inp',
 '1156.inp',
 '1158.inp',
 '1169.inp',
 '1172.inp',
 '1173.inp',
 '1174.inp',
 '1176.inp',
 '1181.inp',
 '1198.inp',
 '1200.inp',
 '1202.inp',
 '1203.inp',
 '1212.inp',
 '1214.inp',
 '1231.inp',
 '1242.inp',
 '1262.inp',
 '1266.inp',
 '1268.inp',
 '1271.inp',
 '1273.inp',
 '1275.inp',
 '1277.inp',
 '1279.inp',
 '1283.inp',
 '1289.inp',
 '1300.inp',
 '1302.inp',
 '1308.inp',
 '1309.inp',
 '1311.inp',
 '1312.inp',
 '1315.inp',
 '1320.inp',
 '1321.inp',
 '1324.inp',
 '1325.inp',
 '1327.inp',
 '1330.inp',
 '1331.inp',
 '1332.inp',
 '1333.inp',
 '1337.inp',
 '1338.inp',
 '1342.inp',
 '1345.inp',
 '1352.inp',
 '1353.inp',
 '1354.inp',
 '237.inp',
 '238.inp',
 '239.inp',
 '240.inp',
 '241.inp',
 '243.inp',
 '244.inp',
 '245.inp',
 '247.inp',
 '250.inp',
 '254.inp',
 '255.inp',
 '257.inp',
 '259.inp',
 '260.inp',
 '263.inp',
 '265.inp',
 '268.inp',
 '273.inp',
 '274.inp',
 '277.inp',
 '278.inp',
 '280.inp',
 '284.inp',
 '295.inp',
 '299.inp',
 '300.inp',
 '303.inp',
 '307.inp',
 '309.inp',
 '313.inp',
 '316.inp',
 '323.inp',
 '325.inp',
 '326.inp',
 '328.inp',
 '331.inp',
 '338.inp',
 '341.inp',
 '342.inp',
 '343.inp',
 '344.inp',
 '347.inp',
 '348.inp',
 '351.inp',
 '360.inp',
 '361.inp',
 '362.inp',
 '363.inp',
 '376.inp',
 '379.inp',
 '380.inp',
 '381.inp',
 '382.inp',
 '383.inp',
 '385.inp',
 '386.inp',
 '387.inp',
 '388.inp',
 '389.inp',
 '390.inp',
 '391.inp',
 '396.inp',
 '399.inp',
 '401.inp',
 '403.inp',
 '404.inp',
 '409.inp',
 '412.inp',
 '415.inp',
 '416.inp',
 '445.inp',
 '446.inp',
 '447.inp',
 '453.inp',
 '462.inp',
 '463.inp',
 '468.inp',
 '469.inp',
 '470.inp',
 '473.inp',
 '474.inp',
 '475.inp',
 '476.inp',
 '477.inp',
 '484.inp',
 '486.inp',
 '489.inp',
 '495.inp',
 '501.inp',
 '507.inp',
 '508.inp',
 '509.inp',
 '511.inp',
 '512.inp',
 '518.inp',
 '522.inp',
 '523.inp',
 '525.inp',
 '526.inp',
 '527.inp',
 '528.inp',
 '529.inp',
 '533.inp',
 '537.inp',
 '540.inp',
 '544.inp',
 '556.inp',
 '558.inp',
 '559.inp',
 '560.inp',
 '561.inp',
 '562.inp',
 '563.inp',
 '566.inp',
 '567.inp',
 '570.inp',
 '571.inp',
 '574.inp',
 '577.inp',
 '583.inp',
 '584.inp',
 '587.inp',
 '588.inp',
 '591.inp',
 '593.inp',
 '597.inp',
 '599.inp',
 '620.inp',
 '621.inp',
 '623.inp',
 '630.inp',
 '635.inp',
 '644.inp',
 '645.inp',
 '646.inp',
 '647.inp',
 '661.inp',
 '663.inp',
 '673.inp',
 '677.inp',
 '678.inp',
 '680.inp',
 '687.inp',
 '689.inp',
 '690.inp',
 '691.inp',
 '692.inp',
 '693.inp',
 '694.inp',
 '695.inp',
 '696.inp',
 '697.inp',
 '705.inp',
 '706.inp',
 '709.inp',
 '710.inp',
 '712.inp',
 '714.inp',
 '715.inp',
 '718.inp',
 '720.inp',
 '721.inp',
 '722.inp',
 '731.inp',
 '735.inp',
 '738.inp',
 '739.inp',
 '741.inp',
 '742.inp',
 '743.inp',
 '744.inp',
 '746.inp',
 '750.inp',
 '753.inp',
 '754.inp',
 '755.inp',
 '756.inp',
 '757.inp',
 '758.inp',
 '759.inp',
 '763.inp',
 '765.inp',
 '766.inp',
 '767.inp',
 '769.inp',
 '770.inp',
 '773.inp',
 '774.inp',
 '775.inp',
 '776.inp',
 '778.inp',
 '779.inp',
 '780.inp',
 '782.inp',
 '783.inp',
 '784.inp',
 '789.inp',
 '790.inp',
 '796.inp',
 '797.inp',
 '798.inp',
 '799.inp',
 '815.inp',
 '818.inp',
 '820.inp',
 '829.inp',
 '830.inp',
 '831.inp',
 '833.inp',
 '835.inp',
 '845.inp',
 '865.inp',
 '871.inp',
 '872.inp',
 '873.inp',
 '880.inp',
 '881.inp',
 '883.inp',
 '885.inp',
 '887.inp',
 '889.inp',
 '894.inp',
 '902.inp',
 '903.inp',
 '904.inp',
 '905.inp',
 '909.inp',
 '910.inp',
 '918.inp',
 '922.inp',
 '930.inp',
 '954.inp',
 '955.inp',
 '958.inp',
 '961.inp',
 '964.inp',
 '979.inp',
 '980.inp',
 '991.inp',
 '993.inp',
 '996.inp',
 '997.inp']

def convert(file_name):
    ############### copies chemkin files to dups folder #############################


    directory ='.'


    #os.system('source activate ct_env') #if on local, this is cantera 2.6 beta
    os.system('source activate cantera_env') #if on discovery, will be cantera 2.5


    #copy folders so i dont screw up the original, and change into the new directory with copies
    os.makedirs('copies', exist_ok=True)

    #copy chem.inp file into dups folder, and will then convert this copy into .cti
    os.command = f'scp {file_name} copies/copy_{file_name}'
    os.system(os.command)
    os.system('scp tran.dat copies/tran.dat')

    #now look in the dups folder
    os.chdir('./copies')









    ############################# converts the dup_chem.inp files to .cti files #############################

    x = 0
    while x == 0: 
    #this is a string of the output when I try to convert this file to a cti file. Will probably produce an error
        output = getoutput( f'ck2cti --input=copy_{file_name} --transport=tran.dat') 
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
                with open(f'copy_{file_name}','r') as f:
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
                with open(f'copy_{file_name}','w+') as f: 
                    for l in data: 
                        f.write(l)
                        x==0

            #if the command generated an error that is not the Duplicate error
            else:
                #if code ever gets to here, just cry
                print('There is another error, see Output')
                print(output)
                x += 1
    os.chdir('../')


for file_name in files: 
    convert(file_name)
                