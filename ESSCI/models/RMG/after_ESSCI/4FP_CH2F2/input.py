 

thermolibs = [
'primaryThermoLibrary',
'Fluorine',
'FFCM1(-)',
'halogens',
'CHOF_G4',
'CHOCl_G4',
'CHOBr_G4',
'CHOFCl_G4',
'CHOFBr_G4',
'CHOFClBr_G4',
'DFT_QCI_thermo',
'2-BTP_G4',
'thermo_DFT_CCSDTF12_BAC',
'SulfurHaynes'
]


database(
thermoLibraries = thermolibs,
reactionLibraries = ['halogens_pdep'],
seedMechanisms= ['FFCM1(-)'],
kineticsDepositories = ['training'],
kineticsFamilies = ['default','halogens','Disproportionation-Y'],
frequenciesLibraries = ['halogens_G4'],
kineticsEstimator = 'rate rules',
)

    
species(
    label = 'N2',
    reactive = False,
    structure = SMILES('N#N')
)
    
    

species(
    label = 'CH2CFCF3',
    reactive = True,
    structure = SMILES('C=C(F)C(F)(F)F')
)      
        
        

species(
    label = '[OH]',
    reactive = True,
    structure = SMILES('[OH]')
)      
        
        

species(
    label = 'O2',
    reactive = True,
    structure = SMILES('[O][O]')
)      
        
        

species(
    label = 'O',
    reactive = True,
    structure = SMILES('O')
)      
        
        

species(
    label = '[H]',
    reactive = True,
    structure = SMILES('[H]')
)      
        
        

species(
    label = '[O]',
    reactive = True,
    structure = SMILES('[O]')
)      
        
        

species(
    label = '[H][H]',
    reactive = True,
    structure = SMILES('[H][H]')
)      
        
        

species(
    label = '[Ar]',
    reactive = True,
    structure = SMILES('[Ar]')
)      
        
        

species(
    label = '[O]O',
    reactive = True,
    structure = SMILES('[O]O')
)      
        
        

species(
    label = '[He]',
    reactive = True,
    structure = SMILES('[He]')
)      
        
        

species(
    label = 'OO',
    reactive = True,
    structure = SMILES('OO')
)      
        
        

species(
    label = '[C-]#[O]',
    reactive = True,
    structure = SMILES('[C-]#[O+]')
)      
        
        

species(
    label = 'O=C=O',
    reactive = True,
    structure = SMILES('O=C=O')
)      
        
        

species(
    label = '[CH]=O',
    reactive = True,
    structure = SMILES('[CH]=O')
)      
        
        

species(
    label = '[C]',
    reactive = True,
    structure = SMILES('[C]')
)      
        
        

species(
    label = '[CH]',
    reactive = True,
    structure = SMILES('[CH]')
)      
        
        

species(
    label = '[CH2]',
    reactive = True,
    structure = SMILES('[CH2]')
)      
        
        

species(
    label = '[CH3]',
    reactive = True,
    structure = SMILES('[CH3]')
)      
        
        

species(
    label = 'C=O',
    reactive = True,
    structure = SMILES('C=O')
)      
        
        

species(
    label = '[CH]=C=O',
    reactive = True,
    structure = SMILES('[CH]=C=O')
)      
        
        

species(
    label = '[C]#C',
    reactive = True,
    structure = SMILES('[C]#C')
)      
        
        

species(
    label = 'C#C',
    reactive = True,
    structure = SMILES('C#C')
)      
        
        

species(
    label = '[C]=C',
    reactive = True,
    structure = SMILES('[C]=C')
)      
        
        

species(
    label = 'CO',
    reactive = True,
    structure = SMILES('CO')
)      
        
        

species(
    label = 'C[O]',
    reactive = True,
    structure = SMILES('C[O]')
)      
        
        

species(
    label = 'C=C=O',
    reactive = True,
    structure = SMILES('C=C=O')
)      
        
        

species(
    label = '[CH]=C',
    reactive = True,
    structure = SMILES('[CH]=C')
)      
        
        

species(
    label = 'C=C',
    reactive = True,
    structure = SMILES('C=C')
)      
        
        

species(
    label = 'C',
    reactive = True,
    structure = SMILES('C')
)      
        
        

species(
    label = 'CC',
    reactive = True,
    structure = SMILES('CC')
)      
        
        

species(
    label = 'C[CH2]',
    reactive = True,
    structure = SMILES('C[CH2]')
)      
        
        

species(
    label = '[CH2]O',
    reactive = True,
    structure = SMILES('[CH2]O')
)      
        
        

species(
    label = 'C[C]=O',
    reactive = True,
    structure = SMILES('C[C]=O')
)      
        
        

species(
    label = '[CH2]C=O',
    reactive = True,
    structure = SMILES('[CH2]C=O')
)      
        
        

species(
    label = 'CC=O',
    reactive = True,
    structure = SMILES('CC=O')
)      
        
        

species(
    label = '[F]',
    reactive = True,
    structure = SMILES('[F]')
)      
        
        

species(
    label = 'F',
    reactive = True,
    structure = SMILES('F')
)      
        
        

species(
    label = '[CH]F',
    reactive = True,
    structure = SMILES('[CH]F')
)      
        
        

species(
    label = 'CH2F2',
    reactive = True,
    structure = SMILES('FCF')
)      
        
        

species(
    label = 'FC(F)F',
    reactive = True,
    structure = SMILES('FC(F)F')
)      
        
        

species(
    label = 'F[C]F',
    reactive = True,
    structure = SMILES('F[C]F')
)      
        
        

species(
    label = 'FC(F)(F)F',
    reactive = True,
    structure = SMILES('FC(F)(F)F')
)      
        
        

species(
    label = 'F[C](F)F',
    reactive = True,
    structure = SMILES('F[C](F)F')
)      
        
        

species(
    label = '[CH2]F',
    reactive = True,
    structure = SMILES('[CH2]F')
)      
        
        

species(
    label = 'O=CF',
    reactive = True,
    structure = SMILES('O=CF')
)      
        
        

species(
    label = '[O]C(F)(F)F',
    reactive = True,
    structure = SMILES('[O]C(F)(F)F')
)      
        
        

species(
    label = 'O=C(F)F',
    reactive = True,
    structure = SMILES('O=C(F)F')
)      
        
        

species(
    label = 'FC=C(F)F',
    reactive = True,
    structure = SMILES('FC=C(F)F')
)      
        
        

species(
    label = 'C=CF',
    reactive = True,
    structure = SMILES('C=CF')
)      
        
        

species(
    label = 'C=C(F)F',
    reactive = True,
    structure = SMILES('C=C(F)F')
)      
        
        

species(
    label = 'C#CF',
    reactive = True,
    structure = SMILES('C#CF')
)      
        
        

species(
    label = 'O=C(F)C(F)(F)F',
    reactive = True,
    structure = SMILES('O=C(F)C(F)(F)F')
)      
        
        

species(
    label = '[CH2]C(=O)F',
    reactive = True,
    structure = SMILES('[CH2]C(=O)F')
)      
        
        

species(
    label = 'C#CC(F)(F)F',
    reactive = True,
    structure = SMILES('C#CC(F)(F)F')
)      
        
        

species(
    label = '[CH]=C(F)C(F)(F)F',
    reactive = True,
    structure = SMILES('[CH]=C(F)C(F)(F)F')
)      
        
        

species(
    label = 'O=C[C](F)C(F)(F)F',
    reactive = True,
    structure = SMILES('O=C[C](F)C(F)(F)F')
)      
        
        

species(
    label = 'C[C](F)C(F)(F)F',
    reactive = True,
    structure = SMILES('C[C](F)C(F)(F)F')
)      
        
        

species(
    label = 'FCC(F)=C(F)F',
    reactive = True,
    structure = SMILES('FCC(F)=C(F)F')
)      
        
        

species(
    label = 'C=C(F)[C](F)F',
    reactive = True,
    structure = SMILES('C=C(F)[C](F)F')
)      
        
        

species(
    label = 'FC[C]C(F)(F)F',
    reactive = True,
    structure = SMILES('FC[C]C(F)(F)F')
)      
        
        

species(
    label = '[CH]C(F)C(F)(F)F',
    reactive = True,
    structure = SMILES('[CH]C(F)C(F)(F)F')
)      
        
        

species(
    label = '[CH2]C(F)(C[C](F)C(F)(F)F)C(F)(F)F',
    reactive = True,
    structure = SMILES('[CH2]C(F)(C[C](F)C(F)(F)F)C(F)(F)F')
)      
        
        

species(
    label = '[CH2]C(F)(F)C(F)(F)F',
    reactive = True,
    structure = SMILES('[CH2]C(F)(F)C(F)(F)F')
)      
        
        

species(
    label = 'FC[C](F)C(F)(F)F',
    reactive = True,
    structure = SMILES('FC[C](F)C(F)(F)F')
)      
        
        

species(
    label = '[CH2]C(F)(O[O])C(F)(F)F',
    reactive = True,
    structure = SMILES('[CH2]C(F)(O[O])C(F)(F)F')
)      
        
        

species(
    label = 'FC=CC(F)(F)F',
    reactive = True,
    structure = SMILES('FC=CC(F)(F)F')
)      
        
        

species(
    label = 'FC(F)=CC(F)F',
    reactive = True,
    structure = SMILES('FC(F)=CC(F)F')
)      
        
        

species(
    label = 'FC=[C]C(F)(F)F',
    reactive = True,
    structure = SMILES('FC=[C]C(F)(F)F')
)      
        
        

species(
    label = 'F[C]CC(F)(F)F',
    reactive = True,
    structure = SMILES('F[C]CC(F)(F)F')
)      
        
        

species(
    label = 'FC=C=C(F)F',
    reactive = True,
    structure = SMILES('FC=C=C(F)F')
)      
        
        

species(
    label = '[CH2]C(F)(F)F',
    reactive = True,
    structure = SMILES('[CH2]C(F)(F)F')
)      
        
        

species(
    label = 'O=O',
    reactive = True,
    structure = SMILES('O=O')
)      
        
        

species(
    label = '[CH2]C([O])(F)C(F)(F)F',
    reactive = True,
    structure = SMILES('[CH2]C([O])(F)C(F)(F)F')
)      
        
        

species(
    label = 'FC(F)(F)C1(F)COO1',
    reactive = True,
    structure = SMILES('FC(F)(F)C1(F)COO1')
)      
        
        

species(
    label = 'FC(F)(F)C1(F)CO1',
    reactive = True,
    structure = SMILES('FC(F)(F)C1(F)CO1')
)      
        
        

species(
    label = 'F[C](F)C[C](F)F',
    reactive = True,
    structure = SMILES('F[C](F)C[C](F)F')
)      
        
        

species(
    label = 'FC(F)[C]C(F)F',
    reactive = True,
    structure = SMILES('FC(F)[C]C(F)F')
)      
        
        

species(
    label = '[CH2]C(F)(C=C(F)C(F)(F)F)C(F)(F)F',
    reactive = True,
    structure = SMILES('[CH2]C(F)(C=C(F)C(F)(F)F)C(F)(F)F')
)      
        
        

species(
    label = 'CC(F)=C(F)F',
    reactive = True,
    structure = SMILES('CC(F)=C(F)F')
)      
        
        

species(
    label = '[O]OCC(F)(F)C(F)(F)F',
    reactive = True,
    structure = SMILES('[O]OCC(F)(F)C(F)(F)F')
)      
        
        

species(
    label = 'C=C(F)C(F)(F)O[O]',
    reactive = True,
    structure = SMILES('C=C(F)C(F)(F)O[O]')
)      
        
        

species(
    label = '[O]OCC(F)=C(F)F',
    reactive = True,
    structure = SMILES('[O]OCC(F)=C(F)F')
)      
        
        

species(
    label = 'C1OO1',
    reactive = True,
    structure = SMILES('C1OO1')
)      
        
        

species(
    label = 'O=C([CH]F)C(F)(F)F',
    reactive = True,
    structure = SMILES('O=C([CH]F)C(F)(F)F')
)      
        
        

species(
    label = 'FC=C(F)C(F)F',
    reactive = True,
    structure = SMILES('FC=C(F)C(F)F')
)      
        
        

species(
    label = 'F[C]1COOC1(F)F',
    reactive = True,
    structure = SMILES('F[C]1COOC1(F)F')
)      
        
        

species(
    label = '[O]OCF',
    reactive = True,
    structure = SMILES('[O]OCF')
)      
        
        

species(
    label = '[O]OC(F)(CF)C(F)(F)F',
    reactive = True,
    structure = SMILES('[O]OC(F)(CF)C(F)(F)F')
)      
        
        

species(
    label = 'F[C](CC(F)(F)F)C(F)(F)F',
    reactive = True,
    structure = SMILES('F[C](CC(F)(F)F)C(F)(F)F')
)      
        
        

species(
    label = '[O]OC(F)(F)F',
    reactive = True,
    structure = SMILES('[O]OC(F)(F)F')
)      
        
        

species(
    label = 'F[CH]C(F)[C](F)F',
    reactive = True,
    structure = SMILES('F[CH]C(F)[C](F)F')
)      
        
        

species(
    label = 'FC=C(F)[C](F)F',
    reactive = True,
    structure = SMILES('FC=C(F)[C](F)F')
)      
        
        

species(
    label = 'FC([CH]C(F)(F)F)C(F)(F)F',
    reactive = True,
    structure = SMILES('FC([CH]C(F)(F)F)C(F)(F)F')
)      
        
        

species(
    label = 'O=CC(F)=C(F)F',
    reactive = True,
    structure = SMILES('O=CC(F)=C(F)F')
)      
        
        

species(
    label = 'CC(F)(O[O])C(F)(F)F',
    reactive = True,
    structure = SMILES('CC(F)(O[O])C(F)(F)F')
)      
        
        

species(
    label = 'OC[C](F)C(F)(F)F',
    reactive = True,
    structure = SMILES('OC[C](F)C(F)(F)F')
)      
        
        

species(
    label = '[O]C[O]',
    reactive = True,
    structure = SMILES('[O]C[O]')
)      
        
        

species(
    label = '[O]C=O',
    reactive = True,
    structure = SMILES('[O]C=O')
)      
        
        

species(
    label = '[O]CC(F)C(F)(F)F',
    reactive = True,
    structure = SMILES('[O]CC(F)C(F)(F)F')
)      
        
        

species(
    label = 'O=CC(F)C(F)(F)F',
    reactive = True,
    structure = SMILES('O=CC(F)C(F)(F)F')
)      
        
        

species(
    label = '[O]OC(F)(CC(F)(F)F)C(F)(F)F',
    reactive = True,
    structure = SMILES('[O]OC(F)(CC(F)(F)F)C(F)(F)F')
)      
        
        

species(
    label = '[O]C(=O)F',
    reactive = True,
    structure = SMILES('[O]C(=O)F')
)      
        
        

species(
    label = 'F[CH][C](F)F',
    reactive = True,
    structure = SMILES('F[CH][C](F)F')
)      
        
        

species(
    label = 'FC1C(F)C1(F)F',
    reactive = True,
    structure = SMILES('FC1C(F)C1(F)F')
)      
        
        

species(
    label = 'O[CH]C(F)C(F)(F)F',
    reactive = True,
    structure = SMILES('O[CH]C(F)C(F)(F)F')
)      
        
        

species(
    label = '[O]C([O])F',
    reactive = True,
    structure = SMILES('[O]C([O])F')
)      
        
        

species(
    label = '[O]OCC(=O)F',
    reactive = True,
    structure = SMILES('[O]OCC(=O)F')
)      
        
        

species(
    label = 'O=C=CF',
    reactive = True,
    structure = SMILES('O=C=CF')
)      
        
        

species(
    label = 'F[C](CC=C(F)C(F)(F)F)C(F)(F)F',
    reactive = True,
    structure = SMILES('F[C](CC=C(F)C(F)(F)F)C(F)(F)F')
)      
        
        

species(
    label = '[O]OC(F)(F)C(F)=CF',
    reactive = True,
    structure = SMILES('[O]OC(F)(F)C(F)=CF')
)      
        
        

species(
    label = '[O]OC(F)C(F)=C(F)F',
    reactive = True,
    structure = SMILES('[O]OC(F)C(F)=C(F)F')
)      
        
        

species(
    label = 'F[C]1C(F)OOC1(F)F',
    reactive = True,
    structure = SMILES('F[C]1C(F)OOC1(F)F')
)      
        
        

species(
    label = 'O=C(F)[CH]F',
    reactive = True,
    structure = SMILES('O=C(F)[CH]F')
)      
        
        

species(
    label = '[O]OC(F)C(=O)C(F)(F)F',
    reactive = True,
    structure = SMILES('[O]OC(F)C(=O)C(F)(F)F')
)      
        
        

species(
    label = '[C]#CF',
    reactive = True,
    structure = SMILES('[C]#CF')
)      
        
        

species(
    label = '[C]=CF',
    reactive = True,
    structure = SMILES('[C]=CF')
)      
        
        

species(
    label = '[CH]=C[O]',
    reactive = True,
    structure = SMILES('[CH]=C[O]')
)      
        
        

species(
    label = 'FC1=CC1(F)F',
    reactive = True,
    structure = SMILES('FC1=CC1(F)F')
)      
        
        

species(
    label = 'F[C](C1CC1(F)C(F)(F)F)C(F)(F)F',
    reactive = True,
    structure = SMILES('F[C](C1CC1(F)C(F)(F)F)C(F)(F)F')
)      
        
        

species(
    label = '[O]OCC(F)(C=C(F)C(F)(F)F)C(F)(F)F',
    reactive = True,
    structure = SMILES('[O]OCC(F)(C=C(F)C(F)(F)F)C(F)(F)F')
)      
        
        

species(
    label = '[O]OC(F)(CC=C(F)C(F)(F)F)C(F)(F)F',
    reactive = True,
    structure = SMILES('[O]OC(F)(CC=C(F)C(F)(F)F)C(F)(F)F')
)      
        
        

species(
    label = 'FC(F)(F)C1(F)[CH]CC(F)(C(F)(F)F)OO1',
    reactive = True,
    structure = SMILES('FC(F)(F)C1(F)[CH]CC(F)(C(F)(F)F)OO1')
)      
        
        

species(
    label = 'F[C](C1CC(F)(C(F)(F)F)OO1)C(F)(F)F',
    reactive = True,
    structure = SMILES('F[C](C1CC(F)(C(F)(F)F)OO1)C(F)(F)F')
)      
        
        

species(
    label = '[O]OC(F)(C1CC(F)(C(F)(F)F)OO1)C(F)(F)F',
    reactive = True,
    structure = SMILES('[O]OC(F)(C1CC(F)(C(F)(F)F)OO1)C(F)(F)F')
)      
        
        

species(
    label = 'FC(F)(F)C1(F)[CH]C(F)(C(F)(F)F)OOC1',
    reactive = True,
    structure = SMILES('FC(F)(F)C1(F)[CH]C(F)(C(F)(F)F)OOC1')
)      
        
        

species(
    label = 'F[C](C1OOCC1(F)C(F)(F)F)C(F)(F)F',
    reactive = True,
    structure = SMILES('F[C](C1OOCC1(F)C(F)(F)F)C(F)(F)F')
)      
        
        

species(
    label = '[O]OC(F)(C=O)C(F)(F)F',
    reactive = True,
    structure = SMILES('[O]OC(F)(C=O)C(F)(F)F')
)      
        
        

species(
    label = '[O]C1OOC1(F)C(F)(F)F',
    reactive = True,
    structure = SMILES('[O]C1OOC1(F)C(F)(F)F')
)      
        
        

species(
    label = 'FC1OO1',
    reactive = True,
    structure = SMILES('FC1OO1')
)      
        
        

species(
    label = 'OC=CF',
    reactive = True,
    structure = SMILES('OC=CF')
)      
        
        

species(
    label = '[O]OC(O)C(F)C(F)(F)F',
    reactive = True,
    structure = SMILES('[O]OC(O)C(F)C(F)(F)F')
)      
        
        

species(
    label = 'FC=COC(F)(F)F',
    reactive = True,
    structure = SMILES('FC=COC(F)(F)F')
)      
        
        

species(
    label = 'O=COO[C](F)C(F)(F)F',
    reactive = True,
    structure = SMILES('O=COO[C](F)C(F)(F)F')
)      
        
        

species(
    label = 'O=[C]C(F)[C](F)F',
    reactive = True,
    structure = SMILES('O=[C]C(F)[C](F)F')
)      
        
        

species(
    label = 'O=C1C(F)C1(F)F',
    reactive = True,
    structure = SMILES('O=C1C(F)C1(F)F')
)      
        
        

species(
    label = 'F[C]1OOC(F)(F)C1F',
    reactive = True,
    structure = SMILES('F[C]1OOC(F)(F)C1F')
)      
        
        

species(
    label = 'O=[C]C(F)(F)[CH]F',
    reactive = True,
    structure = SMILES('O=[C]C(F)(F)[CH]F')
)      
        
        

species(
    label = 'OOC(F)([C]1CC(F)(C(F)(F)F)OO1)C(F)(F)F',
    reactive = True,
    structure = SMILES('OOC(F)([C]1CC(F)(C(F)(F)F)OO1)C(F)(F)F')
)      
        
        

species(
    label = '[O]C(F)(F)C(F)C(=O)F',
    reactive = True,
    structure = SMILES('[O]C(F)(F)C(F)C(=O)F')
)      
        
        

species(
    label = '[O]OC(F)C(=O)F',
    reactive = True,
    structure = SMILES('[O]OC(F)C(=O)F')
)      
        
        

species(
    label = '[O]C(F)(CC(=O)C(F)(OO)C(F)(F)F)C(F)(F)F',
    reactive = True,
    structure = SMILES('[O]C(F)(CC(=O)C(F)(OO)C(F)(F)F)C(F)(F)F')
)      
        
        

species(
    label = 'C=C([O])C(F)(OO)C(F)(F)F',
    reactive = True,
    structure = SMILES('C=C([O])C(F)(OO)C(F)(F)F')
)      
        
        

species(
    label = 'O=C(F)C(=O)CC(F)(F)F',
    reactive = True,
    structure = SMILES('O=C(F)C(=O)CC(F)(F)F')
)      
        
        

species(
    label = 'CC(=O)[C](F)C(F)(F)F',
    reactive = True,
    structure = SMILES('CC(=O)[C](F)C(F)(F)F')
)      
        
        

species(
    label = 'CC(=O)C(F)(O[O])C(F)(F)F',
    reactive = True,
    structure = SMILES('CC(=O)C(F)(O[O])C(F)(F)F')
)      
        
        

species(
    label = '[C]F',
    reactive = True,
    structure = SMILES('[C]F')
)      
        
        

species(
    label = 'O=[C]F',
    reactive = True,
    structure = SMILES('O=[C]F')
)      
        
        

species(
    label = 'FC(F)=C(F)F',
    reactive = True,
    structure = SMILES('FC(F)=C(F)F')
)      
        
        

species(
    label = '[CH2]C(F)F',
    reactive = True,
    structure = SMILES('[CH2]C(F)F')
)      
        
        

species(
    label = 'C[C](F)F',
    reactive = True,
    structure = SMILES('C[C](F)F')
)      
        
        

species(
    label = 'F[CH]C(F)F',
    reactive = True,
    structure = SMILES('F[CH]C(F)F')
)      
        
        

species(
    label = 'FC[C](F)F',
    reactive = True,
    structure = SMILES('FC[C](F)F')
)      
        
        

species(
    label = 'F[C](F)C(F)F',
    reactive = True,
    structure = SMILES('F[C](F)C(F)F')
)      
        
        

species(
    label = 'F[C]=CF',
    reactive = True,
    structure = SMILES('F[C]=CF')
)      
        
        

species(
    label = '[CH]=C(F)F',
    reactive = True,
    structure = SMILES('[CH]=C(F)F')
)      
        
        

species(
    label = 'F[CH]F',
    reactive = True,
    structure = SMILES('F[CH]F')
)      
        
        

species(
    label = '[O]OC(F)F',
    reactive = True,
    structure = SMILES('[O]OC(F)F')
)      
        
        

species(
    label = 'FC(F)C(F)F',
    reactive = True,
    structure = SMILES('FC(F)C(F)F')
)      
        
        

species(
    label = '[O]C(F)F',
    reactive = True,
    structure = SMILES('[O]C(F)F')
)      
        
        

species(
    label = '[O]O[C](F)F',
    reactive = True,
    structure = SMILES('[O]O[C](F)F')
)      
        
        

species(
    label = '[O]C([O])(F)F',
    reactive = True,
    structure = SMILES('[O]C([O])(F)F')
)      
        
        

species(
    label = 'O[C](F)F',
    reactive = True,
    structure = SMILES('O[C](F)F')
)      
        
        

species(
    label = 'FC1(F)OOC1(F)F',
    reactive = True,
    structure = SMILES('FC1(F)OOC1(F)F')
)      
        
        

species(
    label = 'F[C]C(F)F',
    reactive = True,
    structure = SMILES('F[C]C(F)F')
)      
        
        

species(
    label = 'FC1OC1(F)F',
    reactive = True,
    structure = SMILES('FC1OC1(F)F')
)      
        
        

species(
    label = '[O]OC(F)(F)C(F)F',
    reactive = True,
    structure = SMILES('[O]OC(F)(F)C(F)F')
)      
        
        

species(
    label = 'F[C](F)C(F)C(F)F',
    reactive = True,
    structure = SMILES('F[C](F)C(F)C(F)F')
)      
        
        

species(
    label = '[O]C(F)(F)[CH]F',
    reactive = True,
    structure = SMILES('[O]C(F)(F)[CH]F')
)      
        
        

species(
    label = '[O]OCC(F)(F)F',
    reactive = True,
    structure = SMILES('[O]OCC(F)(F)F')
)      
        
        

species(
    label = '[O]OC(F)(F)C=O',
    reactive = True,
    structure = SMILES('[O]OC(F)(F)C=O')
)      
        
        

species(
    label = 'FC1(F)OO1',
    reactive = True,
    structure = SMILES('FC1(F)OO1')
)      
        
        

species(
    label = 'O=C[C](F)F',
    reactive = True,
    structure = SMILES('O=C[C](F)F')
)      
        
        

species(
    label = '[O]OC(F)(F)C(F)C(F)F',
    reactive = True,
    structure = SMILES('[O]OC(F)(F)C(F)C(F)F')
)      
        
        

species(
    label = 'OOC(F)(F)C(F)[C](F)F',
    reactive = True,
    structure = SMILES('OOC(F)(F)C(F)[C](F)F')
)      
        
        

species(
    label = 'F[C](F)CC(F)F',
    reactive = True,
    structure = SMILES('F[C](F)CC(F)F')
)      
        
        

species(
    label = 'FC1(F)CO1',
    reactive = True,
    structure = SMILES('FC1(F)CO1')
)      
        
        

species(
    label = '[O]C1OOC1(F)F',
    reactive = True,
    structure = SMILES('[O]C1OOC1(F)F')
)      
        
        
    
simulator(
    atol = 1e-16,
    rtol = 1e-08,
    sens_atol = 1e-06,
    sens_rtol = 0.0001,
)


generatedSpeciesConstraints(
    allowed=['input species','seed mechanisms','reaction libraries'],
    maximumCarbonAtoms=8,
    maximumOxygenAtoms=6,
    #maximumHeavyAtoms=24,
    maximumRadicalElectrons=2,
    maximumSingletCarbenes=1,
    maximumCarbeneRadicals=0,
    allowSingletO2 = True,
)

options(
    units = "si",
    generateSeedEachIteration = True,
    generateOutputHTML = True,
    generatePlots = True,
    saveSimulationProfiles = True,
    saveEdgeSpecies = False,
    keepIrreversible = True,
    verboseComments = False,
)
    
    
simpleReactor(
        temperature=[(1000,'K'),(2000,'K')],
        pressure= [(1.0,'bar'),(10.0,'bar')],
        nSims=10,
        initialMoleFractions={

        "CH2CFCF3": 0.5,
        "CH2F2": 0.5,
        "O2": 1,
        "N2": 3.76,
        

        },
        # terminationConversion={
        # 'halogen': 0.999,
        # },
        #terminationRateRatio=1e-4,
        #terminationTime=(10,'s'),
        terminationTime=(1,'s'),
        #sensitivity=['halogen','OH'],
        #sensitivityThreshold=0.001,
        )
        
model(
    toleranceMoveToCore = 0.1,
    toleranceInterruptSimulation = 0.1,
    maximumEdgeSpecies = 3e5,
    filterReactions = True,
    filterThreshold = 5e8,
    minCoreSizeForPrune = 50,
    minSpeciesExistIterationsForPrune = 4,
)

pressureDependence(
    method='modified strong collision',
    maximumGrainSize=(0.5,'kcal/mol'),
    minimumNumberOfGrains=250,
    temperatures=(300,2500,'K',8),
    pressures=(0.01,100,'bar',5),
    interpolation=('Chebyshev', 6, 4),
    maximumAtoms=16,
)

    
    
    