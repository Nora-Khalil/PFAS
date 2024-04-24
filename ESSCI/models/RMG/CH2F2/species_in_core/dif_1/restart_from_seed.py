restartFromSeed(path='seed')

# thermolibs = [
# 'primaryThermoLibrary',
# 'FFCM1(-)',
# 'halogens',
# 'CHOF_G4',
# 'CHOCl_G4',
# 'CHOBr_G4',
# 'CHOFCl_G4',
# 'CHOFBr_G4',
# 'CHOFClBr_G4',
# 'DFT_QCI_thermo',
# 'Fluorine',
# '2-BTP_G4',
# 'thermo_DFT_CCSDTF12_BAC',
# 'SulfurHaynes'
# ]

thermolibs = [
'primaryThermoLibrary',
'FFCM1(-)',
'Fluorine',
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
reactionLibraries = ['FFCM1(-)','halogens_pdep'],
seedMechanisms = [],
kineticsDepositories = ['training'],
kineticsFamilies = ['default','halogens','Disproportionation-Y'],
frequenciesLibraries = ['halogens_G4'],
kineticsEstimator = 'rate rules',
)


species(
    label = 'CH2F2',
    reactive = True,
    structure = SMILES('FCF')
)
    
species(
    label = 'O2',
    reactive = True,
    structure = SMILES('[O][O]')
)
    
    
species(
    label = 'N2',
    reactive = False,
    structure = SMILES('N#N')
)

species(
    label = 'H2O',
    reactive = True,
    structure = SMILES('O')
)

species(
    label = 'OH',
    reactive = True,
    structure = SMILES('[OH]')
)

species(
    label = 'H(5)',
    reactive = True,
    structure = SMILES('[H]')
)
    


species(
    label = 'O(6)',
    reactive = True,
    structure = SMILES('[O]')
)
    


species(
    label = 'H2(7)',
    reactive = True,
    structure = SMILES('[H][H]')
)
    


species(
    label = 'HO2(10)',
    reactive = True,
    structure = SMILES('[O]O')
)
    


species(
    label = 'H2O2(11)',
    reactive = True,
    structure = SMILES('OO')
)
    


species(
    label = 'CO(12)',
    reactive = True,
    structure = SMILES('[C-]#[O+]')
)
    


species(
    label = 'CO2(13)',
    reactive = True,
    structure = SMILES('O=C=O')
)
    


species(
    label = 'HCO(14)',
    reactive = True,
    structure = SMILES('[CH]=O')
)
    


species(
    label = 'CH2O(19)',
    reactive = True,
    structure = SMILES('C=O')
)
    


species(
    label = 'F(37)',
    reactive = True,
    structure = SMILES('[F]')
)
    


species(
    label = 'HF(38)',
    reactive = True,
    structure = SMILES('F')
)
    


species(
    label = 'CHF(40)',
    reactive = True,
    structure = SMILES('[CH]F')
)
    


species(
    label = 'CHF3(41)',
    reactive = True,
    structure = SMILES('FC(F)F')
)
    


species(
    label = 'CF2(42)',
    reactive = True,
    structure = SMILES('F[C]F')
)
    


species(
    label = 'CF3(44)',
    reactive = True,
    structure = SMILES('F[C](F)F')
)
    


species(
    label = 'CHFO(46)',
    reactive = True,
    structure = SMILES('O=CF')
)
    


species(
    label = 'CF3O(47)',
    reactive = True,
    structure = SMILES('[O]C(F)(F)F')
)
    


species(
    label = 'CF2O(48)',
    reactive = True,
    structure = SMILES('O=C(F)F')
)
    


species(
    label = 'CF(49)',
    reactive = True,
    structure = SMILES('[C]F')
)
    


species(
    label = 'CFO(50)',
    reactive = True,
    structure = SMILES('O=[C]F')
)
    


species(
    label = 'CHFCF2(54)',
    reactive = True,
    structure = SMILES('FC=C(F)F')
)
    


species(
    label = 'CH2CF2(56)',
    reactive = True,
    structure = SMILES('C=C(F)F')
)
    


species(
    label = 'CF2CF2(60)',
    reactive = True,
    structure = SMILES('FC(F)=C(F)F')
)
    


species(
    label = 'CHF2-CH2(63)',
    reactive = True,
    structure = SMILES('[CH2]C(F)F')
)
    


species(
    label = 'CH3-CF2(64)',
    reactive = True,
    structure = SMILES('C[C](F)F')
)
    


species(
    label = 'CHF2-CHF(66)',
    reactive = True,
    structure = SMILES('F[CH]C(F)F')
)
    


species(
    label = 'CH2F-CF2(67)',
    reactive = True,
    structure = SMILES('FC[C](F)F')
)
    


species(
    label = 'CHF2-CF2(68)',
    reactive = True,
    structure = SMILES('F[C](F)C(F)F')
)
    


species(
    label = 'CHFCF[Z](71)',
    reactive = True,
    structure = SMILES('F[C]=CF')
)
    


species(
    label = 'CF2CH(72)',
    reactive = True,
    structure = SMILES('[CH]=C(F)F')
)
    


species(
    label = 'CHF2(81)',
    reactive = True,
    structure = SMILES('F[CH]F')
)
    


species(
    label = '[O]OC(F)F(128)',
    reactive = True,
    structure = SMILES('[O]OC(F)F')
)
    


species(
    label = 'FC(F)C(F)F(130)',
    reactive = True,
    structure = SMILES('FC(F)C(F)F')
)
    


species(
    label = 'FC1OO1(136)',
    reactive = True,
    structure = SMILES('FC1OO1')
)
    


species(
    label = '[O]C([O])F(138)',
    reactive = True,
    structure = SMILES('[O]C([O])F')
)
    


species(
    label = '[O]C=O(142)',
    reactive = True,
    structure = SMILES('[O]C=O')
)
    


species(
    label = '[O]C(=O)F(143)',
    reactive = True,
    structure = SMILES('[O]C(=O)F')
)
    


species(
    label = '[O]C(F)F(146)',
    reactive = True,
    structure = SMILES('[O]C(F)F')
)
    


species(
    label = '[O]O[C](F)F(147)',
    reactive = True,
    structure = SMILES('[O]O[C](F)F')
)
    


species(
    label = '[O]C([O])(F)F(152)',
    reactive = True,
    structure = SMILES('[O]C([O])(F)F')
)
    


species(
    label = 'O=O(176)',
    reactive = True,
    structure = SMILES('O=O')
)
    


species(
    label = 'O[C](F)F(220)',
    reactive = True,
    structure = SMILES('O[C](F)F')
)
    


species(
    label = 'FC1(F)OOC1(F)F(252)',
    reactive = True,
    structure = SMILES('FC1(F)OOC1(F)F')
)
    


species(
    label = 'F[C]C(F)F(276)',
    reactive = True,
    structure = SMILES('F[C]C(F)F')
)
    


species(
    label = 'FC1OC1(F)F(286)',
    reactive = True,
    structure = SMILES('FC1OC1(F)F')
)
    


species(
    label = '[O]OC(F)(F)C(F)F(406)',
    reactive = True,
    structure = SMILES('[O]OC(F)(F)C(F)F')
)
    


species(
    label = '[O]OC(F)(F)F(416)',
    reactive = True,
    structure = SMILES('[O]OC(F)(F)F')
)
    


species(
    label = 'F[C](F)C(F)C(F)F(444)',
    reactive = True,
    structure = SMILES('F[C](F)C(F)C(F)F')
)
    


species(
    label = '[O]C(F)(F)[CH]F(467)',
    reactive = True,
    structure = SMILES('[O]C(F)(F)[CH]F')
)
    


species(
    label = '[O]OCC(F)(F)F(603)',
    reactive = True,
    structure = SMILES('[O]OCC(F)(F)F')
)
    


species(
    label = '[O]OC(F)(F)C=O(899)',
    reactive = True,
    structure = SMILES('[O]OC(F)(F)C=O')
)
    


species(
    label = 'FC1(F)OO1(932)',
    reactive = True,
    structure = SMILES('FC1(F)OO1')
)
    


species(
    label = '[CH2]C(F)(F)F(962)',
    reactive = True,
    structure = SMILES('[CH2]C(F)(F)F')
)
    


species(
    label = 'O=C[C](F)F(1007)',
    reactive = True,
    structure = SMILES('O=C[C](F)F')
)
    


species(
    label = '[O]OC(F)(F)C(F)C(F)F(1017)',
    reactive = True,
    structure = SMILES('[O]OC(F)(F)C(F)C(F)F')
)
    


species(
    label = 'OOC(F)(F)C(F)[C](F)F(1096)',
    reactive = True,
    structure = SMILES('OOC(F)(F)C(F)[C](F)F')
)
    


species(
    label = 'F[C](F)CC(F)F(1175)',
    reactive = True,
    structure = SMILES('F[C](F)CC(F)F')
)
    


species(
    label = 'FC1(F)CO1(1224)',
    reactive = True,
    structure = SMILES('FC1(F)CO1')
)
    


species(
    label = '[O]C1OOC1(F)F(2126)',
    reactive = True,
    structure = SMILES('[O]C1OOC1(F)F')
)
    
simpleReactor(
        temperature=[(1000,'K'),(2000,'K')],
        pressure=[(1.0,'bar'),(10.0,'bar')],
        nSims=12,
        initialMoleFractions={
        "CH2F2": [0.5,1.0],
        "O2": 1,
        "N2": 3.76,
        },
        terminationConversion={
        'CH2F2': 0.999,
        },
        #terminationRateRatio=1e-4,
        #terminationTime=(10,'s'),
        terminationTime=(1,'s'),
        sensitivity=['CH2F2','OH'],
        sensitivityThreshold=0.001,
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

simulator(
    atol = 1e-16,
    rtol = 1e-08,
    sens_atol = 1e-06,
    sens_rtol = 0.0001,
)

generatedSpeciesConstraints(
    allowed=['input species','seed mechanisms','reaction libraries'],
    maximumCarbonAtoms=4,
    maximumOxygenAtoms=4,
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
    

