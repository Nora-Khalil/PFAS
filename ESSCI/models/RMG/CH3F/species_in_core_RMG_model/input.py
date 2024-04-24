
thermolibs = [
'primaryThermoLibrary',
'FFCM1(-)',
'halogens',
'CHOF_G4',
'CHOCl_G4',
'CHOBr_G4',
'CHOFCl_G4',
'CHOFBr_G4',
'CHOFClBr_G4',
'DFT_QCI_thermo',
'Fluorine',
'2-BTP_G4',
'thermo_DFT_CCSDTF12_BAC',
'SulfurHaynes'
]

thermolibs_Creg = [
'primaryThermoLibrary',
'FFCM1(-)',
'DFT_QCI_thermo',
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
    label = 'CH3F',
    reactive = True,
    structure = SMILES('CF')
)
    
species(
    label = 'O2',
    reactive = True,
    structure = SMILES('[O][O]')
)
    
species(
    label = 'H2O',
    reactive = True,
    structure = SMILES('O')
)
    
species(
    label = 'N2',
    reactive = False,
    structure = SMILES('N#N')
)
    
species(
    label = 'CH4',
    reactive = True,
    structure = SMILES('C')
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
    label = 'H2(8)',
    reactive = True,
    structure = SMILES('[H][H]')
)
    


species(
    label = 'Ar(9)',
    reactive = True,
    structure = SMILES('[Ar]')
)
    


species(
    label = 'He(10)',
    reactive = True,
    structure = SMILES('[He]')
)
    


species(
    label = 'HO2(11)',
    reactive = True,
    structure = SMILES('[O]O')
)
    


species(
    label = 'H2O2(12)',
    reactive = True,
    structure = SMILES('OO')
)
    


species(
    label = 'CO(13)',
    reactive = True,
    structure = SMILES('[C-]#[O+]')
)
    


species(
    label = 'CO2(14)',
    reactive = True,
    structure = SMILES('O=C=O')
)
    


species(
    label = 'HCO(15)',
    reactive = True,
    structure = SMILES('[CH]=O')
)
    


species(
    label = 'C(T)(16)',
    reactive = True,
    structure = SMILES('[C]')
)
    


species(
    label = 'CH(17)',
    reactive = True,
    structure = SMILES('[CH]')
)
    


species(
    label = 'CH2(T)(18)',
    reactive = True,
    structure = SMILES('[CH2]')
)
    


species(
    label = 'CH3(19)',
    reactive = True,
    structure = SMILES('[CH3]')
)
    


species(
    label = 'CH2O(20)',
    reactive = True,
    structure = SMILES('C=O')
)
    


species(
    label = 'HCCO(21)',
    reactive = True,
    structure = SMILES('C#C[O]')
)
    


species(
    label = 'C2H(22)',
    reactive = True,
    structure = SMILES('[C]#C')
)
    


species(
    label = 'C2H2(23)',
    reactive = True,
    structure = SMILES('C#C')
)
    


species(
    label = 'H2CC(24)',
    reactive = True,
    structure = SMILES('[C]=C')
)
    

    


species(
    label = 'CH3OH(26)',
    reactive = True,
    structure = SMILES('CO')
)
    


species(
    label = 'CH3O(27)',
    reactive = True,
    structure = SMILES('C[O]')
)
    


species(
    label = 'CH2CO(28)',
    reactive = True,
    structure = SMILES('C=C=O')
)
    


species(
    label = 'C2H3(29)',
    reactive = True,
    structure = SMILES('[CH]=C')
)
    


species(
    label = 'C2H4(30)',
    reactive = True,
    structure = SMILES('C=C')
)
    


species(
    label = 'C2H6(31)',
    reactive = True,
    structure = SMILES('CC')
)
    


species(
    label = 'C2H5(32)',
    reactive = True,
    structure = SMILES('C[CH2]')
)
    


species(
    label = 'CH2OH(33)',
    reactive = True,
    structure = SMILES('[CH2]O')
)
    


species(
    label = 'CH3CO(34)',
    reactive = True,
    structure = SMILES('C[C]=O')
)
    


species(
    label = 'CH2CHO(35)',
    reactive = True,
    structure = SMILES('[CH2]C=O')
)
    


species(
    label = 'CH3CHO(36)',
    reactive = True,
    structure = SMILES('CC=O')
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
    label = 'CHF(39)',
    reactive = True,
    structure = SMILES('[CH]F')
)
    


species(
    label = 'CH2F2(40)',
    reactive = True,
    structure = SMILES('FCF')
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
    label = 'CH2F(45)',
    reactive = True,
    structure = SMILES('[CH2]F')
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
    label = 'CH2CHF(55)',
    reactive = True,
    structure = SMILES('C=CF')
)
    


species(
    label = 'CH2F-CH2(61)',
    reactive = True,
    structure = SMILES('[CH2]CF')
)
    


species(
    label = 'CH3-CHF(62)',
    reactive = True,
    structure = SMILES('C[CH]F')
)
    


species(
    label = 'CH2CF(69)',
    reactive = True,
    structure = SMILES('C=[C]F')
)
    


species(
    label = 'CHFCH[Z](70)',
    reactive = True,
    structure = SMILES('[CH]=CF')
)
    


species(
    label = 'CH2CFO(78)',
    reactive = True,
    structure = SMILES('[CH2]C(=O)F')
)
    


species(
    label = 'CHF2(81)',
    reactive = True,
    structure = SMILES('F[CH]F')
)
    


species(
    label = '[O]OF(124)',
    reactive = True,
    structure = SMILES('[O]OF')
)
    


species(
    label = 'O=O(128)',
    reactive = True,
    structure = SMILES('O=O')
)
    


species(
    label = '[O]CF(131)',
    reactive = True,
    structure = SMILES('[O]CF')
)
    


species(
    label = 'OF(133)',
    reactive = True,
    structure = SMILES('OF')
)
    


species(
    label = '[O]OCF(139)',
    reactive = True,
    structure = SMILES('[O]OCF')
)
    


species(
    label = 'O=[C]O(145)',
    reactive = True,
    structure = SMILES('O=[C]O')
)
    


species(
    label = '[O]C=O(146)',
    reactive = True,
    structure = SMILES('[O]C=O')
)
    


species(
    label = '[O]C(=O)F(148)',
    reactive = True,
    structure = SMILES('[O]C(=O)F')
)
    


species(
    label = 'C1OO1(156)',
    reactive = True,
    structure = SMILES('C1OO1')
)
    


species(
    label = '[O]C[O](158)',
    reactive = True,
    structure = SMILES('[O]C[O]')
)
    


species(
    label = 'O=CO(160)',
    reactive = True,
    structure = SMILES('O=CO')
)
    


species(
    label = '[O]C([O])F(165)',
    reactive = True,
    structure = SMILES('[O]C([O])F')
)
    


species(
    label = '[CH]=C[O](181)',
    reactive = True,
    structure = SMILES('[CH]=C[O]')
)
    


species(
    label = 'FC1OO1(186)',
    reactive = True,
    structure = SMILES('FC1OO1')
)
    


species(
    label = '[O]C(F)F(189)',
    reactive = True,
    structure = SMILES('[O]C(F)F')
)
    


species(
    label = 'O=C(O)F(220)',
    reactive = True,
    structure = SMILES('O=C(O)F')
)
    


species(
    label = '[O]C(O)(F)F(290)',
    reactive = True,
    structure = SMILES('[O]C(O)(F)F')
)
    


species(
    label = '[O]OC(=O)F(296)',
    reactive = True,
    structure = SMILES('[O]OC(=O)F')
)
    


species(
    label = '[O]C(=O)O(309)',
    reactive = True,
    structure = SMILES('[O]C(=O)O')
)
    


species(
    label = 'FC1CO1(564)',
    reactive = True,
    structure = SMILES('FC1CO1')
)
    


species(
    label = '[O]C1(F)OO1(1309)',
    reactive = True,
    structure = SMILES('[O]C1(F)OO1')
)
    


species(
    label = '[O]OC(F)F(2051)',
    reactive = True,
    structure = SMILES('[O]OC(F)F')
)
    


species(
    label = 'O[C](F)F(2111)',
    reactive = True,
    structure = SMILES('O[C](F)F')
)
    


species(
    label = 'OC(F)(F)F(2239)',
    reactive = True,
    structure = SMILES('OC(F)(F)F')
)
    


species(
    label = '[O]OC(F)(F)F(2608)',
    reactive = True,
    structure = SMILES('[O]OC(F)(F)F')
)



    
simpleReactor(
        temperature=[(1000,'K'),(2000,'K')],
        pressure=[(1.0,'bar'),(10.0,'bar')],
        nSims=12,
        initialMoleFractions={
        "CH3F": [0.5,1.0],
        "O2": 1,
        "N2": 3.76,
        },
        # terminationConversion={
        # 'CH3F': 0.999,
        # },
        #terminationRateRatio=1e-4,
        #terminationTime=(10,'s'),
        terminationTime=(1,'s'),
        sensitivity=['CH3F','OH'],
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
    
