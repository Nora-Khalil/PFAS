
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
'SulfurHaynes',
]

database(
thermoLibraries = thermolibs,
reactionLibraries = ['halogens_pdep'],
seedMechanisms = ['FFCM1(-)'],
kineticsDepositories = ['training'],
kineticsFamilies = ['default','halogens','Disproportionation-Y'],
frequenciesLibraries = ['halogens_G4'],
kineticsEstimator = 'rate rules',
)

#Only put PFPA in core initially, want to see what RMG will place in the core on its own
species(
    label = 'PFPA',
    reactive = True,
    structure = SMILES('OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F')
)

#in N2 bath gas 
species(
    label = 'N2',
    reactive = False,
    structure = SMILES('N#N')
)


simpleReactor(
        temperature=[(1000,'K'),(2000,'K')],
        pressure=(1.0,'bar'),
        nSims=12,
        initialMoleFractions={
        "PFPA": 0.02, #2% mole fraction
        "N2": 0.98,
        },
        terminationConversion={
        "PFPA": 0.95,
        },
        #terminationRateRatio=1e-6,
        terminationTime=(1,'s'),
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
    #method = 'reservoir state',
    maximumGrainSize=(0.5,'kcal/mol'),
    minimumNumberOfGrains=250,
    temperatures=(300,2500,'K',8),
    pressures=(0.01,100,'bar',5),
    interpolation=('Chebyshev', 6, 4),
    maximumAtoms=16,
)

