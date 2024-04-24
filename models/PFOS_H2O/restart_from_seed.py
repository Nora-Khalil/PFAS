restartFromSeed(path='seed')


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

#Only put PFOS, H2O2, H2O in core initially, want to see what RMG will place in the core on its own
species(
    label = 'PFOS',
    reactive = True,
    structure = SMILES('C(C(C(C(C(F)(F)S(=O)(=O)O)(F)F)(F)F)(F)F)(C(C(C(F)(F)F)(F)F)(F)F)(F)F')
)

#not adding H2O2 yet, this is oxidant in HTL + oxidant (OHTL)

# species(
#     label = 'H2O2',
#     reactive = True,
#     structure = SMILES('OO')
# )

species(
    label = 'H2O',
    reactive = True,
    structure = SMILES('O')
)


#in N2 bath gas, needs to be in an inert atmosphere 
species(
    label = 'N2',
    reactive = False,
    structure = SMILES('N#N')
)


simpleReactor(
        temperature=[(523,'K'),(673,'K')], #250 deg C to 400 deg C, typical for HTL
        pressure=[(50,'bar'),(250,'bar')], #typical pressure conditions for HTL
        nSims=12,
        initialMoleFractions={
        "PFOS": [7.2465e-10,1.44994e-9],#1:10 to 1:5 biomass:water by weight
        "N2": 0.5,
        "H2O":0.4999999999275, #even dewatered sludge is usually 80% water, with more including supplementary water
        },
        terminationConversion={
        "PFOS": 0.95,    
        },
        #terminationRateRatio=1e-6,
        terminationTime=(900,'s'), #residence time for HTL is 15 min - 120 min. 
)
    
model(
    toleranceMoveToCore = 0.05,
    toleranceInterruptSimulation = 0.05,
    maximumEdgeSpecies = 3e5,
    filterReactions = True,
    filterThreshold = 5e8,
    minCoreSizeForPrune = 50,
    minSpeciesExistIterationsForPrune = 4,
)

# pressureDependence(
#     method='modified strong collision',
#     #method = 'reservoir state',
#     maximumGrainSize=(0.5,'kcal/mol'),
#     minimumNumberOfGrains=250,
#     temperatures=(300,2500,'K',8),
#     pressures=(0.01,100,'bar',5),
#     interpolation=('Chebyshev', 6, 4),
#     maximumAtoms=16,
# )


simulator(
    atol = 1e-16,
    rtol = 1e-08,
    sens_atol = 1e-06,
    sens_rtol = 0.0001,
)

#left the same as when David ran his refrigerants
generatedSpeciesConstraints(
    allowed=['input species','seed mechanisms','reaction libraries'],
    maximumCarbonAtoms=6, 
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
    saveEdgeSpecies = True,
    keepIrreversible = True,
    verboseComments = False,
)


