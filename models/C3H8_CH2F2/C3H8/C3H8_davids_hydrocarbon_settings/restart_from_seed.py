restartFromSeed(path='seed')

thermolibs_Creg = [
'primaryThermoLibrary',
'FFCM1(-)',
'DFT_QCI_thermo',
]



database(
thermoLibraries = thermolibs_Creg,
reactionLibraries = ['FFCM1(-)'],
seedMechanisms = [],
kineticsDepositories = ['training'],
kineticsFamilies = 'default',
kineticsEstimator = 'rate rules',
)

species(
    label = 'C3H8',
    reactive = True,
    structure = SMILES('CCC')
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
    label='Ar',    # argon
    reactive=False,
    structure=SMILES("[Ar]"),
)

species(
    label = 'N2',
    reactive = False,
    structure = SMILES('N#N')
)

    
# species(
#     label = 'C:',
#     reactive = True,
#     structure = SMILES('[C]')
# )
    
# species(
#     label = 'CH',
#     reactive = True,
#     structure = SMILES('[CH]')
# )
    
simpleReactor(
        temperature=[(1000,'K'),(2000,'K')],
        pressure= [(1.0,'bar'),(10.0,'bar')],
        nSims=12,
        initialMoleFractions={
        "C3H8": [0.1,0.4], #stoic is C3H8: .2, O2: 1, Ar: 15.5
        "O2": 1,
        "Ar": 15.5,
        },
        terminationConversion={
        'C3H8': 0.999,
        },
        terminationRateRatio=1e-8,
        #terminationTime=(10,'s'),
        terminationTime=(1,'s'),
        #sensitivity=['C2H6','OH'],
        #sensitivityThreshold=0.001,
        )
        
model(
    toleranceMoveToCore = 0.1,
    toleranceInterruptSimulation = 0.1,
    maximumEdgeSpecies = 5e5,
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
   # maximumCarbonAtoms=6,
   # maximumOxygenAtoms=4,
   # maximumRadicalElectrons=2,
   # maximumSingletCarbenes=1,
   # maximumCarbeneRadicals=0,
   # allowSingletO2 = True,
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
    

###################################################

# How I got the mole fractions: 
#  
# From Eduardo's presentation, I noticed he didn't use air, but instead used an Ar:93 O2:6 mix
#
#    O2/Ar = .06/.93 = 1/15.5     want O2 to be 1 in ratio

# Stoichiometric equation:

#  C3H8 + 5(O2 + 15.5 Ar) -> 3 CO2 + 4 H2O + 15.5 * 5 Ar 
#
# dividing everything by 5 leaves stoichiometric ratios of C3H8/O2/Ar to be .2/1/15/5
#
#
#
#