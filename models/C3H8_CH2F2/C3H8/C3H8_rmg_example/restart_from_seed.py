restartFromSeed(path='seed')

database(
    thermoLibraries = ['BurkeH2O2','primaryThermoLibrary','thermo_DFT_CCSDTF12_BAC','CBS_QB3_1dHR','DFT_QCI_thermo'],
    reactionLibraries = [],
    seedMechanisms = ['primaryH2O2','ERC-FoundationFuelv0.9'],
    kineticsDepositories = 'default', 
    kineticsFamilies = 'default',
    kineticsEstimator = 'rate rules',
)


species(
    label='C3H8',
    reactive=True,	
    structure=SMILES("CCC"),
)
species(
    label='O2',
    structure=SMILES("[O][O]"),
)
species(
    label='N2',
    reactive=False,
    structure=adjacencyList("""
    1 N u0 p1 c0 {2,T}
	2 N u0 p1 c0 {1,T}
	"""),
)

species(
    label='Ar',    # argon
    reactive=False,
    structure=SMILES("[Ar]"),
)



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

simulator(
    atol=1e-16,
    rtol=1e-8,
)

#model(
#    toleranceMoveToCore=0.3,
#    toleranceInterruptSimulation=0.3,
#    maxNumObjsPerIter=1,
#    terminateAtMaxObjects=True,
#   filterReactions=True,
#    toleranceBranchReactionToCore=0.001,
#    branchingIndex=0.5,
#    branchingRatioMax=1.0,
#)

model(
    toleranceMoveToCore = 0.1,
    toleranceInterruptSimulation = 0.1,
    maximumEdgeSpecies = 5e5,
    filterReactions = True,
    filterThreshold = 5e8,
    minCoreSizeForPrune = 50,
    minSpeciesExistIterationsForPrune = 4,
)

#adding in pressure dependence
pressureDependence(
    method='modified strong collision',
    maximumGrainSize=(0.5,'kcal/mol'),
    minimumNumberOfGrains=250,
    temperatures=(300,2500,'K',8),
    pressures=(0.01,100,'bar',5),
    interpolation=('Chebyshev', 6, 4),
    maximumAtoms=16,
)


options(
    units='si',
    generateSeedEachIteration=False,
    generateOutputHTML=False,
    generatePlots=False,
    saveSimulationProfiles=False,
    verboseComments=False,
    saveEdgeSpecies=False,
    keepIrreversible=False,
)

generatedSpeciesConstraints(
    allowed=['input species','seed mechanisms','reaction libraries'],
    maximumCarbonAtoms=5,
    maximumOxygenAtoms=8,
    maximumNitrogenAtoms=0,
    maximumSiliconAtoms=0,
    maximumSulfurAtoms=0,
    maximumRadicalElectrons=2,
)
