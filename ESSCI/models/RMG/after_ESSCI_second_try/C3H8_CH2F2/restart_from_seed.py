restartFromSeed(path='seed')

#second_try_on_rebase
thermo_libs = [
    'BurkeH2O2',
    'C1_C2_Fluorine', #adding Siddha's as second most trusted 
    'NCSU_C2_C8_PFAS', #adding Westmoreland's thermo as the first most trusted
    'primaryThermoLibrary',
    'FFCM1(-)',
    'CurranPentane_edited',
    'Klippenstein_Glarborg2016',
    'thermo_DFT_CCSDTF12_BAC',
    'DFT_QCI_thermo',
    'CBS_QB3_1dHR',
    'Fluorine',
    'halogens',
    'CHOF_G4',
    'CHOCl_G4',
    'CHOBr_G4',
    'CHOFCl_G4',
    'CHOFBr_G4',
    'CHOFClBr_G4',
    '2-BTP_G4',
    'thermo_DFT_CCSDTF12_BAC',
    'SulfurHaynes'
]

kinetic_libs = [
    'FFCM1(-)',
    'CurranPentane_edited',
    'combustion_core/version5',
    'Klippenstein_Glarborg2016',
    'BurkeH2O2inArHe',
    'BurkeH2O2inN2',
    'halogens_pdep',
]


database(
    thermoLibraries = thermo_libs,
    reactionLibraries = kinetic_libs,
    seedMechanisms = ['BurkeH2O2inN2','BurkeH2O2inArHe', 'C2H4+O_Klipp2017'],  # added BurkeH2O2inAr and BurkeH2O2inN2 for bath gases
    kineticsDepositories = ['training'],
    kineticsFamilies = ['default','halogens','Disproportionation-Y'],
    frequenciesLibraries = ['halogens_G4'],
    kineticsEstimator = 'rate rules',
)

#new F abstraction and H abstraction is retrained. In default and halogens. 

species(
    label = 'CH2F2',
    reactive = True,
    structure = SMILES('FCF')
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

    
simpleReactor(
        temperature=[(1000,'K'),(1850,'K')],
        pressure= [(1.0,'bar'),(10.0,'bar')],
        nSims=12,
        initialMoleFractions={
        "C3H8": [.5, 1.5], #stoic is  C3H8:1, O2: 5,  Ar: 94. Add in from 0.1-4% CH2F2: 0.1-4
        "CH2F2": [0, 5],
        "O2": 5,
        "Ar": 94,
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
    temperatures=(300,2000,'K',8), #due to thermo library max pdep=2500K cannot go to 3000K
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
    saveEdgeSpecies = False,
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
# dividing everything by 5 leaves stoichiometric ratios of C3H8/O2/Ar to be .2/1/15.5, although can be anything.
#

#In cantera, I will be using  C3H8/O2/Ar to be 1/5/94,
#