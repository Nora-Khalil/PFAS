#!/usr/bin/env python
# encoding: utf-8

name = "Functional Group Additivity Values"
shortDesc = ""
longDesc = """

"""
entry(
    index = 0,
    label = "R",
    group = 
"""
1 * R ux
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1,
    label = "C",
    group = 
"""
1 * C u0
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2,
    label = "CJ2_singlet",
    group = 
"""
1 * C u0 p1
""",
    thermo = 'CsJ2_singlet-CsH',
    shortDesc = """Branch for singlet carbenes""",
    longDesc = 
"""

""",
)

entry(
    index = 3,
    label = "CJ2_singlet-F",
    group = 
"""
1 * C u0 p1 {2,S}
2   F u0 p3 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 4,
    label = "CJ2_singlet-FF",
    group = 
"""
1 * C u0 p1 {2,S} {3,S}
2   F u0 p3 {1,S}
3   F u0 p3 {1,S}
""",
    thermo = None,
    shortDesc = """Derived from fluoro-carbene species in CHOF_G4 library""",
    longDesc = 
"""

""",
)

entry(
    index = 5,
    label = "CJ2_singlet-FC",
    group = 
"""
1 * C u0 p1 {2,S} {3,S}
2   F u0 p3 {1,S}
3   C u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.4991,7.56479,8.3499,8.85944,9.50565,9.86451,9.67794],'cal/(mol*K)','+|-',[0.480518,0.552258,0.564543,0.55626,0.482353,0.416025,0.339051]),
        H298 = (28.8678,'kcal/mol','+|-',1.71301),
        S298 = (32.7378,'cal/(mol*K)','+|-',1.12441),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHOF_G4 |         3
""",
)

entry(
    index = 6,
    label = "CJ2_singlet-FCs",
    group = 
"""
1 * C  u0 p1 {2,S} {3,S}
2   F  u0 p3 {1,S}
3   Cs u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHOF_G4 |         6
""",
)

entry(
    index = 7,
    label = "CJ2_singlet-FCO",
    group = 
"""
1 * C   u0 p1 {2,S} {3,S}
2   CO  u0 {1,S} {4,D}
3   F   u0 p3 {1,S}
4   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.8535,8.27139,8.5114,8.65639,8.91107,9.09017,9.27003],'cal/(mol*K)','+|-',[0.483527,0.555716,0.568078,0.559743,0.485373,0.41863,0.341174]),
        H298 = (39.3047,'kcal/mol','+|-',1.72373),
        S298 = (34.1347,'cal/(mol*K)','+|-',1.13145),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHOF_G4 |         3
""",
)

entry(
    index = 8,
    label = "CJ2_singlet-FO",
    group = 
"""
1 * C u0 p1 {2,S} {3,[S,D]}
2   F u0 p3 {1,S}
3   O u0 {1,[S,D]}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.10766,7.30467,8.40938,9.18815,9.90307,10.4767,10.8485],'cal/(mol*K)','+|-',[0.371969,0.427502,0.437012,0.430601,0.373389,0.322045,0.262459]),
        H298 = (-10.8482,'kcal/mol','+|-',1.32604),
        S298 = (28.6471,'cal/(mol*K)','+|-',0.870406),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHOF_G4 |         5
""",
)

entry(
    index = 9,
    label = "CJ2_singlet-Cl",
    group = 
"""
1 * C  u0 p1 {2,S}
2   Cl u0 p3 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 10,
    label = "CJ2_singlet-ClCl",
    group = 
"""
1 * C  u0 p1 {2,S} {3,S}
2   Cl u0 p3 {1,S}
3   Cl u0 p3 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([11.0595,11.9895,12.5825,12.9238,13.3539,13.5813,13.7185],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (54.3567,'kcal/mol','+|-',2.95256),
        S298 = (64.7687,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHOCl_G4 |         1
""",
)

entry(
    index = 11,
    label = "CJ2_singlet-ClC",
    group = 
"""
1 * C  u0 p1 {2,S} {3,S}
2   Cl u0 p3 {1,S}
3   C  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.01184,8.21551,9.18066,9.74126,10.3431,10.6291,10.5688],'cal/(mol*K)','+|-',[0.590806,0.679011,0.694116,0.683932,0.593062,0.511511,0.41687]),
        H298 = (67.962,'kcal/mol','+|-',2.10617),
        S298 = (34.3094,'cal/(mol*K)','+|-',1.38248),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHOCl_G4 |         2
""",
)

entry(
    index = 12,
    label = "CJ2_singlet-ClCs",
    group = 
"""
1 * C  u0 p1 {2,S} {3,S}
2   Cl u0 p3 {1,S}
3   Cs u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.06955,8.70057,8.97448,9.06346,9.20606,9.23058,9.37988],'cal/(mol*K)','+|-',[0.373122,0.428827,0.438366,0.431935,0.374546,0.323043,0.263273]),
        H298 = (72.893,'kcal/mol','+|-',1.33015),
        S298 = (37.4779,'cal/(mol*K)','+|-',0.873103),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHOCl_G4 |         6
""",
)

entry(
    index = 13,
    label = "CJ2_singlet-ClCO",
    group = 
"""
1 * C   u0 p1 {2,S} {3,S}
2   CO  u0 {1,S} {4,D}
3   Cl  u0 p3 {1,S}
4   O2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 14,
    label = "CJ2_singlet-ClO",
    group = 
"""
1 * C  u0 p1 {2,S} {3,[S,D]}
2   Cl u0 p3 {1,S}
3   O  u0 {1,[S,D]}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.41137,8.27666,9.05676,9.63919,10.1058,10.512,10.6901],'cal/(mol*K)','+|-',[0.372153,0.427714,0.437228,0.430814,0.373573,0.322204,0.262589]),
        H298 = (38.7115,'kcal/mol','+|-',1.32669),
        S298 = (31.8433,'cal/(mol*K)','+|-',0.870836),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHOCl_G4 |         5
""",
)

entry(
    index = 15,
    label = "CJ2_singlet-Br",
    group = 
"""
1 * C  u0 p1 {2,S}
2   Br u0 p3 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 16,
    label = "CJ2_singlet-BrBr",
    group = 
"""
1 * C  u0 p1 {2,S} {3,S}
2   Br u0 p3 {1,S}
3   Br u0 p3 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([11.8507,12.5829,12.9374,13.188,13.522,13.693,13.7791],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (80.9121,'kcal/mol','+|-',2.95256),
        S298 = (70.3965,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHOBr_G4 |         1
""",
)

entry(
    index = 17,
    label = "CJ2_singlet-BrC",
    group = 
"""
1 * C  u0 p1 {2,S} {3,S}
2   Br u0 p3 {1,S}
3   C  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.42867,8.57692,9.46091,10.0323,10.5663,10.6839,10.5164],'cal/(mol*K)','+|-',[0.590396,0.67854,0.693634,0.683458,0.59265,0.511155,0.41658]),
        H298 = (80.6611,'kcal/mol','+|-',2.10471),
        S298 = (36.9649,'cal/(mol*K)','+|-',1.38152),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHOBr_G4 |         2
""",
)

entry(
    index = 18,
    label = "CJ2_singlet-BrCs",
    group = 
"""
1 * C  u0 p1 {2,S} {3,S}
2   Br u0 p3 {1,S}
3   Cs u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.33757,8.77681,8.98166,9.08369,9.31807,9.42932,9.70167],'cal/(mol*K)','+|-',[0.416171,0.478303,0.488943,0.48177,0.41776,0.360314,0.293648]),
        H298 = (85.6964,'kcal/mol','+|-',1.48361),
        S298 = (40.0387,'cal/(mol*K)','+|-',0.973838),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHOBr_G4 |         4
""",
)

entry(
    index = 19,
    label = "CJ2_singlet-BrCO",
    group = 
"""
1 * C   u0 p1 {2,S} {3,S}
2   CO  u0 {1,S} {4,D}
3   Br  u0 p3 {1,S}
4   O2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 20,
    label = "CJ2_singlet-BrO",
    group = 
"""
1 * C  u0 p1 {2,S} {3,[S,D]}
2   Br u0 p3 {1,S}
3   O  u0 {1,[S,D]}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.12626,8.76445,9.41893,9.88812,10.1502,10.4411,10.5849],'cal/(mol*K)','+|-',[0.372431,0.428034,0.437555,0.431136,0.373853,0.322445,0.262785]),
        H298 = (52.3298,'kcal/mol','+|-',1.32769),
        S298 = (34.4774,'cal/(mol*K)','+|-',0.871488),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHOBr_G4 |         5
""",
)

entry(
    index = 21,
    label = "CsJ2_singlet-HH",
    group = 
"""
1 * C2s u0 p1 {2,S} {3,S}
2   H   u0 {1,S}
3   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Fitted to DFT_QCI_thermo library""",
    longDesc = 
"""
Fitted to RQCISD(T)/cc-PV(infinity)(Q)Z calculations of:

Goldsmith, C. F.; Magoon, G. R.; Green, W. H., Database of Small Molecule Thermochemistry for Combustion.
J. Phys. Chem. A 2012, 116, 9033-9057.
""",
)

entry(
    index = 22,
    label = "CsJ2_singlet-OsH",
    group = 
"""
1 * C2s u0 p1 {2,S} {3,S}
2   O2s u0 {1,S}
3   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.075,5.312,6.211,6.926,8.355,9.557,10.212],'cal/(mol*K)'),
        H298 = (65.592,'kcal/mol'),
        S298 = (23.749,'cal/(mol*K)'),
    ),
    shortDesc = """Fitted to DFT_QCI_thermo library""",
    longDesc = 
"""
Fitted to RQCISD(T)/cc-PV(infinity)(Q)Z calculations of:

Goldsmith, C. F.; Magoon, G. R.; Green, W. H., Database of Small Molecule Thermochemistry for Combustion.
J. Phys. Chem. A 2012, 116, 9033-9057.
""",
)

entry(
    index = 23,
    label = "CdJ2_singlet-Od",
    group = 
"""
1 * C2d u0 p1 {2,D}
2   O2d u0 p2 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([-1.5,-2.38,-3.32,-4.24,-5.75,-6.88,-8.59],'cal/(mol*K)'),
        H298 = (103.73,'kcal/mol'),
        S298 = (-6.47,'cal/(mol*K)'),
    ),
    shortDesc = """Calculated in relation to formaldehyde from NIST values""",
    longDesc = 
"""

""",
)

entry(
    index = 24,
    label = "CdJ2_singlet-Sd",
    group = 
"""
1 * C2d u0 p1 {2,D}
2   S2d u0 p2 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([-1.97,-2.97,-3.85,-4.6,-5.82,-6.79,-8.44],'cal/(mol*K)'),
        H298 = (143.53,'kcal/mol'),
        S298 = (-6.23,'cal/(mol*K)'),
    ),
    shortDesc = """CBS-QB3 GA 1D-HR Aaron Vandeputte 2009""",
    longDesc = 
"""

""",
)

entry(
    index = 25,
    label = "CdJ2_singlet-(Cdd-Od)",
    group = 
"""
1   Cdd u0 {2,D} {3,D}
2 * C2d u0 p1 {1,D}
3   O2d u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([10.135,11.201,11.749,12.051,12.813,13.581,14.122],'cal/(mol*K)'),
        H298 = (110.367,'kcal/mol'),
        S298 = (53.61,'cal/(mol*K)'),
    ),
    shortDesc = """Fitted to DFT_QCI_thermo library""",
    longDesc = 
"""
Fitted to RQCISD(T)/cc-PV(infinity)(Q)Z calculations of:

Goldsmith, C. F.; Magoon, G. R.; Green, W. H., Database of Small Molecule Thermochemistry for Combustion.
J. Phys. Chem. A 2012, 116, 9033-9057.
""",
)

entry(
    index = 26,
    label = "CsJ2_singlet-CH",
    group = 
"""
1 * C2s u0 p1 {2,S} {3,S}
2   C   u0 {1,S}
3   H   u0 {1,S}
""",
    thermo = 'CsJ2_singlet-CsH',
    shortDesc = """Branch for singlet carbenes single-bonded to one carbon and one hydrogen""",
    longDesc = 
"""

""",
)

entry(
    index = 27,
    label = "CsJ2_singlet-CsH",
    group = 
"""
1 * C2s u0 p1 {2,S} {3,S}
2   Cs  u0 {1,S}
3   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Fitted to DFT_QCI_thermo library""",
    longDesc = 
"""
Fitted to RQCISD(T)/cc-PV(infinity)(Q)Z calculations of:

Goldsmith, C. F.; Magoon, G. R.; Green, W. H., Database of Small Molecule Thermochemistry for Combustion.
J. Phys. Chem. A 2012, 116, 9033-9057.
""",
)

entry(
    index = 28,
    label = "CsJ2_singlet-CtH",
    group = 
"""
1 * C2s u0 p1 {2,S} {3,S}
2   Ct  u0 {1,S}
3   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.675,6.472,6.776,6.9,7.469,8.156,7.58],'cal/(mol*K)'),
        H298 = (88.247,'kcal/mol'),
        S298 = (28.407,'cal/(mol*K)'),
    ),
    shortDesc = """Fitted to DFT_QCI_thermo library""",
    longDesc = 
"""
Fitted to RQCISD(T)/cc-PV(infinity)(Q)Z calculations of:

Goldsmith, C. F.; Magoon, G. R.; Green, W. H., Database of Small Molecule Thermochemistry for Combustion.
J. Phys. Chem. A 2012, 116, 9033-9057.
""",
)

entry(
    index = 29,
    label = "CdJ2_singlet-Cd",
    group = 
"""
1   C   u0 {2,D}
2 * C2d u0 p1 {1,D}
""",
    thermo = 'CdJ2_singlet-Cds',
    shortDesc = """Branch for singlet carbenes double-bonded to one carbon""",
    longDesc = 
"""

""",
)

entry(
    index = 30,
    label = "CdJ2_singlet-Cds",
    group = 
"""
1   Cd  u0 {2,D}
2 * C2d u0 p1 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.372,5.272,4.916,4.506,4.219,4.263,3.97],'cal/(mol*K)'),
        H298 = (92.157,'kcal/mol'),
        S298 = (26.864,'cal/(mol*K)'),
    ),
    shortDesc = """Fitted to DFT_QCI_thermo library""",
    longDesc = 
"""
Fitted to RQCISD(T)/cc-PV(infinity)(Q)Z calculations of:

Goldsmith, C. F.; Magoon, G. R.; Green, W. H., Database of Small Molecule Thermochemistry for Combustion.
J. Phys. Chem. A 2012, 116, 9033-9057.
""",
)

entry(
    index = 31,
    label = "CdJ2_singlet-(Cdd-Cds)",
    group = 
"""
1   Cdd u0 {2,D} {3,D}
2 * C2d u0 p1 {1,D}
3   C   u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.406,3.333,3.175,3.019,3.156,3.468,3.593],'cal/(mol*K)'),
        H298 = (91.676,'kcal/mol'),
        S298 = (26.434,'cal/(mol*K)'),
    ),
    shortDesc = """Fitted to DFT_QCI_thermo library""",
    longDesc = 
"""
Fitted to RQCISD(T)/cc-PV(infinity)(Q)Z calculations of:

Goldsmith, C. F.; Magoon, G. R.; Green, W. H., Database of Small Molecule Thermochemistry for Combustion.
J. Phys. Chem. A 2012, 116, 9033-9057.
""",
)

entry(
    index = 32,
    label = "CsJ2_singlet-(Cds-Cds-Cds-C)C",
    group = 
"""
1   Cd  u0 {2,D} {4,S}
2   Cd  u0 {1,D} {3,S}
3   Cd  u0 {2,S} {5,D}
4 * C2s u0 p1 {1,S} {6,S}
5   C   u0 {3,D}
6   C   u0 {4,S}
""",
    thermo = 'CsJ2_singlet-(Cds-Cds-Cds-Cds)Cs_6_ring',
    shortDesc = """Branch for singlet carbenes delocalized over two conjugated carbon double bonds""",
    longDesc = 
"""

""",
)

entry(
    index = 33,
    label = "CsJ2_singlet-(Cds-Cds-Cds-Cds)Cs_5_ring",
    group = 
"""
1   Cd  u0 {2,D}
2   Cd  u0 {1,D} {3,S} {6,S}
3   Cs  u0 {2,S} {4,S}
4 * C2s u0 p1 {3,S} {5,S}
5   Cd  u0 {4,S} {6,D}
6   Cd  u0 {2,S} {5,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.025,5.12,5.268,4.917,4.799,5.731,5.087],'cal/(mol*K)'),
        H298 = (82.041,'kcal/mol'),
        S298 = (10.325,'cal/(mol*K)'),
    ),
    shortDesc = """Fitted to DFT_QCI_thermo library""",
    longDesc = 
"""
Fitted to RQCISD(T)/cc-PV(infinity)(Q)Z calculations of:

    Goldsmith, C. F.; Magoon, G. R.; Green, W. H., Database of Small Molecule Thermochemistry for Combustion.
    J. Phys. Chem. A 2012, 116, 9033-9057.
""",
)

entry(
    index = 34,
    label = "CsJ2_singlet-(Cds-Cds-Cds-Cds)Cs_6_ring",
    group = 
"""
1   Cd  u0 {2,D} {3,S}
2   Cd  u0 {1,D} {6,S}
3   Cs  u0 {1,S} {4,S}
4 * C2s u0 p1 {3,S} {5,S}
5   Cd  u0 {4,S} {6,D}
6   Cd  u0 {2,S} {5,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.322,5.406,5.554,5.182,5.008,5.927,5.232],'cal/(mol*K)'),
        H298 = (80.263,'kcal/mol'),
        S298 = (12.963,'cal/(mol*K)'),
    ),
    shortDesc = """Fitted to DFT_QCI_thermo library""",
    longDesc = 
"""
Fitted to RQCISD(T)/cc-PV(infinity)(Q)Z calculations of:

    Goldsmith, C. F.; Magoon, G. R.; Green, W. H., Database of Small Molecule Thermochemistry for Combustion.
    J. Phys. Chem. A 2012, 116, 9033-9057.
""",
)

entry(
    index = 35,
    label = "CtJ2_singlet-N5tc",
    group = 
"""
1 * C2tc u0 p1 {2,T}
2   N5tc u0 p0 {1,T}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.81825,6.27186,6.91034,7.64029,8.81025,9.64443,10.709],'cal/(mol*K)','+|-',[1.13738,1.2147,1.25274,1.24384,1.20996,1.18825,1.21246]),
        H298 = (52.0717,'kcal/mol','+|-',4.55876),
        S298 = (33.093,'cal/(mol*K)','+|-',3.49275),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library                 | Number of Species
thermo_DFT_CCSDTF12_BAC |         2
""",
)

entry(
    index = 36,
    label = "CsJ2_singlet-N5tc",
    group = 
"""
1 * C2sc u0 p1 {2,S}
2   N5tc u0 p0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([2.61222,3.80065,4.78208,5.60618,6.82855,7.71398,9.12752],'cal/(mol*K)','+|-',[1.3559,1.44808,1.49342,1.48282,1.44243,1.41654,1.4454]),
        H298 = (48.6834,'kcal/mol','+|-',5.43462),
        S298 = (5.55483,'cal/(mol*K)','+|-',4.1638),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
primaryNS |         1
""",
)

entry(
    index = 37,
    label = "Cbf",
    group = 
"""
1 * Cbf u0
""",
    thermo = 'Cbf-CbCbCbf',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 38,
    label = "Cbf-CbCbCbf",
    group = 
"""
1 * Cbf u0 {2,B} {3,B} {4,B}
2   Cb  u0 {1,B}
3   Cb  u0 {1,B}
4   Cbf u0 {1,B}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.01,3.68,4.2,4.61,5.2,5.7,6.2],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (4.8,'kcal/mol','+|-',0.17),
        S298 = (-5,'cal/(mol*K)','+|-',0.1),
    ),
    shortDesc = """Cbf-CbfCbCb STEIN and FAHR; J. PHYS. CHEM. 1985, 89, 17, 3714""",
    longDesc = 
"""

""",
)

entry(
    index = 39,
    label = "Cbf-CbCbfCbf",
    group = 
"""
1 * Cbf u0 {2,B} {3,B} {4,B}
2   Cb  u0 {1,B}
3   Cbf u0 {1,B}
4   Cbf u0 {1,B}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.01,3.68,4.2,4.61,5.2,5.7,6.2],'cal/(mol*K)','+|-',[0.15,0.15,0.15,0.15,0.15,0.15,0.15]),
        H298 = (3.7,'kcal/mol','+|-',0.3),
        S298 = (-5,'cal/(mol*K)','+|-',0.15),
    ),
    shortDesc = """Cbf-CbfCbfCb STEIN and FAHR; J. PHYS. CHEM. 1985, 89, 17, 3714""",
    longDesc = 
"""

""",
)

entry(
    index = 40,
    label = "Cbf-CbfCbfCbf",
    group = 
"""
1  * Cbf u0 p0 c0 {2,B} {3,B} {6,B}
2    Cbf u0 p0 c0 {1,B} {4,B} {5,B}
3    Cbf u0 p0 c0 {1,B} {8,B} {9,B}
4    Cbf u0 p0 c0 {2,B} {10,B} {11,B}
5    Cbf u0 p0 c0 {2,B} {13,B} {14,B}
6    Cbf u0 p0 c0 {1,B} {15,B} {16,B}
7    C   u0 p0 c0 {8,B} {16,B}
8    C   u0 p0 c0 {3,B} {7,B}
9    C   u0 p0 c0 {3,B} {10,B}
10   C   u0 p0 c0 {4,B} {9,B}
11   C   u0 p0 c0 {4,B} {12,B}
12   C   u0 p0 c0 {11,B} {13,B}
13   C   u0 p0 c0 {5,B} {12,B}
14   C   u0 p0 c0 {5,B} {15,B}
15   C   u0 p0 c0 {6,B} {14,B}
16   C   u0 p0 c0 {6,B} {7,B}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([2,3.11,3.9,4.42,5,5.3,5.7],'cal/(mol*K)','+|-',[0.15,0.15,0.15,0.15,0.15,0.15,0.15]),
        H298 = (1.5,'kcal/mol','+|-',0.3),
        S298 = (1.8,'cal/(mol*K)','+|-',0.15),
    ),
    shortDesc = """Cbf-CbfCbfCbf STEIN and FAHR; J. PHYS. CHEM. 1985, 89, 17, 3714""",
    longDesc = 
"""
The smallest PAH that can have Cbf-CbfCbfCbf is pyrene. Currently the database is restricted
that any group with more three Cbf atoms must have all benzene rings explicitly written out.
Previously, this node would also match one carbon on Benzo[c]phenanthrene and does not now.
Examples from the original source do not include Benzo[c]phenanthrene.
""",
)

entry(
    index = 41,
    label = "Cb",
    group = 
"""
1 * Cb u0
""",
    thermo = 'Cb-Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 42,
    label = "Cb-H",
    group = 
"""
1 * Cb u0 {2,S}
2   H  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.24,4.44,5.46,6.3,7.54,8.41,9.73],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (3.3,'kcal/mol','+|-',0.11),
        S298 = (11.53,'cal/(mol*K)','+|-',0.12),
    ),
    shortDesc = """Cb-H BENSON""",
    longDesc = 
"""

""",
)

entry(
    index = 43,
    label = "Cb-O2s",
    group = 
"""
1 * Cb  u0 {2,S}
2   O2s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.9,5.3,6.2,6.6,6.9,6.9,7.07],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (-0.9,'kcal/mol','+|-',0.16),
        S298 = (-10.2,'cal/(mol*K)','+|-',0.1),
    ),
    shortDesc = """Cb-O BENSON Cp1500=3D Cp1000*(Cp1500/Cp1000: Cb/Cd)""",
    longDesc = 
"""

""",
)

entry(
    index = 44,
    label = "Cb-S",
    group = 
"""
1 * Cb u0 {2,S}
2   S  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.75,3.75,3.47,3.5,5.18,6.15,4.65],'cal/(mol*K)'),
        H298 = (18.76,'kcal/mol'),
        S298 = (-0.62,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 45,
    label = "Cb-C",
    group = 
"""
1 * Cb u0 {2,S}
2   C  u0 {1,S}
""",
    thermo = 'Cb-Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 46,
    label = "Cb-Cs",
    group = 
"""
1 * Cb u0 {2,S}
2   Cs u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([2.67,3.14,3.68,4.15,4.96,5.44,5.98],'cal/(mol*K)','+|-',[0.07,0.07,0.07,0.07,0.07,0.07,0.07]),
        H298 = (5.51,'kcal/mol','+|-',0.13),
        S298 = (-7.69,'cal/(mol*K)','+|-',0.1),
    ),
    shortDesc = """Cb-Cs BENSON""",
    longDesc = 
"""

""",
)

entry(
    index = 47,
    label = "Cb-Cds",
    group = 
"""
1 * Cb         u0 {2,S}
2   [Cd,CO,CS] u0 {1,S}
""",
    thermo = 'Cb-(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 48,
    label = "Cb-(Cds-O2d)",
    group = 
"""
1 * Cb  u0 {2,S}
2   CO  u0 {1,S} {3,D}
3   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.59,3.97,4.38,4.72,5.28,5.61,5.75],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (3.69,'kcal/mol','+|-',0.2),
        S298 = (-7.8,'cal/(mol*K)','+|-',0.1),
    ),
    shortDesc = """Enthalpy from Cb-CO, entropies and heat capacities assigned value of Cb-Cd""",
    longDesc = 
"""

""",
)

entry(
    index = 49,
    label = "Cb-(Cds-Cd)",
    group = 
"""
1 * Cb u0 {2,S}
2   Cd u0 {1,S} {3,D}
3   C  u0 {2,D}
""",
    thermo = 'Cb-(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 50,
    label = "Cb-(Cds-Cds)",
    group = 
"""
1 * Cb u0 {2,S}
2   Cd u0 {1,S} {3,D}
3   Cd u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.59,3.97,4.38,4.72,5.28,5.61,5.75],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (5.69,'kcal/mol','+|-',0.2),
        S298 = (-7.8,'cal/(mol*K)','+|-',0.1),
    ),
    shortDesc = """Cb-Cd STEIN and FAHR; J. PHYS. CHEM. 1985, 89, 17, 3714""",
    longDesc = 
"""

""",
)

entry(
    index = 51,
    label = "Cb-(Cds-Cdd)",
    group = 
"""
1 * Cb  u0 {2,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D}
""",
    thermo = 'Cb-(Cds-Cdd-Cd)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 52,
    label = "Cb-(Cds-Cdd-O2d)",
    group = 
"""
1 * Cb  u0 {2,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {4,D}
4   O2d u0 {3,D}
""",
    thermo = 'Cb-(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 53,
    label = "Cb-(Cds-Cdd-S2d)",
    group = 
"""
1 * Cb  u0 {2,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {4,D}
4   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 54,
    label = "Cb-(Cds-Cdd-Cd)",
    group = 
"""
1 * Cb  u0 {2,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {4,D}
4   C   u0 {3,D}
""",
    thermo = 'Cb-(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 55,
    label = "Cb-Ct",
    group = 
"""
1 * Cb u0 {2,S}
2   Ct u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.59,3.97,4.38,4.72,5.28,5.61,5.75],'cal/(mol*K)','+|-',[0.15,0.15,0.15,0.15,0.15,0.15,0.15]),
        H298 = (5.69,'kcal/mol','+|-',0.3),
        S298 = (-7.8,'cal/(mol*K)','+|-',0.15),
    ),
    shortDesc = """Cb-Ct STEIN and FAHR; J. PHYS. CHEM. 1985, 89, 17, 3714""",
    longDesc = 
"""

""",
)

entry(
    index = 56,
    label = "Cb-(CtN3t)",
    group = 
"""
1 * Cb  u0 {2,S}
2   Ct  u0 {1,S} {3,T}
3   N3t u0 {2,T}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.8,11.2,12.3,13.1,14.2,14.9,16.65],'cal/(mol*K)'),
        H298 = (35.8,'kcal/mol'),
        S298 = (20.5,'cal/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 57,
    label = "Cb-Cb",
    group = 
"""
1 * Cb u0 {2,S}
2   Cb u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.33,4.22,4.89,5.27,5.76,5.95,6.05],'cal/(mol*K)','+|-',[0.15,0.15,0.15,0.15,0.15,0.15,0.15]),
        H298 = (4.96,'kcal/mol','+|-',0.3),
        S298 = (-8.64,'cal/(mol*K)','+|-',0.15),
    ),
    shortDesc = """Cb-Cb BENSON""",
    longDesc = 
"""

""",
)

entry(
    index = 58,
    label = "Cb-F",
    group = 
"""
1 * Cb  u0 {2,S}
2   F1s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.13964,7.28758,8.03688,8.39453,9.21743,9.8138,9.89767],'cal/(mol*K)'),
        H298 = (-43.9283,'kcal/mol'),
        S298 = (15.7765,'cal/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""
Derived from C6H5F in halogens thermo library
""",
)

entry(
    index = 59,
    label = "Cb-Cl",
    group = 
"""
1 * Cb   u0 {2,S}
2   Cl1s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.20014,7.75414,8.28478,8.72226,9.63755,10.1767,9.9642],'cal/(mol*K)'),
        H298 = (-4.0239,'kcal/mol'),
        S298 = (18.5799,'cal/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""
Derived from C6H5Cl in halogens thermo library
""",
)

entry(
    index = 60,
    label = "Cb-Br",
    group = 
"""
1 * Cb   u0 {2,S}
2   Br1s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.62376,8.17957,8.62571,8.95022,9.71948,10.1996,10.0043],'cal/(mol*K)'),
        H298 = (8.50001,'kcal/mol'),
        S298 = (21.3969,'cal/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""
Derived from C6H5Br in halogens thermo library
""",
)

entry(
    index = 61,
    label = "Cb-I",
    group = 
"""
1 * Cb  u0 {2,S}
2   I1s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8,8.9,9.6,9.9,10.3,10.5,10.7],'cal/(mol*K)'),
        H298 = (24,'kcal/mol'),
        S298 = (23.7,'cal/(mol*K)'),
    ),
    shortDesc = """Cb-I BENSON""",
    longDesc = 
"""
Thermochemical Kinetics 2nd Ed., by Sidney Benson (Table A4, p.281)
Cpdata at 1500K was not in the book, Cpdata at 1500K = Cpdata at 1000K + 0.2
""",
)

entry(
    index = 62,
    label = "Cb-N3s",
    group = 
"""
1 * Cb  u0 {2,S}
2   N3s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.95,5.21,5.94,6.32,6.53,6.56,6.635],'cal/(mol*K)'),
        H298 = (-0.5,'kcal/mol'),
        S298 = (-9.69,'cal/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 63,
    label = "Ct",
    group = 
"""
1 * Ct u0
""",
    thermo = 'Ct-CtCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 64,
    label = "CtBrC",
    group = 
"""
1 * Ct u0 {2,S} {3,T}
2   Br u0 {1,S}
3   C  u0 {1,T}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.53226,7.72112,7.96603,8.05651,8.51245,8.85766,9.00327],'cal/(mol*K)','+|-',[0.120214,0.138162,0.141235,0.139163,0.120673,0.10408,0.0848226]),
        H298 = (40.6636,'kcal/mol','+|-',0.428554),
        S298 = (36.0814,'cal/(mol*K)','+|-',0.281301),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library    | Number of Species
CHOBr_G4   |         39
CHOFBr_G4  |         2
CHOClBr_G4 |         2
""",
)

entry(
    index = 65,
    label = "CtCCl",
    group = 
"""
1 * Ct u0 {2,S} {3,T}
2   Cl u0 {1,S}
3   C  u0 {1,T}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.38542,7.64583,7.94497,8.07341,8.56204,8.93336,8.98031],'cal/(mol*K)','+|-',[0.0945988,0.108722,0.111141,0.10951,0.09496,0.0819021,0.0667484]),
        H298 = (30.1852,'kcal/mol','+|-',0.337237),
        S298 = (34.0295,'cal/(mol*K)','+|-',0.221361),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library    | Number of Species
CHOCl_G4   |         40
CHOFCl_G4  |         2
CHOClBr_G4 |         32
""",
)

entry(
    index = 66,
    label = "CtCF",
    group = 
"""
1 * Ct u0 {2,S} {3,T}
2   F  u0 {1,S}
3   C  u0 {1,T}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOF_G4     |         42
CHOFCl_G4   |         34
CHOFClBr_G4 |         12
CHOFBr_G4   |         51
""",
)

entry(
    index = 67,
    label = "Ct-CtH",
    group = 
"""
1 * Ct u0 {2,T} {3,S}
2   Ct u0 {1,T}
3   H  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Ct-H STEIN and FAHR; J. PHYS. CHEM. 1985, 89, 17, 3714""",
    longDesc = 
"""

""",
)

entry(
    index = 68,
    label = "Ct-StH",
    group = 
"""
1 * Ct             u0 {2,T} {3,S}
2   [S4t,S6t,S6td] u0 {1,T}
3   H              u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([-0.3,1.12,1.92,2.33,1.71,1.44,2.24],'cal/(mol*K)'),
        H298 = (98.15,'kcal/mol'),
        S298 = (-9.54,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""

""",
)

entry(
    index = 69,
    label = "Ct-CtOs",
    group = 
"""
1 * Ct  u0 {2,T} {3,S}
2   Ct  u0 {1,T}
3   O2s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Ct-O MELIUS / hc#coh !!!WARNING! Cp1500 value taken as Cp1000""",
    longDesc = 
"""

""",
)

entry(
    index = 70,
    label = "Ct-CtS",
    group = 
"""
1 * Ct u0 {2,T} {3,S}
2   Ct u0 {1,T}
3   S  u0 {1,S}
""",
    thermo = 'Ct-CtS2',
    shortDesc = """CBS-QB3 GA 1D-HR Aaron Vandeputte 2010""",
    longDesc = 
"""

""",
)

entry(
    index = 71,
    label = "Ct-CtS2",
    group = 
"""
1 * Ct  u0 {2,T} {3,S}
2   Ct  u0 {1,T}
3   S2s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.7,3.47,2.94,2.87,4.56,5.68,4.73],'cal/(mol*K)'),
        H298 = (45.23,'kcal/mol'),
        S298 = (14.57,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""

""",
)

entry(
    index = 72,
    label = "Ct-CtS4",
    group = 
"""
1 * Ct                u0 {2,T} {3,S}
2   Ct                u0 {1,T}
3   [S4s,S4d,S4b,S4t] u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([2.19,2.04,1.74,1.81,3.72,4.89,4.48],'cal/(mol*K)'),
        H298 = (56.56,'kcal/mol'),
        S298 = (12.4,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""

""",
)

entry(
    index = 73,
    label = "Ct-CtS6",
    group = 
"""
1 * Ct                      u0 {2,T} {3,S}
2   Ct                      u0 {1,T}
3   [S6s,S6d,S6dd,S6t,S6td] u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.29,3.67,4,4.29,4.74,5.05,5.49],'cal/(mol*K)'),
        H298 = (27.63,'kcal/mol'),
        S298 = (6.32,'cal/(mol*K)'),
    ),
    shortDesc = """CBS-QB3 GA 1D-HR Aaron Vandeputte 2010""",
    longDesc = 
"""

""",
)

entry(
    index = 74,
    label = "Ct-CtC",
    group = 
"""
1 * Ct u0 {2,T} {3,S}
2   Ct u0 {1,T}
3   C  u0 {1,S}
""",
    thermo = 'Ct-CtCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 75,
    label = "Ct-CtCs",
    group = 
"""
1 * Ct u0 {2,T} {3,S}
2   Ct u0 {1,T}
3   Cs u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Ct-Cs STEIN and FAHR; J. PHYS. CHEM. 1985, 89, 17, 3714""",
    longDesc = 
"""

""",
)

entry(
    index = 76,
    label = "Ct-CtCds",
    group = 
"""
1 * Ct         u0 {2,T} {3,S}
2   Ct         u0 {1,T}
3   [Cd,CO,CS] u0 {1,S}
""",
    thermo = 'Ct-Ct(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 77,
    label = "Ct-Ct(Cds-O2d)",
    group = 
"""
1 * Ct  u0 {2,S} {3,T}
2   CO  u0 {1,S} {4,D}
3   Ct  u0 {1,T}
4   O2d u0 {2,D}
""",
    thermo = 'Ct-CtCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 78,
    label = "Ct-Ct(Cds-Cd)",
    group = 
"""
1 * Ct u0 {2,S} {3,T}
2   Cd u0 {1,S} {4,D}
3   Ct u0 {1,T}
4   C  u0 {2,D}
""",
    thermo = 'Ct-Ct(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 79,
    label = "Ct-Ct(Cds-Cds)",
    group = 
"""
1 * Ct u0 {2,S} {3,T}
2   Cd u0 {1,S} {4,D}
3   Ct u0 {1,T}
4   Cd u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([2.57,3.54,3.5,4.92,5.34,5.5,5.8],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (28.2,'kcal/mol','+|-',0.27),
        S298 = (6.43,'cal/(mol*K)','+|-',0.09),
    ),
    shortDesc = """Ct-Cd STEIN and FAHR; J. PHYS. CHEM. 1985, 89, 17, 3714""",
    longDesc = 
"""

""",
)

entry(
    index = 80,
    label = "Ct-Ct(Cds-Cdd)",
    group = 
"""
1 * Ct  u0 {2,S} {3,T}
2   Cd  u0 {1,S} {4,D}
3   Ct  u0 {1,T}
4   Cdd u0 {2,D}
""",
    thermo = 'Ct-Ct(Cds-Cdd-Cd)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 81,
    label = "Ct-Ct(Cds-Cdd-O2d)",
    group = 
"""
1 * Ct  u0 {2,T} {3,S}
2   Ct  u0 {1,T}
3   Cd  u0 {1,S} {4,D}
4   Cdd u0 {3,D} {5,D}
5   O2d u0 {4,D}
""",
    thermo = 'Ct-Ct(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 82,
    label = "Ct-Ct(Cds-Cdd-S2d)",
    group = 
"""
1 * Ct  u0 {2,T} {3,S}
2   Ct  u0 {1,T}
3   Cd  u0 {1,S} {4,D}
4   Cdd u0 {3,D} {5,D}
5   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 83,
    label = "Ct-Ct(Cds-Cdd-Cd)",
    group = 
"""
1 * Ct  u0 {2,T} {3,S}
2   Ct  u0 {1,T}
3   Cd  u0 {1,S} {4,D}
4   Cdd u0 {3,D} {5,D}
5   C   u0 {4,D}
""",
    thermo = 'Ct-Ct(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 84,
    label = "Ct-CtCt",
    group = 
"""
1 * Ct u0 {2,T} {3,S}
2   Ct u0 {1,T}
3   Ct u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.54,4.06,4.4,4.64,5,5.23,5.57],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (25.6,'kcal/mol','+|-',0.27),
        S298 = (5.88,'cal/(mol*K)','+|-',0.09),
    ),
    shortDesc = """Ct-Ct STEIN and FAHR; J. PHYS. CHEM. 1985, 89, 17, 3714""",
    longDesc = 
"""

""",
)

entry(
    index = 85,
    label = "Ct-Ct(CtN3t)",
    group = 
"""
1 * Ct  u0 {2,S} {3,T}
2   Ct  u0 {1,S} {4,T}
3   Ct  u0 {1,T}
4   N3t u0 {2,T}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.857,10.9262,11.6704,12.2046,13.2063,13.9016,14.7678],'cal/(mol*K)','+|-',[0.340749,0.363915,0.375311,0.372646,0.362495,0.355991,0.363243]),
        H298 = (61.4436,'kcal/mol','+|-',1.36577),
        S298 = (34.96,'cal/(mol*K)','+|-',1.0464),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHON_G4  |         4
BurcatNS |         1
""",
)

entry(
    index = 86,
    label = "Ct-CtCb",
    group = 
"""
1 * Ct u0 {2,T} {3,S}
2   Ct u0 {1,T}
3   Cb u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([2.57,3.54,4.5,4.92,5.34,5.5,5.8],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (24.67,'kcal/mol','+|-',0.27),
        S298 = (6.43,'cal/(mol*K)','+|-',0.09),
    ),
    shortDesc = """Ct-Cb STEIN and FAHR; J. PHYS. CHEM. 1985, 89, 17, 3714""",
    longDesc = 
"""

""",
)

entry(
    index = 87,
    label = "Ct-HN",
    group = 
"""
1 * Ct         u0 {2,S} {3,T}
2   H          u0 {1,S}
3   [N3t,N5tc] u0 {1,T}
""",
    thermo = None,
    shortDesc = """Derived from nitrogen species in RMG thermo libraries""",
    longDesc = 
"""

""",
)

entry(
    index = 88,
    label = "Ct-HN5tc",
    group = 
"""
1 * Ct   u0 {2,S} {3,T}
2   H    u0 {1,S}
3   N5tc u0 {1,T}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([13.1965,13.9686,14.707,15.3897,16.5138,17.3711,18.7112],'cal/(mol*K)','+|-',[0.9602,1.02548,1.05759,1.05008,1.02148,1.00315,1.02359]),
        H298 = (39.9057,'kcal/mol','+|-',3.84862),
        S298 = (59.5485,'cal/(mol*K)','+|-',2.94866),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         1
""",
)

entry(
    index = 89,
    label = "Ct-NtN",
    group = 
"""
1 * Ct u0 {2,S} {3,T}
2   N  u0 {1,S}
3   N  u0 {1,T}
""",
    thermo = None,
    shortDesc = """Derived from nitrogen species in RMG thermo libraries""",
    longDesc = 
"""

""",
)

entry(
    index = 90,
    label = "Ct-N5tcN",
    group = 
"""
1 * Ct   u0 {2,T} {3,S}
2   N5tc u0 {1,T}
3   N    u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([11.5403,12.2978,12.6583,12.8343,13.1094,13.3187,13.9193],'cal/(mol*K)','+|-',[0.455745,0.486728,0.50197,0.498405,0.484829,0.47613,0.48583]),
        H298 = (57.745,'kcal/mol','+|-',1.82669),
        S298 = (39.8157,'cal/(mol*K)','+|-',1.39954),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         8
""",
)

entry(
    index = 91,
    label = "Ct-N3tN",
    group = 
"""
1 * Ct  u0 {2,T} {3,S}
2   N3t u0 {1,T}
3   N   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.65905,8.09076,8.17479,8.12017,8.08326,8.11664,8.58483],'cal/(mol*K)','+|-',[0.406998,0.434668,0.448279,0.445096,0.432972,0.425203,0.433866]),
        H298 = (42.4695,'kcal/mol','+|-',1.63131),
        S298 = (31.6607,'cal/(mol*K)','+|-',1.24985),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         14
""",
)

entry(
    index = 92,
    label = "Ct-N3tN3s",
    group = 
"""
1 * Ct  u0 {2,T} {3,S}
2   N3t u0 {1,T}
3   N3s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.05491,8.44107,8.47383,8.35466,8.22932,8.18992,8.52721],'cal/(mol*K)','+|-',[0.356205,0.380422,0.392335,0.389549,0.378938,0.372138,0.37972]),
        H298 = (42.0996,'kcal/mol','+|-',1.42772),
        S298 = (32.0261,'cal/(mol*K)','+|-',1.09387),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHON_G4  |         23
BurcatNS |         1
""",
)

entry(
    index = 93,
    label = "Ct-NtO",
    group = 
"""
1 * Ct u0 {2,S} {3,T}
2   O  u0 {1,S}
3   N  u0 {1,T}
""",
    thermo = None,
    shortDesc = """Derived from nitrogen species in RMG thermo libraries""",
    longDesc = 
"""

""",
)

entry(
    index = 94,
    label = "Ct-N3tO",
    group = 
"""
1 * Ct  u0 {2,S} {3,T}
2   O   u0 {1,S}
3   N3t u0 {1,T}
""",
    thermo = None,
    shortDesc = """Derived from nitrogen species in RMG thermo libraries""",
    longDesc = 
"""

""",
)

entry(
    index = 95,
    label = "Ct-N3tOs",
    group = 
"""
1 * Ct  u0 {2,T} {3,S}
2   N3t u0 {1,T}
3   O2s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.9877,6.39508,6.84559,7.25213,7.83897,8.3159,9.09594],'cal/(mol*K)','+|-',[0.215494,0.230144,0.237352,0.235666,0.229247,0.225133,0.22972]),
        H298 = (37.5554,'kcal/mol','+|-',0.863731),
        S298 = (30.5008,'cal/(mol*K)','+|-',0.661758),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         18
""",
)

entry(
    index = 96,
    label = "Ct-N5tcO",
    group = 
"""
1 * Ct   u0 {2,S} {3,T}
2   O    u0 {1,S}
3   N5tc u0 {1,T}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.99646,10.8243,11.5506,12.1907,13.1875,13.8352,14.8195],'cal/(mol*K)','+|-',[0.482528,0.515332,0.53147,0.527695,0.513321,0.504111,0.514381]),
        H298 = (54.2138,'kcal/mol','+|-',1.93404),
        S298 = (38.0611,'cal/(mol*K)','+|-',1.48179),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         4
""",
)

entry(
    index = 97,
    label = "Ct-NtC",
    group = 
"""
1 * Ct u0 {2,T} {3,S}
2   N  u0 {1,T}
3   C  u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 98,
    label = "Ct-N5tcC",
    group = 
"""
1 * Ct   u0 {2,T} {3,S}
2   N5tc u0 {1,T}
3   C    u0 {1,S}
""",
    thermo = None,
    shortDesc = """Derived from nitrogen species in RMG thermo libraries""",
    longDesc = 
"""

""",
)

entry(
    index = 99,
    label = "Ct-(N5tcO)C",
    group = 
"""
1 * Ct   u0 {2,T} {3,S}
2   N5tc u0 {1,T} {4,S}
3   C    u0 {1,S}
4   O0sc u0 {2,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.75633,10.7132,11.4431,12.2024,13.3001,14.0526,15.1112],'cal/(mol*K)','+|-',[0.3266,0.348803,0.359726,0.357172,0.347443,0.341208,0.34816]),
        H298 = (43.1883,'kcal/mol','+|-',1.30906),
        S298 = (37.3977,'cal/(mol*K)','+|-',1.00295),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         9
""",
)

entry(
    index = 100,
    label = "Ct-N3tC",
    group = 
"""
1 * Ct  u0 {2,T} {3,S}
2   N3t u0 {1,T}
3   C   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.5029,6.68177,7.25713,7.89206,8.47622,8.63402,10.4499],'cal/(mol*K)','+|-',[0.9602,1.02548,1.05759,1.05008,1.02148,1.00315,1.02359]),
        H298 = (37.562,'kcal/mol','+|-',3.84862),
        S298 = (31.9401,'cal/(mol*K)','+|-',2.94866),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
BurcatNS |         1
""",
)

entry(
    index = 101,
    label = "Ct-N3tCb",
    group = 
"""
1 * Ct  u0 {2,T} {3,S}
2   N3t u0 {1,T}
3   Cb  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kcal/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 102,
    label = "Ct-N3tCO",
    group = 
"""
1 * Ct  u0 {2,S} {3,T}
2   CO  u0 {1,S} {4,D}
3   N3t u0 {1,T}
4   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.74108,7.51754,7.99046,8.33222,8.84385,9.20133,10.2819],'cal/(mol*K)','+|-',[0.342801,0.366106,0.377571,0.37489,0.364678,0.358135,0.365431]),
        H298 = (39.8628,'kcal/mol','+|-',1.374),
        S298 = (29.7492,'cal/(mol*K)','+|-',1.0527),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHON_G4  |         4
BurcatNS |         1
""",
)

entry(
    index = 103,
    label = "Ct-N3tCt",
    group = 
"""
1 * Ct  u0 {2,T} {3,S}
2   N3t u0 {1,T}
3   Ct  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.63948,7.283,7.74146,8.10138,8.68031,9.09192,9.65702],'cal/(mol*K)','+|-',[0.43469,0.464242,0.47878,0.475379,0.462431,0.454133,0.463385]),
        H298 = (36.2723,'kcal/mol','+|-',1.7423),
        S298 = (28.9487,'cal/(mol*K)','+|-',1.33488),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         2
""",
)

entry(
    index = 104,
    label = "Ct-CtCt(N3t)",
    group = 
"""
1 * Ct  u0 {2,S} {3,T}
2   Ct  u0 {1,S} {4,T}
3   N3t u0 {1,T}
4   Ct  u0 {2,T}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 105,
    label = "Ct-N3tCd",
    group = 
"""
1 * Ct  u0 {2,T} {3,S}
2   N3t u0 {1,T}
3   Cd  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.66493,6.2676,6.70608,7.35118,8.18292,8.77454,9.47724],'cal/(mol*K)','+|-',[0.21319,0.227684,0.234814,0.233146,0.226796,0.222726,0.227264]),
        H298 = (32.2803,'kcal/mol','+|-',0.854497),
        S298 = (30.9581,'cal/(mol*K)','+|-',0.654684),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHON_G4  |         17
BurcatNS |         2
CHN      |         1
""",
)

entry(
    index = 106,
    label = "Ct-N3tCs",
    group = 
"""
1 * Ct  u0 {2,T} {3,S}
2   N3t u0 {1,T}
3   Cs  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.5919,6.9728,7.37354,7.75887,8.36458,8.78021,9.55871],'cal/(mol*K)','+|-',[0.10658,0.113826,0.117391,0.116557,0.113382,0.111348,0.113616]),
        H298 = (31.6428,'kcal/mol','+|-',0.427189),
        S298 = (29.9778,'cal/(mol*K)','+|-',0.327296),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library                 | Number of Species
CHON_G4                 |         23
CBS_QB3_1dHR            |         1
thermo_DFT_CCSDTF12_BAC |         2
CHN                     |         17
CHON                    |         6
""",
)

entry(
    index = 107,
    label = "Ct-N3tCs(HHH)",
    group = 
"""
1   Cs  u0 {2,S} {3,S} {4,S} {5,S}
2 * Ct  u0 {1,S} {6,T}
3   H   u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   N3t u0 {2,T}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'cal/(mol*K)'),
        H298 = (0,'kcal/mol'),
        S298 = (0,'cal/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 108,
    label = "Ct-CtN",
    group = 
"""
1 * Ct u0 {2,T} {3,S}
2   Ct u0 {1,T}
3   N  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.29661,4.69133,4.58823,4.25036,3.79242,3.60888,3.87325],'cal/(mol*K)','+|-',[1.09628,1.17081,1.20748,1.1989,1.16624,1.14532,1.16865]),
        H298 = (35.9275,'kcal/mol','+|-',4.39405),
        S298 = (8.19713,'cal/(mol*K)','+|-',3.36656),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         1
""",
)

entry(
    index = 109,
    label = "Ct-CtN3sd",
    group = 
"""
1 * Ct        u0 {2,T} {3,S}
2   Ct        u0 {1,T}
3   [N3s,N3d] u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.07117,5.3314,5.20918,4.93063,4.65766,4.48625,4.53719],'cal/(mol*K)','+|-',[0.329224,0.351607,0.362617,0.360042,0.350235,0.343951,0.350958]),
        H298 = (39.3588,'kcal/mol','+|-',1.31958),
        S298 = (8.12398,'cal/(mol*K)','+|-',1.01101),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         58
""",
)

entry(
    index = 110,
    label = "Ct-CtN5sdtc",
    group = 
"""
1 * Ct               u0 {2,T} {3,S}
2   Ct               u0 {1,T}
3   [N5sc,N5dc,N5tc] u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.54752,5.633,5.41461,4.92422,4.55648,4.31883,4.42388],'cal/(mol*K)','+|-',[0.65221,0.696549,0.718362,0.713261,0.693832,0.681382,0.695265]),
        H298 = (49.12,'kcal/mol','+|-',2.61415),
        S298 = (8.50797,'cal/(mol*K)','+|-',2.00286),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         3
""",
)

entry(
    index = 111,
    label = "Ct-CtNOO",
    group = 
"""
1   N5dc u0 {2,S} {3,D} {4,S}
2 * Ct   u0 {1,S} {5,T}
3   O2d  u0 {1,D}
4   O0sc u0 {1,S}
5   Ct   u0 {2,T}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.40428,5.62268,5.50652,5.23364,4.95924,4.75869,4.59335],'cal/(mol*K)','+|-',[0.56617,0.604661,0.623596,0.619168,0.602302,0.591495,0.603545]),
        H298 = (45.3136,'kcal/mol','+|-',2.26929),
        S298 = (6.23492,'cal/(mol*K)','+|-',1.73864),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library                 | Number of Species
BurcatNS                |         1
thermo_DFT_CCSDTF12_BAC |         1
""",
)

entry(
    index = 112,
    label = "Cdd",
    group = 
"""
1 * Cdd u0
""",
    thermo = 'Cdd-CdsCds',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 113,
    label = "Cdd-OdOd",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   O2d u0 {1,D}
3   O2d u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """CHEMKIN DATABASE: S(group) = S(CO2) + Rln(2)""",
    longDesc = 
"""

""",
)

entry(
    index = 114,
    label = "Cdd-OdSd",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   O2d u0 {1,D}
3   S2d u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.81,10.84,11.62,12.18,12.99,13.52,14.14],'cal/(mol*K)'),
        H298 = (-34.78,'kcal/mol'),
        S298 = (55.34,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""

""",
)

entry(
    index = 115,
    label = "Cdd-SdSd",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   S   u0 {1,D}
3   S   u0 {1,D}
""",
    thermo = 'Cdd-S2dS2d',
    shortDesc = """CBS-QB3 GA 1D-HR Aaron Vandeputte 2009""",
    longDesc = 
"""

""",
)

entry(
    index = 116,
    label = "Cdd-S2dS2d",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   S2d u0 {1,D}
3   S2d u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([10.91,11.83,12.49,12.98,13.63,14.01,14.47],'cal/(mol*K)'),
        H298 = (24.5,'kcal/mol'),
        S298 = (58.24,'cal/(mol*K)'),
    ),
    shortDesc = """CBS-QB3 GA 1D-HR Aaron Vandeputte 2009""",
    longDesc = 
"""

""",
)

entry(
    index = 117,
    label = "Cdd-S4S4",
    group = 
"""
1 * Cdd        u0 {2,D} {3,D}
2   [S4d,S4dd] u0 {1,D}
3   [S4d,S4dd] u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([12.44,10.6,10.48,10.8,11.92,12.75,12.7],'cal/(mol*K)'),
        H298 = (54.12,'kcal/mol'),
        S298 = (66.19,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""

""",
)

entry(
    index = 118,
    label = "Cdd-S2S4",
    group = 
"""
1 * Cdd        u0 {2,D} {3,D}
2   S2d        u0 {1,D}
3   [S4d,S4dd] u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([11.88,11.09,11.15,11.47,12.26,12.78,13.01],'cal/(mol*K)'),
        H298 = (60.89,'kcal/mol'),
        S298 = (65.93,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""

""",
)

entry(
    index = 119,
    label = "Cdd-CdOd",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   C   u0 {1,D}
3   O2d u0 {1,D}
""",
    thermo = 'Cdd-CdsOd',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 120,
    label = "Cdd-CdsOd",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   Cd  u0 {1,D}
3   O2d u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """O=C*=C< currently treat the adjacent C as Ck""",
    longDesc = 
"""

""",
)

entry(
    index = 121,
    label = "Cdd-(CdN)Od",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   Cd  u0 {1,D} {4,S}
3   O2d u0 {1,D}
4   N   u0 {2,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.63397,7.16103,7.06937,7.21312,7.51372,7.77288,8.24234],'cal/(mol*K)','+|-',[1.00013,1.06813,1.10158,1.09375,1.06396,1.04487,1.06616]),
        H298 = (6.04835,'kcal/mol','+|-',4.00868),
        S298 = (35.9268,'cal/(mol*K)','+|-',3.0713),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         8
""",
)

entry(
    index = 122,
    label = "Cdd-CddOd",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   Cdd u0 {1,D}
3   O2d u0 {1,D}
""",
    thermo = 'Cdd-(Cdd-Cd)O2d',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 123,
    label = "Cdd-(Cdd-O2d)O2d",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   Cdd u0 {1,D} {4,D}
3   O2d u0 {1,D}
4   O2d u0 {2,D}
""",
    thermo = 'Cdd-CdsOd',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 124,
    label = "Cdd-(Cdd-Cd)O2d",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   Cdd u0 {1,D} {4,D}
3   O2d u0 {1,D}
4   C   u0 {2,D}
""",
    thermo = 'Cdd-CdsOd',
    shortDesc = """O=C*=C= currently not defined. Assigned same value as Ca""",
    longDesc = 
"""

""",
)

entry(
    index = 125,
    label = "Cdd-(CddN)Od",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   Cdd u0 {1,D} {4,D}
3   O2d u0 {1,D}
4   N   u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.15393,6.75455,7.28196,7.70947,8.30474,8.74972,9.39024],'cal/(mol*K)','+|-',[0.498192,0.532061,0.548723,0.544827,0.529986,0.520476,0.53108]),
        H298 = (-7.71295,'kcal/mol','+|-',1.99682),
        S298 = (31.9664,'cal/(mol*K)','+|-',1.52989),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         5
""",
)

entry(
    index = 126,
    label = "Cdd-CdSd",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   C   u0 {1,D}
3   S   u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.88,8.48,8.8,8.99,9.23,9.37,9.58],'cal/(mol*K)'),
        H298 = (40.33,'kcal/mol'),
        S298 = (34.24,'cal/(mol*K)'),
    ),
    shortDesc = """CBS-QB3 GA 1D-HR Aaron Vandeputte 2010""",
    longDesc = 
"""

""",
)

entry(
    index = 127,
    label = "Cdd-CdsSd",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   Cd  u0 {1,D}
3   S   u0 {1,D}
""",
    thermo = 'Cdd-CdsS6d',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 128,
    label = "Cdd-CdsS2d",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   Cd  u0 {1,D}
3   S2d u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.88,8.48,8.8,8.99,9.23,9.37,9.58],'cal/(mol*K)'),
        H298 = (40.33,'kcal/mol'),
        S298 = (34.24,'cal/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 129,
    label = "Cdd-CdsS4d",
    group = 
"""
1 * Cdd        u0 {2,D} {3,D}
2   Cd         u0 {1,D}
3   [S4d,S4dd] u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.88,8.48,8.8,8.99,9.23,9.37,9.58],'cal/(mol*K)'),
        H298 = (40.33,'kcal/mol'),
        S298 = (34.24,'cal/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 130,
    label = "Cdd-CdsS6d",
    group = 
"""
1 * Cdd                   u0 {2,D} {3,D}
2   Cd                    u0 {1,D}
3   [S6d,S6dd,S6ddd,S6td] u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.18,7.47,7.69,7.97,8.76,9.31,9.81],'cal/(mol*K)'),
        H298 = (45.67,'kcal/mol'),
        S298 = (34.13,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""

""",
)

entry(
    index = 131,
    label = "Cdd-CddSd",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   Cdd u0 {1,D}
3   S2d u0 {1,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 132,
    label = "Cdd-(Cdd-S2d)S2d",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   Cdd u0 {1,D} {4,D}
3   S2d u0 {1,D}
4   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 133,
    label = "Cdd-(Cdd-Cd)S2d",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   Cdd u0 {1,D} {4,D}
3   S2d u0 {1,D}
4   C   u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 134,
    label = "Cdd-CdCd",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   C   u0 {1,D}
3   C   u0 {1,D}
""",
    thermo = 'Cdd-CdsCds',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 135,
    label = "Cdd-CddCdd",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   Cdd u0 {1,D}
3   Cdd u0 {1,D}
""",
    thermo = 'Cdd-(Cdd-Cd)(Cdd-Cd)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 136,
    label = "Cdd-(Cdd-O2d)(Cdd-O2d)",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   Cdd u0 {1,D} {4,D}
3   Cdd u0 {1,D} {5,D}
4   O2d u0 {2,D}
5   O2d u0 {3,D}
""",
    thermo = 'Cdd-CdsCds',
    shortDesc = """O=C=C*=C=O, currently not defined. Assigned same value as Ca""",
    longDesc = 
"""

""",
)

entry(
    index = 137,
    label = "Cdd-(Cdd-S2d)(Cdd-S2d)",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   Cdd u0 {1,D} {4,D}
3   Cdd u0 {1,D} {5,D}
4   S2d u0 {2,D}
5   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 138,
    label = "Cdd-(Cdd-O2d)(Cdd-Cd)",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   Cdd u0 {1,D} {4,D}
3   Cdd u0 {1,D} {5,D}
4   O2d u0 {2,D}
5   C   u0 {3,D}
""",
    thermo = 'Cdd-(Cdd-O2d)Cds',
    shortDesc = """O=C=C*=C=C, currently not defined. Assigned same value as Ca""",
    longDesc = 
"""

""",
)

entry(
    index = 139,
    label = "Cdd-(Cdd-S2d)(Cdd-Cd)",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   Cdd u0 {1,D} {4,D}
3   Cdd u0 {1,D} {5,D}
4   S2d u0 {2,D}
5   C   u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 140,
    label = "Cdd-(Cdd-Cd)(Cdd-Cd)",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   Cdd u0 {1,D} {4,D}
3   Cdd u0 {1,D} {5,D}
4   C   u0 {2,D}
5   C   u0 {3,D}
""",
    thermo = 'Cdd-CdsCds',
    shortDesc = """C=C=C*=C=C, currently not defined. Assigned same value as Ca""",
    longDesc = 
"""

""",
)

entry(
    index = 141,
    label = "Cdd-(Cdd-O2d)(Cdd-N3d)",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   Cdd u0 {1,D} {4,D}
3   Cdd u0 {1,D} {5,D}
4   O2d u0 {2,D}
5   N3d u0 {3,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([10.2393,11.1513,11.9342,12.5807,13.4642,14.0388,14.7435],'cal/(mol*K)','+|-',[0.990302,1.05763,1.09075,1.083,1.0535,1.0346,1.05568]),
        H298 = (-1.03148,'kcal/mol','+|-',3.96927),
        S298 = (38.4995,'cal/(mol*K)','+|-',3.04111),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         1
""",
)

entry(
    index = 142,
    label = "Cdd-CddCds",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   Cdd u0 {1,D}
3   Cd  u0 {1,D}
""",
    thermo = 'Cdd-(Cdd-Cd)(Cdd-Cd)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 143,
    label = "Cdd-(Cdd-O2d)Cds",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   Cdd u0 {1,D} {4,D}
3   Cd  u0 {1,D}
4   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.42494,10.6543,11.5241,12.2864,13.2975,13.9637,14.7777],'cal/(mol*K)','+|-',[0.278815,0.320441,0.327569,0.322763,0.279879,0.241393,0.19673]),
        H298 = (17.6301,'kcal/mol','+|-',0.993951),
        S298 = (38.2736,'cal/(mol*K)','+|-',0.677437),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library    | Number of Species
CHOBr_G4   |         2
CHOCl_G4   |         2
CHOF_G4    |         2
CHOFCl_G4  |         1
CHOFBr_G4  |         1
CHOClBr_G4 |         1
""",
)

entry(
    index = 144,
    label = "Cdd-(Cdd-S2d)Cds",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   Cdd u0 {1,D} {4,D}
3   Cd  u0 {1,D}
4   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 145,
    label = "Cdd-(Cdd-Cd)Cds",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   Cdd u0 {1,D} {4,D}
3   Cd  u0 {1,D}
4   C   u0 {2,D}
""",
    thermo = 'Cdd-CdsCds',
    shortDesc = """C=C=C*=C<, currently not defined. Assigned same value as Ca""",
    longDesc = 
"""

""",
)

entry(
    index = 146,
    label = "Cdd-CdsCds",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   Cd  u0 {1,D}
3   Cd  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.9,4.4,4.7,5,5.3,5.5,5.7],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (34.2,'kcal/mol','+|-',0.2),
        S298 = (6,'cal/(mol*K)','+|-',0.1),
    ),
    shortDesc = """Benson's Ca""",
    longDesc = 
"""

""",
)

entry(
    index = 147,
    label = "Cdd-NC",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   N   u0 {1,D}
3   C   u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.93758,5.45406,6.0663,6.84333,8.09172,8.96042,9.82153],'cal/(mol*K)','+|-',[0.624368,0.666815,0.687697,0.682813,0.664214,0.652295,0.665585]),
        H298 = (58.2131,'kcal/mol','+|-',2.50255),
        S298 = (13.6662,'cal/(mol*K)','+|-',1.91736),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         7
""",
)

entry(
    index = 148,
    label = "Cdd-N3dCd",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   N3d u0 {1,D}
3   Cd  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.12203,5.6722,6.15787,6.79214,7.78663,8.54349,9.43468],'cal/(mol*K)','+|-',[0.601617,0.642517,0.662638,0.657933,0.640011,0.628527,0.641332]),
        H298 = (46.9852,'kcal/mol','+|-',2.41137),
        S298 = (12.4036,'cal/(mol*K)','+|-',1.8475),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         44
""",
)

entry(
    index = 149,
    label = "Cdd-N3dCdd",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   N3d u0 {1,D}
3   Cdd u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.98984,5.79563,6.34738,6.99959,7.99243,8.77693,9.83483],'cal/(mol*K)','+|-',[0.727161,0.776597,0.800916,0.795229,0.773567,0.759687,0.775165]),
        H298 = (46.1466,'kcal/mol','+|-',2.91457),
        S298 = (11.3784,'cal/(mol*K)','+|-',2.23303),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         8
""",
)

entry(
    index = 150,
    label = "Cdd-(N3dH)Cdd",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   N3d u0 {1,D} {4,S}
3   Cdd u0 {1,D}
4   H   u0 {2,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.02705,5.81321,6.40931,7.12172,8.15255,8.88491,9.71869],'cal/(mol*K)','+|-',[0.676993,0.723018,0.74566,0.740365,0.720198,0.707275,0.721685]),
        H298 = (45.3269,'kcal/mol','+|-',2.71349),
        S298 = (12.169,'cal/(mol*K)','+|-',2.07897),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         12
""",
)

entry(
    index = 151,
    label = "Cdd-NO",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   N   u0 {1,D}
3   O2d u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.20142,9.97884,10.5705,11.2468,12.264,13.0021,13.8062],'cal/(mol*K)','+|-',[1.08311,1.15674,1.19297,1.18449,1.15223,1.13155,1.15461]),
        H298 = (22.718,'kcal/mol','+|-',4.34125),
        S298 = (40.0573,'cal/(mol*K)','+|-',3.3261),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         1
""",
)

entry(
    index = 152,
    label = "Cdd-N3dOd",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   N3d u0 {1,D}
3   O2d u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.29392,8.14273,8.84411,9.62521,10.7677,11.6742,12.8657],'cal/(mol*K)','+|-',[0.579291,0.618673,0.638047,0.633516,0.61626,0.605202,0.617532]),
        H298 = (-14.5628,'kcal/mol','+|-',2.32188),
        S298 = (38.0639,'cal/(mol*K)','+|-',1.77894),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         20
CHON    |         3
""",
)

entry(
    index = 153,
    label = "Cdd-NN",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   N   u0 {1,D}
3   N   u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.7449,7.86757,8.84375,10.0045,11.6543,12.87,14.1675],'cal/(mol*K)','+|-',[1.06676,1.13928,1.17496,1.16662,1.13484,1.11448,1.13718]),
        H298 = (75.9955,'kcal/mol','+|-',4.27573),
        S298 = (17.9896,'cal/(mol*K)','+|-',3.2759),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         7
""",
)

entry(
    index = 154,
    label = "Cdd-N3dN3d",
    group = 
"""
1 * Cdd u0 {2,D} {3,D}
2   N3d u0 {1,D}
3   N3d u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.99772,8.16411,9.01142,10.083,11.515,12.5823,13.8625],'cal/(mol*K)','+|-',[1.18319,1.26363,1.3032,1.29395,1.2587,1.23611,1.2613]),
        H298 = (48.2205,'kcal/mol','+|-',4.7424),
        S298 = (18.2957,'cal/(mol*K)','+|-',3.63345),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         24
""",
)

entry(
    index = 155,
    label = "Cds",
    group = 
"""
1 * [Cd,CO,CS] u0
""",
    thermo = 'Cds-CdsCsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 156,
    label = "COBrBrO",
    group = 
"""
1 * CO u0 {2,S} {3,S} {4,D}
2   Br u0 {1,S}
3   Br u0 {1,S}
4   O  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([14.8529,16.0692,16.7669,17.3279,18.1456,18.6561,19.2045],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-28.3191,'kcal/mol','+|-',2.95256),
        S298 = (75.0212,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHOBr_G4 |         1
""",
)

entry(
    index = 157,
    label = "CdBrBrC",
    group = 
"""
1 * Cd u0 {2,S} {3,S} {4,D}
2   Br u0 {1,S}
3   Br u0 {1,S}
4   C  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([11.8911,12.6269,13.379,13.8939,14.7842,15.3221,15.7807],'cal/(mol*K)','+|-',[0.119579,0.137432,0.140489,0.138428,0.120036,0.10353,0.0843745]),
        H298 = (18.387,'kcal/mol','+|-',0.42629),
        S298 = (47.5723,'cal/(mol*K)','+|-',0.279815),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOBr_G4    |         57
CHOFClBr_G4 |         2
CHOFBr_G4   |         18
CHOClBr_G4  |         10
""",
)

entry(
    index = 158,
    label = "CdBrBrCdd",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   Br  u0 {1,S}
3   Br  u0 {1,S}
4   Cdd u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([11.9502,12.7019,13.3155,13.7821,14.5777,15.0445,15.4525],'cal/(mol*K)','+|-',[0.155706,0.178952,0.182933,0.180249,0.1563,0.134808,0.109865]),
        H298 = (20.8555,'kcal/mol','+|-',0.555079),
        S298 = (48.4201,'cal/(mol*K)','+|-',0.364352),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOBr_G4    |         15
CHOFClBr_G4 |         1
CHOFBr_G4   |         6
CHOClBr_G4  |         4
""",
)

entry(
    index = 159,
    label = "Cd(Cdd-Od)BrBr",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Br  u0 {1,S}
4   Br  u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([18.3475,19.7117,20.8152,21.7194,23.0359,23.855,24.72],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (14.2101,'kcal/mol','+|-',2.95256),
        S298 = (81.2627,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHOBr_G4 |         1
""",
)

entry(
    index = 160,
    label = "COBrClO",
    group = 
"""
1 * CO u0 {2,S} {3,S} {4,D}
2   Cl u0 {1,S}
3   Br u0 {1,S}
4   O  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([14.251,15.7064,16.5387,17.1489,18.0358,18.5864,19.1679],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-40.5587,'kcal/mol','+|-',2.95256),
        S298 = (72.0675,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library    | Number of Species
CHOClBr_G4 |         1
""",
)

entry(
    index = 161,
    label = "CdBrCCl",
    group = 
"""
1 * Cd u0 {2,S} {3,S} {4,D}
2   Cl u0 {1,S}
3   Br u0 {1,S}
4   C  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([11.3503,12.2686,13.0575,13.6258,14.5833,15.1567,15.6962],'cal/(mol*K)','+|-',[0.0994227,0.114266,0.116808,0.115094,0.0998022,0.0860786,0.0701521]),
        H298 = (7.31075,'kcal/mol','+|-',0.354433),
        S298 = (45.2633,'cal/(mol*K)','+|-',0.232649),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOFClBr_G4 |         14
CHOClBr_G4  |         67
""",
)

entry(
    index = 162,
    label = "Cd(Cdd-O2d)ClBr",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cl  u0 {1,S}
4   Br  u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([18.0253,19.4989,20.6291,21.5582,22.9189,23.7744,24.6953],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (4.1688,'kcal/mol','+|-',2.95256),
        S298 = (78.3906,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library    | Number of Species
CHOClBr_G4 |         1
""",
)

entry(
    index = 163,
    label = "COClClO",
    group = 
"""
1 * CO u0 {2,S} {3,S} {4,D}
2   Cl u0 {1,S}
3   Cl u0 {1,S}
4   O  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([13.6488,15.2632,16.2837,16.9581,17.9245,18.5184,19.1285],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-53.0276,'kcal/mol','+|-',2.95256),
        S298 = (69.1372,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHOCl_G4 |         1
""",
)

entry(
    index = 164,
    label = "CdCClCl",
    group = 
"""
1 * Cd u0 {2,S} {3,S} {4,D}
2   Cl u0 {1,S}
3   Cl u0 {1,S}
4   C  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([11.0157,12.1575,13.0669,13.7142,14.6887,15.2685,15.8766],'cal/(mol*K)','+|-',[0.0976422,0.11222,0.114716,0.113033,0.098015,0.0845371,0.0688958]),
        H298 = (-5.3781,'kcal/mol','+|-',0.348086),
        S298 = (41.8561,'cal/(mol*K)','+|-',0.228482),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library    | Number of Species
CHOCl_G4   |         103
CHOFCl_G4  |         8
CHOClBr_G4 |         37
""",
)

entry(
    index = 165,
    label = "CdCddClCl",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   Cl  u0 {1,S}
3   Cl  u0 {1,S}
4   Cdd u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([11.1947,12.1765,12.8905,13.4327,14.3309,14.8579,15.3141],'cal/(mol*K)','+|-',[0.138643,0.159342,0.162886,0.160497,0.139172,0.120035,0.0978257]),
        H298 = (-1.73724,'kcal/mol','+|-',0.49425),
        S298 = (43.0609,'cal/(mol*K)','+|-',0.324424),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library    | Number of Species
CHOCl_G4   |         19
CHOFCl_G4  |         4
CHOClBr_G4 |         8
""",
)

entry(
    index = 166,
    label = "Cd(Cdd-Od)ClCl",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cl  u0 {1,S}
4   Cl  u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([17.635,19.2806,20.4418,21.3984,22.8049,23.6952,24.6658],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-6.13992,'kcal/mol','+|-',2.95256),
        S298 = (75.5269,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHOCl_G4 |         1
""",
)

entry(
    index = 167,
    label = "COBrFO",
    group = 
"""
1 * CO u0 {2,S} {3,S} {4,D}
2   F  u0 {1,S}
3   Br u0 {1,S}
4   O  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([12.9681,14.5466,15.6413,16.4113,17.5327,18.2426,19.0098],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-85.6424,'kcal/mol','+|-',2.95256),
        S298 = (69.0582,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOFBr_G4 |         1
""",
)

entry(
    index = 168,
    label = "CdBrCF",
    group = 
"""
1 * Cd u0 {2,S} {3,S} {4,D}
2   F  u0 {1,S}
3   Br u0 {1,S}
4   C  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([10.7294,11.8264,12.6955,13.3172,14.3,14.9049,15.5121],'cal/(mol*K)','+|-',[0.0932458,0.107167,0.109551,0.107944,0.0936018,0.0807307,0.0657937]),
        H298 = (-30.8812,'kcal/mol','+|-',0.332413),
        S298 = (42.6519,'cal/(mol*K)','+|-',0.218195),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOFBr_G4 |         97
""",
)

entry(
    index = 169,
    label = "Cd(Cdd-O2d)FBr",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   F   u0 {1,S}
4   Br  u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([17.3679,18.9651,20.1468,21.1261,22.5815,23.5206,24.5882],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-27.1033,'kcal/mol','+|-',2.95256),
        S298 = (75.8337,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOFBr_G4 |         1
""",
)

entry(
    index = 170,
    label = "COClFO",
    group = 
"""
1 * CO u0 {2,S} {3,S} {4,D}
2   F  u0 {1,S}
3   Cl u0 {1,S}
4   O  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([12.4165,14.0483,15.2955,16.1797,17.4077,18.1727,18.9666],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-98.531,'kcal/mol','+|-',2.95256),
        S298 = (66.1682,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOFCl_G4 |         1
""",
)

entry(
    index = 171,
    label = "CdCClF",
    group = 
"""
1 * Cd u0 {2,S} {3,S} {4,D}
2   F  u0 {1,S}
3   Cl u0 {1,S}
4   C  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([10.0601,11.2945,12.2608,12.9567,14.048,14.732,15.4604],'cal/(mol*K)','+|-',[0.0888591,0.102125,0.104397,0.102866,0.0891983,0.0769328,0.0626985]),
        H298 = (-43.1382,'kcal/mol','+|-',0.316775),
        S298 = (39.8406,'cal/(mol*K)','+|-',0.20793),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOFCl_G4   |         62
CHOFClBr_G4 |         45
""",
)

entry(
    index = 172,
    label = "Cd(Cdd-O2d)FCl",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   F   u0 {1,S}
4   Cl  u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([16.8933,18.6957,19.9263,20.9432,22.4563,23.4346,24.5487],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-38.4637,'kcal/mol','+|-',2.95256),
        S298 = (72.8961,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOFCl_G4 |         1
""",
)

entry(
    index = 173,
    label = "COFFO",
    group = 
"""
1 * CO u0 {2,S} {3,S} {4,D}
2   F  u0 {1,S}
3   F  u0 {1,S}
4   O  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHOF_G4 |         1
""",
)

entry(
    index = 174,
    label = "CdCFF",
    group = 
"""
1 * Cd u0 {2,S} {3,S} {4,D}
2   F  u0 {1,S}
3   F  u0 {1,S}
4   C  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOF_G4     |         102
CHOFCl_G4   |         31
CHOFClBr_G4 |         6
CHOFBr_G4   |         81
""",
)

entry(
    index = 175,
    label = "CdCddFF",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   F   u0 {1,S}
3   F   u0 {1,S}
4   Cdd u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.89331,10.2631,11.3446,12.1241,13.3367,14.1187,15.0479],'cal/(mol*K)','+|-',[0.120334,0.1383,0.141376,0.139302,0.120794,0.104183,0.0849071]),
        H298 = (-81.2229,'kcal/mol','+|-',0.428981),
        S298 = (37.6055,'cal/(mol*K)','+|-',0.281581),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOF_G4     |         19
CHOFCl_G4   |         8
CHOFClBr_G4 |         1
CHOFBr_G4   |         16
""",
)

entry(
    index = 176,
    label = "Cd(Cdd-Od)FF",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   F   u0 {1,S}
4   F   u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHOF_G4 |         1
""",
)

entry(
    index = 177,
    label = "COBrHO",
    group = 
"""
1 * CO u0 {2,S} {3,S} {4,D}
2   H  u0 {1,S}
3   Br u0 {1,S}
4   O  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([46.4592,51.6749,56.0267,59.6558,65.388,69.4959,75.2835],'J/(mol*K)'),
        H298 = (-129.956,'kJ/mol'),
        S298 = (271.223,'J/(mol*K)'),
    ),
    shortDesc = """G4 calc""",
    longDesc = 
"""
ODCBr in CHOBr_G4 thermo library
""",
)

entry(
    index = 178,
    label = "CdBrCH",
    group = 
"""
1 * Cd u0 {2,S} {3,S} {4,D}
2   H  u0 {1,S}
3   Br u0 {1,S}
4   C  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.40569,9.50824,10.4846,11.2566,12.4858,13.3329,14.5174],'cal/(mol*K)','+|-',[0.0885984,0.101826,0.104091,0.102564,0.0889366,0.0767071,0.0625145]),
        H298 = (10.1628,'kcal/mol','+|-',0.315846),
        S298 = (37.4419,'cal/(mol*K)','+|-',0.20732),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOBr_G4    |         79
CHOFClBr_G4 |         6
CHOFBr_G4   |         27
CHOClBr_G4  |         18
""",
)

entry(
    index = 179,
    label = "CdBrCddH",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   H   u0 {1,S}
3   Br  u0 {1,S}
4   Cdd u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.15331,9.17945,10.0774,10.8053,12.0717,12.9415,14.1155],'cal/(mol*K)','+|-',[0.137848,0.158428,0.161953,0.159577,0.138374,0.119347,0.0972649]),
        H298 = (12.5562,'kcal/mol','+|-',0.491417),
        S298 = (38.5178,'cal/(mol*K)','+|-',0.322564),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOBr_G4    |         17
CHOFClBr_G4 |         2
CHOFBr_G4   |         6
CHOClBr_G4  |         6
""",
)

entry(
    index = 180,
    label = "Cd(Cdd-Od)BrH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   H   u0 {1,S}
4   Br  u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([14.7956,16.6724,17.9852,19.0614,20.7366,21.9109,23.5205],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (2.25528,'kcal/mol','+|-',2.95256),
        S298 = (70.5685,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHOBr_G4 |         1
""",
)

entry(
    index = 181,
    label = "COClHO",
    group = 
"""
1 * CO u0 {2,S} {3,S} {4,D}
2   H  u0 {1,S}
3   Cl u0 {1,S}
4   O  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([10.7154,12.0444,13.1242,14.0367,15.4773,16.5071,17.9425],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-44.188,'kcal/mol','+|-',2.95256),
        S298 = (61.943,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHOCl_G4 |         1
""",
)

entry(
    index = 182,
    label = "CdCClH",
    group = 
"""
1 * Cd u0 {2,S} {3,S} {4,D}
2   H  u0 {1,S}
3   Cl u0 {1,S}
4   C  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.07436,9.30622,10.3225,11.1302,12.3659,13.2299,14.4837],'cal/(mol*K)','+|-',[0.0731681,0.0840918,0.0859625,0.0847013,0.0734474,0.0633478,0.051627]),
        H298 = (-1.55969,'kcal/mol','+|-',0.260838),
        S298 = (34.7433,'cal/(mol*K)','+|-',0.171213),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library    | Number of Species
CHOCl_G4   |         104
CHOFCl_G4  |         18
CHOClBr_G4 |         87
""",
)

entry(
    index = 183,
    label = "CdCddClH",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   H   u0 {1,S}
3   Cl  u0 {1,S}
4   Cdd u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.84129,8.93976,9.86814,10.6136,11.9064,12.8118,14.0543],'cal/(mol*K)','+|-',[0.124162,0.142698,0.145873,0.143733,0.124636,0.107497,0.0876077]),
        H298 = (1.44247,'kcal/mol','+|-',0.442626),
        S298 = (35.8666,'cal/(mol*K)','+|-',0.290538),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library    | Number of Species
CHOCl_G4   |         19
CHOFCl_G4  |         6
CHOClBr_G4 |         16
""",
)

entry(
    index = 184,
    label = "Cd(Cdd-Od)ClH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   H   u0 {1,S}
4   Cl  u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([14.3677,16.3272,17.7335,18.8545,20.5975,21.8169,23.478],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-7.571,'kcal/mol','+|-',2.95256),
        S298 = (67.7298,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHOCl_G4 |         1
""",
)

entry(
    index = 185,
    label = "COFHO",
    group = 
"""
1 * CO u0 {2,S} {3,S} {4,D}
2   H  u0 {1,S}
3   F  u0 {1,S}
4   O  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHOF_G4 |         1
""",
)

entry(
    index = 186,
    label = "CdCFH",
    group = 
"""
1 * Cd u0 {2,S} {3,S} {4,D}
2   H  u0 {1,S}
3   F  u0 {1,S}
4   C  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOF_G4     |         102
CHOFCl_G4   |         81
CHOFClBr_G4 |         36
CHOFBr_G4   |         141
""",
)

entry(
    index = 187,
    label = "CdCddFH",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   H   u0 {1,S}
3   F   u0 {1,S}
4   Cdd u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.85661,8.05208,9.11854,9.97185,11.4171,12.4518,13.931],'cal/(mol*K)','+|-',[0.105197,0.120902,0.123592,0.121779,0.105598,0.0910778,0.0742264]),
        H298 = (-34.9589,'kcal/mol','+|-',0.375018),
        S298 = (33.0088,'cal/(mol*K)','+|-',0.24616),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOF_G4     |         19
CHOFCl_G4   |         16
CHOFClBr_G4 |         5
CHOFBr_G4   |         22
""",
)

entry(
    index = 188,
    label = "Cd(Cdd-Od)FH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   H   u0 {1,S}
4   F   u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHOF_G4 |         1
""",
)

entry(
    index = 189,
    label = "COBrOO",
    group = 
"""
1 * CO u0 {2,S} {3,S} {4,D}
2   O  u0 {1,S}
3   Br u0 {1,S}
4   O  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([10.6292,12.883,14.4842,15.5841,16.7235,17.1257,16.7967],'cal/(mol*K)','+|-',[0.217606,0.250094,0.255657,0.251906,0.218437,0.1884,0.153542]),
        H298 = (-45.4297,'kcal/mol','+|-',0.775747),
        S298 = (36.1862,'cal/(mol*K)','+|-',0.509197),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOBr_G4    |         5
CHOFClBr_G4 |         1
CHOFBr_G4   |         6
CHOClBr_G4  |         3
""",
)

entry(
    index = 190,
    label = "CdBrCO",
    group = 
"""
1 * Cd u0 {2,S} {3,S} {4,D}
2   O  u0 {1,S}
3   Br u0 {1,S}
4   C  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.84488,8.45135,9.05382,9.58008,10.3948,10.8908,11.2435],'cal/(mol*K)','+|-',[0.0927652,0.106615,0.108986,0.107387,0.0931193,0.0803146,0.0654546]),
        H298 = (12.4381,'kcal/mol','+|-',0.3307),
        S298 = (18.6973,'cal/(mol*K)','+|-',0.21707),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOBr_G4    |         35
CHOFClBr_G4 |         6
CHOFBr_G4   |         43
CHOClBr_G4  |         27
""",
)

entry(
    index = 191,
    label = "CdBrCddO",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   O   u0 {1,S}
3   Br  u0 {1,S}
4   Cdd u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.68277,7.97037,8.51241,9.07043,9.97839,10.5028,10.8887],'cal/(mol*K)','+|-',[0.204332,0.234838,0.240062,0.23654,0.205112,0.176908,0.144176]),
        H298 = (17.3145,'kcal/mol','+|-',0.728427),
        S298 = (19.3282,'cal/(mol*K)','+|-',0.478137),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOBr_G4    |         6
CHOFClBr_G4 |         1
CHOFBr_G4   |         6
CHOClBr_G4  |         4
""",
)

entry(
    index = 192,
    label = "Cd(Cdd-Od)BrO",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   O   u0 {1,S}
4   Br  u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([13.7483,14.938,16.0867,17.1569,18.7784,19.7704,20.6771],'cal/(mol*K)','+|-',[0.585871,0.67334,0.688318,0.67822,0.588108,0.507238,0.413388]),
        H298 = (11.7236,'kcal/mol','+|-',2.08858),
        S298 = (51.8483,'cal/(mol*K)','+|-',1.37094),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHOBr_G4 |         2
""",
)

entry(
    index = 193,
    label = "COClOO",
    group = 
"""
1 * CO u0 {2,S} {3,S} {4,D}
2   O  u0 {1,S}
3   Cl u0 {1,S}
4   O  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([10.0863,12.4243,14.0742,15.1924,16.3723,16.7847,16.6471],'cal/(mol*K)','+|-',[0.264962,0.30452,0.311294,0.306727,0.265973,0.2294,0.186956]),
        H298 = (-57.4862,'kcal/mol','+|-',0.944567),
        S298 = (36.7364,'cal/(mol*K)','+|-',0.62001),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library    | Number of Species
CHOCl_G4   |         6
CHOFCl_G4  |         3
CHOClBr_G4 |         1
""",
)

entry(
    index = 194,
    label = "CdCClO",
    group = 
"""
1 * Cd u0 {2,S} {3,S} {4,D}
2   O  u0 {1,S}
3   Cl u0 {1,S}
4   C  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.35537,8.15155,8.78872,9.34639,10.2196,10.7592,11.1595],'cal/(mol*K)','+|-',[0.0942948,0.108373,0.110783,0.109158,0.0946548,0.0816389,0.0665339]),
        H298 = (0.403727,'kcal/mol','+|-',0.336153),
        S298 = (15.8607,'cal/(mol*K)','+|-',0.22065),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOCl_G4    |         45
CHOFCl_G4   |         26
CHOFClBr_G4 |         9
CHOClBr_G4  |         20
""",
)

entry(
    index = 195,
    label = "CdClCddO",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   O   u0 {1,S}
3   Cl  u0 {1,S}
4   Cdd u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.44892,7.93259,8.51799,9.06198,9.95124,10.4648,10.8514],'cal/(mol*K)','+|-',[0.234382,0.269375,0.275367,0.271327,0.235277,0.202924,0.165379]),
        H298 = (4.66086,'kcal/mol','+|-',0.835553),
        S298 = (16.2031,'cal/(mol*K)','+|-',0.548454),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOCl_G4    |         6
CHOFCl_G4   |         4
CHOFClBr_G4 |         1
CHOClBr_G4  |         2
""",
)

entry(
    index = 196,
    label = "Cd(Cdd-Od)ClO",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   O   u0 {1,S}
4   Cl  u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([13.5302,15.0416,16.3015,17.36,18.9226,19.9083,20.8767],'cal/(mol*K)','+|-',[0.478525,0.549967,0.562201,0.553953,0.480352,0.414299,0.337645]),
        H298 = (-1.12052,'kcal/mol','+|-',1.7059),
        S298 = (49.2622,'cal/(mol*K)','+|-',1.11975),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library    | Number of Species
CHOCl_G4   |         2
CHOClBr_G4 |         1
""",
)

entry(
    index = 197,
    label = "COFOO",
    group = 
"""
1 * CO u0 {2,S} {3,S} {4,D}
2   O  u0 {1,S}
3   F  u0 {1,S}
4   O  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOF_G4   |         6
CHOFCl_G4 |         1
CHOFBr_G4 |         1
""",
)

entry(
    index = 198,
    label = "CdCFO",
    group = 
"""
1 * Cd u0 {2,S} {3,S} {4,D}
2   O  u0 {1,S}
3   F  u0 {1,S}
4   C  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.48124,7.36014,8.11725,8.7856,9.8226,10.4951,11.0622],'cal/(mol*K)','+|-',[0.0932786,0.107205,0.109589,0.107982,0.0936347,0.0807591,0.0658169]),
        H298 = (-41.8157,'kcal/mol','+|-',0.33253),
        S298 = (12.7998,'cal/(mol*K)','+|-',0.218272),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOF_G4     |         46
CHOFCl_G4   |         17
CHOFClBr_G4 |         4
CHOFBr_G4   |         26
""",
)

entry(
    index = 199,
    label = "CdCddFO",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   O   u0 {1,S}
3   F   u0 {1,S}
4   Cdd u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.41997,7.39129,8.18373,8.83083,9.8356,10.4096,10.814],'cal/(mol*K)','+|-',[0.2555,0.293645,0.300178,0.295774,0.256476,0.221208,0.18028]),
        H298 = (-37.2323,'kcal/mol','+|-',0.910837),
        S298 = (13.1568,'cal/(mol*K)','+|-',0.59787),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOF_G4   |         6
CHOFCl_G4 |         2
CHOFBr_G4 |         3
""",
)

entry(
    index = 200,
    label = "Cd(Cdd-Od)FO",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   O   u0 {1,S}
4   F   u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([12.6636,14.2807,15.6019,16.7152,18.3753,19.4704,20.673],'cal/(mol*K)','+|-',[0.47857,0.550019,0.562254,0.554006,0.480397,0.414339,0.337677]),
        H298 = (-42.5571,'kcal/mol','+|-',1.70606),
        S298 = (45.4935,'cal/(mol*K)','+|-',1.11985),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOF_G4   |         2
CHOFBr_G4 |         1
""",
)

entry(
    index = 201,
    label = "COBrCO",
    group = 
"""
1 * CO u0 {2,S} {3,S} {4,D}
2   C  u0 {1,S}
3   Br u0 {1,S}
4   O  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([11.0443,12.2825,13.0421,13.6005,14.3203,14.695,14.6613],'cal/(mol*K)','+|-',[0.15381,0.176773,0.180705,0.178054,0.154397,0.133166,0.108527]),
        H298 = (-29.836,'kcal/mol','+|-',0.548319),
        S298 = (43.5232,'cal/(mol*K)','+|-',0.359914),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOBr_G4    |         9
CHOFClBr_G4 |         2
CHOFBr_G4   |         11
CHOClBr_G4  |         8
""",
)

entry(
    index = 202,
    label = "COBrCsO",
    group = 
"""
1 * CO u0 {2,S} {3,S} {4,D}
2   Cs u0 {1,S}
3   Br u0 {1,S}
4   O  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([10.1797,11.154,11.9627,12.6369,13.6522,14.292,14.9772],'cal/(mol*K)','+|-',[0.115427,0.13266,0.135611,0.133621,0.115867,0.0999347,0.0814445]),
        H298 = (-32.6061,'kcal/mol','+|-',0.411487),
        S298 = (45.2129,'cal/(mol*K)','+|-',0.270098),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOBr_G4    |         15
CHOFClBr_G4 |         7
CHOFBr_G4   |         22
CHOClBr_G4  |         14
""",
)

entry(
    index = 203,
    label = "CdBrCC",
    group = 
"""
1 * Cd u0 {2,S} {3,S} {4,D}
2   C  u0 {1,S}
3   Br u0 {1,S}
4   C  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.23239,9.08202,9.52339,9.96241,10.641,10.9934,11.0227],'cal/(mol*K)','+|-',[0.14949,0.171809,0.175631,0.173054,0.150061,0.129426,0.10548]),
        H298 = (12.0115,'kcal/mol','+|-',0.532921),
        S298 = (17.5739,'cal/(mol*K)','+|-',0.349807),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOBr_G4    |         15
CHOFClBr_G4 |         7
CHOFBr_G4   |         30
CHOClBr_G4  |         20
""",
)

entry(
    index = 204,
    label = "CdBrCsCd",
    group = 
"""
1 * Cd u0 {2,S} {3,S} {4,D}
2   Cs u0 {1,S}
3   Br u0 {1,S}
4   Cd u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.95215,8.65629,9.20499,9.62911,10.3256,10.7009,11.0195],'cal/(mol*K)','+|-',[0.0906525,0.104187,0.106504,0.104942,0.0909986,0.0784855,0.0639639]),
        H298 = (10.3483,'kcal/mol','+|-',0.323168),
        S298 = (17.2098,'cal/(mol*K)','+|-',0.212127),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOBr_G4    |         54
CHOFClBr_G4 |         14
CHOFBr_G4   |         75
CHOClBr_G4  |         41
""",
)

entry(
    index = 205,
    label = "CdBrCtCd",
    group = 
"""
1 * Cd u0 {2,S} {3,S} {4,D}
2   Ct u0 {1,S}
3   Br u0 {1,S}
4   Cd u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.44187,9.03983,10.1224,9.60488,10.2839,10.7362,11.2817],'cal/(mol*K)','+|-',[0.265245,0.304845,0.311627,0.307055,0.266258,0.229645,0.187156]),
        H298 = (10.5743,'kcal/mol','+|-',0.945577),
        S298 = (17.7803,'cal/(mol*K)','+|-',0.620673),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library    | Number of Species
CHOBr_G4   |         6
CHOFBr_G4  |         3
CHOClBr_G4 |         2
""",
)

entry(
    index = 206,
    label = "CdBrCddC",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   C   u0 {1,S}
3   Br  u0 {1,S}
4   Cdd u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.08435,8.57377,9.03308,9.40037,10.0388,10.3618,10.5848],'cal/(mol*K)','+|-',[0.201591,0.231688,0.236842,0.233367,0.202361,0.174534,0.142242]),
        H298 = (13.4805,'kcal/mol','+|-',0.718655),
        S298 = (18.6686,'cal/(mol*K)','+|-',0.471722),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOBr_G4    |         7
CHOFClBr_G4 |         1
CHOFBr_G4   |         7
CHOClBr_G4  |         4
""",
)

entry(
    index = 207,
    label = "Cd(Cdd-Od)CBr",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   C   u0 {1,S}
4   Br  u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([14.5571,15.7125,16.6396,17.3984,18.5526,19.263,20.0317],'cal/(mol*K)','+|-',[0.242221,0.278384,0.284577,0.280402,0.243146,0.209711,0.17091]),
        H298 = (2.0736,'kcal/mol','+|-',0.863499),
        S298 = (51.0559,'cal/(mol*K)','+|-',0.566797),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOBr_G4    |         3
CHOFClBr_G4 |         1
CHOFBr_G4   |         6
CHOClBr_G4  |         3
""",
)

entry(
    index = 208,
    label = "COCClO",
    group = 
"""
1 * CO u0 {2,S} {3,S} {4,D}
2   C  u0 {1,S}
3   Cl u0 {1,S}
4   O  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([10.5914,11.8451,12.6167,13.2224,14.0395,14.4934,14.6394],'cal/(mol*K)','+|-',[0.183857,0.211306,0.216006,0.212837,0.184558,0.15918,0.129728]),
        H298 = (-41.8655,'kcal/mol','+|-',0.655433),
        S298 = (40.7814,'cal/(mol*K)','+|-',0.430224),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library    | Number of Species
CHOCl_G4   |         10
CHOFCl_G4  |         8
CHOClBr_G4 |         1
""",
)

entry(
    index = 209,
    label = "COCsClO",
    group = 
"""
1 * CO u0 {2,S} {3,S} {4,D}
2   Cs u0 {1,S}
3   Cl u0 {1,S}
4   O  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.85506,10.9622,11.8168,12.5161,13.5699,14.244,14.9485],'cal/(mol*K)','+|-',[0.136634,0.157033,0.160526,0.158171,0.137156,0.118296,0.0964083]),
        H298 = (-44.7684,'kcal/mol','+|-',0.487089),
        S298 = (42.0164,'cal/(mol*K)','+|-',0.319723),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library    | Number of Species
CHOCl_G4   |         22
CHOFCl_G4  |         14
CHOClBr_G4 |         6
""",
)

entry(
    index = 210,
    label = "CdCCCl",
    group = 
"""
1 * Cd u0 {2,S} {3,S} {4,D}
2   C  u0 {1,S}
3   Cl u0 {1,S}
4   C  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.74827,8.85784,9.41817,9.93055,10.7147,11.1055,11.0823],'cal/(mol*K)','+|-',[0.155837,0.179102,0.183087,0.180401,0.156432,0.134921,0.109957]),
        H298 = (0.0247896,'kcal/mol','+|-',0.555544),
        S298 = (13.9508,'cal/(mol*K)','+|-',0.364657),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOCl_G4    |         22
CHOFCl_G4   |         16
CHOFClBr_G4 |         9
CHOClBr_G4  |         18
""",
)

entry(
    index = 211,
    label = "CdCsCdCl",
    group = 
"""
1 * Cd u0 {2,S} {3,S} {4,D}
2   Cs u0 {1,S}
3   Cl u0 {1,S}
4   Cd u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.55596,8.32131,8.89446,9.3577,10.1166,10.5297,10.8938],'cal/(mol*K)','+|-',[0.0904726,0.10398,0.106293,0.104733,0.090818,0.0783297,0.063837]),
        H298 = (-1.51023,'kcal/mol','+|-',0.322527),
        S298 = (14.3323,'cal/(mol*K)','+|-',0.211706),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOCl_G4    |         108
CHOFCl_G4   |         34
CHOFClBr_G4 |         18
CHOClBr_G4  |         60
""",
)

entry(
    index = 212,
    label = "CdCtCdCl",
    group = 
"""
1 * Cd u0 {2,S} {3,S} {4,D}
2   Ct u0 {1,S}
3   Cl u0 {1,S}
4   Cd u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.88057,8.51032,9.62127,9.16703,10.0009,10.5339,11.0764],'cal/(mol*K)','+|-',[0.238947,0.274621,0.28073,0.276611,0.239859,0.206877,0.1686]),
        H298 = (-0.800178,'kcal/mol','+|-',0.851826),
        S298 = (14.5227,'cal/(mol*K)','+|-',0.559135),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOCl_G4    |         6
CHOFCl_G4   |         2
CHOFClBr_G4 |         1
CHOClBr_G4  |         4
""",
)

entry(
    index = 213,
    label = "CdCddCCl",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   C   u0 {1,S}
3   Cl  u0 {1,S}
4   Cdd u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.81938,8.3764,8.87298,9.24504,9.91995,10.283,10.5764],'cal/(mol*K)','+|-',[0.184945,0.212556,0.217285,0.214097,0.185651,0.160122,0.130496]),
        H298 = (0.832565,'kcal/mol','+|-',0.659313),
        S298 = (15.2638,'cal/(mol*K)','+|-',0.43277),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOCl_G4    |         12
CHOFCl_G4   |         4
CHOFClBr_G4 |         1
CHOClBr_G4  |         7
""",
)

entry(
    index = 214,
    label = "Cd(Cdd-Od)CCl",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   C   u0 {1,S}
4   Cl  u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([14.2714,15.5562,16.5089,17.285,18.4734,19.2069,19.9873],'cal/(mol*K)','+|-',[0.321664,0.369687,0.377911,0.372366,0.322892,0.278491,0.226964]),
        H298 = (-8.11511,'kcal/mol','+|-',1.1467),
        S298 = (47.7699,'cal/(mol*K)','+|-',0.752692),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOCl_G4  |         4
CHOFCl_G4 |         3
""",
)

entry(
    index = 215,
    label = "COCFO",
    group = 
"""
1 * CO u0 {2,S} {3,S} {4,D}
2   C  u0 {1,S}
3   F  u0 {1,S}
4   O  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.22943,10.5464,11.5362,12.3222,13.3664,13.9897,14.5745],'cal/(mol*K)','+|-',[0.215176,0.247301,0.252803,0.249094,0.215998,0.186296,0.151827]),
        H298 = (-87.6805,'kcal/mol','+|-',0.767085),
        S298 = (38.0895,'cal/(mol*K)','+|-',0.503512),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOF_G4   |         11
CHOFCl_G4 |         1
CHOFBr_G4 |         1
""",
)

entry(
    index = 216,
    label = "COCsFO",
    group = 
"""
1 * CO u0 {2,S} {3,S} {4,D}
2   Cs u0 {1,S}
3   F  u0 {1,S}
4   O  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOF_G4     |         22
CHOFCl_G4   |         6
CHOFClBr_G4 |         2
CHOFBr_G4   |         9
""",
)

entry(
    index = 217,
    label = "CdCCF",
    group = 
"""
1 * Cd u0 {2,S} {3,S} {4,D}
2   C  u0 {1,S}
3   F  u0 {1,S}
4   C  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.78714,8.13666,8.89332,9.53102,10.4675,10.9499,11.0244],'cal/(mol*K)','+|-',[0.142112,0.163329,0.166962,0.164513,0.142655,0.123039,0.100274]),
        H298 = (-39.4225,'kcal/mol','+|-',0.506618),
        S298 = (10.8538,'cal/(mol*K)','+|-',0.332542),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOF_G4     |         22
CHOFCl_G4   |         14
CHOFClBr_G4 |         6
CHOFBr_G4   |         27
""",
)

entry(
    index = 218,
    label = "CdCsCdF",
    group = 
"""
1 * Cd u0 {2,S} {3,S} {4,D}
2   Cs u0 {1,S}
3   F  u0 {1,S}
4   Cd u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOF_G4     |         108
CHOFCl_G4   |         46
CHOFClBr_G4 |         26
CHOFBr_G4   |         97
""",
)

entry(
    index = 219,
    label = "CdCtCdF",
    group = 
"""
1 * Cd u0 {2,S} {3,S} {4,D}
2   Ct u0 {1,S}
3   F  u0 {1,S}
4   Cd u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.81544,7.58852,8.81763,8.44556,9.44655,10.1128,10.8098],'cal/(mol*K)','+|-',[0.210842,0.24232,0.247711,0.244076,0.211647,0.182544,0.148769]),
        H298 = (-38.5774,'kcal/mol','+|-',0.751635),
        S298 = (11.7419,'cal/(mol*K)','+|-',0.49337),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOF_G4     |         6
CHOFCl_G4   |         4
CHOFClBr_G4 |         1
CHOFBr_G4   |         6
""",
)

entry(
    index = 220,
    label = "CdCddCF",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   C   u0 {1,S}
3   F   u0 {1,S}
4   Cdd u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.87064,7.56576,8.17949,8.64795,9.46354,9.93621,10.4467],'cal/(mol*K)','+|-',[0.172186,0.197893,0.202295,0.199327,0.172844,0.149076,0.121494]),
        H298 = (-36.5736,'kcal/mol','+|-',0.613829),
        S298 = (12.3646,'cal/(mol*K)','+|-',0.402915),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOF_G4     |         12
CHOFCl_G4   |         4
CHOFClBr_G4 |         4
CHOFBr_G4   |         10
""",
)

entry(
    index = 221,
    label = "Cd(Cdd-Od)CF",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   C   u0 {1,S}
4   F   u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([13.2731,14.8399,15.9997,16.9109,18.2099,19.0283,20.0356],'cal/(mol*K)','+|-',[0.423039,0.486197,0.497013,0.489721,0.424654,0.366261,0.298494]),
        H298 = (-39.5676,'kcal/mol','+|-',1.5081),
        S298 = (45.0224,'cal/(mol*K)','+|-',0.989911),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHOF_G4 |         4
""",
)

entry(
    index = 222,
    label = "Cds-OdHH",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   H   u0 {1,S}
4   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """CO-HH BENSON !!!WARNING! Cp1500 value taken as Cp1000, S(group) = S(CH2O) + Rln(2)""",
    longDesc = 
"""

""",
)

entry(
    index = 223,
    label = "Cds-OdOsH",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   O2s u0 {1,S}
4   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([25.88,34.56,42.08,48.16,56.57,61.38,65.84],'J/(mol*K)','+|-',[4.01,4.01,4.01,4.01,4.01,4.01,4.01]),
        H298 = (-211.8,'kJ/mol','+|-',3.42),
        S298 = (124.04,'J/(mol*K)','+|-',4.68),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 224,
    label = "CO-SH",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   S   u0 {1,S}
4   H   u0 {1,S}
""",
    thermo = 'CO-S2H',
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2018""",
    longDesc = 
"""
"
""",
)

entry(
    index = 225,
    label = "CO-S2H",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   S2s u0 {1,S}
4   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.27,9.88,10.36,10.77,12.65,14.02,14.06],'cal/(mol*K)'),
        H298 = (-15.2,'kcal/mol'),
        S298 = (41.26,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 226,
    label = "CO-S4H",
    group = 
"""
1 * CO                u0 {2,D} {3,S} {4,S}
2   O2d               u0 {1,D}
3   [S4s,S4d,S4b,S4t] u0 {1,S}
4   H                 u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.51,9.3,9.66,10.35,13.15,14.94,15.19],'cal/(mol*K)'),
        H298 = (-8.53,'kcal/mol'),
        S298 = (39.05,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 227,
    label = "CO-S6H",
    group = 
"""
1 * CO                      u0 {2,D} {3,S} {4,S}
2   O2d                     u0 {1,D}
3   [S6s,S6d,S6dd,S6t,S6td] u0 {1,S}
4   H                       u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.08,7.22,8.24,9.21,11.83,13.5,14.27],'cal/(mol*K)'),
        H298 = (-8.01,'kcal/mol'),
        S298 = (48.01,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 228,
    label = "Cds-OdOsOs",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   O2s u0 {1,S}
4   O2s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([26.17,39.3,48.25,53.88,58.97,59.63,56.09],'J/(mol*K)','+|-',[6,6,6,6,6,6,6]),
        H298 = (-281.4,'kJ/mol','+|-',5.11),
        S298 = (22.66,'J/(mol*K)','+|-',7),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 229,
    label = "CO-CsSs",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   S2s u0 {1,S}
4   Cs  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.12,8.65,9.04,9.38,11.01,11.97,10.97],'cal/(mol*K)'),
        H298 = (-19.07,'kcal/mol'),
        S298 = (20.3,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 230,
    label = "CO-OsSs",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   O2s u0 {1,S}
4   S2s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.87,9.51,10.75,11.62,13.62,14.53,12.86],'cal/(mol*K)'),
        H298 = (-35.59,'kcal/mol'),
        S298 = (16.37,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""

""",
)

entry(
    index = 231,
    label = "Cds-OdCH",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   C   u0 {1,S}
4   H   u0 {1,S}
""",
    thermo = 'Cds-OdCsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 232,
    label = "Cds-OdCsH",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   Cs  u0 {1,S}
4   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 233,
    label = "Cds-OdCdsH",
    group = 
"""
1 * CO      u0 {2,D} {3,S} {4,S}
2   O2d     u0 {1,D}
3   [Cd,CO] u0 {1,S}
4   H       u0 {1,S}
""",
    thermo = 'Cds-O2d(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 234,
    label = "Cds-O2d(Cds-O2d)H",
    group = 
"""
1 * CO  u0 {2,S} {3,D} {4,S}
2   CO  u0 {1,S} {5,D}
3   O2d u0 {1,D}
4   H   u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([29.76,34.63,39.25,43.32,49.57,53.77,59.32],'J/(mol*K)','+|-',[1.7,1.7,1.7,1.7,1.7,1.7,1.7]),
        H298 = (-104.8,'kJ/mol','+|-',1.45),
        S298 = (140.49,'J/(mol*K)','+|-',1.98),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 235,
    label = "Cds-O2d(Cds-Cd)H",
    group = 
"""
1 * CO  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   O2d u0 {1,D}
4   H   u0 {1,S}
5   C   u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([27.31,34,39.42,43.77,50.16,54.55,60.77],'J/(mol*K)','+|-',[4.9,4.9,4.9,4.9,4.9,4.9,4.9]),
        H298 = (-128.3,'kJ/mol','+|-',5.9),
        S298 = (129.26,'J/(mol*K)','+|-',5.71),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 236,
    label = "Cds-O2d(Cds-Cds)H",
    group = 
"""
1 * CO  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   O2d u0 {1,D}
4   H   u0 {1,S}
5   Cd  u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.45,8.77,9.82,10.7,12.14,12.9,12.9],'cal/(mol*K)','+|-',[0.15,0.15,0.15,0.15,0.15,0.15,0.15]),
        H298 = (-30.9,'kcal/mol','+|-',0.3),
        S298 = (33.4,'cal/(mol*K)','+|-',0.15),
    ),
    shortDesc = """CO-CdH Hf BOZZELLI S,Cp =3D CO/C/H-del(cd syst) !!!WARNING! Cp1500 value taken as Cp1000""",
    longDesc = 
"""

""",
)

entry(
    index = 237,
    label = "Cds-O2d(Cds-Cdd)H",
    group = 
"""
1 * CO  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   O2d u0 {1,D}
4   H   u0 {1,S}
5   Cdd u0 {2,D}
""",
    thermo = 'Cds-O2d(Cds-Cdd-Cd)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 238,
    label = "Cds-O2d(Cds-Cdd-O2d)H",
    group = 
"""
1 * CO  u0 {2,S} {4,D} {5,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {6,D}
4   O2d u0 {1,D}
5   H   u0 {1,S}
6   O2d u0 {3,D}
""",
    thermo = 'Cds-O2d(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 239,
    label = "Cds-O2d(Cds-Cdd-Cd)H",
    group = 
"""
1 * CO  u0 {2,S} {4,D} {5,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {6,D}
4   O2d u0 {1,D}
5   H   u0 {1,S}
6   C   u0 {3,D}
""",
    thermo = 'Cds-O2d(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 240,
    label = "Cds-OdCtH",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   Ct  u0 {1,S}
4   H   u0 {1,S}
""",
    thermo = 'Cds-O2d(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 241,
    label = "Cds-OdCbH",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   Cb  u0 {1,S}
4   H   u0 {1,S}
""",
    thermo = 'Cds-O2d(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 242,
    label = "Cds-OdCOs",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   C   u0 {1,S}
4   O2s u0 {1,S}
""",
    thermo = 'Cds-OdCsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 243,
    label = "Cds-OdCsOs",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   Cs  u0 {1,S}
4   O2s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 244,
    label = "Cds-OdCdsOs",
    group = 
"""
1 * CO      u0 {2,D} {3,S} {4,S}
2   O2d     u0 {1,D}
3   [Cd,CO] u0 {1,S}
4   O2s     u0 {1,S}
""",
    thermo = 'Cds-O2d(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 245,
    label = "Cds-O2d(Cds-O2d)O2s",
    group = 
"""
1 * CO  u0 {2,S} {3,D} {4,S}
2   CO  u0 {1,S} {5,D}
3   O2d u0 {1,D}
4   O2s u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([27.18,34.34,39.85,44.13,49.81,52.4,52.33],'J/(mol*K)','+|-',[3.36,3.36,3.36,3.36,3.36,3.36,3.36]),
        H298 = (-196.2,'kJ/mol','+|-',2.86),
        S298 = (39.37,'J/(mol*K)','+|-',3.92),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 246,
    label = "Cds-O2d(Cds-Cd)O2s",
    group = 
"""
1 * CO  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   O2d u0 {1,D}
4   O2s u0 {1,S}
5   C   u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([28.33,37.84,44.54,49.34,55.45,58.73,60.61],'J/(mol*K)','+|-',[7.46,7.46,7.46,7.46,7.46,7.46,7.46]),
        H298 = (-218.6,'kJ/mol','+|-',6.36),
        S298 = (33.44,'J/(mol*K)','+|-',8.7),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 247,
    label = "Cds-O2d(Cds-Cds)O2s",
    group = 
"""
1 * CO  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   O2d u0 {1,D}
4   O2s u0 {1,S}
5   Cd  u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.97,6.7,7.4,8.02,8.87,9.36,9.36],'cal/(mol*K)','+|-',[0.15,0.15,0.15,0.15,0.15,0.15,0.15]),
        H298 = (-32.1,'kcal/mol','+|-',0.3),
        S298 = (14.78,'cal/(mol*K)','+|-',0.15),
    ),
    shortDesc = """CO-OCd RPS + S Coreected !!!WARNING! Cp1500 value taken as Cp1000""",
    longDesc = 
"""

""",
)

entry(
    index = 248,
    label = "Cds-O2d(Cds-Cdd)O2s",
    group = 
"""
1 * CO  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   O2d u0 {1,D}
4   O2s u0 {1,S}
5   Cdd u0 {2,D}
""",
    thermo = 'Cds-O2d(Cds-Cdd-Cd)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 249,
    label = "Cds-O2d(Cds-Cdd-O2d)O2s",
    group = 
"""
1 * CO  u0 {2,S} {4,D} {5,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {6,D}
4   O2d u0 {1,D}
5   O2s u0 {1,S}
6   O2d u0 {3,D}
""",
    thermo = 'Cds-O2d(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 250,
    label = "Cds-O2d(Cds-Cdd-Cd)O2s",
    group = 
"""
1 * CO  u0 {2,S} {4,D} {5,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {6,D}
4   O2d u0 {1,D}
5   O2s u0 {1,S}
6   C   u0 {3,D}
""",
    thermo = 'Cds-O2d(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 251,
    label = "Cds-OdCtOs",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   Ct  u0 {1,S}
4   O2s u0 {1,S}
""",
    thermo = 'Cds-O2d(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 252,
    label = "Cds-OdCbOs",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   Cb  u0 {1,S}
4   O2s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.97,6.7,7.4,8.02,8.87,9.36,9.36],'cal/(mol*K)','+|-',[0.15,0.15,0.15,0.15,0.15,0.15,0.15]),
        H298 = (-36.6,'kcal/mol','+|-',0.3),
        S298 = (14.78,'cal/(mol*K)','+|-',0.15),
    ),
    shortDesc = """CO-OCb RPS + S Coreected !!!WARNING! Cp1500 value taken as Cp1000""",
    longDesc = 
"""

""",
)

entry(
    index = 253,
    label = "Cds-OdCC",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   C   u0 {1,S}
4   C   u0 {1,S}
""",
    thermo = 'Cds-OdCsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 254,
    label = "Cds-OdCsCs",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   Cs  u0 {1,S}
4   Cs  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 255,
    label = "Cds-OdCdsCs",
    group = 
"""
1 * CO      u0 {2,D} {3,S} {4,S}
2   O2d     u0 {1,D}
3   [Cd,CO] u0 {1,S}
4   Cs      u0 {1,S}
""",
    thermo = 'Cds-O2d(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 256,
    label = "Cds-O2d(Cds-O2d)Cs",
    group = 
"""
1 * CO  u0 {2,S} {3,D} {4,S}
2   CO  u0 {1,S} {5,D}
3   O2d u0 {1,D}
4   Cs  u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([26.77,30.83,34.36,37.27,41.27,43.45,45.25],'J/(mol*K)','+|-',[1.28,1.28,1.28,1.28,1.28,1.28,1.28]),
        H298 = (-122,'kJ/mol','+|-',1.09),
        S298 = (57.8,'J/(mol*K)','+|-',1.5),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 257,
    label = "Cds-O2d(Cds-Cd)Cs",
    group = 
"""
1 * CO  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   O2d u0 {1,D}
4   Cs  u0 {1,S}
5   C   u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([25.26,30.66,34.68,37.69,41.62,43.93,46.69],'J/(mol*K)','+|-',[4.9,4.9,4.9,4.9,4.9,4.9,4.9]),
        H298 = (-130.4,'kJ/mol','+|-',4.17),
        S298 = (47.38,'J/(mol*K)','+|-',5.71),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 258,
    label = "Cds-O2d(Cds-Cds)Cs",
    group = 
"""
1 * CO  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   O2d u0 {1,D}
4   Cs  u0 {1,S}
5   Cd  u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.46,6.32,7.17,7.88,9,9.77,9.77],'cal/(mol*K)','+|-',[0.15,0.15,0.15,0.15,0.15,0.15,0.15]),
        H298 = (-30.9,'kcal/mol','+|-',0.3),
        S298 = (14.6,'cal/(mol*K)','+|-',0.15),
    ),
    shortDesc = """CO-CdCs Hf BENSON =3D CO/Cb/C S,Cp !!!WARNING! Cp1500 value taken as Cp1000""",
    longDesc = 
"""

""",
)

entry(
    index = 259,
    label = "Cds-O2d(Cds-Cdd)Cs",
    group = 
"""
1 * CO  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   O2d u0 {1,D}
4   Cs  u0 {1,S}
5   Cdd u0 {2,D}
""",
    thermo = 'Cds-O2d(Cds-Cdd-Cd)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 260,
    label = "Cds-O2d(Cds-Cdd-O2d)Cs",
    group = 
"""
1 * CO  u0 {2,S} {4,D} {5,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {6,D}
4   O2d u0 {1,D}
5   Cs  u0 {1,S}
6   O2d u0 {3,D}
""",
    thermo = 'Cds-O2d(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 261,
    label = "Cds-O2d(Cds-Cdd-Cd)Cs",
    group = 
"""
1 * CO  u0 {2,S} {4,D} {5,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {6,D}
4   O2d u0 {1,D}
5   Cs  u0 {1,S}
6   C   u0 {3,D}
""",
    thermo = 'Cds-O2d(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 262,
    label = "Cds-OdCdsCds",
    group = 
"""
1 * CO      u0 {2,D} {3,S} {4,S}
2   O2d     u0 {1,D}
3   [Cd,CO] u0 {1,S}
4   [Cd,CO] u0 {1,S}
""",
    thermo = 'Cds-O2d(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 263,
    label = "Cds-O2d(Cds-O2d)(Cds-O2d)",
    group = 
"""
1 * CO  u0 {2,S} {3,S} {4,D}
2   CO  u0 {1,S} {5,D}
3   CO  u0 {1,S} {6,D}
4   O2d u0 {1,D}
5   O2d u0 {2,D}
6   O2d u0 {3,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([31.75,33.35,34.1,34.51,35.19,36.06,38.14],'J/(mol*K)','+|-',[2.41,2.41,2.41,2.41,2.41,2.41,2.41]),
        H298 = (-89.3,'kJ/mol','+|-',2.05),
        S298 = (64.51,'J/(mol*K)','+|-',2.81),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 264,
    label = "Cds-O2d(Cds-Cd)(Cds-O2d)",
    group = 
"""
1 * CO  u0 {2,S} {3,S} {4,D}
2   Cd  u0 {1,S} {5,D}
3   CO  u0 {1,S} {6,D}
4   O2d u0 {1,D}
5   C   u0 {2,D}
6   O2d u0 {3,D}
""",
    thermo = 'Cds-O2d(Cds-Cds)(Cds-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 265,
    label = "Cds-O2d(Cds-Cds)(Cds-O2d)",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   Cd  u0 {1,S} {5,D}
4   CO  u0 {1,S} {6,D}
5   Cd  u0 {3,D}
6   O2d u0 {4,D}
""",
    thermo = 'Cds-O2d(Cds-O2d)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 266,
    label = "Cds-O2d(Cds-Cdd)(Cds-O2d)",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   Cd  u0 {1,S} {5,D}
4   CO  u0 {1,S} {6,D}
5   Cdd u0 {3,D}
6   O2d u0 {4,D}
""",
    thermo = 'Cds-O2d(Cds-Cdd-Cd)(Cds-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 267,
    label = "Cds-O2d(Cds-Cdd-O2d)(Cds-O2d)",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   Cd  u0 {1,S} {5,D}
4   CO  u0 {1,S} {7,D}
5   Cdd u0 {3,D} {6,D}
6   O2d u0 {5,D}
7   O2d u0 {4,D}
""",
    thermo = 'Cds-O2d(Cds-Cdd-O2d)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 268,
    label = "Cds-O2d(Cds-Cdd-Cd)(Cds-O2d)",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   Cd  u0 {1,S} {5,D}
4   CO  u0 {1,S} {7,D}
5   Cdd u0 {3,D} {6,D}
6   C   u0 {5,D}
7   O2d u0 {4,D}
""",
    thermo = 'Cds-O2d(Cds-Cds)(Cds-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 269,
    label = "Cds-O2d(Cds-Cd)(Cds-Cd)",
    group = 
"""
1 * CO  u0 {2,S} {3,S} {4,D}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,S} {6,D}
4   O2d u0 {1,D}
5   C   u0 {2,D}
6   C   u0 {3,D}
""",
    thermo = 'Cds-O2d(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 270,
    label = "Cds-O2d(Cds-Cds)(Cds-Cds)",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   Cd  u0 {1,S} {5,D}
4   Cd  u0 {1,S} {6,D}
5   Cd  u0 {3,D}
6   Cd  u0 {4,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.46,6.32,7.17,7.88,9,9.77,9.77],'cal/(mol*K)','+|-',[0.15,0.15,0.15,0.15,0.15,0.15,0.15]),
        H298 = (-30.9,'kcal/mol','+|-',0.3),
        S298 = (14.6,'cal/(mol*K)','+|-',0.15),
    ),
    shortDesc = """CO-CdCd Estimate BOZZELLI !!!WARNING! Cp1500 value taken as Cp1000""",
    longDesc = 
"""

""",
)

entry(
    index = 271,
    label = "Cds-O2d(Cds-Cdd)(Cds-Cds)",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   Cd  u0 {1,S} {5,D}
4   Cd  u0 {1,S} {6,D}
5   Cdd u0 {3,D}
6   Cd  u0 {4,D}
""",
    thermo = 'Cds-O2d(Cds-Cdd-Cd)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 272,
    label = "Cds-O2d(Cds-Cdd-O2d)(Cds-Cds)",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   Cd  u0 {1,S} {5,D}
4   Cd  u0 {1,S} {6,D}
5   Cdd u0 {3,D} {7,D}
6   Cd  u0 {4,D}
7   O2d u0 {5,D}
""",
    thermo = 'Cds-O2d(Cds-Cdd-O2d)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 273,
    label = "Cds-O2d(Cds-Cdd-Cd)(Cds-Cds)",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   Cd  u0 {1,S} {5,D}
4   Cd  u0 {1,S} {6,D}
5   Cdd u0 {3,D} {7,D}
6   Cd  u0 {4,D}
7   C   u0 {5,D}
""",
    thermo = 'Cds-O2d(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 274,
    label = "Cds-O2d(Cds-Cdd)(Cds-Cdd)",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   Cd  u0 {1,S} {5,D}
4   Cd  u0 {1,S} {6,D}
5   Cdd u0 {3,D}
6   Cdd u0 {4,D}
""",
    thermo = 'Cds-O2d(Cds-Cdd-Cd)(Cds-Cdd-Cd)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 275,
    label = "Cds-O2d(Cds-Cdd-O2d)(Cds-Cdd-O2d)",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   Cd  u0 {1,S} {5,D}
4   Cd  u0 {1,S} {6,D}
5   Cdd u0 {3,D} {7,D}
6   Cdd u0 {4,D} {8,D}
7   O2d u0 {5,D}
8   O2d u0 {6,D}
""",
    thermo = 'Cds-O2d(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 276,
    label = "Cds-O2d(Cds-Cdd-Cd)(Cds-Cdd-O2d)",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   Cd  u0 {1,S} {5,D}
4   Cd  u0 {1,S} {6,D}
5   Cdd u0 {3,D} {7,D}
6   Cdd u0 {4,D} {8,D}
7   C   u0 {5,D}
8   O2d u0 {6,D}
""",
    thermo = 'Cds-O2d(Cds-Cdd-O2d)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 277,
    label = "Cds-O2d(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   Cd  u0 {1,S} {5,D}
4   Cd  u0 {1,S} {6,D}
5   Cdd u0 {3,D} {7,D}
6   Cdd u0 {4,D} {8,D}
7   C   u0 {5,D}
8   C   u0 {6,D}
""",
    thermo = 'Cds-O2d(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 278,
    label = "Cds-OdCtCs",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   Ct  u0 {1,S}
4   Cs  u0 {1,S}
""",
    thermo = 'Cds-O2d(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 279,
    label = "Cds-OdCtCds",
    group = 
"""
1 * CO      u0 {2,D} {3,S} {4,S}
2   O2d     u0 {1,D}
3   Ct      u0 {1,S}
4   [Cd,CO] u0 {1,S}
""",
    thermo = 'Cds-OdCt(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 280,
    label = "Cds-OdCt(Cds-O2d)",
    group = 
"""
1 * CO  u0 {2,S} {3,D} {4,S}
2   CO  u0 {1,S} {5,D}
3   O2d u0 {1,D}
4   Ct  u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = 'Cds-O2d(Cds-Cds)(Cds-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 281,
    label = "Cds-OdCt(Cds-Cd)",
    group = 
"""
1 * CO  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   O2d u0 {1,D}
4   Ct  u0 {1,S}
5   C   u0 {2,D}
""",
    thermo = 'Cds-OdCt(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 282,
    label = "Cds-OdCt(Cds-Cds)",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   Ct  u0 {1,S}
4   Cd  u0 {1,S} {5,D}
5   Cd  u0 {4,D}
""",
    thermo = 'Cds-O2d(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 283,
    label = "Cds-OdCt(Cds-Cdd)",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   Ct  u0 {1,S}
4   Cd  u0 {1,S} {5,D}
5   Cdd u0 {4,D}
""",
    thermo = 'Cds-OdCt(Cds-Cdd-Cd)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 284,
    label = "Cds-OdCt(Cds-Cdd-O2d)",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   Ct  u0 {1,S}
4   Cd  u0 {1,S} {5,D}
5   Cdd u0 {4,D} {6,D}
6   O2d u0 {5,D}
""",
    thermo = 'Cds-O2d(Cds-Cdd-O2d)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 285,
    label = "Cds-OdCt(Cds-Cdd-Cd)",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   Ct  u0 {1,S}
4   Cd  u0 {1,S} {5,D}
5   Cdd u0 {4,D} {6,D}
6   C   u0 {5,D}
""",
    thermo = 'Cds-OdCt(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 286,
    label = "Cds-OdCtCt",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   Ct  u0 {1,S}
4   Ct  u0 {1,S}
""",
    thermo = 'Cds-O2d(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 287,
    label = "Cds-OdCbCs",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   Cb  u0 {1,S}
4   Cs  u0 {1,S}
""",
    thermo = 'Cds-O2d(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 288,
    label = "Cds-OdCbCds",
    group = 
"""
1 * CO      u0 {2,D} {3,S} {4,S}
2   O2d     u0 {1,D}
3   Cb      u0 {1,S}
4   [Cd,CO] u0 {1,S}
""",
    thermo = 'Cds-OdCb(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 289,
    label = "Cds-OdCb(Cds-O2d)",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   Cb  u0 {1,S}
4   CO  u0 {1,S} {5,D}
5   O2d u0 {4,D}
""",
    thermo = 'Cds-O2d(Cds-Cds)(Cds-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 290,
    label = "Cds-OdCb(Cds-Cd)",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   Cb  u0 {1,S}
4   Cd  u0 {1,S} {5,D}
5   C   u0 {4,D}
""",
    thermo = 'Cds-OdCb(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 291,
    label = "Cds-OdCb(Cds-Cds)",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   Cb  u0 {1,S}
4   Cd  u0 {1,S} {5,D}
5   Cd  u0 {4,D}
""",
    thermo = 'Cds-O2d(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 292,
    label = "Cds-OdCb(Cds-Cdd)",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   Cb  u0 {1,S}
4   Cd  u0 {1,S} {5,D}
5   Cdd u0 {4,D}
""",
    thermo = 'Cds-OdCb(Cds-Cdd-Cd)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 293,
    label = "Cds-OdCb(Cds-Cdd-O2d)",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   Cb  u0 {1,S}
4   Cd  u0 {1,S} {5,D}
5   Cdd u0 {4,D} {6,D}
6   O2d u0 {5,D}
""",
    thermo = 'Cds-O2d(Cds-Cdd-O2d)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 294,
    label = "Cds-OdCb(Cds-Cdd-Cd)",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   Cb  u0 {1,S}
4   Cd  u0 {1,S} {5,D}
5   Cdd u0 {4,D} {6,D}
6   C   u0 {5,D}
""",
    thermo = 'Cds-OdCb(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 295,
    label = "Cds-OdCbCt",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   Cb  u0 {1,S}
4   Ct  u0 {1,S}
""",
    thermo = 'Cds-OdCt(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 296,
    label = "Cds-OdCbCb",
    group = 
"""
1 * CO  u0 {2,D} {3,S} {4,S}
2   O2d u0 {1,D}
3   Cb  u0 {1,S}
4   Cb  u0 {1,S}
""",
    thermo = 'Cds-O2d(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 297,
    label = "Cds-CdHH",
    group = 
"""
1 * Cd u0 {2,D} {3,S} {4,S}
2   C  u0 {1,D}
3   H  u0 {1,S}
4   H  u0 {1,S}
""",
    thermo = 'Cds-CdsHH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 298,
    label = "Cds-CdsHH",
    group = 
"""
1 * Cd u0 {2,D} {3,S} {4,S}
2   Cd u0 {1,D}
3   H  u0 {1,S}
4   H  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Cd-HH BENSON""",
    longDesc = 
"""

""",
)

entry(
    index = 299,
    label = "Cds-CddHH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D}
3   H   u0 {1,S}
4   H   u0 {1,S}
""",
    thermo = 'Cds-(Cdd-Cd)HH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 300,
    label = "Cds-(Cdd-O2d)HH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   H   u0 {1,S}
4   H   u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """{CCO/H2} RAMAN & GREEN JPCA 2002, 106, 7937-7949""",
    longDesc = 
"""

""",
)

entry(
    index = 301,
    label = "Cds-(Cdd-S2d)HH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   H   u0 {1,S}
4   H   u0 {1,S}
5   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 302,
    label = "Cds-(Cdd-Cd)HH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   H   u0 {1,S}
4   H   u0 {1,S}
5   C   u0 {2,D}
""",
    thermo = 'Cds-CdsHH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 303,
    label = "Cds-CdOsH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   C   u0 {1,D}
3   O2s u0 {1,S}
4   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([18.08,21.17,24.43,27.41,32.22,35.73,40.97],'J/(mol*K)'),
        H298 = (36.4,'kJ/mol'),
        S298 = (33.51,'J/(mol*K)'),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 304,
    label = "Cds-CdsOsH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cd  u0 {1,D}
3   O2s u0 {1,S}
4   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.75,6.46,7.64,8.35,9.1,9.56,10.46],'cal/(mol*K)','+|-',[0.07,0.07,0.07,0.07,0.07,0.07,0.07]),
        H298 = (2.03,'kcal/mol','+|-',0.19),
        S298 = (6.2,'cal/(mol*K)','+|-',0.1),
    ),
    shortDesc = """Cd-OH BOZZELLI Hf vin-oh RADOM + C/Cd/H, S&Cp LAY""",
    longDesc = 
"""

""",
)

entry(
    index = 305,
    label = "Cds-CddOsH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D}
3   O2s u0 {1,S}
4   H   u0 {1,S}
""",
    thermo = 'Cds-(Cdd-Cd)OsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 306,
    label = "Cds-(Cdd-O2d)OsH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   O2s u0 {1,S}
4   H   u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([11.29,13.67,15.1,16.1,17.36,18.25,19.75],'cal/(mol*K)','+|-',[0.07,0.07,0.07,0.07,0.07,0.07,0.07]),
        H298 = (2.11,'kcal/mol','+|-',0.19),
        S298 = (38.17,'cal/(mol*K)','+|-',0.1),
    ),
    shortDesc = """{CCO/O/H} RAMAN & GREEN JPCA 2002, 106, 7937-7949""",
    longDesc = 
"""

""",
)

entry(
    index = 307,
    label = "Cds-(Cdd-Cd)OsH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   O2s u0 {1,S}
4   H   u0 {1,S}
5   C   u0 {2,D}
""",
    thermo = 'Cds-CdsOsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 308,
    label = "Cds-CdSH",
    group = 
"""
1 * Cd u0 {2,D} {3,S} {4,S}
2   C  u0 {1,D}
3   S  u0 {1,S}
4   H  u0 {1,S}
""",
    thermo = 'Cds-CdsSH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 309,
    label = "Cds-CdsSH",
    group = 
"""
1 * Cd u0 {2,D} {3,S} {4,S}
2   Cd u0 {1,D}
3   S  u0 {1,S}
4   H  u0 {1,S}
""",
    thermo = 'Cds-CdsS4H',
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2018""",
    longDesc = 
"""
"
""",
)

entry(
    index = 310,
    label = "Cds-CdsS2H",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cd  u0 {1,D}
3   S2s u0 {1,S}
4   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.72,6.61,6.84,7.2,9.2,10.44,9.73],'cal/(mol*K)'),
        H298 = (18.92,'kcal/mol'),
        S298 = (12.2,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 311,
    label = "Cds-CdsS4H",
    group = 
"""
1 * Cd                u0 {2,D} {3,S} {4,S}
2   Cd                u0 {1,D}
3   [S4s,S4d,S4b,S4t] u0 {1,S}
4   H                 u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.22,7.89,7.54,8.28,10.6,11.9,10.95],'cal/(mol*K)'),
        H298 = (27.56,'kcal/mol'),
        S298 = (5.06,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 312,
    label = "Cds-CdsS6H",
    group = 
"""
1 * Cd                      u0 {2,D} {3,S} {4,S}
2   Cd                      u0 {1,D}
3   [S6s,S6d,S6dd,S6t,S6td] u0 {1,S}
4   H                       u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.03,5.53,6.45,7.26,9.65,11,10.77],'cal/(mol*K)'),
        H298 = (19.86,'kcal/mol'),
        S298 = (13.58,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 313,
    label = "Cds-CddSsH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D}
3   S2s u0 {1,S}
4   H   u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 314,
    label = "Cds-(Cdd-S2d)SsH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   S2s u0 {1,S}
4   H   u0 {1,S}
5   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 315,
    label = "Cds-(Cdd-Cd)SsH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   S2s u0 {1,S}
4   H   u0 {1,S}
5   C   u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 316,
    label = "Cds-CdOsOs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   C   u0 {1,D}
3   O2s u0 {1,S}
4   O2s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([11.34,11.93,14.86,17.95,22.31,24.6,26.92],'J/(mol*K)','+|-',[7.4,7.4,7.4,7.4,7.4,7.4,7.4]),
        H298 = (28.3,'kJ/mol','+|-',6.3),
        S298 = (-42.69,'J/(mol*K)','+|-',8.63),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 317,
    label = "Cds-CdsOsOs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cd  u0 {1,D}
3   O2s u0 {1,S}
4   O2s u0 {1,S}
""",
    thermo = 'Cds-CdsCsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 318,
    label = "Cds-CddOsOs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D}
3   O2s u0 {1,S}
4   O2s u0 {1,S}
""",
    thermo = 'Cds-(Cdd-Cd)OsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 319,
    label = "Cds-(Cdd-O2d)OsOs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   O2s u0 {1,S}
4   O2s u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([11.56,15.58,17.69,18.67,18.78,18.4,18.01],'cal/(mol*K)','+|-',[0.07,0.07,0.07,0.07,0.07,0.07,0.07]),
        H298 = (2.403,'kcal/mol','+|-',0.19),
        S298 = (13.42,'cal/(mol*K)','+|-',0.1),
    ),
    shortDesc = """{CCO/O2} RAMAN & GREEN JPCA 2002, 106, 7937-7949""",
    longDesc = 
"""

""",
)

entry(
    index = 320,
    label = "Cds-(Cdd-Cd)OsOs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   O2s u0 {1,S}
4   O2s u0 {1,S}
5   C   u0 {2,D}
""",
    thermo = 'Cds-CdsOsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 321,
    label = "Cds-CdSsSs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   C   u0 {1,D}
3   S2s u0 {1,S}
4   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 322,
    label = "Cds-CdsSsSs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cd  u0 {1,D}
3   S2s u0 {1,S}
4   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 323,
    label = "Cds-CddSsSs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D}
3   S2s u0 {1,S}
4   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 324,
    label = "Cds-(Cdd-S2d)SsSs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   S2s u0 {1,S}
4   S2s u0 {1,S}
5   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 325,
    label = "Cds-(Cdd-Cd)SsSs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   S2s u0 {1,S}
4   S2s u0 {1,S}
5   C   u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 326,
    label = "Cds-CdCH",
    group = 
"""
1 * Cd u0 {2,D} {3,S} {4,S}
2   C  u0 {1,D}
3   C  u0 {1,S}
4   H  u0 {1,S}
""",
    thermo = 'Cds-CdsCsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 327,
    label = "Cds-CdsCsH",
    group = 
"""
1 * Cd u0 {2,D} {3,S} {4,S}
2   Cd u0 {1,D}
3   Cs u0 {1,S}
4   H  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Cd-CsH BENSON""",
    longDesc = 
"""

""",
)

entry(
    index = 328,
    label = "Cds-CdsCdsH",
    group = 
"""
1 * Cd      u0 {2,D} {3,S} {4,S}
2   Cd      u0 {1,D}
3   [Cd,CO] u0 {1,S}
4   H       u0 {1,S}
""",
    thermo = 'Cds-Cds(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 329,
    label = "Cd-Cd(CO)H",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   CO  u0 {1,S} {5,D}
3   H   u0 {1,S}
4   Cd  u0 {1,D}
5   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([18.08,21.17,24.43,27.41,32.22,35.73,40.97],'J/(mol*K)'),
        H298 = (36.4,'kJ/mol'),
        S298 = (33.51,'J/(mol*K)'),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 330,
    label = "Cd-Cd(CO-O)H",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   CO  u0 {1,S} {5,D} {6,S}
3   H   u0 {1,S}
4   Cd  u0 {1,D}
5   O2d u0 {2,D}
6   O2s u0 {2,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Dummy""",
    longDesc = 
"""
Set to zero to avoid double counting with Cds-O2d(Cds-Cds)O2s
""",
)

entry(
    index = 331,
    label = "Cds-Cds(Cds-Cd)H",
    group = 
"""
1 * Cd u0 {2,S} {3,D} {4,S}
2   Cd u0 {1,S} {5,D}
3   Cd u0 {1,D}
4   H  u0 {1,S}
5   C  u0 {2,D}
""",
    thermo = 'Cds-Cds(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 332,
    label = "Cds-Cds(Cds-Cds)H",
    group = 
"""
1 * Cd u0 {2,S} {3,D} {4,S}
2   Cd u0 {1,S} {5,D}
3   Cd u0 {1,D}
4   H  u0 {1,S}
5   Cd u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.46,5.79,6.75,7.42,8.35,8.99,9.98],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (6.78,'kcal/mol','+|-',0.2),
        S298 = (6.38,'cal/(mol*K)','+|-',0.1),
    ),
    shortDesc = """Cd-CdH BENSON""",
    longDesc = 
"""

""",
)

entry(
    index = 333,
    label = "Cds-Cds(Cds-Cdd)H",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,D}
4   H   u0 {1,S}
5   Cdd u0 {2,D}
""",
    thermo = 'Cds-Cds(Cds-Cdd-Cd)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 334,
    label = "Cd-Cd(CCO)H",
    group = 
"""
1 * Cd  u0 {2,S} {4,S} {5,D}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {6,D}
4   H   u0 {1,S}
5   Cd  u0 {1,D}
6   O2d u0 {3,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([18.08,21.17,24.43,27.41,32.22,35.73,40.97],'J/(mol*K)'),
        H298 = (36.4,'kJ/mol'),
        S298 = (33.51,'J/(mol*K)'),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 335,
    label = "Cds-Cds(Cds-Cdd-S2d)H",
    group = 
"""
1 * Cd  u0 {2,S} {4,D} {5,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {6,D}
4   Cd  u0 {1,D}
5   H   u0 {1,S}
6   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 336,
    label = "Cds-Cds(Cds-Cdd-Cd)H",
    group = 
"""
1 * Cd  u0 {2,S} {4,D} {5,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {6,D}
4   Cd  u0 {1,D}
5   H   u0 {1,S}
6   C   u0 {3,D}
""",
    thermo = 'Cds-Cds(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 337,
    label = "Cds-CdsCtH",
    group = 
"""
1 * Cd u0 {2,D} {3,S} {4,S}
2   Cd u0 {1,D}
3   Ct u0 {1,S}
4   H  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.46,5.79,6.75,7.42,8.35,8.99,9.98],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (6.78,'kcal/mol','+|-',0.2),
        S298 = (6.38,'cal/(mol*K)','+|-',0.1),
    ),
    shortDesc = """Cd-CtH BENSON""",
    longDesc = 
"""

""",
)

entry(
    index = 338,
    label = "Cds-CdsCbH",
    group = 
"""
1 * Cd u0 {2,D} {3,S} {4,S}
2   Cd u0 {1,D}
3   Cb u0 {1,S}
4   H  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.46,5.79,6.75,7.42,8.35,8.99,9.98],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (6.78,'kcal/mol','+|-',0.2),
        S298 = (6.38,'cal/(mol*K)','+|-',0.1),
    ),
    shortDesc = """Cd-CbH BENSON""",
    longDesc = 
"""

""",
)

entry(
    index = 339,
    label = "Cds-(Cds-Os)CbH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cd  u0 {1,D} {5,S}
3   Cb  u0 {1,S}
4   H   u0 {1,S}
5   O2s u0 {2,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.28,6.83,7.245,7.264,8.226,9.901,10.176],'cal/(mol*K)'),
        H298 = (10.329,'kcal/mol'),
        S298 = (2.958,'cal/(mol*K)'),
    ),
    shortDesc = """CBS-QB3""",
    longDesc = 
"""
Fitted to CBS-QB3 calculations for OC=Cc1ccccc1
""",
)

entry(
    index = 340,
    label = "Cds-CddCsH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D}
3   Cs  u0 {1,S}
4   H   u0 {1,S}
""",
    thermo = 'Cds-(Cdd-Cd)CsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 341,
    label = "Cds-(Cdd-O2d)CsH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cs  u0 {1,S}
4   H   u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([43.83,50.1,55.5,60.05,67.09,72.13,79.55],'J/(mol*K)','+|-',[4,4,4,4,4,4,4]),
        H298 = (-17.6,'kJ/mol','+|-',3.41),
        S298 = (169.15,'J/(mol*K)','+|-',4.67),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 342,
    label = "Cds-(Cdd-S2d)CsH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cs  u0 {1,S}
4   H   u0 {1,S}
5   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 343,
    label = "Cds-(Cdd-Cd)CsH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cs  u0 {1,S}
4   H   u0 {1,S}
5   C   u0 {2,D}
""",
    thermo = 'Cds-CdsCsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 344,
    label = "Cds-CddCdsH",
    group = 
"""
1 * Cd      u0 {2,D} {3,S} {4,S}
2   Cdd     u0 {1,D}
3   [Cd,CO] u0 {1,S}
4   H       u0 {1,S}
""",
    thermo = 'Cds-(Cdd-Cd)(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 345,
    label = "Cds-(Cdd-O2d)(Cds-O2d)H",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   CO  u0 {1,S} {6,D}
4   H   u0 {1,S}
5   O2d u0 {2,D}
6   O2d u0 {3,D}
""",
    thermo = 'Cds-(Cdd-O2d)CsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 346,
    label = "Cds-(Cdd-O2d)(Cds-Cd)H",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   H   u0 {1,S}
5   O2d u0 {2,D}
6   C   u0 {3,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([43.67,52.95,59.65,64.67,71.81,76.72,83.92],'J/(mol*K)','+|-',[5.66,5.66,5.66,5.66,5.66,5.66,5.66]),
        H298 = (-36,'kJ/mol','+|-',4.82),
        S298 = (152.19,'J/(mol*K)','+|-',6.6),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 347,
    label = "Cds-(Cdd-O2d)(Cds-Cds)H",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   H   u0 {1,S}
5   O2d u0 {2,D}
6   Cd  u0 {3,D}
""",
    thermo = 'Cds-(Cdd-O2d)CsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 348,
    label = "Cds-(Cdd-O2d)(Cds-Cdd)H",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   H   u0 {1,S}
5   O2d u0 {2,D}
6   Cdd u0 {3,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cdd-Cd)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 349,
    label = "Cds-(Cdd-O2d)(Cds-Cdd-O2d)H",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   H   u0 {1,S}
6   O2d u0 {3,D}
7   O2d u0 {4,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([10.55,12.41,13.82,14.91,16.51,17.62,19.24],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (-4.998,'kcal/mol','+|-',0.2),
        S298 = (39.06,'cal/(mol*K)','+|-',0.1),
    ),
    shortDesc = """{CCO/H/CCO} RAMAN & GREEN JPCA 2002, 106, 7937-7949""",
    longDesc = 
"""

""",
)

entry(
    index = 350,
    label = "Cds-(Cdd-O2d)(Cds-Cdd-Cd)H",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   H   u0 {1,S}
6   O2d u0 {3,D}
7   C   u0 {4,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 351,
    label = "Cds-(Cdd-S2d)(Cds-Cd)H",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   H   u0 {1,S}
5   S2d u0 {2,D}
6   C   u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 352,
    label = "Cds-(Cdd-S2d)(Cds-Cds)H",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   H   u0 {1,S}
5   S2d u0 {2,D}
6   Cd  u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 353,
    label = "Cds-(Cdd-S2d)(Cds-Cdd)H",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   H   u0 {1,S}
5   S2d u0 {2,D}
6   Cdd u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 354,
    label = "Cds-(Cdd-S2d)(Cds-Cdd-S2d)H",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   H   u0 {1,S}
6   S2d u0 {3,D}
7   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 355,
    label = "Cds-(Cdd-S2d)(Cds-Cdd-Cd)H",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   H   u0 {1,S}
6   S2d u0 {3,D}
7   C   u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 356,
    label = "Cds-(Cdd-Cd)(Cds-O2d)H",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   CO  u0 {1,S} {6,D}
4   H   u0 {1,S}
5   C   u0 {2,D}
6   O2d u0 {3,D}
""",
    thermo = 'Cd-Cd(CO)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 357,
    label = "Cds-(Cdd-Cd)(Cds-Cd)H",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   H   u0 {1,S}
5   C   u0 {2,D}
6   C   u0 {3,D}
""",
    thermo = 'Cds-(Cdd-Cd)(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 358,
    label = "Cds-(Cdd-Cd)(Cds-Cds)H",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   H   u0 {1,S}
5   C   u0 {2,D}
6   Cd  u0 {3,D}
""",
    thermo = 'Cds-Cds(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 359,
    label = "Cds-(Cdd-Cd)(Cds-Cdd)H",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   H   u0 {1,S}
5   C   u0 {2,D}
6   Cdd u0 {3,D}
""",
    thermo = 'Cds-(Cdd-Cd)(Cds-Cdd-Cd)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 360,
    label = "Cds-(Cdd-Cd)(Cds-Cdd-O2d)H",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   H   u0 {1,S}
6   C   u0 {3,D}
7   O2d u0 {4,D}
""",
    thermo = 'Cd-Cd(CCO)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 361,
    label = "Cds-(Cdd-Cd)(Cds-Cdd-S2d)H",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   H   u0 {1,S}
6   C   u0 {3,D}
7   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 362,
    label = "Cds-(Cdd-Cd)(Cds-Cdd-Cd)H",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   H   u0 {1,S}
6   C   u0 {3,D}
7   C   u0 {4,D}
""",
    thermo = 'Cds-(Cdd-Cd)(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 363,
    label = "Cds-CddCtH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D}
3   Ct  u0 {1,S}
4   H   u0 {1,S}
""",
    thermo = 'Cds-(Cdd-Cd)CtH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 364,
    label = "Cds-(Cdd-O2d)CtH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Ct  u0 {1,S}
4   H   u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 365,
    label = "Cds-(Cdd-S2d)CtH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Ct  u0 {1,S}
4   H   u0 {1,S}
5   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 366,
    label = "Cds-(Cdd-Cd)CtH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Ct  u0 {1,S}
4   H   u0 {1,S}
5   C   u0 {2,D}
""",
    thermo = 'Cds-CdsCtH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 367,
    label = "Cds-CddCbH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D}
3   Cb  u0 {1,S}
4   H   u0 {1,S}
""",
    thermo = 'Cds-(Cdd-Cd)CbH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 368,
    label = "Cds-(Cdd-O2d)CbH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cb  u0 {1,S}
4   H   u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 369,
    label = "Cds-(Cdd-S2d)CbH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cb  u0 {1,S}
4   H   u0 {1,S}
5   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 370,
    label = "Cds-(Cdd-Cd)CbH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cb  u0 {1,S}
4   H   u0 {1,S}
5   C   u0 {2,D}
""",
    thermo = 'Cds-CdsCbH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 371,
    label = "Cds-(Cdd-Cd)C=SH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   CS  u0 {1,S} {6,D}
4   H   u0 {1,S}
5   C   u0 {2,D}
6   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 372,
    label = "Cds-(Cdd-S2d)C=SH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   CS  u0 {1,S} {6,D}
4   H   u0 {1,S}
5   S2d u0 {2,D}
6   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 373,
    label = "Cds-CdsC=SH",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {4,S}
2   CS  u0 {1,S} {5,D}
3   Cd  u0 {1,D}
4   H   u0 {1,S}
5   S2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.4,5.23,6.07,6.84,8.03,8.83,9.89],'cal/(mol*K)'),
        H298 = (7.8,'kcal/mol'),
        S298 = (8.3,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 374,
    label = "Cds-CdCO",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   C   u0 {1,D}
3   C   u0 {1,S}
4   O2s u0 {1,S}
""",
    thermo = 'Cds-CdsCsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 375,
    label = "Cds-CdsCdsOs",
    group = 
"""
1 * Cd      u0 {2,D} {3,S} {4,S}
2   Cd      u0 {1,D}
3   [Cd,CO] u0 {1,S}
4   O2s     u0 {1,S}
""",
    thermo = 'Cds-Cds(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 376,
    label = "Cds-Cds(Cds-O2d)O2s",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {4,S}
2   CO  u0 {1,S} {5,D}
3   Cd  u0 {1,D}
4   O2s u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.4,5.37,5.93,6.18,6.5,6.62,6.72],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (5.13,'kcal/mol','+|-',0.2),
        S298 = (-14.6,'cal/(mol*K)','+|-',0.1),
    ),
    shortDesc = """Cd-OCO adj BENSON for RADOM c*coh""",
    longDesc = 
"""

""",
)

entry(
    index = 377,
    label = "Cds-Cds(Cds-Cd)O2s",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,D}
4   O2s u0 {1,S}
5   C   u0 {2,D}
""",
    thermo = 'Cds-Cds(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 378,
    label = "Cds-Cds(Cds-Cds)O2s",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,D}
4   O2s u0 {1,S}
5   Cd  u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.4,5.37,5.93,6.18,6.5,6.62,6.72],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (1.5,'kcal/mol','+|-',0.2),
        S298 = (-14.4,'cal/(mol*K)','+|-',0.1),
    ),
    shortDesc = """Cd-OCd jwb need calc""",
    longDesc = 
"""

""",
)

entry(
    index = 379,
    label = "Cds-Cds(Cds-Cdd)O2s",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,D}
4   O2s u0 {1,S}
5   Cdd u0 {2,D}
""",
    thermo = 'Cds-Cds(Cds-Cdd-Cd)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 380,
    label = "Cds-Cds(Cds-Cdd-O2d)O2s",
    group = 
"""
1 * Cd  u0 {2,S} {4,D} {5,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {6,D}
4   Cd  u0 {1,D}
5   O2s u0 {1,S}
6   O2d u0 {3,D}
""",
    thermo = 'Cds-Cds(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 381,
    label = "Cds-Cds(Cds-Cdd-Cd)O2s",
    group = 
"""
1 * Cd  u0 {2,S} {4,D} {5,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {6,D}
4   Cd  u0 {1,D}
5   O2s u0 {1,S}
6   C   u0 {3,D}
""",
    thermo = 'Cds-Cds(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 382,
    label = "Cds-CdsCtOs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cd  u0 {1,D}
3   Ct  u0 {1,S}
4   O2s u0 {1,S}
""",
    thermo = 'Cds-Cds(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 383,
    label = "Cds-CdsCbOs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cd  u0 {1,D}
3   Cb  u0 {1,S}
4   O2s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.4,5.37,5.93,6.18,6.5,6.62,6.72],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (1.5,'kcal/mol','+|-',0.2),
        S298 = (-14.4,'cal/(mol*K)','+|-',0.1),
    ),
    shortDesc = """Cd-OCb jwb need calc""",
    longDesc = 
"""

""",
)

entry(
    index = 384,
    label = "Cds-CddCdsOs",
    group = 
"""
1 * Cd      u0 {2,D} {3,S} {4,S}
2   Cdd     u0 {1,D}
3   [Cd,CO] u0 {1,S}
4   O2s     u0 {1,S}
""",
    thermo = 'Cds-(Cdd-Cd)(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 385,
    label = "Cds-(Cdd-O2d)(Cds-O2d)O2s",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   CO  u0 {1,S} {6,D}
4   O2s u0 {1,S}
5   O2d u0 {2,D}
6   O2d u0 {3,D}
""",
    thermo = 'Cds-(Cdd-O2d)CsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 386,
    label = "Cds-(Cdd-O2d)(Cds-Cd)O2s",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   O2s u0 {1,S}
5   O2d u0 {2,D}
6   C   u0 {3,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 387,
    label = "Cds-(Cdd-O2d)(Cds-Cds)O2s",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   O2s u0 {1,S}
5   O2d u0 {2,D}
6   Cd  u0 {3,D}
""",
    thermo = 'Cds-(Cdd-O2d)CsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 388,
    label = "Cds-(Cdd-O2d)(Cds-Cdd)O2s",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   O2s u0 {1,S}
5   O2d u0 {2,D}
6   Cdd u0 {3,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cdd-Cd)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 389,
    label = "Cds-(Cdd-O2d)(Cds-Cdd-O2d)O2s",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   O2s u0 {1,S}
6   O2d u0 {3,D}
7   O2d u0 {4,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([11.01,12.97,14.17,14.97,15.8,16.26,16.88],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (1.607,'kcal/mol','+|-',0.2),
        S298 = (17.73,'cal/(mol*K)','+|-',0.1),
    ),
    shortDesc = """{CCO/O/CCO} RAMAN & GREEN JPCA 2002, 106, 7937-7949""",
    longDesc = 
"""

""",
)

entry(
    index = 390,
    label = "Cds-(Cdd-O2d)(Cds-Cdd-Cd)O2s",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   O2s u0 {1,S}
6   O2d u0 {3,D}
7   C   u0 {4,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 391,
    label = "Cds-(Cdd-Cd)(Cds-Cd)O2s",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   O2s u0 {1,S}
5   C   u0 {2,D}
6   C   u0 {3,D}
""",
    thermo = 'Cds-(Cdd-Cd)(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 392,
    label = "Cds-(Cdd-Cd)(Cds-Cds)O2s",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   O2s u0 {1,S}
5   C   u0 {2,D}
6   Cd  u0 {3,D}
""",
    thermo = 'Cds-Cds(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 393,
    label = "Cds-(Cdd-Cd)(Cds-Cdd)O2s",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   O2s u0 {1,S}
5   C   u0 {2,D}
6   Cdd u0 {3,D}
""",
    thermo = 'Cds-(Cdd-Cd)(Cds-Cdd-Cd)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 394,
    label = "Cds-(Cdd-Cd)(Cds-Cdd-O2d)O2s",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   O2s u0 {1,S}
6   C   u0 {3,D}
7   O2d u0 {4,D}
""",
    thermo = 'Cds-Cds(Cds-Cdd-O2d)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 395,
    label = "Cds-(Cdd-Cd)(Cds-Cdd-Cd)O2s",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   O2s u0 {1,S}
6   C   u0 {3,D}
7   C   u0 {4,D}
""",
    thermo = 'Cds-(Cdd-Cd)(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 396,
    label = "Cds-CddCtOs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D}
3   Ct  u0 {1,S}
4   O2s u0 {1,S}
""",
    thermo = 'Cds-(Cdd-Cd)CtOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 397,
    label = "Cds-(Cdd-O2d)CtOs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Ct  u0 {1,S}
4   O2s u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 398,
    label = "Cds-(Cdd-Cd)CtOs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Ct  u0 {1,S}
4   O2s u0 {1,S}
5   C   u0 {2,D}
""",
    thermo = 'Cds-CdsCtOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 399,
    label = "Cds-CddCbOs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D}
3   Cb  u0 {1,S}
4   O2s u0 {1,S}
""",
    thermo = 'Cds-(Cdd-Cd)CbOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 400,
    label = "Cds-(Cdd-O2d)CbOs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cb  u0 {1,S}
4   O2s u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 401,
    label = "Cds-(Cdd-Cd)CbOs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cb  u0 {1,S}
4   O2s u0 {1,S}
5   C   u0 {2,D}
""",
    thermo = 'Cds-CdsCbOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 402,
    label = "Cd-CdCsOs",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   Cs  u0 {1,S}
3   O2s u0 {1,S}
4   C   u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([12.79,15.86,19.67,22.91,26.55,27.85,28.45],'J/(mol*K)','+|-',[5.1,5.1,5.1,5.1,5.1,5.1,5.1]),
        H298 = (33,'kJ/mol','+|-',4.34),
        S298 = (-50.89,'J/(mol*K)','+|-',5.94),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 403,
    label = "Cds-CdsCsOs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cd  u0 {1,D}
3   Cs  u0 {1,S}
4   O2s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.59,4.56,5.04,5.3,5.84,6.07,6.16],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (3.03,'kcal/mol','+|-',0.2),
        S298 = (-12.32,'cal/(mol*K)','+|-',0.1),
    ),
    shortDesc = """Cd-OCs BOZZELLI-RADOM vin-oh and del (ccoh-ccohc)""",
    longDesc = 
"""

""",
)

entry(
    index = 404,
    label = "Cds-CddCsOs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D}
3   Cs  u0 {1,S}
4   O2s u0 {1,S}
""",
    thermo = 'Cds-(Cdd-Cd)CsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 405,
    label = "Cds-(Cdd-O2d)CsOs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cs  u0 {1,S}
4   O2s u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([10.91,12.65,13.59,14.22,15,15.48,16.28],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (3.273,'kcal/mol','+|-',0.2),
        S298 = (18.58,'cal/(mol*K)','+|-',0.1),
    ),
    shortDesc = """{CCO/O/C} RAMAN & GREEN JPCA 2002, 106, 7937-7949""",
    longDesc = 
"""

""",
)

entry(
    index = 406,
    label = "Cds-(Cdd-Cd)CsOs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cs  u0 {1,S}
4   O2s u0 {1,S}
5   C   u0 {2,D}
""",
    thermo = 'Cds-CdsCsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 407,
    label = "Cds-CdCS",
    group = 
"""
1 * Cd u0 {2,D} {3,S} {4,S}
2   C  u0 {1,D}
3   C  u0 {1,S}
4   S  u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 408,
    label = "Cds-CdsCsSs",
    group = 
"""
1 * Cd u0 {2,D} {3,S} {4,S}
2   Cd u0 {1,D}
3   Cs u0 {1,S}
4   S  u0 {1,S}
""",
    thermo = 'Cds-CdsCsS2',
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2018""",
    longDesc = 
"""
"
""",
)

entry(
    index = 409,
    label = "Cds-CdsCsS2",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cd  u0 {1,D}
3   Cs  u0 {1,S}
4   S2s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.3,5.84,5.62,5.57,6.97,7.7,6.11],'cal/(mol*K)'),
        H298 = (20.57,'kcal/mol'),
        S298 = (-8.4,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 410,
    label = "Cds-CdsCsS4",
    group = 
"""
1 * Cd                u0 {2,D} {3,S} {4,S}
2   Cd                u0 {1,D}
3   Cs                u0 {1,S}
4   [S4s,S4d,S4b,S4t] u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.74,5.6,4.93,5.57,7.72,8.79,7.32],'cal/(mol*K)'),
        H298 = (27.53,'kcal/mol'),
        S298 = (-12.89,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 411,
    label = "Cds-CdsCsS6",
    group = 
"""
1 * Cd                      u0 {2,D} {3,S} {4,S}
2   Cd                      u0 {1,D}
3   Cs                      u0 {1,S}
4   [S6s,S6d,S6dd,S6t,S6td] u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.01,5.41,5.84,6.2,7.7,8.36,7.24],'cal/(mol*K)'),
        H298 = (21.44,'kcal/mol'),
        S298 = (-2.61,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 412,
    label = "Cds-CdsCdsSs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cd  u0 {1,D}
3   Cd  u0 {1,S}
4   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 413,
    label = "Cds-Cds(Cds-Cd)S2s",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,D}
4   S2s u0 {1,S}
5   C   u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 414,
    label = "Cds-Cds(Cds-Cds)S2s",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,D}
4   S2s u0 {1,S}
5   Cd  u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 415,
    label = "Cds-Cds(Cds-Cdd)S2s",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,D}
4   S2s u0 {1,S}
5   Cdd u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 416,
    label = "Cds-Cds(Cds-Cdd-S2d)S2s",
    group = 
"""
1 * Cd  u0 {2,S} {4,D} {5,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {6,D}
4   Cd  u0 {1,D}
5   S2s u0 {1,S}
6   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 417,
    label = "Cds-Cds(Cds-Cdd-Cd)S2s",
    group = 
"""
1 * Cd  u0 {2,S} {4,D} {5,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {6,D}
4   Cd  u0 {1,D}
5   S2s u0 {1,S}
6   C   u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 418,
    label = "Cds-CdsCtSs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cd  u0 {1,D}
3   Ct  u0 {1,S}
4   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 419,
    label = "Cds-CdsCbSs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cd  u0 {1,D}
3   Cb  u0 {1,S}
4   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 420,
    label = "Cds-CddCsSs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D}
3   Cs  u0 {1,S}
4   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 421,
    label = "Cds-(Cdd-S2d)CsSs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cs  u0 {1,S}
4   S2s u0 {1,S}
5   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 422,
    label = "Cds-(Cdd-Cd)CsSs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cs  u0 {1,S}
4   S2s u0 {1,S}
5   C   u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 423,
    label = "Cds-CddCdsSs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D}
3   Cd  u0 {1,S}
4   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 424,
    label = "Cds-(Cdd-S2d)(Cds-Cd)S2s",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   S2s u0 {1,S}
5   S2d u0 {2,D}
6   C   u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 425,
    label = "Cds-(Cdd-S2d)(Cds-Cds)S2s",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   S2s u0 {1,S}
5   S2d u0 {2,D}
6   Cd  u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 426,
    label = "Cds-(Cdd-S2d)(Cds-Cdd)S2s",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   S2s u0 {1,S}
5   S2d u0 {2,D}
6   Cdd u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 427,
    label = "Cds-(Cdd-S2d)(Cds-Cdd-S2d)S2s",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   S2s u0 {1,S}
6   S2d u0 {3,D}
7   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 428,
    label = "Cds-(Cdd-S2d)(Cds-Cdd-Cd)S2s",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   S2s u0 {1,S}
6   S2d u0 {3,D}
7   C   u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 429,
    label = "Cds-(Cdd-Cd)(Cds-Cd)S2s",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   S2s u0 {1,S}
5   C   u0 {2,D}
6   C   u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 430,
    label = "Cds-(Cdd-Cd)(Cds-Cds)S2s",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   S2s u0 {1,S}
5   C   u0 {2,D}
6   Cd  u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 431,
    label = "Cds-(Cdd-Cd)(Cds-Cdd)S2s",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   S2s u0 {1,S}
5   C   u0 {2,D}
6   Cdd u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 432,
    label = "Cds-(Cdd-Cd)(Cds-Cdd-S2d)S2s",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   S2s u0 {1,S}
6   C   u0 {3,D}
7   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 433,
    label = "Cds-(Cdd-Cd)(Cds-Cdd-Cd)S2s",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   S2s u0 {1,S}
6   C   u0 {3,D}
7   C   u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 434,
    label = "Cds-CddCtSs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D}
3   Ct  u0 {1,S}
4   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 435,
    label = "Cds-(Cdd-S2d)CtSs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Ct  u0 {1,S}
4   S2s u0 {1,S}
5   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 436,
    label = "Cds-(Cdd-Cd)CtSs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Ct  u0 {1,S}
4   S2s u0 {1,S}
5   C   u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 437,
    label = "Cds-CddCbSs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D}
3   Cb  u0 {1,S}
4   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 438,
    label = "Cds-(Cdd-S2d)CbSs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cb  u0 {1,S}
4   S2s u0 {1,S}
5   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 439,
    label = "Cds-(Cdd-Cd)CbSs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cb  u0 {1,S}
4   S2s u0 {1,S}
5   C   u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 440,
    label = "Cds-(Cdd-S2d)C=SSs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   CS  u0 {1,S} {6,D}
4   S2s u0 {1,S}
5   S2d u0 {2,D}
6   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 441,
    label = "Cds-CdsC=SSs",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {4,S}
2   CS  u0 {1,S} {5,D}
3   Cd  u0 {1,D}
4   S2s u0 {1,S}
5   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 442,
    label = "Cds-CdCC",
    group = 
"""
1 * Cd u0 {2,D} {3,S} {4,S}
2   C  u0 {1,D}
3   C  u0 {1,S}
4   C  u0 {1,S}
""",
    thermo = 'Cds-CdsCsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 443,
    label = "Cds-CdsCsCs",
    group = 
"""
1 * Cd u0 {2,D} {3,S} {4,S}
2   Cd u0 {1,D}
3   Cs u0 {1,S}
4   Cs u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.1,4.61,4.99,5.26,5.8,6.08,6.36],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (10.34,'kcal/mol','+|-',0.24),
        S298 = (-12.7,'cal/(mol*K)','+|-',0.12),
    ),
    shortDesc = """Cd-CsCs BENSON""",
    longDesc = 
"""

""",
)

entry(
    index = 444,
    label = "Cds-CdsCdsCs",
    group = 
"""
1 * Cd      u0 {2,D} {3,S} {4,S}
2   Cd      u0 {1,D}
3   [Cd,CO] u0 {1,S}
4   Cs      u0 {1,S}
""",
    thermo = 'Cds-Cds(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 445,
    label = "Cd-CdCs(CO)",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   CO  u0 {1,S} {5,D}
3   Cs  u0 {1,S}
4   Cd  u0 {1,D}
5   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([15.33,16.82,18.64,20.42,23.2,25,27.1],'J/(mol*K)','+|-',[5.66,5.66,5.66,5.66,5.66,5.66,5.66]),
        H298 = (39,'kJ/mol','+|-',4.82),
        S298 = (-51.26,'J/(mol*K)','+|-',6.6),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 446,
    label = "Cds-Cds(Cds-Cd)Cs",
    group = 
"""
1 * Cd u0 {2,S} {3,D} {4,S}
2   Cd u0 {1,S} {5,D}
3   Cd u0 {1,D}
4   Cs u0 {1,S}
5   C  u0 {2,D}
""",
    thermo = 'Cds-Cds(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 447,
    label = "Cds-Cds(Cds-Cds)Cs",
    group = 
"""
1 * Cd u0 {2,S} {3,D} {4,S}
2   Cd u0 {1,S} {5,D}
3   Cd u0 {1,D}
4   Cs u0 {1,S}
5   Cd u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.4,5.37,5.93,6.18,6.5,6.62,6.72],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (8.88,'kcal/mol','+|-',0.24),
        S298 = (-14.6,'cal/(mol*K)','+|-',0.12),
    ),
    shortDesc = """Cd-CdCs BENSON""",
    longDesc = 
"""

""",
)

entry(
    index = 448,
    label = "Cds-Cds(Cds-Cdd)Cs",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,D}
4   Cs  u0 {1,S}
5   Cdd u0 {2,D}
""",
    thermo = 'Cds-Cds(Cds-Cdd-Cd)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 449,
    label = "Cd-CdCs(CCO)",
    group = 
"""
1 * Cd  u0 {2,S} {4,S} {5,D}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {6,D}
4   Cs  u0 {1,S}
5   Cd  u0 {1,D}
6   O2d u0 {3,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([22.68,24.05,24.63,25.07,25.64,25.84,25.7],'J/(mol*K)','+|-',[8,8,8,8,8,8,8]),
        H298 = (41.6,'kJ/mol','+|-',6.82),
        S298 = (-48.01,'J/(mol*K)','+|-',9.33),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 450,
    label = "Cds-Cds(Cds-Cdd-S2d)Cs",
    group = 
"""
1 * Cd  u0 {2,S} {4,D} {5,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {6,D}
4   Cd  u0 {1,D}
5   Cs  u0 {1,S}
6   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 451,
    label = "Cds-Cds(Cds-Cdd-Cd)Cs",
    group = 
"""
1 * Cd  u0 {2,S} {4,D} {5,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {6,D}
4   Cd  u0 {1,D}
5   Cs  u0 {1,S}
6   C   u0 {3,D}
""",
    thermo = 'Cds-Cds(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 452,
    label = "Cds-CdsCdsCds",
    group = 
"""
1 * Cd      u0 {2,D} {3,S} {4,S}
2   Cd      u0 {1,D}
3   [Cd,CO] u0 {1,S}
4   [Cd,CO] u0 {1,S}
""",
    thermo = 'Cds-Cds(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 453,
    label = "Cds-Cds(Cds-O2d)(Cds-O2d)",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   CO  u0 {1,S} {5,D}
3   CO  u0 {1,S} {6,D}
4   Cd  u0 {1,D}
5   O2d u0 {2,D}
6   O2d u0 {3,D}
""",
    thermo = 'Cds-CdsCsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 454,
    label = "Cds-Cds(Cds-O2d)(Cds-Cd)",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   CO  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {5,D}
4   Cd  u0 {1,D}
5   C   u0 {3,D}
6   O2d u0 {2,D}
""",
    thermo = 'Cds-Cds(Cds-O2d)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 455,
    label = "Cds-Cds(Cds-O2d)(Cds-Cds)",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   CO  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {5,D}
4   Cd  u0 {1,D}
5   Cd  u0 {3,D}
6   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.7,6.13,6.87,7.1,7.2,7.16,7.06],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (11.6,'kcal/mol','+|-',0.24),
        S298 = (-16.5,'cal/(mol*K)','+|-',0.12),
    ),
    shortDesc = """Cd-COCd from CD/CD2/ jwb est 6/97""",
    longDesc = 
"""
AG Vandeputte, added 7 kcal/mol to the following value (see phd M Sabbe)
""",
)

entry(
    index = 456,
    label = "Cds-Cds(Cds-O2d)(Cds-Cdd)",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   CO  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {5,D}
4   Cd  u0 {1,D}
5   Cdd u0 {3,D}
6   O2d u0 {2,D}
""",
    thermo = 'Cds-Cds(Cds-O2d)(Cds-Cdd-Cd)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 457,
    label = "Cds-Cds(Cds-O2d)(Cds-Cdd-O2d)",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {5,D}
2   Cd  u0 {1,S} {4,D}
3   CO  u0 {1,S} {6,D}
4   Cdd u0 {2,D} {7,D}
5   Cd  u0 {1,D}
6   O2d u0 {3,D}
7   O2d u0 {4,D}
""",
    thermo = 'Cd-CdCs(CCO)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 458,
    label = "Cds-Cds(Cds-O2d)(Cds-Cdd-Cd)",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {5,D}
2   Cd  u0 {1,S} {4,D}
3   CO  u0 {1,S} {6,D}
4   Cdd u0 {2,D} {7,D}
5   Cd  u0 {1,D}
6   O2d u0 {3,D}
7   C   u0 {4,D}
""",
    thermo = 'Cds-Cds(Cds-O2d)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 459,
    label = "Cds-Cds(Cds-Cd)(Cds-Cd)",
    group = 
"""
1 * Cd u0 {2,S} {3,S} {4,D}
2   Cd u0 {1,S} {5,D}
3   Cd u0 {1,S} {6,D}
4   Cd u0 {1,D}
5   C  u0 {2,D}
6   C  u0 {3,D}
""",
    thermo = 'Cds-Cds(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 460,
    label = "Cds-Cds(Cds-Cds)(Cds-Cds)",
    group = 
"""
1 * Cd u0 {2,S} {3,S} {4,D}
2   Cd u0 {1,S} {5,D}
3   Cd u0 {1,S} {6,D}
4   Cd u0 {1,D}
5   Cd u0 {2,D}
6   Cd u0 {3,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([1.9,2.69,3.5,4.28,5.57,6.21,7.37],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (11.6,'kcal/mol','+|-',0.24),
        S298 = (-15.67,'cal/(mol*K)','+|-',0.12),
    ),
    shortDesc = """Cd-CdCd Hf=3D est S,Cp mopac nov99""",
    longDesc = 
"""
AG Vandeputte, added 7 kcal/mol to the following value (see phd M Sabbe)
""",
)

entry(
    index = 461,
    label = "Cds-Cds(Cds-Cds)(Cds-Cdd)",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,D}
5   Cd  u0 {2,D}
6   Cdd u0 {3,D}
""",
    thermo = 'Cds-Cds(Cds-Cds)(Cds-Cdd-Cd)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 462,
    label = "Cds-Cds(Cds-Cds)(Cds-Cdd-O2d)",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {5,D}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {6,D}
4   Cdd u0 {2,D} {7,D}
5   Cd  u0 {1,D}
6   Cd  u0 {3,D}
7   O2d u0 {4,D}
""",
    thermo = 'Cd-CdCs(CCO)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 463,
    label = "Cds-Cds(Cds-Cds)(Cds-Cdd-S2d)",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {5,D}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {6,D}
4   Cdd u0 {2,D} {7,D}
5   Cd  u0 {1,D}
6   Cd  u0 {3,D}
7   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 464,
    label = "Cds-Cds(Cds-Cds)(Cds-Cdd-Cd)",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {5,D}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {6,D}
4   Cdd u0 {2,D} {7,D}
5   Cd  u0 {1,D}
6   Cd  u0 {3,D}
7   C   u0 {4,D}
""",
    thermo = 'Cds-Cds(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 465,
    label = "Cds-Cds(Cds-Cdd)(Cds-Cdd)",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,D}
5   Cdd u0 {2,D}
6   Cdd u0 {3,D}
""",
    thermo = 'Cds-Cds(Cds-Cdd-Cd)(Cds-Cdd-Cd)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 466,
    label = "Cds-Cds(Cds-Cdd-O2d)(Cds-Cdd-O2d)",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {6,D}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {7,D}
5   Cdd u0 {3,D} {8,D}
6   Cd  u0 {1,D}
7   O2d u0 {4,D}
8   O2d u0 {5,D}
""",
    thermo = 'Cds-Cds(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 467,
    label = "Cds-Cds(Cds-Cdd-O2d)(Cds-Cdd-Cd)",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {6,D}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {7,D}
5   Cdd u0 {3,D} {8,D}
6   Cd  u0 {1,D}
7   O2d u0 {4,D}
8   C   u0 {5,D}
""",
    thermo = 'Cds-Cds(Cds-Cds)(Cds-Cdd-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 468,
    label = "Cds-Cds(Cds-Cdd-S2d)(Cds-Cdd-S2d)",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {6,D}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {7,D}
5   Cdd u0 {3,D} {8,D}
6   Cd  u0 {1,D}
7   S2d u0 {4,D}
8   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 469,
    label = "Cds-Cds(Cds-Cdd-S2d)(Cds-Cdd-Cd)",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {6,D}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {7,D}
5   Cdd u0 {3,D} {8,D}
6   Cd  u0 {1,D}
7   S2d u0 {4,D}
8   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 470,
    label = "Cds-Cds(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {6,D}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {7,D}
5   Cdd u0 {3,D} {8,D}
6   Cd  u0 {1,D}
7   C   u0 {4,D}
8   C   u0 {5,D}
""",
    thermo = 'Cds-Cds(Cds-Cds)(Cds-Cdd-Cd)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 471,
    label = "Cds-CdsCtCs",
    group = 
"""
1 * Cd u0 {2,D} {3,S} {4,S}
2   Cd u0 {1,D}
3   Ct u0 {1,S}
4   Cs u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.5,3.88,4.88,4.18,4.86,5.4,6.01],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (8.11,'kcal/mol','+|-',0.24),
        S298 = (-13.02,'cal/(mol*K)','+|-',0.12),
    ),
    shortDesc = """Cd-CtCs RAMAN & GREEN JPCA 2002, 106, 11141-11149""",
    longDesc = 
"""

""",
)

entry(
    index = 472,
    label = "Cds-CdsCtCds",
    group = 
"""
1 * Cd      u0 {2,D} {3,S} {4,S}
2   Cd      u0 {1,D}
3   Ct      u0 {1,S}
4   [Cd,CO] u0 {1,S}
""",
    thermo = 'Cds-Cds(Cds-Cds)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 473,
    label = "Cds-CdsCt(Cds-O2d)",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {4,S}
2   CO  u0 {1,S} {5,D}
3   Cd  u0 {1,D}
4   Ct  u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = 'Cds-Cds(Cds-O2d)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 474,
    label = "Cds-CdsCt(Cds-Cd)",
    group = 
"""
1 * Cd u0 {2,S} {3,D} {4,S}
2   Cd u0 {1,S} {5,D}
3   Cd u0 {1,D}
4   Ct u0 {1,S}
5   C  u0 {2,D}
""",
    thermo = 'Cds-Cds(Cds-Cds)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 475,
    label = "Cds-Cds(Cds-Cds)Ct",
    group = 
"""
1 * Cd u0 {2,S} {3,D} {4,S}
2   Cd u0 {1,S} {5,D}
3   Cd u0 {1,D}
4   Ct u0 {1,S}
5   Cd u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.89,5.26,5.98,6.37,6.67,6.78,6.89],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (7.54,'kcal/mol','+|-',0.24),
        S298 = (-14.65,'cal/(mol*K)','+|-',0.12),
    ),
    shortDesc = """Cd-CtCd RAMAN & GREEN JPCA 2002, 106, 11141-11149""",
    longDesc = 
"""

""",
)

entry(
    index = 476,
    label = "Cds-Cds(Cds-Cdd)Ct",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,D}
4   Ct  u0 {1,S}
5   Cdd u0 {2,D}
""",
    thermo = 'Cds-Cds(Cds-Cdd-Cd)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 477,
    label = "Cds-Cds(Cds-Cdd-O2d)Ct",
    group = 
"""
1 * Cd  u0 {2,S} {4,D} {5,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {6,D}
4   Cd  u0 {1,D}
5   Ct  u0 {1,S}
6   O2d u0 {3,D}
""",
    thermo = 'Cds-Cds(Cds-Cds)(Cds-Cdd-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 478,
    label = "Cds-Cds(Cds-Cdd-S2d)Ct",
    group = 
"""
1 * Cd  u0 {2,S} {4,D} {5,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {6,D}
4   Cd  u0 {1,D}
5   Ct  u0 {1,S}
6   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 479,
    label = "Cds-Cds(Cds-Cdd-Cd)Ct",
    group = 
"""
1 * Cd  u0 {2,S} {4,D} {5,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {6,D}
4   Cd  u0 {1,D}
5   Ct  u0 {1,S}
6   C   u0 {3,D}
""",
    thermo = 'Cds-Cds(Cds-Cds)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 480,
    label = "Cds-CdsCtCt",
    group = 
"""
1 * Cd u0 {2,D} {3,S} {4,S}
2   Cd u0 {1,D}
3   Ct u0 {1,S}
4   Ct u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.23,4.59,5.41,5.93,6.48,6.74,7.02],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (8.81,'kcal/mol','+|-',0.24),
        S298 = (-13.51,'cal/(mol*K)','+|-',0.12),
    ),
    shortDesc = """Cd-CtCt RAMAN & GREEN JPCA 2002, 106, 11141-11149""",
    longDesc = 
"""

""",
)

entry(
    index = 481,
    label = "Cds-CdsCbCs",
    group = 
"""
1 * Cd u0 {2,D} {3,S} {4,S}
2   Cd u0 {1,D}
3   Cb u0 {1,S}
4   Cs u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.4,5.37,5.93,6.18,6.5,6.62,6.72],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (8.64,'kcal/mol','+|-',0.24),
        S298 = (-14.6,'cal/(mol*K)','+|-',0.12),
    ),
    shortDesc = """Cd-CbCs BENSON""",
    longDesc = 
"""

""",
)

entry(
    index = 482,
    label = "Cds-CdsCbCds",
    group = 
"""
1 * Cd      u0 {2,D} {3,S} {4,S}
2   Cd      u0 {1,D}
3   Cb      u0 {1,S}
4   [Cd,CO] u0 {1,S}
""",
    thermo = 'Cds-Cds(Cds-Cds)Cb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 483,
    label = "Cds-CdsCb(Cds-O2d)",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {4,S}
2   CO  u0 {1,S} {5,D}
3   Cd  u0 {1,D}
4   Cb  u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = 'Cds-Cds(Cds-O2d)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 484,
    label = "Cds-Cds(Cds-Cd)Cb",
    group = 
"""
1 * Cd u0 {2,S} {3,D} {4,S}
2   Cd u0 {1,S} {5,D}
3   Cd u0 {1,D}
4   Cb u0 {1,S}
5   C  u0 {2,D}
""",
    thermo = 'Cds-Cds(Cds-Cds)Cb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 485,
    label = "Cds-Cds(Cds-Cds)Cb",
    group = 
"""
1 * Cd u0 {2,S} {3,D} {4,S}
2   Cd u0 {1,S} {5,D}
3   Cd u0 {1,D}
4   Cb u0 {1,S}
5   Cd u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.7,6.13,6.87,7.1,7.2,7.16,7.06],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (7.18,'kcal/mol','+|-',0.24),
        S298 = (-16.5,'cal/(mol*K)','+|-',0.12),
    ),
    shortDesc = """Cd-CbCd BOZZELLI =3D Cd/Cs/Cb + (Cd/Cs/Cd - Cd/Cs/Cs)""",
    longDesc = 
"""

""",
)

entry(
    index = 486,
    label = "Cds-Cds(Cds-Cdd)Cb",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,D}
4   Cb  u0 {1,S}
5   Cdd u0 {2,D}
""",
    thermo = 'Cds-Cds(Cds-Cdd-Cd)Cb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 487,
    label = "Cds-Cds(Cds-Cdd-O2d)Cb",
    group = 
"""
1 * Cd  u0 {2,S} {4,D} {5,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {6,D}
4   Cd  u0 {1,D}
5   Cb  u0 {1,S}
6   O2d u0 {3,D}
""",
    thermo = 'Cds-Cds(Cds-Cds)(Cds-Cdd-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 488,
    label = "Cds-Cds(Cds-Cdd-S2d)Cb",
    group = 
"""
1 * Cd  u0 {2,S} {4,D} {5,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {6,D}
4   Cd  u0 {1,D}
5   Cb  u0 {1,S}
6   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 489,
    label = "Cds-Cds(Cds-Cdd-Cd)Cb",
    group = 
"""
1 * Cd  u0 {2,S} {4,D} {5,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {6,D}
4   Cd  u0 {1,D}
5   Cb  u0 {1,S}
6   C   u0 {3,D}
""",
    thermo = 'Cds-Cds(Cds-Cds)Cb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 490,
    label = "Cds-CdsCbCt",
    group = 
"""
1 * Cd u0 {2,D} {3,S} {4,S}
2   Cd u0 {1,D}
3   Cb u0 {1,S}
4   Ct u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([2.22,3.14,4.54,4.11,5.06,5.79,6.71],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (6.7,'kcal/mol','+|-',0.24),
        S298 = (-17.04,'cal/(mol*K)','+|-',0.12),
    ),
    shortDesc = """Cd-CbCt Hf=3D est S,Cp mopac nov99""",
    longDesc = 
"""

""",
)

entry(
    index = 491,
    label = "Cds-CdsCbCb",
    group = 
"""
1 * Cd u0 {2,D} {3,S} {4,S}
2   Cd u0 {1,D}
3   Cb u0 {1,S}
4   Cb u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.7,6.13,6.87,7.1,7.2,7.16,7.06],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (8,'kcal/mol','+|-',0.24),
        S298 = (-16.5,'cal/(mol*K)','+|-',0.12),
    ),
    shortDesc = """Cd-CbCb BOZZELLI =3D Cd/Cs/Cb + (Cd/Cs/Cb - Cd/Cs/Cs)""",
    longDesc = 
"""

""",
)

entry(
    index = 492,
    label = "Cds-CddCsCs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D}
3   Cs  u0 {1,S}
4   Cs  u0 {1,S}
""",
    thermo = 'Cds-(Cdd-Cd)CsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 493,
    label = "Cds-(Cdd-O2d)CsCs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cs  u0 {1,S}
4   Cs  u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([42.55,46.42,50,53.24,58.3,61.71,66.01],'J/(mol*K)','+|-',[4.76,4.76,4.76,4.76,4.76,4.76,4.76]),
        H298 = (0.5,'kJ/mol','+|-',4.06),
        S298 = (84.72,'J/(mol*K)','+|-',5.55),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 494,
    label = "Cds-(Cdd-S2d)CsCs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cs  u0 {1,S}
4   Cs  u0 {1,S}
5   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 495,
    label = "Cds-(Cdd-Cd)CsCs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cs  u0 {1,S}
4   Cs  u0 {1,S}
5   C   u0 {2,D}
""",
    thermo = 'Cds-CdsCsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 496,
    label = "Cds-CddCdsCs",
    group = 
"""
1 * Cd      u0 {2,D} {3,S} {4,S}
2   Cdd     u0 {1,D}
3   [Cd,CO] u0 {1,S}
4   Cs      u0 {1,S}
""",
    thermo = 'Cds-(Cdd-Cd)(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 497,
    label = "Cds-(Cdd-O2d)(Cds-O2d)Cs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   CO  u0 {1,S} {6,D}
4   Cs  u0 {1,S}
5   O2d u0 {2,D}
6   O2d u0 {3,D}
""",
    thermo = 'Cds-(Cdd-O2d)CsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 498,
    label = "Cds-(Cdd-O2d)(Cds-Cd)Cs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cs  u0 {1,S}
5   O2d u0 {2,D}
6   C   u0 {3,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 499,
    label = "Cds-(Cdd-O2d)(Cds-Cds)Cs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cs  u0 {1,S}
5   O2d u0 {2,D}
6   Cd  u0 {3,D}
""",
    thermo = 'Cds-(Cdd-O2d)CsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 500,
    label = "Cds-(Cdd-O2d)(Cds-Cdd)Cs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cs  u0 {1,S}
5   O2d u0 {2,D}
6   Cdd u0 {3,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cdd-Cd)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 501,
    label = "Cds-(Cdd-O2d)(Cds-Cdd-O2d)Cs",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   Cs  u0 {1,S}
6   O2d u0 {3,D}
7   O2d u0 {4,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([10.1,11.24,12.12,12.84,14,14.75,15.72],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (-2.07,'kcal/mol','+|-',0.24),
        S298 = (19.65,'cal/(mol*K)','+|-',0.12),
    ),
    shortDesc = """{CCO/C/CCO} RAMAN & GREEN JPCA 2002, 106, 7937-7949""",
    longDesc = 
"""

""",
)

entry(
    index = 502,
    label = "Cds-(Cdd-O2d)(Cds-Cdd-Cd)Cs",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   Cs  u0 {1,S}
6   O2d u0 {3,D}
7   C   u0 {4,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 503,
    label = "Cds-(Cdd-S2d)(Cds-Cd)Cs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cs  u0 {1,S}
5   S2d u0 {2,D}
6   C   u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 504,
    label = "Cds-(Cdd-S2d)(Cds-Cds)Cs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cs  u0 {1,S}
5   S2d u0 {2,D}
6   Cd  u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 505,
    label = "Cds-(Cdd-S2d)(Cds-Cdd)Cs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cs  u0 {1,S}
5   S2d u0 {2,D}
6   Cdd u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 506,
    label = "Cds-(Cdd-S2d)(Cds-Cdd-S2d)Cs",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   Cs  u0 {1,S}
6   S2d u0 {3,D}
7   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 507,
    label = "Cds-(Cdd-S2d)(Cds-Cdd-Cd)Cs",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   Cs  u0 {1,S}
6   S2d u0 {3,D}
7   C   u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 508,
    label = "Cds-(Cdd-Cd)(Cds-Cd)Cs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cs  u0 {1,S}
5   C   u0 {2,D}
6   C   u0 {3,D}
""",
    thermo = 'Cds-(Cdd-Cd)(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 509,
    label = "Cds-(Cdd-Cd)(Cds-Cds)Cs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cs  u0 {1,S}
5   C   u0 {2,D}
6   Cd  u0 {3,D}
""",
    thermo = 'Cds-Cds(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 510,
    label = "Cds-(Cdd-Cd)(Cds-Cdd)Cs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cs  u0 {1,S}
5   C   u0 {2,D}
6   Cdd u0 {3,D}
""",
    thermo = 'Cds-(Cdd-Cd)(Cds-Cdd-Cd)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 511,
    label = "Cds-(Cdd-Cd)(Cds-Cdd-O2d)Cs",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   Cs  u0 {1,S}
6   C   u0 {3,D}
7   O2d u0 {4,D}
""",
    thermo = 'Cd-CdCs(CCO)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 512,
    label = "Cds-(Cdd-Cd)(Cds-Cdd-S2d)Cs",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   Cs  u0 {1,S}
6   C   u0 {3,D}
7   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 513,
    label = "Cds-(Cdd-Cd)(Cds-Cdd-Cd)Cs",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   Cs  u0 {1,S}
6   C   u0 {3,D}
7   C   u0 {4,D}
""",
    thermo = 'Cds-(Cdd-Cd)(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 514,
    label = "Cds-CddCdsCds",
    group = 
"""
1 * Cd      u0 {2,D} {3,S} {4,S}
2   Cdd     u0 {1,D}
3   [Cd,CO] u0 {1,S}
4   [Cd,CO] u0 {1,S}
""",
    thermo = 'Cds-(Cdd-Cd)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 515,
    label = "Cds-(Cdd-O2d)(Cds-O2d)(Cds-O2d)",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   CO  u0 {1,S} {6,D}
4   CO  u0 {1,S} {7,D}
5   O2d u0 {2,D}
6   O2d u0 {3,D}
7   O2d u0 {4,D}
""",
    thermo = 'Cds-(Cdd-O2d)CsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 516,
    label = "Cds-(Cdd-O2d)(Cds-Cd)(Cds-O2d)",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   CO  u0 {1,S} {7,D}
5   O2d u0 {2,D}
6   C   u0 {3,D}
7   O2d u0 {4,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cds)(Cds-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 517,
    label = "Cds-(Cdd-O2d)(Cds-Cds)(Cds-O2d)",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   CO  u0 {1,S} {7,D}
5   O2d u0 {2,D}
6   Cd  u0 {3,D}
7   O2d u0 {4,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-O2d)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 518,
    label = "Cds-(Cdd-O2d)(Cds-Cdd)(Cds-O2d)",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   CO  u0 {1,S} {7,D}
5   O2d u0 {2,D}
6   Cdd u0 {3,D}
7   O2d u0 {4,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cdd-Cd)(Cds-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 519,
    label = "Cds-(Cdd-O2d)(Cds-Cdd-O2d)(Cds-O2d)",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   Cdd u0 {1,D} {6,D}
4   CO  u0 {1,S} {7,D}
5   Cdd u0 {2,D} {8,D}
6   O2d u0 {3,D}
7   O2d u0 {4,D}
8   O2d u0 {5,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cdd-O2d)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 520,
    label = "Cds-(Cdd-O2d)(Cds-Cdd-Cd)(Cds-O2d)",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   Cdd u0 {1,D} {6,D}
4   CO  u0 {1,S} {7,D}
5   Cdd u0 {2,D} {8,D}
6   O2d u0 {3,D}
7   O2d u0 {4,D}
8   C   u0 {5,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cds)(Cds-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 521,
    label = "Cds-(Cdd-O2d)(Cds-Cd)(Cds-Cd)",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   O2d u0 {2,D}
6   C   u0 {3,D}
7   C   u0 {4,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 522,
    label = "Cds-(Cdd-O2d)(Cds-Cds)(Cds-Cds)",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   O2d u0 {2,D}
6   Cd  u0 {3,D}
7   Cd  u0 {4,D}
""",
    thermo = 'Cds-(Cdd-O2d)CsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 523,
    label = "Cds-(Cdd-O2d)(Cds-Cdd)(Cds-Cds)",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   O2d u0 {2,D}
6   Cdd u0 {3,D}
7   Cd  u0 {4,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cdd-Cd)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 524,
    label = "Cds-(Cdd-O2d)(Cds-Cdd-O2d)(Cds-Cds)",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   Cdd u0 {1,D} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Cdd u0 {2,D} {8,D}
6   O2d u0 {3,D}
7   Cd  u0 {4,D}
8   O2d u0 {5,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cdd-O2d)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 525,
    label = "Cds-(Cdd-O2d)(Cds-Cdd-Cd)(Cds-Cds)",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   Cdd u0 {1,D} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Cdd u0 {2,D} {8,D}
6   O2d u0 {3,D}
7   Cd  u0 {4,D}
8   C   u0 {5,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 526,
    label = "Cds-(Cdd-O2d)(Cds-Cdd)(Cds-Cdd)",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   O2d u0 {2,D}
6   Cdd u0 {3,D}
7   Cdd u0 {4,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cdd-Cd)(Cds-Cdd-Cd)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 527,
    label = "Cds-(Cdd-O2d)(Cds-Cdd-O2d)(Cds-Cdd-O2d)",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cdd u0 {1,D} {7,D}
5   Cdd u0 {2,D} {8,D}
6   Cdd u0 {3,D} {9,D}
7   O2d u0 {4,D}
8   O2d u0 {5,D}
9   O2d u0 {6,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 528,
    label = "Cds-(Cdd-O2d)(Cds-Cdd-O2d)(Cds-Cdd-Cd)",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cdd u0 {1,D} {7,D}
5   Cdd u0 {2,D} {8,D}
6   Cdd u0 {3,D} {9,D}
7   O2d u0 {4,D}
8   O2d u0 {5,D}
9   C   u0 {6,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cdd-O2d)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 529,
    label = "Cds-(Cdd-O2d)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cdd u0 {1,D} {7,D}
5   Cdd u0 {2,D} {8,D}
6   Cdd u0 {3,D} {9,D}
7   O2d u0 {4,D}
8   C   u0 {5,D}
9   C   u0 {6,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 530,
    label = "Cds-(Cdd-Cd)(Cds-O2d)(Cds-O2d)",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   CO  u0 {1,S} {6,D}
4   CO  u0 {1,S} {7,D}
5   C   u0 {2,D}
6   O2d u0 {3,D}
7   O2d u0 {4,D}
""",
    thermo = 'Cds-Cds(Cds-O2d)(Cds-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 531,
    label = "Cds-(Cdd-Cd)(Cds-O2d)(Cds-Cd)",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   CO  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {6,D}
5   C   u0 {2,D}
6   C   u0 {4,D}
7   O2d u0 {3,D}
""",
    thermo = 'Cds-(Cdd-Cd)(Cds-O2d)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 532,
    label = "Cds-(Cdd-Cd)(Cds-O2d)(Cds-Cds)",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   CO  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {6,D}
5   C   u0 {2,D}
6   Cd  u0 {4,D}
7   O2d u0 {3,D}
""",
    thermo = 'Cds-Cds(Cds-O2d)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 533,
    label = "Cds-(Cdd-Cd)(Cds-O2d)(Cds-Cdd)",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   CO  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {6,D}
5   C   u0 {2,D}
6   Cdd u0 {4,D}
7   O2d u0 {3,D}
""",
    thermo = 'Cds-(Cdd-Cd)(Cds-O2d)(Cds-Cdd-Cd)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 534,
    label = "Cds-(Cdd-Cd)(Cds-O2d)(Cds-Cdd-O2d)",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   Cdd u0 {1,D} {6,D}
4   CO  u0 {1,S} {7,D}
5   Cdd u0 {2,D} {8,D}
6   C   u0 {3,D}
7   O2d u0 {4,D}
8   O2d u0 {5,D}
""",
    thermo = 'Cds-Cds(Cds-O2d)(Cds-Cdd-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 535,
    label = "Cds-(Cdd-Cd)(Cds-O2d)(Cds-Cdd-Cd)",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   Cdd u0 {1,D} {6,D}
4   CO  u0 {1,S} {7,D}
5   Cdd u0 {2,D} {8,D}
6   C   u0 {3,D}
7   O2d u0 {4,D}
8   C   u0 {5,D}
""",
    thermo = 'Cds-(Cdd-Cd)(Cds-O2d)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 536,
    label = "Cds-(Cdd-S2d)(Cds-Cd)(Cds-Cd)",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   S2d u0 {2,D}
6   C   u0 {3,D}
7   C   u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 537,
    label = "Cds-(Cdd-S2d)(Cds-Cds)(Cds-Cds)",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   S2d u0 {2,D}
6   Cd  u0 {3,D}
7   Cd  u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 538,
    label = "Cds-(Cdd-S2d)(Cds-Cdd)(Cds-Cds)",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   S2d u0 {2,D}
6   Cdd u0 {3,D}
7   Cd  u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 539,
    label = "Cds-(Cdd-S2d)(Cds-Cdd-S2d)(Cds-Cds)",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   Cdd u0 {1,D} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Cdd u0 {2,D} {8,D}
6   S2d u0 {3,D}
7   Cd  u0 {4,D}
8   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 540,
    label = "Cds-(Cdd-S2d)(Cds-Cdd-Cd)(Cds-Cds)",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   Cdd u0 {1,D} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Cdd u0 {2,D} {8,D}
6   S2d u0 {3,D}
7   Cd  u0 {4,D}
8   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 541,
    label = "Cds-(Cdd-S2d)(Cds-Cdd)(Cds-Cdd)",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   S2d u0 {2,D}
6   Cdd u0 {3,D}
7   Cdd u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 542,
    label = "Cds-(Cdd-S2d)(Cds-Cdd-S2d)(Cds-Cdd-S2d)",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cdd u0 {1,D} {7,D}
5   Cdd u0 {2,D} {8,D}
6   Cdd u0 {3,D} {9,D}
7   S2d u0 {4,D}
8   S2d u0 {5,D}
9   S2d u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 543,
    label = "Cds-(Cdd-S2d)(Cds-Cdd-S2d)(Cds-Cdd-Cd)",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cdd u0 {1,D} {7,D}
5   Cdd u0 {2,D} {8,D}
6   Cdd u0 {3,D} {9,D}
7   S2d u0 {4,D}
8   S2d u0 {5,D}
9   C   u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 544,
    label = "Cds-(Cdd-S2d)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cdd u0 {1,D} {7,D}
5   Cdd u0 {2,D} {8,D}
6   Cdd u0 {3,D} {9,D}
7   S2d u0 {4,D}
8   C   u0 {5,D}
9   C   u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 545,
    label = "Cds-(Cdd-Cd)(Cds-Cd)(Cds-Cd)",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   C   u0 {2,D}
6   C   u0 {3,D}
7   C   u0 {4,D}
""",
    thermo = 'Cds-(Cdd-Cd)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 546,
    label = "Cds-(Cdd-Cd)(Cds-Cds)(Cds-Cds)",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   C   u0 {2,D}
6   Cd  u0 {3,D}
7   Cd  u0 {4,D}
""",
    thermo = 'Cds-Cds(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 547,
    label = "Cds-(Cdd-Cd)(Cds-Cdd)(Cds-Cds)",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   C   u0 {2,D}
6   Cdd u0 {3,D}
7   Cd  u0 {4,D}
""",
    thermo = 'Cds-(Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 548,
    label = "Cds-(Cdd-Cd)(Cds-Cdd-O2d)(Cds-Cds)",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   Cdd u0 {1,D} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Cdd u0 {2,D} {8,D}
6   C   u0 {3,D}
7   Cd  u0 {4,D}
8   O2d u0 {5,D}
""",
    thermo = 'Cds-Cds(Cds-Cds)(Cds-Cdd-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 549,
    label = "Cds-(Cdd-Cd)(Cds-Cdd-S2d)(Cds-Cds)",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   Cdd u0 {1,D} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Cdd u0 {2,D} {8,D}
6   C   u0 {3,D}
7   Cd  u0 {4,D}
8   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 550,
    label = "Cds-(Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cds)",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   Cdd u0 {1,D} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Cdd u0 {2,D} {8,D}
6   C   u0 {3,D}
7   Cd  u0 {4,D}
8   C   u0 {5,D}
""",
    thermo = 'Cds-(Cdd-Cd)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 551,
    label = "Cds-(Cdd-Cd)(Cds-Cdd)(Cds-Cdd)",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   C   u0 {2,D}
6   Cdd u0 {3,D}
7   Cdd u0 {4,D}
""",
    thermo = 'Cds-(Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 552,
    label = "Cds-(Cdd-Cd)(Cds-Cdd-O2d)(Cds-Cdd-O2d)",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cdd u0 {1,D} {7,D}
5   Cdd u0 {2,D} {8,D}
6   Cdd u0 {3,D} {9,D}
7   C   u0 {4,D}
8   O2d u0 {5,D}
9   O2d u0 {6,D}
""",
    thermo = 'Cds-Cds(Cds-Cdd-O2d)(Cds-Cdd-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 553,
    label = "Cds-(Cdd-Cd)(Cds-Cdd-O2d)(Cds-Cdd-Cd)",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cdd u0 {1,D} {7,D}
5   Cdd u0 {2,D} {8,D}
6   Cdd u0 {3,D} {9,D}
7   C   u0 {4,D}
8   O2d u0 {5,D}
9   C   u0 {6,D}
""",
    thermo = 'Cds-(Cdd-Cd)(Cds-Cdd-O2d)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 554,
    label = "Cds-(Cdd-Cd)(Cds-Cdd-S2d)(Cds-Cdd-S2d)",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cdd u0 {1,D} {7,D}
5   Cdd u0 {2,D} {8,D}
6   Cdd u0 {3,D} {9,D}
7   C   u0 {4,D}
8   S2d u0 {5,D}
9   S2d u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 555,
    label = "Cds-(Cdd-Cd)(Cds-Cdd-S2d)(Cds-Cdd-Cd)",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cdd u0 {1,D} {7,D}
5   Cdd u0 {2,D} {8,D}
6   Cdd u0 {3,D} {9,D}
7   C   u0 {4,D}
8   S2d u0 {5,D}
9   C   u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 556,
    label = "Cds-(Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cdd u0 {1,D} {7,D}
5   Cdd u0 {2,D} {8,D}
6   Cdd u0 {3,D} {9,D}
7   C   u0 {4,D}
8   C   u0 {5,D}
9   C   u0 {6,D}
""",
    thermo = 'Cds-(Cdd-Cd)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 557,
    label = "Cds-CddCtCs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D}
3   Ct  u0 {1,S}
4   Cs  u0 {1,S}
""",
    thermo = 'Cds-(Cdd-Cd)CtCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 558,
    label = "Cds-(Cdd-O2d)CtCs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Ct  u0 {1,S}
4   Cs  u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 559,
    label = "Cds-(Cdd-S2d)CtCs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Ct  u0 {1,S}
4   Cs  u0 {1,S}
5   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 560,
    label = "Cds-(Cdd-Cd)CtCs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Ct  u0 {1,S}
4   Cs  u0 {1,S}
5   C   u0 {2,D}
""",
    thermo = 'Cds-CdsCtCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 561,
    label = "Cds-CddCtCds",
    group = 
"""
1 * Cd      u0 {2,D} {3,S} {4,S}
2   Cdd     u0 {1,D}
3   Ct      u0 {1,S}
4   [Cd,CO] u0 {1,S}
""",
    thermo = 'Cds-(Cdd-Cd)(Cds-Cds)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 562,
    label = "Cds-(Cdd-O2d)(Cds-O2d)Ct",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   CO  u0 {1,S} {6,D}
4   Ct  u0 {1,S}
5   O2d u0 {2,D}
6   O2d u0 {3,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cds)(Cds-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 563,
    label = "Cds-(Cdd-O2d)(Cds-Cd)Ct",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Ct  u0 {1,S}
5   O2d u0 {2,D}
6   C   u0 {3,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cds)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 564,
    label = "Cds-(Cdd-O2d)(Cds-Cds)Ct",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Ct  u0 {1,S}
5   O2d u0 {2,D}
6   Cd  u0 {3,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 565,
    label = "Cds-(Cdd-O2d)(Cds-Cdd)Ct",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Ct  u0 {1,S}
5   O2d u0 {2,D}
6   Cdd u0 {3,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cdd-Cd)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 566,
    label = "Cds-(Cdd-O2d)(Cds-Cdd-O2d)Ct",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   Ct  u0 {1,S}
6   O2d u0 {3,D}
7   O2d u0 {4,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cdd-O2d)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 567,
    label = "Cds-(Cdd-O2d)(Cds-Cdd-Cd)Ct",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   Ct  u0 {1,S}
6   O2d u0 {3,D}
7   C   u0 {4,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cds)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 568,
    label = "Cds-(Cdd-S2d)(Cds-Cd)Ct",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Ct  u0 {1,S}
5   S2d u0 {2,D}
6   C   u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 569,
    label = "Cds-(Cdd-S2d)(Cds-Cds)Ct",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Ct  u0 {1,S}
5   S2d u0 {2,D}
6   Cd  u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 570,
    label = "Cds-(Cdd-S2d)(Cds-Cdd)Ct",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Ct  u0 {1,S}
5   S2d u0 {2,D}
6   Cdd u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 571,
    label = "Cds-(Cdd-S2d)(Cds-Cdd-S2d)Ct",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   Ct  u0 {1,S}
6   S2d u0 {3,D}
7   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 572,
    label = "Cds-(Cdd-S2d)(Cds-Cdd-Cd)Ct",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   Ct  u0 {1,S}
6   S2d u0 {3,D}
7   C   u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 573,
    label = "Cds-(Cdd-Cd)(Cds-Cd)Ct",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Ct  u0 {1,S}
5   C   u0 {2,D}
6   C   u0 {3,D}
""",
    thermo = 'Cds-(Cdd-Cd)(Cds-Cds)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 574,
    label = "Cds-(Cdd-Cd)(Cds-Cds)Ct",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Ct  u0 {1,S}
5   C   u0 {2,D}
6   Cd  u0 {3,D}
""",
    thermo = 'Cds-Cds(Cds-Cds)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 575,
    label = "Cds-(Cdd-Cd)(Cds-Cdd)Ct",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Ct  u0 {1,S}
5   C   u0 {2,D}
6   Cdd u0 {3,D}
""",
    thermo = 'Cds-(Cdd-Cd)(Cds-Cdd-Cd)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 576,
    label = "Cds-(Cdd-Cd)(Cds-Cdd-O2d)Ct",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   Ct  u0 {1,S}
6   C   u0 {3,D}
7   O2d u0 {4,D}
""",
    thermo = 'Cds-Cds(Cds-Cdd-O2d)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 577,
    label = "Cds-(Cdd-Cd)(Cds-Cdd-S2d)Ct",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   Ct  u0 {1,S}
6   C   u0 {3,D}
7   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 578,
    label = "Cds-(Cdd-Cd)(Cds-Cdd-Cd)Ct",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   Ct  u0 {1,S}
6   C   u0 {3,D}
7   C   u0 {4,D}
""",
    thermo = 'Cds-(Cdd-Cd)(Cds-Cds)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 579,
    label = "Cds-CddCtCt",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D}
3   Ct  u0 {1,S}
4   Ct  u0 {1,S}
""",
    thermo = 'Cds-(Cdd-Cd)CtCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 580,
    label = "Cds-(Cdd-O2d)CtCt",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Ct  u0 {1,S}
4   Ct  u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 581,
    label = "Cds-(Cdd-S2d)CtCt",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Ct  u0 {1,S}
4   Ct  u0 {1,S}
5   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 582,
    label = "Cds-(Cdd-Cd)CtCt",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Ct  u0 {1,S}
4   Ct  u0 {1,S}
5   C   u0 {2,D}
""",
    thermo = 'Cds-CdsCtCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 583,
    label = "Cds-CddCbCs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D}
3   Cb  u0 {1,S}
4   Cs  u0 {1,S}
""",
    thermo = 'Cds-(Cdd-Cd)CbCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 584,
    label = "Cds-(Cdd-O2d)CbCs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cb  u0 {1,S}
4   Cs  u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 585,
    label = "Cds-(Cdd-S2d)CbCs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cb  u0 {1,S}
4   Cs  u0 {1,S}
5   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 586,
    label = "Cds-(Cdd-Cd)CbCs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cb  u0 {1,S}
4   Cs  u0 {1,S}
5   C   u0 {2,D}
""",
    thermo = 'Cds-CdsCbCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 587,
    label = "Cds-CddCbCds",
    group = 
"""
1 * Cd      u0 {2,D} {3,S} {4,S}
2   Cdd     u0 {1,D}
3   Cb      u0 {1,S}
4   [Cd,CO] u0 {1,S}
""",
    thermo = 'Cds-(Cdd-Cd)(Cds-Cds)Cb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 588,
    label = "Cds-(Cdd-O2d)(Cds-O2d)Cb",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   CO  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   O2d u0 {2,D}
6   O2d u0 {3,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cds)(Cds-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 589,
    label = "Cds-(Cdd-O2d)(Cds-Cd)Cb",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   O2d u0 {2,D}
6   C   u0 {3,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cds)Cb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 590,
    label = "Cds-(Cdd-O2d)(Cds-Cds)Cb",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   O2d u0 {2,D}
6   Cd  u0 {3,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 591,
    label = "Cds-(Cdd-O2d)(Cds-Cdd)Cb",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   O2d u0 {2,D}
6   Cdd u0 {3,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cdd-Cd)Cb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 592,
    label = "Cds-(Cdd-O2d)(Cds-Cdd-O2d)Cb",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   Cb  u0 {1,S}
6   O2d u0 {3,D}
7   O2d u0 {4,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cdd-O2d)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 593,
    label = "Cds-(Cdd-O2d)(Cds-Cdd-Cd)Cb",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   Cb  u0 {1,S}
6   O2d u0 {3,D}
7   C   u0 {4,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cds)Cb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 594,
    label = "Cds-(Cdd-S2d)(Cds-Cd)Cb",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   S2d u0 {2,D}
6   C   u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 595,
    label = "Cds-(Cdd-S2d)(Cds-Cds)Cb",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   S2d u0 {2,D}
6   Cd  u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 596,
    label = "Cds-(Cdd-S2d)(Cds-Cdd)Cb",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   S2d u0 {2,D}
6   Cdd u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 597,
    label = "Cds-(Cdd-S2d)(Cds-Cdd-S2d)Cb",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   Cb  u0 {1,S}
6   S2d u0 {3,D}
7   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 598,
    label = "Cds-(Cdd-S2d)(Cds-Cdd-Cd)Cb",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   Cb  u0 {1,S}
6   S2d u0 {3,D}
7   C   u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 599,
    label = "Cds-(Cdd-Cd)(Cds-Cd)Cb",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   C   u0 {2,D}
6   C   u0 {3,D}
""",
    thermo = 'Cds-(Cdd-Cd)(Cds-Cds)Cb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 600,
    label = "Cds-(Cdd-Cd)(Cds-Cds)Cb",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   C   u0 {2,D}
6   Cd  u0 {3,D}
""",
    thermo = 'Cds-Cds(Cds-Cds)Cb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 601,
    label = "Cds-(Cdd-Cd)(Cds-Cdd)Cb",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   C   u0 {2,D}
6   Cdd u0 {3,D}
""",
    thermo = 'Cds-(Cdd-Cd)(Cds-Cdd-Cd)Cb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 602,
    label = "Cds-(Cdd-Cd)(Cds-Cdd-O2d)Cb",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   Cb  u0 {1,S}
6   C   u0 {3,D}
7   O2d u0 {4,D}
""",
    thermo = 'Cds-Cds(Cds-Cdd-O2d)Cb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 603,
    label = "Cds-(Cdd-Cd)(Cds-Cdd-S2d)Cb",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   Cb  u0 {1,S}
6   C   u0 {3,D}
7   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 604,
    label = "Cds-(Cdd-Cd)(Cds-Cdd-Cd)Cb",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {5,S}
2   Cd  u0 {1,S} {4,D}
3   Cdd u0 {1,D} {6,D}
4   Cdd u0 {2,D} {7,D}
5   Cb  u0 {1,S}
6   C   u0 {3,D}
7   C   u0 {4,D}
""",
    thermo = 'Cds-(Cdd-Cd)(Cds-Cds)Cb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 605,
    label = "Cds-CddCbCt",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D}
3   Cb  u0 {1,S}
4   Ct  u0 {1,S}
""",
    thermo = 'Cds-(Cdd-Cd)CbCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 606,
    label = "Cds-(Cdd-O2d)CbCt",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cb  u0 {1,S}
4   Ct  u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cds)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 607,
    label = "Cds-(Cdd-S2d)CbCt",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cb  u0 {1,S}
4   Ct  u0 {1,S}
5   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 608,
    label = "Cds-(Cdd-Cd)CbCt",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cb  u0 {1,S}
4   Ct  u0 {1,S}
5   C   u0 {2,D}
""",
    thermo = 'Cds-CdsCbCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 609,
    label = "Cds-CddCbCb",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D}
3   Cb  u0 {1,S}
4   Cb  u0 {1,S}
""",
    thermo = 'Cds-(Cdd-Cd)CbCb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 610,
    label = "Cds-(Cdd-O2d)CbCb",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cb  u0 {1,S}
4   Cb  u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = 'Cds-(Cdd-O2d)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 611,
    label = "Cds-(Cdd-S2d)CbCb",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cb  u0 {1,S}
4   Cb  u0 {1,S}
5   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 612,
    label = "Cds-(Cdd-Cd)CbCb",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cb  u0 {1,S}
4   Cb  u0 {1,S}
5   C   u0 {2,D}
""",
    thermo = 'Cds-CdsCbCb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 613,
    label = "Cds-CdsC=SC=S",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   CS  u0 {1,S} {5,D}
3   CS  u0 {1,S} {6,D}
4   Cd  u0 {1,D}
5   S2d u0 {2,D}
6   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 614,
    label = "Cds-(Cdd-Cd)C=S(Cds-Cd)",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   CS  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {6,D}
5   C   u0 {2,D}
6   C   u0 {4,D}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 615,
    label = "Cds-(Cdd-Cd)C=S(Cds-Cds)",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   CS  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {6,D}
5   C   u0 {2,D}
6   Cd  u0 {4,D}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 616,
    label = "Cds-(Cdd-Cd)C=S(Cds-Cdd)",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   CS  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {6,D}
5   C   u0 {2,D}
6   Cdd u0 {4,D}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 617,
    label = "Cds-(Cdd-Cd)C=S(Cds-Cdd-Cd)",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   Cdd u0 {1,D} {6,D}
4   CS  u0 {1,S} {7,D}
5   Cdd u0 {2,D} {8,D}
6   C   u0 {3,D}
7   S2d u0 {4,D}
8   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 618,
    label = "Cds-(Cdd-Cd)C=S(Cds-Cdd-S2d)",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   Cdd u0 {1,D} {6,D}
4   CS  u0 {1,S} {7,D}
5   Cdd u0 {2,D} {8,D}
6   C   u0 {3,D}
7   S2d u0 {4,D}
8   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 619,
    label = "Cds-(Cdd-S2d)C=SCs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   CS  u0 {1,S} {6,D}
4   Cs  u0 {1,S}
5   S2d u0 {2,D}
6   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 620,
    label = "Cds-(Cdd-S2d)C=SCt",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   CS  u0 {1,S} {6,D}
4   Ct  u0 {1,S}
5   S2d u0 {2,D}
6   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 621,
    label = "Cds-(Cdd-S2d)C=SCb",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   CS  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   S2d u0 {2,D}
6   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 622,
    label = "Cds-(Cdd-Cd)C=SC=S",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   CS  u0 {1,S} {6,D}
4   CS  u0 {1,S} {7,D}
5   C   u0 {2,D}
6   S2d u0 {3,D}
7   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 623,
    label = "Cds-(Cdd-S2d)(Cds-Cd)C=S",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   CS  u0 {1,S} {7,D}
5   S2d u0 {2,D}
6   C   u0 {3,D}
7   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 624,
    label = "Cds-(Cdd-S2d)(Cds-Cds)C=S",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   CS  u0 {1,S} {7,D}
5   S2d u0 {2,D}
6   Cd  u0 {3,D}
7   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 625,
    label = "Cds-(Cdd-S2d)(Cds-Cdd)C=S",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   Cd  u0 {1,S} {6,D}
4   CS  u0 {1,S} {7,D}
5   S2d u0 {2,D}
6   Cdd u0 {3,D}
7   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 626,
    label = "Cds-(Cdd-S2d)(Cds-Cdd-S2d)C=S",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   Cdd u0 {1,D} {6,D}
4   CS  u0 {1,S} {7,D}
5   Cdd u0 {2,D} {8,D}
6   S2d u0 {3,D}
7   S2d u0 {4,D}
8   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 627,
    label = "Cds-(Cdd-S2d)(Cds-Cdd-Cd)C=S",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   Cdd u0 {1,D} {6,D}
4   CS  u0 {1,S} {7,D}
5   Cdd u0 {2,D} {8,D}
6   S2d u0 {3,D}
7   S2d u0 {4,D}
8   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 628,
    label = "Cds-CdsCbC=S",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {4,S}
2   CS  u0 {1,S} {5,D}
3   Cd  u0 {1,D}
4   Cb  u0 {1,S}
5   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 629,
    label = "Cds-CdsCtC=S",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {4,S}
2   CS  u0 {1,S} {5,D}
3   Cd  u0 {1,D}
4   Ct  u0 {1,S}
5   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 630,
    label = "Cds-CdsC=SCs",
    group = 
"""
1 * Cd  u0 {2,S} {3,D} {4,S}
2   CS  u0 {1,S} {5,D}
3   Cd  u0 {1,D}
4   Cs  u0 {1,S}
5   S2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.21,4.51,4.77,4.99,5.4,5.66,5.98],'cal/(mol*K)'),
        H298 = (9.24,'kcal/mol'),
        S298 = (-11.25,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 631,
    label = "Cds-CdsC=S(Cds-Cd)",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   CS  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {5,D}
4   Cd  u0 {1,D}
5   C   u0 {3,D}
6   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 632,
    label = "Cds-CdsC=S(Cds-Cds)",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   CS  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {5,D}
4   Cd  u0 {1,D}
5   Cd  u0 {3,D}
6   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 633,
    label = "Cds-CdsC=S(Cds-Cdd)",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   CS  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {5,D}
4   Cd  u0 {1,D}
5   Cdd u0 {3,D}
6   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 634,
    label = "Cds-CdsC=S(Cds-Cdd-Cd)",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {5,D}
2   Cd  u0 {1,S} {4,D}
3   CS  u0 {1,S} {6,D}
4   Cdd u0 {2,D} {7,D}
5   Cd  u0 {1,D}
6   S2d u0 {3,D}
7   C   u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 635,
    label = "Cds-CdsC=S(Cds-Cdd-S2d)",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {5,D}
2   Cd  u0 {1,S} {4,D}
3   CS  u0 {1,S} {6,D}
4   Cdd u0 {2,D} {7,D}
5   Cd  u0 {1,D}
6   S2d u0 {3,D}
7   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 636,
    label = "Cds-(Cdd-S2d)C=SC=S",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   CS  u0 {1,S} {6,D}
4   CS  u0 {1,S} {7,D}
5   S2d u0 {2,D}
6   S2d u0 {3,D}
7   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 637,
    label = "C=S-SsSs",
    group = 
"""
1 * CS  u0 {2,D} {3,S} {4,S}
2   S2d u0 {1,D}
3   S2s u0 {1,S}
4   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 638,
    label = "C=S-CH",
    group = 
"""
1 * CS u0 {2,D} {3,S} {4,S}
2   S  u0 {1,D}
3   C  u0 {1,S}
4   H  u0 {1,S}
""",
    thermo = 'C=S-CsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 639,
    label = "C=S-CsH",
    group = 
"""
1 * CS u0 {2,D} {3,S} {4,S}
2   S  u0 {1,D}
3   Cs u0 {1,S}
4   H  u0 {1,S}
""",
    thermo = 'C=S2-CsH',
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2018""",
    longDesc = 
"""
"
""",
)

entry(
    index = 640,
    label = "C=S2-CsH",
    group = 
"""
1 * CS  u0 {2,D} {3,S} {4,S}
2   S2d u0 {1,D}
3   Cs  u0 {1,S}
4   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([10.09,9.4,9.55,10.13,11.37,12.36,13.87],'cal/(mol*K)'),
        H298 = (24.14,'kcal/mol'),
        S298 = (36.84,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 641,
    label = "C=S4-CsH",
    group = 
"""
1 * CS         u0 {2,D} {3,S} {4,S}
2   [S4d,S4dd] u0 {1,D}
3   Cs         u0 {1,S}
4   H          u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.55,7.96,8.45,9.22,10.82,11.81,12.94],'cal/(mol*K)'),
        H298 = (16.51,'kcal/mol'),
        S298 = (38.28,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 642,
    label = "C=S-CdsH",
    group = 
"""
1 * CS      u0 {2,D} {3,S} {4,S}
2   S2d     u0 {1,D}
3   [Cd,Cb] u0 {1,S}
4   H       u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.41,9.18,10.55,11.53,12.8,13.57,14.62],'cal/(mol*K)'),
        H298 = (24.85,'kcal/mol'),
        S298 = (33.97,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 643,
    label = "C=S-(Cds-Cd)H",
    group = 
"""
1 * CS          u0 {2,S} {3,D} {4,S}
2   Cd          u0 {1,S} {5,D}
3   S2d         u0 {1,D}
4   H           u0 {1,S}
5   [Cd,Cdd,CO] u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 644,
    label = "C=S-(Cds-Cdd)H",
    group = 
"""
1 * CS  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   S2d u0 {1,D}
4   H   u0 {1,S}
5   Cdd u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 645,
    label = "C=S-(Cds-Cdd-Cd)H",
    group = 
"""
1 * CS  u0 {2,S} {4,D} {5,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {6,D}
4   S2d u0 {1,D}
5   H   u0 {1,S}
6   C   u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 646,
    label = "C=S-(Cds-Cdd-S2d)H",
    group = 
"""
1 * CS  u0 {2,S} {4,D} {5,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {6,D}
4   S2d u0 {1,D}
5   H   u0 {1,S}
6   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 647,
    label = "C=S-(Cds-Cds)H",
    group = 
"""
1 * CS  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   S2d u0 {1,D}
4   H   u0 {1,S}
5   Cd  u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 648,
    label = "C=S-CtH",
    group = 
"""
1 * CS  u0 {2,D} {3,S} {4,S}
2   S2d u0 {1,D}
3   Ct  u0 {1,S}
4   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.07,8.97,10.47,10.17,11.38,12.33,13.6],'cal/(mol*K)'),
        H298 = (30.53,'kcal/mol'),
        S298 = (36.94,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 649,
    label = "C=S-C=SH",
    group = 
"""
1 * CS  u0 {2,S} {3,D} {4,S}
2   CS  u0 {1,S} {5,D}
3   S2d u0 {1,D}
4   H   u0 {1,S}
5   S2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.79,9.18,10.41,11.42,12.82,13.64,14.54],'cal/(mol*K)'),
        H298 = (26.96,'kcal/mol'),
        S298 = (35.65,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 650,
    label = "C=S-CC",
    group = 
"""
1 * CS  u0 {2,D} {3,S} {4,S}
2   S2d u0 {1,D}
3   C   u0 {1,S}
4   C   u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 651,
    label = "C=S-CbCds",
    group = 
"""
1 * CS  u0 {2,D} {3,S} {4,S}
2   S2d u0 {1,D}
3   Cb  u0 {1,S}
4   Cd  u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 652,
    label = "C=S-Cb(Cds-Cd)",
    group = 
"""
1 * CS          u0 {2,S} {3,D} {4,S}
2   Cd          u0 {1,S} {5,D}
3   S2d         u0 {1,D}
4   Cb          u0 {1,S}
5   [Cd,Cdd,CO] u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 653,
    label = "C=S-Cb(Cds-Cds)",
    group = 
"""
1 * CS  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   S2d u0 {1,D}
4   Cb  u0 {1,S}
5   Cd  u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 654,
    label = "C=S-Cb(Cds-Cdd)",
    group = 
"""
1 * CS  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   S2d u0 {1,D}
4   Cb  u0 {1,S}
5   Cdd u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 655,
    label = "C=S-Cb(Cds-Cdd-S2d)",
    group = 
"""
1 * CS  u0 {2,S} {4,D} {5,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {6,D}
4   S2d u0 {1,D}
5   Cb  u0 {1,S}
6   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 656,
    label = "C=S-Cb(Cds-Cdd-Cd)",
    group = 
"""
1 * CS  u0 {2,S} {4,D} {5,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {6,D}
4   S2d u0 {1,D}
5   Cb  u0 {1,S}
6   C   u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 657,
    label = "C=S-CtCt",
    group = 
"""
1 * CS  u0 {2,D} {3,S} {4,S}
2   S2d u0 {1,D}
3   Ct  u0 {1,S}
4   Ct  u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 658,
    label = "C=S-CbCb",
    group = 
"""
1 * CS  u0 {2,D} {3,S} {4,S}
2   S2d u0 {1,D}
3   Cb  u0 {1,S}
4   Cb  u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 659,
    label = "C=S-CdsCds",
    group = 
"""
1 * CS  u0 {2,D} {3,S} {4,S}
2   S2d u0 {1,D}
3   Cd  u0 {1,S}
4   Cd  u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 660,
    label = "C=S-(Cds-Cd)(Cds-Cd)",
    group = 
"""
1 * CS          u0 {2,S} {3,S} {4,D}
2   Cd          u0 {1,S} {5,D}
3   Cd          u0 {1,S} {6,D}
4   S2d         u0 {1,D}
5   [Cd,Cdd,CO] u0 {2,D}
6   [Cd,Cdd,CO] u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 661,
    label = "C=S-(Cds-Cdd)(Cds-Cds)",
    group = 
"""
1 * CS  u0 {2,S} {3,S} {4,D}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,S} {6,D}
4   S2d u0 {1,D}
5   Cdd u0 {2,D}
6   Cd  u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 662,
    label = "C=S-(Cds-Cdd-Cd)(Cds-Cds)",
    group = 
"""
1 * CS  u0 {2,S} {3,S} {5,D}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {6,D}
4   Cdd u0 {2,D} {7,D}
5   S2d u0 {1,D}
6   Cd  u0 {3,D}
7   C   u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 663,
    label = "C=S-(Cds-Cdd-S2d)(Cds-Cds)",
    group = 
"""
1 * CS  u0 {2,S} {3,S} {5,D}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {6,D}
4   Cdd u0 {2,D} {7,D}
5   S2d u0 {1,D}
6   Cd  u0 {3,D}
7   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 664,
    label = "C=S-(Cds-Cds)(Cds-Cds)",
    group = 
"""
1 * CS  u0 {2,S} {3,S} {4,D}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,S} {6,D}
4   S2d u0 {1,D}
5   Cd  u0 {2,D}
6   Cd  u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 665,
    label = "C=S-(Cds-Cdd)(Cds-Cdd)",
    group = 
"""
1 * CS  u0 {2,S} {3,S} {4,D}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,S} {6,D}
4   S2d u0 {1,D}
5   Cdd u0 {2,D}
6   Cdd u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 666,
    label = "C=S-(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group = 
"""
1 * CS  u0 {2,S} {3,S} {6,D}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {7,D}
5   Cdd u0 {3,D} {8,D}
6   S2d u0 {1,D}
7   C   u0 {4,D}
8   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 667,
    label = "C=S-(Cds-Cdd-S2d)(Cds-Cdd-S2d)",
    group = 
"""
1 * CS  u0 {2,S} {3,S} {6,D}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {7,D}
5   Cdd u0 {3,D} {8,D}
6   S2d u0 {1,D}
7   S2d u0 {4,D}
8   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 668,
    label = "C=S-(Cds-Cdd-Cd)(Cds-Cdd-S2d)",
    group = 
"""
1 * CS  u0 {2,S} {3,S} {6,D}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {7,D}
5   Cdd u0 {3,D} {8,D}
6   S2d u0 {1,D}
7   C   u0 {4,D}
8   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 669,
    label = "C=S-CtCds",
    group = 
"""
1 * CS  u0 {2,D} {3,S} {4,S}
2   S2d u0 {1,D}
3   Ct  u0 {1,S}
4   Cd  u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 670,
    label = "C=S-Ct(Cds-Cd)",
    group = 
"""
1 * CS          u0 {2,S} {3,D} {4,S}
2   Cd          u0 {1,S} {5,D}
3   S2d         u0 {1,D}
4   Ct          u0 {1,S}
5   [Cd,Cdd,CO] u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 671,
    label = "C=S-Ct(Cds-Cds)",
    group = 
"""
1 * CS  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   S2d u0 {1,D}
4   Ct  u0 {1,S}
5   Cd  u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 672,
    label = "C=S-Ct(Cds-Cdd)",
    group = 
"""
1 * CS  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   S2d u0 {1,D}
4   Ct  u0 {1,S}
5   Cdd u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 673,
    label = "C=S-Ct(Cds-Cdd-Cd)",
    group = 
"""
1 * CS  u0 {2,S} {4,D} {5,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {6,D}
4   S2d u0 {1,D}
5   Ct  u0 {1,S}
6   C   u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 674,
    label = "C=S-Ct(Cds-Cdd-S2d)",
    group = 
"""
1 * CS  u0 {2,S} {4,D} {5,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {6,D}
4   S2d u0 {1,D}
5   Ct  u0 {1,S}
6   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 675,
    label = "C=S-CbCt",
    group = 
"""
1 * CS  u0 {2,D} {3,S} {4,S}
2   S2d u0 {1,D}
3   Cb  u0 {1,S}
4   Ct  u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 676,
    label = "C=S-CsCs",
    group = 
"""
1 * CS  u0 {2,D} {3,S} {4,S}
2   S2d u0 {1,D}
3   Cs  u0 {1,S}
4   Cs  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([10.64,8.16,7.45,7.69,8.72,9.5,10.5],'cal/(mol*K)'),
        H298 = (20.9,'kcal/mol'),
        S298 = (16.55,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 677,
    label = "C=S-CdsCs",
    group = 
"""
1 * CS  u0 {2,D} {3,S} {4,S}
2   S2d u0 {1,D}
3   Cd  u0 {1,S}
4   Cs  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.58,9.38,9.53,9.87,10.47,10.79,11.17],'cal/(mol*K)'),
        H298 = (23.84,'kcal/mol'),
        S298 = (12.34,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""

""",
)

entry(
    index = 678,
    label = "C=S-(Cds-Cd)Cs",
    group = 
"""
1 * CS          u0 {2,S} {3,D} {4,S}
2   Cd          u0 {1,S} {5,D}
3   S2d         u0 {1,D}
4   Cs          u0 {1,S}
5   [Cd,Cdd,CO] u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 679,
    label = "C=S-(Cds-Cds)Cs",
    group = 
"""
1 * CS  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   S2d u0 {1,D}
4   Cs  u0 {1,S}
5   Cd  u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 680,
    label = "C=S-(Cds-Cdd)Cs",
    group = 
"""
1 * CS  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   S2d u0 {1,D}
4   Cs  u0 {1,S}
5   Cdd u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 681,
    label = "C=S-(Cds-Cdd-S2d)Cs",
    group = 
"""
1 * CS  u0 {2,S} {4,D} {5,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {6,D}
4   S2d u0 {1,D}
5   Cs  u0 {1,S}
6   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 682,
    label = "C=S-(Cds-Cdd-Cd)Cs",
    group = 
"""
1 * CS  u0 {2,S} {4,D} {5,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {6,D}
4   S2d u0 {1,D}
5   Cs  u0 {1,S}
6   C   u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 683,
    label = "C=S-CtCs",
    group = 
"""
1 * CS  u0 {2,D} {3,S} {4,S}
2   S2d u0 {1,D}
3   Ct  u0 {1,S}
4   Cs  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.44,8.3,8.72,7.99,8.85,9.55,10.38],'cal/(mol*K)'),
        H298 = (26.63,'kcal/mol'),
        S298 = (16.54,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""

""",
)

entry(
    index = 684,
    label = "C=S-CbCs",
    group = 
"""
1 * CS  u0 {2,D} {3,S} {4,S}
2   S2d u0 {1,D}
3   Cb  u0 {1,S}
4   Cs  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.94,8.73,9,9.46,10.17,10.53,10.9],'cal/(mol*K)'),
        H298 = (23.58,'kcal/mol'),
        S298 = (13.65,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 685,
    label = "C=S-C=SCs",
    group = 
"""
1 * CS  u0 {2,S} {3,D} {4,S}
2   CS  u0 {1,S} {5,D}
3   S2d u0 {1,D}
4   Cs  u0 {1,S}
5   S2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.9,8.29,8.4,8.85,9.71,10.22,10.71],'cal/(mol*K)'),
        H298 = (24.34,'kcal/mol'),
        S298 = (15.85,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 686,
    label = "C=S-CtC=S",
    group = 
"""
1 * CS  u0 {2,S} {3,D} {4,S}
2   CS  u0 {1,S} {5,D}
3   S2d u0 {1,D}
4   Ct  u0 {1,S}
5   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 687,
    label = "C=S-(Cds-Cd)C=S",
    group = 
"""
1 * CS  u0 {2,S} {3,D} {4,S}
2   CS  u0 {1,S} {5,D}
3   S2d u0 {1,D}
4   Cd  u0 {1,S}
5   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 688,
    label = "C=S-(Cds-Cdd)C=S",
    group = 
"""
1 * CS  u0 {2,S} {3,S} {4,D}
2   Cd  u0 {1,S} {5,D}
3   CS  u0 {1,S} {6,D}
4   S2d u0 {1,D}
5   Cdd u0 {2,D}
6   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 689,
    label = "C=S-(Cds-Cdd-Cd)C=S",
    group = 
"""
1 * CS  u0 {2,S} {3,S} {5,D}
2   Cd  u0 {1,S} {4,D}
3   CS  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {6,D}
5   S2d u0 {1,D}
6   C   u0 {4,D}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 690,
    label = "C=S-(Cds-Cdd-S2d)C=S",
    group = 
"""
1 * CS  u0 {2,S} {3,S} {5,D}
2   Cd  u0 {1,S} {4,D}
3   CS  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {6,D}
5   S2d u0 {1,D}
6   S2d u0 {4,D}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 691,
    label = "C=S-(Cds-Cds)C=S",
    group = 
"""
1 * CS  u0 {2,S} {3,S} {4,D}
2   Cd  u0 {1,S} {5,D}
3   CS  u0 {1,S} {6,D}
4   S2d u0 {1,D}
5   Cd  u0 {2,D}
6   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 692,
    label = "C=S-C=SC=S",
    group = 
"""
1 * CS  u0 {2,S} {3,S} {4,D}
2   CS  u0 {1,S} {5,D}
3   CS  u0 {1,S} {6,D}
4   S2d u0 {1,D}
5   S2d u0 {2,D}
6   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 693,
    label = "C=S-CbC=S",
    group = 
"""
1 * CS  u0 {2,S} {3,D} {4,S}
2   CS  u0 {1,S} {5,D}
3   S2d u0 {1,D}
4   Cb  u0 {1,S}
5   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 694,
    label = "C=S-HH",
    group = 
"""
1 * CS u0 {2,D} {3,S} {4,S}
2   S  u0 {1,D}
3   H  u0 {1,S}
4   H  u0 {1,S}
""",
    thermo = 'C=S2d-HH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 695,
    label = "C=S2d-HH",
    group = 
"""
1 * CS  u0 {2,D} {3,S} {4,S}
2   S2d u0 {1,D}
3   H   u0 {1,S}
4   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.08,10.35,11.52,12.5,14.08,15.25,17.14],'cal/(mol*K)'),
        H298 = (27.7,'kcal/mol'),
        S298 = (56.5,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 696,
    label = "C=S4d-HH",
    group = 
"""
1 * CS         u0 {2,D} {3,S} {4,S}
2   [S4d,S4dd] u0 {1,D}
3   H          u0 {1,S}
4   H          u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.82,9.96,10.83,11.76,13.57,14.92,16.57],'cal/(mol*K)'),
        H298 = (16.47,'kcal/mol'),
        S298 = (57.73,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 697,
    label = "C=S6dd-HH",
    group = 
"""
1 * CS   u0 {2,D} {3,S} {4,S}
2   S6dd u0 {1,D}
3   H    u0 {1,S}
4   H    u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.79,11.27,12.37,13.22,14.74,15.89,17.4],'cal/(mol*K)'),
        H298 = (14.2,'kcal/mol'),
        S298 = (54.75,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 698,
    label = "C=S6ddd-HH",
    group = 
"""
1 * CS           u0 {2,D} {3,S} {4,S}
2   [S6ddd,S6td] u0 {1,D}
3   H            u0 {1,S}
4   H            u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([10.85,10.47,11.15,12.02,14.01,15.34,17.31],'cal/(mol*K)'),
        H298 = (10.97,'kcal/mol'),
        S298 = (52.97,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 699,
    label = "C=S-SH",
    group = 
"""
1 * CS u0 {2,D} {3,S} {4,S}
2   S  u0 {1,D}
3   S  u0 {1,S}
4   H  u0 {1,S}
""",
    thermo = 'C=S-S2H',
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2018""",
    longDesc = 
"""
"
""",
)

entry(
    index = 700,
    label = "C=S-S2H",
    group = 
"""
1 * CS  u0 {2,D} {3,S} {4,S}
2   S2d u0 {1,D}
3   S2s u0 {1,S}
4   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.85,10.7,11.13,11.7,14.22,15.76,15.35],'cal/(mol*K)'),
        H298 = (30.11,'kcal/mol'),
        S298 = (39.21,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 701,
    label = "C=S-S4H",
    group = 
"""
1 * CS                u0 {2,D} {3,S} {4,S}
2   S2d               u0 {1,D}
3   [S4s,S4d,S4b,S4t] u0 {1,S}
4   H                 u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.98,9.29,9.55,10.22,13.01,14.82,15.15],'cal/(mol*K)'),
        H298 = (45.35,'kcal/mol'),
        S298 = (42.42,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 702,
    label = "C=S-S6H",
    group = 
"""
1 * CS                      u0 {2,D} {3,S} {4,S}
2   S2d                     u0 {1,D}
3   [S6s,S6d,S6dd,S6t,S6td] u0 {1,S}
4   H                       u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.5,9.98,11.08,11.93,13.14,13.95,15.08],'cal/(mol*K)'),
        H298 = (21.69,'kcal/mol'),
        S298 = (34.23,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2018""",
    longDesc = 
"""
"
""",
)

entry(
    index = 703,
    label = "C=S6-S2H",
    group = 
"""
1 * CS                    u0 {2,D} {3,S} {4,S}
2   [S6d,S6dd,S6ddd,S6td] u0 {1,D}
3   S2s                   u0 {1,S}
4   H                     u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([12.67,15.3,16,16.78,20.51,22.83,21.04],'cal/(mol*K)'),
        H298 = (32.31,'kcal/mol'),
        S298 = (64.82,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 704,
    label = "C=S-CSs",
    group = 
"""
1 * CS u0 {2,D} {3,S} {4,S}
2   S  u0 {1,D}
3   C  u0 {1,S}
4   S  u0 {1,S}
""",
    thermo = 'C=S-CsSs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 705,
    label = "C=S-CbSs",
    group = 
"""
1 * CS  u0 {2,D} {3,S} {4,S}
2   S2d u0 {1,D}
3   Cb  u0 {1,S}
4   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 706,
    label = "C=S-CdsSs",
    group = 
"""
1 * CS  u0 {2,D} {3,S} {4,S}
2   S2d u0 {1,D}
3   Cd  u0 {1,S}
4   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 707,
    label = "C=S-(Cds-Cd)S2s",
    group = 
"""
1 * CS          u0 {2,S} {3,D} {4,S}
2   Cd          u0 {1,S} {5,D}
3   S2d         u0 {1,D}
4   S2s         u0 {1,S}
5   [Cd,Cdd,CO] u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 708,
    label = "C=S-(Cds-Cds)S2s",
    group = 
"""
1 * CS  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   S2d u0 {1,D}
4   S2s u0 {1,S}
5   Cd  u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 709,
    label = "C=S-(Cds-Cdd)S2s",
    group = 
"""
1 * CS  u0 {2,S} {3,D} {4,S}
2   Cd  u0 {1,S} {5,D}
3   S2d u0 {1,D}
4   S2s u0 {1,S}
5   Cdd u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 710,
    label = "C=S-(Cds-Cdd-Cd)S2s",
    group = 
"""
1 * CS  u0 {2,S} {4,D} {5,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {6,D}
4   S2d u0 {1,D}
5   S2s u0 {1,S}
6   C   u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 711,
    label = "C=S-(Cds-Cdd-S2d)S2s",
    group = 
"""
1 * CS  u0 {2,S} {4,D} {5,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {6,D}
4   S2d u0 {1,D}
5   S2s u0 {1,S}
6   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 712,
    label = "C=S-S(CO)",
    group = 
"""
1 * CS  u0 {2,D} {3,S} {4,S}
2   S2d u0 {1,D}
3   CO  u0 {1,S}
4   S   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.4,8.38,9.16,9.8,10.72,11.25,11.66],'cal/(mol*K)'),
        H298 = (21.35,'kcal/mol'),
        S298 = (14.52,'cal/(mol*K)'),
    ),
    shortDesc = """CBS-QB3 GA 1D-HR Aaron Vandeputte 2010""",
    longDesc = 
"""

""",
)

entry(
    index = 713,
    label = "C=S-CtSs",
    group = 
"""
1 * CS  u0 {2,D} {3,S} {4,S}
2   S2d u0 {1,D}
3   Ct  u0 {1,S}
4   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 714,
    label = "C=S-CsSs",
    group = 
"""
1 * CS  u0 {2,D} {3,S} {4,S}
2   S2d u0 {1,D}
3   Cs  u0 {1,S}
4   S2s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([10.76,9.54,8.97,9.22,11.6,13.02,12.13],'cal/(mol*K)'),
        H298 = (26.79,'kcal/mol'),
        S298 = (18.82,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""

""",
)

entry(
    index = 715,
    label = "C=S-C=SSs",
    group = 
"""
1 * CS  u0 {2,S} {3,D} {4,S}
2   CS  u0 {1,S} {5,D}
3   S2d u0 {1,D}
4   S2s u0 {1,S}
5   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 716,
    label = "Cds-CdIH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   C   u0 {1,D}
3   I1s u0 {1,S}
4   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.8,10,10.9,11.6,12.6,13.3,14],'cal/(mol*K)'),
        H298 = (24.5,'kcal/mol'),
        S298 = (40.5,'cal/(mol*K)'),
    ),
    shortDesc = """Cd-(I)(H) BENSON""",
    longDesc = 
"""
Thermochemical Kinetics 2nd Ed., by Sidney Benson (Table A4, p.280)
Cpdata at 1500K was not in the book, Cpdata at 1500K = Cpdata at 1000K + 0.7
""",
)

entry(
    index = 717,
    label = "C=S-OsH",
    group = 
"""
1 * CS  u0 {2,D} {3,S} {4,S}
2   S   u0 {1,D}
3   O2s u0 {1,S}
4   H   u0 {1,S}
""",
    thermo = 'C=S2-OsH',
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2018""",
    longDesc = 
"""
"
""",
)

entry(
    index = 718,
    label = "C=S2-OsH",
    group = 
"""
1 * CS  u0 {2,D} {3,S} {4,S}
2   S2d u0 {1,D}
3   O2s u0 {1,S}
4   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.1,10.53,11.43,12.21,13.48,14.27,15.01],'cal/(mol*K)'),
        H298 = (19.05,'kcal/mol'),
        S298 = (34.45,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 719,
    label = "C=S4-OsH",
    group = 
"""
1 * CS         u0 {2,D} {3,S} {4,S}
2   [S4d,S4dd] u0 {1,D}
3   O2s        u0 {1,S}
4   H          u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.67,8.18,9.52,10.9,13.34,14.97,16.45],'cal/(mol*K)'),
        H298 = (9.61,'kcal/mol'),
        S298 = (32.61,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 720,
    label = "C=S-CsOs",
    group = 
"""
1 * CS  u0 {2,D} {3,S} {4,S}
2   S2d u0 {1,D}
3   O2s u0 {1,S}
4   Cs  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([10.55,9.84,9.67,10.07,11.21,11.91,12.2],'cal/(mol*K)'),
        H298 = (11.72,'kcal/mol'),
        S298 = (12.22,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""

""",
)

entry(
    index = 721,
    label = "C=S-OsOs",
    group = 
"""
1 * CS  u0 {2,D} {3,S} {4,S}
2   S2d u0 {1,D}
3   O2s u0 {1,S}
4   O2s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([12.25,13.98,14.47,14.69,14.74,14.19,12.39],'cal/(mol*K)'),
        H298 = (9.69,'kcal/mol'),
        S298 = (11.28,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""

""",
)

entry(
    index = 722,
    label = "C=S-OsS",
    group = 
"""
1 * CS  u0 {2,D} {3,S} {4,S}
2   S2d u0 {1,D}
3   O2s u0 {1,S}
4   S   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.39,10.42,10.78,11.43,14.16,15.66,14.34],'cal/(mol*K)'),
        H298 = (36.94,'kcal/mol'),
        S298 = (16.85,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""

""",
)

entry(
    index = 723,
    label = "Cd-HHN",
    group = 
"""
1 * Cd u0 {2,S} {3,S} {4,D}
2   H  u0 {1,S}
3   H  u0 {1,S}
4   N  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.05807,7.45336,8.74103,9.99047,11.9119,13.3251,15.1939],'cal/(mol*K)','+|-',[0.625601,0.668132,0.689055,0.684162,0.665526,0.653584,0.6669]),
        H298 = (20.4503,'kcal/mol','+|-',2.5075),
        S298 = (21.9548,'cal/(mol*K)','+|-',1.92115),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library                 | Number of Species
CHON_G4                 |         5
thermo_DFT_CCSDTF12_BAC |         1
""",
)

entry(
    index = 724,
    label = "Cd-N3dHH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   N3d u0 {1,D}
3   H   u0 {1,S}
4   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.01377,8.53539,9.95906,11.4181,13.7456,15.5167,18.0328],'cal/(mol*K)','+|-',[0.46982,0.501761,0.517474,0.513799,0.499803,0.490835,0.500835]),
        H298 = (25.2069,'kcal/mol','+|-',1.88311),
        S298 = (33.1036,'cal/(mol*K)','+|-',1.44276),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHON_G4  |         95
BurcatNS |         1
""",
)

entry(
    index = 725,
    label = "CO-HNO",
    group = 
"""
1 * CO  u0 {2,S} {3,S} {4,D}
2   N   u0 {1,S}
3   H   u0 {1,S}
4   O2d u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.97298,9.90797,10.5369,10.9322,11.5881,12.1148,13.2149],'cal/(mol*K)','+|-',[0.41316,0.441248,0.455066,0.451834,0.439527,0.43164,0.440434]),
        H298 = (-31.0551,'kcal/mol','+|-',1.656),
        S298 = (33.9112,'cal/(mol*K)','+|-',1.26877),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         13
""",
)

entry(
    index = 726,
    label = "Cds-OdN3sH",
    group = 
"""
1 * CO  u0 {2,S} {3,S} {4,D}
2   N3s u0 {1,S}
3   H   u0 {1,S}
4   O2d u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.13543,9.0286,9.74035,10.2805,11.2193,11.9291,13.347],'cal/(mol*K)','+|-',[0.379499,0.405299,0.417991,0.415023,0.403718,0.396474,0.404551]),
        H298 = (-38.3062,'kcal/mol','+|-',1.52108),
        S298 = (32.7398,'cal/(mol*K)','+|-',1.1654),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         21
CHON    |         1
""",
)

entry(
    index = 727,
    label = "CO-CNO",
    group = 
"""
1 * CO  u0 {2,S} {3,S} {4,D}
2   C   u0 {1,S}
3   N   u0 {1,S}
4   O2d u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.38216,8.1686,8.64646,8.93052,9.40138,9.75912,9.89671],'cal/(mol*K)','+|-',[0.554163,0.591837,0.610371,0.606036,0.589529,0.57895,0.590746]),
        H298 = (-35.1496,'kcal/mol','+|-',2.22116),
        S298 = (14.9763,'cal/(mol*K)','+|-',1.70177),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         7
CHON    |         1
""",
)

entry(
    index = 728,
    label = "Cds-OdN3sCs",
    group = 
"""
1 * CO  u0 {2,S} {3,S} {4,D}
2   N3s u0 {1,S}
3   Cs  u0 {1,S}
4   O2d u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.33194,7.33141,8.10367,8.61099,9.34269,9.83821,10.8027],'cal/(mol*K)','+|-',[0.5554,0.593158,0.611733,0.607389,0.590844,0.580243,0.592064]),
        H298 = (-38.3664,'kcal/mol','+|-',2.22612),
        S298 = (12.3176,'cal/(mol*K)','+|-',1.70557),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHON_G4  |         5
BurcatNS |         1
CHON     |         1
""",
)

entry(
    index = 729,
    label = "Cd-HNN",
    group = 
"""
1 * Cd u0 {2,S} {3,S} {4,D}
2   N  u0 {1,S}
3   H  u0 {1,S}
4   N  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.42223,8.6295,9.48229,10.2587,11.4762,12.4219,13.9459],'cal/(mol*K)','+|-',[0.412497,0.44054,0.454336,0.45111,0.438822,0.430948,0.439728]),
        H298 = (22.5918,'kcal/mol','+|-',1.65335),
        S298 = (12.9049,'cal/(mol*K)','+|-',1.26673),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         75
""",
)

entry(
    index = 730,
    label = "Cd-NNN",
    group = 
"""
1 * Cd u0 {2,S} {3,S} {4,D}
2   N  u0 {1,S}
3   N  u0 {1,S}
4   N  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.49931,9.63454,9.8904,9.91133,9.76016,9.6238,9.75725],'cal/(mol*K)','+|-',[0.6255,0.668024,0.688944,0.684051,0.665418,0.653478,0.666792]),
        H298 = (21.4944,'kcal/mol','+|-',2.50709),
        S298 = (-8.34113,'cal/(mol*K)','+|-',1.92084),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHON_G4  |         8
BurcatNS |         2
""",
)

entry(
    index = 731,
    label = "CO-NNOd",
    group = 
"""
1 * CO  u0 {2,S} {3,S} {4,D}
2   N   u0 {1,S}
3   N   u0 {1,S}
4   O2d u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.99486,9.97228,10.2296,10.049,9.8567,9.73768,10.3846],'cal/(mol*K)','+|-',[1.19903,1.28055,1.32065,1.31127,1.27555,1.25266,1.27818]),
        H298 = (-42.4746,'kcal/mol','+|-',4.80588),
        S298 = (11.8277,'cal/(mol*K)','+|-',3.68209),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         1
""",
)

entry(
    index = 732,
    label = "CO-N3sN3sOd",
    group = 
"""
1 * CO  u0 {2,S} {3,S} {4,D}
2   N3s u0 {1,S}
3   N3s u0 {1,S}
4   O2d u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.52948,9.3404,9.67084,9.61718,9.58742,9.58444,10.4295],'cal/(mol*K)','+|-',[1.00029,1.06829,1.10174,1.09392,1.06412,1.04503,1.06632]),
        H298 = (-39.1244,'kcal/mol','+|-',4.00929),
        S298 = (12.2427,'cal/(mol*K)','+|-',3.07176),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         3
""",
)

entry(
    index = 733,
    label = "CO-NN3dOd",
    group = 
"""
1 * CO  u0 {2,S} {3,S} {4,D}
2   N   u0 {1,S}
3   N3d u0 {1,S}
4   O2d u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.28111,9.82546,9.93617,9.84018,9.71743,9.61164,10.032],'cal/(mol*K)','+|-',[0.944701,1.00893,1.04052,1.03313,1.00499,0.986957,1.00706]),
        H298 = (-36.737,'kcal/mol','+|-',3.78649),
        S298 = (14.0051,'cal/(mol*K)','+|-',2.90107),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         3
""",
)

entry(
    index = 734,
    label = "CO-NOO",
    group = 
"""
1 * CO         u0 {2,S} {3,S} {4,D}
2   [O2s,O0sc] u0 {1,S}
3   N          u0 {1,S}
4   O2d        u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.29633,8.75338,9.70517,10.2833,10.9302,11.1733,11.3548],'cal/(mol*K)','+|-',[0.48014,0.512782,0.52884,0.525084,0.510781,0.501616,0.511836]),
        H298 = (-50.016,'kcal/mol','+|-',1.92447),
        S298 = (9.28759,'cal/(mol*K)','+|-',1.47445),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHON_G4  |         8
BurcatNS |         1
""",
)

entry(
    index = 735,
    label = "Cd-HNO",
    group = 
"""
1 * Cd         u0 {2,S} {3,S} {4,D}
2   [O2s,O0sc] u0 {1,S}
3   H          u0 {1,S}
4   N          u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.23211,6.70772,8.30334,9.92021,12.3839,14.1715,16.2809],'cal/(mol*K)','+|-',[0.711305,0.759662,0.783451,0.777888,0.756699,0.743121,0.758261]),
        H298 = (13.2007,'kcal/mol','+|-',2.85101),
        S298 = (11.6454,'cal/(mol*K)','+|-',2.18434),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         7
""",
)

entry(
    index = 736,
    label = "Cd-HN3dO",
    group = 
"""
1 * Cd         u0 {2,S} {3,S} {4,D}
2   [O2s,O0sc] u0 {1,S}
3   H          u0 {1,S}
4   N3d        u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.7762,7.33959,8.86424,10.2721,12.2602,13.6883,15.302],'cal/(mol*K)','+|-',[0.506326,0.540749,0.557683,0.553722,0.538639,0.528974,0.539751]),
        H298 = (12.2349,'kcal/mol','+|-',2.02943),
        S298 = (12.1327,'cal/(mol*K)','+|-',1.55487),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         22
""",
)

entry(
    index = 737,
    label = "Cd-HNdOH",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   O2s u0 {1,S} {5,S}
3   H   u0 {1,S}
4   N3d u0 {1,D}
5   H   u0 {2,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.20445,6.0074,7.82167,9.59567,12.2128,14.0887,16.186],'cal/(mol*K)','+|-',[0.680025,0.726256,0.748999,0.74368,0.723423,0.710442,0.724916]),
        H298 = (11.1074,'kcal/mol','+|-',2.72564),
        S298 = (11.8829,'cal/(mol*K)','+|-',2.08828),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         16
""",
)

entry(
    index = 738,
    label = "Cd-NNO",
    group = 
"""
1 * Cd         u0 {2,S} {3,S} {4,D}
2   [O2s,O0sc] u0 {1,S}
3   N          u0 {1,S}
4   N          u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.0612,7.40677,8.29289,9.04073,10.0482,10.7668,11.4898],'cal/(mol*K)','+|-',[0.623592,0.665986,0.686842,0.681964,0.663388,0.651485,0.664758]),
        H298 = (11.0621,'kcal/mol','+|-',2.49944),
        S298 = (-7.56197,'cal/(mol*K)','+|-',1.91498),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         12
""",
)

entry(
    index = 739,
    label = "Cd-OONd",
    group = 
"""
1 * Cd         u0 {2,S} {3,S} {4,D}
2   [O2s,O0sc] u0 {1,S}
3   [O2s,O0sc] u0 {1,S}
4   N          u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.24913,5.30554,6.45718,7.79575,9.72488,11.1836,12.7793],'cal/(mol*K)','+|-',[1.37523,1.46872,1.51472,1.50396,1.46299,1.43674,1.46602]),
        H298 = (16.1606,'kcal/mol','+|-',5.51212),
        S298 = (-5.52981,'cal/(mol*K)','+|-',4.22318),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         1
""",
)

entry(
    index = 740,
    label = "Cd-OON3d",
    group = 
"""
1 * Cd         u0 {2,S} {3,S} {4,D}
2   [O2s,O0sc] u0 {1,S}
3   [O2s,O0sc] u0 {1,S}
4   N3d        u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.04453,6.21563,8.16977,9.81647,11.8446,13.1079,13.6818],'cal/(mol*K)','+|-',[0.943655,1.00781,1.03937,1.03199,1.00388,0.985864,1.00595]),
        H298 = (-2.62386,'kcal/mol','+|-',3.7823),
        S298 = (-9.01087,'cal/(mol*K)','+|-',2.89786),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         5
""",
)

entry(
    index = 741,
    label = "Cd-CHN",
    group = 
"""
1 * Cd u0 {2,S} {3,S} {4,D}
2   C  u0 {1,S}
3   H  u0 {1,S}
4   N  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.28312,7.46692,8.83544,9.27198,11.012,12.3642,14.1459],'cal/(mol*K)','+|-',[0.55812,0.596063,0.61473,0.610364,0.593738,0.583085,0.594964]),
        H298 = (25.2828,'kcal/mol','+|-',2.23702),
        S298 = (13.7486,'cal/(mol*K)','+|-',1.71392),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         10
""",
)

entry(
    index = 742,
    label = "Cd-HN(CO)",
    group = 
"""
1 * Cd  u0 {2,S} {3,S} {4,D}
2   CO  u0 {1,S} {5,D}
3   H   u0 {1,S}
4   N   u0 {1,D}
5   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.55426,7.24124,8.76127,10.2406,12.4042,14.0882,16.8936],'cal/(mol*K)','+|-',[0.609275,0.650696,0.671073,0.666308,0.648158,0.636528,0.649496]),
        H298 = (21.9277,'kcal/mol','+|-',2.44206),
        S298 = (10.3118,'cal/(mol*K)','+|-',1.87102),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         6
""",
)

entry(
    index = 743,
    label = "Cd-N3dCsH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   N3d u0 {1,D}
3   Cs  u0 {1,S}
4   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.40837,7.42116,8.38401,9.49169,11.3085,12.6402,14.4255],'cal/(mol*K)','+|-',[0.476626,0.509029,0.52497,0.521241,0.507043,0.497945,0.50809]),
        H298 = (24.5509,'kcal/mol','+|-',1.91038),
        S298 = (12.7353,'cal/(mol*K)','+|-',1.46366),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         48
""",
)

entry(
    index = 744,
    label = "Cd-N3dCdH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   N3d u0 {1,D}
3   Cd  u0 {1,S}
4   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.16769,7.73454,9.06038,10.418,12.4592,13.9115,15.5672],'cal/(mol*K)','+|-',[0.504593,0.538897,0.555773,0.551826,0.536795,0.527163,0.537903]),
        H298 = (22.4063,'kcal/mol','+|-',2.02248),
        S298 = (13.6393,'cal/(mol*K)','+|-',1.54955),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         15
""",
)

entry(
    index = 745,
    label = "Cd-N5dcCH",
    group = 
"""
1 * Cd           u0 {2,D} {3,S} {4,S}
2   [N5dc,N5ddc] u0 {1,D}
3   C            u0 {1,S}
4   H            u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.69843,6.76448,7.96237,8.98225,10.8756,12.2054,13.9638],'cal/(mol*K)','+|-',[0.622626,0.664955,0.685778,0.680908,0.662361,0.650476,0.663728]),
        H298 = (19.4758,'kcal/mol','+|-',2.49557),
        S298 = (13.7131,'cal/(mol*K)','+|-',1.91201),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         7
""",
)

entry(
    index = 746,
    label = "Cd-CNNd",
    group = 
"""
1 * Cd u0 {2,S} {3,S} {4,D}
2   C  u0 {1,S}
3   N  u0 {1,S}
4   N  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.98644,8.08457,8.80038,8.78615,9.34172,9.87758,10.992],'cal/(mol*K)','+|-',[0.714584,0.763165,0.787064,0.781474,0.760187,0.746547,0.761757]),
        H298 = (23.8071,'kcal/mol','+|-',2.86415),
        S298 = (-7.97221,'cal/(mol*K)','+|-',2.19441),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         3
""",
)

entry(
    index = 747,
    label = "Cd-CsNNd",
    group = 
"""
1 * Cd u0 {2,S} {3,S} {4,D}
2   Cs u0 {1,S}
3   N  u0 {1,S}
4   N  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.53151,7.43274,7.86822,8.29934,9.00201,9.48175,10.2332],'cal/(mol*K)','+|-',[0.490589,0.523941,0.540348,0.536511,0.521897,0.512532,0.522974]),
        H298 = (21.9592,'kcal/mol','+|-',1.96635),
        S298 = (-7.34892,'cal/(mol*K)','+|-',1.50654),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         14
""",
)

entry(
    index = 748,
    label = "Cd-CdNNd",
    group = 
"""
1 * Cd u0 {2,S} {3,S} {4,D}
2   Cd u0 {1,S}
3   N  u0 {1,S}
4   N  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.88269,7.76195,8.02405,8.35404,8.93654,9.42673,10.2211],'cal/(mol*K)','+|-',[0.822198,0.878094,0.905592,0.899161,0.874669,0.858974,0.876475]),
        H298 = (22.6159,'kcal/mol','+|-',3.29548),
        S298 = (-8.48726,'cal/(mol*K)','+|-',2.52488),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         2
""",
)

entry(
    index = 749,
    label = "Cd-NNCd",
    group = 
"""
1 * Cd       u0 {2,D} {3,S} {4,S}
2   [Cd,Cdd] u0 {1,D}
3   N        u0 {1,S}
4   N        u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.78458,7.61118,7.54963,7.03601,6.21522,5.56294,5.0242],'cal/(mol*K)','+|-',[0.63743,0.680765,0.702084,0.697098,0.67811,0.665942,0.67951]),
        H298 = (16.0964,'kcal/mol','+|-',2.55491),
        S298 = (-14.6581,'cal/(mol*K)','+|-',1.95748),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         10
""",
)

entry(
    index = 750,
    label = "Cd-NNCdd",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D}
3   N   u0 {1,S}
4   N   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.94023,9.5681,9.89459,9.1468,7.71132,6.40607,4.62419],'cal/(mol*K)','+|-',[1.11741,1.19337,1.23075,1.22201,1.18872,1.16739,1.19117]),
        H298 = (20.438,'kcal/mol','+|-',4.47873),
        S298 = (-13.8985,'cal/(mol*K)','+|-',3.43143),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         1
""",
)

entry(
    index = 751,
    label = "Cd-NN(CddOd)",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   N   u0 {1,S}
4   N   u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.30487,7.06782,7.42805,7.21595,6.82604,6.50871,6.53978],'cal/(mol*K)','+|-',[1.38503,1.47919,1.52551,1.51468,1.47342,1.44698,1.47646]),
        H298 = (10.3677,'kcal/mol','+|-',5.5514),
        S298 = (-18.0332,'cal/(mol*K)','+|-',4.25327),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         1
""",
)

entry(
    index = 752,
    label = "Cd-CNO",
    group = 
"""
1 * Cd         u0 {2,S} {3,S} {4,D}
2   C          u0 {1,S}
3   [O2s,O0sc] u0 {1,S}
4   N          u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.60622,4.91974,6.34774,7.67042,9.76141,11.2514,12.7732],'cal/(mol*K)','+|-',[0.659354,0.704179,0.726231,0.721074,0.701432,0.688846,0.702881]),
        H298 = (10.0293,'kcal/mol','+|-',2.64278),
        S298 = (-7.67445,'cal/(mol*K)','+|-',2.0248),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         13
""",
)

entry(
    index = 753,
    label = "Cd-CCN",
    group = 
"""
1 * Cd u0 {2,S} {3,S} {4,D}
2   C  u0 {1,S}
3   C  u0 {1,S}
4   N  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.11779,5.85494,6.52268,7.31523,8.71633,9.59261,10.4457],'cal/(mol*K)','+|-',[0.67601,0.721968,0.744577,0.739289,0.719152,0.706248,0.720636]),
        H298 = (24.379,'kcal/mol','+|-',2.70954),
        S298 = (-6.5901,'cal/(mol*K)','+|-',2.07595),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         4
""",
)

entry(
    index = 754,
    label = "Cd-N3dCsCs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   N3d u0 {1,D}
3   Cs  u0 {1,S}
4   Cs  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.29597,6.02398,6.70143,7.53182,8.92895,9.87331,10.8465],'cal/(mol*K)','+|-',[0.583453,0.623118,0.642631,0.638068,0.620687,0.60955,0.621969]),
        H298 = (23.8057,'kcal/mol','+|-',2.33856),
        S298 = (-7.34045,'cal/(mol*K)','+|-',1.79172),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         7
""",
)

entry(
    index = 755,
    label = "Cds-CNH",
    group = 
"""
1 * Cd u0 {2,D} {3,S} {4,S}
2   C  u0 {1,D}
3   N  u0 {1,S}
4   H  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.00927,7.08731,7.70433,8.02559,8.49362,8.80073,9.40865],'cal/(mol*K)','+|-',[0.361591,0.386173,0.398267,0.395438,0.384667,0.377765,0.385461]),
        H298 = (8.73352,'kcal/mol','+|-',1.44931),
        S298 = (6.21287,'cal/(mol*K)','+|-',1.1104),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library        | Number of Species
CHON_G4        |         18
BurcatNS       |         1
NitrogenCurran |         1
""",
)

entry(
    index = 756,
    label = "Cd-CddNH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D}
3   N   u0 {1,S}
4   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.21149,7.2651,7.81284,8.05333,8.46122,8.75557,9.52088],'cal/(mol*K)','+|-',[0.353783,0.377835,0.389667,0.3869,0.376361,0.369608,0.377138]),
        H298 = (12.3279,'kcal/mol','+|-',1.41801),
        S298 = (8.74965,'cal/(mol*K)','+|-',1.08643),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         19
""",
)

entry(
    index = 757,
    label = "Cd-(CddOd)NH",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   Cdd u0 {1,D} {5,D}
3   N   u0 {1,S}
4   H   u0 {1,S}
5   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.11075,7.17491,8.09895,8.48174,9.00965,9.3474,10.0153],'cal/(mol*K)','+|-',[1.0432,1.11412,1.14901,1.14085,1.10978,1.08987,1.11207]),
        H298 = (-4.12589,'kcal/mol','+|-',4.18131),
        S298 = (3.36251,'cal/(mol*K)','+|-',3.20356),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         6
""",
)

entry(
    index = 758,
    label = "Cd-CdHN3s",
    group = 
"""
1 * Cd  u0 {2,D} {5,S} {6,S}
2   Cd  u0 {1,D} {3,S} {4,S}
3   R   u0 {2,S}
4   R   u0 {2,S}
5   H   u0 {1,S}
6   N3s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.70913,6.39709,6.84953,7.16841,7.77264,8.23152,9.17025],'cal/(mol*K)','+|-',[0.312893,0.334165,0.344629,0.342182,0.332861,0.326889,0.333548]),
        H298 = (13.1531,'kcal/mol','+|-',1.25412),
        S298 = (7.14708,'cal/(mol*K)','+|-',0.960859),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         49
""",
)

entry(
    index = 759,
    label = "Cd-CdHN1s",
    group = 
"""
1 * Cd  u0 {2,D} {5,S} {6,S}
2   Cd  u0 {1,D} {3,S} {4,S}
3   R   u0 {2,S}
4   R   u0 {2,S}
5   H   u0 {1,S}
6   N1s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.87795,4.61405,5.16405,5.59583,6.21763,6.62765,7.27944],'cal/(mol*K)','+|-',[0.48034,0.512995,0.52906,0.525303,0.510994,0.501825,0.512049]),
        H298 = (47.0566,'kcal/mol','+|-',1.92527),
        S298 = (17.4527,'cal/(mol*K)','+|-',1.47507),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
NOx2018 |         1
""",
)

entry(
    index = 760,
    label = "Cds-CCN",
    group = 
"""
1 * Cd u0 {2,D} {3,S} {4,S}
2   C  u0 {1,D}
3   C  u0 {1,S}
4   N  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.22719,5.90843,6.06409,5.93438,5.79871,5.61067,5.409],'cal/(mol*K)','+|-',[0.443429,0.473575,0.488406,0.484937,0.471728,0.463263,0.472702]),
        H298 = (13.4126,'kcal/mol','+|-',1.77733),
        S298 = (-13.184,'cal/(mol*K)','+|-',1.36172),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         8
""",
)

entry(
    index = 761,
    label = "Cd-CdCsN3s",
    group = 
"""
1 * Cd  u0 {2,D} {5,S} {6,S}
2   Cd  u0 {1,D} {3,S} {4,S}
3   R   u0 {2,S}
4   R   u0 {2,S}
5   Cs  u0 {1,S}
6   N3s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.81749,5.23011,5.3326,5.35724,5.55247,5.63631,5.895],'cal/(mol*K)','+|-',[0.414582,0.442767,0.456632,0.45339,0.44104,0.433126,0.44195]),
        H298 = (13.7301,'kcal/mol','+|-',1.6617),
        S298 = (-14.2927,'cal/(mol*K)','+|-',1.27313),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         10
""",
)

entry(
    index = 762,
    label = "Cs",
    group = 
"""
1 * Cs u0
""",
    thermo = 'Cs-CsCsCsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 763,
    label = "CsBrBrBrBr",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Br u0 {1,S}
3   Br u0 {1,S}
4   Br u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([22.2,23.1653,23.8955,24.4568,25.1755,25.5076,25.5753],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (24.6832,'kcal/mol','+|-',2.95256),
        S298 = (90.5794,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHOBr_G4 |         1
""",
)

entry(
    index = 764,
    label = "CsBrBrBrCl",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cl u0 {1,S}
3   Br u0 {1,S}
4   Br u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([21.5643,23.1885,24.0245,24.435,24.881,25.2672,25.6222],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (12.7167,'kcal/mol','+|-',2.95256),
        S298 = (87.733,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library    | Number of Species
CHOClBr_G4 |         1
""",
)

entry(
    index = 765,
    label = "CsBrBrClCl",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cl u0 {1,S}
3   Cl u0 {1,S}
4   Br u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([21.3493,22.6417,23.4642,24.1032,24.9405,25.3519,25.5121],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (0.684689,'kcal/mol','+|-',2.95256),
        S298 = (84.8768,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library    | Number of Species
CHOClBr_G4 |         1
""",
)

entry(
    index = 766,
    label = "CsBrClClCl",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cl u0 {1,S}
3   Cl u0 {1,S}
4   Cl u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([20.8351,22.4252,23.2785,23.9436,24.8249,25.2707,25.483],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-11.415,'kcal/mol','+|-',2.95256),
        S298 = (82.0248,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library    | Number of Species
CHOClBr_G4 |         1
""",
)

entry(
    index = 767,
    label = "CsClClClCl",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cl u0 {1,S}
3   Cl u0 {1,S}
4   Cl u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([20.2225,22.1656,23.0904,23.7842,24.7107,25.1888,25.4468],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-23.6447,'kcal/mol','+|-',2.95256),
        S298 = (79.122,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHOCl_G4 |         1
""",
)

entry(
    index = 768,
    label = "CsBrBrBrF",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   F  u0 {1,S}
3   Br u0 {1,S}
4   Br u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([20.5813,21.9256,22.8586,23.5934,24.5831,25.1034,25.4039],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-30.6932,'kcal/mol','+|-',2.95256),
        S298 = (85.0524,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOFBr_G4 |         1
""",
)

entry(
    index = 769,
    label = "CsBrBrClF",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   F  u0 {1,S}
3   Cl u0 {1,S}
4   Br u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([20.0652,21.6695,22.6414,23.41,24.4547,25.0156,25.371],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-43.4525,'kcal/mol','+|-',2.95256),
        S298 = (82.1635,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOFClBr_G4 |         1
""",
)

entry(
    index = 770,
    label = "CsBrClClF",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   F  u0 {1,S}
3   Cl u0 {1,S}
4   Cl u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([19.4819,21.4138,22.4447,23.2423,24.3344,24.9305,25.3363],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-56.2348,'kcal/mol','+|-',2.95256),
        S298 = (79.2811,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOFClBr_G4 |         1
""",
)

entry(
    index = 771,
    label = "CsClClClF",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   F  u0 {1,S}
3   Cl u0 {1,S}
4   Cl u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([18.8318,21.0799,22.248,23.0786,24.22,24.8486,25.2953],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-69.1504,'kcal/mol','+|-',2.95256),
        S298 = (76.3671,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOFCl_G4 |         1
""",
)

entry(
    index = 772,
    label = "CsBrBrFF",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   F  u0 {1,S}
3   F  u0 {1,S}
4   Br u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([18.6639,20.5895,21.7368,22.6455,23.9155,24.6395,25.2097],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-90.989,'kcal/mol','+|-',2.95256),
        S298 = (79.3068,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOFBr_G4 |         1
""",
)

entry(
    index = 773,
    label = "CsBrClFF",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   F  u0 {1,S}
3   F  u0 {1,S}
4   Cl u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([18.028,20.262,21.5329,22.4736,23.7937,24.553,25.1701],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-104.398,'kcal/mol','+|-',2.95256),
        S298 = (76.397,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOFClBr_G4 |         1
""",
)

entry(
    index = 774,
    label = "CsClClFF",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   F  u0 {1,S}
3   F  u0 {1,S}
4   Cl u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([17.3593,19.8372,21.3215,22.305,23.6817,24.4741,25.1243],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-117.942,'kcal/mol','+|-',2.95256),
        S298 = (73.468,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOFCl_G4 |         1
""",
)

entry(
    index = 775,
    label = "CsBrFFF",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   F  u0 {1,S}
3   F  u0 {1,S}
4   F  u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([16.5248,18.9875,20.5814,21.6768,23.2336,24.1601,24.9908],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-155.373,'kcal/mol','+|-',2.95256),
        S298 = (73.3976,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOFBr_G4 |         1
""",
)

entry(
    index = 776,
    label = "CsClFFF",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   F  u0 {1,S}
3   F  u0 {1,S}
4   F  u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([15.8883,18.4822,20.2974,21.4942,23.1292,24.0935,24.9411],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-169.46,'kcal/mol','+|-',2.95256),
        S298 = (70.4738,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOFCl_G4 |         1
""",
)

entry(
    index = 777,
    label = "CsFFFF",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   F  u0 {1,S}
3   F  u0 {1,S}
4   F  u0 {1,S}
5   F  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHOF_G4 |         1
""",
)

entry(
    index = 778,
    label = "CsBrBrBrH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   H  u0 {1,S}
3   Br u0 {1,S}
4   Br u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([17.3414,18.8145,19.956,20.8998,22.2982,23.1982,24.2369],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (11.7505,'kcal/mol','+|-',2.95256),
        S298 = (81.4264,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHOBr_G4 |         1
""",
)

entry(
    index = 779,
    label = "CsBrBrClH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   H  u0 {1,S}
3   Cl u0 {1,S}
4   Br u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([16.9025,18.5486,19.7307,20.7085,22.1589,23.0942,24.1796],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-0.0714357,'kcal/mol','+|-',2.95256),
        S298 = (78.6122,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library    | Number of Species
CHOClBr_G4 |         1
""",
)

entry(
    index = 780,
    label = "CsBrClClH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   H  u0 {1,S}
3   Cl u0 {1,S}
4   Cl u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([16.4103,18.2942,19.5216,20.5204,22.0122,22.9858,24.1415],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-12.0198,'kcal/mol','+|-',2.95256),
        S298 = (75.7978,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library    | Number of Species
CHOClBr_G4 |         1
""",
)

entry(
    index = 781,
    label = "CsClClClH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   H  u0 {1,S}
3   Cl u0 {1,S}
4   Cl u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([15.8415,17.951,19.3099,20.338,21.8768,22.8853,24.0947],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-24.1628,'kcal/mol','+|-',2.95256),
        S298 = (72.9466,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHOCl_G4 |         1
""",
)

entry(
    index = 782,
    label = "CsBrBrFH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   H  u0 {1,S}
3   F  u0 {1,S}
4   Br u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([15.764,17.5897,18.9005,19.9895,21.63,22.7153,24.0264],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-41.743,'kcal/mol','+|-',2.95256),
        S298 = (75.8694,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOFBr_G4 |         1
""",
)

entry(
    index = 783,
    label = "CsBrClFH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   H  u0 {1,S}
3   F  u0 {1,S}
4   Cl u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([15.2097,17.2605,18.6731,19.7903,21.4797,22.6048,23.981],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-54.6156,'kcal/mol','+|-',2.95256),
        S298 = (73.0213,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOFClBr_G4 |         1
""",
)

entry(
    index = 784,
    label = "CsClClFH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   H  u0 {1,S}
3   F  u0 {1,S}
4   Cl u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([14.6377,16.8328,18.4205,19.5917,21.347,22.5116,23.9284],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-67.6634,'kcal/mol','+|-',2.95256),
        S298 = (70.1509,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOFCl_G4 |         1
""",
)

entry(
    index = 785,
    label = "CsBrFFH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   H  u0 {1,S}
3   F  u0 {1,S}
4   F  u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([13.9378,16.071,17.7129,18.9804,20.9017,22.1953,23.7984],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-101.149,'kcal/mol','+|-',2.95256),
        S298 = (70.1269,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOFBr_G4 |         1
""",
)

entry(
    index = 786,
    label = "CsClFFH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   H  u0 {1,S}
3   F  u0 {1,S}
4   F  u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([13.4197,15.5913,17.3673,18.7412,20.7724,22.1199,23.7446],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-115.092,'kcal/mol','+|-',2.95256),
        S298 = (67.2325,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOFCl_G4 |         1
""",
)

entry(
    index = 787,
    label = "CsFFFH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   H  u0 {1,S}
3   F  u0 {1,S}
4   F  u0 {1,S}
5   F  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHOF_G4 |         1
""",
)

entry(
    index = 788,
    label = "CsBrBrHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   H  u0 {1,S}
3   H  u0 {1,S}
4   Br u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([13.1955,15.0566,16.5168,17.7072,19.5902,20.9458,22.8923],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (1.10843,'kcal/mol','+|-',2.95256),
        S298 = (71.5649,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHOBr_G4 |         1
""",
)

entry(
    index = 789,
    label = "CsBrClHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   H  u0 {1,S}
3   H  u0 {1,S}
4   Cl u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([12.7529,14.6412,16.1975,17.4575,19.4322,20.8422,22.839],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-10.4491,'kcal/mol','+|-',2.95256),
        S298 = (68.7931,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library    | Number of Species
CHOClBr_G4 |         1
""",
)

entry(
    index = 790,
    label = "CsClClHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   H  u0 {1,S}
3   H  u0 {1,S}
4   Cl u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([12.3459,14.2143,15.8271,17.1717,19.2681,20.7457,22.7888],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-22.2439,'kcal/mol','+|-',2.95256),
        S298 = (65.9908,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHOCl_G4 |         1
""",
)

entry(
    index = 791,
    label = "CsBrFHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   H  u0 {1,S}
3   H  u0 {1,S}
4   F  u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([11.7798,13.6441,15.278,16.6424,18.812,20.3959,22.6805],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-50.061,'kcal/mol','+|-',2.95256),
        S298 = (66.038,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOFBr_G4 |         1
""",
)

entry(
    index = 792,
    label = "CsClFHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   H  u0 {1,S}
3   H  u0 {1,S}
4   F  u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([11.3221,13.2296,14.9546,16.3944,18.635,20.272,22.631],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-62.8507,'kcal/mol','+|-',2.95256),
        S298 = (63.2032,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOFCl_G4 |         1
""",
)

entry(
    index = 793,
    label = "CsFFHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   H  u0 {1,S}
3   H  u0 {1,S}
4   F  u0 {1,S}
5   F  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHOF_G4 |         1
""",
)

entry(
    index = 794,
    label = "CsBrHHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   H  u0 {1,S}
3   H  u0 {1,S}
4   H  u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([10.1894,11.855,13.416,14.7962,17.0822,18.8397,21.6182],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-8.40784,'kcal/mol','+|-',2.95256),
        S298 = (60.9239,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHOBr_G4 |         1
""",
)

entry(
    index = 795,
    label = "CsClHHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   H  u0 {1,S}
3   H  u0 {1,S}
4   H  u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.83687,11.4168,13.0205,14.4951,16.8804,18.6969,21.5672],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-19.5024,'kcal/mol','+|-',2.95256),
        S298 = (58.1653,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHOCl_G4 |         1
""",
)

entry(
    index = 796,
    label = "CsFHHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   H  u0 {1,S}
3   H  u0 {1,S}
4   H  u0 {1,S}
5   F  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHOF_G4 |         1
""",
)

entry(
    index = 797,
    label = "CsBrBrBrO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   O  u0 {1,S}
3   Br u0 {1,S}
4   Br u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([17.3633,18.4532,19.2425,19.8333,20.7539,21.1785,21.3004],'cal/(mol*K)','+|-',[0.151432,0.17404,0.177912,0.175302,0.15201,0.131108,0.10685]),
        H298 = (14.148,'kcal/mol','+|-',0.539843),
        S298 = (59.717,'cal/(mol*K)','+|-',0.354351),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library    | Number of Species
CHOBr_G4   |         20
CHOFBr_G4  |         9
CHOClBr_G4 |         2
""",
)

entry(
    index = 798,
    label = "CsBrBrClO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   O  u0 {1,S}
3   Cl u0 {1,S}
4   Br u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([17.1996,18.3305,19.0709,19.6222,20.5007,20.9242,21.0633],'cal/(mol*K)','+|-',[0.20875,0.239916,0.245253,0.241655,0.209547,0.180733,0.147293]),
        H298 = (1.59819,'kcal/mol','+|-',0.744177),
        S298 = (55.4811,'cal/(mol*K)','+|-',0.488475),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOFClBr_G4 |         2
CHOClBr_G4  |         14
""",
)

entry(
    index = 799,
    label = "CsBrClClO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   O  u0 {1,S}
3   Cl u0 {1,S}
4   Cl u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([16.7744,18.0576,18.8404,19.4251,20.3525,20.8157,21.0078],'cal/(mol*K)','+|-',[0.20875,0.239916,0.245253,0.241655,0.209547,0.180733,0.147293]),
        H298 = (-10.6892,'kcal/mol','+|-',0.744177),
        S298 = (52.9253,'cal/(mol*K)','+|-',0.488475),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOFClBr_G4 |         2
CHOClBr_G4  |         14
""",
)

entry(
    index = 800,
    label = "CsClClClO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   O  u0 {1,S}
3   Cl u0 {1,S}
4   Cl u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([16.0862,17.8162,18.8603,19.5755,20.5155,20.9764,21.2635],'cal/(mol*K)','+|-',[0.120908,0.13896,0.142051,0.139967,0.12137,0.104681,0.0853123]),
        H298 = (-22.6274,'kcal/mol','+|-',0.431028),
        S298 = (50.891,'cal/(mol*K)','+|-',0.282925),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library    | Number of Species
CHOCl_G4   |         39
CHOFCl_G4  |         1
CHOClBr_G4 |         3
""",
)

entry(
    index = 801,
    label = "CsBrBrFO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   O  u0 {1,S}
3   F  u0 {1,S}
4   Br u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([16.5841,17.747,18.5454,19.1455,20.0919,20.5791,20.8904],'cal/(mol*K)','+|-',[0.167283,0.192258,0.196535,0.193652,0.167922,0.144831,0.118034]),
        H298 = (-44.0282,'kcal/mol','+|-',0.596351),
        S298 = (54.0422,'cal/(mol*K)','+|-',0.391443),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOFBr_G4 |         25
""",
)

entry(
    index = 802,
    label = "CsBrClFO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   O  u0 {1,S}
3   F  u0 {1,S}
4   Cl u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([15.7639,17.1399,18.0433,18.7404,19.8482,20.4387,20.8633],'cal/(mol*K)','+|-',[0.222443,0.255652,0.261339,0.257505,0.223292,0.192587,0.156954]),
        H298 = (-57.2202,'kcal/mol','+|-',0.792989),
        S298 = (50.2375,'cal/(mol*K)','+|-',0.520515),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOFClBr_G4 |         14
""",
)

entry(
    index = 803,
    label = "CsClClFO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   O  u0 {1,S}
3   F  u0 {1,S}
4   Cl u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([15.1762,16.9105,17.9699,18.7582,19.8987,20.4579,20.7724],'cal/(mol*K)','+|-',[0.214755,0.246817,0.252308,0.248606,0.215575,0.185932,0.15153]),
        H298 = (-70.2117,'kcal/mol','+|-',0.765585),
        S298 = (47.5564,'cal/(mol*K)','+|-',0.502527),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOFCl_G4   |         12
CHOFClBr_G4 |         3
""",
)

entry(
    index = 804,
    label = "CsBrFFO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   O  u0 {1,S}
3   F  u0 {1,S}
4   F  u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([14.6601,16.2223,17.2998,18.1427,19.3657,19.9912,20.5078],'cal/(mol*K)','+|-',[0.170802,0.196302,0.200669,0.197725,0.171454,0.147878,0.120517]),
        H298 = (-106.939,'kcal/mol','+|-',0.608895),
        S298 = (48.821,'cal/(mol*K)','+|-',0.399677),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOFBr_G4 |         24
""",
)

entry(
    index = 805,
    label = "CsClFFO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   O  u0 {1,S}
3   F  u0 {1,S}
4   F  u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([13.9139,15.8223,17.0385,17.9331,19.2329,19.9163,20.452],'cal/(mol*K)','+|-',[0.214755,0.246817,0.252308,0.248606,0.215575,0.185932,0.15153]),
        H298 = (-120.15,'kcal/mol','+|-',0.765585),
        S298 = (44.5499,'cal/(mol*K)','+|-',0.502527),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOFCl_G4   |         12
CHOFClBr_G4 |         3
""",
)

entry(
    index = 806,
    label = "CsFFFO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   O  u0 {1,S}
3   F  u0 {1,S}
4   F  u0 {1,S}
5   F  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOF_G4   |         41
CHOFCl_G4 |         3
CHOFBr_G4 |         8
""",
)

entry(
    index = 807,
    label = "CsBrBrHO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   O  u0 {1,S}
3   H  u0 {1,S}
4   Br u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([13.4368,15.0833,16.3509,17.3216,18.7127,19.4564,20.094],'cal/(mol*K)','+|-',[0.0991652,0.11397,0.116505,0.114796,0.0995438,0.0858556,0.0699704]),
        H298 = (2.51501,'kcal/mol','+|-',0.353516),
        S298 = (49.3094,'cal/(mol*K)','+|-',0.232046),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOBr_G4    |         31
CHOFClBr_G4 |         1
CHOFBr_G4   |         24
CHOClBr_G4  |         9
""",
)

entry(
    index = 808,
    label = "CsBrClHO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   O  u0 {1,S}
3   H  u0 {1,S}
4   Cl u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([13.2538,14.8788,16.0494,16.9491,18.3321,19.1356,20.0602],'cal/(mol*K)','+|-',[0.142782,0.164099,0.167749,0.165288,0.143327,0.123618,0.100746]),
        H298 = (-9.48474,'kcal/mol','+|-',0.509005),
        S298 = (45.3805,'cal/(mol*K)','+|-',0.334109),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOFClBr_G4 |         9
CHOClBr_G4  |         25
""",
)

entry(
    index = 809,
    label = "CsClClHO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   O  u0 {1,S}
3   H  u0 {1,S}
4   Cl u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([12.6157,14.4743,15.8317,16.8555,18.3112,19.1416,20.0642],'cal/(mol*K)','+|-',[0.103929,0.119445,0.122102,0.12031,0.104325,0.0899798,0.0733315]),
        H298 = (-22.3052,'kcal/mol','+|-',0.370497),
        S298 = (42.7759,'cal/(mol*K)','+|-',0.243193),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library    | Number of Species
CHOCl_G4   |         41
CHOFCl_G4  |         9
CHOClBr_G4 |         9
""",
)

entry(
    index = 810,
    label = "CsBrFHO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   O  u0 {1,S}
3   H  u0 {1,S}
4   F  u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([12.3599,14.1063,15.3815,16.377,17.8719,18.747,19.7499],'cal/(mol*K)','+|-',[0.139186,0.159965,0.163524,0.161125,0.139717,0.120505,0.0982086]),
        H298 = (-53.9333,'kcal/mol','+|-',0.496185),
        S298 = (43.7068,'cal/(mol*K)','+|-',0.325694),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOFBr_G4 |         36
""",
)

entry(
    index = 811,
    label = "CsClFHO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   O  u0 {1,S}
3   H  u0 {1,S}
4   F  u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([11.8079,13.6445,14.9954,16.0418,17.6237,18.5886,19.733],'cal/(mol*K)','+|-',[0.142991,0.164339,0.167995,0.16553,0.143537,0.1238,0.100894]),
        H298 = (-67.2347,'kcal/mol','+|-',0.509752),
        S298 = (39.9529,'cal/(mol*K)','+|-',0.334599),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOFCl_G4   |         25
CHOFClBr_G4 |         9
""",
)

entry(
    index = 812,
    label = "CsFFHO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   O  u0 {1,S}
3   H  u0 {1,S}
4   F  u0 {1,S}
5   F  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOF_G4     |         41
CHOFCl_G4   |         9
CHOFClBr_G4 |         1
CHOFBr_G4   |         18
""",
)

entry(
    index = 813,
    label = "CsBrHHO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   O  u0 {1,S}
3   H  u0 {1,S}
4   H  u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.45313,11.3574,12.8899,14.0962,15.9092,17.092,18.7982],'cal/(mol*K)','+|-',[0.0786121,0.0903486,0.0923584,0.0910034,0.0789122,0.0680611,0.0554683]),
        H298 = (-5.79178,'kcal/mol','+|-',0.280246),
        S298 = (39.4417,'cal/(mol*K)','+|-',0.183952),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOBr_G4    |         38
CHOFClBr_G4 |         10
CHOFBr_G4   |         37
CHOClBr_G4  |         24
""",
)

entry(
    index = 814,
    label = "CsClHHO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   O  u0 {1,S}
3   H  u0 {1,S}
4   H  u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.1541,11.1017,12.6495,13.8529,15.6712,16.8708,18.6265],'cal/(mol*K)','+|-',[0.0891754,0.102489,0.104769,0.103232,0.0895159,0.0772067,0.0629217]),
        H298 = (-17.9509,'kcal/mol','+|-',0.317903),
        S298 = (36.5727,'cal/(mol*K)','+|-',0.20867),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library    | Number of Species
CHOCl_G4   |         41
CHOFCl_G4  |         24
CHOClBr_G4 |         18
""",
)

entry(
    index = 815,
    label = "CsFHHO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   O  u0 {1,S}
3   H  u0 {1,S}
4   H  u0 {1,S}
5   F  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOF_G4     |         41
CHOFCl_G4   |         18
CHOFClBr_G4 |         7
CHOFBr_G4   |         24
""",
)

entry(
    index = 816,
    label = "CsBrBrOO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   O  u0 {1,S}
3   O  u0 {1,S}
4   Br u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([14.4317,16.2897,17.2922,17.8377,18.2847,18.3798,17.7485],'cal/(mol*K)','+|-',[0.208809,0.239983,0.245322,0.241723,0.209606,0.180783,0.147334]),
        H298 = (-3.60237,'kcal/mol','+|-',0.744386),
        S298 = (25.836,'cal/(mol*K)','+|-',0.488612),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library    | Number of Species
CHOBr_G4   |         7
CHOFBr_G4  |         6
CHOClBr_G4 |         3
""",
)

entry(
    index = 817,
    label = "CsBrClOO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   O  u0 {1,S}
3   O  u0 {1,S}
4   Cl u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([13.8971,15.8102,16.833,17.4294,17.9904,18.1243,17.5525],'cal/(mol*K)','+|-',[0.262923,0.302176,0.308898,0.304366,0.263926,0.227634,0.185517]),
        H298 = (-16.1987,'kcal/mol','+|-',0.937297),
        S298 = (22.9283,'cal/(mol*K)','+|-',0.615238),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOFClBr_G4 |         3
CHOClBr_G4  |         7
""",
)

entry(
    index = 818,
    label = "CsClClOO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   O  u0 {1,S}
3   O  u0 {1,S}
4   Cl u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([13.6488,15.9971,17.2143,17.9121,18.4889,18.5564,17.7902],'cal/(mol*K)','+|-',[0.197881,0.227424,0.232483,0.229072,0.198637,0.171322,0.139624]),
        H298 = (-28.0095,'kcal/mol','+|-',0.70543),
        S298 = (19.6958,'cal/(mol*K)','+|-',0.463041),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOCl_G4    |         11
CHOFCl_G4   |         3
CHOFClBr_G4 |         1
CHOClBr_G4  |         3
""",
)

entry(
    index = 819,
    label = "CsBrFOO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   O  u0 {1,S}
3   O  u0 {1,S}
4   F  u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([12.4113,14.4445,15.6693,16.4469,17.2319,17.4979,17.1277],'cal/(mol*K)','+|-',[0.277597,0.319041,0.326138,0.321353,0.278657,0.240339,0.195871]),
        H298 = (-64.6008,'kcal/mol','+|-',0.989609),
        S298 = (20.3961,'cal/(mol*K)','+|-',0.649576),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOFBr_G4 |         9
""",
)

entry(
    index = 820,
    label = "CsClFOO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   O  u0 {1,S}
3   O  u0 {1,S}
4   F  u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([11.8322,14.3201,15.7602,16.6404,17.5165,17.8137,17.414],'cal/(mol*K)','+|-',[0.26301,0.302277,0.309001,0.304468,0.264015,0.22771,0.185579]),
        H298 = (-77.4106,'kcal/mol','+|-',0.93761),
        S298 = (16.4445,'cal/(mol*K)','+|-',0.615444),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOFCl_G4   |         7
CHOFClBr_G4 |         3
""",
)

entry(
    index = 821,
    label = "CsFFOO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   O  u0 {1,S}
3   O  u0 {1,S}
4   F  u0 {1,S}
5   F  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([10.5929,13.2366,14.7874,15.7616,16.7963,17.2268,17.0373],'cal/(mol*K)','+|-',[0.198271,0.227872,0.232941,0.229523,0.199028,0.17166,0.139899]),
        H298 = (-127.083,'kcal/mol','+|-',0.706818),
        S298 = (14.3851,'cal/(mol*K)','+|-',0.463953),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOF_G4   |         11
CHOFCl_G4 |         3
CHOFBr_G4 |         4
""",
)

entry(
    index = 822,
    label = "CsBrHOO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   O  u0 {1,S}
3   O  u0 {1,S}
4   H  u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([11.2611,13.8442,15.1607,15.907,16.6399,16.9367,16.7074],'cal/(mol*K)','+|-',[0.161597,0.185723,0.189855,0.187069,0.162214,0.139908,0.114022]),
        H298 = (-12.7152,'kcal/mol','+|-',0.576081),
        S298 = (13.7418,'cal/(mol*K)','+|-',0.378137),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOBr_G4    |         9
CHOFClBr_G4 |         1
CHOFBr_G4   |         11
CHOClBr_G4  |         6
""",
)

entry(
    index = 823,
    label = "CsClHOO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   O  u0 {1,S}
3   O  u0 {1,S}
4   H  u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([10.7034,13.3907,14.7394,15.4767,16.2209,16.5293,16.4058],'cal/(mol*K)','+|-',[0.175566,0.201777,0.206266,0.20324,0.176236,0.152002,0.123878]),
        H298 = (-25.2043,'kcal/mol','+|-',0.625878),
        S298 = (11.4175,'cal/(mol*K)','+|-',0.410824),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOCl_G4    |         11
CHOFCl_G4   |         6
CHOFClBr_G4 |         2
CHOClBr_G4  |         4
""",
)

entry(
    index = 824,
    label = "CsFHOO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   O  u0 {1,S}
3   O  u0 {1,S}
4   H  u0 {1,S}
5   F  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.8292,11.2462,12.7308,13.7012,14.8047,15.4098,15.8185],'cal/(mol*K)','+|-',[0.188677,0.216846,0.22167,0.218418,0.189398,0.163354,0.13313]),
        H298 = (-71.8804,'kcal/mol','+|-',0.672619),
        S298 = (9.72474,'cal/(mol*K)','+|-',0.441505),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOF_G4   |         11
CHOFCl_G4 |         4
CHOFBr_G4 |         5
""",
)

entry(
    index = 825,
    label = "CsBrOOO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   O  u0 {1,S}
3   O  u0 {1,S}
4   O  u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([11.9169,15.3667,16.5941,16.8867,16.7723,16.117,14.086],'cal/(mol*K)'),
        H298 = (-21.5933,'kcal/mol'),
        S298 = (-8.38973,'cal/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
library:CHOBr_G4 label:OC(O)(O)Br smiles:OC(O)(O)Br H298:-140.04 kcal/mol
""",
)

entry(
    index = 826,
    label = "CsClOOO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   O  u0 {1,S}
3   O  u0 {1,S}
4   O  u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([12.328,14.7653,15.296,15.4042,15.1769,14.6381,13.2612],'cal/(mol*K)'),
        H298 = (-34.3228,'kcal/mol'),
        S298 = (-9.60223,'cal/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
library:CHOCl_G4 label:OC(O)(O)Cl smiles:OC(O)(O)Cl H298:-152.77 kcal/mol
""",
)

entry(
    index = 827,
    label = "CsFOOO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   O  u0 {1,S}
3   O  u0 {1,S}
4   O  u0 {1,S}
5   F  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.97097,11.1967,12.478,13.3169,14.153,14.2521,13.3172],'cal/(mol*K)'),
        H298 = (-86.4118,'kcal/mol'),
        S298 = (-11.6373,'cal/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
library:CHOF_G4 label:OC(O)(O)F smiles:OC(O)(O)F H298:-204.86 kcal/mol
""",
)

entry(
    index = 828,
    label = "CsBrBrBrC",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   Br u0 {1,S}
4   Br u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = 'CsBrBrBrCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 829,
    label = "CsBrBrBrCs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cs u0 {1,S}
3   Br u0 {1,S}
4   Br u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([17.4122,18.53,19.363,20.053,21.0386,21.5434,21.503],'cal/(mol*K)','+|-',[0.12646,0.145341,0.148574,0.146394,0.126943,0.109487,0.0892298]),
        H298 = (12.1483,'kcal/mol','+|-',0.450821),
        S298 = (57.8864,'cal/(mol*K)','+|-',0.295917),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOBr_G4    |         34
CHOFClBr_G4 |         1
CHOFBr_G4   |         26
CHOClBr_G4  |         8
""",
)

entry(
    index = 830,
    label = "CsBrBrBrCd",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S}
3   Br u0 {1,S}
4   Br u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([17.6507,18.7216,19.4984,20.0983,20.8697,21.2687,21.3253],'cal/(mol*K)','+|-',[0.163824,0.188282,0.192471,0.189647,0.164449,0.141836,0.115593]),
        H298 = (14.9465,'kcal/mol','+|-',0.584018),
        S298 = (58.6923,'cal/(mol*K)','+|-',0.383347),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library    | Number of Species
CHOBr_G4   |         13
CHOFBr_G4  |         14
CHOClBr_G4 |         5
""",
)

entry(
    index = 831,
    label = "CsBrBrBrCt",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Ct u0 {1,S}
3   Br u0 {1,S}
4   Br u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([17.8789,19.0202,19.82,20.4219,21.1404,21.4729,20.3221],'cal/(mol*K)','+|-',[0.293773,0.337632,0.345143,0.340079,0.294894,0.254344,0.207285]),
        H298 = (19.5485,'kcal/mol','+|-',1.04728),
        S298 = (61.1521,'cal/(mol*K)','+|-',0.687428),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library    | Number of Species
CHOBr_G4   |         5
CHOFBr_G4  |         2
CHOClBr_G4 |         1
""",
)

entry(
    index = 832,
    label = "CsBrBrBrCO",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   Br  u0 {1,S}
4   Br  u0 {1,S}
5   Br  u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([18.3381,19.3802,19.9425,20.3137,20.7951,21.0537,21.1577],'cal/(mol*K)','+|-',[0.371968,0.427501,0.437011,0.430599,0.373388,0.322044,0.262458]),
        H298 = (19.4336,'kcal/mol','+|-',1.32603),
        S298 = (60.0329,'cal/(mol*K)','+|-',0.870403),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOBr_G4  |         4
CHOFBr_G4 |         1
""",
)

entry(
    index = 833,
    label = "CsC2sBrBrBr",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 p1 {1,S}
3   Br u0 {1,S}
4   Br u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = 'CsBrBrBrCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 834,
    label = "CsBrBrCCl",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   Cl u0 {1,S}
4   Br u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([17.1505,18.421,19.2811,19.9264,20.8319,21.339,21.3908],'cal/(mol*K)','+|-',[0.109646,0.126016,0.128819,0.126929,0.110065,0.0949301,0.0773659]),
        H298 = (1.88513,'kcal/mol','+|-',0.39088),
        S298 = (54.6003,'cal/(mol*K)','+|-',0.256572),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOFClBr_G4 |         14
CHOClBr_G4  |         50
""",
)

entry(
    index = 835,
    label = "CsBrCClCl",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   Cl u0 {1,S}
4   Cl u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([16.6796,18.0784,19.0151,19.7259,20.7023,21.2349,21.2948],'cal/(mol*K)','+|-',[0.109646,0.126016,0.128819,0.126929,0.110065,0.0949301,0.0773659]),
        H298 = (-9.91184,'kcal/mol','+|-',0.39088),
        S298 = (52.0802,'cal/(mol*K)','+|-',0.256572),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOFClBr_G4 |         14
CHOClBr_G4  |         50
""",
)

entry(
    index = 836,
    label = "CsCClClCl",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   Cl u0 {1,S}
4   Cl u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = 'CsClClClCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 837,
    label = "CsClClClCs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cs u0 {1,S}
3   Cl u0 {1,S}
4   Cl u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([16.0977,17.6655,18.7188,19.5381,20.6723,21.2955,21.5677],'cal/(mol*K)','+|-',[0.103232,0.118644,0.121283,0.119503,0.103626,0.0893762,0.0728396]),
        H298 = (-23.5535,'kcal/mol','+|-',0.368012),
        S298 = (49.8015,'cal/(mol*K)','+|-',0.241561),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library                  | Number of Species
CHOCl_G4                 |         159
Chlorinated_Hydrocarbons |         2
CHOFCl_G4                |         5
CHOClBr_G4               |         8
""",
)

entry(
    index = 838,
    label = "CsClClClCd",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S}
3   Cl u0 {1,S}
4   Cl u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([16.5485,18.1615,19.1699,19.9271,20.8583,21.3585,21.568],'cal/(mol*K)','+|-',[0.132198,0.151935,0.155315,0.153036,0.132703,0.114455,0.0932784]),
        H298 = (-19.1204,'kcal/mol','+|-',0.471276),
        S298 = (50.1525,'cal/(mol*K)','+|-',0.309344),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library    | Number of Species
CHOCl_G4   |         56
CHOFCl_G4  |         2
CHOClBr_G4 |         3
""",
)

entry(
    index = 839,
    label = "CsClClClCt",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Ct u0 {1,S}
3   Cl u0 {1,S}
4   Cl u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([16.4022,17.977,18.9309,19.5531,20.45,20.9133,19.9295],'cal/(mol*K)','+|-',[0.293513,0.337333,0.344837,0.339778,0.294633,0.254119,0.207101]),
        H298 = (-13.3593,'kcal/mol','+|-',1.04635),
        S298 = (52.6026,'cal/(mol*K)','+|-',0.686818),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOCl_G4  |         7
CHOFCl_G4 |         1
""",
)

entry(
    index = 840,
    label = "CsClClClCO",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   Cl  u0 {1,S}
4   Cl  u0 {1,S}
5   Cl  u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([17.0668,18.4444,19.2517,19.8568,20.6651,21.0835,21.284],'cal/(mol*K)','+|-',[0.24002,0.275854,0.281991,0.277854,0.240937,0.207806,0.169357]),
        H298 = (-14.0361,'kcal/mol','+|-',0.855652),
        S298 = (50.1884,'cal/(mol*K)','+|-',0.561647),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library    | Number of Species
CHOCl_G4   |         8
CHOClBr_G4 |         1
""",
)

entry(
    index = 841,
    label = "CsC2sClClCl",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 p1 {1,S}
3   Cl u0 {1,S}
4   Cl u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([16.314,17.5106,18.3821,19.1881,20.2325,20.8086,20.941],'cal/(mol*K)','+|-',[0.908394,1.04401,1.06724,1.05158,0.911862,0.786473,0.640958]),
        H298 = (-18.1,'kcal/mol','+|-',3.23835),
        S298 = (52.6739,'cal/(mol*K)','+|-',2.12564),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHOCl_G4 |         1
""",
)

entry(
    index = 842,
    label = "CsBrBrCF",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   F  u0 {1,S}
4   Br u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([16.7544,18.0281,18.8845,19.5378,20.4598,20.9579,20.9593],'cal/(mol*K)','+|-',[0.0902607,0.103736,0.106044,0.104488,0.0906053,0.0781463,0.0636875]),
        H298 = (-40.6379,'kcal/mol','+|-',0.321772),
        S298 = (52.7126,'cal/(mol*K)','+|-',0.21121),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOFBr_G4 |         105
""",
)

entry(
    index = 843,
    label = "CsBrCClF",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   F  u0 {1,S}
4   Cl u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([15.9075,17.4004,18.4071,19.1786,20.2529,20.8563,21.042],'cal/(mol*K)','+|-',[0.122479,0.140765,0.143897,0.141785,0.122947,0.106041,0.0864209]),
        H298 = (-53.176,'kcal/mol','+|-',0.436629),
        S298 = (49.0834,'cal/(mol*K)','+|-',0.286602),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOFClBr_G4 |         50
""",
)

entry(
    index = 844,
    label = "CsCClClF",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   F  u0 {1,S}
4   Cl u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([15.4764,17.0706,18.1014,18.9072,20.0563,20.7001,20.8521],'cal/(mol*K)','+|-',[0.122072,0.140296,0.143417,0.141313,0.122538,0.105688,0.0861331]),
        H298 = (-66.044,'kcal/mol','+|-',0.435175),
        S298 = (46.9087,'cal/(mol*K)','+|-',0.285647),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOFCl_G4   |         39
CHOFClBr_G4 |         12
""",
)

entry(
    index = 845,
    label = "CsBrCFF",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   F  u0 {1,S}
4   F  u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([14.9331,16.4747,17.5504,18.3961,19.6055,20.2969,20.5681],'cal/(mol*K)','+|-',[0.0897499,0.103149,0.105444,0.103897,0.0900926,0.0777041,0.0633271]),
        H298 = (-99.8076,'kcal/mol','+|-',0.319951),
        S298 = (48.1267,'cal/(mol*K)','+|-',0.210015),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOFBr_G4 |         106
""",
)

entry(
    index = 846,
    label = "CsCClFF",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   F  u0 {1,S}
4   F  u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([14.0665,15.9019,17.1353,18.0838,19.4211,20.2032,20.5969],'cal/(mol*K)','+|-',[0.122072,0.140296,0.143417,0.141313,0.122538,0.105688,0.0861331]),
        H298 = (-113.473,'kcal/mol','+|-',0.435175),
        S298 = (44.5773,'cal/(mol*K)','+|-',0.285647),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOFCl_G4   |         39
CHOFClBr_G4 |         12
""",
)

entry(
    index = 847,
    label = "CsCFFF",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   F  u0 {1,S}
4   F  u0 {1,S}
5   F  u0 {1,S}
""",
    thermo = 'CsCsFFF',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 848,
    label = "CsCsFFF",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cs u0 {1,S}
3   F  u0 {1,S}
4   F  u0 {1,S}
5   F  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOF_G4     |         163
CHOFCl_G4   |         7
CHOFClBr_G4 |         1
CHOFBr_G4   |         39
""",
)

entry(
    index = 849,
    label = "CsCdFFF",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S}
3   F  u0 {1,S}
4   F  u0 {1,S}
5   F  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOF_G4   |         56
CHOFBr_G4 |         9
""",
)

entry(
    index = 850,
    label = "CsCtFFF",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Ct u0 {1,S}
3   F  u0 {1,S}
4   F  u0 {1,S}
5   F  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOF_G4   |         7
CHOFBr_G4 |         2
""",
)

entry(
    index = 851,
    label = "CsCOFFF",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   F   u0 {1,S}
4   F   u0 {1,S}
5   F   u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOF_G4   |         8
CHOFCl_G4 |         1
CHOFBr_G4 |         3
""",
)

entry(
    index = 852,
    label = "CsC2sFFF",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 p1 {1,S}
3   F  u0 {1,S}
4   F  u0 {1,S}
5   F  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHOF_G4 |         1
""",
)

entry(
    index = 853,
    label = "CsBrBrCH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   H  u0 {1,S}
4   Br u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([13.6551,15.0942,16.256,17.22,18.5476,19.3769,19.9874],'cal/(mol*K)','+|-',[0.0828955,0.0952715,0.0973908,0.0959619,0.0832119,0.0717696,0.0584906]),
        H298 = (4.79193,'kcal/mol','+|-',0.295515),
        S298 = (49.1847,'cal/(mol*K)','+|-',0.193975),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOBr_G4    |         46
CHOFClBr_G4 |         4
CHOFBr_G4   |         44
CHOClBr_G4  |         21
""",
)

entry(
    index = 854,
    label = "CsBrBrCsH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cs u0 {1,S}
3   H  u0 {1,S}
4   Br u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([13.8055,15.4986,16.769,17.7972,19.2211,20.013,20.525],'cal/(mol*K)','+|-',[0.065305,0.0750549,0.0767245,0.0755988,0.0655544,0.0565401,0.0460789]),
        H298 = (2.56116,'kcal/mol','+|-',0.232807),
        S298 = (47.8869,'cal/(mol*K)','+|-',0.152814),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOBr_G4    |         75
CHOFClBr_G4 |         10
CHOFBr_G4   |         74
CHOClBr_G4  |         33
""",
)

entry(
    index = 855,
    label = "CsBrCClH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   H  u0 {1,S}
4   Cl u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([13.2524,14.8939,16.126,17.1406,18.5972,19.4825,20.2423],'cal/(mol*K)','+|-',[0.0673062,0.0773548,0.0790755,0.0779154,0.0675632,0.0582726,0.0474909]),
        H298 = (-7.96045,'kcal/mol','+|-',0.239941),
        S298 = (45.2769,'cal/(mol*K)','+|-',0.157496),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOFClBr_G4 |         53
CHOClBr_G4  |         112
""",
)

entry(
    index = 856,
    label = "CsCClClH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   H  u0 {1,S}
4   Cl u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([12.6667,14.3235,15.6018,16.6377,18.1077,19.0613,19.8792],'cal/(mol*K)','+|-',[0.0821716,0.0944396,0.0965404,0.095124,0.0824854,0.0711429,0.0579799]),
        H298 = (-17.6924,'kcal/mol','+|-',0.292935),
        S298 = (43.9815,'cal/(mol*K)','+|-',0.192281),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library    | Number of Species
CHOCl_G4   |         72
CHOFCl_G4  |         17
CHOClBr_G4 |         18
""",
)

entry(
    index = 857,
    label = "CsCsClClH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cs u0 {1,S}
3   H  u0 {1,S}
4   Cl u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([12.9467,14.8218,16.1741,17.2659,18.7964,19.6715,20.4338],'cal/(mol*K)','+|-',[0.0600193,0.06898,0.0705145,0.0694799,0.0602485,0.0519638,0.0423493]),
        H298 = (-20.9422,'kcal/mol','+|-',0.213964),
        S298 = (42.474,'cal/(mol*K)','+|-',0.140445),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library                  | Number of Species
CHOCl_G4                 |         162
Chlorinated_Hydrocarbons |         2
CHOFCl_G4                |         26
CHOClBr_G4               |         41
""",
)

entry(
    index = 858,
    label = "CsBrCFH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   H  u0 {1,S}
4   F  u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([12.6116,14.3531,15.6497,16.7106,18.2306,19.1515,19.9835],'cal/(mol*K)','+|-',[0.0647582,0.0744264,0.076082,0.0749658,0.0650054,0.0560666,0.045693]),
        H298 = (-48.2752,'kcal/mol','+|-',0.230858),
        S298 = (42.6139,'cal/(mol*K)','+|-',0.151534),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOFBr_G4 |         183
""",
)

entry(
    index = 859,
    label = "CsCClFH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   H  u0 {1,S}
4   F  u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([12.0304,13.7733,15.0641,16.1288,17.7187,18.7304,19.7378],'cal/(mol*K)','+|-',[0.0676875,0.077793,0.0795235,0.0783568,0.0679459,0.0586028,0.0477599]),
        H298 = (-60.7464,'kcal/mol','+|-',0.2413),
        S298 = (39.5409,'cal/(mol*K)','+|-',0.158389),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOFCl_G4   |         105
CHOFClBr_G4 |         58
""",
)

entry(
    index = 860,
    label = "CsCFFH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   H  u0 {1,S}
4   F  u0 {1,S}
5   F  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOF_G4     |         71
CHOFCl_G4   |         14
CHOFClBr_G4 |         2
CHOFBr_G4   |         36
""",
)

entry(
    index = 861,
    label = "CsCsFFH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cs u0 {1,S}
3   H  u0 {1,S}
4   F  u0 {1,S}
5   F  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOF_G4     |         163
CHOFCl_G4   |         38
CHOFClBr_G4 |         9
CHOFBr_G4   |         103
""",
)

entry(
    index = 862,
    label = "CsBrCHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   H  u0 {1,S}
4   H  u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.60714,11.3502,12.7826,14.0053,15.799,17.0943,18.7376],'cal/(mol*K)','+|-',[0.0611201,0.0702451,0.0718077,0.0707542,0.0613535,0.0529168,0.043126]),
        H298 = (-3.61451,'kcal/mol','+|-',0.217888),
        S298 = (39.548,'cal/(mol*K)','+|-',0.143021),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOBr_G4    |         63
CHOFClBr_G4 |         20
CHOFBr_G4   |         76
CHOClBr_G4  |         48
""",
)

entry(
    index = 863,
    label = "CsBrCsHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cs u0 {1,S}
3   H  u0 {1,S}
4   H  u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([10.1467,11.9455,13.3867,14.5914,16.4048,17.6368,19.1929],'cal/(mol*K)','+|-',[0.0484643,0.0556998,0.0569389,0.0561035,0.0486493,0.0419596,0.0341961]),
        H298 = (-5.65225,'kcal/mol','+|-',0.172771),
        S298 = (37.9981,'cal/(mol*K)','+|-',0.113406),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOBr_G4    |         114
CHOFClBr_G4 |         38
CHOFBr_G4   |         148
CHOClBr_G4  |         76
""",
)

entry(
    index = 864,
    label = "CsCClHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   H  u0 {1,S}
4   H  u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.2051,10.9015,12.3424,13.5855,15.4688,16.8304,18.5266],'cal/(mol*K)','+|-',[0.0691911,0.0795211,0.0812901,0.0800975,0.0694553,0.0599046,0.0488209]),
        H298 = (-14.2491,'kcal/mol','+|-',0.246661),
        S298 = (36.8183,'cal/(mol*K)','+|-',0.161907),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOCl_G4    |         72
CHOFCl_G4   |         45
CHOFClBr_G4 |         2
CHOClBr_G4  |         39
""",
)

entry(
    index = 865,
    label = "CsClCsHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cs u0 {1,S}
3   H  u0 {1,S}
4   H  u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.73272,11.5646,13.0312,14.2592,16.118,17.3667,19.0044],'cal/(mol*K)','+|-',[0.0497673,0.0571974,0.0584697,0.0576119,0.0499573,0.0430878,0.0351156]),
        H298 = (-16.8453,'kcal/mol','+|-',0.177416),
        S298 = (35.5899,'cal/(mol*K)','+|-',0.116455),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library                  | Number of Species
CHOCl_G4                 |         157
Chlorinated_Hydrocarbons |         2
CHOFCl_G4                |         74
CHOFClBr_G4              |         4
CHOClBr_G4               |         108
""",
)

entry(
    index = 866,
    label = "CsCFHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   H  u0 {1,S}
4   H  u0 {1,S}
5   F  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOF_G4     |         72
CHOFCl_G4   |         36
CHOFClBr_G4 |         13
CHOFBr_G4   |         64
""",
)

entry(
    index = 867,
    label = "CsCsFHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cs u0 {1,S}
3   H  u0 {1,S}
4   H  u0 {1,S}
5   F  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOF_G4     |         163
CHOFCl_G4   |         103
CHOFClBr_G4 |         49
CHOFBr_G4   |         195
""",
)

entry(
    index = 868,
    label = "CsBrBrCO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   O  u0 {1,S}
4   Br u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([14.5433,16.1629,17.0347,17.6132,18.2123,18.3306,17.6446],'cal/(mol*K)','+|-',[0.0864582,0.0993661,0.101576,0.100086,0.0867883,0.0748541,0.0610044]),
        H298 = (1.69082,'kcal/mol','+|-',0.308216),
        S298 = (24.6289,'cal/(mol*K)','+|-',0.202312),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOBr_G4    |         39
CHOFClBr_G4 |         1
CHOFBr_G4   |         50
CHOClBr_G4  |         16
""",
)

entry(
    index = 869,
    label = "CsBrCClO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   O  u0 {1,S}
4   Cl u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([13.9494,15.837,16.8196,17.4768,18.1622,18.2947,17.5929],'cal/(mol*K)','+|-',[0.120931,0.138986,0.142077,0.139993,0.121393,0.1047,0.0853284]),
        H298 = (-10.9098,'kcal/mol','+|-',0.431109),
        S298 = (21.1226,'cal/(mol*K)','+|-',0.282978),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOFClBr_G4 |         16
CHOClBr_G4  |         33
""",
)

entry(
    index = 870,
    label = "CsCClClO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   O  u0 {1,S}
4   Cl u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([13.389,15.3514,16.3795,17.07,17.8052,17.9899,17.3758],'cal/(mol*K)','+|-',[0.0843206,0.0969094,0.0990652,0.0976118,0.0846426,0.0730035,0.0594962]),
        H298 = (-22.8011,'kcal/mol','+|-',0.300596),
        S298 = (19.5947,'cal/(mol*K)','+|-',0.19731),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOCl_G4    |         86
CHOFCl_G4   |         14
CHOFClBr_G4 |         2
CHOClBr_G4  |         12
""",
)

entry(
    index = 871,
    label = "CsBrCFO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   O  u0 {1,S}
4   F  u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([12.904,14.8838,16.0095,16.7781,17.6155,17.8602,17.2953],'cal/(mol*K)','+|-',[0.112322,0.129092,0.131963,0.130027,0.112751,0.0972469,0.079254]),
        H298 = (-55.5978,'kcal/mol','+|-',0.40042),
        S298 = (19.4721,'cal/(mol*K)','+|-',0.262834),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOFBr_G4 |         58
""",
)

entry(
    index = 872,
    label = "CsCClFO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   O  u0 {1,S}
4   F  u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([12.3519,14.5607,15.7159,16.496,17.3994,17.6879,17.1349],'cal/(mol*K)','+|-',[0.128194,0.147333,0.150611,0.148401,0.128684,0.110989,0.0904532]),
        H298 = (-69.0461,'kcal/mol','+|-',0.457002),
        S298 = (15.9389,'cal/(mol*K)','+|-',0.299974),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOFCl_G4   |         31
CHOFClBr_G4 |         12
""",
)

entry(
    index = 873,
    label = "CsCFFO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   O  u0 {1,S}
4   F  u0 {1,S}
5   F  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOF_G4   |         87
CHOFCl_G4 |         11
CHOFBr_G4 |         30
""",
)

entry(
    index = 874,
    label = "CsBrCHO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   O  u0 {1,S}
4   H  u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([10.4045,12.3718,13.621,14.5309,15.6574,16.1899,16.3217],'cal/(mol*K)','+|-',[0.0615776,0.0707709,0.0723452,0.0712838,0.0618127,0.0533129,0.0434488]),
        H298 = (-5.43336,'kcal/mol','+|-',0.219519),
        S298 = (14.939,'cal/(mol*K)','+|-',0.144091),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOBr_G4    |         62
CHOFClBr_G4 |         14
CHOFBr_G4   |         93
CHOClBr_G4  |         49
""",
)

entry(
    index = 875,
    label = "CsCClHO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   O  u0 {1,S}
4   H  u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.81745,11.9188,13.2374,14.1848,15.3747,15.9694,16.18],'cal/(mol*K)','+|-',[0.0668073,0.0767814,0.0784894,0.0773379,0.0670624,0.0578407,0.0471389]),
        H298 = (-17.7685,'kcal/mol','+|-',0.238163),
        S298 = (12.3477,'cal/(mol*K)','+|-',0.156329),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOCl_G4    |         87
CHOFCl_G4   |         47
CHOFClBr_G4 |         15
CHOClBr_G4  |         29
""",
)

entry(
    index = 876,
    label = "CsCFHO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   O  u0 {1,S}
4   H  u0 {1,S}
5   F  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.22865,10.3914,11.9148,13.0267,14.4607,15.2575,15.841],'cal/(mol*K)','+|-',[0.0681607,0.0783369,0.0800795,0.0789046,0.068421,0.0590125,0.0480939]),
        H298 = (-61.8215,'kcal/mol','+|-',0.242987),
        S298 = (10.1098,'cal/(mol*K)','+|-',0.159496),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOF_G4     |         86
CHOFCl_G4   |         29
CHOFClBr_G4 |         3
CHOFBr_G4   |         49
""",
)

entry(
    index = 877,
    label = "CsBrCOO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   O  u0 {1,S}
4   O  u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([11.187,13.424,14.4479,15.03,15.4427,15.2289,13.8402],'cal/(mol*K)','+|-',[0.1817,0.208828,0.213473,0.210341,0.182394,0.157313,0.128207]),
        H298 = (-11.9517,'kcal/mol','+|-',0.647747),
        S298 = (-8.10363,'cal/(mol*K)','+|-',0.425178),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOBr_G4    |         7
CHOFClBr_G4 |         1
CHOFBr_G4   |         10
CHOClBr_G4  |         4
""",
)

entry(
    index = 878,
    label = "CsCClOO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   O  u0 {1,S}
4   O  u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([10.2021,12.7685,13.9537,14.568,14.9918,14.8264,13.5524],'cal/(mol*K)','+|-',[0.18667,0.214539,0.219311,0.216094,0.187383,0.161616,0.131713]),
        H298 = (-24.8739,'kcal/mol','+|-',0.665463),
        S298 = (-10.4212,'cal/(mol*K)','+|-',0.436807),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOCl_G4    |         12
CHOFCl_G4   |         4
CHOFClBr_G4 |         1
CHOClBr_G4  |         4
""",
)

entry(
    index = 879,
    label = "CsCFOO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   O  u0 {1,S}
4   O  u0 {1,S}
5   F  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOF_G4     |         12
CHOFCl_G4   |         4
CHOFClBr_G4 |         1
CHOFBr_G4   |         7
""",
)

entry(
    index = 880,
    label = "CsBrBrCC",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   C  u0 {1,S}
4   Br u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([13.9848,15.2596,15.9896,16.5405,17.1724,17.4386,17.0386],'cal/(mol*K)','+|-',[0.149337,0.171633,0.175451,0.172877,0.149907,0.129294,0.105372]),
        H298 = (5.51692,'kcal/mol','+|-',0.532375),
        S298 = (26.7511,'cal/(mol*K)','+|-',0.349449),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOBr_G4    |         12
CHOFClBr_G4 |         2
CHOFBr_G4   |         14
CHOClBr_G4  |         9
""",
)

entry(
    index = 881,
    label = "CsBrBrCsCs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cs u0 {1,S}
3   Cs u0 {1,S}
4   Br u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([14.0149,15.386,16.1965,16.8103,17.6137,17.882,17.5188],'cal/(mol*K)','+|-',[0.112863,0.129714,0.132599,0.130654,0.113294,0.0977154,0.0796359]),
        H298 = (2.55179,'kcal/mol','+|-',0.402349),
        S298 = (25.4383,'cal/(mol*K)','+|-',0.2641),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOBr_G4    |         19
CHOFClBr_G4 |         6
CHOFBr_G4   |         29
CHOClBr_G4  |         15
""",
)

entry(
    index = 882,
    label = "CsBrCCCl",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   C  u0 {1,S}
4   Cl u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([13.5696,14.9243,15.699,16.2857,17.0413,17.3248,16.8986],'cal/(mol*K)','+|-',[0.13211,0.151834,0.155211,0.152934,0.132615,0.114379,0.0932163]),
        H298 = (-7.8729,'kcal/mol','+|-',0.470962),
        S298 = (23.0113,'cal/(mol*K)','+|-',0.309138),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOFClBr_G4 |         13
CHOClBr_G4  |         30
""",
)

entry(
    index = 883,
    label = "CsCCClCl",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   C  u0 {1,S}
4   Cl u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([13.152,14.64,15.4351,16.0703,16.8266,17.1288,16.5686],'cal/(mol*K)','+|-',[0.135387,0.155599,0.159061,0.156727,0.135904,0.117216,0.0955281]),
        H298 = (-17.4834,'kcal/mol','+|-',0.482642),
        S298 = (21.6844,'cal/(mol*K)','+|-',0.316804),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOCl_G4    |         41
CHOFCl_G4   |         4
CHOFClBr_G4 |         1
CHOClBr_G4  |         8
""",
)

entry(
    index = 884,
    label = "CsClClCsCs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cs u0 {1,S}
3   Cs u0 {1,S}
4   Cl u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([13.4719,15.0106,15.7943,16.3897,17.1827,17.4287,16.9428],'cal/(mol*K)','+|-',[0.117396,0.134923,0.137924,0.1359,0.117844,0.101639,0.0828339]),
        H298 = (-20.6328,'kcal/mol','+|-',0.418506),
        S298 = (20.4615,'cal/(mol*K)','+|-',0.274706),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOCl_G4    |         73
CHOFCl_G4   |         8
CHOFClBr_G4 |         1
CHOClBr_G4  |         9
""",
)

entry(
    index = 885,
    label = "CsBrCCF",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   C  u0 {1,S}
4   F  u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([13.2262,14.6469,15.3962,15.9577,16.6734,16.9381,16.5545],'cal/(mol*K)','+|-',[0.118265,0.135922,0.138946,0.136907,0.118717,0.102392,0.0834474]),
        H298 = (-49.3413,'kcal/mol','+|-',0.421606),
        S298 = (20.5507,'cal/(mol*K)','+|-',0.276741),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOFBr_G4 |         57
""",
)

entry(
    index = 886,
    label = "CsCCClF",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   C  u0 {1,S}
4   F  u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([12.3841,13.8916,14.6902,15.3272,16.2184,16.6085,16.3437],'cal/(mol*K)','+|-',[0.130721,0.150238,0.15358,0.151326,0.13122,0.113176,0.0922363]),
        H298 = (-61.7874,'kcal/mol','+|-',0.466011),
        S298 = (17.5229,'cal/(mol*K)','+|-',0.305888),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOFCl_G4   |         27
CHOFClBr_G4 |         17
""",
)

entry(
    index = 887,
    label = "CsCCFF",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   C  u0 {1,S}
4   F  u0 {1,S}
5   F  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOF_G4     |         41
CHOFCl_G4   |         3
CHOFClBr_G4 |         2
CHOFBr_G4   |         11
""",
)

entry(
    index = 888,
    label = "CsCsCsFF",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cs u0 {1,S}
3   Cs u0 {1,S}
4   F  u0 {1,S}
5   F  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOF_G4     |         76
CHOFCl_G4   |         9
CHOFClBr_G4 |         1
CHOFBr_G4   |         33
""",
)

entry(
    index = 889,
    label = "CsBrCCH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   C  u0 {1,S}
4   H  u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.96866,11.3913,12.4029,13.2398,14.3431,15.0075,15.4414],'cal/(mol*K)','+|-',[0.117444,0.134978,0.13798,0.135956,0.117892,0.101681,0.0828677]),
        H298 = (-2.24303,'kcal/mol','+|-',0.418677),
        S298 = (17.3826,'cal/(mol*K)','+|-',0.274818),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOBr_G4    |         24
CHOFClBr_G4 |         3
CHOFBr_G4   |         28
CHOClBr_G4  |         14
""",
)

entry(
    index = 890,
    label = "CsBrCsCsH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cs u0 {1,S}
3   Cs u0 {1,S}
4   H  u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.94605,11.396,12.471,13.3333,14.579,15.2575,15.7926],'cal/(mol*K)','+|-',[0.0826673,0.0950093,0.0971228,0.0956979,0.082983,0.0715721,0.0583296]),
        H298 = (-3.83905,'kcal/mol','+|-',0.294702),
        S298 = (16.3063,'cal/(mol*K)','+|-',0.193441),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOBr_G4    |         37
CHOFClBr_G4 |         10
CHOFBr_G4   |         68
CHOClBr_G4  |         29
""",
)

entry(
    index = 891,
    label = "CsCCClH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   C  u0 {1,S}
4   H  u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.16325,10.6765,11.7764,12.6948,13.9744,14.7845,15.3096],'cal/(mol*K)','+|-',[0.109363,0.125691,0.128487,0.126602,0.109781,0.0946852,0.0771663]),
        H298 = (-13.1538,'kcal/mol','+|-',0.389872),
        S298 = (14.4562,'cal/(mol*K)','+|-',0.25591),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOCl_G4    |         41
CHOFCl_G4   |         14
CHOFClBr_G4 |         6
CHOClBr_G4  |         19
""",
)

entry(
    index = 892,
    label = "CsCsCsClH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cs u0 {1,S}
3   Cs u0 {1,S}
4   H  u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.37603,10.8757,11.964,12.8431,14.1535,14.8999,15.5761],'cal/(mol*K)','+|-',[0.0832924,0.0957276,0.0978571,0.0964214,0.0836103,0.0721132,0.0587706]),
        H298 = (-15.2695,'kcal/mol','+|-',0.29693),
        S298 = (13.8611,'cal/(mol*K)','+|-',0.194904),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library                  | Number of Species
CHOCl_G4                 |         74
Chlorinated_Hydrocarbons |         3
CHOFCl_G4                |         29
CHOFClBr_G4              |         12
CHOClBr_G4               |         33
""",
)

entry(
    index = 893,
    label = "CsCCFH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   C  u0 {1,S}
4   H  u0 {1,S}
5   F  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.01653,9.5087,10.6909,11.7025,13.1613,14.1292,14.9033],'cal/(mol*K)','+|-',[0.1039,0.119412,0.122068,0.120277,0.104296,0.0899548,0.0733111]),
        H298 = (-52.0569,'kcal/mol','+|-',0.370394),
        S298 = (12.0517,'cal/(mol*K)','+|-',0.243125),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOF_G4     |         41
CHOFCl_G4   |         11
CHOFClBr_G4 |         9
CHOFBr_G4   |         26
""",
)

entry(
    index = 894,
    label = "CsCsCsFH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cs u0 {1,S}
3   Cs u0 {1,S}
4   H  u0 {1,S}
5   F  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.49154,9.97409,11.0744,11.9677,13.3429,14.2021,15.1414],'cal/(mol*K)','+|-',[0.0757731,0.0870857,0.0890229,0.0877169,0.0760624,0.0656031,0.0534651]),
        H298 = (-55.1509,'kcal/mol','+|-',0.270125),
        S298 = (11.7666,'cal/(mol*K)','+|-',0.177309),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOF_G4     |         76
CHOFCl_G4   |         32
CHOFClBr_G4 |         11
CHOFBr_G4   |         71
""",
)

entry(
    index = 895,
    label = "CsBrCCO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   C  u0 {1,S}
4   O  u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([10.4839,12.3943,13.2831,13.8604,14.4559,14.4541,13.4233],'cal/(mol*K)','+|-',[0.164112,0.188613,0.192809,0.18998,0.164738,0.142085,0.115796]),
        H298 = (-6.75758,'kcal/mol','+|-',0.585045),
        S298 = (-7.44993,'cal/(mol*K)','+|-',0.384021),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOBr_G4    |         8
CHOFClBr_G4 |         1
CHOFBr_G4   |         14
CHOClBr_G4  |         6
""",
)

entry(
    index = 896,
    label = "CsCCClO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   C  u0 {1,S}
4   O  u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.95921,11.9351,12.8264,13.4089,14.0484,14.0917,13.1085],'cal/(mol*K)','+|-',[0.158522,0.182188,0.186241,0.183509,0.159127,0.137246,0.111852]),
        H298 = (-18.6266,'kcal/mol','+|-',0.565117),
        S298 = (-9.52281,'cal/(mol*K)','+|-',0.37094),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOCl_G4    |         19
CHOFCl_G4   |         5
CHOFClBr_G4 |         2
CHOClBr_G4  |         7
""",
)

entry(
    index = 897,
    label = "CsCCFO",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   C  u0 {1,S}
4   O  u0 {1,S}
5   F  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOF_G4     |         20
CHOFCl_G4   |         7
CHOFClBr_G4 |         2
CHOFBr_G4   |         18
""",
)

entry(
    index = 898,
    label = "CsBrCCC",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   C  u0 {1,S}
4   C  u0 {1,S}
5   Br u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.54855,10.52,11.0952,11.5673,12.3491,12.6698,12.4641],'cal/(mol*K)','+|-',[0.185715,0.213442,0.21819,0.214989,0.186424,0.160789,0.13104]),
        H298 = (-2.86688,'kcal/mol','+|-',0.662059),
        S298 = (-4.67826,'cal/(mol*K)','+|-',0.434573),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOBr_G4    |         6
CHOFClBr_G4 |         1
CHOFBr_G4   |         10
CHOClBr_G4  |         7
""",
)

entry(
    index = 899,
    label = "CsCCCCl",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   C  u0 {1,S}
4   C  u0 {1,S}
5   Cl u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.21509,10.4368,11.1091,11.6359,12.4374,12.6709,12.2742],'cal/(mol*K)','+|-',[0.167756,0.192801,0.19709,0.194198,0.168396,0.14524,0.118367]),
        H298 = (-13.7218,'kcal/mol','+|-',0.598034),
        S298 = (-7.03497,'cal/(mol*K)','+|-',0.392547),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library                  | Number of Species
CHOCl_G4                 |         20
Chlorinated_Hydrocarbons |         4
CHOFCl_G4                |         4
CHOFClBr_G4              |         1
CHOClBr_G4               |         5
""",
)

entry(
    index = 900,
    label = "CsCCCF",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   C  u0 {1,S}
4   C  u0 {1,S}
5   F  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.38048,9.51432,10.1781,10.6966,11.4614,11.7498,11.6067],'cal/(mol*K)','+|-',[0.154098,0.177104,0.181044,0.178388,0.154686,0.133415,0.108731]),
        H298 = (-55.7617,'kcal/mol','+|-',0.549345),
        S298 = (-8.39095,'cal/(mol*K)','+|-',0.360588),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOF_G4     |         20
CHOFCl_G4   |         5
CHOFClBr_G4 |         2
CHOFBr_G4   |         15
""",
)

entry(
    index = 901,
    label = "Cs-HHHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   H  u0 {1,S}
3   H  u0 {1,S}
4   H  u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """CHEMKIN DATABASE S(group) = S(CH4) + Rln(12)""",
    longDesc = 
"""

""",
)

entry(
    index = 902,
    label = "Cs-CHHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   H  u0 {1,S}
4   H  u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = 'Cs-CsHHH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 903,
    label = "Cs-CsHHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cs u0 {1,S}
3   H  u0 {1,S}
4   H  u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Cs-CsHHH BENSON""",
    longDesc = 
"""

""",
)

entry(
    index = 904,
    label = "Cs-CdsHHH",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   [Cd,CO] u0 {1,S}
3   H       u0 {1,S}
4   H       u0 {1,S}
5   H       u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)HHH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 905,
    label = "Cs-(Cds-O2d)HHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   H   u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 906,
    label = "Cs-(Cds-Cd)HHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   H  u0 {1,S}
4   H  u0 {1,S}
5   H  u0 {1,S}
6   C  u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)HHH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 907,
    label = "Cs-(Cds-Cds)HHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   H  u0 {1,S}
4   H  u0 {1,S}
5   H  u0 {1,S}
6   Cd u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Cs-CdHHH BENSON (Assigned Cs-CsHHH)""",
    longDesc = 
"""

""",
)

entry(
    index = 908,
    label = "Cs-(Cds-Cdd)HHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   H   u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)HHH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 909,
    label = "Cs-(Cds-Cdd-O2d)HHH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   H   u0 {1,S}
7   O2d u0 {3,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([25.31,32.07,38.44,44.06,53.36,60.63,72.47],'J/(mol*K)'),
        H298 = (-42.9,'kJ/mol'),
        S298 = (127.12,'J/(mol*K)'),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 910,
    label = "Cs-(Cds-Cdd-S2d)HHH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   H   u0 {1,S}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 911,
    label = "Cs-(Cds-Cdd-Cd)HHH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   H   u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)HHH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 912,
    label = "Cs-CtHHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Ct u0 {1,S}
3   H  u0 {1,S}
4   H  u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.19,7.84,9.4,10.79,13.02,14.77,17.58],'cal/(mol*K)','+|-',[0.08,0.08,0.08,0.08,0.08,0.08,0.08]),
        H298 = (-10.2,'kcal/mol','+|-',0.15),
        S298 = (30.41,'cal/(mol*K)','+|-',0.08),
    ),
    shortDesc = """Cs-CtHHH BENSON (Assigned Cs-CsHHH)""",
    longDesc = 
"""

""",
)

entry(
    index = 913,
    label = "Cs-(CtN3t)HHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Ct  u0 {1,S} {6,T}
3   H   u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   N3t u0 {2,T}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([12.2791,14.3896,16.359,18.1113,21.0145,23.2453,26.7496],'cal/(mol*K)','+|-',[0.9602,1.02548,1.05759,1.05008,1.02148,1.00315,1.02359]),
        H298 = (17.323,'kcal/mol','+|-',3.84862),
        S298 = (59.9823,'cal/(mol*K)','+|-',2.94866),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         1
""",
)

entry(
    index = 914,
    label = "Cs-CbHHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cb u0 {1,S}
3   H  u0 {1,S}
4   H  u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.19,7.84,9.4,10.79,13.02,14.77,17.58],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (-10.2,'kcal/mol','+|-',0.18),
        S298 = (30.41,'cal/(mol*K)','+|-',0.14),
    ),
    shortDesc = """Cs-CbHHH BENSON (Assigned Cs-CsHHH)""",
    longDesc = 
"""

""",
)

entry(
    index = 915,
    label = "Cs-C=SHHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   H   u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   S2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.99,7.23,9.48,11.01,13.13,14.7,17.27],'cal/(mol*K)'),
        H298 = (-7.1,'kcal/mol'),
        S298 = (31.12,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 916,
    label = "Cs-OsHHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   O2s u0 {1,S}
3   H   u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 917,
    label = "Cs-OsOsHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   O2s u0 {1,S}
3   O2s u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([27.35,43.68,53.55,58.15,60.86,61.66,63.53],'J/(mol*K)','+|-',[5.77,5.77,5.77,5.77,5.77,5.77,5.77]),
        H298 = (-67.5,'kJ/mol','+|-',4.92),
        S298 = (17.89,'J/(mol*K)','+|-',6.74),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 918,
    label = "Cs-OsOsOsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   O2s u0 {1,S}
3   O2s u0 {1,S}
4   O2s u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.54,6,7.17,8.05,9.31,10.05,10.05],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (-21.23,'kcal/mol','+|-',0.2),
        S298 = (-12.07,'cal/(mol*K)','+|-',0.1),
    ),
    shortDesc = """Cs-OOOH BOZZELLI del C/C2/O - C/C3/O, series !!!WARNING! Cp1500 value taken as Cp1000""",
    longDesc = 
"""

""",
)

entry(
    index = 919,
    label = "Cs-OsSHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   O2s u0 {1,S}
3   S   u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = 'Cs-OsS2HH',
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2018""",
    longDesc = 
"""
"
""",
)

entry(
    index = 920,
    label = "Cs-OsS2HH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   O2s u0 {1,S}
3   S2s u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.61,10.2,11.36,12.27,14.72,16.15,15.88],'cal/(mol*K)'),
        H298 = (3.32,'kcal/mol'),
        S298 = (11.26,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 921,
    label = "Cs-OsS4HH",
    group = 
"""
1 * Cs                u0 {2,S} {3,S} {4,S} {5,S}
2   O2s               u0 {1,S}
3   [S4s,S4d,S4b,S4t] u0 {1,S}
4   H                 u0 {1,S}
5   H                 u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.97,7.92,9.33,10.67,14.39,16.61,16.95],'cal/(mol*K)'),
        H298 = (5.62,'kcal/mol'),
        S298 = (8.3,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 922,
    label = "Cs-OsSSH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   O2s u0 {1,S}
3   S   u0 {1,S}
4   S   u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = 'Cs-OsS2S2H',
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2018""",
    longDesc = 
"""
"
""",
)

entry(
    index = 923,
    label = "Cs-OsS2S2H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   O2s u0 {1,S}
3   S2s u0 {1,S}
4   S2s u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.48,8.15,9.41,10.4,11.71,12.83,14.46],'cal/(mol*K)'),
        H298 = (-10.34,'kcal/mol'),
        S298 = (6.83,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2018""",
    longDesc = 
"""
"
""",
)

entry(
    index = 924,
    label = "Cs-OsS4S2H",
    group = 
"""
1 * Cs                u0 {2,S} {3,S} {4,S} {5,S}
2   O2s               u0 {1,S}
3   S2s               u0 {1,S}
4   [S4s,S4d,S4b,S4t] u0 {1,S}
5   H                 u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.61,9,10.51,11.92,16.57,18.75,16.38],'cal/(mol*K)'),
        H298 = (19.45,'kcal/mol'),
        S298 = (-12.67,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 925,
    label = "Cs-OsOsSH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   O2s u0 {1,S}
3   O2s u0 {1,S}
4   S   u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = 'Cs-OsOsS2H',
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2018""",
    longDesc = 
"""
"
""",
)

entry(
    index = 926,
    label = "Cs-OsOsS2H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   O2s u0 {1,S}
3   O2s u0 {1,S}
4   S2s u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.28,9.26,10.65,11.64,14.06,15.03,13.66],'cal/(mol*K)'),
        H298 = (-3.58,'kcal/mol'),
        S298 = (-10.1,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 927,
    label = "Cs-OsOsS4H",
    group = 
"""
1 * Cs                u0 {2,S} {3,S} {4,S} {5,S}
2   O2s               u0 {1,S}
3   O2s               u0 {1,S}
4   [S4s,S4d,S4b,S4t] u0 {1,S}
5   H                 u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.65,7.55,8.69,9.85,13.15,14.98,14.49],'cal/(mol*K)'),
        H298 = (-3.81,'kcal/mol'),
        S298 = (-12.44,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 928,
    label = "Cs-SsHHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   S  u0 {1,S}
3   H  u0 {1,S}
4   H  u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = 'Cs-S2sHHH',
    shortDesc = """CBS-QB3 GA 1D-HR Aaron Vandeputte 2010""",
    longDesc = 
"""

""",
)

entry(
    index = 929,
    label = "Cs-S2sHHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   S2s u0 {1,S}
3   H   u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.68,8.89,10.23,11.52,14.91,17.19,18.48],'cal/(mol*K)'),
        H298 = (1.96,'kcal/mol'),
        S298 = (35.84,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 930,
    label = "Cs-S4HHH",
    group = 
"""
1 * Cs                u0 {2,S} {3,S} {4,S} {5,S}
2   [S4s,S4d,S4b,S4t] u0 {1,S}
3   H                 u0 {1,S}
4   H                 u0 {1,S}
5   H                 u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.48,6.85,8.1,9.45,13.61,16.46,18.29],'cal/(mol*K)'),
        H298 = (5.53,'kcal/mol'),
        S298 = (33.83,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 931,
    label = "Cs-S6HHH",
    group = 
"""
1 * Cs                      u0 {2,S} {3,S} {4,S} {5,S}
2   [S6s,S6d,S6dd,S6t,S6td] u0 {1,S}
3   H                       u0 {1,S}
4   H                       u0 {1,S}
5   H                       u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.93,7.19,8.82,10.34,13.83,16.2,18.09],'cal/(mol*K)'),
        H298 = (-0.84,'kcal/mol'),
        S298 = (41.29,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 932,
    label = "Cs-SsSsHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   S  u0 {1,S}
3   S  u0 {1,S}
4   H  u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.04,11.14,12.08,12.95,16.92,19,17.34],'cal/(mol*K)'),
        H298 = (19.1,'kcal/mol'),
        S298 = (16.46,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 933,
    label = "Cs-SsSsSsH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   S  u0 {1,S}
3   S  u0 {1,S}
4   S  u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.14,12.16,12.79,13.27,17.74,19.66,15.75],'cal/(mol*K)'),
        H298 = (36.14,'kcal/mol'),
        S298 = (-0.63,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 934,
    label = "Cs-CCHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   C  u0 {1,S}
4   H  u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = 'Cs-CsCsHH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 935,
    label = "Cs-CsCsHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cs u0 {1,S}
3   Cs u0 {1,S}
4   H  u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Cs-CsCsHH BENSON""",
    longDesc = 
"""

""",
)

entry(
    index = 936,
    label = "Cs-CdsCsHH",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   [Cd,CO] u0 {1,S}
3   Cs      u0 {1,S}
4   H       u0 {1,S}
5   H       u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CsHH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 937,
    label = "Cs-(Cds-O2d)CsHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([26.91,30.8,34.98,38.91,45.56,50.73,58.93],'J/(mol*K)','+|-',[1.53,1.53,1.53,1.53,1.53,1.53,1.53]),
        H298 = (-21.5,'kJ/mol','+|-',1.3),
        S298 = (40.32,'J/(mol*K)','+|-',1.78),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 938,
    label = "Cs-(Cds-Cd)CsHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cs u0 {1,S}
4   H  u0 {1,S}
5   H  u0 {1,S}
6   C  u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)CsHH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 939,
    label = "Cs-(Cds-Cds)CsHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cs u0 {1,S}
4   H  u0 {1,S}
5   H  u0 {1,S}
6   Cd u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.12,6.86,8.32,9.49,11.22,12.48,14.36],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (-4.76,'kcal/mol','+|-',0.16),
        S298 = (9.8,'cal/(mol*K)','+|-',0.1),
    ),
    shortDesc = """Cs-CdCsHH BENSON""",
    longDesc = 
"""

""",
)

entry(
    index = 940,
    label = "Cs-(Cds-Cdd)CsHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)CsHH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 941,
    label = "Cs-(Cds-Cdd-O2d)CsHH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cs  u0 {1,S}
5   H   u0 {1,S}
6   H   u0 {1,S}
7   O2d u0 {3,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.35,6.83,8.25,9.45,11.19,12.46,14.34],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (-5.723,'kcal/mol','+|-',0.16),
        S298 = (9.37,'cal/(mol*K)','+|-',0.1),
    ),
    shortDesc = """{C/C/H2/CCO} RAMAN & GREEN JPCA 2002, 106, 7937-7949""",
    longDesc = 
"""

""",
)

entry(
    index = 942,
    label = "Cs-(Cds-Cdd-S2d)CsHH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cs  u0 {1,S}
5   H   u0 {1,S}
6   H   u0 {1,S}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 943,
    label = "Cs-(Cds-Cdd-Cd)CsHH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cs  u0 {1,S}
5   H   u0 {1,S}
6   H   u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)CsHH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 944,
    label = "Cs-CdsCdsHH",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   [Cd,CO] u0 {1,S}
3   [Cd,CO] u0 {1,S}
4   H       u0 {1,S}
5   H       u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)HH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 945,
    label = "Cs-(Cds-O2d)(Cds-O2d)HH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   CO  u0 {1,S} {7,D}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   O2d u0 {2,D}
7   O2d u0 {3,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([27.77,32.81,37.1,40.67,46.39,50.85,58.25],'J/(mol*K)','+|-',[4.19,4.19,4.19,4.19,4.19,4.19,4.19]),
        H298 = (-10,'kJ/mol','+|-',3.57),
        S298 = (40.1,'J/(mol*K)','+|-',4.88),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 946,
    label = "Cs-(Cds-O2d)(Cds-Cd)HH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   C   u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([24.94,31.41,36.47,40.49,46.72,51.49,59.29],'J/(mol*K)','+|-',[3.34,3.34,3.34,3.34,3.34,3.34,3.34]),
        H298 = (-16.9,'kJ/mol','+|-',2.85),
        S298 = (40.18,'J/(mol*K)','+|-',3.9),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 947,
    label = "Cs-(Cds-O2d)(Cds-Cds)HH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   Cd  u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.75,7.11,8.92,10.32,12.16,13.61,13.61],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (-3.8,'kcal/mol','+|-',0.16),
        S298 = (6.31,'cal/(mol*K)','+|-',0.1),
    ),
    shortDesc = """Cs-COCdHH BENSON Hf, Mopac =3D S,Cp nov99 !!!WARNING! Cp1500 value taken as Cp1000""",
    longDesc = 
"""

""",
)

entry(
    index = 948,
    label = "Cs-(Cds-O2d)(Cds-Cdd)HH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-Cd)HH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 949,
    label = "Cs-(Cds-O2d)(Cds-Cdd-O2d)HH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CO  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   H   u0 {1,S}
6   H   u0 {1,S}
7   O2d u0 {3,D}
8   O2d u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)CsHH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 950,
    label = "Cs-(Cds-O2d)(Cds-Cdd-Cd)HH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CO  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   H   u0 {1,S}
6   H   u0 {1,S}
7   O2d u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)HH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 951,
    label = "Cs-(Cds-Cd)(Cds-Cd)HH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cd u0 {1,S} {7,D}
4   H  u0 {1,S}
5   H  u0 {1,S}
6   C  u0 {2,D}
7   C  u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)HH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 952,
    label = "Cs-(Cds-Cds)(Cds-Cds)HH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cd u0 {1,S} {7,D}
4   H  u0 {1,S}
5   H  u0 {1,S}
6   Cd u0 {2,D}
7   Cd u0 {3,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.7,6.8,8.4,9.6,11.3,12.6,14.4],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (-4.29,'kcal/mol','+|-',0.16),
        S298 = (10.2,'cal/(mol*K)','+|-',0.1),
    ),
    shortDesc = """Cs-CdCdHH BENSON""",
    longDesc = 
"""

""",
)

entry(
    index = 953,
    label = "Cs-(Cds-Cdd)(Cds-Cds)HH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {2,D}
7   Cd  u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cds)HH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 954,
    label = "Cs-Cd(CCO)HH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   H   u0 {1,S}
6   H   u0 {1,S}
7   Cd  u0 {3,D}
8   O2d u0 {4,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([25.85,31.99,37.06,41.14,47.42,52.15,59.73],'J/(mol*K)','+|-',[6.93,6.93,6.93,6.93,6.93,6.93,6.93]),
        H298 = (-22.2,'kJ/mol','+|-',5.9),
        S298 = (37.92,'J/(mol*K)','+|-',8.08),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 955,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cds)HH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   H   u0 {1,S}
6   H   u0 {1,S}
7   Cd  u0 {3,D}
8   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 956,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cds)HH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   H   u0 {1,S}
6   H   u0 {1,S}
7   Cd  u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)HH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 957,
    label = "Cs-(Cds-Cdd)(Cds-Cdd)HH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {2,D}
7   Cdd u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)HH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 958,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)HH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   H   u0 {1,S}
7   H   u0 {1,S}
8   O2d u0 {4,D}
9   O2d u0 {5,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.68,8.28,9.58,10.61,12.04,13.13,14.87],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (-5.301,'kcal/mol','+|-',0.16),
        S298 = (7.18,'cal/(mol*K)','+|-',0.1),
    ),
    shortDesc = """{C/H2/CCO2} RAMAN & GREEN JPCA 2002, 106, 7937-7949""",
    longDesc = 
"""

""",
)

entry(
    index = 959,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-Cd)HH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   H   u0 {1,S}
7   H   u0 {1,S}
8   O2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-Cd(CCO)HH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 960,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-S2d)HH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   H   u0 {1,S}
7   H   u0 {1,S}
8   S2d u0 {4,D}
9   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 961,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-Cd)HH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   H   u0 {1,S}
7   H   u0 {1,S}
8   S2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 962,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)HH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   H   u0 {1,S}
7   H   u0 {1,S}
8   C   u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)HH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 963,
    label = "Cs-CtCsHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Ct u0 {1,S}
3   Cs u0 {1,S}
4   H  u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.95,6.56,7.93,9.08,10.86,12.19,14.2],'cal/(mol*K)','+|-',[0.08,0.08,0.08,0.08,0.08,0.08,0.08]),
        H298 = (-4.73,'kcal/mol','+|-',0.28),
        S298 = (10.3,'cal/(mol*K)','+|-',0.07),
    ),
    shortDesc = """Cs-CtCsHH BENSON""",
    longDesc = 
"""

""",
)

entry(
    index = 964,
    label = "Cs-CtCdsHH",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Ct      u0 {1,S}
3   [Cd,CO] u0 {1,S}
4   H       u0 {1,S}
5   H       u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CtHH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 965,
    label = "Cs-(Cds-O2d)CtHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.85,6.22,8.01,9.43,11.29,12.76,12.76],'cal/(mol*K)','+|-',[0.08,0.08,0.08,0.08,0.08,0.08,0.08]),
        H298 = (-5.4,'kcal/mol','+|-',0.28),
        S298 = (7.68,'cal/(mol*K)','+|-',0.07),
    ),
    shortDesc = """Cs-COCtHH BENSON Hf, Mopac =3D S,Cp nov99 !!!WARNING! Cp1500 value taken as Cp1000""",
    longDesc = 
"""

""",
)

entry(
    index = 966,
    label = "Cs-(Cds-Cd)CtHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Ct u0 {1,S}
4   H  u0 {1,S}
5   H  u0 {1,S}
6   C  u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)CtHH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 967,
    label = "Cs-(Cds-Cds)CtHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Ct u0 {1,S}
4   H  u0 {1,S}
5   H  u0 {1,S}
6   Cd u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.4,6.33,7.9,9.16,10.93,12.29,13.43],'cal/(mol*K)','+|-',[0.08,0.08,0.08,0.08,0.08,0.08,0.08]),
        H298 = (-3.49,'kcal/mol','+|-',0.28),
        S298 = (9.31,'cal/(mol*K)','+|-',0.07),
    ),
    shortDesc = """Cs-CtCdHH RAMAN & GREEN JPCA 2002, 106, 11141-11149""",
    longDesc = 
"""

""",
)

entry(
    index = 968,
    label = "Cs-(Cds-Cdd)CtHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)CtHH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 969,
    label = "Cs-(Cds-Cdd-O2d)CtHH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Ct  u0 {1,S}
5   H   u0 {1,S}
6   H   u0 {1,S}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-Cd(CCO)HH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 970,
    label = "Cs-(Cds-Cdd-S2d)CtHH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Ct  u0 {1,S}
5   H   u0 {1,S}
6   H   u0 {1,S}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 971,
    label = "Cs-(Cds-Cdd-Cd)CtHH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Ct  u0 {1,S}
5   H   u0 {1,S}
6   H   u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)CtHH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 972,
    label = "Cs-CtCtHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Ct u0 {1,S}
3   Ct u0 {1,S}
4   H  u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4,6.07,7.71,9.03,10.88,12.3,12.48],'cal/(mol*K)','+|-',[0.08,0.08,0.08,0.08,0.08,0.08,0.08]),
        H298 = (-0.82,'kcal/mol','+|-',0.28),
        S298 = (10.04,'cal/(mol*K)','+|-',0.07),
    ),
    shortDesc = """Cs-CtCtHH RAMAN & GREEN JPCA 2002, 106, 11141-11149""",
    longDesc = 
"""

""",
)

entry(
    index = 973,
    label = "Cs-CbCsHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cb u0 {1,S}
3   Cs u0 {1,S}
4   H  u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.84,7.61,8.98,10.01,11.49,12.54,13.76],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (-4.86,'kcal/mol','+|-',0.2),
        S298 = (9.34,'cal/(mol*K)','+|-',0.19),
    ),
    shortDesc = """Cs-CbCsHH BENSON""",
    longDesc = 
"""

""",
)

entry(
    index = 974,
    label = "Cs-CbCdsHH",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Cb      u0 {1,S}
3   [Cd,CO] u0 {1,S}
4   H       u0 {1,S}
5   H       u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CbHH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 975,
    label = "Cs-(Cds-O2d)CbHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.38,7.59,9.25,10.51,12.19,13.52,13.52],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (-5.4,'kcal/mol','+|-',0.2),
        S298 = (5.89,'cal/(mol*K)','+|-',0.19),
    ),
    shortDesc = """Cs-COCbHH BENSON Hf, Mopac =3D S,Cp nov99 !!!WARNING! Cp1500 value taken as Cp1000""",
    longDesc = 
"""

""",
)

entry(
    index = 976,
    label = "Cs-(Cds-Cd)CbHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cb u0 {1,S}
4   H  u0 {1,S}
5   H  u0 {1,S}
6   C  u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)CbHH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 977,
    label = "Cs-(Cds-Cds)CbHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cb u0 {1,S}
4   H  u0 {1,S}
5   H  u0 {1,S}
6   Cd u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.51,6.76,8.61,10.01,11.97,13.4,15.47],'cal/(mol*K)','+|-',[0.2,0.2,0.2,0.2,0.2,0.2,0.2]),
        H298 = (-4.29,'kcal/mol','+|-',0.2),
        S298 = (2,'cal/(mol*K)','+|-',0.19),
    ),
    shortDesc = """Cs-CbCdHH Hf=Stein S,Cp=3D mopac nov99""",
    longDesc = 
"""

""",
)

entry(
    index = 978,
    label = "Cs-(Cds-Cdd)CbHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)CbHH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 979,
    label = "Cs-(Cds-Cdd-O2d)CbHH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   H   u0 {1,S}
6   H   u0 {1,S}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-Cd(CCO)HH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 980,
    label = "Cs-(Cds-Cdd-S2d)CbHH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   H   u0 {1,S}
6   H   u0 {1,S}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 981,
    label = "Cs-(Cds-Cdd-Cd)CbHH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   H   u0 {1,S}
6   H   u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)CbHH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 982,
    label = "Cs-CbCtHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cb u0 {1,S}
3   Ct u0 {1,S}
4   H  u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.28,6.43,8.16,9.5,11.36,12.74,13.7],'cal/(mol*K)','+|-',[0.08,0.08,0.08,0.08,0.08,0.08,0.08]),
        H298 = (-4.29,'kcal/mol','+|-',0.28),
        S298 = (9.84,'cal/(mol*K)','+|-',0.07),
    ),
    shortDesc = """Cs-CbCtHH Hf=Stein S,Cp=3D mopac nov99""",
    longDesc = 
"""

""",
)

entry(
    index = 983,
    label = "Cs-CbCbHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cb u0 {1,S}
3   Cb u0 {1,S}
4   H  u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.67,7.7,9.31,10.52,12.21,13.47,15.11],'cal/(mol*K)','+|-',[0.2,0.2,0.2,0.2,0.2,0.2,0.2]),
        H298 = (-4.29,'kcal/mol','+|-',0.2),
        S298 = (8.07,'cal/(mol*K)','+|-',0.19),
    ),
    shortDesc = """Cs-CbCbHH Hf=3Dbsn/Cs/Cd2/H2 S,Cp=3D mopac nov99""",
    longDesc = 
"""

""",
)

entry(
    index = 984,
    label = "Cs-C=SCtHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 985,
    label = "Cs-C=SCsHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   S2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.04,6.19,8.2,9.42,10.99,12.05,13.72],'cal/(mol*K)'),
        H298 = (-1.72,'kcal/mol'),
        S298 = (10.53,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 986,
    label = "Cs-C=S(Cds-Cd)HH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   C   u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 987,
    label = "Cs-C=S(Cds-Cdd)HH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 988,
    label = "Cs-C=S(Cds-Cdd-Cd)HH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CS  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   H   u0 {1,S}
6   H   u0 {1,S}
7   S2d u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 989,
    label = "Cs-C=S(Cds-Cdd-S2d)HH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CS  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   H   u0 {1,S}
6   H   u0 {1,S}
7   S2d u0 {3,D}
8   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 990,
    label = "Cs-C=S(Cds-Cds)HH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   Cd  u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 991,
    label = "Cs-C=SC=SHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   CS  u0 {1,S} {7,D}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   S2d u0 {2,D}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 992,
    label = "Cs-C=SCbHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 993,
    label = "Cs-CCCH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   C  u0 {1,S}
4   C  u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = 'Cs-CsCsCsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 994,
    label = "Cs-CsCsCsH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cs u0 {1,S}
3   Cs u0 {1,S}
4   Cs u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.54,6,7.17,8.05,9.31,10.05,11.17],'cal/(mol*K)','+|-',[0.07,0.07,0.07,0.07,0.07,0.07,0.07]),
        H298 = (-1.9,'kcal/mol','+|-',0.15),
        S298 = (-12.07,'cal/(mol*K)','+|-',0.07),
    ),
    shortDesc = """Cs-CsCsCsH BENSON""",
    longDesc = 
"""

""",
)

entry(
    index = 995,
    label = "Cs-CdsCsCsH",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   [Cd,CO] u0 {1,S}
3   Cs      u0 {1,S}
4   Cs      u0 {1,S}
5   H       u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CsCsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 996,
    label = "Cs-(Cds-O2d)CsCsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   Cs  u0 {1,S}
5   H   u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([23.68,27.86,31.26,34,38.07,41,45.46],'J/(mol*K)','+|-',[3.34,3.34,3.34,3.34,3.34,3.34,3.34]),
        H298 = (-5.4,'kJ/mol','+|-',2.85),
        S298 = (-47.41,'J/(mol*K)','+|-',3.9),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 997,
    label = "Cs-(Cds-Cd)CsCsH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cs u0 {1,S}
4   Cs u0 {1,S}
5   H  u0 {1,S}
6   C  u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)CsCsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 998,
    label = "Cs-(Cds-Cds)CsCsH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cs u0 {1,S}
4   Cs u0 {1,S}
5   H  u0 {1,S}
6   Cd u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.16,5.91,7.34,8.19,9.46,10.19,11.28],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (-1.48,'kcal/mol','+|-',0.27),
        S298 = (-11.69,'cal/(mol*K)','+|-',0.15),
    ),
    shortDesc = """Cs-CdCsCsH BENSON""",
    longDesc = 
"""

""",
)

entry(
    index = 999,
    label = "Cs-(Cds-Cdd)CsCsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   Cs  u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)CsCsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1000,
    label = "Cs-(Cds-Cdd-O2d)CsCsH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cs  u0 {1,S}
5   Cs  u0 {1,S}
6   H   u0 {1,S}
7   O2d u0 {3,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([21.23,27.55,32.36,35.85,40.37,43.16,46.94],'J/(mol*K)','+|-',[6.93,6.93,6.93,6.93,6.93,6.93,6.93]),
        H298 = (-11.1,'kJ/mol','+|-',5.9),
        S298 = (-47.59,'J/(mol*K)','+|-',8.08),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 1001,
    label = "Cs-(Cds-Cdd-S2d)CsCsH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cs  u0 {1,S}
5   Cs  u0 {1,S}
6   H   u0 {1,S}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1002,
    label = "Cs-(Cds-Cdd-Cd)CsCsH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cs  u0 {1,S}
5   Cs  u0 {1,S}
6   H   u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.16,5.91,7.34,8.19,9.46,10.19,11.28],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (-1.48,'kcal/mol','+|-',0.27),
        S298 = (-11.69,'cal/(mol*K)','+|-',0.15),
    ),
    shortDesc = """Cs-CdCsCsH BENSON""",
    longDesc = 
"""

""",
)

entry(
    index = 1003,
    label = "Cs-CtCsCsH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Ct u0 {1,S}
3   Cs u0 {1,S}
4   Cs u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.99,5.61,6.85,7.78,9.1,9.9,11.12],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (-1.72,'kcal/mol','+|-',0.27),
        S298 = (-11.19,'cal/(mol*K)','+|-',0.15),
    ),
    shortDesc = """Cs-CtCsCsH BENSON""",
    longDesc = 
"""

""",
)

entry(
    index = 1004,
    label = "Cs-CbCsCsH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cb u0 {1,S}
3   Cs u0 {1,S}
4   Cs u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.88,6.66,7.9,8.75,9.73,10.25,10.68],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (-0.98,'kcal/mol','+|-',0.27),
        S298 = (-12.15,'cal/(mol*K)','+|-',0.15),
    ),
    shortDesc = """Cs-CbCsCsH BENSON""",
    longDesc = 
"""

""",
)

entry(
    index = 1005,
    label = "Cs-CdsCdsCsH",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   [Cd,CO] u0 {1,S}
3   [Cd,CO] u0 {1,S}
4   Cs      u0 {1,S}
5   H       u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1006,
    label = "Cs-(Cds-O2d)(Cds-O2d)CsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   CO  u0 {1,S} {7,D}
4   Cs  u0 {1,S}
5   H   u0 {1,S}
6   O2d u0 {2,D}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-CsCsCsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1007,
    label = "Cs-(Cds-O2d)(Cds-Cd)CsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cs  u0 {1,S}
5   H   u0 {1,S}
6   C   u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([29.32,32.99,35.49,37.28,39.75,41.6,44.96],'J/(mol*K)','+|-',[3.34,3.34,3.34,3.34,3.34,3.34,3.34]),
        H298 = (-2.2,'kJ/mol','+|-',2.85),
        S298 = (-50.47,'J/(mol*K)','+|-',3.9),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 1008,
    label = "Cs-(Cds-O2d)(Cds-Cds)CsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cs  u0 {1,S}
5   H   u0 {1,S}
6   Cd  u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)CsCsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1009,
    label = "Cs-(Cds-O2d)(Cds-Cdd)CsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cs  u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-Cd)CsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1010,
    label = "Cs-(Cds-O2d)(Cds-Cdd-O2d)CsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CO  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cs  u0 {1,S}
6   H   u0 {1,S}
7   O2d u0 {3,D}
8   O2d u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)CsCsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1011,
    label = "Cs-(Cds-O2d)(Cds-Cdd-Cd)CsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CO  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cs  u0 {1,S}
6   H   u0 {1,S}
7   O2d u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)CsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1012,
    label = "Cs-(Cds-Cd)(Cds-Cd)CsH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cd u0 {1,S} {7,D}
4   Cs u0 {1,S}
5   H  u0 {1,S}
6   C  u0 {2,D}
7   C  u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1013,
    label = "Cs-(Cds-Cds)(Cds-Cds)CsH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cd u0 {1,S} {7,D}
4   Cs u0 {1,S}
5   H  u0 {1,S}
6   Cd u0 {2,D}
7   Cd u0 {3,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.28,6.54,7.67,8.48,9.45,10.18,11.24],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (-1.1,'kcal/mol','+|-',0.27),
        S298 = (-13.03,'cal/(mol*K)','+|-',0.15),
    ),
    shortDesc = """Cs-CdCdCsH RAMAN & GREEN JPCA 2002, 106, 11141-11149""",
    longDesc = 
"""

""",
)

entry(
    index = 1014,
    label = "Cs-(Cds-Cdd)(Cds-Cds)CsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cs  u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {2,D}
7   Cd  u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cds)CsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1015,
    label = "Cs-CsCd(CCO)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cs  u0 {1,S}
6   H   u0 {1,S}
7   Cd  u0 {3,D}
8   O2d u0 {4,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([24.45,31.59,36.01,38.8,42.13,44.21,47.25],'J/(mol*K)','+|-',[6.93,6.93,6.93,6.93,6.93,6.93,6.93]),
        H298 = (-10.4,'kJ/mol','+|-',5.9),
        S298 = (-54.03,'J/(mol*K)','+|-',8.08),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 1016,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cds)CsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cs  u0 {1,S}
6   H   u0 {1,S}
7   Cd  u0 {3,D}
8   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1017,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cds)CsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cs  u0 {1,S}
6   H   u0 {1,S}
7   Cd  u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1018,
    label = "Cs-(Cds-Cdd)(Cds-Cdd)CsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cs  u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {2,D}
7   Cdd u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1019,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)CsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cs  u0 {1,S}
7   H   u0 {1,S}
8   O2d u0 {4,D}
9   O2d u0 {5,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([21.19,28,33.91,38.75,46.07,51.36,59.45],'J/(mol*K)','+|-',[3.46,3.46,3.46,3.46,3.46,3.46,3.46]),
        H298 = (-21.1,'kJ/mol','+|-',2.95),
        S298 = (40.95,'J/(mol*K)','+|-',4.04),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 1020,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-Cd)CsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cs  u0 {1,S}
7   H   u0 {1,S}
8   O2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-CsCd(CCO)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1021,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-S2d)CsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cs  u0 {1,S}
7   H   u0 {1,S}
8   S2d u0 {4,D}
9   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1022,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-Cd)CsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cs  u0 {1,S}
7   H   u0 {1,S}
8   S2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1023,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cs  u0 {1,S}
7   H   u0 {1,S}
8   C   u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1024,
    label = "Cs-CtCdsCsH",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Ct      u0 {1,S}
3   [Cd,CO] u0 {1,S}
4   Cs      u0 {1,S}
5   H       u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CtCsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1025,
    label = "Cs-(Cds-O2d)CtCsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   Cs  u0 {1,S}
5   H   u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)CsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1026,
    label = "Cs-(Cds-Cd)CtCsH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Ct u0 {1,S}
4   Cs u0 {1,S}
5   H  u0 {1,S}
6   C  u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)CtCsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1027,
    label = "Cs-(Cds-Cds)CtCsH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Ct u0 {1,S}
4   Cs u0 {1,S}
5   H  u0 {1,S}
6   Cd u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.55,7.21,8.39,9.17,10,10.61,10.51],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (-6.9,'kcal/mol','+|-',0.27),
        S298 = (-13.48,'cal/(mol*K)','+|-',0.15),
    ),
    shortDesc = """Cs-CtCdCsH RAMAN & GREEN JPCA 2002, 106, 11141-11149""",
    longDesc = 
"""

""",
)

entry(
    index = 1028,
    label = "Cs-(Cds-Cdd)CtCsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   Cs  u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)CtCsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1029,
    label = "Cs-(Cds-Cdd-O2d)CtCsH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Ct  u0 {1,S}
5   Cs  u0 {1,S}
6   H   u0 {1,S}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-CsCd(CCO)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1030,
    label = "Cs-(Cds-Cdd-S2d)CtCsH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Ct  u0 {1,S}
5   Cs  u0 {1,S}
6   H   u0 {1,S}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1031,
    label = "Cs-(Cds-Cdd-Cd)CtCsH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Ct  u0 {1,S}
5   Cs  u0 {1,S}
6   H   u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)CtCsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1032,
    label = "Cs-CbCdsCsH",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Cb      u0 {1,S}
3   [Cd,CO] u0 {1,S}
4   Cs      u0 {1,S}
5   H       u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CbCsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1033,
    label = "Cs-(Cds-O2d)CbCsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Cs  u0 {1,S}
5   H   u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)CsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1034,
    label = "Cs-(Cds-Cd)CbCsH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cb u0 {1,S}
4   Cs u0 {1,S}
5   H  u0 {1,S}
6   C  u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)CbCsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1035,
    label = "Cs-(Cds-Cds)CbCsH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cb u0 {1,S}
4   Cs u0 {1,S}
5   H  u0 {1,S}
6   Cd u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.5,6.57,8.07,8.89,9.88,10.39,10.79],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (-1.56,'kcal/mol','+|-',0.27),
        S298 = (-11.77,'cal/(mol*K)','+|-',0.15),
    ),
    shortDesc = """Cs-CbCdCsH BOZZELLI =3D Cs/Cs2/Cd/H + (Cs/Cs2/Cb/H - Cs/Cs3/H)""",
    longDesc = 
"""

""",
)

entry(
    index = 1036,
    label = "Cs-(Cds-Cdd)CbCsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Cs  u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)CbCsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1037,
    label = "Cs-(Cds-Cdd-O2d)CbCsH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Cs  u0 {1,S}
6   H   u0 {1,S}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-CsCd(CCO)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1038,
    label = "Cs-(Cds-Cdd-Cd)CbCsH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Cs  u0 {1,S}
6   H   u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)CbCsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1039,
    label = "Cs-CtCtCsH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Ct u0 {1,S}
3   Ct u0 {1,S}
4   Cs u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.27,5.32,6.9,8.03,9.33,10.21,9.38],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (1.72,'kcal/mol','+|-',0.27),
        S298 = (-11.61,'cal/(mol*K)','+|-',0.15),
    ),
    shortDesc = """Cs-CtCtCsH RAMAN & GREEN JPCA 2002, 106, 11141-11149""",
    longDesc = 
"""

""",
)

entry(
    index = 1040,
    label = "Cs-CbCtCsH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cb u0 {1,S}
3   Ct u0 {1,S}
4   Cs u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.33,6.27,7.58,8.48,9.52,10.1,10.63],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (-1.55,'kcal/mol','+|-',0.27),
        S298 = (-11.65,'cal/(mol*K)','+|-',0.15),
    ),
    shortDesc = """Cs-CbCtCsH BOZZELLI =3D Cs/Cs2/Cb/H + (Cs/Cs2/Ct/H - Cs/Cs3/H)""",
    longDesc = 
"""

""",
)

entry(
    index = 1041,
    label = "Cs-CbCbCsH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cb u0 {1,S}
3   Cb u0 {1,S}
4   Cs u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.22,7.32,8.63,8.45,10.15,10.45,10.89],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (-1.06,'kcal/mol','+|-',0.27),
        S298 = (-12.23,'cal/(mol*K)','+|-',0.15),
    ),
    shortDesc = """Cs-CbCbCsCs BOZZELLI =3D Cs/Cs2/Cb/H + (Cs/Cs2/Cb/H - Cs/Cs3/H)""",
    longDesc = 
"""

""",
)

entry(
    index = 1042,
    label = "Cs-CdsCdsCdsH",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   [Cd,CO] u0 {1,S}
3   [Cd,CO] u0 {1,S}
4   [Cd,CO] u0 {1,S}
5   H       u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1043,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-O2d)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   CO  u0 {1,S} {7,D}
4   CO  u0 {1,S} {8,D}
5   H   u0 {1,S}
6   O2d u0 {2,D}
7   O2d u0 {3,D}
8   O2d u0 {4,D}
""",
    thermo = 'Cs-CsCsCsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1044,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-Cd)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   CO  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {6,D}
5   H   u0 {1,S}
6   C   u0 {4,D}
7   O2d u0 {2,D}
8   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1045,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-Cds)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   CO  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {6,D}
5   H   u0 {1,S}
6   Cd  u0 {4,D}
7   O2d u0 {2,D}
8   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)CsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1046,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-Cdd)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   CO  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {6,D}
5   H   u0 {1,S}
6   Cdd u0 {4,D}
7   O2d u0 {2,D}
8   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-Cdd-Cd)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1047,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-Cdd-O2d)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CO  u0 {1,S} {7,D}
4   CO  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   H   u0 {1,S}
7   O2d u0 {3,D}
8   O2d u0 {4,D}
9   O2d u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)CsCsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1048,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-Cdd-Cd)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CO  u0 {1,S} {7,D}
4   CO  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   H   u0 {1,S}
7   O2d u0 {3,D}
8   O2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1049,
    label = "Cs-(Cds-O2d)(Cds-Cd)(Cds-Cd)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   H   u0 {1,S}
6   C   u0 {3,D}
7   C   u0 {4,D}
8   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([29.26,34.41,37.4,39.22,41.43,43.04,46.12],'J/(mol*K)','+|-',[3.34,3.34,3.34,3.34,3.34,3.34,3.34]),
        H298 = (2.9,'kJ/mol','+|-',2.85),
        S298 = (-53.2,'J/(mol*K)','+|-',3.9),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 1050,
    label = "Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   H   u0 {1,S}
6   Cd  u0 {3,D}
7   Cd  u0 {4,D}
8   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)CsCsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1051,
    label = "Cs-(Cds-O2d)(Cds-Cdd)(Cds-Cds)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   H   u0 {1,S}
6   Cdd u0 {3,D}
7   Cd  u0 {4,D}
8   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-Cd)(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1052,
    label = "Cs-(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cds)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CO  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {7,D}
5   Cdd u0 {2,D} {9,D}
6   H   u0 {1,S}
7   Cd  u0 {4,D}
8   O2d u0 {3,D}
9   O2d u0 {5,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-O2d)CsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1053,
    label = "Cs-(Cds-O2d)(Cds-Cdd-Cd)(Cds-Cds)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CO  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {7,D}
5   Cdd u0 {2,D} {9,D}
6   H   u0 {1,S}
7   Cd  u0 {4,D}
8   O2d u0 {3,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1054,
    label = "Cs-(Cds-O2d)(Cds-Cdd)(Cds-Cdd)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   H   u0 {1,S}
6   Cdd u0 {3,D}
7   Cdd u0 {4,D}
8   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-Cd)(Cds-Cdd-Cd)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1055,
    label = "Cs-(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cdd-O2d)H",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    CO  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    H   u0 {1,S}
8    O2d u0 {4,D}
9    O2d u0 {5,D}
10   O2d u0 {6,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)CsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1056,
    label = "Cs-(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cdd-Cd)H",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    CO  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    H   u0 {1,S}
8    O2d u0 {4,D}
9    O2d u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1057,
    label = "Cs-(Cds-O2d)(Cds-Cdd-Cd)(Cds-Cdd-Cd)H",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    CO  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    H   u0 {1,S}
8    O2d u0 {4,D}
9    C   u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1058,
    label = "Cs-(Cds-Cd)(Cds-Cd)(Cds-Cd)H",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cd u0 {1,S} {7,D}
4   Cd u0 {1,S} {8,D}
5   H  u0 {1,S}
6   C  u0 {2,D}
7   C  u0 {3,D}
8   C  u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1059,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cd u0 {1,S} {7,D}
4   Cd u0 {1,S} {8,D}
5   H  u0 {1,S}
6   Cd u0 {2,D}
7   Cd u0 {3,D}
8   Cd u0 {4,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.51,5.96,7.13,7.98,9.06,9.9,11.23],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (0.41,'kcal/mol','+|-',0.27),
        S298 = (-11.82,'cal/(mol*K)','+|-',0.15),
    ),
    shortDesc = """Cs-CdCdCdH RAMAN & GREEN JPC 2002""",
    longDesc = 
"""

""",
)

entry(
    index = 1060,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   H   u0 {1,S}
6   Cd  u0 {2,D}
7   Cd  u0 {3,D}
8   Cdd u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1061,
    label = "Cs-CdCd(CCO)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   H   u0 {1,S}
7   Cd  u0 {3,D}
8   Cd  u0 {4,D}
9   O2d u0 {5,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([27.62,35.4,39.24,41.25,43.4,44.87,47.43],'J/(mol*K)','+|-',[6.93,6.93,6.93,6.93,6.93,6.93,6.93]),
        H298 = (-6.8,'kJ/mol','+|-',5.9),
        S298 = (-55.37,'J/(mol*K)','+|-',8.08),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 1062,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-S2d)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   H   u0 {1,S}
7   Cd  u0 {3,D}
8   Cd  u0 {4,D}
9   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1063,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   H   u0 {1,S}
7   Cd  u0 {3,D}
8   Cd  u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1064,
    label = "Cs-(Cds-Cds)(Cds-Cdd)(Cds-Cdd)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   H   u0 {1,S}
6   Cd  u0 {2,D}
7   Cdd u0 {3,D}
8   Cdd u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1065,
    label = "Cs-(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-O2d)H",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    H   u0 {1,S}
8    Cd  u0 {4,D}
9    O2d u0 {5,D}
10   O2d u0 {6,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)CsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1066,
    label = "Cs-(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-Cd)H",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    H   u0 {1,S}
8    Cd  u0 {4,D}
9    O2d u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = 'Cs-CdCd(CCO)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1067,
    label = "Cs-(Cds-Cds)(Cds-Cdd-S2d)(Cds-Cdd-S2d)H",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    H   u0 {1,S}
8    Cd  u0 {4,D}
9    S2d u0 {5,D}
10   S2d u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1068,
    label = "Cs-(Cds-Cds)(Cds-Cdd-S2d)(Cds-Cdd-Cd)H",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    H   u0 {1,S}
8    Cd  u0 {4,D}
9    S2d u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1069,
    label = "Cs-(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)H",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    H   u0 {1,S}
8    Cd  u0 {4,D}
9    C   u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1070,
    label = "Cs-(Cds-Cdd)(Cds-Cdd)(Cds-Cdd)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   H   u0 {1,S}
6   Cdd u0 {2,D}
7   Cdd u0 {3,D}
8   Cdd u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1071,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)(Cds-Cdd-O2d)H",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    H   u0 {1,S}
9    O2d u0 {5,D}
10   O2d u0 {6,D}
11   O2d u0 {7,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1072,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)(Cds-Cdd-Cd)H",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    H   u0 {1,S}
9    O2d u0 {5,D}
10   O2d u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-O2d)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1073,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-Cd)(Cds-Cdd-Cd)H",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    H   u0 {1,S}
9    O2d u0 {5,D}
10   C   u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = 'Cs-CdCd(CCO)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1074,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-S2d)(Cds-Cdd-S2d)H",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    H   u0 {1,S}
9    S2d u0 {5,D}
10   S2d u0 {6,D}
11   S2d u0 {7,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1075,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-S2d)(Cds-Cdd-Cd)H",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    H   u0 {1,S}
9    S2d u0 {5,D}
10   S2d u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1076,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-Cd)(Cds-Cdd-Cd)H",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    H   u0 {1,S}
9    S2d u0 {5,D}
10   C   u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1077,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)H",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    H   u0 {1,S}
9    C   u0 {5,D}
10   C   u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1078,
    label = "Cs-CtCdsCdsH",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Ct      u0 {1,S}
3   [Cd,CO] u0 {1,S}
4   [Cd,CO] u0 {1,S}
5   H       u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CtH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1079,
    label = "Cs-(Cds-O2d)(Cds-O2d)CtH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   CO  u0 {1,S} {7,D}
4   Ct  u0 {1,S}
5   H   u0 {1,S}
6   O2d u0 {2,D}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1080,
    label = "Cs-(Cds-O2d)(Cds-Cd)CtH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Ct  u0 {1,S}
5   H   u0 {1,S}
6   C   u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)CtH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1081,
    label = "Cs-(Cds-O2d)(Cds-Cds)CtH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Ct  u0 {1,S}
5   H   u0 {1,S}
6   Cd  u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1082,
    label = "Cs-(Cds-O2d)(Cds-Cdd)CtH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Ct  u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-Cd)CtH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1083,
    label = "Cs-(Cds-O2d)(Cds-Cdd-O2d)CtH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CO  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Ct  u0 {1,S}
6   H   u0 {1,S}
7   O2d u0 {3,D}
8   O2d u0 {4,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1084,
    label = "Cs-(Cds-O2d)(Cds-Cdd-Cd)CtH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CO  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Ct  u0 {1,S}
6   H   u0 {1,S}
7   O2d u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)CtH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1085,
    label = "Cs-(Cds-Cd)(Cds-Cd)CtH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cd u0 {1,S} {7,D}
4   Ct u0 {1,S}
5   H  u0 {1,S}
6   C  u0 {2,D}
7   C  u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CtH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1086,
    label = "Cs-(Cds-Cds)(Cds-Cds)CtH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cd u0 {1,S} {7,D}
4   Ct u0 {1,S}
5   H  u0 {1,S}
6   Cd u0 {2,D}
7   Cd u0 {3,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.68,7.85,8.62,9.16,9.81,10.42,10.49],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (1.88,'kcal/mol','+|-',0.27),
        S298 = (-13.75,'cal/(mol*K)','+|-',0.15),
    ),
    shortDesc = """Cs-CtCdCdH RAMAN & GREEN JPCA 2002, 106, 11141-11149""",
    longDesc = 
"""

""",
)

entry(
    index = 1087,
    label = "Cs-(Cds-Cdd)(Cds-Cds)CtH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Ct  u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {2,D}
7   Cd  u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cds)CtH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1088,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cds)CtH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Ct  u0 {1,S}
6   H   u0 {1,S}
7   Cd  u0 {3,D}
8   O2d u0 {4,D}
""",
    thermo = 'Cs-CdCd(CCO)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1089,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cds)CtH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Ct  u0 {1,S}
6   H   u0 {1,S}
7   Cd  u0 {3,D}
8   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1090,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cds)CtH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Ct  u0 {1,S}
6   H   u0 {1,S}
7   Cd  u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CtH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1091,
    label = "Cs-(Cds-Cdd)(Cds-Cdd)CtH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Ct  u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {2,D}
7   Cdd u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CtH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1092,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)CtH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Ct  u0 {1,S}
7   H   u0 {1,S}
8   O2d u0 {4,D}
9   O2d u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-O2d)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1093,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-Cd)CtH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Ct  u0 {1,S}
7   H   u0 {1,S}
8   O2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)(Cds-Cds)CtH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1094,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-S2d)CtH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Ct  u0 {1,S}
7   H   u0 {1,S}
8   S2d u0 {4,D}
9   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1095,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-Cd)CtH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Ct  u0 {1,S}
7   H   u0 {1,S}
8   S2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1096,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CtH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Ct  u0 {1,S}
7   H   u0 {1,S}
8   C   u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CtH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1097,
    label = "Cs-CbCdsCdsH",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Cb      u0 {1,S}
3   [Cd,CO] u0 {1,S}
4   [Cd,CO] u0 {1,S}
5   H       u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CbH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1098,
    label = "Cs-(Cds-O2d)(Cds-O2d)CbH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   CO  u0 {1,S} {7,D}
4   Cb  u0 {1,S}
5   H   u0 {1,S}
6   O2d u0 {2,D}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1099,
    label = "Cs-(Cds-O2d)(Cds-Cd)CbH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   H   u0 {1,S}
6   C   u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)CbH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1100,
    label = "Cs-(Cds-O2d)(Cds-Cds)CbH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   H   u0 {1,S}
6   Cd  u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1101,
    label = "Cs-(Cds-O2d)(Cds-Cdd)CbH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-Cd)CbH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1102,
    label = "Cs-(Cds-O2d)(Cds-Cdd-O2d)CbH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CO  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   H   u0 {1,S}
7   O2d u0 {3,D}
8   O2d u0 {4,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1103,
    label = "Cs-(Cds-O2d)(Cds-Cdd-Cd)CbH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CO  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   H   u0 {1,S}
7   O2d u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)CbH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1104,
    label = "Cs-(Cds-Cd)(Cds-Cd)CbH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cd u0 {1,S} {7,D}
4   Cb u0 {1,S}
5   H  u0 {1,S}
6   C  u0 {2,D}
7   C  u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CbH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1105,
    label = "Cs-(Cds-Cds)(Cds-Cds)CbH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cd u0 {1,S} {7,D}
4   Cb u0 {1,S}
5   H  u0 {1,S}
6   Cd u0 {2,D}
7   Cd u0 {3,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.12,6.51,8.24,9,10.03,10.53,10.89],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (-1.39,'kcal/mol','+|-',0.27),
        S298 = (-11.39,'cal/(mol*K)','+|-',0.15),
    ),
    shortDesc = """Cs-CbCdCdH BOZZELLI =3D Cs/Cs/Cd2/H + (Cs/Cs2/Cb/H - Cs/Cs3/H)""",
    longDesc = 
"""

""",
)

entry(
    index = 1106,
    label = "Cs-(Cds-Cdd)(Cds-Cds)CbH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cb  u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {2,D}
7   Cd  u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cds)CbH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1107,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cds)CbH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   H   u0 {1,S}
7   Cd  u0 {3,D}
8   O2d u0 {4,D}
""",
    thermo = 'Cs-CdCd(CCO)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1108,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cds)CbH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   H   u0 {1,S}
7   Cd  u0 {3,D}
8   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1109,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cds)CbH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   H   u0 {1,S}
7   Cd  u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CbH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1110,
    label = "Cs-(Cds-Cdd)(Cds-Cdd)CbH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cb  u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {2,D}
7   Cdd u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CbH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1111,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)CbH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cb  u0 {1,S}
7   H   u0 {1,S}
8   O2d u0 {4,D}
9   O2d u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-O2d)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1112,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-Cd)CbH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cb  u0 {1,S}
7   H   u0 {1,S}
8   O2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)(Cds-Cds)CbH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1113,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-S2d)CbH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cb  u0 {1,S}
7   H   u0 {1,S}
8   S2d u0 {4,D}
9   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1114,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-Cd)CbH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cb  u0 {1,S}
7   H   u0 {1,S}
8   S2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1115,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CbH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cb  u0 {1,S}
7   H   u0 {1,S}
8   C   u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CbH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1116,
    label = "Cs-CtCtCdsH",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Ct      u0 {1,S}
3   Ct      u0 {1,S}
4   [Cd,CO] u0 {1,S}
5   H       u0 {1,S}
""",
    thermo = 'Cs-CtCt(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1117,
    label = "Cs-CtCt(Cds-O2d)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   Ct  u0 {1,S}
5   H   u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1118,
    label = "Cs-CtCt(Cds-Cd)H",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Ct u0 {1,S}
4   Ct u0 {1,S}
5   H  u0 {1,S}
6   C  u0 {2,D}
""",
    thermo = 'Cs-CtCt(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1119,
    label = "Cs-CtCt(Cds-Cds)H",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Ct u0 {1,S}
4   Ct u0 {1,S}
5   H  u0 {1,S}
6   Cd u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.58,5.68,7.11,8.12,9.27,10.13,9.44],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (4.73,'kcal/mol','+|-',0.27),
        S298 = (-11.46,'cal/(mol*K)','+|-',0.15),
    ),
    shortDesc = """Cs-CtCtCdH RAMAN & GREEN JPCA 2002, 106, 11141-11149""",
    longDesc = 
"""

""",
)

entry(
    index = 1120,
    label = "Cs-CtCt(Cds-Cdd)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   Ct  u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = 'Cs-CtCt(Cds-Cdd-Cd)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1121,
    label = "Cs-CtCt(Cds-Cdd-O2d)H",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Ct  u0 {1,S}
5   Ct  u0 {1,S}
6   H   u0 {1,S}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-CdCd(CCO)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1122,
    label = "Cs-CtCt(Cds-Cdd-S2d)H",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Ct  u0 {1,S}
5   Ct  u0 {1,S}
6   H   u0 {1,S}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1123,
    label = "Cs-CtCt(Cds-Cdd-Cd)H",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Ct  u0 {1,S}
5   Ct  u0 {1,S}
6   H   u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = 'Cs-CtCt(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1124,
    label = "Cs-CbCtCdsH",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Cb      u0 {1,S}
3   Ct      u0 {1,S}
4   [Cd,CO] u0 {1,S}
5   H       u0 {1,S}
""",
    thermo = 'Cs-CbCt(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1125,
    label = "Cs-CbCt(Cds-O2d)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Ct  u0 {1,S}
5   H   u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)CtH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1126,
    label = "Cs-CbCt(Cds-Cd)H",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cb u0 {1,S}
4   Ct u0 {1,S}
5   H  u0 {1,S}
6   C  u0 {2,D}
""",
    thermo = 'Cs-CbCt(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1127,
    label = "Cs-CbCt(Cds-Cds)H",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cb u0 {1,S}
4   Ct u0 {1,S}
5   H  u0 {1,S}
6   Cd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CtH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1128,
    label = "Cs-CbCt(Cds-Cdd)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Ct  u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = 'Cs-CbCt(Cds-Cdd-Cd)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1129,
    label = "Cs-CbCt(Cds-Cdd-O2d)H",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Ct  u0 {1,S}
6   H   u0 {1,S}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)(Cds-Cds)CtH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1130,
    label = "Cs-CbCt(Cds-Cdd-S2d)H",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Ct  u0 {1,S}
6   H   u0 {1,S}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1131,
    label = "Cs-CbCt(Cds-Cdd-Cd)H",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Ct  u0 {1,S}
6   H   u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = 'Cs-CbCt(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1132,
    label = "Cs-CbCbCdsH",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Cb      u0 {1,S}
3   Cb      u0 {1,S}
4   [Cd,CO] u0 {1,S}
5   H       u0 {1,S}
""",
    thermo = 'Cs-CbCb(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1133,
    label = "Cs-CbCb(Cds-O2d)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Cb  u0 {1,S}
5   H   u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1134,
    label = "Cs-CbCb(Cds-Cd)H",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cb u0 {1,S}
4   Cb u0 {1,S}
5   H  u0 {1,S}
6   C  u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1135,
    label = "Cs-CbCb(Cds-Cds)H",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cb u0 {1,S}
4   Cb u0 {1,S}
5   H  u0 {1,S}
6   Cd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1136,
    label = "Cs-CbCb(Cds-Cdd)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Cb  u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = 'Cs-CbCb(Cds-Cdd-Cd)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1137,
    label = "Cs-CbCb(Cds-Cdd-O2d)H",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Cb  u0 {1,S}
6   H   u0 {1,S}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-CdCd(CCO)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1138,
    label = "Cs-CbCb(Cds-Cdd-S2d)H",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Cb  u0 {1,S}
6   H   u0 {1,S}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1139,
    label = "Cs-CbCb(Cds-Cdd-Cd)H",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Cb  u0 {1,S}
6   H   u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = 'Cs-CbCb(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1140,
    label = "Cs-CtCtCtH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Ct u0 {1,S}
3   Ct u0 {1,S}
4   Ct u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.03,5.27,6.78,7.88,9.14,10.08,8.47],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (10.11,'kcal/mol','+|-',0.27),
        S298 = (-10.46,'cal/(mol*K)','+|-',0.15),
    ),
    shortDesc = """Cs-CtCtCtH RAMAN & GREEN JPCA 2002, 106, 11141-11149""",
    longDesc = 
"""

""",
)

entry(
    index = 1141,
    label = "Cs-CbCtCtH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cb u0 {1,S}
3   Ct u0 {1,S}
4   Ct u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = 'Cs-CtCt(Cds-Cds)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1142,
    label = "Cs-CbCbCtH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cb u0 {1,S}
3   Cb u0 {1,S}
4   Ct u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CtH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1143,
    label = "Cs-CbCbCbH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cb u0 {1,S}
3   Cb u0 {1,S}
4   Cb u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.56,7.98,9.36,10.15,10.57,10.65,9.7],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (-0.34,'kcal/mol','+|-',0.27),
        S298 = (-12.31,'cal/(mol*K)','+|-',0.15),
    ),
    shortDesc = """Cs-CbCbCbH BOZZELLI =3D Cs/Cs/Cb2/H + (Cs/Cs2/Cb/H - Cs/Cs3/H)""",
    longDesc = 
"""

""",
)

entry(
    index = 1144,
    label = "Cs-C=SC=SCbH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   CS  u0 {1,S} {7,D}
4   Cb  u0 {1,S}
5   H   u0 {1,S}
6   S2d u0 {2,D}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1145,
    label = "Cs-C=S(Cds-Cd)(Cds-Cd)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   H   u0 {1,S}
6   C   u0 {3,D}
7   C   u0 {4,D}
8   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1146,
    label = "Cs-C=S(Cds-Cdd)(Cds-Cds)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   H   u0 {1,S}
6   Cdd u0 {3,D}
7   Cd  u0 {4,D}
8   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1147,
    label = "Cs-C=S(Cds-Cdd-Cd)(Cds-Cds)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CS  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {7,D}
5   Cdd u0 {2,D} {9,D}
6   H   u0 {1,S}
7   Cd  u0 {4,D}
8   S2d u0 {3,D}
9   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1148,
    label = "Cs-C=S(Cds-Cdd-S2d)(Cds-Cds)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CS  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {7,D}
5   Cdd u0 {2,D} {9,D}
6   H   u0 {1,S}
7   Cd  u0 {4,D}
8   S2d u0 {3,D}
9   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1149,
    label = "Cs-C=S(Cds-Cds)(Cds-Cds)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   H   u0 {1,S}
6   Cd  u0 {3,D}
7   Cd  u0 {4,D}
8   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1150,
    label = "Cs-C=S(Cds-Cdd)(Cds-Cdd)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   H   u0 {1,S}
6   Cdd u0 {3,D}
7   Cdd u0 {4,D}
8   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1151,
    label = "Cs-C=S(Cds-Cdd-Cd)(Cds-Cdd-Cd)H",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    CS  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    H   u0 {1,S}
8    S2d u0 {4,D}
9    C   u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1152,
    label = "Cs-C=S(Cds-Cdd-S2d)(Cds-Cdd-S2d)H",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    CS  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    H   u0 {1,S}
8    S2d u0 {4,D}
9    S2d u0 {5,D}
10   S2d u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1153,
    label = "Cs-C=S(Cds-Cdd-S2d)(Cds-Cdd-Cd)H",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    CS  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    H   u0 {1,S}
8    S2d u0 {4,D}
9    S2d u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1154,
    label = "Cs-C=S(Cds-Cd)CtH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Ct  u0 {1,S}
5   H   u0 {1,S}
6   C   u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1155,
    label = "Cs-C=S(Cds-Cdd)CtH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Ct  u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1156,
    label = "Cs-C=S(Cds-Cdd-S2d)CtH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CS  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Ct  u0 {1,S}
6   H   u0 {1,S}
7   S2d u0 {3,D}
8   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1157,
    label = "Cs-C=S(Cds-Cdd-Cd)CtH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CS  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Ct  u0 {1,S}
6   H   u0 {1,S}
7   S2d u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1158,
    label = "Cs-C=S(Cds-Cds)CtH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Ct  u0 {1,S}
5   H   u0 {1,S}
6   Cd  u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1159,
    label = "Cs-C=SC=SCtH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   CS  u0 {1,S} {7,D}
4   Ct  u0 {1,S}
5   H   u0 {1,S}
6   S2d u0 {2,D}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1160,
    label = "Cs-C=SCtCsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   Cs  u0 {1,S}
5   H   u0 {1,S}
6   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1161,
    label = "Cs-C=SC=SCsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   CS  u0 {1,S} {7,D}
4   Cs  u0 {1,S}
5   H   u0 {1,S}
6   S2d u0 {2,D}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1162,
    label = "Cs-C=S(Cds-Cd)CbH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   H   u0 {1,S}
6   C   u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1163,
    label = "Cs-C=S(Cds-Cds)CbH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   H   u0 {1,S}
6   Cd  u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1164,
    label = "Cs-C=S(Cds-Cdd)CbH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1165,
    label = "Cs-C=S(Cds-Cdd-S2d)CbH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CS  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   H   u0 {1,S}
7   S2d u0 {3,D}
8   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1166,
    label = "Cs-C=S(Cds-Cdd-Cd)CbH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CS  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   H   u0 {1,S}
7   S2d u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1167,
    label = "Cs-C=S(Cds-Cd)CsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cs  u0 {1,S}
5   H   u0 {1,S}
6   C   u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1168,
    label = "Cs-C=S(Cds-Cds)CsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cs  u0 {1,S}
5   H   u0 {1,S}
6   Cd  u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1169,
    label = "Cs-C=S(Cds-Cdd)CsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cs  u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1170,
    label = "Cs-C=S(Cds-Cdd-Cd)CsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CS  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cs  u0 {1,S}
6   H   u0 {1,S}
7   S2d u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1171,
    label = "Cs-C=S(Cds-Cdd-S2d)CsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CS  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cs  u0 {1,S}
6   H   u0 {1,S}
7   S2d u0 {3,D}
8   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1172,
    label = "Cs-CbCtC=SH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Ct  u0 {1,S}
5   H   u0 {1,S}
6   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1173,
    label = "Cs-C=SC=SC=SH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   CS  u0 {1,S} {7,D}
4   CS  u0 {1,S} {8,D}
5   H   u0 {1,S}
6   S2d u0 {2,D}
7   S2d u0 {3,D}
8   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1174,
    label = "Cs-C=SCsCsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   Cs  u0 {1,S}
5   H   u0 {1,S}
6   S2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([2.36,5.4,7.24,8.24,9.35,9.91,10.56],'cal/(mol*K)'),
        H298 = (2.29,'kcal/mol'),
        S298 = (-10.76,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 1175,
    label = "Cs-CtCtC=SH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   Ct  u0 {1,S}
5   H   u0 {1,S}
6   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1176,
    label = "Cs-CbCbC=SH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Cb  u0 {1,S}
5   H   u0 {1,S}
6   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1177,
    label = "Cs-C=SC=S(Cds-Cd)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   CS  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {6,D}
5   H   u0 {1,S}
6   C   u0 {4,D}
7   S2d u0 {2,D}
8   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1178,
    label = "Cs-C=SC=S(Cds-Cds)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   CS  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {6,D}
5   H   u0 {1,S}
6   Cd  u0 {4,D}
7   S2d u0 {2,D}
8   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1179,
    label = "Cs-C=SC=S(Cds-Cdd)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   CS  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {6,D}
5   H   u0 {1,S}
6   Cdd u0 {4,D}
7   S2d u0 {2,D}
8   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1180,
    label = "Cs-C=SC=S(Cds-Cdd-S2d)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CS  u0 {1,S} {7,D}
4   CS  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   H   u0 {1,S}
7   S2d u0 {3,D}
8   S2d u0 {4,D}
9   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1181,
    label = "Cs-C=SC=S(Cds-Cdd-Cd)H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CS  u0 {1,S} {7,D}
4   CS  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   H   u0 {1,S}
7   S2d u0 {3,D}
8   S2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1182,
    label = "Cs-CCCC",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   C  u0 {1,S}
4   C  u0 {1,S}
5   C  u0 {1,S}
""",
    thermo = 'Cs-CsCsCsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1183,
    label = "Cs-CsCsCsCs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cs u0 {1,S}
3   Cs u0 {1,S}
4   Cs u0 {1,S}
5   Cs u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.37,6.13,7.36,8.12,8.77,8.76,8.12],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (0.5,'kcal/mol','+|-',0.27),
        S298 = (-35.1,'cal/(mol*K)','+|-',0.15),
    ),
    shortDesc = """Cs-CsCsCsCs BENSON""",
    longDesc = 
"""

""",
)

entry(
    index = 1184,
    label = "Cs-CdsCsCsCs",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   [Cd,CO] u0 {1,S}
3   Cs      u0 {1,S}
4   Cs      u0 {1,S}
5   Cs      u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CsCsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1185,
    label = "Cs-(Cds-O2d)CsCsCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   Cs  u0 {1,S}
5   Cs  u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([22.68,27.48,30.12,31.51,32.36,32.39,32.42],'J/(mol*K)','+|-',[3.34,3.34,3.34,3.34,3.34,3.34,3.34]),
        H298 = (4.6,'kJ/mol','+|-',2.85),
        S298 = (-140.94,'J/(mol*K)','+|-',3.9),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 1186,
    label = "Cs-(Cds-Cd)CsCsCs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cs u0 {1,S}
4   Cs u0 {1,S}
5   Cs u0 {1,S}
6   C  u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)CsCsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1187,
    label = "Cs-(Cds-Cds)CsCsCs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cs u0 {1,S}
4   Cs u0 {1,S}
5   Cs u0 {1,S}
6   Cd u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.99,6.04,7.43,8.26,8.92,8.96,8.23],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (1.68,'kcal/mol','+|-',0.27),
        S298 = (-34.72,'cal/(mol*K)','+|-',0.15),
    ),
    shortDesc = """Cs-CdCsCsCs BENSON""",
    longDesc = 
"""

""",
)

entry(
    index = 1188,
    label = "Cs-(Cds-Cdd)CsCsCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   Cs  u0 {1,S}
5   Cs  u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)CsCsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1189,
    label = "Cs-(Cds-Cdd-O2d)CsCsCs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cs  u0 {1,S}
5   Cs  u0 {1,S}
6   Cs  u0 {1,S}
7   O2d u0 {3,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([20.63,27.65,31.98,34.41,36.16,36.25,35.2],'J/(mol*K)','+|-',[6.93,6.93,6.93,6.93,6.93,6.93,6.93]),
        H298 = (-4.5,'kJ/mol','+|-',5.9),
        S298 = (-144.08,'J/(mol*K)','+|-',8.08),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 1190,
    label = "Cs-(Cds-Cdd-S2d)CsCsCs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cs  u0 {1,S}
5   Cs  u0 {1,S}
6   Cs  u0 {1,S}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1191,
    label = "Cs-(Cds-Cdd-Cd)CsCsCs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cs  u0 {1,S}
5   Cs  u0 {1,S}
6   Cs  u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)CsCsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1192,
    label = "Cs-CtCsCsCs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Ct u0 {1,S}
3   Cs u0 {1,S}
4   Cs u0 {1,S}
5   Cs u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.37,6.79,8.09,8.78,9.19,8.96,7.63],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (2.81,'kcal/mol','+|-',0.27),
        S298 = (-35.18,'cal/(mol*K)','+|-',0.15),
    ),
    shortDesc = """Cs-CtCsCsCs BENSON""",
    longDesc = 
"""

""",
)

entry(
    index = 1193,
    label = "Cs-CbCsCsCs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cb u0 {1,S}
3   Cs u0 {1,S}
4   Cs u0 {1,S}
5   Cs u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.37,6.79,8.09,8.78,9.19,8.96,7.63],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (2.81,'kcal/mol','+|-',0.26),
        S298 = (-35.18,'cal/(mol*K)','+|-',0.13),
    ),
    shortDesc = """Cs-CbCsCsCs BENSON""",
    longDesc = 
"""

""",
)

entry(
    index = 1194,
    label = "Cs-CdsCdsCsCs",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   [Cd,CO] u0 {1,S}
3   [Cd,CO] u0 {1,S}
4   Cs      u0 {1,S}
5   Cs      u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1195,
    label = "Cs-(Cds-O2d)(Cds-O2d)CsCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   CO  u0 {1,S} {7,D}
4   Cs  u0 {1,S}
5   Cs  u0 {1,S}
6   O2d u0 {2,D}
7   O2d u0 {3,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([33.76,33.42,32.6,31.91,31.01,30.55,30.35],'J/(mol*K)','+|-',[5.08,5.08,5.08,5.08,5.08,5.08,5.08]),
        H298 = (14.9,'kJ/mol','+|-',4.33),
        S298 = (-146.69,'J/(mol*K)','+|-',5.92),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 1196,
    label = "Cs-(Cds-O2d)(Cds-Cd)CsCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cs  u0 {1,S}
5   Cs  u0 {1,S}
6   C   u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([26.01,30.13,32.44,33.51,33.75,33.26,32.55],'J/(mol*K)','+|-',[3.34,3.34,3.34,3.34,3.34,3.34,3.34]),
        H298 = (9.8,'kJ/mol','+|-',2.85),
        S298 = (-146.74,'J/(mol*K)','+|-',3.9),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 1197,
    label = "Cs-(Cds-O2d)(Cds-Cds)CsCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cs  u0 {1,S}
5   Cs  u0 {1,S}
6   Cd  u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)CsCsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1198,
    label = "Cs-(Cds-O2d)(Cds-Cdd)CsCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cs  u0 {1,S}
5   Cs  u0 {1,S}
6   Cdd u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-Cd)CsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1199,
    label = "Cs-(Cds-O2d)(Cds-Cdd-O2d)CsCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CO  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cs  u0 {1,S}
6   Cs  u0 {1,S}
7   O2d u0 {3,D}
8   O2d u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)CsCsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1200,
    label = "Cs-(Cds-O2d)(Cds-Cdd-Cd)CsCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CO  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cs  u0 {1,S}
6   Cs  u0 {1,S}
7   O2d u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)CsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1201,
    label = "Cs-(Cds-Cd)(Cds-Cd)CsCs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cd u0 {1,S} {7,D}
4   Cs u0 {1,S}
5   Cs u0 {1,S}
6   C  u0 {2,D}
7   C  u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1202,
    label = "Cs-(Cds-Cds)(Cds-Cds)CsCs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cd u0 {1,S} {7,D}
4   Cs u0 {1,S}
5   Cs u0 {1,S}
6   Cd u0 {2,D}
7   Cd u0 {3,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.99,6.04,7.43,8.26,8.92,8.96,8.23],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (1.68,'kcal/mol','+|-',0.26),
        S298 = (-34.72,'cal/(mol*K)','+|-',0.13),
    ),
    shortDesc = """Cs-CdCdCsCs BENSON""",
    longDesc = 
"""

""",
)

entry(
    index = 1203,
    label = "Cs-(Cds-Cdd)(Cds-Cds)CsCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cs  u0 {1,S}
5   Cs  u0 {1,S}
6   Cdd u0 {2,D}
7   Cd  u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cds)CsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1204,
    label = "Cs-CsCsCd(CCO)",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cs  u0 {1,S}
6   Cs  u0 {1,S}
7   Cd  u0 {3,D}
8   O2d u0 {4,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([25.48,31.89,35.19,36.68,37.19,36.66,34.96],'J/(mol*K)','+|-',[6.93,6.93,6.93,6.93,6.93,6.93,6.93]),
        H298 = (2.9,'kJ/mol','+|-',5.9),
        S298 = (-144.6,'J/(mol*K)','+|-',8.08),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 1205,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cds)CsCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cs  u0 {1,S}
6   Cs  u0 {1,S}
7   Cd  u0 {3,D}
8   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1206,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cds)CsCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cs  u0 {1,S}
6   Cs  u0 {1,S}
7   Cd  u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1207,
    label = "Cs-(Cds-Cdd)(Cds-Cdd)CsCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cs  u0 {1,S}
5   Cs  u0 {1,S}
6   Cdd u0 {2,D}
7   Cdd u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1208,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)CsCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cs  u0 {1,S}
7   Cs  u0 {1,S}
8   O2d u0 {4,D}
9   O2d u0 {5,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.73,8.1,9.02,9.53,9.66,9.52,8.93],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (-2.987,'kcal/mol','+|-',0.26),
        S298 = (-36.46,'cal/(mol*K)','+|-',0.13),
    ),
    shortDesc = """{C/C2/CCO2} RAMAN & GREEN JPCA 2002, 106, 7937-7949""",
    longDesc = 
"""

""",
)

entry(
    index = 1209,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-Cd)CsCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cs  u0 {1,S}
7   Cs  u0 {1,S}
8   O2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-CsCsCd(CCO)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1210,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-S2d)CsCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cs  u0 {1,S}
7   Cs  u0 {1,S}
8   S2d u0 {4,D}
9   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1211,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-Cd)CsCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cs  u0 {1,S}
7   Cs  u0 {1,S}
8   S2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1212,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CsCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cs  u0 {1,S}
7   Cs  u0 {1,S}
8   C   u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1213,
    label = "Cs-CtCdsCsCs",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Ct      u0 {1,S}
3   [Cd,CO] u0 {1,S}
4   Cs      u0 {1,S}
5   Cs      u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CtCsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1214,
    label = "Cs-(Cds-O2d)CtCsCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   Cs  u0 {1,S}
5   Cs  u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)CsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1215,
    label = "Cs-(Cds-Cd)CtCsCs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Ct u0 {1,S}
4   Cs u0 {1,S}
5   Cs u0 {1,S}
6   C  u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)CtCsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1216,
    label = "Cs-(Cds-Cds)CtCsCs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Ct u0 {1,S}
4   Cs u0 {1,S}
5   Cs u0 {1,S}
6   Cd u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.99,6.7,8.16,8.92,9.34,9.16,7.14],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (2.99,'kcal/mol','+|-',0.26),
        S298 = (-34.8,'cal/(mol*K)','+|-',0.13),
    ),
    shortDesc = """Cs-CtCdCsCs BOZZELLI =3D Cs/Cs3/Cd + (Cs/Cs3/Ct - Cs/Cs4)""",
    longDesc = 
"""

""",
)

entry(
    index = 1217,
    label = "Cs-(Cds-Cdd)CtCsCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   Cs  u0 {1,S}
5   Cs  u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)CtCsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1218,
    label = "Cs-(Cds-Cdd-O2d)CtCsCs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Ct  u0 {1,S}
5   Cs  u0 {1,S}
6   Cs  u0 {1,S}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-CsCsCd(CCO)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1219,
    label = "Cs-(Cds-Cdd-S2d)CtCsCs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Ct  u0 {1,S}
5   Cs  u0 {1,S}
6   Cs  u0 {1,S}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1220,
    label = "Cs-(Cds-Cdd-Cd)CtCsCs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Ct  u0 {1,S}
5   Cs  u0 {1,S}
6   Cs  u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)CtCsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1221,
    label = "Cs-CbCdsCsCs",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Cb      u0 {1,S}
3   [Cd,CO] u0 {1,S}
4   Cs      u0 {1,S}
5   Cs      u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CbCsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1222,
    label = "Cs-(Cds-O2d)CbCsCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Cs  u0 {1,S}
5   Cs  u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)CsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1223,
    label = "Cs-(Cds-Cd)CbCsCs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cb u0 {1,S}
4   Cs u0 {1,S}
5   Cs u0 {1,S}
6   C  u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)CbCsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1224,
    label = "Cs-(Cds-Cds)CbCsCs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cb u0 {1,S}
4   Cs u0 {1,S}
5   Cs u0 {1,S}
6   Cd u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.99,6.7,8.16,8.92,9.34,9.16,7.14],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (2.99,'kcal/mol','+|-',0.26),
        S298 = (-34.8,'cal/(mol*K)','+|-',0.13),
    ),
    shortDesc = """Cs-CbCdCsCs BOZZELLI =3D Cs/Cs3/Cb + (Cs/Cs3/Cd - Cs/Cs4)""",
    longDesc = 
"""

""",
)

entry(
    index = 1225,
    label = "Cs-(Cds-Cdd)CbCsCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Cs  u0 {1,S}
5   Cs  u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)CbCsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1226,
    label = "Cs-(Cds-Cdd-O2d)CbCsCs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Cs  u0 {1,S}
6   Cs  u0 {1,S}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-CsCsCd(CCO)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1227,
    label = "Cs-(Cds-Cdd-S2d)CbCsCs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Cs  u0 {1,S}
6   Cs  u0 {1,S}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1228,
    label = "Cs-(Cds-Cdd-Cd)CbCsCs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Cs  u0 {1,S}
6   Cs  u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)CbCsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1229,
    label = "Cs-CtCtCsCs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Ct u0 {1,S}
3   Ct u0 {1,S}
4   Cs u0 {1,S}
5   Cs u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.57,5.98,7.51,8.37,9,9.02,8.34],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (1.16,'kcal/mol','+|-',0.26),
        S298 = (-35.26,'cal/(mol*K)','+|-',0.13),
    ),
    shortDesc = """Cs-CtCtCsCs BOZZELLI =3D Cs/Cs3/Ct + (Cs/Cs3/Ct - Cs/Cs4)""",
    longDesc = 
"""

""",
)

entry(
    index = 1230,
    label = "Cs-CbCtCsCs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cb u0 {1,S}
3   Ct u0 {1,S}
4   Cs u0 {1,S}
5   Cs u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.57,5.98,7.51,8.37,9,9.02,8.34],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (1.16,'kcal/mol','+|-',0.26),
        S298 = (-35.26,'cal/(mol*K)','+|-',0.13),
    ),
    shortDesc = """Cs-CbCtCsCs BOZZELLI =3D Cs/Cs3/Cb + (Cs/Cs3/Ct - Cs/Cs4)""",
    longDesc = 
"""

""",
)

entry(
    index = 1231,
    label = "Cs-CbCbCsCs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cb u0 {1,S}
3   Cb u0 {1,S}
4   Cs u0 {1,S}
5   Cs u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.57,5.98,7.51,8.37,9,9.02,8.34],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (1.16,'kcal/mol','+|-',0.26),
        S298 = (-35.26,'cal/(mol*K)','+|-',0.13),
    ),
    shortDesc = """Cs-CbCbCsCs BENSON""",
    longDesc = 
"""

""",
)

entry(
    index = 1232,
    label = "Cs-CdsCdsCdsCs",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   [Cd,CO] u0 {1,S}
3   [Cd,CO] u0 {1,S}
4   [Cd,CO] u0 {1,S}
5   Cs      u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1233,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-O2d)Cs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   CO  u0 {1,S} {7,D}
4   CO  u0 {1,S} {8,D}
5   Cs  u0 {1,S}
6   O2d u0 {2,D}
7   O2d u0 {3,D}
8   O2d u0 {4,D}
""",
    thermo = 'Cs-CsCsCsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1234,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-Cd)Cs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   CO  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {6,D}
5   Cs  u0 {1,S}
6   C   u0 {4,D}
7   O2d u0 {2,D}
8   O2d u0 {3,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([35.99,39.53,39.94,39.09,36.71,34.8,32.51],'J/(mol*K)','+|-',[5.08,5.08,5.08,5.08,5.08,5.08,5.08]),
        H298 = (19.9,'kJ/mol','+|-',4.33),
        S298 = (-150.69,'J/(mol*K)','+|-',5.92),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 1235,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-Cds)Cs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   CO  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {6,D}
5   Cs  u0 {1,S}
6   Cd  u0 {4,D}
7   O2d u0 {2,D}
8   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)CsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1236,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-Cdd)Cs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   CO  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {6,D}
5   Cs  u0 {1,S}
6   Cdd u0 {4,D}
7   O2d u0 {2,D}
8   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-Cdd-Cd)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1237,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-Cdd-O2d)Cs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CO  u0 {1,S} {7,D}
4   CO  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   Cs  u0 {1,S}
7   O2d u0 {3,D}
8   O2d u0 {4,D}
9   O2d u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)CsCsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1238,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-Cdd-Cd)Cs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CO  u0 {1,S} {7,D}
4   CO  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   Cs  u0 {1,S}
7   O2d u0 {3,D}
8   O2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1239,
    label = "Cs-(Cds-O2d)(Cds-Cd)(Cds-Cd)Cs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Cs  u0 {1,S}
6   C   u0 {3,D}
7   C   u0 {4,D}
8   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1240,
    label = "Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)Cs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Cs  u0 {1,S}
6   Cd  u0 {3,D}
7   Cd  u0 {4,D}
8   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)CsCsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1241,
    label = "Cs-(Cds-O2d)(Cds-Cdd)(Cds-Cds)Cs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Cs  u0 {1,S}
6   Cdd u0 {3,D}
7   Cd  u0 {4,D}
8   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-Cd)(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1242,
    label = "Cs-(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cds)Cs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CO  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {7,D}
5   Cdd u0 {2,D} {9,D}
6   Cs  u0 {1,S}
7   Cd  u0 {4,D}
8   O2d u0 {3,D}
9   O2d u0 {5,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-O2d)CsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1243,
    label = "Cs-(Cds-O2d)(Cds-Cdd-Cd)(Cds-Cds)Cs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CO  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {7,D}
5   Cdd u0 {2,D} {9,D}
6   Cs  u0 {1,S}
7   Cd  u0 {4,D}
8   O2d u0 {3,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1244,
    label = "Cs-(Cds-O2d)(Cds-Cdd)(Cds-Cdd)Cs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Cs  u0 {1,S}
6   Cdd u0 {3,D}
7   Cdd u0 {4,D}
8   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1245,
    label = "Cs-(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cdd-O2d)Cs",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    CO  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cs  u0 {1,S}
8    O2d u0 {4,D}
9    O2d u0 {5,D}
10   O2d u0 {6,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)CsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1246,
    label = "Cs-(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cdd-Cd)Cs",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    CO  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cs  u0 {1,S}
8    O2d u0 {4,D}
9    O2d u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1247,
    label = "Cs-(Cds-O2d)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cs",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    CO  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cs  u0 {1,S}
8    O2d u0 {4,D}
9    C   u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1248,
    label = "Cs-(Cds-Cd)(Cds-Cd)(Cds-Cd)Cs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cd u0 {1,S} {7,D}
4   Cd u0 {1,S} {8,D}
5   Cs u0 {1,S}
6   C  u0 {2,D}
7   C  u0 {3,D}
8   C  u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1249,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Cs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cd u0 {1,S} {7,D}
4   Cd u0 {1,S} {8,D}
5   Cs u0 {1,S}
6   Cd u0 {2,D}
7   Cd u0 {3,D}
8   Cd u0 {4,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.32,5.86,7.57,8.54,9.22,9.36,8.45],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (2.54,'kcal/mol','+|-',0.26),
        S298 = (-33.96,'cal/(mol*K)','+|-',0.13),
    ),
    shortDesc = """Cs-CdCdCdCs BOZZELLI =3D Cs/Cs2/Cd2 + (Cs/Cs3/Cd - Cs/Cs4)""",
    longDesc = 
"""

""",
)

entry(
    index = 1250,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd)Cs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   Cs  u0 {1,S}
6   Cd  u0 {2,D}
7   Cd  u0 {3,D}
8   Cdd u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1251,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)Cs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   Cs  u0 {1,S}
7   Cd  u0 {3,D}
8   Cd  u0 {4,D}
9   O2d u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)CsCsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1252,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-S2d)Cs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   Cs  u0 {1,S}
7   Cd  u0 {3,D}
8   Cd  u0 {4,D}
9   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1253,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)Cs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   Cs  u0 {1,S}
7   Cd  u0 {3,D}
8   Cd  u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1254,
    label = "Cs-(Cds-Cds)(Cds-Cdd)(Cds-Cdd)Cs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   Cs  u0 {1,S}
6   Cd  u0 {2,D}
7   Cdd u0 {3,D}
8   Cdd u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1255,
    label = "Cs-(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-O2d)Cs",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cs  u0 {1,S}
8    Cd  u0 {4,D}
9    O2d u0 {5,D}
10   O2d u0 {6,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)CsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1256,
    label = "Cs-(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-Cd)Cs",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cs  u0 {1,S}
8    Cd  u0 {4,D}
9    O2d u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1257,
    label = "Cs-(Cds-Cds)(Cds-Cdd-S2d)(Cds-Cdd-S2d)Cs",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cs  u0 {1,S}
8    Cd  u0 {4,D}
9    S2d u0 {5,D}
10   S2d u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1258,
    label = "Cs-(Cds-Cds)(Cds-Cdd-S2d)(Cds-Cdd-Cd)Cs",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cs  u0 {1,S}
8    Cd  u0 {4,D}
9    S2d u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1259,
    label = "Cs-(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cs",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cs  u0 {1,S}
8    Cd  u0 {4,D}
9    C   u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1260,
    label = "Cs-(Cds-Cdd)(Cds-Cdd)(Cds-Cdd)Cs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   Cs  u0 {1,S}
6   Cdd u0 {2,D}
7   Cdd u0 {3,D}
8   Cdd u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1261,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)(Cds-Cdd-O2d)Cs",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    Cs  u0 {1,S}
9    O2d u0 {5,D}
10   O2d u0 {6,D}
11   O2d u0 {7,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1262,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)(Cds-Cdd-Cd)Cs",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    Cs  u0 {1,S}
9    O2d u0 {5,D}
10   O2d u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-O2d)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1263,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cs",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    Cs  u0 {1,S}
9    O2d u0 {5,D}
10   C   u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1264,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-S2d)(Cds-Cdd-S2d)Cs",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    Cs  u0 {1,S}
9    S2d u0 {5,D}
10   S2d u0 {6,D}
11   S2d u0 {7,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1265,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-S2d)(Cds-Cdd-Cd)Cs",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    Cs  u0 {1,S}
9    S2d u0 {5,D}
10   S2d u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1266,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cs",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    Cs  u0 {1,S}
9    S2d u0 {5,D}
10   C   u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1267,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cs",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    Cs  u0 {1,S}
9    C   u0 {5,D}
10   C   u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1268,
    label = "Cs-CtCdsCdsCs",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Ct      u0 {1,S}
3   [Cd,CO] u0 {1,S}
4   [Cd,CO] u0 {1,S}
5   Cs      u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CtCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1269,
    label = "Cs-(Cds-O2d)(Cds-O2d)CtCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   CO  u0 {1,S} {7,D}
4   Ct  u0 {1,S}
5   Cs  u0 {1,S}
6   O2d u0 {2,D}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1270,
    label = "Cs-(Cds-O2d)(Cds-Cd)CtCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Ct  u0 {1,S}
5   Cs  u0 {1,S}
6   C   u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)CtCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1271,
    label = "Cs-(Cds-O2d)(Cds-Cds)CtCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Ct  u0 {1,S}
5   Cs  u0 {1,S}
6   Cd  u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1272,
    label = "Cs-(Cds-O2d)(Cds-Cdd)CtCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Ct  u0 {1,S}
5   Cs  u0 {1,S}
6   Cdd u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-Cd)CtCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1273,
    label = "Cs-(Cds-O2d)(Cds-Cdd-O2d)CtCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CO  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Ct  u0 {1,S}
6   Cs  u0 {1,S}
7   O2d u0 {3,D}
8   O2d u0 {4,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1274,
    label = "Cs-(Cds-O2d)(Cds-Cdd-Cd)CtCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CO  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Ct  u0 {1,S}
6   Cs  u0 {1,S}
7   O2d u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)CtCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1275,
    label = "Cs-(Cds-Cd)(Cds-Cd)CtCs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cd u0 {1,S} {7,D}
4   Ct u0 {1,S}
5   Cs u0 {1,S}
6   C  u0 {2,D}
7   C  u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CtCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1276,
    label = "Cs-(Cds-Cds)(Cds-Cds)CtCs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cd u0 {1,S} {7,D}
4   Ct u0 {1,S}
5   Cs u0 {1,S}
6   Cd u0 {2,D}
7   Cd u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1277,
    label = "Cs-(Cds-Cdd)(Cds-Cds)CtCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Ct  u0 {1,S}
5   Cs  u0 {1,S}
6   Cdd u0 {2,D}
7   Cd  u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cds)CtCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1278,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cds)CtCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Ct  u0 {1,S}
6   Cs  u0 {1,S}
7   Cd  u0 {3,D}
8   O2d u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1279,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cds)CtCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Ct  u0 {1,S}
6   Cs  u0 {1,S}
7   Cd  u0 {3,D}
8   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1280,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cds)CtCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Ct  u0 {1,S}
6   Cs  u0 {1,S}
7   Cd  u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CtCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1281,
    label = "Cs-(Cds-Cdd)(Cds-Cdd)CtCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Ct  u0 {1,S}
5   Cs  u0 {1,S}
6   Cdd u0 {2,D}
7   Cdd u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CtCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1282,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)CtCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Ct  u0 {1,S}
7   Cs  u0 {1,S}
8   O2d u0 {4,D}
9   O2d u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-O2d)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1283,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-Cd)CtCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Ct  u0 {1,S}
7   Cs  u0 {1,S}
8   O2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)(Cds-Cds)CtCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1284,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-S2d)CtCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Ct  u0 {1,S}
7   Cs  u0 {1,S}
8   S2d u0 {4,D}
9   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1285,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-Cd)CtCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Ct  u0 {1,S}
7   Cs  u0 {1,S}
8   S2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1286,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CtCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Ct  u0 {1,S}
7   Cs  u0 {1,S}
8   C   u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CtCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1287,
    label = "Cs-CbCdsCdsCs",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Cb      u0 {1,S}
3   [Cd,CO] u0 {1,S}
4   [Cd,CO] u0 {1,S}
5   Cs      u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CbCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1288,
    label = "Cs-(Cds-O2d)(Cds-O2d)CbCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   CO  u0 {1,S} {7,D}
4   Cb  u0 {1,S}
5   Cs  u0 {1,S}
6   O2d u0 {2,D}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1289,
    label = "Cs-(Cds-O2d)(Cds-Cd)CbCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   Cs  u0 {1,S}
6   C   u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)CbCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1290,
    label = "Cs-(Cds-O2d)(Cds-Cds)CbCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   Cs  u0 {1,S}
6   Cd  u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1291,
    label = "Cs-(Cds-O2d)(Cds-Cdd)CbCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   Cs  u0 {1,S}
6   Cdd u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-Cd)CbCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1292,
    label = "Cs-(Cds-O2d)(Cds-Cdd-O2d)CbCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CO  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   Cs  u0 {1,S}
7   O2d u0 {3,D}
8   O2d u0 {4,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1293,
    label = "Cs-(Cds-O2d)(Cds-Cdd-Cd)CbCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CO  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   Cs  u0 {1,S}
7   O2d u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)CbCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1294,
    label = "Cs-(Cds-Cd)(Cds-Cd)CbCs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cd u0 {1,S} {7,D}
4   Cb u0 {1,S}
5   Cs u0 {1,S}
6   C  u0 {2,D}
7   C  u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CbCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1295,
    label = "Cs-(Cds-Cds)(Cds-Cds)CbCs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cd u0 {1,S} {7,D}
4   Cb u0 {1,S}
5   Cs u0 {1,S}
6   Cd u0 {2,D}
7   Cd u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1296,
    label = "Cs-(Cds-Cdd)(Cds-Cds)CbCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cb  u0 {1,S}
5   Cs  u0 {1,S}
6   Cdd u0 {2,D}
7   Cd  u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cds)CbCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1297,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cds)CbCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   Cs  u0 {1,S}
7   Cd  u0 {3,D}
8   O2d u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1298,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cds)CbCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   Cs  u0 {1,S}
7   Cd  u0 {3,D}
8   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1299,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cds)CbCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   Cs  u0 {1,S}
7   Cd  u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CbCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1300,
    label = "Cs-(Cds-Cdd)(Cds-Cdd)CbCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cb  u0 {1,S}
5   Cs  u0 {1,S}
6   Cdd u0 {2,D}
7   Cdd u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CbCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1301,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)CbCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cb  u0 {1,S}
7   Cs  u0 {1,S}
8   O2d u0 {4,D}
9   O2d u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-O2d)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1302,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-Cd)CbCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cb  u0 {1,S}
7   Cs  u0 {1,S}
8   O2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)(Cds-Cds)CbCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1303,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-S2d)CbCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cb  u0 {1,S}
7   Cs  u0 {1,S}
8   S2d u0 {4,D}
9   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1304,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-Cd)CbCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cb  u0 {1,S}
7   Cs  u0 {1,S}
8   S2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1305,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CbCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cb  u0 {1,S}
7   Cs  u0 {1,S}
8   C   u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CbCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1306,
    label = "Cs-CtCtCdsCs",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Ct      u0 {1,S}
3   Ct      u0 {1,S}
4   [Cd,CO] u0 {1,S}
5   Cs      u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CtCtCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1307,
    label = "Cs-(Cds-O2d)CtCtCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   Ct  u0 {1,S}
5   Cs  u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1308,
    label = "Cs-(Cds-Cd)CtCtCs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Ct u0 {1,S}
4   Ct u0 {1,S}
5   Cs u0 {1,S}
6   C  u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)CtCtCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1309,
    label = "Cs-(Cds-Cds)CtCtCs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Ct u0 {1,S}
4   Ct u0 {1,S}
5   Cs u0 {1,S}
6   Cd u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.99,7.36,8.89,9.58,9.76,9.16,7.25],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (5.1,'kcal/mol','+|-',0.26),
        S298 = (-34.88,'cal/(mol*K)','+|-',0.13),
    ),
    shortDesc = """Cs-CtCtCdCs BOZZELLI =3D Cs/Cd2/Cs2 + (Cs/Cs3/Ct - Cs/Cs4)""",
    longDesc = 
"""

""",
)

entry(
    index = 1310,
    label = "Cs-(Cds-Cdd)CtCtCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   Ct  u0 {1,S}
5   Cs  u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)CtCtCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1311,
    label = "Cs-(Cds-Cdd-O2d)CtCtCs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Ct  u0 {1,S}
5   Ct  u0 {1,S}
6   Cs  u0 {1,S}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1312,
    label = "Cs-(Cds-Cdd-S2d)CtCtCs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Ct  u0 {1,S}
5   Ct  u0 {1,S}
6   Cs  u0 {1,S}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1313,
    label = "Cs-(Cds-Cdd-Cd)CtCtCs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Ct  u0 {1,S}
5   Ct  u0 {1,S}
6   Cs  u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)CtCtCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1314,
    label = "Cs-CbCtCdsCs",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Cb      u0 {1,S}
3   Ct      u0 {1,S}
4   [Cd,CO] u0 {1,S}
5   Cs      u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CbCtCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1315,
    label = "Cs-(Cds-O2d)CbCtCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Ct  u0 {1,S}
5   Cs  u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)CtCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1316,
    label = "Cs-(Cds-Cd)CbCtCs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cb u0 {1,S}
4   Ct u0 {1,S}
5   Cs u0 {1,S}
6   C  u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)CbCtCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1317,
    label = "Cs-(Cds-Cds)CbCtCs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cb u0 {1,S}
4   Ct u0 {1,S}
5   Cs u0 {1,S}
6   Cd u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.99,7.36,8.89,9.58,9.76,9.16,7.25],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (5.1,'kcal/mol','+|-',0.26),
        S298 = (-34.88,'cal/(mol*K)','+|-',0.13),
    ),
    shortDesc = """Cs-CbCtCdCs BOZZELLI =3D Cs/Cb/Cd/Cs2 + (Cs/Cs3/Ct - Cs/Cs4)""",
    longDesc = 
"""

""",
)

entry(
    index = 1318,
    label = "Cs-(Cds-Cdd)CbCtCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Ct  u0 {1,S}
5   Cs  u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)CbCtCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1319,
    label = "Cs-(Cds-Cdd-O2d)CbCtCs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Ct  u0 {1,S}
6   Cs  u0 {1,S}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)(Cds-Cds)CtCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1320,
    label = "Cs-(Cds-Cdd-S2d)CbCtCs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Ct  u0 {1,S}
6   Cs  u0 {1,S}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1321,
    label = "Cs-(Cds-Cdd-Cd)CbCtCs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Ct  u0 {1,S}
6   Cs  u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.99,7.36,8.89,9.58,9.76,9.16,7.25],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (5.1,'kcal/mol','+|-',0.26),
        S298 = (-34.88,'cal/(mol*K)','+|-',0.13),
    ),
    shortDesc = """Cs-CbCtCdCs BOZZELLI =3D Cs/Cb/Cd/Cs2 + (Cs/Cs3/Ct - Cs/Cs4)""",
    longDesc = 
"""

""",
)

entry(
    index = 1322,
    label = "Cs-CbCbCdsCs",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Cb      u0 {1,S}
3   Cb      u0 {1,S}
4   [Cd,CO] u0 {1,S}
5   Cs      u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CbCbCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1323,
    label = "Cs-(Cds-O2d)CbCbCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Cb  u0 {1,S}
5   Cs  u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1324,
    label = "Cs-(Cds-Cd)CbCbCs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cb u0 {1,S}
4   Cb u0 {1,S}
5   Cs u0 {1,S}
6   C  u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)CbCbCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1325,
    label = "Cs-(Cds-Cds)CbCbCs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cb u0 {1,S}
4   Cb u0 {1,S}
5   Cs u0 {1,S}
6   Cd u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.99,7.36,8.89,9.58,9.76,9.16,7.25],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (5.1,'kcal/mol','+|-',0.26),
        S298 = (-34.88,'cal/(mol*K)','+|-',0.13),
    ),
    shortDesc = """Cs-CbCbCdCs BOZZELLI =3D Cs/Cs2/Cb2 + (Cs/Cs3/Cd - Cs/Cs4)""",
    longDesc = 
"""

""",
)

entry(
    index = 1326,
    label = "Cs-(Cds-Cdd)CbCbCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Cb  u0 {1,S}
5   Cs  u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)CbCbCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1327,
    label = "Cs-(Cds-Cdd-O2d)CbCbCs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Cb  u0 {1,S}
6   Cs  u0 {1,S}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1328,
    label = "Cs-(Cds-Cdd-S2d)CbCbCs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Cb  u0 {1,S}
6   Cs  u0 {1,S}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1329,
    label = "Cs-(Cds-Cdd-Cd)CbCbCs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Cb  u0 {1,S}
6   Cs  u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)CbCbCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1330,
    label = "Cs-CtCtCtCs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Ct u0 {1,S}
3   Ct u0 {1,S}
4   Ct u0 {1,S}
5   Cs u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.37,8.11,9.55,10.1,10.03,9.36,6.65],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (6.23,'kcal/mol','+|-',0.26),
        S298 = (-35.34,'cal/(mol*K)','+|-',0.13),
    ),
    shortDesc = """Cs-CtCtCtCs BOZZELLI =3D Cs/Cs2/Ct2 + (Cs/Cs3/Ct - Cs/Cs4)""",
    longDesc = 
"""

""",
)

entry(
    index = 1331,
    label = "Cs-CbCtCtCs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cb u0 {1,S}
3   Ct u0 {1,S}
4   Ct u0 {1,S}
5   Cs u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.37,8.11,9.55,10.1,10.03,9.36,6.65],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (6.23,'kcal/mol','+|-',0.26),
        S298 = (-35.34,'cal/(mol*K)','+|-',0.13),
    ),
    shortDesc = """Cs-CbCtCtCs BOZZELLI =3D Cs/Cs2/Cb/Ct + (Cs/Cs3/Ct - Cs/Cs4)""",
    longDesc = 
"""

""",
)

entry(
    index = 1332,
    label = "Cs-CbCbCtCs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cb u0 {1,S}
3   Cb u0 {1,S}
4   Ct u0 {1,S}
5   Cs u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.37,8.11,9.55,10.1,10.03,9.36,6.65],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (6.43,'kcal/mol','+|-',0.26),
        S298 = (-35.34,'cal/(mol*K)','+|-',0.13),
    ),
    shortDesc = """Cs-CbCbCtCs BOZZELLI =3D Cs/Cs2/Cb2 + (Cs/Cs3/Ct - Cs/Cs4)""",
    longDesc = 
"""

""",
)

entry(
    index = 1333,
    label = "Cs-CbCbCbCs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cb u0 {1,S}
3   Cb u0 {1,S}
4   Cb u0 {1,S}
5   Cs u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.37,8.11,9.55,10.1,10.03,9.36,6.65],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (6.23,'kcal/mol','+|-',0.26),
        S298 = (-35.34,'cal/(mol*K)','+|-',0.13),
    ),
    shortDesc = """Cs-CbCbCbCs BOZZELLI =3D Cs/Cs2/Cb2 + (Cs/Cs3/Cb - Cs/Cs4)""",
    longDesc = 
"""

""",
)

entry(
    index = 1334,
    label = "Cs-CdsCdsCdsCds",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   [Cd,CO] u0 {1,S}
3   [Cd,CO] u0 {1,S}
4   [Cd,CO] u0 {1,S}
5   [Cd,CO] u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1335,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-O2d)(Cds-O2d)",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   CO  u0 {1,S} {7,D}
4   CO  u0 {1,S} {8,D}
5   CO  u0 {1,S} {9,D}
6   O2d u0 {2,D}
7   O2d u0 {3,D}
8   O2d u0 {4,D}
9   O2d u0 {5,D}
""",
    thermo = 'Cs-CsCsCsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1336,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-O2d)(Cds-Cd)",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   CO  u0 {1,S} {8,D}
4   CO  u0 {1,S} {9,D}
5   Cd  u0 {1,S} {6,D}
6   C   u0 {5,D}
7   O2d u0 {2,D}
8   O2d u0 {3,D}
9   O2d u0 {4,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-O2d)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1337,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-O2d)(Cds-Cds)",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   CO  u0 {1,S} {8,D}
4   CO  u0 {1,S} {9,D}
5   Cd  u0 {1,S} {6,D}
6   Cd  u0 {5,D}
7   O2d u0 {2,D}
8   O2d u0 {3,D}
9   O2d u0 {4,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-O2d)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1338,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-O2d)(Cds-Cdd)",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   CO  u0 {1,S} {8,D}
4   CO  u0 {1,S} {9,D}
5   Cd  u0 {1,S} {6,D}
6   Cdd u0 {5,D}
7   O2d u0 {2,D}
8   O2d u0 {3,D}
9   O2d u0 {4,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-O2d)(Cds-Cdd-Cd)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1339,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-O2d)(Cds-Cdd-O2d)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    CO  u0 {1,S} {7,D}
4    CO  u0 {1,S} {8,D}
5    CO  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    O2d u0 {3,D}
8    O2d u0 {4,D}
9    O2d u0 {5,D}
10   O2d u0 {6,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)CsCsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1340,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-O2d)(Cds-Cdd-Cd)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    CO  u0 {1,S} {7,D}
4    CO  u0 {1,S} {8,D}
5    CO  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    O2d u0 {3,D}
8    O2d u0 {4,D}
9    O2d u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-O2d)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1341,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-Cd)(Cds-Cd)",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {8,D}
3   CO  u0 {1,S} {9,D}
4   Cd  u0 {1,S} {6,D}
5   Cd  u0 {1,S} {7,D}
6   C   u0 {4,D}
7   C   u0 {5,D}
8   O2d u0 {2,D}
9   O2d u0 {3,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([42.49,50.96,52.27,50.54,45.33,41.1,35.7],'J/(mol*K)','+|-',[5.08,5.08,5.08,5.08,5.08,5.08,5.08]),
        H298 = (25.2,'kJ/mol','+|-',4.33),
        S298 = (-168.67,'J/(mol*K)','+|-',5.92),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 1342,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-Cds)(Cds-Cds)",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {8,D}
3   CO  u0 {1,S} {9,D}
4   Cd  u0 {1,S} {6,D}
5   Cd  u0 {1,S} {7,D}
6   Cd  u0 {4,D}
7   Cd  u0 {5,D}
8   O2d u0 {2,D}
9   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)CsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1343,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-Cdd)(Cds-Cds)",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {8,D}
3   CO  u0 {1,S} {9,D}
4   Cd  u0 {1,S} {6,D}
5   Cd  u0 {1,S} {7,D}
6   Cdd u0 {4,D}
7   Cd  u0 {5,D}
8   O2d u0 {2,D}
9   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-Cdd-Cd)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1344,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cds)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    CO  u0 {1,S} {8,D}
4    CO  u0 {1,S} {9,D}
5    Cd  u0 {1,S} {7,D}
6    Cdd u0 {2,D} {10,D}
7    Cd  u0 {5,D}
8    O2d u0 {3,D}
9    O2d u0 {4,D}
10   O2d u0 {6,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-Cdd-O2d)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1345,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-Cdd-Cd)(Cds-Cds)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    CO  u0 {1,S} {8,D}
4    CO  u0 {1,S} {9,D}
5    Cd  u0 {1,S} {7,D}
6    Cdd u0 {2,D} {10,D}
7    Cd  u0 {5,D}
8    O2d u0 {3,D}
9    O2d u0 {4,D}
10   C   u0 {6,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1346,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-Cdd)(Cds-Cdd)",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {8,D}
3   CO  u0 {1,S} {9,D}
4   Cd  u0 {1,S} {6,D}
5   Cd  u0 {1,S} {7,D}
6   Cdd u0 {4,D}
7   Cdd u0 {5,D}
8   O2d u0 {2,D}
9   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-Cdd-Cd)(Cds-Cdd-Cd)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1347,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cdd-O2d)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    CO  u0 {1,S} {8,D}
5    CO  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    O2d u0 {4,D}
9    O2d u0 {5,D}
10   O2d u0 {6,D}
11   O2d u0 {7,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)CsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1348,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cdd-Cd)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    CO  u0 {1,S} {8,D}
5    CO  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    O2d u0 {4,D}
9    O2d u0 {5,D}
10   O2d u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1349,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    CO  u0 {1,S} {8,D}
5    CO  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    O2d u0 {4,D}
9    O2d u0 {5,D}
10   C   u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1350,
    label = "Cs-(Cds-O2d)(Cds-Cd)(Cds-Cd)(Cds-Cd)",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {9,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Cd  u0 {1,S} {8,D}
6   C   u0 {3,D}
7   C   u0 {4,D}
8   C   u0 {5,D}
9   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1351,
    label = "Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)(Cds-Cds)",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {9,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Cd  u0 {1,S} {8,D}
6   Cd  u0 {3,D}
7   Cd  u0 {4,D}
8   Cd  u0 {5,D}
9   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)CsCsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1352,
    label = "Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)(Cds-Cdd)",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {9,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Cd  u0 {1,S} {8,D}
6   Cd  u0 {3,D}
7   Cd  u0 {4,D}
8   Cdd u0 {5,D}
9   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1353,
    label = "Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    CO  u0 {1,S} {9,D}
4    Cd  u0 {1,S} {7,D}
5    Cd  u0 {1,S} {8,D}
6    Cdd u0 {2,D} {10,D}
7    Cd  u0 {4,D}
8    Cd  u0 {5,D}
9    O2d u0 {3,D}
10   O2d u0 {6,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-O2d)CsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1354,
    label = "Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    CO  u0 {1,S} {9,D}
4    Cd  u0 {1,S} {7,D}
5    Cd  u0 {1,S} {8,D}
6    Cdd u0 {2,D} {10,D}
7    Cd  u0 {4,D}
8    Cd  u0 {5,D}
9    O2d u0 {3,D}
10   C   u0 {6,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1355,
    label = "Cs-(Cds-O2d)(Cds-Cds)(Cds-Cdd)(Cds-Cdd)",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {9,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Cd  u0 {1,S} {8,D}
6   Cd  u0 {3,D}
7   Cdd u0 {4,D}
8   Cdd u0 {5,D}
9   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1356,
    label = "Cs-(Cds-O2d)(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-O2d)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    CO  u0 {1,S} {9,D}
5    Cd  u0 {1,S} {8,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    Cd  u0 {5,D}
9    O2d u0 {4,D}
10   O2d u0 {6,D}
11   O2d u0 {7,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cdd-O2d)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1357,
    label = "Cs-(Cds-O2d)(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-Cd)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    CO  u0 {1,S} {9,D}
5    Cd  u0 {1,S} {8,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    Cd  u0 {5,D}
9    O2d u0 {4,D}
10   O2d u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1358,
    label = "Cs-(Cds-O2d)(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    CO  u0 {1,S} {9,D}
5    Cd  u0 {1,S} {8,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    Cd  u0 {5,D}
9    O2d u0 {4,D}
10   C   u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1359,
    label = "Cs-(Cds-O2d)(Cds-Cdd)(Cds-Cdd)(Cds-Cdd)",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {9,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Cd  u0 {1,S} {8,D}
6   Cdd u0 {3,D}
7   Cdd u0 {4,D}
8   Cdd u0 {5,D}
9   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1360,
    label = "Cs-(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cdd-O2d)(Cds-Cdd-O2d)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    Cd  u0 {1,S} {8,D}
5    CO  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    Cdd u0 {4,D} {12,D}
9    O2d u0 {5,D}
10   O2d u0 {6,D}
11   O2d u0 {7,D}
12   O2d u0 {8,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)(Cds-Cdd-O2d)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1361,
    label = "Cs-(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cdd-O2d)(Cds-Cdd-Cd)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    Cd  u0 {1,S} {8,D}
5    CO  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    Cdd u0 {4,D} {12,D}
9    O2d u0 {5,D}
10   O2d u0 {6,D}
11   O2d u0 {7,D}
12   C   u0 {8,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1362,
    label = "Cs-(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    Cd  u0 {1,S} {8,D}
5    CO  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    Cdd u0 {4,D} {12,D}
9    O2d u0 {5,D}
10   O2d u0 {6,D}
11   C   u0 {7,D}
12   C   u0 {8,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1363,
    label = "Cs-(Cds-O2d)(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    Cd  u0 {1,S} {8,D}
5    CO  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    Cdd u0 {4,D} {12,D}
9    O2d u0 {5,D}
10   C   u0 {6,D}
11   C   u0 {7,D}
12   C   u0 {8,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1364,
    label = "Cs-(Cds-Cd)(Cds-Cd)(Cds-Cd)(Cds-Cd)",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cd u0 {1,S} {7,D}
4   Cd u0 {1,S} {8,D}
5   Cd u0 {1,S} {9,D}
6   C  u0 {2,D}
7   C  u0 {3,D}
8   C  u0 {4,D}
9   C  u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1365,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cds)",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cd u0 {1,S} {7,D}
4   Cd u0 {1,S} {8,D}
5   Cd u0 {1,S} {9,D}
6   Cd u0 {2,D}
7   Cd u0 {3,D}
8   Cd u0 {4,D}
9   Cd u0 {5,D}
""",
    thermo = 'Cs-CsCsCsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1366,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cdd)",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   Cd  u0 {1,S} {9,D}
6   Cd  u0 {2,D}
7   Cd  u0 {3,D}
8   Cd  u0 {4,D}
9   Cdd u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1367,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    Cd  u0 {1,S} {8,D}
5    Cd  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cd  u0 {3,D}
8    Cd  u0 {4,D}
9    Cd  u0 {5,D}
10   O2d u0 {6,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)CsCsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1368,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cdd-S2d)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    Cd  u0 {1,S} {8,D}
5    Cd  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cd  u0 {3,D}
8    Cd  u0 {4,D}
9    Cd  u0 {5,D}
10   S2d u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1369,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    Cd  u0 {1,S} {8,D}
5    Cd  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cd  u0 {3,D}
8    Cd  u0 {4,D}
9    Cd  u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1370,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd)(Cds-Cdd)",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   Cd  u0 {1,S} {9,D}
6   Cd  u0 {2,D}
7   Cd  u0 {3,D}
8   Cdd u0 {4,D}
9   Cdd u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1371,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-O2d)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    Cd  u0 {1,S} {8,D}
5    Cd  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    Cd  u0 {4,D}
9    Cd  u0 {5,D}
10   O2d u0 {6,D}
11   O2d u0 {7,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)CsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1372,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-Cd)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    Cd  u0 {1,S} {8,D}
5    Cd  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    Cd  u0 {4,D}
9    Cd  u0 {5,D}
10   O2d u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1373,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-S2d)(Cds-Cdd-S2d)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    Cd  u0 {1,S} {8,D}
5    Cd  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    Cd  u0 {4,D}
9    Cd  u0 {5,D}
10   S2d u0 {6,D}
11   S2d u0 {7,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1374,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-S2d)(Cds-Cdd-Cd)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    Cd  u0 {1,S} {8,D}
5    Cd  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    Cd  u0 {4,D}
9    Cd  u0 {5,D}
10   S2d u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1375,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    Cd  u0 {1,S} {8,D}
5    Cd  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    Cd  u0 {4,D}
9    Cd  u0 {5,D}
10   C   u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1376,
    label = "Cs-(Cds-Cds)(Cds-Cdd)(Cds-Cdd)(Cds-Cdd)",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   Cd  u0 {1,S} {9,D}
6   Cd  u0 {2,D}
7   Cdd u0 {3,D}
8   Cdd u0 {4,D}
9   Cdd u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1377,
    label = "Cs-(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-O2d)(Cds-Cdd-O2d)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    Cd  u0 {1,S} {8,D}
5    Cd  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    Cdd u0 {4,D} {12,D}
9    Cd  u0 {5,D}
10   O2d u0 {6,D}
11   O2d u0 {7,D}
12   O2d u0 {8,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)(Cds-Cdd-O2d)Cs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1378,
    label = "Cs-(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-O2d)(Cds-Cdd-Cd)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    Cd  u0 {1,S} {8,D}
5    Cd  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    Cdd u0 {4,D} {12,D}
9    Cd  u0 {5,D}
10   O2d u0 {6,D}
11   O2d u0 {7,D}
12   C   u0 {8,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1379,
    label = "Cs-(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    Cd  u0 {1,S} {8,D}
5    Cd  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    Cdd u0 {4,D} {12,D}
9    Cd  u0 {5,D}
10   O2d u0 {6,D}
11   C   u0 {7,D}
12   C   u0 {8,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1380,
    label = "Cs-(Cds-Cds)(Cds-Cdd-S2d)(Cds-Cdd-S2d)(Cds-Cdd-S2d)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    Cd  u0 {1,S} {8,D}
5    Cd  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    Cdd u0 {4,D} {12,D}
9    Cd  u0 {5,D}
10   S2d u0 {6,D}
11   S2d u0 {7,D}
12   S2d u0 {8,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1381,
    label = "Cs-(Cds-Cds)(Cds-Cdd-S2d)(Cds-Cdd-S2d)(Cds-Cdd-Cd)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    Cd  u0 {1,S} {8,D}
5    Cd  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    Cdd u0 {4,D} {12,D}
9    Cd  u0 {5,D}
10   S2d u0 {6,D}
11   S2d u0 {7,D}
12   C   u0 {8,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1382,
    label = "Cs-(Cds-Cds)(Cds-Cdd-S2d)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    Cd  u0 {1,S} {8,D}
5    Cd  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    Cdd u0 {4,D} {12,D}
9    Cd  u0 {5,D}
10   S2d u0 {6,D}
11   C   u0 {7,D}
12   C   u0 {8,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1383,
    label = "Cs-(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    Cd  u0 {1,S} {8,D}
5    Cd  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    Cdd u0 {4,D} {12,D}
9    Cd  u0 {5,D}
10   C   u0 {6,D}
11   C   u0 {7,D}
12   C   u0 {8,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1384,
    label = "Cs-(Cds-Cdd)(Cds-Cdd)(Cds-Cdd)(Cds-Cdd)",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   Cd  u0 {1,S} {9,D}
6   Cdd u0 {2,D}
7   Cdd u0 {3,D}
8   Cdd u0 {4,D}
9   Cdd u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1385,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)(Cds-Cdd-O2d)(Cds-Cdd-O2d)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    Cd  u0 {1,S} {8,D}
5    Cd  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    Cdd u0 {4,D} {12,D}
9    Cdd u0 {5,D} {13,D}
10   O2d u0 {6,D}
11   O2d u0 {7,D}
12   O2d u0 {8,D}
13   O2d u0 {9,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1386,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)(Cds-Cdd-O2d)(Cds-Cdd-Cd)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    Cd  u0 {1,S} {8,D}
5    Cd  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    Cdd u0 {4,D} {12,D}
9    Cdd u0 {5,D} {13,D}
10   O2d u0 {6,D}
11   O2d u0 {7,D}
12   O2d u0 {8,D}
13   C   u0 {9,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-O2d)(Cds-Cdd-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1387,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    Cd  u0 {1,S} {8,D}
5    Cd  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    Cdd u0 {4,D} {12,D}
9    Cdd u0 {5,D} {13,D}
10   O2d u0 {6,D}
11   O2d u0 {7,D}
12   C   u0 {8,D}
13   C   u0 {9,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1388,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    Cd  u0 {1,S} {8,D}
5    Cd  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    Cdd u0 {4,D} {12,D}
9    Cdd u0 {5,D} {13,D}
10   O2d u0 {6,D}
11   C   u0 {7,D}
12   C   u0 {8,D}
13   C   u0 {9,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1389,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-S2d)(Cds-Cdd-S2d)(Cds-Cdd-S2d)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    Cd  u0 {1,S} {8,D}
5    Cd  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    Cdd u0 {4,D} {12,D}
9    Cdd u0 {5,D} {13,D}
10   S2d u0 {6,D}
11   S2d u0 {7,D}
12   S2d u0 {8,D}
13   S2d u0 {9,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1390,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-S2d)(Cds-Cdd-S2d)(Cds-Cdd-Cd)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    Cd  u0 {1,S} {8,D}
5    Cd  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    Cdd u0 {4,D} {12,D}
9    Cdd u0 {5,D} {13,D}
10   S2d u0 {6,D}
11   S2d u0 {7,D}
12   S2d u0 {8,D}
13   C   u0 {9,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1391,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-S2d)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    Cd  u0 {1,S} {8,D}
5    Cd  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    Cdd u0 {4,D} {12,D}
9    Cdd u0 {5,D} {13,D}
10   S2d u0 {6,D}
11   S2d u0 {7,D}
12   C   u0 {8,D}
13   C   u0 {9,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1392,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    Cd  u0 {1,S} {8,D}
5    Cd  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    Cdd u0 {4,D} {12,D}
9    Cdd u0 {5,D} {13,D}
10   S2d u0 {6,D}
11   C   u0 {7,D}
12   C   u0 {8,D}
13   C   u0 {9,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1393,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    Cd  u0 {1,S} {8,D}
5    Cd  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    Cdd u0 {4,D} {12,D}
9    Cdd u0 {5,D} {13,D}
10   C   u0 {6,D}
11   C   u0 {7,D}
12   C   u0 {8,D}
13   C   u0 {9,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1394,
    label = "Cs-CtCdsCdsCds",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Ct      u0 {1,S}
3   [Cd,CO] u0 {1,S}
4   [Cd,CO] u0 {1,S}
5   [Cd,CO] u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1395,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-O2d)Ct",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   CO  u0 {1,S} {7,D}
4   CO  u0 {1,S} {8,D}
5   Ct  u0 {1,S}
6   O2d u0 {2,D}
7   O2d u0 {3,D}
8   O2d u0 {4,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-O2d)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1396,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-Cd)Ct",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   CO  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {6,D}
5   Ct  u0 {1,S}
6   C   u0 {4,D}
7   O2d u0 {2,D}
8   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-Cds)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1397,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-Cds)Ct",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   CO  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {6,D}
5   Ct  u0 {1,S}
6   Cd  u0 {4,D}
7   O2d u0 {2,D}
8   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1398,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-Cdd)Ct",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   CO  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {6,D}
5   Ct  u0 {1,S}
6   Cdd u0 {4,D}
7   O2d u0 {2,D}
8   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-Cdd-Cd)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1399,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-Cdd-O2d)Ct",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CO  u0 {1,S} {7,D}
4   CO  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   Ct  u0 {1,S}
7   O2d u0 {3,D}
8   O2d u0 {4,D}
9   O2d u0 {5,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1400,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-Cdd-Cd)Ct",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CO  u0 {1,S} {7,D}
4   CO  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   Ct  u0 {1,S}
7   O2d u0 {3,D}
8   O2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-Cds)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1401,
    label = "Cs-(Cds-O2d)(Cds-Cd)(Cds-Cd)Ct",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Ct  u0 {1,S}
6   C   u0 {3,D}
7   C   u0 {4,D}
8   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1402,
    label = "Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)Ct",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Ct  u0 {1,S}
6   Cd  u0 {3,D}
7   Cd  u0 {4,D}
8   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1403,
    label = "Cs-(Cds-O2d)(Cds-Cdd)(Cds-Cds)Ct",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Ct  u0 {1,S}
6   Cdd u0 {3,D}
7   Cd  u0 {4,D}
8   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-Cd)(Cds-Cds)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1404,
    label = "Cs-(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cds)Ct",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CO  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {7,D}
5   Cdd u0 {2,D} {9,D}
6   Ct  u0 {1,S}
7   Cd  u0 {4,D}
8   O2d u0 {3,D}
9   O2d u0 {5,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1405,
    label = "Cs-(Cds-O2d)(Cds-Cdd-Cd)(Cds-Cds)Ct",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CO  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {7,D}
5   Cdd u0 {2,D} {9,D}
6   Ct  u0 {1,S}
7   Cd  u0 {4,D}
8   O2d u0 {3,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1406,
    label = "Cs-(Cds-O2d)(Cds-Cdd)(Cds-Cdd)Ct",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Ct  u0 {1,S}
6   Cdd u0 {3,D}
7   Cdd u0 {4,D}
8   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1407,
    label = "Cs-(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cdd-O2d)Ct",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    CO  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Ct  u0 {1,S}
8    O2d u0 {4,D}
9    O2d u0 {5,D}
10   O2d u0 {6,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1408,
    label = "Cs-(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cdd-Cd)Ct",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    CO  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Ct  u0 {1,S}
8    O2d u0 {4,D}
9    O2d u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cds)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1409,
    label = "Cs-(Cds-O2d)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Ct",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    CO  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Ct  u0 {1,S}
8    O2d u0 {4,D}
9    C   u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1410,
    label = "Cs-(Cds-Cd)(Cds-Cd)(Cds-Cd)Ct",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cd u0 {1,S} {7,D}
4   Cd u0 {1,S} {8,D}
5   Ct u0 {1,S}
6   C  u0 {2,D}
7   C  u0 {3,D}
8   C  u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1411,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Ct",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cd u0 {1,S} {7,D}
4   Cd u0 {1,S} {8,D}
5   Ct u0 {1,S}
6   Cd u0 {2,D}
7   Cd u0 {3,D}
8   Cd u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1412,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd)Ct",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   Ct  u0 {1,S}
6   Cd  u0 {2,D}
7   Cd  u0 {3,D}
8   Cdd u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1413,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)Ct",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   Ct  u0 {1,S}
7   Cd  u0 {3,D}
8   Cd  u0 {4,D}
9   O2d u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1414,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-S2d)Ct",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   Ct  u0 {1,S}
7   Cd  u0 {3,D}
8   Cd  u0 {4,D}
9   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1415,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)Ct",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   Ct  u0 {1,S}
7   Cd  u0 {3,D}
8   Cd  u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1416,
    label = "Cs-(Cds-Cds)(Cds-Cdd)(Cds-Cdd)Ct",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   Ct  u0 {1,S}
6   Cd  u0 {2,D}
7   Cdd u0 {3,D}
8   Cdd u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1417,
    label = "Cs-(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-O2d)Ct",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Ct  u0 {1,S}
8    Cd  u0 {4,D}
9    O2d u0 {5,D}
10   O2d u0 {6,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1418,
    label = "Cs-(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-Cd)Ct",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Ct  u0 {1,S}
8    Cd  u0 {4,D}
9    O2d u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1419,
    label = "Cs-(Cds-Cds)(Cds-Cdd-S2d)(Cds-Cdd-S2d)Ct",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Ct  u0 {1,S}
8    Cd  u0 {4,D}
9    S2d u0 {5,D}
10   S2d u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1420,
    label = "Cs-(Cds-Cds)(Cds-Cdd-S2d)(Cds-Cdd-Cd)Ct",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Ct  u0 {1,S}
8    Cd  u0 {4,D}
9    S2d u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1421,
    label = "Cs-(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Ct",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Ct  u0 {1,S}
8    Cd  u0 {4,D}
9    C   u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1422,
    label = "Cs-(Cds-Cdd)(Cds-Cdd)(Cds-Cdd)Ct",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   Ct  u0 {1,S}
6   Cdd u0 {2,D}
7   Cdd u0 {3,D}
8   Cdd u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1423,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)(Cds-Cdd-O2d)Ct",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    Ct  u0 {1,S}
9    O2d u0 {5,D}
10   O2d u0 {6,D}
11   O2d u0 {7,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-O2d)(Cds-Cdd-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1424,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)(Cds-Cdd-Cd)Ct",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    Ct  u0 {1,S}
9    O2d u0 {5,D}
10   O2d u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-O2d)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1425,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Ct",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    Ct  u0 {1,S}
9    O2d u0 {5,D}
10   C   u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1426,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-S2d)(Cds-Cdd-S2d)Ct",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    Ct  u0 {1,S}
9    S2d u0 {5,D}
10   S2d u0 {6,D}
11   S2d u0 {7,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1427,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-S2d)(Cds-Cdd-Cd)Ct",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    Ct  u0 {1,S}
9    S2d u0 {5,D}
10   S2d u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1428,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Ct",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    Ct  u0 {1,S}
9    S2d u0 {5,D}
10   C   u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1429,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Ct",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    Ct  u0 {1,S}
9    C   u0 {5,D}
10   C   u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1430,
    label = "Cs-CbCdsCdsCds",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Cb      u0 {1,S}
3   [Cd,CO] u0 {1,S}
4   [Cd,CO] u0 {1,S}
5   [Cd,CO] u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Cb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1431,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-O2d)Cb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   CO  u0 {1,S} {7,D}
4   CO  u0 {1,S} {8,D}
5   Cb  u0 {1,S}
6   O2d u0 {2,D}
7   O2d u0 {3,D}
8   O2d u0 {4,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-O2d)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1432,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-Cd)Cb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   CO  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {6,D}
5   Cb  u0 {1,S}
6   C   u0 {4,D}
7   O2d u0 {2,D}
8   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-Cds)Cb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1433,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-Cds)Cb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   CO  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {6,D}
5   Cb  u0 {1,S}
6   Cd  u0 {4,D}
7   O2d u0 {2,D}
8   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1434,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-Cdd)Cb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   CO  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {6,D}
5   Cb  u0 {1,S}
6   Cdd u0 {4,D}
7   O2d u0 {2,D}
8   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-Cdd-Cd)Cb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1435,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-Cdd-O2d)Cb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CO  u0 {1,S} {7,D}
4   CO  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   Cb  u0 {1,S}
7   O2d u0 {3,D}
8   O2d u0 {4,D}
9   O2d u0 {5,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1436,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-Cdd-Cd)Cb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CO  u0 {1,S} {7,D}
4   CO  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   Cb  u0 {1,S}
7   O2d u0 {3,D}
8   O2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-Cds)Cb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1437,
    label = "Cs-(Cds-O2d)(Cds-Cd)(Cds-Cd)Cb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Cb  u0 {1,S}
6   C   u0 {3,D}
7   C   u0 {4,D}
8   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)Cb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1438,
    label = "Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)Cb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Cb  u0 {1,S}
6   Cd  u0 {3,D}
7   Cd  u0 {4,D}
8   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1439,
    label = "Cs-(Cds-O2d)(Cds-Cdd)(Cds-Cds)Cb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Cb  u0 {1,S}
6   Cdd u0 {3,D}
7   Cd  u0 {4,D}
8   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-Cd)(Cds-Cds)Cb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1440,
    label = "Cs-(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cds)Cb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CO  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {7,D}
5   Cdd u0 {2,D} {9,D}
6   Cb  u0 {1,S}
7   Cd  u0 {4,D}
8   O2d u0 {3,D}
9   O2d u0 {5,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1441,
    label = "Cs-(Cds-O2d)(Cds-Cdd-Cd)(Cds-Cds)Cb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CO  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {7,D}
5   Cdd u0 {2,D} {9,D}
6   Cb  u0 {1,S}
7   Cd  u0 {4,D}
8   O2d u0 {3,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)Cb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1442,
    label = "Cs-(Cds-O2d)(Cds-Cdd)(Cds-Cdd)Cb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Cb  u0 {1,S}
6   Cdd u0 {3,D}
7   Cdd u0 {4,D}
8   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1443,
    label = "Cs-(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cdd-O2d)Cb",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    CO  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cb  u0 {1,S}
8    O2d u0 {4,D}
9    O2d u0 {5,D}
10   O2d u0 {6,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1444,
    label = "Cs-(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cdd-Cd)Cb",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    CO  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cb  u0 {1,S}
8    O2d u0 {4,D}
9    O2d u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cds)Cb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1445,
    label = "Cs-(Cds-O2d)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cb",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    CO  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cb  u0 {1,S}
8    O2d u0 {4,D}
9    C   u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)Cb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1446,
    label = "Cs-(Cds-Cd)(Cds-Cd)(Cds-Cd)Cb",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cd u0 {1,S} {7,D}
4   Cd u0 {1,S} {8,D}
5   Cb u0 {1,S}
6   C  u0 {2,D}
7   C  u0 {3,D}
8   C  u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Cb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1447,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Cb",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cd u0 {1,S} {7,D}
4   Cd u0 {1,S} {8,D}
5   Cb u0 {1,S}
6   Cd u0 {2,D}
7   Cd u0 {3,D}
8   Cd u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1448,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd)Cb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   Cb  u0 {1,S}
6   Cd  u0 {2,D}
7   Cd  u0 {3,D}
8   Cdd u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)Cb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1449,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)Cb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   Cb  u0 {1,S}
7   Cd  u0 {3,D}
8   Cd  u0 {4,D}
9   O2d u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1450,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-S2d)Cb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   Cb  u0 {1,S}
7   Cd  u0 {3,D}
8   Cd  u0 {4,D}
9   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1451,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)Cb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   Cb  u0 {1,S}
7   Cd  u0 {3,D}
8   Cd  u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Cb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1452,
    label = "Cs-(Cds-Cds)(Cds-Cdd)(Cds-Cdd)Cb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   Cb  u0 {1,S}
6   Cd  u0 {2,D}
7   Cdd u0 {3,D}
8   Cdd u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1453,
    label = "Cs-(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-O2d)Cb",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cb  u0 {1,S}
8    Cd  u0 {4,D}
9    O2d u0 {5,D}
10   O2d u0 {6,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1454,
    label = "Cs-(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-Cd)Cb",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cb  u0 {1,S}
8    Cd  u0 {4,D}
9    O2d u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)Cb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1455,
    label = "Cs-(Cds-Cds)(Cds-Cdd-S2d)(Cds-Cdd-S2d)Cb",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cb  u0 {1,S}
8    Cd  u0 {4,D}
9    S2d u0 {5,D}
10   S2d u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1456,
    label = "Cs-(Cds-Cds)(Cds-Cdd-S2d)(Cds-Cdd-Cd)Cb",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cb  u0 {1,S}
8    Cd  u0 {4,D}
9    S2d u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1457,
    label = "Cs-(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cb",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cb  u0 {1,S}
8    Cd  u0 {4,D}
9    C   u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Cb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1458,
    label = "Cs-(Cds-Cdd)(Cds-Cdd)(Cds-Cdd)Cb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   Cb  u0 {1,S}
6   Cdd u0 {2,D}
7   Cdd u0 {3,D}
8   Cdd u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1459,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)(Cds-Cdd-O2d)Cb",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    Cb  u0 {1,S}
9    O2d u0 {5,D}
10   O2d u0 {6,D}
11   O2d u0 {7,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-O2d)(Cds-Cdd-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1460,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)(Cds-Cdd-Cd)Cb",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    Cb  u0 {1,S}
9    O2d u0 {5,D}
10   O2d u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-O2d)Cb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1461,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cb",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    Cb  u0 {1,S}
9    O2d u0 {5,D}
10   C   u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)Cb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1462,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-S2d)(Cds-Cdd-S2d)Cb",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    Cb  u0 {1,S}
9    S2d u0 {5,D}
10   S2d u0 {6,D}
11   S2d u0 {7,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1463,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-S2d)(Cds-Cdd-Cd)Cb",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    Cb  u0 {1,S}
9    S2d u0 {5,D}
10   S2d u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1464,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cb",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    Cb  u0 {1,S}
9    S2d u0 {5,D}
10   C   u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1465,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cb",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    Cb  u0 {1,S}
9    C   u0 {5,D}
10   C   u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Cb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1466,
    label = "Cs-CtCtCdsCds",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Ct      u0 {1,S}
3   Ct      u0 {1,S}
4   [Cd,CO] u0 {1,S}
5   [Cd,CO] u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CtCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1467,
    label = "Cs-(Cds-O2d)(Cds-O2d)CtCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   CO  u0 {1,S} {7,D}
4   Ct  u0 {1,S}
5   Ct  u0 {1,S}
6   O2d u0 {2,D}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1468,
    label = "Cs-(Cds-O2d)(Cds-Cd)CtCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Ct  u0 {1,S}
5   Ct  u0 {1,S}
6   C   u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)CtCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1469,
    label = "Cs-(Cds-O2d)(Cds-Cds)CtCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Ct  u0 {1,S}
5   Ct  u0 {1,S}
6   Cd  u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1470,
    label = "Cs-(Cds-O2d)(Cds-Cdd)CtCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Ct  u0 {1,S}
5   Ct  u0 {1,S}
6   Cdd u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-Cd)CtCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1471,
    label = "Cs-(Cds-O2d)(Cds-Cdd-O2d)CtCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CO  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Ct  u0 {1,S}
6   Ct  u0 {1,S}
7   O2d u0 {3,D}
8   O2d u0 {4,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1472,
    label = "Cs-(Cds-O2d)(Cds-Cdd-Cd)CtCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CO  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Ct  u0 {1,S}
6   Ct  u0 {1,S}
7   O2d u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)CtCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1473,
    label = "Cs-(Cds-Cd)(Cds-Cd)CtCt",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cd u0 {1,S} {7,D}
4   Ct u0 {1,S}
5   Ct u0 {1,S}
6   C  u0 {2,D}
7   C  u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CtCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1474,
    label = "Cs-(Cds-Cds)(Cds-Cds)CtCt",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cd u0 {1,S} {7,D}
4   Ct u0 {1,S}
5   Ct u0 {1,S}
6   Cd u0 {2,D}
7   Cd u0 {3,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.61,7.3,8.97,9.69,9.84,9.42,7.36],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (5.48,'kcal/mol','+|-',0.26),
        S298 = (-34.5,'cal/(mol*K)','+|-',0.13),
    ),
    shortDesc = """Cs-CtCtCdCd BOZZELLI =3D Cs/Cs/Cd/Ct2 + (Cs/Cs3/Cd - Cs/Cs4)""",
    longDesc = 
"""

""",
)

entry(
    index = 1475,
    label = "Cs-(Cds-Cdd)(Cds-Cds)CtCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Ct  u0 {1,S}
5   Ct  u0 {1,S}
6   Cdd u0 {2,D}
7   Cd  u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cds)CtCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1476,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cds)CtCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Ct  u0 {1,S}
6   Ct  u0 {1,S}
7   Cd  u0 {3,D}
8   O2d u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1477,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cds)CtCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Ct  u0 {1,S}
6   Ct  u0 {1,S}
7   Cd  u0 {3,D}
8   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1478,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cds)CtCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Ct  u0 {1,S}
6   Ct  u0 {1,S}
7   Cd  u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CtCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1479,
    label = "Cs-(Cds-Cdd)(Cds-Cdd)CtCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Ct  u0 {1,S}
5   Ct  u0 {1,S}
6   Cdd u0 {2,D}
7   Cdd u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CtCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1480,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)CtCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Ct  u0 {1,S}
7   Ct  u0 {1,S}
8   O2d u0 {4,D}
9   O2d u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1481,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-Cd)CtCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Ct  u0 {1,S}
7   Ct  u0 {1,S}
8   O2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)(Cds-Cds)CtCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1482,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-S2d)CtCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Ct  u0 {1,S}
7   Ct  u0 {1,S}
8   S2d u0 {4,D}
9   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1483,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-Cd)CtCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Ct  u0 {1,S}
7   Ct  u0 {1,S}
8   S2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1484,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CtCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Ct  u0 {1,S}
7   Ct  u0 {1,S}
8   C   u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CtCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1485,
    label = "Cs-CbCtCdsCds",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Cb      u0 {1,S}
3   Ct      u0 {1,S}
4   [Cd,CO] u0 {1,S}
5   [Cd,CO] u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CbCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1486,
    label = "Cs-(Cds-O2d)(Cds-O2d)CbCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   CO  u0 {1,S} {7,D}
4   Cb  u0 {1,S}
5   Ct  u0 {1,S}
6   O2d u0 {2,D}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-Cds)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1487,
    label = "Cs-(Cds-O2d)(Cds-Cd)CbCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   Ct  u0 {1,S}
6   C   u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)CbCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1488,
    label = "Cs-(Cds-O2d)(Cds-Cds)CbCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   Ct  u0 {1,S}
6   Cd  u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1489,
    label = "Cs-(Cds-O2d)(Cds-Cdd)CbCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   Ct  u0 {1,S}
6   Cdd u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-Cd)CbCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1490,
    label = "Cs-(Cds-O2d)(Cds-Cdd-O2d)CbCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CO  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   Ct  u0 {1,S}
7   O2d u0 {3,D}
8   O2d u0 {4,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cds)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1491,
    label = "Cs-(Cds-O2d)(Cds-Cdd-Cd)CbCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CO  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   Ct  u0 {1,S}
7   O2d u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)CbCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1492,
    label = "Cs-(Cds-Cd)(Cds-Cd)CbCt",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cd u0 {1,S} {7,D}
4   Cb u0 {1,S}
5   Ct u0 {1,S}
6   C  u0 {2,D}
7   C  u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CbCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1493,
    label = "Cs-(Cds-Cds)(Cds-Cds)CbCt",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cd u0 {1,S} {7,D}
4   Cb u0 {1,S}
5   Ct u0 {1,S}
6   Cd u0 {2,D}
7   Cd u0 {3,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.61,7.3,8.97,9.69,9.84,9.42,7.36],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (5.48,'kcal/mol','+|-',0.26),
        S298 = (-34.5,'cal/(mol*K)','+|-',0.13),
    ),
    shortDesc = """Cs-CbCtCdCd BOZZELLI =3D Cs/Cs/Cb/Cd2 + (Cs/Cs3/Ct - Cs/Cs4)""",
    longDesc = 
"""

""",
)

entry(
    index = 1494,
    label = "Cs-(Cds-Cdd)(Cds-Cds)CbCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cb  u0 {1,S}
5   Ct  u0 {1,S}
6   Cdd u0 {2,D}
7   Cd  u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cds)CbCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1495,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cds)CbCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   Ct  u0 {1,S}
7   Cd  u0 {3,D}
8   O2d u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1496,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cds)CbCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   Ct  u0 {1,S}
7   Cd  u0 {3,D}
8   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1497,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cds)CbCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   Ct  u0 {1,S}
7   Cd  u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CbCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1498,
    label = "Cs-(Cds-Cdd)(Cds-Cdd)CbCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cb  u0 {1,S}
5   Ct  u0 {1,S}
6   Cdd u0 {2,D}
7   Cdd u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CbCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1499,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)CbCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cb  u0 {1,S}
7   Ct  u0 {1,S}
8   O2d u0 {4,D}
9   O2d u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-O2d)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1500,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-Cd)CbCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cb  u0 {1,S}
7   Ct  u0 {1,S}
8   O2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)(Cds-Cds)CbCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1501,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-S2d)CbCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cb  u0 {1,S}
7   Ct  u0 {1,S}
8   S2d u0 {4,D}
9   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1502,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-Cd)CbCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cb  u0 {1,S}
7   Ct  u0 {1,S}
8   S2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1503,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CbCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cb  u0 {1,S}
7   Ct  u0 {1,S}
8   C   u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CbCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1504,
    label = "Cs-CbCbCdsCds",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Cb      u0 {1,S}
3   Cb      u0 {1,S}
4   [Cd,CO] u0 {1,S}
5   [Cd,CO] u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CbCb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1505,
    label = "Cs-(Cds-O2d)(Cds-O2d)CbCb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   CO  u0 {1,S} {7,D}
4   Cb  u0 {1,S}
5   Cb  u0 {1,S}
6   O2d u0 {2,D}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1506,
    label = "Cs-(Cds-O2d)(Cds-Cd)CbCb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   Cb  u0 {1,S}
6   C   u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)CbCb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1507,
    label = "Cs-(Cds-O2d)(Cds-Cds)CbCb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   Cb  u0 {1,S}
6   Cd  u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1508,
    label = "Cs-(Cds-O2d)(Cds-Cdd)CbCb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   Cb  u0 {1,S}
6   Cdd u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-Cd)CbCb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1509,
    label = "Cs-(Cds-O2d)(Cds-Cdd-O2d)CbCb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CO  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   Cb  u0 {1,S}
7   O2d u0 {3,D}
8   O2d u0 {4,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1510,
    label = "Cs-(Cds-O2d)(Cds-Cdd-Cd)CbCb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CO  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   Cb  u0 {1,S}
7   O2d u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)CbCb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1511,
    label = "Cs-(Cds-Cd)(Cds-Cd)CbCb",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cd u0 {1,S} {7,D}
4   Cb u0 {1,S}
5   Cb u0 {1,S}
6   C  u0 {2,D}
7   C  u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CbCb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1512,
    label = "Cs-(Cds-Cds)(Cds-Cds)CbCb",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cd u0 {1,S} {7,D}
4   Cb u0 {1,S}
5   Cb u0 {1,S}
6   Cd u0 {2,D}
7   Cd u0 {3,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.61,7.3,8.97,9.69,9.84,9.42,7.36],'cal/(mol*K)','+|-',[0.13,0.13,0.13,0.13,0.13,0.13,0.13]),
        H298 = (5.48,'kcal/mol','+|-',0.26),
        S298 = (-34.5,'cal/(mol*K)','+|-',0.13),
    ),
    shortDesc = """Cs-CbCbCdCd BOZZELLI =3D Cs/Cs/Cb2/Cd + (Cs/Cs3/Cd - Cs/Cs4)""",
    longDesc = 
"""

""",
)

entry(
    index = 1513,
    label = "Cs-(Cds-Cdd)(Cds-Cds)CbCb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cb  u0 {1,S}
5   Cb  u0 {1,S}
6   Cdd u0 {2,D}
7   Cd  u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cds)CbCb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1514,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cds)CbCb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   Cb  u0 {1,S}
7   Cd  u0 {3,D}
8   O2d u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1515,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cds)CbCb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   Cb  u0 {1,S}
7   Cd  u0 {3,D}
8   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1516,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cds)CbCb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   Cb  u0 {1,S}
7   Cd  u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CbCb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1517,
    label = "Cs-(Cds-Cdd)(Cds-Cdd)CbCb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cb  u0 {1,S}
5   Cb  u0 {1,S}
6   Cdd u0 {2,D}
7   Cdd u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CbCb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1518,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)CbCb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cb  u0 {1,S}
7   Cb  u0 {1,S}
8   O2d u0 {4,D}
9   O2d u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1519,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-Cd)CbCb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cb  u0 {1,S}
7   Cb  u0 {1,S}
8   O2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)(Cds-Cds)CbCb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1520,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-S2d)CbCb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cb  u0 {1,S}
7   Cb  u0 {1,S}
8   S2d u0 {4,D}
9   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1521,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-Cd)CbCb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cb  u0 {1,S}
7   Cb  u0 {1,S}
8   S2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1522,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CbCb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cb  u0 {1,S}
7   Cb  u0 {1,S}
8   C   u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CbCb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1523,
    label = "Cs-CtCtCtCds",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Ct      u0 {1,S}
3   Ct      u0 {1,S}
4   Ct      u0 {1,S}
5   [Cd,CO] u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CtCtCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1524,
    label = "Cs-(Cds-O2d)CtCtCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   Ct  u0 {1,S}
5   Ct  u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1525,
    label = "Cs-(Cds-Cd)CtCtCt",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Ct u0 {1,S}
4   Ct u0 {1,S}
5   Ct u0 {1,S}
6   C  u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1526,
    label = "Cs-(Cds-Cds)CtCtCt",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Ct u0 {1,S}
4   Ct u0 {1,S}
5   Ct u0 {1,S}
6   Cd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1527,
    label = "Cs-(Cds-Cdd)CtCtCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   Ct  u0 {1,S}
5   Ct  u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)CtCtCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1528,
    label = "Cs-(Cds-Cdd-O2d)CtCtCt",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Ct  u0 {1,S}
5   Ct  u0 {1,S}
6   Ct  u0 {1,S}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1529,
    label = "Cs-(Cds-Cdd-S2d)CtCtCt",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Ct  u0 {1,S}
5   Ct  u0 {1,S}
6   Ct  u0 {1,S}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1530,
    label = "Cs-(Cds-Cdd-Cd)CtCtCt",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Ct  u0 {1,S}
5   Ct  u0 {1,S}
6   Ct  u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)CtCtCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1531,
    label = "Cs-CbCtCtCds",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Cb      u0 {1,S}
3   Ct      u0 {1,S}
4   Ct      u0 {1,S}
5   [Cd,CO] u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CbCtCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1532,
    label = "Cs-(Cds-O2d)CbCtCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Ct  u0 {1,S}
5   Ct  u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)CtCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1533,
    label = "Cs-(Cds-Cd)CbCtCt",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cb u0 {1,S}
4   Ct u0 {1,S}
5   Ct u0 {1,S}
6   C  u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)CbCtCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1534,
    label = "Cs-(Cds-Cds)CbCtCt",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cb u0 {1,S}
4   Ct u0 {1,S}
5   Ct u0 {1,S}
6   Cd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CtCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1535,
    label = "Cs-(Cds-Cdd)CbCtCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Ct  u0 {1,S}
5   Ct  u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)CbCtCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1536,
    label = "Cs-(Cds-Cdd-O2d)CbCtCt",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Ct  u0 {1,S}
6   Ct  u0 {1,S}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)(Cds-Cds)CtCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1537,
    label = "Cs-(Cds-Cdd-S2d)CbCtCt",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Ct  u0 {1,S}
6   Ct  u0 {1,S}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1538,
    label = "Cs-(Cds-Cdd-Cd)CbCtCt",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Ct  u0 {1,S}
6   Ct  u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)CbCtCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1539,
    label = "Cs-CbCbCtCds",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Cb      u0 {1,S}
3   Cb      u0 {1,S}
4   Ct      u0 {1,S}
5   [Cd,CO] u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CbCbCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1540,
    label = "Cs-(Cds-O2d)CbCbCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Cb  u0 {1,S}
5   Ct  u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1541,
    label = "Cs-(Cds-Cd)CbCbCt",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cb u0 {1,S}
4   Cb u0 {1,S}
5   Ct u0 {1,S}
6   C  u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)CbCbCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1542,
    label = "Cs-(Cds-Cds)CbCbCt",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cb u0 {1,S}
4   Cb u0 {1,S}
5   Ct u0 {1,S}
6   Cd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1543,
    label = "Cs-(Cds-Cdd)CbCbCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Cb  u0 {1,S}
5   Ct  u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)CbCbCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1544,
    label = "Cs-(Cds-Cdd-O2d)CbCbCt",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Cb  u0 {1,S}
6   Ct  u0 {1,S}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1545,
    label = "Cs-(Cds-Cdd-S2d)CbCbCt",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Cb  u0 {1,S}
6   Ct  u0 {1,S}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1546,
    label = "Cs-(Cds-Cdd-Cd)CbCbCt",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Cb  u0 {1,S}
6   Ct  u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)CbCbCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1547,
    label = "Cs-CbCbCbCds",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Cb      u0 {1,S}
3   Cb      u0 {1,S}
4   Cb      u0 {1,S}
5   [Cd,CO] u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CbCbCb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1548,
    label = "Cs-(Cds-O2d)CbCbCb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Cb  u0 {1,S}
5   Cb  u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1549,
    label = "Cs-(Cds-Cd)CbCbCb",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cb u0 {1,S}
4   Cb u0 {1,S}
5   Cb u0 {1,S}
6   C  u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)CbCbCb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1550,
    label = "Cs-(Cds-Cds)CbCbCb",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cd u0 {1,S} {6,D}
3   Cb u0 {1,S}
4   Cb u0 {1,S}
5   Cb u0 {1,S}
6   Cd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1551,
    label = "Cs-(Cds-Cdd)CbCbCb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Cb  u0 {1,S}
5   Cb  u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)CbCbCb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1552,
    label = "Cs-(Cds-Cdd-O2d)CbCbCb",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Cb  u0 {1,S}
6   Cb  u0 {1,S}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1553,
    label = "Cs-(Cds-Cdd-S2d)CbCbCb",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Cb  u0 {1,S}
6   Cb  u0 {1,S}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1554,
    label = "Cs-(Cds-Cdd-Cd)CbCbCb",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Cb  u0 {1,S}
6   Cb  u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)CbCbCb',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1555,
    label = "Cs-CtCtCtCt",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Ct u0 {1,S}
3   Ct u0 {1,S}
4   Ct u0 {1,S}
5   Ct u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1556,
    label = "Cs-CbCtCtCt",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cb u0 {1,S}
3   Ct u0 {1,S}
4   Ct u0 {1,S}
5   Ct u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CtCtCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1557,
    label = "Cs-CbCbCtCt",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cb u0 {1,S}
3   Cb u0 {1,S}
4   Ct u0 {1,S}
5   Ct u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CtCt',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1558,
    label = "Cs-CbCbCbCt",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cb u0 {1,S}
3   Cb u0 {1,S}
4   Cb u0 {1,S}
5   Ct u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Ct',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1559,
    label = "Cs-CbCbCbCb",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cb u0 {1,S}
3   Cb u0 {1,S}
4   Cb u0 {1,S}
5   Cb u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cds)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1560,
    label = "Cs-C=SCbCtCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Ct  u0 {1,S}
5   Ct  u0 {1,S}
6   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1561,
    label = "Cs-C=S(Cds-Cd)(Cds-Cd)(Cds-Cd)",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {9,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Cd  u0 {1,S} {8,D}
6   C   u0 {3,D}
7   C   u0 {4,D}
8   C   u0 {5,D}
9   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1562,
    label = "Cs-C=S(Cds-Cds)(Cds-Cds)(Cds-Cdd)",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {9,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Cd  u0 {1,S} {8,D}
6   Cd  u0 {3,D}
7   Cd  u0 {4,D}
8   Cdd u0 {5,D}
9   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1563,
    label = "Cs-C=S(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    CS  u0 {1,S} {9,D}
4    Cd  u0 {1,S} {7,D}
5    Cd  u0 {1,S} {8,D}
6    Cdd u0 {2,D} {10,D}
7    Cd  u0 {4,D}
8    Cd  u0 {5,D}
9    S2d u0 {3,D}
10   C   u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1564,
    label = "Cs-C=S(Cds-Cds)(Cds-Cds)(Cds-Cdd-S2d)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    CS  u0 {1,S} {9,D}
4    Cd  u0 {1,S} {7,D}
5    Cd  u0 {1,S} {8,D}
6    Cdd u0 {2,D} {10,D}
7    Cd  u0 {4,D}
8    Cd  u0 {5,D}
9    S2d u0 {3,D}
10   S2d u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1565,
    label = "Cs-C=S(Cds-Cdd)(Cds-Cdd)(Cds-Cdd)",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {9,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Cd  u0 {1,S} {8,D}
6   Cdd u0 {3,D}
7   Cdd u0 {4,D}
8   Cdd u0 {5,D}
9   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1566,
    label = "Cs-C=S(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    Cd  u0 {1,S} {8,D}
5    CS  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    Cdd u0 {4,D} {12,D}
9    S2d u0 {5,D}
10   C   u0 {6,D}
11   C   u0 {7,D}
12   C   u0 {8,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1567,
    label = "Cs-C=S(Cds-Cdd-S2d)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    Cd  u0 {1,S} {8,D}
5    CS  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    Cdd u0 {4,D} {12,D}
9    S2d u0 {5,D}
10   S2d u0 {6,D}
11   C   u0 {7,D}
12   C   u0 {8,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1568,
    label = "Cs-C=S(Cds-Cdd-S2d)(Cds-Cdd-S2d)(Cds-Cdd-S2d)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    Cd  u0 {1,S} {8,D}
5    CS  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    Cdd u0 {4,D} {12,D}
9    S2d u0 {5,D}
10   S2d u0 {6,D}
11   S2d u0 {7,D}
12   S2d u0 {8,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1569,
    label = "Cs-C=S(Cds-Cdd-S2d)(Cds-Cdd-S2d)(Cds-Cdd-Cd)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    Cd  u0 {1,S} {8,D}
5    CS  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    Cdd u0 {4,D} {12,D}
9    S2d u0 {5,D}
10   S2d u0 {6,D}
11   S2d u0 {7,D}
12   C   u0 {8,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1570,
    label = "Cs-C=S(Cds-Cds)(Cds-Cds)(Cds-Cds)",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {9,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Cd  u0 {1,S} {8,D}
6   Cd  u0 {3,D}
7   Cd  u0 {4,D}
8   Cd  u0 {5,D}
9   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1571,
    label = "Cs-C=S(Cds-Cds)(Cds-Cdd)(Cds-Cdd)",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {9,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Cd  u0 {1,S} {8,D}
6   Cd  u0 {3,D}
7   Cdd u0 {4,D}
8   Cdd u0 {5,D}
9   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1572,
    label = "Cs-C=S(Cds-Cds)(Cds-Cdd-S2d)(Cds-Cdd-S2d)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    CS  u0 {1,S} {9,D}
5    Cd  u0 {1,S} {8,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    Cd  u0 {5,D}
9    S2d u0 {4,D}
10   S2d u0 {6,D}
11   S2d u0 {7,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1573,
    label = "Cs-C=S(Cds-Cds)(Cds-Cdd-S2d)(Cds-Cdd-Cd)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    CS  u0 {1,S} {9,D}
5    Cd  u0 {1,S} {8,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    Cd  u0 {5,D}
9    S2d u0 {4,D}
10   S2d u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1574,
    label = "Cs-C=S(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    CS  u0 {1,S} {9,D}
5    Cd  u0 {1,S} {8,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    Cd  u0 {5,D}
9    S2d u0 {4,D}
10   C   u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1575,
    label = "Cs-C=S(Cds-Cd)CtCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Ct  u0 {1,S}
5   Ct  u0 {1,S}
6   C   u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1576,
    label = "Cs-C=S(Cds-Cds)CtCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Ct  u0 {1,S}
5   Ct  u0 {1,S}
6   Cd  u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1577,
    label = "Cs-C=S(Cds-Cdd)CtCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Ct  u0 {1,S}
5   Ct  u0 {1,S}
6   Cdd u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1578,
    label = "Cs-C=S(Cds-Cdd-S2d)CtCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CS  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Ct  u0 {1,S}
6   Ct  u0 {1,S}
7   S2d u0 {3,D}
8   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1579,
    label = "Cs-C=S(Cds-Cdd-Cd)CtCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CS  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Ct  u0 {1,S}
6   Ct  u0 {1,S}
7   S2d u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1580,
    label = "Cs-C=S(Cds-Cd)CtCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Ct  u0 {1,S}
5   Cs  u0 {1,S}
6   C   u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1581,
    label = "Cs-C=S(Cds-Cds)CtCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Ct  u0 {1,S}
5   Cs  u0 {1,S}
6   Cd  u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1582,
    label = "Cs-C=S(Cds-Cdd)CtCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Ct  u0 {1,S}
5   Cs  u0 {1,S}
6   Cdd u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1583,
    label = "Cs-C=S(Cds-Cdd-S2d)CtCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CS  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Ct  u0 {1,S}
6   Cs  u0 {1,S}
7   S2d u0 {3,D}
8   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1584,
    label = "Cs-C=S(Cds-Cdd-Cd)CtCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CS  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Ct  u0 {1,S}
6   Cs  u0 {1,S}
7   S2d u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1585,
    label = "Cs-C=SCbCbCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Cb  u0 {1,S}
5   Ct  u0 {1,S}
6   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1586,
    label = "Cs-C=SCbCsCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Cs  u0 {1,S}
5   Cs  u0 {1,S}
6   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1587,
    label = "Cs-C=SCbCbCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Cb  u0 {1,S}
5   Cs  u0 {1,S}
6   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1588,
    label = "Cs-C=SCtCtCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   Ct  u0 {1,S}
5   Ct  u0 {1,S}
6   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1589,
    label = "Cs-C=S(Cds-Cd)(Cds-Cd)Cs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Cs  u0 {1,S}
6   C   u0 {3,D}
7   C   u0 {4,D}
8   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1590,
    label = "Cs-C=S(Cds-Cdd)(Cds-Cdd)Cs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Cs  u0 {1,S}
6   Cdd u0 {3,D}
7   Cdd u0 {4,D}
8   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1591,
    label = "Cs-C=S(Cds-Cdd-S2d)(Cds-Cdd-Cd)Cs",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    CS  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cs  u0 {1,S}
8    S2d u0 {4,D}
9    S2d u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1592,
    label = "Cs-C=S(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cs",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    CS  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cs  u0 {1,S}
8    S2d u0 {4,D}
9    C   u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1593,
    label = "Cs-C=S(Cds-Cdd-S2d)(Cds-Cdd-S2d)Cs",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    CS  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cs  u0 {1,S}
8    S2d u0 {4,D}
9    S2d u0 {5,D}
10   S2d u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1594,
    label = "Cs-C=S(Cds-Cds)(Cds-Cds)Cs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Cs  u0 {1,S}
6   Cd  u0 {3,D}
7   Cd  u0 {4,D}
8   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1595,
    label = "Cs-C=S(Cds-Cdd)(Cds-Cds)Cs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Cs  u0 {1,S}
6   Cdd u0 {3,D}
7   Cd  u0 {4,D}
8   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1596,
    label = "Cs-C=S(Cds-Cdd-S2d)(Cds-Cds)Cs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CS  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {7,D}
5   Cdd u0 {2,D} {9,D}
6   Cs  u0 {1,S}
7   Cd  u0 {4,D}
8   S2d u0 {3,D}
9   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1597,
    label = "Cs-C=S(Cds-Cdd-Cd)(Cds-Cds)Cs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CS  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {7,D}
5   Cdd u0 {2,D} {9,D}
6   Cs  u0 {1,S}
7   Cd  u0 {4,D}
8   S2d u0 {3,D}
9   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1598,
    label = "Cs-C=SC=SCtCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   CS  u0 {1,S} {7,D}
4   Ct  u0 {1,S}
5   Ct  u0 {1,S}
6   S2d u0 {2,D}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1599,
    label = "Cs-C=SCsCsCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   Cs  u0 {1,S}
5   Cs  u0 {1,S}
6   S2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([2.49,5.54,7.04,7.55,7.81,7.65,7.27],'cal/(mol*K)'),
        H298 = (4.38,'kcal/mol'),
        S298 = (-33.21,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 1600,
    label = "Cs-C=SCtCtCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   Ct  u0 {1,S}
5   Cs  u0 {1,S}
6   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1601,
    label = "Cs-C=SC=SC=SCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   CS  u0 {1,S} {7,D}
4   CS  u0 {1,S} {8,D}
5   Ct  u0 {1,S}
6   S2d u0 {2,D}
7   S2d u0 {3,D}
8   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1602,
    label = "Cs-C=SC=SC=SCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   CS  u0 {1,S} {7,D}
4   CS  u0 {1,S} {8,D}
5   Cs  u0 {1,S}
6   S2d u0 {2,D}
7   S2d u0 {3,D}
8   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1603,
    label = "Cs-C=SC=SC=SC=S",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   CS  u0 {1,S} {7,D}
4   CS  u0 {1,S} {8,D}
5   CS  u0 {1,S} {9,D}
6   S2d u0 {2,D}
7   S2d u0 {3,D}
8   S2d u0 {4,D}
9   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1604,
    label = "Cs-C=SCtCsCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   Cs  u0 {1,S}
5   Cs  u0 {1,S}
6   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1605,
    label = "Cs-C=SC=SC=SCb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   CS  u0 {1,S} {7,D}
4   CS  u0 {1,S} {8,D}
5   Cb  u0 {1,S}
6   S2d u0 {2,D}
7   S2d u0 {3,D}
8   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1606,
    label = "Cs-C=SC=SC=S(Cds-Cd)",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   CS  u0 {1,S} {8,D}
4   CS  u0 {1,S} {9,D}
5   Cd  u0 {1,S} {6,D}
6   C   u0 {5,D}
7   S2d u0 {2,D}
8   S2d u0 {3,D}
9   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1607,
    label = "Cs-C=SC=SC=S(Cds-Cdd)",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   CS  u0 {1,S} {8,D}
4   CS  u0 {1,S} {9,D}
5   Cd  u0 {1,S} {6,D}
6   Cdd u0 {5,D}
7   S2d u0 {2,D}
8   S2d u0 {3,D}
9   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1608,
    label = "Cs-C=SC=SC=S(Cds-Cdd-Cd)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    CS  u0 {1,S} {7,D}
4    CS  u0 {1,S} {8,D}
5    CS  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    S2d u0 {3,D}
8    S2d u0 {4,D}
9    S2d u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1609,
    label = "Cs-C=SC=SC=S(Cds-Cdd-S2d)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    CS  u0 {1,S} {7,D}
4    CS  u0 {1,S} {8,D}
5    CS  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    S2d u0 {3,D}
8    S2d u0 {4,D}
9    S2d u0 {5,D}
10   S2d u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1610,
    label = "Cs-C=SC=SC=S(Cds-Cds)",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   CS  u0 {1,S} {8,D}
4   CS  u0 {1,S} {9,D}
5   Cd  u0 {1,S} {6,D}
6   Cd  u0 {5,D}
7   S2d u0 {2,D}
8   S2d u0 {3,D}
9   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1611,
    label = "Cs-C=S(Cds-Cd)(Cds-Cd)Ct",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Ct  u0 {1,S}
6   C   u0 {3,D}
7   C   u0 {4,D}
8   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1612,
    label = "Cs-C=S(Cds-Cdd)(Cds-Cdd)Ct",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Ct  u0 {1,S}
6   Cdd u0 {3,D}
7   Cdd u0 {4,D}
8   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1613,
    label = "Cs-C=S(Cds-Cdd-Cd)(Cds-Cdd-Cd)Ct",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    CS  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Ct  u0 {1,S}
8    S2d u0 {4,D}
9    C   u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1614,
    label = "Cs-C=S(Cds-Cdd-S2d)(Cds-Cdd-S2d)Ct",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    CS  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Ct  u0 {1,S}
8    S2d u0 {4,D}
9    S2d u0 {5,D}
10   S2d u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1615,
    label = "Cs-C=S(Cds-Cdd-S2d)(Cds-Cdd-Cd)Ct",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    CS  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Ct  u0 {1,S}
8    S2d u0 {4,D}
9    S2d u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1616,
    label = "Cs-C=S(Cds-Cds)(Cds-Cds)Ct",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Ct  u0 {1,S}
6   Cd  u0 {3,D}
7   Cd  u0 {4,D}
8   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1617,
    label = "Cs-C=S(Cds-Cdd)(Cds-Cds)Ct",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Ct  u0 {1,S}
6   Cdd u0 {3,D}
7   Cd  u0 {4,D}
8   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1618,
    label = "Cs-C=S(Cds-Cdd-Cd)(Cds-Cds)Ct",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CS  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {7,D}
5   Cdd u0 {2,D} {9,D}
6   Ct  u0 {1,S}
7   Cd  u0 {4,D}
8   S2d u0 {3,D}
9   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1619,
    label = "Cs-C=S(Cds-Cdd-S2d)(Cds-Cds)Ct",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CS  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {7,D}
5   Cdd u0 {2,D} {9,D}
6   Ct  u0 {1,S}
7   Cd  u0 {4,D}
8   S2d u0 {3,D}
9   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1620,
    label = "Cs-C=SC=SCtCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   CS  u0 {1,S} {7,D}
4   Ct  u0 {1,S}
5   Cs  u0 {1,S}
6   S2d u0 {2,D}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1621,
    label = "Cs-C=SC=SCbCb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   CS  u0 {1,S} {7,D}
4   Cb  u0 {1,S}
5   Cb  u0 {1,S}
6   S2d u0 {2,D}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1622,
    label = "Cs-C=S(Cds-Cd)CsCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cs  u0 {1,S}
5   Cs  u0 {1,S}
6   C   u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1623,
    label = "Cs-C=S(Cds-Cds)CsCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cs  u0 {1,S}
5   Cs  u0 {1,S}
6   Cd  u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1624,
    label = "Cs-C=S(Cds-Cdd)CsCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cs  u0 {1,S}
5   Cs  u0 {1,S}
6   Cdd u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1625,
    label = "Cs-C=S(Cds-Cdd-Cd)CsCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CS  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cs  u0 {1,S}
6   Cs  u0 {1,S}
7   S2d u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1626,
    label = "Cs-C=S(Cds-Cdd-S2d)CsCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CS  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cs  u0 {1,S}
6   Cs  u0 {1,S}
7   S2d u0 {3,D}
8   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1627,
    label = "Cs-C=SC=SCbCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   CS  u0 {1,S} {7,D}
4   Cb  u0 {1,S}
5   Ct  u0 {1,S}
6   S2d u0 {2,D}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1628,
    label = "Cs-C=S(Cds-Cd)CbCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   Ct  u0 {1,S}
6   C   u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1629,
    label = "Cs-C=S(Cds-Cds)CbCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   Ct  u0 {1,S}
6   Cd  u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1630,
    label = "Cs-C=S(Cds-Cdd)CbCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   Ct  u0 {1,S}
6   Cdd u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1631,
    label = "Cs-C=S(Cds-Cdd-S2d)CbCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CS  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   Ct  u0 {1,S}
7   S2d u0 {3,D}
8   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1632,
    label = "Cs-C=S(Cds-Cdd-Cd)CbCt",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CS  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   Ct  u0 {1,S}
7   S2d u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1633,
    label = "Cs-C=SC=SCsCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   CS  u0 {1,S} {7,D}
4   Cs  u0 {1,S}
5   Cs  u0 {1,S}
6   S2d u0 {2,D}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1634,
    label = "Cs-C=S(Cds-Cd)CbCb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   Cb  u0 {1,S}
6   C   u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1635,
    label = "Cs-C=S(Cds-Cds)CbCb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   Cb  u0 {1,S}
6   Cd  u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1636,
    label = "Cs-C=S(Cds-Cdd)CbCb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   Cb  u0 {1,S}
6   Cdd u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1637,
    label = "Cs-C=S(Cds-Cdd-S2d)CbCb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CS  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   Cb  u0 {1,S}
7   S2d u0 {3,D}
8   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1638,
    label = "Cs-C=S(Cds-Cdd-Cd)CbCb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CS  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   Cb  u0 {1,S}
7   S2d u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1639,
    label = "Cs-C=SC=S(Cds-Cd)Ct",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   CS  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {6,D}
5   Ct  u0 {1,S}
6   C   u0 {4,D}
7   S2d u0 {2,D}
8   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1640,
    label = "Cs-C=SC=S(Cds-Cds)Ct",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   CS  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {6,D}
5   Ct  u0 {1,S}
6   Cd  u0 {4,D}
7   S2d u0 {2,D}
8   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1641,
    label = "Cs-C=SC=S(Cds-Cdd)Ct",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   CS  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {6,D}
5   Ct  u0 {1,S}
6   Cdd u0 {4,D}
7   S2d u0 {2,D}
8   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1642,
    label = "Cs-C=SC=S(Cds-Cdd-Cd)Ct",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CS  u0 {1,S} {7,D}
4   CS  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   Ct  u0 {1,S}
7   S2d u0 {3,D}
8   S2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1643,
    label = "Cs-C=SC=S(Cds-Cdd-S2d)Ct",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CS  u0 {1,S} {7,D}
4   CS  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   Ct  u0 {1,S}
7   S2d u0 {3,D}
8   S2d u0 {4,D}
9   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1644,
    label = "Cs-C=SC=S(Cds-Cd)Cs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   CS  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {6,D}
5   Cs  u0 {1,S}
6   C   u0 {4,D}
7   S2d u0 {2,D}
8   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1645,
    label = "Cs-C=SC=S(Cds-Cds)Cs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   CS  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {6,D}
5   Cs  u0 {1,S}
6   Cd  u0 {4,D}
7   S2d u0 {2,D}
8   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1646,
    label = "Cs-C=SC=S(Cds-Cdd)Cs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   CS  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {6,D}
5   Cs  u0 {1,S}
6   Cdd u0 {4,D}
7   S2d u0 {2,D}
8   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1647,
    label = "Cs-C=SC=S(Cds-Cdd-S2d)Cs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CS  u0 {1,S} {7,D}
4   CS  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   Cs  u0 {1,S}
7   S2d u0 {3,D}
8   S2d u0 {4,D}
9   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1648,
    label = "Cs-C=SC=S(Cds-Cdd-Cd)Cs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CS  u0 {1,S} {7,D}
4   CS  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   Cs  u0 {1,S}
7   S2d u0 {3,D}
8   S2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1649,
    label = "Cs-C=SC=S(Cds-Cd)(Cds-Cd)",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {8,D}
3   CS  u0 {1,S} {9,D}
4   Cd  u0 {1,S} {6,D}
5   Cd  u0 {1,S} {7,D}
6   C   u0 {4,D}
7   C   u0 {5,D}
8   S2d u0 {2,D}
9   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1650,
    label = "Cs-C=SC=S(Cds-Cdd)(Cds-Cds)",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {8,D}
3   CS  u0 {1,S} {9,D}
4   Cd  u0 {1,S} {6,D}
5   Cd  u0 {1,S} {7,D}
6   Cdd u0 {4,D}
7   Cd  u0 {5,D}
8   S2d u0 {2,D}
9   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1651,
    label = "Cs-C=SC=S(Cds-Cdd-S2d)(Cds-Cds)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    CS  u0 {1,S} {8,D}
4    CS  u0 {1,S} {9,D}
5    Cd  u0 {1,S} {7,D}
6    Cdd u0 {2,D} {10,D}
7    Cd  u0 {5,D}
8    S2d u0 {3,D}
9    S2d u0 {4,D}
10   S2d u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1652,
    label = "Cs-C=SC=S(Cds-Cdd-Cd)(Cds-Cds)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    CS  u0 {1,S} {8,D}
4    CS  u0 {1,S} {9,D}
5    Cd  u0 {1,S} {7,D}
6    Cdd u0 {2,D} {10,D}
7    Cd  u0 {5,D}
8    S2d u0 {3,D}
9    S2d u0 {4,D}
10   C   u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1653,
    label = "Cs-C=SC=S(Cds-Cdd)(Cds-Cdd)",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {8,D}
3   CS  u0 {1,S} {9,D}
4   Cd  u0 {1,S} {6,D}
5   Cd  u0 {1,S} {7,D}
6   Cdd u0 {4,D}
7   Cdd u0 {5,D}
8   S2d u0 {2,D}
9   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1654,
    label = "Cs-C=SC=S(Cds-Cdd-S2d)(Cds-Cdd-S2d)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    CS  u0 {1,S} {8,D}
5    CS  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    S2d u0 {4,D}
9    S2d u0 {5,D}
10   S2d u0 {6,D}
11   S2d u0 {7,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1655,
    label = "Cs-C=SC=S(Cds-Cdd-S2d)(Cds-Cdd-Cd)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    CS  u0 {1,S} {8,D}
5    CS  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    S2d u0 {4,D}
9    S2d u0 {5,D}
10   S2d u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1656,
    label = "Cs-C=SC=S(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2    Cd  u0 {1,S} {6,D}
3    Cd  u0 {1,S} {7,D}
4    CS  u0 {1,S} {8,D}
5    CS  u0 {1,S} {9,D}
6    Cdd u0 {2,D} {10,D}
7    Cdd u0 {3,D} {11,D}
8    S2d u0 {4,D}
9    S2d u0 {5,D}
10   C   u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1657,
    label = "Cs-C=SC=S(Cds-Cds)(Cds-Cds)",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {8,D}
3   CS  u0 {1,S} {9,D}
4   Cd  u0 {1,S} {6,D}
5   Cd  u0 {1,S} {7,D}
6   Cd  u0 {4,D}
7   Cd  u0 {5,D}
8   S2d u0 {2,D}
9   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1658,
    label = "Cs-C=SC=S(Cds-Cd)Cb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   CS  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {6,D}
5   Cb  u0 {1,S}
6   C   u0 {4,D}
7   S2d u0 {2,D}
8   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1659,
    label = "Cs-C=SC=S(Cds-Cdd)Cb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   CS  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {6,D}
5   Cb  u0 {1,S}
6   Cdd u0 {4,D}
7   S2d u0 {2,D}
8   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1660,
    label = "Cs-C=SC=S(Cds-Cdd-S2d)Cb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CS  u0 {1,S} {7,D}
4   CS  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   Cb  u0 {1,S}
7   S2d u0 {3,D}
8   S2d u0 {4,D}
9   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1661,
    label = "Cs-C=SC=S(Cds-Cdd-Cd)Cb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CS  u0 {1,S} {7,D}
4   CS  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   Cb  u0 {1,S}
7   S2d u0 {3,D}
8   S2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1662,
    label = "Cs-C=SC=S(Cds-Cds)Cb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   CS  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {6,D}
5   Cb  u0 {1,S}
6   Cd  u0 {4,D}
7   S2d u0 {2,D}
8   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1663,
    label = "Cs-C=SCbCtCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Ct  u0 {1,S}
5   Cs  u0 {1,S}
6   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1664,
    label = "Cs-C=S(Cds-Cd)CbCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   Cs  u0 {1,S}
6   C   u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1665,
    label = "Cs-C=S(Cds-Cds)CbCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   Cs  u0 {1,S}
6   Cd  u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1666,
    label = "Cs-C=S(Cds-Cdd)CbCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   Cs  u0 {1,S}
6   Cdd u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1667,
    label = "Cs-C=S(Cds-Cdd-S2d)CbCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CS  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   Cs  u0 {1,S}
7   S2d u0 {3,D}
8   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1668,
    label = "Cs-C=S(Cds-Cdd-Cd)CbCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CS  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   Cs  u0 {1,S}
7   S2d u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1669,
    label = "Cs-C=S(Cds-Cd)(Cds-Cd)Cb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Cb  u0 {1,S}
6   C   u0 {3,D}
7   C   u0 {4,D}
8   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1670,
    label = "Cs-C=S(Cds-Cdd)(Cds-Cdd)Cb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Cb  u0 {1,S}
6   Cdd u0 {3,D}
7   Cdd u0 {4,D}
8   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1671,
    label = "Cs-C=S(Cds-Cdd-S2d)(Cds-Cdd-Cd)Cb",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    CS  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cb  u0 {1,S}
8    S2d u0 {4,D}
9    S2d u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1672,
    label = "Cs-C=S(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cb",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    CS  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cb  u0 {1,S}
8    S2d u0 {4,D}
9    C   u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1673,
    label = "Cs-C=S(Cds-Cdd-S2d)(Cds-Cdd-S2d)Cb",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    CS  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cb  u0 {1,S}
8    S2d u0 {4,D}
9    S2d u0 {5,D}
10   S2d u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1674,
    label = "Cs-C=S(Cds-Cds)(Cds-Cds)Cb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Cb  u0 {1,S}
6   Cd  u0 {3,D}
7   Cd  u0 {4,D}
8   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1675,
    label = "Cs-C=S(Cds-Cdd)(Cds-Cds)Cb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   Cb  u0 {1,S}
6   Cdd u0 {3,D}
7   Cd  u0 {4,D}
8   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1676,
    label = "Cs-C=S(Cds-Cdd-S2d)(Cds-Cds)Cb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CS  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {7,D}
5   Cdd u0 {2,D} {9,D}
6   Cb  u0 {1,S}
7   Cd  u0 {4,D}
8   S2d u0 {3,D}
9   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1677,
    label = "Cs-C=S(Cds-Cdd-Cd)(Cds-Cds)Cb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CS  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {7,D}
5   Cdd u0 {2,D} {9,D}
6   Cb  u0 {1,S}
7   Cd  u0 {4,D}
8   S2d u0 {3,D}
9   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1678,
    label = "Cs-C=SCbCbCb",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Cb  u0 {1,S}
5   Cb  u0 {1,S}
6   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1679,
    label = "Cs-C=SC=SCbCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   CS  u0 {1,S} {7,D}
4   Cb  u0 {1,S}
5   Cs  u0 {1,S}
6   S2d u0 {2,D}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1680,
    label = "Cs-CCCOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   C   u0 {1,S}
3   C   u0 {1,S}
4   C   u0 {1,S}
5   O2s u0 {1,S}
""",
    thermo = 'Cs-CsCsCsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1681,
    label = "Cs-CsCsCsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cs  u0 {1,S}
3   Cs  u0 {1,S}
4   Cs  u0 {1,S}
5   O2s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([23.99,31.2,34.89,36.47,36.78,36.05,34.4],'J/(mol*K)','+|-',[3.81,3.81,3.81,3.81,3.81,3.81,3.81]),
        H298 = (-20.3,'kJ/mol','+|-',3.24),
        S298 = (-144.38,'J/(mol*K)','+|-',4.44),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 1682,
    label = "Cs-CdsCsCsOs",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   [Cd,CO] u0 {1,S}
3   Cs      u0 {1,S}
4   Cs      u0 {1,S}
5   O2s     u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CsCsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1683,
    label = "Cs-(Cds-O2d)CsCsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   Cs  u0 {1,S}
5   O2s u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([28.15,35.17,38.11,38.72,37.49,35.88,33.45],'J/(mol*K)','+|-',[5.16,5.16,5.16,5.16,5.16,5.16,5.16]),
        H298 = (-10.9,'kJ/mol','+|-',4.39),
        S298 = (-148.7,'J/(mol*K)','+|-',6.02),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 1684,
    label = "Cs-(Cds-Cd)CsCsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   Cs  u0 {1,S}
5   O2s u0 {1,S}
6   C   u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([29.24,37.61,40.84,41.46,40.06,38.2,35.08],'J/(mol*K)','+|-',[3.81,3.81,3.81,3.81,3.81,3.81,3.81]),
        H298 = (-14.6,'kJ/mol','+|-',3.24),
        S298 = (-153.23,'J/(mol*K)','+|-',4.44),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 1685,
    label = "Cs-(Cds-Cds)CsCsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   Cs  u0 {1,S}
5   O2s u0 {1,S}
6   Cd  u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.63,6.79,7.95,8.4,8.8,8.44,8.44],'cal/(mol*K)','+|-',[0.2,0.2,0.2,0.2,0.2,0.2,0.2]),
        H298 = (-6.6,'kcal/mol','+|-',0.4),
        S298 = (-32.56,'cal/(mol*K)','+|-',0.2),
    ),
    shortDesc = """Cs-OCdCsCs BOZZELLI C/C3/O - (C/C3/H - C/Cb/C2/H), Hf-1 !!!WARNING! Cp1500 value taken as Cp1000""",
    longDesc = 
"""

""",
)

entry(
    index = 1686,
    label = "Cs-(Cds-Cdd)CsCsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   Cs  u0 {1,S}
5   O2s u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)CsCsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1687,
    label = "Cs-(Cds-Cdd-O2d)CsCsOs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cs  u0 {1,S}
5   Cs  u0 {1,S}
6   O2s u0 {1,S}
7   O2d u0 {3,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.39,9.66,10.03,10.07,9.64,9.26,8.74],'cal/(mol*K)','+|-',[0.2,0.2,0.2,0.2,0.2,0.2,0.2]),
        H298 = (-9.725,'kcal/mol','+|-',0.4),
        S298 = (-36.5,'cal/(mol*K)','+|-',0.2),
    ),
    shortDesc = """{C/CCO/O/C2} RAMAN & GREEN JPCA 2002, 106, 7937-7949""",
    longDesc = 
"""

""",
)

entry(
    index = 1688,
    label = "Cs-(Cds-Cdd-Cd)CsCsOs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cs  u0 {1,S}
5   Cs  u0 {1,S}
6   O2s u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)CsCsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1689,
    label = "Cs-OsCtCsCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   O2s u0 {1,S}
3   Ct  u0 {1,S}
4   Cs  u0 {1,S}
5   Cs  u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CsCsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1690,
    label = "Cs-CbCsCsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   Cs  u0 {1,S}
4   Cs  u0 {1,S}
5   O2s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.63,6.79,7.95,8.4,8.8,8.44,8.44],'cal/(mol*K)','+|-',[0.2,0.2,0.2,0.2,0.2,0.2,0.2]),
        H298 = (-6.6,'kcal/mol','+|-',0.4),
        S298 = (-32.56,'cal/(mol*K)','+|-',0.2),
    ),
    shortDesc = """Cs-OCbCsCs BOZZELLI C/C3/O - (C/C3/H - C/Cb/C2/H), Hf-1 !!!WARNING! Cp1500 value taken as Cp1000""",
    longDesc = 
"""

""",
)

entry(
    index = 1691,
    label = "Cs-CdsCdsCsOs",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   [Cd,CO] u0 {1,S}
3   [Cd,CO] u0 {1,S}
4   Cs      u0 {1,S}
5   O2s     u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1692,
    label = "Cs-(Cds-O2d)(Cds-O2d)CsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   CO  u0 {1,S} {7,D}
4   Cs  u0 {1,S}
5   O2s u0 {1,S}
6   O2d u0 {2,D}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-CsCsCsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1693,
    label = "Cs-(Cds-O2d)(Cds-Cd)CsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cs  u0 {1,S}
5   O2s u0 {1,S}
6   C   u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([33.75,42.15,45.09,44.95,41.74,38.55,34.46],'J/(mol*K)','+|-',[4.3,4.3,4.3,4.3,4.3,4.3,4.3]),
        H298 = (-3.9,'kJ/mol','+|-',3.66),
        S298 = (-158.3,'J/(mol*K)','+|-',5.02),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 1694,
    label = "Cs-(Cds-O2d)(Cds-Cds)CsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cs  u0 {1,S}
5   O2s u0 {1,S}
6   Cd  u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)CsCsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1695,
    label = "Cs-(Cds-O2d)(Cds-Cdd)CsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cs  u0 {1,S}
5   O2s u0 {1,S}
6   Cdd u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-Cd)CsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1696,
    label = "Cs-(Cds-O2d)(Cds-Cdd-O2d)CsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CO  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cs  u0 {1,S}
6   O2s u0 {1,S}
7   O2d u0 {3,D}
8   O2d u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)CsCsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1697,
    label = "Cs-(Cds-O2d)(Cds-Cdd-Cd)CsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CO  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cs  u0 {1,S}
6   O2s u0 {1,S}
7   O2d u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)CsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1698,
    label = "Cs-(Cds-Cd)(Cds-Cd)CsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cs  u0 {1,S}
5   O2s u0 {1,S}
6   C   u0 {2,D}
7   C   u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1699,
    label = "Cs-(Cds-Cds)(Cds-Cds)CsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cs  u0 {1,S}
5   O2s u0 {1,S}
6   Cd  u0 {2,D}
7   Cd  u0 {3,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.61,5.98,7.51,8.37,9,9.02,8.34],'cal/(mol*K)','+|-',[0.2,0.2,0.2,0.2,0.2,0.2,0.2]),
        H298 = (-8.01,'kcal/mol','+|-',0.4),
        S298 = (-34.34,'cal/(mol*K)','+|-',0.2),
    ),
    shortDesc = """Cs-OCdCdCs Hf jwb 697 S,Cp from C/Cd2/C2""",
    longDesc = 
"""

""",
)

entry(
    index = 1700,
    label = "Cs-(Cds-Cdd)(Cds-Cds)CsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cs  u0 {1,S}
5   O2s u0 {1,S}
6   Cdd u0 {2,D}
7   Cd  u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cds)CsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1701,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cds)CsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cs  u0 {1,S}
6   O2s u0 {1,S}
7   Cd  u0 {3,D}
8   O2d u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)CsCsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1702,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cds)CsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cs  u0 {1,S}
6   O2s u0 {1,S}
7   Cd  u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1703,
    label = "Cs-(Cds-Cdd)(Cds-Cdd)CsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cs  u0 {1,S}
5   O2s u0 {1,S}
6   Cdd u0 {2,D}
7   Cdd u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1704,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)CsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cs  u0 {1,S}
7   O2s u0 {1,S}
8   O2d u0 {4,D}
9   O2d u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1705,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-Cd)CsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cs  u0 {1,S}
7   O2s u0 {1,S}
8   O2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)(Cds-Cds)CsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1706,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cs  u0 {1,S}
7   O2s u0 {1,S}
8   C   u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1707,
    label = "Cs-CtCdsCsOs",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Ct      u0 {1,S}
3   [Cd,CO] u0 {1,S}
4   Cs      u0 {1,S}
5   O2s     u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CtCsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1708,
    label = "Cs-(Cds-O2d)CtCsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   Cs  u0 {1,S}
5   O2s u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)CsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1709,
    label = "Cs-(Cds-Cd)CtCsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   Cs  u0 {1,S}
5   O2s u0 {1,S}
6   C   u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)CtCsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1710,
    label = "Cs-(Cds-Cds)CtCsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   Cs  u0 {1,S}
5   O2s u0 {1,S}
6   Cd  u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1711,
    label = "Cs-(Cds-Cdd)CtCsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   Cs  u0 {1,S}
5   O2s u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)CtCsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1712,
    label = "Cs-(Cds-Cdd-O2d)CtCsOs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Ct  u0 {1,S}
5   Cs  u0 {1,S}
6   O2s u0 {1,S}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)(Cds-Cds)CsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1713,
    label = "Cs-(Cds-Cdd-Cd)CtCsOs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Ct  u0 {1,S}
5   Cs  u0 {1,S}
6   O2s u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)CtCsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1714,
    label = "Cs-CbCdsCsOs",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Cb      u0 {1,S}
3   [Cd,CO] u0 {1,S}
4   Cs      u0 {1,S}
5   O2s     u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CbCsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1715,
    label = "Cs-(Cds-O2d)CbCsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Cs  u0 {1,S}
5   O2s u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)CsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1716,
    label = "Cs-(Cds-Cd)CbCsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Cs  u0 {1,S}
5   O2s u0 {1,S}
6   C   u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)CbCsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1717,
    label = "Cs-(Cds-Cds)CbCsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Cs  u0 {1,S}
5   O2s u0 {1,S}
6   Cd  u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1718,
    label = "Cs-(Cds-Cdd)CbCsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Cs  u0 {1,S}
5   O2s u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)CbCsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1719,
    label = "Cs-(Cds-Cdd-O2d)CbCsOs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Cs  u0 {1,S}
6   O2s u0 {1,S}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)(Cds-Cds)CsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1720,
    label = "Cs-(Cds-Cdd-Cd)CbCsOs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Cs  u0 {1,S}
6   O2s u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)CbCsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1721,
    label = "Cs-CtCtCsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Ct  u0 {1,S}
3   Ct  u0 {1,S}
4   Cs  u0 {1,S}
5   O2s u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1722,
    label = "Cs-CbCtCsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   Ct  u0 {1,S}
4   Cs  u0 {1,S}
5   O2s u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CtCsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1723,
    label = "Cs-CbCbCsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   Cb  u0 {1,S}
4   Cs  u0 {1,S}
5   O2s u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1724,
    label = "Cs-CdsCdsCdsOs",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   [Cd,CO] u0 {1,S}
3   [Cd,CO] u0 {1,S}
4   [Cd,CO] u0 {1,S}
5   O2s     u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1725,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-O2d)O2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   CO  u0 {1,S} {7,D}
4   CO  u0 {1,S} {8,D}
5   O2s u0 {1,S}
6   O2d u0 {2,D}
7   O2d u0 {3,D}
8   O2d u0 {4,D}
""",
    thermo = 'Cs-CsCsCsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1726,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-Cd)O2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   CO  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {6,D}
5   O2s u0 {1,S}
6   C   u0 {4,D}
7   O2d u0 {2,D}
8   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1727,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-Cds)O2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   CO  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {6,D}
5   O2s u0 {1,S}
6   Cd  u0 {4,D}
7   O2d u0 {2,D}
8   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)CsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1728,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-Cdd)O2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   CO  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {6,D}
5   O2s u0 {1,S}
6   Cdd u0 {4,D}
7   O2d u0 {2,D}
8   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-Cdd-Cd)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1729,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-Cdd-O2d)O2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CO  u0 {1,S} {7,D}
4   CO  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   O2s u0 {1,S}
7   O2d u0 {3,D}
8   O2d u0 {4,D}
9   O2d u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)CsCsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1730,
    label = "Cs-(Cds-O2d)(Cds-O2d)(Cds-Cdd-Cd)O2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CO  u0 {1,S} {7,D}
4   CO  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   O2s u0 {1,S}
7   O2d u0 {3,D}
8   O2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1731,
    label = "Cs-(Cds-O2d)(Cds-Cd)(Cds-Cd)O2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   O2s u0 {1,S}
6   C   u0 {3,D}
7   C   u0 {4,D}
8   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([36.85,46.04,49,48.85,45.61,42.23,37.25],'J/(mol*K)','+|-',[4.09,4.09,4.09,4.09,4.09,4.09,4.09]),
        H298 = (3,'kJ/mol','+|-',3.49),
        S298 = (-160.69,'J/(mol*K)','+|-',4.77),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 1732,
    label = "Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)O2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   O2s u0 {1,S}
6   Cd  u0 {3,D}
7   Cd  u0 {4,D}
8   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)CsCsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1733,
    label = "Cs-(Cds-O2d)(Cds-Cdd)(Cds-Cds)O2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   O2s u0 {1,S}
6   Cdd u0 {3,D}
7   Cd  u0 {4,D}
8   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-Cd)(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1734,
    label = "Cs-(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cds)O2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CO  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {7,D}
5   Cdd u0 {2,D} {9,D}
6   O2s u0 {1,S}
7   Cd  u0 {4,D}
8   O2d u0 {3,D}
9   O2d u0 {5,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-O2d)CsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1735,
    label = "Cs-(Cds-O2d)(Cds-Cdd-Cd)(Cds-Cds)O2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CO  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {7,D}
5   Cdd u0 {2,D} {9,D}
6   O2s u0 {1,S}
7   Cd  u0 {4,D}
8   O2d u0 {3,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1736,
    label = "Cs-(Cds-O2d)(Cds-Cdd)(Cds-Cdd)O2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   O2s u0 {1,S}
6   Cdd u0 {3,D}
7   Cdd u0 {4,D}
8   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-Cd)(Cds-Cdd-Cd)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1737,
    label = "Cs-(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cdd-O2d)O2s",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    CO  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    O2s u0 {1,S}
8    O2d u0 {4,D}
9    O2d u0 {5,D}
10   O2d u0 {6,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)CsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1738,
    label = "Cs-(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cdd-Cd)O2s",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    CO  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    O2s u0 {1,S}
8    O2d u0 {4,D}
9    O2d u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1739,
    label = "Cs-(Cds-O2d)(Cds-Cdd-Cd)(Cds-Cdd-Cd)O2s",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    CO  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    O2s u0 {1,S}
8    O2d u0 {4,D}
9    C   u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1740,
    label = "Cs-(Cds-Cd)(Cds-Cd)(Cds-Cd)O2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   O2s u0 {1,S}
6   C   u0 {2,D}
7   C   u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1741,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)O2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   O2s u0 {1,S}
6   Cd  u0 {2,D}
7   Cd  u0 {3,D}
8   Cd  u0 {4,D}
""",
    thermo = 'Cs-CsCsCsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1742,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd)O2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   O2s u0 {1,S}
6   Cd  u0 {2,D}
7   Cd  u0 {3,D}
8   Cdd u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1743,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)O2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   O2s u0 {1,S}
7   Cd  u0 {3,D}
8   Cd  u0 {4,D}
9   O2d u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)CsCsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1744,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)O2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   O2s u0 {1,S}
7   Cd  u0 {3,D}
8   Cd  u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1745,
    label = "Cs-(Cds-Cds)(Cds-Cdd)(Cds-Cdd)O2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   O2s u0 {1,S}
6   Cd  u0 {2,D}
7   Cdd u0 {3,D}
8   Cdd u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1746,
    label = "Cs-(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-O2d)O2s",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    O2s u0 {1,S}
8    Cd  u0 {4,D}
9    O2d u0 {5,D}
10   O2d u0 {6,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)CsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1747,
    label = "Cs-(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-Cd)O2s",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    O2s u0 {1,S}
8    Cd  u0 {4,D}
9    O2d u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1748,
    label = "Cs-(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)O2s",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    O2s u0 {1,S}
8    Cd  u0 {4,D}
9    C   u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1749,
    label = "Cs-(Cds-Cdd)(Cds-Cdd)(Cds-Cdd)O2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   O2s u0 {1,S}
6   Cdd u0 {2,D}
7   Cdd u0 {3,D}
8   Cdd u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1750,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)(Cds-Cdd-O2d)O2s",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    O2s u0 {1,S}
9    O2d u0 {5,D}
10   O2d u0 {6,D}
11   O2d u0 {7,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1751,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)(Cds-Cdd-Cd)O2s",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    O2s u0 {1,S}
9    O2d u0 {5,D}
10   O2d u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-O2d)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1752,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-Cd)(Cds-Cdd-Cd)O2s",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    O2s u0 {1,S}
9    O2d u0 {5,D}
10   C   u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1753,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)O2s",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    O2s u0 {1,S}
9    C   u0 {5,D}
10   C   u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1754,
    label = "Cs-CtCdsCdsOs",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Ct      u0 {1,S}
3   [Cd,CO] u0 {1,S}
4   [Cd,CO] u0 {1,S}
5   O2s     u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CtOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1755,
    label = "Cs-(Cds-O2d)(Cds-O2d)CtOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   CO  u0 {1,S} {7,D}
4   Ct  u0 {1,S}
5   O2s u0 {1,S}
6   O2d u0 {2,D}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1756,
    label = "Cs-(Cds-O2d)(Cds-Cd)CtOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Ct  u0 {1,S}
5   O2s u0 {1,S}
6   C   u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)CtOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1757,
    label = "Cs-(Cds-O2d)(Cds-Cds)CtOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Ct  u0 {1,S}
5   O2s u0 {1,S}
6   Cd  u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1758,
    label = "Cs-(Cds-O2d)(Cds-Cdd)CtOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Ct  u0 {1,S}
5   O2s u0 {1,S}
6   Cdd u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-Cd)CtOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1759,
    label = "Cs-(Cds-O2d)(Cds-Cdd-O2d)CtOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CO  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Ct  u0 {1,S}
6   O2s u0 {1,S}
7   O2d u0 {3,D}
8   O2d u0 {4,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1760,
    label = "Cs-(Cds-O2d)(Cds-Cdd-Cd)CtOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CO  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Ct  u0 {1,S}
6   O2s u0 {1,S}
7   O2d u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)CtOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1761,
    label = "Cs-(Cds-Cd)(Cds-Cd)CtOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Ct  u0 {1,S}
5   O2s u0 {1,S}
6   C   u0 {2,D}
7   C   u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CtOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1762,
    label = "Cs-(Cds-Cds)(Cds-Cds)CtOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Ct  u0 {1,S}
5   O2s u0 {1,S}
6   Cd  u0 {2,D}
7   Cd  u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1763,
    label = "Cs-(Cds-Cdd)(Cds-Cds)CtOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Ct  u0 {1,S}
5   O2s u0 {1,S}
6   Cdd u0 {2,D}
7   Cd  u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cds)CtOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1764,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cds)CtOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Ct  u0 {1,S}
6   O2s u0 {1,S}
7   Cd  u0 {3,D}
8   O2d u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1765,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cds)CtOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Ct  u0 {1,S}
6   O2s u0 {1,S}
7   Cd  u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CtOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1766,
    label = "Cs-(Cds-Cdd)(Cds-Cdd)CtOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Ct  u0 {1,S}
5   O2s u0 {1,S}
6   Cdd u0 {2,D}
7   Cdd u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CtOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1767,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)CtOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Ct  u0 {1,S}
7   O2s u0 {1,S}
8   O2d u0 {4,D}
9   O2d u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-O2d)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1768,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-Cd)CtOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Ct  u0 {1,S}
7   O2s u0 {1,S}
8   O2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)(Cds-Cds)CtOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1769,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CtOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Ct  u0 {1,S}
7   O2s u0 {1,S}
8   C   u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CtOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1770,
    label = "Cs-CbCdsCdsOs",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Cb      u0 {1,S}
3   [Cd,CO] u0 {1,S}
4   [Cd,CO] u0 {1,S}
5   O2s     u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CbOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1771,
    label = "Cs-(Cds-O2d)(Cds-O2d)CbOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   CO  u0 {1,S} {7,D}
4   Cb  u0 {1,S}
5   O2s u0 {1,S}
6   O2d u0 {2,D}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-O2d)(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1772,
    label = "Cs-(Cds-O2d)(Cds-Cd)CbOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   O2s u0 {1,S}
6   C   u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)CbOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1773,
    label = "Cs-(Cds-O2d)(Cds-Cds)CbOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   O2s u0 {1,S}
6   Cd  u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1774,
    label = "Cs-(Cds-O2d)(Cds-Cdd)CbOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   O2s u0 {1,S}
6   Cdd u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-Cd)CbOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1775,
    label = "Cs-(Cds-O2d)(Cds-Cdd-O2d)CbOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CO  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   O2s u0 {1,S}
7   O2d u0 {3,D}
8   O2d u0 {4,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-O2d)(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1776,
    label = "Cs-(Cds-O2d)(Cds-Cdd-Cd)CbOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CO  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   O2s u0 {1,S}
7   O2d u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)CbOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1777,
    label = "Cs-(Cds-Cd)(Cds-Cd)CbOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cb  u0 {1,S}
5   O2s u0 {1,S}
6   C   u0 {2,D}
7   C   u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CbOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1778,
    label = "Cs-(Cds-Cds)(Cds-Cds)CbOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cb  u0 {1,S}
5   O2s u0 {1,S}
6   Cd  u0 {2,D}
7   Cd  u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1779,
    label = "Cs-(Cds-Cdd)(Cds-Cds)CbOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cb  u0 {1,S}
5   O2s u0 {1,S}
6   Cdd u0 {2,D}
7   Cd  u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cds)CbOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1780,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cds)CbOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   O2s u0 {1,S}
7   Cd  u0 {3,D}
8   O2d u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1781,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cds)CbOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   O2s u0 {1,S}
7   Cd  u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CbOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1782,
    label = "Cs-(Cds-Cdd)(Cds-Cdd)CbOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cb  u0 {1,S}
5   O2s u0 {1,S}
6   Cdd u0 {2,D}
7   Cdd u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CbOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1783,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)CbOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cb  u0 {1,S}
7   O2s u0 {1,S}
8   O2d u0 {4,D}
9   O2d u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cdd-O2d)(Cds-Cdd-O2d)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1784,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-Cd)CbOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cb  u0 {1,S}
7   O2s u0 {1,S}
8   O2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)(Cds-Cds)CbOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1785,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CbOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cb  u0 {1,S}
7   O2s u0 {1,S}
8   C   u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CbOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1786,
    label = "Cs-CtCtCdsOs",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Ct      u0 {1,S}
3   Ct      u0 {1,S}
4   [Cd,CO] u0 {1,S}
5   O2s     u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CtCtOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1787,
    label = "Cs-(Cds-O2d)CtCtOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   Ct  u0 {1,S}
5   O2s u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1788,
    label = "Cs-(Cds-Cd)CtCtOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   Ct  u0 {1,S}
5   O2s u0 {1,S}
6   C   u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)CtCtOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1789,
    label = "Cs-(Cds-Cds)CtCtOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   Ct  u0 {1,S}
5   O2s u0 {1,S}
6   Cd  u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1790,
    label = "Cs-(Cds-Cdd)CtCtOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   Ct  u0 {1,S}
5   O2s u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)CtCtOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1791,
    label = "Cs-(Cds-Cdd-O2d)CtCtOs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Ct  u0 {1,S}
5   Ct  u0 {1,S}
6   O2s u0 {1,S}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1792,
    label = "Cs-(Cds-Cdd-Cd)CtCtOs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Ct  u0 {1,S}
5   Ct  u0 {1,S}
6   O2s u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)CtCtOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1793,
    label = "Cs-CbCtCdsOs",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Cb      u0 {1,S}
3   Ct      u0 {1,S}
4   [Cd,CO] u0 {1,S}
5   O2s     u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CbCtOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1794,
    label = "Cs-(Cds-O2d)CbCtOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Ct  u0 {1,S}
5   O2s u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)CtOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1795,
    label = "Cs-(Cds-Cd)CbCtOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Ct  u0 {1,S}
5   O2s u0 {1,S}
6   C   u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)CbCtOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1796,
    label = "Cs-(Cds-Cds)CbCtOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Ct  u0 {1,S}
5   O2s u0 {1,S}
6   Cd  u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CtOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1797,
    label = "Cs-(Cds-Cdd)CbCtOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Ct  u0 {1,S}
5   O2s u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)CbCtOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1798,
    label = "Cs-(Cds-Cdd-O2d)CbCtOs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Ct  u0 {1,S}
6   O2s u0 {1,S}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)(Cds-Cds)CtOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1799,
    label = "Cs-(Cds-Cdd-Cd)CbCtOs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Ct  u0 {1,S}
6   O2s u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)CbCtOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1800,
    label = "Cs-CbCbCdsOs",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Cb      u0 {1,S}
3   Cb      u0 {1,S}
4   [Cd,CO] u0 {1,S}
5   O2s     u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CbCbOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1801,
    label = "Cs-(Cds-O2d)CbCbOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Cb  u0 {1,S}
5   O2s u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1802,
    label = "Cs-(Cds-Cd)CbCbOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Cb  u0 {1,S}
5   O2s u0 {1,S}
6   C   u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)CbCbOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1803,
    label = "Cs-(Cds-Cds)CbCbOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Cb  u0 {1,S}
5   O2s u0 {1,S}
6   Cd  u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1804,
    label = "Cs-(Cds-Cdd)CbCbOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Cb  u0 {1,S}
5   O2s u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)CbCbOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1805,
    label = "Cs-(Cds-Cdd-O2d)CbCbOs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Cb  u0 {1,S}
6   O2s u0 {1,S}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-O2d)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1806,
    label = "Cs-(Cds-Cdd-Cd)CbCbOs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Cb  u0 {1,S}
6   O2s u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)CbCbOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1807,
    label = "Cs-CtCtCtOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Ct  u0 {1,S}
3   Ct  u0 {1,S}
4   Ct  u0 {1,S}
5   O2s u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1808,
    label = "Cs-CbCtCtOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   Ct  u0 {1,S}
4   Ct  u0 {1,S}
5   O2s u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CtCtOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1809,
    label = "Cs-CbCbCtOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   Cb  u0 {1,S}
4   Ct  u0 {1,S}
5   O2s u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)CtOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1810,
    label = "Cs-CbCbCbOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   Cb  u0 {1,S}
4   Cb  u0 {1,S}
5   O2s u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)O2s',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1811,
    label = "Cs-CCOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   C   u0 {1,S}
3   C   u0 {1,S}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
""",
    thermo = 'Cs-CsCsOsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1812,
    label = "Cs-CsCsOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cs  u0 {1,S}
3   Cs  u0 {1,S}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([27.88,43.75,51.85,54,50.77,45.94,38.31],'J/(mol*K)','+|-',[5.77,5.77,5.77,5.77,5.77,5.77,5.77]),
        H298 = (-69.2,'kJ/mol','+|-',4.92),
        S298 = (-163.77,'J/(mol*K)','+|-',6.74),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 1813,
    label = "Cs-CdsCsOsOs",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   [Cd,CO] u0 {1,S}
3   Cs      u0 {1,S}
4   O2s     u0 {1,S}
5   O2s     u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CsOsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1814,
    label = "Cs-(Cds-O2d)CsOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = 'Cs-CsCsOsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1815,
    label = "Cs-(Cds-Cd)CsOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
6   C   u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([27.95,42.92,51.33,54.81,53.92,49.73,41.11],'J/(mol*K)','+|-',[5.77,5.77,5.77,5.77,5.77,5.77,5.77]),
        H298 = (-62.8,'kJ/mol','+|-',4.92),
        S298 = (-170.44,'J/(mol*K)','+|-',6.74),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 1816,
    label = "Cs-(Cds-Cds)CsOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
6   Cd  u0 {2,D}
""",
    thermo = 'Cs-CsCsOsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1817,
    label = "Cs-(Cds-Cdd)CsOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)CsOsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1818,
    label = "Cs-(Cds-Cdd-O2d)CsOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cs  u0 {1,S}
5   O2s u0 {1,S}
6   O2s u0 {1,S}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)CsOsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1819,
    label = "Cs-(Cds-Cdd-Cd)CsOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cs  u0 {1,S}
5   O2s u0 {1,S}
6   O2s u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)CsOsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1820,
    label = "Cs-CdsCdsOsOs",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   [Cd,CO] u0 {1,S}
3   [Cd,CO] u0 {1,S}
4   O2s     u0 {1,S}
5   O2s     u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)OsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1821,
    label = "Cs-(Cds-O2d)(Cds-O2d)OsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   CO  u0 {1,S} {7,D}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
6   O2d u0 {2,D}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-CsCsOsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1822,
    label = "Cs-(Cds-O2d)(Cds-Cd)OsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
6   C   u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)OsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1823,
    label = "Cs-(Cds-O2d)(Cds-Cds)OsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
6   Cd  u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)CsOsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1824,
    label = "Cs-(Cds-O2d)(Cds-Cdd)OsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
6   Cdd u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-Cd)OsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1825,
    label = "Cs-(Cds-O2d)(Cds-Cdd-O2d)OsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CO  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   O2s u0 {1,S}
6   O2s u0 {1,S}
7   O2d u0 {3,D}
8   O2d u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)CsOsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1826,
    label = "Cs-(Cds-O2d)(Cds-Cdd-Cd)OsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CO  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   O2s u0 {1,S}
6   O2s u0 {1,S}
7   O2d u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)OsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1827,
    label = "Cs-(Cds-Cd)(Cds-Cd)OsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
6   C   u0 {2,D}
7   C   u0 {3,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([30.08,45.85,54.7,58.39,57.78,53.65,44.31],'J/(mol*K)','+|-',[5.77,5.77,5.77,5.77,5.77,5.77,5.77]),
        H298 = (-55.7,'kJ/mol','+|-',4.92),
        S298 = (-179.76,'J/(mol*K)','+|-',6.74),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 1828,
    label = "Cs-(Cds-Cds)(Cds-Cds)OsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
6   Cd  u0 {2,D}
7   Cd  u0 {3,D}
""",
    thermo = 'Cs-CsCsOsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1829,
    label = "Cs-(Cds-Cdd)(Cds-Cds)OsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
6   Cdd u0 {2,D}
7   Cd  u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cds)OsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1830,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cds)OsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   O2s u0 {1,S}
6   O2s u0 {1,S}
7   Cd  u0 {3,D}
8   O2d u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)CsOsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1831,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cds)OsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   O2s u0 {1,S}
6   O2s u0 {1,S}
7   Cd  u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)OsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1832,
    label = "Cs-(Cds-Cdd)(Cds-Cdd)OsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
6   Cdd u0 {2,D}
7   Cdd u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)OsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1833,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)OsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   O2s u0 {1,S}
7   O2s u0 {1,S}
8   O2d u0 {4,D}
9   O2d u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)OsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1834,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-Cd)OsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   O2s u0 {1,S}
7   O2s u0 {1,S}
8   O2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)(Cds-Cds)OsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1835,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)OsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   O2s u0 {1,S}
7   O2s u0 {1,S}
8   C   u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)OsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1836,
    label = "Cs-CtCsOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Ct  u0 {1,S}
3   Cs  u0 {1,S}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CsOsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1837,
    label = "Cs-CtCdsOsOs",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Ct      u0 {1,S}
3   [Cd,CO] u0 {1,S}
4   O2s     u0 {1,S}
5   O2s     u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CtOsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1838,
    label = "Cs-(Cds-O2d)CtOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)OsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1839,
    label = "Cs-(Cds-Cd)CtOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
6   C   u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)CtOsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1840,
    label = "Cs-(Cds-Cds)CtOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
6   Cd  u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)OsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1841,
    label = "Cs-(Cds-Cdd)CtOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)CtOsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1842,
    label = "Cs-(Cds-Cdd-O2d)CtOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Ct  u0 {1,S}
5   O2s u0 {1,S}
6   O2s u0 {1,S}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)(Cds-Cds)OsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1843,
    label = "Cs-(Cds-Cdd-Cd)CtOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Ct  u0 {1,S}
5   O2s u0 {1,S}
6   O2s u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)CtOsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1844,
    label = "Cs-CtCtOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Ct  u0 {1,S}
3   Ct  u0 {1,S}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)OsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1845,
    label = "Cs-CbCsOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   Cs  u0 {1,S}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CsOsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1846,
    label = "Cs-CbCdsOsOs",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Cb      u0 {1,S}
3   [Cd,CO] u0 {1,S}
4   O2s     u0 {1,S}
5   O2s     u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CbOsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1847,
    label = "Cs-(Cds-O2d)CbOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)OsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1848,
    label = "Cs-(Cds-Cd)CbOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
6   C   u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)CbOsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1849,
    label = "Cs-(Cds-Cds)CbOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
6   Cd  u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)OsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1850,
    label = "Cs-(Cds-Cdd)CbOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)CbOsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1851,
    label = "Cs-(Cds-Cdd-O2d)CbOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   O2s u0 {1,S}
6   O2s u0 {1,S}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)(Cds-Cds)OsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1852,
    label = "Cs-(Cds-Cdd-Cd)CbOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   O2s u0 {1,S}
6   O2s u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)CbOsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1853,
    label = "Cs-CbCtOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   Ct  u0 {1,S}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CtOsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1854,
    label = "Cs-CbCbOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   Cb  u0 {1,S}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)OsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1855,
    label = "Cs-COsOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   C   u0 {1,S}
3   O2s u0 {1,S}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
""",
    thermo = 'Cs-CsOsOsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1856,
    label = "Cs-CsOsOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cs  u0 {1,S}
3   O2s u0 {1,S}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.33,6.19,7.25,7.7,8.2,8.24,8.24],'cal/(mol*K)','+|-',[0.2,0.2,0.2,0.2,0.2,0.2,0.2]),
        H298 = (-19,'kcal/mol','+|-',0.4),
        S298 = (-33.56,'cal/(mol*K)','+|-',0.2),
    ),
    shortDesc = """Cs-OOOCs BOZZELLI est !!!WARNING! Cp1500 value taken as Cp1000""",
    longDesc = 
"""

""",
)

entry(
    index = 1857,
    label = "Cs-CdsOsOsOs",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   [Cd,CO] u0 {1,S}
3   O2s     u0 {1,S}
4   O2s     u0 {1,S}
5   O2s     u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)OsOsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1858,
    label = "Cs-(Cds-O2d)OsOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   O2s u0 {1,S}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = 'Cs-CsOsOsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1859,
    label = "Cs-(Cds-Cd)OsOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   O2s u0 {1,S}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
6   C   u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)OsOsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1860,
    label = "Cs-(Cds-Cds)OsOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   O2s u0 {1,S}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
6   Cd  u0 {2,D}
""",
    thermo = 'Cs-CsOsOsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1861,
    label = "Cs-(Cds-Cdd)OsOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   O2s u0 {1,S}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)OsOsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1862,
    label = "Cs-(Cds-Cdd-O2d)OsOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
6   O2s u0 {1,S}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)OsOsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1863,
    label = "Cs-(Cds-Cdd-Cd)OsOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
6   O2s u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)OsOsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1864,
    label = "Cs-CtOsOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Ct  u0 {1,S}
3   O2s u0 {1,S}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)OsOsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1865,
    label = "Cs-CbOsOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   O2s u0 {1,S}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)OsOsOs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1866,
    label = "Cs-OsOsOsOs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   O2s u0 {1,S}
3   O2s u0 {1,S}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.33,6.13,7.25,7.7,8.2,8.24,8.24],'cal/(mol*K)','+|-',[0.2,0.2,0.2,0.2,0.2,0.2,0.2]),
        H298 = (-23,'kcal/mol','+|-',0.4),
        S298 = (-35.56,'cal/(mol*K)','+|-',0.2),
    ),
    shortDesc = """Cs-OOOO BOZZELLI est !!!WARNING! Cp1500 value taken as Cp1000""",
    longDesc = 
"""

""",
)

entry(
    index = 1867,
    label = "Cs-COsOsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   C   u0 {1,S}
3   O2s u0 {1,S}
4   O2s u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = 'Cs-CsOsOsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1868,
    label = "Cs-CsOsOsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cs  u0 {1,S}
3   O2s u0 {1,S}
4   O2s u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.25,7.1,8.81,9.55,10.31,11.05,11.05],'cal/(mol*K)','+|-',[0.12,0.12,0.12,0.12,0.12,0.12,0.12]),
        H298 = (-16,'kcal/mol','+|-',0.24),
        S298 = (-12.07,'cal/(mol*K)','+|-',0.12),
    ),
    shortDesc = """Cs-OOCsH BENSON Hf, BOZZELLI C/C3/H - C/C2/O/H !!!WARNING! Cp1500 value taken as Cp1000""",
    longDesc = 
"""

""",
)

entry(
    index = 1869,
    label = "Cs-CdsOsOsH",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   [Cd,CO] u0 {1,S}
3   O2s     u0 {1,S}
4   O2s     u0 {1,S}
5   H       u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)OsOsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1870,
    label = "Cs-(Cds-O2d)OsOsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   O2s u0 {1,S}
4   O2s u0 {1,S}
5   H   u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = 'Cs-CsOsOsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1871,
    label = "Cs-(Cds-Cd)OsOsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   O2s u0 {1,S}
4   O2s u0 {1,S}
5   H   u0 {1,S}
6   C   u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)OsOsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1872,
    label = "Cs-(Cds-Cds)OsOsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   O2s u0 {1,S}
4   O2s u0 {1,S}
5   H   u0 {1,S}
6   Cd  u0 {2,D}
""",
    thermo = 'Cs-CsOsOsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1873,
    label = "Cs-(Cds-Cdd)OsOsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   O2s u0 {1,S}
4   O2s u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)OsOsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1874,
    label = "Cs-(Cds-Cdd-O2d)OsOsH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
6   H   u0 {1,S}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)OsOsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1875,
    label = "Cs-(Cds-Cdd-Cd)OsOsH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   O2s u0 {1,S}
5   O2s u0 {1,S}
6   H   u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)OsOsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1876,
    label = "Cs-CtOsOsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Ct  u0 {1,S}
3   O2s u0 {1,S}
4   O2s u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)OsOsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1877,
    label = "Cs-CbOsOsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   O2s u0 {1,S}
4   O2s u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)OsOsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1878,
    label = "Cs-COsSH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   C   u0 {1,S}
3   O2s u0 {1,S}
4   S   u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = 'Cs-CsOsSH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1879,
    label = "Cs-CsOsSH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cs  u0 {1,S}
3   O2s u0 {1,S}
4   S   u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = 'Cs-CsOsS2H',
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2018""",
    longDesc = 
"""
"
""",
)

entry(
    index = 1880,
    label = "Cs-CsOsS2H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cs  u0 {1,S}
3   O2s u0 {1,S}
4   S2s u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.76,10.63,11.62,12.1,13.03,14.22,12.83],'cal/(mol*K)'),
        H298 = (4.17,'kcal/mol'),
        S298 = (-13.3,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 1881,
    label = "Cs-CsOsS4H",
    group = 
"""
1 * Cs                u0 {2,S} {3,S} {4,S} {5,S}
2   Cs                u0 {1,S}
3   O2s               u0 {1,S}
4   [S4s,S4d,S4b,S4t] u0 {1,S}
5   H                 u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.67,8.88,10.44,11.92,15.43,17.07,15.18],'cal/(mol*K)'),
        H298 = (5.25,'kcal/mol'),
        S298 = (-20.4,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 1882,
    label = "Cs-CdsOsSsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S}
3   O2s u0 {1,S}
4   S2s u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.88,11.94,12.92,13.33,14.67,15.02,13.48],'cal/(mol*K)'),
        H298 = (5.37,'kcal/mol'),
        S298 = (-11.83,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 1883,
    label = "Cs-CtOsSsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Ct  u0 {1,S}
3   O2s u0 {1,S}
4   S2s u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1884,
    label = "Cs-CbOsSsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   O2s u0 {1,S}
4   S2s u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.94,11.75,12.59,12.98,14.4,14.78,12.9],'cal/(mol*K)'),
        H298 = (4.61,'kcal/mol'),
        S298 = (-13.38,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 1885,
    label = "Cs-CCOsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   C   u0 {1,S}
3   C   u0 {1,S}
4   O2s u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1886,
    label = "Cs-CsCsOsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cs  u0 {1,S}
3   Cs  u0 {1,S}
4   O2s u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.94,10.89,11.45,11.37,11.9,11.9,9.51],'cal/(mol*K)'),
        H298 = (6.54,'kcal/mol'),
        S298 = (-38.36,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 1887,
    label = "Cs-COsOsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   C   u0 {1,S}
3   O2s u0 {1,S}
4   O2s u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = 'Cs-CsOsOsSs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1888,
    label = "Cs-CsOsOsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cs  u0 {1,S}
3   O2s u0 {1,S}
4   O2s u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.58,8.97,9.76,10.23,11.93,12.31,10.02],'cal/(mol*K)'),
        H298 = (-5.27,'kcal/mol'),
        S298 = (-33.54,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 1889,
    label = "Cs-CCOsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   C   u0 {1,S}
3   C   u0 {1,S}
4   O2s u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = 'Cs-CsCsOsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1890,
    label = "Cs-CsCsOsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cs  u0 {1,S}
3   Cs  u0 {1,S}
4   O2s u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([21.99,29.03,34.22,37.78,41.96,44.27,47.11],'J/(mol*K)','+|-',[3.32,3.32,3.32,3.32,3.32,3.32,3.32]),
        H298 = (-25.1,'kJ/mol','+|-',2.83),
        S298 = (-52.05,'J/(mol*K)','+|-',3.88),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 1891,
    label = "Cs-CdsCsOsH",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   [Cd,CO] u0 {1,S}
3   Cs      u0 {1,S}
4   O2s     u0 {1,S}
5   H       u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CsOsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1892,
    label = "Cs-(Cds-O2d)CsOsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   O2s u0 {1,S}
5   H   u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.47,6.82,8.45,9.17,10.24,10.8,11.02],'cal/(mol*K)','+|-',[0.12,0.12,0.12,0.12,0.12,0.12,0.12]),
        H298 = (-6,'kcal/mol','+|-',0.24),
        S298 = (-11.1,'cal/(mol*K)','+|-',0.12),
    ),
    shortDesc = """Cs-OCOCsH BOZZELLI""",
    longDesc = 
"""

""",
)

entry(
    index = 1893,
    label = "Cs-(Cds-Cd)CsOsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   O2s u0 {1,S}
5   H   u0 {1,S}
6   C   u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([29.84,38.86,43.83,46.37,48.34,49.06,49.94],'J/(mol*K)','+|-',[3.74,3.74,3.74,3.74,3.74,3.74,3.74]),
        H298 = (-24,'kJ/mol','+|-',3.19),
        S298 = (-61.06,'J/(mol*K)','+|-',4.36),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 1894,
    label = "Cs-(Cds-Cds)CsOsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   O2s u0 {1,S}
5   H   u0 {1,S}
6   Cd  u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.47,6.82,8.45,9.17,10.24,10.8,11.02],'cal/(mol*K)','+|-',[0.12,0.12,0.12,0.12,0.12,0.12,0.12]),
        H298 = (-6,'kcal/mol','+|-',0.24),
        S298 = (-11.1,'cal/(mol*K)','+|-',0.12),
    ),
    shortDesc = """Cs-OCdCsH BOZZELLI""",
    longDesc = 
"""

""",
)

entry(
    index = 1895,
    label = "Cs-(Cds-Cdd)CsOsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   O2s u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)CsOsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1896,
    label = "Cs-(Cds-Cdd-O2d)CsOsH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cs  u0 {1,S}
5   O2s u0 {1,S}
6   H   u0 {1,S}
7   O2d u0 {3,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.2,8.49,9.33,9.92,10.5,10.92,11.71],'cal/(mol*K)','+|-',[0.12,0.12,0.12,0.12,0.12,0.12,0.12]),
        H298 = (-8.37,'kcal/mol','+|-',0.24),
        S298 = (-13.04,'cal/(mol*K)','+|-',0.12),
    ),
    shortDesc = """{C/CCO/O/C/H} RAMAN & GREEN JPCA 2002, 106, 7937-7949""",
    longDesc = 
"""

""",
)

entry(
    index = 1897,
    label = "Cs-(Cds-Cdd-Cd)CsOsH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cs  u0 {1,S}
5   O2s u0 {1,S}
6   H   u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)CsOsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1898,
    label = "Cs-CdsCdsOsH",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   [Cd,CO] u0 {1,S}
3   [Cd,CO] u0 {1,S}
4   O2s     u0 {1,S}
5   H       u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)OsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1899,
    label = "Cs-(Cds-O2d)(Cds-O2d)OsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   CO  u0 {1,S} {7,D}
4   O2s u0 {1,S}
5   H   u0 {1,S}
6   O2d u0 {2,D}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-CsCsOsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1900,
    label = "Cs-(Cds-O2d)(Cds-Cd)OsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   O2s u0 {1,S}
5   H   u0 {1,S}
6   C   u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)OsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1901,
    label = "Cs-(Cds-O2d)(Cds-Cds)OsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   O2s u0 {1,S}
5   H   u0 {1,S}
6   Cd  u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)CsOsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1902,
    label = "Cs-(Cds-O2d)(Cds-Cdd)OsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   O2s u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {3,D}
7   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cdd-Cd)OsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1903,
    label = "Cs-(Cds-O2d)(Cds-Cdd-O2d)OsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CO  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   O2s u0 {1,S}
6   H   u0 {1,S}
7   O2d u0 {3,D}
8   O2d u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)CsOsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1904,
    label = "Cs-(Cds-O2d)(Cds-Cdd-Cd)OsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CO  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   O2s u0 {1,S}
6   H   u0 {1,S}
7   O2d u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)OsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1905,
    label = "Cs-(Cds-Cd)(Cds-Cd)OsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   O2s u0 {1,S}
5   H   u0 {1,S}
6   C   u0 {2,D}
7   C   u0 {3,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([29.82,38.47,43.27,45.7,47.5,48.09,48.78],'J/(mol*K)','+|-',[3.64,3.64,3.64,3.64,3.64,3.64,3.64]),
        H298 = (-17.4,'kJ/mol','+|-',3.1),
        S298 = (-64.14,'J/(mol*K)','+|-',4.24),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 1906,
    label = "Cs-(Cds-Cds)(Cds-Cds)OsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   O2s u0 {1,S}
5   H   u0 {1,S}
6   Cd  u0 {2,D}
7   Cd  u0 {3,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.21,6.6,8.26,9.05,10.23,10.86,11.04],'cal/(mol*K)','+|-',[0.12,0.12,0.12,0.12,0.12,0.12,0.12]),
        H298 = (-6.67,'kcal/mol','+|-',0.24),
        S298 = (-10.42,'cal/(mol*K)','+|-',0.12),
    ),
    shortDesc = """Cs-OCdCdH BOZZELLI""",
    longDesc = 
"""

""",
)

entry(
    index = 1907,
    label = "Cs-(Cds-Cdd)(Cds-Cds)OsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   O2s u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {2,D}
7   Cd  u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cds)OsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1908,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cds)OsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   O2s u0 {1,S}
6   H   u0 {1,S}
7   Cd  u0 {3,D}
8   O2d u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)CsOsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1909,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cds)OsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   O2s u0 {1,S}
6   H   u0 {1,S}
7   Cd  u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)OsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1910,
    label = "Cs-(Cds-Cdd)(Cds-Cdd)OsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   O2s u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {2,D}
7   Cdd u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)OsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1911,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-O2d)OsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   O2s u0 {1,S}
7   H   u0 {1,S}
8   O2d u0 {4,D}
9   O2d u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)OsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1912,
    label = "Cs-(Cds-Cdd-O2d)(Cds-Cdd-Cd)OsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   O2s u0 {1,S}
7   H   u0 {1,S}
8   O2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)(Cds-Cds)OsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1913,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)OsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   O2s u0 {1,S}
7   H   u0 {1,S}
8   C   u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)OsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1914,
    label = "Cs-CtCsOsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Ct  u0 {1,S}
3   Cs  u0 {1,S}
4   O2s u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CsOsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1915,
    label = "Cs-CtCdsOsH",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Ct      u0 {1,S}
3   [Cd,CO] u0 {1,S}
4   O2s     u0 {1,S}
5   H       u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CtOsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1916,
    label = "Cs-(Cds-O2d)CtOsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   O2s u0 {1,S}
5   H   u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)OsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1917,
    label = "Cs-(Cds-Cd)CtOsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   O2s u0 {1,S}
5   H   u0 {1,S}
6   C   u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)CtOsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1918,
    label = "Cs-(Cds-Cds)CtOsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   O2s u0 {1,S}
5   H   u0 {1,S}
6   Cd  u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)OsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1919,
    label = "Cs-(Cds-Cdd)CtOsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   O2s u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)CtOsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1920,
    label = "Cs-(Cds-Cdd-O2d)CtOsH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Ct  u0 {1,S}
5   O2s u0 {1,S}
6   H   u0 {1,S}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)(Cds-Cds)OsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1921,
    label = "Cs-(Cds-Cdd-Cd)CtOsH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Ct  u0 {1,S}
5   O2s u0 {1,S}
6   H   u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)CtOsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1922,
    label = "Cs-CtCtOsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Ct  u0 {1,S}
3   Ct  u0 {1,S}
4   O2s u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)OsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1923,
    label = "Cs-CbCsOsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   Cs  u0 {1,S}
4   O2s u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.47,6.82,8.45,9.17,10.24,10.8,11.02],'cal/(mol*K)','+|-',[0.12,0.12,0.12,0.12,0.12,0.12,0.12]),
        H298 = (-6,'kcal/mol','+|-',0.24),
        S298 = (-11.1,'cal/(mol*K)','+|-',0.12),
    ),
    shortDesc = """Cs-OCbCsH BOZZELLI =3D C/Cd/C/H/O Jul 91""",
    longDesc = 
"""

""",
)

entry(
    index = 1924,
    label = "Cs-CbCdsOsH",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   Cb      u0 {1,S}
3   [Cd,CO] u0 {1,S}
4   O2s     u0 {1,S}
5   H       u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CbOsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1925,
    label = "Cs-(Cds-O2d)CbOsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   O2s u0 {1,S}
5   H   u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = 'Cs-(Cds-O2d)(Cds-Cds)OsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1926,
    label = "Cs-(Cds-Cd)CbOsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   O2s u0 {1,S}
5   H   u0 {1,S}
6   C   u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)CbOsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1927,
    label = "Cs-(Cds-Cds)CbOsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   O2s u0 {1,S}
5   H   u0 {1,S}
6   Cd  u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)OsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1928,
    label = "Cs-(Cds-Cdd)CbOsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   O2s u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)CbOsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1929,
    label = "Cs-(Cds-Cdd-O2d)CbOsH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   O2s u0 {1,S}
6   H   u0 {1,S}
7   O2d u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cdd-O2d)(Cds-Cds)OsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1930,
    label = "Cs-(Cds-Cdd-Cd)CbOsH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   O2s u0 {1,S}
6   H   u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)CbOsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1931,
    label = "Cs-CbCtOsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   Ct  u0 {1,S}
4   O2s u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)CtOsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1932,
    label = "Cs-CbCbOsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   Cb  u0 {1,S}
4   O2s u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)(Cds-Cds)OsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1933,
    label = "Cs-COsHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   C   u0 {1,S}
3   O2s u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = 'Cs-CsOsHH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1934,
    label = "Cs-CsOsHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cs  u0 {1,S}
3   O2s u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 1935,
    label = "Cs-CdsOsHH",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   [Cd,CO] u0 {1,S}
3   O2s     u0 {1,S}
4   H       u0 {1,S}
5   H       u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)OsHH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1936,
    label = "Cs-(Cds-O2d)OsHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   O2s u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([26.75,34.37,40.77,45.37,51.2,54.96,60.79],'J/(mol*K)','+|-',[4.34,4.34,4.34,4.34,4.34,4.34,4.34]),
        H298 = (-19.8,'kJ/mol','+|-',3.7),
        S298 = (31.54,'J/(mol*K)','+|-',5.06),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 1937,
    label = "Cs-(Cds-Cd)OsHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   O2s u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   C   u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([28.42,35.65,40.62,44.31,49.79,53.92,60.6],'J/(mol*K)','+|-',[3.38,3.38,3.38,3.38,3.38,3.38,3.38]),
        H298 = (-26.6,'kJ/mol','+|-',2.88),
        S298 = (34.59,'J/(mol*K)','+|-',3.95),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 1938,
    label = "Cs-(Cds-Cds)OsHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   O2s u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   Cd  u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.12,6.86,8.32,9.49,11.22,12.48,14.4],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (-6.76,'kcal/mol','+|-',0.2),
        S298 = (9.8,'cal/(mol*K)','+|-',0.1),
    ),
    shortDesc = """Cs-OCdHH BOZZELLI Hf PEDLEY c*ccoh C/C/Cd/H2""",
    longDesc = 
"""

""",
)

entry(
    index = 1939,
    label = "Cs-(Cds-Cdd)OsHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   O2s u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cdd-Cd)OsHH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1940,
    label = "Cs-(Cds-Cdd-O2d)OsHH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   O2s u0 {1,S}
5   H   u0 {1,S}
6   H   u0 {1,S}
7   O2d u0 {3,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.15,8.67,9.75,10.65,11.93,12.97,14.86],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (-8.68,'kcal/mol','+|-',0.2),
        S298 = (8.43,'cal/(mol*K)','+|-',0.1),
    ),
    shortDesc = """{C/CCO/O/H2} RAMAN & GREEN JPCA 2002, 106, 7937-7949""",
    longDesc = 
"""

""",
)

entry(
    index = 1941,
    label = "Cs-(Cds-Cdd-Cd)OsHH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   O2s u0 {1,S}
5   H   u0 {1,S}
6   H   u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = 'Cs-(Cds-Cds)OsHH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1942,
    label = "Cs-CtOsHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Ct  u0 {1,S}
3   O2s u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.12,6.86,8.32,9.49,11.22,12.48,14.4],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (-6.76,'kcal/mol','+|-',0.2),
        S298 = (9.8,'cal/(mol*K)','+|-',0.1),
    ),
    shortDesc = """Cs-OCtHH BOZZELLI assigned C/Cd/H2/O""",
    longDesc = 
"""

""",
)

entry(
    index = 1943,
    label = "Cs-CbOsHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   O2s u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cds)OsHH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1944,
    label = "Cs-CCCS",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   C  u0 {1,S}
4   C  u0 {1,S}
5   S  u0 {1,S}
""",
    thermo = 'Cs-CsCsCsS',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1945,
    label = "Cs-CsCsCsS",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cs u0 {1,S}
3   Cs u0 {1,S}
4   Cs u0 {1,S}
5   S  u0 {1,S}
""",
    thermo = 'Cs-CsCsCsS2',
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2018""",
    longDesc = 
"""
"
""",
)

entry(
    index = 1946,
    label = "Cs-CsCsCsS2",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cs  u0 {1,S}
3   Cs  u0 {1,S}
4   Cs  u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.97,6.93,7.58,8.09,10.09,10.99,9.6],'cal/(mol*K)'),
        H298 = (10.75,'kcal/mol'),
        S298 = (-31.6,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 1947,
    label = "Cs-CsCsCsS4",
    group = 
"""
1 * Cs                u0 {2,S} {3,S} {4,S} {5,S}
2   Cs                u0 {1,S}
3   Cs                u0 {1,S}
4   Cs                u0 {1,S}
5   [S4s,S4d,S4b,S4t] u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.57,5.49,5.72,6.11,8.54,9.77,8.79],'cal/(mol*K)'),
        H298 = (14.82,'kcal/mol'),
        S298 = (-29.47,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 1948,
    label = "Cs-CdsCsCsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S}
3   Cs  u0 {1,S}
4   Cs  u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1949,
    label = "Cs-(Cds-Cd)CsCsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   Cs  u0 {1,S}
5   S2s u0 {1,S}
6   C   u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1950,
    label = "Cs-(Cds-Cds)CsCsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   Cs  u0 {1,S}
5   S2s u0 {1,S}
6   Cd  u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1951,
    label = "Cs-(Cds-Cdd)CsCsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   Cs  u0 {1,S}
5   S2s u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1952,
    label = "Cs-(Cds-Cdd-S2d)CsCsSs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cs  u0 {1,S}
5   Cs  u0 {1,S}
6   S2s u0 {1,S}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1953,
    label = "Cs-(Cds-Cdd-Cd)CsCsSs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cs  u0 {1,S}
5   Cs  u0 {1,S}
6   S2s u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1954,
    label = "Cs-SsCtCsCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   S2s u0 {1,S}
3   Ct  u0 {1,S}
4   Cs  u0 {1,S}
5   Cs  u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1955,
    label = "Cs-CbCsCsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   Cs  u0 {1,S}
4   Cs  u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1956,
    label = "Cs-CdsCdsCsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S}
3   Cd  u0 {1,S}
4   Cs  u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1957,
    label = "Cs-(Cds-Cd)(Cds-Cd)CsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cs  u0 {1,S}
5   S2s u0 {1,S}
6   C   u0 {2,D}
7   C   u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1958,
    label = "Cs-(Cds-Cds)(Cds-Cds)CsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cs  u0 {1,S}
5   S2s u0 {1,S}
6   Cd  u0 {2,D}
7   Cd  u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1959,
    label = "Cs-(Cds-Cdd)(Cds-Cds)CsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cs  u0 {1,S}
5   S2s u0 {1,S}
6   Cdd u0 {2,D}
7   Cd  u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1960,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cds)CsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cs  u0 {1,S}
6   S2s u0 {1,S}
7   Cd  u0 {3,D}
8   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1961,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cds)CsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cs  u0 {1,S}
6   S2s u0 {1,S}
7   Cd  u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1962,
    label = "Cs-(Cds-Cdd)(Cds-Cdd)CsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cs  u0 {1,S}
5   S2s u0 {1,S}
6   Cdd u0 {2,D}
7   Cdd u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1963,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-S2d)CsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cs  u0 {1,S}
7   S2s u0 {1,S}
8   S2d u0 {4,D}
9   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1964,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-Cd)CsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cs  u0 {1,S}
7   S2s u0 {1,S}
8   S2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1965,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cs  u0 {1,S}
7   S2s u0 {1,S}
8   C   u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1966,
    label = "Cs-CtCdsCsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Ct  u0 {1,S}
3   Cd  u0 {1,S}
4   Cs  u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1967,
    label = "Cs-(Cds-Cd)CtCsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   Cs  u0 {1,S}
5   S2s u0 {1,S}
6   C   u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1968,
    label = "Cs-(Cds-Cds)CtCsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   Cs  u0 {1,S}
5   S2s u0 {1,S}
6   Cd  u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1969,
    label = "Cs-(Cds-Cdd)CtCsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   Cs  u0 {1,S}
5   S2s u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1970,
    label = "Cs-(Cds-Cdd-S2d)CtCsSs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Ct  u0 {1,S}
5   Cs  u0 {1,S}
6   S2s u0 {1,S}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1971,
    label = "Cs-(Cds-Cdd-Cd)CtCsSs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Ct  u0 {1,S}
5   Cs  u0 {1,S}
6   S2s u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1972,
    label = "Cs-CbCdsCsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   Cd  u0 {1,S}
4   Cs  u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1973,
    label = "Cs-(Cds-Cd)CbCsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Cs  u0 {1,S}
5   S2s u0 {1,S}
6   C   u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1974,
    label = "Cs-(Cds-Cds)CbCsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Cs  u0 {1,S}
5   S2s u0 {1,S}
6   Cd  u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1975,
    label = "Cs-(Cds-Cdd)CbCsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Cs  u0 {1,S}
5   S2s u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1976,
    label = "Cs-(Cds-Cdd-S2d)CbCsSs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Cs  u0 {1,S}
6   S2s u0 {1,S}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1977,
    label = "Cs-(Cds-Cdd-Cd)CbCsSs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Cs  u0 {1,S}
6   S2s u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1978,
    label = "Cs-CtCtCsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Ct  u0 {1,S}
3   Ct  u0 {1,S}
4   Cs  u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1979,
    label = "Cs-CbCtCsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   Ct  u0 {1,S}
4   Cs  u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1980,
    label = "Cs-CbCbCsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   Cb  u0 {1,S}
4   Cs  u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1981,
    label = "Cs-CdsCdsCdsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S}
3   Cd  u0 {1,S}
4   Cd  u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1982,
    label = "Cs-(Cds-Cd)(Cds-Cd)(Cds-Cd)S2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   S2s u0 {1,S}
6   C   u0 {2,D}
7   C   u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1983,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)S2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   S2s u0 {1,S}
6   Cd  u0 {2,D}
7   Cd  u0 {3,D}
8   Cd  u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1984,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd)S2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   S2s u0 {1,S}
6   Cd  u0 {2,D}
7   Cd  u0 {3,D}
8   Cdd u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1985,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-S2d)S2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   S2s u0 {1,S}
7   Cd  u0 {3,D}
8   Cd  u0 {4,D}
9   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1986,
    label = "Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)S2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   S2s u0 {1,S}
7   Cd  u0 {3,D}
8   Cd  u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1987,
    label = "Cs-(Cds-Cds)(Cds-Cdd)(Cds-Cdd)S2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   S2s u0 {1,S}
6   Cd  u0 {2,D}
7   Cdd u0 {3,D}
8   Cdd u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1988,
    label = "Cs-(Cds-Cds)(Cds-Cdd-S2d)(Cds-Cdd-S2d)S2s",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    S2s u0 {1,S}
8    Cd  u0 {4,D}
9    S2d u0 {5,D}
10   S2d u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1989,
    label = "Cs-(Cds-Cds)(Cds-Cdd-S2d)(Cds-Cdd-Cd)S2s",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    S2s u0 {1,S}
8    Cd  u0 {4,D}
9    S2d u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1990,
    label = "Cs-(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)S2s",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    S2s u0 {1,S}
8    Cd  u0 {4,D}
9    C   u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1991,
    label = "Cs-(Cds-Cdd)(Cds-Cdd)(Cds-Cdd)S2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cd  u0 {1,S} {8,D}
5   S2s u0 {1,S}
6   Cdd u0 {2,D}
7   Cdd u0 {3,D}
8   Cdd u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1992,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-S2d)(Cds-Cdd-S2d)S2s",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    S2s u0 {1,S}
9    S2d u0 {5,D}
10   S2d u0 {6,D}
11   S2d u0 {7,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1993,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-S2d)(Cds-Cdd-Cd)S2s",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    S2s u0 {1,S}
9    S2d u0 {5,D}
10   S2d u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1994,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-Cd)(Cds-Cdd-Cd)S2s",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    S2s u0 {1,S}
9    S2d u0 {5,D}
10   C   u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1995,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)S2s",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {8,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    Cd  u0 {1,S} {7,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    Cdd u0 {4,D} {11,D}
8    S2s u0 {1,S}
9    C   u0 {5,D}
10   C   u0 {6,D}
11   C   u0 {7,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1996,
    label = "Cs-CtCdsCdsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Ct  u0 {1,S}
3   Cd  u0 {1,S}
4   Cd  u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1997,
    label = "Cs-(Cds-Cd)(Cds-Cd)CtSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Ct  u0 {1,S}
5   S2s u0 {1,S}
6   C   u0 {2,D}
7   C   u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1998,
    label = "Cs-(Cds-Cds)(Cds-Cds)CtSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Ct  u0 {1,S}
5   S2s u0 {1,S}
6   Cd  u0 {2,D}
7   Cd  u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1999,
    label = "Cs-(Cds-Cdd)(Cds-Cds)CtSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Ct  u0 {1,S}
5   S2s u0 {1,S}
6   Cdd u0 {2,D}
7   Cd  u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2000,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cds)CtSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Ct  u0 {1,S}
6   S2s u0 {1,S}
7   Cd  u0 {3,D}
8   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2001,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cds)CtSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Ct  u0 {1,S}
6   S2s u0 {1,S}
7   Cd  u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2002,
    label = "Cs-(Cds-Cdd)(Cds-Cdd)CtSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Ct  u0 {1,S}
5   S2s u0 {1,S}
6   Cdd u0 {2,D}
7   Cdd u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2003,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-S2d)CtSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Ct  u0 {1,S}
7   S2s u0 {1,S}
8   S2d u0 {4,D}
9   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2004,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-Cd)CtSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Ct  u0 {1,S}
7   S2s u0 {1,S}
8   S2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2005,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CtSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Ct  u0 {1,S}
7   S2s u0 {1,S}
8   C   u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2006,
    label = "Cs-CbCdsCdsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   Cd  u0 {1,S}
4   Cd  u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2007,
    label = "Cs-(Cds-Cd)(Cds-Cd)CbSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cb  u0 {1,S}
5   S2s u0 {1,S}
6   C   u0 {2,D}
7   C   u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2008,
    label = "Cs-(Cds-Cds)(Cds-Cds)CbSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cb  u0 {1,S}
5   S2s u0 {1,S}
6   Cd  u0 {2,D}
7   Cd  u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2009,
    label = "Cs-(Cds-Cdd)(Cds-Cds)CbSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cb  u0 {1,S}
5   S2s u0 {1,S}
6   Cdd u0 {2,D}
7   Cd  u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2010,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cds)CbSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   S2s u0 {1,S}
7   Cd  u0 {3,D}
8   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2011,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cds)CbSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   S2s u0 {1,S}
7   Cd  u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2012,
    label = "Cs-(Cds-Cdd)(Cds-Cdd)CbSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   Cb  u0 {1,S}
5   S2s u0 {1,S}
6   Cdd u0 {2,D}
7   Cdd u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2013,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-S2d)CbSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cb  u0 {1,S}
7   S2s u0 {1,S}
8   S2d u0 {4,D}
9   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2014,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-Cd)CbSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cb  u0 {1,S}
7   S2s u0 {1,S}
8   S2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2015,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CbSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   Cb  u0 {1,S}
7   S2s u0 {1,S}
8   C   u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2016,
    label = "Cs-CtCtCdsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Ct  u0 {1,S}
3   Ct  u0 {1,S}
4   Cd  u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2017,
    label = "Cs-(Cds-Cd)CtCtSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   Ct  u0 {1,S}
5   S2s u0 {1,S}
6   C   u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2018,
    label = "Cs-(Cds-Cds)CtCtSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   Ct  u0 {1,S}
5   S2s u0 {1,S}
6   Cd  u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2019,
    label = "Cs-(Cds-Cdd)CtCtSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   Ct  u0 {1,S}
5   S2s u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2020,
    label = "Cs-(Cds-Cdd-S2d)CtCtSs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Ct  u0 {1,S}
5   Ct  u0 {1,S}
6   S2s u0 {1,S}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2021,
    label = "Cs-(Cds-Cdd-Cd)CtCtSs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Ct  u0 {1,S}
5   Ct  u0 {1,S}
6   S2s u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2022,
    label = "Cs-CbCtCdsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   Ct  u0 {1,S}
4   Cd  u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2023,
    label = "Cs-(Cds-Cd)CbCtSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Ct  u0 {1,S}
5   S2s u0 {1,S}
6   C   u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2024,
    label = "Cs-(Cds-Cds)CbCtSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Ct  u0 {1,S}
5   S2s u0 {1,S}
6   Cd  u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2025,
    label = "Cs-(Cds-Cdd)CbCtSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Ct  u0 {1,S}
5   S2s u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2026,
    label = "Cs-(Cds-Cdd-S2d)CbCtSs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Ct  u0 {1,S}
6   S2s u0 {1,S}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2027,
    label = "Cs-(Cds-Cdd-Cd)CbCtSs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Ct  u0 {1,S}
6   S2s u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2028,
    label = "Cs-CbCbCdsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   Cb  u0 {1,S}
4   Cd  u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2029,
    label = "Cs-(Cds-Cd)CbCbSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Cb  u0 {1,S}
5   S2s u0 {1,S}
6   C   u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2030,
    label = "Cs-(Cds-Cds)CbCbSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Cb  u0 {1,S}
5   S2s u0 {1,S}
6   Cd  u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2031,
    label = "Cs-(Cds-Cdd)CbCbSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Cb  u0 {1,S}
5   S2s u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2032,
    label = "Cs-(Cds-Cdd-S2d)CbCbSs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Cb  u0 {1,S}
6   S2s u0 {1,S}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2033,
    label = "Cs-(Cds-Cdd-Cd)CbCbSs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   Cb  u0 {1,S}
6   S2s u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2034,
    label = "Cs-CtCtCtSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Ct  u0 {1,S}
3   Ct  u0 {1,S}
4   Ct  u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2035,
    label = "Cs-CbCtCtSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   Ct  u0 {1,S}
4   Ct  u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2036,
    label = "Cs-CbCbCtSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   Cb  u0 {1,S}
4   Ct  u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2037,
    label = "Cs-CbCbCbSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   Cb  u0 {1,S}
4   Cb  u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2038,
    label = "Cs-C=SCbCsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Cs  u0 {1,S}
5   S2s u0 {1,S}
6   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2039,
    label = "Cs-C=SCsCsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   Cs  u0 {1,S}
5   S2s u0 {1,S}
6   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2040,
    label = "Cs-C=S(Cds-Cd)(Cds-Cd)S2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   S2s u0 {1,S}
6   C   u0 {3,D}
7   C   u0 {4,D}
8   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2041,
    label = "Cs-C=S(Cds-Cdd)(Cds-Cdd)S2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   S2s u0 {1,S}
6   Cdd u0 {3,D}
7   Cdd u0 {4,D}
8   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2042,
    label = "Cs-C=S(Cds-Cdd-Cd)(Cds-Cdd-Cd)S2s",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    CS  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    S2s u0 {1,S}
8    S2d u0 {4,D}
9    C   u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2043,
    label = "Cs-C=S(Cds-Cdd-S2d)(Cds-Cdd-Cd)S2s",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    CS  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    S2s u0 {1,S}
8    S2d u0 {4,D}
9    S2d u0 {5,D}
10   C   u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2044,
    label = "Cs-C=S(Cds-Cdd-S2d)(Cds-Cdd-S2d)S2s",
    group = 
"""
1  * Cs  u0 {2,S} {3,S} {4,S} {7,S}
2    Cd  u0 {1,S} {5,D}
3    Cd  u0 {1,S} {6,D}
4    CS  u0 {1,S} {8,D}
5    Cdd u0 {2,D} {9,D}
6    Cdd u0 {3,D} {10,D}
7    S2s u0 {1,S}
8    S2d u0 {4,D}
9    S2d u0 {5,D}
10   S2d u0 {6,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2045,
    label = "Cs-C=S(Cds-Cdd)(Cds-Cds)S2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   S2s u0 {1,S}
6   Cdd u0 {3,D}
7   Cd  u0 {4,D}
8   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2046,
    label = "Cs-C=S(Cds-Cdd-Cd)(Cds-Cds)S2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CS  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {7,D}
5   Cdd u0 {2,D} {9,D}
6   S2s u0 {1,S}
7   Cd  u0 {4,D}
8   S2d u0 {3,D}
9   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2047,
    label = "Cs-C=S(Cds-Cdd-S2d)(Cds-Cds)S2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CS  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {7,D}
5   Cdd u0 {2,D} {9,D}
6   S2s u0 {1,S}
7   Cd  u0 {4,D}
8   S2d u0 {3,D}
9   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2048,
    label = "Cs-C=S(Cds-Cds)(Cds-Cds)S2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {8,D}
3   Cd  u0 {1,S} {6,D}
4   Cd  u0 {1,S} {7,D}
5   S2s u0 {1,S}
6   Cd  u0 {3,D}
7   Cd  u0 {4,D}
8   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2049,
    label = "Cs-C=S(Cds-Cd)CtSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Ct  u0 {1,S}
5   S2s u0 {1,S}
6   C   u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2050,
    label = "Cs-C=S(Cds-Cds)CtSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Ct  u0 {1,S}
5   S2s u0 {1,S}
6   Cd  u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2051,
    label = "Cs-C=S(Cds-Cdd)CtSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Ct  u0 {1,S}
5   S2s u0 {1,S}
6   Cdd u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2052,
    label = "Cs-C=S(Cds-Cdd-S2d)CtSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CS  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Ct  u0 {1,S}
6   S2s u0 {1,S}
7   S2d u0 {3,D}
8   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2053,
    label = "Cs-C=S(Cds-Cdd-Cd)CtSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CS  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Ct  u0 {1,S}
6   S2s u0 {1,S}
7   S2d u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2054,
    label = "Cs-C=SCtCsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   Cs  u0 {1,S}
5   S2s u0 {1,S}
6   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2055,
    label = "Cs-C=SC=SC=SSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   CS  u0 {1,S} {7,D}
4   CS  u0 {1,S} {8,D}
5   S2s u0 {1,S}
6   S2d u0 {2,D}
7   S2d u0 {3,D}
8   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2056,
    label = "Cs-C=SC=S(Cds-Cd)S2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   CS  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {6,D}
5   S2s u0 {1,S}
6   C   u0 {4,D}
7   S2d u0 {2,D}
8   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2057,
    label = "Cs-C=SC=S(Cds-Cds)S2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   CS  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {6,D}
5   S2s u0 {1,S}
6   Cd  u0 {4,D}
7   S2d u0 {2,D}
8   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2058,
    label = "Cs-C=SC=S(Cds-Cdd)S2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   CS  u0 {1,S} {8,D}
4   Cd  u0 {1,S} {6,D}
5   S2s u0 {1,S}
6   Cdd u0 {4,D}
7   S2d u0 {2,D}
8   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2059,
    label = "Cs-C=SC=S(Cds-Cdd-S2d)S2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CS  u0 {1,S} {7,D}
4   CS  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   S2s u0 {1,S}
7   S2d u0 {3,D}
8   S2d u0 {4,D}
9   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2060,
    label = "Cs-C=SC=S(Cds-Cdd-Cd)S2s",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {6,S}
2   Cd  u0 {1,S} {5,D}
3   CS  u0 {1,S} {7,D}
4   CS  u0 {1,S} {8,D}
5   Cdd u0 {2,D} {9,D}
6   S2s u0 {1,S}
7   S2d u0 {3,D}
8   S2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2061,
    label = "Cs-C=SCbCbSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Cb  u0 {1,S}
5   S2s u0 {1,S}
6   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2062,
    label = "Cs-C=SC=SCbSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   CS  u0 {1,S} {7,D}
4   Cb  u0 {1,S}
5   S2s u0 {1,S}
6   S2d u0 {2,D}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2063,
    label = "Cs-C=SC=SCsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   CS  u0 {1,S} {7,D}
4   Cs  u0 {1,S}
5   S2s u0 {1,S}
6   S2d u0 {2,D}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2064,
    label = "Cs-C=SCtCtSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   Ct  u0 {1,S}
5   S2s u0 {1,S}
6   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2065,
    label = "Cs-C=S(Cds-Cd)CbSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   S2s u0 {1,S}
6   C   u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2066,
    label = "Cs-C=S(Cds-Cdd)CbSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   S2s u0 {1,S}
6   Cdd u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2067,
    label = "Cs-C=S(Cds-Cdd-Cd)CbSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CS  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   S2s u0 {1,S}
7   S2d u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2068,
    label = "Cs-C=S(Cds-Cdd-S2d)CbSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CS  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cb  u0 {1,S}
6   S2s u0 {1,S}
7   S2d u0 {3,D}
8   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2069,
    label = "Cs-C=S(Cds-Cds)CbSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cb  u0 {1,S}
5   S2s u0 {1,S}
6   Cd  u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2070,
    label = "Cs-C=SCbCtSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   Ct  u0 {1,S}
5   S2s u0 {1,S}
6   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2071,
    label = "Cs-C=SC=SCtSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   CS  u0 {1,S} {7,D}
4   Ct  u0 {1,S}
5   S2s u0 {1,S}
6   S2d u0 {2,D}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2072,
    label = "Cs-C=S(Cds-Cd)CsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cs  u0 {1,S}
5   S2s u0 {1,S}
6   C   u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2073,
    label = "Cs-C=S(Cds-Cds)CsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cs  u0 {1,S}
5   S2s u0 {1,S}
6   Cd  u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2074,
    label = "Cs-C=S(Cds-Cdd)CsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   Cs  u0 {1,S}
5   S2s u0 {1,S}
6   Cdd u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2075,
    label = "Cs-C=S(Cds-Cdd-S2d)CsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CS  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cs  u0 {1,S}
6   S2s u0 {1,S}
7   S2d u0 {3,D}
8   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2076,
    label = "Cs-C=S(Cds-Cdd-Cd)CsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CS  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   Cs  u0 {1,S}
6   S2s u0 {1,S}
7   S2d u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2077,
    label = "Cs-CCSS",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   C  u0 {1,S}
4   S  u0 {1,S}
5   S  u0 {1,S}
""",
    thermo = 'Cs-CsCsSS',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2078,
    label = "Cs-CsCsSS",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cs u0 {1,S}
3   Cs u0 {1,S}
4   S  u0 {1,S}
5   S  u0 {1,S}
""",
    thermo = 'Cs-CsCsS2S2',
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2018""",
    longDesc = 
"""
"
""",
)

entry(
    index = 2079,
    label = "Cs-CsCsS2S2",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cs  u0 {1,S}
3   Cs  u0 {1,S}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.83,10.39,10.63,10.76,13.38,14.21,10.84],'cal/(mol*K)'),
        H298 = (24.52,'kcal/mol'),
        S298 = (-26.78,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 2080,
    label = "Cs-CsCsS6S2",
    group = 
"""
1 * Cs                      u0 {2,S} {3,S} {4,S} {5,S}
2   Cs                      u0 {1,S}
3   Cs                      u0 {1,S}
4   [S6s,S6d,S6dd,S6t,S6td] u0 {1,S}
5   S2s                     u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.61,8.47,9.47,10.23,13.38,14.44,10.88],'cal/(mol*K)'),
        H298 = (20.26,'kcal/mol'),
        S298 = (-21.75,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 2081,
    label = "Cs-CdsCsSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S}
3   Cs  u0 {1,S}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2082,
    label = "Cs-(Cds-Cd)CsSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
6   C   u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2083,
    label = "Cs-(Cds-Cds)CsSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
6   Cd  u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2084,
    label = "Cs-(Cds-Cdd)CsSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2085,
    label = "Cs-(Cds-Cdd-S2d)CsSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cs  u0 {1,S}
5   S2s u0 {1,S}
6   S2s u0 {1,S}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2086,
    label = "Cs-(Cds-Cdd-Cd)CsSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cs  u0 {1,S}
5   S2s u0 {1,S}
6   S2s u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2087,
    label = "Cs-CdsCdsSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S}
3   Cd  u0 {1,S}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2088,
    label = "Cs-(Cds-Cd)(Cds-Cd)SsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
6   C   u0 {2,D}
7   C   u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2089,
    label = "Cs-(Cds-Cds)(Cds-Cds)SsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
6   Cd  u0 {2,D}
7   Cd  u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2090,
    label = "Cs-(Cds-Cdd)(Cds-Cds)SsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
6   Cdd u0 {2,D}
7   Cd  u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2091,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cds)SsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   S2s u0 {1,S}
6   S2s u0 {1,S}
7   Cd  u0 {3,D}
8   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2092,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cds)SsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   S2s u0 {1,S}
6   S2s u0 {1,S}
7   Cd  u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2093,
    label = "Cs-(Cds-Cdd)(Cds-Cdd)SsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
6   Cdd u0 {2,D}
7   Cdd u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2094,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-S2d)SsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   S2s u0 {1,S}
7   S2s u0 {1,S}
8   S2d u0 {4,D}
9   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2095,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-Cd)SsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   S2s u0 {1,S}
7   S2s u0 {1,S}
8   S2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2096,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)SsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   S2s u0 {1,S}
7   S2s u0 {1,S}
8   C   u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2097,
    label = "Cs-CtCsSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Ct  u0 {1,S}
3   Cs  u0 {1,S}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2098,
    label = "Cs-CtCdsSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Ct  u0 {1,S}
3   Cd  u0 {1,S}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2099,
    label = "Cs-(Cds-Cd)CtSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
6   C   u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2100,
    label = "Cs-(Cds-Cds)CtSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
6   Cd  u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2101,
    label = "Cs-(Cds-Cdd)CtSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2102,
    label = "Cs-(Cds-Cdd-S2d)CtSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Ct  u0 {1,S}
5   S2s u0 {1,S}
6   S2s u0 {1,S}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2103,
    label = "Cs-(Cds-Cdd-Cd)CtSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Ct  u0 {1,S}
5   S2s u0 {1,S}
6   S2s u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2104,
    label = "Cs-CtCtSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Ct  u0 {1,S}
3   Ct  u0 {1,S}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2105,
    label = "Cs-CbCsSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   Cs  u0 {1,S}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2106,
    label = "Cs-CbCdsSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   Cd  u0 {1,S}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2107,
    label = "Cs-(Cds-Cd)CbSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
6   C   u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2108,
    label = "Cs-(Cds-Cds)CbSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
6   Cd  u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2109,
    label = "Cs-(Cds-Cdd)CbSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2110,
    label = "Cs-(Cds-Cdd-S2d)CbSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   S2s u0 {1,S}
6   S2s u0 {1,S}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2111,
    label = "Cs-(Cds-Cdd-Cd)CbSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   S2s u0 {1,S}
6   S2s u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2112,
    label = "Cs-CbCtSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   Ct  u0 {1,S}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2113,
    label = "Cs-CbCbSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   Cb  u0 {1,S}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2114,
    label = "Cs-C=SCsSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
6   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2115,
    label = "Cs-C=S(Cds-Cd)SsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
6   C   u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2116,
    label = "Cs-C=S(Cds-Cdd)SsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
6   Cdd u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2117,
    label = "Cs-C=S(Cds-Cdd-Cd)SsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CS  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   S2s u0 {1,S}
6   S2s u0 {1,S}
7   S2d u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2118,
    label = "Cs-C=S(Cds-Cdd-S2d)SsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CS  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   S2s u0 {1,S}
6   S2s u0 {1,S}
7   S2d u0 {3,D}
8   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2119,
    label = "Cs-C=S(Cds-Cds)SsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
6   Cd  u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2120,
    label = "Cs-C=SC=SSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   CS  u0 {1,S} {7,D}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
6   S2d u0 {2,D}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2121,
    label = "Cs-C=SCbSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
6   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2122,
    label = "Cs-C=SCtSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
6   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2123,
    label = "Cs-CSsSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   C   u0 {1,S}
3   S2s u0 {1,S}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2124,
    label = "Cs-CsSsSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cs  u0 {1,S}
3   S2s u0 {1,S}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.76,12.55,12.86,12.95,16.67,17.86,12.77],'cal/(mol*K)'),
        H298 = (37.76,'kcal/mol'),
        S298 = (-23.39,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 2125,
    label = "Cs-CdsSsSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S}
3   S2s u0 {1,S}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2126,
    label = "Cs-(Cds-Cd)SsSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   S2s u0 {1,S}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
6   C   u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2127,
    label = "Cs-(Cds-Cds)SsSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   S2s u0 {1,S}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
6   Cd  u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2128,
    label = "Cs-(Cds-Cdd)SsSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   S2s u0 {1,S}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2129,
    label = "Cs-(Cds-Cdd-S2d)SsSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
6   S2s u0 {1,S}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2130,
    label = "Cs-(Cds-Cdd-Cd)SsSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
6   S2s u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2131,
    label = "Cs-CtSsSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Ct  u0 {1,S}
3   S2s u0 {1,S}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2132,
    label = "Cs-CbSsSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   S2s u0 {1,S}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2133,
    label = "Cs-C=SSsSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   S2s u0 {1,S}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
6   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2134,
    label = "Cs-SsSsSsSs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   S2s u0 {1,S}
3   S2s u0 {1,S}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.33,15.09,16.62,17.94,24.67,27.95,23.34],'cal/(mol*K)'),
        H298 = (33.21,'kcal/mol'),
        S298 = (-20.45,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 2135,
    label = "Cs-CSSH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   S  u0 {1,S}
4   S  u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = 'Cs-CsSSH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2136,
    label = "Cs-CsSSH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cs u0 {1,S}
3   S  u0 {1,S}
4   S  u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = 'Cs-CsS2S2H',
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2018""",
    longDesc = 
"""
"
""",
)

entry(
    index = 2137,
    label = "Cs-CsS2S2H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cs  u0 {1,S}
3   S2s u0 {1,S}
4   S2s u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.47,10.83,11.17,11.54,14.74,16.21,13.93],'cal/(mol*K)'),
        H298 = (22.59,'kcal/mol'),
        S298 = (-4.77,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 2138,
    label = "Cs-CsS4S2H",
    group = 
"""
1 * Cs                u0 {2,S} {3,S} {4,S} {5,S}
2   Cs                u0 {1,S}
3   [S4s,S4d,S4b,S4t] u0 {1,S}
4   S2s               u0 {1,S}
5   H                 u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.81,11.54,12.42,13.49,17.89,19.74,16.4],'cal/(mol*K)'),
        H298 = (18.39,'kcal/mol'),
        S298 = (-14.28,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 2139,
    label = "Cs-CsS6S2H",
    group = 
"""
1 * Cs                      u0 {2,S} {3,S} {4,S} {5,S}
2   Cs                      u0 {1,S}
3   [S6s,S6d,S6dd,S6t,S6td] u0 {1,S}
4   S2s                     u0 {1,S}
5   H                       u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.79,9.79,12.04,13.48,17.54,19.23,16.19],'cal/(mol*K)'),
        H298 = (18.74,'kcal/mol'),
        S298 = (-6.26,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 2140,
    label = "Cs-CdsSsSsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S}
3   S2s u0 {1,S}
4   S2s u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2141,
    label = "Cs-(Cds-Cd)SsSsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   S2s u0 {1,S}
4   S2s u0 {1,S}
5   H   u0 {1,S}
6   C   u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2142,
    label = "Cs-(Cds-Cds)SsSsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   S2s u0 {1,S}
4   S2s u0 {1,S}
5   H   u0 {1,S}
6   Cd  u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2143,
    label = "Cs-(Cds-Cdd)SsSsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   S2s u0 {1,S}
4   S2s u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2144,
    label = "Cs-(Cds-Cdd-S2d)SsSsH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
6   H   u0 {1,S}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2145,
    label = "Cs-(Cds-Cdd-Cd)SsSsH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   S2s u0 {1,S}
5   S2s u0 {1,S}
6   H   u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2146,
    label = "Cs-CtSsSsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Ct  u0 {1,S}
3   S2s u0 {1,S}
4   S2s u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2147,
    label = "Cs-CbSsSsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   S2s u0 {1,S}
4   S2s u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2148,
    label = "Cs-C=SSsSsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   S2s u0 {1,S}
4   S2s u0 {1,S}
5   H   u0 {1,S}
6   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2149,
    label = "Cs-CCSH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   C  u0 {1,S}
4   S  u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = 'Cs-CsCsSH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2150,
    label = "Cs-CsCsSH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cs u0 {1,S}
3   Cs u0 {1,S}
4   S  u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = 'Cs-CsCsS2H',
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2018""",
    longDesc = 
"""
"
""",
)

entry(
    index = 2151,
    label = "Cs-CsCsS2H",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cs  u0 {1,S}
3   Cs  u0 {1,S}
4   S2s u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.41,7.56,8.52,9.25,11.57,12.8,12.07],'cal/(mol*K)'),
        H298 = (11.07,'kcal/mol'),
        S298 = (-7.17,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 2152,
    label = "Cs-CsCsS4H",
    group = 
"""
1 * Cs                u0 {2,S} {3,S} {4,S} {5,S}
2   Cs                u0 {1,S}
3   Cs                u0 {1,S}
4   [S4s,S4d,S4b,S4t] u0 {1,S}
5   H                 u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.45,5.07,5.65,6.45,9.65,11.6,12.52],'cal/(mol*K)'),
        H298 = (12.82,'kcal/mol'),
        S298 = (-7.24,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 2153,
    label = "Cs-CsCsS6H",
    group = 
"""
1 * Cs                      u0 {2,S} {3,S} {4,S} {5,S}
2   Cs                      u0 {1,S}
3   Cs                      u0 {1,S}
4   [S6s,S6d,S6dd,S6t,S6td] u0 {1,S}
5   H                       u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.72,5.38,6.37,7.25,9.88,11.39,11.46],'cal/(mol*K)'),
        H298 = (5.98,'kcal/mol'),
        S298 = (-5.67,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 2154,
    label = "Cs-CdsCsSH",
    group = 
"""
1 * Cs      u0 {2,S} {3,S} {4,S} {5,S}
2   [Cd,CO] u0 {1,S}
3   Cs      u0 {1,S}
4   S       u0 {1,S}
5   H       u0 {1,S}
""",
    thermo = 'Cs-(Cds-Cd)CsSsH',
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2018""",
    longDesc = 
"""
"
""",
)

entry(
    index = 2155,
    label = "Cs-CdsCsS4H",
    group = 
"""
1 * Cs                u0 {2,S} {3,S} {4,S} {5,S}
2   [Cd,CO]           u0 {1,S}
3   Cs                u0 {1,S}
4   [S4s,S4d,S4b,S4t] u0 {1,S}
5   H                 u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.04343,9.7511,9.94209,10.4998,13.2832,14.9311,14.2721],'cal/(mol*K)'),
        H298 = (5.13823,'kcal/mol'),
        S298 = (-1.82552,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 2156,
    label = "Cs-(Cds-Cd)CsSsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   S2s u0 {1,S}
5   H   u0 {1,S}
6   C   u0 {2,D}
""",
    thermo = 'Cs-(Cds-Cds)CsSsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2157,
    label = "Cs-(Cds-Cds)CsSsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   S2s u0 {1,S}
5   H   u0 {1,S}
6   Cd  u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.97,10.76,11.55,11.83,13.17,13.89,12.78],'cal/(mol*K)'),
        H298 = (10.35,'kcal/mol'),
        S298 = (-9.71,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""

""",
)

entry(
    index = 2158,
    label = "Cs-(Cds-Cdd)CsSsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   S2s u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2159,
    label = "Cs-(Cds-Cdd-S2d)CsSsH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cs  u0 {1,S}
5   S2s u0 {1,S}
6   H   u0 {1,S}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2160,
    label = "Cs-(Cds-Cdd-Cd)CsSsH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cs  u0 {1,S}
5   S2s u0 {1,S}
6   H   u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2161,
    label = "Cs-CdsCdsSsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S}
3   Cd  u0 {1,S}
4   S2s u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2162,
    label = "Cs-(Cds-Cd)(Cds-Cd)SsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   S2s u0 {1,S}
5   H   u0 {1,S}
6   C   u0 {2,D}
7   C   u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2163,
    label = "Cs-(Cds-Cds)(Cds-Cds)SsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   S2s u0 {1,S}
5   H   u0 {1,S}
6   Cd  u0 {2,D}
7   Cd  u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2164,
    label = "Cs-(Cds-Cdd)(Cds-Cds)SsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   S2s u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {2,D}
7   Cd  u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2165,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cds)SsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   S2s u0 {1,S}
6   H   u0 {1,S}
7   Cd  u0 {3,D}
8   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2166,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cds)SsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   S2s u0 {1,S}
6   H   u0 {1,S}
7   Cd  u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2167,
    label = "Cs-(Cds-Cdd)(Cds-Cdd)SsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cd  u0 {1,S} {7,D}
4   S2s u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {2,D}
7   Cdd u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2168,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-S2d)SsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   S2s u0 {1,S}
7   H   u0 {1,S}
8   S2d u0 {4,D}
9   S2d u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2169,
    label = "Cs-(Cds-Cdd-S2d)(Cds-Cdd-Cd)SsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   S2s u0 {1,S}
7   H   u0 {1,S}
8   S2d u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2170,
    label = "Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)SsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {6,S} {7,S}
2   Cd  u0 {1,S} {4,D}
3   Cd  u0 {1,S} {5,D}
4   Cdd u0 {2,D} {8,D}
5   Cdd u0 {3,D} {9,D}
6   S2s u0 {1,S}
7   H   u0 {1,S}
8   C   u0 {4,D}
9   C   u0 {5,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2171,
    label = "Cs-CtCsSsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Ct  u0 {1,S}
3   Cs  u0 {1,S}
4   S2s u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.62,7.79,8.7,9.44,11.69,12.98,11.46],'cal/(mol*K)'),
        H298 = (13.55,'kcal/mol'),
        S298 = (-6.16,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 2172,
    label = "Cs-CtCdsSsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Ct  u0 {1,S}
3   Cd  u0 {1,S}
4   S2s u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2173,
    label = "Cs-(Cds-Cd)CtSsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   S2s u0 {1,S}
5   H   u0 {1,S}
6   C   u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2174,
    label = "Cs-(Cds-Cds)CtSsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   S2s u0 {1,S}
5   H   u0 {1,S}
6   Cd  u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2175,
    label = "Cs-(Cds-Cdd)CtSsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   S2s u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2176,
    label = "Cs-(Cds-Cdd-S2d)CtSsH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Ct  u0 {1,S}
5   S2s u0 {1,S}
6   H   u0 {1,S}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2177,
    label = "Cs-(Cds-Cdd-Cd)CtSsH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Ct  u0 {1,S}
5   S2s u0 {1,S}
6   H   u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2178,
    label = "Cs-CtCtSsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Ct  u0 {1,S}
3   Ct  u0 {1,S}
4   S2s u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2179,
    label = "Cs-CbCsSsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   Cs  u0 {1,S}
4   S2s u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.44,8.88,9.92,10.57,12.5,13.34,12.17],'cal/(mol*K)'),
        H298 = (11.58,'kcal/mol'),
        S298 = (-9,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 2180,
    label = "Cs-CbCdsSsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   Cd  u0 {1,S}
4   S2s u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2181,
    label = "Cs-(Cds-Cd)CbSsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   S2s u0 {1,S}
5   H   u0 {1,S}
6   C   u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2182,
    label = "Cs-(Cds-Cds)CbSsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   S2s u0 {1,S}
5   H   u0 {1,S}
6   Cd  u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2183,
    label = "Cs-(Cds-Cdd)CbSsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   S2s u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2184,
    label = "Cs-(Cds-Cdd-S2d)CbSsH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   S2s u0 {1,S}
6   H   u0 {1,S}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2185,
    label = "Cs-(Cds-Cdd-Cd)CbSsH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   Cb  u0 {1,S}
5   S2s u0 {1,S}
6   H   u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2186,
    label = "Cs-CbCtSsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   Ct  u0 {1,S}
4   S2s u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2187,
    label = "Cs-CbCbSsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   Cb  u0 {1,S}
4   S2s u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2188,
    label = "Cs-C=SCbSsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   Cb  u0 {1,S}
4   S2s u0 {1,S}
5   H   u0 {1,S}
6   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2189,
    label = "Cs-C=SC=SSsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   CS  u0 {1,S} {7,D}
4   S2s u0 {1,S}
5   H   u0 {1,S}
6   S2d u0 {2,D}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2190,
    label = "Cs-C=SCsSsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   S2s u0 {1,S}
5   H   u0 {1,S}
6   S2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.76,10.49,12.16,12.69,14.03,14.39,12.81],'cal/(mol*K)'),
        H298 = (12.61,'kcal/mol'),
        S298 = (-10.23,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 2191,
    label = "Cs-C=SCtSsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   Ct  u0 {1,S}
4   S2s u0 {1,S}
5   H   u0 {1,S}
6   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2192,
    label = "Cs-C=S(Cds-Cd)SsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   S2s u0 {1,S}
5   H   u0 {1,S}
6   C   u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2193,
    label = "Cs-C=S(Cds-Cdd)SsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   S2s u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2194,
    label = "Cs-C=S(Cds-Cdd-Cd)SsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CS  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   S2s u0 {1,S}
6   H   u0 {1,S}
7   S2d u0 {3,D}
8   C   u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2195,
    label = "Cs-C=S(Cds-Cdd-S2d)SsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {5,S} {6,S}
2   Cd  u0 {1,S} {4,D}
3   CS  u0 {1,S} {7,D}
4   Cdd u0 {2,D} {8,D}
5   S2s u0 {1,S}
6   H   u0 {1,S}
7   S2d u0 {3,D}
8   S2d u0 {4,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2196,
    label = "Cs-C=S(Cds-Cds)SsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {7,D}
3   Cd  u0 {1,S} {6,D}
4   S2s u0 {1,S}
5   H   u0 {1,S}
6   Cd  u0 {3,D}
7   S2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2197,
    label = "Cs-CSHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   S  u0 {1,S}
4   H  u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = 'Cs-CsSHH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2198,
    label = "Cs-CsSHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   Cs u0 {1,S}
3   S  u0 {1,S}
4   H  u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = 'Cs-CsS2HH',
    shortDesc = """""",
    longDesc = 
"""
"
""",
)

entry(
    index = 2199,
    label = "Cs-CsS2HH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cs  u0 {1,S}
3   S2s u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.96,8.17,9.34,10.35,13.25,15.03,15.38],'cal/(mol*K)'),
        H298 = (6.95,'kcal/mol'),
        S298 = (14.5,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 2200,
    label = "Cs-CsS4HH",
    group = 
"""
1 * Cs                u0 {2,S} {3,S} {4,S} {5,S}
2   Cs                u0 {1,S}
3   [S4s,S4d,S4b,S4t] u0 {1,S}
4   H                 u0 {1,S}
5   H                 u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.12,6.09,6.9,7.99,11.61,14.01,15.09],'cal/(mol*K)'),
        H298 = (9.13,'kcal/mol'),
        S298 = (12.35,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 2201,
    label = "Cs-CsS6HH",
    group = 
"""
1 * Cs                      u0 {2,S} {3,S} {4,S} {5,S}
2   Cs                      u0 {1,S}
3   [S6s,S6d,S6dd,S6t,S6td] u0 {1,S}
4   H                       u0 {1,S}
5   H                       u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.91,5.81,7.26,8.55,11.77,13.8,14.82],'cal/(mol*K)'),
        H298 = (3.84,'kcal/mol'),
        S298 = (19.66,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 2202,
    label = "Cs-CdsSsHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S}
3   S2s u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.96,10.86,11.45,12.01,14.21,15.74,15.89],'cal/(mol*K)'),
        H298 = (7.47,'kcal/mol'),
        S298 = (12.31,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 2203,
    label = "Cs-(Cds-Cd)SsHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   S2s u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   C   u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2204,
    label = "Cs-(Cds-Cds)SsHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   S2s u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   Cd  u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2205,
    label = "Cs-(Cds-Cdd)SsHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   S2s u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   Cdd u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2206,
    label = "Cs-(Cds-Cdd-S2d)SsHH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   S2s u0 {1,S}
5   H   u0 {1,S}
6   H   u0 {1,S}
7   S2d u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2207,
    label = "Cs-(Cds-Cdd-Cd)SsHH",
    group = 
"""
1 * Cs  u0 {2,S} {4,S} {5,S} {6,S}
2   Cd  u0 {1,S} {3,D}
3   Cdd u0 {2,D} {7,D}
4   S2s u0 {1,S}
5   H   u0 {1,S}
6   H   u0 {1,S}
7   C   u0 {3,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2208,
    label = "Cs-CtSsHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Ct  u0 {1,S}
3   S2s u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.95,8.29,9.45,10.53,13.39,15.28,14.76],'cal/(mol*K)'),
        H298 = (10.19,'kcal/mol'),
        S298 = (15.23,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 2209,
    label = "Cs-CbSsHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cb  u0 {1,S}
3   S2s u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.54,8.82,10.03,11.08,13.8,15.37,15.33],'cal/(mol*K)'),
        H298 = (6.16,'kcal/mol'),
        S298 = (12.86,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 2210,
    label = "Cs-C=SSsHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CS  u0 {1,S} {6,D}
3   S2s u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   S2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.7,11.06,12.6,13.3,15.32,16.4,15.98],'cal/(mol*K)'),
        H298 = (10.03,'kcal/mol'),
        S298 = (11.36,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 2211,
    label = "Cs-CIHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   C   u0 {1,S}
3   I1s u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.2,11,12.9,13.9,15.8,17.2,18.6],'cal/(mol*K)'),
        H298 = (8,'kcal/mol'),
        S298 = (43,'cal/(mol*K)'),
    ),
    shortDesc = """C-(I)(H)2(C) BENSON""",
    longDesc = 
"""
Thermochemical Kinetics 2nd Ed., by Sidney Benson (Table A4, p.280)
Cpdata at 1500K = Cpdata at 1000K + 1.4
""",
)

entry(
    index = 2212,
    label = "Cs-CIIH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   C   u0 {1,S}
3   I1s u0 {1,S}
4   I1s u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([12.2,13.4,16.4,17,17.7,18.4,19.1],'cal/(mol*K)'),
        H298 = (26,'kcal/mol'),
        S298 = (54.6,'cal/(mol*K)'),
    ),
    shortDesc = """C-(I)2(C)(H) BENSON""",
    longDesc = 
"""
Thermochemical Kinetics 2nd Ed., by Sidney Benson (Table A4, p.280)
Cpdata from 600 to 1500K estimated (base on entry 2088)
""",
)

entry(
    index = 2213,
    label = "Cs-CCIH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   C   u0 {1,S}
3   C   u0 {1,S}
4   I1s u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.2,10.9,12.2,13,14.2,14.8,15.4],'cal/(mol*K)'),
        H298 = (10.5,'kcal/mol'),
        S298 = (21.3,'cal/(mol*K)'),
    ),
    shortDesc = """C-(I)(H)(C)2 BENSON""",
    longDesc = 
"""
Thermochemical Kinetics 2nd Ed., by Sidney Benson (Table A4, p.280)
Cpdata at 1500K = Cpdata at 1000K + 0.6
""",
)

entry(
    index = 2214,
    label = "Cs-CCCI",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   C   u0 {1,S}
3   C   u0 {1,S}
4   C   u0 {1,S}
5   I1s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.7,11.4,12.7,13.9,14.7,15.3,15.9],'cal/(mol*K)'),
        H298 = (13,'kcal/mol'),
        S298 = (0,'cal/(mol*K)'),
    ),
    shortDesc = """C-(I)(C)3  BENSON""",
    longDesc = 
"""
Thermochemical Kinetics 2nd Ed., by Sidney Benson (Table A4, p.280)
Cpdata from 400 to 1500K estimated (base on entry 2092)
""",
)

entry(
    index = 2215,
    label = "Cs-HHNN",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   N  u0 {1,S}
3   N  u0 {1,S}
4   H  u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.40223,9.82386,10.6253,10.9822,11.6154,12.046,13.1963],'cal/(mol*K)','+|-',[0.634015,0.677118,0.698323,0.693364,0.674477,0.662375,0.67587]),
        H298 = (-3.98748,'kcal/mol','+|-',2.54122),
        S298 = (8.38663,'cal/(mol*K)','+|-',1.94699),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library                 | Number of Species
CHON_G4                 |         62
thermo_DFT_CCSDTF12_BAC |         1
""",
)

entry(
    index = 2216,
    label = "Cs-HNNN",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   N  u0 {1,S}
3   N  u0 {1,S}
4   N  u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.71645,11.4005,11.7334,11.2695,10.4835,9.76382,9.30964],'cal/(mol*K)','+|-',[1.01212,1.08093,1.11478,1.10686,1.07671,1.05739,1.07893]),
        H298 = (-3.90546,'kcal/mol','+|-',4.05672),
        S298 = (-13.0866,'cal/(mol*K)','+|-',3.1081),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHON_G4  |         7
BurcatNS |         1
""",
)

entry(
    index = 2217,
    label = "Cs-NNNN",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   N  u0 {1,S}
3   N  u0 {1,S}
4   N  u0 {1,S}
5   N  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([10.9928,11.6289,10.9865,9.75775,8.38313,7.18791,5.4612],'cal/(mol*K)','+|-',[1.75746,1.87694,1.93572,1.92197,1.86962,1.83607,1.87348]),
        H298 = (29.2456,'kcal/mol','+|-',7.04415),
        S298 = (-46.5692,'cal/(mol*K)','+|-',5.39696),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library                 | Number of Species
thermo_DFT_CCSDTF12_BAC |         1
""",
)

entry(
    index = 2218,
    label = "Cs-N3sdN3sdN3sdN3sd",
    group = 
"""
1 * Cs        u0 {2,S} {3,S} {4,S} {5,S}
2   [N3s,N3d] u0 {1,S}
3   [N3s,N3d] u0 {1,S}
4   [N3s,N3d] u0 {1,S}
5   [N3s,N3d] u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.67635,11.0833,11.8998,11.2529,9.74518,8.2474,6.55188],'cal/(mol*K)','+|-',[1.6183,1.72832,1.78245,1.76979,1.72158,1.69069,1.72514]),
        H298 = (-9.71968,'kcal/mol','+|-',6.48639),
        S298 = (-36.82,'cal/(mol*K)','+|-',4.96963),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         1
""",
)

entry(
    index = 2219,
    label = "Cs-HHNO",
    group = 
"""
1 * Cs         u0 {2,S} {3,S} {4,S} {5,S}
2   [O2s,O0sc] u0 {1,S}
3   N          u0 {1,S}
4   H          u0 {1,S}
5   H          u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.93262,8.88412,10.1865,11.0276,12.3018,13.0993,14.3444],'cal/(mol*K)','+|-',[0.332453,0.355054,0.366173,0.363573,0.353669,0.347323,0.354399]),
        H298 = (-9.28074,'kcal/mol','+|-',1.33252),
        S298 = (6.33483,'cal/(mol*K)','+|-',1.02092),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         64
""",
)

entry(
    index = 2220,
    label = "Cs-HHO(NO)",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   N   u0 {1,S} {6,D}
3   O2s u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.33372,9.50589,11.53,12.6936,14.2205,14.9164,15.5717],'cal/(mol*K)','+|-',[0.679153,0.725324,0.748039,0.742726,0.722495,0.709531,0.723987]),
        H298 = (-9.71782,'kcal/mol','+|-',2.72214),
        S298 = (1.91765,'cal/(mol*K)','+|-',2.0856),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         3
""",
)

entry(
    index = 2221,
    label = "Cs-HNNO",
    group = 
"""
1 * Cs         u0 {2,S} {3,S} {4,S} {5,S}
2   [O2s,O0sc] u0 {1,S}
3   N          u0 {1,S}
4   N          u0 {1,S}
5   H          u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.00005,10.2827,11.365,11.5557,11.5141,11.1548,10.7162],'cal/(mol*K)','+|-',[0.716702,0.765426,0.789396,0.78379,0.76244,0.748759,0.764014]),
        H298 = (-11.932,'kcal/mol','+|-',2.87264),
        S298 = (-17.367,'cal/(mol*K)','+|-',2.20091),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         9
""",
)

entry(
    index = 2222,
    label = "Cs-NNNO",
    group = 
"""
1 * Cs         u0 {2,S} {3,S} {4,S} {5,S}
2   [O2s,O0sc] u0 {1,S}
3   N          u0 {1,S}
4   N          u0 {1,S}
5   N          u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([10.0712,13.0847,13.9882,13.1906,11.1778,9.1615,6.15878],'cal/(mol*K)','+|-',[1.36986,1.46298,1.5088,1.49808,1.45728,1.43113,1.46029]),
        H298 = (-15.7945,'kcal/mol','+|-',5.49057),
        S298 = (-33.7236,'cal/(mol*K)','+|-',4.20667),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         1
""",
)

entry(
    index = 2223,
    label = "Cs-HNOO",
    group = 
"""
1 * Cs         u0 {2,S} {3,S} {4,S} {5,S}
2   [O2s,O0sc] u0 {1,S}
3   [O2s,O0sc] u0 {1,S}
4   N          u0 {1,S}
5   H          u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.31027,8.8625,10.3772,11.1993,12.1525,12.515,12.6237],'cal/(mol*K)','+|-',[0.481704,0.514452,0.530563,0.526795,0.512445,0.50325,0.513503]),
        H298 = (-19.4136,'kcal/mol','+|-',1.93074),
        S298 = (-19.4679,'cal/(mol*K)','+|-',1.47926),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         7
""",
)

entry(
    index = 2224,
    label = "Cs-HOO(NO)",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   N   u0 {1,S} {6,D}
3   O2s u0 {1,S}
4   O2s u0 {1,S}
5   H   u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.00727,8.3137,9.73942,10.6081,11.6874,12.1759,12.5726],'cal/(mol*K)','+|-',[1.03564,1.10605,1.14069,1.13259,1.10174,1.08197,1.10401]),
        H298 = (-15.2237,'kcal/mol','+|-',4.15101),
        S298 = (-20.5641,'cal/(mol*K)','+|-',3.18035),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         1
""",
)

entry(
    index = 2225,
    label = "Cs-NNOO",
    group = 
"""
1 * Cs         u0 {2,S} {3,S} {4,S} {5,S}
2   [O2s,O0sc] u0 {1,S}
3   [O2s,O0sc] u0 {1,S}
4   N          u0 {1,S}
5   N          u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.80657,9.33604,10.2399,10.324,10.1484,9.64457,8.79236],'cal/(mol*K)','+|-',[1.16026,1.23914,1.27795,1.26887,1.23431,1.21216,1.23686]),
        H298 = (-23.6222,'kcal/mol','+|-',4.6505),
        S298 = (-39.4234,'cal/(mol*K)','+|-',3.56304),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         1
""",
)

entry(
    index = 2226,
    label = "Cs-CHHN",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   N  u0 {1,S}
4   H  u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.4382,8.81296,9.71675,10.3283,11.3268,12.0934,13.3355],'cal/(mol*K)','+|-',[0.352863,0.376852,0.388654,0.385894,0.375382,0.368647,0.376157]),
        H298 = (-2.1676,'kcal/mol','+|-',1.41433),
        S298 = (8.99137,'cal/(mol*K)','+|-',1.0836),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHON_G4  |         39
BurcatNS |         1
""",
)

entry(
    index = 2227,
    label = "Cs-N5dcCsHH",
    group = 
"""
1 * Cs   u0 {2,S} {3,S} {4,S} {5,S}
2   N5dc u0 {1,S}
3   Cs   u0 {1,S}
4   H    u0 {1,S}
5   H    u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'cal/(mol*K)'),
        H298 = (0,'kcal/mol'),
        S298 = (0,'cal/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2228,
    label = "Cs-(N5dcOdO0sc)CsHH",
    group = 
"""
1 * Cs   u0 {2,S} {3,S} {4,S} {5,S}
2   N5dc u0 {1,S} {6,D} {7,S}
3   Cs   u0 {1,S}
4   H    u0 {1,S}
5   H    u0 {1,S}
6   O2d  u0 {2,D}
7   O0sc u0 {2,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.11334,7.44242,8.5783,9.407,10.6291,11.4782,12.6045],'cal/(mol*K)','+|-',[0.417652,0.446046,0.460014,0.456748,0.444306,0.436334,0.445223]),
        H298 = (-10.195,'kcal/mol','+|-',1.67401),
        S298 = (8.12434,'cal/(mol*K)','+|-',1.28256),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library                 | Number of Species
BurcatNS                |         2
thermo_DFT_CCSDTF12_BAC |         2
CHON                    |         14
""",
)

entry(
    index = 2229,
    label = "Cs-N3dCsHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   N3d u0 {1,S}
3   Cs  u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.7856,6.93263,7.88967,8.58557,9.80276,10.8102,12.6648],'cal/(mol*K)','+|-',[0.46392,0.495459,0.510974,0.507346,0.493526,0.48467,0.494545]),
        H298 = (-3.74819,'kcal/mol','+|-',1.85946),
        S298 = (8.21727,'cal/(mol*K)','+|-',1.42464),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHON_G4  |         6
BurcatNS |         1
CHON     |         2
""",
)

entry(
    index = 2230,
    label = "Cs-(N3dCd)CsHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   N3d u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   Cd  u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.82076,8.35083,9.60684,10.4892,11.993,12.9367,14.3857],'cal/(mol*K)','+|-',[0.584556,0.624296,0.643846,0.639274,0.621861,0.610702,0.623145]),
        H298 = (-2.22016,'kcal/mol','+|-',2.34298),
        S298 = (9.49076,'cal/(mol*K)','+|-',1.7951),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         7
""",
)

entry(
    index = 2231,
    label = "Cs-(N3dN3d)CsHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   N3d u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   N3d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.48953,11.0042,12.2879,13.473,15.4766,16.922,19.0554],'cal/(mol*K)','+|-',[0.37976,0.405577,0.418278,0.415308,0.403995,0.396746,0.404829]),
        H298 = (22.42,'kcal/mol','+|-',1.52213),
        S298 = (16.1014,'cal/(mol*K)','+|-',1.1662),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         7
""",
)

entry(
    index = 2232,
    label = "Cs-N3dCOHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   N3d u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.06269,8.55996,9.68665,10.5695,12.1001,13.1433,14.8482],'cal/(mol*K)','+|-',[0.764133,0.816081,0.841638,0.835661,0.812898,0.798312,0.814576]),
        H298 = (-2.15801,'kcal/mol','+|-',3.06275),
        S298 = (5.43806,'cal/(mol*K)','+|-',2.34656),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         2
""",
)

entry(
    index = 2233,
    label = "Cs-N3dCdHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   N3d u0 {1,S}
3   Cd  u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.41722,7.75772,8.82571,9.59371,10.8786,11.865,13.5253],'cal/(mol*K)','+|-',[0.589366,0.629433,0.649145,0.644535,0.626978,0.615728,0.628273]),
        H298 = (-3.41688,'kcal/mol','+|-',2.36226),
        S298 = (9.79685,'cal/(mol*K)','+|-',1.80988),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         4
""",
)

entry(
    index = 2234,
    label = "Cs-N3dCtHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   N3d u0 {1,S}
3   Ct  u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.60067,8.07316,9.0705,9.78654,10.9804,11.8935,13.0227],'cal/(mol*K)','+|-',[0.551417,0.588905,0.607347,0.603034,0.586607,0.576082,0.587819]),
        H298 = (-3.48194,'kcal/mol','+|-',2.21016),
        S298 = (10.3307,'cal/(mol*K)','+|-',1.69334),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         5
""",
)

entry(
    index = 2235,
    label = "Cs-NCsHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   N  u0 {1,S}
3   Cs u0 {1,S}
4   H  u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.45231,8.69775,9.47901,10.0354,11.1199,11.9231,13.4181],'cal/(mol*K)','+|-',[0.586948,0.626852,0.646482,0.641891,0.624406,0.613202,0.625695]),
        H298 = (-6.25629,'kcal/mol','+|-',2.35257),
        S298 = (7.54475,'cal/(mol*K)','+|-',1.80245),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         4
""",
)

entry(
    index = 2236,
    label = "Cs-N3sCsHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   N3s u0 {1,S}
3   Cs  u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.72324,7.93026,8.97178,9.74895,10.9321,11.7571,13.1277],'cal/(mol*K)','+|-',[0.323601,0.3456,0.356423,0.353892,0.344252,0.338075,0.344963]),
        H298 = (-2.61142,'kcal/mol','+|-',1.29704),
        S298 = (9.83938,'cal/(mol*K)','+|-',0.993741),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHON_G4  |         66
BurcatNS |         1
CHN      |         43
CHON     |         9
""",
)

entry(
    index = 2237,
    label = "Cs-CHNN",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   N  u0 {1,S}
4   N  u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.02857,10.2939,10.5221,10.2699,9.84078,9.49988,9.28489],'cal/(mol*K)','+|-',[0.809959,0.865024,0.892113,0.885777,0.861649,0.846188,0.863428]),
        H298 = (0.393621,'kcal/mol','+|-',3.24643),
        S298 = (-12.1951,'cal/(mol*K)','+|-',2.48729),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         4
""",
)

entry(
    index = 2238,
    label = "Cs-NNCsH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   N  u0 {1,S}
3   N  u0 {1,S}
4   Cs u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.24166,9.33161,9.70182,9.63398,9.58224,9.4462,9.60581],'cal/(mol*K)','+|-',[0.693821,0.74099,0.764194,0.758767,0.738099,0.724855,0.739623]),
        H298 = (-2.52448,'kcal/mol','+|-',2.78093),
        S298 = (-13.9066,'cal/(mol*K)','+|-',2.13065),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHON_G4  |         11
BurcatNS |         1
""",
)

entry(
    index = 2239,
    label = "Cs-CNNN",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   N  u0 {1,S}
4   N  u0 {1,S}
5   N  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([13.1708,14.1044,13.7412,12.7575,11.4387,10.1513,8.20207],'cal/(mol*K)','+|-',[1.18634,1.26699,1.30667,1.29739,1.26205,1.2394,1.26466]),
        H298 = (19.9029,'kcal/mol','+|-',4.75502),
        S298 = (-38.6482,'cal/(mol*K)','+|-',3.64312),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHON_G4  |         1
BurcatNS |         1
""",
)

entry(
    index = 2240,
    label = "Cs-CN3dsN3dsN3ds",
    group = 
"""
1 * Cs        u0 {2,S} {3,S} {4,S} {5,S}
2   C         u0 {1,S}
3   [N3s,N3d] u0 {1,S}
4   [N3s,N3d] u0 {1,S}
5   [N3s,N3d] u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.64506,11.964,12.4292,11.52,9.69687,7.91876,5.4402],'cal/(mol*K)','+|-',[1.36986,1.46298,1.5088,1.49808,1.45728,1.43113,1.46029]),
        H298 = (-4.1318,'kcal/mol','+|-',5.49057),
        S298 = (-32.8429,'cal/(mol*K)','+|-',4.20667),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         1
""",
)

entry(
    index = 2241,
    label = "Cs-CHNO",
    group = 
"""
1 * Cs         u0 {2,S} {3,S} {4,S} {5,S}
2   C          u0 {1,S}
3   [O2s,O0sc] u0 {1,S}
4   N          u0 {1,S}
5   H          u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.05675,8.788,9.7105,10.1547,10.7165,10.858,10.8637],'cal/(mol*K)','+|-',[0.398794,0.425905,0.439243,0.436124,0.424244,0.416632,0.42512]),
        H298 = (-8.01295,'kcal/mol','+|-',1.59842),
        S298 = (-15.8961,'cal/(mol*K)','+|-',1.22465),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         16
""",
)

entry(
    index = 2242,
    label = "Cs-CNNO",
    group = 
"""
1 * Cs         u0 {2,S} {3,S} {4,S} {5,S}
2   C          u0 {1,S}
3   [O2s,O0sc] u0 {1,S}
4   N          u0 {1,S}
5   N          u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.39017,11.288,12.36,12.0313,10.9456,9.62469,7.37044],'cal/(mol*K)','+|-',[1.16026,1.23914,1.27795,1.26887,1.23431,1.21216,1.23686]),
        H298 = (-13.077,'kcal/mol','+|-',4.6505),
        S298 = (-38.7231,'cal/(mol*K)','+|-',3.56304),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         1
""",
)

entry(
    index = 2243,
    label = "Cs-CNOO",
    group = 
"""
1 * Cs         u0 {2,S} {3,S} {4,S} {5,S}
2   C          u0 {1,S}
3   [O2s,O0sc] u0 {1,S}
4   [O2s,O0sc] u0 {1,S}
5   N          u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.33474,8.87704,9.54314,9.7458,9.77737,9.36489,8.20634],'cal/(mol*K)','+|-',[1.01392,1.08285,1.11676,1.10883,1.07863,1.05928,1.08086]),
        H298 = (-18.9312,'kcal/mol','+|-',4.06395),
        S298 = (-37.0406,'cal/(mol*K)','+|-',3.11364),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         1
""",
)

entry(
    index = 2244,
    label = "Cs-CCHN",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   C  u0 {1,S}
4   N  u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = None,
    shortDesc = """Derived from nitrogen species in RMG thermo libraries""",
    longDesc = 
"""

""",
)

entry(
    index = 2245,
    label = "Cs-N3dCsCsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   N3d u0 {1,S}
3   Cs  u0 {1,S}
4   Cs  u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.38415,7.44907,7.94256,8.21287,8.80238,9.12379,9.71106],'cal/(mol*K)','+|-',[0.762029,0.813835,0.839321,0.83336,0.81066,0.796114,0.812334]),
        H298 = (-3.34035,'kcal/mol','+|-',3.05432),
        S298 = (-13.2562,'cal/(mol*K)','+|-',2.34011),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         2
""",
)

entry(
    index = 2246,
    label = "Cs-(N3dN3d)CsCsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   N3d u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   Cs  u0 {1,S}
5   H   u0 {1,S}
6   N3d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.59191,10.6166,11.2121,11.7476,12.8403,13.6955,15.1241],'cal/(mol*K)','+|-',[0.973259,1.03943,1.07198,1.06436,1.03537,1.01679,1.03751]),
        H298 = (24.3702,'kcal/mol','+|-',3.90096),
        S298 = (-4.24004,'cal/(mol*K)','+|-',2.98877),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         1
""",
)

entry(
    index = 2247,
    label = "Cs-NCsCsH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   N  u0 {1,S}
3   Cs u0 {1,S}
4   Cs u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.41794,7.88289,8.14773,8.22253,8.52655,8.69213,9.11098],'cal/(mol*K)','+|-',[1.017,1.08614,1.12016,1.1122,1.08191,1.06249,1.08414]),
        H298 = (-2.96805,'kcal/mol','+|-',4.07629),
        S298 = (-10.3844,'cal/(mol*K)','+|-',3.1231),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         1
""",
)

entry(
    index = 2248,
    label = "Cs-N5dcCsCsH",
    group = 
"""
1 * Cs   u0 {2,S} {3,S} {4,S} {5,S}
2   N5dc u0 {1,S}
3   Cs   u0 {1,S}
4   Cs   u0 {1,S}
5   H    u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.47689,6.71666,7.71095,8.31512,9.0022,9.21721,9.06393],'cal/(mol*K)','+|-',[0.46892,0.500799,0.516482,0.512814,0.498845,0.489894,0.499875]),
        H298 = (-6.42174,'kcal/mol','+|-',1.8795),
        S298 = (-14.4091,'cal/(mol*K)','+|-',1.44),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library| Number of Species
CHON |         11
""",
)

entry(
    index = 2249,
    label = "Cs-N3sCsCsH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   N3s u0 {1,S}
3   Cs  u0 {1,S}
4   Cs  u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.10597,7.1538,8.11993,8.72777,9.33531,9.58782,10.0013],'cal/(mol*K)','+|-',[0.361372,0.385939,0.398025,0.395198,0.384434,0.377536,0.385227]),
        H298 = (1.16673,'kcal/mol','+|-',1.44843),
        S298 = (-12.3648,'cal/(mol*K)','+|-',1.10973),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         9
CHN     |         21
CHON    |         1
""",
)

entry(
    index = 2250,
    label = "Cs-NCsCdtH",
    group = 
"""
1 * Cs         u0 {2,S} {3,S} {4,S} {5,S}
2   N          u0 {1,S}
3   Cs         u0 {1,S}
4   [Cd,Ct,CO] u0 {1,S}
5   H          u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.24456,8.58036,9.00649,9.13736,9.37785,9.48483,9.41277],'cal/(mol*K)','+|-',[0.644515,0.688332,0.709887,0.704846,0.685647,0.673344,0.687062]),
        H298 = (-0.00838218,'kcal/mol','+|-',2.58331),
        S298 = (-12.7763,'cal/(mol*K)','+|-',1.97923),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         3
""",
)

entry(
    index = 2251,
    label = "Cs-CsN3sH(Cds-O2d)",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   CO  u0 {1,S} {6,D}
3   N3s u0 {1,S}
4   Cs  u0 {1,S}
5   H   u0 {1,S}
6   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.0659,7.07261,7.70277,7.95645,7.97087,7.8627,7.82994],'cal/(mol*K)','+|-',[0.591967,0.632211,0.652009,0.647379,0.629745,0.618445,0.631045]),
        H298 = (-1.55176,'kcal/mol','+|-',2.37269),
        S298 = (-14.8221,'cal/(mol*K)','+|-',1.81786),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHON_G4  |         1
BurcatNS |         1
CHON     |         2
""",
)

entry(
    index = 2252,
    label = "Cs-CsN3sH(Cds-N3d)",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   Cd  u0 {1,S} {6,D}
3   N3s u0 {1,S}
4   Cs  u0 {1,S}
5   H   u0 {1,S}
6   N3d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.23069,7.98096,8.25792,8.28513,8.43756,8.66653,9.46624],'cal/(mol*K)','+|-',[1.02433,1.09397,1.12823,1.12022,1.0897,1.07015,1.09195]),
        H298 = (0.175692,'kcal/mol','+|-',4.10567),
        S298 = (-10.9785,'cal/(mol*K)','+|-',3.14561),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         1
""",
)

entry(
    index = 2253,
    label = "Cs-CCNN",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   C  u0 {1,S}
4   N  u0 {1,S}
5   N  u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2254,
    label = "Cs-NNCsCs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   N  u0 {1,S}
3   N  u0 {1,S}
4   Cs u0 {1,S}
5   Cs u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.60462,8.44952,8.2066,7.60439,6.74705,5.95177,5.16843],'cal/(mol*K)','+|-',[1.16026,1.23914,1.27795,1.26887,1.23431,1.21216,1.23686]),
        H298 = (-0.637045,'kcal/mol','+|-',4.6505),
        S298 = (-30.36,'cal/(mol*K)','+|-',3.56304),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         1
""",
)

entry(
    index = 2255,
    label = "Cs-N5dcN5dcCsCs",
    group = 
"""
1 * Cs   u0 {2,S} {3,S} {4,S} {5,S}
2   N5dc u0 {1,S}
3   N5dc u0 {1,S}
4   Cs   u0 {1,S}
5   Cs   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.4106,8.10942,8.05312,7.56577,6.71762,6.01715,5.11996],'cal/(mol*K)','+|-',[1.20982,1.29206,1.33253,1.32306,1.28702,1.26393,1.28968]),
        H298 = (6.47103,'kcal/mol','+|-',4.84911),
        S298 = (-35.9117,'cal/(mol*K)','+|-',3.71521),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
BurcatNS |         1
""",
)

entry(
    index = 2256,
    label = "Cs-CCNO",
    group = 
"""
1 * Cs         u0 {2,S} {3,S} {4,S} {5,S}
2   C          u0 {1,S}
3   C          u0 {1,S}
4   [O2s,O0sc] u0 {1,S}
5   N          u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.67747,7.67187,7.96768,7.92763,7.80769,7.41834,6.6028],'cal/(mol*K)','+|-',[1.01392,1.08285,1.11676,1.10883,1.07863,1.05928,1.08086]),
        H298 = (-7.93776,'kcal/mol','+|-',4.06395),
        S298 = (-32.7147,'cal/(mol*K)','+|-',3.11364),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         1
""",
)

entry(
    index = 2257,
    label = "Cs-CCCN",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   C  u0 {1,S}
3   C  u0 {1,S}
4   C  u0 {1,S}
5   N  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.17023,6.87008,7.07447,7.05012,7.1538,7.11075,6.80941],'cal/(mol*K)','+|-',[0.639504,0.68298,0.704368,0.699366,0.680316,0.668109,0.68172]),
        H298 = (-3.30313,'kcal/mol','+|-',2.56322),
        S298 = (-31.1732,'cal/(mol*K)','+|-',1.96384),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library      | Number of Species
CBS_QB3_1dHR |         1
""",
)

entry(
    index = 2258,
    label = "Cs-N5dcCsCsCs",
    group = 
"""
1 * Cs   u0 {2,S} {3,S} {4,S} {5,S}
2   N5dc u0 {1,S}
3   Cs   u0 {1,S}
4   Cs   u0 {1,S}
5   Cs   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.86276,6.41576,7.47448,7.85991,7.71684,7.18634,6.03706],'cal/(mol*K)','+|-',[0.56617,0.604661,0.623596,0.619168,0.602302,0.591495,0.603545]),
        H298 = (-3.26773,'kcal/mol','+|-',2.26929),
        S298 = (-32.7727,'cal/(mol*K)','+|-',1.73864),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library| Number of Species
CHON |         5
""",
)

entry(
    index = 2259,
    label = "Cs-N3dCsCsCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   N3d u0 {1,S}
3   Cs  u0 {1,S}
4   Cs  u0 {1,S}
5   Cs  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([18.239,23.5129,27.0374,28.7404,30.6922,34.3289,37.5281],'cal/(mol*K)','+|-',[1.04482,1.11586,1.1508,1.14263,1.1115,1.09156,1.1138]),
        H298 = (13.1507,'kcal/mol','+|-',4.1878),
        S298 = (-30.0664,'cal/(mol*K)','+|-',3.20854),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library| Number of Species
CHON |         1
""",
)

entry(
    index = 2260,
    label = "Cs-(N3dN3d)CsCsCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   N3d u0 {1,S} {6,D}
3   Cs  u0 {1,S}
4   Cs  u0 {1,S}
5   Cs  u0 {1,S}
6   N3d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'cal/(mol*K)'),
        H298 = (-1.9,'kcal/mol'),
        S298 = (-34.7,'cal/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2261,
    label = "Cs-NCsCsCs",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   N  u0 {1,S}
3   Cs u0 {1,S}
4   Cs u0 {1,S}
5   Cs u0 {1,S}
""",
    thermo = 'Cs-N3sCsCsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2262,
    label = "Cs-N3sCsCsCs",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   N3s u0 {1,S}
3   Cs  u0 {1,S}
4   Cs  u0 {1,S}
5   Cs  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.61725,5.86697,6.92168,7.33987,7.06736,6.54724,5.79192],'cal/(mol*K)','+|-',[0.41391,0.442049,0.455892,0.452655,0.440325,0.432424,0.441234]),
        H298 = (5.37237,'kcal/mol','+|-',1.65901),
        S298 = (-37.0622,'cal/(mol*K)','+|-',1.27107),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         2
CHN     |         9
CHON    |         3
""",
)

entry(
    index = 2263,
    label = "Cs-NCCtCt",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   N  u0 {1,S}
3   C  u0 {1,S}
4   Ct u0 {1,S}
5   Ct u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2264,
    label = "Cs-NHHH",
    group = 
"""
1 * Cs u0 {2,S} {3,S} {4,S} {5,S}
2   N  u0 {1,S}
3   H  u0 {1,S}
4   H  u0 {1,S}
5   H  u0 {1,S}
""",
    thermo = 'Cs-N3sHHH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2265,
    label = "Cs-N3dHHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   N3d u0 {1,S}
3   H   u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.09633,8.39151,9.56074,10.5507,12.3136,13.6946,16.2755],'cal/(mol*K)','+|-',[0.676658,0.72266,0.745291,0.739998,0.719841,0.706925,0.721327]),
        H298 = (-9.07535,'kcal/mol','+|-',2.71214),
        S298 = (30.8015,'cal/(mol*K)','+|-',2.07794),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library                 | Number of Species
CHON_G4                 |         2
thermo_DFT_CCSDTF12_BAC |         1
""",
)

entry(
    index = 2266,
    label = "Cs-(N3dCd)HHH",
    group = 
"""
1 * Cs       u0 {2,S} {3,S} {4,S} {5,S}
2   N3d      u0 {1,S} {6,D}
3   H        u0 {1,S}
4   H        u0 {1,S}
5   H        u0 {1,S}
6   [Cd,Cdd] u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.17214,8.5373,9.9686,11.176,13.2783,14.8793,17.6244],'cal/(mol*K)','+|-',[0.434588,0.464133,0.478668,0.475269,0.462323,0.454027,0.463277]),
        H298 = (-4.75119,'kcal/mol','+|-',1.74189),
        S298 = (31.251,'cal/(mol*K)','+|-',1.33457),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         33
""",
)

entry(
    index = 2267,
    label = "Cs-(N3dN3d)HHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   N3d u0 {1,S} {6,D}
3   H   u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
6   N3d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.23039,10.961,12.6783,14.2805,16.9915,19.0887,22.379],'cal/(mol*K)','+|-',[0.22935,0.244942,0.252613,0.250819,0.243986,0.239609,0.24449]),
        H298 = (17.9013,'kcal/mol','+|-',0.919266),
        S298 = (37.2817,'cal/(mol*K)','+|-',0.704307),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         17
""",
)

entry(
    index = 2268,
    label = "Cs-N3sHHH",
    group = 
"""
1 * Cs  u0 {2,S} {3,S} {4,S} {5,S}
2   N3s u0 {1,S}
3   H   u0 {1,S}
4   H   u0 {1,S}
5   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.62523,8.77722,9.86496,10.8099,12.5403,13.9615,16.6573],'cal/(mol*K)','+|-',[0.333741,0.35643,0.367592,0.364982,0.35504,0.348669,0.355773]),
        H298 = (-5.6792,'kcal/mol','+|-',1.33768),
        S298 = (31.5835,'cal/(mol*K)','+|-',1.02488),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library                 | Number of Species
CHON_G4                 |         159
BurcatNS                |         1
thermo_DFT_CCSDTF12_BAC |         1
CHN                     |         21
CHON                    |         3
""",
)

entry(
    index = 2269,
    label = "Cs-N5sdtcHHH",
    group = 
"""
1 * Cs               u0 {2,S} {3,S} {4,S} {5,S}
2   [N5sc,N5dc,N5tc] u0 {1,S}
3   H                u0 {1,S}
4   H                u0 {1,S}
5   H                u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.90205,8.3084,9.4894,10.4687,12.2375,13.6625,16.1572],'cal/(mol*K)','+|-',[0.613584,0.655298,0.675819,0.67102,0.652742,0.641029,0.654089]),
        H298 = (-10.27,'kcal/mol','+|-',2.45933),
        S298 = (27.8734,'cal/(mol*K)','+|-',1.88425),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library                 | Number of Species
CHON_G4                 |         7
thermo_DFT_CCSDTF12_BAC |         2
""",
)

entry(
    index = 2270,
    label = "O",
    group = 
"""
1 * O u0
""",
    thermo = 'O2s-CsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2271,
    label = "O0sc",
    group = 
"""
1 * O0sc u0
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'cal/(mol*K)'),
        H298 = (0,'kcal/mol'),
        S298 = (0,'cal/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""
This should be 0 since we account for this atom with
neighbor at the center
""",
)

entry(
    index = 2272,
    label = "Oa(S)",
    group = 
"""
1 * O u0 p3 c0
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.968,4.968,4.968,4.968,4.968,4.968,4.968],'cal/(mol*K)'),
        H298 = (104.81,'kcal/mol'),
        S298 = (34.25,'cal/(mol*K)'),
    ),
    shortDesc = """PrimaryTHermoLibrary""",
    longDesc = 
"""
H298: ATcT version 1.110
level of theory energy: CCSD(T)F12A/cc-pVQZ-F12//CCSD(T)/cc-pVQZ
level of theory frequency: B3LYP/6-311++g(d,p)//B3LYP/6-311++g(d,p)
""",
)

entry(
    index = 2273,
    label = "O2d",
    group = 
"""
1 * O2d u0
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2274,
    label = "O2d-Cd",
    group = 
"""
1 * O2d u0 {2,D}
2   CO  u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """In this case the C is treated as the central atom""",
    longDesc = 
"""

""",
)

entry(
    index = 2275,
    label = "O2d-O2d",
    group = 
"""
1 * O2d u0 {2,D}
2   O2d u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.5,3.575,3.685,3.8,3.99,4.12,4.29],'cal/(mol*K)'),
        H298 = (14.01,'kcal/mol'),
        S298 = (24.085,'cal/(mol*K)'),
    ),
    shortDesc = """A. Vandeputte""",
    longDesc = 
"""

""",
)

entry(
    index = 2276,
    label = "O2d-Sd",
    group = 
"""
1 * O2d u0 {2,D}
2   S   ux {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
Inferred from a least squares fit from 40 species mostly calculated at cbsqb3, 4/2017, Ryan Gillis
""",
)

entry(
    index = 2277,
    label = "O2d-N3d",
    group = 
"""
1 * O2d u0 {2,D}
2   N3d u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'cal/(mol*K)'),
        H298 = (0,'kcal/mol'),
        S298 = (0,'cal/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2278,
    label = "O2d-N5dc",
    group = 
"""
1 * O2d  u0 {2,D}
2   N5dc u0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'cal/(mol*K)'),
        H298 = (0,'kcal/mol'),
        S298 = (0,'cal/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2279,
    label = "O2s",
    group = 
"""
1 * O2s u0
""",
    thermo = 'O2s-(Cds-Cd)(Cds-Cd)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2280,
    label = "O2sBrBr",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   Br  u0 {1,S}
3   Br  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([12.0459,12.7195,13.0348,13.2663,13.5717,13.7245,13.7921],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (24.5169,'kcal/mol','+|-',2.95256),
        S298 = (70.8159,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHOBr_G4 |         1
""",
)

entry(
    index = 2281,
    label = "O2sBrCl",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   Cl  u0 {1,S}
3   Br  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([11.7044,12.5173,12.9194,13.1709,13.5059,13.6779,13.769],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (21.5156,'kcal/mol','+|-',2.95256),
        S298 = (68.1197,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library    | Number of Species
CHOClBr_G4 |         1
""",
)

entry(
    index = 2282,
    label = "O2sClCl",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   Cl  u0 {1,S}
3   Cl  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([11.3756,12.2778,12.7915,13.0774,13.4492,13.6409,13.7462],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (18.5768,'kcal/mol','+|-',2.95256),
        S298 = (65.4078,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHOCl_G4 |         1
""",
)

entry(
    index = 2283,
    label = "O2sBrF",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   F   u0 {1,S}
3   Br  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([11.0349,11.9306,12.4937,12.8325,13.2825,13.5298,13.7028],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (16.8651,'kcal/mol','+|-',2.95256),
        S298 = (65.6233,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOFBr_G4 |         1
""",
)

entry(
    index = 2284,
    label = "O2sClF",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   F   u0 {1,S}
3   Cl  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([10.7218,11.6511,12.3053,12.7133,13.2176,13.4926,13.6787],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (13.066,'kcal/mol','+|-',2.95256),
        S298 = (62.8723,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOFCl_G4 |         1
""",
)

entry(
    index = 2285,
    label = "O2sFF",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   F   u0 {1,S}
3   F   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHOF_G4 |         1
""",
)

entry(
    index = 2286,
    label = "O2sBrH",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   H   u0 {1,S}
3   Br  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.13406,9.76541,10.2406,10.6068,11.2023,11.6552,12.389],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-14.711,'kcal/mol','+|-',2.95256),
        S298 = (59.1556,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHOBr_G4 |         1
""",
)

entry(
    index = 2287,
    label = "O2sClH",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   H   u0 {1,S}
3   Cl  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.89789,9.51248,10.0281,10.431,11.072,11.5614,12.3509],'cal/(mol*K)','+|-',[0.828227,0.951879,0.973053,0.958777,0.831389,0.717066,0.584393]),
        H298 = (-17.7584,'kcal/mol','+|-',2.95256),
        S298 = (56.499,'cal/(mol*K)','+|-',1.93805),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library  | Number of Species
CHOCl_G4 |         1
""",
)

entry(
    index = 2288,
    label = "O2sFH",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   H   u0 {1,S}
3   F   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHOF_G4 |         1
""",
)

entry(
    index = 2289,
    label = "O2sBrO",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   O   u0 {1,S}
3   Br  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([9.31398,9.6066,9.69709,9.68989,9.52101,9.2884,9.23091],'cal/(mol*K)','+|-',[0.0892288,0.10255,0.104832,0.103294,0.0895695,0.0772529,0.0629594]),
        H298 = (19.0958,'kcal/mol','+|-',0.318093),
        S298 = (38.41,'cal/(mol*K)','+|-',0.208795),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOBr_G4    |         23
CHOFClBr_G4 |         9
CHOFBr_G4   |         33
CHOClBr_G4  |         21
""",
)

entry(
    index = 2290,
    label = "O2sClO",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   O   u0 {1,S}
3   Cl  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.99675,9.33066,9.43256,9.46127,9.41254,9.23488,9.16751],'cal/(mol*K)','+|-',[0.117276,0.134785,0.137783,0.135762,0.117724,0.101536,0.0827492]),
        H298 = (15.7173,'kcal/mol','+|-',0.418079),
        S298 = (35.7869,'cal/(mol*K)','+|-',0.274425),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library    | Number of Species
CHOCl_G4   |         26
CHOFCl_G4  |         21
CHOClBr_G4 |         1
""",
)

entry(
    index = 2291,
    label = "O2sFO",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   O   u0 {1,S}
3   F   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library   | Number of Species
CHOF_G4   |         26
CHOFCl_G4 |         1
CHOFBr_G4 |         1
""",
)

entry(
    index = 2292,
    label = "O2sBrC",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   C   u0 {1,S}
3   Br  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.25839,8.58754,8.7708,8.94016,9.23695,9.39925,9.50708],'cal/(mol*K)','+|-',[0.0325225,0.037378,0.0382094,0.0376489,0.0326466,0.0281575,0.0229477]),
        H298 = (-0.215186,'kcal/mol','+|-',0.11594),
        S298 = (38.794,'cal/(mol*K)','+|-',0.0761025),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOBr_G4    |         159
CHOFClBr_G4 |         78
CHOFBr_G4   |         307
CHOClBr_G4  |         174
""",
)

entry(
    index = 2293,
    label = "O2sCCl",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   C   u0 {1,S}
3   Cl  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.95062,8.32428,8.52309,8.72182,9.07901,9.28258,9.43221],'cal/(mol*K)','+|-',[0.0413581,0.0475327,0.0485901,0.0478772,0.041516,0.0358072,0.029182]),
        H298 = (-2.55835,'kcal/mol','+|-',0.147438),
        S298 = (35.8183,'cal/(mol*K)','+|-',0.0967778),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOCl_G4    |         214
CHOFCl_G4   |         168
CHOFClBr_G4 |         3
CHOClBr_G4  |         33
""",
)

entry(
    index = 2294,
    label = "O2sCF",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   C   u0 {1,S}
3   F   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.40519,7.85409,8.12951,8.39159,8.89407,9.19261,9.46834],'cal/(mol*K)','+|-',[0.0468839,0.0538835,0.0550822,0.0542741,0.0470629,0.0405914,0.0330811]),
        H298 = (-6.87509,'kcal/mol','+|-',0.167137),
        S298 = (33.6636,'cal/(mol*K)','+|-',0.109708),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library     | Number of Species
CHOF_G4     |         215
CHOFCl_G4   |         34
CHOFClBr_G4 |         12
CHOFBr_G4   |         55
""",
)

entry(
    index = 2295,
    label = "O2s-HH",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   H   u0 {1,S}
3   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """O-HH WATER. !!!Using NIST value for H2O, S(group) = S(H2O) + Rln(2)""",
    longDesc = 
"""

""",
)

entry(
    index = 2296,
    label = "O2s-OsH",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   O2s u0 {1,S}
3   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """O-OH SANDIA 1/2*H2O2""",
    longDesc = 
"""

""",
)

entry(
    index = 2297,
    label = "O2s-(Os-CdOd)H",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   O2s u0 {1,S} {4,S}
3   H   u0 {1,S}
4   CO  u0 {2,S} {5,D}
5   O2d u0 {4,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.53,6.31,5.48,5.01,4.14,4.11,5.27],'cal/(mol*K)'),
        H298 = (0,'kcal/mol'),
        S298 = (27.83,'cal/(mol*K)','+|-',0.07),
    ),
    shortDesc = """H298 set to 0 to avoid double counting with O2s-O2s(Cds-O2d)""",
    longDesc = 
"""
Cpdata fit to OCHOOH in Thermo library: DFT_QCI_thermo
S298 taken from O2s-OsH
""",
)

entry(
    index = 2298,
    label = "O2s-OsOs",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   O2s u0 {1,S}
3   O2s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([2.2,3.64,4.2,4.34,4.62,4.9,4.9],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (8.85,'kcal/mol','+|-',0.16),
        S298 = (9.4,'cal/(mol*K)','+|-',0.1),
    ),
    shortDesc = """O-OO LAY 1997=20 !!!WARNING! Cp1500 value taken as Cp1000""",
    longDesc = 
"""

""",
)

entry(
    index = 2299,
    label = "O2s-SsOs",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   S   u0 {1,S}
3   O2s u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.22,3.81,4.08,4.22,4.78,4.9,4.88],'cal/(mol*K)'),
        H298 = (-3.35,'kcal/mol'),
        S298 = (8.69,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""

""",
)

entry(
    index = 2300,
    label = "O2s-CH",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   C   u0 {1,S}
3   H   u0 {1,S}
""",
    thermo = 'O2s-CsH',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2301,
    label = "O2s-CtH",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   Ct  u0 {1,S}
3   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """O-CtH BENSON (Assigned O-CsH)""",
    longDesc = 
"""

""",
)

entry(
    index = 2302,
    label = "O2s-CdsH",
    group = 
"""
1 * O2s     u0 {2,S} {3,S}
2   [Cd,CO] u0 {1,S}
3   H       u0 {1,S}
""",
    thermo = 'O2s-(Cds-Cd)H',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2303,
    label = "O2s-(Cds-O2d)H",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   CO  u0 {1,S} {4,D}
3   H   u0 {1,S}
4   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 2304,
    label = "O2s-(Cds-Cd)H",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   Cd  u0 {1,S} {4,D}
3   H   u0 {1,S}
4   C   u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([24.6,30.3,32.52,33.15,33.29,33.55,34.97],'J/(mol*K)','+|-',[4.18,4.18,4.18,4.18,4.18,4.18,4.18]),
        H298 = (-188.1,'kJ/mol','+|-',3.56),
        S298 = (106.3,'J/(mol*K)','+|-',4.87),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 2305,
    label = "O2s-(Cds-Nd)H",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   Cd  u0 {1,S} {4,D}
3   H   u0 {1,S}
4   N   u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([5.36464,6.0985,6.4769,6.70326,7.12478,7.39134,7.89159],'cal/(mol*K)','+|-',[0.412633,0.440686,0.454486,0.451259,0.438967,0.43109,0.439873]),
        H298 = (-38.2558,'kcal/mol','+|-',1.65389),
        S298 = (26.3761,'cal/(mol*K)','+|-',1.26715),
    ),
    shortDesc = """Derived from RMG Thermo Libraries""",
    longDesc = 
"""
Fitted using sklearn Ridge regression with alpha = 1e-06
Library | Number of Species
CHON_G4 |         46
""",
)

entry(
    index = 2306,
    label = "O2s-CsH",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   Cs  u0 {1,S}
3   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 2307,
    label = "O2s-CbH",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   Cb  u0 {1,S}
3   H   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.3,4.5,4.82,5.23,6.02,6.61,7.44],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (-37.9,'kcal/mol','+|-',0.16),
        S298 = (29.1,'cal/(mol*K)','+|-',0.1),
    ),
    shortDesc = """O-CbH BENSON (Assigned O-CsH)""",
    longDesc = 
"""

""",
)

entry(
    index = 2308,
    label = "O2s-CSH",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   CS  u0 {1,S} {4,D}
3   H   u0 {1,S}
4   S2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([2.4,3.15,4.19,5.05,6.21,7.04,8.41],'cal/(mol*K)'),
        H298 = (-47.58,'kcal/mol'),
        S298 = (27.77,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 2309,
    label = "O2s-OsC",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   O2s u0 {1,S}
3   C   u0 {1,S}
""",
    thermo = 'O2s-OsCs',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2310,
    label = "O2s-OsCt",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   O2s u0 {1,S}
3   Ct  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.9,4.31,4.6,4.84,5.32,5.8,5.8],'cal/(mol*K)','+|-',[0.15,0.15,0.15,0.15,0.15,0.15,0.15]),
        H298 = (7,'kcal/mol','+|-',0.3),
        S298 = (10.8,'cal/(mol*K)','+|-',0.15),
    ),
    shortDesc = """O-OCb Hf JWB plot S,Cp assigned O/O/Cd !!!WARNING! Cp1500 value taken as Cp1000""",
    longDesc = 
"""

""",
)

entry(
    index = 2311,
    label = "O2s-OsCds",
    group = 
"""
1 * O2s     u0 {2,S} {3,S}
2   O2s     u0 {1,S}
3   [Cd,CO] u0 {1,S}
""",
    thermo = 'O2s-O2s(Cds-Cd)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2312,
    label = "O2s-O2s(Cds-O2d)",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   CO  u0 {1,S} {4,D}
3   O2s u0 {1,S}
4   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.53,5.02,5.79,6.08,6.54,6.49,6.49],'cal/(mol*K)','+|-',[0.15,0.15,0.15,0.15,0.15,0.15,0.15]),
        H298 = (-20.22,'kcal/mol','+|-',0.3),
        S298 = (9.11,'cal/(mol*K)','+|-',0.15),
    ),
    shortDesc = """O-OCO jwl cbsQ 99 cqcho=20 !!!WARNING! Cp1500 value taken as Cp1000""",
    longDesc = 
"""
David - increased group value by 3 kcal/mol for more accurate estimates of: 
SMILES            H298 ref            before             after
O=COO       -69.1 (DFT_QCI)           -73.8              -70.8
O=C(C)OO    -84.1 (Klipp_Glar)        -86.5              -83.5
""",
)

entry(
    index = 2313,
    label = "O2s-O2s(Cds-Cd)",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   Cd  u0 {1,S} {4,D}
3   O2s u0 {1,S}
4   C   u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.5,3.87,3.95,4.15,4.73,4.89,4.89],'cal/(mol*K)','+|-',[0.15,0.15,0.15,0.15,0.15,0.15,0.15]),
        H298 = (1.64,'kcal/mol','+|-',0.3),
        S298 = (10.12,'cal/(mol*K)','+|-',0.15),
    ),
    shortDesc = """O-OCd WESTMORELAND S,Cp LAY'9405 !!!WARNING! Cp1500 value taken as Cp1000""",
    longDesc = 
"""

""",
)

entry(
    index = 2314,
    label = "O2s-OsCs",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   O2s u0 {1,S}
3   Cs  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.9,4.31,4.6,4.84,5.32,5.8,5.8],'cal/(mol*K)','+|-',[0.15,0.15,0.15,0.15,0.15,0.15,0.15]),
        H298 = (-5.4,'kcal/mol','+|-',0.3),
        S298 = (8.54,'cal/(mol*K)','+|-',0.15),
    ),
    shortDesc = """O-OCs LAY 1997 !!!WARNING! Cp1500 value taken as Cp1000""",
    longDesc = 
"""

""",
)

entry(
    index = 2315,
    label = "O2s-OsCb",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   O2s u0 {1,S}
3   Cb  u0 {1,S}
""",
    thermo = 'O2s-O2s(Cds-Cd)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2316,
    label = "O2s-CC",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   C   u0 {1,S}
3   C   u0 {1,S}
""",
    thermo = 'O2s-(Cds-Cd)(Cds-Cd)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2317,
    label = "O2s-CtCt",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   Ct  u0 {1,S}
3   Ct  u0 {1,S}
""",
    thermo = 'O2s-(Cds-Cd)(Cds-Cd)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2318,
    label = "O2s-CtCds",
    group = 
"""
1 * O2s     u0 {2,S} {3,S}
2   Ct      u0 {1,S}
3   [Cd,CO] u0 {1,S}
""",
    thermo = 'O2s-(Cds-Cd)(Cds-Cd)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2319,
    label = "O2s-Ct(Cds-O2d)",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   CO  u0 {1,S} {4,D}
3   Ct  u0 {1,S}
4   O2d u0 {2,D}
""",
    thermo = 'O2s-(Cds-Cd)(Cds-Cd)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2320,
    label = "O2s-Ct(Cds-Cd)",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   Cd  u0 {1,S} {4,D}
3   Ct  u0 {1,S}
4   C   u0 {2,D}
""",
    thermo = 'O2s-(Cds-Cd)(Cds-Cd)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2321,
    label = "O2s-CtCs",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   Ct  u0 {1,S}
3   Cs  u0 {1,S}
""",
    thermo = 'O2s-Cs(Cds-Cd)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2322,
    label = "O2s-CtCb",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   Ct  u0 {1,S}
3   Cb  u0 {1,S}
""",
    thermo = 'O2s-(Cds-Cd)(Cds-Cd)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2323,
    label = "O2s-CdsCds",
    group = 
"""
1 * O2s     u0 {2,S} {3,S}
2   [Cd,CO] u0 {1,S}
3   [Cd,CO] u0 {1,S}
""",
    thermo = 'O2s-(Cds-Cd)(Cds-Cd)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2324,
    label = "O2s-(Cds-O2d)(Cds-O2d)",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   CO  u0 {1,S} {4,D}
3   CO  u0 {1,S} {5,D}
4   O2d u0 {2,D}
5   O2d u0 {3,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([18.4,11.55,6.97,3.72,-0.53,-2.57,-1.41],'J/(mol*K)','+|-',[6.51,6.51,6.51,6.51,6.51,6.51,6.51]),
        H298 = (-46.4,'kJ/mol','+|-',5.54),
        S298 = (80.8,'J/(mol*K)','+|-',7.59),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 2325,
    label = "O2s-(Cds-O2d)(Cds-Cd)",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   CO  u0 {1,S} {5,D}
3   Cd  u0 {1,S} {4,D}
4   C   u0 {3,D}
5   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([20.02,19.61,18.5,17.71,17.02,16.49,15.33],'J/(mol*K)','+|-',[8.17,8.17,8.17,8.17,8.17,8.17,8.17]),
        H298 = (-100.6,'kJ/mol','+|-',6.96),
        S298 = (38.43,'J/(mol*K)','+|-',9.53),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 2326,
    label = "O2s-(Cds-Cd)(Cds-Cd)",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   Cd  u0 {1,S}
3   Cd  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.4,3.7,3.7,3.8,4.4,4.6,4.8],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (-19.61,'kcal/mol','+|-',0.19),
        S298 = (10,'cal/(mol*K)','+|-',0.1),
    ),
    shortDesc = """O-CdCd BOZZELLI""",
    longDesc = 
"""

""",
)

entry(
    index = 2327,
    label = "O2s-CdsCs",
    group = 
"""
1 * O2s     u0 {2,S} {3,S}
2   [Cd,CO] u0 {1,S}
3   Cs      u0 {1,S}
""",
    thermo = 'O2s-Cs(Cds-Cd)',
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2328,
    label = "O2s-Cs(Cds-O2d)",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   CO  u0 {1,S} {4,D}
3   Cs  u0 {1,S}
4   O2d u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 2329,
    label = "O2s-Cs(Cds-Cd)",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   Cd  u0 {1,S} {4,D}
3   Cs  u0 {1,S}
4   C   u0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([19.07,23.32,25.26,25.92,25.5,24.52,22.72],'J/(mol*K)','+|-',[3.47,3.47,3.47,3.47,3.47,3.47,3.47]),
        H298 = (-123.9,'kJ/mol','+|-',2.96),
        S298 = (18.91,'J/(mol*K)','+|-',4.05),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 2330,
    label = "O2s-CdsCb",
    group = 
"""
1 * O2s     u0 {2,S} {3,S}
2   [Cd,CO] u0 {1,S}
3   Cb      u0 {1,S}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2331,
    label = "O2s-Cb(Cds-O2d)",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   CO  u0 {1,S} {4,D}
3   Cb  u0 {1,S}
4   O2d u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2332,
    label = "O2s-Cb(Cds-Cd)",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   Cd  u0 {1,S} {4,D}
3   Cb  u0 {1,S}
4   C   u0 {2,D}
""",
    thermo = None,
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2333,
    label = "O2s-CsCs",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   Cs  u0 {1,S}
3   Cs  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([14.7,13.4,13.58,14.54,16.71,18.29,20.17],'J/(mol*K)','+|-',[2.44,2.44,2.44,2.44,2.44,2.44,2.44]),
        H298 = (-98.6,'kJ/mol','+|-',2.08),
        S298 = (38.61,'J/(mol*K)','+|-',2.85),
    ),
    shortDesc = """\Derived from CBS-QB3 calculation with 1DHR treatment""",
    longDesc = 
"""
Derived using calculations at B3LYP/6-311G(d,p)/CBS-QB3 level of theory. 1DH-rotors
optimized at the B3LYP/6-31G(d).Paraskevas et al, Chem. Eur. J. 2013, 19, 16431-16452,
DOI: 10.1002/chem.201301381
""",
)

entry(
    index = 2334,
    label = "O2s-CsCb",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   Cs  u0 {1,S}
3   Cb  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.4,3.7,3.7,3.8,4.4,4.6,4.6],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (-22.6,'kcal/mol','+|-',0.19),
        S298 = (9.7,'cal/(mol*K)','+|-',0.1),
    ),
    shortDesc = """O-CbCs REID, PRAUSNITZ and SHERWOOD !!!WARNING! Cp1500 value taken as Cp1000""",
    longDesc = 
"""

""",
)

entry(
    index = 2335,
    label = "O2s-CbCb",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   Cb  u0 {1,S}
3   Cb  u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([1.19,-0.24,-0.72,-0.51,0.43,1.36,1.75],'cal/(mol*K)','+|-',[0.1,0.1,0.1,0.1,0.1,0.1,0.1]),
        H298 = (-18.77,'kcal/mol','+|-',0.19),
        S298 = (13.59,'cal/(mol*K)','+|-',0.1),
    ),
    shortDesc = """O-CbCb CHERN 1/97 Hf PEDLEY, Mopac""",
    longDesc = 
"""

""",
)

entry(
    index = 2336,
    label = "O2s-Cs(Cds-S2d)",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   Cs  u0 {1,S}
3   CS  u0 {1,S} {4,D}
4   S2d u0 {3,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([1.05,1.22,1.84,2.44,3.31,3.94,4.86],'cal/(mol*K)'),
        H298 = (-30.58,'kcal/mol'),
        S298 = (5.73,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""

""",
)

entry(
    index = 2337,
    label = "O2s-CS",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   S   ux {1,S}
3   C   ux {1,S}
""",
    thermo = 'O2s-CS4',
    shortDesc = """Sulfur/Oxygen Extension, Ryan Gillis""",
    longDesc = 
"""
"
""",
)

entry(
    index = 2338,
    label = "O2s-CS2",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   S2s ux {1,S}
3   Cs  ux {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([2.78,3.48,3.83,4.16,5.02,5.4,5.86],'cal/(mol*K)'),
        H298 = (-22.57,'kcal/mol'),
        S298 = (7.17,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 2339,
    label = "O2s-CS4",
    group = 
"""
1 * O2s               u0 {2,S} {3,S}
2   [S4s,S4d,S4b,S4t] ux {1,S}
3   C                 ux {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 2340,
    label = "O2s-CS6",
    group = 
"""
1 * O2s                     u0 {2,S} {3,S}
2   [S6s,S6d,S6dd,S6t,S6td] ux {1,S}
3   C                       ux {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 2341,
    label = "O2s-SH",
    group = 
"""
1 * O2s u0 {2,S} {3,S}
2   S   ux {1,S}
3   H   ux {1,S}
""",
    thermo = 'O2s-S_nonDeH',
    shortDesc = """Sulfur/Oxygen Extension, Ryan Gillis""",
    longDesc = 
"""
"
""",
)

entry(
    index = 2342,
    label = "O2s-S_nonDeH",
    group = 
"""
1 * O2s           u0 {2,S} {3,S}
2   [S2s,S4s,S6s] ux {1,S}
3   H             ux {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([3.88,4.91,5.45,5.86,6.66,7.13,7.81],'cal/(mol*K)'),
        H298 = (-36.34,'kcal/mol'),
        S298 = (29.09,'cal/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 2343,
    label = "O2s-S_DeH",
    group = 
"""
1 * O2s            u0 {2,S} {3,S}
2   [S4d,S6d,S6dd] ux {1,S}
3   H              ux {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'J/(mol*K)'),
        H298 = (0,'kJ/mol'),
        S298 = (0,'J/(mol*K)'),
    ),
    shortDesc = """RMG-type entries for Sulfur Groups, based on quantum calculations perfomred by Vandeputte (2011), Gillis, Class (2013), and Bozzelli, refit by Ryan Gillis in 2019""",
    longDesc = 
"""
"
""",
)

entry(
    index = 2344,
    label = "O2s-N",
    group = 
"""
1 * O2s u0 {2,S}
2   N   u0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0,0,0,0,0,0,0],'cal/(mol*K)'),
        H298 = (0,'kcal/mol'),
        S298 = (0,'cal/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2345,
    label = "O2s-N5tc",
    group = 
"""
1 * O2s  u0 {2,S}
2   N5tc u0 {1,S}
""",
    thermo = None,
    shortDesc = """Derived from nitrogen species in RMG thermo libraries""",
    longDesc = 
"""

""",
)

entry(
    index = 2346,
    label = "O2s-N5tcH",
    group = 
"""
1 * O2s  u0 {2,S} {3,S}
2   N5tc u0 {1,S}
3   H    u0 {1,S}
""",
    thermo = ThermoData(
    ),
"""