species(
    label = 'F[C](F)C(F)(F)C(F)(F)[C](F)F(619)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {12,S}
8  F u0 p3 c0 {12,S}
9  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
10 C u0 p0 c0 {3,S} {4,S} {9,S} {12,S}
11 C u1 p0 c0 {5,S} {6,S} {9,S}
12 C u1 p0 c0 {7,S} {8,S} {10,S}
"""),
    E0 = (-1338.27,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([181,249,286,344,486,552,511,665,553,637,1169,1241,1190,1306,146,234,414,562,504,606,1176,1296,1354,1460,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.360733,'amu*angstrom^2'), symmetry=1, barrier=(8.29397,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.360642,'amu*angstrom^2'), symmetry=1, barrier=(8.29187,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.360749,'amu*angstrom^2'), symmetry=1, barrier=(8.29434,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (200.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.48109,0.138783,-0.000270705,2.51251e-07,-8.7455e-11,-160777,35.2172], Tmin=(100,'K'), Tmax=(858.308,'K')), NASAPolynomial(coeffs=[13.079,0.0327531,-1.86904e-05,3.68687e-09,-2.53489e-13,-161870,-24.6107], Tmin=(858.308,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1338.27,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsF1sF1s)(F1s)(F1s)) + radical(Csj(Cs-CsF1sF1s)(F1s)(F1s))"""),
)

species(
    label = 'CF2CF2(60)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {6,S}
4 F u0 p3 c0 {6,S}
5 C u0 p0 c0 {1,S} {2,S} {6,D}
6 C u0 p0 c0 {3,S} {4,S} {5,D}
"""),
    E0 = (-688.535,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(99.9936,'amu')),
        NonlinearRotor(inertia=([91.6969,155.94,247.638],'amu*angstrom^2'), symmetry=4),
        HarmonicOscillator(frequencies=([198.559,203.927,398.069,428.925,524.86,551.48,555.021,803.576,1208.21,1359.53,1361.83,1918.98],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.015,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2113.54,'J/mol'), sigma=(4.647,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.86796,0.00896674,9.20781e-05,-2.50451e-07,1.90637e-10,-82806.8,9.05323], Tmin=(10,'K'), Tmax=(470.327,'K')), NASAPolynomial(coeffs=[5.38595,0.0191802,-1.42425e-05,4.78608e-09,-5.97116e-13,-83205.3,0.155967], Tmin=(470.327,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-688.535,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""FC(F)DC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CF2(42)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u2 p0 c0 {1,S} {2,S}
"""),
    E0 = (33.7272,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([700.388,993.445,1307.31],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (50.0074,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2405.34,'J/mol'), sigma=(4.22615,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=375.71 K, Pc=72.31 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28591,0.0107608,-1.05382e-05,4.89881e-09,-8.86384e-13,4216.69,13.1348], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[5.33121,0.00197748,-9.60248e-07,2.10704e-10,-1.5954e-14,3366.49,-2.56367], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(33.7272,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CF2(T)""", comment="""Thermo library: halogens"""),
)

species(
    label = 'F[C](F)C(F)(F)[C](F)F(519)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {8,S}
4 F u0 p3 c0 {8,S}
5 F u0 p3 c0 {9,S}
6 F u0 p3 c0 {9,S}
7 C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8 C u1 p0 c0 {3,S} {4,S} {7,S}
9 C u1 p0 c0 {5,S} {6,S} {7,S}
"""),
    E0 = (-917.879,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,146,234,414,562,504,606,1176,1296,1354,1460,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.339856,'amu*angstrom^2'), symmetry=1, barrier=(7.81396,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.647923,'amu*angstrom^2'), symmetry=1, barrier=(14.897,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (150.022,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.242198,0.0921598,-0.000169907,1.51452e-07,-5.14312e-11,-110269,27.3581], Tmin=(100,'K'), Tmax=(842.398,'K')), NASAPolynomial(coeffs=[12.0022,0.0196248,-1.10215e-05,2.18602e-09,-1.51603e-13,-111658,-23.8442], Tmin=(842.398,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-917.879,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsF1sF1s)(F1s)(F1s)) + radical(Csj(Cs-CsF1sF1s)(F1s)(F1s))"""),
)

species(
    label = 'FC1(F)C(F)(F)C(F)(F)C1(F)F(781)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {12,S}
8  F u0 p3 c0 {12,S}
9  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
10 C u0 p0 c0 {3,S} {4,S} {9,S} {12,S}
11 C u0 p0 c0 {5,S} {6,S} {9,S} {12,S}
12 C u0 p0 c0 {7,S} {8,S} {10,S} {11,S}
"""),
    E0 = (-1617.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (200.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.353055,0.0990718,-0.000133766,8.91252e-08,-2.32203e-11,-194395,26.6178], Tmin=(100,'K'), Tmax=(943.085,'K')), NASAPolynomial(coeffs=[17.8964,0.0216675,-1.06511e-05,2.09379e-09,-1.48983e-13,-197837,-60.3554], Tmin=(943.085,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1617.57,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + group(CsCsCsFF) + group(CsCsCsFF) + group(CsCsCsFF) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)2) + ring(Cs-Cs-Cs(F)(F)-Cs)"""),
)

species(
    label = 'F(37)',
    structure = adjacencyList("""multiplicity 2
1 F u1 p3 c0
"""),
    E0 = (72.8916,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (18.9984,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(665.158,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.90371,-0.000635296,2.64735e-07,7.69063e-11,-5.45254e-14,8672.27,2.70828], Tmin=(298,'K'), Tmax=(1400,'K')), NASAPolynomial(coeffs=[2.65117,-0.00014013,5.19236e-08,-8.84954e-12,5.9028e-16,8758.29,4.07857], Tmin=(1400,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(72.8916,'kJ/mol'), Cp0=(20.7862,'J/mol/K'), CpInf=(20.7862,'J/mol/K'), label="""F""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'F[C](F)C(F)(F)C(F)=C(F)F(759)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {11,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {8,S} {11,D}
10 C u1 p0 c0 {4,S} {5,S} {8,S}
11 C u0 p0 c0 {6,S} {7,S} {9,D}
"""),
    E0 = (-1164.9,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,323,467,575,827,1418,190,488,555,1236,1407,182,240,577,636,1210,1413,180],'cm^-1')),
        HinderedRotor(inertia=(0.263211,'amu*angstrom^2'), symmetry=1, barrier=(6.05174,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.263679,'amu*angstrom^2'), symmetry=1, barrier=(6.0625,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (181.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.1053,0.0841355,-0.000138223,1.21252e-07,-4.32074e-11,-140102,16.8264], Tmin=(10,'K'), Tmax=(670.261,'K')), NASAPolynomial(coeffs=[11.627,0.0332798,-2.44112e-05,8.05121e-09,-9.84913e-13,-141244,-20.8764], Tmin=(670.261,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1164.9,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), label="""F[C](F)C(F)(F)C(F)DC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[C](F)[C](F)F(230)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {6,S}
4 F u0 p3 c0 {6,S}
5 C u1 p0 c0 {1,S} {2,S} {6,S}
6 C u1 p0 c0 {3,S} {4,S} {5,S}
"""),
    E0 = (-493.512,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([146,234,414,562,504,606,1176,1296,1354,1460,180],'cm^-1')),
        HinderedRotor(inertia=(0.200092,'amu*angstrom^2'), symmetry=1, barrier=(4.6005,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.015,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.04125,0.0516239,-0.000104942,1.02887e-07,-3.74503e-11,-59293.5,18.0917], Tmin=(100,'K'), Tmax=(852.08,'K')), NASAPolynomial(coeffs=[5.53325,0.0162398,-9.21972e-06,1.83691e-09,-1.27509e-13,-59199.1,5.8493], Tmin=(852.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-493.512,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCsF1sF1s) + radical(CsCsF1sF1s)"""),
)

species(
    label = 'F2(106)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 F u0 p3 c0 {1,S}
"""),
    E0 = (-8.80492,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(37.9968,'amu')),
        LinearRotor(inertia=(18.2987,'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([1076.6],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (37.9968,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66955,0.00529418,-6.47091e-06,3.91833e-09,-9.38451e-13,-980.801,7.82659], Tmin=(298,'K'), Tmax=(1150,'K')), NASAPolynomial(coeffs=[4.05774,0.000600034,-2.19218e-07,4.31508e-11,-3.12588e-15,-1324.39,0.863214], Tmin=(1150,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-8.80492,'kJ/mol'), Cp0=(29.1007,'J/mol/K'), CpInf=(37.4151,'J/mol/K'), label="""F2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'F[C](F)C(F)=C(F)[C](F)F(760)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {1,S} {8,D} {9,S}
8  C u0 p0 c0 {2,S} {7,D} {10,S}
9  C u1 p0 c0 {3,S} {4,S} {7,S}
10 C u1 p0 c0 {5,S} {6,S} {8,S}
"""),
    E0 = (-827.516,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([206,336,431,607,515,611,528,696,1312,1446,124,198,231,363,414,566,474,694,695,865,1255,1461],'cm^-1')),
        HinderedRotor(inertia=(0.31341,'amu*angstrom^2'), symmetry=1, barrier=(7.20592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.307831,'amu*angstrom^2'), symmetry=1, barrier=(7.07765,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (162.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.322626,0.10945,-0.000211377,1.9481e-07,-6.70884e-11,-99385.2,30.3248], Tmin=(100,'K'), Tmax=(872.086,'K')), NASAPolynomial(coeffs=[10.9888,0.0261641,-1.41082e-05,2.71576e-09,-1.83654e-13,-100164,-15.8515], Tmin=(872.086,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-827.516,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cd-F1sCd)(F1s)(F1s)) + radical(Csj(Cd-F1sCd)(F1s)(F1s))"""),
)

species(
    label = 'CF2(87)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p1 c0 {1,S} {2,S}
"""),
    E0 = (-203.712,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([192,594,627],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (50.0074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28591,0.0107608,-1.05382e-05,4.89881e-09,-8.86384e-13,-24340.7,13.1348], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[5.33121,0.00197748,-9.60248e-07,2.10704e-10,-1.5954e-14,-25190.9,-2.56367], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-203.712,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CF2""", comment="""Thermo library: halogens"""),
)

species(
    label = 'F[C](F)C(F)(F)[C](F)C(F)(F)F(782)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {12,S}
8  F u0 p3 c0 {12,S}
9  C u0 p0 c0 {1,S} {2,S} {11,S} {12,S}
10 C u0 p0 c0 {3,S} {4,S} {5,S} {11,S}
11 C u1 p0 c0 {6,S} {9,S} {10,S}
12 C u1 p0 c0 {7,S} {8,S} {9,S}
"""),
    E0 = (-1370.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (200.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.08741,0.123416,-0.000219793,1.92987e-07,-6.52804e-11,-164618,35.7193], Tmin=(100,'K'), Tmax=(825.849,'K')), NASAPolynomial(coeffs=[15.0709,0.0280712,-1.55921e-05,3.10018e-09,-2.16368e-13,-166705,-35.6167], Tmin=(825.849,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1370.14,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sF1s)(F1s)) + radical(Csj(Cs-CsF1sF1s)(F1s)(F1s))"""),
)

species(
    label = 'F[C](F)[C](F)C(F)(F)C(F)(F)F(783)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {12,S}
8  F u0 p3 c0 {12,S}
9  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
10 C u0 p0 c0 {3,S} {4,S} {5,S} {9,S}
11 C u1 p0 c0 {6,S} {9,S} {12,S}
12 C u1 p0 c0 {7,S} {8,S} {11,S}
"""),
    E0 = (-1372.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (200.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.991738,0.131854,-0.000268896,2.63703e-07,-9.59654e-11,-164936,36.3531], Tmin=(100,'K'), Tmax=(854.241,'K')), NASAPolynomial(coeffs=[7.63369,0.041947,-2.40745e-05,4.78183e-09,-3.30955e-13,-164603,6.67431], Tmin=(854.241,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1372.67,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(F1s)) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
)

species(
    label = 'N2',
    structure = adjacencyList("""1 N u0 p1 c0 {2,T}
2 N u0 p1 c0 {1,T}
"""),
    E0 = (-8.64289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0137,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(810.913,'J/mol'), sigma=(3.621,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53101,-0.000123661,-5.02999e-07,2.43531e-09,-1.40881e-12,-1046.98,2.96747], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.95258,0.0013969,-4.92632e-07,7.8601e-11,-4.60755e-15,-923.949,5.87189], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-8.64289,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""N2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'Ne',
    structure = adjacencyList("""1 Ne u0 p4 c0
"""),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (20.1801,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1235.53,'J/mol'), sigma=(3.758e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,3.35532], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,3.35532], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ne""", comment="""Thermo library: primaryThermoLibrary"""),
)

transitionState(
    label = 'TS1',
    E0 = (-613.708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-159.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-605.424,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-321.056,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-454.498,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-262.463,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-103.35,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-397.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-447.893,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-403.485,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['F[C](F)C(F)(F)C(F)(F)[C](F)F(619)'],
    products = ['CF2CF2(60)', 'CF2CF2(60)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CF2(42)', 'F[C](F)C(F)(F)[C](F)F(519)'],
    products = ['F[C](F)C(F)(F)C(F)(F)[C](F)F(619)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(4.0899e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['F[C](F)C(F)(F)C(F)(F)[C](F)F(619)'],
    products = ['FC1(F)C(F)(F)C(F)(F)C1(F)F(781)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_noH;Cpri_rad_out_noH]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['F(37)', 'F[C](F)C(F)(F)C(F)=C(F)F(759)'],
    products = ['F[C](F)C(F)(F)C(F)(F)[C](F)F(619)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(46.388,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CF2CF2(60)', 'F[C](F)[C](F)F(230)'],
    products = ['F[C](F)C(F)(F)C(F)(F)[C](F)F(619)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.61741e-05,'m^3/(mol*s)'), n=3.0961, Ea=(2.98668,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction6',
    reactants = ['F[C](F)[C](F)F(230)', 'F[C](F)[C](F)F(230)'],
    products = ['F[C](F)C(F)(F)C(F)(F)[C](F)F(619)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(5.26262e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0
Ea raised from -2.8 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['F2(106)', 'F[C](F)C(F)=C(F)[C](F)F(760)'],
    products = ['F[C](F)C(F)(F)C(F)(F)[C](F)F(619)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(8.40976,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CF2(87)', 'F[C](F)C(F)(F)[C](F)F(519)'],
    products = ['F[C](F)C(F)(F)C(F)(F)[C](F)F(619)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1573.45,'m^3/(mol*s)'), n=1.25031, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F[C](F)C(F)(F)C(F)(F)[C](F)F(619)'],
    products = ['F[C](F)C(F)(F)[C](F)C(F)(F)F(782)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.03052e+12,'s^-1'), n=0.18834, Ea=(165.815,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F[C](F)C(F)(F)C(F)(F)[C](F)F(619)'],
    products = ['F[C](F)[C](F)C(F)(F)C(F)(F)F(783)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(0.0372321,'s^-1'), n=4.16824, Ea=(210.224,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 4.0"""),
)

network(
    label = 'PDepNetwork #304',
    isomers = [
        'F[C](F)C(F)(F)C(F)(F)[C](F)F(619)',
    ],
    reactants = [
        ('CF2CF2(60)', 'CF2CF2(60)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #304',
    Tmin = (300,'K'),
    Tmax = (2500,'K'),
    Tcount = 8,
    Tlist = ([302.558,324.028,372.925,464.512,632.697,950.724,1545.17,2335.46],'K'),
    Pmin = (0.01,'bar'),
    Pmax = (100,'bar'),
    Pcount = 5,
    Plist = ([0.0125282,0.0667467,1,14.982,79.8202],'bar'),
    maximumGrainSize = (0.5,'kcal/mol'),
    minimumGrainCount = 250,
    method = 'modified strong collision',
    interpolationModel = ('Chebyshev', 6, 4),
    activeKRotor = True,
    activeJRotor = True,
    rmgmode = True,
)

