species(
    label = 'OOC(F)(F)C(F)[C](F)F(1096)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {7,S} {9,S}
7  O u0 p2 c0 {6,S} {12,S}
8  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
9  C u0 p0 c0 {2,S} {3,S} {6,S} {8,S}
10 C u1 p0 c0 {4,S} {5,S} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-1050.73,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,259,529,569,1128,1321,1390,3140,223,363,546,575,694,1179,1410,190,488,555,1236,1407,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.381526,'amu*angstrom^2'), symmetry=1, barrier=(8.77202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.76383,'amu*angstrom^2'), symmetry=1, barrier=(40.5539,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.76189,'amu*angstrom^2'), symmetry=1, barrier=(40.5094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.628843,'amu*angstrom^2'), symmetry=1, barrier=(14.4583,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (165.039,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3545.48,'J/mol'), sigma=(5.86849,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=553.80 K, Pc=39.81 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.482266,0.107352,-0.000178926,1.529e-07,-5.1546e-11,-126221,31.1515], Tmin=(100,'K'), Tmax=(777.733,'K')), NASAPolynomial(coeffs=[13.6254,0.0287898,-1.58252e-05,3.16502e-09,-2.23437e-13,-128233,-32.196], Tmin=(777.733,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1050.73,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCsCsFH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
)

species(
    label = 'CF2(89)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p1 c0 {1,S} {2,S}
"""),
    E0 = (-196.899,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([192,594,627],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (50.0074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28591,0.0107608,-1.05382e-05,4.89881e-09,-8.86384e-13,-23521.2,13.1348], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[5.33121,0.00197748,-9.60248e-07,2.10704e-10,-1.5954e-14,-24371.4,-2.56367], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-196.899,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CF2""", comment="""Thermo library: Fluorine"""),
)

species(
    label = 'OOC(F)[C](F)F(436)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {7,S}
4 O u0 p2 c0 {5,S} {6,S}
5 O u0 p2 c0 {4,S} {9,S}
6 C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
7 C u1 p0 c0 {2,S} {3,S} {6,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-587.25,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,487,638,688,1119,1325,1387,3149,190,488,555,1236,1407,180],'cm^-1')),
        HinderedRotor(inertia=(0.665616,'amu*angstrom^2'), symmetry=1, barrier=(15.3038,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.6675,'amu*angstrom^2'), symmetry=1, barrier=(15.3471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.79135,'amu*angstrom^2'), symmetry=1, barrier=(41.1867,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.54011,0.0496725,-6.0877e-05,3.93291e-08,-1.05092e-11,-70629.7,12.2903], Tmin=(10,'K'), Tmax=(857.261,'K')), NASAPolynomial(coeffs=[8.85859,0.0248563,-1.74546e-05,5.56076e-09,-6.61428e-13,-71541.6,-12.5492], Tmin=(857.261,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-587.25,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), label="""OOC(F)[C](F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OOC([C](F)F)C(F)(F)F(499)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {7,S} {8,S}
7  O u0 p2 c0 {6,S} {12,S}
8  C u0 p0 c0 {6,S} {9,S} {10,S} {11,S}
9  C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
10 C u1 p0 c0 {4,S} {5,S} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-1063.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (165.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.313577,0.102841,-0.000162209,1.30706e-07,-4.17427e-11,-127716,30.2371], Tmin=(100,'K'), Tmax=(768.726,'K')), NASAPolynomial(coeffs=[14.1869,0.0273873,-1.49734e-05,3.01421e-09,-2.14681e-13,-129945,-35.9051], Tmin=(768.726,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1063.12,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-CsOs) + longDistanceInteraction_noncyclic(Cs(F)3-R-Cs(F)2) + group(CsCsFFH) + radical(Csj(Cs-O2sCsH)(F1s)(F1s))"""),
)

species(
    label = 'OOC(F)(F)C(F)(F)[CH]F(500)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {7,S} {9,S}
7  O u0 p2 c0 {6,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {6,S} {8,S}
10 C u1 p0 c0 {5,S} {8,S} {11,S}
11 H u0 p0 c0 {10,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-1051.34,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,215,315,519,588,595,1205,1248,223,363,546,575,694,1179,1410,334,575,1197,1424,3202,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.297682,'amu*angstrom^2'), symmetry=1, barrier=(6.8443,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.74376,'amu*angstrom^2'), symmetry=1, barrier=(40.0924,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.29766,'amu*angstrom^2'), symmetry=1, barrier=(6.84379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.74403,'amu*angstrom^2'), symmetry=1, barrier=(40.0986,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (165.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.65814,0.113645,-0.000200749,1.78076e-07,-6.10334e-11,-126290,30.9793], Tmin=(100,'K'), Tmax=(826.408,'K')), NASAPolynomial(coeffs=[12.92,0.0300356,-1.65221e-05,3.26637e-09,-2.27482e-13,-127923,-28.2425], Tmin=(826.408,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1051.34,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(CsCsFHH) + radical(Csj(Cs-CsF1sF1s)(F1s)(H))"""),
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
    label = 'OOC(F)(F)C(F)[C]F(501)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8  C u0 p0 c0 {3,S} {7,S} {9,S} {10,S}
9  C u2 p0 c0 {4,S} {8,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-665.244,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([250,417,511,1155,1315,1456,3119,90,150,250,247,323,377,431,433,1065,1274,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (146.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.175858,0.0912693,-0.000137822,1.05653e-07,-3.21643e-11,-79879.2,26.2343], Tmin=(100,'K'), Tmax=(804.596,'K')), NASAPolynomial(coeffs=[13.5695,0.024681,-1.36772e-05,2.7852e-09,-2.00411e-13,-82034.4,-35.4695], Tmin=(804.596,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-665.244,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(CsCCl_triplet)"""),
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
    label = 'OOC(F)(F)[CH]F(437)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {7,S}
4 O u0 p2 c0 {5,S} {6,S}
5 O u0 p2 c0 {4,S} {9,S}
6 C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
7 C u1 p0 c0 {3,S} {6,S} {8,S}
8 H u0 p0 c0 {7,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-608.958,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,253,525,597,667,842,1178,1324,334,575,1197,1424,3202,180],'cm^-1')),
        HinderedRotor(inertia=(0.556187,'amu*angstrom^2'), symmetry=1, barrier=(12.7878,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31631,'amu*angstrom^2'), symmetry=1, barrier=(30.2646,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.556382,'amu*angstrom^2'), symmetry=1, barrier=(12.7923,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.42439,0.0503,-4.92602e-05,4.41425e-09,1.4757e-11,-73236,12.7569], Tmin=(10,'K'), Tmax=(609.328,'K')), NASAPolynomial(coeffs=[9.56374,0.0235609,-1.6825e-05,5.45817e-09,-6.59694e-13,-74236,-15.8867], Tmin=(609.328,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-608.958,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), label="""OOC(F)(F)[CH]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OH(4)',
    structure = adjacencyList("""multiplicity 2
1 O u1 p2 c0 {2,S}
2 H u0 p0 c0 {1,S}
"""),
    E0 = (28.3945,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3668.68],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (17.0074,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51457,2.92733e-05,-5.3215e-07,1.01947e-09,-3.85939e-13,3414.25,2.10435], Tmin=(100,'K'), Tmax=(1145.76,'K')), NASAPolynomial(coeffs=[3.07194,0.00060402,-1.39806e-08,-2.13441e-11,2.48061e-15,3579.39,4.57801], Tmin=(1145.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.3945,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH(D)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'FC1C(F)(F)OC1(F)F(341)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  O u0 p2 c0 {8,S} {9,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
9  C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-1203.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (148.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.78842,0.0142565,0.000167223,-4.24302e-07,3.09636e-10,-144772,12.9607], Tmin=(10,'K'), Tmax=(476.767,'K')), NASAPolynomial(coeffs=[4.91851,0.0408334,-2.9839e-05,9.88299e-09,-1.217e-12,-145290,4.0479], Tmin=(476.767,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1203.73,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), label="""FC1C(F)(F)OC1(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OOC(F)(F)C=C(F)F(502)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8  C u0 p0 c0 {7,S} {9,D} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {8,D}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-935.437,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,275,321,533,585,746,850,1103,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,180],'cm^-1')),
        HinderedRotor(inertia=(0.250264,'amu*angstrom^2'), symmetry=1, barrier=(5.75407,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.991048,'amu*angstrom^2'), symmetry=1, barrier=(22.7861,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.30447,'amu*angstrom^2'), symmetry=1, barrier=(52.9843,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.496325,0.0839304,-0.000124176,9.63844e-08,-3.00496e-11,-112387,25.4601], Tmin=(100,'K'), Tmax=(783.188,'K')), NASAPolynomial(coeffs=[11.722,0.0265947,-1.43595e-05,2.90233e-09,-2.08158e-13,-114145,-25.9535], Tmin=(783.188,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-935.437,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCFFO) + group(Cds-CdsCsH) + group(CdCFF)"""),
)

species(
    label = 'H(5)',
    structure = adjacencyList("""multiplicity 2
1 H u1 p0 c0
"""),
    E0 = (211.805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (1.00797,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1205.6,'J/mol'), sigma=(2.05,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.05029e-14,-1.20495e-16,5.23432e-20,-6.99814e-24,25474.2,-0.444973], Tmin=(100,'K'), Tmax=(4230.24,'K')), NASAPolynomial(coeffs=[2.5,3.58763e-09,-1.22015e-12,1.84118e-16,-1.04e-20,25474.2,-0.444948], Tmin=(4230.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.805,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""H""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'OOC(F)(F)C(F)=C(F)F(504)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {7,S} {8,S}
7  O u0 p2 c0 {6,S} {11,S}
8  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
9  C u0 p0 c0 {3,S} {8,S} {10,D}
10 C u0 p0 c0 {4,S} {5,S} {9,D}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-1087.92,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,275,321,533,585,746,850,1103,323,467,575,827,1418,182,240,577,636,1210,1413,180],'cm^-1')),
        HinderedRotor(inertia=(0.325481,'amu*angstrom^2'), symmetry=1, barrier=(7.48344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.618986,'amu*angstrom^2'), symmetry=1, barrier=(14.2317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.31086,'amu*angstrom^2'), symmetry=1, barrier=(53.1313,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (164.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.350879,0.104868,-0.000181285,1.57358e-07,-5.31601e-11,-130699,28.2332], Tmin=(100,'K'), Tmax=(811.02,'K')), NASAPolynomial(coeffs=[13.3298,0.026011,-1.43834e-05,2.85788e-09,-2.00023e-13,-132544,-32.5945], Tmin=(811.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1087.92,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CdCsCdF) + group(CdCFF) + longDistanceInteraction_noncyclic(Cds(F)2=Cds(F))"""),
)

species(
    label = 'OOC(F)(F)[CH][C](F)F(505)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8  C u1 p0 c0 {7,S} {9,S} {10,S}
9  C u1 p0 c0 {3,S} {4,S} {8,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-678.647,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,253,525,597,667,842,1178,1324,3025,407.5,1350,352.5,190,488,555,1236,1407,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.222306,'amu*angstrom^2'), symmetry=1, barrier=(5.11125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.222489,'amu*angstrom^2'), symmetry=1, barrier=(5.11547,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.8486,'amu*angstrom^2'), symmetry=1, barrier=(42.503,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31555,'amu*angstrom^2'), symmetry=1, barrier=(30.2472,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (146.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.036242,0.0993823,-0.000177104,1.60414e-07,-5.62348e-11,-81487.3,31.0779], Tmin=(100,'K'), Tmax=(818.333,'K')), NASAPolynomial(coeffs=[10.7393,0.0291871,-1.63135e-05,3.25442e-09,-2.28138e-13,-82664.1,-15.162], Tmin=(818.333,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-678.647,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(CsCFFO) + group(CsCsFFH) + radical(CCJCOOH) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = '[O]C(F)(F)C(F)[C](F)F(338)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
9  C u1 p0 c0 {4,S} {5,S} {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-871.425,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,351,323,533,609,664,892,1120,1201,190,488,555,1236,1407,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.831114,'amu*angstrom^2'), symmetry=1, barrier=(19.1089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.519034,'amu*angstrom^2'), symmetry=1, barrier=(11.9336,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (148.031,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3397,'J/mol'), sigma=(5.67045,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=530.60 K, Pc=42.28 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.141903,0.0954747,-0.000174281,1.59093e-07,-5.58175e-11,-104680,26.5838], Tmin=(100,'K'), Tmax=(822.779,'K')), NASAPolynomial(coeffs=[10.65,0.026186,-1.47764e-05,2.96419e-09,-2.08139e-13,-105793,-18.3179], Tmin=(822.779,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-871.425,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sF1s)) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
)

species(
    label = 'HO2(10)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,S} {3,S}
2 O u1 p2 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (2.49012,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1064.4,1465.7,3224.93],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (33.0068,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.3018,-0.00474912,2.11583e-05,-2.42764e-08,9.29225e-12,264.018,3.71666], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.17229,0.00188118,-3.46277e-07,1.94658e-11,1.76257e-16,31.0207,2.95768], Tmin=(1000,'K'), Tmax=(5000,'K'))], Tmin=(200,'K'), Tmax=(5000,'K'), E0=(2.49012,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HO2""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'F[C](F)C(F)[C](F)F(272)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {7,S}
4 F u0 p3 c0 {8,S}
5 F u0 p3 c0 {8,S}
6 C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7 C u1 p0 c0 {2,S} {3,S} {6,S}
8 C u1 p0 c0 {4,S} {5,S} {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-711.656,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,146,234,414,562,504,606,1176,1296,1354,1460,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.130631,'amu*angstrom^2'), symmetry=1, barrier=(3.00347,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130993,'amu*angstrom^2'), symmetry=1, barrier=(3.01179,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03372,0.0738643,-0.000133097,1.21982e-07,-4.27798e-11,-85493.9,26.1783], Tmin=(100,'K'), Tmax=(838.788,'K')), NASAPolynomial(coeffs=[8.34888,0.0225509,-1.19533e-05,2.34634e-09,-1.62743e-13,-86143.1,-4.38197], Tmin=(838.788,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-711.656,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsF1sH)(F1s)(F1s)) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
)

species(
    label = '[O]OC(F)(F)C(F)[C](F)F(489)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {7,S} {9,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
9  C u0 p0 c0 {2,S} {3,S} {6,S} {8,S}
10 C u1 p0 c0 {4,S} {5,S} {8,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-898.725,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,259,529,569,1128,1321,1390,3140,223,363,546,575,694,1179,1410,190,488,555,1236,1407,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.384579,'amu*angstrom^2'), symmetry=1, barrier=(8.84222,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31388,'amu*angstrom^2'), symmetry=1, barrier=(30.2087,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.383968,'amu*angstrom^2'), symmetry=1, barrier=(8.82819,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (164.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.270775,0.104398,-0.000185892,1.65006e-07,-5.63342e-11,-107948,31.1495], Tmin=(100,'K'), Tmax=(835.656,'K')), NASAPolynomial(coeffs=[12.2173,0.0268082,-1.46426e-05,2.87794e-09,-1.99394e-13,-109413,-23.1345], Tmin=(835.656,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-898.725,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCsCsFH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + radical(ROOJ) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
)

species(
    label = 'OOC(F)(F)[C](F)[C](F)F(507)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {7,S} {8,S}
7  O u0 p2 c0 {6,S} {11,S}
8  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
9  C u1 p0 c0 {3,S} {8,S} {10,S}
10 C u1 p0 c0 {4,S} {5,S} {9,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-852.707,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,253,525,597,667,842,1178,1324,212,367,445,1450,190,488,555,1236,1407,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.234022,'amu*angstrom^2'), symmetry=1, barrier=(5.38062,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09069,'amu*angstrom^2'), symmetry=1, barrier=(25.077,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.43013,'amu*angstrom^2'), symmetry=1, barrier=(55.8734,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.234104,'amu*angstrom^2'), symmetry=1, barrier=(5.38251,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (164.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.361666,0.105707,-0.000186131,1.6302e-07,-5.51798e-11,-102409,33.6799], Tmin=(100,'K'), Tmax=(825.46,'K')), NASAPolynomial(coeffs=[13.1689,0.0256058,-1.41615e-05,2.80066e-09,-1.94964e-13,-104148,-26.0022], Tmin=(825.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-852.707,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCsCsFH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + radical(Csj(Cs-F1sF1sO2s)(Cs-F1sF1sH)(F1s)) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
)

species(
    label = 'OOC(F)(F)C(F)[C]F-2(510)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {7,S}
6  O u0 p2 c0 {5,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8  C u0 p0 c0 {3,S} {7,S} {9,S} {10,S}
9  C u0 p1 c0 {4,S} {8,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-670.714,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,223,363,546,575,694,1179,1410,353,444,1253,3145,617,898,1187,294.081,294.088,294.091,2612.01],'cm^-1')),
        HinderedRotor(inertia=(0.103661,'amu*angstrom^2'), symmetry=1, barrier=(6.36424,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.558799,'amu*angstrom^2'), symmetry=1, barrier=(34.3058,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.558744,'amu*angstrom^2'), symmetry=1, barrier=(34.3057,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.558971,'amu*angstrom^2'), symmetry=1, barrier=(34.3056,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.374964,0.0868527,-0.000129181,9.90581e-08,-3.03692e-11,-80544.2,27.2972], Tmin=(100,'K'), Tmax=(797.215,'K')), NASAPolynomial(coeffs=[12.5509,0.0257517,-1.41999e-05,2.89282e-09,-2.08403e-13,-82485.3,-28.6837], Tmin=(797.215,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-670.714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCCFH) + group(CJ2_singlet-FCs)"""),
)

species(
    label = 'OOC(F)(F)[C](F)C(F)F(495)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {7,S} {8,S}
7  O u0 p2 c0 {6,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {10,S} {11,S}
10 C u1 p0 c0 {5,S} {8,S} {9,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-1054.08,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,253,525,597,667,842,1178,1324,522,611,926,1093,1137,1374,1416,3112,212,367,445,1450,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.584564,'amu*angstrom^2'), symmetry=1, barrier=(13.4403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.583758,'amu*angstrom^2'), symmetry=1, barrier=(13.4217,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06855,'amu*angstrom^2'), symmetry=1, barrier=(24.5681,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.2591,'amu*angstrom^2'), symmetry=1, barrier=(51.9412,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (165.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.190018,0.0990841,-0.000148994,1.13505e-07,-3.418e-11,-126632,31.7767], Tmin=(100,'K'), Tmax=(814.905,'K')), NASAPolynomial(coeffs=[14.8016,0.0255016,-1.35591e-05,2.71399e-09,-1.93331e-13,-129076,-37.4818], Tmin=(814.905,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1054.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCsCsFH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + radical(Csj(Cs-F1sF1sO2s)(Cs-F1sF1sH)(F1s))"""),
)

species(
    label = '[O]OC(F)(F)C(F)C(F)F(1017)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {7,S} {9,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
9  C u0 p0 c0 {2,S} {3,S} {6,S} {8,S}
10 C u0 p0 c0 {4,S} {5,S} {8,S} {12,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-1100.1,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,250,417,511,1155,1315,1456,3119,223,363,546,575,694,1179,1410,235,523,627,1123,1142,1372,1406,3097,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.567912,'amu*angstrom^2'), symmetry=1, barrier=(13.0574,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.56688,'amu*angstrom^2'), symmetry=1, barrier=(13.0337,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49654,'amu*angstrom^2'), symmetry=1, barrier=(34.4085,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (165.039,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3545.48,'J/mol'), sigma=(5.86849,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=553.80 K, Pc=39.81 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.158567,0.0985478,-0.000151832,1.20119e-07,-3.76472e-11,-132168,29.4551], Tmin=(100,'K'), Tmax=(783.329,'K')), NASAPolynomial(coeffs=[13.9446,0.026534,-1.39383e-05,2.76652e-09,-1.95666e-13,-134378,-35.1413], Tmin=(783.329,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1100.1,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCsCsFH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + radical(ROOJ)"""),
)

species(
    label = '[O]C(F)(F)C(F)C(O)(F)F(511)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {9,S} {12,S}
7  O u1 p2 c0 {10,S}
8  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
9  C u0 p0 c0 {2,S} {3,S} {6,S} {8,S}
10 C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-1273.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (165.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.107071,0.0963347,-0.000140793,1.02921e-07,-2.96078e-11,-153029,29.7], Tmin=(100,'K'), Tmax=(853.515,'K')), NASAPolynomial(coeffs=[15.5147,0.0231244,-1.21328e-05,2.42864e-09,-1.73477e-13,-155696,-43.1922], Tmin=(853.515,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1273.54,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCsCsFH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + radical(O2sj(Cs-CsF1sF1s))"""),
)

species(
    label = 'OOC(F)(F)[CH]C(F)(F)F(512)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  O u0 p2 c0 {7,S} {8,S}
7  O u0 p2 c0 {6,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {5,S} {10,S}
10 C u1 p0 c0 {8,S} {9,S} {11,S}
11 H u0 p0 c0 {10,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-1126.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (165.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.239521,0.101418,-0.000161336,1.32605e-07,-4.33714e-11,-135299,31.0685], Tmin=(100,'K'), Tmax=(749.969,'K')), NASAPolynomial(coeffs=[13.4102,0.0286174,-1.57302e-05,3.17293e-09,-2.26168e-13,-137346,-30.8562], Tmin=(749.969,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1126.14,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(CsCFFO) + group(CsCsFFF) + radical(CCJCOOH)"""),
)

species(
    label = 'O=C(F)C(F)C(F)(F)F(513)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {9,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
9  C u0 p0 c0 {5,S} {6,D} {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-1271.85,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (148.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.45138,0.0554271,-6.36074e-05,3.79523e-08,-9.28041e-12,-152967,13.2558], Tmin=(10,'K'), Tmax=(938.495,'K')), NASAPolynomial(coeffs=[10.1572,0.0268458,-1.79258e-05,5.50192e-09,-6.36127e-13,-154226,-18.6702], Tmin=(938.495,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1271.85,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""ODC(F)C(F)C(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    E0 = (-51.528,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-218.834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-361.125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-59.4571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-42.3348,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-441.978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-266.979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-341.605,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-72.8598,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-301.503,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-176.27,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-154.025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-108.007,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-64.927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-272.961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-384.335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-489.326,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-407.795,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-371.294,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-373.136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction27',
    reactants = ['CF2(89)', 'OOC(F)[C](F)F(436)'],
    products = ['OOC(F)(F)C(F)[C](F)F(1096)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.40908e-07,'m^3/(mol*s)'), n=3.45227, Ea=(199.725,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s',), comment="""Estimated from node CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction28',
    reactants = ['OOC(F)(F)C(F)[C](F)F(1096)'],
    products = ['OOC([C](F)F)C(F)(F)F(499)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction29',
    reactants = ['OOC(F)(F)C(F)(F)[CH]F(500)'],
    products = ['OOC(F)(F)C(F)[C](F)F(1096)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction30',
    reactants = ['F(37)', 'OOC(F)(F)C(F)[C]F(501)'],
    products = ['OOC(F)(F)C(F)[C](F)F(1096)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [H/Val7_rad;Birad] for rate rule [Val7_rad;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['CF2(42)', 'OOC(F)(F)[CH]F(437)'],
    products = ['OOC(F)(F)C(F)[C](F)F(1096)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_sec_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction32',
    reactants = ['OOC(F)(F)C(F)[C](F)F(1096)'],
    products = ['OH(4)', 'FC1C(F)(F)OC1(F)F(341)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.31e+11,'s^-1','*|/',1.74), n=0, Ea=(75.8559,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R3OO_SS;C_ter_rad_intra;OOH]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['F(37)', 'OOC(F)(F)C=C(F)F(502)'],
    products = ['OOC(F)(F)C(F)[C](F)F(1096)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(62.6707,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction34',
    reactants = ['H(5)', 'OOC(F)(F)C(F)=C(F)F(504)'],
    products = ['OOC(F)(F)C(F)[C](F)F(1096)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(56.8734,'m^3/(mol*s)'), n=1.75834, Ea=(1.61913,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.2648472359903826, var=0.02782886759889551, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_4R!H->C_Sp-4C-1COS_Ext-1COS-R_N-6R!H-inRing_N-7R!H-inRing_Ext-4C-R_N-Sp-8R!H#4C',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_4R!H->C_Sp-4C-1COS_Ext-1COS-R_N-6R!H-inRing_N-7R!H-inRing_Ext-4C-R_N-Sp-8R!H#4C"""),
)

reaction(
    label = 'reaction35',
    reactants = ['F(37)', 'OOC(F)(F)[CH][C](F)F(505)'],
    products = ['OOC(F)(F)C(F)[C](F)F(1096)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -28.9 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction36',
    reactants = ['OH(4)', '[O]C(F)(F)C(F)[C](F)F(338)'],
    products = ['OOC(F)(F)C(F)[C](F)F(1096)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(8.63278,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction37',
    reactants = ['HO2(10)', 'F[C](F)C(F)[C](F)F(272)'],
    products = ['OOC(F)(F)C(F)[C](F)F(1096)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.6844e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction38',
    reactants = ['H(5)', '[O]OC(F)(F)C(F)[C](F)F(489)'],
    products = ['OOC(F)(F)C(F)[C](F)F(1096)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(43950,'m^3/(mol*s)'), n=1, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction39',
    reactants = ['H(5)', 'OOC(F)(F)[C](F)[C](F)F(507)'],
    products = ['OOC(F)(F)C(F)[C](F)F(1096)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C
Ea raised from -0.4 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction40',
    reactants = ['F(37)', 'OOC(F)(F)C(F)[C]F-2(510)'],
    products = ['OOC(F)(F)C(F)[C](F)F(1096)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(518506,'m^3/(mol*s)'), n=0.4717, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->H_N-2Br1sCl1sF1s->Cl1s_3BrCClFINOPSSi->F',), comment="""Estimated from node Root_N-3R->H_N-2Br1sCl1sF1s->Cl1s_3BrCClFINOPSSi->F"""),
)

reaction(
    label = 'reaction41',
    reactants = ['CF2(89)', 'OOC(F)(F)[CH]F(437)'],
    products = ['OOC(F)(F)C(F)[C](F)F(1096)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(786.723,'m^3/(mol*s)'), n=1.25031, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl"""),
)

reaction(
    label = 'reaction42',
    reactants = ['OOC(F)(F)C(F)[C](F)F(1096)'],
    products = ['OOC(F)(F)[C](F)C(F)F(495)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.70914e+10,'s^-1'), n=0.735063, Ea=(133.5,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_noH;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['OOC(F)(F)C(F)[C](F)F(1096)'],
    products = ['[O]OC(F)(F)C(F)C(F)F(1017)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2044.71,'s^-1'), n=2.0375, Ea=(28.5087,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;C_rad_out_single;XH_out] for rate rule [R5H_SSSS;C_rad_out_noH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['OOC(F)(F)C(F)[C](F)F(1096)'],
    products = ['[O]C(F)(F)C(F)C(O)(F)F(511)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(4.79e+10,'s^-1'), n=0, Ea=(110.039,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R3OOH_SS;C_rad_out_noH]
Euclidian distance = 0
family: intra_OH_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['OOC(F)(F)C(F)[C](F)F(1096)'],
    products = ['OOC(F)(F)[CH]C(F)(F)F(512)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(146.54,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction45',
    reactants = ['OOC(F)(F)C(F)[C](F)F(1096)'],
    products = ['OH(4)', 'O=C(F)C(F)C(F)(F)F(513)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(144.699,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

network(
    label = 'PDepNetwork #231',
    isomers = [
        'OOC(F)(F)C(F)[C](F)F(1096)',
    ],
    reactants = [
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #231',
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

