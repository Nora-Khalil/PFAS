species(
    label = 'F[C](F)C(F)C(F)[C](F)F(476)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {7,S} {10,S} {12,S}
9  C u1 p0 c0 {3,S} {4,S} {7,S}
10 C u1 p0 c0 {5,S} {6,S} {8,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-938.204,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([207,311,436,622,447,691,1074,1182,1276,1366,1311,1469,3045,3235,146,234,414,562,504,606,1176,1296,1354,1460,278.796,278.835,278.917],'cm^-1')),
        HinderedRotor(inertia=(0.227219,'amu*angstrom^2'), symmetry=1, barrier=(12.5354,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.226902,'amu*angstrom^2'), symmetry=1, barrier=(12.5249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.680775,'amu*angstrom^2'), symmetry=1, barrier=(37.6045,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (164.049,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.120225,0.0979875,-0.000148241,1.14506e-07,-3.50311e-11,-112698,30.7704], Tmin=(100,'K'), Tmax=(801.779,'K')), NASAPolynomial(coeffs=[14.2718,0.026185,-1.39062e-05,2.80487e-09,-2.00889e-13,-115006,-35.4826], Tmin=(801.779,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-938.204,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsF1sH)(F1s)(F1s)) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
)

species(
    label = 'CHFCF2(54)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 C u0 p0 c0 {1,S} {5,D} {6,S}
5 C u0 p0 c0 {2,S} {3,S} {4,D}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-505.825,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(82.003,'amu')),
        NonlinearRotor(inertia=([47.283,130.848,178.131],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([226.38,309.801,488.171,585.738,628.102,776.692,950.625,1195.27,1295.73,1386.74,1841.68,3239.84],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0245,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2091.09,'J/mol'), sigma=(4.442,'angstroms'), dipoleMoment=(1.4,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35521,0.0323113,-3.57904e-05,2.09758e-08,-5.1422e-12,-60645.8,19.2938], Tmin=(298,'K'), Tmax=(1100,'K')), NASAPolynomial(coeffs=[9.29782,0.00669066,-2.71841e-06,4.9759e-10,-3.35867e-14,-62705.1,-20.9391], Tmin=(1100,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-505.825,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(133.032,'J/mol/K'), label="""CHFCF2""", comment="""Thermo library: Fluorine"""),
)

species(
    label = 'F[CH]C(F)(F)C(F)[C](F)F(475)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
9  C u1 p0 c0 {4,S} {5,S} {7,S}
10 C u1 p0 c0 {6,S} {8,S} {12,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-934.787,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,215,315,519,588,595,1205,1248,190,488,555,1236,1407,334,575,1197,1424,3202,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.209985,'amu*angstrom^2'), symmetry=1, barrier=(4.82797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.209966,'amu*angstrom^2'), symmetry=1, barrier=(4.82754,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.605753,'amu*angstrom^2'), symmetry=1, barrier=(13.9274,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (164.049,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.543761,0.112742,-0.000205623,1.86499e-07,-6.44413e-11,-112278,31.8882], Tmin=(100,'K'), Tmax=(847.347,'K')), NASAPolynomial(coeffs=[11.6602,0.0304619,-1.62979e-05,3.18556e-09,-2.19646e-13,-113460,-19.7416], Tmin=(847.347,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-934.787,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsF1sH)(F1s)(F1s)) + radical(Csj(Cs-CsF1sF1s)(F1s)(H))"""),
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
    label = 'F[CH]C(F)[C](F)F(270)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {7,S}
4 F u0 p3 c0 {7,S}
5 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6 C u1 p0 c0 {2,S} {5,S} {9,S}
7 C u1 p0 c0 {3,S} {4,S} {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-506.234,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,334,575,1197,1424,3202,190,488,555,1236,1407,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.165162,'amu*angstrom^2'), symmetry=1, barrier=(3.7974,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165272,'amu*angstrom^2'), symmetry=1, barrier=(3.79992,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2384,0.0729047,-0.000142738,1.40492e-07,-5.15286e-11,-60798.1,23.0028], Tmin=(100,'K'), Tmax=(853.913,'K')), NASAPolynomial(coeffs=[5.03346,0.0276127,-1.48438e-05,2.90818e-09,-2.00628e-13,-60443.1,11.1669], Tmin=(853.913,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-506.234,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sCsH)(F1s)(H)) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
)

species(
    label = 'FC1C(F)C(F)(F)C1(F)F(591)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {7,S} {10,S} {12,S}
9  C u0 p0 c0 {3,S} {4,S} {7,S} {10,S}
10 C u0 p0 c0 {5,S} {6,S} {8,S} {9,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-1187.41,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (164.049,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.582478,0.072053,-7.34827e-05,3.74311e-08,-7.5451e-12,-142687,24.9602], Tmin=(100,'K'), Tmax=(1201.72,'K')), NASAPolynomial(coeffs=[16.0582,0.020541,-9.18452e-06,1.76083e-09,-1.24406e-13,-146406,-52.5448], Tmin=(1201.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1187.41,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + group(CsCsCsFH) + group(CsCsCsFF) + group(CsCsCsFF) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + ring(Cs-Cs-Cs(F)(F)-Cs)"""),
)

species(
    label = 'FC(F)=C(F)C(F)C(F)F(592)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {12,S}
9  C u0 p0 c0 {4,S} {7,S} {10,D}
10 C u0 p0 c0 {5,S} {6,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-1167.83,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (164.049,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.24763,0.0658696,-3.75703e-05,-6.86567e-08,8.14977e-11,-140453,14.6047], Tmin=(10,'K'), Tmax=(509.31,'K')), NASAPolynomial(coeffs=[9.65531,0.0390391,-2.77436e-05,9.05216e-09,-1.10452e-12,-141411,-14.9769], Tmin=(509.31,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1167.83,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), label="""FC(F)DC(F)C(F)C(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'F[C](F)C(F)C=C(F)F(593)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {6,S} {9,D} {11,S}
8  C u1 p0 c0 {2,S} {3,S} {6,S}
9  C u0 p0 c0 {4,S} {5,S} {7,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-812.438,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,3010,987.5,1337.5,450,1655,190,488,555,1236,1407,182,240,577,636,1210,1413,180],'cm^-1')),
        HinderedRotor(inertia=(0.587091,'amu*angstrom^2'), symmetry=1, barrier=(13.4984,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.33265,'amu*angstrom^2'), symmetry=1, barrier=(7.64829,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (145.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.32842,0.0664288,-9.68022e-05,8.22097e-08,-2.96234e-11,-97713,14.7494], Tmin=(10,'K'), Tmax=(647.847,'K')), NASAPolynomial(coeffs=[8.37817,0.0352502,-2.4613e-05,7.92368e-09,-9.57009e-13,-98367.3,-7.42059], Tmin=(647.847,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-812.438,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), label="""F[C](F)C(F)CDC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[CH][C](F)F(185)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 C u1 p0 c0 {1,S} {5,S} {6,S}
5 C u1 p0 c0 {2,S} {3,S} {4,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-282.469,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([334,575,1197,1424,3202,190,488,555,1236,1407,1399.55],'cm^-1')),
        HinderedRotor(inertia=(0.462103,'amu*angstrom^2'), symmetry=1, barrier=(10.6247,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.0245,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.7193,0.0310345,-4.23424e-05,3.17512e-08,-9.73807e-12,-33929.6,16.1924], Tmin=(100,'K'), Tmax=(789.997,'K')), NASAPolynomial(coeffs=[6.47706,0.0120087,-6.21892e-06,1.26852e-09,-9.20647e-14,-34523.4,-1.05094], Tmin=(789.997,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-282.469,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: Fluorine + radical(CsCsF1sH) + radical(CsCsF1sF1s)"""),
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
    label = 'F[C](F)C(F)C(F)=C(F)F(594)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {7,S} {10,D}
9  C u1 p0 c0 {3,S} {4,S} {7,S}
10 C u0 p0 c0 {5,S} {6,S} {8,D}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-963.659,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,323,467,575,827,1418,190,488,555,1236,1407,182,240,577,636,1210,1413,180],'cm^-1')),
        HinderedRotor(inertia=(0.348376,'amu*angstrom^2'), symmetry=1, barrier=(8.00986,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.196883,'amu*angstrom^2'), symmetry=1, barrier=(4.52672,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (163.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.26073,0.0765157,-0.000129888,1.29378e-07,-5.42634e-11,-115903,15.6322], Tmin=(10,'K'), Tmax=(564.614,'K')), NASAPolynomial(coeffs=[8.65564,0.0382955,-2.83493e-05,9.48603e-09,-1.17739e-12,-116513,-7.3113], Tmin=(564.614,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-963.659,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), label="""F[C](F)C(F)C(F)DC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F2(108)',
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
    label = 'F[C](F)C=C[C](F)F(595)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {6,D} {7,S} {9,S}
6  C u0 p0 c0 {5,D} {8,S} {10,S}
7  C u1 p0 c0 {1,S} {2,S} {5,S}
8  C u1 p0 c0 {3,S} {4,S} {6,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-514.446,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,124,198,231,363,414,566,474,694,695,865,1255,1461],'cm^-1')),
        HinderedRotor(inertia=(0.00185187,'amu*angstrom^2'), symmetry=1, barrier=(0.131987,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.752473,'amu*angstrom^2'), symmetry=1, barrier=(54.4753,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58262,0.0580861,-6.82311e-05,4.39241e-08,-1.17277e-11,-61790.7,24.0263], Tmin=(100,'K'), Tmax=(896.17,'K')), NASAPolynomial(coeffs=[9.02887,0.0248507,-1.2603e-05,2.54267e-09,-1.83904e-13,-63125.4,-11.0814], Tmin=(896.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-514.446,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cd-CdH)(F1s)(F1s)) + radical(Csj(Cd-CdH)(F1s)(F1s))"""),
)

species(
    label = 'HF(38)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 H u0 p0 c0 {1,S}
"""),
    E0 = (-281.113,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(20.0062,'amu')),
        LinearRotor(inertia=(0.809097,'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([4113.43],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (20.0064,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(2743.78,'J/mol'), sigma=(3.148,'angstroms'), dipoleMoment=(1.92,'De'), polarizability=(2.46,'angstroms^3'), rotrelaxcollnum=1.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.43657,0.000486021,-1.2524e-06,1.36475e-09,-4.09574e-13,-33800.1,1.20682], Tmin=(298,'K'), Tmax=(1250,'K')), NASAPolynomial(coeffs=[2.7813,0.00103959,-2.41735e-07,2.68416e-11,-1.09766e-15,-33504.2,5.0197], Tmin=(1250,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-281.113,'kJ/mol'), Cp0=(29.1007,'J/mol/K'), CpInf=(37.4151,'J/mol/K'), label="""HF""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'F[C](F)C=C(F)[C](F)F(596)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {7,D} {8,S} {10,S}
7  C u0 p0 c0 {1,S} {6,D} {9,S}
8  C u1 p0 c0 {2,S} {3,S} {6,S}
9  C u1 p0 c0 {4,S} {5,S} {7,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-684.449,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,271,519,563,612,1379,124,198,231,363,414,566,474,694,695,865,1255,1461],'cm^-1')),
        HinderedRotor(inertia=(0.00994083,'amu*angstrom^2'), symmetry=1, barrier=(8.1003,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.364,'amu*angstrom^2'), symmetry=1, barrier=(54.3529,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.664737,0.0802216,-0.000129315,1.11182e-07,-3.827e-11,-82206.7,27.7969], Tmin=(100,'K'), Tmax=(745.782,'K')), NASAPolynomial(coeffs=[10.293,0.0256991,-1.38572e-05,2.79206e-09,-1.98964e-13,-83562.7,-15.2924], Tmin=(745.782,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-684.449,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cd-CdH)(F1s)(F1s)) + radical(Csj(Cd-F1sCd)(F1s)(F1s))"""),
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
    label = 'F[C](F)C(F)[C](F)C(F)F(597)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {9,S} {12,S}
9  C u1 p0 c0 {4,S} {7,S} {8,S}
10 C u1 p0 c0 {5,S} {6,S} {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-941.163,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (164.049,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.108247,0.10003,-0.000169592,1.51007e-07,-5.27778e-11,-113057,32.8344], Tmin=(100,'K'), Tmax=(799.67,'K')), NASAPolynomial(coeffs=[10.9215,0.0317849,-1.7056e-05,3.39759e-09,-2.39287e-13,-114403,-15.2985], Tmin=(799.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-941.163,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sCsH)(Cs-F1sF1sH)(F1s)) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
)

species(
    label = 'F[C](F)[C](F)C(F)C(F)F(598)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {12,S}
9  C u1 p0 c0 {4,S} {7,S} {10,S}
10 C u1 p0 c0 {5,S} {6,S} {9,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-941.163,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (164.049,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.108247,0.10003,-0.000169592,1.51007e-07,-5.27778e-11,-113057,32.8344], Tmin=(100,'K'), Tmax=(799.67,'K')), NASAPolynomial(coeffs=[10.9215,0.0317849,-1.7056e-05,3.39759e-09,-2.39287e-13,-114403,-15.2985], Tmin=(799.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-941.163,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sCsH)(Cs-F1sF1sH)(F1s)) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
)

species(
    label = 'F[C](F)C(F)[CH]C(F)(F)F(599)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {4,S} {9,S}
9  C u1 p0 c0 {7,S} {8,S} {12,S}
10 C u1 p0 c0 {5,S} {6,S} {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-1004.83,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (164.049,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.038437,0.0960047,-0.000159621,1.40809e-07,-4.90771e-11,-120719,30.6516], Tmin=(100,'K'), Tmax=(788.777,'K')), NASAPolynomial(coeffs=[10.7977,0.0311132,-1.65742e-05,3.30394e-09,-2.33214e-13,-122095,-16.6657], Tmin=(788.777,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1004.83,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sCsH)(Cs-F1sF1sF1s)(H)) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
)

species(
    label = 'F[C](F)[CH]C(F)C(F)(F)F(600)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
9  C u1 p0 c0 {7,S} {10,S} {12,S}
10 C u1 p0 c0 {5,S} {6,S} {9,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-994.766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (164.049,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0806256,0.0956932,-0.000161047,1.44446e-07,-5.10032e-11,-119511,32.0239], Tmin=(100,'K'), Tmax=(796.736,'K')), NASAPolynomial(coeffs=[10.0212,0.0325917,-1.74048e-05,3.46742e-09,-2.44461e-13,-120676,-11.0465], Tmin=(796.736,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-994.766,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sCsH)(Cs-F1sF1sH)(H)) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
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
    E0 = (-385.36,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-224.625,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (80.337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-377.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-321.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-127.841,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-231.894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-198.156,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-12.0939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (47.2897,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-135.981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-150.289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-251.861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-263.003,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-235.019,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-185.435,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['F[C](F)C(F)C(F)[C](F)F(476)'],
    products = ['CHFCF2(54)', 'CHFCF2(54)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['F[CH]C(F)(F)C(F)[C](F)F(475)'],
    products = ['F[C](F)C(F)C(F)[C](F)F(476)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CF2(42)', 'F[CH]C(F)[C](F)F(270)'],
    products = ['F[C](F)C(F)C(F)[C](F)F(476)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_sec_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['F[C](F)C(F)C(F)[C](F)F(476)'],
    products = ['FC1C(F)C(F)(F)C1(F)F(591)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_noH;Cpri_rad_out_noH]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['F[C](F)C(F)C(F)[C](F)F(476)'],
    products = ['FC(F)=C(F)C(F)C(F)F(592)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['F(37)', 'F[C](F)C(F)C=C(F)F(593)'],
    products = ['F[C](F)C(F)C(F)[C](F)F(476)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(58.8616,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CHFCF2(54)', 'F[CH][C](F)F(185)'],
    products = ['F[C](F)C(F)C(F)[C](F)F(476)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(3.55601,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(5)', 'F[C](F)C(F)C(F)=C(F)F(594)'],
    products = ['F[C](F)C(F)C(F)[C](F)F(476)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(56.8734,'m^3/(mol*s)'), n=1.75834, Ea=(0.854637,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.2648472359903826, var=0.02782886759889551, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_4R!H->C_Sp-4C-1COS_Ext-1COS-R_N-6R!H-inRing_N-7R!H-inRing_Ext-4C-R_N-Sp-8R!H#4C',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_4R!H->C_Sp-4C-1COS_Ext-1COS-R_N-6R!H-inRing_N-7R!H-inRing_Ext-4C-R_N-Sp-8R!H#4C"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F[CH][C](F)F(185)', 'F[CH][C](F)F(185)'],
    products = ['F[C](F)C(F)C(F)[C](F)F(476)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.31566e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -9.3 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction10',
    reactants = ['F2(108)', 'F[C](F)C=C[C](F)F(595)'],
    products = ['F[C](F)C(F)C(F)[C](F)F(476)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(17.6971,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction11',
    reactants = ['HF(38)', 'F[C](F)C=C(F)[C](F)F(596)'],
    products = ['F[C](F)C(F)C(F)[C](F)F(476)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(276.738,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CF2(89)', 'F[CH]C(F)[C](F)F(270)'],
    products = ['F[C](F)C(F)C(F)[C](F)F(476)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(786.723,'m^3/(mol*s)'), n=1.25031, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F[C](F)C(F)C(F)[C](F)F(476)'],
    products = ['F[C](F)C(F)[C](F)C(F)F(597)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(5.41827e+10,'s^-1'), n=0.735063, Ea=(133.5,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_noH;Cs_H_out_noH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F[C](F)C(F)C(F)[C](F)F(476)'],
    products = ['F[C](F)[C](F)C(F)C(F)F(598)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.01948e+07,'s^-1'), n=1.33047, Ea=(122.357,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_noH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F[C](F)C(F)C(F)[C](F)F(476)'],
    products = ['F[C](F)C(F)[CH]C(F)(F)F(599)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.01526e+12,'s^-1'), n=0.18834, Ea=(150.341,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F[C](F)C(F)C(F)[C](F)F(476)'],
    products = ['F[C](F)[CH]C(F)C(F)(F)F(600)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(199.925,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

network(
    label = 'PDepNetwork #220',
    isomers = [
        'F[C](F)C(F)C(F)[C](F)F(476)',
    ],
    reactants = [
        ('CHFCF2(54)', 'CHFCF2(54)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #220',
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

