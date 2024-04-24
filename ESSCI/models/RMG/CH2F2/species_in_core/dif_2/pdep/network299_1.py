species(
    label = 'F[C](F)C(F)C(F)(F)[C](F)F(614)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {11,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {8,S} {11,S} {12,S}
10 C u1 p0 c0 {4,S} {5,S} {8,S}
11 C u1 p0 c0 {6,S} {7,S} {9,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-1133.72,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,259,529,569,1128,1321,1390,3140,146,234,414,562,504,606,1176,1296,1354,1460,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.302756,'amu*angstrom^2'), symmetry=1, barrier=(6.96095,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.59378,'amu*angstrom^2'), symmetry=1, barrier=(13.6522,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.300765,'amu*angstrom^2'), symmetry=1, barrier=(6.91519,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (182.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.966503,0.123611,-0.000230582,2.09333e-07,-7.19372e-11,-136190,33.9702], Tmin=(100,'K'), Tmax=(854.928,'K')), NASAPolynomial(coeffs=[12.8964,0.0301518,-1.64274e-05,3.20953e-09,-2.20359e-13,-137515,-24.625], Tmin=(854.928,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1133.72,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsF1sF1s)(F1s)(F1s)) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
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
    label = 'CHFCF2(54)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 C u0 p0 c0 {1,S} {5,D} {6,S}
5 C u0 p0 c0 {2,S} {3,S} {4,D}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-511.455,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93201,0.00413887,7.03233e-05,-1.49827e-07,9.42397e-11,-61509.2,9.50863], Tmin=(10,'K'), Tmax=(540.182,'K')), NASAPolynomial(coeffs=[3.82032,0.0197002,-1.38027e-05,4.49179e-09,-5.49698e-13,-61712.1,7.9889], Tmin=(540.182,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-511.455,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""FCDC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[CH]C(F)(F)C(F)(F)[C](F)F(613)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  F u0 p3 c0 {11,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {8,S} {11,S}
10 C u1 p0 c0 {5,S} {6,S} {8,S}
11 C u1 p0 c0 {7,S} {9,S} {12,S}
12 H u0 p0 c0 {11,S}
"""),
    E0 = (-1135.15,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (182.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.881098,0.125064,-0.000242529,2.31895e-07,-8.37412e-11,-136368,33.3204], Tmin=(100,'K'), Tmax=(837.489,'K')), NASAPolynomial(coeffs=[9.59408,0.0385328,-2.21736e-05,4.4463e-09,-3.10929e-13,-136843,-7.7165], Tmin=(837.489,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1135.15,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsF1sF1s)(F1s)(F1s)) + radical(Csj(Cs-CsF1sF1s)(F1s)(H))"""),
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
    label = 'F[C](F)C(F)[C](F)F(459)',
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
    label = 'F[CH]C(F)(F)[C](F)F(473)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {7,S}
4 F u0 p3 c0 {8,S}
5 F u0 p3 c0 {8,S}
6 C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7 C u1 p0 c0 {3,S} {6,S} {9,S}
8 C u1 p0 c0 {4,S} {5,S} {6,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-716.093,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,334,575,1197,1424,3202,190,488,555,1236,1407,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.279061,'amu*angstrom^2'), symmetry=1, barrier=(6.41616,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.279095,'amu*angstrom^2'), symmetry=1, barrier=(6.41695,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.615081,0.086223,-0.000166487,1.56139e-07,-5.49764e-11,-86015.5,25.0614], Tmin=(100,'K'), Tmax=(855.934,'K')), NASAPolynomial(coeffs=[8.68409,0.0235774,-1.30006e-05,2.55445e-09,-1.75869e-13,-86483.3,-7.27556], Tmin=(855.934,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-716.093,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsF1sF1s)(F1s)(H)) + radical(Csj(Cs-CsF1sF1s)(F1s)(F1s))"""),
)

species(
    label = 'FC1C(F)(F)C(F)(F)C1(F)F(756)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {11,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {8,S} {11,S} {12,S}
10 C u0 p0 c0 {4,S} {5,S} {8,S} {11,S}
11 C u0 p0 c0 {6,S} {7,S} {9,S} {10,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-1409.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (182.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.178052,0.0849522,-0.000101446,6.04112e-08,-1.4174e-11,-169336,26.9574], Tmin=(100,'K'), Tmax=(1040.5,'K')), NASAPolynomial(coeffs=[16.6219,0.0217364,-1.03125e-05,2.0194e-09,-1.44136e-13,-172758,-53.0271], Tmin=(1040.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1409.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + group(CsCsCsFH) + group(CsCsCsFF) + group(CsCsCsFF) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)2) + ring(Cs-Cs-Cs(F)(F)-Cs)"""),
)

species(
    label = 'FC(F)=C(F)C(F)(F)C(F)F(757)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {11,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {8,S} {12,S}
10 C u0 p0 c0 {5,S} {8,S} {11,D}
11 C u0 p0 c0 {6,S} {7,S} {10,D}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-1367.72,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (182.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.09548,0.0816355,-0.000115091,8.59847e-08,-2.61827e-11,-164494,15.5498], Tmin=(10,'K'), Tmax=(775.449,'K')), NASAPolynomial(coeffs=[12.2392,0.0344694,-2.38551e-05,7.54717e-09,-8.94865e-13,-165912,-26.2379], Tmin=(775.449,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1367.72,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), label="""FC(F)DC(F)C(F)(F)C(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'F[C](F)C(F)C(F)=C(F)F(487)',
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
    label = 'F[CH][C](F)F(183)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 C u1 p0 c0 {1,S} {5,S} {6,S}
5 C u1 p0 c0 {2,S} {3,S} {4,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-285.254,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([334,575,1197,1424,3202,190,488,555,1236,1407,1765.02],'cm^-1')),
        HinderedRotor(inertia=(0.510429,'amu*angstrom^2'), symmetry=1, barrier=(11.7358,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.0245,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.54626,0.0360262,-6.18601e-05,5.6897e-08,-2.03341e-11,-34259.7,17.1201], Tmin=(100,'K'), Tmax=(816.913,'K')), NASAPolynomial(coeffs=[5.79757,0.0130515,-6.72074e-06,1.3276e-09,-9.30775e-14,-34555.5,3.53266], Tmin=(816.913,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-285.254,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCsF1sH) + radical(CsCsF1sF1s)"""),
)

species(
    label = 'F[C](F)C(F)(F)C=C(F)F(758)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {7,S} {10,D} {11,S}
9  C u1 p0 c0 {3,S} {4,S} {7,S}
10 C u0 p0 c0 {5,S} {6,S} {8,D}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-1020.45,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,3010,987.5,1337.5,450,1655,190,488,555,1236,1407,182,240,577,636,1210,1413,180],'cm^-1')),
        HinderedRotor(inertia=(0.387794,'amu*angstrom^2'), symmetry=1, barrier=(8.91614,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.388716,'amu*angstrom^2'), symmetry=1, barrier=(8.93736,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (163.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.26277,0.0742697,-0.00011564,1.01609e-07,-3.72755e-11,-122731,15.4668], Tmin=(10,'K'), Tmax=(642.128,'K')), NASAPolynomial(coeffs=[9.4219,0.0359026,-2.60146e-05,8.55839e-09,-1.04798e-12,-123522,-11.5191], Tmin=(642.128,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1020.45,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), label="""F[C](F)C(F)(F)CDC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'F[C](F)C=C(F)[C](F)F(489)',
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
    label = 'F[C](F)C(F)(F)[C](F)C(F)F(761)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {11,S}
8  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {10,S} {12,S}
10 C u1 p0 c0 {5,S} {8,S} {9,S}
11 C u1 p0 c0 {6,S} {7,S} {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-1142.07,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (182.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.966108,0.12524,-0.000238853,2.20465e-07,-7.65475e-11,-137196,34.7299], Tmin=(100,'K'), Tmax=(859.711,'K')), NASAPolynomial(coeffs=[11.8908,0.0318842,-1.74548e-05,3.40648e-09,-2.33315e-13,-138167,-18.1459], Tmin=(859.711,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1142.07,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(F1s)) + radical(Csj(Cs-CsF1sF1s)(F1s)(F1s))"""),
)

species(
    label = 'F[C](F)[C](F)C(F)(F)C(F)F(762)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {11,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {8,S} {12,S}
10 C u1 p0 c0 {5,S} {8,S} {11,S}
11 C u1 p0 c0 {6,S} {7,S} {10,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-1143.37,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (182.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.82895,0.122787,-0.000235807,2.20309e-07,-7.73106e-11,-137357,34.9678], Tmin=(100,'K'), Tmax=(858.357,'K')), NASAPolynomial(coeffs=[10.7178,0.0336945,-1.8456e-05,3.60814e-09,-2.47593e-13,-138040,-11.4037], Tmin=(858.357,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1143.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(F1s)) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
)

species(
    label = 'F[C](F)C(F)[C](F)C(F)(F)F(763)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {11,S}
8  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
9  C u0 p0 c0 {2,S} {3,S} {4,S} {10,S}
10 C u1 p0 c0 {5,S} {8,S} {9,S}
11 C u1 p0 c0 {6,S} {7,S} {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-1170.27,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (182.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.311869,0.108038,-0.000196086,1.83228e-07,-6.589e-11,-140608,34.7535], Tmin=(100,'K'), Tmax=(821.089,'K')), NASAPolynomial(coeffs=[9.4414,0.0360834,-1.99865e-05,3.9953e-09,-2.8059e-13,-141386,-5.36068], Tmin=(821.089,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1170.27,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sCsH)(Cs-F1sF1sF1s)(F1s)) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
)

species(
    label = 'F[C](F)[C](F)C(F)C(F)(F)F(764)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {11,S}
8  C u0 p0 c0 {1,S} {9,S} {10,S} {12,S}
9  C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
10 C u1 p0 c0 {5,S} {8,S} {11,S}
11 C u1 p0 c0 {6,S} {7,S} {10,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-1165.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (182.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.408504,0.115146,-0.000225956,2.18599e-07,-7.88731e-11,-140041,34.639], Tmin=(100,'K'), Tmax=(858.09,'K')), NASAPolynomial(coeffs=[7.28623,0.0385887,-2.10013e-05,4.10656e-09,-2.82235e-13,-139863,7.42263], Tmin=(858.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1165.54,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sCsH)(Cs-F1sF1sH)(F1s)) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
)

species(
    label = 'F[C](F)[CH]C(F)(F)C(F)(F)F(765)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {11,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
10 C u1 p0 c0 {8,S} {11,S} {12,S}
11 C u1 p0 c0 {6,S} {7,S} {10,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-1206.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (182.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.277209,0.107841,-0.000197542,1.86336e-07,-6.74444e-11,-144970,34.5226], Tmin=(100,'K'), Tmax=(823.262,'K')), NASAPolynomial(coeffs=[8.85117,0.0371063,-2.05909e-05,4.11611e-09,-2.88954e-13,-145579,-2.31122], Tmin=(823.262,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1206.52,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H)) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = 'F[C](F)C(F)(F)[CH]C(F)(F)F(766)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {11,S}
8  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {5,S} {10,S}
10 C u1 p0 c0 {8,S} {9,S} {12,S}
11 C u1 p0 c0 {6,S} {7,S} {8,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-1211.4,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (182.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.734328,0.120253,-0.000229184,2.1629e-07,-7.7315e-11,-145543,33.9767], Tmin=(100,'K'), Tmax=(839.186,'K')), NASAPolynomial(coeffs=[10.0782,0.0360751,-2.0378e-05,4.0626e-09,-2.83268e-13,-146208,-9.444], Tmin=(839.186,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1211.4,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sF1s)(H)) + radical(Csj(Cs-CsF1sF1s)(F1s)(F1s))"""),
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
    E0 = (-480.652,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-320.717,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-24.8633,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-29.3004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-472.368,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-417.252,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-190.576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-318.083,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-231.512,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-345.532,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-298.809,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-125.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-25.608,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-241.608,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-262.303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-266.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-347.153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-358.295,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-317.174,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-269.316,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-288.774,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-335.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['F[C](F)C(F)C(F)(F)[C](F)F(614)'],
    products = ['CF2CF2(60)', 'CHFCF2(54)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['F[C](F)C(F)C(F)(F)[C](F)F(614)'],
    products = ['F[CH]C(F)(F)C(F)(F)[C](F)F(613)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CF2(42)', 'F[C](F)C(F)[C](F)F(459)'],
    products = ['F[C](F)C(F)C(F)(F)[C](F)F(614)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4.0899e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['CF2(42)', 'F[CH]C(F)(F)[C](F)F(473)'],
    products = ['F[C](F)C(F)C(F)(F)[C](F)F(614)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_sec_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['F[C](F)C(F)C(F)(F)[C](F)F(614)'],
    products = ['FC1C(F)(F)C(F)(F)C1(F)F(756)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_noH;Cpri_rad_out_noH]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['F[C](F)C(F)C(F)(F)[C](F)F(614)'],
    products = ['FC(F)=C(F)C(F)(F)C(F)F(757)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F(37)', 'F[C](F)C(F)C(F)=C(F)F(487)'],
    products = ['F[C](F)C(F)C(F)(F)[C](F)F(614)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(47.127,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CF2CF2(60)', 'F[CH][C](F)F(183)'],
    products = ['F[C](F)C(F)C(F)(F)[C](F)F(614)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(8.08706e-06,'m^3/(mol*s)'), n=3.0961, Ea=(2.63964,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F(37)', 'F[C](F)C(F)(F)C=C(F)F(758)'],
    products = ['F[C](F)C(F)C(F)(F)[C](F)F(614)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(62.9765,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F[C](F)[C](F)F(230)', 'CHFCF2(54)'],
    products = ['F[C](F)C(F)C(F)(F)[C](F)F(614)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(8.08706e-06,'m^3/(mol*s)'), n=3.0961, Ea=(6.36902,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(5)', 'F[C](F)C(F)(F)C(F)=C(F)F(759)'],
    products = ['F[C](F)C(F)C(F)(F)[C](F)F(614)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(56.8734,'m^3/(mol*s)'), n=1.75834, Ea=(1.21837,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.2648472359903826, var=0.02782886759889551, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_4R!H->C_Sp-4C-1COS_Ext-1COS-R_N-6R!H-inRing_N-7R!H-inRing_Ext-4C-R_N-Sp-8R!H#4C',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_4R!H->C_Sp-4C-1COS_Ext-1COS-R_N-6R!H-inRing_N-7R!H-inRing_Ext-4C-R_N-Sp-8R!H#4C"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F[C](F)[C](F)F(230)', 'F[CH][C](F)F(183)'],
    products = ['F[C](F)C(F)C(F)(F)[C](F)F(614)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(5.26262e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0
Ea raised from -3.9 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['F2(106)', 'F[C](F)C=C(F)[C](F)F(489)'],
    products = ['F[C](F)C(F)C(F)(F)[C](F)F(614)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(14.581,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction14',
    reactants = ['HF(38)', 'F[C](F)C(F)=C(F)[C](F)F(760)'],
    products = ['F[C](F)C(F)C(F)(F)[C](F)F(614)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(8.28222,'m^3/(mol*s)'), n=1.29695, Ea=(213.956,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CF2(87)', 'F[C](F)C(F)[C](F)F(459)'],
    products = ['F[C](F)C(F)C(F)(F)[C](F)F(614)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1573.45,'m^3/(mol*s)'), n=1.25031, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction16',
    reactants = ['CF2(87)', 'F[CH]C(F)(F)[C](F)F(473)'],
    products = ['F[C](F)C(F)C(F)(F)[C](F)F(614)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(786.723,'m^3/(mol*s)'), n=1.25031, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F[C](F)C(F)C(F)(F)[C](F)F(614)'],
    products = ['F[C](F)C(F)(F)[C](F)C(F)F(761)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.70914e+10,'s^-1'), n=0.735063, Ea=(133.5,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_noH;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['F[C](F)C(F)C(F)(F)[C](F)F(614)'],
    products = ['F[C](F)[C](F)C(F)(F)C(F)F(762)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.50974e+07,'s^-1'), n=1.33047, Ea=(122.357,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_noH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['F[C](F)C(F)C(F)(F)[C](F)F(614)'],
    products = ['F[C](F)C(F)[C](F)C(F)(F)F(763)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.01526e+12,'s^-1'), n=0.18834, Ea=(163.478,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction20',
    reactants = ['F[C](F)C(F)C(F)(F)[C](F)F(614)'],
    products = ['F[C](F)[C](F)C(F)C(F)(F)F(764)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(211.336,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction21',
    reactants = ['F[C](F)C(F)C(F)(F)[C](F)F(614)'],
    products = ['F[C](F)[CH]C(F)(F)C(F)(F)F(765)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(191.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction22',
    reactants = ['F[C](F)C(F)C(F)(F)[C](F)F(614)'],
    products = ['F[C](F)C(F)(F)[CH]C(F)(F)F(766)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(145.597,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

network(
    label = 'PDepNetwork #299',
    isomers = [
        'F[C](F)C(F)C(F)(F)[C](F)F(614)',
    ],
    reactants = [
        ('CF2CF2(60)', 'CHFCF2(54)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #299',
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

