species(
    label = 'F[CH]C(F)(F)C(F)(F)[C](F)F(672)',
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
    modes = [
        HarmonicOscillator(frequencies=([181,249,286,344,486,552,511,665,553,637,1169,1241,1190,1306,190,488,555,1236,1407,334,575,1197,1424,3202,180,180,1755.07],'cm^-1')),
        HinderedRotor(inertia=(0.408626,'amu*angstrom^2'), symmetry=1, barrier=(9.39512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.407956,'amu*angstrom^2'), symmetry=1, barrier=(9.37971,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.95469,'amu*angstrom^2'), symmetry=1, barrier=(44.9422,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (182.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.881098,0.125064,-0.000242529,2.31895e-07,-8.37412e-11,-136368,33.3204], Tmin=(100,'K'), Tmax=(837.489,'K')), NASAPolynomial(coeffs=[9.59408,0.0385328,-2.21736e-05,4.4463e-09,-3.10929e-13,-136843,-7.7165], Tmin=(837.489,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1135.15,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsF1sF1s)(F1s)(F1s)) + radical(Csj(Cs-CsF1sF1s)(F1s)(H))"""),
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
    E0 = (-675.664,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.59735,0.0339059,-4.12239e-05,2.60507e-08,-6.78805e-12,-81178.9,12.7972], Tmin=(298,'K'), Tmax=(1100,'K')), NASAPolynomial(coeffs=[11.3243,0.00499631,-2.11782e-06,3.98619e-10,-2.75039e-14,-83426.6,-31.2702], Tmin=(1100,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-675.664,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(133.032,'J/mol/K'), label="""CF2CF2""", comment="""Thermo library: Fluorine"""),
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
    label = 'F[C](F)C(F)C(F)(F)[C](F)F(673)',
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(2970.8,'J/mol'), sigma=(5.0835,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=464.03 K, Pc=51.31 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.966503,0.123611,-0.000230582,2.09333e-07,-7.19372e-11,-136190,33.9702], Tmin=(100,'K'), Tmax=(854.928,'K')), NASAPolynomial(coeffs=[12.8964,0.0301518,-1.64274e-05,3.20953e-09,-2.20359e-13,-137515,-24.625], Tmin=(854.928,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1133.72,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsF1sF1s)(F1s)(F1s)) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
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
    label = 'F[CH]C(F)(F)[C](F)F(394)',
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
    label = 'CHF(40)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u2 p0 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (218.467,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([952.939,954.11,4000],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (32.017,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2402.79,'J/mol'), sigma=(4.33176,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=375.31 K, Pc=67.08 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.92903,-0.00113288,1.33435e-05,-1.59096e-08,5.65026e-12,26280.4,4.5427], Tmin=(100,'K'), Tmax=(1017.48,'K')), NASAPolynomial(coeffs=[5.22566,0.00101071,-4.91442e-07,1.4946e-10,-1.404e-14,25641.7,-3.57711], Tmin=(1017.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(218.467,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo library: Fluorine + radical(CH2_triplet)"""),
)

species(
    label = 'F[C](F)C(F)(F)[C](F)F(617)',
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
    label = 'FC1C(F)(F)C(F)(F)C1(F)F(769)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.178052,0.0849522,-0.000101446,6.04112e-08,-1.4174e-11,-169336,26.9574], Tmin=(100,'K'), Tmax=(1040.5,'K')), NASAPolynomial(coeffs=[16.6219,0.0217364,-1.03125e-05,2.0194e-09,-1.44136e-13,-172758,-53.0271], Tmin=(1040.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1409.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + group(CsCsCsFH) + group(CsCsCsFF) + group(CsCsCsFF) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)2) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)2) + ring(Cs-Cs-Cs(F)(F)-Cs)"""),
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
    label = 'F[CH]C(F)(F)C(F)=C(F)F(605)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {3,S} {7,S} {10,D}
9  C u1 p0 c0 {4,S} {7,S} {11,S}
10 C u0 p0 c0 {5,S} {6,S} {8,D}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-961.307,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,323,467,575,827,1418,334,575,1197,1424,3202,182,240,577,636,1210,1413,180],'cm^-1')),
        HinderedRotor(inertia=(0.191665,'amu*angstrom^2'), symmetry=1, barrier=(4.40677,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.190372,'amu*angstrom^2'), symmetry=1, barrier=(4.37702,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (163.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.12252,0.0816027,-0.000148704,1.53156e-07,-6.46569e-11,-115617,16.0595], Tmin=(10,'K'), Tmax=(569.52,'K')), NASAPolynomial(coeffs=[9.82035,0.0345609,-2.4806e-05,8.12442e-09,-9.92988e-13,-116380,-12.483], Tmin=(569.52,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-961.307,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), label="""F[CH]C(F)(F)C(F)DC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'FC=C(F)C(F)(F)[C](F)F(780)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {10,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {3,S} {7,S} {10,D}
9  C u1 p0 c0 {4,S} {5,S} {7,S}
10 C u0 p0 c0 {6,S} {8,D} {11,S}
11 H u0 p0 c0 {10,S}
"""),
    E0 = (-983.778,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,323,467,575,827,1418,190,488,555,1236,1407,194,682,905,1196,1383,3221,180],'cm^-1')),
        HinderedRotor(inertia=(0.340983,'amu*angstrom^2'), symmetry=1, barrier=(7.83987,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.340749,'amu*angstrom^2'), symmetry=1, barrier=(7.83449,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (163.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.16636,0.0761475,-0.000114894,9.24254e-08,-3.02088e-11,-118317,15.5698], Tmin=(10,'K'), Tmax=(728.023,'K')), NASAPolynomial(coeffs=[11.4231,0.030782,-2.1424e-05,6.83233e-09,-8.16457e-13,-119519,-21.6434], Tmin=(728.023,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-983.778,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), label="""FCDC(F)C(F)(F)[C](F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[C](F)[C](F)F(315)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {6,S}
4 F u0 p3 c0 {6,S}
5 C u1 p0 c0 {1,S} {2,S} {6,S}
6 C u1 p0 c0 {3,S} {4,S} {5,S}
"""),
    E0 = (-489.838,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([146,234,414,562,504,606,1176,1296,1354,1460,1701.96],'cm^-1')),
        HinderedRotor(inertia=(0.556297,'amu*angstrom^2'), symmetry=1, barrier=(12.7904,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.015,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.42114,0.0380951,-5.64738e-05,4.33909e-08,-1.33257e-11,-58860.1,18.3256], Tmin=(100,'K'), Tmax=(795.463,'K')), NASAPolynomial(coeffs=[7.71944,0.011453,-6.23605e-06,1.28829e-09,-9.38394e-14,-59703.1,-6.02335], Tmin=(795.463,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-489.838,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: Fluorine + radical(CsCsF1sF1s) + radical(CsCsF1sF1s)"""),
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
    label = 'F[CH]C(F)=C(F)[C](F)F(608)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {7,D} {8,S}
7  C u0 p0 c0 {2,S} {6,D} {9,S}
8  C u1 p0 c0 {3,S} {6,S} {10,S}
9  C u1 p0 c0 {4,S} {5,S} {7,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-643.189,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([206,336,431,607,515,611,528,696,1312,1446,234,589,736,816,1240,3237,161,297,490,584,780,1358],'cm^-1')),
        HinderedRotor(inertia=(0.0104319,'amu*angstrom^2'), symmetry=1, barrier=(58.8204,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146401,'amu*angstrom^2'), symmetry=1, barrier=(58.8197,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (144.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.849592,0.075931,-0.000113259,9.05043e-08,-2.92264e-11,-77250.4,26.7006], Tmin=(100,'K'), Tmax=(755.307,'K')), NASAPolynomial(coeffs=[10.2903,0.0259276,-1.39424e-05,2.83166e-09,-2.03783e-13,-78676.4,-16.195], Tmin=(755.307,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-643.189,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cd-CdF1s)(F1s)(H)) + radical(Csj(Cd-F1sCd)(F1s)(F1s))"""),
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
    label = 'CHF(88)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 C u0 p1 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (153.289,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1169.21,1416.41,2978.06],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (32.017,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2178.39,'J/mol'), sigma=(4.123,'angstroms'), dipoleMoment=(1.8,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.34484,0.00235461,1.93983e-06,-2.65251e-09,7.91169e-13,18514,7.05286], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[4.48366,0.00174964,-5.0479e-07,1.08953e-10,-9.87898e-15,17958.1,0.289222], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(153.289,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CHF""", comment="""Thermo library: Fluorine"""),
)

species(
    label = 'F[CH]C(F)(F)[C](F)C(F)(F)F(781)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {10,S}
7  F u0 p3 c0 {11,S}
8  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {5,S} {10,S}
10 C u1 p0 c0 {6,S} {8,S} {9,S}
11 C u1 p0 c0 {7,S} {8,S} {12,S}
12 H u0 p0 c0 {11,S}
"""),
    E0 = (-1176.22,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (182.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.76134,0.119069,-0.000221186,2.04318e-07,-7.19857e-11,-141309,33.3151], Tmin=(100,'K'), Tmax=(835.098,'K')), NASAPolynomial(coeffs=[11.3999,0.0336093,-1.88114e-05,3.7453e-09,-2.61271e-13,-142391,-17.4838], Tmin=(835.098,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1176.22,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sF1s)(F1s)) + radical(Csj(Cs-CsF1sF1s)(F1s)(H))"""),
)

species(
    label = 'F[C](F)[C](F)C(F)(F)C(F)F(775)',
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
    label = 'F[CH][C](F)C(F)(F)C(F)(F)F(782)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {10,S}
7  F u0 p3 c0 {11,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
10 C u1 p0 c0 {6,S} {8,S} {11,S}
11 C u1 p0 c0 {7,S} {10,S} {12,S}
12 H u0 p0 c0 {11,S}
"""),
    E0 = (-1180.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (182.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.707981,0.11772,-0.000217911,2.01522e-07,-7.11858e-11,-141784,33.3052], Tmin=(100,'K'), Tmax=(832.711,'K')), NASAPolynomial(coeffs=[11.165,0.0339949,-1.90108e-05,3.7886e-09,-2.64629e-13,-142836,-16.2449], Tmin=(832.711,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1180.16,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sF1sCs)(Cs-F1sHH)(F1s)) + radical(Csj(Cs-F1sCsH)(F1s)(H))"""),
)

species(
    label = 'F[C](F)C(F)(F)[C](F)C(F)F(774)',
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
    E0 = (-529.447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-368.079,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-76.6623,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-93.7086,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-521.163,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-236.631,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-351.098,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-253.217,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-384.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-166.604,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-36.078,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-307.288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-158.887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-367.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-306.089,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-324.125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-350.646,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['F[CH]C(F)(F)C(F)(F)[C](F)F(672)'],
    products = ['CF2CF2(60)', 'CHFCF2(54)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['F[C](F)C(F)C(F)(F)[C](F)F(673)'],
    products = ['F[CH]C(F)(F)C(F)(F)[C](F)F(672)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CF2(42)', 'F[CH]C(F)(F)[C](F)F(394)'],
    products = ['F[CH]C(F)(F)C(F)(F)[C](F)F(672)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['CHF(40)', 'F[C](F)C(F)(F)[C](F)F(617)'],
    products = ['F[CH]C(F)(F)C(F)(F)[C](F)F(672)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4.0899e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['F[CH]C(F)(F)C(F)(F)[C](F)F(672)'],
    products = ['FC1C(F)(F)C(F)(F)C1(F)F(769)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_noH;Cpri_rad_out_1H]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['F(37)', 'F[CH]C(F)(F)C(F)=C(F)F(605)'],
    products = ['F[CH]C(F)(F)C(F)(F)[C](F)F(672)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(46.0815,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CF2CF2(60)', 'F[CH][C](F)F(185)'],
    products = ['F[CH]C(F)(F)C(F)(F)[C](F)F(672)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(8.08706e-06,'m^3/(mol*s)'), n=3.0961, Ea=(1.33149,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F(37)', 'FC=C(F)C(F)(F)[C](F)F(780)'],
    products = ['F[CH]C(F)(F)C(F)(F)[C](F)F(672)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(51.9658,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F[C](F)[C](F)F(315)', 'CHFCF2(54)'],
    products = ['F[CH]C(F)(F)C(F)(F)[C](F)F(672)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(8.08706e-06,'m^3/(mol*s)'), n=3.0961, Ea=(5.04301,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F[C](F)[C](F)F(315)', 'F[CH][C](F)F(185)'],
    products = ['F[CH]C(F)(F)C(F)(F)[C](F)F(672)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.26262e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0
Ea raised from -5.7 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction11',
    reactants = ['F2(108)', 'F[CH]C(F)=C(F)[C](F)F(608)'],
    products = ['F[CH]C(F)(F)C(F)(F)[C](F)F(672)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(10.2122,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CF2(89)', 'F[CH]C(F)(F)[C](F)F(394)'],
    products = ['F[CH]C(F)(F)C(F)(F)[C](F)F(672)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(786.723,'m^3/(mol*s)'), n=1.25031, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CHF(88)', 'F[C](F)C(F)(F)[C](F)F(617)'],
    products = ['F[CH]C(F)(F)C(F)(F)[C](F)F(672)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1573.45,'m^3/(mol*s)'), n=1.25031, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F[CH]C(F)(F)C(F)(F)[C](F)F(672)'],
    products = ['F[CH]C(F)(F)[C](F)C(F)(F)F(781)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.01526e+12,'s^-1'), n=0.18834, Ea=(162.03,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F[CH]C(F)(F)C(F)(F)[C](F)F(672)'],
    products = ['F[C](F)[C](F)C(F)(F)C(F)F(775)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(223.358,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F[CH]C(F)(F)C(F)(F)[C](F)F(672)'],
    products = ['F[CH][C](F)C(F)(F)C(F)(F)F(782)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(205.323,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F[CH]C(F)(F)C(F)(F)[C](F)F(672)'],
    products = ['F[C](F)C(F)(F)[C](F)C(F)F(774)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.01526e+12,'s^-1'), n=0.18834, Ea=(178.802,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

network(
    label = 'PDepNetwork #279',
    isomers = [
        'F[CH]C(F)(F)C(F)(F)[C](F)F(672)',
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
    label = 'PDepNetwork #279',
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

