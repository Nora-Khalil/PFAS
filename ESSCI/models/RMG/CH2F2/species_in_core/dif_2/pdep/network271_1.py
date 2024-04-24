species(
    label = '[O]C(F)C(F)(F)[CH]F(590)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  O u1 p2 c0 {7,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {3,S} {5,S} {6,S} {9,S}
8  C u1 p0 c0 {4,S} {6,S} {10,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-678.425,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,391,562,707,872,1109,1210,1289,3137,334,575,1197,1424,3202,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.56927,'amu*angstrom^2'), symmetry=1, barrier=(36.0807,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.57028,'amu*angstrom^2'), symmetry=1, barrier=(36.1038,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16197,0.065254,-7.29524e-05,4.01401e-08,-8.81858e-12,-81495.8,23.0346], Tmin=(100,'K'), Tmax=(1095.55,'K')), NASAPolynomial(coeffs=[13.6313,0.0197269,-1.06178e-05,2.20804e-09,-1.6265e-13,-84228,-38.2605], Tmin=(1095.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-678.425,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sH)) + radical(Csj(Cs-CsF1sF1s)(F1s)(H))"""),
)

species(
    label = 'CHFO(46)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 H u0 p0 c0 {3,S}
"""),
    E0 = (-394.294,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(48.0011,'amu')),
        NonlinearRotor(inertia=([5.46358,42.889,48.3526],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([674.139,1050.3,1121.98,1395.32,1906.92,3070.52],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (48.0164,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2914.22,'J/mol'), sigma=(4.906,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05425,-0.00375329,3.18277e-05,-4.07297e-08,1.68956e-11,-47423,6.56837], Tmin=(10,'K'), Tmax=(745.315,'K')), NASAPolynomial(coeffs=[2.27714,0.0104751,-6.24845e-06,1.77292e-09,-1.93455e-13,-47288.4,13.7455], Tmin=(745.315,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-394.294,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""ODCF""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]C(F)C(F)[C](F)F(591)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u1 p2 c0 {7,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {2,S} {5,S} {6,S} {10,S}
8  C u1 p0 c0 {3,S} {4,S} {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-671.557,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,391,562,707,872,1109,1210,1289,3137,190,488,555,1236,1407,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.330503,'amu*angstrom^2'), symmetry=1, barrier=(7.59891,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.53103,'amu*angstrom^2'), symmetry=1, barrier=(35.2014,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.041,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3403.26,'J/mol'), sigma=(5.47543,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=531.58 K, Pc=47.04 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11081,0.0681287,-9.1777e-05,6.42432e-08,-1.79738e-11,-80669.6,24.598], Tmin=(100,'K'), Tmax=(871.53,'K')), NASAPolynomial(coeffs=[11.3938,0.0209344,-1.05518e-05,2.11197e-09,-1.51666e-13,-82462.1,-23.5976], Tmin=(871.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-671.557,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sH)) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
)

species(
    label = 'CHF(40)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u2 p0 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (214.928,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([787.277,1376.72,4000],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (32.017,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2402.79,'J/mol'), sigma=(4.33176,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=375.31 K, Pc=67.08 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93333,-0.00026336,8.89187e-06,-1.03032e-08,3.50809e-12,25853.7,4.33729], Tmin=(100,'K'), Tmax=(1056.12,'K')), NASAPolynomial(coeffs=[4.72426,0.00164131,-7.73114e-07,1.90987e-10,-1.59925e-14,25413.4,-0.815528], Tmin=(1056.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.928,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CH2_triplet)"""),
)

species(
    label = '[O]C(F)[C](F)F(386)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {6,S}
4 O u1 p2 c0 {5,S}
5 C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
6 C u1 p0 c0 {2,S} {3,S} {5,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-445.931,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([391,562,707,872,1109,1210,1289,3137,190,488,555,1236,1407,180],'cm^-1')),
        HinderedRotor(inertia=(0.813996,'amu*angstrom^2'), symmetry=1, barrier=(18.7154,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0239,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.14839,0.0397111,-4.32916e-05,2.26821e-08,-4.63769e-12,-53565.5,18.7797], Tmin=(100,'K'), Tmax=(1191.22,'K')), NASAPolynomial(coeffs=[11.353,0.00880293,-4.37161e-06,9.00449e-10,-6.63946e-14,-55758.5,-27.2377], Tmin=(1191.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-445.931,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sH)) + radical(Csj(Cs-F1sO2sH)(F1s)(F1s))"""),
)

species(
    label = 'O(6)',
    structure = adjacencyList("""multiplicity 3
1 O u2 p2 c0
"""),
    E0 = (243.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (15.9994,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.05029e-14,-1.20495e-16,5.23432e-20,-6.99814e-24,29230.2,5.12616], Tmin=(100,'K'), Tmax=(4230.24,'K')), NASAPolynomial(coeffs=[2.5,3.58763e-09,-1.22015e-12,1.84118e-16,-1.04e-20,29230.2,5.12619], Tmin=(4230.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.034,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'F[CH]C(F)(F)[CH]F(494)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {6,S}
4 F u0 p3 c0 {7,S}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6 C u1 p0 c0 {3,S} {5,S} {8,S}
7 C u1 p0 c0 {4,S} {5,S} {9,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-512.117,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,262,406,528,622,1148,1246,1368,1480,3164,3240,296.521,296.532],'cm^-1')),
        HinderedRotor(inertia=(0.140083,'amu*angstrom^2'), symmetry=1, barrier=(8.73776,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.140016,'amu*angstrom^2'), symmetry=1, barrier=(8.73754,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.959027,0.0727694,-0.000123276,1.05179e-07,-3.48434e-11,-61489.6,22.3742], Tmin=(100,'K'), Tmax=(820.048,'K')), NASAPolynomial(coeffs=[10.8922,0.0172827,-8.91382e-06,1.74493e-09,-1.21312e-13,-62882.2,-22.1345], Tmin=(820.048,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-512.117,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsF1sF1s)(F1s)(H)) + radical(Csj(Cs-CsF1sF1s)(F1s)(H))"""),
)

species(
    label = 'FC1OC(F)C1(F)F(595)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {7,S} {8,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {3,S} {5,S} {6,S} {9,S}
8  C u0 p0 c0 {4,S} {5,S} {6,S} {10,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-948.542,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.79771,0.0212497,5.86956e-05,-1.24581e-07,6.69727e-11,-114081,11.7873], Tmin=(10,'K'), Tmax=(670.753,'K')), NASAPolynomial(coeffs=[5.21407,0.035376,-2.33744e-05,7.16023e-09,-8.2925e-13,-114779,1.73458], Tmin=(670.753,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-948.542,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), label="""FC1OC(F)C1(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C(F)C(F)(F)CF(698)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {3,S} {6,S} {9,S} {10,S}
8  C u0 p0 c0 {4,S} {5,D} {6,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-1032.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.4409,0.0561789,-8.40319e-05,7.82471e-08,-3.16671e-11,-124168,12.659], Tmin=(10,'K'), Tmax=(575.789,'K')), NASAPolynomial(coeffs=[6.82319,0.0326822,-2.28201e-05,7.37404e-09,-8.94985e-13,-124558,-1.79148], Tmin=(575.789,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1032.39,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""ODC(F)C(F)(F)CF""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]C(F)C(F)=CF(699)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {7,S}
4 O u1 p2 c0 {5,S}
5 C u0 p0 c0 {1,S} {4,S} {6,S} {8,S}
6 C u0 p0 c0 {2,S} {5,S} {7,D}
7 C u0 p0 c0 {3,S} {6,D} {9,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-494.125,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([361,769,965,1078,1132,1246,3247,323,467,575,827,1418,194,682,905,1196,1383,3221,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.165178,'amu*angstrom^2'), symmetry=1, barrier=(3.79776,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.70591,0.0282874,0.00014724,-7.81679e-07,1.06059e-09,-59432.1,11.9956], Tmin=(10,'K'), Tmax=(275.514,'K')), NASAPolynomial(coeffs=[5.56266,0.0302683,-2.1093e-05,6.86002e-09,-8.41495e-13,-59644.2,3.43828], Tmin=(275.514,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-494.125,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""[O]C(F)C(F)DCF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O][CH]F(168)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 O u1 p2 c0 {3,S}
3 C u1 p0 c0 {1,S} {2,S} {4,S}
4 H u0 p0 c0 {3,S}
"""),
    E0 = (-19.6796,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([580,1155,1237,1373,3147,180],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (48.0164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53592,0.00799825,-5.22271e-07,-4.56612e-09,2.08752e-12,-2348.33,9.31392], Tmin=(100,'K'), Tmax=(1079.18,'K')), NASAPolynomial(coeffs=[5.97187,0.0038822,-1.62979e-06,3.36437e-10,-2.54188e-14,-3160.18,-3.94923], Tmin=(1079.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-19.6796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsOJ) + radical(CsF1sHO2s)"""),
)

species(
    label = 'O=CC(F)(F)[CH]F(700)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {6,S}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6 C u1 p0 c0 {3,S} {5,S} {8,S}
7 C u0 p0 c0 {4,D} {5,S} {9,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-576.995,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([249,341,381,651,928,1217,1251,334,575,1197,1424,3202,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.307692,'amu*angstrom^2'), symmetry=1, barrier=(7.07444,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.283425,'amu*angstrom^2'), symmetry=1, barrier=(6.5165,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.32413,0.0699665,-0.000233296,4.73652e-07,-3.72445e-10,-69396.2,12.1051], Tmin=(10,'K'), Tmax=(369.236,'K')), NASAPolynomial(coeffs=[6.12547,0.0284953,-1.96315e-05,6.28355e-09,-7.57912e-13,-69527.3,2.40812], Tmin=(369.236,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-576.995,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""ODCC(F)(F)[CH]F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'O=C(F)C(F)(F)[CH]F(475)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {7,S}
4 F u0 p3 c0 {8,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7 C u1 p0 c0 {3,S} {6,S} {9,S}
8 C u0 p0 c0 {4,S} {5,D} {6,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-831.698,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([249,341,381,651,928,1217,1251,334,575,1197,1424,3202,486,617,768,1157,1926,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.331337,'amu*angstrom^2'), symmetry=1, barrier=(7.6181,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.33233,'amu*angstrom^2'), symmetry=1, barrier=(7.64092,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (129.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.09251,0.0725567,-0.000154812,1.7269e-07,-7.45513e-11,-100025,13.7795], Tmin=(10,'K'), Tmax=(650.241,'K')), NASAPolynomial(coeffs=[9.67075,0.0218463,-1.42006e-05,4.29846e-09,-4.94478e-13,-100664,-13.4601], Tmin=(650.241,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-831.698,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""ODC(F)C(F)(F)[CH]F""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]C=C(F)[CH]F(701)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 O u1 p2 c0 {6,S}
4 C u0 p0 c0 {1,S} {5,S} {6,D}
5 C u1 p0 c0 {2,S} {4,S} {7,S}
6 C u0 p0 c0 {3,S} {4,D} {8,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-255.526,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([271,519,563,612,1379,234,589,736,816,1240,3237,3010,987.5,1337.5,450,1655,534.248],'cm^-1')),
        HinderedRotor(inertia=(0.23444,'amu*angstrom^2'), symmetry=1, barrier=(47.4817,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0441,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28043,0.0363294,-2.79393e-05,9.98028e-09,-1.41201e-12,-30669.9,17.4481], Tmin=(100,'K'), Tmax=(1642.15,'K')), NASAPolynomial(coeffs=[11.8659,0.0129809,-6.61185e-06,1.32195e-09,-9.38741e-14,-33818,-33.5506], Tmin=(1642.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-255.526,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(C=COJ) + radical(Csj(Cd-CdF1s)(F1s)(H))"""),
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
    label = 'O=C(F)[C](F)[CH]F(476)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {7,S}
4 O u0 p2 c0 {7,D}
5 C u1 p0 c0 {1,S} {6,S} {7,S}
6 C u1 p0 c0 {2,S} {5,S} {8,S}
7 C u0 p0 c0 {3,S} {4,D} {5,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-469.624,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([234,347,1316,1464,334,575,1197,1424,3202,611,648,830,1210,1753,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.139167,'amu*angstrom^2'), symmetry=1, barrier=(3.19972,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66661,0.0568091,-9.44748e-05,8.42658e-08,-2.96633e-11,-56404,21.8015], Tmin=(100,'K'), Tmax=(793.064,'K')), NASAPolynomial(coeffs=[7.65033,0.0194496,-1.02342e-05,2.03671e-09,-1.43696e-13,-57127.3,-4.25576], Tmin=(793.064,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-469.624,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCOCsF1s) + radical(Csj(Cs-F1sCOH)(F1s)(H))"""),
)

species(
    label = 'CHF(86)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 C u0 p1 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (138.756,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1169.21,1416.41,2978.06],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (32.017,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2178.39,'J/mol'), sigma=(4.123,'angstroms'), dipoleMoment=(1.8,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.34484,0.00235461,1.93983e-06,-2.65251e-09,7.91169e-13,16766.1,7.05286], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[4.48366,0.00174964,-5.0479e-07,1.08953e-10,-9.87898e-15,16210.2,0.289222], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(138.756,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CHF""", comment="""Thermo library: halogens"""),
)

species(
    label = 'O[C](F)C(F)(F)[CH]F(702)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {7,S} {10,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u1 p0 c0 {3,S} {5,S} {6,S}
8  C u1 p0 c0 {4,S} {6,S} {9,S}
9  H u0 p0 c0 {8,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-711.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.446803,0.0870978,-0.000153262,1.38702e-07,-4.89797e-11,-85414.5,25.6018], Tmin=(100,'K'), Tmax=(797.416,'K')), NASAPolynomial(coeffs=[10.1236,0.0259627,-1.4571e-05,2.94521e-09,-2.08752e-13,-86557.4,-16.3818], Tmin=(797.416,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-711.168,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sF1sCs)(F1s)(O2s-H)) + radical(Csj(Cs-CsF1sF1s)(F1s)(H))"""),
)

species(
    label = '[O][C](F)C(F)(F)CF(703)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  O u1 p2 c0 {8,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {3,S} {6,S} {9,S} {10,S}
8  C u1 p0 c0 {4,S} {5,S} {6,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-683.835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.928365,0.0741272,-0.000109605,8.60433e-08,-2.72597e-11,-82141.8,23.6598], Tmin=(100,'K'), Tmax=(769.517,'K')), NASAPolynomial(coeffs=[10.4163,0.0248071,-1.34644e-05,2.74986e-09,-1.98675e-13,-83602,-19.6279], Tmin=(769.517,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-683.835,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sH)) + radical(Csj(Cs-F1sF1sCs)(F1s)(O2s-H))"""),
)

species(
    label = 'F[CH][C](F)C(F)OF(704)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {6,S}
6  C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
7  C u1 p0 c0 {2,S} {6,S} {8,S}
8  C u1 p0 c0 {3,S} {7,S} {10,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-360.73,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,487,638,688,1119,1325,1387,3149,212,367,445,1450,334,575,1197,1424,3202,180,180,281.704],'cm^-1')),
        HinderedRotor(inertia=(0.577359,'amu*angstrom^2'), symmetry=1, barrier=(13.2746,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106519,'amu*angstrom^2'), symmetry=1, barrier=(2.44907,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.2332,'amu*angstrom^2'), symmetry=1, barrier=(51.3457,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.495464,0.0901781,-0.000172527,1.66312e-07,-6.06334e-11,-43272.3,25.839], Tmin=(100,'K'), Tmax=(836.348,'K')), NASAPolynomial(coeffs=[7.11536,0.0311033,-1.74091e-05,3.47288e-09,-2.42678e-13,-43420.8,0.816789], Tmin=(836.348,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-360.73,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sO2sH)(Cs-F1sHH)(F1s)) + radical(Csj(Cs-F1sCsH)(F1s)(H))"""),
)

species(
    label = '[O]C(F)[C](F)C(F)F(690)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  O u1 p2 c0 {6,S}
6  C u0 p0 c0 {1,S} {5,S} {8,S} {10,S}
7  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
8  C u1 p0 c0 {4,S} {6,S} {7,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-671.097,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([391,562,707,872,1109,1210,1289,3137,522,611,926,1093,1137,1374,1416,3112,212,367,445,1450,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.449187,'amu*angstrom^2'), symmetry=1, barrier=(10.3277,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.623147,'amu*angstrom^2'), symmetry=1, barrier=(14.3274,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.865885,0.0757694,-0.000124873,1.09125e-07,-3.7501e-11,-80608,25.3878], Tmin=(100,'K'), Tmax=(807.061,'K')), NASAPolynomial(coeffs=[9.43139,0.0243063,-1.24778e-05,2.44833e-09,-1.71127e-13,-81697.2,-12.2818], Tmin=(807.061,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-671.097,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sH)) + radical(Csj(Cs-F1sO2sH)(Cs-F1sF1sH)(F1s))"""),
)

species(
    label = 'F[CH]C(F)(F)[CH]OF(705)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {7,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u1 p0 c0 {5,S} {6,S} {9,S}
8  C u1 p0 c0 {3,S} {6,S} {10,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-377.957,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,215,315,519,588,595,1205,1248,3025,407.5,1350,352.5,334,575,1197,1424,3202,241.742,241.903,242.027],'cm^-1')),
        HinderedRotor(inertia=(0.256352,'amu*angstrom^2'), symmetry=1, barrier=(10.6626,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.486073,'amu*angstrom^2'), symmetry=1, barrier=(20.1667,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.485965,'amu*angstrom^2'), symmetry=1, barrier=(20.1698,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0125793,0.0957835,-0.000166503,1.41616e-07,-4.67017e-11,-45320.4,26.3716], Tmin=(100,'K'), Tmax=(800.589,'K')), NASAPolynomial(coeffs=[14.2033,0.0186248,-1.04485e-05,2.10041e-09,-1.47923e-13,-47400.2,-37.8231], Tmin=(800.589,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-377.957,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCsJO) + radical(Csj(Cs-CsF1sF1s)(F1s)(H))"""),
)

species(
    label = '[O][CH]C(F)(F)C(F)F(706)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {7,S}
5  O u1 p2 c0 {8,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {3,S} {4,S} {6,S} {9,S}
8  C u1 p0 c0 {5,S} {6,S} {10,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-695.109,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.927933,0.0763913,-0.000130951,1.22832e-07,-4.53742e-11,-83500.3,24.4879], Tmin=(100,'K'), Tmax=(789.46,'K')), NASAPolynomial(coeffs=[7.09905,0.0307339,-1.68592e-05,3.39841e-09,-2.41351e-13,-84026.2,-0.985275], Tmin=(789.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-695.109,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCOJ) + radical(Csj(Cs-F1sF1sCs)(O2s-H)(H))"""),
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
    E0 = (-320.019,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-153.216,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (127.403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (89.323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-311.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-256.619,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-19.3998,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-172.251,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-131.144,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-274.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-235.254,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (53.4729,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (111.778,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-128.393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (51.2312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-194.559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-164.226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (97.0298,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-133.806,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (60.6783,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-100.388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(F)C(F)(F)[CH]F(590)'],
    products = ['CHFO(46)', 'CHFCF2(54)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]C(F)C(F)[C](F)F(591)'],
    products = ['[O]C(F)C(F)(F)[CH]F(590)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CHF(40)', '[O]C(F)[C](F)F(386)'],
    products = ['[O]C(F)C(F)(F)[CH]F(590)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(6)', 'F[CH]C(F)(F)[CH]F(494)'],
    products = ['[O]C(F)C(F)(F)[CH]F(590)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3335.46,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_sec_rad;O_birad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]C(F)C(F)(F)[CH]F(590)'],
    products = ['FC1OC(F)C1(F)F(595)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_1H;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]C(F)C(F)(F)[CH]F(590)'],
    products = ['O=C(F)C(F)(F)CF(698)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F(37)', '[O]C(F)C(F)=CF(699)'],
    products = ['[O]C(F)C(F)(F)[CH]F(590)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(43.4274,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O][CH]F(168)', 'CHFCF2(54)'],
    products = ['[O]C(F)C(F)(F)[CH]F(590)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(9.4589e-08,'m^3/(mol*s)'), n=3.53001, Ea=(0.47777,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.4944253016374622, var=1.7828810760479818, Tref=1000.0, N=135, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing_Ext-1R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F(37)', 'O=CC(F)(F)[CH]F(700)'],
    products = ['[O]C(F)C(F)(F)[CH]F(590)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.08261e+20,'m^3/(mol*s)'), n=-5.07836, Ea=(14.5529,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CHFO(46)', 'F[CH][C](F)F(183)'],
    products = ['[O]C(F)C(F)(F)[CH]F(590)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(46.5733,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(5)', 'O=C(F)C(F)(F)[CH]F(475)'],
    products = ['[O]C(F)C(F)(F)[CH]F(590)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(15.9,'m^3/(mol*s)'), n=1.84, Ea=(26.2332,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_N-Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_N-Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O][CH]F(168)', 'F[CH][C](F)F(183)'],
    products = ['[O]C(F)C(F)(F)[CH]F(590)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F2(106)', '[O]C=C(F)[CH]F(701)'],
    products = ['[O]C(F)C(F)(F)[CH]F(590)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(17.7024,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction14',
    reactants = ['HF(38)', 'O=C(F)[C](F)[CH]F(476)'],
    products = ['[O]C(F)C(F)(F)[CH]F(590)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.14111,'m^3/(mol*s)'), n=1.29695, Ea=(263.938,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CHF(86)', '[O]C(F)[C](F)F(386)'],
    products = ['[O]C(F)C(F)(F)[CH]F(590)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(786.723,'m^3/(mol*s)'), n=1.25031, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]C(F)C(F)(F)[CH]F(590)'],
    products = ['O[C](F)C(F)(F)[CH]F(702)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.66008e+10,'s^-1'), n=0.776642, Ea=(125.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;Cs_H_out_noH] + [R2H_S;O_rad_out;Cs_H_out] for rate rule [R2H_S;O_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]C(F)C(F)(F)[CH]F(590)'],
    products = ['[O][C](F)C(F)(F)CF(703)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(91.367,'s^-1'), n=3.04268, Ea=(155.793,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_1H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['F[CH][C](F)C(F)OF(704)'],
    products = ['[O]C(F)C(F)(F)[CH]F(590)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(99.3541,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O]C(F)[C](F)C(F)F(690)'],
    products = ['[O]C(F)C(F)(F)[CH]F(590)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(8.71296e+12,'s^-1'), n=0.369741, Ea=(178.885,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_N-4R!H->C_Ext-2R!H-R_5R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_N-4R!H->C_Ext-2R!H-R_5R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction20',
    reactants = ['F[CH]C(F)(F)[CH]OF(705)'],
    products = ['[O]C(F)C(F)(F)[CH]F(590)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(80.2295,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]C(F)C(F)(F)[CH]F(590)'],
    products = ['[O][CH]C(F)(F)C(F)F(706)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(219.632,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

network(
    label = 'PDepNetwork #271',
    isomers = [
        '[O]C(F)C(F)(F)[CH]F(590)',
    ],
    reactants = [
        ('CHFO(46)', 'CHFCF2(54)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #271',
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

