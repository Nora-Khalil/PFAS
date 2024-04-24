species(
    label = '[O]C(F)C(F)(F)[C](F)F(596)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {3,S} {6,S} {7,S} {10,S}
9  C u1 p0 c0 {4,S} {5,S} {7,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-878.656,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,391,562,707,872,1109,1210,1289,3137,190,488,555,1236,1407,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.422139,'amu*angstrom^2'), symmetry=1, barrier=(9.70581,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.423151,'amu*angstrom^2'), symmetry=1, barrier=(9.72907,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (148.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.112869,0.097307,-0.000182447,1.67383e-07,-5.81171e-11,-105549,27.2421], Tmin=(100,'K'), Tmax=(851.553,'K')), NASAPolynomial(coeffs=[10.4133,0.025247,-1.38098e-05,2.71007e-09,-1.86769e-13,-106445,-15.756], Tmin=(851.553,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-878.656,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sH)) + radical(Csj(Cs-CsF1sF1s)(F1s)(F1s))"""),
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
    label = 'FC1OC(F)(F)C1(F)F(358)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  O u0 p2 c0 {8,S} {9,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {3,S} {6,S} {7,S} {10,S}
9  C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-1180.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (148.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.71256,0.0257385,7.45682e-05,-1.86832e-07,1.16264e-10,-141941,13.1956], Tmin=(10,'K'), Tmax=(589.954,'K')), NASAPolynomial(coeffs=[6.25413,0.0368921,-2.59639e-05,8.33009e-09,-9.99094e-13,-142735,-1.91138], Tmin=(589.954,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1180.2,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), label="""FC1OC(F)(F)C1(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C(F)C(F)(F)C(F)F(750)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {9,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {3,S} {4,S} {7,S} {10,S}
9  C u0 p0 c0 {5,S} {6,D} {7,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-1236.58,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (148.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.46542,0.0583387,-7.64623e-05,5.64267e-08,-1.77795e-11,-148727,13.233], Tmin=(10,'K'), Tmax=(723.416,'K')), NASAPolynomial(coeffs=[8.07976,0.0328244,-2.35587e-05,7.67318e-09,-9.31136e-13,-149395,-7.53457], Tmin=(723.416,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1236.58,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""ODC(F)C(F)(F)C(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]C(F)C(F)=C(F)F(686)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {8,S}
4 F u0 p3 c0 {8,S}
5 O u1 p2 c0 {6,S}
6 C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
7 C u0 p0 c0 {2,S} {6,S} {8,D}
8 C u0 p0 c0 {3,S} {4,S} {7,D}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-688.092,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([361,769,965,1078,1132,1246,3247,323,467,575,827,1418,182,240,577,636,1210,1413,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.490926,'amu*angstrom^2'), symmetry=1, barrier=(11.2874,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (129.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.67943,0.0236114,0.000161177,-5.32832e-07,4.64216e-10,-82755.4,13.1547], Tmin=(10,'K'), Tmax=(426.943,'K')), NASAPolynomial(coeffs=[7.98521,0.0295583,-2.23418e-05,7.66806e-09,-9.73587e-13,-83544.9,-8.89409], Tmin=(426.943,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-688.092,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""[O]C(F)C(F)DC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'O=CC(F)(F)[C](F)F(562)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {7,S}
4 F u0 p3 c0 {7,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7 C u1 p0 c0 {3,S} {4,S} {6,S}
8 C u0 p0 c0 {5,D} {6,S} {9,S}
9 H u0 p0 c0 {8,S}
"""),
    E0 = (-781.834,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([249,341,381,651,928,1217,1251,190,488,555,1236,1407,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.458696,'amu*angstrom^2'), symmetry=1, barrier=(10.5463,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.458927,'amu*angstrom^2'), symmetry=1, barrier=(10.5516,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (129.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.46705,0.0472046,-1.13744e-05,-1.21811e-07,1.37968e-10,-94032.2,14.3018], Tmin=(10,'K'), Tmax=(439.005,'K')), NASAPolynomial(coeffs=[8.05744,0.0275797,-2.01745e-05,6.74516e-09,-8.39508e-13,-94649.2,-6.50166], Tmin=(439.005,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-781.834,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""ODCC(F)(F)[C](F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'O=C(F)C(F)(F)[C](F)F(520)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {8,S}
4 F u0 p3 c0 {8,S}
5 F u0 p3 c0 {9,S}
6 O u0 p2 c0 {9,D}
7 C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8 C u1 p0 c0 {3,S} {4,S} {7,S}
9 C u0 p0 c0 {5,S} {6,D} {7,S}
"""),
    E0 = (-1032.75,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([249,341,381,651,928,1217,1251,190,488,555,1236,1407,486,617,768,1157,1926,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.332332,'amu*angstrom^2'), symmetry=1, barrier=(7.64097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.333042,'amu*angstrom^2'), symmetry=1, barrier=(7.65728,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (147.023,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.2563,0.0709694,-0.000139664,1.49914e-07,-6.49112e-11,-124211,14.8198], Tmin=(10,'K'), Tmax=(559.52,'K')), NASAPolynomial(coeffs=[9.54024,0.0260455,-1.92285e-05,6.41639e-09,-7.94562e-13,-124914,-11.8477], Tmin=(559.52,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1032.75,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""ODC(F)C(F)(F)[C](F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]C=C(F)[C](F)F(689)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {6,S}
4 O u1 p2 c0 {7,S}
5 C u0 p0 c0 {1,S} {6,S} {7,D}
6 C u1 p0 c0 {2,S} {3,S} {5,S}
7 C u0 p0 c0 {4,S} {5,D} {8,S}
8 H u0 p0 c0 {7,S}
"""),
    E0 = (-438.735,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([271,519,563,612,1379,161,297,490,584,780,1358,3010,987.5,1337.5,450,1655,291.478],'cm^-1')),
        HinderedRotor(inertia=(0.00194419,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06308,0.0417084,-3.67637e-05,1.54384e-08,-2.58165e-12,-52697.5,21.5788], Tmin=(100,'K'), Tmax=(1411.13,'K')), NASAPolynomial(coeffs=[11.8789,0.0138844,-7.18735e-06,1.46542e-09,-1.06151e-13,-55467.7,-29.157], Tmin=(1411.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-438.735,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(C=COJ) + radical(Csj(Cd-F1sCd)(F1s)(F1s))"""),
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
    label = 'O=C(F)[C](F)[C](F)F(465)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {7,S}
4 F u0 p3 c0 {8,S}
5 O u0 p2 c0 {8,D}
6 C u1 p0 c0 {1,S} {7,S} {8,S}
7 C u1 p0 c0 {2,S} {3,S} {6,S}
8 C u0 p0 c0 {4,S} {5,D} {6,S}
"""),
    E0 = (-679.265,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([234,347,1316,1464,190,488,555,1236,1407,611,648,830,1210,1753,213.146,213.181],'cm^-1')),
        HinderedRotor(inertia=(0.310447,'amu*angstrom^2'), symmetry=1, barrier=(10.0173,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.44559,'amu*angstrom^2'), symmetry=1, barrier=(46.6913,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (128.025,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28617,0.0654735,-0.000109958,9.59754e-08,-3.31717e-11,-81604.5,24.1632], Tmin=(100,'K'), Tmax=(762.412,'K')), NASAPolynomial(coeffs=[9.40916,0.0188768,-1.04531e-05,2.12042e-09,-1.51252e-13,-82727.5,-12.0636], Tmin=(762.412,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-679.265,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCOCsF1s) + radical(Csj(Cs-F1sCOH)(F1s)(F1s))"""),
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
    label = 'O[C](F)C(F)(F)[C](F)F(751)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  O u0 p2 c0 {8,S} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u1 p0 c0 {3,S} {6,S} {7,S}
9  C u1 p0 c0 {4,S} {5,S} {7,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-911.399,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (148.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.246233,0.114742,-0.00024608,2.41853e-07,-8.65585e-11,-109483,28.5437], Tmin=(100,'K'), Tmax=(879.667,'K')), NASAPolynomial(coeffs=[6.97424,0.0313541,-1.76824e-05,3.42696e-09,-2.31103e-13,-108798,5.75267], Tmin=(879.667,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-911.399,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sF1sCs)(F1s)(O2s-H)) + radical(Csj(Cs-CsF1sF1s)(F1s)(F1s))"""),
)

species(
    label = '[O][C](F)C(F)(F)C(F)F(479)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {9,S}
6  O u1 p2 c0 {9,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {3,S} {4,S} {7,S} {10,S}
9  C u1 p0 c0 {5,S} {6,S} {7,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-887.497,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (148.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0974746,0.104473,-0.000217372,2.12569e-07,-7.61301e-11,-106618,27.6412], Tmin=(100,'K'), Tmax=(876.931,'K')), NASAPolynomial(coeffs=[6.55909,0.0311779,-1.70432e-05,3.28878e-09,-2.22206e-13,-106067,6.92335], Tmin=(876.931,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-887.497,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sH)) + radical(Csj(Cs-F1sF1sCs)(F1s)(O2s-H))"""),
)

species(
    label = 'FOC(F)[C](F)[C](F)F(752)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {7,S}
7  C u0 p0 c0 {1,S} {6,S} {8,S} {10,S}
8  C u1 p0 c0 {2,S} {7,S} {9,S}
9  C u1 p0 c0 {3,S} {4,S} {8,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-554.238,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,487,638,688,1119,1325,1387,3149,212,367,445,1450,190,488,555,1236,1407,180,180,1429.91],'cm^-1')),
        HinderedRotor(inertia=(0.343221,'amu*angstrom^2'), symmetry=1, barrier=(7.89132,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.344694,'amu*angstrom^2'), symmetry=1, barrier=(7.92519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.20156,'amu*angstrom^2'), symmetry=1, barrier=(50.6183,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (148.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.374325,0.094179,-0.00018371,1.79044e-07,-6.56836e-11,-66542.9,29.0937], Tmin=(100,'K'), Tmax=(838.028,'K')), NASAPolynomial(coeffs=[6.75477,0.0328339,-1.86157e-05,3.72302e-09,-2.60254e-13,-66527.6,5.91095], Tmin=(838.028,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-554.238,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sO2sH)(Cs-F1sF1sH)(F1s)) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
)

species(
    label = '[O]C(F)[C](F)C(F)(F)F(753)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {9,S}
6  O u1 p2 c0 {7,S}
7  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
8  C u0 p0 c0 {1,S} {3,S} {4,S} {9,S}
9  C u1 p0 c0 {5,S} {7,S} {8,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-905.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (148.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.282175,0.0911879,-0.000162421,1.46335e-07,-5.10419e-11,-108731,26.3982], Tmin=(100,'K'), Tmax=(813.654,'K')), NASAPolynomial(coeffs=[10.7357,0.0252575,-1.40707e-05,2.82167e-09,-1.98567e-13,-109951,-18.9207], Tmin=(813.654,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-905.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sH)) + radical(Csj(Cs-F1sO2sH)(Cs-F1sF1sF1s)(F1s))"""),
)

species(
    label = 'FO[CH]C(F)(F)[C](F)F(754)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {8,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u1 p0 c0 {6,S} {7,S} {10,S}
9  C u1 p0 c0 {3,S} {4,S} {7,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-582.669,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,215,315,519,588,595,1205,1248,3025,407.5,1350,352.5,190,488,555,1236,1407,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.978346,'amu*angstrom^2'), symmetry=1, barrier=(22.4941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.976445,'amu*angstrom^2'), symmetry=1, barrier=(22.4504,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.978419,'amu*angstrom^2'), symmetry=1, barrier=(22.4958,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (148.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.329329,0.102221,-0.000169483,1.34165e-07,-4.10935e-11,-69929.5,27.734], Tmin=(100,'K'), Tmax=(807.206,'K')), NASAPolynomial(coeffs=[17.0601,0.0160518,-9.3617e-06,1.92485e-09,-1.38086e-13,-72736.9,-52.4359], Tmin=(807.206,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-582.669,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCsJO) + radical(Csj(Cs-CsF1sF1s)(F1s)(F1s))"""),
)

species(
    label = '[O][CH]C(F)(F)C(F)(F)F(755)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {8,S}
6  O u1 p2 c0 {9,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {3,S} {4,S} {5,S} {7,S}
9  C u1 p0 c0 {6,S} {7,S} {10,S}
10 H u0 p0 c0 {9,S}
"""),
    E0 = (-926.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (148.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.878782,0.0808758,-0.000149802,1.48704e-07,-5.64864e-11,-111292,26.4715], Tmin=(100,'K'), Tmax=(815,'K')), NASAPolynomial(coeffs=[4.80511,0.0359989,-2.00777e-05,4.03873e-09,-2.85206e-13,-111081,13.5498], Tmin=(815,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-926.168,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCOJ) + radical(Csj(Cs-F1sF1sCs)(O2s-H)(H))"""),
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
    E0 = (-396.983,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (69.4696,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (8.61423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-388.699,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-333.583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-91.0381,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-226.541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-211.304,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-354.897,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-311.836,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-31.5185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (49.9812,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-208.575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-167.97,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-271.523,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-274.626,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (25.5853,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-228.017,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-18.9174,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-193.475,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(F)C(F)(F)[C](F)F(596)'],
    products = ['CHFO(46)', 'CF2CF2(60)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CF2(42)', '[O]C(F)[C](F)F(386)'],
    products = ['[O]C(F)C(F)(F)[C](F)F(596)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['O(6)', 'F[CH]C(F)(F)[C](F)F(473)'],
    products = ['[O]C(F)C(F)(F)[C](F)F(596)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_sec_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]C(F)C(F)(F)[C](F)F(596)'],
    products = ['FC1OC(F)(F)C1(F)F(358)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_noH;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]C(F)C(F)(F)[C](F)F(596)'],
    products = ['O=C(F)C(F)(F)C(F)F(750)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['F(37)', '[O]C(F)C(F)=C(F)F(686)'],
    products = ['[O]C(F)C(F)(F)[C](F)F(596)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(42.4891,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O][CH]F(168)', 'CF2CF2(60)'],
    products = ['[O]C(F)C(F)(F)[C](F)F(596)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.89178e-07,'m^3/(mol*s)'), n=3.53001, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.4944253016374622, var=1.7828810760479818, Tref=1000.0, N=135, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing_Ext-1R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing_Ext-1R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F(37)', 'O=CC(F)(F)[C](F)F(562)'],
    products = ['[O]C(F)C(F)(F)[C](F)F(596)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.08261e+20,'m^3/(mol*s)'), n=-5.07836, Ea=(15.9653,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CHFO(46)', 'F[C](F)[C](F)F(230)'],
    products = ['[O]C(F)C(F)(F)[C](F)F(596)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(8.08706e-06,'m^3/(mol*s)'), n=3.0961, Ea=(51.2365,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(5)', 'O=C(F)C(F)(F)[C](F)F(520)'],
    products = ['[O]C(F)C(F)(F)[C](F)F(596)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(15.9,'m^3/(mol*s)'), n=1.84, Ea=(27.4411,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_N-Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_N-Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O][CH]F(168)', 'F[C](F)[C](F)F(230)'],
    products = ['[O]C(F)C(F)(F)[C](F)F(596)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.808e+07,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F2(106)', '[O]C=C(F)[C](F)F(689)'],
    products = ['[O]C(F)C(F)(F)[C](F)F(596)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(15.8483,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction13',
    reactants = ['HF(38)', 'O=C(F)[C](F)[C](F)F(465)'],
    products = ['[O]C(F)C(F)(F)[C](F)F(596)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.14111,'m^3/(mol*s)'), n=1.29695, Ea=(270.13,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['CF2(87)', '[O]C(F)[C](F)F(386)'],
    products = ['[O]C(F)C(F)(F)[C](F)F(596)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(786.723,'m^3/(mol*s)'), n=1.25031, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O]C(F)C(F)(F)[C](F)F(596)'],
    products = ['O[C](F)C(F)(F)[C](F)F(751)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.66008e+10,'s^-1'), n=0.776642, Ea=(125.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;Cs_H_out_noH] + [R2H_S;O_rad_out;Cs_H_out] for rate rule [R2H_S;O_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]C(F)C(F)(F)[C](F)F(596)'],
    products = ['[O][C](F)C(F)(F)C(F)F(479)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.50974e+07,'s^-1'), n=1.33047, Ea=(122.357,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_noH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['FOC(F)[C](F)[C](F)F(752)'],
    products = ['[O]C(F)C(F)(F)[C](F)F(596)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(98.1506,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]C(F)C(F)(F)[C](F)F(596)'],
    products = ['[O]C(F)[C](F)C(F)(F)F(753)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.01526e+12,'s^-1'), n=0.18834, Ea=(168.965,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction19',
    reactants = ['FO[CH]C(F)(F)[C](F)F(754)'],
    products = ['[O]C(F)C(F)(F)[C](F)F(596)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(82.079,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O]C(F)C(F)(F)[C](F)F(596)'],
    products = ['[O][CH]C(F)(F)C(F)(F)F(755)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(203.508,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

network(
    label = 'PDepNetwork #277',
    isomers = [
        '[O]C(F)C(F)(F)[C](F)F(596)',
    ],
    reactants = [
        ('CHFO(46)', 'CF2CF2(60)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #277',
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

