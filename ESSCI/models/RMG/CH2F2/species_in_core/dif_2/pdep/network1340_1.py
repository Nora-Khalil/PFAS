species(
    label = 'O=C[C](F)OC[C](F)F(3715)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {6,S} {7,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {4,S} {8,S} {10,S} {11,S}
7  C u1 p0 c0 {3,S} {4,S} {9,S}
8  C u1 p0 c0 {1,S} {2,S} {6,S}
9  C u0 p0 c0 {5,D} {7,S} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-637.259,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,280,501,1494,1531,190,488,555,1236,1407,2782.5,750,1395,475,1775,1000,250.873,250.873,250.875,886.505],'cm^-1')),
        HinderedRotor(inertia=(0.119164,'amu*angstrom^2'), symmetry=1, barrier=(5.32208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.119167,'amu*angstrom^2'), symmetry=1, barrier=(5.32208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12394,'amu*angstrom^2'), symmetry=1, barrier=(50.1977,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12395,'amu*angstrom^2'), symmetry=1, barrier=(50.1977,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.465973,0.0871166,-0.000143504,1.30944e-07,-4.71509e-11,-76526.4,28.7803], Tmin=(100,'K'), Tmax=(807.79,'K')), NASAPolynomial(coeffs=[7.65037,0.035717,-1.86751e-05,3.67257e-09,-2.57005e-13,-77170.8,-1.15107], Tmin=(807.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-637.259,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsOsHH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFFH) + group(Cds-OdCsH) + radical(CsCOF1sO2s) + radical(Csj(Cs-O2sHH)(F1s)(F1s))"""),
)

species(
    label = 'O=CC(=O)F(557)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,D}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {2,D} {5,S} {6,S}
5 C u0 p0 c0 {1,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-483.116,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,286,619,818,1246,1924],'cm^-1')),
        HinderedRotor(inertia=(0.753127,'amu*angstrom^2'), symmetry=1, barrier=(17.3159,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (76.0265,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3426.83,'J/mol'), sigma=(5.15159,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=535.26 K, Pc=56.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88796,0.00794544,8.17071e-05,-2.34537e-07,1.90378e-10,-58102.5,9.57665], Tmin=(10,'K'), Tmax=(438.71,'K')), NASAPolynomial(coeffs=[4.97623,0.0166174,-1.15194e-05,3.74172e-09,-4.59031e-13,-58377,3.18363], Tmin=(438.71,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-483.116,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""ODCC(DO)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CH2CF2(56)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {4,D} {5,S} {6,S}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-361.616,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(64.0125,'amu')),
        NonlinearRotor(inertia=([45.7027,48.2614,93.9642],'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([437.293,557.015,653.832,726.079,816.319,956,966.438,1345.56,1413.22,1792.31,3202.97,3303.55],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2091.09,'J/mol'), sigma=(4.442,'angstroms'), dipoleMoment=(1.4,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.10281,-0.0101072,0.000121983,-2.28108e-07,1.37933e-10,-43490.6,7.77929], Tmin=(10,'K'), Tmax=(534.293,'K')), NASAPolynomial(coeffs=[2.52167,0.0198841,-1.31824e-05,4.13929e-09,-4.93215e-13,-43580.8,11.9914], Tmin=(534.293,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-361.616,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""CDC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C[C]F-2(2705)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {2,D} {4,S} {5,S}
4 C u2 p0 c0 {1,S} {3,S}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (22.7053,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,163,1167],'cm^-1')),
        HinderedRotor(inertia=(0.0337629,'amu*angstrom^2'), symmetry=1, barrier=(13.6228,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.0271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.34131,0.0165482,-1.72266e-05,1.26793e-08,-4.54539e-12,2752.64,10.5202], Tmin=(100,'K'), Tmax=(638.12,'K')), NASAPolynomial(coeffs=[4.07864,0.0119262,-6.36172e-06,1.32812e-09,-9.81566e-14,2658.54,7.29432], Tmin=(638.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(22.7053,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCCl_triplet)"""),
)

species(
    label = '[O]C[C](F)F(728)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u1 p2 c0 {4,S}
4 C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
5 C u1 p0 c0 {1,S} {2,S} {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-233.591,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,190,488,555,1236,1407,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0983613,'amu*angstrom^2'), symmetry=1, barrier=(2.26152,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3214.57,'J/mol'), sigma=(5.31812,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=502.11 K, Pc=48.49 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.3035,0.0417665,-6.88725e-05,6.30488e-08,-2.27983e-11,-28037.7,16.0779], Tmin=(100,'K'), Tmax=(796.61,'K')), NASAPolynomial(coeffs=[5.91256,0.0167111,-8.63847e-06,1.7145e-09,-1.20926e-13,-28392.7,0.867685], Tmin=(796.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-233.591,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCOJ) + radical(Csj(Cs-O2sHH)(F1s)(F1s))"""),
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
    label = '[CH2]O[C](F)C=O(3929)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {4,S} {6,S}
3 O u0 p2 c0 {5,D}
4 C u1 p0 c0 {1,S} {2,S} {5,S}
5 C u0 p0 c0 {3,D} {4,S} {7,S}
6 C u1 p0 c0 {2,S} {8,S} {9,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-203.967,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([280,501,1494,1531,2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,991.468,1447.98],'cm^-1')),
        HinderedRotor(inertia=(2.39779,'amu*angstrom^2'), symmetry=1, barrier=(55.1298,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.345491,'amu*angstrom^2'), symmetry=1, barrier=(7.94353,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.40057,'amu*angstrom^2'), symmetry=1, barrier=(55.1938,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (90.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.4615,0.0387779,-3.31098e-05,5.68874e-09,8.9303e-12,-24480.9,18.0671], Tmin=(100,'K'), Tmax=(533.911,'K')), NASAPolynomial(coeffs=[5.17615,0.02491,-1.2326e-05,2.43414e-09,-1.73509e-13,-24863,5.81043], Tmin=(533.911,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-203.967,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-OsHHH) + group(Cds-OdCsH) + radical(CsCOF1sO2s) + radical(CsJOCH3)"""),
)

species(
    label = 'O=CC1(F)OCC1(F)F(3717)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {9,S}
7  C u0 p0 c0 {2,S} {3,S} {6,S} {8,S}
8  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
9  C u0 p0 c0 {5,D} {6,S} {12,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-887.214,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44094,0.058535,-4.99851e-05,2.16982e-08,-3.908e-12,-106617,23.807], Tmin=(100,'K'), Tmax=(1280.42,'K')), NASAPolynomial(coeffs=[11.4362,0.0273098,-1.3405e-05,2.65221e-09,-1.89287e-13,-109177,-26.8852], Tmin=(1280.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-887.214,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCFO) + group(CsCsCsFF) + group(Cs-CsOsHH) + group(Cds-OdCsH) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + ring(O2s-Cs-Cs-Cs(F))"""),
)

species(
    label = 'O=CC(F)OC=C(F)F(3930)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {6,S} {7,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {4,S} {8,S} {10,S}
7  C u0 p0 c0 {4,S} {9,D} {11,S}
8  C u0 p0 c0 {5,D} {6,S} {12,S}
9  C u0 p0 c0 {2,S} {3,S} {7,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-838.142,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.306834,0.0757829,-7.93372e-05,3.99764e-08,-7.84901e-12,-100668,26.8139], Tmin=(100,'K'), Tmax=(1242.77,'K')), NASAPolynomial(coeffs=[18.7748,0.0163413,-7.592e-06,1.48953e-09,-1.06821e-13,-105258,-66.2972], Tmin=(1242.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-838.142,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs)"""),
)

species(
    label = 'O=C=C(F)OCC(F)F(3720)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {12,S}
8  C u0 p0 c0 {3,S} {4,S} {9,D}
9  C u0 p0 c0 {5,D} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-811.372,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.464644,0.0810234,-0.000100419,6.29263e-08,-1.56675e-11,-97460.8,26.0198], Tmin=(100,'K'), Tmax=(977.916,'K')), NASAPolynomial(coeffs=[14.6433,0.0230277,-1.14604e-05,2.28134e-09,-1.63856e-13,-100234,-42.0675], Tmin=(977.916,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-811.372,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + missing(O2d-Cdd) + group(Cs-CsOsHH) + group(CsCsFFH) + group(Cd(Cdd-Od)FO) + missing(Cdd-CdO2d)"""),
)

species(
    label = 'F[C](F)COC1(F)[CH]O1(3931)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {6,S} {7,S}
5  O u0 p2 c0 {6,S} {8,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
7  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
8  C u1 p0 c0 {5,S} {6,S} {12,S}
9  C u1 p0 c0 {2,S} {3,S} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-491.189,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.640135,0.080201,-0.000106615,7.61885e-08,-2.21651e-11,-58961,27.3146], Tmin=(100,'K'), Tmax=(833.175,'K')), NASAPolynomial(coeffs=[11.2162,0.0294266,-1.52048e-05,3.04645e-09,-2.18471e-13,-60723.3,-21.7785], Tmin=(833.175,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-491.189,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(CsCFOO) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(CsCsFFH) + ring(Cs(F)(O2)-O2s-Cs) + radical(Csj(Cs-F1sO2sO2s)(O2s-Cs)(H)_ring) + radical(Csj(Cs-O2sHH)(F1s)(F1s))"""),
)

species(
    label = 'F[C]1[CH]OC(F)(F)CO1(3932)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {7,S} {9,S}
6  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
8  C u1 p0 c0 {3,S} {4,S} {9,S}
9  C u1 p0 c0 {5,S} {8,S} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-650.099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.102683,0.0814856,-0.000100001,6.24827e-08,-1.43807e-11,-78033.4,22.885], Tmin=(100,'K'), Tmax=(1284.82,'K')), NASAPolynomial(coeffs=[15.4967,0.0144443,-1.6128e-07,-5.14439e-10,5.50483e-14,-80516.9,-50.3481], Tmin=(1284.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-650.099,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(Cs-CsOsHH) + group(CsCFFO) + group(CsCFHO) + group(Cs-CsOsHH) + ring(1,4-Dioxane) + radical(CsCsF1sO2s) + radical(CCsJOCs)"""),
)

species(
    label = '[O]C1[C](F)OCC1(F)F(3933)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {8,S} {9,S}
5  O u1 p2 c0 {7,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
8  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
9  C u1 p0 c0 {3,S} {4,S} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-616.071,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.771054,0.0670394,-6.50951e-05,2.99817e-08,-4.17504e-12,-73976.3,21.9015], Tmin=(100,'K'), Tmax=(909.352,'K')), NASAPolynomial(coeffs=[14.0978,0.0199849,-6.55693e-06,1.05357e-09,-6.757e-14,-76878.3,-43.7557], Tmin=(909.352,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-616.071,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCsCsFF) + group(CsCFHO) + group(Cs-CsOsHH) + ring(Tetrahydrofuran) + radical(CC(C)OJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[CH2][C](F)F(185)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 C u1 p0 c0 {4,S} {5,S} {6,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-109.048,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,190,488,555,1236,1407],'cm^-1')),
        HinderedRotor(inertia=(0.00258864,'amu*angstrom^2'), symmetry=1, barrier=(7.63529,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.25819,0.0186259,-1.56413e-05,8.223e-09,-2.04577e-12,-13090.7,13.5552], Tmin=(100,'K'), Tmax=(887.666,'K')), NASAPolynomial(coeffs=[4.4701,0.0131647,-6.4128e-06,1.29206e-09,-9.37476e-14,-13305.9,7.85287], Tmin=(887.666,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-109.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Cs_P) + radical(CsCsF1sF1s)"""),
)

species(
    label = '[O][CH]C(=O)F(438)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u1 p2 c0 {4,S}
3 O u0 p2 c0 {5,D}
4 C u1 p0 c0 {2,S} {5,S} {6,S}
5 C u0 p0 c0 {1,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-214.477,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,611,648,830,1210,1753,380.326,381.637],'cm^-1')),
        HinderedRotor(inertia=(0.483418,'amu*angstrom^2'), symmetry=1, barrier=(49.7921,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (76.0265,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.80089,0.0290438,-3.0283e-05,1.66613e-08,-3.81625e-12,-25754.7,13.8243], Tmin=(100,'K'), Tmax=(1028.25,'K')), NASAPolynomial(coeffs=[6.95398,0.0128879,-6.71484e-06,1.38085e-09,-1.01081e-13,-26608.7,-6.32768], Tmin=(1028.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-214.477,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(C=OCOJ) + radical(OCJC=O)"""),
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
    label = 'O=C[C](F)OC=C(F)F(3740)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {7,S}
5  O u0 p2 c0 {9,D}
6  C u1 p0 c0 {1,S} {4,S} {9,S}
7  C u0 p0 c0 {4,S} {8,D} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {7,D}
9  C u0 p0 c0 {5,D} {6,S} {11,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-698.916,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([280,501,1494,1531,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,2782.5,750,1395,475,1775,1000,252.504,252.664,252.73],'cm^-1')),
        HinderedRotor(inertia=(0.855143,'amu*angstrom^2'), symmetry=1, barrier=(38.7734,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.860405,'amu*angstrom^2'), symmetry=1, barrier=(38.7698,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.860148,'amu*angstrom^2'), symmetry=1, barrier=(38.7567,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (139.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.840523,0.0688158,-7.30507e-05,3.84531e-08,-8.03782e-12,-83945.8,25.1163], Tmin=(100,'K'), Tmax=(1156.16,'K')), NASAPolynomial(coeffs=[14.9572,0.0199766,-9.68736e-06,1.91687e-09,-1.37558e-13,-87210.1,-45.0369], Tmin=(1156.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-698.916,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(CsCOF1sO2s)"""),
)

species(
    label = 'O=C=C(F)OC[C](F)F(3934)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
7  C u1 p0 c0 {1,S} {2,S} {6,S}
8  C u0 p0 c0 {3,S} {4,S} {9,D}
9  C u0 p0 c0 {5,D} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-608.905,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,190,488,555,1236,1407,197,221,431,657,2120,512.5,787.5,431.247,431.247,431.307,431.33,431.368],'cm^-1')),
        HinderedRotor(inertia=(0.011321,'amu*angstrom^2'), symmetry=1, barrier=(1.49426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0890994,'amu*angstrom^2'), symmetry=1, barrier=(11.7642,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.231284,'amu*angstrom^2'), symmetry=1, barrier=(30.539,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (139.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.361292,0.0852977,-0.000122535,8.85695e-08,-2.52592e-11,-73107.9,27.2302], Tmin=(100,'K'), Tmax=(859.948,'K')), NASAPolynomial(coeffs=[14.0883,0.0214479,-1.11634e-05,2.23014e-09,-1.59269e-13,-75468.8,-36.9234], Tmin=(859.948,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-608.905,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + missing(O2d-Cdd) + group(Cs-CsOsHH) + group(CsCsFFH) + group(Cd(Cdd-Od)FO) + missing(Cdd-CdO2d) + radical(Csj(Cs-O2sHH)(F1s)(F1s))"""),
)

species(
    label = 'O=C[C]F(628)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {2,D} {4,S} {5,S}
4 C u0 p1 c0 {1,S} {3,S}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (35.6539,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,262,1290],'cm^-1')),
        HinderedRotor(inertia=(0.407026,'amu*angstrom^2'), symmetry=1, barrier=(9.35834,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (60.0271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.79853,0.0176142,-2.7066e-05,3.1192e-08,-1.61007e-11,4288.74,9.0865], Tmin=(10,'K'), Tmax=(518.444,'K')), NASAPolynomial(coeffs=[4.38817,0.0119962,-7.71993e-06,2.33921e-09,-2.7033e-13,4241.97,6.76768], Tmin=(518.444,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(35.6539,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""ODC[C]F""", comment="""Thermo library: CHOF_G4"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(897.963,'J/mol'), sigma=(3.977,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28591,0.0107608,-1.05382e-05,4.89881e-09,-8.86384e-13,-24340.7,13.1348], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[5.33121,0.00197748,-9.60248e-07,2.10704e-10,-1.5954e-14,-25190.9,-2.56367], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-203.712,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CF2""", comment="""Thermo library: halogens"""),
)

species(
    label = 'O=C[C](F)O[CH]C(F)F(3935)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
7  C u1 p0 c0 {4,S} {6,S} {11,S}
8  C u1 p0 c0 {3,S} {4,S} {9,S}
9  C u0 p0 c0 {5,D} {8,S} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-653.705,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.362751,0.0889366,-0.000145595,1.30401e-07,-4.6109e-11,-78500.2,28.9282], Tmin=(100,'K'), Tmax=(809.421,'K')), NASAPolynomial(coeffs=[8.61739,0.0340183,-1.76443e-05,3.45625e-09,-2.41276e-13,-79373.7,-6.29198], Tmin=(809.421,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-653.705,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsOsHH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFFH) + group(Cds-OdCsH) + radical(Csj(Cs-F1sF1sH)(O2s-Cs)(H)) + radical(CsCOF1sO2s)"""),
)

species(
    label = 'O=[C]C(F)OC[C](F)F(3936)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {7,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {4,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {4,S} {9,S} {12,S}
8  C u1 p0 c0 {2,S} {3,S} {6,S}
9  C u1 p0 c0 {5,D} {7,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-617.944,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,355,410,600,1181,1341,1420,3056,190,488,555,1236,1407,1855,455,950,245.302,245.31,245.326,245.328],'cm^-1')),
        HinderedRotor(inertia=(0.169142,'amu*angstrom^2'), symmetry=1, barrier=(7.22218,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169158,'amu*angstrom^2'), symmetry=1, barrier=(7.22197,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.607843,'amu*angstrom^2'), symmetry=1, barrier=(25.9544,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169167,'amu*angstrom^2'), symmetry=1, barrier=(7.22278,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.133505,0.102242,-0.000181991,1.64014e-07,-5.65557e-11,-74183.5,30.3926], Tmin=(100,'K'), Tmax=(849.596,'K')), NASAPolynomial(coeffs=[10.5092,0.0301072,-1.57415e-05,3.04278e-09,-2.08709e-13,-75196.9,-14.5387], Tmin=(849.596,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-617.944,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsOsHH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFFH) + group(Cds-OdCsH) + radical(Csj(Cs-O2sHH)(F1s)(F1s)) + radical(COj(Cs-F1sO2sH)(O2d))"""),
)

species(
    label = 'O=CC(F)O[CH][C](F)F(3937)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {6,S} {7,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {4,S} {8,S} {10,S}
7  C u1 p0 c0 {4,S} {9,S} {11,S}
8  C u0 p0 c0 {5,D} {6,S} {12,S}
9  C u1 p0 c0 {2,S} {3,S} {7,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-590.464,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([232,360,932,1127,1349,1365,3045,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,190,488,555,1236,1407,268.543,268.545,268.556,2381.83],'cm^-1')),
        HinderedRotor(inertia=(0.223096,'amu*angstrom^2'), symmetry=1, barrier=(11.4177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223156,'amu*angstrom^2'), symmetry=1, barrier=(11.4175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.622535,'amu*angstrom^2'), symmetry=1, barrier=(31.8624,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.62247,'amu*angstrom^2'), symmetry=1, barrier=(31.8625,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00181795,0.0969439,-0.000162796,1.43161e-07,-4.93871e-11,-70881.1,30.8446], Tmin=(100,'K'), Tmax=(808.742,'K')), NASAPolynomial(coeffs=[10.9524,0.0303161,-1.60972e-05,3.17262e-09,-2.21819e-13,-72244.6,-17.1407], Tmin=(808.742,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-590.464,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsOsHH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFFH) + group(Cds-OdCsH) + radical(Csj(Cs-F1sF1sH)(O2s-Cs)(H)) + radical(Csj(Cs-O2sHH)(F1s)(F1s))"""),
)

species(
    label = 'O=[C][C](F)OCC(F)F(3938)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {12,S}
8  C u1 p0 c0 {3,S} {4,S} {9,S}
9  C u1 p0 c0 {5,D} {8,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-681.185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.229778,0.0942036,-0.00016466,1.51052e-07,-5.31725e-11,-81802.7,28.4681], Tmin=(100,'K'), Tmax=(848.101,'K')), NASAPolynomial(coeffs=[8.17147,0.0338143,-1.72916e-05,3.32715e-09,-2.2823e-13,-82325,-3.67522], Tmin=(848.101,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-681.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsOsHH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFFH) + group(Cds-OdCsH) + radical(CsCOF1sO2s) + radical(COj(Cs-F1sO2sH)(O2d))"""),
)

species(
    label = '[O]C(C[C](F)F)C(=O)F(3712)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u1 p2 c0 {6,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {4,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {6,S} {8,S} {11,S} {12,S}
8  C u1 p0 c0 {2,S} {3,S} {7,S}
9  C u0 p0 c0 {1,S} {5,D} {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-617.045,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,190,488,555,1236,1407,486,617,768,1157,1926,356.176,356.212,1686.23,4000],'cm^-1')),
        HinderedRotor(inertia=(0.115527,'amu*angstrom^2'), symmetry=1, barrier=(10.4011,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115522,'amu*angstrom^2'), symmetry=1, barrier=(10.4012,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.319391,'amu*angstrom^2'), symmetry=1, barrier=(28.7534,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.916191,0.0736061,-9.63505e-05,6.92066e-08,-2.03687e-11,-74107.4,32.3273], Tmin=(100,'K'), Tmax=(822.146,'K')), NASAPolynomial(coeffs=[10.1294,0.0287836,-1.45764e-05,2.90056e-09,-2.07275e-13,-75622.4,-10.3171], Tmin=(822.146,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-617.045,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(CsCsFFH) + group(COCsFO) + radical(C=OCOJ) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
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
    label = '[CH]=C(F)OC[C](F)F(3939)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {5,S} {7,S}
5  C u0 p0 c0 {4,S} {6,S} {9,S} {10,S}
6  C u1 p0 c0 {1,S} {2,S} {5,S}
7  C u0 p0 c0 {3,S} {4,S} {8,D}
8  C u1 p0 c0 {7,D} {11,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-316.136,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,190,488,555,1236,1407,293,496,537,1218,3120,650,792.5,1650,344.641,344.657,344.664,344.671],'cm^-1')),
        HinderedRotor(inertia=(0.263674,'amu*angstrom^2'), symmetry=1, barrier=(22.2296,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.213087,'amu*angstrom^2'), symmetry=1, barrier=(17.9699,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146642,'amu*angstrom^2'), symmetry=1, barrier=(12.3646,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.494684,0.0812539,-0.000111985,7.71242e-08,-2.09093e-11,-37899.7,25.3638], Tmin=(100,'K'), Tmax=(904.299,'K')), NASAPolynomial(coeffs=[14.3808,0.0198311,-1.01003e-05,2.01234e-09,-1.4397e-13,-40411.1,-40.2316], Tmin=(904.299,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-316.136,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(CsCsFFH) + group(CdCFO) + group(Cds-CdsHH) + radical(Csj(Cs-O2sHH)(F1s)(F1s)) + radical(Cdj(Cd-F1sO2s)(H))"""),
)

species(
    label = 'FC1=COC(F)(F)CO1(3940)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {7,S} {9,S}
6  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
8  C u0 p0 c0 {3,S} {4,S} {9,D}
9  C u0 p0 c0 {5,S} {8,D} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-890.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.158996,0.0688455,-6.59666e-05,3.11005e-08,-5.45096e-12,-106985,19.757], Tmin=(100,'K'), Tmax=(1657.38,'K')), NASAPolynomial(coeffs=[18.0788,0.0109685,-1.04052e-06,-6.15344e-11,1.06991e-14,-111127,-71.7007], Tmin=(1657.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-890.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(CsCFFO) + group(CdCFO) + group(Cds-CdsOsH) + longDistanceInteraction_cyclic(Cd(F)=CdOs) + ring(23dihydro14dioxin)"""),
)

species(
    label = 'OC=C(F)OC=C(F)F(3941)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {6,S} {7,S}
5  O u0 p2 c0 {8,S} {12,S}
6  C u0 p0 c0 {1,S} {4,S} {8,D}
7  C u0 p0 c0 {4,S} {9,D} {10,S}
8  C u0 p0 c0 {5,S} {6,D} {11,S}
9  C u0 p0 c0 {2,S} {3,S} {7,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-754.788,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.573552,0.0786514,-9.51427e-05,5.86771e-08,-1.44516e-11,-90659.3,27.765], Tmin=(100,'K'), Tmax=(985.789,'K')), NASAPolynomial(coeffs=[14.0612,0.0239233,-1.18673e-05,2.36009e-09,-1.69451e-13,-93318.5,-37.1121], Tmin=(985.789,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-754.788,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs)"""),
)

species(
    label = '[O][CH]C1(F)OCC1(F)F(3942)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {6,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u1 p2 c0 {9,S}
6  C u0 p0 c0 {3,S} {4,S} {7,S} {9,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
8  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
9  C u1 p0 c0 {5,S} {6,S} {12,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-560.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.575262,0.0843162,-0.000136958,1.24903e-07,-4.51455e-11,-67241.4,27.1147], Tmin=(100,'K'), Tmax=(803.295,'K')), NASAPolynomial(coeffs=[7.3528,0.0357359,-1.85487e-05,3.64835e-09,-2.55644e-13,-67851.8,-1.12], Tmin=(803.295,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-560.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCCFO) + group(CsCsCsFF) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + ring(O2s-Cs-Cs-Cs(F)) + radical(CCOJ) + radical(CCsJOH) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
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
    label = 'O=C=[C]OC[C](F)F(3943)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  O u0 p2 c0 {5,S} {7,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {3,S} {6,S} {9,S} {10,S}
6  C u1 p0 c0 {1,S} {2,S} {5,S}
7  C u1 p0 c0 {3,S} {8,D}
8  C u0 p0 c0 {4,D} {7,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-180.954,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,190,488,555,1236,1407,1685,370,2120,512.5,787.5,301.354,301.911,302.437,302.669],'cm^-1')),
        HinderedRotor(inertia=(0.242736,'amu*angstrom^2'), symmetry=1, barrier=(15.7103,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.240374,'amu*angstrom^2'), symmetry=1, barrier=(15.7083,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.243015,'amu*angstrom^2'), symmetry=1, barrier=(15.711,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.583971,0.0783902,-0.0001165,8.49456e-08,-2.40364e-11,-21643.6,27.3076], Tmin=(100,'K'), Tmax=(872.884,'K')), NASAPolynomial(coeffs=[14.4848,0.0146909,-7.0382e-06,1.34599e-09,-9.34131e-14,-24070.4,-37.8663], Tmin=(872.884,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-180.954,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + missing(O2d-Cdd) + group(Cs-CsOsHH) + group(CsCsFFH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(Csj(Cs-O2sHH)(F1s)(F1s)) + radical(C=CJO)"""),
)

species(
    label = 'O[C]=C(F)OC[C](F)F(3944)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {9,S} {12,S}
6  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
7  C u1 p0 c0 {1,S} {2,S} {6,S}
8  C u0 p0 c0 {3,S} {4,S} {9,D}
9  C u1 p0 c0 {5,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-520.842,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,190,488,555,1236,1407,293,496,537,1218,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0899483,0.0949598,-0.000138332,9.99336e-08,-2.82505e-11,-62499.9,30.9445], Tmin=(100,'K'), Tmax=(870.431,'K')), NASAPolynomial(coeffs=[16.041,0.0208304,-1.05841e-05,2.09005e-09,-1.48182e-13,-65308,-44.6394], Tmin=(870.431,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-520.842,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(CsCsFFH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + radical(Csj(Cs-O2sHH)(F1s)(F1s)) + radical(C=CJO)"""),
)

species(
    label = 'OC=C(F)O[CH][C](F)F(3945)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {6,S} {7,S}
5  O u0 p2 c0 {8,S} {12,S}
6  C u0 p0 c0 {1,S} {4,S} {8,D}
7  C u1 p0 c0 {4,S} {9,S} {10,S}
8  C u0 p0 c0 {5,S} {6,D} {11,S}
9  C u1 p0 c0 {2,S} {3,S} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-566.658,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,326,540,652,719,1357,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,190,488,555,1236,1407,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.926212,'amu*angstrom^2'), symmetry=1, barrier=(21.2954,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.922088,'amu*angstrom^2'), symmetry=1, barrier=(21.2006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.923851,'amu*angstrom^2'), symmetry=1, barrier=(21.2411,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.925463,'amu*angstrom^2'), symmetry=1, barrier=(21.2782,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.817931,0.104506,-0.000145562,9.3883e-08,-2.30324e-11,-67978.3,30.5865], Tmin=(100,'K'), Tmax=(1012.75,'K')), NASAPolynomial(coeffs=[23.3506,0.0090489,-4.18059e-06,8.15388e-10,-5.8532e-14,-72873.7,-86.3189], Tmin=(1012.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-566.658,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(CsCsFFH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + radical(CCsJOC(O)) + radical(Csj(Cs-O2sHH)(F1s)(F1s))"""),
)

species(
    label = 'FOC=[C]OC[C](F)F(3946)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {5,S}
4  O u0 p2 c0 {6,S} {9,S}
5  O u0 p2 c0 {3,S} {8,S}
6  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
7  C u1 p0 c0 {1,S} {2,S} {6,S}
8  C u0 p0 c0 {5,S} {9,D} {12,S}
9  C u1 p0 c0 {4,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-204.623,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,2850,1437.5,1250,1305,750,350,190,488,555,1236,1407,3010,987.5,1337.5,450,1655,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.267695,0.0937368,-0.000126562,8.19135e-08,-2.03957e-11,-24456.6,32.0249], Tmin=(100,'K'), Tmax=(993.622,'K')), NASAPolynomial(coeffs=[19.5315,0.0140312,-6.23507e-06,1.1802e-09,-8.26804e-14,-28391.1,-63.3678], Tmin=(993.622,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-204.623,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2sCF) + group(Cs-CsOsHH) + group(CsCsFFH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(Csj(Cs-O2sHH)(F1s)(F1s)) + radical(C=CJO)"""),
)

species(
    label = '[O]C=[C]OCC(F)(F)F(3947)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {6,S} {9,S}
5  O u1 p2 c0 {8,S}
6  C u0 p0 c0 {4,S} {7,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
8  C u0 p0 c0 {5,S} {9,D} {12,S}
9  C u1 p0 c0 {4,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-658.461,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.612656,0.086274,-9.84942e-05,5.19002e-08,-1.02439e-11,-79015.5,31.1709], Tmin=(100,'K'), Tmax=(1357.25,'K')), NASAPolynomial(coeffs=[24.9026,0.00425627,-3.12172e-07,-2.83182e-11,3.15889e-15,-85313.4,-97.4052], Tmin=(1357.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-658.461,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-CsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJO)"""),
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
    E0 = (-252.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (173.525,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (214.171,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-244.564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-163.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-227.875,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-69.3933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-182.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-190.891,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-170.355,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-180.531,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-101.757,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-12.6891,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (60.8857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (186.474,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-20.9377,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-240.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-96.585,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-48.484,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-189.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-80.9429,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (311.309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-245.317,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-227.875,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-175.619,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (205.331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (26.3484,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-87.6349,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (247.92,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-64.6591,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=C[C](F)OC[C](F)F(3715)'],
    products = ['O=CC(=O)F(557)', 'CH2CF2(56)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O=C[C]F-2(2705)', '[O]C[C](F)F(728)'],
    products = ['O=C[C](F)OC[C](F)F(3715)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['CF2(42)', '[CH2]O[C](F)C=O(3929)'],
    products = ['O=C[C](F)OC[C](F)F(3715)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/O;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['O=C[C](F)OC[C](F)F(3715)'],
    products = ['O=CC1(F)OCC1(F)F(3717)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_noH;Cpri_rad_out_noH]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O=C[C](F)OC[C](F)F(3715)'],
    products = ['O=CC(F)OC=C(F)F(3930)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O=C[C](F)OC[C](F)F(3715)'],
    products = ['O=C=C(F)OCC(F)F(3720)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O=C[C](F)OC[C](F)F(3715)'],
    products = ['F[C](F)COC1(F)[CH]O1(3931)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(9.27213e+10,'s^-1'), n=0.543712, Ea=(183.455,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra_cs] for rate rule [R3_CO;carbonyl_intra_H;radadd_intra_cs]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['O=C[C](F)OC[C](F)F(3715)'],
    products = ['F[C]1[CH]OC(F)(F)CO1(3932)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(8.58001e+07,'s^-1'), n=0.730566, Ea=(70.179,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra_cs] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_cs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['O=C[C](F)OC[C](F)F(3715)'],
    products = ['[O]C1[C](F)OCC1(F)F(3933)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.54917e+06,'s^-1'), n=1.13913, Ea=(61.9573,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;multiplebond_intra;radadd_intra_cs] for rate rule [R6;carbonylbond_intra_H;radadd_intra_cs]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['O=CC(=O)F(557)', '[CH2][C](F)F(185)'],
    products = ['O=C[C](F)OC[C](F)F(3715)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(37.3985,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O][CH]C(=O)F(438)', 'CH2CF2(56)'],
    products = ['O=C[C](F)OC[C](F)F(3715)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.30416e-16,'m^3/(mol*s)'), n=6.08142, Ea=(11.151,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.42494241855262277, var=2.081425725945606, Tref=1000.0, N=380, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(5)', 'O=C[C](F)OC=C(F)F(3740)'],
    products = ['O=C[C](F)OC[C](F)F(3715)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.09572,'m^3/(mol*s)'), n=2.03221, Ea=(0.943906,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.053734944950677085, var=3.1017726857793777, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_N-4R!H->C',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_N-4R!H->C"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(5)', 'O=C=C(F)OC[C](F)F(3934)'],
    products = ['O=C[C](F)OC[C](F)F(3715)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.09572,'m^3/(mol*s)'), n=2.03221, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.053734944950677085, var=3.1017726857793777, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_N-4R!H->C',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_N-4R!H->C"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O][CH]C(=O)F(438)', '[CH2][C](F)F(185)'],
    products = ['O=C[C](F)OC[C](F)F(3715)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O=C[C]F(628)', '[O]C[C](F)F(728)'],
    products = ['O=C[C](F)OC[C](F)F(3715)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(128827,'m^3/(mol*s)'), n=0.469398, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.009898522871762284, var=0.3372166703721302, Tref=1000.0, N=6, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['CF2(87)', '[CH2]O[C](F)C=O(3929)'],
    products = ['O=C[C](F)OC[C](F)F(3715)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(2.3306,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O=C[C](F)OC[C](F)F(3715)'],
    products = ['O=C[C](F)O[CH]C(F)F(3935)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.245e+06,'s^-1'), n=1.223, Ea=(12.1085,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 442 used for R2H_S;C_rad_out_noH;Cs_H_out_H/NonDeO
Exact match found for rate rule [R2H_S;C_rad_out_noH;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['O=[C]C(F)OC[C](F)F(3936)'],
    products = ['O=C[C](F)OC[C](F)F(3715)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(6.23078e+06,'s^-1'), n=1.82328, Ea=(136.948,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_noH] for rate rule [R2H_S;CO_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['O=CC(F)O[CH][C](F)F(3937)'],
    products = ['O=C[C](F)OC[C](F)F(3715)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.22107e+09,'s^-1'), n=1.12, Ea=(157.569,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_O;Y_rad_out;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['O=C[C](F)OC[C](F)F(3715)'],
    products = ['O=[C][C](F)OCC(F)F(3938)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(57774.9,'s^-1'), n=1.76812, Ea=(63.6725,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_single;XH_out] for rate rule [R5HJ_3;C_rad_out_noH;CO_H_out]
Euclidian distance = 1.7320508075688772
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]C(C[C](F)F)C(=O)F(3712)'],
    products = ['O=C[C](F)OC[C](F)F(3715)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(151.691,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction22',
    reactants = ['O(6)', '[CH]=C(F)OC[C](F)F(3939)'],
    products = ['O=C[C](F)OC[C](F)F(3715)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['O=C[C](F)OC[C](F)F(3715)'],
    products = ['FC1=COC(F)(F)CO1(3940)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;C_rad_out_single;Ypri_rad_out] for rate rule [R6_SSSDS;C_rad_out_noH;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['O=C[C](F)OC[C](F)F(3715)'],
    products = ['OC=C(F)OC=C(F)F(3941)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['O=C[C](F)OC[C](F)F(3715)'],
    products = ['[O][CH]C1(F)OCC1(F)F(3942)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.06771e+06,'s^-1'), n=1.35044, Ea=(77.2295,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 76.8 to 77.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction26',
    reactants = ['HF(38)', 'O=C=[C]OC[C](F)F(3943)'],
    products = ['O=C[C](F)OC[C](F)F(3715)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(282.988,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction27',
    reactants = ['O[C]=C(F)OC[C](F)F(3944)'],
    products = ['O=C[C](F)OC[C](F)F(3715)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['OC=C(F)O[CH][C](F)F(3945)'],
    products = ['O=C[C](F)OC[C](F)F(3715)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(722272,'s^-1'), n=1.6737, Ea=(94.6126,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSMS;Y_rad_out;XH_out] for rate rule [R5H_SSMS;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['FOC=[C]OC[C](F)F(3946)'],
    products = ['O=C[C](F)OC[C](F)F(3715)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(68.1331,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction30',
    reactants = ['O=C[C](F)OC[C](F)F(3715)'],
    products = ['[O]C=[C]OCC(F)(F)F(3947)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(188.189,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

network(
    label = 'PDepNetwork #1340',
    isomers = [
        'O=C[C](F)OC[C](F)F(3715)',
    ],
    reactants = [
        ('O=CC(=O)F(557)', 'CH2CF2(56)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1340',
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

