species(
    label = 'O=C[C](F)OC(F)(F)[CH]F(3697)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {9,S}
8  C u1 p0 c0 {4,S} {5,S} {10,S}
9  C u1 p0 c0 {3,S} {7,S} {11,S}
10 C u0 p0 c0 {6,D} {8,S} {12,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-858.086,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,525,597,667,842,1178,1324,280,501,1494,1531,334,575,1197,1424,3202,2782.5,750,1395,475,1775,1000,283.683,283.894,283.925,284.128],'cm^-1')),
        HinderedRotor(inertia=(0.809694,'amu*angstrom^2'), symmetry=1, barrier=(46.0572,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.807218,'amu*angstrom^2'), symmetry=1, barrier=(46.0683,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.126303,'amu*angstrom^2'), symmetry=1, barrier=(7.2236,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.804963,'amu*angstrom^2'), symmetry=1, barrier=(46.0622,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0114602,0.0960629,-0.000150679,1.27097e-07,-4.31363e-11,-103067,29.8835], Tmin=(100,'K'), Tmax=(744.227,'K')), NASAPolynomial(coeffs=[11.4987,0.0319834,-1.70601e-05,3.40261e-09,-2.41089e-13,-104719,-21.8339], Tmin=(744.227,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-858.086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFHH) + group(Cds-OdCsH) + radical(CsCOF1sO2s) + radical(Csj(Cs-F1sF1sO2s)(F1s)(H))"""),
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
    label = '[O]C(F)(F)[CH]F(474)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {6,S}
4 O u1 p2 c0 {5,S}
5 C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6 C u1 p0 c0 {3,S} {5,S} {7,S}
7 H u0 p0 c0 {6,S}
"""),
    E0 = (-445.565,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,525,597,667,842,1178,1324,334,575,1197,1424,3202,180,395.38],'cm^-1')),
        HinderedRotor(inertia=(0.422093,'amu*angstrom^2'), symmetry=1, barrier=(9.70474,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0239,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3211.42,'J/mol'), sigma=(5.51205,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=501.62 K, Pc=43.51 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.99253,0.0473334,-6.8592e-05,5.03378e-08,-1.45679e-11,-53519.7,18.0995], Tmin=(100,'K'), Tmax=(847.727,'K')), NASAPolynomial(coeffs=[9.47106,0.0120471,-6.15724e-06,1.23976e-09,-8.91188e-14,-54787.7,-16.745], Tmin=(847.727,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-445.565,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sF1s)) + radical(Csj(Cs-F1sF1sO2s)(F1s)(H))"""),
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
    label = 'O=C[C](F)O[C](F)F(2641)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {8,S}
3 F u0 p3 c0 {8,S}
4 O u0 p2 c0 {6,S} {8,S}
5 O u0 p2 c0 {7,D}
6 C u1 p0 c0 {1,S} {4,S} {7,S}
7 C u0 p0 c0 {5,D} {6,S} {9,S}
8 C u1 p0 c0 {2,S} {3,S} {4,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-626.807,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([280,501,1494,1531,2782.5,750,1395,475,1775,1000,493,600,700,1144,1293,180,180,2179.8],'cm^-1')),
        HinderedRotor(inertia=(2.38293,'amu*angstrom^2'), symmetry=1, barrier=(54.7882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177256,'amu*angstrom^2'), symmetry=1, barrier=(4.07545,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.38224,'amu*angstrom^2'), symmetry=1, barrier=(54.7724,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45319,0.0645948,-0.000112688,1.07893e-07,-4.00386e-11,-75304.1,24.5605], Tmin=(100,'K'), Tmax=(819.346,'K')), NASAPolynomial(coeffs=[5.28345,0.0286641,-1.5362e-05,3.03579e-09,-2.12433e-13,-75353.4,10.3746], Tmin=(819.346,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-626.807,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsFFHO) + group(Cds-OdCsH) + radical(CsCOF1sO2s) + radical(Csj(F1s)(F1s)(O2s-Cs))"""),
)

species(
    label = 'O=CC1(F)OC(F)(F)C1F(3705)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {10,S}
8  C u0 p0 c0 {2,S} {7,S} {9,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
10 C u0 p0 c0 {6,D} {7,S} {12,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-1103.65,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.813295,0.0659006,-5.87059e-05,2.53468e-08,-4.35136e-12,-132620,26.3534], Tmin=(100,'K'), Tmax=(1387.25,'K')), NASAPolynomial(coeffs=[16.4054,0.0209418,-1.00926e-05,1.98466e-09,-1.41167e-13,-136946,-53.9733], Tmin=(1387.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1103.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCFO) + group(CsCsCsFH) + group(CsCFFO) + group(Cds-OdCsH) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + ring(Cs-Cs(F)(F)-O2s-Cs)"""),
)

species(
    label = 'O=C=C(F)OC(F)(F)CF(3707)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
9  C u0 p0 c0 {4,S} {5,S} {10,D}
10 C u0 p0 c0 {6,D} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-1030.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0326317,0.0886408,-0.000109226,6.53784e-08,-1.53218e-11,-123786,27.6251], Tmin=(100,'K'), Tmax=(1043.54,'K')), NASAPolynomial(coeffs=[18.0256,0.0196714,-1.00877e-05,2.04317e-09,-1.48528e-13,-127541,-59.9471], Tmin=(1043.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1030.39,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + missing(O2d-Cdd) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + group(Cd(Cdd-Od)FO) + missing(Cdd-CdO2d)"""),
)

species(
    label = 'F[CH]C(F)(F)OC1(F)[CH]O1(3833)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {7,S} {9,S}
7  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
8  C u0 p0 c0 {2,S} {3,S} {5,S} {10,S}
9  C u1 p0 c0 {6,S} {7,S} {11,S}
10 C u1 p0 c0 {4,S} {8,S} {12,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-712.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.167548,0.0891923,-0.000114397,7.36415e-08,-1.88532e-11,-85501.7,28.3936], Tmin=(100,'K'), Tmax=(950.782,'K')), NASAPolynomial(coeffs=[15.4163,0.0250406,-1.31896e-05,2.67765e-09,-1.94078e-13,-88401.4,-44.4033], Tmin=(950.782,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-712.015,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(CsCFOO) + group(Cs-CsOsHH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + ring(Cs(F)(O2)-O2s-Cs) + radical(Csj(Cs-F1sO2sO2s)(O2s-Cs)(H)_ring) + radical(Csj(Cs-F1sF1sO2s)(F1s)(H))"""),
)

species(
    label = 'F[C]1[CH]OC(F)C(F)(F)O1(3773)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {8,S} {10,S}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
9  C u1 p0 c0 {5,S} {10,S} {12,S}
10 C u1 p0 c0 {4,S} {6,S} {9,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-852.203,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.569354,0.0886992,-0.000108854,6.61266e-08,-1.47519e-11,-102322,25.0515], Tmin=(100,'K'), Tmax=(1317.18,'K')), NASAPolynomial(coeffs=[18.3551,0.0117311,1.00203e-06,-7.13332e-10,6.73346e-14,-105615,-65.0409], Tmin=(1317.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-852.203,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(CsCFHO) + group(CsCFFO) + group(Cs-CsOsHH) + group(CsCFHO) + ring(1,4-Dioxane) + radical(CCsJOCs) + radical(CsCsF1sO2s) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = '[O]C1[C](F)OC(F)(F)C1F(3817)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {9,S} {10,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {11,S}
8  C u0 p0 c0 {6,S} {7,S} {10,S} {12,S}
9  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
10 C u1 p0 c0 {4,S} {5,S} {8,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-820.922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0890399,0.0757171,-8.06823e-05,4.33337e-08,-8.92938e-12,-98584.5,25.2877], Tmin=(100,'K'), Tmax=(1289.62,'K')), NASAPolynomial(coeffs=[18.1318,0.0146477,-3.71085e-06,4.72944e-10,-2.53586e-14,-102814,-64.7008], Tmin=(1289.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-820.922,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCFHO) + group(CsCFFO) + ring(Tetrahydrofuran) + radical(CC(C)OJ) + radical(CsCsF1sO2s) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F))"""),
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
    label = 'O=C[C](F)OC(F)=CF(3834)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {7,S}
5  O u0 p2 c0 {9,D}
6  C u1 p0 c0 {2,S} {4,S} {9,S}
7  C u0 p0 c0 {1,S} {4,S} {8,D}
8  C u0 p0 c0 {3,S} {7,D} {10,S}
9  C u0 p0 c0 {5,D} {6,S} {11,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-693.596,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([280,501,1494,1531,326,540,652,719,1357,194,682,905,1196,1383,3221,2782.5,750,1395,475,1775,1000,189.799,189.819,954.607],'cm^-1')),
        HinderedRotor(inertia=(2.2032,'amu*angstrom^2'), symmetry=1, barrier=(56.3212,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.20284,'amu*angstrom^2'), symmetry=1, barrier=(56.3214,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.803376,'amu*angstrom^2'), symmetry=1, barrier=(20.5446,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (139.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.96309,0.0715946,-8.63769e-05,5.46166e-08,-1.40398e-11,-83314.9,24.3318], Tmin=(100,'K'), Tmax=(937.759,'K')), NASAPolynomial(coeffs=[11.6709,0.02592,-1.33172e-05,2.67699e-09,-1.92993e-13,-85323.2,-26.6392], Tmin=(937.759,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-693.596,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-OdCsH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(CsCOF1sO2s)"""),
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
    label = 'O=C=C(F)OC(F)(F)[CH]F(3835)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8  C u1 p0 c0 {3,S} {7,S} {11,S}
9  C u0 p0 c0 {4,S} {5,S} {10,D}
10 C u0 p0 c0 {6,D} {9,D}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-829.731,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,525,597,667,842,1178,1324,334,575,1197,1424,3202,197,221,431,657,2120,512.5,787.5,187.465,517.491,517.586,517.729,517.793],'cm^-1')),
        HinderedRotor(inertia=(0.217987,'amu*angstrom^2'), symmetry=1, barrier=(41.4349,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.879882,'amu*angstrom^2'), symmetry=1, barrier=(20.2302,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.879961,'amu*angstrom^2'), symmetry=1, barrier=(20.232,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (157.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.150178,0.0947414,-0.000131879,8.80541e-08,-2.28363e-11,-99647,28.4491], Tmin=(100,'K'), Tmax=(949.14,'K')), NASAPolynomial(coeffs=[18.2736,0.017098,-9.17405e-06,1.86836e-09,-1.35518e-13,-103144,-59.4737], Tmin=(949.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-829.731,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + missing(O2d-Cdd) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + group(Cd(Cdd-Od)FO) + missing(Cdd-CdO2d) + radical(Csj(Cs-F1sF1sO2s)(F1s)(H))"""),
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
    label = 'O=[C]C(F)OC(F)(F)[CH]F(3836)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {9,S}
8  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
9  C u1 p0 c0 {4,S} {7,S} {12,S}
10 C u1 p0 c0 {6,D} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-838.77,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,525,597,667,842,1178,1324,355,410,600,1181,1341,1420,3056,334,575,1197,1424,3202,1855,455,950,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.46827,'amu*angstrom^2'), symmetry=1, barrier=(10.7664,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.469711,'amu*angstrom^2'), symmetry=1, barrier=(10.7996,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.91701,'amu*angstrom^2'), symmetry=1, barrier=(44.0758,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.925273,'amu*angstrom^2'), symmetry=1, barrier=(21.2739,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.668372,0.111931,-0.000192073,1.64367e-07,-5.45084e-11,-100722,31.6976], Tmin=(100,'K'), Tmax=(826.704,'K')), NASAPolynomial(coeffs=[14.5443,0.0260356,-1.39228e-05,2.72313e-09,-1.88571e-13,-102817,-36.2597], Tmin=(826.704,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-838.77,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFHH) + group(Cds-OdCsH) + radical(Csj(Cs-F1sF1sO2s)(F1s)(H)) + radical(COj(Cs-F1sO2sH)(O2d))"""),
)

species(
    label = 'O=[C][C](F)OC(F)(F)CF(3837)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8  C u0 p0 c0 {3,S} {7,S} {11,S} {12,S}
9  C u1 p0 c0 {4,S} {5,S} {10,S}
10 C u1 p0 c0 {6,D} {9,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-900.204,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.17258,0.101425,-0.000171922,1.514e-07,-5.19513e-11,-108129,29.9705], Tmin=(100,'K'), Tmax=(824.547,'K')), NASAPolynomial(coeffs=[11.1824,0.0311121,-1.6306e-05,3.18218e-09,-2.20745e-13,-109483,-19.4809], Tmin=(824.547,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-900.204,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFHH) + group(Cds-OdCsH) + radical(CsCOF1sO2s) + radical(COj(Cs-F1sO2sH)(O2d))"""),
)

species(
    label = 'O=CC(F)(F)O[C](F)[CH]F(3838)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {9,S}
8  C u1 p0 c0 {3,S} {5,S} {10,S}
9  C u0 p0 c0 {6,D} {7,S} {11,S}
10 C u1 p0 c0 {4,S} {8,S} {12,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-797.785,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([251,367,519,700,855,1175,1303,395,473,707,1436,2782.5,750,1395,475,1775,1000,334,575,1197,1424,3202,220.78,220.78,220.78,220.78],'cm^-1')),
        HinderedRotor(inertia=(0.26594,'amu*angstrom^2'), symmetry=1, barrier=(9.19878,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.265939,'amu*angstrom^2'), symmetry=1, barrier=(9.19878,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.26594,'amu*angstrom^2'), symmetry=1, barrier=(9.19878,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17187,'amu*angstrom^2'), symmetry=1, barrier=(40.5347,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.655011,0.113973,-0.000204031,1.81116e-07,-6.16026e-11,-95794.9,32.4389], Tmin=(100,'K'), Tmax=(843.832,'K')), NASAPolynomial(coeffs=[12.9676,0.0286532,-1.54887e-05,3.02427e-09,-2.08403e-13,-97355.3,-26.5928], Tmin=(843.832,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-797.785,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-OdCsH) + radical(Csj(Cs-F1sHH)(F1s)(O2s-Cs)) + radical(Csj(Cs-F1sO2sH)(F1s)(H))"""),
)

species(
    label = 'O=C[C](F)O[C](F)C(F)F(3743)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {8,S} {9,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {11,S}
8  C u1 p0 c0 {3,S} {5,S} {7,S}
9  C u1 p0 c0 {4,S} {5,S} {10,S}
10 C u0 p0 c0 {6,D} {9,S} {12,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-842.103,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([522,611,926,1093,1137,1374,1416,3112,395,473,707,1436,280,501,1494,1531,2782.5,750,1395,475,1775,1000,180,180,180,1758.2],'cm^-1')),
        HinderedRotor(inertia=(1.9197,'amu*angstrom^2'), symmetry=1, barrier=(44.1377,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.91903,'amu*angstrom^2'), symmetry=1, barrier=(44.1223,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.338939,'amu*angstrom^2'), symmetry=1, barrier=(7.79288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.91976,'amu*angstrom^2'), symmetry=1, barrier=(44.1391,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26922,0.073472,-5.96711e-05,-5.12477e-08,8.99447e-11,-101197,28.2543], Tmin=(100,'K'), Tmax=(464.257,'K')), NASAPolynomial(coeffs=[8.3166,0.0369608,-1.99209e-05,3.98967e-09,-2.83096e-13,-102112,-3.14741], Tmin=(464.257,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-842.103,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCFHO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-OdCsH) + radical(Csj(Cs-F1sF1sH)(F1s)(O2s-Cs)) + radical(CsCOF1sO2s)"""),
)

species(
    label = '[O]C(C(=O)F)C(F)(F)[CH]F(3694)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {9,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {10,D}
7  C u0 p0 c0 {5,S} {8,S} {10,S} {11,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {9,S}
9  C u1 p0 c0 {4,S} {8,S} {12,S}
10 C u0 p0 c0 {3,S} {6,D} {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-821.193,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,215,315,519,588,595,1205,1248,334,575,1197,1424,3202,486,617,768,1157,1926,292.442,292.526,292.939,2043.72],'cm^-1')),
        HinderedRotor(inertia=(0.0019618,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128037,'amu*angstrom^2'), symmetry=1, barrier=(7.77456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.175243,'amu*angstrom^2'), symmetry=1, barrier=(10.6533,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0284937,0.0966005,-0.000158527,1.35867e-07,-4.59936e-11,-98629.3,34.8096], Tmin=(100,'K'), Tmax=(790.963,'K')), NASAPolynomial(coeffs=[11.9716,0.0284296,-1.50516e-05,2.97451e-09,-2.08728e-13,-100293,-18.7899], Tmin=(790.963,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-821.193,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + group(COCsFO) + radical(C=OCOJ) + radical(Csj(Cs-CsF1sF1s)(F1s)(H))"""),
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
    label = '[CH]=C(F)OC(F)(F)[CH]F(3839)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {6,S} {8,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
7  C u1 p0 c0 {3,S} {6,S} {10,S}
8  C u0 p0 c0 {4,S} {5,S} {9,D}
9  C u1 p0 c0 {8,D} {11,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-536.962,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,525,597,667,842,1178,1324,334,575,1197,1424,3202,293,496,537,1218,3120,650,792.5,1650,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.14084,'amu*angstrom^2'), symmetry=1, barrier=(26.2302,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14056,'amu*angstrom^2'), symmetry=1, barrier=(26.2236,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14033,'amu*angstrom^2'), symmetry=1, barrier=(26.2184,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0642501,0.0912874,-0.000123532,7.96713e-08,-1.98883e-11,-64436.7,26.7509], Tmin=(100,'K'), Tmax=(986.954,'K')), NASAPolynomial(coeffs=[18.6955,0.0152574,-7.98043e-06,1.61946e-09,-1.17625e-13,-68139.8,-63.5078], Tmin=(986.954,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-536.962,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + group(CdCFO) + group(Cds-CdsHH) + radical(Csj(Cs-F1sF1sO2s)(F1s)(H)) + radical(Cdj(Cd-F1sO2s)(H))"""),
)

species(
    label = 'FC1=COC(F)C(F)(F)O1(3784)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {8,S} {10,S}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
9  C u0 p0 c0 {5,S} {10,D} {12,S}
10 C u0 p0 c0 {4,S} {6,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-1093.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.673729,0.0765628,-7.63139e-05,3.63484e-08,-6.37989e-12,-131271,22.1005], Tmin=(100,'K'), Tmax=(1661.41,'K')), NASAPolynomial(coeffs=[20.8532,0.00832769,1.05063e-07,-2.59821e-10,2.31336e-14,-136159,-85.868], Tmin=(1661.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1093.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(CsCFHO) + group(CsCFFO) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cd(F)=CdOs) + ring(23dihydro14dioxin)"""),
)

species(
    label = '[O][CH]C1(F)OC(F)(F)C1F(3840)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u1 p2 c0 {10,S}
7  C u0 p0 c0 {2,S} {5,S} {8,S} {10,S}
8  C u0 p0 c0 {1,S} {7,S} {9,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
10 C u1 p0 c0 {6,S} {7,S} {12,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-776.465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (158.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.41935,0.0859866,-0.000125082,1.00736e-07,-3.31504e-11,-93265,27.977], Tmin=(100,'K'), Tmax=(739.171,'K')), NASAPolynomial(coeffs=[10.2453,0.0328177,-1.71945e-05,3.43828e-09,-2.45163e-13,-94717.7,-16.4587], Tmin=(739.171,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-776.465,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCCFO) + group(CsCsCsFH) + group(CsCFFO) + group(Cs-CsOsHH) + ring(Cs-Cs(F)(F)-O2s-Cs) + radical(CCOJ) + radical(CCsJOH) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F))"""),
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
    label = 'O=C=[C]OC(F)(F)[CH]F(3841)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
7  C u1 p0 c0 {3,S} {6,S} {10,S}
8  C u1 p0 c0 {4,S} {9,D}
9  C u0 p0 c0 {5,D} {8,D}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-401.78,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,525,597,667,842,1178,1324,334,575,1197,1424,3202,1685,370,2120,512.5,787.5,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.03476,'amu*angstrom^2'), symmetry=1, barrier=(23.7911,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03426,'amu*angstrom^2'), symmetry=1, barrier=(23.7797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03145,'amu*angstrom^2'), symmetry=1, barrier=(23.7151,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (138.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0430687,0.0881736,-0.000127002,8.59142e-08,-2.22505e-11,-48181.4,28.6325], Tmin=(100,'K'), Tmax=(957.24,'K')), NASAPolynomial(coeffs=[18.6645,0.0103589,-5.063e-06,9.88246e-10,-7.00422e-14,-51746.4,-60.3908], Tmin=(957.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-401.78,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + missing(O2d-Cdd) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(Csj(Cs-F1sF1sO2s)(F1s)(H)) + radical(C=CJO)"""),
)

species(
    label = 'O[C]=C(F)OC(F)(F)[CH]F(3842)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {10,S} {12,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8  C u1 p0 c0 {3,S} {7,S} {11,S}
9  C u0 p0 c0 {4,S} {5,S} {10,D}
10 C u1 p0 c0 {6,S} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-741.668,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,253,525,597,667,842,1178,1324,334,575,1197,1424,3202,293,496,537,1218,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.613085,0.104538,-0.000148135,1.00009e-07,-2.60828e-11,-89038.5,32.2054], Tmin=(100,'K'), Tmax=(946.517,'K')), NASAPolynomial(coeffs=[20.2204,0.0164942,-8.60429e-06,1.73081e-09,-1.24662e-13,-92982.3,-67.1586], Tmin=(946.517,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-741.668,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + radical(Csj(Cs-F1sF1sO2s)(F1s)(H)) + radical(C=CJO)"""),
)

species(
    label = 'FC=C(F)O[C](F)[CH]OF(3843)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {6,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {4,S} {9,S}
7  C u1 p0 c0 {2,S} {5,S} {9,S}
8  C u0 p0 c0 {1,S} {5,S} {10,D}
9  C u1 p0 c0 {6,S} {7,S} {11,S}
10 C u0 p0 c0 {3,S} {8,D} {12,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-410.336,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,395,473,707,1436,326,540,652,719,1357,3025,407.5,1350,352.5,194,682,905,1196,1383,3221,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.848005,0.114456,-0.000186717,1.49264e-07,-4.65311e-11,-49184.7,32.3236], Tmin=(100,'K'), Tmax=(791.572,'K')), NASAPolynomial(coeffs=[17.3486,0.0225013,-1.24591e-05,2.49726e-09,-1.76623e-13,-52065.4,-51.2107], Tmin=(791.572,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-410.336,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2sCF) + group(CsCFHO) + group(Cs-CsOsHH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(CsCsF1sO2s) + radical(CCsJO)"""),
)

species(
    label = 'F[CH]C(F)(F)O[C]=COF(3844)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {6,S}
5  O u0 p2 c0 {7,S} {10,S}
6  O u0 p2 c0 {4,S} {9,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8  C u1 p0 c0 {3,S} {7,S} {11,S}
9  C u0 p0 c0 {6,S} {10,D} {12,S}
10 C u1 p0 c0 {5,S} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-425.45,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,253,525,597,667,842,1178,1324,334,575,1197,1424,3202,3010,987.5,1337.5,450,1655,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.903718,0.104671,-0.000141217,8.84495e-08,-2.10674e-11,-50990.3,33.6889], Tmin=(100,'K'), Tmax=(1043.28,'K')), NASAPolynomial(coeffs=[23.9907,0.00922401,-3.98603e-06,7.57833e-10,-5.3955e-14,-56184.7,-87.4671], Tmin=(1043.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-425.45,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2sCF) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(Csj(Cs-F1sF1sO2s)(F1s)(H)) + radical(C=CJO)"""),
)

species(
    label = '[O]C=[C]OC(F)(F)C(F)F(3845)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {7,S} {10,S}
6  O u1 p2 c0 {9,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8  C u0 p0 c0 {3,S} {4,S} {7,S} {11,S}
9  C u0 p0 c0 {6,S} {10,D} {12,S}
10 C u1 p0 c0 {5,S} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-851.979,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([223,363,546,575,694,1179,1410,235,523,627,1123,1142,1372,1406,3097,3010,987.5,1337.5,450,1655,1685,370,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.1592,'amu*angstrom^2'), symmetry=1, barrier=(26.6522,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15845,'amu*angstrom^2'), symmetry=1, barrier=(26.6351,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15753,'amu*angstrom^2'), symmetry=1, barrier=(26.6139,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (158.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.02429,0.0957926,-0.000114568,6.24362e-08,-1.27324e-11,-102276,34.0847], Tmin=(100,'K'), Tmax=(1295.03,'K')), NASAPolynomial(coeffs=[27.3087,0.00248766,2.12779e-07,-1.0527e-10,7.61107e-15,-109129,-108.056], Tmin=(1295.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-851.979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)2) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJO)"""),
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
    E0 = (-355.495,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (79.7306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (90.7113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-347.211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-330.522,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-172.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-285.316,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-293.538,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-238.624,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-69.2876,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-223.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-115.336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (2.86013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (92.6791,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (14.5392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-199.232,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-328.613,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-97.5245,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-137.611,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-175.374,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (208.663,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-347.964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-273.874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (102.685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-76.2982,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (137.016,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (145.274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-153.431,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=C[C](F)OC(F)(F)[CH]F(3697)'],
    products = ['CHFCF2(54)', 'O=CC(=O)F(557)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O=C[C]F-2(2705)', '[O]C(F)(F)[CH]F(474)'],
    products = ['O=C[C](F)OC(F)(F)[CH]F(3697)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['CHF(40)', 'O=C[C](F)O[C](F)F(2641)'],
    products = ['O=C[C](F)OC(F)(F)[CH]F(3697)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['O=C[C](F)OC(F)(F)[CH]F(3697)'],
    products = ['O=CC1(F)OC(F)(F)C1F(3705)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_noH;Cpri_rad_out_1H]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O=C[C](F)OC(F)(F)[CH]F(3697)'],
    products = ['O=C=C(F)OC(F)(F)CF(3707)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O=C[C](F)OC(F)(F)[CH]F(3697)'],
    products = ['F[CH]C(F)(F)OC1(F)[CH]O1(3833)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(9.27213e+10,'s^-1'), n=0.543712, Ea=(183.455,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra_cs] for rate rule [R3_CO;carbonyl_intra_H;radadd_intra_cs]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O=C[C](F)OC(F)(F)[CH]F(3697)'],
    products = ['F[C]1[CH]OC(F)C(F)(F)O1(3773)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(8.58001e+07,'s^-1'), n=0.730566, Ea=(70.179,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra_cs] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_cs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['O=C[C](F)OC(F)(F)[CH]F(3697)'],
    products = ['[O]C1[C](F)OC(F)(F)C1F(3817)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.54917e+06,'s^-1'), n=1.13913, Ea=(61.9573,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;multiplebond_intra;radadd_intra_cs] for rate rule [R6;carbonylbond_intra_H;radadd_intra_cs]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F[CH][C](F)F(183)', 'O=CC(=O)F(557)'],
    products = ['O=C[C](F)OC(F)(F)[CH]F(3697)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(9.07578e-06,'m^3/(mol*s)'), n=3.04336, Ea=(27.1552,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3757377757886876, var=2.242054186761003, Tref=1000.0, N=1042, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F(37)', 'O=C[C](F)OC(F)=CF(3834)'],
    products = ['O=C[C](F)OC(F)(F)[CH]F(3697)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(48.8258,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CHFCF2(54)', '[O][CH]C(=O)F(438)'],
    products = ['O=C[C](F)OC(F)(F)[CH]F(3697)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.0603e-20,'m^3/(mol*s)'), n=7.37188, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.4161614592809219, var=2.5944177025948685, Tref=1000.0, N=250, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R_Ext-1R!H-R_N-6R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R_Ext-1R!H-R_N-6R!H-inRing"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(5)', 'O=C=C(F)OC(F)(F)[CH]F(3835)'],
    products = ['O=C[C](F)OC(F)(F)[CH]F(3697)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.09572,'m^3/(mol*s)'), n=2.03221, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.053734944950677085, var=3.1017726857793777, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_N-4R!H->C',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_N-4R!H->C"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F[CH][C](F)F(183)', '[O][CH]C(=O)F(438)'],
    products = ['O=C[C](F)OC(F)(F)[CH]F(3697)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['O=C[C]F(628)', '[O]C(F)(F)[CH]F(474)'],
    products = ['O=C[C](F)OC(F)(F)[CH]F(3697)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(128827,'m^3/(mol*s)'), n=0.469398, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.009898522871762284, var=0.3372166703721302, Tref=1000.0, N=6, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CHF(86)', 'O=C[C](F)O[C](F)F(2641)'],
    products = ['O=C[C](F)OC(F)(F)[CH]F(3697)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(786.723,'m^3/(mol*s)'), n=1.25031, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl"""),
)

reaction(
    label = 'reaction16',
    reactants = ['O=[C]C(F)OC(F)(F)[CH]F(3836)'],
    products = ['O=C[C](F)OC(F)(F)[CH]F(3697)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(6.23078e+06,'s^-1'), n=1.82328, Ea=(136.948,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_noH] for rate rule [R2H_S;CO_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O=C[C](F)OC(F)(F)[CH]F(3697)'],
    products = ['O=[C][C](F)OC(F)(F)CF(3837)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.0378492,'s^-1'), n=3.26, Ea=(26.8822,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_1H;XH_out] for rate rule [R5HJ_3;C_rad_out_1H;CO_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['O=CC(F)(F)O[C](F)[CH]F(3838)'],
    products = ['O=C[C](F)OC(F)(F)[CH]F(3697)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(197.67,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction19',
    reactants = ['O=C[C](F)O[C](F)C(F)F(3743)'],
    products = ['O=C[C](F)OC(F)(F)[CH]F(3697)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(8.5608e+11,'s^-1'), n=0.253035, Ea=(201.901,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_N-4R!H->C_Ext-2R!H-R_N-5R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_N-4R!H->C_Ext-2R!H-R_N-5R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O]C(C(=O)F)C(F)(F)[CH]F(3694)'],
    products = ['O=C[C](F)OC(F)(F)[CH]F(3697)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(143.228,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction21',
    reactants = ['O(6)', '[CH]=C(F)OC(F)(F)[CH]F(3839)'],
    products = ['O=C[C](F)OC(F)(F)[CH]F(3697)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['O=C[C](F)OC(F)(F)[CH]F(3697)'],
    products = ['FC1=COC(F)C(F)(F)O1(3784)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;C_rad_out_1H;Ypri_rad_out] for rate rule [R6_SSSDS;C_rad_out_1H;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['O=C[C](F)OC(F)(F)[CH]F(3697)'],
    products = ['[O][CH]C1(F)OC(F)(F)C1F(3840)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.06771e+06,'s^-1'), n=1.35044, Ea=(81.621,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 80.3 to 81.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction24',
    reactants = ['HF(38)', 'O=C=[C]OC(F)(F)[CH]F(3841)'],
    products = ['O=C[C](F)OC(F)(F)[CH]F(3697)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(282.988,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction25',
    reactants = ['O[C]=C(F)OC(F)(F)[CH]F(3842)'],
    products = ['O=C[C](F)OC(F)(F)[CH]F(3697)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['FC=C(F)O[C](F)[CH]OF(3843)'],
    products = ['O=C[C](F)OC(F)(F)[CH]F(3697)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(44.7612,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction27',
    reactants = ['F[CH]C(F)(F)O[C]=COF(3844)'],
    products = ['O=C[C](F)OC(F)(F)[CH]F(3697)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(68.1331,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[O]C=[C]OC(F)(F)C(F)F(3845)'],
    products = ['O=C[C](F)OC(F)(F)[CH]F(3697)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(9.39758e+07,'s^-1'), n=1.42748, Ea=(195.957,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F
Multiplied by reaction path degeneracy 2.0"""),
)

network(
    label = 'PDepNetwork #1322',
    isomers = [
        'O=C[C](F)OC(F)(F)[CH]F(3697)',
    ],
    reactants = [
        ('CHFCF2(54)', 'O=CC(=O)F(557)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1322',
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

