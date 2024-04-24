species(
    label = 'C=C([O])C(F)(F)[CH]F(3306)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {7,S}
4  O u1 p2 c0 {6,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u0 p0 c0 {4,S} {5,S} {8,D}
7  C u1 p0 c0 {3,S} {5,S} {9,S}
8  C u0 p0 c0 {6,D} {10,S} {11,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-467.791,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,350,440,435,1725,334,575,1197,1424,3202,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0899896,'amu*angstrom^2'), symmetry=1, barrier=(2.06904,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.439401,'amu*angstrom^2'), symmetry=1, barrier=(10.1027,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.911817,0.0744209,-0.000105464,7.64294e-08,-2.04545e-11,-56157.1,25.95], Tmin=(100,'K'), Tmax=(651.816,'K')), NASAPolynomial(coeffs=[10.4193,0.0250454,-1.2478e-05,2.43461e-09,-1.70891e-13,-57587,-17.3104], Tmin=(651.816,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-467.791,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Csj(Cs-F1sF1sCd)(F1s)(H))"""),
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
    label = 'CH2CO(75)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,D}
2 C u0 p0 c0 {3,D} {4,S} {5,S}
3 C u0 p0 c0 {1,D} {2,D}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (-60.8183,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2120,512.5,787.5],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0366,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13241,0.0181319,-1.74093e-05,9.35336e-09,-2.01725e-12,-7148.09,13.3808], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.75871,0.00635124,-2.25955e-06,3.62322e-10,-2.15856e-14,-8085.33,-4.9649], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-60.8183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CH2CO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'C=C([O])C(F)[C](F)F(3309)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  O u1 p2 c0 {6,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {4,S} {5,S} {8,D}
7  C u1 p0 c0 {2,S} {3,S} {5,S}
8  C u0 p0 c0 {6,D} {10,S} {11,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-460.809,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,350,440,435,1725,190,488,555,1236,1407,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0119486,'amu*angstrom^2'), symmetry=1, barrier=(8.49497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.837418,'amu*angstrom^2'), symmetry=1, barrier=(19.2539,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3544.65,'J/mol'), sigma=(5.68142,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=553.67 K, Pc=43.86 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08051,0.0692867,-9.12708e-05,6.48539e-08,-1.8729e-11,-55321.9,25.3045], Tmin=(100,'K'), Tmax=(840.016,'K')), NASAPolynomial(coeffs=[10.3176,0.025303,-1.27332e-05,2.52618e-09,-1.80223e-13,-56873.8,-17.6493], Tmin=(840.016,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-460.809,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Csj(Cs-F1sCdH)(F1s)(F1s))"""),
)

species(
    label = '[O]C(CF)=C(F)[CH]F(3484)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u1 p2 c0 {6,S}
5  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
6  C u0 p0 c0 {4,S} {5,S} {7,D}
7  C u0 p0 c0 {2,S} {6,D} {8,S}
8  C u1 p0 c0 {3,S} {7,S} {11,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-464.999,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([212,931,1104,1251,1325,1474,3102,3155,350,440,435,1725,271,519,563,612,1379,234,589,736,816,1240,3237,180,1776.45],'cm^-1')),
        HinderedRotor(inertia=(2.23717,'amu*angstrom^2'), symmetry=1, barrier=(51.4369,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.214211,'amu*angstrom^2'), symmetry=1, barrier=(4.92513,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.900111,0.0748152,-0.000115949,1.00966e-07,-3.54093e-11,-55821.3,24.8929], Tmin=(100,'K'), Tmax=(781.106,'K')), NASAPolynomial(coeffs=[8.28294,0.0296849,-1.52196e-05,2.9912e-09,-2.10108e-13,-56751.2,-7.47106], Tmin=(781.106,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-464.999,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(Cds-CdsCsOs) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(C=C(C)OJ) + radical(Csj(Cd-CdF1s)(F1s)(H))"""),
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
    label = '[CH2]C(=O)[C](F)F(3383)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {3,D} {5,S} {6,S}
5 C u1 p0 c0 {4,S} {7,S} {8,S}
6 C u1 p0 c0 {1,S} {2,S} {4,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-279.314,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,3000,3100,440,815,1455,1000,179,346,818,1406,1524,241.327],'cm^-1')),
        HinderedRotor(inertia=(0.0200138,'amu*angstrom^2'), symmetry=1, barrier=(12.2066,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.28776,'amu*angstrom^2'), symmetry=1, barrier=(53.2065,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0441,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00969,0.0471597,-6.16119e-05,4.28745e-08,-1.20388e-11,-33525,18.4741], Tmin=(100,'K'), Tmax=(865.241,'K')), NASAPolynomial(coeffs=[8.69208,0.016267,-8.05515e-06,1.60882e-09,-1.15551e-13,-34681.4,-12.7972], Tmin=(865.241,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-279.314,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(CO-CsO2d)(F1s)(F1s)) + radical(C2JC=O)"""),
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
    label = 'C=[C]C(F)(F)[CH]F(2834)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  F u0 p3 c0 {4,S}
3  F u0 p3 c0 {5,S}
4  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
5  C u1 p0 c0 {3,S} {4,S} {8,S}
6  C u0 p0 c0 {7,D} {9,S} {10,S}
7  C u1 p0 c0 {4,S} {6,D}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-177.914,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([136,307,446,511,682,757,1180,1185,334,575,1197,1424,3202,2950,3100,1380,975,1025,1650,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.17789,'amu*angstrom^2'), symmetry=1, barrier=(4.09005,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.45853,'amu*angstrom^2'), symmetry=1, barrier=(10.5425,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33505,0.0636796,-9.6896e-05,8.19008e-08,-2.77841e-11,-21307,24.785], Tmin=(100,'K'), Tmax=(783.68,'K')), NASAPolynomial(coeffs=[8.42788,0.0232294,-1.13423e-05,2.20521e-09,-1.54239e-13,-22288.3,-6.87274], Tmin=(783.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-177.914,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sF1sCd)(F1s)(H)) + radical(Cdj(Cs-F1sF1sCs)(Cd-HH))"""),
)

species(
    label = 'C=C1OC(F)C1(F)F(3312)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {6,S}
4  O u0 p2 c0 {6,S} {7,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u0 p0 c0 {3,S} {4,S} {5,S} {9,S}
7  C u0 p0 c0 {4,S} {5,S} {8,D}
8  C u0 p0 c0 {7,D} {10,S} {11,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-684.923,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06487,0.0489146,-4.46949e-06,-4.48364e-08,2.53218e-11,-82257,16.407], Tmin=(100,'K'), Tmax=(919.08,'K')), NASAPolynomial(coeffs=[20.6392,0.00466784,9.20703e-07,-2.74983e-10,1.58734e-14,-87584.4,-85.7837], Tmin=(919.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-684.923,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(CsCCFF) + group(CsCFHO) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + ring(2methyleneoxetane)"""),
)

species(
    label = 'F[CH]C(F)(F)[C]1CO1(3485)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {7,S}
5  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
6  C u0 p0 c0 {4,S} {7,S} {9,S} {10,S}
7  C u1 p0 c0 {4,S} {5,S} {6,S}
8  C u1 p0 c0 {3,S} {5,S} {11,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-358.27,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1262,0.072413,-0.000121247,1.14564e-07,-4.2176e-11,-42995.3,25.2707], Tmin=(100,'K'), Tmax=(824.909,'K')), NASAPolynomial(coeffs=[5.27638,0.033777,-1.73308e-05,3.37798e-09,-2.34861e-13,-43050.2,9.86491], Tmin=(824.909,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-358.27,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsOsHH) + group(CsCsFHH) + ring(Cs(C-FF)-Cs-O2s) + radical(C2CsJO) + radical(Csj(Cs-CsF1sF1s)(F1s)(H))"""),
)

species(
    label = '[O][C]1CC(F)C1(F)F(3486)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  O u1 p2 c0 {8,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
7  C u0 p0 c0 {5,S} {8,S} {10,S} {11,S}
8  C u1 p0 c0 {4,S} {6,S} {7,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-380.358,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34909,0.055724,-4.97687e-05,2.2376e-08,-4.0464e-12,-45649,22.6051], Tmin=(100,'K'), Tmax=(1314.73,'K')), NASAPolynomial(coeffs=[13.0731,0.0200542,-9.07233e-06,1.73982e-09,-1.22371e-13,-48731.8,-37.1646], Tmin=(1314.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-380.358,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCsCsFF) + group(Cs-CsCsHH) + group(CsCsCsFH) + ring(Cs-Cs-Cs(F)(F)-Cs) + radical(CC(C)OJ) + radical(C2CsJOH) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F))"""),
)

species(
    label = '[CH2]C1([O])C(F)C1(F)F(3447)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  O u1 p2 c0 {5,S}
5  C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
6  C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
7  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
8  C u1 p0 c0 {5,S} {10,S} {11,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-288.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.93071,0.0689546,-8.03213e-05,4.78776e-08,-1.13671e-11,-34617.4,22.4169], Tmin=(100,'K'), Tmax=(1023.56,'K')), NASAPolynomial(coeffs=[13.2601,0.0207716,-9.70995e-06,1.88652e-09,-1.33926e-13,-37141.4,-37.3526], Tmin=(1023.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-288.735,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(CsCsCsFH) + group(CsCsCsFF) + group(Cs-CsHHH) + ring(Cs(F)(F)-Cs(C)-Cs) + radical(CC(C)2OJ) + radical(CJC(C)2O) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(3ring-Cs(F)2-Cs(F))"""),
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
    label = 'C=C([O])C(F)=CF(3487)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  O u1 p2 c0 {4,S}
4  C u0 p0 c0 {3,S} {5,S} {6,D}
5  C u0 p0 c0 {1,S} {4,S} {7,D}
6  C u0 p0 c0 {4,D} {8,S} {9,S}
7  C u0 p0 c0 {2,S} {5,D} {10,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-344.129,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,280,518,736,852,873,2950,3100,1380,975,1025,1650,194,682,905,1196,1383,3221,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.522041,'amu*angstrom^2'), symmetry=1, barrier=(12.0028,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28217,0.0600102,-7.10476e-05,4.27072e-08,-1.01261e-11,-41291.4,19.8397], Tmin=(100,'K'), Tmax=(1030.73,'K')), NASAPolynomial(coeffs=[12.6043,0.0160714,-7.10391e-06,1.34861e-09,-9.46481e-14,-43625.4,-35.1256], Tmin=(1030.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-344.129,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(CdCCF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(Cds-CdsHH) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2][C]=O(2901)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,D}
2 C u1 p0 c0 {3,S} {4,S} {5,S}
3 C u1 p0 c0 {1,D} {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (160.185,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,539.612,539.669],'cm^-1')),
        HinderedRotor(inertia=(0.000578908,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0366,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2980.68,'J/mol'), sigma=(5.03063,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=465.58 K, Pc=53.12 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39564,0.0101365,2.30772e-06,-8.97606e-09,3.68258e-12,19290.3,10.0702], Tmin=(100,'K'), Tmax=(1068.89,'K')), NASAPolynomial(coeffs=[6.35051,0.00638958,-2.69372e-06,5.4222e-10,-4.02484e-14,18240.9,-6.33577], Tmin=(1068.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CJC=O) + radical(CsCJ=O)"""),
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
    label = '[CH]=C(O)C(F)(F)[CH]F(3488)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {6,S} {10,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u0 p0 c0 {4,S} {5,S} {8,D}
7  C u1 p0 c0 {3,S} {5,S} {9,S}
8  C u1 p0 c0 {6,D} {11,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-358.5,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,248,333,466,604,684,796,1061,1199,350,440,435,1725,334,575,1197,1424,3202,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.756632,'amu*angstrom^2'), symmetry=1, barrier=(17.3965,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.755687,'amu*angstrom^2'), symmetry=1, barrier=(17.3747,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.756906,'amu*angstrom^2'), symmetry=1, barrier=(17.4028,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.431713,0.082632,-0.000118578,8.49038e-08,-2.38206e-11,-42992.7,26.5619], Tmin=(100,'K'), Tmax=(876.611,'K')), NASAPolynomial(coeffs=[14.4212,0.0187972,-9.34718e-06,1.833e-09,-1.29575e-13,-45445.3,-39.0869], Tmin=(876.611,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-358.5,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Csj(Cs-F1sF1sCd)(F1s)(H)) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])C(F)(F)CF(3489)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {6,S}
4  O u1 p2 c0 {7,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u0 p0 c0 {3,S} {5,S} {9,S} {10,S}
7  C u0 p0 c0 {4,S} {5,S} {8,D}
8  C u1 p0 c0 {7,D} {11,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-421.801,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([274,345,380,539,705,1166,1213,528,1116,1182,1331,1402,1494,3075,3110,350,440,435,1725,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.654295,'amu*angstrom^2'), symmetry=1, barrier=(15.0435,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.6532,'amu*angstrom^2'), symmetry=1, barrier=(15.0184,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.522879,0.0826699,-0.00013147,1.09949e-07,-3.63464e-11,-50611.7,24.6929], Tmin=(100,'K'), Tmax=(801.711,'K')), NASAPolynomial(coeffs=[11.1111,0.0244051,-1.22848e-05,2.38114e-09,-1.65383e-13,-52134.7,-22.9592], Tmin=(801.711,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-421.801,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C(OF)=C(F)[CH]F(3490)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {5,S}
5  C u0 p0 c0 {4,S} {6,D} {7,S}
6  C u0 p0 c0 {1,S} {5,D} {8,S}
7  C u1 p0 c0 {5,S} {9,S} {10,S}
8  C u1 p0 c0 {2,S} {6,S} {11,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-130.322,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,350,440,435,1725,271,519,563,612,1379,3000,3100,440,815,1455,1000,234,589,736,816,1240,3237,180],'cm^-1')),
        HinderedRotor(inertia=(0.00265767,'amu*angstrom^2'), symmetry=1, barrier=(8.83827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.02482,'amu*angstrom^2'), symmetry=1, barrier=(46.5545,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.013983,'amu*angstrom^2'), symmetry=1, barrier=(46.4969,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.692582,0.0801384,-0.000129846,1.13692e-07,-3.92539e-11,-15562.2,26.1553], Tmin=(100,'K'), Tmax=(815.848,'K')), NASAPolynomial(coeffs=[8.99691,0.0283914,-1.44217e-05,2.7992e-09,-1.94316e-13,-16550,-9.96801], Tmin=(815.848,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-130.322,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-Cds)HHH) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cds(F)) + group(Cds-CdsCsOs) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(C=C(O)CJ) + radical(Csj(Cd-CdF1s)(F1s)(H))"""),
)

species(
    label = 'C=C([O])[C](F)C(F)F(3450)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {6,S}
4  O u1 p2 c0 {7,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
6  C u1 p0 c0 {3,S} {5,S} {7,S}
7  C u0 p0 c0 {4,S} {6,S} {8,D}
8  C u0 p0 c0 {7,D} {10,S} {11,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-526.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.896753,0.0691355,-8.05878e-05,4.74619e-08,-1.10761e-11,-63172.3,24.3264], Tmin=(100,'K'), Tmax=(1043.16,'K')), NASAPolynomial(coeffs=[13.8625,0.0194175,-9.09553e-06,1.77174e-09,-1.26068e-13,-65877.3,-38.7737], Tmin=(1043.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-526.168,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CdOs) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(CsCdCsF1s)"""),
)

species(
    label = 'O=[C]CC(F)(F)[CH]F(3308)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
7  C u1 p0 c0 {3,S} {5,S} {11,S}
8  C u1 p0 c0 {4,D} {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-481.091,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,2750,2850,1437.5,1250,1305,750,350,334,575,1197,1424,3202,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.240731,'amu*angstrom^2'), symmetry=1, barrier=(5.53488,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.240973,'amu*angstrom^2'), symmetry=1, barrier=(5.54043,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.240847,'amu*angstrom^2'), symmetry=1, barrier=(5.53756,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3296.48,'J/mol'), sigma=(5.67221,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=514.90 K, Pc=40.99 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.512626,0.0889179,-0.000164637,1.544e-07,-5.46255e-11,-57748,27.1661], Tmin=(100,'K'), Tmax=(860.916,'K')), NASAPolynomial(coeffs=[7.42999,0.0300879,-1.56325e-05,3.00437e-09,-2.04871e-13,-57949.9,0.574332], Tmin=(860.916,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-481.091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFHH) + group(Cds-OdCsH) + radical(Csj(Cs-CsF1sF1s)(F1s)(H)) + radical(CCCJ=O)"""),
)

species(
    label = 'CH2(T)(66)',
    structure = adjacencyList("""multiplicity 3
1 C u2 p0 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (381.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1066.91,2790.97,3622.38],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.01192,-0.000154978,3.26297e-06,-2.40422e-09,5.69496e-13,45867.7,0.533201], Tmin=(100,'K'), Tmax=(1104.65,'K')), NASAPolynomial(coeffs=[3.14983,0.00296674,-9.76056e-07,1.54115e-10,-9.50339e-15,46058.1,4.77808], Tmin=(1104.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'O=[C]C(F)(F)[CH]F(446)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {6,S}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6 C u1 p0 c0 {3,S} {5,S} {8,S}
7 C u1 p0 c0 {4,D} {5,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-418.438,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([169,350,401,508,725,1243,1223,334,575,1197,1424,3202,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.433749,'amu*angstrom^2'), symmetry=1, barrier=(9.97275,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159145,'amu*angstrom^2'), symmetry=1, barrier=(3.65906,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39354,0.06634,-0.000124963,1.18098e-07,-4.2238e-11,-50241.3,22.4052], Tmin=(100,'K'), Tmax=(843.468,'K')), NASAPolynomial(coeffs=[7.00178,0.0211417,-1.15012e-05,2.2718e-09,-1.57747e-13,-50525.6,0.226044], Tmin=(843.468,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-418.438,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sF1sCO)(F1s)(H)) + radical(COj(Cs-F1sF1sCs)(O2d))"""),
)

species(
    label = 'O=C1CC(F)C1(F)F(3416)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {5,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
8  C u0 p0 c0 {4,D} {6,S} {7,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-698.637,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48067,0.0479342,-3.27595e-05,9.04961e-09,-4.44664e-13,-83929.7,22.2405], Tmin=(100,'K'), Tmax=(1224.84,'K')), NASAPolynomial(coeffs=[12.9842,0.0191411,-8.2436e-06,1.55448e-09,-1.08587e-13,-87405.9,-38.2772], Tmin=(1224.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-698.637,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + group(Cs-(Cds-O2d)CsHH) + group(CsCCFF) + group(Cds-OdCsCs) + longDistanceInteraction_cyclic(Cs(F)2-CO) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + ring(Cyclobutanone)"""),
)

species(
    label = '[CH2][C]1OC(F)C1(F)F(3491)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {6,S}
4  O u0 p2 c0 {6,S} {7,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u0 p0 c0 {3,S} {4,S} {5,S} {9,S}
7  C u1 p0 c0 {4,S} {5,S} {8,S}
8  C u1 p0 c0 {7,S} {10,S} {11,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-396.853,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.903204,0.0752949,-0.000121209,1.06654e-07,-3.69216e-11,-47625.8,24.2984], Tmin=(100,'K'), Tmax=(825.615,'K')), NASAPolynomial(coeffs=[8.21028,0.0280692,-1.39259e-05,2.67958e-09,-1.85065e-13,-48429.3,-7.11327], Tmin=(825.615,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-396.853,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(CsCsCsFF) + group(CsCFHO) + group(Cs-CsHHH) + ring(O2s-Cs-Cs-Cs(F)) + radical(C2CsJOCs) + radical(CJC(C)OC) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
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
    E0 = (-221.353,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-54.436,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (88.6634,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (182.052,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (311.558,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-213.069,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (9.86338,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-96.1886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-42.2973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (35.1026,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-74.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-92.5557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (121.37,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (105.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (80.0674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-131.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (210.093,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-67.1729,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-64.0344,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (209.37,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-213.069,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (20.8171,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C([O])C(F)(F)[CH]F(3306)'],
    products = ['CHFCF2(54)', 'CH2CO(75)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=C([O])C(F)[C](F)F(3309)'],
    products = ['C=C([O])C(F)(F)[CH]F(3306)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]C(CF)=C(F)[CH]F(3484)'],
    products = ['C=C([O])C(F)(F)[CH]F(3306)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(307.225,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CHF(40)', '[CH2]C(=O)[C](F)F(3383)'],
    products = ['C=C([O])C(F)(F)[CH]F(3306)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(6)', 'C=[C]C(F)(F)[CH]F(2834)'],
    products = ['C=C([O])C(F)(F)[CH]F(3306)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=C([O])C(F)(F)[CH]F(3306)'],
    products = ['C=C1OC(F)C1(F)F(3312)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_1H;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=C([O])C(F)(F)[CH]F(3306)'],
    products = ['F[CH]C(F)(F)[C]1CO1(3485)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R3_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=C([O])C(F)(F)[CH]F(3306)'],
    products = ['[O][C]1CC(F)C1(F)F(3486)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.94431e+07,'s^-1'), n=1.16299, Ea=(125.164,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=C([O])C(F)(F)[CH]F(3306)'],
    products = ['[CH2]C1([O])C(F)C1(F)F(3447)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.17979e+10,'s^-1'), n=0.520585, Ea=(179.056,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 178.5 to 179.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction10',
    reactants = ['F(37)', 'C=C([O])C(F)=CF(3487)'],
    products = ['C=C([O])C(F)(F)[CH]F(3306)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(59.9017,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CHFCF2(54)', '[CH2][C]=O(2901)'],
    products = ['C=C([O])C(F)(F)[CH]F(3306)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1256.55,'m^3/(mol*s)'), n=0.85398, Ea=(30.6947,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.03938838252749272, var=0.8044325964401412, Tref=1000.0, N=80, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-1R!H-R_N-5R!H-inRing_Ext-1R!H-R_Ext-2C-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_2R!H->C_Ext-1R!H-R_N-5R!H-inRing_Ext-1R!H-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F[CH][C](F)F(183)', 'CH2CO(75)'],
    products = ['C=C([O])C(F)(F)[CH]F(3306)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(7.07808,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F[CH][C](F)F(183)', '[CH2][C]=O(2901)'],
    products = ['C=C([O])C(F)(F)[CH]F(3306)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['CHF(86)', '[CH2]C(=O)[C](F)F(3383)'],
    products = ['C=C([O])C(F)(F)[CH]F(3306)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(786.723,'m^3/(mol*s)'), n=1.25031, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=C(O)C(F)(F)[CH]F(3488)'],
    products = ['C=C([O])C(F)(F)[CH]F(3306)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=C([O])C(F)(F)CF(3489)'],
    products = ['C=C([O])C(F)(F)[CH]F(3306)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_1H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(OF)=C(F)[CH]F(3490)'],
    products = ['C=C([O])C(F)(F)[CH]F(3306)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(93.9767,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=C([O])C(F)(F)[CH]F(3306)'],
    products = ['C=C([O])[C](F)C(F)F(3450)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.01526e+12,'s^-1'), n=0.18834, Ea=(154.18,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=C([O])C(F)(F)[CH]F(3306)'],
    products = ['O=[C]CC(F)(F)[CH]F(3308)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction20',
    reactants = ['CH2(T)(66)', 'O=[C]C(F)(F)[CH]F(446)'],
    products = ['C=C([O])C(F)(F)[CH]F(3306)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=C([O])C(F)(F)[CH]F(3306)'],
    products = ['O=C1CC(F)C1(F)F(3416)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_1H;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=C([O])C(F)(F)[CH]F(3306)'],
    products = ['[CH2][C]1OC(F)C1(F)F(3491)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(242.17,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_cs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

network(
    label = 'PDepNetwork #1265',
    isomers = [
        'C=C([O])C(F)(F)[CH]F(3306)',
    ],
    reactants = [
        ('CHFCF2(54)', 'CH2CO(75)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1265',
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

