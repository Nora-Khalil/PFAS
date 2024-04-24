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
    label = 'CO(12)',
    structure = adjacencyList("""1 O u0 p1 c+1 {2,T}
2 C u0 p1 c-1 {1,T}
"""),
    E0 = (-118.741,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2193.04],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.01,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(815.652,'J/mol'), sigma=(3.65,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.95,'angstroms^3'), rotrelaxcollnum=1.8, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.56838,-0.000852125,2.48918e-06,-1.56331e-09,3.13595e-13,-14284.3,3.57912], Tmin=(100,'K'), Tmax=(1571.64,'K')), NASAPolynomial(coeffs=[2.91307,0.00164658,-6.88615e-07,1.21037e-10,-7.84018e-15,-14180.9,6.71045], Tmin=(1571.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-118.741,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""CO""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[O]C(F)C[C](F)F(799)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  O u1 p2 c0 {6,S}
5  C u0 p0 c0 {6,S} {7,S} {8,S} {9,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {10,S}
7  C u1 p0 c0 {2,S} {3,S} {5,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-488.826,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,391,562,707,872,1109,1210,1289,3137,190,488,555,1236,1407,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.1756,'amu*angstrom^2'), symmetry=1, barrier=(4.03738,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.175281,'amu*angstrom^2'), symmetry=1, barrier=(4.03006,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.05,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3383.02,'J/mol'), sigma=(5.55661,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=528.42 K, Pc=44.74 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67671,0.0551558,-6.70751e-05,4.51592e-08,-1.25048e-11,-58712,22.8831], Tmin=(100,'K'), Tmax=(870.329,'K')), NASAPolynomial(coeffs=[8.76219,0.0225906,-1.09485e-05,2.16584e-09,-1.54866e-13,-59945.3,-10.3161], Tmin=(870.329,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-488.826,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sH)) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = '[CH2]C(F)(F)C([O])C(=O)F(3708)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  O u1 p2 c0 {6,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {4,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
8  C u1 p0 c0 {7,S} {11,S} {12,S}
9  C u0 p0 c0 {3,S} {5,D} {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-639.423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.628122,0.0789636,-0.000101574,6.74518e-08,-1.79464e-11,-76787.5,31.5643], Tmin=(100,'K'), Tmax=(913.811,'K')), NASAPolynomial(coeffs=[13.0164,0.0247368,-1.25618e-05,2.51341e-09,-1.80616e-13,-79051.6,-27.0853], Tmin=(913.811,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-639.423,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCsCsFF) + group(Cs-CsHHH) + group(COCsFO) + radical(C=OCOJ) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
)

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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3461.84,'J/mol'), sigma=(5.59194,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=540.73 K, Pc=44.92 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.465973,0.0871166,-0.000143504,1.30944e-07,-4.71509e-11,-76526.4,28.7803], Tmin=(100,'K'), Tmax=(807.79,'K')), NASAPolynomial(coeffs=[7.65037,0.035717,-1.86751e-05,3.67257e-09,-2.57005e-13,-77170.8,-1.15107], Tmin=(807.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-637.259,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsOsHH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFFH) + group(Cds-OdCsH) + radical(CsCOF1sO2s) + radical(Csj(Cs-O2sHH)(F1s)(F1s))"""),
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
    label = 'O=C(F)[CH]C[C](F)F(2720)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {6,S} {7,S} {9,S} {10,S}
6  C u1 p0 c0 {5,S} {8,S} {11,S}
7  C u1 p0 c0 {1,S} {2,S} {5,S}
8  C u0 p0 c0 {3,S} {4,D} {6,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-528.767,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,190,488,555,1236,1407,611,648,830,1210,1753,180,180,619.636],'cm^-1')),
        HinderedRotor(inertia=(0.152173,'amu*angstrom^2'), symmetry=1, barrier=(3.49876,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0128451,'amu*angstrom^2'), symmetry=1, barrier=(3.49955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.20637,'amu*angstrom^2'), symmetry=1, barrier=(50.7288,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25145,0.066115,-8.60905e-05,6.34542e-08,-1.94495e-11,-63502.1,25.9142], Tmin=(100,'K'), Tmax=(786.37,'K')), NASAPolynomial(coeffs=[8.60831,0.028693,-1.47081e-05,2.93777e-09,-2.10279e-13,-64659.2,-7.81031], Tmin=(786.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-528.767,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFFH) + group(COCsFO) + radical(CCJC=O) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
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
    label = '[CH2]C([O])C(=O)F(3980)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u1 p2 c0 {4,S}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
5 C u1 p0 c0 {4,S} {8,S} {9,S}
6 C u0 p0 c0 {1,S} {3,D} {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-171.875,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,486,617,768,1157,1926,323.978,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00160297,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.369234,'amu*angstrom^2'), symmetry=1, barrier=(27.3591,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (90.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94289,0.0408922,-3.76385e-05,1.73638e-08,-3.16312e-12,-20594.2,23.9303], Tmin=(100,'K'), Tmax=(1328.26,'K')), NASAPolynomial(coeffs=[11.5998,0.0118109,-4.79712e-06,8.80376e-10,-6.067e-14,-23159.5,-25.3999], Tmin=(1328.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-171.875,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(COCsFO) + radical(C=OCOJ) + radical(CJCO)"""),
)

species(
    label = 'O=C(F)C1CC(F)(F)O1(3718)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {4,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {6,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
9  C u0 p0 c0 {3,S} {5,D} {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-949.364,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34604,0.045971,-8.08113e-06,-2.35823e-08,1.19598e-11,-114076,27.0192], Tmin=(100,'K'), Tmax=(1040.67,'K')), NASAPolynomial(coeffs=[15.3666,0.0191753,-8.51293e-06,1.71332e-09,-1.27352e-13,-118461,-48.2301], Tmin=(1040.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-949.364,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(CsCFFO) + group(COCsFO) + ring(Cs-Cs(F)(F)-O2s-Cs)"""),
)

species(
    label = 'O=C(F)C(=O)CC(F)F(3981)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {8,D}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {7,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {12,S}
8  C u0 p0 c0 {4,D} {6,S} {9,S}
9  C u0 p0 c0 {3,S} {5,D} {8,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-985.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.510792,0.0840477,-0.000130791,1.13222e-07,-3.95953e-11,-118368,27.0805], Tmin=(100,'K'), Tmax=(763.507,'K')), NASAPolynomial(coeffs=[9.30013,0.0319935,-1.67234e-05,3.31754e-09,-2.34437e-13,-119535,-11.8048], Tmin=(763.507,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-985.152,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(CsCsFFH) + group(Cds-O2d(Cds-O2d)Cs) + group(COCFO)"""),
)

species(
    label = 'O=C(F)C(O)C=C(F)F(3982)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {6,S} {12,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {4,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {6,S} {9,D} {11,S}
8  C u0 p0 c0 {1,S} {5,D} {6,S}
9  C u0 p0 c0 {2,S} {3,S} {7,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {4,S}
"""),
    E0 = (-917.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01694,0.0661093,-6.81811e-05,3.56633e-08,-7.49514e-12,-110201,29.7534], Tmin=(100,'K'), Tmax=(1142.27,'K')), NASAPolynomial(coeffs=[13.5203,0.0223252,-1.06848e-05,2.10661e-09,-1.50826e-13,-113057,-32.2311], Tmin=(1142.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-917.152,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(COCsFO) + group(CdCFF)"""),
)

species(
    label = 'F[C](F)CC1OO[C]1F(3983)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {5,S} {6,S}
5  O u0 p2 c0 {4,S} {8,S}
6  C u0 p0 c0 {4,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {6,S} {9,S} {11,S} {12,S}
8  C u1 p0 c0 {1,S} {5,S} {6,S}
9  C u1 p0 c0 {2,S} {3,S} {7,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-313.954,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.770618,0.0777602,-0.000117048,1.01066e-07,-3.55423e-11,-37650.1,29.7124], Tmin=(100,'K'), Tmax=(764.907,'K')), NASAPolynomial(coeffs=[8.27489,0.0325247,-1.65882e-05,3.26697e-09,-2.30236e-13,-38622.8,-3.33422], Tmin=(764.907,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-313.954,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCFHO) + group(CsCsFFH) + ring(Cs-Cs(F)-O2s-O2s) + radical(CsCsF1sO2s) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = '[O]C1CC(F)(F)O[C]1F(3984)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {8,S} {9,S}
5  O u1 p2 c0 {7,S}
6  C u0 p0 c0 {7,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {5,S} {6,S} {9,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
9  C u1 p0 c0 {3,S} {4,S} {7,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-632.243,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.295119,0.0660536,-6.14204e-05,2.9053e-08,-5.25222e-12,-75894.6,24.3394], Tmin=(100,'K'), Tmax=(1531.88,'K')), NASAPolynomial(coeffs=[16.8532,0.01482,-3.42189e-06,4.04211e-10,-2.06113e-14,-80029.2,-59.5428], Tmin=(1531.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-632.243,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCFHO) + group(CsCFFO) + ring(Tetrahydrofuran) + radical(CC(C)OJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]C1(F)OC1C[C](F)F(3985)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u1 p2 c0 {8,S}
6  C u0 p0 c0 {4,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {6,S} {9,S} {11,S} {12,S}
8  C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
9  C u1 p0 c0 {2,S} {3,S} {7,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-494.814,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12613,0.069045,-8.15887e-05,5.43937e-08,-1.52194e-11,-59414,26.7716], Tmin=(100,'K'), Tmax=(853.876,'K')), NASAPolynomial(coeffs=[9.09372,0.0317195,-1.60171e-05,3.19677e-09,-2.29381e-13,-60774.6,-10.4085], Tmin=(853.876,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-494.814,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCFOO) + group(CsCsFFH) + ring(Cs(F)(O2)-O2s-Cs) + radical(O2sj(Cs-F1sO2sCs)) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = '[O]C1CC(F)(F)C1([O])F(3986)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u1 p2 c0 {6,S}
5  O u1 p2 c0 {9,S}
6  C u0 p0 c0 {4,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {6,S} {8,S} {11,S} {12,S}
8  C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
9  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-527.069,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.807852,0.0639167,-5.7167e-05,2.51147e-08,-4.36495e-12,-63271.7,26.4336], Tmin=(100,'K'), Tmax=(1382.19,'K')), NASAPolynomial(coeffs=[16.3287,0.0189997,-8.42115e-06,1.60317e-09,-1.12321e-13,-67562.2,-53.4689], Tmin=(1382.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-527.069,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCCFO) + group(Cs-CsCsHH) + group(CsCsCsFF) + ring(Cs-Cs-Cs(F)(F)-Cs) + radical(CC(C)OJ) + radical(O2sj(Cs-F1sCsCs)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'O=C(C[C](F)F)[C](O)F(3987)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {8,S} {12,S}
5  O u0 p2 c0 {7,D}
6  C u0 p0 c0 {7,S} {9,S} {10,S} {11,S}
7  C u0 p0 c0 {5,D} {6,S} {8,S}
8  C u1 p0 c0 {3,S} {4,S} {7,S}
9  C u1 p0 c0 {1,S} {2,S} {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {4,S}
"""),
    E0 = (-701.893,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,280,501,1494,1531,190,488,555,1236,1407,180,180,1292.48],'cm^-1')),
        HinderedRotor(inertia=(0.327997,'amu*angstrom^2'), symmetry=1, barrier=(7.54129,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.327606,'amu*angstrom^2'), symmetry=1, barrier=(7.53231,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.327462,'amu*angstrom^2'), symmetry=1, barrier=(7.52898,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.42588,'amu*angstrom^2'), symmetry=1, barrier=(55.7757,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.151806,0.0981826,-0.000180103,1.69131e-07,-6.00044e-11,-84292.7,28.739], Tmin=(100,'K'), Tmax=(860.832,'K')), NASAPolynomial(coeffs=[7.38147,0.0348599,-1.79608e-05,3.44263e-09,-2.34591e-13,-84435.9,1.34132], Tmin=(860.832,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-701.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsHH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFFH) + group(Cds-OdCsCs) + radical(CsCOF1sO2s) + radical(Csj(Cs-COHH)(F1s)(F1s))"""),
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
    label = 'CFO(50)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u1 p0 c0 {1,S} {2,D}
"""),
    E0 = (-190.359,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(46.9933,'amu')),
        NonlinearRotor(inertia=([2.65864,43.9175,46.5761],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([636.046,1085.22,1918.24],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (47.0084,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(7150.45,'J/mol'), sigma=(4,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02994,-0.00208576,2.20825e-05,-3.30628e-08,1.5841e-11,-22895,6.85532], Tmin=(10,'K'), Tmax=(668.186,'K')), NASAPolynomial(coeffs=[3.39014,0.00554951,-3.6e-06,1.08417e-09,-1.23809e-13,-22894.5,9.04838], Tmin=(668.186,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-190.359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""OD[C]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=CC[C](F)F(911)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
5 C u1 p0 c0 {1,S} {2,S} {4,S}
6 C u0 p0 c0 {3,D} {4,S} {9,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-411.653,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,190,488,555,1236,1407,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.173403,'amu*angstrom^2'), symmetry=1, barrier=(3.98689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172994,'amu*angstrom^2'), symmetry=1, barrier=(3.97748,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.83124,0.0223284,0.000421225,-3.79918e-06,9.54663e-09,-49508.1,10.3772], Tmin=(10,'K'), Tmax=(147.606,'K')), NASAPolynomial(coeffs=[5.19566,0.0253007,-1.49258e-05,4.17585e-09,-4.48631e-13,-49591.9,4.93097], Tmin=(147.606,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-411.653,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""ODCC[C](F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'O=C(F)C(=O)C[C](F)F(3988)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {7,D}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {7,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {4,D} {6,S} {9,S}
8  C u1 p0 c0 {1,S} {2,S} {6,S}
9  C u0 p0 c0 {3,S} {5,D} {7,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-783.833,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,190,488,555,1236,1407,286,619,818,1246,1924,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.139809,'amu*angstrom^2'), symmetry=1, barrier=(3.21448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139768,'amu*angstrom^2'), symmetry=1, barrier=(3.21354,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.562915,'amu*angstrom^2'), symmetry=1, barrier=(12.9425,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (139.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.328558,0.091472,-0.000163408,1.49716e-07,-5.26133e-11,-94151.5,27.6462], Tmin=(100,'K'), Tmax=(840.119,'K')), NASAPolynomial(coeffs=[8.97219,0.0295347,-1.57142e-05,3.06951e-09,-2.12237e-13,-94870.4,-8.18348], Tmin=(840.119,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-783.833,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(CsCsFFH) + group(Cds-O2d(Cds-O2d)Cs) + group(COCFO) + radical(Csj(Cs-COHH)(F1s)(F1s))"""),
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
    label = '[O]C(C=C(F)F)C(=O)F(3822)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  O u1 p2 c0 {6,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {4,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {6,S} {9,D} {11,S}
8  C u0 p0 c0 {1,S} {5,D} {6,S}
9  C u0 p0 c0 {2,S} {3,S} {7,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-673.419,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,486,617,768,1157,1926,182,240,577,636,1210,1413,308.156,310.011,3749.21],'cm^-1')),
        HinderedRotor(inertia=(0.00909944,'amu*angstrom^2'), symmetry=1, barrier=(11.5746,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.361323,'amu*angstrom^2'), symmetry=1, barrier=(24.1993,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (139.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29124,0.0627057,-6.99733e-05,4.08506e-08,-9.69873e-12,-80898.7,29.5558], Tmin=(100,'K'), Tmax=(1011.32,'K')), NASAPolynomial(coeffs=[11.2662,0.0232531,-1.14578e-05,2.27744e-09,-1.63535e-13,-82916.3,-18.68], Tmin=(1011.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-673.419,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(COCsFO) + group(CdCFF) + radical(C=OCOJ)"""),
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
    label = 'O=[C]C(=O)C[C](F)F(3973)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {6,D}
4  O u0 p2 c0 {8,D}
5  C u0 p0 c0 {6,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {3,D} {5,S} {8,S}
7  C u1 p0 c0 {1,S} {2,S} {5,S}
8  C u1 p0 c0 {4,D} {6,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-360.377,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,190,488,555,1236,1407,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.268422,'amu*angstrom^2'), symmetry=1, barrier=(6.17154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.268537,'amu*angstrom^2'), symmetry=1, barrier=(6.17419,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.557671,'amu*angstrom^2'), symmetry=1, barrier=(12.822,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.65718,0.0847461,-0.000157724,1.46584e-07,-5.13876e-11,-43233.7,26.2791], Tmin=(100,'K'), Tmax=(861.901,'K')), NASAPolynomial(coeffs=[8.03356,0.0262847,-1.38156e-05,2.65892e-09,-1.81187e-13,-43605.3,-2.99084], Tmin=(861.901,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-360.377,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(CsCsFFH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(Csj(Cs-COHH)(F1s)(F1s)) + radical(CCCJ=O)"""),
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
    label = 'O=C(F)[C](O)C[C](F)F(3989)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {7,S} {12,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {7,S} {8,S} {10,S} {11,S}
7  C u1 p0 c0 {4,S} {6,S} {9,S}
8  C u1 p0 c0 {1,S} {2,S} {6,S}
9  C u0 p0 c0 {3,S} {5,D} {7,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {4,S}
"""),
    E0 = (-684.15,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.330035,0.0877846,-0.000139784,1.18751e-07,-4.01604e-11,-82158.8,33.0581], Tmin=(100,'K'), Tmax=(782.342,'K')), NASAPolynomial(coeffs=[10.9151,0.0280601,-1.45279e-05,2.85773e-09,-2.00473e-13,-83643.5,-14.3145], Tmin=(782.342,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-684.15,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCsFFH) + group(COCsFO) + radical(C2CsJOH) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = '[O]C([CH]C(F)F)C(=O)F(3990)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  O u1 p2 c0 {6,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {11,S}
8  C u1 p0 c0 {6,S} {7,S} {12,S}
9  C u0 p0 c0 {3,S} {5,D} {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-617.684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.997574,0.0665333,-6.92665e-05,3.65242e-08,-7.72648e-12,-74182.5,33.1574], Tmin=(100,'K'), Tmax=(1136.31,'K')), NASAPolynomial(coeffs=[13.6302,0.0220644,-1.05647e-05,2.08409e-09,-1.49289e-13,-77053.4,-29.4018], Tmin=(1136.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-617.684,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(CsCsFFH) + group(COCsFO) + radical(C=OCOJ) + radical(CCJCO)"""),
)

species(
    label = 'O=C(F)C(O)[CH][C](F)F(3991)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {6,S} {12,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {4,S} {7,S} {8,S} {10,S}
7  C u1 p0 c0 {6,S} {9,S} {11,S}
8  C u0 p0 c0 {1,S} {5,D} {6,S}
9  C u1 p0 c0 {2,S} {3,S} {7,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {4,S}
"""),
    E0 = (-660.876,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.954603,0.0703887,-8.09684e-05,4.80216e-08,-1.14937e-11,-79378.1,33.9007], Tmin=(100,'K'), Tmax=(1007.28,'K')), NASAPolynomial(coeffs=[12.6089,0.0241097,-1.20532e-05,2.41124e-09,-1.73764e-13,-81726,-22.4092], Tmin=(1007.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-660.876,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(CsCsFFH) + group(COCsFO) + radical(CCJCO) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = '[O]C(F)=C([O])CC(F)F(3992)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  O u1 p2 c0 {8,S}
5  O u1 p2 c0 {9,S}
6  C u0 p0 c0 {7,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {12,S}
8  C u0 p0 c0 {4,S} {6,S} {9,D}
9  C u0 p0 c0 {3,S} {5,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-729.359,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.678997,0.0793043,-0.000107518,7.91412e-08,-2.37298e-11,-87607.6,28.5847], Tmin=(100,'K'), Tmax=(809.375,'K')), NASAPolynomial(coeffs=[10.7716,0.0294275,-1.50855e-05,3.00934e-09,-2.15074e-13,-89241.4,-17.9723], Tmin=(809.375,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-729.359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFFH) + group(Cds-CdsCsOs) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = 'O=[C]C(C[C](F)F)OF(3993)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {6,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {4,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {6,S} {8,S} {11,S} {12,S}
8  C u1 p0 c0 {1,S} {2,S} {7,S}
9  C u1 p0 c0 {5,D} {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-312.316,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,190,488,555,1236,1407,1855,455,950,350.096,351.691,351.858],'cm^-1')),
        HinderedRotor(inertia=(0.0019321,'amu*angstrom^2'), symmetry=1, barrier=(21.9372,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130634,'amu*angstrom^2'), symmetry=1, barrier=(11.483,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131234,'amu*angstrom^2'), symmetry=1, barrier=(11.486,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.316552,'amu*angstrom^2'), symmetry=1, barrier=(27.6012,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.565757,0.0808113,-0.000108026,7.47881e-08,-2.07237e-11,-37443.9,32.8698], Tmin=(100,'K'), Tmax=(879.18,'K')), NASAPolynomial(coeffs=[12.8398,0.0249698,-1.2756e-05,2.54856e-09,-1.8265e-13,-39602.2,-24.7654], Tmin=(879.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-312.316,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(CsCsFFH) + group(Cds-OdCsH) + radical(Csj(Cs-CsHH)(F1s)(F1s)) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C([C]=O)CC(F)(F)F(3994)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u1 p2 c0 {7,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {7,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {4,S} {6,S} {9,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
9  C u1 p0 c0 {5,D} {7,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-650.588,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (140.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.850827,0.0705732,-8.13898e-05,4.79943e-08,-1.12693e-11,-78135.3,31.6961], Tmin=(100,'K'), Tmax=(1034.58,'K')), NASAPolynomial(coeffs=[13.6054,0.0212607,-9.89395e-06,1.92385e-09,-1.36737e-13,-80774.4,-30.2708], Tmin=(1034.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-650.588,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCsFFF) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCCJ=O)"""),
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
    E0 = (-210.334,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (179.303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (35.2673,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-58.6427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (120.978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (268.563,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-202.049,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-146.933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-146.933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (92.9072,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-156.841,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-88.1034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-120.358,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-30.466,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-140.459,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-107.133,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-139.109,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-141.051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-28.4405,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (83.1859,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (40.1519,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (31.1231,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-56.4254,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-29.8344,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-135.053,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-87.9766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (198.095,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (15.3644,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(C[C](F)F)C(=O)F(3712)'],
    products = ['O=CC(=O)F(557)', 'CH2CF2(56)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(12)', '[O]C(F)C[C](F)F(799)'],
    products = ['[O]C(C[C](F)F)C(=O)F(3712)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.0026956,'m^3/(mol*s)'), n=2.93313, Ea=(380.159,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_2Br1sCl1sF1sH->F1s',), comment="""Estimated from node Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_2Br1sCl1sF1sH->F1s"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]C(C[C](F)F)C(=O)F(3712)'],
    products = ['[CH2]C(F)(F)C([O])C(=O)F(3708)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]C(C[C](F)F)C(=O)F(3712)'],
    products = ['O=C[C](F)OC[C](F)F(3715)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(151.691,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(6)', 'O=C(F)[CH]C[C](F)F(2720)'],
    products = ['[O]C(C[C](F)F)C(=O)F(3712)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/OneDeC;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['CF2(42)', '[CH2]C([O])C(=O)F(3980)'],
    products = ['[O]C(C[C](F)F)C(=O)F(3712)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C(C[C](F)F)C(=O)F(3712)'],
    products = ['O=C(F)C1CC(F)(F)O1(3718)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;O_rad;Cpri_rad_out_noH]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C(C[C](F)F)C(=O)F(3712)'],
    products = ['O=C(F)C(=O)CC(F)F(3981)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C(C[C](F)F)C(=O)F(3712)'],
    products = ['O=C(F)C(O)C=C(F)F(3982)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C(C[C](F)F)C(=O)F(3712)'],
    products = ['F[C](F)CC1OO[C]1F(3983)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(303.241,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]C(C[C](F)F)C(=O)F(3712)'],
    products = ['[O]C1CC(F)(F)O[C]1F(3984)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.47116e+08,'s^-1'), n=0.669085, Ea=(53.4924,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs] for rate rule [R5_SS_CO;carbonyl_intra;radadd_intra_cs]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C(C[C](F)F)C(=O)F(3712)'],
    products = ['[O]C1(F)OC1C[C](F)F(3985)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.785e+11,'s^-1'), n=0.342, Ea=(122.23,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonylbond_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic
Ea raised from 121.9 to 122.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]C(C[C](F)F)C(=O)F(3712)'],
    products = ['[O]C1CC(F)(F)C1([O])F(3986)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(902977,'s^-1'), n=1.63829, Ea=(89.9757,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs] for rate rule [R5_SS_CO;carbonylbond_intra;radadd_intra_cs]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic
Ea raised from 88.5 to 90.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction14',
    reactants = ['O=C(C[C](F)F)[C](O)F(3987)'],
    products = ['[O]C(C[C](F)F)C(=O)F(3712)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(264.716,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_Ext-1R!H-R',), comment="""Estimated from node Root_3R!H->C_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O=CC(=O)F(557)', '[CH2][C](F)F(185)'],
    products = ['[O]C(C[C](F)F)C(=O)F(3712)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(44.9942,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['CFO(50)', 'O=CC[C](F)F(911)'],
    products = ['[O]C(C[C](F)F)C(=O)F(3712)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=-1.07934e-11, Ea=(88.1675,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(5)', 'O=C(F)C(=O)C[C](F)F(3988)'],
    products = ['[O]C(C[C](F)F)C(=O)F(3712)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(39.7,'m^3/(mol*s)'), n=1.88, Ea=(26.2085,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O][CH]C(=O)F(438)', 'CH2CF2(56)'],
    products = ['[O]C(C[C](F)F)C(=O)F(3712)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(28.3314,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(5)', '[O]C(C=C(F)F)C(=O)F(3822)'],
    products = ['[O]C(C[C](F)F)C(=O)F(3712)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(306600,'m^3/(mol*s)'), n=0.481, Ea=(26.4631,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_4R!H->C_Sp-4C-1COS_N-6R!H-inRing_Ext-4C-R_Ext-4C-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_4R!H->C_Sp-4C-1COS_N-6R!H-inRing_Ext-4C-R_Ext-4C-R"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O][CH]C(=O)F(438)', '[CH2][C](F)F(185)'],
    products = ['[O]C(C[C](F)F)C(=O)F(3712)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction21',
    reactants = ['HF(38)', 'O=[C]C(=O)C[C](F)F(3973)'],
    products = ['[O]C(C[C](F)F)C(=O)F(3712)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(274.931,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CF2(87)', '[CH2]C([O])C(=O)F(3980)'],
    products = ['[O]C(C[C](F)F)C(=O)F(3712)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[O]C(C[C](F)F)C(=O)F(3712)'],
    products = ['O=C(F)[C](O)C[C](F)F(3989)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.70223e+09,'s^-1'), n=1.15155, Ea=(153.908,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_OneDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_CO]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]C(C[C](F)F)C(=O)F(3712)'],
    products = ['[O]C([CH]C(F)F)C(=O)F(3990)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.76836e+09,'s^-1'), n=1.1815, Ea=(180.499,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_noH;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O]C(C[C](F)F)C(=O)F(3712)'],
    products = ['O=C(F)C(O)[CH][C](F)F(3991)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(223829,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O]C(C[C](F)F)C(=O)F(3712)'],
    products = ['[O]C(F)=C([O])CC(F)F(3992)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.50974e+07,'s^-1'), n=1.33047, Ea=(122.357,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_noH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['O=[C]C(C[C](F)F)OF(3993)'],
    products = ['[O]C(C[C](F)F)C(=O)F(3712)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(103.7,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[O]C(C[C](F)F)C(=O)F(3712)'],
    products = ['[O]C([C]=O)CC(F)(F)F(3994)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(0.00363316,'s^-1'), n=4.43046, Ea=(225.698,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R4F_Ext-3R!H-R',), comment="""Estimated from node R4F_Ext-3R!H-R"""),
)

network(
    label = 'PDepNetwork #1337',
    isomers = [
        '[O]C(C[C](F)F)C(=O)F(3712)',
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
    label = 'PDepNetwork #1337',
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

