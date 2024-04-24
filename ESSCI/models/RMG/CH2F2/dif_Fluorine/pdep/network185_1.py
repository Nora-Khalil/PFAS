species(
    label = 'O=C(F)OC(O)F(570)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {6,S} {7,S}
4 O u0 p2 c0 {6,S} {9,S}
5 O u0 p2 c0 {7,D}
6 C u0 p0 c0 {1,S} {3,S} {4,S} {8,S}
7 C u0 p0 c0 {2,S} {3,S} {5,D}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (-1020.47,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,509,613,660,1171,1360,1414,3084,482,664,788,1296,1923,221.051,225.396,1066.89],'cm^-1')),
        HinderedRotor(inertia=(0.268854,'amu*angstrom^2'), symmetry=1, barrier=(9.36951,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.809722,'amu*angstrom^2'), symmetry=1, barrier=(28.8895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.605702,'amu*angstrom^2'), symmetry=1, barrier=(21.5102,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71583,0.0500006,-5.30385e-05,2.7775e-08,-5.77562e-12,-122651,21.0141], Tmin=(100,'K'), Tmax=(1161.05,'K')), NASAPolynomial(coeffs=[12.0233,0.0144896,-7.16047e-06,1.43206e-09,-1.03373e-13,-125045,-30.2527], Tmin=(1161.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1020.47,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(CsFHOO) + group(COFOO)"""),
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
    label = 'CO2(13)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,D}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,D} {2,D}
"""),
    E0 = (-403.138,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([664.558,664.681,1368.23,2264.78],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0094,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(2028.74,'J/mol'), sigma=(3.763,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(2.65,'angstroms^3'), rotrelaxcollnum=2.1, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.35681,0.00898413,-7.12206e-06,2.4573e-09,-1.42885e-13,-48372,9.9009], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.63651,0.00274146,-9.95898e-07,1.60387e-10,-9.16199e-15,-49024.9,-1.9349], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-403.138,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), label="""CO2""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'CHFO(46)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 H u0 p0 c0 {3,S}
"""),
    E0 = (-386.87,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.1185,0.010774,-5.51582e-06,5.81519e-10,2.17042e-13,-46356.7,14.6524], Tmin=(298,'K'), Tmax=(1150,'K')), NASAPolynomial(coeffs=[4.81925,0.00505418,-2.01421e-06,3.63591e-10,-2.44724e-14,-47263.1,0.0972251], Tmin=(1150,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-386.87,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(83.1447,'J/mol/K'), label="""CHFO""", comment="""Thermo library: Fluorine"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.56838,-0.000852125,2.48918e-06,-1.56331e-09,3.13595e-13,-14284.3,3.57912], Tmin=(100,'K'), Tmax=(1571.64,'K')), NASAPolynomial(coeffs=[2.91307,0.00164658,-6.88615e-07,1.21037e-10,-7.84018e-15,-14180.9,6.71045], Tmin=(1571.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-118.741,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""CO""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'OC(F)OF(505)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {5,S} {7,S}
4 O u0 p2 c0 {2,S} {5,S}
5 C u0 p0 c0 {1,S} {3,S} {4,S} {6,S}
6 H u0 p0 c0 {5,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (-514.434,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,509,613,660,1171,1360,1414,3084,180],'cm^-1')),
        HinderedRotor(inertia=(0.884565,'amu*angstrom^2'), symmetry=1, barrier=(20.3379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.883872,'amu*angstrom^2'), symmetry=1, barrier=(20.3219,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0222,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.85504,0.00891939,0.00010045,-2.42206e-07,1.6199e-10,-61863,9.86619], Tmin=(10,'K'), Tmax=(542.673,'K')), NASAPolynomial(coeffs=[6.62744,0.019591,-1.50291e-05,5.28637e-09,-6.8652e-13,-62622,-6.03462], Tmin=(542.673,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-514.434,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), label="""OC(F)OF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]C(F)O[C](O)F(612)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {6,S} {7,S}
4 O u0 p2 c0 {7,S} {9,S}
5 O u1 p2 c0 {6,S}
6 C u0 p0 c0 {1,S} {3,S} {5,S} {8,S}
7 C u1 p0 c0 {2,S} {3,S} {4,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (-593.5,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,451,553,637,1069,1180,1265,1301,3056,482,586,761,1411,180,180,937.787],'cm^-1')),
        HinderedRotor(inertia=(0.69705,'amu*angstrom^2'), symmetry=1, barrier=(16.0265,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115694,'amu*angstrom^2'), symmetry=1, barrier=(2.66002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.58218,'amu*angstrom^2'), symmetry=1, barrier=(36.3774,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49865,0.0520881,-5.60822e-05,2.93266e-08,-5.98402e-12,-71289,24.5693], Tmin=(100,'K'), Tmax=(1197.09,'K')), NASAPolynomial(coeffs=[13.6369,0.0115287,-5.25952e-06,1.02299e-09,-7.3064e-14,-74195.1,-36.1741], Tmin=(1197.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-593.5,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(CsFHOO) + group(CsFHOO) + radical(O2sj(Cs-F1sO2sH)) + radical(Csj(F1s)(O2s-Cs)(O2s-H))"""),
)

species(
    label = 'O[C](F)O[C](O)F(613)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {6,S} {7,S}
4 O u0 p2 c0 {6,S} {9,S}
5 O u0 p2 c0 {7,S} {8,S}
6 C u1 p0 c0 {1,S} {3,S} {4,S}
7 C u1 p0 c0 {2,S} {3,S} {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (-644.2,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,428,536,497,675,664,858,1354,1468,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.123224,'amu*angstrom^2'), symmetry=1, barrier=(2.83317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.56166,'amu*angstrom^2'), symmetry=1, barrier=(35.9055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0722142,'amu*angstrom^2'), symmetry=1, barrier=(20.895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.072666,'amu*angstrom^2'), symmetry=1, barrier=(20.9069,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15203,0.0574332,-6.55017e-05,3.51406e-08,-7.22206e-12,-77372.4,25.4561], Tmin=(100,'K'), Tmax=(1201.99,'K')), NASAPolynomial(coeffs=[16.1312,0.00758541,-3.29505e-06,6.38484e-10,-4.60247e-14,-80973.3,-49.5655], Tmin=(1201.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-644.2,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(191.233,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(CsFHOO) + group(CsFHOO) + radical(Csj(F1s)(O2s-Cs)(O2s-H)) + radical(Csj(F1s)(O2s-Cs)(O2s-H))"""),
)

species(
    label = '[O]C(F)OC([O])F(614)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {6,S} {7,S}
4 O u1 p2 c0 {6,S}
5 O u1 p2 c0 {7,S}
6 C u0 p0 c0 {1,S} {3,S} {4,S} {8,S}
7 C u0 p0 c0 {2,S} {3,S} {5,S} {9,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-542.8,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([421,481,502,604,564,710,1025,1113,1111,1249,1195,1335,1232,1370,3006,3106,180,306.222,2404.16],'cm^-1')),
        HinderedRotor(inertia=(0.530355,'amu*angstrom^2'), symmetry=1, barrier=(12.1939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.56241,'amu*angstrom^2'), symmetry=1, barrier=(35.9228,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.8439,0.0467562,-4.66981e-05,2.35465e-08,-4.7566e-12,-65205.6,22.3015], Tmin=(100,'K'), Tmax=(1188.57,'K')), NASAPolynomial(coeffs=[11.1362,0.0154838,-7.23105e-06,1.40922e-09,-1.0025e-13,-67414.5,-24.1335], Tmin=(1188.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-542.8,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(CsFHOO) + group(CsFHOO) + radical(O2sj(Cs-F1sO2sH)) + radical(O2sj(Cs-F1sO2sH))"""),
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
    label = 'O=C(F)O[CH]O(615)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {5,S} {6,S}
3 O u0 p2 c0 {5,S} {8,S}
4 O u0 p2 c0 {6,D}
5 C u1 p0 c0 {2,S} {3,S} {7,S}
6 C u0 p0 c0 {1,S} {2,S} {4,D}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (-595.669,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,482,664,788,1296,1923,448.271,756.52,758.831],'cm^-1')),
        HinderedRotor(inertia=(0.162635,'amu*angstrom^2'), symmetry=1, barrier=(3.73931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.701108,'amu*angstrom^2'), symmetry=1, barrier=(16.1199,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04421,'amu*angstrom^2'), symmetry=1, barrier=(24.0083,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.0338,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81287,0.0464932,-5.47881e-05,3.13845e-08,-6.97396e-12,-71562.3,19.9356], Tmin=(100,'K'), Tmax=(1106.44,'K')), NASAPolynomial(coeffs=[12.1838,0.00899982,-3.9579e-06,7.5739e-10,-5.37144e-14,-73857.2,-31.1472], Tmin=(1106.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-595.669,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(Cs-OsOsHH) + group(COFOO) + radical(OCJO)"""),
)

species(
    label = 'O=[C]OC(O)F(616)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {5,S} {6,S}
3 O u0 p2 c0 {5,S} {8,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
6 C u1 p0 c0 {2,S} {4,D}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (-601.414,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,509,613,660,1171,1360,1414,3084,1855,455,950,180,289.42],'cm^-1')),
        HinderedRotor(inertia=(0.631935,'amu*angstrom^2'), symmetry=1, barrier=(14.5294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.63795,'amu*angstrom^2'), symmetry=1, barrier=(60.6516,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162094,'amu*angstrom^2'), symmetry=1, barrier=(38.6475,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.0338,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.08758,0.0435144,-4.86055e-05,2.79276e-08,-6.46486e-12,-72265.7,17.6225], Tmin=(100,'K'), Tmax=(1041.31,'K')), NASAPolynomial(coeffs=[9.56882,0.0147766,-7.20898e-06,1.42478e-09,-1.02019e-13,-73823.8,-18.773], Tmin=(1041.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-601.414,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(CsFHOO) + group(Cds-OdOsH) + radical((O)CJOC)"""),
)

species(
    label = '[O]C(=O)F(143)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 O u1 p2 c0 {4,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
"""),
    E0 = (-381.43,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(62.9882,'amu')),
        NonlinearRotor(inertia=([36.4335,44.7291,81.1626],'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([509.533,537.43,765.261,1011.27,1202.3,1550.29],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (63.0078,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3607.92,'J/mol'), sigma=(5.19854,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=563.55 K, Pc=58.27 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05952,-0.0056915,7.13967e-05,-1.29651e-07,7.44832e-11,-45874.2,8.21459], Tmin=(10,'K'), Tmax=(575.482,'K')), NASAPolynomial(coeffs=[3.42951,0.0118543,-8.65614e-06,2.84335e-09,-3.46285e-13,-46019.7,9.01155], Tmin=(575.482,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-381.43,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""[O]C(DO)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O[CH]F(180)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {3,S} {5,S}
3 C u1 p0 c0 {1,S} {2,S} {4,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (-244.274,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,580,1155,1237,1373,3147],'cm^-1')),
        HinderedRotor(inertia=(0.780979,'amu*angstrom^2'), symmetry=1, barrier=(17.9562,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (49.0244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.95258,0.00277528,4.2927e-05,-8.89919e-08,5.36995e-11,-29376.8,7.35238], Tmin=(10,'K'), Tmax=(577.676,'K')), NASAPolynomial(coeffs=[4.36429,0.0110503,-7.44974e-06,2.48543e-09,-3.1757e-13,-29610,3.98525], Tmin=(577.676,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-244.274,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""O[CH]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CFO(50)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u1 p0 c0 {1,S} {2,D}
"""),
    E0 = (-181.965,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(46.9933,'amu')),
        NonlinearRotor(inertia=([2.65864,43.9175,46.5761],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([636.046,1085.22,1918.24],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (47.0084,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.94061,0.00747092,-6.21833e-06,2.62916e-09,-4.73222e-13,-21792.3,11.1689], Tmin=(298,'K'), Tmax=(1250,'K')), NASAPolynomial(coeffs=[5.14996,0.00184938,-7.40232e-07,1.32785e-10,-8.8306e-15,-22488.4,-0.49683], Tmin=(1250,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-181.964,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CFO""", comment="""Thermo library: Fluorine"""),
)

species(
    label = '[O]C(O)F(415)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {4,S} {6,S}
3 O u1 p2 c0 {4,S}
4 C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
5 H u0 p0 c0 {4,S}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (-401.555,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,451,553,637,1069,1180,1265,1301,3056],'cm^-1')),
        HinderedRotor(inertia=(0.412376,'amu*angstrom^2'), symmetry=1, barrier=(9.48133,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (65.0238,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.91831,0.00720244,6.35725e-05,-2.08276e-07,2.06569e-10,-48296.9,9.02236], Tmin=(10,'K'), Tmax=(327.664,'K')), NASAPolynomial(coeffs=[3.62426,0.0171767,-1.13157e-05,3.55881e-09,-4.26934e-13,-48311.9,9.58992], Tmin=(327.664,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-401.555,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""[O]C(O)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'O=C(F)O[CH]F(472)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {5,S} {6,S}
4 O u0 p2 c0 {6,D}
5 C u1 p0 c0 {1,S} {3,S} {7,S}
6 C u0 p0 c0 {2,S} {3,S} {4,D}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-604.807,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([580,1155,1237,1373,3147,482,664,788,1296,1923,180,492.863,647.158],'cm^-1')),
        HinderedRotor(inertia=(0.740594,'amu*angstrom^2'), symmetry=1, barrier=(17.0277,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.58964,'amu*angstrom^2'), symmetry=1, barrier=(36.549,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.0249,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.75173,0.0249972,-2.08089e-06,-2.63495e-08,1.65145e-11,-72740.3,11.2485], Tmin=(10,'K'), Tmax=(755.828,'K')), NASAPolynomial(coeffs=[7.32436,0.018139,-1.23821e-05,3.82767e-09,-4.43211e-13,-73624.5,-7.26399], Tmin=(755.828,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-604.807,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), label="""ODC(F)O[CH]F""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]C(F)OC(=O)F(617)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {6,S} {7,S}
4 O u1 p2 c0 {6,S}
5 O u0 p2 c0 {7,D}
6 C u0 p0 c0 {1,S} {3,S} {4,S} {8,S}
7 C u0 p0 c0 {2,S} {3,S} {5,D}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-764.934,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([451,553,637,1069,1180,1265,1301,3056,482,664,788,1296,1923,180,580.509,2104.54],'cm^-1')),
        HinderedRotor(inertia=(1.30332,'amu*angstrom^2'), symmetry=1, barrier=(29.9658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30324,'amu*angstrom^2'), symmetry=1, barrier=(29.9641,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.024,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91906,0.0405969,-3.6544e-05,1.56493e-08,-2.62167e-12,-91921.2,21.8246], Tmin=(100,'K'), Tmax=(1440.43,'K')), NASAPolynomial(coeffs=[12.9256,0.0100323,-4.7156e-06,9.1841e-10,-6.49979e-14,-95092,-35.2926], Tmin=(1440.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-764.934,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(CsFHOO) + group(COFOO) + radical(O2sj(Cs-F1sO2sH))"""),
)

species(
    label = 'O=C(F)O[C](O)F(389)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {6,S} {7,S}
4 O u0 p2 c0 {6,S} {8,S}
5 O u0 p2 c0 {7,D}
6 C u1 p0 c0 {1,S} {3,S} {4,S}
7 C u0 p0 c0 {2,S} {3,S} {5,D}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-814.278,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,482,586,761,1411,482,664,788,1296,1923,191.359,1054.98,1905.36],'cm^-1')),
        HinderedRotor(inertia=(0.591423,'amu*angstrom^2'), symmetry=1, barrier=(14.7891,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11416,'amu*angstrom^2'), symmetry=1, barrier=(27.4996,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.55459,'amu*angstrom^2'), symmetry=1, barrier=(40.1069,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.024,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00527,0.0434887,-4.46043e-05,2.23764e-08,-4.46435e-12,-97863,22.1567], Tmin=(100,'K'), Tmax=(1205,'K')), NASAPolynomial(coeffs=[11.2161,0.012913,-6.54266e-06,1.31844e-09,-9.54184e-14,-100083,-23.9976], Tmin=(1205,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-814.278,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(CsFHOO) + group(COFOO) + radical(CsF1sO2sO2s)"""),
)

species(
    label = 'O=C(O)F(141)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (-625.647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (64.0158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.96655,0.00185137,4.47909e-05,-7.96564e-08,4.23622e-11,-75247.4,7.80392], Tmin=(10,'K'), Tmax=(614.343,'K')), NASAPolynomial(coeffs=[2.83475,0.0173246,-1.27762e-05,4.2862e-09,-5.35304e-13,-75261.2,11.4681], Tmin=(614.343,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-625.647,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""ODC(O)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=COC(=O)F(337)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {5,S} {6,S}
3 O u0 p2 c0 {5,D}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {2,S} {3,D} {7,S}
6 C u0 p0 c0 {1,S} {2,S} {4,D}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-707.29,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,482,664,788,1296,1923,249.355,4000],'cm^-1')),
        HinderedRotor(inertia=(0.299229,'amu*angstrom^2'), symmetry=1, barrier=(13.2007,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.299283,'amu*angstrom^2'), symmetry=1, barrier=(13.2008,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (92.0259,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.52888,0.0335236,-3.63657e-05,1.95715e-08,-4.22042e-12,-85015.4,19.1008], Tmin=(100,'K'), Tmax=(1113.61,'K')), NASAPolynomial(coeffs=[8.89044,0.010673,-5.58599e-06,1.14479e-09,-8.36428e-14,-86432.2,-12.2745], Tmin=(1113.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-707.29,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-O2d)) + group(Cds-OdOsH) + group(COFOO)"""),
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
    E0 = (-525.674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-142.588,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-131.414,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-237.146,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-119.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-124.092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-129.837,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-227.019,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-184.835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-177.727,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-154.444,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-203.788,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-463.795,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-436.577,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=C(F)OC(O)F(570)'],
    products = ['HF(38)', 'CO2(13)', 'CHFO(46)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5.65824e+14,'s^-1'), n=-0.436116, Ea=(96.1059,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.9927400314419011, var=19.60363129417426, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R!H->C_5Br1sCl1sF1sH->F1s_Ext-2C-R_N-7R!H->C',), comment="""Estimated from node Root_N-1R!H->C_5Br1sCl1sF1sH->F1s_Ext-2C-R_N-7R!H->C"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(12)', 'OC(F)OF(505)'],
    products = ['O=C(F)OC(O)F(570)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.000742267,'m^3/(mol*s)'), n=2.83796, Ea=(91.9008,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.2632766423807829, var=145.9379134136138, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-1COCbCdCsCtHNOSSidSis->Cs',), comment="""Estimated from node Root_N-1COCbCdCsCtHNOSSidSis->Cs"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]C(F)O[C](O)F(612)'],
    products = ['O=C(F)OC(O)F(570)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O[C](F)O[C](O)F(613)'],
    products = ['O=C(F)OC(O)F(570)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]C(F)OC([O])F(614)'],
    products = ['O=C(F)OC(O)F(570)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['F(37)', 'O=C(F)O[CH]O(615)'],
    products = ['O=C(F)OC(O)F(570)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.52887e+10,'m^3/(mol*s)'), n=-0.421056, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.030895812897821735, var=3.1393582276389975, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F(37)', 'O=[C]OC(O)F(616)'],
    products = ['O=C(F)OC(O)F(570)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.52887e+10,'m^3/(mol*s)'), n=-0.421056, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.030895812897821735, var=3.1393582276389975, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C(=O)F(143)', 'O[CH]F(180)'],
    products = ['O=C(F)OC(O)F(570)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.47663e+07,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CFO(50)', '[O]C(O)F(415)'],
    products = ['O=C(F)OC(O)F(570)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction10',
    reactants = ['OH(4)', 'O=C(F)O[CH]F(472)'],
    products = ['O=C(F)OC(O)F(570)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.7e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(5)', '[O]C(F)OC(=O)F(617)'],
    products = ['O=C(F)OC(O)F(570)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(5)', 'O=C(F)O[C](O)F(389)'],
    products = ['O=C(F)OC(O)F(570)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.1766,'m^3/(mol*s)'), n=1.94174, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_N-4R!H->C_Ext-2C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_N-4R!H->C_Ext-2C-R
Ea raised from -1.4 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['O=C(F)OC(O)F(570)'],
    products = ['CHFO(46)', 'O=C(O)F(141)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.82035e+06,'s^-1'), n=1.67734, Ea=(157.985,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.2285837836968013, var=2.96975919141528, Tref=1000.0, N=18, data_mean=0.0, correlation='Root_N-1R!H->C_N-5R!H->N_Ext-2R!H-R',), comment="""Estimated from node Root_N-1R!H->C_N-5R!H->N_Ext-2R!H-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['HF(38)', 'O=COC(=O)F(337)'],
    products = ['O=C(F)OC(O)F(570)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.109156,'m^3/(mol*s)'), n=1.86531, Ea=(153.141,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd',), comment="""Estimated from node HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd"""),
)

network(
    label = 'PDepNetwork #185',
    isomers = [
        'O=C(F)OC(O)F(570)',
    ],
    reactants = [
        ('HF(38)', 'CO2(13)', 'CHFO(46)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #185',
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

