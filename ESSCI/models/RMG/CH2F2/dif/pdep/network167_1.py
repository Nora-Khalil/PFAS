species(
    label = 'OC(F)OC(F)F(352)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {6,S} {7,S}
5  O u0 p2 c0 {6,S} {10,S}
6  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
7  C u0 p0 c0 {2,S} {3,S} {4,S} {9,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-1084.91,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,509,613,660,1171,1360,1414,3084,519,593,830,1115,1166,1389,1413,3114,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.480687,'amu*angstrom^2'), symmetry=1, barrier=(11.0519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.481017,'amu*angstrom^2'), symmetry=1, barrier=(11.0595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.84579,'amu*angstrom^2'), symmetry=1, barrier=(19.4464,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.70499,0.0220185,0.000153714,-4.94209e-07,4.27634e-10,-130484,12.9624], Tmin=(10,'K'), Tmax=(423.094,'K')), NASAPolynomial(coeffs=[6.89189,0.0318825,-2.30467e-05,7.72902e-09,-9.68067e-13,-131112,-3.90157], Tmin=(423.094,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1084.91,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), label="""OC(F)OC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'CF2(42)',
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
    label = 'OC(O)F(608)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {4,S} {6,S}
3 O u0 p2 c0 {4,S} {7,S}
4 C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
5 H u0 p0 c0 {4,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (-655.938,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,509,613,660,1171,1360,1414,3084],'cm^-1')),
        HinderedRotor(inertia=(0.75594,'amu*angstrom^2'), symmetry=1, barrier=(17.3806,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.757577,'amu*angstrom^2'), symmetry=1, barrier=(17.4182,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (66.0318,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.89583,0.00622994,8.13087e-05,-1.81346e-07,1.14908e-10,-78884.9,8.27912], Tmin=(10,'K'), Tmax=(563.154,'K')), NASAPolynomial(coeffs=[5.5716,0.0179261,-1.27023e-05,4.35646e-09,-5.63739e-13,-79447.9,-2.16576], Tmin=(563.154,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-655.938,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), label="""OC(O)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CHF(40)',
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
    label = 'OC(F)OF(468)',
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
    label = 'O[CH]OC(F)F(609)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {5,S} {6,S}
4 O u0 p2 c0 {6,S} {9,S}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
6 C u1 p0 c0 {3,S} {4,S} {8,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (-658.714,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,519,593,830,1115,1166,1389,1413,3114,3025,407.5,1350,352.5,180,180,1066.88],'cm^-1')),
        HinderedRotor(inertia=(1.1037,'amu*angstrom^2'), symmetry=1, barrier=(25.3763,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109568,'amu*angstrom^2'), symmetry=1, barrier=(2.51919,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.10282,'amu*angstrom^2'), symmetry=1, barrier=(25.3559,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.0408,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.68464,0.0356052,-2.37181e-05,2.31672e-09,2.31171e-12,-79226,11.5802], Tmin=(10,'K'), Tmax=(968.036,'K')), NASAPolynomial(coeffs=[9.88261,0.0188063,-1.13417e-05,3.19646e-09,-3.43897e-13,-80838.9,-20.2525], Tmin=(968.036,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-658.714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), label="""O[CH]OC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OC(F)O[CH]F(610)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {5,S} {6,S}
4 O u0 p2 c0 {5,S} {9,S}
5 C u0 p0 c0 {1,S} {3,S} {4,S} {7,S}
6 C u1 p0 c0 {2,S} {3,S} {8,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (-647.333,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,509,613,660,1171,1360,1414,3084,580,1155,1237,1373,3147,180,180,1249.15],'cm^-1')),
        HinderedRotor(inertia=(0.764793,'amu*angstrom^2'), symmetry=1, barrier=(17.5841,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.10323,'amu*angstrom^2'), symmetry=1, barrier=(25.3654,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.0408,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.72124,0.0316456,-1.0776e-05,-1.25521e-08,8.2614e-12,-77857,11.0532], Tmin=(10,'K'), Tmax=(887.994,'K')), NASAPolynomial(coeffs=[8.44216,0.0217748,-1.33504e-05,3.83124e-09,-4.19414e-13,-79144.7,-13.6914], Tmin=(887.994,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-647.333,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""OC(F)O[CH]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]C(F)F(170)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u1 p2 c0 {4,S}
4 C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
5 H u0 p0 c0 {4,S}
"""),
    E0 = (-426.241,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(66.9995,'amu')),
        NonlinearRotor(inertia=([47.2408,48.1768,88.2639],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([480.226,502.831,615.399,942.972,1118.44,1159.86,1274.53,1324.03,2842.12],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.0148,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.07235,-0.00674902,8.67846e-05,-1.53143e-07,8.6196e-11,-51263.6,9.19392], Tmin=(10,'K'), Tmax=(577.889,'K')), NASAPolynomial(coeffs=[2.80054,0.016572,-1.1432e-05,3.6345e-09,-4.33761e-13,-51359,12.5348], Tmin=(577.889,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-426.241,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""[O]C(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O[CH]F(290)',
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
    label = 'CHF2(81)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u1 p0 c0 {1,S} {2,S} {4,S}
4 H u0 p0 c0 {3,S}
"""),
    E0 = (-256.71,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(51.0046,'amu')),
        NonlinearRotor(inertia=([7.43413,45.9439,52.5803],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([549.125,1005.77,1195.1,1212.61,1359.42,3085.19],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (51.0154,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2178.39,'J/mol'), sigma=(4.123,'angstroms'), dipoleMoment=(1.8,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05476,-0.0040567,3.90133e-05,-5.51349e-08,2.50461e-11,-30875.2,7.58714], Tmin=(10,'K'), Tmax=(697.139,'K')), NASAPolynomial(coeffs=[2.58942,0.0108145,-6.89144e-06,2.06262e-09,-2.34597e-13,-30827.9,13.0014], Tmin=(697.139,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-256.71,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""F[CH]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]C(O)F(296)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51457,2.92734e-05,-5.32151e-07,1.01948e-09,-3.85939e-13,3414.25,2.10435], Tmin=(100,'K'), Tmax=(1145.76,'K')), NASAPolynomial(coeffs=[3.07194,0.00060402,-1.39806e-08,-2.13441e-11,2.48061e-15,3579.39,4.57801], Tmin=(1145.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.3945,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH(D)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'F[CH]OC(F)F(299)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {6,S}
4 O u0 p2 c0 {5,S} {6,S}
5 C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
6 C u1 p0 c0 {3,S} {4,S} {8,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-670.593,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([519,593,830,1115,1166,1389,1413,3114,580,1155,1237,1373,3147,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.35387,'amu*angstrom^2'), symmetry=1, barrier=(31.1282,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00193046,'amu*angstrom^2'), symmetry=1, barrier=(10.7431,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.0318,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.5602,0.0466294,-0.000151773,3.47811e-07,-3.06675e-10,-80653.3,12.3867], Tmin=(10,'K'), Tmax=(351.878,'K')), NASAPolynomial(coeffs=[4.39046,0.0262196,-1.79984e-05,5.74934e-09,-6.92539e-13,-80643.8,10.2136], Tmin=(351.878,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-670.593,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""F[CH]OC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,1.49298e-14,-2.05327e-17,9.23927e-21,-1.28064e-24,25474.2,-0.444973], Tmin=(100,'K'), Tmax=(3959.16,'K')), NASAPolynomial(coeffs=[2.5,1.19009e-10,-4.23987e-14,6.68965e-18,-3.94352e-22,25474.2,-0.444972], Tmin=(3959.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.805,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""H""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[O]C(F)OC(F)F(530)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {7,S}
4 O u0 p2 c0 {6,S} {7,S}
5 O u1 p2 c0 {7,S}
6 C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
7 C u0 p0 c0 {3,S} {4,S} {5,S} {9,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-832.392,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([519,593,830,1115,1166,1389,1413,3114,451,553,637,1069,1180,1265,1301,3056,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.312684,'amu*angstrom^2'), symmetry=1, barrier=(7.18923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.724027,'amu*angstrom^2'), symmetry=1, barrier=(16.6468,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.59814,0.0349147,4.09671e-06,-7.24873e-08,5.72024e-11,-100110,13.1403], Tmin=(10,'K'), Tmax=(558.625,'K')), NASAPolynomial(coeffs=[6.82145,0.0277041,-1.91547e-05,6.11581e-09,-7.33603e-13,-100718,-2.74978], Tmin=(558.625,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-832.392,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""[O]C(F)OC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O[C](F)OC(F)F(601)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {7,S}
4 O u0 p2 c0 {6,S} {7,S}
5 O u0 p2 c0 {7,S} {9,S}
6 C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
7 C u1 p0 c0 {3,S} {4,S} {5,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-878.876,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,519,593,830,1115,1166,1389,1413,3114,482,586,761,1411,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.3577,'amu*angstrom^2'), symmetry=1, barrier=(31.2162,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.506329,'amu*angstrom^2'), symmetry=1, barrier=(11.6415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.506285,'amu*angstrom^2'), symmetry=1, barrier=(11.6405,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.4959,0.0441543,-2.72259e-05,-2.74232e-08,3.17005e-11,-105700,12.7297], Tmin=(10,'K'), Tmax=(594.27,'K')), NASAPolynomial(coeffs=[8.65821,0.0248921,-1.76919e-05,5.72895e-09,-6.92084e-13,-106587,-11.7895], Tmin=(594.27,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-878.876,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), label="""O[C](F)OC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OC(F)O[C](F)F(604)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {7,S}
4 O u0 p2 c0 {6,S} {7,S}
5 O u0 p2 c0 {6,S} {9,S}
6 C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
7 C u1 p0 c0 {2,S} {3,S} {4,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-874.243,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,509,613,660,1171,1360,1414,3084,493,600,700,1144,1293,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.257643,'amu*angstrom^2'), symmetry=1, barrier=(5.92372,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.258409,'amu*angstrom^2'), symmetry=1, barrier=(5.94133,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.60461,'amu*angstrom^2'), symmetry=1, barrier=(13.9012,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.41486,0.0505296,-6.87282e-05,4.96794e-08,-1.45838e-11,-105143,13.5918], Tmin=(10,'K'), Tmax=(806.995,'K')), NASAPolynomial(coeffs=[9.3964,0.0208812,-1.36192e-05,4.15332e-09,-4.80183e-13,-106108,-13.983], Tmin=(806.995,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-874.243,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), label="""OC(F)O[C](F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=COC(F)F(162)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {5,S} {6,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
6 C u0 p0 c0 {3,S} {4,D} {8,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-806.64,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([519,593,830,1115,1166,1389,1413,3114,2782.5,750,1395,475,1775,1000,253.197,974.364],'cm^-1')),
        HinderedRotor(inertia=(0.819988,'amu*angstrom^2'), symmetry=1, barrier=(37.4431,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.821166,'amu*angstrom^2'), symmetry=1, barrier=(37.4316,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0328,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.92173,0.00488721,0.000104944,-2.18912e-07,1.38545e-10,-97016.1,12.6043], Tmin=(10,'K'), Tmax=(511.764,'K')), NASAPolynomial(coeffs=[2.26171,0.0333089,-2.36362e-05,7.56645e-09,-9.03623e-13,-97048.5,17.5246], Tmin=(511.764,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-806.64,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""ODCOC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    E0 = (-588.937,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-416.012,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (2.64959,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-217.488,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-206.107,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-302.181,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-289.932,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-273.865,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-252.253,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-298.737,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-294.104,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-548.033,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['OC(F)OC(F)F(352)'],
    products = ['HF(38)', 'CHFO(46)', 'CHFO(46)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.13165e+15,'s^-1'), n=-0.436116, Ea=(127.644,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.9927400314419011, var=19.60363129417426, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R!H->C_5Br1sCl1sF1sH->F1s_Ext-2C-R_N-7R!H->C',), comment="""Estimated from node Root_N-1R!H->C_5Br1sCl1sF1sH->F1s_Ext-2C-R_N-7R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CF2(42)', 'OC(O)F(608)'],
    products = ['OC(F)OC(F)F(352)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.5279e-12,'m^3/(mol*s)'), n=5.02686, Ea=(75.3048,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='OH_N-2Br1sCl1sF1sHI1s->H_2Br1sCl1sF1sI1s->F1s',), comment="""Estimated from node OH_N-2Br1sCl1sF1sHI1s->H_2Br1sCl1sF1sI1s->F1s
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CHF(40)', 'OC(F)OF(468)'],
    products = ['OC(F)OC(F)F(352)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.13916e+09,'m^3/(mol*s)'), n=-1.00842, Ea=(9.99325,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17406617270594374, var=16.167420070870673, Tref=1000.0, N=4, data_mean=0.0, correlation='OHOY',), comment="""Estimated from node OHOY"""),
)

reaction(
    label = 'reaction4',
    reactants = ['F(37)', 'O[CH]OC(F)F(609)'],
    products = ['OC(F)OC(F)F(352)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.52887e+10,'m^3/(mol*s)'), n=-0.421056, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.030895812897821735, var=3.1393582276389975, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction5',
    reactants = ['F(37)', 'OC(F)O[CH]F(610)'],
    products = ['OC(F)OC(F)F(352)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]C(F)F(170)', 'O[CH]F(290)'],
    products = ['OC(F)OC(F)F(352)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CHF2(81)', '[O]C(O)F(296)'],
    products = ['OC(F)OC(F)F(352)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction8',
    reactants = ['OH(4)', 'F[CH]OC(F)F(299)'],
    products = ['OC(F)OC(F)F(352)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.7e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(5)', '[O]C(F)OC(F)F(530)'],
    products = ['OC(F)OC(F)F(352)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(5)', 'O[C](F)OC(F)F(601)'],
    products = ['OC(F)OC(F)F(352)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.1766,'m^3/(mol*s)'), n=1.94174, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_N-4R!H->C_Ext-2C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_N-4R!H->C_Ext-2C-R
Ea raised from -1.6 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(5)', 'OC(F)O[C](F)F(604)'],
    products = ['OC(F)OC(F)F(352)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.1766,'m^3/(mol*s)'), n=1.94174, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_N-4R!H->C_Ext-2C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_N-4R!H->C_Ext-2C-R
Ea raised from -2.1 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction12',
    reactants = ['HF(38)', 'O=COC(F)F(162)'],
    products = ['OC(F)OC(F)F(352)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(0.109156,'m^3/(mol*s)'), n=1.86531, Ea=(171.386,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd',), comment="""Estimated from node HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd"""),
)

network(
    label = 'PDepNetwork #167',
    isomers = [
        'OC(F)OC(F)F(352)',
    ],
    reactants = [
        ('HF(38)', 'CHFO(46)', 'CHFO(46)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #167',
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

