species(
    label = '[O]OOC([O])F(140)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {3,S} {6,S}
3 O u0 p2 c0 {2,S} {5,S}
4 O u1 p2 c0 {6,S}
5 O u1 p2 c0 {3,S}
6 C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
7 H u0 p0 c0 {6,S}
"""),
    E0 = (-131.798,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([451,553,637,1069,1180,1265,1301,3056,180,342.147,521.907,522.255,3457.88],'cm^-1')),
        HinderedRotor(inertia=(0.211554,'amu*angstrom^2'), symmetry=1, barrier=(41.0052,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0149407,'amu*angstrom^2'), symmetry=1, barrier=(2.89298,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.02817,0.0366543,-3.84973e-05,1.93692e-08,-3.69751e-12,-15774.9,22.6682], Tmin=(100,'K'), Tmax=(1392.19,'K')), NASAPolynomial(coeffs=[12.3025,0.00456309,-1.15059e-06,1.58656e-10,-9.57176e-15,-18386.5,-29.4042], Tmin=(1392.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-131.798,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsOs) + group(O2s-OsH) + group(CsFHOO) + radical(O2sj(Cs-F1sO2sH)) + radical(ROOJ)"""),
)

species(
    label = 'O2(2)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {2,S}
2 O u1 p2 c0 {1,S}
"""),
    E0 = (-8.62683,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1487.4],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.6,'angstroms^3'), rotrelaxcollnum=3.8, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53732,-0.00121571,5.31618e-06,-4.89443e-09,1.45845e-12,-1038.59,4.68368], Tmin=(100,'K'), Tmax=(1074.56,'K')), NASAPolynomial(coeffs=[3.15382,0.00167804,-7.69971e-07,1.51275e-10,-1.08782e-14,-1040.82,6.16754], Tmin=(1074.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.62683,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'FC1OO1(129)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {4,S}
3 O u0 p2 c0 {2,S} {4,S}
4 C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
5 H u0 p0 c0 {4,S}
"""),
    E0 = (-250.373,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,616.963,617.097,617.511,915.899,1417.41,1417.45,1417.51],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (64.0158,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2964.22,'J/mol'), sigma=(4.92574,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=463.00 K, Pc=56.28 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.07975,-0.00670713,7.61072e-05,-1.23288e-07,6.38153e-11,-30112,8.13301], Tmin=(10,'K'), Tmax=(624.171,'K')), NASAPolynomial(coeffs=[2.49608,0.0165225,-1.11531e-05,3.48846e-09,-4.10791e-13,-30169.1,12.9858], Tmin=(624.171,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-250.373,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""FC1OO1""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]C([O])F(131)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u1 p2 c0 {4,S}
3 O u1 p2 c0 {4,S}
4 C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
5 H u0 p0 c0 {4,S}
"""),
    E0 = (-145.018,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([451,553,637,1069,1180,1265,1301,3056,1064.74],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (64.0158,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3649.77,'J/mol'), sigma=(5.8153,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=570.09 K, Pc=42.11 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.47854,0.0027869,3.40158e-05,-4.87418e-08,1.9469e-11,-17414.9,12.3961], Tmin=(100,'K'), Tmax=(944.743,'K')), NASAPolynomial(coeffs=[8.92813,0.00190223,1.90731e-07,-1.26307e-11,-4.31643e-15,-19434.8,-18.8259], Tmin=(944.743,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-145.018,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-F1sO2sH)) + radical(O2sj(Cs-F1sO2sH))"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,1.49298e-14,-2.05327e-17,9.23927e-21,-1.28064e-24,29230.2,5.12616], Tmin=(100,'K'), Tmax=(3959.16,'K')), NASAPolynomial(coeffs=[2.5,1.19009e-10,-4.23987e-14,6.68965e-18,-3.94352e-22,29230.2,5.12617], Tmin=(3959.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.034,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[O]OC([O])F(211)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u1 p2 c0 {5,S}
4 O u1 p2 c0 {2,S}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
6 H u0 p0 c0 {5,S}
"""),
    E0 = (-167.775,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,451,553,637,1069,1180,1265,1301,3056],'cm^-1')),
        HinderedRotor(inertia=(0.491221,'amu*angstrom^2'), symmetry=1, barrier=(11.2941,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.73707,0.0289956,-3.52662e-05,2.3045e-08,-6.07486e-12,-20134.2,16.0103], Tmin=(100,'K'), Tmax=(921.435,'K')), NASAPolynomial(coeffs=[7.07196,0.0101777,-4.63276e-06,8.81524e-10,-6.15894e-14,-20933.1,-4.54834], Tmin=(921.435,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-167.775,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(CsFHOO) + radical(O2sj(Cs-F1sO2sH)) + radical(ROOJ)"""),
)

species(
    label = '[O]OO[CH]F(212)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {3,S} {5,S}
3 O u0 p2 c0 {2,S} {4,S}
4 O u1 p2 c0 {3,S}
5 C u1 p0 c0 {1,S} {2,S} {6,S}
6 H u0 p0 c0 {5,S}
"""),
    E0 = (25.1845,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([580,1155,1237,1373,3147,368.017,368.133,368.218,1887.73,4000],'cm^-1')),
        HinderedRotor(inertia=(0.165211,'amu*angstrom^2'), symmetry=1, barrier=(15.9052,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.478144,'amu*angstrom^2'), symmetry=1, barrier=(45.9789,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.63192,0.0302498,-3.83392e-05,2.46662e-08,-6.21158e-12,3078.15,18.3779], Tmin=(100,'K'), Tmax=(976.553,'K')), NASAPolynomial(coeffs=[8.24696,0.00725037,-3.01159e-06,5.48905e-10,-3.74905e-14,1981.47,-8.57817], Tmin=(976.553,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(25.1845,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsH) + group(CsFHHO) + radical(ROOJ) + radical(CsF1sHO2s)"""),
)

species(
    label = 'FC1OOOO1(213)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {4,S} {6,S}
3 O u0 p2 c0 {5,S} {6,S}
4 O u0 p2 c0 {2,S} {5,S}
5 O u0 p2 c0 {3,S} {4,S}
6 C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
7 H u0 p0 c0 {6,S}
"""),
    E0 = (-257.182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.14412,-0.00774365,0.000117949,-1.60457e-07,6.41312e-11,-30876.4,17.1342], Tmin=(100,'K'), Tmax=(930.029,'K')), NASAPolynomial(coeffs=[18.6843,-0.00632619,5.57741e-06,-9.95418e-10,5.43257e-14,-36718.9,-72.5807], Tmin=(930.029,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-257.182,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsOs) + group(CsFHOO) + ring(Cyclopentane)"""),
)

species(
    label = 'FC1OOO1(214)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u0 p2 c0 {4,S} {5,S}
4 O u0 p2 c0 {2,S} {3,S}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
6 H u0 p0 c0 {5,S}
"""),
    E0 = (-211.392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0152,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.17664,-0.000295691,7.50019e-05,-1.03847e-07,4.09613e-11,-25378.2,12.0172], Tmin=(100,'K'), Tmax=(950.174,'K')), NASAPolynomial(coeffs=[14.3675,-0.000485961,1.23112e-06,-1.16925e-10,-4.84091e-15,-29622.9,-52.5456], Tmin=(950.174,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-211.392,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsFHOO) + ring(Cyclobutane)"""),
)

species(
    label = 'O=C(F)OOO(215)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {3,S} {6,S}
3 O u0 p2 c0 {2,S} {4,S}
4 O u0 p2 c0 {3,S} {7,S}
5 O u0 p2 c0 {6,D}
6 C u0 p0 c0 {1,S} {2,S} {5,D}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-563.958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.14582,0.0236363,2.83453e-05,-6.862e-08,3.1757e-11,-67745.9,20.7459], Tmin=(100,'K'), Tmax=(932.857,'K')), NASAPolynomial(coeffs=[19.0565,-0.00499796,3.83551e-06,-6.83612e-10,3.8107e-14,-72810,-69.8957], Tmin=(932.857,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-563.958,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(145.503,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-OsOs) + group(O2s-OsH) + group(COFOO)"""),
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
    label = '[O]OOC=O(216)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,S} {5,S}
2 O u0 p2 c0 {1,S} {4,S}
3 O u0 p2 c0 {5,D}
4 O u1 p2 c0 {2,S}
5 C u0 p0 c0 {1,S} {3,D} {6,S}
6 H u0 p0 c0 {5,S}
"""),
    E0 = (-189.35,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,511.104,511.448,512.158,512.166],'cm^-1')),
        HinderedRotor(inertia=(0.000644707,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.246701,'amu*angstrom^2'), symmetry=1, barrier=(45.934,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (77.0162,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.00133,0.00571179,5.63158e-05,-8.74819e-08,3.66576e-11,-22722.6,17.7464], Tmin=(100,'K'), Tmax=(925.97,'K')), NASAPolynomial(coeffs=[14.8314,-0.00354144,3.51145e-06,-6.55266e-10,3.76131e-14,-26707.6,-48.1048], Tmin=(925.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-189.35,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-OsOs) + group(O2s-OsH) + group(Cds-OdOsH) + radical(ROOJ)"""),
)

species(
    label = '[O]O[O](152)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {3,S}
2 O u1 p2 c0 {1,S}
3 O u1 p2 c0 {1,S}
"""),
    E0 = (192.544,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([767.866,3898.78,3898.78],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (47.9982,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.47471,0.00468947,-6.33886e-06,3.99267e-09,-7.67832e-13,23182.7,11.3223], Tmin=(100,'K'), Tmax=(1873.31,'K')), NASAPolynomial(coeffs=[2.72432,0.000671973,1.37793e-06,-3.54951e-10,2.60883e-14,24449.9,18.0454], Tmin=(1873.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.544,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsOs) + group(O2s-OsH) + group(O2s-OsH) + radical(ROOJ) + radical(ROOJ)"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05425,-0.00375329,3.18277e-05,-4.07297e-08,1.68956e-11,-47423,6.56837], Tmin=(10,'K'), Tmax=(745.315,'K')), NASAPolynomial(coeffs=[2.27714,0.0104751,-6.24845e-06,1.77292e-09,-1.93455e-13,-47288.4,13.7455], Tmin=(745.315,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-394.294,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""ODCF""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]OOC(=O)F(146)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {3,S} {6,S}
3 O u0 p2 c0 {2,S} {5,S}
4 O u0 p2 c0 {6,D}
5 O u1 p2 c0 {3,S}
6 C u0 p0 c0 {1,S} {2,S} {4,D}
"""),
    E0 = (-411.953,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([482,664,788,1296,1923,589.236,589.237,589.241,589.241,589.242],'cm^-1')),
        HinderedRotor(inertia=(0.116089,'amu*angstrom^2'), symmetry=1, barrier=(28.6019,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116086,'amu*angstrom^2'), symmetry=1, barrier=(28.6018,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.0066,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3633.72,'J/mol'), sigma=(5.52188,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=567.58 K, Pc=48.97 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.42965,0.0198099,2.44594e-05,-6.03764e-08,2.84581e-11,-49476.3,20.4862], Tmin=(100,'K'), Tmax=(920.985,'K')), NASAPolynomial(coeffs=[17.2534,-0.00628772,4.61037e-06,-8.72837e-10,5.39354e-14,-53830.4,-58.6248], Tmin=(920.985,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-411.953,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-OsOs) + group(O2s-OsH) + group(COFOO) + radical(ROOJ)"""),
)

species(
    label = '[O][CH]F(127)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53592,0.00799823,-5.22203e-07,-4.56621e-09,2.08755e-12,-2348.33,9.31391], Tmin=(100,'K'), Tmax=(1079.18,'K')), NASAPolynomial(coeffs=[5.97186,0.00388221,-1.6298e-06,3.36439e-10,-2.54189e-14,-3160.17,-3.94918], Tmin=(1079.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-19.6796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsOJ) + radical(CsF1sHO2s)"""),
)

species(
    label = '[O]OO[C](O)F(217)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {4,S} {6,S}
3 O u0 p2 c0 {6,S} {7,S}
4 O u0 p2 c0 {2,S} {5,S}
5 O u1 p2 c0 {4,S}
6 C u1 p0 c0 {1,S} {2,S} {3,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (-181.142,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06048,0.0401423,-4.84397e-05,2.82256e-08,-6.31286e-12,-21714.3,23.196], Tmin=(100,'K'), Tmax=(1107.16,'K')), NASAPolynomial(coeffs=[11.5068,0.006014,-2.20152e-06,3.83369e-10,-2.5937e-14,-23806,-23.3384], Tmin=(1107.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-181.142,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(145.503,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsOs) + group(O2s-OsH) + group(CsFHOO) + radical(ROOJ) + radical(CsF1sO2sO2s)"""),
)

species(
    label = '[O][C](F)OOO(218)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {3,S} {6,S}
3 O u0 p2 c0 {2,S} {4,S}
4 O u0 p2 c0 {3,S} {7,S}
5 O u1 p2 c0 {6,S}
6 C u1 p0 c0 {1,S} {2,S} {5,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-77.6159,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([482,586,761,1411,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2569,0.0313724,-1.72658e-05,-5.58388e-09,5.63491e-12,-9266.05,23.2678], Tmin=(100,'K'), Tmax=(982.895,'K')), NASAPolynomial(coeffs=[12.4199,0.00573497,-2.13408e-06,4.26878e-10,-3.32573e-14,-12023.3,-29.4508], Tmin=(982.895,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-77.6159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(145.503,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsOs) + group(O2s-OsH) + group(CsFHOO) + radical(O2sj(Cs-F1sO2sH)) + radical(CsF1sO2sO2s)"""),
)

species(
    label = '[O]OO[CH]OF(219)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {4,S} {6,S}
3 O u0 p2 c0 {1,S} {6,S}
4 O u0 p2 c0 {2,S} {5,S}
5 O u1 p2 c0 {4,S}
6 C u1 p0 c0 {2,S} {3,S} {7,S}
7 H u0 p0 c0 {6,S}
"""),
    E0 = (171.975,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3025,407.5,1350,352.5,180,180,595.147,595.728,597.579,4000],'cm^-1')),
        HinderedRotor(inertia=(0.05851,'amu*angstrom^2'), symmetry=1, barrier=(14.7801,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184685,'amu*angstrom^2'), symmetry=1, barrier=(46.697,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106071,'amu*angstrom^2'), symmetry=1, barrier=(14.7859,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51139,0.054079,-8.01787e-05,5.51085e-08,-1.43249e-11,20774,22.4271], Tmin=(100,'K'), Tmax=(959.669,'K')), NASAPolynomial(coeffs=[13.6581,0.00345379,-1.05523e-06,1.46678e-10,-8.00866e-15,18442.4,-35.6748], Tmin=(959.669,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(171.975,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(145.503,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-OsOs) + group(O2s-OsH) + group(Cs-OsOsHH) + radical(ROOJ) + radical(OCJO)"""),
)

species(
    label = '[O][CH]OOOF(220)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {6,S}
3 O u0 p2 c0 {2,S} {4,S}
4 O u0 p2 c0 {1,S} {3,S}
5 O u1 p2 c0 {6,S}
6 C u1 p0 c0 {2,S} {5,S} {7,S}
7 H u0 p0 c0 {6,S}
"""),
    E0 = (208.161,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([277,555,632,3025,407.5,1350,352.5,506.162,506.166,506.169,506.171,506.176],'cm^-1')),
        HinderedRotor(inertia=(0.000657964,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000657974,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00065795,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9762,0.051469,-9.3811e-05,9.00724e-08,-3.32842e-11,25102.1,21.683], Tmin=(100,'K'), Tmax=(818.287,'K')), NASAPolynomial(coeffs=[5.53001,0.0202904,-1.13487e-05,2.27018e-09,-1.59595e-13,24982.7,8.07509], Tmin=(818.287,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(208.161,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(145.503,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsOs) + group(O2sFO) + group(Cs-OsOsHH) + radical(OCOJ) + radical(OCJO)"""),
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
    E0 = (2.6308,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (111.998,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (304.957,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-87.9467,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (70.0066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-70.0862,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-13.2659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-43.4071,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-72.063,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-79.5182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (209.603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (18.9125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (43.0118,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (288.167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (314.659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]OOC([O])F(140)'],
    products = ['O2(2)', 'FC1OO1(129)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(4.11041e+14,'s^-1'), n=-0.372113, Ea=(97.6903,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO_S;Y_rad_intra;OO_intra]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O(6)', '[O]OC([O])F(211)'],
    products = ['[O]OOC([O])F(140)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['O(6)', '[O]OO[CH]F(212)'],
    products = ['[O]OOC([O])F(140)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_sec_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]OOC([O])F(140)'],
    products = ['FC1OOOO1(213)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SSSS;Y_rad_out;Ypri_rad_out] for rate rule [R5_SSSS;O_rad;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]OOC([O])F(140)'],
    products = ['O(6)', 'FC1OOO1(214)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.67697e+10,'s^-1'), n=0.54, Ea=(165.066,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnOO;Y_rad_intra;OOJ] + [R3OO;Y_rad_intra;OO] for rate rule [R3OO;Y_rad_intra;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]OOC([O])F(140)'],
    products = ['O=C(F)OOO(215)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F(37)', '[O]OOC=O(216)'],
    products = ['[O]OOC([O])F(140)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.08261e+20,'m^3/(mol*s)'), n=-5.07836, Ea=(66.4544,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]O[O](152)', 'CHFO(46)'],
    products = ['[O]OOC([O])F(140)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.35412e-07,'m^3/(mol*s)'), n=3.42791, Ea=(121.605,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.9589805767432776, var=2.714361599831737, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_N-3R->C_2R!H->O_Ext-1R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_N-3R->C_2R!H->O_Ext-1R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(5)', '[O]OOC(=O)F(146)'],
    products = ['[O]OOC([O])F(140)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(15.9,'m^3/(mol*s)'), n=1.84, Ea=(91.3469,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_N-Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_N-Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction10',
    reactants = ['O2(2)', '[O]C([O])F(131)'],
    products = ['[O]OOC([O])F(140)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(9.96604e+08,'m^3/(mol*s)'), n=0.0472428, Ea=(37.3882,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]O[O](152)', '[O][CH]F(127)'],
    products = ['[O]OOC([O])F(140)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.01e+06,'m^3/(mol*s)'), n=-1.3397e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Sp-4R!H-2C_Ext-2C-R_N-5R!H->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Sp-4R!H-2C_Ext-2C-R_N-5R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]OOC([O])F(140)'],
    products = ['[O]OO[C](O)F(217)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.15e+14,'s^-1','+|-',2), n=-0.27, Ea=(113.972,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1800,'K'), comment="""Estimated using template [R2H_S;O_rad_out;Cs_H_out] for rate rule [R2H_S;O_rad_out;Cs_H_out_OOH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O][C](F)OOO(218)'],
    products = ['[O]OOC([O])F(140)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;Y_rad_out;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]OO[CH]OF(219)'],
    products = ['[O]OOC([O])F(140)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(79.4537,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O][CH]OOOF(220)'],
    products = ['[O]OOC([O])F(140)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(69.7596,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

network(
    label = 'PDepNetwork #10',
    isomers = [
        '[O]OOC([O])F(140)',
    ],
    reactants = [
        ('O2(2)', 'FC1OO1(129)'),
        ('O2(2)', '[O]C([O])F(131)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #10',
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

