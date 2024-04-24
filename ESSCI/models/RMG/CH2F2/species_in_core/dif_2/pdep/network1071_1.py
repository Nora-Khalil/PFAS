species(
    label = 'FC1(F)C=C1(2805)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
4 C u0 p0 c0 {3,S} {5,D} {6,S}
5 C u0 p0 c0 {3,S} {4,D} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-162.834,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,180,180,367.544,1003.92,1003.93,1003.97,1003.97,1003.98,1003.98,1003.99,1836.32],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (76.0447,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.94756,0.00289645,7.1868e-05,-1.26533e-07,6.69252e-11,-19583.6,7.86706], Tmin=(10,'K'), Tmax=(615.385,'K')), NASAPolynomial(coeffs=[2.03451,0.0280179,-2.02892e-05,6.80462e-09,-8.52836e-13,-19588.4,14.2158], Tmin=(615.385,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-162.834,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""FC1(F)CDC1""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'C2H2(70)',
    structure = adjacencyList("""1 C u0 p0 c0 {2,T} {3,S}
2 C u0 p0 c0 {1,T} {4,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (217.784,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,584.389,584.389,2772.01],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (26.0372,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.80868,0.0233616,-3.55172e-05,2.80153e-08,-8.50075e-12,26429,13.9397], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.65878,0.00488397,-1.60829e-06,2.46975e-10,-1.38606e-14,25759.4,-3.99838], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(217.784,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(87.302,'J/(mol*K)'), label="""C2H2""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'FC1=CC1F(2066)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {5,S}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4 C u0 p0 c0 {3,S} {5,D} {7,S}
5 C u0 p0 c0 {2,S} {3,S} {4,D}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-116.215,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,362.965,364.521,365.785,371.642,1108.64,1109.45,1111.66,1112.61,1113.37,1113.72,2091.2],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (76.0447,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.94882,0.00290669,7.19373e-05,-1.30708e-07,7.19211e-11,-13976,9.81734], Tmin=(10,'K'), Tmax=(588.625,'K')), NASAPolynomial(coeffs=[2.23287,0.026471,-1.8446e-05,6.01458e-09,-7.38973e-13,-13980.2,15.4347], Tmin=(588.625,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-116.215,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""FC1DCC1F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[CH]C=C(F)F(2812)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {4,D} {5,S} {6,S}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
5 C u2 p0 c0 {3,S} {7,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-20.8714,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,413.662,414.021,415.836],'cm^-1')),
        HinderedRotor(inertia=(0.418708,'amu*angstrom^2'), symmetry=1, barrier=(51.5185,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (76.0447,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.51816,0.0287758,-1.6516e-05,4.47079e-09,-4.78588e-13,-2453.53,15.7021], Tmin=(100,'K'), Tmax=(2137.89,'K')), NASAPolynomial(coeffs=[11.4201,0.0121203,-4.83013e-06,8.26775e-10,-5.24689e-14,-6259.86,-34.0088], Tmin=(2137.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-20.8714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(AllylJ2_triplet)"""),
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
    label = 'F[C]1C=C1(2884)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 C u1 p0 c0 {1,S} {3,S} {4,S}
3 C u0 p0 c0 {2,S} {4,D} {5,S}
4 C u0 p0 c0 {2,S} {3,D} {6,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (272.445,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,568.162,580.286,581.34,585.263,1088.39,1096.8,1098.06,2095.96],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.0463,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.09484,-0.00895101,0.000112671,-2.02633e-07,1.16837e-10,32769.4,8.63914], Tmin=(10,'K'), Tmax=(563.668,'K')), NASAPolynomial(coeffs=[2.58143,0.0201187,-1.34664e-05,4.24522e-09,-5.06175e-13,32648.8,12.4898], Tmin=(563.668,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(272.445,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""F[C]1CDC1""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'FC1(F)[C]=C1(2885)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
4 C u0 p0 c0 {3,S} {5,D} {6,S}
5 C u1 p0 c0 {3,S} {4,D}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (71.299,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,180,180,532.916,941.038,941.198,941.403,941.744,941.763,941.801,1783.81],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (75.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.94042,0.00332616,6.51123e-05,-1.21761e-07,6.69658e-11,8577.56,9.31645], Tmin=(10,'K'), Tmax=(613.804,'K')), NASAPolynomial(coeffs=[3.31003,0.0220561,-1.6392e-05,5.57215e-09,-7.03772e-13,8379.51,9.80628], Tmin=(613.804,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(71.299,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""FC1(F)[C]DC1""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC1(F)[C]C1(2886)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
4 C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
5 C u0 p1 c0 {3,S} {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (95.833,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,417.319,417.32,417.321,417.321,417.322,417.325,706.916,1339.22,1339.22,1339.23,1339.23],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (76.0447,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.90624,0.00588841,9.08554e-05,-2.04946e-07,1.35751e-10,11532.3,9.69899], Tmin=(10,'K'), Tmax=(517.5,'K')), NASAPolynomial(coeffs=[4.11155,0.0236024,-1.64341e-05,5.33953e-09,-6.53013e-13,11252.6,6.34668], Tmin=(517.5,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(95.833,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""FC1(F)[C]C1""", comment="""Thermo library: 2-BTP_G4"""),
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
    label = 'FC1=C=C1(2887)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 C u0 p0 c0 {3,S} {4,D} {5,S}
3 C u0 p0 c0 {1,S} {2,S} {4,D}
4 C u0 p0 c0 {2,D} {3,D}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (310.491,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,180,180,667.844,963.626,963.643,963.646,963.647],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.0383,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.04038,0.016802,-5.40047e-06,-8.25831e-09,5.21795e-12,37381.9,7.19607], Tmin=(100,'K'), Tmax=(953.46,'K')), NASAPolynomial(coeffs=[8.69708,0.00353864,-1.00245e-06,1.8125e-10,-1.41143e-14,35827.4,-22.3198], Tmin=(953.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(310.491,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(CdCddCF) + group(Cdd-CdsCds) + ring(Cyclopropadiene)"""),
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
    E0 = (-33.529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (121.713,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-63.2434,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (297.579,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (235.346,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (72.5886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (159.658,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['CF2(87)', 'C2H2(70)'],
    products = ['FC1(F)C=C1(2805)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.77e+15,'cm^3/(mol*s)'), n=-0.662, Ea=(0.157737,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(2000,'K'), comment="""Estimated using template [carbene;mb_tb_unsub] for rate rule [CF2;mb_tb_unsub]
Euclidian distance = 1.0
family: 1+2_Cycloaddition"""),
)

reaction(
    label = 'reaction2',
    reactants = ['FC1=CC1F(2066)'],
    products = ['FC1(F)C=C1(2805)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(285.686,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]C=C(F)F(2812)'],
    products = ['FC1(F)C=C1(2805)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.86291e+13,'s^-1'), n=0.0122167, Ea=(5.3862,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_single;Ypri_rad_out] + [R3;Y_rad_out;Ypri_rad_out] for rate rule [R3_SD;C_rad_out_noH;CdsinglepriH_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['F(37)', 'F[C]1C=C1(2884)'],
    products = ['FC1(F)C=C1(2805)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.71905e+10,'m^3/(mol*s)'), n=-0.966851, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.07464897564796032, var=0.299509236916578, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(5)', 'FC1(F)[C]=C1(2885)'],
    products = ['FC1(F)C=C1(2805)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.37186e+06,'m^3/(mol*s)'), n=0.39882, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_2BrCClFHNO-inRing_Ext-2BrCClFHNO-R_Ext-3R!H-R_Sp-4R!H-3R!H_Sp-3R!H-2BrCClFHNO_4R!H-inRing_Ext-3R!H-R',), comment="""Estimated from node Root_1R->H_N-2R->S_2BrCClFHNO-inRing_Ext-2BrCClFHNO-R_Ext-3R!H-R_Sp-4R!H-3R!H_Sp-3R!H-2BrCClFHNO_4R!H-inRing_Ext-3R!H-R"""),
)

reaction(
    label = 'reaction6',
    reactants = ['FC1(F)[C]C1(2886)'],
    products = ['FC1(F)C=C1(2805)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4.44322e+17,'s^-1'), n=-1.50998, Ea=(24.5138,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.5018375041037424, var=16.310967916345923, Tref=1000.0, N=4, data_mean=0.0, correlation='CCH_Ext-3C-R_N-4R!H->Br_N-4CClFINOPSSi->F_1C-inRing',), comment="""Estimated from node CCH_Ext-3C-R_N-4R!H->Br_N-4CClFINOPSSi->F_1C-inRing
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction7',
    reactants = ['HF(38)', 'FC1=C=C1(2887)'],
    products = ['FC1(F)C=C1(2805)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(178.039,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

network(
    label = 'PDepNetwork #1071',
    isomers = [
        'FC1(F)C=C1(2805)',
    ],
    reactants = [
        ('CF2(87)', 'C2H2(70)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1071',
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

