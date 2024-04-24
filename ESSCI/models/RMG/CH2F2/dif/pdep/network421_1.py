species(
    label = 'CC(F)(F)F(1133)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-767.338,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,193,295,551,588,656,1146,1192,1350],'cm^-1')),
        HinderedRotor(inertia=(0.272102,'amu*angstrom^2'), symmetry=1, barrier=(6.25615,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0404,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.92157,0.00475571,9.14975e-05,-1.88599e-07,1.16734e-10,-92284.5,8.23158], Tmin=(10,'K'), Tmax=(539.223,'K')), NASAPolynomial(coeffs=[3.23895,0.026964,-1.79736e-05,5.70921e-09,-6.91212e-13,-92460.1,8.79206], Tmin=(539.223,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-767.338,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""CC(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[CH2]C(F)(F)F(1125)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 F u0 p3 c0 {4,S}
4 C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
5 C u1 p0 c0 {4,S} {6,S} {7,S}
6 H u0 p0 c0 {5,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-547.919,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,522,585,610,849,1160,1215,1402,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.0709201,'amu*angstrom^2'), symmetry=1, barrier=(7.86184,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0324,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2403.71,'J/mol'), sigma=(4.911,'angstroms'), dipoleMoment=(2.3,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.90955,0.00544702,8.47571e-05,-1.82077e-07,1.1401e-10,-65893.7,10.0072], Tmin=(10,'K'), Tmax=(552.079,'K')), NASAPolynomial(coeffs=[4.35096,0.0223417,-1.57381e-05,5.20006e-09,-6.46882e-13,-66248.6,5.36671], Tmin=(552.079,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-547.919,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""[CH2]C(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CHF3(41)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 F u0 p3 c0 {4,S}
4 C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
5 H u0 p0 c0 {4,S}
"""),
    E0 = (-707.987,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(70.003,'amu')),
        NonlinearRotor(inertia=([48.8559,48.8561,89.1073],'amu*angstrom^2'), symmetry=3),
        HarmonicOscillator(frequencies=([506.193,506.206,703.466,1151.07,1179.84,1179.87,1414.02,1414.02,3086.64],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (70.0138,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2178.39,'J/mol'), sigma=(4.123,'angstroms'), dipoleMoment=(1.8,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.06655,-0.00569699,7.11567e-05,-1.16193e-07,6.06504e-11,-85150.3,7.48265], Tmin=(10,'K'), Tmax=(615.111,'K')), NASAPolynomial(coeffs=[2.45498,0.0163704,-1.09138e-05,3.38161e-09,-3.95687e-13,-85171.3,12.6925], Tmin=(615.111,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-707.987,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""FC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CH2(S)(24)',
    structure = adjacencyList("""1 C u0 p1 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (419.091,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1369.93,2896.01,2896.03],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.10264,-0.00144069,5.4507e-06,-3.58003e-09,7.56196e-13,50400.6,-0.411767], Tmin=(100,'K'), Tmax=(1442.35,'K')), NASAPolynomial(coeffs=[2.62647,0.00394764,-1.49925e-06,2.54541e-10,-1.62957e-14,50691.8,6.78384], Tmin=(1442.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(419.091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: primaryThermoLibrary"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(897.963,'J/mol'), sigma=(3.977,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28591,0.0107608,-1.05382e-05,4.89881e-09,-8.86384e-13,-24340.7,13.1348], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[5.33121,0.00197748,-9.60248e-07,2.10704e-10,-1.5954e-14,-25190.9,-2.56367], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-203.712,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CF2""", comment="""Thermo library: halogens"""),
)

species(
    label = 'CH3F(39)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (-246.136,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(34.0219,'amu')),
        NonlinearRotor(inertia=([3.21737,19.5168,19.5168],'amu*angstrom^2'), symmetry=3),
        HarmonicOscillator(frequencies=([1097.73,1206.94,1206.96,1493.65,1493.66,1509.77,3011.42,3086.11,3086.13],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (34.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05191,-0.00340134,1.7341e-05,2.49456e-09,-1.5179e-11,-29605.1,3.91759], Tmin=(10,'K'), Tmax=(534.708,'K')), NASAPolynomial(coeffs=[0.487074,0.0143029,-7.17972e-06,1.7168e-09,-1.57925e-13,-29095.7,20.0824], Tmin=(534.708,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-246.136,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CF""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'CH3-CF2(64)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (-315.516,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,190,488,555,1236,1407],'cm^-1')),
        HinderedRotor(inertia=(0.203516,'amu*angstrom^2'), symmetry=1, barrier=(4.67924,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (65.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.97026,0.00248557,9.36596e-05,-2.62618e-07,2.48555e-10,-37947.6,9.23011], Tmin=(10,'K'), Tmax=(267.682,'K')), NASAPolynomial(coeffs=[2.69162,0.021592,-1.34044e-05,4.02185e-09,-4.66099e-13,-37879.1,13.7136], Tmin=(267.682,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-315.516,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""C[C](F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CF3(44)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 F u0 p3 c0 {4,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
"""),
    E0 = (-483.19,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(68.9952,'amu')),
        NonlinearRotor(inertia=([46.5949,46.5949,90.0934],'amu*angstrom^2'), symmetry=3),
        HarmonicOscillator(frequencies=([504.804,504.824,700.362,1092.66,1279.94,1279.95],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05587,-0.0054938,7.30062e-05,-1.34908e-07,7.857e-11,-58113,8.13986], Tmin=(10,'K'), Tmax=(570.231,'K')), NASAPolynomial(coeffs=[3.53924,0.0118582,-8.75006e-06,2.89362e-09,-3.54041e-13,-58277.3,8.38507], Tmin=(570.231,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-483.19,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""F[C](F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CH3(18)',
    structure = adjacencyList("""multiplicity 2
1 C u1 p0 c0 {2,S} {3,S} {4,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
"""),
    E0 = (136.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([604.263,1333.71,1492.19,2836.77,2836.77,3806.92],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (15.0346,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.65718,0.0021266,5.45839e-06,-6.6181e-09,2.46571e-12,16422.7,1.67354], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.97812,0.00579785,-1.97558e-06,3.07298e-10,-1.79174e-14,16509.5,4.72248], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(136.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH3""", comment="""Thermo library: FFCM1(-)"""),
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
    E0 = (-250.465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-59.4636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-25.191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-13.192,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-96.5667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-75.8879,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['CC(F)(F)F(1133)'],
    products = ['HF(38)', 'CH2CF2(56)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1e+14,'s^-1'), n=0, Ea=(68700,'cal/mol'), T0=(1,'K'), comment="""Matched reaction 2 C2H3F3 <=> CH2CF2 + HF in XY_Addition_MultipleBond/training
This reaction matched rate rule [HF_N-3COCdCddCtO2d->Ct_3CdO2d->Cd_Ext-4COCdCddCtO2d-R_5R!H->F_Ext-4COCdCddCtO2d-R]
family: XY_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CHF3(41)', 'CH2(S)(24)'],
    products = ['CC(F)(F)F(1133)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.67081e+07,'m^3/(mol*s)'), n=-0.25099, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05783420702249951, var=0.11014630093021238, Tref=1000.0, N=12, data_mean=0.0, correlation='CH_Ext-4CbCdCsCt-R',), comment="""Estimated from node CH_Ext-4CbCdCsCt-R"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CF2(42)', 'CH3F(39)'],
    products = ['CC(F)(F)F(1133)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(0.6158,'cm^3/(mol*s)'), n=3.54561, Ea=(195.225,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), Tmax=(2500,'K'), comment="""Fitted to 50 data points; dA = *|/ 1.20788, dn = +|- 0.0248129, dEa = +|- 0.135031 kJ/molMatched reaction 50 CF2_r1 + CH3F <=> C2H3F3-2 in 1,2_Insertion_carbene/training
This reaction matched rate rule [CY_N-2Br1sCl1sF1sHI1s->H_N-5Br1sCl1sF1sI1s->Cl1s_5Br1sF1s->F1s]
family: 1,2_Insertion_carbene"""),
)

reaction(
    label = 'reaction4',
    reactants = ['F(37)', 'CH3-CF2(64)'],
    products = ['CC(F)(F)F(1133)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CF3(44)', 'CH3(18)'],
    products = ['CC(F)(F)F(1133)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.78e+33,'cm^3/(mol*s)'), n=-6.64, Ea=(5020,'cal/mol'), T0=(1,'K'), comment="""Matched reaction 192 CH3 + CF3 <=> C2H3F3-2 in R_Recombination/training
This reaction matched rate rule [Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_2CF->C_Ext-1C-R_Ext-1C-R]
family: R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(5)', '[CH2]C(F)(F)F(1125)'],
    products = ['CC(F)(F)F(1133)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7.27e+42,'cm^3/(mol*s)'), n=-9.86, Ea=(7360,'cal/mol'), T0=(1,'K'), comment="""Matched reaction 185 C2H2F3-3 + H <=> C2H3F3-2 in R_Recombination/training
This reaction matched rate rule [Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_N-4R!H->C_Sp-3C-2C_Ext-3C-R_Ext-3C-R]
family: R_Recombination"""),
)

network(
    label = 'PDepNetwork #421',
    isomers = [
        'CC(F)(F)F(1133)',
    ],
    reactants = [
        ('HF(38)', 'CH2CF2(56)'),
        ('H(5)', '[CH2]C(F)(F)F(1125)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #421',
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

