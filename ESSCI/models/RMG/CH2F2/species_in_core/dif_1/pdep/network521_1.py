species(
    label = 'OC(F)(F)C(F)(F)C(F)(F)F(1363)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {11,S}
8  O u0 p2 c0 {10,S} {12,S}
9  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
10 C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
11 C u0 p0 c0 {5,S} {6,S} {7,S} {9,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-1771.58,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,222,329,445,522,589,1214,1475,223,363,546,575,694,1179,1410,193,295,551,588,656,1146,1192,1350,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.124592,'amu*angstrom^2'), symmetry=1, barrier=(2.86462,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123854,'amu*angstrom^2'), symmetry=1, barrier=(2.84765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.718217,'amu*angstrom^2'), symmetry=1, barrier=(16.5132,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (186.028,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.14326,0.075296,-8.20782e-05,3.27491e-08,1.11663e-12,-213065,15.125], Tmin=(10,'K'), Tmax=(700.031,'K')), NASAPolynomial(coeffs=[13.3418,0.0322632,-2.25284e-05,7.14035e-09,-8.45426e-13,-214866,-33.1072], Tmin=(700.031,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1771.58,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), label="""OC(F)(F)C(F)(F)C(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'CF2O(48)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
"""),
    E0 = (-650.702,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(65.9917,'amu')),
        NonlinearRotor(inertia=([42.7382,43.0674,85.8056],'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([584.127,619.74,793.429,989.231,1281.66,1989.03],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (66.0068,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2914.22,'J/mol'), sigma=(4.906,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17614,0.021087,-2.39216e-05,1.39949e-08,-3.37605e-12,-77951.3,19.1134], Tmin=(298,'K'), Tmax=(1150,'K')), NASAPolynomial(coeffs=[6.79131,0.00338305,-1.43012e-06,2.71226e-10,-1.90524e-14,-79454,-9.48312], Tmin=(1150,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-650.702,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(83.1447,'J/mol/K'), label="""CF2O""", comment="""Thermo library: Fluorine"""),
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
    E0 = (-675.664,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.59735,0.0339059,-4.12239e-05,2.60507e-08,-6.78805e-12,-81178.9,12.7972], Tmin=(298,'K'), Tmax=(1100,'K')), NASAPolynomial(coeffs=[11.3243,0.00499631,-2.11782e-06,3.98619e-10,-2.75039e-14,-83426.6,-31.2702], Tmin=(1100,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-675.664,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(133.032,'J/mol/K'), label="""CF2CF2""", comment="""Thermo library: Fluorine"""),
)

species(
    label = 'CF2(89)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p1 c0 {1,S} {2,S}
"""),
    E0 = (-196.899,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([192,594,627],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (50.0074,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(897.963,'J/mol'), sigma=(3.977,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28591,0.0107608,-1.05382e-05,4.89881e-09,-8.86384e-13,-23521.2,13.1348], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[5.33121,0.00197748,-9.60248e-07,2.10704e-10,-1.5954e-14,-24371.4,-2.56367], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-196.899,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CF2""", comment="""Thermo library: Fluorine"""),
)

species(
    label = 'OC(F)(F)C(F)(F)F(1019)',
    structure = adjacencyList("""1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {8,S}
4 F u0 p3 c0 {8,S}
5 F u0 p3 c0 {8,S}
6 O u0 p2 c0 {7,S} {9,S}
7 C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
8 C u0 p0 c0 {3,S} {4,S} {5,S} {7,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-1350.04,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,223,363,546,575,694,1179,1410,193,295,551,588,656,1146,1192,1350,180],'cm^-1')),
        HinderedRotor(inertia=(0.138227,'amu*angstrom^2'), symmetry=1, barrier=(3.17811,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.407754,'amu*angstrom^2'), symmetry=1, barrier=(9.37507,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.021,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.57082,0.0364322,3.32484e-05,-1.72022e-07,1.4359e-10,-162369,12.2295], Tmin=(10,'K'), Tmax=(492.269,'K')), NASAPolynomial(coeffs=[7.65639,0.0287673,-2.11978e-05,7.0786e-09,-8.77171e-13,-163080,-7.72772], Tmin=(492.269,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1350.04,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""OC(F)(F)C(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CF4(90)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 F u0 p3 c0 {5,S}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {4,S}
"""),
    E0 = (-946.658,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(87.9936,'amu')),
        NonlinearRotor(inertia=([88.4665,88.4666,88.467],'amu*angstrom^2'), symmetry=12),
        HarmonicOscillator(frequencies=([429.709,429.73,625.518,625.532,625.546,914.019,1289.73,1289.74,1289.75],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (88.0043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.329574,0.0332696,-3.88267e-05,2.20593e-08,-4.94837e-12,-113504,21.2017], Tmin=(298,'K'), Tmax=(1250,'K')), NASAPolynomial(coeffs=[9.80832,0.00328269,-1.34045e-06,2.43738e-10,-1.62862e-14,-116025,-27.0002], Tmin=(1250,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-946.655,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(108.088,'J/mol/K'), label="""CF4""", comment="""Thermo library: Fluorine"""),
)

species(
    label = 'OC(F)(F)[C]F(576)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {6,S}
4 O u0 p2 c0 {5,S} {7,S}
5 C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6 C u0 p1 c0 {3,S} {5,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-539.334,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,617,898,1187],'cm^-1')),
        HinderedRotor(inertia=(1.38402,'amu*angstrom^2'), symmetry=1, barrier=(31.8212,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.39593,'amu*angstrom^2'), symmetry=1, barrier=(32.0951,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0239,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03035,0.046278,-6.6989e-05,4.89416e-08,-1.41201e-11,-64798.6,17.5213], Tmin=(100,'K'), Tmax=(850.047,'K')), NASAPolynomial(coeffs=[9.35662,0.011803,-6.15344e-06,1.22947e-09,-8.77629e-14,-66044.1,-16.6334], Tmin=(850.047,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-539.334,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCFFO) + group(CJ2_singlet-FCs)"""),
)

species(
    label = 'OC(F)(F)F(254)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 O u0 p2 c0 {5,S} {6,S}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {4,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-924.469,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,431,527,613,668,835,1199,1245,1322],'cm^-1')),
        HinderedRotor(inertia=(0.246742,'amu*angstrom^2'), symmetry=1, barrier=(8.25205,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0132,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3135.61,'J/mol'), sigma=(5.37414,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=489.78 K, Pc=45.84 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93145,0.00398703,6.6773e-05,-1.35516e-07,8.01823e-11,-111184,9.48547], Tmin=(10,'K'), Tmax=(580.755,'K')), NASAPolynomial(coeffs=[4.01947,0.0193555,-1.41814e-05,4.77794e-09,-5.99414e-13,-111464,6.78935], Tmin=(580.755,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-924.469,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""OC(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[C]C(F)(F)F(307)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 F u0 p3 c0 {6,S}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
6 C u0 p1 c0 {4,S} {5,S}
"""),
    E0 = (-531.432,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([261,346,529,1231,1306,617,898,1187,310.428,1028.96,1029.05],'cm^-1')),
        HinderedRotor(inertia=(0.085726,'amu*angstrom^2'), symmetry=1, barrier=(5.82341,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.015,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88935,0.00712406,8.594e-05,-2.10572e-07,1.46542e-10,-63912.6,11.0368], Tmin=(10,'K'), Tmax=(506.569,'K')), NASAPolynomial(coeffs=[4.91752,0.0201688,-1.53534e-05,5.20649e-09,-6.50075e-13,-64288.3,4.0956], Tmin=(506.569,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-531.432,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""F[C]C(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OF(155)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (-95.2653,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(36.0011,'amu')),
        NonlinearRotor(inertia=([0.860315,18.4105,19.2708],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([1005.07,1417.2,3730.42],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (36.0058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02848,-0.00192233,1.33972e-05,-1.61993e-08,6.45714e-12,-11458,4.36077], Tmin=(10,'K'), Tmax=(768.674,'K')), NASAPolynomial(coeffs=[3.2293,0.00415811,-2.21832e-06,5.96331e-10,-6.32051e-14,-11391.9,7.63678], Tmin=(768.674,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-95.2653,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""OF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[C]C(F)(F)C(F)(F)F(1640)',
    structure = adjacencyList("""1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {8,S}
4 F u0 p3 c0 {8,S}
5 F u0 p3 c0 {8,S}
6 F u0 p3 c0 {9,S}
7 C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8 C u0 p0 c0 {3,S} {4,S} {5,S} {7,S}
9 C u0 p1 c0 {6,S} {7,S}
"""),
    E0 = (-983.983,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([347,453,1141,1468,193,295,551,588,656,1146,1192,1350,617,898,1187,201.873,201.895,201.916,2724.26],'cm^-1')),
        HinderedRotor(inertia=(0.396095,'amu*angstrom^2'), symmetry=1, barrier=(11.4611,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23123,'amu*angstrom^2'), symmetry=1, barrier=(35.616,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (150.022,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.877471,0.0745126,-0.000116664,9.28379e-08,-2.92724e-11,-118239,24.0019], Tmin=(100,'K'), Tmax=(778.336,'K')), NASAPolynomial(coeffs=[11.5611,0.0196054,-1.0843e-05,2.19546e-09,-1.57019e-13,-119902,-24.8627], Tmin=(778.336,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-983.983,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFF) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)2) + group(CJ2_singlet-FCs)"""),
)

species(
    label = 'OC(F)(C(F)(F)F)C(F)(F)F(1641)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {11,S}
8  O u0 p2 c0 {9,S} {12,S}
9  C u0 p0 c0 {1,S} {8,S} {10,S} {11,S}
10 C u0 p0 c0 {2,S} {3,S} {4,S} {9,S}
11 C u0 p0 c0 {5,S} {6,S} {7,S} {9,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-1777.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (186.028,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39813,0.0481848,0.000138213,-5.18744e-07,4.43172e-10,-213786,14.3013], Tmin=(10,'K'), Tmax=(455.196,'K')), NASAPolynomial(coeffs=[10.6569,0.0422826,-3.30801e-05,1.14821e-08,-1.46102e-12,-215046,-21.5923], Tmin=(455.196,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1777.54,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), label="""OC(F)(C(F)(F)F)C(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC(F)=C(F)C(F)(F)F(879)',
    structure = adjacencyList("""1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {7,S}
4 F u0 p3 c0 {8,S}
5 F u0 p3 c0 {9,S}
6 F u0 p3 c0 {9,S}
7 C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
8 C u0 p0 c0 {4,S} {7,S} {9,D}
9 C u0 p0 c0 {5,S} {6,S} {8,D}
"""),
    E0 = (-1171.65,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([219,296,586,564,718,793,1177,1228,323,467,575,827,1418,182,240,577,636,1210,1413,2018.52],'cm^-1')),
        HinderedRotor(inertia=(0.370476,'amu*angstrom^2'), symmetry=1, barrier=(8.51796,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (150.022,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.48756,0.0445189,-1.69994e-05,-4.62355e-08,4.1721e-11,-140912,13.4084], Tmin=(10,'K'), Tmax=(600.481,'K')), NASAPolynomial(coeffs=[8.64988,0.0279311,-2.00276e-05,6.49183e-09,-7.83507e-13,-141852,-11.5361], Tmin=(600.481,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1171.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""FC(F)DC(F)C(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'OC(F)(F)[C](F)C(F)(F)F(1279)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {10,S}
7  O u0 p2 c0 {8,S} {11,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {5,S} {10,S}
10 C u1 p0 c0 {6,S} {8,S} {9,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-1367.81,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,253,525,597,667,842,1178,1324,253,522,585,610,849,1160,1215,1402,212,367,445,1450,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.209374,'amu*angstrom^2'), symmetry=1, barrier=(4.81393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.209341,'amu*angstrom^2'), symmetry=1, barrier=(4.81317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.209356,'amu*angstrom^2'), symmetry=1, barrier=(4.81351,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (167.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.41756,0.0717297,0.00020306,-2.74937e-06,6.58302e-09,-164509,15.6032], Tmin=(10,'K'), Tmax=(175.741,'K')), NASAPolynomial(coeffs=[7.27623,0.0415986,-3.22097e-05,1.11931e-08,-1.43327e-12,-164734,1.16173], Tmin=(175.741,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1367.81,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), label="""OC(F)(F)[C](F)C(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O[C](F)C(F)(F)C(F)(F)F(1642)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {10,S}
7  O u0 p2 c0 {10,S} {11,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
10 C u1 p0 c0 {6,S} {7,S} {8,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-1350.94,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,215,315,519,588,595,1205,1248,193,295,551,588,656,1146,1192,1350,395,473,707,1436,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.44904,'amu*angstrom^2'), symmetry=1, barrier=(10.3243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.449708,'amu*angstrom^2'), symmetry=1, barrier=(10.3397,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.448414,'amu*angstrom^2'), symmetry=1, barrier=(10.3099,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (167.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31865,0.0573393,0.000102545,-5.75278e-07,6.01017e-10,-162484,14.4733], Tmin=(10,'K'), Tmax=(394.689,'K')), NASAPolynomial(coeffs=[11.2112,0.0343767,-2.69063e-05,9.43945e-09,-1.21594e-12,-163551,-21.8929], Tmin=(394.689,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1350.94,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), label="""O[C](F)C(F)(F)C(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OC(F)(F)C(F)(F)[C](F)F(1024)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  O u0 p2 c0 {9,S} {11,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
10 C u1 p0 c0 {5,S} {6,S} {8,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-1339.81,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,215,315,519,588,595,1205,1248,223,363,546,575,694,1179,1410,190,488,555,1236,1407,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.304132,'amu*angstrom^2'), symmetry=1, barrier=(6.99259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.303192,'amu*angstrom^2'), symmetry=1, barrier=(6.97097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.302103,'amu*angstrom^2'), symmetry=1, barrier=(6.94593,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (167.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.08592,0.0848458,-0.000149394,1.39989e-07,-5.27882e-11,-161139,15.0603], Tmin=(10,'K'), Tmax=(639.355,'K')), NASAPolynomial(coeffs=[11.7585,0.0305876,-2.20983e-05,7.25571e-09,-8.86893e-13,-162248,-22.9006], Tmin=(639.355,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1339.81,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), label="""OC(F)(F)C(F)(F)[C](F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'F[C](F)C(F)(F)C(F)(F)F(1643)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {10,S}
7  F u0 p3 c0 {10,S}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
10 C u1 p0 c0 {6,S} {7,S} {8,S}
"""),
    E0 = (-1354.15,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,193,295,551,588,656,1146,1192,1350,190,488,555,1236,1407,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.250645,'amu*angstrom^2'), symmetry=1, barrier=(5.76281,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.25059,'amu*angstrom^2'), symmetry=1, barrier=(5.76156,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (169.021,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.10717,0.0723747,-0.00010609,7.89753e-08,-2.34473e-11,-162859,14.936], Tmin=(10,'K'), Tmax=(806.883,'K')), NASAPolynomial(coeffs=[12.7827,0.0244098,-1.69226e-05,5.30323e-09,-6.21197e-13,-164421,-29.6667], Tmin=(806.883,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1354.15,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""F[C](F)C(F)(F)C(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]C(F)(F)C(F)(F)C(F)(F)F(1644)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {11,S}
6  F u0 p3 c0 {11,S}
7  F u0 p3 c0 {11,S}
8  O u1 p2 c0 {10,S}
9  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
10 C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
11 C u0 p0 c0 {5,S} {6,S} {7,S} {9,S}
"""),
    E0 = (-1509.76,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([222,329,445,522,589,1214,1475,351,323,533,609,664,892,1120,1201,193,295,551,588,656,1146,1192,1350,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.31308,'amu*angstrom^2'), symmetry=1, barrier=(7.19831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.809423,'amu*angstrom^2'), symmetry=1, barrier=(18.6102,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (185.02,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.32864,0.0564558,5.21026e-05,-2.89456e-07,2.4747e-10,-181577,15.4201], Tmin=(10,'K'), Tmax=(490.723,'K')), NASAPolynomial(coeffs=[11.1876,0.0378079,-2.97084e-05,1.02701e-08,-1.29932e-12,-182895,-22.4714], Tmin=(490.723,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1509.76,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), label="""[O]C(F)(F)C(F)(F)C(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O[C](F)F(220)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {4,S} {5,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (-470.293,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,493,600,700,1144,1293],'cm^-1')),
        HinderedRotor(inertia=(0.307072,'amu*angstrom^2'), symmetry=1, barrier=(7.06019,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.0148,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3063.56,'J/mol'), sigma=(4.99758,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=478.52 K, Pc=55.69 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.94331,0.003562,5.16413e-05,-1.18327e-07,7.89485e-11,-56559.5,10.5484], Tmin=(10,'K'), Tmax=(519.044,'K')), NASAPolynomial(coeffs=[4.2538,0.012858,-9.00307e-06,2.95251e-09,-3.63989e-13,-56749.2,7.73732], Tmin=(519.044,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-470.293,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""O[C](F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[C](F)C(F)(F)F(424)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {6,S}
4 F u0 p3 c0 {7,S}
5 F u0 p3 c0 {7,S}
6 C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
7 C u1 p0 c0 {4,S} {5,S} {6,S}
"""),
    E0 = (-909.392,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,522,585,610,849,1160,1215,1402,190,488,555,1236,1407,508.62],'cm^-1')),
        HinderedRotor(inertia=(0.492968,'amu*angstrom^2'), symmetry=1, barrier=(11.3343,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (119.013,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2688.9,'J/mol'), sigma=(4.9,'angstroms'), dipoleMoment=(1.5,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57479,0.0461494,-5.69635e-05,3.48343e-08,-8.52354e-12,-109277,20.4809], Tmin=(298,'K'), Tmax=(1150,'K')), NASAPolynomial(coeffs=[13.9595,0.00491266,-2.10538e-06,3.99134e-10,-2.7649e-14,-112424,-41.9106], Tmin=(1150,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-909.395,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(153.818,'J/mol/K'), label="""CF3-CF2""", comment="""Thermo library: Fluorine"""),
)

species(
    label = 'CF3(44)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 F u0 p3 c0 {4,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
"""),
    E0 = (-479.453,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(68.9952,'amu')),
        NonlinearRotor(inertia=([46.5949,46.5949,90.0934],'amu*angstrom^2'), symmetry=3),
        HarmonicOscillator(frequencies=([504.804,504.824,700.362,1092.66,1279.94,1279.95],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0058,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1006.05,'J/mol'), sigma=(4.32,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35729,0.0215758,-2.3902e-05,1.29611e-08,-2.80173e-12,-57391.3,18.6474], Tmin=(298,'K'), Tmax=(1250,'K')), NASAPolynomial(coeffs=[7.73167,0.00227115,-8.90496e-07,1.52938e-10,-9.50838e-15,-59145.7,-14.0202], Tmin=(1250,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-479.453,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(83.1447,'J/mol/K'), label="""CF3""", comment="""Thermo library: Fluorine"""),
)

species(
    label = 'OC(F)(F)[C](F)F(586)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {7,S}
4 F u0 p3 c0 {7,S}
5 O u0 p2 c0 {6,S} {8,S}
6 C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
7 C u1 p0 c0 {3,S} {4,S} {6,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-907.917,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,253,525,597,667,842,1178,1324,190,488,555,1236,1407,180],'cm^-1')),
        HinderedRotor(inertia=(0.391062,'amu*angstrom^2'), symmetry=1, barrier=(8.99127,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.390952,'amu*angstrom^2'), symmetry=1, barrier=(8.98875,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (117.022,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.56174,0.0389036,1.94218e-05,-2.22218e-07,2.51043e-10,-109199,12.6539], Tmin=(10,'K'), Tmax=(390.957,'K')), NASAPolynomial(coeffs=[7.58136,0.0231046,-1.7134e-05,5.81753e-09,-7.33941e-13,-109707,-5.43885], Tmin=(390.957,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-907.917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""OC(F)(F)[C](F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C(F)C(F)(F)C(F)(F)F(1645)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {10,S}
7  O u0 p2 c0 {10,D}
8  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
10 C u0 p0 c0 {6,S} {7,D} {8,S}
"""),
    E0 = (-1466.96,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,348,427,405,1245,1236,1280,193,295,551,588,656,1146,1192,1350,486,617,768,1157,1926,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.174012,'amu*angstrom^2'), symmetry=1, barrier=(4.00087,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.383664,'amu*angstrom^2'), symmetry=1, barrier=(8.8212,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (166.022,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36994,0.0682513,-0.000108793,9.95316e-08,-3.84993e-11,-176437,13.9257], Tmin=(10,'K'), Tmax=(605.776,'K')), NASAPolynomial(coeffs=[8.40417,0.0350098,-2.64819e-05,8.94648e-09,-1.11532e-12,-177047,-7.83814], Tmin=(605.776,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1466.96,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""ODC(F)C(F)(F)C(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F2(108)',
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
    label = 'OC(F)=C(F)C(F)(F)F(1281)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {9,S}
6  O u0 p2 c0 {9,S} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
8  C u0 p0 c0 {4,S} {7,S} {9,D}
9  C u0 p0 c0 {5,S} {6,S} {8,D}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-1175.06,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,219,296,586,564,718,793,1177,1228,323,467,575,827,1418,326,540,652,719,1357,180],'cm^-1')),
        HinderedRotor(inertia=(0.371381,'amu*angstrom^2'), symmetry=1, barrier=(8.53878,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.370822,'amu*angstrom^2'), symmetry=1, barrier=(8.52593,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (148.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.62611,0.0272005,0.000179856,-5.83927e-07,4.96534e-10,-141325,13.4436], Tmin=(10,'K'), Tmax=(438.608,'K')), NASAPolynomial(coeffs=[8.85209,0.0334496,-2.5878e-05,9.00461e-09,-1.15195e-12,-142302,-13.373], Tmin=(438.608,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1175.06,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""OC(F)DC(F)C(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OC(F)(F)C(F)=C(F)F(382)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  O u0 p2 c0 {7,S} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
8  C u0 p0 c0 {3,S} {7,S} {9,D}
9  C u0 p0 c0 {4,S} {5,S} {8,D}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-1154.51,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,275,321,533,585,746,850,1103,323,467,575,827,1418,182,240,577,636,1210,1413,180],'cm^-1')),
        HinderedRotor(inertia=(0.102433,'amu*angstrom^2'), symmetry=1, barrier=(2.35514,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0987329,'amu*angstrom^2'), symmetry=1, barrier=(2.27006,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (148.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.56187,0.035635,0.000138624,-5.6967e-07,5.63035e-10,-138858,13.4425], Tmin=(10,'K'), Tmax=(387.981,'K')), NASAPolynomial(coeffs=[8.14435,0.0331468,-2.4791e-05,8.45085e-09,-1.06863e-12,-139551,-8.66771], Tmin=(387.981,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1154.51,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""OC(F)(F)C(F)DC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    E0 = (-759.301,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-468.123,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-541.763,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-522.415,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-332.398,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-737.061,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-299.601,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-535.333,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-559.192,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-542.532,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-531.401,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-590.235,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-562.436,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-644.167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-651.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-855.419,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-448.351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-427.8,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['OC(F)(F)C(F)(F)C(F)(F)F(1363)'],
    products = ['HF(38)', 'CF2O(48)', 'CF2CF2(60)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5.95626e+10,'s^-1'), n=0.978818, Ea=(276.759,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R!H->C_N-5Br1sCl1sF1sH->H_5Br1sCl1sF1s->F1s_Ext-2C-R_N-7R!H->O_Sp-7BrCClFINPSSi-2C_Ext-2C-R',), comment="""Estimated from node Root_1R!H->C_N-5Br1sCl1sF1sH->H_5Br1sCl1sF1s->F1s_Ext-2C-R_N-7R!H->O_Sp-7BrCClFINPSSi-2C_Ext-2C-R
Multiplied by reaction path degeneracy 3.0"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CF2(89)', 'OC(F)(F)C(F)(F)F(1019)'],
    products = ['OC(F)(F)C(F)(F)C(F)(F)F(1363)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.22791e-07,'m^3/(mol*s)'), n=3.59921, Ea=(343.295,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CC_N-2Br1sCl1sF1sHI1s->H',), comment="""Estimated from node CC_N-2Br1sCl1sF1sHI1s->H"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CF4(90)', 'OC(F)(F)[C]F(576)'],
    products = ['OC(F)(F)C(F)(F)C(F)(F)F(1363)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.40188e+54,'m^3/(mol*s)'), n=-13.541, Ea=(208.711,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.513119107657762, var=99.27123869380007, Tref=1000.0, N=63, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction4',
    reactants = ['OC(F)(F)F(254)', 'F[C]C(F)(F)F(307)'],
    products = ['OC(F)(F)C(F)(F)C(F)(F)F(1363)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.05141e+54,'m^3/(mol*s)'), n=-13.541, Ea=(197.967,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.513119107657762, var=99.27123869380007, Tref=1000.0, N=63, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 3.0"""),
)

reaction(
    label = 'reaction5',
    reactants = ['OF(155)', 'F[C]C(F)(F)C(F)(F)F(1640)'],
    products = ['OC(F)(F)C(F)(F)C(F)(F)F(1363)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.13916e+09,'m^3/(mol*s)'), n=-1.00842, Ea=(11.3304,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17406617270594374, var=16.167420070870673, Tref=1000.0, N=4, data_mean=0.0, correlation='OHOY',), comment="""Estimated from node OHOY"""),
)

reaction(
    label = 'reaction6',
    reactants = ['OC(F)(F)C(F)(F)C(F)(F)F(1363)'],
    products = ['OC(F)(C(F)(F)F)C(F)(F)F(1641)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction7',
    reactants = ['OF(155)', 'FC(F)=C(F)C(F)(F)F(879)'],
    products = ['OC(F)(F)C(F)(F)C(F)(F)F(1363)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.000286045,'m^3/(mol*s)'), n=2.81783, Ea=(231.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OH] for rate rule [Cd/disub_Cd/disub;H_OH]
Euclidian distance = 1.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction8',
    reactants = ['OC(F)(F)F(254)', 'CF2CF2(60)'],
    products = ['OC(F)(F)C(F)(F)C(F)(F)F(1363)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.58e-11,'m^3/(mol*s)'), n=3.97, Ea=(329.281,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;R_OH] for rate rule [Cd/disub_Cd/disub;R_OH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F(37)', 'OC(F)(F)[C](F)C(F)(F)F(1279)'],
    products = ['OC(F)(F)C(F)(F)C(F)(F)F(1363)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+06,'m^3/(mol*s)'), n=0, Ea=(0.208499,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F(37)', 'O[C](F)C(F)(F)C(F)(F)F(1642)'],
    products = ['OC(F)(F)C(F)(F)C(F)(F)F(1363)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C
Ea raised from -0.9 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', 'OC(F)(F)C(F)(F)[C](F)F(1024)'],
    products = ['OC(F)(F)C(F)(F)C(F)(F)F(1363)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C"""),
)

reaction(
    label = 'reaction12',
    reactants = ['OH(4)', 'F[C](F)C(F)(F)C(F)(F)F(1643)'],
    products = ['OC(F)(F)C(F)(F)C(F)(F)F(1363)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.7e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(5)', '[O]C(F)(F)C(F)(F)C(F)(F)F(1644)'],
    products = ['OC(F)(F)C(F)(F)C(F)(F)F(1363)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction14',
    reactants = ['O[C](F)F(220)', 'F[C](F)C(F)(F)F(424)'],
    products = ['OC(F)(F)C(F)(F)C(F)(F)F(1363)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -13.8 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['CF3(44)', 'OC(F)(F)[C](F)F(586)'],
    products = ['OC(F)(F)C(F)(F)C(F)(F)F(1363)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -12.4 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['HF(38)', 'O=C(F)C(F)(F)C(F)(F)F(1645)'],
    products = ['OC(F)(F)C(F)(F)C(F)(F)F(1363)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.109156,'m^3/(mol*s)'), n=1.86531, Ea=(157.139,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd',), comment="""Estimated from node HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F2(108)', 'OC(F)=C(F)C(F)(F)F(1281)'],
    products = ['OC(F)(F)C(F)(F)C(F)(F)F(1363)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction18',
    reactants = ['F2(108)', 'OC(F)(F)C(F)=C(F)F(382)'],
    products = ['OC(F)(F)C(F)(F)C(F)(F)F(1363)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

network(
    label = 'PDepNetwork #521',
    isomers = [
        'OC(F)(F)C(F)(F)C(F)(F)F(1363)',
    ],
    reactants = [
        ('HF(38)', 'CF2O(48)', 'CF2CF2(60)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #521',
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

