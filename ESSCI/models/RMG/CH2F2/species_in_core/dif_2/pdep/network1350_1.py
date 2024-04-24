species(
    label = 'O=C(F)C(OF)C(O)(F)F(3725)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u0 p2 c0 {9,S} {12,S}
7  O u0 p2 c0 {10,D}
8  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
9  C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
10 C u0 p0 c0 {3,S} {7,D} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-1112.93,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,486,617,768,1157,1926,180,832.159,832.293],'cm^-1')),
        HinderedRotor(inertia=(0.222196,'amu*angstrom^2'), symmetry=1, barrier=(5.10873,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.755382,'amu*angstrom^2'), symmetry=1, barrier=(17.3677,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.755393,'amu*angstrom^2'), symmetry=1, barrier=(17.368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.755377,'amu*angstrom^2'), symmetry=1, barrier=(17.3676,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (162.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.238049,0.0863469,-0.000109184,6.77682e-08,-1.65912e-11,-133722,31.899], Tmin=(100,'K'), Tmax=(996.118,'K')), NASAPolynomial(coeffs=[16.4045,0.0214287,-1.14265e-05,2.34262e-09,-1.71014e-13,-136943,-46.0321], Tmin=(996.118,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1112.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCFFO) + group(COCsFO)"""),
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
    E0 = (-618.61,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.06775,-0.00614966,6.73615e-05,-1.17623e-07,6.56735e-11,-74400.6,7.68563], Tmin=(10,'K'), Tmax=(584.627,'K')), NASAPolynomial(coeffs=[3.15981,0.0116963,-8.27581e-06,2.66621e-09,-3.20585e-13,-74493.3,9.87819], Tmin=(584.627,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-618.61,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""ODC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'OC(F)(F)C(F)OF(1196)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {7,S}
6  O u0 p2 c0 {8,S} {10,S}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
8  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-954.131,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,261,493,600,1152,1365,1422,3097,223,363,546,575,694,1179,1410,216.058,1602.83],'cm^-1')),
        HinderedRotor(inertia=(0.391102,'amu*angstrom^2'), symmetry=1, barrier=(13.0237,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.496566,'amu*angstrom^2'), symmetry=1, barrier=(16.5278,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20014,'amu*angstrom^2'), symmetry=1, barrier=(39.2298,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.03,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3380.85,'J/mol'), sigma=(5.63196,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=528.08 K, Pc=42.94 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.50095,0.0463791,-2.66712e-05,-1.51866e-08,1.62125e-11,-114752,13.6584], Tmin=(10,'K'), Tmax=(730.091,'K')), NASAPolynomial(coeffs=[9.74283,0.0273007,-1.85374e-05,5.75107e-09,-6.69776e-13,-116066,-17.2511], Tmin=(730.091,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-954.131,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), label="""OC(F)(F)C(F)OF""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'O=C(F)C(O)OF(4126)',
    structure = adjacencyList("""1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {6,S} {9,S}
4 O u0 p2 c0 {2,S} {6,S}
5 O u0 p2 c0 {7,D}
6 C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
7 C u0 p0 c0 {1,S} {5,D} {6,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {3,S}
"""),
    E0 = (-657.296,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,1380,1390,370,380,2900,435,486,617,768,1157,1926,222.623,232.319],'cm^-1')),
        HinderedRotor(inertia=(0.00361791,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118904,'amu*angstrom^2'), symmetry=1, barrier=(3.86076,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.632893,'amu*angstrom^2'), symmetry=1, barrier=(21.3695,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70374,0.0545929,-7.11828e-05,4.86336e-08,-1.34331e-11,-78975.3,23.2397], Tmin=(100,'K'), Tmax=(877.769,'K')), NASAPolynomial(coeffs=[9.59361,0.0186384,-9.74064e-06,1.9679e-09,-1.42004e-13,-80360.4,-13.7956], Tmin=(877.769,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-657.296,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(Cs-CsOsOsH) + group(COCsFO)"""),
)

species(
    label = 'OF(153)',
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
    label = 'O=C(F)C([C]F)OF(4159)',
    structure = adjacencyList("""1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {8,S}
3 F u0 p3 c0 {4,S}
4 O u0 p2 c0 {3,S} {6,S}
5 O u0 p2 c0 {7,D}
6 C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
7 C u0 p0 c0 {1,S} {5,D} {6,S}
8 C u0 p1 c0 {2,S} {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-315.836,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,486,617,768,1157,1926,617,898,1187,180,186.661],'cm^-1')),
        HinderedRotor(inertia=(0.0049749,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.386151,'amu*angstrom^2'), symmetry=1, barrier=(8.87838,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.04582,'amu*angstrom^2'), symmetry=1, barrier=(47.0374,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29171,0.0644042,-9.37966e-05,7.09998e-08,-2.14832e-11,-37893.1,25.7045], Tmin=(100,'K'), Tmax=(807.9,'K')), NASAPolynomial(coeffs=[10.3782,0.0194121,-1.02538e-05,2.0552e-09,-1.46801e-13,-39361.2,-16.1933], Tmin=(807.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-315.836,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsOsH) + group(COCsFO) + group(CJ2_singlet-FCs)"""),
)

species(
    label = 'O=C(F)C(F)C(O)(F)OF(4142)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {6,S}
5  O u0 p2 c0 {9,S} {12,S}
6  O u0 p2 c0 {4,S} {9,S}
7  O u0 p2 c0 {10,D}
8  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
9  C u0 p0 c0 {2,S} {5,S} {6,S} {8,S}
10 C u0 p0 c0 {3,S} {7,D} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-1087.85,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,206,363,518,1175,1317,1481,3059,375,361,584,565,722,1474,486,617,768,1157,1926,209.678,210.58,210.962],'cm^-1')),
        HinderedRotor(inertia=(0.275595,'amu*angstrom^2'), symmetry=1, barrier=(8.64878,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.778759,'amu*angstrom^2'), symmetry=1, barrier=(24.5018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.273024,'amu*angstrom^2'), symmetry=1, barrier=(8.65283,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.995559,'amu*angstrom^2'), symmetry=1, barrier=(31.0404,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (162.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.200543,0.10183,-0.000159527,1.2346e-07,-3.59393e-11,-130696,30.061], Tmin=(100,'K'), Tmax=(651.189,'K')), NASAPolynomial(coeffs=[13.7749,0.0275819,-1.5213e-05,3.06577e-09,-2.1816e-13,-132762,-33.2554], Tmin=(651.189,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1087.85,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(COCsFO)"""),
)

species(
    label = 'FOF(716)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 O u0 p2 c0 {1,S} {2,S}
"""),
    E0 = (16.0257,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(53.9917,'amu')),
        NonlinearRotor(inertia=([8.34135,45.8552,54.1965],'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([486.983,915.713,1034.58],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (53.9962,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.03449,-0.00336126,4.20188e-05,-7.82672e-08,4.57369e-11,1928.28,6.37844], Tmin=(10,'K'), Tmax=(574.663,'K')), NASAPolynomial(coeffs=[3.91346,0.00598694,-4.58407e-06,1.5533e-09,-1.93094e-13,1801.75,5.67331], Tmin=(574.663,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(16.0257,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""FOF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C=CC(O)(F)F(3369)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {5,S} {9,S}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
6 C u0 p0 c0 {5,S} {7,D} {8,S}
7 C u0 p0 c0 {4,D} {6,D}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {3,S}
"""),
    E0 = (-695.014,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,275,321,533,585,746,850,1103,3010,987.5,1337.5,450,1655,2120,512.5,787.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.415431,'amu*angstrom^2'), symmetry=1, barrier=(9.55158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.415597,'amu*angstrom^2'), symmetry=1, barrier=(9.55539,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56816,0.0574135,-8.03592e-05,5.90387e-08,-1.73616e-11,-83506.8,20.2138], Tmin=(100,'K'), Tmax=(830.473,'K')), NASAPolynomial(coeffs=[9.76572,0.0179293,-9.04214e-06,1.78791e-09,-1.27016e-13,-84868.4,-17.8118], Tmin=(830.473,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-695.014,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + missing(O2d-Cdd) + group(CsCFFO) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d)"""),
)

species(
    label = 'O=C(F)C=C(O)F(4160)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {6,S} {9,S}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {6,D} {7,S} {8,S}
6 C u0 p0 c0 {1,S} {3,S} {5,D}
7 C u0 p0 c0 {2,S} {4,D} {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {3,S}
"""),
    E0 = (-714.052,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,326,540,652,719,1357,255,533,799,832,1228,180],'cm^-1')),
        HinderedRotor(inertia=(0.469109,'amu*angstrom^2'), symmetry=1, barrier=(10.7857,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.525757,'amu*angstrom^2'), symmetry=1, barrier=(12.0882,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56119,0.05468,-6.62762e-05,4.03038e-08,-9.67269e-12,-85793.7,19.8412], Tmin=(100,'K'), Tmax=(1017.19,'K')), NASAPolynomial(coeffs=[11.8138,0.0143621,-6.82078e-06,1.3363e-09,-9.53486e-14,-87879.4,-29.7966], Tmin=(1017.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-714.052,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cd-Cd(CO)H) + group(CdCFO) + group(COCFO)"""),
)

species(
    label = 'H2O(3)',
    structure = adjacencyList("""1 O u0 p2 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (-251.755,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1601.58,3620.23,4000],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (18.0153,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4759.22,'J/mol'), sigma=(2.605,'angstroms'), dipoleMoment=(1.844,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=4.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05764,-0.000787936,2.90877e-06,-1.47519e-09,2.12843e-13,-30281.6,-0.311364], Tmin=(100,'K'), Tmax=(1130.24,'K')), NASAPolynomial(coeffs=[2.84325,0.00275109,-7.81031e-07,1.07244e-10,-5.79392e-15,-29958.6,5.91042], Tmin=(1130.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-251.755,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""H2O""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'O=C(F)C(OF)=C(F)F(4161)',
    structure = adjacencyList("""1 F u0 p3 c0 {8,S}
2 F u0 p3 c0 {8,S}
3 F u0 p3 c0 {9,S}
4 F u0 p3 c0 {5,S}
5 O u0 p2 c0 {4,S} {7,S}
6 O u0 p2 c0 {9,D}
7 C u0 p0 c0 {5,S} {8,D} {9,S}
8 C u0 p0 c0 {1,S} {2,S} {7,D}
9 C u0 p0 c0 {3,S} {6,D} {7,S}
"""),
    E0 = (-726.202,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,350,440,435,1725,182,240,577,636,1210,1413,255,533,799,832,1228,180,1276.3],'cm^-1')),
        HinderedRotor(inertia=(0.143159,'amu*angstrom^2'), symmetry=1, barrier=(3.2915,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144311,'amu*angstrom^2'), symmetry=1, barrier=(3.318,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (144.024,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.739981,0.0819755,-0.000151572,1.41482e-07,-5.04932e-11,-87234.4,25.1831], Tmin=(100,'K'), Tmax=(832.075,'K')), NASAPolynomial(coeffs=[8.25926,0.0258929,-1.45329e-05,2.89047e-09,-2.01732e-13,-87795.6,-5.56381], Tmin=(832.075,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-726.202,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cds-Cds(Cds-O2d)O2s) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + group(COCFO)"""),
)

species(
    label = 'OC(F)(F)C=C(F)OOF(4162)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {7,S}
5  O u0 p2 c0 {8,S} {12,S}
6  O u0 p2 c0 {7,S} {10,S}
7  O u0 p2 c0 {4,S} {6,S}
8  C u0 p0 c0 {1,S} {2,S} {5,S} {9,S}
9  C u0 p0 c0 {8,S} {10,D} {11,S}
10 C u0 p0 c0 {3,S} {6,S} {9,D}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-788.322,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3615,1310,387.5,850,1000,277,555,632,275,321,533,585,746,850,1103,3010,987.5,1337.5,450,1655,326,540,652,719,1357],'cm^-1')),
        HinderedRotor(inertia=(0.634981,'amu*angstrom^2'), symmetry=1, barrier=(14.5995,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.635802,'amu*angstrom^2'), symmetry=1, barrier=(14.6183,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (162.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.212802,0.104201,-0.00018502,1.69078e-07,-5.98418e-11,-94672.7,32.0567], Tmin=(100,'K'), Tmax=(818.173,'K')), NASAPolynomial(coeffs=[10.2926,0.0330455,-1.82763e-05,3.63991e-09,-2.55156e-13,-95729.2,-12.4685], Tmin=(818.173,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-788.322,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2sFO) + group(CsCFFO) + group(Cds-CdsCsH) + group(CdCFO)"""),
)

species(
    label = 'OC(F)(F)OC(F)=COF(4163)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {7,S}
5  O u0 p2 c0 {8,S} {9,S}
6  O u0 p2 c0 {8,S} {12,S}
7  O u0 p2 c0 {4,S} {10,S}
8  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
9  C u0 p0 c0 {3,S} {5,S} {10,D}
10 C u0 p0 c0 {7,S} {9,D} {11,S}
11 H u0 p0 c0 {10,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-1017.27,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,231,791,435,565,619,662,854,1178,1396,326,540,652,719,1357,3010,987.5,1337.5,450,1655,261.632,261.696,261.91,262.177],'cm^-1')),
        HinderedRotor(inertia=(0.266834,'amu*angstrom^2'), symmetry=1, barrier=(12.9842,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.266144,'amu*angstrom^2'), symmetry=1, barrier=(12.9823,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.680824,'amu*angstrom^2'), symmetry=1, barrier=(33.1621,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.53446,'amu*angstrom^2'), symmetry=1, barrier=(26.059,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (162.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.147309,0.0958631,-0.000130556,8.69866e-08,-2.27331e-11,-122204,29.1851], Tmin=(100,'K'), Tmax=(937.86,'K')), NASAPolynomial(coeffs=[17.3113,0.0214038,-1.14706e-05,2.33873e-09,-1.69726e-13,-125479,-53.9232], Tmin=(937.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1017.27,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2sCF) + group(CsFFOO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH)"""),
)

species(
    label = '[O]C(F)[C](OF)C(O)(F)F(4164)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {6,S}
5  O u0 p2 c0 {8,S} {12,S}
6  O u0 p2 c0 {4,S} {10,S}
7  O u1 p2 c0 {9,S}
8  C u0 p0 c0 {1,S} {2,S} {5,S} {10,S}
9  C u0 p0 c0 {3,S} {7,S} {10,S} {11,S}
10 C u1 p0 c0 {6,S} {8,S} {9,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-757.588,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,253,525,597,667,842,1178,1324,391,562,707,872,1109,1210,1289,3137,360,370,350,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.305258,'amu*angstrom^2'), symmetry=1, barrier=(7.01849,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.305789,'amu*angstrom^2'), symmetry=1, barrier=(7.03068,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.305636,'amu*angstrom^2'), symmetry=1, barrier=(7.02717,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.305196,'amu*angstrom^2'), symmetry=1, barrier=(7.01705,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (162.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.614444,0.113059,-0.000203544,1.80639e-07,-6.12706e-11,-90961.7,33.7276], Tmin=(100,'K'), Tmax=(847.961,'K')), NASAPolynomial(coeffs=[13.034,0.0275358,-1.48612e-05,2.89521e-09,-1.99003e-13,-92516.3,-25.3859], Tmin=(847.961,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-757.588,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCFHO) + radical(O2sj(Cs-CsF1sH)) + radical(C2CsJO)"""),
)

species(
    label = '[O]C(F)C(OF)C([O])(F)F(4165)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u1 p2 c0 {9,S}
7  O u1 p2 c0 {10,S}
8  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
9  C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
10 C u0 p0 c0 {3,S} {7,S} {8,S} {12,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-677.725,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,351,323,533,609,664,892,1120,1201,391,562,707,872,1109,1210,1289,3137,202.707,202.72,202.77],'cm^-1')),
        HinderedRotor(inertia=(0.262352,'amu*angstrom^2'), symmetry=1, barrier=(7.65632,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.262518,'amu*angstrom^2'), symmetry=1, barrier=(7.6562,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.262513,'amu*angstrom^2'), symmetry=1, barrier=(7.65612,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (162.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.287545,0.10158,-0.000162241,1.31852e-07,-4.22517e-11,-81363.9,32.7298], Tmin=(100,'K'), Tmax=(767.813,'K')), NASAPolynomial(coeffs=[14.3312,0.0254206,-1.34551e-05,2.66485e-09,-1.87628e-13,-83608.8,-33.9347], Tmin=(767.813,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-677.725,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCFHO) + radical(O2sj(Cs-CsF1sF1s)) + radical(O2sj(Cs-CsF1sH))"""),
)

species(
    label = '[O]C(F)(F)C(OF)[C](O)F(4166)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u0 p2 c0 {10,S} {12,S}
7  O u1 p2 c0 {9,S}
8  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
9  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
10 C u1 p0 c0 {3,S} {6,S} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-709.744,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,1380,1390,370,380,2900,435,351,323,533,609,664,892,1120,1201,395,473,707,1436,192.846,192.901,192.923],'cm^-1')),
        HinderedRotor(inertia=(0.264042,'amu*angstrom^2'), symmetry=1, barrier=(6.97022,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.264021,'amu*angstrom^2'), symmetry=1, barrier=(6.97026,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.263947,'amu*angstrom^2'), symmetry=1, barrier=(6.96991,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.54576,'amu*angstrom^2'), symmetry=1, barrier=(14.4065,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (162.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.715978,0.114691,-0.000204685,1.79796e-07,-6.06157e-11,-85203.2,34.1603], Tmin=(100,'K'), Tmax=(839.348,'K')), NASAPolynomial(coeffs=[13.9545,0.0265283,-1.45182e-05,2.84896e-09,-1.9689e-13,-87023.1,-30.2176], Tmin=(839.348,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-709.744,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCFHO) + radical(O2sj(Cs-CsF1sF1s)) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'OC(F)=C(OF)C(O)(F)F(4167)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {6,S}
5  O u0 p2 c0 {8,S} {11,S}
6  O u0 p2 c0 {4,S} {9,S}
7  O u0 p2 c0 {10,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {5,S} {9,S}
9  C u0 p0 c0 {6,S} {8,S} {10,D}
10 C u0 p0 c0 {3,S} {7,S} {9,D}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-1017.17,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,231,791,275,321,533,585,746,850,1103,350,440,435,1725,326,540,652,719,1357,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.227511,'amu*angstrom^2'), symmetry=1, barrier=(5.23092,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.22749,'amu*angstrom^2'), symmetry=1, barrier=(5.23045,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.227518,'amu*angstrom^2'), symmetry=1, barrier=(5.23109,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.74186,'amu*angstrom^2'), symmetry=1, barrier=(40.0488,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (162.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.650737,0.116568,-0.000217699,2.00787e-07,-7.03824e-11,-122184,30.3309], Tmin=(100,'K'), Tmax=(844.077,'K')), NASAPolynomial(coeffs=[11.1424,0.0324758,-1.81363e-05,3.58053e-09,-2.47995e-13,-123170,-18.6133], Tmin=(844.077,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1017.17,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(O2s-(Cds-Cd)H) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(Cds-CdsCsOs) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs)"""),
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
    label = 'O=C(F)C(OF)[C](O)F(4168)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {7,S}
5  O u0 p2 c0 {8,S} {11,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
8  C u1 p0 c0 {1,S} {5,S} {7,S}
9  C u0 p0 c0 {2,S} {6,D} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-683.672,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,1380,1390,370,380,2900,435,395,473,707,1436,486,617,768,1157,1926,316.372,326.086,329.648],'cm^-1')),
        HinderedRotor(inertia=(0.0645098,'amu*angstrom^2'), symmetry=1, barrier=(5.11265,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0718809,'amu*angstrom^2'), symmetry=1, barrier=(5.19245,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180683,'amu*angstrom^2'), symmetry=1, barrier=(13.0472,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.31618,'amu*angstrom^2'), symmetry=1, barrier=(23.1728,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (143.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.540458,0.0820424,-0.00011889,8.7601e-08,-2.56336e-11,-82107.5,31.4106], Tmin=(100,'K'), Tmax=(836.262,'K')), NASAPolynomial(coeffs=[12.9922,0.0224845,-1.20637e-05,2.44055e-09,-1.75452e-13,-84190.1,-26.4356], Tmin=(836.262,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-683.672,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCFHO) + group(COCsFO) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'O=[C]C(OF)C(O)(F)F(4169)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {7,S}
5  O u0 p2 c0 {8,S} {11,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
9  C u1 p0 c0 {6,D} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-698.98,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,1855,455,950,292.993,296.862],'cm^-1')),
        HinderedRotor(inertia=(0.241162,'amu*angstrom^2'), symmetry=1, barrier=(14.8693,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00190718,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.241974,'amu*angstrom^2'), symmetry=1, barrier=(14.877,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.241405,'amu*angstrom^2'), symmetry=1, barrier=(14.87,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (143.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.439524,0.0802829,-0.000105369,6.68343e-08,-1.65105e-11,-83941.2,31.0653], Tmin=(100,'K'), Tmax=(994.195,'K')), NASAPolynomial(coeffs=[16.4528,0.0158557,-8.16431e-06,1.65246e-09,-1.19927e-13,-87125.2,-46.0963], Tmin=(994.195,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-698.98,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCFFO) + group(Cds-OdCsH) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C(C(=O)F)C(O)(F)F(4170)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {8,S} {11,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
9  C u0 p0 c0 {3,S} {6,D} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-1003.71,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,486,617,768,1157,1926,355.602,355.602,355.608],'cm^-1')),
        HinderedRotor(inertia=(0.00133309,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0689286,'amu*angstrom^2'), symmetry=1, barrier=(6.1847,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0689296,'amu*angstrom^2'), symmetry=1, barrier=(6.18487,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (143.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.858235,0.0722455,-9.06484e-05,5.70748e-08,-1.42506e-11,-120608,30.2797], Tmin=(100,'K'), Tmax=(975.736,'K')), NASAPolynomial(coeffs=[13.6431,0.0198348,-1.00782e-05,2.02621e-09,-1.46333e-13,-123103,-31.0859], Tmin=(975.736,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1003.71,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCFFO) + group(COCsFO) + radical(C=OCOJ)"""),
)

species(
    label = '[O]F(154)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 O u1 p2 c0 {1,S}
"""),
    E0 = (102.686,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(34.9933,'amu')),
        LinearRotor(inertia=(15.6231,'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([1151.69],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (34.9978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51779,-0.00115859,7.7708e-06,-9.55983e-09,3.74809e-12,12350.1,5.42233], Tmin=(10,'K'), Tmax=(823.913,'K')), NASAPolynomial(coeffs=[3.16214,0.00226288,-1.54389e-06,4.73852e-10,-5.4015e-14,12351.2,6.72012], Tmin=(823.913,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(102.686,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""[O]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C(F)[CH]C(O)(F)F(4171)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {10,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
7  C u1 p0 c0 {6,S} {8,S} {9,S}
8  C u0 p0 c0 {3,S} {5,D} {7,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (-877.657,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,253,525,597,667,842,1178,1324,3025,407.5,1350,352.5,611,648,830,1210,1753,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.204998,'amu*angstrom^2'), symmetry=1, barrier=(4.7133,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20498,'amu*angstrom^2'), symmetry=1, barrier=(4.7129,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204974,'amu*angstrom^2'), symmetry=1, barrier=(4.71276,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19834,0.067655,-9.94818e-05,7.61687e-08,-2.25237e-11,-105463,25.9667], Tmin=(100,'K'), Tmax=(657.9,'K')), NASAPolynomial(coeffs=[9.65114,0.0227443,-1.18653e-05,2.36036e-09,-1.6749e-13,-106715,-12.3401], Tmin=(657.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-877.657,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsHH) + group(CsCFFO) + group(COCsFO) + radical(CCJCO)"""),
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
    label = 'O=C(F)C(OF)[C](F)F(1764)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {7,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
8  C u1 p0 c0 {2,S} {3,S} {7,S}
9  C u0 p0 c0 {1,S} {6,D} {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-700.682,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,190,488,555,1236,1407,486,617,768,1157,1926,212.217,212.396,691.702],'cm^-1')),
        HinderedRotor(inertia=(0.522,'amu*angstrom^2'), symmetry=1, barrier=(16.7259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00374253,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.699945,'amu*angstrom^2'), symmetry=1, barrier=(22.3809,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (145.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.853722,0.0744802,-0.000101767,6.95982e-08,-1.89229e-11,-84163.8,29.2715], Tmin=(100,'K'), Tmax=(896.883,'K')), NASAPolynomial(coeffs=[12.9879,0.0203625,-1.12562e-05,2.31905e-09,-1.69008e-13,-86340.3,-27.948], Tmin=(896.883,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-700.682,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-O2d)CsOsH) + group(CsCsFFH) + group(COCsFO) + radical(CsCsF1sF1s)"""),
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
    label = '[O]C(F)(F)C(OF)C(=O)F(4172)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u1 p2 c0 {9,S}
7  O u0 p2 c0 {10,D}
8  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
9  C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
10 C u0 p0 c0 {3,S} {7,D} {8,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-852.965,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,351,323,533,609,664,892,1120,1201,486,617,768,1157,1926,180,180,1018.62],'cm^-1')),
        HinderedRotor(inertia=(0.0890477,'amu*angstrom^2'), symmetry=1, barrier=(2.04738,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.673973,'amu*angstrom^2'), symmetry=1, barrier=(15.496,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.671758,'amu*angstrom^2'), symmetry=1, barrier=(15.445,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (161.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.469572,0.0821205,-0.000106978,6.85788e-08,-1.73781e-11,-102465,31.6593], Tmin=(100,'K'), Tmax=(962.237,'K')), NASAPolynomial(coeffs=[15.2246,0.0207886,-1.13759e-05,2.34777e-09,-1.71733e-13,-105304,-38.9581], Tmin=(962.237,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-852.965,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCFFO) + group(COCsFO) + radical(O2sj(Cs-CsF1sF1s))"""),
)

species(
    label = 'O[C](F)F(232)',
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
    label = 'O=C(F)[CH]OF(2999)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {3,S}
3 O u0 p2 c0 {2,S} {5,S}
4 O u0 p2 c0 {6,D}
5 C u1 p0 c0 {3,S} {6,S} {7,S}
6 C u0 p0 c0 {1,S} {4,D} {5,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-323.826,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3025,407.5,1350,352.5,611,648,830,1210,1753,180,783.373],'cm^-1')),
        HinderedRotor(inertia=(0.41351,'amu*angstrom^2'), symmetry=1, barrier=(9.50742,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.48155,'amu*angstrom^2'), symmetry=1, barrier=(57.0558,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.0249,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.01266,0.0470476,-7.14345e-05,5.69319e-08,-1.80482e-11,-38878.9,18.3743], Tmin=(100,'K'), Tmax=(773.881,'K')), NASAPolynomial(coeffs=[8.45194,0.0137648,-6.92359e-06,1.35885e-09,-9.56141e-14,-39875.5,-11.0409], Tmin=(773.881,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-323.826,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(OCJC=O)"""),
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
    label = 'OC(F)(F)[CH]OF(1227)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {5,S}
4 O u0 p2 c0 {6,S} {9,S}
5 O u0 p2 c0 {3,S} {7,S}
6 C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
7 C u1 p0 c0 {5,S} {6,S} {8,S}
8 H u0 p0 c0 {7,S}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (-571.479,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,253,525,597,667,842,1178,1324,3025,407.5,1350,352.5,200.971,201.785],'cm^-1')),
        HinderedRotor(inertia=(0.733687,'amu*angstrom^2'), symmetry=1, barrier=(20.6575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.331115,'amu*angstrom^2'), symmetry=1, barrier=(9.56013,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.591119,'amu*angstrom^2'), symmetry=1, barrier=(17.0053,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.995105,0.0707693,-0.000107109,7.97536e-08,-2.32357e-11,-68629,22.2258], Tmin=(100,'K'), Tmax=(844.494,'K')), NASAPolynomial(coeffs=[12.7536,0.0150762,-8.18967e-06,1.66667e-09,-1.20036e-13,-70615,-32.5153], Tmin=(844.494,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-571.479,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCsJO)"""),
)

species(
    label = 'O=C(F)[C](OF)C(O)(F)F(4173)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {6,S}
5  O u0 p2 c0 {8,S} {11,S}
6  O u0 p2 c0 {4,S} {9,S}
7  O u0 p2 c0 {10,D}
8  C u0 p0 c0 {1,S} {2,S} {5,S} {9,S}
9  C u1 p0 c0 {6,S} {8,S} {10,S}
10 C u0 p0 c0 {3,S} {7,D} {9,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-932.828,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,253,525,597,667,842,1178,1324,360,370,350,611,648,830,1210,1753,219.626,219.626,219.626],'cm^-1')),
        HinderedRotor(inertia=(0.147128,'amu*angstrom^2'), symmetry=1, barrier=(5.03606,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147128,'amu*angstrom^2'), symmetry=1, barrier=(5.03606,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147129,'amu*angstrom^2'), symmetry=1, barrier=(5.03606,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.575106,'amu*angstrom^2'), symmetry=1, barrier=(19.6853,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (161.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00635981,0.095242,-0.00015428,1.25801e-07,-4.04485e-11,-112056,33.1444], Tmin=(100,'K'), Tmax=(765.034,'K')), NASAPolynomial(coeffs=[13.7858,0.0231798,-1.29572e-05,2.62246e-09,-1.86977e-13,-114164,-29.64], Tmin=(765.034,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-932.828,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCFFO) + group(COCsFO) + radical(C2CsJO)"""),
)

species(
    label = 'OC(F)=COF(1234)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {5,S} {8,S}
4 O u0 p2 c0 {2,S} {6,S}
5 C u0 p0 c0 {1,S} {3,S} {6,D}
6 C u0 p0 c0 {4,S} {5,D} {7,S}
7 H u0 p0 c0 {6,S}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (-376.042,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,231,791,326,540,652,719,1357,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(0.995017,'amu*angstrom^2'), symmetry=1, barrier=(22.8774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.996095,'amu*angstrom^2'), symmetry=1, barrier=(22.9022,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0328,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72754,0.0499513,-6.12171e-05,3.67588e-08,-8.61435e-12,-45145.5,17.7213], Tmin=(100,'K'), Tmax=(1046.64,'K')), NASAPolynomial(coeffs=[11.9808,0.0107666,-5.06019e-06,9.89728e-10,-7.07151e-14,-47291.8,-32.2121], Tmin=(1046.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-376.042,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2sCF) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH)"""),
)

species(
    label = 'O=C(F)C(OF)C(=O)F(4174)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {7,S}
5  O u0 p2 c0 {8,D}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {1,S} {5,D} {7,S}
9  C u0 p0 c0 {2,S} {6,D} {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-828.49,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,430,542,558,676,687,849,1103,1211,1906,1946,286.154,292.228,2887.39],'cm^-1')),
        HinderedRotor(inertia=(0.153044,'amu*angstrom^2'), symmetry=1, barrier=(9.25421,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.565526,'amu*angstrom^2'), symmetry=1, barrier=(33.9413,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.564632,'amu*angstrom^2'), symmetry=1, barrier=(33.9445,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.201,0.0670854,-8.89421e-05,6.21929e-08,-1.7662e-11,-99548.4,26.9311], Tmin=(100,'K'), Tmax=(852.866,'K')), NASAPolynomial(coeffs=[10.4463,0.0237235,-1.26766e-05,2.5766e-09,-1.86352e-13,-101125,-16.2004], Tmin=(852.866,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-828.49,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-CsCsOsH) + group(COCsFO) + group(COCsFO)"""),
)

species(
    label = 'O=C(F)C(=O)C(O)(F)F(3107)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {7,S} {10,S}
5  O u0 p2 c0 {8,D}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
8  C u0 p0 c0 {5,D} {7,S} {9,S}
9  C u0 p0 c0 {3,S} {6,D} {8,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (-1146.12,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,251,367,519,700,855,1175,1303,375,552.5,462.5,1710,286,619,818,1246,1924,243.592,243.607],'cm^-1')),
        HinderedRotor(inertia=(0.168168,'amu*angstrom^2'), symmetry=1, barrier=(7.08401,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168245,'amu*angstrom^2'), symmetry=1, barrier=(7.08495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.583608,'amu*angstrom^2'), symmetry=1, barrier=(24.5817,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.615836,0.0804929,-0.000126656,1.01562e-07,-3.22185e-11,-137730,26.3539], Tmin=(100,'K'), Tmax=(774.354,'K')), NASAPolynomial(coeffs=[12.1398,0.020962,-1.13334e-05,2.27231e-09,-1.61471e-13,-139515,-26.2952], Tmin=(774.354,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1146.12,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cds-O2d(Cds-O2d)Cs) + group(COCFO)"""),
)

species(
    label = 'O=C(F)C(OF)=C(O)F(4175)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {7,S}
5  O u0 p2 c0 {8,S} {10,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {4,S} {8,D} {9,S}
8  C u0 p0 c0 {1,S} {5,S} {7,D}
9  C u0 p0 c0 {2,S} {6,D} {7,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-736.01,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,3615,1277.5,1000,350,440,435,1725,326,540,652,719,1357,255,533,799,832,1228,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.13582,'amu*angstrom^2'), symmetry=1, barrier=(3.12278,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134791,'amu*angstrom^2'), symmetry=1, barrier=(3.09912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.957395,'amu*angstrom^2'), symmetry=1, barrier=(22.0124,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.636,0.0841522,-0.000126193,7.84646e-08,-8.77389e-12,-88410.3,24.966], Tmin=(100,'K'), Tmax=(575.517,'K')), NASAPolynomial(coeffs=[12.309,0.0216634,-1.19126e-05,2.36613e-09,-1.6573e-13,-90062.6,-27.5825], Tmin=(575.517,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-736.01,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-O2d)O2s) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(COCFO)"""),
)

species(
    label = 'O=C=C(OF)C(O)(F)F(4089)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {5,S}
4  O u0 p2 c0 {7,S} {10,S}
5  O u0 p2 c0 {3,S} {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
8  C u0 p0 c0 {5,S} {7,S} {9,D}
9  C u0 p0 c0 {6,D} {8,D}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (-676.909,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,231,791,275,321,533,585,746,850,1103,350,440,435,1725,2120,512.5,787.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.316404,'amu*angstrom^2'), symmetry=1, barrier=(7.27474,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.99251,'amu*angstrom^2'), symmetry=1, barrier=(91.7956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.318254,'amu*angstrom^2'), symmetry=1, barrier=(7.31727,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0190321,0.0999273,-0.000186827,1.68471e-07,-5.73275e-11,-81279.6,26.1512], Tmin=(100,'K'), Tmax=(865.186,'K')), NASAPolynomial(coeffs=[11.4267,0.0232193,-1.25886e-05,2.42887e-09,-1.64946e-13,-82369.7,-22.2646], Tmin=(865.186,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-676.909,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + missing(O2d-Cdd) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(Cds-(Cdd-O2d)CsOs) + missing(Cdd-CdO2d)"""),
)

species(
    label = 'O=C=COF(1390)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {4,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {2,S} {5,D} {6,S}
5 C u0 p0 c0 {3,D} {4,D}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-81.3175,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,3010,987.5,1337.5,450,1655,2120,512.5,787.5,499.481],'cm^-1')),
        HinderedRotor(inertia=(1.57797,'amu*angstrom^2'), symmetry=1, barrier=(36.2806,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (76.0265,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.89204,0.00775819,8.52992e-05,-2.45346e-07,2.02193e-10,-9777.47,10.2957], Tmin=(10,'K'), Tmax=(427.3,'K')), NASAPolynomial(coeffs=[4.71619,0.0177834,-1.21691e-05,3.88442e-09,-4.70166e-13,-10009.9,5.12526], Tmin=(427.3,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-81.3175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""ODCDCOF""", comment="""Thermo library: CHOF_G4"""),
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
    E0 = (-569.953,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-203.842,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-151.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (103.407,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-284.816,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (8.60448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (37.8028,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-242.127,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-183.454,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-396.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-231.251,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-148.716,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-180.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-337.948,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-106.745,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-122.052,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-426.78,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-270.935,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-168.251,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-137.125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-283.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-257.802,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-216.988,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-399.405,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-438.41,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-476.329,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-332.171,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-261.409,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-351.23,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=C(F)C(OF)C(O)(F)F(3725)'],
    products = ['HF(38)', 'CF2O(48)', 'O=CC(=O)F(557)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(0.00285451,'s^-1'), n=4.62568, Ea=(38.9447,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(12)', 'OC(F)(F)C(F)OF(1196)'],
    products = ['O=C(F)C(OF)C(O)(F)F(3725)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.0026956,'m^3/(mol*s)'), n=2.93313, Ea=(364.994,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_2Br1sCl1sF1sH->F1s',), comment="""Estimated from node Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_2Br1sCl1sF1sH->F1s"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CF2(87)', 'O=C(F)C(O)OF(4126)'],
    products = ['O=C(F)C(OF)C(O)(F)F(3725)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.40908e-07,'m^3/(mol*s)'), n=3.45227, Ea=(205.693,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s',), comment="""Estimated from node CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction4',
    reactants = ['OF(153)', 'O=C(F)C([C]F)OF(4159)'],
    products = ['O=C(F)C(OF)C(O)(F)F(3725)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.13916e+09,'m^3/(mol*s)'), n=-1.00842, Ea=(10.4717,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17406617270594374, var=16.167420070870673, Tref=1000.0, N=4, data_mean=0.0, correlation='OHOY',), comment="""Estimated from node OHOY"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O=C(F)C(F)C(O)(F)OF(4142)'],
    products = ['O=C(F)C(OF)C(O)(F)F(3725)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction6',
    reactants = ['FOF(716)', 'O=C=CC(O)(F)F(3369)'],
    products = ['O=C(F)C(OF)C(O)(F)F(3725)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.75748e-06,'m^3/(mol*s)'), n=3.49511, Ea=(183.557,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [doublebond;H_OH] for rate rule [Cdd_Cd_HNd;H_OH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction7',
    reactants = ['FOF(716)', 'O=C(F)C=C(O)F(4160)'],
    products = ['O=C(F)C(OF)C(O)(F)F(3725)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00057209,'m^3/(mol*s)'), n=2.81783, Ea=(231.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OH] for rate rule [Cd/disub_Cd/monosub;H_OH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H2O(3)', 'O=C(F)C(OF)=C(F)F(4161)'],
    products = ['O=C(F)C(OF)C(O)(F)F(3725)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.00057209,'m^3/(mol*s)'), n=2.81783, Ea=(231.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OH] for rate rule [Cd/disub_Cd/disub;H_OH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction9',
    reactants = ['OC(F)(F)C=C(F)OOF(4162)'],
    products = ['O=C(F)C(OF)C(O)(F)F(3725)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(100.833,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction10',
    reactants = ['OC(F)(F)OC(F)=COF(4163)'],
    products = ['O=C(F)C(OF)C(O)(F)F(3725)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(116.647,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]C(F)[C](OF)C(O)(F)F(4164)'],
    products = ['O=C(F)C(OF)C(O)(F)F(3725)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C(F)C(OF)C([O])(F)F(4165)'],
    products = ['O=C(F)C(OF)C(O)(F)F(3725)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]C(F)(F)C(OF)[C](O)F(4166)'],
    products = ['O=C(F)C(OF)C(O)(F)F(3725)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['OC(F)=C(OF)C(O)(F)F(4167)'],
    products = ['O=C(F)C(OF)C(O)(F)F(3725)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(175.189,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_Ext-1R!H-R',), comment="""Estimated from node Root_3R!H->C_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F(37)', 'O=C(F)C(OF)[C](O)F(4168)'],
    products = ['O=C(F)C(OF)C(O)(F)F(3725)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F(37)', 'O=[C]C(OF)C(O)(F)F(4169)'],
    products = ['O=C(F)C(OF)C(O)(F)F(3725)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.52887e+10,'m^3/(mol*s)'), n=-0.421056, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.030895812897821735, var=3.1393582276389975, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F(37)', '[O]C(C(=O)F)C(O)(F)F(4170)'],
    products = ['O=C(F)C(OF)C(O)(F)F(3725)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.55564e+07,'m^3/(mol*s)'), n=-5.63145e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.7121440562946592, var=2.9508506800589083, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]F(154)', 'O=C(F)[CH]C(O)(F)F(4171)'],
    products = ['O=C(F)C(OF)C(O)(F)F(3725)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction19',
    reactants = ['OH(4)', 'O=C(F)C(OF)[C](F)F(1764)'],
    products = ['O=C(F)C(OF)C(O)(F)F(3725)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.7e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(5)', '[O]C(F)(F)C(OF)C(=O)F(4172)'],
    products = ['O=C(F)C(OF)C(O)(F)F(3725)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction21',
    reactants = ['O[C](F)F(232)', 'O=C(F)[CH]OF(2999)'],
    products = ['O=C(F)C(OF)C(O)(F)F(3725)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(6.86081,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CFO(50)', 'OC(F)(F)[CH]OF(1227)'],
    products = ['O=C(F)C(OF)C(O)(F)F(3725)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -3.5 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(5)', 'O=C(F)[C](OF)C(O)(F)F(4173)'],
    products = ['O=C(F)C(OF)C(O)(F)F(3725)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.7979e+07,'m^3/(mol*s)'), n=0.240345, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_N-Sp-3C=2C_Ext-4R!H-R_Ext-2C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_N-Sp-3C=2C_Ext-4R!H-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction24',
    reactants = ['O=C(F)C(OF)C(O)(F)F(3725)'],
    products = ['CF2O(48)', 'OC(F)=COF(1234)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(351.483,'s^-1'), n=2.72887, Ea=(209.492,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.13055717705123002, var=19.959654271360805, Tref=1000.0, N=68, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction25',
    reactants = ['HF(38)', 'O=C(F)C(OF)C(=O)F(4174)'],
    products = ['O=C(F)C(OF)C(O)(F)F(3725)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.218312,'m^3/(mol*s)'), n=1.86531, Ea=(167.157,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd',), comment="""Estimated from node HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction26',
    reactants = ['HF(38)', 'O=C(F)C(=O)C(O)(F)F(3107)'],
    products = ['O=C(F)C(OF)C(O)(F)F(3725)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(446.869,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction27',
    reactants = ['HF(38)', 'O=C(F)C(OF)=C(O)F(4175)'],
    products = ['O=C(F)C(OF)C(O)(F)F(3725)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.14111,'m^3/(mol*s)'), n=1.29695, Ea=(180.916,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction28',
    reactants = ['HF(38)', 'O=C=C(OF)C(O)(F)F(4089)'],
    products = ['O=C(F)C(OF)C(O)(F)F(3725)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(192.578,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction29',
    reactants = ['O=C(F)C(OF)C(O)(F)F(3725)'],
    products = ['HF(38)', 'CF2O(48)', 'O=C=COF(1390)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.98542e+10,'s^-1'), n=0.978818, Ea=(257.668,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R!H->C_N-5Br1sCl1sF1sH->H_5Br1sCl1sF1s->F1s_Ext-2C-R_N-7R!H->O_Sp-7BrCClFINPSSi-2C_Ext-2C-R',), comment="""Estimated from node Root_1R!H->C_N-5Br1sCl1sF1sH->H_5Br1sCl1sF1s->F1s_Ext-2C-R_N-7R!H->O_Sp-7BrCClFINPSSi-2C_Ext-2C-R"""),
)

network(
    label = 'PDepNetwork #1350',
    isomers = [
        'O=C(F)C(OF)C(O)(F)F(3725)',
    ],
    reactants = [
        ('HF(38)', 'CF2O(48)', 'O=CC(=O)F(557)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1350',
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

