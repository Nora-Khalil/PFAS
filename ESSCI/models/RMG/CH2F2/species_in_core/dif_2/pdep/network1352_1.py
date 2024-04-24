species(
    label = 'O=C(F)C(O)OC(F)(F)F(3727)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {8,S} {9,S}
6  O u0 p2 c0 {8,S} {12,S}
7  O u0 p2 c0 {10,D}
8  C u0 p0 c0 {5,S} {6,S} {10,S} {11,S}
9  C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
10 C u0 p0 c0 {4,S} {7,D} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-1451.19,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,431,527,613,668,835,1199,1245,1322,486,617,768,1157,1926,264.045,264.059,1972.53,1972.54],'cm^-1')),
        HinderedRotor(inertia=(0.195059,'amu*angstrom^2'), symmetry=1, barrier=(9.65146,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.534789,'amu*angstrom^2'), symmetry=1, barrier=(26.4606,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.195056,'amu*angstrom^2'), symmetry=1, barrier=(9.65164,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.195063,'amu*angstrom^2'), symmetry=1, barrier=(9.65147,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (162.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.605277,0.0812782,-0.000112568,8.31628e-08,-2.49005e-11,-174421,31.5443], Tmin=(100,'K'), Tmax=(811.713,'K')), NASAPolynomial(coeffs=[11.3203,0.0284752,-1.49893e-05,3.0199e-09,-2.16928e-13,-176161,-17.9143], Tmin=(811.713,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1451.19,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsFFFO) + group(COCsFO)"""),
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
    label = 'OC(F)OC(F)(F)F(1198)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {7,S} {10,S}
7  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
8  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-1320.24,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,509,613,660,1171,1360,1414,3084,431,527,613,668,835,1199,1245,1322,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.35669,'amu*angstrom^2'), symmetry=1, barrier=(8.20102,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.35659,'amu*angstrom^2'), symmetry=1, barrier=(8.1987,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.688734,'amu*angstrom^2'), symmetry=1, barrier=(15.8353,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.03,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3451.5,'J/mol'), sigma=(5.81127,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=539.12 K, Pc=39.91 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.66954,0.0244471,0.000162466,-5.2726e-07,4.53462e-10,-158788,13.5872], Tmin=(10,'K'), Tmax=(429.783,'K')), NASAPolynomial(coeffs=[7.67544,0.0325672,-2.4337e-05,8.30398e-09,-1.05007e-12,-159552,-7.2345], Tmin=(429.783,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1320.24,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), label="""OC(F)OC(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'FOC(F)(F)F(1140)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {6,S}
4 F u0 p3 c0 {5,S}
5 O u0 p2 c0 {4,S} {6,S}
6 C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
"""),
    E0 = (-757.394,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,431,527,613,668,835,1199,1245,1322,180],'cm^-1')),
        HinderedRotor(inertia=(0.86978,'amu*angstrom^2'), symmetry=1, barrier=(19.9979,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.004,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88199,0.0071794,8.93998e-05,-2.06387e-07,1.33773e-10,-91085.4,10.4612], Tmin=(10,'K'), Tmax=(550.265,'K')), NASAPolynomial(coeffs=[5.53551,0.0205748,-1.6396e-05,5.72482e-09,-7.29448e-13,-91652.1,-0.0246061], Tmin=(550.265,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-757.394,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""FOC(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C=CO(1378)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,S} {6,S}
2 O u0 p2 c0 {4,D}
3 C u0 p0 c0 {1,S} {4,D} {5,S}
4 C u0 p0 c0 {2,D} {3,D}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {1,S}
"""),
    E0 = (-165.709,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,2120,512.5,787.5],'cm^-1')),
        HinderedRotor(inertia=(0.80813,'amu*angstrom^2'), symmetry=1, barrier=(18.5805,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.86771,0.0228707,-2.2339e-05,1.13145e-08,-2.26468e-12,-19887.6,11.2463], Tmin=(100,'K'), Tmax=(1217.45,'K')), NASAPolynomial(coeffs=[7.78054,0.0067292,-2.45109e-06,4.23933e-10,-2.82988e-14,-21083.9,-13.4219], Tmin=(1217.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-165.709,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""hydroxyketene""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = 'O=C=COC(F)(F)F(4269)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {6,S}
4 O u0 p2 c0 {6,S} {7,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {2,S} {3,S} {4,S}
7 C u0 p0 c0 {4,S} {8,D} {9,S}
8 C u0 p0 c0 {5,D} {7,D}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-855.407,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([431,527,613,668,835,1199,1245,1322,3010,987.5,1337.5,450,1655,2120,512.5,787.5,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.976985,'amu*angstrom^2'), symmetry=1, barrier=(22.4628,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.973166,'amu*angstrom^2'), symmetry=1, barrier=(22.375,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.399127,0.0648266,-7.52066e-05,3.93818e-08,-7.60131e-12,-102740,23.6366], Tmin=(100,'K'), Tmax=(1445.82,'K')), NASAPolynomial(coeffs=[20.5574,-5.52391e-05,1.56012e-06,-3.74349e-10,2.67438e-14,-107616,-77.7535], Tmin=(1445.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-855.407,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + missing(O2d-Cdd) + group(CsFFFO) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d)"""),
)

species(
    label = 'OC=C(F)OOC(F)(F)F(4270)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {9,S}
7  O u0 p2 c0 {10,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
9  C u0 p0 c0 {4,S} {6,S} {10,D}
10 C u0 p0 c0 {7,S} {9,D} {11,S}
11 H u0 p0 c0 {10,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-1089.56,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1277.5,1000,431,527,613,668,835,1199,1245,1322,326,540,652,719,1357,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(0.583821,'amu*angstrom^2'), symmetry=1, barrier=(13.4232,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.584406,'amu*angstrom^2'), symmetry=1, barrier=(13.4366,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.18806,'amu*angstrom^2'), symmetry=1, barrier=(50.3078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.755434,'amu*angstrom^2'), symmetry=1, barrier=(17.3689,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (162.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.23853,0.0882295,-0.000116481,7.76137e-08,-2.05905e-11,-130913,30.104], Tmin=(100,'K'), Tmax=(918.304,'K')), NASAPolynomial(coeffs=[14.7405,0.0250594,-1.32936e-05,2.70043e-09,-1.95542e-13,-133576,-38.6234], Tmin=(918.304,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1089.56,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(CsFFFO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH)"""),
)

species(
    label = 'OOC(F)=COC(F)(F)F(4271)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {8,S} {9,S}
6  O u0 p2 c0 {7,S} {10,S}
7  O u0 p2 c0 {6,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
9  C u0 p0 c0 {5,S} {10,D} {11,S}
10 C u0 p0 c0 {4,S} {6,S} {9,D}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-1070.8,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,431,527,613,668,835,1199,1245,1322,3010,987.5,1337.5,450,1655,326,540,652,719,1357,234.09,234.12,234.252],'cm^-1')),
        HinderedRotor(inertia=(0.259556,'amu*angstrom^2'), symmetry=1, barrier=(10.1008,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.259473,'amu*angstrom^2'), symmetry=1, barrier=(10.0991,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.07712,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.07648,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (162.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.250897,0.0870057,-0.000111858,7.17403e-08,-1.82482e-11,-128656,29.4367], Tmin=(100,'K'), Tmax=(957.874,'K')), NASAPolynomial(coeffs=[15.4593,0.0234972,-1.24064e-05,2.52397e-09,-1.83213e-13,-131570,-43.2805], Tmin=(957.874,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1070.8,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsFFFO) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs)"""),
)

species(
    label = '[O]C(F)[C](O)OC(F)(F)F(4272)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {9,S} {10,S}
6  O u0 p2 c0 {10,S} {12,S}
7  O u1 p2 c0 {8,S}
8  C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
9  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
10 C u1 p0 c0 {5,S} {6,S} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-1069.9,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,391,562,707,872,1109,1210,1289,3137,431,527,613,668,835,1199,1245,1322,360,370,350,207.479,207.48,207.482,2167.6],'cm^-1')),
        HinderedRotor(inertia=(0.00391604,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.226148,'amu*angstrom^2'), symmetry=1, barrier=(6.90834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.226147,'amu*angstrom^2'), symmetry=1, barrier=(6.90834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.226148,'amu*angstrom^2'), symmetry=1, barrier=(6.90834,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (162.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.259896,0.0909693,-0.000136033,9.95168e-08,-2.59137e-11,-128553,33.4295], Tmin=(100,'K'), Tmax=(634.001,'K')), NASAPolynomial(coeffs=[12.4436,0.0267643,-1.40898e-05,2.79533e-09,-1.97399e-13,-130353,-21.8047], Tmin=(634.001,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1069.9,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCFHO) + group(CsFFFO) + radical(O2sj(Cs-CsF1sH)) + radical(Cs_P)"""),
)

species(
    label = '[O]C(F)C([O])OC(F)(F)F(4273)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {8,S} {10,S}
6  O u1 p2 c0 {8,S}
7  O u1 p2 c0 {9,S}
8  C u0 p0 c0 {5,S} {6,S} {9,S} {11,S}
9  C u0 p0 c0 {1,S} {7,S} {8,S} {12,S}
10 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-1049.44,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,391,562,707,872,1109,1210,1289,3137,431,527,613,668,835,1199,1245,1322,180,180,180,1815.18,1815.83],'cm^-1')),
        HinderedRotor(inertia=(0.26598,'amu*angstrom^2'), symmetry=1, barrier=(6.1154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.265886,'amu*angstrom^2'), symmetry=1, barrier=(6.11325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.266411,'amu*angstrom^2'), symmetry=1, barrier=(6.12531,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (162.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.30429,0.0914622,-0.000155677,1.43149e-07,-5.15758e-11,-126096,33.4045], Tmin=(100,'K'), Tmax=(808.574,'K')), NASAPolynomial(coeffs=[8.21141,0.0349978,-1.87466e-05,3.71555e-09,-2.60786e-13,-126807,0.444069], Tmin=(808.574,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1049.44,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCFHO) + group(CsFFFO) + radical(CCOJ) + radical(O2sj(Cs-CsF1sH))"""),
)

species(
    label = '[O]C(OC(F)(F)F)[C](O)F(4274)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {8,S} {9,S}
6  O u0 p2 c0 {10,S} {12,S}
7  O u1 p2 c0 {8,S}
8  C u0 p0 c0 {5,S} {7,S} {10,S} {11,S}
9  C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
10 C u1 p0 c0 {4,S} {6,S} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-1081.46,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,431,527,613,668,835,1199,1245,1322,395,473,707,1436,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (162.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0502537,0.103598,-0.000194168,1.85074e-07,-6.69139e-11,-129938,34.5764], Tmin=(100,'K'), Tmax=(841.099,'K')), NASAPolynomial(coeffs=[7.67155,0.0364026,-1.99895e-05,3.94368e-09,-2.73795e-13,-130159,5.06701], Tmin=(841.099,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1081.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCFHO) + group(CsFFFO) + radical(CCOJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'OC(F)=C(O)OC(F)(F)F(4275)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {8,S} {9,S}
6  O u0 p2 c0 {9,S} {12,S}
7  O u0 p2 c0 {10,S} {11,S}
8  C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
9  C u0 p0 c0 {5,S} {6,S} {10,D}
10 C u0 p0 c0 {4,S} {7,S} {9,D}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-1352.43,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,431,527,613,668,835,1199,1245,1322,350,440,435,1725,326,540,652,719,1357,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.725012,'amu*angstrom^2'), symmetry=1, barrier=(16.6695,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.725134,'amu*angstrom^2'), symmetry=1, barrier=(16.6723,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.72523,'amu*angstrom^2'), symmetry=1, barrier=(16.6745,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.724991,'amu*angstrom^2'), symmetry=1, barrier=(16.669,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (162.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.620161,0.100791,-0.000138656,8.93134e-08,-2.19769e-11,-162493,28.2954], Tmin=(100,'K'), Tmax=(1007.9,'K')), NASAPolynomial(coeffs=[21.9853,0.011077,-5.13994e-06,9.99497e-10,-7.1372e-14,-167049,-80.9406], Tmin=(1007.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1352.43,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(CsFFFO) + group(Cds-CdsCsCs) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs)"""),
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
    label = 'O=C(F)C(O)O[C](F)F(4276)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {7,S} {9,S}
5  O u0 p2 c0 {7,S} {11,S}
6  O u0 p2 c0 {8,D}
7  C u0 p0 c0 {4,S} {5,S} {8,S} {10,S}
8  C u0 p0 c0 {1,S} {6,D} {7,S}
9  C u1 p0 c0 {2,S} {3,S} {4,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-1006.03,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,486,617,768,1157,1926,493,600,700,1144,1293,233.27,233.27,1905.7,1905.7],'cm^-1')),
        HinderedRotor(inertia=(0.207919,'amu*angstrom^2'), symmetry=1, barrier=(8.02856,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.701275,'amu*angstrom^2'), symmetry=1, barrier=(27.0789,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.207917,'amu*angstrom^2'), symmetry=1, barrier=(8.02855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20792,'amu*angstrom^2'), symmetry=1, barrier=(8.02856,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (143.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.656279,0.0816698,-0.00013689,1.22816e-07,-4.34661e-11,-120885,31.8666], Tmin=(100,'K'), Tmax=(797.057,'K')), NASAPolynomial(coeffs=[8.92721,0.0286385,-1.5402e-05,3.06298e-09,-2.1563e-13,-121837,-3.86347], Tmin=(797.057,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1006.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsFFHO) + group(COCsFO) + radical(Csj(F1s)(F1s)(O2s-Cs))"""),
)

species(
    label = 'O=[C]C(O)OC(F)(F)F(4277)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u0 p2 c0 {7,S} {11,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
8  C u0 p0 c0 {1,S} {2,S} {3,S} {4,S}
9  C u1 p0 c0 {6,D} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-1042.53,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,431,527,613,668,835,1199,1245,1322,1855,455,950,259.208,1934.11,1935.34],'cm^-1')),
        HinderedRotor(inertia=(0.213757,'amu*angstrom^2'), symmetry=1, barrier=(10.2911,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.21561,'amu*angstrom^2'), symmetry=1, barrier=(10.3174,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.212886,'amu*angstrom^2'), symmetry=1, barrier=(10.2888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.557692,'amu*angstrom^2'), symmetry=1, barrier=(26.6588,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (143.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01387,0.0716055,-9.96298e-05,7.46685e-08,-2.27488e-11,-125285,29.8289], Tmin=(100,'K'), Tmax=(797.452,'K')), NASAPolynomial(coeffs=[10.1323,0.0258651,-1.35878e-05,2.73384e-09,-1.96155e-13,-126739,-12.098], Tmin=(797.452,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1042.53,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsFFFO) + group(Cds-OdCsH) + radical(CsCJ=O)"""),
)

species(
    label = 'CF3O(90)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 O u1 p2 c0 {5,S}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {4,S}
"""),
    E0 = (-649.398,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(84.9901,'amu')),
        NonlinearRotor(inertia=([83.2476,86.3025,90.6062],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([268.861,397.819,573.092,594.815,616.809,900.817,1178.75,1236.18,1267.81],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93566,0.0036847,5.89034e-05,-1.19115e-07,6.93035e-11,-78101,9.85551], Tmin=(10,'K'), Tmax=(596.518,'K')), NASAPolynomial(coeffs=[4.22715,0.0170524,-1.32403e-05,4.57332e-09,-5.80483e-13,-78408.4,6.31482], Tmin=(596.518,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-649.398,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""[O]C(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C(F)[CH]O(2997)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,S} {7,S}
3 O u0 p2 c0 {5,D}
4 C u1 p0 c0 {2,S} {5,S} {6,S}
5 C u0 p0 c0 {1,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {2,S}
"""),
    E0 = (-466.175,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,611,648,830,1210,1753,634.681],'cm^-1')),
        HinderedRotor(inertia=(0.178812,'amu*angstrom^2'), symmetry=1, barrier=(51.1156,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.178794,'amu*angstrom^2'), symmetry=1, barrier=(51.1155,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (77.0344,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.9486,0.00327696,8.62766e-05,-1.78856e-07,1.15095e-10,-56066.9,9.87278], Tmin=(10,'K'), Tmax=(486.336,'K')), NASAPolynomial(coeffs=[2.12805,0.0278262,-1.89743e-05,5.9068e-09,-6.93439e-13,-56003.1,16.1793], Tmin=(486.336,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-466.175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""ODC(F)[CH]O""", comment="""Thermo library: CHOF_G4"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(1006.05,'J/mol'), sigma=(4.32,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05587,-0.0054938,7.30062e-05,-1.34908e-07,7.857e-11,-58113,8.13986], Tmin=(10,'K'), Tmax=(570.231,'K')), NASAPolynomial(coeffs=[3.53924,0.0118582,-8.75006e-06,2.89362e-09,-3.54041e-13,-58277.3,8.38507], Tmin=(570.231,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-483.19,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""F[C](F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]C(O)C(=O)F(4278)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {5,S} {8,S}
3 O u1 p2 c0 {5,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
6 C u0 p0 c0 {1,S} {4,D} {5,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {2,S}
"""),
    E0 = (-548.071,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,486,617,768,1157,1926,180,3056.5],'cm^-1')),
        HinderedRotor(inertia=(0.0139113,'amu*angstrom^2'), symmetry=1, barrier=(6.66147,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.529262,'amu*angstrom^2'), symmetry=1, barrier=(12.1688,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.0338,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26506,0.0412616,-5.57075e-05,4.24648e-08,-1.329e-11,-65858.1,21.8266], Tmin=(100,'K'), Tmax=(775.28,'K')), NASAPolynomial(coeffs=[7.02534,0.0167011,-8.18786e-06,1.60206e-09,-1.13152e-13,-66596.2,0.0727119], Tmin=(775.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-548.071,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(COCsFO) + radical(C=OCOJ)"""),
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
    label = 'O=C(F)[CH]OC(F)(F)F(4279)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
8  C u1 p0 c0 {5,S} {9,S} {10,S}
9  C u0 p0 c0 {4,S} {6,D} {8,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-1056.08,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([431,527,613,668,835,1199,1245,1322,3025,407.5,1350,352.5,611,648,830,1210,1753,263.697,263.699,263.701,263.703],'cm^-1')),
        HinderedRotor(inertia=(0.00242422,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128582,'amu*angstrom^2'), symmetry=1, barrier=(6.34498,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.308459,'amu*angstrom^2'), symmetry=1, barrier=(15.2209,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (145.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.777496,0.0746179,-0.000103428,7.10688e-08,-1.91788e-11,-126904,26.5052], Tmin=(100,'K'), Tmax=(909.112,'K')), NASAPolynomial(coeffs=[13.7905,0.0173623,-8.95924e-06,1.7935e-09,-1.28646e-13,-129271,-35.0351], Tmin=(909.112,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1056.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)OsHH) + group(CsFFFO) + group(COCsFO) + radical(CCsJOCs)"""),
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
    label = '[O]C(OC(F)(F)F)C(=O)F(4280)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {8,S} {9,S}
6  O u1 p2 c0 {8,S}
7  O u0 p2 c0 {10,D}
8  C u0 p0 c0 {5,S} {6,S} {10,S} {11,S}
9  C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
10 C u0 p0 c0 {4,S} {7,D} {8,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-1207.46,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,431,527,613,668,835,1199,1245,1322,486,617,768,1157,1926,286.56,291.147,291.969,1787.06,1787.09],'cm^-1')),
        HinderedRotor(inertia=(0.125906,'amu*angstrom^2'), symmetry=1, barrier=(7.31143,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120206,'amu*angstrom^2'), symmetry=1, barrier=(7.32788,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125822,'amu*angstrom^2'), symmetry=1, barrier=(7.34491,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (161.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.603542,0.0814958,-0.000128823,1.09896e-07,-3.7678e-11,-145108,32.3132], Tmin=(100,'K'), Tmax=(752.916,'K')), NASAPolynomial(coeffs=[10.1371,0.0275162,-1.46455e-05,2.92185e-09,-2.07011e-13,-146449,-10.3479], Tmin=(752.916,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1207.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsFFFO) + group(COCsFO) + radical(C=OCOJ)"""),
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
    label = 'O[CH]OC(F)(F)F(1238)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {6,S}
4 O u0 p2 c0 {6,S} {7,S}
5 O u0 p2 c0 {7,S} {9,S}
6 C u0 p0 c0 {1,S} {2,S} {3,S} {4,S}
7 C u1 p0 c0 {4,S} {5,S} {8,S}
8 H u0 p0 c0 {7,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-897.696,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,431,527,613,668,835,1199,1245,1322,3025,407.5,1350,352.5,242.293,242.294,242.303],'cm^-1')),
        HinderedRotor(inertia=(0.00287308,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161932,'amu*angstrom^2'), symmetry=1, barrier=(6.74249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.594783,'amu*angstrom^2'), symmetry=1, barrier=(24.7601,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.73696,0.0184023,0.000148993,-4.39586e-07,3.50946e-10,-107965,13.0776], Tmin=(10,'K'), Tmax=(458.766,'K')), NASAPolynomial(coeffs=[7.32095,0.028945,-2.21207e-05,7.63896e-09,-9.71498e-13,-108734,-6.21355], Tmin=(458.766,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-897.696,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), label="""O[CH]OC(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C(F)[C](O)OC(F)(F)F(4281)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {8,S} {9,S}
6  O u0 p2 c0 {9,S} {11,S}
7  O u0 p2 c0 {10,D}
8  C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
9  C u1 p0 c0 {5,S} {6,S} {10,S}
10 C u0 p0 c0 {4,S} {7,D} {9,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-1245.94,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,431,527,613,668,835,1199,1245,1322,360,370,350,611,648,830,1210,1753,180,180,1051.6,1052.21],'cm^-1')),
        HinderedRotor(inertia=(0.246162,'amu*angstrom^2'), symmetry=1, barrier=(5.65976,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.247377,'amu*angstrom^2'), symmetry=1, barrier=(5.68768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.246741,'amu*angstrom^2'), symmetry=1, barrier=(5.67306,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.247762,'amu*angstrom^2'), symmetry=1, barrier=(5.69654,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (161.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.497783,0.0838824,-0.000130405,1.06383e-07,-3.47348e-11,-149733,32.6075], Tmin=(100,'K'), Tmax=(750.027,'K')), NASAPolynomial(coeffs=[11.4157,0.0256412,-1.38983e-05,2.79908e-09,-1.99705e-13,-151370,-16.9221], Tmin=(750.027,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1245.94,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsFFFO) + group(COCsFO) + radical(Cs_P)"""),
)

species(
    label = 'O=C=C(O)OC(F)(F)F(4230)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u0 p2 c0 {8,S} {10,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {4,S}
8  C u0 p0 c0 {4,S} {5,S} {9,D}
9  C u0 p0 c0 {6,D} {8,D}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-1044.96,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,431,527,613,668,835,1199,1245,1322,350,440,435,1725,2120,512.5,787.5,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.23945,'amu*angstrom^2'), symmetry=1, barrier=(28.4975,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23976,'amu*angstrom^2'), symmetry=1, barrier=(28.5046,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23876,'amu*angstrom^2'), symmetry=1, barrier=(28.4816,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.309427,0.0694013,-2.19745e-05,-6.18528e-08,4.03598e-11,-125500,25.0751], Tmin=(100,'K'), Tmax=(901.026,'K')), NASAPolynomial(coeffs=[35.8953,-0.0217549,1.39609e-05,-2.74747e-09,1.83546e-13,-134848,-161.489], Tmin=(901.026,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1044.96,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(CsFFFO) + group(Cds-(Cdd-O2d)OsOs) + missing(Cdd-CdO2d)"""),
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
    E0 = (-645.833,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-407.028,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-186.497,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-87.7412,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-113.458,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-345.692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-331.985,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-393.945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-332.387,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-402.834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-524.731,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-279.479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-315.976,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-461.915,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-377.604,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-374.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-341.994,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-434.397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-380.48,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-466.619,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=C(F)C(O)OC(F)(F)F(3727)'],
    products = ['HF(38)', 'CF2O(48)', 'O=CC(=O)F(557)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.58825e+11,'s^-1'), n=0.745914, Ea=(151.699,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R!H->C_5Br1sCl1sF1sH->F1s_Ext-2C-R_7R!H->C',), comment="""Estimated from node Root_N-1R!H->C_5Br1sCl1sF1sH->F1s_Ext-2C-R_7R!H->C
Multiplied by reaction path degeneracy 3.0"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(12)', 'OC(F)OC(F)(F)F(1198)'],
    products = ['O=C(F)C(O)OC(F)(F)F(3727)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.0026956,'m^3/(mol*s)'), n=2.93313, Ea=(378.297,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_2Br1sCl1sF1sH->F1s',), comment="""Estimated from node Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_2Br1sCl1sF1sH->F1s"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CF2(87)', 'O=C(F)C(O)OF(4126)'],
    products = ['O=C(F)C(O)OC(F)(F)F(3727)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.13916e+09,'m^3/(mol*s)'), n=-1.00842, Ea=(20.8544,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17406617270594374, var=16.167420070870673, Tref=1000.0, N=4, data_mean=0.0, correlation='OHOY',), comment="""Estimated from node OHOY"""),
)

reaction(
    label = 'reaction4',
    reactants = ['FOC(F)(F)F(1140)', 'O=C=CO(1378)'],
    products = ['O=C(F)C(O)OC(F)(F)F(3727)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.52781e-06,'m^3/(mol*s)'), n=3.38614, Ea=(181.703,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [doublebond;H_OR] for rate rule [Cdd_Cd_HNd;H_OR]
Euclidian distance = 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction5',
    reactants = ['OF(153)', 'O=C=COC(F)(F)F(4269)'],
    products = ['O=C(F)C(O)OC(F)(F)F(3727)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.37874e-06,'m^3/(mol*s)'), n=3.49511, Ea=(183.557,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [doublebond;H_OH] for rate rule [Cdd_Cd_HNd;H_OH]
Euclidian distance = 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction6',
    reactants = ['OC=C(F)OOC(F)(F)F(4270)'],
    products = ['O=C(F)C(O)OC(F)(F)F(3727)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(90.2106,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction7',
    reactants = ['OOC(F)=COC(F)(F)F(4271)'],
    products = ['O=C(F)C(O)OC(F)(F)F(3727)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(85.1575,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C(F)[C](O)OC(F)(F)F(4272)'],
    products = ['O=C(F)C(O)OC(F)(F)F(3727)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C(F)C([O])OC(F)(F)F(4273)'],
    products = ['O=C(F)C(O)OC(F)(F)F(3727)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C(OC(F)(F)F)[C](O)F(4274)'],
    products = ['O=C(F)C(O)OC(F)(F)F(3727)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['OC(F)=C(O)OC(F)(F)F(4275)'],
    products = ['O=C(F)C(O)OC(F)(F)F(3727)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(174.043,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_Ext-1R!H-R',), comment="""Estimated from node Root_3R!H->C_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', 'O=C(F)C(O)O[C](F)F(4276)'],
    products = ['O=C(F)C(O)OC(F)(F)F(3727)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', 'O=[C]C(O)OC(F)(F)F(4277)'],
    products = ['O=C(F)C(O)OC(F)(F)F(3727)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.52887e+10,'m^3/(mol*s)'), n=-0.421056, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.030895812897821735, var=3.1393582276389975, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['CF3O(90)', 'O=C(F)[CH]O(2997)'],
    products = ['O=C(F)C(O)OC(F)(F)F(3727)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CF3(44)', '[O]C(O)C(=O)F(4278)'],
    products = ['O=C(F)C(O)OC(F)(F)F(3727)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['OH(4)', 'O=C(F)[CH]OC(F)(F)F(4279)'],
    products = ['O=C(F)C(O)OC(F)(F)F(3727)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.7e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(5)', '[O]C(OC(F)(F)F)C(=O)F(4280)'],
    products = ['O=C(F)C(O)OC(F)(F)F(3727)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CFO(50)', 'O[CH]OC(F)(F)F(1238)'],
    products = ['O=C(F)C(O)OC(F)(F)F(3727)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -7.0 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(5)', 'O=C(F)[C](O)OC(F)(F)F(4281)'],
    products = ['O=C(F)C(O)OC(F)(F)F(3727)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.7979e+07,'m^3/(mol*s)'), n=0.240345, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_N-Sp-3C=2C_Ext-4R!H-R_Ext-2C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_N-Sp-3C=2C_Ext-4R!H-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction20',
    reactants = ['HF(38)', 'O=C=C(O)OC(F)(F)F(4230)'],
    products = ['O=C(F)C(O)OC(F)(F)F(3727)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(205.792,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

network(
    label = 'PDepNetwork #1352',
    isomers = [
        'O=C(F)C(O)OC(F)(F)F(3727)',
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
    label = 'PDepNetwork #1352',
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

