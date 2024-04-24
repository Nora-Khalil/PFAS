species(
    label = 'O=C(F)C(OF)OC(F)F(3723)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {6,S}
5  O u0 p2 c0 {8,S} {9,S}
6  O u0 p2 c0 {4,S} {8,S}
7  O u0 p2 c0 {10,D}
8  C u0 p0 c0 {5,S} {6,S} {10,S} {11,S}
9  C u0 p0 c0 {1,S} {2,S} {5,S} {12,S}
10 C u0 p0 c0 {3,S} {7,D} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-1083.33,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,519,593,830,1115,1166,1389,1413,3114,486,617,768,1157,1926,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (162.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.415584,0.0881,-0.000128215,9.41058e-08,-2.48627e-11,-130174,31.7826], Tmin=(100,'K'), Tmax=(612.319,'K')), NASAPolynomial(coeffs=[10.8474,0.0313391,-1.7058e-05,3.44904e-09,-2.46933e-13,-131665,-15.1709], Tmin=(612.319,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1083.33,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2sCF) + group(Cs-CsOsOsH) + group(CsFFHO) + group(COCsFO)"""),
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
    label = 'FOC(F)OC(F)F(1197)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {6,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {4,S} {7,S}
7  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
8  C u0 p0 c0 {2,S} {3,S} {5,S} {10,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-929.974,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,509,613,660,1171,1360,1414,3084,519,593,830,1115,1166,1389,1413,3114,195.503,195.523,195.541,195.617],'cm^-1')),
        HinderedRotor(inertia=(0.338861,'amu*angstrom^2'), symmetry=1, barrier=(9.20677,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.339302,'amu*angstrom^2'), symmetry=1, barrier=(9.20711,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12692,'amu*angstrom^2'), symmetry=1, barrier=(30.5824,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.03,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2940.86,'J/mol'), sigma=(4.81151,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=459.36 K, Pc=59.91 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.54699,0.0379871,5.58382e-05,-2.3765e-07,1.97646e-10,-111847,13.6596], Tmin=(10,'K'), Tmax=(477.99,'K')), NASAPolynomial(coeffs=[7.81606,0.0328879,-2.42688e-05,8.12383e-09,-1.00921e-12,-112605,-7.44468], Tmin=(477.99,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-929.974,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), label="""FOC(F)OC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'O=C(F)C(OF)OF(4127)',
    structure = adjacencyList("""1 F u0 p3 c0 {8,S}
2 F u0 p3 c0 {4,S}
3 F u0 p3 c0 {5,S}
4 O u0 p2 c0 {2,S} {7,S}
5 O u0 p2 c0 {3,S} {7,S}
6 O u0 p2 c0 {8,D}
7 C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
8 C u0 p0 c0 {1,S} {6,D} {7,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-522.788,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([504,610,1067,1155,1380,1390,370,380,2900,435,486,617,768,1157,1926,212.031,212.067,2173.9],'cm^-1')),
        HinderedRotor(inertia=(0.838196,'amu*angstrom^2'), symmetry=1, barrier=(26.7608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.409006,'amu*angstrom^2'), symmetry=1, barrier=(13.0522,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.838312,'amu*angstrom^2'), symmetry=1, barrier=(26.7618,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.023,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20189,0.0672313,-9.8739e-05,7.44878e-08,-2.24964e-11,-62781.3,24.5204], Tmin=(100,'K'), Tmax=(808.13,'K')), NASAPolynomial(coeffs=[10.725,0.0200941,-1.12452e-05,2.30918e-09,-1.67271e-13,-64320.4,-19.3945], Tmin=(808.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-522.788,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2sCF) + group(Cs-CsOsOsH) + group(COCsFO)"""),
)

species(
    label = 'FOC(F)F(253)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {4,S}
4 O u0 p2 c0 {3,S} {5,S}
5 C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6 H u0 p0 c0 {5,S}
"""),
    E0 = (-530.939,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,519,593,830,1115,1166,1389,1413,3114,180],'cm^-1')),
        HinderedRotor(inertia=(0.880111,'amu*angstrom^2'), symmetry=1, barrier=(20.2355,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0132,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88456,0.0076465,8.48623e-05,-2.22012e-07,1.63539e-10,-63851.5,9.5466], Tmin=(10,'K'), Tmax=(485.118,'K')), NASAPolynomial(coeffs=[5.31172,0.0175524,-1.27817e-05,4.268e-09,-5.32145e-13,-64245,1.06521], Tmin=(485.118,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-530.939,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""FOC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'O=C=COC(F)F(3356)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {5,S} {6,S}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
6 C u0 p0 c0 {3,S} {7,D} {9,S}
7 C u0 p0 c0 {4,D} {6,D}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-622.055,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([519,593,830,1115,1166,1389,1413,3114,3010,987.5,1337.5,450,1655,2120,512.5,787.5,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.06567,'amu*angstrom^2'), symmetry=1, barrier=(24.5018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06488,'amu*angstrom^2'), symmetry=1, barrier=(24.4836,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.704139,0.0594455,-6.71594e-05,3.47841e-08,-6.67293e-12,-74686.1,21.9085], Tmin=(100,'K'), Tmax=(1449.04,'K')), NASAPolynomial(coeffs=[18.7553,0.00171777,7.74896e-07,-2.32532e-10,1.74277e-14,-79088.3,-69.0115], Tmin=(1449.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-622.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + missing(O2d-Cdd) + group(CsFFHO) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d)"""),
)

species(
    label = 'FOC=C(F)OOC(F)F(4128)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {7,S}
5  O u0 p2 c0 {6,S} {8,S}
6  O u0 p2 c0 {5,S} {9,S}
7  O u0 p2 c0 {4,S} {10,S}
8  C u0 p0 c0 {1,S} {2,S} {5,S} {11,S}
9  C u0 p0 c0 {3,S} {6,S} {10,D}
10 C u0 p0 c0 {7,S} {9,D} {12,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-698.402,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,231,791,519,593,830,1115,1166,1389,1413,3114,326,540,652,719,1357,3010,987.5,1337.5,450,1655,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.216328,'amu*angstrom^2'), symmetry=1, barrier=(4.9738,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.215374,'amu*angstrom^2'), symmetry=1, barrier=(4.95186,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.1691,'amu*angstrom^2'), symmetry=1, barrier=(49.8719,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (162.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0205449,0.100154,-0.000175444,1.63648e-07,-5.95754e-11,-83865,32.8424], Tmin=(100,'K'), Tmax=(805.054,'K')), NASAPolynomial(coeffs=[8.48831,0.037515,-2.07947e-05,4.16552e-09,-2.93811e-13,-84575.2,-2.26473], Tmin=(805.054,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-698.402,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2sCF) + group(CsFFHO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH)"""),
)

species(
    label = 'FOOC(F)=COC(F)F(4129)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {7,S}
5  O u0 p2 c0 {8,S} {9,S}
6  O u0 p2 c0 {7,S} {10,S}
7  O u0 p2 c0 {4,S} {6,S}
8  C u0 p0 c0 {1,S} {2,S} {5,S} {11,S}
9  C u0 p0 c0 {5,S} {10,D} {12,S}
10 C u0 p0 c0 {3,S} {6,S} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-742.148,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,277,555,632,519,593,830,1115,1166,1389,1413,3114,3010,987.5,1337.5,450,1655,326,540,652,719,1357],'cm^-1')),
        HinderedRotor(inertia=(1.15801,'amu*angstrom^2'), symmetry=1, barrier=(26.6249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15789,'amu*angstrom^2'), symmetry=1, barrier=(26.6221,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15815,'amu*angstrom^2'), symmetry=1, barrier=(26.6282,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.15736,'amu*angstrom^2'), symmetry=1, barrier=(26.61,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (162.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.131559,0.0919816,-0.000128651,9.16722e-08,-2.60711e-11,-89126.5,30.6177], Tmin=(100,'K'), Tmax=(857.704,'K')), NASAPolynomial(coeffs=[14.1259,0.0267185,-1.45166e-05,2.96074e-09,-2.14224e-13,-91527.1,-34.7491], Tmin=(857.704,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-742.148,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2sFO) + group(CsFFHO) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs)"""),
)

species(
    label = '[O]C(F)[C](OF)OC(F)F(4130)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {6,S}
5  O u0 p2 c0 {9,S} {10,S}
6  O u0 p2 c0 {4,S} {10,S}
7  O u1 p2 c0 {8,S}
8  C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
9  C u0 p0 c0 {2,S} {3,S} {5,S} {12,S}
10 C u1 p0 c0 {5,S} {6,S} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-702.044,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,391,562,707,872,1109,1210,1289,3137,519,593,830,1115,1166,1389,1413,3114,360,370,350,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (162.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.235598,0.102412,-0.000174054,1.52353e-07,-5.22596e-11,-84292.8,34.7088], Tmin=(100,'K'), Tmax=(802.037,'K')), NASAPolynomial(coeffs=[12.1351,0.0292927,-1.594e-05,3.16832e-09,-2.22453e-13,-85909.8,-19.9535], Tmin=(802.037,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-702.044,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2sCF) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCFHO) + group(CsFFHO) + radical(O2sj(Cs-CsF1sH)) + radical(Cs_P)"""),
)

species(
    label = '[O]C(F)C(OF)O[C](F)F(4131)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {6,S}
5  O u0 p2 c0 {8,S} {10,S}
6  O u0 p2 c0 {4,S} {8,S}
7  O u1 p2 c0 {9,S}
8  C u0 p0 c0 {5,S} {6,S} {9,S} {11,S}
9  C u0 p0 c0 {1,S} {7,S} {8,S} {12,S}
10 C u1 p0 c0 {2,S} {3,S} {5,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-695.481,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,391,562,707,872,1109,1210,1289,3137,493,600,700,1144,1293,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (162.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.233063,0.103671,-0.000181123,1.62429e-07,-5.6595e-11,-83504.8,35.1706], Tmin=(100,'K'), Tmax=(816.848,'K')), NASAPolynomial(coeffs=[11.1857,0.0310161,-1.69675e-05,3.36701e-09,-2.35642e-13,-84811.8,-14.19], Tmin=(816.848,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-695.481,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2sCF) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCFHO) + group(CsFFHO) + radical(O2sj(Cs-CsF1sH)) + radical(Csj(F1s)(F1s)(O2s-Cs))"""),
)

species(
    label = 'O[C](F)C(OF)O[C](F)F(4132)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {6,S}
5  O u0 p2 c0 {8,S} {10,S}
6  O u0 p2 c0 {4,S} {8,S}
7  O u0 p2 c0 {9,S} {12,S}
8  C u0 p0 c0 {5,S} {6,S} {9,S} {11,S}
9  C u1 p0 c0 {1,S} {7,S} {8,S}
10 C u1 p0 c0 {2,S} {3,S} {5,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-727.501,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,1380,1390,370,380,2900,435,395,473,707,1436,493,600,700,1144,1293,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (162.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.585769,0.115787,-0.000219569,2.0435e-07,-7.19695e-11,-87347.3,36.3359], Tmin=(100,'K'), Tmax=(847.388,'K')), NASAPolynomial(coeffs=[10.6249,0.0324579,-1.82323e-05,3.6004e-09,-2.49095e-13,-88155.4,-9.45028], Tmin=(847.388,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-727.501,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2sCF) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCFHO) + group(CsFFHO) + radical(CsCsF1sO2s) + radical(Csj(F1s)(F1s)(O2s-Cs))"""),
)

species(
    label = 'OC(F)=C(OF)OC(F)F(4133)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {6,S}
5  O u0 p2 c0 {8,S} {9,S}
6  O u0 p2 c0 {4,S} {9,S}
7  O u0 p2 c0 {10,S} {12,S}
8  C u0 p0 c0 {1,S} {2,S} {5,S} {11,S}
9  C u0 p0 c0 {5,S} {6,S} {10,D}
10 C u0 p0 c0 {3,S} {7,S} {9,D}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-961.274,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,3615,1277.5,1000,519,593,830,1115,1166,1389,1413,3114,350,440,435,1725,326,540,652,719,1357,269.317,269.634,269.746,269.976],'cm^-1')),
        HinderedRotor(inertia=(0.195967,'amu*angstrom^2'), symmetry=1, barrier=(10.122,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.195572,'amu*angstrom^2'), symmetry=1, barrier=(10.116,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.634785,'amu*angstrom^2'), symmetry=1, barrier=(32.78,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.195406,'amu*angstrom^2'), symmetry=1, barrier=(10.1182,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (162.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.577129,0.108845,-0.000182537,1.5332e-07,-5.03276e-11,-115458,29.9696], Tmin=(100,'K'), Tmax=(798.032,'K')), NASAPolynomial(coeffs=[14.9853,0.0248582,-1.34285e-05,2.65461e-09,-1.85655e-13,-117751,-40.4055], Tmin=(798.032,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-961.274,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2sCF) + group(O2s-(Cds-Cd)H) + group(CsFFHO) + group(Cds-CdsCsCs) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs)"""),
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
    label = 'O=C(F)C(OF)O[CH]F(4134)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {5,S}
4  O u0 p2 c0 {7,S} {9,S}
5  O u0 p2 c0 {3,S} {7,S}
6  O u0 p2 c0 {8,D}
7  C u0 p0 c0 {4,S} {5,S} {8,S} {10,S}
8  C u0 p0 c0 {1,S} {6,D} {7,S}
9  C u1 p0 c0 {2,S} {4,S} {11,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-652.202,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,486,617,768,1157,1926,580,1155,1237,1373,3147,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (143.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.455658,0.0861628,-0.000144185,1.28294e-07,-4.52254e-11,-78322,30.9006], Tmin=(100,'K'), Tmax=(783.157,'K')), NASAPolynomial(coeffs=[9.77128,0.0288623,-1.58178e-05,3.17098e-09,-2.24405e-13,-79483,-9.86171], Tmin=(783.157,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-652.202,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2sCF) + group(Cs-CsOsOsH) + group(CsFHHO) + group(COCsFO) + radical(Csj(F1s)(O2s-Cs)(H))"""),
)

species(
    label = 'O=[C]C(OF)OC(F)F(4135)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {5,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u0 p2 c0 {3,S} {7,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
8  C u0 p0 c0 {1,S} {2,S} {4,S} {11,S}
9  C u1 p0 c0 {6,D} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-674.666,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,519,593,830,1115,1166,1389,1413,3114,1855,455,950,272.815,272.816,272.817,2444.34],'cm^-1')),
        HinderedRotor(inertia=(0.141682,'amu*angstrom^2'), symmetry=1, barrier=(7.48318,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.539875,'amu*angstrom^2'), symmetry=1, barrier=(28.5142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141684,'amu*angstrom^2'), symmetry=1, barrier=(7.48319,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.539873,'amu*angstrom^2'), symmetry=1, barrier=(28.5142,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (143.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.619107,0.0814671,-0.000129672,1.11992e-07,-3.90167e-11,-81028.6,30.7683], Tmin=(100,'K'), Tmax=(741.616,'K')), NASAPolynomial(coeffs=[9.82381,0.0284081,-1.54532e-05,3.11163e-09,-2.21688e-13,-82300,-10.2549], Tmin=(741.616,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-674.666,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2sCF) + group(Cs-CsOsOsH) + group(CsFFHO) + group(Cds-OdCsH) + radical(CsCJ=O)"""),
)

species(
    label = '[O]C(OC(F)F)C(=O)F(4136)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
8  C u0 p0 c0 {1,S} {2,S} {4,S} {11,S}
9  C u0 p0 c0 {3,S} {6,D} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-974.104,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,519,593,830,1115,1166,1389,1413,3114,486,617,768,1157,1926,202.924,203.02,203.517,203.538,1919.7],'cm^-1')),
        HinderedRotor(inertia=(0.00407577,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160327,'amu*angstrom^2'), symmetry=1, barrier=(4.70065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160041,'amu*angstrom^2'), symmetry=1, barrier=(4.69354,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (143.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.82418,0.0771509,-0.000124591,1.10502e-07,-3.90682e-11,-117050,30.8851], Tmin=(100,'K'), Tmax=(787.839,'K')), NASAPolynomial(coeffs=[8.47479,0.0290178,-1.52609e-05,3.02133e-09,-2.12679e-13,-117968,-2.37062], Tmin=(787.839,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-974.104,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsFFHO) + group(COCsFO) + radical(C=OCOJ)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3063.56,'J/mol'), sigma=(4.99758,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=478.52 K, Pc=55.69 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.07235,-0.00674902,8.67846e-05,-1.53143e-07,8.6196e-11,-51263.6,9.19392], Tmin=(10,'K'), Tmax=(577.889,'K')), NASAPolynomial(coeffs=[2.80054,0.016572,-1.1432e-05,3.6345e-09,-4.33761e-13,-51359,12.5348], Tmin=(577.889,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-426.241,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""[O]C(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]C(OF)C(=O)F(4137)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {3,S}
3 O u0 p2 c0 {2,S} {6,S}
4 O u1 p2 c0 {6,S}
5 O u0 p2 c0 {7,D}
6 C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
7 C u0 p0 c0 {1,S} {5,D} {6,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-413.563,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,486,617,768,1157,1926,335.977,337.932,1490.25],'cm^-1')),
        HinderedRotor(inertia=(0.0952546,'amu*angstrom^2'), symmetry=1, barrier=(7.67996,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0945978,'amu*angstrom^2'), symmetry=1, barrier=(7.67891,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.024,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88882,0.0520439,-7.43838e-05,5.15874e-08,-1.1605e-11,-49669.4,23.3701], Tmin=(100,'K'), Tmax=(601.933,'K')), NASAPolynomial(coeffs=[8.20196,0.0180785,-9.64678e-06,1.93247e-09,-1.37512e-13,-50574.2,-5.08449], Tmin=(601.933,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-413.563,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(COCsFO) + radical(C=OCOJ)"""),
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
    label = 'O=C(F)[CH]OC(F)F(4138)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {7,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {2,S} {4,S} {9,S}
7  C u1 p0 c0 {4,S} {8,S} {10,S}
8  C u0 p0 c0 {3,S} {5,D} {7,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-822.73,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([519,593,830,1115,1166,1389,1413,3114,3025,407.5,1350,352.5,611,648,830,1210,1753,307.668,307.668,307.669,307.669],'cm^-1')),
        HinderedRotor(inertia=(0.114928,'amu*angstrom^2'), symmetry=1, barrier=(7.72009,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114928,'amu*angstrom^2'), symmetry=1, barrier=(7.72009,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.284281,'amu*angstrom^2'), symmetry=1, barrier=(19.0962,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0065,0.0701373,-9.85558e-05,7.0649e-08,-2.00716e-11,-98847.6,25.0495], Tmin=(100,'K'), Tmax=(861.897,'K')), NASAPolynomial(coeffs=[12.0049,0.0190916,-9.71393e-06,1.92732e-09,-1.37255e-13,-100743,-26.3765], Tmin=(861.897,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-822.73,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)OsHH) + group(CsFFHO) + group(COCsFO) + radical(CCsJOCs)"""),
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
    label = 'FO[CH]OC(F)F(1209)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {5,S}
4 O u0 p2 c0 {6,S} {7,S}
5 O u0 p2 c0 {3,S} {7,S}
6 C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
7 C u1 p0 c0 {4,S} {5,S} {9,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-529.581,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,519,593,830,1115,1166,1389,1413,3114,3025,407.5,1350,352.5,180,180,180,1634.79],'cm^-1')),
        HinderedRotor(inertia=(0.196505,'amu*angstrom^2'), symmetry=1, barrier=(4.51805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.195814,'amu*angstrom^2'), symmetry=1, barrier=(4.50216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197526,'amu*angstrom^2'), symmetry=1, barrier=(4.5415,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33122,0.077336,-0.00017024,1.81274e-07,-6.95234e-11,-63615.7,22.1955], Tmin=(100,'K'), Tmax=(861.764,'K')), NASAPolynomial(coeffs=[0.18225,0.0370208,-2.06104e-05,4.05141e-09,-2.78857e-13,-61722.7,37.4022], Tmin=(861.764,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-529.581,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(OCJO)"""),
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
    label = 'O=C(F)[C](OF)OC(F)F(4139)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {6,S}
5  O u0 p2 c0 {8,S} {9,S}
6  O u0 p2 c0 {4,S} {9,S}
7  O u0 p2 c0 {10,D}
8  C u0 p0 c0 {1,S} {2,S} {5,S} {11,S}
9  C u1 p0 c0 {5,S} {6,S} {10,S}
10 C u0 p0 c0 {3,S} {7,D} {9,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-878.083,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,519,593,830,1115,1166,1389,1413,3114,360,370,350,611,648,830,1210,1753,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (161.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.161814,0.0929669,-0.000157266,1.38727e-07,-4.83867e-11,-105479,33.3409], Tmin=(100,'K'), Tmax=(777.735,'K')), NASAPolynomial(coeffs=[11.0618,0.0282657,-1.58123e-05,3.18866e-09,-2.26233e-13,-106913,-14.8251], Tmin=(777.735,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-878.083,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2sCF) + group(Cs-CsOsOsH) + group(CsFFHO) + group(COCsFO) + radical(Cs_P)"""),
)

species(
    label = 'O=C(F)C(OF)O[C](F)F(4140)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {6,S}
5  O u0 p2 c0 {8,S} {10,S}
6  O u0 p2 c0 {4,S} {8,S}
7  O u0 p2 c0 {9,D}
8  C u0 p0 c0 {5,S} {6,S} {9,S} {11,S}
9  C u0 p0 c0 {1,S} {7,D} {8,S}
10 C u1 p0 c0 {2,S} {3,S} {5,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-871.52,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,486,617,768,1157,1926,493,600,700,1144,1293,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (161.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.157379,0.0943141,-0.00016467,1.49261e-07,-5.29192e-11,-104691,33.8274], Tmin=(100,'K'), Tmax=(799.427,'K')), NASAPolynomial(coeffs=[10.1395,0.02994,-1.68102e-05,3.38014e-09,-2.3881e-13,-105826,-9.21229], Tmin=(799.427,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-871.52,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2sCF) + group(Cs-CsOsOsH) + group(CsFFHO) + group(COCsFO) + radical(Csj(F1s)(F1s)(O2s-Cs))"""),
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
    label = 'O=C(F)C(=O)OC(F)F(3055)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u0 p2 c0 {8,D}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {2,S} {4,S} {10,S}
8  C u0 p0 c0 {4,S} {5,D} {9,S}
9  C u0 p0 c0 {3,S} {6,D} {8,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-1174.43,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([519,593,830,1115,1166,1389,1413,3114,286,619,818,1246,1924,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37312,0.064167,-8.43745e-05,6.1786e-08,-1.88993e-11,-141163,25.1585], Tmin=(100,'K'), Tmax=(785.358,'K')), NASAPolynomial(coeffs=[8.4759,0.0279904,-1.52773e-05,3.13039e-09,-2.27321e-13,-142278,-7.3921], Tmin=(785.358,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1174.43,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(CsFFHO) + group(Cds-O2d(Cds-O2d)O2s) + group(COCFO)"""),
)

species(
    label = 'O=C=C(OF)OC(F)F(3928)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {5,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u0 p2 c0 {3,S} {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {2,S} {4,S} {10,S}
8  C u0 p0 c0 {4,S} {5,S} {9,D}
9  C u0 p0 c0 {6,D} {8,D}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-653.798,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,519,593,830,1115,1166,1389,1413,3114,350,440,435,1725,2120,512.5,787.5,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.28844,'amu*angstrom^2'), symmetry=1, barrier=(29.6238,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.288,'amu*angstrom^2'), symmetry=1, barrier=(29.6136,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.28647,'amu*angstrom^2'), symmetry=1, barrier=(29.5786,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.35572,0.0900945,-0.000109108,5.68516e-08,-1.06522e-11,-78417,30.6716], Tmin=(100,'K'), Tmax=(1552.58,'K')), NASAPolynomial(coeffs=[29.0521,-0.00885957,6.40951e-06,-1.30214e-09,8.87364e-14,-85374.7,-121.405], Tmin=(1552.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-653.798,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2sCF) + missing(O2d-Cdd) + group(CsFFHO) + group(Cds-(Cdd-O2d)OsOs) + missing(Cdd-CdO2d)"""),
)

species(
    label = 'F2(106)',
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
    E0 = (-585.474,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-215.321,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-319.182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (92.4202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (35.2762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (43.356,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-148.558,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-180.139,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-213.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-204.679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-236.698,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-331.136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-113.481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-135.945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-435.384,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-284.238,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-204.444,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-254.215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-254.111,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-200.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-193.886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-389.933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-502.094,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-272.926,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-358.099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=C(F)C(OF)OC(F)F(3723)'],
    products = ['HF(38)', 'CF2O(48)', 'O=CC(=O)F(557)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(0.00285451,'s^-1'), n=4.62568, Ea=(32.0264,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(12)', 'FOC(F)OC(F)F(1197)'],
    products = ['O=C(F)C(OF)OC(F)F(3723)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.0026956,'m^3/(mol*s)'), n=2.93313, Ea=(367.565,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_2Br1sCl1sF1sH->F1s',), comment="""Estimated from node Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_2Br1sCl1sF1sH->F1s"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CF2(87)', 'O=C(F)C(O)OF(4126)'],
    products = ['O=C(F)C(OF)OC(F)F(3723)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.76395e-12,'m^3/(mol*s)'), n=5.02686, Ea=(75.9973,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='OH_N-2Br1sCl1sF1sHI1s->H_2Br1sCl1sF1sI1s->F1s',), comment="""Estimated from node OH_N-2Br1sCl1sF1sHI1s->H_2Br1sCl1sF1sI1s->F1s"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CHF(86)', 'O=C(F)C(OF)OF(4127)'],
    products = ['O=C(F)C(OF)OC(F)F(3723)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4.27832e+09,'m^3/(mol*s)'), n=-1.00842, Ea=(10.6239,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17406617270594374, var=16.167420070870673, Tref=1000.0, N=4, data_mean=0.0, correlation='OHOY',), comment="""Estimated from node OHOY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction5',
    reactants = ['FOC(F)F(253)', 'O=C=COF(1390)'],
    products = ['O=C(F)C(OF)OC(F)F(3723)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.52781e-06,'m^3/(mol*s)'), n=3.38614, Ea=(181.703,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [doublebond;H_OR] for rate rule [Cdd_Cd_HNd;H_OR]
Euclidian distance = 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction6',
    reactants = ['FOF(716)', 'O=C=COC(F)F(3356)'],
    products = ['O=C(F)C(OF)OC(F)F(3723)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.75748e-06,'m^3/(mol*s)'), n=3.49511, Ea=(183.557,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [doublebond;H_OH] for rate rule [Cdd_Cd_HNd;H_OH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction7',
    reactants = ['FOC=C(F)OOC(F)F(4128)'],
    products = ['O=C(F)C(OF)OC(F)F(3723)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(84.0147,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction8',
    reactants = ['FOOC(F)=COC(F)F(4129)'],
    products = ['O=C(F)C(OF)OC(F)F(3723)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(96.1799,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C(F)[C](OF)OC(F)F(4130)'],
    products = ['O=C(F)C(OF)OC(F)F(3723)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C(F)C(OF)O[C](F)F(4131)'],
    products = ['O=C(F)C(OF)OC(F)F(3723)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['O[C](F)C(OF)O[C](F)F(4132)'],
    products = ['O=C(F)C(OF)OC(F)F(3723)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['OC(F)=C(OF)OC(F)F(4133)'],
    products = ['O=C(F)C(OF)OC(F)F(3723)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(164.309,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_Ext-1R!H-R',), comment="""Estimated from node Root_3R!H->C_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', 'O=C(F)C(OF)O[CH]F(4134)'],
    products = ['O=C(F)C(OF)OC(F)F(3723)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F(37)', 'O=[C]C(OF)OC(F)F(4135)'],
    products = ['O=C(F)C(OF)OC(F)F(3723)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.52887e+10,'m^3/(mol*s)'), n=-0.421056, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.030895812897821735, var=3.1393582276389975, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F(37)', '[O]C(OC(F)F)C(=O)F(4136)'],
    products = ['O=C(F)C(OF)OC(F)F(3723)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.55564e+07,'m^3/(mol*s)'), n=-5.63145e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.7121440562946592, var=2.9508506800589083, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]C(F)F(170)', 'O=C(F)[CH]OF(2999)'],
    products = ['O=C(F)C(OF)OC(F)F(3723)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction17',
    reactants = ['CHF2(81)', '[O]C(OF)C(=O)F(4137)'],
    products = ['O=C(F)C(OF)OC(F)F(3723)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]F(154)', 'O=C(F)[CH]OC(F)F(4138)'],
    products = ['O=C(F)C(OF)OC(F)F(3723)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CFO(50)', 'FO[CH]OC(F)F(1209)'],
    products = ['O=C(F)C(OF)OC(F)F(3723)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -7.2 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(5)', 'O=C(F)[C](OF)OC(F)F(4139)'],
    products = ['O=C(F)C(OF)OC(F)F(3723)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.7979e+07,'m^3/(mol*s)'), n=0.240345, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_N-Sp-3C=2C_Ext-4R!H-R_Ext-2C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_N-Sp-3C=2C_Ext-4R!H-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(5)', 'O=C(F)C(OF)O[C](F)F(4140)'],
    products = ['O=C(F)C(OF)OC(F)F(3723)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.1766,'m^3/(mol*s)'), n=1.94174, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_N-4R!H->C_Ext-2C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_N-4R!H->C_Ext-2C-R
Ea raised from -2.1 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['O=C(F)C(OF)OC(F)F(3723)'],
    products = ['CF2O(48)', 'OC(F)=COF(1234)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.75934e+14,'s^-1'), n=-0.752179, Ea=(227.567,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.09940871819346306, var=3.010161187748991, Tref=1000.0, N=31, data_mean=0.0, correlation='Root_1R!H->C',), comment="""Estimated from node Root_1R!H->C"""),
)

reaction(
    label = 'reaction23',
    reactants = ['HF(38)', 'O=C(F)C(=O)OC(F)F(3055)'],
    products = ['O=C(F)C(OF)OC(F)F(3723)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(487.621,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction24',
    reactants = ['HF(38)', 'O=C=C(OF)OC(F)F(3928)'],
    products = ['O=C(F)C(OF)OC(F)F(3723)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(196.156,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction25',
    reactants = ['O=C(F)C(OF)OC(F)F(3723)'],
    products = ['F2(106)', 'CHFO(46)', 'O=CC(=O)F(557)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.00570902,'s^-1'), n=4.62568, Ea=(259.402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

network(
    label = 'PDepNetwork #1348',
    isomers = [
        'O=C(F)C(OF)OC(F)F(3723)',
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
    label = 'PDepNetwork #1348',
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

