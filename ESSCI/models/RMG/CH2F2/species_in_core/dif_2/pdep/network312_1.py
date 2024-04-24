species(
    label = 'OC(F)(F)C(F)(F)OF(625)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {7,S}
6  O u0 p2 c0 {8,S} {10,S}
7  O u0 p2 c0 {5,S} {9,S}
8  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
9  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-1168.68,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,183,263,327,399,503,589,506,644,608,780,1145,1213,1339,1481,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.367915,'amu*angstrom^2'), symmetry=1, barrier=(8.45909,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07194,'amu*angstrom^2'), symmetry=1, barrier=(24.6461,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.369212,'amu*angstrom^2'), symmetry=1, barrier=(8.48892,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.02,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.44447,0.0457561,8.58585e-05,-4.03399e-07,3.7423e-10,-140559,13.6644], Tmin=(10,'K'), Tmax=(436.067,'K')), NASAPolynomial(coeffs=[9.9846,0.0318864,-2.50848e-05,8.7639e-09,-1.12164e-12,-141568,-17.488], Tmin=(436.067,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1168.68,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), label="""OC(F)(F)C(F)(F)OF""", comment="""Thermo library: CHOF_G4"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28591,0.0107608,-1.05382e-05,4.89881e-09,-8.86384e-13,-24340.7,13.1348], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[5.33121,0.00197748,-9.60248e-07,2.10704e-10,-1.5954e-14,-25190.9,-2.56367], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-203.712,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CF2""", comment="""Thermo library: halogens"""),
)

species(
    label = 'OC(F)(F)OF(707)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {5,S}
4 O u0 p2 c0 {6,S} {7,S}
5 O u0 p2 c0 {3,S} {6,S}
6 C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-744.768,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,435,565,619,662,854,1178,1396,180],'cm^-1')),
        HinderedRotor(inertia=(0.698652,'amu*angstrom^2'), symmetry=1, barrier=(16.0634,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.700175,'amu*angstrom^2'), symmetry=1, barrier=(16.0984,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.013,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.79123,0.0143291,0.000118735,-3.5357e-07,2.81483e-10,-89567.8,10.2299], Tmin=(10,'K'), Tmax=(466.2,'K')), NASAPolynomial(coeffs=[7.41097,0.0196271,-1.52851e-05,5.35113e-09,-6.88488e-13,-90300.3,-8.70804], Tmin=(466.2,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-744.768,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), label="""OC(F)(F)OF""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'F[C]C(F)(F)OF(714)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {7,S}
4 F u0 p3 c0 {5,S}
5 O u0 p2 c0 {4,S} {6,S}
6 C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
7 C u0 p1 c0 {3,S} {6,S}
"""),
    E0 = (-404.826,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2850,1437.5,1250,1305,750,350,617,898,1187,180],'cm^-1')),
        HinderedRotor(inertia=(1.52162,'amu*angstrom^2'), symmetry=1, barrier=(34.985,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.52282,'amu*angstrom^2'), symmetry=1, barrier=(35.0125,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.014,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52555,0.0589818,-9.49275e-05,7.55253e-08,-2.36061e-11,-48604.5,19.5038], Tmin=(100,'K'), Tmax=(787.137,'K')), NASAPolynomial(coeffs=[10.5434,0.0131539,-7.59295e-06,1.55453e-09,-1.11631e-13,-50024.1,-21.8434], Tmin=(787.137,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-404.826,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCFFO) + group(CJ2_singlet-FCs)"""),
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
    label = 'OC(F)(F)[C]F(715)',
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
    label = 'OC(F)(OF)C(F)(F)F(717)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {7,S}
6  O u0 p2 c0 {8,S} {10,S}
7  O u0 p2 c0 {5,S} {8,S}
8  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
9  C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-1170.8,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,375,361,584,565,722,1474,193,295,551,588,656,1146,1192,1350,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.534747,'amu*angstrom^2'), symmetry=1, barrier=(12.2949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.53381,'amu*angstrom^2'), symmetry=1, barrier=(12.2733,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.732223,'amu*angstrom^2'), symmetry=1, barrier=(16.8353,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (152.02,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.5443,0.0339888,0.000196267,-7.10619e-07,6.46552e-10,-140811,12.8483], Tmin=(10,'K'), Tmax=(420.079,'K')), NASAPolynomial(coeffs=[10.7309,0.0306776,-2.44361e-05,8.65695e-09,-1.12194e-12,-141990,-22.4287], Tmin=(420.079,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1170.8,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), label="""OC(F)(OF)C(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FOC(F)=C(F)F(718)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {7,S}
4 F u0 p3 c0 {5,S}
5 O u0 p2 c0 {4,S} {6,S}
6 C u0 p0 c0 {1,S} {5,S} {7,D}
7 C u0 p0 c0 {2,S} {3,S} {6,D}
"""),
    E0 = (-566.291,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,326,540,652,719,1357,182,240,577,636,1210,1413,640.155],'cm^-1')),
        HinderedRotor(inertia=(0.101302,'amu*angstrom^2'), symmetry=1, barrier=(2.32914,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.014,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.70077,0.0244422,6.10847e-05,-2.43384e-07,2.18666e-10,-68108.8,12.2443], Tmin=(10,'K'), Tmax=(428.07,'K')), NASAPolynomial(coeffs=[6.23355,0.0239046,-1.80789e-05,6.12494e-09,-7.67187e-13,-68537.5,-0.301166], Tmin=(428.07,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-566.291,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""FOC(F)DC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OC(F)=C(F)F(719)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {6,S}
4 O u0 p2 c0 {5,S} {7,S}
5 C u0 p0 c0 {1,S} {4,S} {6,D}
6 C u0 p0 c0 {2,S} {3,S} {5,D}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-684.209,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,326,540,652,719,1357,182,240,577,636,1210,1413],'cm^-1')),
        HinderedRotor(inertia=(0.0938502,'amu*angstrom^2'), symmetry=1, barrier=(2.1578,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0239,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.81959,0.0136405,9.8061e-05,-3.14486e-07,2.74877e-10,-82290.6,11.3581], Tmin=(10,'K'), Tmax=(414.43,'K')), NASAPolynomial(coeffs=[5.45368,0.0214222,-1.53545e-05,5.09528e-09,-6.32816e-13,-82628.3,2.47368], Tmin=(414.43,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-684.209,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""OC(F)DC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'O[C](F)C(F)(F)OF(720)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {8,S}
4 F u0 p3 c0 {5,S}
5 O u0 p2 c0 {4,S} {7,S}
6 O u0 p2 c0 {8,S} {9,S}
7 C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8 C u1 p0 c0 {3,S} {6,S} {7,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-745.312,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,253,525,597,667,842,1178,1324,395,473,707,1436,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.489157,'amu*angstrom^2'), symmetry=1, barrier=(11.2467,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.488514,'amu*angstrom^2'), symmetry=1, barrier=(11.2319,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.48755,'amu*angstrom^2'), symmetry=1, barrier=(11.2097,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (133.022,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.120383,0.110758,-0.000240513,2.31767e-07,-8.07934e-11,-89510.9,24.1496], Tmin=(100,'K'), Tmax=(898.753,'K')), NASAPolynomial(coeffs=[8.53213,0.0233375,-1.29778e-05,2.43593e-09,-1.58548e-13,-89090.8,-5.68039], Tmin=(898.753,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-745.312,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sF1sO2s)(F1s)(O2s-H))"""),
)

species(
    label = 'OC(F)(F)[C](F)OF(721)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {8,S}
4 F u0 p3 c0 {6,S}
5 O u0 p2 c0 {7,S} {9,S}
6 O u0 p2 c0 {4,S} {8,S}
7 C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8 C u1 p0 c0 {3,S} {6,S} {7,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-758.32,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,253,525,597,667,842,1178,1324,395,473,707,1436,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.259049,'amu*angstrom^2'), symmetry=1, barrier=(5.95604,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.912516,'amu*angstrom^2'), symmetry=1, barrier=(20.9805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.259084,'amu*angstrom^2'), symmetry=1, barrier=(5.95686,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (133.022,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.954528,0.0726762,-0.000114932,8.92724e-08,-2.63918e-11,-91100.4,25.0667], Tmin=(100,'K'), Tmax=(691.276,'K')), NASAPolynomial(coeffs=[11.8807,0.0165805,-8.676e-06,1.71411e-09,-1.20476e-13,-92781.3,-24.8434], Tmin=(691.276,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-758.32,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]C(F)(F)C(O)(F)F(722)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {8,S}
4 F u0 p3 c0 {8,S}
5 O u0 p2 c0 {7,S} {9,S}
6 O u1 p2 c0 {8,S}
7 C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8 C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-1072.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,223,363,546,575,694,1179,1410,351,323,533,609,664,892,1120,1201,180],'cm^-1')),
        HinderedRotor(inertia=(0.171742,'amu*angstrom^2'), symmetry=1, barrier=(3.94869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171702,'amu*angstrom^2'), symmetry=1, barrier=(3.94777,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (133.022,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24997,0.0653559,-9.69483e-05,7.3886e-08,-2.2371e-11,-128882,21.2313], Tmin=(100,'K'), Tmax=(809.036,'K')), NASAPolynomial(coeffs=[10.7752,0.0182623,-9.63497e-06,1.93837e-09,-1.38754e-13,-130423,-22.704], Tmin=(809.036,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1072.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sF1s))"""),
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
    label = 'FOC(F)(F)[C](F)F(723)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {8,S}
4 F u0 p3 c0 {8,S}
5 F u0 p3 c0 {6,S}
6 O u0 p2 c0 {5,S} {7,S}
7 C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
8 C u1 p0 c0 {3,S} {4,S} {7,S}
"""),
    E0 = (-744.472,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,253,525,597,667,842,1178,1324,190,488,555,1236,1407,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.822843,'amu*angstrom^2'), symmetry=1, barrier=(18.9188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.822102,'amu*angstrom^2'), symmetry=1, barrier=(18.9017,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (135.013,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.62228,0.0269041,0.0001746,-6.00025e-07,5.19722e-10,-89530.3,13.0155], Tmin=(10,'K'), Tmax=(443.604,'K')), NASAPolynomial(coeffs=[11.1479,0.0214546,-1.80038e-05,6.57499e-09,-8.70046e-13,-90812.1,-24.0954], Tmin=(443.604,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-744.472,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""FOC(F)(F)[C](F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]C(F)(F)C(F)(F)OF(724)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {8,S}
2 F u0 p3 c0 {8,S}
3 F u0 p3 c0 {9,S}
4 F u0 p3 c0 {9,S}
5 F u0 p3 c0 {6,S}
6 O u0 p2 c0 {5,S} {8,S}
7 O u1 p2 c0 {9,S}
8 C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
9 C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
"""),
    E0 = (-898.289,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,223,363,546,575,694,1179,1410,351,323,533,609,664,892,1120,1201,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.824484,'amu*angstrom^2'), symmetry=1, barrier=(18.9565,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.824523,'amu*angstrom^2'), symmetry=1, barrier=(18.9574,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (151.012,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.5958,0.0288532,0.000181869,-5.96026e-07,4.99985e-10,-108035,14.3218], Tmin=(10,'K'), Tmax=(450.82,'K')), NASAPolynomial(coeffs=[10.2473,0.0305115,-2.55295e-05,9.21288e-09,-1.20115e-12,-109251,-19.3069], Tmin=(450.82,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-898.289,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""[O]C(F)(F)C(F)(F)OF""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'OC(F)(F)[C](F)F(725)',
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
    label = 'FO[C](F)F(192)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {4,S}
4 O u0 p2 c0 {3,S} {5,S}
5 C u1 p0 c0 {1,S} {2,S} {4,S}
"""),
    E0 = (-317.517,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,493,600,700,1144,1293,180],'cm^-1')),
        HinderedRotor(inertia=(0.523179,'amu*angstrom^2'), symmetry=1, barrier=(12.0289,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66912,0.0316874,-4.80467e-05,3.64259e-08,-1.08808e-11,-38142.7,14.3743], Tmin=(100,'K'), Tmax=(821.793,'K')), NASAPolynomial(coeffs=[7.60337,0.00767068,-4.20998e-06,8.64397e-10,-6.26024e-14,-38953.7,-8.46216], Tmin=(821.793,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-317.517,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsF1sF1sO2s)"""),
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
    label = 'O=C(F)C(F)(F)OF(726)',
    structure = adjacencyList("""1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {8,S}
4 F u0 p3 c0 {5,S}
5 O u0 p2 c0 {4,S} {7,S}
6 O u0 p2 c0 {8,D}
7 C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8 C u0 p0 c0 {3,S} {6,D} {7,S}
"""),
    E0 = (-861.217,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,251,367,519,700,855,1175,1303,486,617,768,1157,1926,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.253857,'amu*angstrom^2'), symmetry=1, barrier=(5.83667,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.665343,'amu*angstrom^2'), symmetry=1, barrier=(15.2976,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.014,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.50848,0.0492077,-6.50416e-05,4.38118e-08,-1.19234e-11,-103579,13.1827], Tmin=(10,'K'), Tmax=(861.743,'K')), NASAPolynomial(coeffs=[9.78147,0.0200902,-1.43581e-05,4.60177e-09,-5.48256e-13,-104661,-16.1475], Tmin=(861.743,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-861.217,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""ODC(F)C(F)(F)OF""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'O=C(F)C(O)(F)F(435)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {7,S}
4 O u0 p2 c0 {6,S} {8,S}
5 O u0 p2 c0 {7,D}
6 C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
7 C u0 p0 c0 {3,S} {5,D} {6,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-1029.36,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,251,367,519,700,855,1175,1303,486,617,768,1157,1926,180],'cm^-1')),
        HinderedRotor(inertia=(1.01043,'amu*angstrom^2'), symmetry=1, barrier=(23.2317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01274,'amu*angstrom^2'), symmetry=1, barrier=(23.2848,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.023,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.77229,0.0142182,0.000141818,-3.61218e-07,2.49514e-10,-123789,11.5845], Tmin=(10,'K'), Tmax=(533.644,'K')), NASAPolynomial(coeffs=[8.70117,0.0245747,-2.02509e-05,7.35176e-09,-9.67885e-13,-124989,-15.4095], Tmin=(533.644,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1029.36,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""ODC(F)C(O)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OC(F)=C(F)OF(727)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {5,S}
4 O u0 p2 c0 {6,S} {8,S}
5 O u0 p2 c0 {3,S} {7,S}
6 C u0 p0 c0 {1,S} {4,S} {7,D}
7 C u0 p0 c0 {2,S} {5,S} {6,D}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-596.147,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,231,791,262,390,483,597,572,732,631,807,1275,1439,419.165],'cm^-1')),
        HinderedRotor(inertia=(1.25663,'amu*angstrom^2'), symmetry=1, barrier=(28.8923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.023,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.75579,0.0170646,0.000138563,-4.03061e-07,3.18173e-10,-71696.6,11.6484], Tmin=(10,'K'), Tmax=(461.657,'K')), NASAPolynomial(coeffs=[6.83474,0.0284575,-2.21515e-05,7.64998e-09,-9.69279e-13,-72386.5,-5.21985], Tmin=(461.657,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-596.147,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""OC(F)DC(F)OF""", comment="""Thermo library: CHOF_G4"""),
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
    E0 = (-713.441,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-295.211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-52.9147,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-73.9557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-437.808,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (4.22621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-2.39994,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-238.431,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-251.439,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-565.492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-282.088,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-252.495,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-371.242,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-353.821,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-552.234,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-517.902,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-166.569,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['OC(F)(F)C(F)(F)OF(625)'],
    products = ['HF(38)', 'CF2O(48)', 'CF2O(48)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(0.00285451,'s^-1'), n=4.62568, Ea=(21.2505,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CF2(87)', 'OC(F)(F)OF(707)'],
    products = ['OC(F)(F)C(F)(F)OF(625)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.81816e-07,'m^3/(mol*s)'), n=3.45227, Ea=(219.28,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s',), comment="""Estimated from node CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction3',
    reactants = ['OF(153)', 'F[C]C(F)(F)OF(714)'],
    products = ['OC(F)(F)C(F)(F)OF(625)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.13916e+09,'m^3/(mol*s)'), n=-1.00842, Ea=(13.1878,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17406617270594374, var=16.167420070870673, Tref=1000.0, N=4, data_mean=0.0, correlation='OHOY',), comment="""Estimated from node OHOY"""),
)

reaction(
    label = 'reaction4',
    reactants = ['FOF(716)', 'OC(F)(F)[C]F(715)'],
    products = ['OC(F)(F)C(F)(F)OF(625)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4.27832e+09,'m^3/(mol*s)'), n=-1.00842, Ea=(15.3636,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17406617270594374, var=16.167420070870673, Tref=1000.0, N=4, data_mean=0.0, correlation='OHOY',), comment="""Estimated from node OHOY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction5',
    reactants = ['OC(F)(OF)C(F)(F)F(717)'],
    products = ['OC(F)(F)C(F)(F)OF(625)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(6e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
Multiplied by reaction path degeneracy 6.0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction6',
    reactants = ['OF(153)', 'FOC(F)=C(F)F(718)'],
    products = ['OC(F)(F)C(F)(F)OF(625)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.000286045,'m^3/(mol*s)'), n=2.81783, Ea=(231.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OH] for rate rule [Cd/disub_Cd/disub;H_OH]
Euclidian distance = 1.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction7',
    reactants = ['FOF(716)', 'OC(F)=C(F)F(719)'],
    products = ['OC(F)(F)C(F)(F)OF(625)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00057209,'m^3/(mol*s)'), n=2.81783, Ea=(231.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OH] for rate rule [Cd/disub_Cd/disub;H_OH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F(37)', 'O[C](F)C(F)(F)OF(720)'],
    products = ['OC(F)(F)C(F)(F)OF(625)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C
Ea raised from -1.2 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction9',
    reactants = ['F(37)', 'OC(F)(F)[C](F)OF(721)'],
    products = ['OC(F)(F)C(F)(F)OF(625)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C
Ea raised from -0.1 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction10',
    reactants = ['F(37)', '[O]C(F)(F)C(O)(F)F(722)'],
    products = ['OC(F)(F)C(F)(F)OF(625)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.55564e+07,'m^3/(mol*s)'), n=-5.63145e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.7121440562946592, var=2.9508506800589083, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C"""),
)

reaction(
    label = 'reaction11',
    reactants = ['OH(4)', 'FOC(F)(F)[C](F)F(723)'],
    products = ['OC(F)(F)C(F)(F)OF(625)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.7e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(5)', '[O]C(F)(F)C(F)(F)OF(724)'],
    products = ['OC(F)(F)C(F)(F)OF(625)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]F(154)', 'OC(F)(F)[C](F)F(725)'],
    products = ['OC(F)(F)C(F)(F)OF(625)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['FO[C](F)F(192)', 'O[C](F)F(232)'],
    products = ['OC(F)(F)C(F)(F)OF(625)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -11.1 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['HF(38)', 'O=C(F)C(F)(F)OF(726)'],
    products = ['OC(F)(F)C(F)(F)OF(625)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(0.109156,'m^3/(mol*s)'), n=1.86531, Ea=(156.107,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd',), comment="""Estimated from node HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F2(106)', 'O=C(F)C(O)(F)F(435)'],
    products = ['OC(F)(F)C(F)(F)OF(625)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(86.2704,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F2(106)', 'OC(F)=C(F)OF(727)'],
    products = ['OC(F)(F)C(F)(F)OF(625)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(4.39281,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

network(
    label = 'PDepNetwork #312',
    isomers = [
        'OC(F)(F)C(F)(F)OF(625)',
    ],
    reactants = [
        ('HF(38)', 'CF2O(48)', 'CF2O(48)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #312',
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

