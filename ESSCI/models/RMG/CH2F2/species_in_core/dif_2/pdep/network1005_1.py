species(
    label = '[O]OC(F)(F)C(OF)C(O)(F)F(2481)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {11,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {10,S}
7  O u0 p2 c0 {11,S} {14,S}
8  O u0 p2 c0 {9,S} {12,S}
9  O u1 p2 c0 {8,S}
10 C u0 p0 c0 {6,S} {11,S} {12,S} {13,S}
11 C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
12 C u0 p0 c0 {3,S} {4,S} {8,S} {10,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (-1175.46,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,492.5,1135,1000,1380,1390,370,380,2900,435,183,263,327,399,503,589,506,644,608,780,1145,1213,1339,1481,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (197.037,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.7598,0.137651,-0.000236652,2.01471e-07,-6.66302e-11,-141179,38.244], Tmin=(100,'K'), Tmax=(811.848,'K')), NASAPolynomial(coeffs=[17.7566,0.0300013,-1.65237e-05,3.27347e-09,-2.28531e-13,-143969,-49.5112], Tmin=(811.848,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1175.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCFFO) + radical(ROOJ)"""),
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
    label = '[O]OC(F)(F)C=O(898)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {5,S} {6,S}
4 O u0 p2 c0 {7,D}
5 O u1 p2 c0 {3,S}
6 C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
7 C u0 p0 c0 {4,D} {6,S} {8,S}
8 H u0 p0 c0 {7,S}
"""),
    E0 = (-551.124,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,251,367,519,700,855,1175,1303,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.588668,'amu*angstrom^2'), symmetry=1, barrier=(13.5346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.588151,'amu*angstrom^2'), symmetry=1, barrier=(13.5227,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.024,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3731.44,'J/mol'), sigma=(6.05191,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=582.84 K, Pc=38.2 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57192,0.0584795,-9.64588e-05,8.34312e-08,-2.83592e-11,-66202.3,20.9945], Tmin=(100,'K'), Tmax=(811.135,'K')), NASAPolynomial(coeffs=[8.46086,0.0179336,-9.32178e-06,1.82203e-09,-1.26844e-13,-67103.6,-9.46559], Tmin=(811.135,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-551.124,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cds-OdCsH) + radical(ROOJ)"""),
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
    label = '[O]OC(F)(F)C(O)OF(2502)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {5,S}
4  O u0 p2 c0 {8,S} {11,S}
5  O u0 p2 c0 {3,S} {8,S}
6  O u0 p2 c0 {7,S} {9,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
9  C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-719.028,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,492.5,1135,1000,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,209.257,212.821],'cm^-1')),
        HinderedRotor(inertia=(0.269547,'amu*angstrom^2'), symmetry=1, barrier=(8.73881,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.798602,'amu*angstrom^2'), symmetry=1, barrier=(23.6481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.290217,'amu*angstrom^2'), symmetry=1, barrier=(8.73514,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.288593,'amu*angstrom^2'), symmetry=1, barrier=(8.74919,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (147.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0749845,0.0977123,-0.000167461,1.43624e-07,-4.79255e-11,-86340.1,30.3757], Tmin=(100,'K'), Tmax=(812.553,'K')), NASAPolynomial(coeffs=[13.1933,0.0231068,-1.25893e-05,2.48961e-09,-1.73797e-13,-88189.7,-28.9952], Tmin=(812.553,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-719.028,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2sCF) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(CsCFFO) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(OF)C(O)(F)F(2545)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {9,S} {11,S}
7  O u1 p2 c0 {5,S}
8  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
9  C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-719.028,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,492.5,1135,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,208.379,209.043],'cm^-1')),
        HinderedRotor(inertia=(0.282513,'amu*angstrom^2'), symmetry=1, barrier=(8.74025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.282713,'amu*angstrom^2'), symmetry=1, barrier=(8.73949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.282444,'amu*angstrom^2'), symmetry=1, barrier=(8.74344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.763512,'amu*angstrom^2'), symmetry=1, barrier=(23.648,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (147.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0749258,0.0977115,-0.000167458,1.43619e-07,-4.7923e-11,-86340.1,30.3755], Tmin=(100,'K'), Tmax=(812.581,'K')), NASAPolynomial(coeffs=[13.1932,0.023107,-1.25894e-05,2.48963e-09,-1.73799e-13,-88189.7,-28.9946], Tmin=(812.581,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-719.028,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(CsCFFO) + radical(ROOJ)"""),
)

species(
    label = '[O]OF(152)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {3,S}
3 O u1 p2 c0 {2,S}
"""),
    E0 = (13.5955,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(50.9882,'amu')),
        NonlinearRotor(inertia=([6.28174,46.6771,52.9589],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([449.023,679.882,1504.5],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (50.9972,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.63698,0.00939662,-1.79756e-05,1.61481e-08,-5.19278e-12,1646.95,7.91768], Tmin=(100,'K'), Tmax=(982.37,'K')), NASAPolynomial(coeffs=[4.37142,0.00223351,-6.66892e-07,7.82214e-11,-2.87138e-15,1703.99,5.41223], Tmin=(982.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(13.5955,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""OOF""", comment="""Thermo library: halogens"""),
)

species(
    label = 'OC(F)(F)C([C]F)OF-2(1767)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {5,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {3,S} {7,S}
6  O u0 p2 c0 {8,S} {11,S}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
9  C u0 p1 c0 {4,S} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-601.077,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,617,898,1187,256.034,256.074],'cm^-1')),
        HinderedRotor(inertia=(0.210678,'amu*angstrom^2'), symmetry=1, barrier=(9.80607,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.210797,'amu*angstrom^2'), symmetry=1, barrier=(9.80607,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.210785,'amu*angstrom^2'), symmetry=1, barrier=(9.80604,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.682531,'amu*angstrom^2'), symmetry=1, barrier=(31.7592,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.110076,0.0918322,-0.000145105,1.15039e-07,-3.58139e-11,-72158.6,29.1875], Tmin=(100,'K'), Tmax=(790.906,'K')), NASAPolynomial(coeffs=[14.0613,0.0212733,-1.12842e-05,2.23834e-09,-1.57875e-13,-74365.4,-34.8461], Tmin=(790.906,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-601.077,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CJ2_singlet-FCs)"""),
)

species(
    label = '[O]OC(F)(F)C(F)C(O)(F)OF(2523)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {7,S}
6  O u0 p2 c0 {11,S} {14,S}
7  O u0 p2 c0 {5,S} {11,S}
8  O u0 p2 c0 {9,S} {12,S}
9  O u1 p2 c0 {8,S}
10 C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
11 C u0 p0 c0 {2,S} {6,S} {7,S} {10,S}
12 C u0 p0 c0 {3,S} {4,S} {8,S} {10,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-1154.46,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,492.5,1135,1000,250,417,511,1155,1315,1456,3119,375,361,584,565,722,1474,223,363,546,575,694,1179,1410,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (197.037,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.89131,0.140624,-0.000241984,2.05428e-07,-6.76802e-11,-138648,37.6301], Tmin=(100,'K'), Tmax=(813.597,'K')), NASAPolynomial(coeffs=[18.3264,0.0297477,-1.64032e-05,3.24655e-09,-2.26401e-13,-141558,-53.4034], Tmin=(813.597,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1154.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2sCF) + group(O2s-OsH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)(OF)C(F)C(O)(F)F(2546)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {11,S}
7  O u0 p2 c0 {9,S} {11,S}
8  O u0 p2 c0 {12,S} {14,S}
9  O u1 p2 c0 {7,S}
10 C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
11 C u0 p0 c0 {2,S} {6,S} {7,S} {10,S}
12 C u0 p0 c0 {3,S} {4,S} {8,S} {10,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-1154.46,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,492.5,1135,1000,3615,1277.5,1000,250,417,511,1155,1315,1456,3119,375,361,584,565,722,1474,223,363,546,575,694,1179,1410,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (197.037,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.89131,0.140624,-0.000241984,2.05428e-07,-6.76802e-11,-138648,37.6301], Tmin=(100,'K'), Tmax=(813.597,'K')), NASAPolynomial(coeffs=[18.3264,0.0297477,-1.64032e-05,3.24655e-09,-2.26401e-13,-141558,-53.4034], Tmin=(813.597,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1154.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-CsH) + group(O2s-OsH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + radical(ROOJ)"""),
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
    label = '[O]OC(F)(F)C=C(O)F(2547)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {6,S} {7,S}
5  O u0 p2 c0 {9,S} {11,S}
6  O u1 p2 c0 {4,S}
7  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
8  C u0 p0 c0 {7,S} {9,D} {10,S}
9  C u0 p0 c0 {3,S} {5,S} {8,D}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-784.608,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,275,321,533,585,746,850,1103,3010,987.5,1337.5,450,1655,326,540,652,719,1357,180],'cm^-1')),
        HinderedRotor(inertia=(0.622698,'amu*angstrom^2'), symmetry=1, barrier=(14.317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.622414,'amu*angstrom^2'), symmetry=1, barrier=(14.3105,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.623609,'amu*angstrom^2'), symmetry=1, barrier=(14.338,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (143.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.11369,0.0913458,-0.000143804,1.13666e-07,-3.52197e-11,-94232.1,26.9915], Tmin=(100,'K'), Tmax=(795.433,'K')), NASAPolynomial(coeffs=[14.1415,0.0207881,-1.07181e-05,2.09961e-09,-1.46947e-13,-96463.2,-37.4707], Tmin=(795.433,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-784.608,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(CsCFFO) + group(Cds-CdsCsH) + group(CdCFO) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)=CC(O)(F)F(2548)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {7,S} {11,S}
5  O u0 p2 c0 {6,S} {9,S}
6  O u1 p2 c0 {5,S}
7  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
8  C u0 p0 c0 {7,S} {9,D} {10,S}
9  C u0 p0 c0 {3,S} {5,S} {8,D}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-731.618,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,275,321,533,585,746,850,1103,3010,987.5,1337.5,450,1655,326,540,652,719,1357,180],'cm^-1')),
        HinderedRotor(inertia=(0.226307,'amu*angstrom^2'), symmetry=1, barrier=(5.20324,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.227411,'amu*angstrom^2'), symmetry=1, barrier=(5.22864,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.226265,'amu*angstrom^2'), symmetry=1, barrier=(5.20227,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (143.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.397587,0.0910777,-0.000167283,1.56148e-07,-5.5393e-11,-87875,29.2443], Tmin=(100,'K'), Tmax=(848.087,'K')), NASAPolynomial(coeffs=[8.09784,0.030182,-1.6108e-05,3.14045e-09,-2.16431e-13,-88297.3,-1.42524], Tmin=(848.087,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-731.618,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsCFFO) + group(Cds-CdsCsH) + group(CdCFO) + radical(ROOJ)"""),
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
    label = '[O]OC(F)(F)C(OF)=C(F)F(2549)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {7,S}
6  O u0 p2 c0 {8,S} {9,S}
7  O u0 p2 c0 {5,S} {10,S}
8  O u1 p2 c0 {6,S}
9  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
10 C u0 p0 c0 {7,S} {9,S} {11,D}
11 C u0 p0 c0 {3,S} {4,S} {10,D}
"""),
    E0 = (-783.855,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,231,791,275,321,533,585,746,850,1103,350,440,435,1725,182,240,577,636,1210,1413,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.168764,'amu*angstrom^2'), symmetry=1, barrier=(3.88021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168277,'amu*angstrom^2'), symmetry=1, barrier=(3.86902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168449,'amu*angstrom^2'), symmetry=1, barrier=(3.87297,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (179.022,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.588004,0.122068,-0.000250944,2.4724e-07,-9.0124e-11,-94131.1,32.2879], Tmin=(100,'K'), Tmax=(858.18,'K')), NASAPolynomial(coeffs=[6.66795,0.0399394,-2.29563e-05,4.53716e-09,-3.12721e-13,-93597.6,8.75607], Tmin=(858.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-783.855,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-OsH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(Cds-CdsCsOs) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(O)(F)F(327)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {5,S} {6,S}
4 O u0 p2 c0 {6,S} {7,S}
5 O u1 p2 c0 {3,S}
6 C u0 p0 c0 {1,S} {2,S} {3,S} {4,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-655.428,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,435,565,619,662,854,1178,1396],'cm^-1')),
        HinderedRotor(inertia=(0.627674,'amu*angstrom^2'), symmetry=1, barrier=(14.4315,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.623938,'amu*angstrom^2'), symmetry=1, barrier=(14.3456,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.0136,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3791.54,'J/mol'), sigma=(6.20894,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=592.23 K, Pc=35.94 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.89632,0.0496154,-7.99025e-05,6.55003e-08,-2.10556e-11,-78757.2,18.1468], Tmin=(100,'K'), Tmax=(804.92,'K')), NASAPolynomial(coeffs=[9.02522,0.0121076,-6.12679e-06,1.18408e-09,-8.19294e-14,-79837.5,-14.2802], Tmin=(804.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-655.428,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(CsFFOO) + radical(ROOJ)"""),
)

species(
    label = 'FOC=C(F)F(1226)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {4,S}
4 O u0 p2 c0 {3,S} {5,S}
5 C u0 p0 c0 {4,S} {6,D} {7,S}
6 C u0 p0 c0 {1,S} {2,S} {5,D}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-380.308,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,1091.78],'cm^-1')),
        HinderedRotor(inertia=(0.655419,'amu*angstrom^2'), symmetry=1, barrier=(29.8096,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0239,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.79505,0.0174471,4.03379e-05,-1.08299e-07,6.89768e-11,-45737.6,11.4304], Tmin=(10,'K'), Tmax=(583.249,'K')), NASAPolynomial(coeffs=[5.4314,0.0228414,-1.62698e-05,5.25091e-09,-6.31965e-13,-46211.1,1.99525], Tmin=(583.249,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-380.308,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""FOCDC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    collisionModel = TransportData(shapeIndex=0, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.05029e-14,-1.20495e-16,5.23432e-20,-6.99814e-24,29230.2,5.12616], Tmin=(100,'K'), Tmax=(4230.24,'K')), NASAPolynomial(coeffs=[2.5,3.58763e-09,-1.22015e-12,1.84118e-16,-1.04e-20,29230.2,5.12619], Tmin=(4230.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.034,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[O]C(F)(F)C(OF)C(O)(F)F(2550)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {9,S}
7  O u0 p2 c0 {10,S} {13,S}
8  O u1 p2 c0 {11,S}
9  C u0 p0 c0 {6,S} {10,S} {11,S} {12,S}
10 C u0 p0 c0 {1,S} {2,S} {7,S} {9,S}
11 C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-1139.01,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,351,323,533,609,664,892,1120,1201,265.267,265.267,265.268],'cm^-1')),
        HinderedRotor(inertia=(0.173432,'amu*angstrom^2'), symmetry=1, barrier=(8.66011,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.173431,'amu*angstrom^2'), symmetry=1, barrier=(8.6601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.361069,'amu*angstrom^2'), symmetry=1, barrier=(18.0295,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.173431,'amu*angstrom^2'), symmetry=1, barrier=(8.6601,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (181.038,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.00996,0.118741,-0.000193389,1.57016e-07,-4.99867e-11,-136818,34.6288], Tmin=(100,'K'), Tmax=(774.064,'K')), NASAPolynomial(coeffs=[16.8557,0.0264106,-1.44524e-05,2.89148e-09,-2.04453e-13,-139584,-46.9858], Tmin=(774.064,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1139.01,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCFFO) + radical(O2sj(Cs-CsF1sF1s))"""),
)

species(
    label = 'HO2(10)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,S} {3,S}
2 O u1 p2 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (2.49012,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1064.4,1465.7,3224.93],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (33.0068,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.3018,-0.00474912,2.11583e-05,-2.42764e-08,9.29225e-12,264.018,3.71666], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.17229,0.00188118,-3.46277e-07,1.94658e-11,1.76257e-16,31.0207,2.95768], Tmin=(1000,'K'), Tmax=(5000,'K'))], Tmin=(200,'K'), Tmax=(5000,'K'), E0=(2.49012,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HO2""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'OC(F)(F)C(OF)=C(F)F(1758)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {7,S}
6  O u0 p2 c0 {8,S} {11,S}
7  O u0 p2 c0 {5,S} {9,S}
8  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
9  C u0 p0 c0 {7,S} {8,S} {10,D}
10 C u0 p0 c0 {3,S} {4,S} {9,D}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-1007.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,231,791,275,321,533,585,746,850,1103,350,440,435,1725,182,240,577,636,1210,1413,180,1145.53],'cm^-1')),
        HinderedRotor(inertia=(0.224111,'amu*angstrom^2'), symmetry=1, barrier=(5.15274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.226473,'amu*angstrom^2'), symmetry=1, barrier=(5.20706,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.221901,'amu*angstrom^2'), symmetry=1, barrier=(5.10194,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (164.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.151901,0.108442,-0.000214064,2.08221e-07,-7.58015e-11,-121025,29.202], Tmin=(100,'K'), Tmax=(849.09,'K')), NASAPolynomial(coeffs=[7.03244,0.0368309,-2.08393e-05,4.12631e-09,-2.85898e-13,-120884,3.73368], Tmin=(849.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1007.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(Cds-CdsCsOs) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs)"""),
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
    label = '[O]OC(F)(F)C(OF)[C](O)F(2551)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {8,S} {10,S}
7  O u0 p2 c0 {11,S} {13,S}
8  O u1 p2 c0 {6,S}
9  C u0 p0 c0 {5,S} {10,S} {11,S} {12,S}
10 C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
11 C u1 p0 c0 {3,S} {7,S} {9,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-746.203,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,492.5,1135,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,395,473,707,1436,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (178.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.38142,0.13249,-0.000243476,2.175e-07,-7.39005e-11,-89567.2,37.4799], Tmin=(100,'K'), Tmax=(851.872,'K')), NASAPolynomial(coeffs=[14.6475,0.0304965,-1.68175e-05,3.28668e-09,-2.25707e-13,-91328.3,-31.5883], Tmin=(851.872,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-746.203,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCFHO) + radical(ROOJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]O[C](F)C(OF)C(O)(F)F(2552)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {10,S} {13,S}
7  O u0 p2 c0 {8,S} {11,S}
8  O u1 p2 c0 {7,S}
9  C u0 p0 c0 {5,S} {10,S} {11,S} {12,S}
10 C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
11 C u1 p0 c0 {3,S} {7,S} {9,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-746.203,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,492.5,1135,1000,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,395,473,707,1436,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (178.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.38141,0.13249,-0.000243475,2.175e-07,-7.39003e-11,-89567.2,37.4799], Tmin=(100,'K'), Tmax=(851.872,'K')), NASAPolynomial(coeffs=[14.6475,0.0304965,-1.68175e-05,3.28669e-09,-2.25708e-13,-91328.3,-31.5882], Tmin=(851.872,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-746.203,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(CsCFHO) + group(CsCFFO) + radical(ROOJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]OC(F)(F)C([O])C(O)(F)F(2553)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {10,S} {13,S}
6  O u0 p2 c0 {8,S} {11,S}
7  O u1 p2 c0 {9,S}
8  O u1 p2 c0 {6,S}
9  C u0 p0 c0 {7,S} {10,S} {11,S} {12,S}
10 C u0 p0 c0 {1,S} {2,S} {5,S} {9,S}
11 C u0 p0 c0 {3,S} {4,S} {6,S} {9,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (-1079.61,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,1380,1390,370,380,2900,435,183,263,327,399,503,589,506,644,608,780,1145,1213,1339,1481,249.132,249.14,249.141],'cm^-1')),
        HinderedRotor(inertia=(0.199917,'amu*angstrom^2'), symmetry=1, barrier=(8.8059,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.199902,'amu*angstrom^2'), symmetry=1, barrier=(8.80584,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.537568,'amu*angstrom^2'), symmetry=1, barrier=(23.6784,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.199905,'amu*angstrom^2'), symmetry=1, barrier=(8.80584,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (178.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.09178,0.120939,-0.000202725,1.69742e-07,-5.55141e-11,-129672,35.4845], Tmin=(100,'K'), Tmax=(797.519,'K')), NASAPolynomial(coeffs=[16.4184,0.0270142,-1.4593e-05,2.88409e-09,-2.01651e-13,-132271,-43.8136], Tmin=(797.519,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1079.61,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCFFO) + radical(CC(C)OJ) + radical(ROOJ)"""),
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
    label = '[O]OC(F)(F)[CH]C(O)(F)F(2554)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {8,S} {12,S}
6  O u0 p2 c0 {7,S} {9,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {1,S} {2,S} {5,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {6,S} {10,S}
10 C u1 p0 c0 {8,S} {9,S} {11,S}
11 H u0 p0 c0 {10,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-937.086,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,202,304,469,581,520,674,579,755,732,952,1136,1220,1240,1408,3025,407.5,1350,352.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.409556,'amu*angstrom^2'), symmetry=1, barrier=(9.4165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.411756,'amu*angstrom^2'), symmetry=1, barrier=(9.46708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.41101,'amu*angstrom^2'), symmetry=1, barrier=(9.44994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.602518,'amu*angstrom^2'), symmetry=1, barrier=(13.8531,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (162.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.537232,0.109842,-0.00019177,1.67516e-07,-5.66054e-11,-112552,33.6305], Tmin=(100,'K'), Tmax=(826.59,'K')), NASAPolynomial(coeffs=[13.4299,0.0272957,-1.48328e-05,2.92078e-09,-2.02987e-13,-114350,-28.0018], Tmin=(826.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-937.086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(CsCFFO) + group(CsCFFO) + radical(ROOJ) + radical(CCJCOOH)"""),
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
    label = '[O]OC(F)(F)C(OF)[C](F)F(2555)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {9,S}
7  O u0 p2 c0 {8,S} {10,S}
8  O u1 p2 c0 {7,S}
9  C u0 p0 c0 {6,S} {10,S} {11,S} {12,S}
10 C u0 p0 c0 {1,S} {2,S} {7,S} {9,S}
11 C u1 p0 c0 {3,S} {4,S} {9,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-763.692,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,492.5,1135,1000,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,190,488,555,1236,1407,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.429287,'amu*angstrom^2'), symmetry=1, barrier=(9.87014,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.42951,'amu*angstrom^2'), symmetry=1, barrier=(9.87529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.429394,'amu*angstrom^2'), symmetry=1, barrier=(9.87262,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49595,'amu*angstrom^2'), symmetry=1, barrier=(34.3947,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (180.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.10235,0.124782,-0.000226907,2.01562e-07,-6.85845e-11,-91679.3,35.5747], Tmin=(100,'K'), Tmax=(836.517,'K')), NASAPolynomial(coeffs=[14.5898,0.0284886,-1.61195e-05,3.19491e-09,-2.21831e-13,-93560.9,-32.8841], Tmin=(836.517,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-763.692,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCsFFH) + radical(ROOJ) + radical(Csj(Cs-O2sCsH)(F1s)(F1s))"""),
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
    label = '[O]OC(F)(F)C(OF)C([O])(F)F(2556)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {11,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {10,S}
7  O u0 p2 c0 {9,S} {11,S}
8  O u1 p2 c0 {12,S}
9  O u1 p2 c0 {7,S}
10 C u0 p0 c0 {6,S} {11,S} {12,S} {13,S}
11 C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
12 C u0 p0 c0 {3,S} {4,S} {8,S} {10,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-915.496,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,492.5,1135,1000,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,351,323,533,609,664,892,1120,1201,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.383008,'amu*angstrom^2'), symmetry=1, barrier=(8.8061,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.383012,'amu*angstrom^2'), symmetry=1, barrier=(8.80621,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.383291,'amu*angstrom^2'), symmetry=1, barrier=(8.81262,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.966854,'amu*angstrom^2'), symmetry=1, barrier=(22.2299,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (196.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.54995,0.133724,-0.000235685,2.04119e-07,-6.82808e-11,-109920,38.0791], Tmin=(100,'K'), Tmax=(823.106,'K')), NASAPolynomial(coeffs=[16.7624,0.0290287,-1.6274e-05,3.23031e-09,-2.25159e-13,-112403,-43.471], Tmin=(823.106,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-915.496,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCFFO) + radical(O2sj(Cs-CsF1sF1s)) + radical(ROOJ)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53732,-0.00121572,5.3162e-06,-4.89446e-09,1.45846e-12,-1038.59,4.68368], Tmin=(100,'K'), Tmax=(1074.55,'K')), NASAPolynomial(coeffs=[3.15382,0.00167804,-7.69974e-07,1.51275e-10,-1.08782e-14,-1040.82,6.16756], Tmin=(1074.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.62683,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'OC(F)(F)C(OF)[C](F)F(1668)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {8,S}
7  O u0 p2 c0 {9,S} {12,S}
8  C u0 p0 c0 {6,S} {9,S} {10,S} {11,S}
9  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
10 C u1 p0 c0 {3,S} {4,S} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-987.202,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,190,488,555,1236,1407,296.884,296.884,296.884],'cm^-1')),
        HinderedRotor(inertia=(0.161204,'amu*angstrom^2'), symmetry=1, barrier=(10.0827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161204,'amu*angstrom^2'), symmetry=1, barrier=(10.0827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161204,'amu*angstrom^2'), symmetry=1, barrier=(10.0827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.402056,'amu*angstrom^2'), symmetry=1, barrier=(25.147,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (165.039,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3545.48,'J/mol'), sigma=(5.86849,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=553.80 K, Pc=39.81 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.635244,0.110766,-0.000188551,1.60495e-07,-5.3346e-11,-118574,32.3793], Tmin=(100,'K'), Tmax=(797.711,'K')), NASAPolynomial(coeffs=[14.842,0.0255805,-1.4122e-05,2.81298e-09,-1.97453e-13,-120803,-37.2803], Tmin=(797.711,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-987.202,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCsFFH) + radical(Csj(Cs-O2sCsH)(F1s)(F1s))"""),
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
    label = '[O]OC(F)(F)[CH]OF(1523)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {5,S}
4 O u0 p2 c0 {6,S} {7,S}
5 O u0 p2 c0 {3,S} {8,S}
6 O u1 p2 c0 {4,S}
7 C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
8 C u1 p0 c0 {5,S} {7,S} {9,S}
9 H u0 p0 c0 {8,S}
"""),
    E0 = (-335.075,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,557,1111,253,525,597,667,842,1178,1324,3025,407.5,1350,352.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.729116,'amu*angstrom^2'), symmetry=1, barrier=(16.7638,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.728358,'amu*angstrom^2'), symmetry=1, barrier=(16.7464,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.72859,'amu*angstrom^2'), symmetry=1, barrier=(16.7517,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.023,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.353357,0.0868688,-0.000151317,1.27876e-07,-4.17026e-11,-40175.2,25.8821], Tmin=(100,'K'), Tmax=(824.175,'K')), NASAPolynomial(coeffs=[13.3943,0.0161223,-8.99081e-06,1.7756e-09,-1.23257e-13,-42071.6,-32.9752], Tmin=(824.175,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-335.075,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-OsH) + group(CsCFFO) + group(Cs-CsOsHH) + radical(ROOJ) + radical(CCsJO)"""),
)

species(
    label = '[O]O[C](F)F(171)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {4,S} {5,S}
4 O u1 p2 c0 {3,S}
5 C u1 p0 c0 {1,S} {2,S} {3,S}
"""),
    E0 = (-195.642,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,493,600,700,1144,1293],'cm^-1')),
        HinderedRotor(inertia=(0.711545,'amu*angstrom^2'), symmetry=1, barrier=(16.3598,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.0062,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3217.18,'J/mol'), sigma=(5.19785,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=502.52 K, Pc=51.98 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.8245,0.013883,6.3822e-05,-2.51617e-07,2.43597e-10,-23530.5,10.5618], Tmin=(10,'K'), Tmax=(394.515,'K')), NASAPolynomial(coeffs=[5.8912,0.0131528,-1.02966e-05,3.57092e-09,-4.54892e-13,-23851,0.518718], Tmin=(394.515,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-195.642,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""[O]O[C](F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]OC(F)(F)[C](OF)C(O)(F)F(2557)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {8,S}
6  O u0 p2 c0 {10,S} {13,S}
7  O u0 p2 c0 {9,S} {11,S}
8  O u0 p2 c0 {5,S} {12,S}
9  O u1 p2 c0 {7,S}
10 C u0 p0 c0 {1,S} {2,S} {6,S} {12,S}
11 C u0 p0 c0 {3,S} {4,S} {7,S} {12,S}
12 C u1 p0 c0 {8,S} {10,S} {11,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-995.359,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,557,1111,202,304,469,581,520,674,579,755,732,952,1136,1220,1240,1408,360,370,350,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (196.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.80162,0.144221,-0.000273068,2.47062e-07,-8.44354e-11,-119521,38.8131], Tmin=(100,'K'), Tmax=(856.815,'K')), NASAPolynomial(coeffs=[15.2644,0.0315073,-1.78991e-05,3.51407e-09,-2.41071e-13,-121233,-33.8055], Tmin=(856.815,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-995.359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCFFO) + radical(ROOJ) + radical(C2CsJO)"""),
)

species(
    label = '[O]OC(F)(F)C(OF)C(=O)F(2558)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {8,S} {10,S}
7  O u0 p2 c0 {11,D}
8  O u1 p2 c0 {6,S}
9  C u0 p0 c0 {5,S} {10,S} {11,S} {12,S}
10 C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
11 C u0 p0 c0 {3,S} {7,D} {9,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-889.423,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,492.5,1135,1000,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,486,617,768,1157,1926,310.809,310.869,310.874],'cm^-1')),
        HinderedRotor(inertia=(0.115404,'amu*angstrom^2'), symmetry=1, barrier=(7.9189,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.233085,'amu*angstrom^2'), symmetry=1, barrier=(15.9676,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115491,'amu*angstrom^2'), symmetry=1, barrier=(7.92048,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.454947,'amu*angstrom^2'), symmetry=1, barrier=(31.175,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (177.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.230265,0.10026,-0.000146675,1.07139e-07,-3.09309e-11,-106827,35.1067], Tmin=(100,'K'), Tmax=(848.579,'K')), NASAPolynomial(coeffs=[15.702,0.0251604,-1.39243e-05,2.84728e-09,-2.05865e-13,-109531,-39.1412], Tmin=(848.579,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-889.423,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCFFO) + group(COCsFO) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)(F)C(=O)C(O)(F)F(2559)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {9,S} {12,S}
6  O u0 p2 c0 {8,S} {10,S}
7  O u0 p2 c0 {11,D}
8  O u1 p2 c0 {6,S}
9  C u0 p0 c0 {1,S} {2,S} {5,S} {11,S}
10 C u0 p0 c0 {3,S} {4,S} {6,S} {11,S}
11 C u0 p0 c0 {7,D} {9,S} {10,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-1202.72,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,218,284,327,407,463,575,605,795,835,875,1140,1210,1216,1390,375,552.5,462.5,1710,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.441235,'amu*angstrom^2'), symmetry=1, barrier=(10.1449,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.441138,'amu*angstrom^2'), symmetry=1, barrier=(10.1426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.441215,'amu*angstrom^2'), symmetry=1, barrier=(10.1444,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.44097,'amu*angstrom^2'), symmetry=1, barrier=(10.1388,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (177.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.902432,0.11971,-0.000216294,1.90663e-07,-6.41383e-11,-144489,34.1643], Tmin=(100,'K'), Tmax=(850.726,'K')), NASAPolynomial(coeffs=[14.2576,0.02688,-1.46203e-05,2.84614e-09,-1.95203e-13,-146288,-31.9398], Tmin=(850.726,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1202.72,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cds-OdCsCs) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)(F)C(OF)=C(O)F(2560)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {6,S}
5  O u0 p2 c0 {8,S} {9,S}
6  O u0 p2 c0 {4,S} {10,S}
7  O u0 p2 c0 {11,S} {12,S}
8  O u1 p2 c0 {5,S}
9  C u0 p0 c0 {1,S} {2,S} {5,S} {10,S}
10 C u0 p0 c0 {6,S} {9,S} {11,D}
11 C u0 p0 c0 {3,S} {7,S} {10,D}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-793.663,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,231,791,3615,1277.5,1000,275,321,533,585,746,850,1103,350,440,435,1725,326,540,652,719,1357,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.301667,'amu*angstrom^2'), symmetry=1, barrier=(6.93591,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.301407,'amu*angstrom^2'), symmetry=1, barrier=(6.92995,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.301594,'amu*angstrom^2'), symmetry=1, barrier=(6.93424,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.41572,'amu*angstrom^2'), symmetry=1, barrier=(55.5421,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (177.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.0905,0.13024,-0.000254754,2.40051e-07,-8.48162e-11,-95289.9,33.4297], Tmin=(100,'K'), Tmax=(856.097,'K')), NASAPolynomial(coeffs=[10.79,0.0355628,-2.02406e-05,3.9883e-09,-2.74557e-13,-95888.8,-13.6585], Tmin=(856.097,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-793.663,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(Cds-CdsCsOs) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)=C(OF)C(O)(F)F(2561)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {6,S}
5  O u0 p2 c0 {9,S} {12,S}
6  O u0 p2 c0 {4,S} {10,S}
7  O u0 p2 c0 {8,S} {11,S}
8  O u1 p2 c0 {7,S}
9  C u0 p0 c0 {1,S} {2,S} {5,S} {10,S}
10 C u0 p0 c0 {6,S} {9,S} {11,D}
11 C u0 p0 c0 {3,S} {7,S} {10,D}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-740.674,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,231,791,492.5,1135,1000,275,321,533,585,746,850,1103,350,440,435,1725,326,540,652,719,1357,180,2373.91],'cm^-1')),
        HinderedRotor(inertia=(0.264676,'amu*angstrom^2'), symmetry=1, barrier=(6.08541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.264608,'amu*angstrom^2'), symmetry=1, barrier=(6.08386,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.264755,'amu*angstrom^2'), symmetry=1, barrier=(6.08723,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.264745,'amu*angstrom^2'), symmetry=1, barrier=(6.087,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (177.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.710855,0.128729,-0.000273318,2.75261e-07,-1.01447e-10,-88936.9,35.3462], Tmin=(100,'K'), Tmax=(864.173,'K')), NASAPolynomial(coeffs=[4.49741,0.0454046,-2.58994e-05,5.09454e-09,-3.49584e-13,-87625.9,23.7731], Tmin=(864.173,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-740.674,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2sCF) + group(O2s-OsH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(Cds-CdsCsOs) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)=COF(629)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {5,S} {6,S}
4 O u0 p2 c0 {2,S} {7,S}
5 O u1 p2 c0 {3,S}
6 C u0 p0 c0 {1,S} {3,S} {7,D}
7 C u0 p0 c0 {4,S} {6,D} {8,S}
8 H u0 p0 c0 {7,S}
"""),
    E0 = (-99.5423,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,231,791,326,540,652,719,1357,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(0.213515,'amu*angstrom^2'), symmetry=1, barrier=(4.90913,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.213911,'amu*angstrom^2'), symmetry=1, barrier=(4.91824,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.024,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57387,0.0632844,-0.000121339,1.17857e-07,-4.29602e-11,-11894.4,23.0683], Tmin=(100,'K'), Tmax=(852.006,'K')), NASAPolynomial(coeffs=[5.2683,0.0238131,-1.28928e-05,2.52049e-09,-1.73717e-13,-11720.8,10.5496], Tmin=(852.006,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-99.5423,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2sCF) + group(O2s-OsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + radical(ROOJ)"""),
)

species(
    label = 'OOC(F)(F)[C](OF)C(O)(F)F(2562)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {11,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {8,S}
6  O u0 p2 c0 {9,S} {10,S}
7  O u0 p2 c0 {11,S} {13,S}
8  O u0 p2 c0 {5,S} {12,S}
9  O u0 p2 c0 {6,S} {14,S}
10 C u0 p0 c0 {3,S} {4,S} {6,S} {12,S}
11 C u0 p0 c0 {1,S} {2,S} {7,S} {12,S}
12 C u1 p0 c0 {8,S} {10,S} {11,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-1147.36,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (197.037,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.04926,0.147628,-0.000267804,2.37275e-07,-8.06511e-11,-137792,38.9428], Tmin=(100,'K'), Tmax=(832.816,'K')), NASAPolynomial(coeffs=[16.8183,0.0332282,-1.89258e-05,3.76335e-09,-2.61916e-13,-140110,-43.6793], Tmin=(832.816,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1147.36,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCFFO) + radical(C2CsJO)"""),
)

species(
    label = '[O]C(F)(F)C(OF)C(F)(F)OO(2563)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {11,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {7,S}
6  O u0 p2 c0 {8,S} {11,S}
7  O u0 p2 c0 {5,S} {10,S}
8  O u0 p2 c0 {6,S} {14,S}
9  O u1 p2 c0 {12,S}
10 C u0 p0 c0 {7,S} {11,S} {12,S} {13,S}
11 C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
12 C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-1067.5,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,557,1111,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,351,323,533,609,664,892,1120,1201,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (197.037,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.57905,0.134142,-0.00021767,1.73607e-07,-5.32223e-11,-128200,37.4492], Tmin=(100,'K'), Tmax=(680.539,'K')), NASAPolynomial(coeffs=[17.9247,0.0314699,-1.77397e-05,3.58751e-09,-2.55225e-13,-131132,-51.1742], Tmin=(680.539,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1067.5,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCFFO) + radical(O2sj(Cs-CsF1sF1s))"""),
)

species(
    label = '[O]C(C(O)(F)F)C(F)(F)OOF(2564)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {12,S}
2  F u0 p3 c0 {12,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {8,S}
6  O u0 p2 c0 {8,S} {11,S}
7  O u0 p2 c0 {12,S} {14,S}
8  O u0 p2 c0 {5,S} {6,S}
9  O u1 p2 c0 {10,S}
10 C u0 p0 c0 {9,S} {11,S} {12,S} {13,S}
11 C u0 p0 c0 {3,S} {4,S} {6,S} {10,S}
12 C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (-1136.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3615,1277.5,1000,277,555,632,1380,1390,370,380,2900,435,183,263,327,399,503,589,506,644,608,780,1145,1213,1339,1481],'cm^-1')),
        HinderedRotor(inertia=(1.06912,'amu*angstrom^2'), symmetry=1, barrier=(24.5811,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06917,'amu*angstrom^2'), symmetry=1, barrier=(24.5824,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06919,'amu*angstrom^2'), symmetry=1, barrier=(24.5827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06942,'amu*angstrom^2'), symmetry=1, barrier=(24.588,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06979,'amu*angstrom^2'), symmetry=1, barrier=(24.5966,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (197.037,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.54671,0.132005,-0.000212111,1.69927e-07,-5.35372e-11,-136477,37.7529], Tmin=(100,'K'), Tmax=(780.888,'K')), NASAPolynomial(coeffs=[18.2681,0.0305065,-1.71421e-05,3.47672e-09,-2.48309e-13,-139571,-52.9413], Tmin=(780.888,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1136.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2sFO) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCFFO) + radical(CC(C)OJ)"""),
)

species(
    label = 'O[C](F)C(OF)C(F)(F)OOF(2565)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {11,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {7,S}
5  F u0 p3 c0 {9,S}
6  O u0 p2 c0 {9,S} {11,S}
7  O u0 p2 c0 {4,S} {10,S}
8  O u0 p2 c0 {12,S} {14,S}
9  O u0 p2 c0 {5,S} {6,S}
10 C u0 p0 c0 {7,S} {11,S} {12,S} {13,S}
11 C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
12 C u1 p0 c0 {3,S} {8,S} {10,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-802.907,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,557,1111,3615,1277.5,1000,277,555,632,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,395,473,707,1436],'cm^-1')),
        HinderedRotor(inertia=(0.865839,'amu*angstrom^2'), symmetry=1, barrier=(19.9073,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.868848,'amu*angstrom^2'), symmetry=1, barrier=(19.9765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.869184,'amu*angstrom^2'), symmetry=1, barrier=(19.9842,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.33292,'amu*angstrom^2'), symmetry=1, barrier=(30.6465,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.869436,'amu*angstrom^2'), symmetry=1, barrier=(19.99,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.32614,'amu*angstrom^2'), symmetry=1, barrier=(30.4905,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (197.037,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.98825,0.145563,-0.000260994,2.30062e-07,-7.81423e-11,-96365,40.28], Tmin=(100,'K'), Tmax=(824.874,'K')), NASAPolynomial(coeffs=[16.8466,0.0333527,-1.89816e-05,3.78516e-09,-2.64352e-13,-98762.1,-42.6563], Tmin=(824.874,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-802.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-CsH) + group(O2sFO) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCFHO) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'OC(F)(F)C(OF)[C](F)OOF(2566)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {11,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {6,S}
5  F u0 p3 c0 {9,S}
6  O u0 p2 c0 {4,S} {10,S}
7  O u0 p2 c0 {11,S} {14,S}
8  O u0 p2 c0 {9,S} {12,S}
9  O u0 p2 c0 {5,S} {8,S}
10 C u0 p0 c0 {6,S} {11,S} {12,S} {13,S}
11 C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
12 C u1 p0 c0 {3,S} {8,S} {10,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (-802.907,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,3615,1310,387.5,850,1000,277,555,632,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,395,473,707,1436],'cm^-1')),
        HinderedRotor(inertia=(0.920903,'amu*angstrom^2'), symmetry=1, barrier=(21.1734,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.91835,'amu*angstrom^2'), symmetry=1, barrier=(21.1147,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.91869,'amu*angstrom^2'), symmetry=1, barrier=(21.1225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.55041,'amu*angstrom^2'), symmetry=1, barrier=(35.647,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.919765,'amu*angstrom^2'), symmetry=1, barrier=(21.1472,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.918381,'amu*angstrom^2'), symmetry=1, barrier=(21.1154,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (197.037,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.98856,0.145567,-0.000261011,2.30087e-07,-7.81551e-11,-96365,40.2811], Tmin=(100,'K'), Tmax=(824.805,'K')), NASAPolynomial(coeffs=[16.8472,0.0333517,-1.8981e-05,3.78501e-09,-2.64339e-13,-98762.3,-42.6596], Tmin=(824.805,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-802.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-CsH) + group(O2sFO) + group(Cs-CsCsOsH) + group(CsCFHO) + group(CsCFFO) + radical(CsCsF1sO2s)"""),
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
    E0 = (-595.038,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-173.671,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-173.671,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-22.8532,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-312.086,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-312.086,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (6.58614,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (59.576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-260.441,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-163.081,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-352.596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-457.299,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-129.936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-129.936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-463.344,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-291.025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-191.923,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-160.316,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-452.454,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-261.993,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-223.745,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-240.179,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-460.001,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-496.679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-352.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-320.869,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-347.971,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-528.241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-454.583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-413.596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (-198.057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-175.744,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]OC(F)(F)C(OF)C(O)(F)F(2481)'],
    products = ['HF(38)', 'CF2O(48)', '[O]OC(F)(F)C=O(898)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(0.00285451,'s^-1'), n=4.62568, Ea=(37.0509,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CF2(87)', '[O]OC(F)(F)C(O)OF(2502)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2481)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.40908e-07,'m^3/(mol*s)'), n=3.45227, Ea=(205.694,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s',), comment="""Estimated from node CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CF2(87)', '[O]OC(OF)C(O)(F)F(2545)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2481)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.40908e-07,'m^3/(mol*s)'), n=3.45227, Ea=(205.694,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s',), comment="""Estimated from node CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]OF(152)', 'OC(F)(F)C([C]F)OF-2(1767)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2481)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.13916e+09,'m^3/(mol*s)'), n=-1.00842, Ea=(21.2535,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17406617270594374, var=16.167420070870673, Tref=1000.0, N=4, data_mean=0.0, correlation='OHOY',), comment="""Estimated from node OHOY"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]OC(F)(F)C(F)C(O)(F)OF(2523)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2481)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]OC(F)(OF)C(F)C(O)(F)F(2546)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2481)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction7',
    reactants = ['FOF(716)', '[O]OC(F)(F)C=C(O)F(2547)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2481)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00057209,'m^3/(mol*s)'), n=2.81783, Ea=(231.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OH] for rate rule [Cd/disub_Cd/monosub;H_OH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction8',
    reactants = ['FOF(716)', '[O]OC(F)=CC(O)(F)F(2548)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2481)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.00057209,'m^3/(mol*s)'), n=2.81783, Ea=(231.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OH] for rate rule [Cd/disub_Cd/monosub;H_OH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H2O(3)', '[O]OC(F)(F)C(OF)=C(F)F(2549)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2481)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.00057209,'m^3/(mol*s)'), n=2.81783, Ea=(231.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OH] for rate rule [Cd/disub_Cd/disub;H_OH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]OC(O)(F)F(327)', 'FOC=C(F)F(1226)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2481)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.79e-11,'m^3/(mol*s)'), n=3.97, Ea=(329.281,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;R_OH] for rate rule [Cd/monosub_Cd/disub;R_OH]
Euclidian distance = 1.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction11',
    reactants = ['O(6)', '[O]C(F)(F)C(OF)C(O)(F)F(2550)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2481)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]OC(F)(F)C(OF)C(O)(F)F(2481)'],
    products = ['HO2(10)', 'OC(F)(F)C(OF)=C(F)F(1758)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.00545e+13,'s^-1'), n=-0.00666667, Ea=(174.79,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO_NdNd]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', '[O]OC(F)(F)C(OF)[C](O)F(2551)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2481)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F(37)', '[O]O[C](F)C(OF)C(O)(F)F(2552)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2481)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F(37)', '[O]OC(F)(F)C([O])C(O)(F)F(2553)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2481)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.49151e+08,'m^3/(mol*s)'), n=0.0472428, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]F(154)', '[O]OC(F)(F)[CH]C(O)(F)F(2554)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2481)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['OH(4)', '[O]OC(F)(F)C(OF)[C](F)F(2555)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2481)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(5)', '[O]OC(F)(F)C(OF)C([O])(F)F(2556)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2481)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction19',
    reactants = ['O2(2)', 'OC(F)(F)C(OF)[C](F)F(1668)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2481)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.6844e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction20',
    reactants = ['O[C](F)F(232)', '[O]OC(F)(F)[CH]OF(1523)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2481)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]O[C](F)F(171)', 'OC(F)(F)[CH]OF(1227)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2481)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(5)', '[O]OC(F)(F)[C](OF)C(O)(F)F(2557)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2481)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(43950,'m^3/(mol*s)'), n=1, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction23',
    reactants = ['HF(38)', '[O]OC(F)(F)C(OF)C(=O)F(2558)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2481)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(0.109156,'m^3/(mol*s)'), n=1.86531, Ea=(167.161,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd',), comment="""Estimated from node HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd"""),
)

reaction(
    label = 'reaction24',
    reactants = ['HF(38)', '[O]OC(F)(F)C(=O)C(O)(F)F(2559)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2481)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(443.781,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction25',
    reactants = ['HF(38)', '[O]OC(F)(F)C(OF)=C(O)F(2560)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2481)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.14111,'m^3/(mol*s)'), n=1.29695, Ea=(179.265,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction26',
    reactants = ['HF(38)', '[O]OC(F)=C(OF)C(O)(F)F(2561)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2481)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.14111,'m^3/(mol*s)'), n=1.29695, Ea=(157.543,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O]OC(F)(F)C(OF)C(O)(F)F(2481)'],
    products = ['HF(38)', 'CF2O(48)', '[O]OC(F)=COF(629)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.97084e+10,'s^-1'), n=0.978818, Ea=(284.118,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R!H->C_N-5Br1sCl1sF1sH->H_5Br1sCl1sF1s->F1s_Ext-2C-R_N-7R!H->O_Sp-7BrCClFINPSSi-2C_Ext-2C-R',), comment="""Estimated from node Root_1R!H->C_N-5Br1sCl1sF1sH->H_5Br1sCl1sF1s->F1s_Ext-2C-R_N-7R!H->O_Sp-7BrCClFINPSSi-2C_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[O]OC(F)(F)C(OF)C(O)(F)F(2481)'],
    products = ['OOC(F)(F)[C](OF)C(O)(F)F(2562)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(5.29e+09,'s^-1'), n=0.75, Ea=(103.847,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 272 used for R4H_SSS_OCs;O_rad_out;Cs_H_out_NDMustO
Exact match found for rate rule [R4H_SSS_OCs;O_rad_out;Cs_H_out_NDMustO]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[O]C(F)(F)C(OF)C(F)(F)OO(2563)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2481)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(123903,'s^-1'), n=1.46258, Ea=(69.5427,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_SSSSS;O_rad_out;XH_out] for rate rule [R6H_SSSSS;O_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[O]C(C(O)(F)F)C(F)(F)OOF(2564)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2481)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(179.343,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction31',
    reactants = ['O[C](F)C(OF)C(F)(F)OOF(2565)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2481)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(61.4745,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction32',
    reactants = ['OC(F)(F)C(OF)[C](F)OOF(2566)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2481)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(83.7872,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

network(
    label = 'PDepNetwork #1005',
    isomers = [
        '[O]OC(F)(F)C(OF)C(O)(F)F(2481)',
    ],
    reactants = [
        ('HF(38)', 'CF2O(48)', '[O]OC(F)(F)C=O(898)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1005',
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

