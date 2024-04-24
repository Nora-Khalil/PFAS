species(
    label = '[O]OC(F)(F)C(OF)C(O)(F)F(2294)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.75988,0.137652,-0.000236656,2.01478e-07,-6.66339e-11,-141179,38.2443], Tmin=(100,'K'), Tmax=(811.815,'K')), NASAPolynomial(coeffs=[17.7568,0.030001,-1.65235e-05,3.27343e-09,-2.28527e-13,-143969,-49.512], Tmin=(811.815,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1175.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCFFO) + radical(ROOJ)"""),
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
    label = '[O]OC(F)(F)C=O(622)',
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
        HinderedRotor(inertia=(0.589456,'amu*angstrom^2'), symmetry=1, barrier=(13.5528,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.587337,'amu*angstrom^2'), symmetry=1, barrier=(13.504,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.024,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3731.44,'J/mol'), sigma=(6.05191,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=582.84 K, Pc=38.2 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57191,0.0584796,-9.64595e-05,8.34323e-08,-2.83598e-11,-66202.3,20.9945], Tmin=(100,'K'), Tmax=(811.126,'K')), NASAPolynomial(coeffs=[8.46089,0.0179335,-9.32175e-06,1.82202e-09,-1.26843e-13,-67103.6,-9.46572], Tmin=(811.126,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-551.124,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cds-OdCsH) + radical(ROOJ)"""),
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
    label = '[O]OC(F)(F)C(O)OF(2315)',
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
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,492.5,1135,1000,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,208.606,208.901],'cm^-1')),
        HinderedRotor(inertia=(0.280497,'amu*angstrom^2'), symmetry=1, barrier=(8.739,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.771927,'amu*angstrom^2'), symmetry=1, barrier=(23.6479,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.282685,'amu*angstrom^2'), symmetry=1, barrier=(8.73862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.281621,'amu*angstrom^2'), symmetry=1, barrier=(8.74555,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (147.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0748993,0.0977111,-0.000167456,1.43617e-07,-4.79219e-11,-86340.1,30.3754], Tmin=(100,'K'), Tmax=(812.593,'K')), NASAPolynomial(coeffs=[13.1932,0.0231071,-1.25895e-05,2.48965e-09,-1.738e-13,-88189.6,-28.9944], Tmin=(812.593,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-719.028,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2sCF) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(CsCFFO) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(OF)C(O)(F)F(2358)',
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
        HarmonicOscillator(frequencies=([557,1111,492.5,1135,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,208.786,208.786],'cm^-1')),
        HinderedRotor(inertia=(0.282575,'amu*angstrom^2'), symmetry=1, barrier=(8.74106,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.764473,'amu*angstrom^2'), symmetry=1, barrier=(23.6479,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.282576,'amu*angstrom^2'), symmetry=1, barrier=(8.74106,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.282576,'amu*angstrom^2'), symmetry=1, barrier=(8.74106,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (147.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0748599,0.0977106,-0.000167454,1.43614e-07,-4.79202e-11,-86340.1,30.3753], Tmin=(100,'K'), Tmax=(812.612,'K')), NASAPolynomial(coeffs=[13.1931,0.0231072,-1.25895e-05,2.48967e-09,-1.73802e-13,-88189.6,-28.994], Tmin=(812.612,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-719.028,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(CsCFFO) + radical(ROOJ)"""),
)

species(
    label = '[O]OF(124)',
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3046.52,'J/mol'), sigma=(4.95848,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=475.86 K, Pc=56.7 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.63697,0.00939667,-1.79758e-05,1.61484e-08,-5.19288e-12,1646.96,7.91769], Tmin=(100,'K'), Tmax=(982.363,'K')), NASAPolynomial(coeffs=[4.37144,0.00223349,-6.66876e-07,7.82176e-11,-2.87106e-15,1703.99,5.41214], Tmin=(982.363,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(13.5955,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""OOF""", comment="""Thermo library: halogens"""),
)

species(
    label = 'OC(F)(F)C([C]F)OF-2(2118)',
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
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,617,898,1187,255.765,256.288],'cm^-1')),
        HinderedRotor(inertia=(0.211652,'amu*angstrom^2'), symmetry=1, barrier=(9.80798,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208983,'amu*angstrom^2'), symmetry=1, barrier=(9.80644,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.211579,'amu*angstrom^2'), symmetry=1, barrier=(9.80375,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.682456,'amu*angstrom^2'), symmetry=1, barrier=(31.7592,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.11007,0.0918323,-0.000145105,1.15039e-07,-3.5814e-11,-72158.6,29.1875], Tmin=(100,'K'), Tmax=(790.858,'K')), NASAPolynomial(coeffs=[14.0613,0.0212733,-1.12842e-05,2.23834e-09,-1.57875e-13,-74365.4,-34.8461], Tmin=(790.858,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-601.077,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CJ2_singlet-FCs)"""),
)

species(
    label = '[O]OC(F)(F)C(F)C(O)(F)OF(2336)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.89147,0.140627,-0.000241993,2.05442e-07,-6.76871e-11,-138648,37.6307], Tmin=(100,'K'), Tmax=(813.533,'K')), NASAPolynomial(coeffs=[18.3267,0.0297471,-1.64029e-05,3.24647e-09,-2.26394e-13,-141558,-53.405], Tmin=(813.533,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1154.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2sCF) + group(O2s-OsH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)(OF)C(F)C(O)(F)F(2359)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.89147,0.140627,-0.000241993,2.05442e-07,-6.76871e-11,-138648,37.6307], Tmin=(100,'K'), Tmax=(813.533,'K')), NASAPolynomial(coeffs=[18.3267,0.0297471,-1.64029e-05,3.24647e-09,-2.26394e-13,-141558,-53.405], Tmin=(813.533,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1154.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-CsH) + group(O2s-OsH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + radical(ROOJ)"""),
)

species(
    label = 'FOF(428)',
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
    label = '[O]OC(F)(F)C=C(O)F(2360)',
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
        HinderedRotor(inertia=(0.62098,'amu*angstrom^2'), symmetry=1, barrier=(14.2776,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.624206,'amu*angstrom^2'), symmetry=1, barrier=(14.3517,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.623558,'amu*angstrom^2'), symmetry=1, barrier=(14.3368,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (143.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.113666,0.0913461,-0.000143805,1.13668e-07,-3.52202e-11,-94232.1,26.9916], Tmin=(100,'K'), Tmax=(795.139,'K')), NASAPolynomial(coeffs=[14.1414,0.0207882,-1.07182e-05,2.09962e-09,-1.46948e-13,-96463.2,-37.4705], Tmin=(795.139,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-784.608,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(CsCFFO) + group(Cds-CdsCsH) + group(CdCFO) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)=CC(O)(F)F(2361)',
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
        HinderedRotor(inertia=(0.226824,'amu*angstrom^2'), symmetry=1, barrier=(5.21514,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.226268,'amu*angstrom^2'), symmetry=1, barrier=(5.20234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.226917,'amu*angstrom^2'), symmetry=1, barrier=(5.21727,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (143.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.397651,0.0910768,-0.00016728,1.56143e-07,-5.53906e-11,-87875,29.2441], Tmin=(100,'K'), Tmax=(848.099,'K')), NASAPolynomial(coeffs=[8.09771,0.0301822,-1.61081e-05,3.14049e-09,-2.16434e-13,-88297.2,-1.42451], Tmin=(848.099,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-731.618,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsCFFO) + group(Cds-CdsCsH) + group(CdCFO) + radical(ROOJ)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05764,-0.000787936,2.90878e-06,-1.47519e-09,2.12843e-13,-30281.6,-0.311364], Tmin=(100,'K'), Tmax=(1130.24,'K')), NASAPolynomial(coeffs=[2.84325,0.00275109,-7.81032e-07,1.07244e-10,-5.79392e-15,-29958.6,5.91042], Tmin=(1130.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-251.755,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""H2O""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[O]OC(F)(F)C(OF)=C(F)F(2362)',
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
        HinderedRotor(inertia=(0.168389,'amu*angstrom^2'), symmetry=1, barrier=(3.8716,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168454,'amu*angstrom^2'), symmetry=1, barrier=(3.87308,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168773,'amu*angstrom^2'), symmetry=1, barrier=(3.88042,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (179.022,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.587968,0.122067,-0.000250942,2.47238e-07,-9.01227e-11,-94131.1,32.2878], Tmin=(100,'K'), Tmax=(858.183,'K')), NASAPolynomial(coeffs=[6.66788,0.0399395,-2.29564e-05,4.53718e-09,-3.12722e-13,-93597.6,8.7565], Tmin=(858.183,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-783.855,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-OsH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(Cds-CdsCsOs) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(O)(F)F(283)',
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
        HinderedRotor(inertia=(0.625973,'amu*angstrom^2'), symmetry=1, barrier=(14.3923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.625191,'amu*angstrom^2'), symmetry=1, barrier=(14.3744,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.0136,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3791.54,'J/mol'), sigma=(6.20894,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=592.23 K, Pc=35.94 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.89633,0.0496153,-7.9902e-05,6.54995e-08,-2.10552e-11,-78757.2,18.1467], Tmin=(100,'K'), Tmax=(804.939,'K')), NASAPolynomial(coeffs=[9.02521,0.0121076,-6.12681e-06,1.18408e-09,-8.19298e-14,-79837.4,-14.2801], Tmin=(804.939,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-655.428,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(CsFFOO) + radical(ROOJ)"""),
)

species(
    label = 'FOC=C(F)F(565)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,1.49298e-14,-2.05327e-17,9.23927e-21,-1.28064e-24,29230.2,5.12616], Tmin=(100,'K'), Tmax=(3959.16,'K')), NASAPolynomial(coeffs=[2.5,1.19009e-10,-4.23987e-14,6.68965e-18,-3.94352e-22,29230.2,5.12617], Tmin=(3959.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.034,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[O]C(F)(F)C(OF)C(O)(F)F(2363)',
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
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,351,323,533,609,664,892,1120,1201,265.092,265.229,265.324],'cm^-1')),
        HinderedRotor(inertia=(0.173379,'amu*angstrom^2'), symmetry=1, barrier=(8.6613,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.173415,'amu*angstrom^2'), symmetry=1, barrier=(8.66005,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.360792,'amu*angstrom^2'), symmetry=1, barrier=(18.0296,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.173432,'amu*angstrom^2'), symmetry=1, barrier=(8.65894,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (181.038,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.00996,0.118741,-0.000193389,1.57016e-07,-4.99867e-11,-136818,34.6288], Tmin=(100,'K'), Tmax=(774.014,'K')), NASAPolynomial(coeffs=[16.8557,0.0264106,-1.44524e-05,2.89148e-09,-2.04453e-13,-139584,-46.9857], Tmin=(774.014,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1139.01,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCFFO) + radical(O2sj(Cs-CsF1sF1s))"""),
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
    label = 'OC(F)(F)C(OF)=C(F)F(2109)',
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
        HarmonicOscillator(frequencies=([3615,1277.5,1000,231,791,275,321,533,585,746,850,1103,350,440,435,1725,182,240,577,636,1210,1413,180,1146.38],'cm^-1')),
        HinderedRotor(inertia=(0.224426,'amu*angstrom^2'), symmetry=1, barrier=(5.15999,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.224206,'amu*angstrom^2'), symmetry=1, barrier=(5.15493,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.224208,'amu*angstrom^2'), symmetry=1, barrier=(5.15498,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (164.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.152126,0.108445,-0.000214076,2.08238e-07,-7.58099e-11,-121025,29.2028], Tmin=(100,'K'), Tmax=(849.061,'K')), NASAPolynomial(coeffs=[7.03291,0.0368301,-2.08388e-05,4.12619e-09,-2.85888e-13,-120884,3.73106], Tmin=(849.061,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1007.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(Cds-CdsCsOs) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs)"""),
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
    label = '[O]OC(F)(F)C(OF)[C](O)F(2364)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.38141,0.13249,-0.000243475,2.17499e-07,-7.39001e-11,-89567.2,37.4799], Tmin=(100,'K'), Tmax=(851.873,'K')), NASAPolynomial(coeffs=[14.6475,0.0304966,-1.68176e-05,3.28669e-09,-2.25708e-13,-91328.3,-31.5881], Tmin=(851.873,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-746.203,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCFHO) + radical(ROOJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]O[C](F)C(OF)C(O)(F)F(2365)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.38142,0.13249,-0.000243476,2.175e-07,-7.39007e-11,-89567.2,37.48], Tmin=(100,'K'), Tmax=(851.871,'K')), NASAPolynomial(coeffs=[14.6475,0.0304965,-1.68175e-05,3.28668e-09,-2.25707e-13,-91328.3,-31.5883], Tmin=(851.871,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-746.203,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(CsCFHO) + group(CsCFFO) + radical(ROOJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]OC(F)(F)C([O])C(O)(F)F(2366)',
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
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,1380,1390,370,380,2900,435,183,263,327,399,503,589,506,644,608,780,1145,1213,1339,1481,249.142,249.145,249.152],'cm^-1')),
        HinderedRotor(inertia=(0.199916,'amu*angstrom^2'), symmetry=1, barrier=(8.80577,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.199921,'amu*angstrom^2'), symmetry=1, barrier=(8.80586,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.537608,'amu*angstrom^2'), symmetry=1, barrier=(23.6784,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.199917,'amu*angstrom^2'), symmetry=1, barrier=(8.80595,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (178.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.09175,0.120938,-0.000202724,1.6974e-07,-5.55129e-11,-129672,35.4844], Tmin=(100,'K'), Tmax=(797.538,'K')), NASAPolynomial(coeffs=[16.4183,0.0270142,-1.45931e-05,2.8841e-09,-2.01652e-13,-132271,-43.8134], Tmin=(797.538,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1079.61,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCFFO) + radical(CC(C)OJ) + radical(ROOJ)"""),
)

species(
    label = '[O]F(228)',
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
    label = '[O]OC(F)(F)[CH]C(O)(F)F(2367)',
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
        HinderedRotor(inertia=(0.409615,'amu*angstrom^2'), symmetry=1, barrier=(9.41787,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.410448,'amu*angstrom^2'), symmetry=1, barrier=(9.437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.412258,'amu*angstrom^2'), symmetry=1, barrier=(9.47862,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.602536,'amu*angstrom^2'), symmetry=1, barrier=(13.8535,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (162.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.537156,0.109841,-0.000191766,1.6751e-07,-5.66024e-11,-112552,33.6302], Tmin=(100,'K'), Tmax=(826.613,'K')), NASAPolynomial(coeffs=[13.4297,0.027296,-1.48329e-05,2.92082e-09,-2.02991e-13,-114350,-28.001], Tmin=(826.613,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-937.086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(CsCFFO) + group(CsCFFO) + radical(ROOJ) + radical(CCJCOOH)"""),
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
    label = '[O]OC(F)(F)C(OF)[C](F)F(2368)',
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
        HinderedRotor(inertia=(0.429259,'amu*angstrom^2'), symmetry=1, barrier=(9.86952,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.430149,'amu*angstrom^2'), symmetry=1, barrier=(9.88996,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.428783,'amu*angstrom^2'), symmetry=1, barrier=(9.85857,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49595,'amu*angstrom^2'), symmetry=1, barrier=(34.3947,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (180.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.10236,0.124782,-0.000226908,2.01563e-07,-6.8585e-11,-91679.3,35.5747], Tmin=(100,'K'), Tmax=(836.514,'K')), NASAPolynomial(coeffs=[14.5898,0.0284886,-1.61195e-05,3.1949e-09,-2.21831e-13,-93560.9,-32.8842], Tmin=(836.514,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-763.692,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCsFFH) + radical(ROOJ) + radical(Csj(Cs-O2sCsH)(F1s)(F1s))"""),
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
    label = '[O]OC(F)(F)C(OF)C([O])(F)F(2369)',
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
        HinderedRotor(inertia=(0.383211,'amu*angstrom^2'), symmetry=1, barrier=(8.81078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.383008,'amu*angstrom^2'), symmetry=1, barrier=(8.8061,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.383087,'amu*angstrom^2'), symmetry=1, barrier=(8.80793,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.96686,'amu*angstrom^2'), symmetry=1, barrier=(22.23,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (196.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.55004,0.133725,-0.00023569,2.04127e-07,-6.82846e-11,-109920,38.0794], Tmin=(100,'K'), Tmax=(823.078,'K')), NASAPolynomial(coeffs=[16.7626,0.0290284,-1.62738e-05,3.23027e-09,-2.25155e-13,-112403,-43.472], Tmin=(823.078,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-915.496,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCFFO) + radical(O2sj(Cs-CsF1sF1s)) + radical(ROOJ)"""),
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
    label = 'OC(F)(F)C(OF)[C](F)F(1927)',
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
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,190,488,555,1236,1407,296.769,296.877,297.035],'cm^-1')),
        HinderedRotor(inertia=(0.161424,'amu*angstrom^2'), symmetry=1, barrier=(10.0805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161255,'amu*angstrom^2'), symmetry=1, barrier=(10.0861,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161019,'amu*angstrom^2'), symmetry=1, barrier=(10.0815,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.401912,'amu*angstrom^2'), symmetry=1, barrier=(25.1471,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (165.039,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3545.48,'J/mol'), sigma=(5.86849,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=553.80 K, Pc=39.81 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.635297,0.110767,-0.000188554,1.60499e-07,-5.33484e-11,-118574,32.3795], Tmin=(100,'K'), Tmax=(797.681,'K')), NASAPolynomial(coeffs=[14.8421,0.0255803,-1.41219e-05,2.81296e-09,-1.97451e-13,-120803,-37.2808], Tmin=(797.681,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-987.202,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCsFFH) + radical(Csj(Cs-O2sCsH)(F1s)(F1s))"""),
)

species(
    label = 'O[C](F)F(239)',
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
    label = '[O]OC(F)(F)[CH]OF(2248)',
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
        HinderedRotor(inertia=(0.728895,'amu*angstrom^2'), symmetry=1, barrier=(16.7587,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.728505,'amu*angstrom^2'), symmetry=1, barrier=(16.7498,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.72866,'amu*angstrom^2'), symmetry=1, barrier=(16.7533,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.023,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.3534,0.0868683,-0.000151314,1.27872e-07,-4.17009e-11,-40175.2,25.8819], Tmin=(100,'K'), Tmax=(824.202,'K')), NASAPolynomial(coeffs=[13.3942,0.0161224,-8.9909e-06,1.77562e-09,-1.23259e-13,-42071.6,-32.9747], Tmin=(824.202,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-335.075,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-OsH) + group(CsCFFO) + group(Cs-CsOsHH) + radical(ROOJ) + radical(CCsJO)"""),
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
    label = 'OC(F)(F)[CH]OF(568)',
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
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,253,525,597,667,842,1178,1324,3025,407.5,1350,352.5,201.033,201.254],'cm^-1')),
        HinderedRotor(inertia=(0.326313,'amu*angstrom^2'), symmetry=1, barrier=(9.56002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.602307,'amu*angstrom^2'), symmetry=1, barrier=(17.0059,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.727334,'amu*angstrom^2'), symmetry=1, barrier=(20.6568,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.995122,0.0707691,-0.000107108,7.97529e-08,-2.32354e-11,-68629,22.2258], Tmin=(100,'K'), Tmax=(844.561,'K')), NASAPolynomial(coeffs=[12.7537,0.0150761,-8.18963e-06,1.66666e-09,-1.20035e-13,-70615,-32.5155], Tmin=(844.561,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-571.479,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCsJO)"""),
)

species(
    label = '[O]OC(F)(F)[C](OF)C(O)(F)F(2370)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.80168,0.144221,-0.000273071,2.47067e-07,-8.44375e-11,-119521,38.8133], Tmin=(100,'K'), Tmax=(856.806,'K')), NASAPolynomial(coeffs=[15.2645,0.031507,-1.7899e-05,3.51404e-09,-2.41068e-13,-121233,-33.8062], Tmin=(856.806,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-995.359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCFFO) + radical(ROOJ) + radical(C2CsJO)"""),
)

species(
    label = '[O]OC(F)(F)C(OF)C(=O)F(2371)',
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
        HarmonicOscillator(frequencies=([557,1111,492.5,1135,1000,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,486,617,768,1157,1926,310.572,311.054,311.288],'cm^-1')),
        HinderedRotor(inertia=(0.115607,'amu*angstrom^2'), symmetry=1, barrier=(7.91976,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.232943,'amu*angstrom^2'), symmetry=1, barrier=(15.9676,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115306,'amu*angstrom^2'), symmetry=1, barrier=(7.91963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.456208,'amu*angstrom^2'), symmetry=1, barrier=(31.1749,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (177.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.230258,0.10026,-0.000146674,1.07138e-07,-3.09308e-11,-106827,35.1067], Tmin=(100,'K'), Tmax=(848.595,'K')), NASAPolynomial(coeffs=[15.702,0.0251603,-1.39243e-05,2.84727e-09,-2.05865e-13,-109531,-39.1413], Tmin=(848.595,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-889.423,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCFFO) + group(COCsFO) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)(F)C(=O)C(O)(F)F(2372)',
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
        HinderedRotor(inertia=(0.440281,'amu*angstrom^2'), symmetry=1, barrier=(10.1229,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.441274,'amu*angstrom^2'), symmetry=1, barrier=(10.1458,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.440912,'amu*angstrom^2'), symmetry=1, barrier=(10.1374,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.442132,'amu*angstrom^2'), symmetry=1, barrier=(10.1655,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (177.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.902338,0.119709,-0.00021629,1.90656e-07,-6.41348e-11,-144489,34.164], Tmin=(100,'K'), Tmax=(850.747,'K')), NASAPolynomial(coeffs=[14.2574,0.0268804,-1.46205e-05,2.84619e-09,-1.95208e-13,-146288,-31.9387], Tmin=(850.747,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1202.72,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cds-OdCsCs) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)(F)C(OF)=C(O)F(2373)',
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
        HinderedRotor(inertia=(0.301747,'amu*angstrom^2'), symmetry=1, barrier=(6.93775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.301156,'amu*angstrom^2'), symmetry=1, barrier=(6.92417,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.301618,'amu*angstrom^2'), symmetry=1, barrier=(6.9348,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.41164,'amu*angstrom^2'), symmetry=1, barrier=(55.4483,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (177.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.0902,0.130236,-0.000254738,2.40028e-07,-8.48051e-11,-95289.9,33.4287], Tmin=(100,'K'), Tmax=(856.135,'K')), NASAPolynomial(coeffs=[10.7894,0.0355639,-2.02413e-05,3.98847e-09,-2.74572e-13,-95888.5,-13.6549], Tmin=(856.135,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-793.663,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(Cds-CdsCsOs) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)=C(OF)C(O)(F)F(2374)',
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
        HarmonicOscillator(frequencies=([3615,1277.5,1000,231,791,492.5,1135,1000,275,321,533,585,746,850,1103,350,440,435,1725,326,540,652,719,1357,180,2374.05],'cm^-1')),
        HinderedRotor(inertia=(0.265325,'amu*angstrom^2'), symmetry=1, barrier=(6.10035,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.264866,'amu*angstrom^2'), symmetry=1, barrier=(6.08979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.264405,'amu*angstrom^2'), symmetry=1, barrier=(6.07918,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.263984,'amu*angstrom^2'), symmetry=1, barrier=(6.06952,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (177.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.710536,0.128725,-0.000273302,2.75238e-07,-1.01435e-10,-88936.9,35.3451], Tmin=(100,'K'), Tmax=(864.2,'K')), NASAPolynomial(coeffs=[4.4967,0.0454059,-2.59001e-05,5.09473e-09,-3.496e-13,-87625.6,23.7771], Tmin=(864.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-740.674,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-Cd)) + group(O2sCF) + group(O2s-OsH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(Cds-CdsCsOs) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)=COF(2182)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57388,0.0632842,-0.000121338,1.17856e-07,-4.29597e-11,-11894.4,23.0683], Tmin=(100,'K'), Tmax=(852.009,'K')), NASAPolynomial(coeffs=[5.26827,0.0238132,-1.28928e-05,2.5205e-09,-1.73717e-13,-11720.8,10.5497], Tmin=(852.009,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-99.5423,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2sCF) + group(O2s-OsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + radical(ROOJ)"""),
)

species(
    label = 'OOC(F)(F)[C](OF)C(O)(F)F(2375)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.04952,0.147632,-0.000267818,2.37296e-07,-8.06616e-11,-137792,38.9438], Tmin=(100,'K'), Tmax=(832.764,'K')), NASAPolynomial(coeffs=[16.8188,0.0332273,-1.89253e-05,3.76321e-09,-2.61905e-13,-140110,-43.6821], Tmin=(832.764,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1147.36,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCFFO) + radical(C2CsJO)"""),
)

species(
    label = '[O]C(F)(F)C(OF)C(F)(F)OO(2376)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.5788,0.134138,-0.000217653,1.73578e-07,-5.32052e-11,-128200,37.4483], Tmin=(100,'K'), Tmax=(680.474,'K')), NASAPolynomial(coeffs=[17.9244,0.0314703,-1.77399e-05,3.58758e-09,-2.55231e-13,-131132,-51.173], Tmin=(680.474,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1067.5,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCFFO) + radical(O2sj(Cs-CsF1sF1s))"""),
)

species(
    label = '[O]C(C(O)(F)F)C(F)(F)OOF(2377)',
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
        HinderedRotor(inertia=(1.06925,'amu*angstrom^2'), symmetry=1, barrier=(24.5842,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0693,'amu*angstrom^2'), symmetry=1, barrier=(24.5853,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06941,'amu*angstrom^2'), symmetry=1, barrier=(24.5879,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06936,'amu*angstrom^2'), symmetry=1, barrier=(24.5866,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06938,'amu*angstrom^2'), symmetry=1, barrier=(24.5873,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (197.037,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.54672,0.132005,-0.000212111,1.69928e-07,-5.35374e-11,-136477,37.753], Tmin=(100,'K'), Tmax=(780.842,'K')), NASAPolynomial(coeffs=[18.2681,0.0305065,-1.71421e-05,3.47672e-09,-2.4831e-13,-139571,-52.9413], Tmin=(780.842,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1136.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2sFO) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCFFO) + radical(CC(C)OJ)"""),
)

species(
    label = 'O[C](F)C(OF)C(F)(F)OOF(2378)',
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
        HinderedRotor(inertia=(1.32971,'amu*angstrom^2'), symmetry=1, barrier=(30.5728,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.868294,'amu*angstrom^2'), symmetry=1, barrier=(19.9638,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.32955,'amu*angstrom^2'), symmetry=1, barrier=(30.569,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.86834,'amu*angstrom^2'), symmetry=1, barrier=(19.9649,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.868312,'amu*angstrom^2'), symmetry=1, barrier=(19.9642,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.868328,'amu*angstrom^2'), symmetry=1, barrier=(19.9646,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (197.037,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.98839,0.145565,-0.000261002,2.30073e-07,-7.81479e-11,-96365,40.2805], Tmin=(100,'K'), Tmax=(824.843,'K')), NASAPolynomial(coeffs=[16.8469,0.0333522,-1.89814e-05,3.78509e-09,-2.64346e-13,-98762.2,-42.6577], Tmin=(824.843,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-802.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-CsH) + group(O2sFO) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCFHO) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'OC(F)(F)C(OF)[C](F)OOF(2379)',
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
        HinderedRotor(inertia=(1.32984,'amu*angstrom^2'), symmetry=1, barrier=(30.5757,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.868309,'amu*angstrom^2'), symmetry=1, barrier=(19.9641,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.868295,'amu*angstrom^2'), symmetry=1, barrier=(19.9638,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.32944,'amu*angstrom^2'), symmetry=1, barrier=(30.5665,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.868331,'amu*angstrom^2'), symmetry=1, barrier=(19.9646,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.868326,'amu*angstrom^2'), symmetry=1, barrier=(19.9645,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (197.037,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.98834,0.145565,-0.000260999,2.30069e-07,-7.81461e-11,-96365,40.2804], Tmin=(100,'K'), Tmax=(824.853,'K')), NASAPolynomial(coeffs=[16.8468,0.0333524,-1.89815e-05,3.78512e-09,-2.64348e-13,-98762.2,-42.6573], Tmin=(824.853,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-802.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-CsH) + group(O2sFO) + group(Cs-CsCsOsH) + group(CsCFHO) + group(CsCFFO) + radical(CsCsF1sO2s)"""),
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
    reactants = ['[O]OC(F)(F)C(OF)C(O)(F)F(2294)'],
    products = ['HF(38)', 'CF2O(48)', '[O]OC(F)(F)C=O(622)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(0.00285451,'s^-1'), n=4.62568, Ea=(37.0509,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CF2(42)', '[O]OC(F)(F)C(O)OF(2315)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2294)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.40908e-07,'m^3/(mol*s)'), n=3.45227, Ea=(205.694,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s',), comment="""Estimated from node CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CF2(42)', '[O]OC(OF)C(O)(F)F(2358)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2294)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.40908e-07,'m^3/(mol*s)'), n=3.45227, Ea=(205.694,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s',), comment="""Estimated from node CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]OF(124)', 'OC(F)(F)C([C]F)OF-2(2118)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2294)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.13916e+09,'m^3/(mol*s)'), n=-1.00842, Ea=(21.2535,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17406617270594374, var=16.167420070870673, Tref=1000.0, N=4, data_mean=0.0, correlation='OHOY',), comment="""Estimated from node OHOY"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]OC(F)(F)C(F)C(O)(F)OF(2336)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2294)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]OC(F)(OF)C(F)C(O)(F)F(2359)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2294)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction7',
    reactants = ['FOF(428)', '[O]OC(F)(F)C=C(O)F(2360)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2294)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.00057209,'m^3/(mol*s)'), n=2.81783, Ea=(231.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OH] for rate rule [Cd/disub_Cd/monosub;H_OH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction8',
    reactants = ['FOF(428)', '[O]OC(F)=CC(O)(F)F(2361)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2294)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.00057209,'m^3/(mol*s)'), n=2.81783, Ea=(231.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OH] for rate rule [Cd/disub_Cd/monosub;H_OH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H2O(3)', '[O]OC(F)(F)C(OF)=C(F)F(2362)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2294)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.00057209,'m^3/(mol*s)'), n=2.81783, Ea=(231.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OH] for rate rule [Cd/disub_Cd/disub;H_OH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]OC(O)(F)F(283)', 'FOC=C(F)F(565)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2294)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.79e-11,'m^3/(mol*s)'), n=3.97, Ea=(329.281,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;R_OH] for rate rule [Cd/monosub_Cd/disub;R_OH]
Euclidian distance = 1.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction11',
    reactants = ['O(6)', '[O]C(F)(F)C(OF)C(O)(F)F(2363)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2294)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]OC(F)(F)C(OF)C(O)(F)F(2294)'],
    products = ['HO2(10)', 'OC(F)(F)C(OF)=C(F)F(2109)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.00545e+13,'s^-1'), n=-0.00666667, Ea=(174.79,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO_NdNd]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', '[O]OC(F)(F)C(OF)[C](O)F(2364)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2294)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F(37)', '[O]O[C](F)C(OF)C(O)(F)F(2365)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2294)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F(37)', '[O]OC(F)(F)C([O])C(O)(F)F(2366)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2294)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.49151e+08,'m^3/(mol*s)'), n=0.0472428, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]F(228)', '[O]OC(F)(F)[CH]C(O)(F)F(2367)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2294)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['OH(4)', '[O]OC(F)(F)C(OF)[C](F)F(2368)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2294)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(5)', '[O]OC(F)(F)C(OF)C([O])(F)F(2369)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2294)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction19',
    reactants = ['O2(2)', 'OC(F)(F)C(OF)[C](F)F(1927)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2294)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.6844e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction20',
    reactants = ['O[C](F)F(239)', '[O]OC(F)(F)[CH]OF(2248)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2294)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]O[C](F)F(171)', 'OC(F)(F)[CH]OF(568)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2294)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(5)', '[O]OC(F)(F)[C](OF)C(O)(F)F(2370)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2294)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(43950,'m^3/(mol*s)'), n=1, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction23',
    reactants = ['HF(38)', '[O]OC(F)(F)C(OF)C(=O)F(2371)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2294)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(0.109156,'m^3/(mol*s)'), n=1.86531, Ea=(167.161,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd',), comment="""Estimated from node HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd"""),
)

reaction(
    label = 'reaction24',
    reactants = ['HF(38)', '[O]OC(F)(F)C(=O)C(O)(F)F(2372)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2294)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(443.781,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction25',
    reactants = ['HF(38)', '[O]OC(F)(F)C(OF)=C(O)F(2373)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2294)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.14111,'m^3/(mol*s)'), n=1.29695, Ea=(179.265,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction26',
    reactants = ['HF(38)', '[O]OC(F)=C(OF)C(O)(F)F(2374)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2294)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.14111,'m^3/(mol*s)'), n=1.29695, Ea=(157.543,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O]OC(F)(F)C(OF)C(O)(F)F(2294)'],
    products = ['HF(38)', 'CF2O(48)', '[O]OC(F)=COF(2182)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.97084e+10,'s^-1'), n=0.978818, Ea=(284.118,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R!H->C_N-5Br1sCl1sF1sH->H_5Br1sCl1sF1s->F1s_Ext-2C-R_N-7R!H->O_Sp-7BrCClFINPSSi-2C_Ext-2C-R',), comment="""Estimated from node Root_1R!H->C_N-5Br1sCl1sF1sH->H_5Br1sCl1sF1s->F1s_Ext-2C-R_N-7R!H->O_Sp-7BrCClFINPSSi-2C_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[O]OC(F)(F)C(OF)C(O)(F)F(2294)'],
    products = ['OOC(F)(F)[C](OF)C(O)(F)F(2375)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(5.29e+09,'s^-1'), n=0.75, Ea=(103.847,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 272 used for R4H_SSS_OCs;O_rad_out;Cs_H_out_NDMustO
Exact match found for rate rule [R4H_SSS_OCs;O_rad_out;Cs_H_out_NDMustO]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[O]C(F)(F)C(OF)C(F)(F)OO(2376)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2294)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(123903,'s^-1'), n=1.46258, Ea=(69.5427,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_SSSSS;O_rad_out;XH_out] for rate rule [R6H_SSSSS;O_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[O]C(C(O)(F)F)C(F)(F)OOF(2377)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2294)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(179.343,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction31',
    reactants = ['O[C](F)C(OF)C(F)(F)OOF(2378)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2294)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(61.4745,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction32',
    reactants = ['OC(F)(F)C(OF)[C](F)OOF(2379)'],
    products = ['[O]OC(F)(F)C(OF)C(O)(F)F(2294)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(83.7872,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

network(
    label = 'PDepNetwork #906',
    isomers = [
        '[O]OC(F)(F)C(OF)C(O)(F)F(2294)',
    ],
    reactants = [
        ('HF(38)', 'CF2O(48)', '[O]OC(F)(F)C=O(622)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #906',
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

