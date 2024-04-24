species(
    label = '[O]OC(F)(F)C(O)C(F)(F)OF(2483)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {11,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {7,S}
6  O u0 p2 c0 {10,S} {14,S}
7  O u0 p2 c0 {5,S} {11,S}
8  O u0 p2 c0 {9,S} {12,S}
9  O u1 p2 c0 {8,S}
10 C u0 p0 c0 {6,S} {11,S} {12,S} {13,S}
11 C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
12 C u0 p0 c0 {3,S} {4,S} {8,S} {10,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-1175.46,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,492.5,1135,1000,1380,1390,370,380,2900,435,183,263,327,399,503,589,506,644,608,780,1145,1213,1339,1481,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.75989,0.137652,-0.000236657,2.01479e-07,-6.66341e-11,-141179,38.2443], Tmin=(100,'K'), Tmax=(811.812,'K')), NASAPolynomial(coeffs=[17.7568,0.030001,-1.65235e-05,3.27343e-09,-2.28527e-13,-143969,-49.5121], Tmin=(811.812,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1175.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2sCF) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCFFO) + radical(ROOJ)"""),
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
    label = '[O]OC(O)C(F)(F)OF(2522)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {6,S}
4  O u0 p2 c0 {8,S} {11,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {3,S} {9,S}
7  O u1 p2 c0 {5,S}
8  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
9  C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-719.028,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,557,1111,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,209.257,212.821],'cm^-1')),
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
    label = 'OC([C]F)C(F)(F)OF-2(1753)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {6,S}
5  O u0 p2 c0 {7,S} {11,S}
6  O u0 p2 c0 {4,S} {8,S}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
9  C u0 p1 c0 {3,S} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-601.077,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,617,898,1187,256.061,256.061],'cm^-1')),
        HinderedRotor(inertia=(0.210756,'amu*angstrom^2'), symmetry=1, barrier=(9.80606,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.210757,'amu*angstrom^2'), symmetry=1, barrier=(9.80606,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.210755,'amu*angstrom^2'), symmetry=1, barrier=(9.80606,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.682582,'amu*angstrom^2'), symmetry=1, barrier=(31.7592,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.110076,0.0918322,-0.000145105,1.15039e-07,-3.58139e-11,-72158.6,29.1875], Tmin=(100,'K'), Tmax=(790.906,'K')), NASAPolynomial(coeffs=[14.0613,0.0212733,-1.12842e-05,2.23834e-09,-1.57875e-13,-74365.4,-34.8461], Tmin=(790.906,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-601.077,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CJ2_singlet-FCs)"""),
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
    label = '[O]OC(O)(F)C(F)C(F)(F)OF(2524)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {8,S}
6  O u0 p2 c0 {11,S} {14,S}
7  O u0 p2 c0 {9,S} {11,S}
8  O u0 p2 c0 {5,S} {12,S}
9  O u1 p2 c0 {7,S}
10 C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
11 C u0 p0 c0 {2,S} {6,S} {7,S} {10,S}
12 C u0 p0 c0 {3,S} {4,S} {8,S} {10,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-1154.46,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,557,1111,250,417,511,1155,1315,1456,3119,375,361,584,565,722,1474,223,363,546,575,694,1179,1410,200,800,1600],'cm^-1')),
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
    label = '[O]OC(F)(F)C=C(F)OF(2525)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {6,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {4,S} {10,S}
7  O u1 p2 c0 {5,S}
8  C u0 p0 c0 {1,S} {2,S} {5,S} {9,S}
9  C u0 p0 c0 {8,S} {10,D} {11,S}
10 C u0 p0 c0 {3,S} {6,S} {9,D}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-626.802,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,231,791,275,321,533,585,746,850,1103,3010,987.5,1337.5,450,1655,326,540,652,719,1357,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.280064,'amu*angstrom^2'), symmetry=1, barrier=(6.43922,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.280019,'amu*angstrom^2'), symmetry=1, barrier=(6.43818,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.280004,'amu*angstrom^2'), symmetry=1, barrier=(6.43785,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (161.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.15137,0.10489,-0.000196404,1.83418e-07,-6.50326e-11,-75250.5,30.3997], Tmin=(100,'K'), Tmax=(843.909,'K')), NASAPolynomial(coeffs=[9.39642,0.0320088,-1.77583e-05,3.50175e-09,-2.42577e-13,-75878.2,-8.21401], Tmin=(843.909,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-626.802,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-OsH) + group(CsCFFO) + group(Cds-CdsCsH) + group(CdCFO) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)=CC(F)(F)OF(2526)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u0 p2 c0 {7,S} {10,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {1,S} {2,S} {5,S} {9,S}
9  C u0 p0 c0 {8,S} {10,D} {11,S}
10 C u0 p0 c0 {3,S} {6,S} {9,D}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-597.111,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,492.5,1135,1000,275,321,533,585,746,850,1103,3010,987.5,1337.5,450,1655,326,540,652,719,1357,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.278573,'amu*angstrom^2'), symmetry=1, barrier=(6.40494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.278487,'amu*angstrom^2'), symmetry=1, barrier=(6.40297,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.279091,'amu*angstrom^2'), symmetry=1, barrier=(6.41685,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (161.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0985843,0.103684,-0.000194902,1.82328e-07,-6.47017e-11,-71681.3,31.1956], Tmin=(100,'K'), Tmax=(844.018,'K')), NASAPolynomial(coeffs=[9.31029,0.0314835,-1.75164e-05,3.45768e-09,-2.39618e-13,-72286.1,-6.77537], Tmin=(844.018,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-597.111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsCFFO) + group(Cds-CdsCsH) + group(CdCFO) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)(F)C(O)=C(F)F(2527)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {9,S} {11,S}
7  O u1 p2 c0 {5,S}
8  C u0 p0 c0 {1,S} {2,S} {5,S} {9,S}
9  C u0 p0 c0 {6,S} {8,S} {10,D}
10 C u0 p0 c0 {3,S} {4,S} {9,D}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-941.662,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,275,321,533,585,746,850,1103,350,440,435,1725,182,240,577,636,1210,1413,180],'cm^-1')),
        HinderedRotor(inertia=(0.378803,'amu*angstrom^2'), symmetry=1, barrier=(8.70943,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.376328,'amu*angstrom^2'), symmetry=1, barrier=(8.65253,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.8388,'amu*angstrom^2'), symmetry=1, barrier=(42.2777,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (161.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.413294,0.1097,-0.00020302,1.84449e-07,-6.37261e-11,-113109,29.1969], Tmin=(100,'K'), Tmax=(848.257,'K')), NASAPolynomial(coeffs=[11.6376,0.0283136,-1.56726e-05,3.07572e-09,-2.12061e-13,-114270,-21.7505], Tmin=(848.257,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-941.662,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(Cds-CdsCsOs) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + radical(ROOJ)"""),
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
    label = '[O]C(F)(F)C(O)C(F)(F)OF(2528)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {7,S}
6  O u0 p2 c0 {9,S} {13,S}
7  O u0 p2 c0 {5,S} {10,S}
8  O u1 p2 c0 {11,S}
9  C u0 p0 c0 {6,S} {10,S} {11,S} {12,S}
10 C u0 p0 c0 {1,S} {2,S} {7,S} {9,S}
11 C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-1139.01,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,351,323,533,609,664,892,1120,1201,265.267,265.268,265.268],'cm^-1')),
        HinderedRotor(inertia=(0.173432,'amu*angstrom^2'), symmetry=1, barrier=(8.6601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.173431,'amu*angstrom^2'), symmetry=1, barrier=(8.66011,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.36107,'amu*angstrom^2'), symmetry=1, barrier=(18.0295,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.173432,'amu*angstrom^2'), symmetry=1, barrier=(8.66011,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (181.038,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.00996,0.118741,-0.000193389,1.57016e-07,-4.99867e-11,-136818,34.6288], Tmin=(100,'K'), Tmax=(774.064,'K')), NASAPolynomial(coeffs=[16.8557,0.0264106,-1.44524e-05,2.89148e-09,-2.04453e-13,-139584,-46.9858], Tmin=(774.064,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1139.01,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCFFO) + radical(O2sj(Cs-CsF1sF1s))"""),
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
    label = 'OC(=C(F)F)C(F)(F)OF(1746)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {8,S}
7  O u0 p2 c0 {9,S} {11,S}
8  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
9  C u0 p0 c0 {7,S} {8,S} {10,D}
10 C u0 p0 c0 {3,S} {4,S} {9,D}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-1030.66,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,275,321,533,585,746,850,1103,350,440,435,1725,182,240,577,636,1210,1413,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.567983,'amu*angstrom^2'), symmetry=1, barrier=(13.059,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.300147,'amu*angstrom^2'), symmetry=1, barrier=(6.90097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.90287,'amu*angstrom^2'), symmetry=1, barrier=(43.7507,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (164.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.459655,0.108511,-0.000193135,1.7078e-07,-5.83616e-11,-123810,28.0137], Tmin=(100,'K'), Tmax=(820.389,'K')), NASAPolynomial(coeffs=[13.1572,0.0266084,-1.50244e-05,2.99668e-09,-2.09655e-13,-125522,-31.8026], Tmin=(820.389,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1030.66,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-Cd)H) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(Cds-CdsCsOs) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs)"""),
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
    label = '[O]OC(F)(F)C(O)[C](F)OF(2529)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {7,S}
5  O u0 p2 c0 {9,S} {13,S}
6  O u0 p2 c0 {8,S} {10,S}
7  O u0 p2 c0 {4,S} {11,S}
8  O u1 p2 c0 {6,S}
9  C u0 p0 c0 {5,S} {10,S} {11,S} {12,S}
10 C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
11 C u1 p0 c0 {3,S} {7,S} {9,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (-746.203,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,557,1111,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,395,473,707,1436,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.38128,0.132488,-0.000243469,2.17489e-07,-7.38953e-11,-89567.2,37.4795], Tmin=(100,'K'), Tmax=(851.897,'K')), NASAPolynomial(coeffs=[14.6472,0.0304971,-1.68179e-05,3.28676e-09,-2.25714e-13,-91328.2,-31.5866], Tmin=(851.897,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-746.203,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2sCF) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCFHO) + radical(ROOJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]O[C](F)C(O)C(F)(F)OF(2530)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {6,S}
5  O u0 p2 c0 {9,S} {13,S}
6  O u0 p2 c0 {4,S} {10,S}
7  O u0 p2 c0 {8,S} {11,S}
8  O u1 p2 c0 {7,S}
9  C u0 p0 c0 {5,S} {10,S} {11,S} {12,S}
10 C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
11 C u1 p0 c0 {3,S} {7,S} {9,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (-746.203,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,492.5,1135,1000,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,395,473,707,1436,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.38144,0.13249,-0.000243477,2.17501e-07,-7.39012e-11,-89567.2,37.48], Tmin=(100,'K'), Tmax=(851.868,'K')), NASAPolynomial(coeffs=[14.6476,0.0304965,-1.68175e-05,3.28667e-09,-2.25706e-13,-91328.3,-31.5885], Tmin=(851.868,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-746.203,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2sCF) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(CsCFHO) + group(CsCFFO) + radical(ROOJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]OC(F)(F)C(O)C([O])(F)F(2531)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {9,S} {13,S}
6  O u0 p2 c0 {8,S} {10,S}
7  O u1 p2 c0 {11,S}
8  O u1 p2 c0 {6,S}
9  C u0 p0 c0 {5,S} {10,S} {11,S} {12,S}
10 C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
11 C u0 p0 c0 {3,S} {4,S} {7,S} {9,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {5,S}
"""),
    E0 = (-1050,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,351,323,533,609,664,892,1120,1201,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.386202,'amu*angstrom^2'), symmetry=1, barrier=(8.87953,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.386259,'amu*angstrom^2'), symmetry=1, barrier=(8.88084,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.386056,'amu*angstrom^2'), symmetry=1, barrier=(8.8762,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.386091,'amu*angstrom^2'), symmetry=1, barrier=(8.87699,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (178.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.05418,0.121121,-0.000208074,1.77933e-07,-5.8959e-11,-126114,36.1293], Tmin=(100,'K'), Tmax=(825.245,'K')), NASAPolynomial(coeffs=[15.5556,0.0277172,-1.48597e-05,2.91166e-09,-2.01853e-13,-128416,-38.1528], Tmin=(825.245,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1050,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCFFO) + radical(O2sj(Cs-CsF1sF1s)) + radical(ROOJ)"""),
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
    label = '[O]OC(F)(F)[CH]C(F)(F)OF(2532)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {9,S}
7  O u0 p2 c0 {8,S} {10,S}
8  O u1 p2 c0 {7,S}
9  C u0 p0 c0 {1,S} {2,S} {6,S} {11,S}
10 C u0 p0 c0 {3,S} {4,S} {7,S} {11,S}
11 C u1 p0 c0 {9,S} {10,S} {12,S}
12 H u0 p0 c0 {11,S}
"""),
    E0 = (-802.578,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,492.5,1135,1000,202,304,469,581,520,674,579,755,732,952,1136,1220,1240,1408,3025,407.5,1350,352.5,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.446269,'amu*angstrom^2'), symmetry=1, barrier=(10.2606,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.446238,'amu*angstrom^2'), symmetry=1, barrier=(10.2599,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.446726,'amu*angstrom^2'), symmetry=1, barrier=(10.2711,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.51764,'amu*angstrom^2'), symmetry=1, barrier=(34.8936,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (180.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.03324,0.122448,-0.000219394,1.93721e-07,-6.59361e-11,-96357.9,35.5812], Tmin=(100,'K'), Tmax=(824.401,'K')), NASAPolynomial(coeffs=[14.6372,0.0286062,-1.62464e-05,3.23928e-09,-2.2628e-13,-98336.5,-33.3234], Tmin=(824.401,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-802.578,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-OsH) + group(Cs-CsCsHH) + group(CsCFFO) + group(CsCFFO) + radical(ROOJ) + radical(CCJCOOH)"""),
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
    label = '[O]OC(F)(F)C([O])C(F)(F)OF(2533)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {11,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {6,S}
6  O u0 p2 c0 {5,S} {11,S}
7  O u0 p2 c0 {9,S} {12,S}
8  O u1 p2 c0 {10,S}
9  O u1 p2 c0 {7,S}
10 C u0 p0 c0 {8,S} {11,S} {12,S} {13,S}
11 C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
12 C u0 p0 c0 {3,S} {4,S} {7,S} {10,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-945.103,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,492.5,1135,1000,1380,1390,370,380,2900,435,183,263,327,399,503,589,506,644,608,780,1145,1213,1339,1481,233.038,233.038,233.039,233.039],'cm^-1')),
        HinderedRotor(inertia=(0.231962,'amu*angstrom^2'), symmetry=1, barrier=(8.93921,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.231962,'amu*angstrom^2'), symmetry=1, barrier=(8.93921,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.890046,'amu*angstrom^2'), symmetry=1, barrier=(34.3,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.432105,'amu*angstrom^2'), symmetry=1, barrier=(16.6521,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (196.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.59206,0.133602,-0.000230582,1.96306e-07,-6.5029e-11,-113478,37.4501], Tmin=(100,'K'), Tmax=(800.492,'K')), NASAPolynomial(coeffs=[17.6328,0.0283117,-1.59988e-05,3.20066e-09,-2.24779e-13,-116261,-49.1746], Tmin=(800.492,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-945.103,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2sCF) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCFFO) + radical(CC(C)OJ) + radical(ROOJ)"""),
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
    label = '[O]OC(F)(F)C(O)[C](F)F(2534)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {8,S} {12,S}
6  O u0 p2 c0 {7,S} {9,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {5,S} {9,S} {10,S} {11,S}
9  C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
10 C u1 p0 c0 {3,S} {4,S} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-898.2,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,190,488,555,1236,1407,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.357072,'amu*angstrom^2'), symmetry=1, barrier=(8.20978,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.357969,'amu*angstrom^2'), symmetry=1, barrier=(8.2304,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.357494,'amu*angstrom^2'), symmetry=1, barrier=(8.2195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.800197,'amu*angstrom^2'), symmetry=1, barrier=(18.3981,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (162.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.607135,0.112187,-0.000199332,1.75439e-07,-5.92994e-11,-107873,33.6268], Tmin=(100,'K'), Tmax=(840.284,'K')), NASAPolynomial(coeffs=[13.3816,0.0271796,-1.47067e-05,2.8766e-09,-1.98554e-13,-109574,-27.5576], Tmin=(840.284,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-898.2,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCsFFH) + radical(ROOJ) + radical(Csj(Cs-O2sCsH)(F1s)(F1s))"""),
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
    label = 'OC([C](F)F)C(F)(F)OF(1670)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {7,S}
6  O u0 p2 c0 {8,S} {12,S}
7  O u0 p2 c0 {5,S} {9,S}
8  C u0 p0 c0 {6,S} {9,S} {10,S} {11,S}
9  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
10 C u1 p0 c0 {3,S} {4,S} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-987.202,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,190,488,555,1236,1407,296.864,296.886,296.891],'cm^-1')),
        HinderedRotor(inertia=(0.161147,'amu*angstrom^2'), symmetry=1, barrier=(10.0835,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161192,'amu*angstrom^2'), symmetry=1, barrier=(10.0817,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161229,'amu*angstrom^2'), symmetry=1, barrier=(10.0828,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.40214,'amu*angstrom^2'), symmetry=1, barrier=(25.147,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (165.039,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3545.48,'J/mol'), sigma=(5.86849,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=553.80 K, Pc=39.81 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.635244,0.110766,-0.000188551,1.60495e-07,-5.3346e-11,-118574,32.3793], Tmin=(100,'K'), Tmax=(797.711,'K')), NASAPolynomial(coeffs=[14.842,0.0255805,-1.4122e-05,2.81298e-09,-1.97453e-13,-120803,-37.2803], Tmin=(797.711,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-987.202,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCsFFH) + radical(Csj(Cs-O2sCsH)(F1s)(F1s))"""),
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
    label = '[O]OC(F)(F)[CH]O(1521)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {5,S} {6,S}
4 O u0 p2 c0 {7,S} {9,S}
5 O u1 p2 c0 {3,S}
6 C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
7 C u1 p0 c0 {4,S} {6,S} {8,S}
8 H u0 p0 c0 {7,S}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (-475.94,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,253,525,597,667,842,1178,1324,3025,407.5,1350,352.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.294388,'amu*angstrom^2'), symmetry=1, barrier=(6.76857,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.295162,'amu*angstrom^2'), symmetry=1, barrier=(6.78636,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.294623,'amu*angstrom^2'), symmetry=1, barrier=(6.77395,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.882327,0.0842581,-0.00017616,1.72311e-07,-6.13484e-11,-57144.9,24.0213], Tmin=(100,'K'), Tmax=(889.397,'K')), NASAPolynomial(coeffs=[5.46318,0.0257653,-1.36059e-05,2.56448e-09,-1.70124e-13,-56461.1,10.8832], Tmin=(889.397,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-475.94,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(CsCFFO) + group(Cs-CsOsHH) + radical(ROOJ) + radical(Csj(Cs-F1sF1sO2s)(O2s-H)(H))"""),
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
    label = 'O[CH]C(F)(F)OF(1218)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {4,S}
4 O u0 p2 c0 {3,S} {6,S}
5 O u0 p2 c0 {7,S} {9,S}
6 C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
7 C u1 p0 c0 {5,S} {6,S} {8,S}
8 H u0 p0 c0 {7,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-552.736,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,253,525,597,667,842,1178,1324,3025,407.5,1350,352.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.403572,'amu*angstrom^2'), symmetry=1, barrier=(9.27892,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.402425,'amu*angstrom^2'), symmetry=1, barrier=(9.25255,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.402962,'amu*angstrom^2'), symmetry=1, barrier=(9.2649,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.610852,0.0911412,-0.000191878,1.8838e-07,-6.76233e-11,-66372.5,21.9075], Tmin=(100,'K'), Tmax=(876.473,'K')), NASAPolynomial(coeffs=[6.12535,0.0268665,-1.49479e-05,2.8952e-09,-1.95849e-13,-65837,4.5996], Tmin=(876.473,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-552.736,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sF1sO2s)(O2s-H)(H))"""),
)

species(
    label = '[O]OC(F)(F)[C](O)C(F)(F)OF(2535)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {7,S}
6  O u0 p2 c0 {9,S} {10,S}
7  O u0 p2 c0 {5,S} {11,S}
8  O u0 p2 c0 {12,S} {13,S}
9  O u1 p2 c0 {6,S}
10 C u0 p0 c0 {1,S} {2,S} {6,S} {12,S}
11 C u0 p0 c0 {3,S} {4,S} {7,S} {12,S}
12 C u1 p0 c0 {8,S} {10,S} {11,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-998.836,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,557,1111,3615,1277.5,1000,202,304,469,581,520,674,579,755,732,952,1136,1220,1240,1408,360,370,350,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.96308,0.146984,-0.000275857,2.46481e-07,-8.3363e-11,-119933,38.3964], Tmin=(100,'K'), Tmax=(855.86,'K')), NASAPolynomial(coeffs=[16.6592,0.0295892,-1.68982e-05,3.31982e-09,-2.27726e-13,-122008,-42.051], Tmin=(855.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-998.836,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2sCF) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCFFO) + radical(ROOJ) + radical(C2CsJOH)"""),
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
    label = '[O]OC(F)(F)C(O)C(=O)F(2536)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {8,S} {12,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {10,D}
7  O u1 p2 c0 {5,S}
8  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
9  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
10 C u0 p0 c0 {3,S} {6,D} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {4,S}
"""),
    E0 = (-1023.93,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,486,617,768,1157,1926,255.612,257.402],'cm^-1')),
        HinderedRotor(inertia=(0.00255448,'amu*angstrom^2'), symmetry=1, barrier=(0.119871,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.207042,'amu*angstrom^2'), symmetry=1, barrier=(9.70537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.206235,'amu*angstrom^2'), symmetry=1, barrier=(9.71053,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.43411,'amu*angstrom^2'), symmetry=1, barrier=(20.4248,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (159.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.247498,0.0879293,-0.000120311,8.30254e-08,-2.27112e-11,-123020,33.2177], Tmin=(100,'K'), Tmax=(893.53,'K')), NASAPolynomial(coeffs=[14.6091,0.0236369,-1.23798e-05,2.49643e-09,-1.79796e-13,-125586,-34.4519], Tmin=(893.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1023.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCFFO) + group(COCsFO) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)(F)C(O)=C(F)OF(2537)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {7,S}
5  O u0 p2 c0 {8,S} {9,S}
6  O u0 p2 c0 {10,S} {12,S}
7  O u0 p2 c0 {4,S} {11,S}
8  O u1 p2 c0 {5,S}
9  C u0 p0 c0 {1,S} {2,S} {5,S} {10,S}
10 C u0 p0 c0 {6,S} {9,S} {11,D}
11 C u0 p0 c0 {3,S} {7,S} {10,D}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-793.663,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,231,791,275,321,533,585,746,850,1103,350,440,435,1725,326,540,652,719,1357,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.301301,'amu*angstrom^2'), symmetry=1, barrier=(6.92751,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.301611,'amu*angstrom^2'), symmetry=1, barrier=(6.93464,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.301601,'amu*angstrom^2'), symmetry=1, barrier=(6.9344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.41005,'amu*angstrom^2'), symmetry=1, barrier=(55.4118,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (177.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.0905,0.13024,-0.000254754,2.40051e-07,-8.48162e-11,-95289.9,33.4297], Tmin=(100,'K'), Tmax=(856.097,'K')), NASAPolynomial(coeffs=[10.79,0.0355628,-2.02406e-05,3.9883e-09,-2.74557e-13,-95888.8,-13.6585], Tmin=(856.097,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-793.663,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2sCF) + group(O2s-OsH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(Cds-CdsCsOs) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)=C(O)C(F)(F)OF(2538)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u0 p2 c0 {10,S} {12,S}
7  O u0 p2 c0 {8,S} {11,S}
8  O u1 p2 c0 {7,S}
9  C u0 p0 c0 {1,S} {2,S} {5,S} {10,S}
10 C u0 p0 c0 {6,S} {9,S} {11,D}
11 C u0 p0 c0 {3,S} {7,S} {10,D}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-763.972,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,492.5,1135,1000,275,321,533,585,746,850,1103,350,440,435,1725,326,540,652,719,1357,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.3735,'amu*angstrom^2'), symmetry=1, barrier=(8.5875,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.373161,'amu*angstrom^2'), symmetry=1, barrier=(8.5797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.373207,'amu*angstrom^2'), symmetry=1, barrier=(8.58076,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.373298,'amu*angstrom^2'), symmetry=1, barrier=(8.58286,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (177.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.03726,0.129028,-0.000253228,2.38927e-07,-8.44693e-11,-91720.7,34.224], Tmin=(100,'K'), Tmax=(856.243,'K')), NASAPolynomial(coeffs=[10.7028,0.0350395,-1.99998e-05,3.94452e-09,-2.71623e-13,-92296.2,-12.2136], Tmin=(856.243,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-763.972,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(Cds-CdsCsOs) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)=CO(2539)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {5,S} {8,S}
3 O u0 p2 c0 {4,S} {6,S}
4 O u1 p2 c0 {3,S}
5 C u0 p0 c0 {2,S} {6,D} {7,S}
6 C u0 p0 c0 {1,S} {3,S} {5,D}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {2,S}
"""),
    E0 = (-257.348,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (93.0338,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.82782,0.0498996,-6.94357e-05,4.92149e-08,-1.37138e-11,-30875.5,19.6982], Tmin=(100,'K'), Tmax=(882.076,'K')), NASAPolynomial(coeffs=[10.0829,0.0124643,-5.77428e-06,1.09905e-09,-7.64361e-14,-32331.8,-19.0916], Tmin=(882.076,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-257.348,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + radical(ROOJ)"""),
)

species(
    label = 'OOC(F)(F)[C](O)C(F)(F)OF(2540)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {11,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {7,S}
6  O u0 p2 c0 {9,S} {10,S}
7  O u0 p2 c0 {5,S} {11,S}
8  O u0 p2 c0 {12,S} {13,S}
9  O u0 p2 c0 {6,S} {14,S}
10 C u0 p0 c0 {3,S} {4,S} {6,S} {12,S}
11 C u0 p0 c0 {1,S} {2,S} {7,S} {12,S}
12 C u1 p0 c0 {8,S} {10,S} {11,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-1150.84,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (197.037,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.20738,0.150349,-0.000270433,2.36471e-07,-7.94772e-11,-138204,38.5143], Tmin=(100,'K'), Tmax=(829.321,'K')), NASAPolynomial(coeffs=[18.2018,0.0313304,-1.7937e-05,3.57201e-09,-2.48817e-13,-140881,-51.8614], Tmin=(829.321,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1150.84,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2sCF) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCFFO) + radical(C2CsJOH)"""),
)

species(
    label = '[O]C(C(F)(F)OO)C(F)(F)OF(2541)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {12,S}
2  F u0 p3 c0 {12,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {7,S}
6  O u0 p2 c0 {8,S} {11,S}
7  O u0 p2 c0 {5,S} {12,S}
8  O u0 p2 c0 {6,S} {14,S}
9  O u1 p2 c0 {10,S}
10 C u0 p0 c0 {9,S} {11,S} {12,S} {13,S}
11 C u0 p0 c0 {3,S} {4,S} {6,S} {10,S}
12 C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-1097.11,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,557,1111,1380,1390,370,380,2900,435,183,263,327,399,503,589,506,644,608,780,1145,1213,1339,1481,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.64519,0.134447,-0.000214965,1.70816e-07,-5.33799e-11,-131758,36.8987], Tmin=(100,'K'), Tmax=(787.063,'K')), NASAPolynomial(coeffs=[18.74,0.0308442,-1.7516e-05,3.56978e-09,-2.55819e-13,-134967,-56.5669], Tmin=(787.063,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1097.11,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2sCF) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCFFO) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]C(F)(F)C(O)C(F)(F)OOF(2542)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {11,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {8,S}
6  O u0 p2 c0 {8,S} {11,S}
7  O u0 p2 c0 {10,S} {14,S}
8  O u0 p2 c0 {5,S} {6,S}
9  O u1 p2 c0 {12,S}
10 C u0 p0 c0 {7,S} {11,S} {12,S} {13,S}
11 C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
12 C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (-1106.71,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3615,1277.5,1000,277,555,632,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,351,323,533,609,664,892,1120,1201],'cm^-1')),
        HinderedRotor(inertia=(0.885285,'amu*angstrom^2'), symmetry=1, barrier=(20.3544,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.891138,'amu*angstrom^2'), symmetry=1, barrier=(20.489,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.889441,'amu*angstrom^2'), symmetry=1, barrier=(20.45,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.891407,'amu*angstrom^2'), symmetry=1, barrier=(20.4952,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (197.037,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.37836,0.130224,-0.000208013,1.60564e-07,-4.60358e-11,-132924,37.952], Tmin=(100,'K'), Tmax=(653.549,'K')), NASAPolynomial(coeffs=[17.4033,0.0312283,-1.74266e-05,3.50977e-09,-2.4905e-13,-135719,-47.2773], Tmin=(653.549,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1106.71,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2sFO) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCFFO) + radical(O2sj(Cs-CsF1sF1s))"""),
)

species(
    label = 'OC([C](F)OF)C(F)(F)OOF(2543)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {11,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {9,S}
6  O u0 p2 c0 {9,S} {11,S}
7  O u0 p2 c0 {10,S} {14,S}
8  O u0 p2 c0 {4,S} {12,S}
9  O u0 p2 c0 {5,S} {6,S}
10 C u0 p0 c0 {7,S} {11,S} {12,S} {13,S}
11 C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
12 C u1 p0 c0 {3,S} {8,S} {10,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (-802.907,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3615,1277.5,1000,557,1111,277,555,632,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,395,473,707,1436],'cm^-1')),
        HinderedRotor(inertia=(0.869328,'amu*angstrom^2'), symmetry=1, barrier=(19.9876,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.867979,'amu*angstrom^2'), symmetry=1, barrier=(19.9565,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.868202,'amu*angstrom^2'), symmetry=1, barrier=(19.9617,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.867866,'amu*angstrom^2'), symmetry=1, barrier=(19.9539,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.33234,'amu*angstrom^2'), symmetry=1, barrier=(30.633,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.32677,'amu*angstrom^2'), symmetry=1, barrier=(30.505,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (197.037,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.98856,0.145567,-0.000261011,2.30087e-07,-7.81551e-11,-96365,40.2811], Tmin=(100,'K'), Tmax=(824.805,'K')), NASAPolynomial(coeffs=[16.8472,0.0333517,-1.8981e-05,3.78501e-09,-2.64339e-13,-98762.3,-42.6596], Tmin=(824.805,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-802.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2sCF) + group(O2sFO) + group(Cs-CsCsOsH) + group(CsCFFO) + group(CsCFHO) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'OC([C](F)OOF)C(F)(F)OF(2544)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {11,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {7,S}
5  F u0 p3 c0 {9,S}
6  O u0 p2 c0 {10,S} {14,S}
7  O u0 p2 c0 {4,S} {11,S}
8  O u0 p2 c0 {9,S} {12,S}
9  O u0 p2 c0 {5,S} {8,S}
10 C u0 p0 c0 {6,S} {11,S} {12,S} {13,S}
11 C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
12 C u1 p0 c0 {3,S} {8,S} {10,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {6,S}
"""),
    E0 = (-802.907,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,3615,1310,387.5,850,1000,277,555,632,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,395,473,707,1436],'cm^-1')),
        HinderedRotor(inertia=(0.919892,'amu*angstrom^2'), symmetry=1, barrier=(21.1501,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.918493,'amu*angstrom^2'), symmetry=1, barrier=(21.118,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.917904,'amu*angstrom^2'), symmetry=1, barrier=(21.1044,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.920163,'amu*angstrom^2'), symmetry=1, barrier=(21.1564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.919541,'amu*angstrom^2'), symmetry=1, barrier=(21.1421,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.55063,'amu*angstrom^2'), symmetry=1, barrier=(35.652,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (197.037,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.98856,0.145567,-0.000261011,2.30087e-07,-7.81551e-11,-96365,40.2811], Tmin=(100,'K'), Tmax=(824.805,'K')), NASAPolynomial(coeffs=[16.8472,0.0333517,-1.8981e-05,3.78501e-09,-2.64339e-13,-98762.3,-42.6596], Tmin=(824.805,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-802.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2sCF) + group(O2sFO) + group(Cs-CsCsOsH) + group(CsCFHO) + group(CsCFFO) + radical(CsCsF1sO2s)"""),
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
    E0 = (-606.359,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-184.992,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-184.992,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-34.1745,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-323.407,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-323.407,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (41.7799,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (71.4714,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-273.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-363.918,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-491.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-141.257,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-141.257,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-445.058,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-242.13,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-201.244,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-263.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-463.775,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-261.403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-216.323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-254.977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-418.546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-363.458,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-346.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-316.469,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-539.563,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-506.271,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (-408.881,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (-209.378,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (-187.066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]OC(F)(F)C(O)C(F)(F)OF(2483)'],
    products = ['HF(38)', 'CF2O(48)', '[O]OC(F)(F)C=O(898)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(0.00285451,'s^-1'), n=4.62568, Ea=(37.0509,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CF2(87)', '[O]OC(F)(F)C(O)OF(2502)'],
    products = ['[O]OC(F)(F)C(O)C(F)(F)OF(2483)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.40908e-07,'m^3/(mol*s)'), n=3.45227, Ea=(205.694,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s',), comment="""Estimated from node CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CF2(87)', '[O]OC(O)C(F)(F)OF(2522)'],
    products = ['[O]OC(F)(F)C(O)C(F)(F)OF(2483)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.40908e-07,'m^3/(mol*s)'), n=3.45227, Ea=(205.694,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s',), comment="""Estimated from node CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]OF(152)', 'OC([C]F)C(F)(F)OF-2(1753)'],
    products = ['[O]OC(F)(F)C(O)C(F)(F)OF(2483)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.13916e+09,'m^3/(mol*s)'), n=-1.00842, Ea=(21.2535,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17406617270594374, var=16.167420070870673, Tref=1000.0, N=4, data_mean=0.0, correlation='OHOY',), comment="""Estimated from node OHOY"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]OC(F)(F)C(F)C(O)(F)OF(2523)'],
    products = ['[O]OC(F)(F)C(O)C(F)(F)OF(2483)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]OC(O)(F)C(F)C(F)(F)OF(2524)'],
    products = ['[O]OC(F)(F)C(O)C(F)(F)OF(2483)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction7',
    reactants = ['OF(153)', '[O]OC(F)(F)C=C(F)OF(2525)'],
    products = ['[O]OC(F)(F)C(O)C(F)(F)OF(2483)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.000286045,'m^3/(mol*s)'), n=2.81783, Ea=(231.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OH] for rate rule [Cd/disub_Cd/monosub;H_OH]
Euclidian distance = 1.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction8',
    reactants = ['OF(153)', '[O]OC(F)=CC(F)(F)OF(2526)'],
    products = ['[O]OC(F)(F)C(O)C(F)(F)OF(2483)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.000286045,'m^3/(mol*s)'), n=2.81783, Ea=(231.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OH] for rate rule [Cd/disub_Cd/monosub;H_OH]
Euclidian distance = 1.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction9',
    reactants = ['OF(153)', '[O]OC(F)(F)C(O)=C(F)F(2527)'],
    products = ['[O]OC(F)(F)C(O)C(F)(F)OF(2483)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.000286045,'m^3/(mol*s)'), n=2.81783, Ea=(231.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OH] for rate rule [Cd/disub_Cd/disub;H_OH]
Euclidian distance = 1.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction10',
    reactants = ['O(6)', '[O]C(F)(F)C(O)C(F)(F)OF(2528)'],
    products = ['[O]OC(F)(F)C(O)C(F)(F)OF(2483)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]OC(F)(F)C(O)C(F)(F)OF(2483)'],
    products = ['HO2(10)', 'OC(=C(F)F)C(F)(F)OF(1746)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.00545e+13,'s^-1'), n=-0.00666667, Ea=(151.89,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO_NdNd]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', '[O]OC(F)(F)C(O)[C](F)OF(2529)'],
    products = ['[O]OC(F)(F)C(O)C(F)(F)OF(2483)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', '[O]O[C](F)C(O)C(F)(F)OF(2530)'],
    products = ['[O]OC(F)(F)C(O)C(F)(F)OF(2483)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F(37)', '[O]OC(F)(F)C(O)C([O])(F)F(2531)'],
    products = ['[O]OC(F)(F)C(O)C(F)(F)OF(2483)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.49151e+08,'m^3/(mol*s)'), n=0.0472428, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C"""),
)

reaction(
    label = 'reaction15',
    reactants = ['OH(4)', '[O]OC(F)(F)[CH]C(F)(F)OF(2532)'],
    products = ['[O]OC(F)(F)C(O)C(F)(F)OF(2483)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(5)', '[O]OC(F)(F)C([O])C(F)(F)OF(2533)'],
    products = ['[O]OC(F)(F)C(O)C(F)(F)OF(2483)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]F(154)', '[O]OC(F)(F)C(O)[C](F)F(2534)'],
    products = ['[O]OC(F)(F)C(O)C(F)(F)OF(2483)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['O2(2)', 'OC([C](F)F)C(F)(F)OF(1670)'],
    products = ['[O]OC(F)(F)C(O)C(F)(F)OF(2483)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7.6844e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction19',
    reactants = ['FO[C](F)F(192)', '[O]OC(F)(F)[CH]O(1521)'],
    products = ['[O]OC(F)(F)C(O)C(F)(F)OF(2483)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O]O[C](F)F(171)', 'O[CH]C(F)(F)OF(1218)'],
    products = ['[O]OC(F)(F)C(O)C(F)(F)OF(2483)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(5)', '[O]OC(F)(F)[C](O)C(F)(F)OF(2535)'],
    products = ['[O]OC(F)(F)C(O)C(F)(F)OF(2483)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(43950,'m^3/(mol*s)'), n=1, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction22',
    reactants = ['F2(106)', '[O]OC(F)(F)C(O)C(=O)F(2536)'],
    products = ['[O]OC(F)(F)C(O)C(F)(F)OF(2483)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(82.1367,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction23',
    reactants = ['HF(38)', '[O]OC(F)(F)C(O)=C(F)OF(2537)'],
    products = ['[O]OC(F)(F)C(O)C(F)(F)OF(2483)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.14111,'m^3/(mol*s)'), n=1.29695, Ea=(179.265,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction24',
    reactants = ['HF(38)', '[O]OC(F)=C(O)C(F)(F)OF(2538)'],
    products = ['[O]OC(F)(F)C(O)C(F)(F)OF(2483)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.14111,'m^3/(mol*s)'), n=1.29695, Ea=(166.831,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O]OC(F)(F)C(O)C(F)(F)OF(2483)'],
    products = ['F2(106)', 'CF2O(48)', '[O]OC(F)=CO(2539)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.00570902,'s^-1'), n=4.62568, Ea=(326.941,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O]OC(F)(F)C(O)C(F)(F)OF(2483)'],
    products = ['OOC(F)(F)[C](O)C(F)(F)OF(2540)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(5.29e+09,'s^-1'), n=0.75, Ea=(103.847,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 272 used for R4H_SSS_OCs;O_rad_out;Cs_H_out_NDMustO
Exact match found for rate rule [R4H_SSS_OCs;O_rad_out;Cs_H_out_NDMustO]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O]C(C(F)(F)OO)C(F)(F)OF(2541)'],
    products = ['[O]OC(F)(F)C(O)C(F)(F)OF(2483)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.37227e+06,'s^-1'), n=1.56745, Ea=(58.7826,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;O_rad_out;XH_out] for rate rule [R5H_SSSS;O_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[O]C(F)(F)C(O)C(F)(F)OOF(2542)'],
    products = ['[O]OC(F)(F)C(O)C(F)(F)OF(2483)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(165.773,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction29',
    reactants = ['OC([C](F)OF)C(F)(F)OOF(2543)'],
    products = ['[O]OC(F)(F)C(O)C(F)(F)OF(2483)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(61.4745,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction30',
    reactants = ['OC([C](F)OOF)C(F)(F)OF(2544)'],
    products = ['[O]OC(F)(F)C(O)C(F)(F)OF(2483)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(83.7872,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

network(
    label = 'PDepNetwork #1007',
    isomers = [
        '[O]OC(F)(F)C(O)C(F)(F)OF(2483)',
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
    label = 'PDepNetwork #1007',
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

