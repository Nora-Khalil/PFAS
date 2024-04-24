species(
    label = '[O]OC(F)(F)C(F)OC(O)(F)F(2298)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {10,S} {12,S}
7  O u0 p2 c0 {9,S} {11,S}
8  O u0 p2 c0 {12,S} {14,S}
9  O u1 p2 c0 {7,S}
10 C u0 p0 c0 {1,S} {6,S} {11,S} {13,S}
11 C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
12 C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-1498.29,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,261,493,600,1152,1365,1422,3097,223,363,546,575,694,1179,1410,435,565,619,662,854,1178,1396,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.25402,0.124997,-0.000196786,1.56624e-07,-4.92369e-11,-180022,36.3531], Tmin=(100,'K'), Tmax=(781.552,'K')), NASAPolynomial(coeffs=[17.0233,0.0314524,-1.72471e-05,3.47571e-09,-2.4773e-13,-182879,-47.3194], Tmin=(781.552,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1498.29,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(CsCFHO) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsFFOO) + radical(ROOJ)"""),
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
    label = 'CHF(40)',
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
    label = '[O]OC(F)(F)OC(O)(F)F(2380)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {9,S} {10,S}
6  O u0 p2 c0 {9,S} {11,S}
7  O u0 p2 c0 {8,S} {10,S}
8  O u1 p2 c0 {7,S}
9  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
10 C u0 p0 c0 {3,S} {4,S} {5,S} {7,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-1294.6,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,408,462,496,634,573,665,597,727,750,958,1142,1214,1310,1482,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.489275,'amu*angstrom^2'), symmetry=1, barrier=(11.2494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.490682,'amu*angstrom^2'), symmetry=1, barrier=(11.2817,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.489779,'amu*angstrom^2'), symmetry=1, barrier=(11.261,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.95371,'amu*angstrom^2'), symmetry=1, barrier=(44.9197,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (165.02,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.234967,0.101107,-0.000170906,1.45041e-07,-4.81806e-11,-155559,30.2196], Tmin=(100,'K'), Tmax=(792.179,'K')), NASAPolynomial(coeffs=[13.8601,0.023868,-1.3164e-05,2.62281e-09,-1.84278e-13,-157602,-33.2958], Tmin=(792.179,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1294.6,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(CsFFOO) + group(CsFFOO) + radical(ROOJ)"""),
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
    label = '[O]OC(F)OC(O)(F)F(2381)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {8,S} {9,S}
5  O u0 p2 c0 {8,S} {11,S}
6  O u0 p2 c0 {7,S} {9,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
9  C u0 p0 c0 {3,S} {4,S} {6,S} {10,S}
10 H u0 p0 c0 {9,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-1062.47,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,435,565,619,662,854,1178,1396,509,613,660,1171,1360,1414,3084,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.492846,'amu*angstrom^2'), symmetry=1, barrier=(11.3315,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.492828,'amu*angstrom^2'), symmetry=1, barrier=(11.3311,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.492794,'amu*angstrom^2'), symmetry=1, barrier=(11.3303,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.51385,'amu*angstrom^2'), symmetry=1, barrier=(34.8064,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (147.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0634882,0.0938564,-0.000156469,1.32259e-07,-4.37632e-11,-127651,28.4909], Tmin=(100,'K'), Tmax=(803.475,'K')), NASAPolynomial(coeffs=[12.8965,0.0233291,-1.24064e-05,2.44041e-09,-1.70291e-13,-129499,-29.2792], Tmin=(803.475,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1062.47,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(CsFHOO) + group(CsFFOO) + radical(ROOJ)"""),
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
    label = 'OC(F)(F)OC(F)[C]F-2(2167)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {7,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
8  C u0 p0 c0 {3,S} {5,S} {9,S} {10,S}
9  C u0 p1 c0 {4,S} {8,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-943.796,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,435,565,619,662,854,1178,1396,2750,2850,1437.5,1250,1305,750,350,617,898,1187,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.62462,'amu*angstrom^2'), symmetry=1, barrier=(37.3533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.86002,'amu*angstrom^2'), symmetry=1, barrier=(19.7736,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.860433,'amu*angstrom^2'), symmetry=1, barrier=(19.783,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.860529,'amu*angstrom^2'), symmetry=1, barrier=(19.7853,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.540237,0.0823862,-0.000117616,8.61308e-08,-2.51871e-11,-113393,27.0964], Tmin=(100,'K'), Tmax=(834.902,'K')), NASAPolynomial(coeffs=[12.686,0.0241965,-1.3072e-05,2.65378e-09,-1.91273e-13,-115422,-29.3085], Tmin=(834.902,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-943.796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsFFOO) + group(CsCFHO) + group(CJ2_singlet-FCs)"""),
)

species(
    label = '[O]OC(F)(OC(O)(F)F)C(F)F(2382)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {10,S} {12,S}
7  O u0 p2 c0 {9,S} {10,S}
8  O u0 p2 c0 {12,S} {14,S}
9  O u1 p2 c0 {7,S}
10 C u0 p0 c0 {1,S} {6,S} {7,S} {11,S}
11 C u0 p0 c0 {2,S} {3,S} {10,S} {13,S}
12 C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
13 H u0 p0 c0 {11,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-1496.84,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,375,361,584,565,722,1474,235,523,627,1123,1142,1372,1406,3097,435,565,619,662,854,1178,1396,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.66148,0.135576,-0.000231818,1.98295e-07,-6.61114e-11,-179835,37.1898], Tmin=(100,'K'), Tmax=(806.257,'K')), NASAPolynomial(coeffs=[16.9459,0.0318176,-1.74912e-05,3.47135e-09,-2.42943e-13,-182464,-46.2664], Tmin=(806.257,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1496.84,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(CsCFOO) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsFFOO) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(OC(O)(F)F)C(F)(F)F(2383)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {11,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {10,S} {12,S}
7  O u0 p2 c0 {9,S} {10,S}
8  O u0 p2 c0 {12,S} {14,S}
9  O u1 p2 c0 {7,S}
10 C u0 p0 c0 {6,S} {7,S} {11,S} {13,S}
11 C u0 p0 c0 {1,S} {2,S} {3,S} {10,S}
12 C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-1517.83,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (197.037,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.32597,0.127574,-0.000214406,1.83473e-07,-6.17003e-11,-182371,37.5663], Tmin=(100,'K'), Tmax=(787.123,'K')), NASAPolynomial(coeffs=[15.6271,0.0329948,-1.81086e-05,3.61388e-09,-2.54523e-13,-184779,-38.5059], Tmin=(787.123,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1517.83,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-CsOs) + group(CsFFOO) + radical(ROOJ)"""),
)

species(
    label = 'OC(F)(F)OF(421)',
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
    label = '[O]OC(F)=CF(659)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {4,S} {5,S}
4 O u1 p2 c0 {3,S}
5 C u0 p0 c0 {1,S} {3,S} {6,D}
6 C u0 p0 c0 {2,S} {5,D} {7,S}
7 H u0 p0 c0 {6,S}
"""),
    E0 = (-249.52,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,326,540,652,719,1357,194,682,905,1196,1383,3221],'cm^-1')),
        HinderedRotor(inertia=(0.358734,'amu*angstrom^2'), symmetry=1, barrier=(8.24799,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.0249,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.75962,0.0203814,3.21936e-05,-1.12324e-07,8.28e-11,-30007.8,12.0802], Tmin=(10,'K'), Tmax=(521.492,'K')), NASAPolynomial(coeffs=[5.63612,0.021414,-1.5147e-05,4.91851e-09,-5.9756e-13,-30413.3,2.2378], Tmin=(521.492,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-249.52,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""[O]OC(F)DCF""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]C(F)(F)C(F)OC(O)(F)F(2384)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {9,S} {11,S}
7  O u0 p2 c0 {11,S} {13,S}
8  O u1 p2 c0 {10,S}
9  C u0 p0 c0 {1,S} {6,S} {10,S} {12,S}
10 C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
11 C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-1461.83,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,261,493,600,1152,1365,1422,3097,351,323,533,609,664,892,1120,1201,435,565,619,662,854,1178,1396,238.34,238.951,243.2,247.961],'cm^-1')),
        HinderedRotor(inertia=(0.162674,'amu*angstrom^2'), symmetry=1, barrier=(6.55035,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.383102,'amu*angstrom^2'), symmetry=1, barrier=(15.8142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162535,'amu*angstrom^2'), symmetry=1, barrier=(6.56592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.84019,'amu*angstrom^2'), symmetry=1, barrier=(34.7297,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (181.038,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.45667,0.105521,-0.000151499,1.09406e-07,-3.13088e-11,-175664,32.5679], Tmin=(100,'K'), Tmax=(855.101,'K')), NASAPolynomial(coeffs=[16.1628,0.0277778,-1.51205e-05,3.07936e-09,-2.22379e-13,-178506,-45.0096], Tmin=(855.101,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1461.83,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(CsCFHO) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsFFOO) + radical(O2sj(Cs-CsF1sF1s))"""),
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
    label = 'OC(F)(F)OC(F)=C(F)F(2162)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {8,S} {9,S}
7  O u0 p2 c0 {8,S} {11,S}
8  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
9  C u0 p0 c0 {3,S} {6,S} {10,D}
10 C u0 p0 c0 {4,S} {5,S} {9,D}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-1346.53,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,435,565,619,662,854,1178,1396,326,540,652,719,1357,182,240,577,636,1210,1413,226.407,226.444,226.508],'cm^-1')),
        HinderedRotor(inertia=(0.339398,'amu*angstrom^2'), symmetry=1, barrier=(12.3481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.339323,'amu*angstrom^2'), symmetry=1, barrier=(12.3484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.803125,'amu*angstrom^2'), symmetry=1, barrier=(29.2157,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (164.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00896679,0.0924813,-0.000131397,9.06755e-08,-2.44343e-11,-161810,27.1707], Tmin=(100,'K'), Tmax=(912.139,'K')), NASAPolynomial(coeffs=[16.8193,0.0187633,-1.01695e-05,2.07253e-09,-1.50028e-13,-164877,-52.3836], Tmin=(912.139,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1346.53,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(CsFFOO) + group(CdCFO) + group(CdCFF) + longDistanceInteraction_noncyclic(Cds(F)2=Cds(F))"""),
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
    label = '[O]OC(F)(F)[CH]OC(O)(F)F(2385)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {10,S} {11,S}
6  O u0 p2 c0 {8,S} {9,S}
7  O u0 p2 c0 {10,S} {13,S}
8  O u1 p2 c0 {6,S}
9  C u0 p0 c0 {1,S} {2,S} {6,S} {11,S}
10 C u0 p0 c0 {3,S} {4,S} {5,S} {7,S}
11 C u1 p0 c0 {5,S} {9,S} {12,S}
12 H u0 p0 c0 {11,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-1112.02,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,253,525,597,667,842,1178,1324,435,565,619,662,854,1178,1396,3025,407.5,1350,352.5,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.42326,0.129686,-0.000226213,1.92367e-07,-6.29602e-11,-133559,35.7545], Tmin=(100,'K'), Tmax=(837.361,'K')), NASAPolynomial(coeffs=[17.3468,0.0255825,-1.38586e-05,2.70462e-09,-1.86297e-13,-136197,-48.4453], Tmin=(837.361,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1112.02,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsHH) + group(CsCFFO) + group(CsFFOO) + radical(ROOJ) + radical(CCsJOCs)"""),
)

species(
    label = '[O]O[C](F)C(F)OC(O)(F)F(2386)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {9,S} {10,S}
6  O u0 p2 c0 {10,S} {13,S}
7  O u0 p2 c0 {8,S} {11,S}
8  O u1 p2 c0 {7,S}
9  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
10 C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
11 C u1 p0 c0 {4,S} {7,S} {9,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-1074.25,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,487,638,688,1119,1325,1387,3149,435,565,619,662,854,1178,1396,395,473,707,1436,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.02927,0.121864,-0.00021207,1.86023e-07,-6.33824e-11,-129032,36.0284], Tmin=(100,'K'), Tmax=(816.522,'K')), NASAPolynomial(coeffs=[14.1567,0.0316967,-1.74484e-05,3.4607e-09,-2.41878e-13,-130986,-30.9368], Tmin=(816.522,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1074.25,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsFFOO) + radical(ROOJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]OC(F)(F)C(F)O[C](O)F(2387)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {9,S} {11,S}
6  O u0 p2 c0 {8,S} {10,S}
7  O u0 p2 c0 {11,S} {13,S}
8  O u1 p2 c0 {6,S}
9  C u0 p0 c0 {1,S} {5,S} {10,S} {12,S}
10 C u0 p0 c0 {2,S} {3,S} {6,S} {9,S}
11 C u1 p0 c0 {4,S} {5,S} {7,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-1061.34,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,261,493,600,1152,1365,1422,3097,223,363,546,575,694,1179,1410,482,586,761,1411,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.578445,0.107532,-0.000153624,1.09152e-07,-3.0572e-11,-127491,35.1488], Tmin=(100,'K'), Tmax=(875.115,'K')), NASAPolynomial(coeffs=[17.2283,0.0261372,-1.41044e-05,2.86213e-09,-2.06374e-13,-130607,-48.3827], Tmin=(875.115,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1061.34,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(CsCFHO) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsFHOO) + radical(ROOJ) + radical(Csj(F1s)(O2s-Cs)(O2s-H))"""),
)

species(
    label = '[O]C(O)(F)F(155)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {5,S} {6,S}
4 O u1 p2 c0 {5,S}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {4,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-638.905,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,375,486,586,623,710,938,1161,1294],'cm^-1')),
        HinderedRotor(inertia=(0.596909,'amu*angstrom^2'), symmetry=1, barrier=(13.7241,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0142,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3646.55,'J/mol'), sigma=(6.01309,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=569.58 K, Pc=38.06 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.90705,0.00545186,7.4862e-05,-1.61389e-07,9.87578e-11,-76836.8,10.4258], Tmin=(10,'K'), Tmax=(579.229,'K')), NASAPolynomial(coeffs=[5.15887,0.0182822,-1.39769e-05,4.86853e-09,-6.26572e-13,-77342.1,1.96031], Tmin=(579.229,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-638.905,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""[O]C(O)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]OC(F)(F)[CH]F(657)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {7,S}
4 O u0 p2 c0 {5,S} {6,S}
5 O u1 p2 c0 {4,S}
6 C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
7 C u1 p0 c0 {3,S} {6,S} {8,S}
8 H u0 p0 c0 {7,S}
"""),
    E0 = (-455.824,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,253,525,597,667,842,1178,1324,334,575,1197,1424,3202,180],'cm^-1')),
        HinderedRotor(inertia=(0.482095,'amu*angstrom^2'), symmetry=1, barrier=(11.0843,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.481426,'amu*angstrom^2'), symmetry=1, barrier=(11.0689,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.023,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3362.45,'J/mol'), sigma=(5.71223,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=525.21 K, Pc=40.93 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05144,0.0759341,-0.00015029,1.41518e-07,-4.95006e-11,-54727.4,21.7766], Tmin=(100,'K'), Tmax=(871.493,'K')), NASAPolynomial(coeffs=[7.94351,0.0195367,-1.05965e-05,2.05046e-09,-1.39101e-13,-54988.3,-5.13026], Tmin=(871.493,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-455.824,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(ROOJ) + radical(Csj(Cs-F1sF1sO2s)(F1s)(H))"""),
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
    label = '[O]OC(F)(F)C([O])F(2249)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {8,S}
4 O u0 p2 c0 {6,S} {7,S}
5 O u1 p2 c0 {8,S}
6 O u1 p2 c0 {4,S}
7 C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
8 C u0 p0 c0 {3,S} {5,S} {7,S} {9,S}
9 H u0 p0 c0 {8,S}
"""),
    E0 = (-632.549,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,223,363,546,575,694,1179,1410,391,562,707,872,1109,1210,1289,3137,180],'cm^-1')),
        HinderedRotor(inertia=(0.501054,'amu*angstrom^2'), symmetry=1, barrier=(11.5202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.500171,'amu*angstrom^2'), symmetry=1, barrier=(11.4999,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.023,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.96838,0.0716145,-0.00011463,9.37918e-08,-3.02314e-11,-75973.6,24.8435], Tmin=(100,'K'), Tmax=(780.676,'K')), NASAPolynomial(coeffs=[11.1816,0.0180434,-9.31463e-06,1.82045e-09,-1.27012e-13,-77530.4,-21.6587], Tmin=(780.676,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-632.549,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCFHO) + radical(O2sj(Cs-CsF1sH)) + radical(ROOJ)"""),
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
    label = 'OC(F)(F)OC(F)[C](F)F(1931)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {8,S} {9,S}
7  O u0 p2 c0 {9,S} {12,S}
8  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
9  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
10 C u1 p0 c0 {4,S} {5,S} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-1304.03,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,487,638,688,1119,1325,1387,3149,435,565,619,662,854,1178,1396,190,488,555,1236,1407,240.477,240.794,240.803,241.448],'cm^-1')),
        HinderedRotor(inertia=(0.229666,'amu*angstrom^2'), symmetry=1, barrier=(9.36423,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.227841,'amu*angstrom^2'), symmetry=1, barrier=(9.37175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.791851,'amu*angstrom^2'), symmetry=1, barrier=(32.3014,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.79063,'amu*angstrom^2'), symmetry=1, barrier=(32.302,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (165.039,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3545.48,'J/mol'), sigma=(5.86849,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=553.80 K, Pc=39.81 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0498013,0.0943451,-0.000135239,9.89087e-08,-2.88804e-11,-156702,30.2576], Tmin=(100,'K'), Tmax=(835.87,'K')), NASAPolynomial(coeffs=[14.0433,0.0273854,-1.50875e-05,3.08683e-09,-2.23448e-13,-159042,-34.7454], Tmin=(835.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1304.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCFHO) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsFFOO) + radical(Csj(Cs-F1sO2sH)(F1s)(F1s))"""),
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
    label = '[O]OC(F)(F)C(F)O[C](F)F(2388)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {9,S} {11,S}
7  O u0 p2 c0 {8,S} {10,S}
8  O u1 p2 c0 {7,S}
9  C u0 p0 c0 {1,S} {6,S} {10,S} {12,S}
10 C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
11 C u1 p0 c0 {4,S} {5,S} {6,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-1073.35,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,261,493,600,1152,1365,1422,3097,223,363,546,575,694,1179,1410,493,600,700,1144,1293,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.177579,'amu*angstrom^2'), symmetry=1, barrier=(4.0829,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.409636,'amu*angstrom^2'), symmetry=1, barrier=(9.41833,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174293,'amu*angstrom^2'), symmetry=1, barrier=(4.00734,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.04033,'amu*angstrom^2'), symmetry=1, barrier=(46.9112,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (180.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.629141,0.112755,-0.000197601,1.75461e-07,-6.05649e-11,-128938,34.7247], Tmin=(100,'K'), Tmax=(812.495,'K')), NASAPolynomial(coeffs=[12.773,0.030914,-1.72274e-05,3.43578e-09,-2.40999e-13,-130592,-23.9285], Tmin=(812.495,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1073.35,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(CsCFHO) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsFFHO) + radical(ROOJ) + radical(Csj(F1s)(F1s)(O2s-Cs))"""),
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
    label = '[O]OC(F)(F)C(F)OC([O])(F)F(2389)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {10,S} {12,S}
7  O u0 p2 c0 {9,S} {11,S}
8  O u1 p2 c0 {12,S}
9  O u1 p2 c0 {7,S}
10 C u0 p0 c0 {1,S} {6,S} {11,S} {13,S}
11 C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
12 C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-1231.78,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,261,493,600,1152,1365,1422,3097,223,363,546,575,694,1179,1410,375,486,586,623,710,938,1161,1294,249.166,249.168,249.169,249.169],'cm^-1')),
        HinderedRotor(inertia=(0.552746,'amu*angstrom^2'), symmetry=1, barrier=(24.3526,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169569,'amu*angstrom^2'), symmetry=1, barrier=(7.47046,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.995582,'amu*angstrom^2'), symmetry=1, barrier=(43.8619,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169565,'amu*angstrom^2'), symmetry=1, barrier=(7.47046,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (196.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.79929,0.113761,-0.000167989,1.23861e-07,-3.60746e-11,-147983,35.6282], Tmin=(100,'K'), Tmax=(841.557,'K')), NASAPolynomial(coeffs=[17.1815,0.0282999,-1.56669e-05,3.19816e-09,-2.30843e-13,-151010,-48.0177], Tmin=(841.557,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1231.78,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(CsCFHO) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsFFOO) + radical(O2sj(Cs-F1sF1sO2s)) + radical(ROOJ)"""),
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
    label = 'OC(F)(F)O[CH]F(602)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {7,S}
4 O u0 p2 c0 {6,S} {7,S}
5 O u0 p2 c0 {6,S} {9,S}
6 C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
7 C u1 p0 c0 {3,S} {4,S} {8,S}
8 H u0 p0 c0 {7,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-885.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,435,565,619,662,854,1178,1396,580,1155,1237,1373,3147,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.703585,'amu*angstrom^2'), symmetry=1, barrier=(16.1768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.703232,'amu*angstrom^2'), symmetry=1, barrier=(16.1687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.703342,'amu*angstrom^2'), symmetry=1, barrier=(16.1712,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.66975,0.0232301,0.000172568,-5.59225e-07,4.71644e-10,-106454,12.1225], Tmin=(10,'K'), Tmax=(448.854,'K')), NASAPolynomial(coeffs=[10.0657,0.023206,-1.77497e-05,6.23801e-09,-8.1149e-13,-107602,-20.0041], Tmin=(448.854,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-885.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), label="""OC(F)(F)O[CH]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]OC(F)(F)[C](F)OC(O)(F)F(2390)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {11,S} {12,S}
7  O u0 p2 c0 {9,S} {10,S}
8  O u0 p2 c0 {11,S} {13,S}
9  O u1 p2 c0 {7,S}
10 C u0 p0 c0 {1,S} {2,S} {7,S} {12,S}
11 C u0 p0 c0 {3,S} {4,S} {6,S} {8,S}
12 C u1 p0 c0 {5,S} {6,S} {10,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-1303.74,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,253,525,597,667,842,1178,1324,435,565,619,662,854,1178,1396,395,473,707,1436,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.55586,0.135742,-0.000246118,2.18225e-07,-7.41485e-11,-156616,38.2231], Tmin=(100,'K'), Tmax=(836.75,'K')), NASAPolynomial(coeffs=[15.5709,0.0310382,-1.74921e-05,3.46225e-09,-2.40242e-13,-158683,-36.5753], Tmin=(836.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1303.74,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(CsCFHO) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsFFOO) + radical(ROOJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]OC(F)(F)C(F)OC(=O)F(2391)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {9,S} {11,S}
6  O u0 p2 c0 {8,S} {10,S}
7  O u0 p2 c0 {11,D}
8  O u1 p2 c0 {6,S}
9  C u0 p0 c0 {1,S} {5,S} {10,S} {12,S}
10 C u0 p0 c0 {2,S} {3,S} {6,S} {9,S}
11 C u0 p0 c0 {4,S} {5,S} {7,D}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-1232.77,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,261,493,600,1152,1365,1422,3097,223,363,546,575,694,1179,1410,482,664,788,1296,1923,322.926,322.936,322.964,323.737],'cm^-1')),
        HinderedRotor(inertia=(0.223428,'amu*angstrom^2'), symmetry=1, barrier=(16.5352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.53616,'amu*angstrom^2'), symmetry=1, barrier=(39.664,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.535852,'amu*angstrom^2'), symmetry=1, barrier=(39.6632,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0781574,'amu*angstrom^2'), symmetry=1, barrier=(5.80405,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (177.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0509277,0.0936368,-0.000125988,8.54386e-08,-2.31281e-11,-148132,31.6508], Tmin=(100,'K'), Tmax=(899.103,'K')), NASAPolynomial(coeffs=[15.0196,0.0270436,-1.48893e-05,3.0619e-09,-2.22968e-13,-150824,-38.9724], Tmin=(899.103,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1232.77,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-OsCs) + group(O2s-OsH) + group(CsCFHO) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(COFOO) + radical(ROOJ)"""),
)

species(
    label = 'F2(77)',
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
    label = '[O]OC(F)=COC(O)(F)F(2392)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {10,S}
4  O u0 p2 c0 {8,S} {9,S}
5  O u0 p2 c0 {8,S} {12,S}
6  O u0 p2 c0 {7,S} {10,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
9  C u0 p0 c0 {4,S} {10,D} {11,S}
10 C u0 p0 c0 {3,S} {6,S} {9,D}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {5,S}
"""),
    E0 = (-898.579,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,435,565,619,662,854,1178,1396,3010,987.5,1337.5,450,1655,326,540,652,719,1357,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.92746,'amu*angstrom^2'), symmetry=1, barrier=(21.3241,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.516425,'amu*angstrom^2'), symmetry=1, barrier=(11.8736,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.517273,'amu*angstrom^2'), symmetry=1, barrier=(11.8931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.516294,'amu*angstrom^2'), symmetry=1, barrier=(11.8706,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (159.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.129349,0.0968285,-0.000142543,1.04634e-07,-3.01872e-11,-107931,31.4541], Tmin=(100,'K'), Tmax=(851.531,'K')), NASAPolynomial(coeffs=[15.6547,0.0226864,-1.19436e-05,2.3908e-09,-1.70681e-13,-110619,-42.1587], Tmin=(851.531,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-898.579,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + group(CsFFOO) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(F)=C(F)OC(O)(F)F(2393)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {9,S} {10,S}
6  O u0 p2 c0 {9,S} {12,S}
7  O u0 p2 c0 {8,S} {11,S}
8  O u1 p2 c0 {7,S}
9  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
10 C u0 p0 c0 {3,S} {5,S} {11,D}
11 C u0 p0 c0 {4,S} {7,S} {10,D}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-1085.66,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,435,565,619,662,854,1178,1396,262,390,483,597,572,732,631,807,1275,1439,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.456881,'amu*angstrom^2'), symmetry=1, barrier=(10.5046,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.376119,'amu*angstrom^2'), symmetry=1, barrier=(8.64772,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.380652,'amu*angstrom^2'), symmetry=1, barrier=(8.75195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.379156,'amu*angstrom^2'), symmetry=1, barrier=(8.71755,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (177.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.763409,0.115186,-0.000202194,1.76453e-07,-5.95902e-11,-130413,34.4321], Tmin=(100,'K'), Tmax=(821.635,'K')), NASAPolynomial(coeffs=[14.2536,0.0274191,-1.52005e-05,3.01344e-09,-2.10184e-13,-132386,-32.0547], Tmin=(821.635,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1085.66,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + group(CsFFOO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFO) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + radical(ROOJ)"""),
)

species(
    label = 'OOC(F)(F)[C](F)OC(O)(F)F(2394)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {11,S} {12,S}
7  O u0 p2 c0 {9,S} {10,S}
8  O u0 p2 c0 {11,S} {13,S}
9  O u0 p2 c0 {7,S} {14,S}
10 C u0 p0 c0 {1,S} {2,S} {7,S} {12,S}
11 C u0 p0 c0 {3,S} {4,S} {6,S} {8,S}
12 C u1 p0 c0 {5,S} {6,S} {10,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-1455.74,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3615,1277.5,1000,253,525,597,667,842,1178,1324,435,565,619,662,854,1178,1396,395,473,707,1436,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.77876,0.138843,-0.000239725,2.06943e-07,-6.97438e-11,-174888,38.2652], Tmin=(100,'K'), Tmax=(796.291,'K')), NASAPolynomial(coeffs=[17.0159,0.0329531,-1.86344e-05,3.73952e-09,-2.63451e-13,-177518,-45.8419], Tmin=(796.291,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1455.74,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(CsCFHO) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsFFOO) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]C(F)(F)OC(F)C(F)(F)OO(2395)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {10,S} {12,S}
7  O u0 p2 c0 {8,S} {11,S}
8  O u0 p2 c0 {7,S} {14,S}
9  O u1 p2 c0 {12,S}
10 C u0 p0 c0 {1,S} {6,S} {11,S} {13,S}
11 C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
12 C u0 p0 c0 {4,S} {5,S} {6,S} {9,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-1383.78,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,261,493,600,1152,1365,1422,3097,223,363,546,575,694,1179,1410,375,486,586,623,710,938,1161,1294,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.898906,0.115351,-0.000155961,1.04526e-07,-2.77676e-11,-166261,35.2314], Tmin=(100,'K'), Tmax=(918.324,'K')), NASAPolynomial(coeffs=[18.6657,0.0301284,-1.67503e-05,3.4598e-09,-2.52645e-13,-169854,-57.4887], Tmin=(918.324,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1383.78,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(CsCFHO) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsFFOO) + radical(O2sj(Cs-F1sF1sO2s))"""),
)

species(
    label = 'OC(F)(F)O[CH]C(F)(F)OOF(2396)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {9,S}
6  O u0 p2 c0 {11,S} {12,S}
7  O u0 p2 c0 {9,S} {10,S}
8  O u0 p2 c0 {11,S} {14,S}
9  O u0 p2 c0 {5,S} {7,S}
10 C u0 p0 c0 {1,S} {2,S} {7,S} {12,S}
11 C u0 p0 c0 {3,S} {4,S} {6,S} {8,S}
12 C u1 p0 c0 {6,S} {10,S} {13,S}
13 H u0 p0 c0 {12,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-1168.72,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3615,1277.5,1000,277,555,632,253,525,597,667,842,1178,1324,435,565,619,662,854,1178,1396,3025,407.5,1350,352.5,200],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.73318,0.138575,-0.00022512,1.73048e-07,-4.88023e-11,-140370,37.5286], Tmin=(100,'K'), Tmax=(656.507,'K')), NASAPolynomial(coeffs=[19.2328,0.0290251,-1.63844e-05,3.29278e-09,-2.32652e-13,-143515,-57.7833], Tmin=(656.507,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1168.72,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2sFO) + group(Cs-CsOsHH) + group(CsCFFO) + group(CsFFOO) + radical(CCsJOCs)"""),
)

species(
    label = 'OC(F)(F)OC(F)[C](F)OOF(2397)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {9,S}
6  O u0 p2 c0 {10,S} {11,S}
7  O u0 p2 c0 {9,S} {12,S}
8  O u0 p2 c0 {11,S} {14,S}
9  O u0 p2 c0 {5,S} {7,S}
10 C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
11 C u0 p0 c0 {2,S} {3,S} {6,S} {8,S}
12 C u1 p0 c0 {4,S} {7,S} {10,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-1130.95,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3615,1277.5,1000,277,555,632,487,638,688,1119,1325,1387,3149,435,565,619,662,854,1178,1396,395,473,707,1436,200],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.17046,0.128103,-0.000197431,1.39564e-07,-3.07937e-11,-135850,37.2325], Tmin=(100,'K'), Tmax=(596.837,'K')), NASAPolynomial(coeffs=[16.0266,0.035177,-2.00007e-05,4.05601e-09,-2.88883e-13,-138300,-40.1892], Tmin=(596.837,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1130.95,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2sFO) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsFFOO) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'O[C](F)OC(F)C(F)(F)OOF(2398)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {9,S}
6  O u0 p2 c0 {10,S} {12,S}
7  O u0 p2 c0 {9,S} {11,S}
8  O u0 p2 c0 {12,S} {14,S}
9  O u0 p2 c0 {5,S} {7,S}
10 C u0 p0 c0 {1,S} {6,S} {11,S} {13,S}
11 C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
12 C u1 p0 c0 {4,S} {6,S} {8,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-1118.04,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3615,1277.5,1000,277,555,632,261,493,600,1152,1365,1422,3097,223,363,546,575,694,1179,1410,482,586,761,1411,200],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.05275,0.118936,-0.000164732,1.12392e-07,-3.02843e-11,-134294,37.48], Tmin=(100,'K'), Tmax=(907.219,'K')), NASAPolynomial(coeffs=[19.294,0.0292249,-1.64036e-05,3.39278e-09,-2.477e-13,-137986,-58.7001], Tmin=(907.219,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1118.04,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2sFO) + group(CsCFHO) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsFHOO) + radical(Csj(F1s)(O2s-Cs)(O2s-H))"""),
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
    E0 = (-678.083,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-318.795,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-349.595,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-203.65,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-494.923,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-496.372,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-70.4962,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-515.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-636.718,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-336.206,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-298.438,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-285.527,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-391.809,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-399.923,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-609.732,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-342.032,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-317.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-377.91,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-389.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-633.85,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-204.464,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-497.353,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-668.933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-608.65,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (-393.719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (-342.802,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (-355.214,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]OC(F)(F)C(F)OC(O)(F)F(2298)'],
    products = ['HF(38)', 'CF2O(48)', '[O]OC(F)(F)C=O(622)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(6.14547e+10,'s^-1'), n=0.631481, Ea=(117.288,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R!H->C_5Br1sCl1sF1sH->F1s_Ext-2C-R_N-7R!H->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R!H->C_5Br1sCl1sF1sH->F1s_Ext-2C-R_N-7R!H->C_Ext-2C-R"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CHF(40)', '[O]OC(F)(F)OC(O)(F)F(2380)'],
    products = ['[O]OC(F)(F)C(F)OC(O)(F)F(2298)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.12553e-07,'m^3/(mol*s)'), n=3.34134, Ea=(134.127,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CO_2Br1sCl1sF1sHI1s->F1s_N-3Br1sCCl1sF1sHI1s->F1s',), comment="""Estimated from node CO_2Br1sCl1sF1sHI1s->F1s_N-3Br1sCCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CF2(42)', '[O]OC(F)OC(O)(F)F(2381)'],
    products = ['[O]OC(F)(F)C(F)OC(O)(F)F(2298)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.40908e-07,'m^3/(mol*s)'), n=3.45227, Ea=(213.672,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s',), comment="""Estimated from node CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]OF(124)', 'OC(F)(F)OC(F)[C]F-2(2167)'],
    products = ['[O]OC(F)(F)C(F)OC(O)(F)F(2298)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.13916e+09,'m^3/(mol*s)'), n=-1.00842, Ea=(23.6309,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17406617270594374, var=16.167420070870673, Tref=1000.0, N=4, data_mean=0.0, correlation='OHOY',), comment="""Estimated from node OHOY"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]OC(F)(OC(O)(F)F)C(F)F(2382)'],
    products = ['[O]OC(F)(F)C(F)OC(O)(F)F(2298)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]OC(F)(F)C(F)OC(O)(F)F(2298)'],
    products = ['[O]OC(OC(O)(F)F)C(F)(F)F(2383)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction7',
    reactants = ['OC(F)(F)OF(421)', '[O]OC(F)=CF(659)'],
    products = ['[O]OC(F)(F)C(F)OC(O)(F)F(2298)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.000170623,'m^3/(mol*s)'), n=2.84144, Ea=(220.872,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OR] for rate rule [Cd/disub_Cd/monosub;H_OR]
Euclidian distance = 1.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction8',
    reactants = ['O(6)', '[O]C(F)(F)C(F)OC(O)(F)F(2384)'],
    products = ['[O]OC(F)(F)C(F)OC(O)(F)F(2298)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]OC(F)(F)C(F)OC(O)(F)F(2298)'],
    products = ['HO2(10)', 'OC(F)(F)OC(F)=C(F)F(2162)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(9.58174e+10,'s^-1'), n=0.573333, Ea=(158.653,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F(37)', '[O]OC(F)(F)[CH]OC(O)(F)F(2385)'],
    products = ['[O]OC(F)(F)C(F)OC(O)(F)F(2298)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', '[O]O[C](F)C(F)OC(O)(F)F(2386)'],
    products = ['[O]OC(F)(F)C(F)OC(O)(F)F(2298)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', '[O]OC(F)(F)C(F)O[C](O)F(2387)'],
    products = ['[O]OC(F)(F)C(F)OC(O)(F)F(2298)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]C(O)(F)F(155)', '[O]OC(F)(F)[CH]F(657)'],
    products = ['[O]OC(F)(F)C(F)OC(O)(F)F(2298)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['O[C](F)F(239)', '[O]OC(F)(F)C([O])F(2249)'],
    products = ['[O]OC(F)(F)C(F)OC(O)(F)F(2298)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.505e+06,'m^3/(mol*s)'), n=-1.3397e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Sp-4R!H-2C_Ext-2C-R_N-5R!H->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Sp-4R!H-2C_Ext-2C-R_N-5R!H->C"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O2(2)', 'OC(F)(F)OC(F)[C](F)F(1931)'],
    products = ['[O]OC(F)(F)C(F)OC(O)(F)F(2298)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.6844e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction16',
    reactants = ['OH(4)', '[O]OC(F)(F)C(F)O[C](F)F(2388)'],
    products = ['[O]OC(F)(F)C(F)OC(O)(F)F(2298)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(5)', '[O]OC(F)(F)C(F)OC([O])(F)F(2389)'],
    products = ['[O]OC(F)(F)C(F)OC(O)(F)F(2298)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]O[C](F)F(171)', 'OC(F)(F)O[CH]F(602)'],
    products = ['[O]OC(F)(F)C(F)OC(O)(F)F(2298)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(5)', '[O]OC(F)(F)[C](F)OC(O)(F)F(2390)'],
    products = ['[O]OC(F)(F)C(F)OC(O)(F)F(2298)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(43950,'m^3/(mol*s)'), n=1, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction20',
    reactants = ['HF(38)', '[O]OC(F)(F)C(F)OC(=O)F(2391)'],
    products = ['[O]OC(F)(F)C(F)OC(O)(F)F(2298)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.109156,'m^3/(mol*s)'), n=1.86531, Ea=(177.116,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd',), comment="""Estimated from node HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd"""),
)

reaction(
    label = 'reaction21',
    reactants = ['F2(77)', '[O]OC(F)=COC(O)(F)F(2392)'],
    products = ['[O]OC(F)(F)C(F)OC(O)(F)F(2298)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction22',
    reactants = ['HF(38)', '[O]OC(F)=C(F)OC(O)(F)F(2393)'],
    products = ['[O]OC(F)(F)C(F)OC(O)(F)F(2298)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.14111,'m^3/(mol*s)'), n=1.29695, Ea=(166.499,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction23',
    reactants = ['OOC(F)(F)[C](F)OC(O)(F)F(2394)'],
    products = ['[O]OC(F)(F)C(F)OC(O)(F)F(2298)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_single;O_H_out] for rate rule [R4H_SSS;C_rad_out_noH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]C(F)(F)OC(F)C(F)(F)OO(2395)'],
    products = ['[O]OC(F)(F)C(F)OC(O)(F)F(2298)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(627096,'s^-1'), n=1.03067, Ea=(72.2138,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7H;O_rad_out;XH_out] for rate rule [R7H;O_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['OC(F)(F)O[CH]C(F)(F)OOF(2396)'],
    products = ['[O]OC(F)(F)C(F)OC(O)(F)F(2298)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(72.0826,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction26',
    reactants = ['OC(F)(F)OC(F)[C](F)OOF(2397)'],
    products = ['[O]OC(F)(F)C(F)OC(O)(F)F(2298)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(85.232,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction27',
    reactants = ['O[C](F)OC(F)C(F)(F)OOF(2398)'],
    products = ['[O]OC(F)(F)C(F)OC(O)(F)F(2298)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(59.9083,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

network(
    label = 'PDepNetwork #910',
    isomers = [
        '[O]OC(F)(F)C(F)OC(O)(F)F(2298)',
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
    label = 'PDepNetwork #910',
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

