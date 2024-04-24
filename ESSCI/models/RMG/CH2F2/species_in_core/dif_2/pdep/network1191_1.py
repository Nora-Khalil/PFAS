species(
    label = 'O=C=C(F)OC(O)(F)F(3015)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u0 p2 c0 {7,S} {10,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
8  C u0 p0 c0 {3,S} {4,S} {9,D}
9  C u0 p0 c0 {6,D} {8,D}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-1023.4,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,435,565,619,662,854,1178,1396,197,221,431,657,2120,512.5,787.5,180,180,180,1009.18],'cm^-1')),
        HinderedRotor(inertia=(0.692969,'amu*angstrom^2'), symmetry=1, barrier=(15.9327,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.19309,'amu*angstrom^2'), symmetry=1, barrier=(15.9329,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0417858,'amu*angstrom^2'), symmetry=1, barrier=(30.1977,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.51135,0.0778693,-0.000102604,6.4522e-08,-1.57216e-11,-122961,25.3783], Tmin=(100,'K'), Tmax=(1010.05,'K')), NASAPolynomial(coeffs=[16.767,0.0134939,-7.00219e-06,1.42155e-09,-1.03435e-13,-126245,-53.2081], Tmin=(1010.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1023.4,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + missing(O2d-Cdd) + group(CsFFOO) + group(Cd(Cdd-Od)FO) + missing(Cdd-CdO2d)"""),
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
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.615836,0.0804929,-0.000126656,1.01562e-07,-3.22185e-11,-137730,26.3539], Tmin=(100,'K'), Tmax=(774.354,'K')), NASAPolynomial(coeffs=[12.1398,0.020962,-1.13334e-05,2.27231e-09,-1.61471e-13,-139515,-26.2952], Tmin=(774.354,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1146.12,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cds-O2d(Cds-O2d)Cs) + group(COCFO)"""),
)

species(
    label = '[O]C(F)(F)OC(F)[C]=O(3069)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u1 p2 c0 {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
9  C u1 p0 c0 {6,D} {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-765.924,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([355,410,600,1181,1341,1420,3056,375,486,586,623,710,938,1161,1294,1855,455,950,198.179,198.199,198.212],'cm^-1')),
        HinderedRotor(inertia=(0.265676,'amu*angstrom^2'), symmetry=1, barrier=(7.39079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.973355,'amu*angstrom^2'), symmetry=1, barrier=(27.1364,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.265189,'amu*angstrom^2'), symmetry=1, barrier=(7.39127,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.460062,0.0836366,-0.000133177,1.06851e-07,-3.36959e-11,-91997.3,27.8612], Tmin=(100,'K'), Tmax=(780.752,'K')), NASAPolynomial(coeffs=[12.9286,0.0197592,-1.0459e-05,2.06892e-09,-1.45551e-13,-93944.4,-29.207], Tmin=(780.752,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-765.924,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsFFOO) + group(Cds-OdCsH) + radical(O2sj(Cs-F1sF1sO2s)) + radical(COj(Cs-F1sO2sH)(O2d))"""),
)

species(
    label = '[O]C(F)(F)O[C](F)C=O(3008)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
8  C u1 p0 c0 {3,S} {4,S} {9,S}
9  C u0 p0 c0 {6,D} {8,S} {10,S}
10 H u0 p0 c0 {9,S}
"""),
    E0 = (-785.24,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,486,586,623,710,938,1161,1294,280,501,1494,1531,2782.5,750,1395,475,1775,1000,180,180,830.588],'cm^-1')),
        HinderedRotor(inertia=(2.35272,'amu*angstrom^2'), symmetry=1, barrier=(54.0936,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124557,'amu*angstrom^2'), symmetry=1, barrier=(2.8638,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.35097,'amu*angstrom^2'), symmetry=1, barrier=(54.0533,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3887.01,'J/mol'), sigma=(6.28312,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=607.14 K, Pc=35.56 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2617,0.0658876,-8.42604e-05,5.81342e-08,-1.64939e-11,-94348.8,25.5385], Tmin=(100,'K'), Tmax=(849.293,'K')), NASAPolynomial(coeffs=[9.73811,0.0259652,-1.37499e-05,2.78548e-09,-2.01193e-13,-95788.6,-13.9706], Tmin=(849.293,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-785.24,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsFFOO) + group(Cds-OdCsH) + radical(O2sj(Cs-F1sF1sO2s)) + radical(CsCOF1sO2s)"""),
)

species(
    label = '[O]C(F)(F)OC(F)=[C]O(3076)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u0 p2 c0 {9,S} {10,S}
6  O u1 p2 c0 {7,S}
7  C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
8  C u0 p0 c0 {3,S} {4,S} {9,D}
9  C u1 p0 c0 {5,S} {8,D}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-668.822,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,375,486,586,623,710,938,1161,1294,293,496,537,1218,1685,370,180,205.451,531.256,531.695],'cm^-1')),
        HinderedRotor(inertia=(0.812248,'amu*angstrom^2'), symmetry=1, barrier=(18.6752,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0934343,'amu*angstrom^2'), symmetry=1, barrier=(18.6805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.812267,'amu*angstrom^2'), symmetry=1, barrier=(18.6756,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.244792,0.0796647,-0.000102294,6.10071e-08,-1.3892e-11,-80302.8,29.3243], Tmin=(100,'K'), Tmax=(1087.73,'K')), NASAPolynomial(coeffs=[19.6062,0.00846499,-4.10831e-06,8.28999e-10,-6.08643e-14,-84514.8,-65.7115], Tmin=(1087.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-668.822,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(CsFFOO) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + radical(O2sj(Cs-F1sF1sO2s)) + radical(C=CJO)"""),
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
    label = 'O=C=C(F)O[C](O)F(3253)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {6,S} {7,S}
4 O u0 p2 c0 {6,S} {9,S}
5 O u0 p2 c0 {8,D}
6 C u1 p0 c0 {1,S} {3,S} {4,S}
7 C u0 p0 c0 {2,S} {3,S} {8,D}
8 C u0 p0 c0 {5,D} {7,D}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (-585.088,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,482,586,761,1411,197,221,431,657,2120,512.5,787.5,180,1032.59,1033.53,1033.66],'cm^-1')),
        HinderedRotor(inertia=(0.694378,'amu*angstrom^2'), symmetry=1, barrier=(15.9651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.694492,'amu*angstrom^2'), symmetry=1, barrier=(15.9677,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.399687,'amu*angstrom^2'), symmetry=1, barrier=(41.1481,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04837,0.0647491,-8.221e-05,4.99608e-08,-1.17585e-11,-70263.3,24.9722], Tmin=(100,'K'), Tmax=(1045.62,'K')), NASAPolynomial(coeffs=[15.0067,0.0113513,-5.60758e-06,1.12038e-09,-8.10293e-14,-73182.3,-42.9915], Tmin=(1045.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-585.088,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + missing(O2d-Cdd) + group(CsFHOO) + group(Cd(Cdd-Od)FO) + missing(Cdd-CdO2d) + radical(CsF1sO2sO2s)"""),
)

species(
    label = 'O=C=[C]OC(O)(F)F(3254)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {6,S} {7,S}
4 O u0 p2 c0 {6,S} {9,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {2,S} {3,S} {4,S}
7 C u1 p0 c0 {3,S} {8,D}
8 C u0 p0 c0 {5,D} {7,D}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (-595.447,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,435,565,619,662,854,1178,1396,1685,370,2120,512.5,787.5,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.813631,'amu*angstrom^2'), symmetry=1, barrier=(18.707,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.813882,'amu*angstrom^2'), symmetry=1, barrier=(18.7128,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.813664,'amu*angstrom^2'), symmetry=1, barrier=(18.7077,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.687323,0.0714899,-9.83206e-05,6.30752e-08,-1.5406e-11,-71495.1,25.6247], Tmin=(100,'K'), Tmax=(1017.21,'K')), NASAPolynomial(coeffs=[17.1397,0.00679322,-2.91662e-06,5.48067e-10,-3.85481e-14,-74842.2,-54.029], Tmin=(1017.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-595.447,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + missing(O2d-Cdd) + group(CsFFOO) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(C=CJO)"""),
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
    label = 'O=[C]C(=O)F(997)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {4,D}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {2,D} {5,S}
5 C u1 p0 c0 {3,D} {4,S}
"""),
    E0 = (-325.34,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([233,496,705,1150,2014,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.06462,'amu*angstrom^2'), symmetry=1, barrier=(24.4776,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (75.0185,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.90371,0.00648742,6.66758e-05,-1.78696e-07,1.33913e-10,-39126.2,10.0186], Tmin=(10,'K'), Tmax=(477.098,'K')), NASAPolynomial(coeffs=[4.98548,0.0142778,-1.08248e-05,3.66779e-09,-4.58866e-13,-39421.3,3.58924], Tmin=(477.098,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-325.34,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""OD[C]C(DO)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]C(O)(F)F(323)',
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
    label = 'O=C=[C]F(965)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {4,D}
3 C u1 p0 c0 {1,S} {4,D}
4 C u0 p0 c0 {2,D} {3,D}
"""),
    E0 = (33.6712,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(58.9933,'amu')),
        NonlinearRotor(inertia=([3.58812,117.154,120.742],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([280.764,365.894,563.98,866.433,1429.66,2068.95],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (59.0191,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.94745,0.00353157,4.23087e-05,-1.10022e-07,8.19717e-11,4052.13,8.19715], Tmin=(10,'K'), Tmax=(472.411,'K')), NASAPolynomial(coeffs=[4.41104,0.00925722,-6.51505e-06,2.12227e-09,-2.60232e-13,3900.64,5.16846], Tmin=(472.411,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(33.6712,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""ODCD[C]F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'O=C=C(F)O[C](F)F(3255)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {7,S}
4 O u0 p2 c0 {6,S} {7,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {4,S} {8,D}
7 C u1 p0 c0 {2,S} {3,S} {4,S}
8 C u0 p0 c0 {5,D} {6,D}
"""),
    E0 = (-595.86,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([197,221,431,657,493,600,700,1144,1293,2120,512.5,787.5,383.47,383.478,383.479,383.48],'cm^-1')),
        HinderedRotor(inertia=(0.00114644,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00114634,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (125.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51138,0.0549321,-6.58003e-05,3.76386e-08,-8.40599e-12,-71575.9,23.5316], Tmin=(100,'K'), Tmax=(1093.97,'K')), NASAPolynomial(coeffs=[13.4173,0.0113984,-6.10782e-06,1.26129e-09,-9.26963e-14,-74180.8,-34.9767], Tmin=(1093.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-595.86,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + missing(O2d-Cdd) + group(CsFFHO) + group(Cd(Cdd-Od)FO) + missing(Cdd-CdO2d) + radical(Csj(F1s)(F1s)(O2s-Cd))"""),
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
    label = '[O]C(F)(F)OC(F)=C=O(3068)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {8,S}
4 O u0 p2 c0 {7,S} {8,S}
5 O u1 p2 c0 {7,S}
6 O u0 p2 c0 {9,D}
7 C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
8 C u0 p0 c0 {3,S} {4,S} {9,D}
9 C u0 p0 c0 {6,D} {8,D}
"""),
    E0 = (-756.885,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,486,586,623,710,938,1161,1294,197,221,431,657,2120,512.5,787.5,180,180,605.291,1032.51],'cm^-1')),
        HinderedRotor(inertia=(0.0298094,'amu*angstrom^2'), symmetry=1, barrier=(22.5229,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.980101,'amu*angstrom^2'), symmetry=1, barrier=(22.5344,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.025,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.693978,0.0700419,-8.6691e-05,4.99381e-08,-1.10312e-11,-90910.8,25.6163], Tmin=(100,'K'), Tmax=(1117.55,'K')), NASAPolynomial(coeffs=[17.7928,0.0088409,-4.54605e-06,9.35242e-10,-6.91167e-14,-94732.6,-58.7761], Tmin=(1117.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-756.885,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + missing(O2d-Cdd) + group(CsFFOO) + group(Cd(Cdd-Od)FO) + missing(Cdd-CdO2d) + radical(O2sj(Cs-F1sF1sO2s))"""),
)

species(
    label = 'O=C(F)[C-]=[O+]C(O)(F)F(3256)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p1 c+1 {7,S} {9,D}
5  O u0 p2 c0 {7,S} {10,S}
6  O u0 p2 c0 {8,D}
7  C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
8  C u0 p0 c0 {3,S} {6,D} {9,S}
9  C u0 p1 c-1 {4,D} {8,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-779.65,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.75996,'amu*angstrom^2'), symmetry=1, barrier=(40.4649,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.74693,'amu*angstrom^2'), symmetry=1, barrier=(40.1655,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.76127,'amu*angstrom^2'), symmetry=1, barrier=(40.4951,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.282343,0.0927131,-0.000168356,1.54637e-07,-5.44044e-11,-93647,38.0415], Tmin=(100,'K'), Tmax=(837.233,'K')), NASAPolynomial(coeffs=[9.3528,0.0284084,-1.55779e-05,3.06794e-09,-2.12879e-13,-94430.9,0.282204], Tmin=(837.233,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-779.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsFFOO) + group(COCFO) + group(CsJ2_singlet-CsH)"""),
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
    label = 'O=C=C(F)OC(=O)F(3257)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {6,S} {7,S}
4 O u0 p2 c0 {7,D}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {3,S} {8,D}
7 C u0 p0 c0 {2,S} {3,S} {4,D}
8 C u0 p0 c0 {5,D} {6,D}
"""),
    E0 = (-732.326,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([197,221,431,657,482,664,788,1296,1923,2120,512.5,787.5,180,180,898.987,899.048],'cm^-1')),
        HinderedRotor(inertia=(0.021069,'amu*angstrom^2'), symmetry=1, barrier=(4.74746,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31462,'amu*angstrom^2'), symmetry=1, barrier=(30.2258,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (122.027,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45625,0.0612953,-9.1904e-05,7.0404e-08,-2.1563e-11,-87991.6,21.8662], Tmin=(100,'K'), Tmax=(797.474,'K')), NASAPolynomial(coeffs=[10.1152,0.0178631,-1.021e-05,2.10936e-09,-1.5316e-13,-89372.7,-17.9485], Tmin=(797.474,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-732.326,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + missing(O2d-Cdd) + group(Cd(Cdd-Od)FO) + group(COFOO) + missing(Cdd-CdO2d)"""),
)

species(
    label = 'O=C=C=O(2992)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,D}
2 O u0 p2 c0 {4,D}
3 C u0 p0 c0 {1,D} {4,D}
4 C u0 p0 c0 {2,D} {3,D}
"""),
    E0 = (63.731,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2110,2130,495,530,650,925,180],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.0201,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5355,0.0240331,-4.26141e-05,3.87328e-08,-1.34308e-11,7697.1,25.5769], Tmin=(100,'K'), Tmax=(857.391,'K')), NASAPolynomial(coeffs=[4.78325,0.0077353,-3.93445e-06,7.52105e-10,-5.12481e-14,7525.26,16.3242], Tmin=(857.391,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(63.731,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(87.302,'J/(mol*K)'), label="""OCCO(S)""", comment="""Thermo library: primaryThermoLibrary"""),
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
    label = 'OC(F)(F)[O+]=[C-]F(3258)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {7,S}
4 O u0 p1 c+1 {6,S} {7,D}
5 O u0 p2 c0 {6,S} {8,S}
6 C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
7 C u0 p1 c-1 {3,S} {4,D}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-860.791,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.023,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.84435,0.0503585,-6.10027e-05,3.72391e-08,-9.11806e-12,-103454,18.2547], Tmin=(100,'K'), Tmax=(987.51,'K')), NASAPolynomial(coeffs=[10.3981,0.0157122,-8.37843e-06,1.71405e-09,-1.24874e-13,-105144,-22.9051], Tmin=(987.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-860.791,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsFFOO) + group(CJ2_singlet-FO)"""),
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
    E0 = (-408.773,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-421.106,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-244.975,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-264.29,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-147.872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-16.2205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-26.5787,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-293.064,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-109.257,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-71.4894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-49.1039,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-167.035,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-353.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-313.748,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-460.952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=C=C(F)OC(O)(F)F(3015)'],
    products = ['CF2O(48)', 'O=CC(=O)F(557)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(4.02942e+10,'s^-1'), n=0.375329, Ea=(118.648,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.06684469846581474, var=24.965000021544366, Tref=1000.0, N=36, data_mean=0.0, correlation='Root_N-1R!H->C_N-5R!H->N',), comment="""Estimated from node Root_N-1R!H->C_N-5R!H->N"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O=C=C(F)OC(O)(F)F(3015)'],
    products = ['O=C(F)C(=O)C(O)(F)F(3107)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(106.315,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]C(F)(F)OC(F)[C]=O(3069)'],
    products = ['O=C=C(F)OC(O)(F)F(3015)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]C(F)(F)O[C](F)C=O(3008)'],
    products = ['O=C=C(F)OC(O)(F)F(3015)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]C(F)(F)OC(F)=[C]O(3076)'],
    products = ['O=C=C(F)OC(O)(F)F(3015)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['F(37)', 'O=C=C(F)O[C](O)F(3253)'],
    products = ['O=C=C(F)OC(O)(F)F(3015)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F(37)', 'O=C=[C]OC(O)(F)F(3254)'],
    products = ['O=C=C(F)OC(O)(F)F(3015)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.52887e+10,'m^3/(mol*s)'), n=-0.421056, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.030895812897821735, var=3.1393582276389975, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_Ext-5R!H-R"""),
)

reaction(
    label = 'reaction8',
    reactants = ['O[C](F)F(232)', 'O=[C]C(=O)F(997)'],
    products = ['O=C=C(F)OC(O)(F)F(3015)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(6.59276,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C(O)(F)F(323)', 'O=C=[C]F(965)'],
    products = ['O=C=C(F)OC(O)(F)F(3015)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction10',
    reactants = ['OH(4)', 'O=C=C(F)O[C](F)F(3255)'],
    products = ['O=C=C(F)OC(O)(F)F(3015)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.7e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(5)', '[O]C(F)(F)OC(F)=C=O(3068)'],
    products = ['O=C=C(F)OC(O)(F)F(3015)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction12',
    reactants = ['O=C(F)[C-]=[O+]C(O)(F)F(3256)'],
    products = ['O=C=C(F)OC(O)(F)F(3015)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.91033e+10,'s^-1'), n=0.827, Ea=(116.639,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CCY',), comment="""Estimated from node CCY"""),
)

reaction(
    label = 'reaction13',
    reactants = ['HF(38)', 'O=C=C(F)OC(=O)F(3257)'],
    products = ['O=C=C(F)OC(O)(F)F(3015)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.109156,'m^3/(mol*s)'), n=1.86531, Ea=(163.71,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd',), comment="""Estimated from node HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd"""),
)

reaction(
    label = 'reaction14',
    reactants = ['O=C=C(F)OC(O)(F)F(3015)'],
    products = ['HF(38)', 'CF2O(48)', 'O=C=C=O(2992)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(6.14547e+10,'s^-1'), n=0.631481, Ea=(213.673,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R!H->C_5Br1sCl1sF1sH->F1s_Ext-2C-R_N-7R!H->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R!H->C_5Br1sCl1sF1sH->F1s_Ext-2C-R_N-7R!H->C_Ext-2C-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O=C=C(F)OC(O)(F)F(3015)'],
    products = ['CO(12)', 'OC(F)(F)[O+]=[C-]F(3258)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.15699e+14,'s^-1'), n=0.0573689, Ea=(66.4688,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='RFC=C=O',), comment="""Estimated from node RFC=C=O"""),
)

network(
    label = 'PDepNetwork #1191',
    isomers = [
        'O=C=C(F)OC(O)(F)F(3015)',
    ],
    reactants = [
        ('CF2O(48)', 'O=CC(=O)F(557)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1191',
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

