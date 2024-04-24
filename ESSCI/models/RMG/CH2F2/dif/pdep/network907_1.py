species(
    label = '[O]OC(F)(F)C(O)OC(F)(F)F(2295)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {11,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {10,S} {12,S}
7  O u0 p2 c0 {10,S} {14,S}
8  O u0 p2 c0 {9,S} {11,S}
9  O u1 p2 c0 {8,S}
10 C u0 p0 c0 {6,S} {7,S} {11,S} {13,S}
11 C u0 p0 c0 {1,S} {2,S} {8,S} {10,S}
12 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (-1512.92,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,431,527,613,668,835,1199,1245,1322,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.19867,0.124775,-0.000210575,1.81055e-07,-6.09623e-11,-181785,38.7659], Tmin=(100,'K'), Tmax=(801.66,'K')), NASAPolynomial(coeffs=[15.0434,0.0327122,-1.7695e-05,3.5061e-09,-2.45665e-13,-184035,-33.7933], Tmin=(801.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1512.92,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(CsCFFO) + group(CsFFFO) + radical(ROOJ)"""),
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
    label = '[O]OC(O)OC(F)(F)F(2399)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {8,S} {9,S}
5  O u0 p2 c0 {8,S} {11,S}
6  O u0 p2 c0 {7,S} {8,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {4,S} {5,S} {6,S} {10,S}
9  C u0 p0 c0 {1,S} {2,S} {3,S} {4,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-1036.58,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,1380,1390,370,380,2900,435,431,527,613,668,835,1199,1245,1322,180,180,3177.58],'cm^-1')),
        HinderedRotor(inertia=(0.224506,'amu*angstrom^2'), symmetry=1, barrier=(5.16183,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.224598,'amu*angstrom^2'), symmetry=1, barrier=(5.16394,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.224525,'amu*angstrom^2'), symmetry=1, barrier=(5.16226,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.224546,'amu*angstrom^2'), symmetry=1, barrier=(5.16274,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (147.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.622562,0.0878831,-0.00016727,1.62797e-07,-5.96843e-11,-124563,29.8607], Tmin=(100,'K'), Tmax=(845.474,'K')), NASAPolynomial(coeffs=[5.71108,0.0339082,-1.84616e-05,3.63009e-09,-2.51453e-13,-124355,12.4858], Tmin=(845.474,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1036.58,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-OsOsOsH) + group(CsFFFO) + radical(ROOJ)"""),
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
    label = 'OC([C]F)OC(F)(F)F-2(2178)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {7,S} {11,S}
7  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
8  C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
9  C u0 p1 c0 {4,S} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-938.534,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,431,527,613,668,835,1199,1245,1322,617,898,1187,245.245,245.976,3091.74],'cm^-1')),
        HinderedRotor(inertia=(0.180134,'amu*angstrom^2'), symmetry=1, barrier=(7.83901,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180261,'amu*angstrom^2'), symmetry=1, barrier=(7.84225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.181369,'amu*angstrom^2'), symmetry=1, barrier=(7.85496,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.788521,'amu*angstrom^2'), symmetry=1, barrier=(34.1569,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.68999,0.0787018,-0.000117964,9.29505e-08,-2.92783e-11,-112766,29.644], Tmin=(100,'K'), Tmax=(777.163,'K')), NASAPolynomial(coeffs=[11.3058,0.0240621,-1.25032e-05,2.48272e-09,-1.76015e-13,-114416,-18.8949], Tmin=(777.163,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-938.534,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsFFFO) + group(CJ2_singlet-FCs)"""),
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
    label = '[O]OC(F)(OC(F)(F)F)C(O)F(2400)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {10,S} {12,S}
7  O u0 p2 c0 {9,S} {10,S}
8  O u0 p2 c0 {11,S} {14,S}
9  O u1 p2 c0 {7,S}
10 C u0 p0 c0 {1,S} {6,S} {7,S} {11,S}
11 C u0 p0 c0 {2,S} {8,S} {10,S} {13,S}
12 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
13 H u0 p0 c0 {11,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-1497.47,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,375,361,584,565,722,1474,261,493,600,1152,1365,1422,3097,431,527,613,668,835,1199,1245,1322,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.979765,0.121596,-0.000186714,1.36351e-07,-3.44534e-11,-179937,35.95], Tmin=(100,'K'), Tmax=(626.408,'K')), NASAPolynomial(coeffs=[15.871,0.032487,-1.76159e-05,3.51524e-09,-2.48364e-13,-182410,-40.3594], Tmin=(626.408,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1497.47,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsFFFO) + radical(ROOJ)"""),
)

species(
    label = '[O]OC(O)(F)C(F)OC(F)(F)F(2401)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {11,S} {12,S}
7  O u0 p2 c0 {10,S} {14,S}
8  O u0 p2 c0 {9,S} {10,S}
9  O u1 p2 c0 {8,S}
10 C u0 p0 c0 {1,S} {7,S} {8,S} {11,S}
11 C u0 p0 c0 {2,S} {6,S} {10,S} {13,S}
12 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
13 H u0 p0 c0 {11,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (-1497.47,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,375,361,584,565,722,1474,261,493,600,1152,1365,1422,3097,431,527,613,668,835,1199,1245,1322,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.979765,0.121596,-0.000186714,1.36351e-07,-3.44534e-11,-179937,35.95], Tmin=(100,'K'), Tmax=(626.408,'K')), NASAPolynomial(coeffs=[15.871,0.032487,-1.76159e-05,3.51524e-09,-2.48364e-13,-182410,-40.3594], Tmin=(626.408,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1497.47,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsFFFO) + radical(ROOJ)"""),
)

species(
    label = 'FOC(F)(F)F(875)',
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
    label = '[O]OC(F)=CO(2352)',
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
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,3010,987.5,1337.5,450,1655,326,540,652,719,1357],'cm^-1')),
        HinderedRotor(inertia=(0.618904,'amu*angstrom^2'), symmetry=1, barrier=(14.2298,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.620672,'amu*angstrom^2'), symmetry=1, barrier=(14.2705,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.0338,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.82782,0.0498996,-6.94358e-05,4.9215e-08,-1.37139e-11,-30875.5,19.6982], Tmin=(100,'K'), Tmax=(882.061,'K')), NASAPolynomial(coeffs=[10.0828,0.0124643,-5.77428e-06,1.09905e-09,-7.64362e-14,-32331.8,-19.0915], Tmin=(882.061,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-257.348,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + radical(ROOJ)"""),
)

species(
    label = 'OF(173)',
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
    label = '[O]OC(F)=COC(F)(F)F(2402)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {8,S} {9,S}
6  O u0 p2 c0 {7,S} {10,S}
7  O u1 p2 c0 {6,S}
8  C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
9  C u0 p0 c0 {5,S} {10,D} {11,S}
10 C u0 p0 c0 {4,S} {6,S} {9,D}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-918.795,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,431,527,613,668,835,1199,1245,1322,3010,987.5,1337.5,450,1655,326,540,652,719,1357,199.034,199.034,199.036],'cm^-1')),
        HinderedRotor(inertia=(0.321751,'amu*angstrom^2'), symmetry=1, barrier=(9.04478,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.321749,'amu*angstrom^2'), symmetry=1, barrier=(9.04479,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.321754,'amu*angstrom^2'), symmetry=1, barrier=(9.0448,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (161.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.38221,0.0850083,-0.000122301,8.87847e-08,-2.54684e-11,-110380,29.7221], Tmin=(100,'K'), Tmax=(854.573,'K')), NASAPolynomial(coeffs=[13.877,0.021842,-1.14255e-05,2.28701e-09,-1.63488e-13,-112686,-33.2615], Tmin=(854.573,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-918.795,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(CsFFFO) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(ROOJ)"""),
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
    label = '[O]C(F)(F)C(O)OC(F)(F)F(2403)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {9,S} {11,S}
7  O u0 p2 c0 {9,S} {13,S}
8  O u1 p2 c0 {10,S}
9  C u0 p0 c0 {6,S} {7,S} {10,S} {12,S}
10 C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
11 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {7,S}
"""),
    E0 = (-1476.46,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,351,323,533,609,664,892,1120,1201,431,527,613,668,835,1199,1245,1322,180,180,180,4000],'cm^-1')),
        HinderedRotor(inertia=(0.25146,'amu*angstrom^2'), symmetry=1, barrier=(5.78156,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.251706,'amu*angstrom^2'), symmetry=1, barrier=(5.78721,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.253176,'amu*angstrom^2'), symmetry=1, barrier=(5.82101,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11755,'amu*angstrom^2'), symmetry=1, barrier=(25.6948,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (181.038,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.435672,0.105686,-0.000166564,1.35427e-07,-4.37148e-11,-177425,35.1049], Tmin=(100,'K'), Tmax=(760.48,'K')), NASAPolynomial(coeffs=[14.1078,0.0291858,-1.5663e-05,3.13381e-09,-2.22418e-13,-179637,-31.0758], Tmin=(760.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1476.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCFFO) + group(CsFFFO) + radical(O2sj(Cs-CsF1sF1s))"""),
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
    label = 'OC(OC(F)(F)F)=C(F)F(2174)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {8,S} {9,S}
7  O u0 p2 c0 {9,S} {11,S}
8  C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
9  C u0 p0 c0 {6,S} {7,S} {10,D}
10 C u0 p0 c0 {4,S} {5,S} {9,D}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-1342.62,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,431,527,613,668,835,1199,1245,1322,350,440,435,1725,182,240,577,636,1210,1413,299.647,299.659,299.669],'cm^-1')),
        HinderedRotor(inertia=(0.283622,'amu*angstrom^2'), symmetry=1, barrier=(18.0778,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197407,'amu*angstrom^2'), symmetry=1, barrier=(12.5771,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197352,'amu*angstrom^2'), symmetry=1, barrier=(12.5769,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (164.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00562612,0.0910516,-0.000128833,8.79306e-08,-2.32769e-11,-161339,26.7187], Tmin=(100,'K'), Tmax=(931.326,'K')), NASAPolynomial(coeffs=[17.4298,0.0162173,-8.30754e-06,1.6571e-09,-1.18682e-13,-164585,-56.1036], Tmin=(931.326,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1342.62,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(CsFFFO) + group(Cds-CdsCsCs) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs)"""),
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
    label = '[O]O[C](F)C(O)OC(F)(F)F(2404)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {9,S} {10,S}
6  O u0 p2 c0 {9,S} {13,S}
7  O u0 p2 c0 {8,S} {11,S}
8  O u1 p2 c0 {7,S}
9  C u0 p0 c0 {5,S} {6,S} {11,S} {12,S}
10 C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
11 C u1 p0 c0 {4,S} {7,S} {9,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-1083.66,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,1380,1390,370,380,2900,435,431,527,613,668,835,1199,1245,1322,395,473,707,1436,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.824841,0.119668,-0.000217583,1.97292e-07,-6.82949e-11,-130173,38.0181], Tmin=(100,'K'), Tmax=(845.553,'K')), NASAPolynomial(coeffs=[11.9634,0.0331556,-1.7958e-05,3.51187e-09,-2.42213e-13,-131406,-16.0328], Tmin=(845.553,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1083.66,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(CsCFHO) + group(CsFFFO) + radical(ROOJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]OC(F)(F)C(O)O[C](F)F(2405)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {9,S} {11,S}
6  O u0 p2 c0 {9,S} {13,S}
7  O u0 p2 c0 {8,S} {10,S}
8  O u1 p2 c0 {7,S}
9  C u0 p0 c0 {5,S} {6,S} {10,S} {12,S}
10 C u0 p0 c0 {1,S} {2,S} {7,S} {9,S}
11 C u1 p0 c0 {3,S} {4,S} {5,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {6,S}
"""),
    E0 = (-1067.76,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,493,600,700,1144,1293,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.932215,0.122323,-0.000223364,2.0307e-07,-7.05864e-11,-128258,38.3338], Tmin=(100,'K'), Tmax=(839.487,'K')), NASAPolynomial(coeffs=[12.2504,0.0336032,-1.85478e-05,3.65679e-09,-2.53525e-13,-129558,-17.5212], Tmin=(839.487,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1067.76,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(CsCFFO) + group(CsFFHO) + radical(ROOJ) + radical(Csj(F1s)(F1s)(O2s-Cs))"""),
)

species(
    label = 'CF3O(47)',
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
    label = '[O]OC(F)(F)[CH]O(2246)',
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
        HinderedRotor(inertia=(0.294626,'amu*angstrom^2'), symmetry=1, barrier=(6.77403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.294504,'amu*angstrom^2'), symmetry=1, barrier=(6.77124,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.294738,'amu*angstrom^2'), symmetry=1, barrier=(6.7766,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.882035,0.0842618,-0.000176174,1.72331e-07,-6.1358e-11,-57144.9,24.0223], Tmin=(100,'K'), Tmax=(889.356,'K')), NASAPolynomial(coeffs=[5.46392,0.025764,-1.36051e-05,2.56429e-09,-1.70108e-13,-56461.4,10.8791], Tmin=(889.356,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-475.94,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(CsCFFO) + group(Cs-CsOsHH) + radical(ROOJ) + radical(Csj(Cs-F1sF1sO2s)(O2s-H)(H))"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05587,-0.0054938,7.30062e-05,-1.34908e-07,7.857e-11,-58113,8.13986], Tmin=(10,'K'), Tmax=(570.231,'K')), NASAPolynomial(coeffs=[3.53924,0.0118582,-8.75006e-06,2.89362e-09,-3.54041e-13,-58277.3,8.38507], Tmin=(570.231,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-483.19,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""F[C](F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]OC(F)(F)C([O])O(2406)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  O u0 p2 c0 {7,S} {10,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u1 p2 c0 {7,S}
6  O u1 p2 c0 {4,S}
7  C u0 p0 c0 {3,S} {5,S} {8,S} {9,S}
8  C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-627.831,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.182567,'amu*angstrom^2'), symmetry=1, barrier=(4.19758,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.181469,'amu*angstrom^2'), symmetry=1, barrier=(4.17234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.18316,'amu*angstrom^2'), symmetry=1, barrier=(4.21121,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (128.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.739866,0.0829089,-0.000154239,1.45526e-07,-5.21261e-11,-75404,28.1387], Tmin=(100,'K'), Tmax=(843.976,'K')), NASAPolynomial(coeffs=[7.39812,0.02803,-1.52516e-05,2.99622e-09,-2.07426e-13,-75697.2,2.06674], Tmin=(843.976,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-627.831,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(CsCFFO) + radical(CCOJ) + radical(ROOJ)"""),
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
    label = '[O]OC(F)(F)[CH]OC(F)(F)F(2407)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {10,S} {11,S}
7  O u0 p2 c0 {8,S} {9,S}
8  O u1 p2 c0 {7,S}
9  C u0 p0 c0 {1,S} {2,S} {7,S} {11,S}
10 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
11 C u1 p0 c0 {6,S} {9,S} {12,S}
12 H u0 p0 c0 {11,S}
"""),
    E0 = (-1132.23,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,253,525,597,667,842,1178,1324,431,527,613,668,835,1199,1245,1322,3025,407.5,1350,352.5,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.31637,'amu*angstrom^2'), symmetry=1, barrier=(7.27397,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.317251,'amu*angstrom^2'), symmetry=1, barrier=(7.29423,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.317519,'amu*angstrom^2'), symmetry=1, barrier=(7.30038,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.658279,'amu*angstrom^2'), symmetry=1, barrier=(15.1351,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (180.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.909235,0.117832,-0.000205828,1.76292e-07,-5.81259e-11,-136009,34.014], Tmin=(100,'K'), Tmax=(838.33,'K')), NASAPolynomial(coeffs=[15.5605,0.0247537,-1.33501e-05,2.60318e-09,-1.79305e-13,-138261,-39.501], Tmin=(838.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1132.23,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsHH) + group(CsCFFO) + group(CsFFFO) + radical(ROOJ) + radical(CCsJOCs)"""),
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
    label = '[O]OC(F)(F)C([O])OC(F)(F)F(2408)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {11,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {10,S} {12,S}
7  O u0 p2 c0 {9,S} {11,S}
8  O u1 p2 c0 {10,S}
9  O u1 p2 c0 {7,S}
10 C u0 p0 c0 {6,S} {8,S} {11,S} {13,S}
11 C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
12 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
13 H u0 p0 c0 {10,S}
"""),
    E0 = (-1287.22,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,431,527,613,668,835,1199,1245,1322,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (196.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.891652,0.122727,-0.000225543,2.09948e-07,-7.48445e-11,-154655,38.5213], Tmin=(100,'K'), Tmax=(830.634,'K')), NASAPolynomial(coeffs=[10.4983,0.0388692,-2.17251e-05,4.3201e-09,-3.01647e-13,-155546,-8.29136], Tmin=(830.634,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1287.22,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(CsCFFO) + group(CsFFFO) + radical(CCOJ) + radical(ROOJ)"""),
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
    label = 'OC(OC(F)(F)F)[C](F)F(1928)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  O u0 p2 c0 {8,S} {9,S}
7  O u0 p2 c0 {8,S} {12,S}
8  C u0 p0 c0 {6,S} {7,S} {10,S} {11,S}
9  C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
10 C u1 p0 c0 {4,S} {5,S} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-1318.9,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,431,527,613,668,835,1199,1245,1322,190,488,555,1236,1407,259.944,259.952,259.961,1400.56],'cm^-1')),
        HinderedRotor(inertia=(0.129607,'amu*angstrom^2'), symmetry=1, barrier=(6.21592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129619,'amu*angstrom^2'), symmetry=1, barrier=(6.21599,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12964,'amu*angstrom^2'), symmetry=1, barrier=(6.21604,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129631,'amu*angstrom^2'), symmetry=1, barrier=(6.21582,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (165.039,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3613.01,'J/mol'), sigma=(6.04626,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=564.34 K, Pc=37.09 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00990972,0.097211,-0.000165662,1.47174e-07,-5.11217e-11,-158493,33.0602], Tmin=(100,'K'), Tmax=(811.379,'K')), NASAPolynomial(coeffs=[10.7262,0.0304019,-1.6308e-05,3.22386e-09,-2.256e-13,-159772,-13.5641], Tmin=(811.379,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1318.9,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(CsCsFFH) + group(CsFFFO) + radical(Csj(Cs-O2sO2sH)(F1s)(F1s))"""),
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
    label = 'O[CH]OC(F)(F)F(603)',
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
        HarmonicOscillator(frequencies=([3615,1277.5,1000,431,527,613,668,835,1199,1245,1322,3025,407.5,1350,352.5,242.257,242.257,242.265],'cm^-1')),
        HinderedRotor(inertia=(0.00287244,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161892,'amu*angstrom^2'), symmetry=1, barrier=(6.74246,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.594529,'amu*angstrom^2'), symmetry=1, barrier=(24.7601,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.73696,0.0184023,0.000148993,-4.39586e-07,3.50946e-10,-107965,13.0776], Tmin=(10,'K'), Tmax=(458.766,'K')), NASAPolynomial(coeffs=[7.32095,0.028945,-2.21207e-05,7.63896e-09,-9.71498e-13,-108734,-6.21355], Tmin=(458.766,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-897.696,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), label="""O[CH]OC(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]OC(F)(F)[C](O)OC(F)(F)F(2409)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {11,S} {12,S}
7  O u0 p2 c0 {9,S} {10,S}
8  O u0 p2 c0 {12,S} {13,S}
9  O u1 p2 c0 {7,S}
10 C u0 p0 c0 {1,S} {2,S} {7,S} {12,S}
11 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
12 C u1 p0 c0 {6,S} {8,S} {10,S}
13 H u0 p0 c0 {8,S}
"""),
    E0 = (-1307.67,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,253,525,597,667,842,1178,1324,431,527,613,668,835,1199,1245,1322,360,370,350,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.19159,0.125877,-0.000222355,1.95025e-07,-6.60935e-11,-157101,39.4273], Tmin=(100,'K'), Tmax=(823.859,'K')), NASAPolynomial(coeffs=[14.9776,0.0301694,-1.67791e-05,3.32792e-09,-2.32059e-13,-159182,-31.9036], Tmin=(823.859,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1307.67,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(CsCFFO) + group(CsFFFO) + radical(ROOJ) + radical(Cs_P)"""),
)

species(
    label = '[O]OC(F)=C(O)OC(F)(F)F(2410)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {11,S}
5  O u0 p2 c0 {9,S} {10,S}
6  O u0 p2 c0 {10,S} {12,S}
7  O u0 p2 c0 {8,S} {11,S}
8  O u1 p2 c0 {7,S}
9  C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
10 C u0 p0 c0 {5,S} {6,S} {11,D}
11 C u0 p0 c0 {4,S} {7,S} {10,D}
12 H u0 p0 c0 {6,S}
"""),
    E0 = (-1075.93,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,492.5,1135,1000,431,527,613,668,835,1199,1245,1322,350,440,435,1725,326,540,652,719,1357,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.355675,'amu*angstrom^2'), symmetry=1, barrier=(8.17767,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.355741,'amu*angstrom^2'), symmetry=1, barrier=(8.17919,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.355444,'amu*angstrom^2'), symmetry=1, barrier=(8.17236,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.355657,'amu*angstrom^2'), symmetry=1, barrier=(8.17726,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (177.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.748486,0.113817,-0.000197676,1.68955e-07,-5.57012e-11,-129243,33.552], Tmin=(100,'K'), Tmax=(834.004,'K')), NASAPolynomial(coeffs=[15.1931,0.0242633,-1.30548e-05,2.55001e-09,-1.76031e-13,-131446,-37.7329], Tmin=(834.004,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1075.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(CsFFFO) + group(Cds-CdsCsCs) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(ROOJ)"""),
)

species(
    label = 'OOC(F)(F)[C](O)OC(F)(F)F(2411)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {10,S}
2  F u0 p3 c0 {10,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {11,S}
5  F u0 p3 c0 {11,S}
6  O u0 p2 c0 {11,S} {12,S}
7  O u0 p2 c0 {9,S} {10,S}
8  O u0 p2 c0 {12,S} {13,S}
9  O u0 p2 c0 {7,S} {14,S}
10 C u0 p0 c0 {1,S} {2,S} {7,S} {12,S}
11 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
12 C u1 p0 c0 {6,S} {8,S} {10,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {9,S}
"""),
    E0 = (-1459.68,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (197.037,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.10436,0.124535,-0.000195775,1.48288e-07,-4.06602e-11,-175387,38.4013], Tmin=(100,'K'), Tmax=(633.318,'K')), NASAPolynomial(coeffs=[16.1291,0.0326395,-1.82663e-05,3.69109e-09,-2.62678e-13,-177909,-39.5519], Tmin=(633.318,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1459.68,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(CsCFFO) + group(CsFFFO) + radical(Cs_P)"""),
)

species(
    label = '[O]C(OC(F)(F)F)C(F)(F)OO(2412)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {11,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {12,S}
6  O u0 p2 c0 {10,S} {12,S}
7  O u0 p2 c0 {8,S} {11,S}
8  O u0 p2 c0 {7,S} {14,S}
9  O u1 p2 c0 {10,S}
10 C u0 p0 c0 {6,S} {9,S} {11,S} {13,S}
11 C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
12 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-1439.22,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,431,527,613,668,835,1199,1245,1322,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.12783,0.126005,-0.000219883,1.99812e-07,-7.10366e-11,-172926,38.6098], Tmin=(100,'K'), Tmax=(799.547,'K')), NASAPolynomial(coeffs=[11.9606,0.0407522,-2.2848e-05,4.59257e-09,-3.24447e-13,-174387,-17.654], Tmin=(799.547,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1439.22,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(CsCFFO) + group(CsFFFO) + radical(CCOJ)"""),
)

species(
    label = 'OC(OC(F)(F)F)[C](F)OOF(2413)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {11,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {11,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {9,S}
6  O u0 p2 c0 {10,S} {11,S}
7  O u0 p2 c0 {10,S} {14,S}
8  O u0 p2 c0 {9,S} {12,S}
9  O u0 p2 c0 {5,S} {8,S}
10 C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
11 C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
12 C u1 p0 c0 {4,S} {8,S} {10,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {7,S}
"""),
    E0 = (-1140.36,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3615,1310,387.5,850,1000,277,555,632,1380,1390,370,380,2900,435,431,527,613,668,835,1199,1245,1322,395,473,707,1436,200],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.43099,0.132737,-0.000235107,2.09912e-07,-7.25954e-11,-136971,40.8156], Tmin=(100,'K'), Tmax=(817.806,'K')), NASAPolynomial(coeffs=[14.1468,0.0360395,-2.01385e-05,4.01428e-09,-2.81188e-13,-138834,-27.0129], Tmin=(817.806,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1140.36,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2sFO) + group(Cs-CsOsOsH) + group(CsCFHO) + group(CsFFFO) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'OC(O[C](F)F)C(F)(F)OOF(2414)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {11,S}
2  F u0 p3 c0 {11,S}
3  F u0 p3 c0 {12,S}
4  F u0 p3 c0 {12,S}
5  F u0 p3 c0 {9,S}
6  O u0 p2 c0 {10,S} {12,S}
7  O u0 p2 c0 {9,S} {11,S}
8  O u0 p2 c0 {10,S} {14,S}
9  O u0 p2 c0 {5,S} {7,S}
10 C u0 p0 c0 {6,S} {8,S} {11,S} {13,S}
11 C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
12 C u1 p0 c0 {3,S} {4,S} {6,S}
13 H u0 p0 c0 {10,S}
14 H u0 p0 c0 {8,S}
"""),
    E0 = (-1124.46,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,3615,1277.5,1000,277,555,632,1380,1390,370,380,2900,435,223,363,546,575,694,1179,1410,493,600,700,1144,1293,200],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.53623,0.135367,-0.000240814,2.15624e-07,-7.48819e-11,-135055,41.1236], Tmin=(100,'K'), Tmax=(811.686,'K')), NASAPolynomial(coeffs=[14.4167,0.0365173,-2.07462e-05,4.16353e-09,-2.92865e-13,-136979,-28.4059], Tmin=(811.686,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1124.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(307.635,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2sFO) + group(Cs-CsOsOsH) + group(CsCFFO) + group(CsFFHO) + radical(Csj(F1s)(F1s)(O2s-Cs))"""),
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
    E0 = (-663.547,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-340.215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-201.414,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-199.615,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-496.202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-496.202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-91.5989,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-79.9956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-531.157,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-633.211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-308.497,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-292.597,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-423.067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-408.75,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-401.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-373.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-625.259,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-391.066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-393.598,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-498.285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-706.802,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-678.166,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-354.305,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-364.555,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]OC(F)(F)C(O)OC(F)(F)F(2295)'],
    products = ['HF(38)', 'CF2O(48)', '[O]OC(F)(F)C=O(622)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.58825e+11,'s^-1'), n=0.745914, Ea=(147.103,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R!H->C_5Br1sCl1sF1sH->F1s_Ext-2C-R_7R!H->C',), comment="""Estimated from node Root_N-1R!H->C_5Br1sCl1sF1sH->F1s_Ext-2C-R_7R!H->C
Multiplied by reaction path degeneracy 3.0"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CF2(42)', '[O]OC(O)OC(F)(F)F(2399)'],
    products = ['[O]OC(F)(F)C(O)OC(F)(F)F(2295)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.40908e-07,'m^3/(mol*s)'), n=3.45227, Ea=(197.805,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s',), comment="""Estimated from node CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]OF(124)', 'OC([C]F)OC(F)(F)F-2(2178)'],
    products = ['[O]OC(F)(F)C(O)OC(F)(F)F(2295)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.13916e+09,'m^3/(mol*s)'), n=-1.00842, Ea=(21.2535,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17406617270594374, var=16.167420070870673, Tref=1000.0, N=4, data_mean=0.0, correlation='OHOY',), comment="""Estimated from node OHOY"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CF2(42)', '[O]OC(F)(F)C(O)OF(2315)'],
    products = ['[O]OC(F)(F)C(O)OC(F)(F)F(2295)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.13916e+09,'m^3/(mol*s)'), n=-1.00842, Ea=(20.8544,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17406617270594374, var=16.167420070870673, Tref=1000.0, N=4, data_mean=0.0, correlation='OHOY',), comment="""Estimated from node OHOY"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]OC(F)(OC(F)(F)F)C(O)F(2400)'],
    products = ['[O]OC(F)(F)C(O)OC(F)(F)F(2295)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]OC(O)(F)C(F)OC(F)(F)F(2401)'],
    products = ['[O]OC(F)(F)C(O)OC(F)(F)F(2295)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction7',
    reactants = ['FOC(F)(F)F(875)', '[O]OC(F)=CO(2352)'],
    products = ['[O]OC(F)(F)C(O)OC(F)(F)F(2295)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.000170623,'m^3/(mol*s)'), n=2.84144, Ea=(220.872,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OR] for rate rule [Cd/disub_Cd/monosub;H_OR]
Euclidian distance = 1.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction8',
    reactants = ['OF(173)', '[O]OC(F)=COC(F)(F)F(2402)'],
    products = ['[O]OC(F)(F)C(O)OC(F)(F)F(2295)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(0.000286045,'m^3/(mol*s)'), n=2.81783, Ea=(231.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OH] for rate rule [Cd/disub_Cd/monosub;H_OH]
Euclidian distance = 1.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction9',
    reactants = ['O(6)', '[O]C(F)(F)C(O)OC(F)(F)F(2403)'],
    products = ['[O]OC(F)(F)C(O)OC(F)(F)F(2295)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]OC(F)(F)C(O)OC(F)(F)F(2295)'],
    products = ['HO2(10)', 'OC(OC(F)(F)F)=C(F)F(2174)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.00545e+13,'s^-1'), n=-0.00666667, Ea=(177.438,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO_NdNd]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', '[O]O[C](F)C(O)OC(F)(F)F(2404)'],
    products = ['[O]OC(F)(F)C(O)OC(F)(F)F(2295)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', '[O]OC(F)(F)C(O)O[C](F)F(2405)'],
    products = ['[O]OC(F)(F)C(O)OC(F)(F)F(2295)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CF3O(47)', '[O]OC(F)(F)[CH]O(2246)'],
    products = ['[O]OC(F)(F)C(O)OC(F)(F)F(2295)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['CF3(44)', '[O]OC(F)(F)C([O])O(2406)'],
    products = ['[O]OC(F)(F)C(O)OC(F)(F)F(2295)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.505e+06,'m^3/(mol*s)'), n=-1.3397e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Sp-4R!H-2C_Ext-2C-R_N-5R!H->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Sp-4R!H-2C_Ext-2C-R_N-5R!H->C"""),
)

reaction(
    label = 'reaction15',
    reactants = ['OH(4)', '[O]OC(F)(F)[CH]OC(F)(F)F(2407)'],
    products = ['[O]OC(F)(F)C(O)OC(F)(F)F(2295)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(5)', '[O]OC(F)(F)C([O])OC(F)(F)F(2408)'],
    products = ['[O]OC(F)(F)C(O)OC(F)(F)F(2295)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O2(2)', 'OC(OC(F)(F)F)[C](F)F(1928)'],
    products = ['[O]OC(F)(F)C(O)OC(F)(F)F(2295)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.6844e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]O[C](F)F(171)', 'O[CH]OC(F)(F)F(603)'],
    products = ['[O]OC(F)(F)C(O)OC(F)(F)F(2295)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(5)', '[O]OC(F)(F)[C](O)OC(F)(F)F(2409)'],
    products = ['[O]OC(F)(F)C(O)OC(F)(F)F(2295)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(43950,'m^3/(mol*s)'), n=1, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction20',
    reactants = ['HF(38)', '[O]OC(F)=C(O)OC(F)(F)F(2410)'],
    products = ['[O]OC(F)(F)C(O)OC(F)(F)F(2295)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.14111,'m^3/(mol*s)'), n=1.29695, Ea=(156.49,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]OC(F)(F)C(O)OC(F)(F)F(2295)'],
    products = ['OOC(F)(F)[C](O)OC(F)(F)F(2411)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(5.29e+09,'s^-1'), n=0.75, Ea=(103.847,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 272 used for R4H_SSS_OCs;O_rad_out;Cs_H_out_NDMustO
Exact match found for rate rule [R4H_SSS_OCs;O_rad_out;Cs_H_out_NDMustO]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O]C(OC(F)(F)F)C(F)(F)OO(2412)'],
    products = ['[O]OC(F)(F)C(O)OC(F)(F)F(2295)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.37227e+06,'s^-1'), n=1.56745, Ea=(58.7826,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;O_rad_out;XH_out] for rate rule [R5H_SSSS;O_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['OC(OC(F)(F)F)[C](F)OOF(2413)'],
    products = ['[O]OC(F)(F)C(O)OC(F)(F)F(2295)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(83.7872,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction24',
    reactants = ['OC(O[C](F)F)C(F)(F)OOF(2414)'],
    products = ['[O]OC(F)(F)C(O)OC(F)(F)F(2295)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(57.6368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

network(
    label = 'PDepNetwork #907',
    isomers = [
        '[O]OC(F)(F)C(O)OC(F)(F)F(2295)',
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
    label = 'PDepNetwork #907',
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

