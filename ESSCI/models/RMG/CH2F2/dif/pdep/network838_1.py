species(
    label = '[O]C(O)(F)C(F)(F)OF(1911)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {6,S}
5  O u0 p2 c0 {8,S} {10,S}
6  O u0 p2 c0 {4,S} {9,S}
7  O u1 p2 c0 {8,S}
8  C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
9  C u0 p0 c0 {2,S} {3,S} {6,S} {8,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-897.905,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,373,336,502,477,708,1220,1477,223,363,546,575,694,1179,1410,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.505746,'amu*angstrom^2'), symmetry=1, barrier=(11.6281,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.26897,'amu*angstrom^2'), symmetry=1, barrier=(29.1761,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.507306,'amu*angstrom^2'), symmetry=1, barrier=(11.664,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (149.021,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.273524,0.0906263,-0.000145211,1.09145e-07,-2.91319e-11,-107867,26.6765], Tmin=(100,'K'), Tmax=(639.097,'K')), NASAPolynomial(coeffs=[13.6731,0.0201645,-1.1292e-05,2.2656e-09,-1.59991e-13,-109853,-34.1111], Tmin=(639.097,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-897.905,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2sCF) + group(CsCFOO) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + radical(O2sj(Cs-F1sO2sCs))"""),
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
    label = '[O]C(=O)F(136)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 O u1 p2 c0 {4,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
"""),
    E0 = (-381.43,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(62.9882,'amu')),
        NonlinearRotor(inertia=([36.4335,44.7291,81.1626],'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([509.533,537.43,765.261,1011.27,1202.3,1550.29],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (63.0078,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3607.92,'J/mol'), sigma=(5.19854,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=563.55 K, Pc=58.27 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05952,-0.0056915,7.13967e-05,-1.29651e-07,7.44832e-11,-45874.2,8.21459], Tmin=(10,'K'), Tmax=(575.482,'K')), NASAPolynomial(coeffs=[3.42951,0.0118543,-8.65614e-06,2.84335e-09,-3.46285e-13,-46019.7,9.01155], Tmin=(575.482,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-381.43,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""[O]C(DO)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]C(O)(F)OF(1944)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {6,S} {7,S}
4 O u0 p2 c0 {2,S} {6,S}
5 O u1 p2 c0 {6,S}
6 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (-514.261,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,1380,1390,370,380,2900,435,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.22494,'amu*angstrom^2'), symmetry=1, barrier=(51.1557,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06143,'amu*angstrom^2'), symmetry=1, barrier=(51.1411,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.0136,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93745,0.0611311,-0.000138632,1.52053e-07,-5.99323e-11,-61792.3,17.696], Tmin=(100,'K'), Tmax=(848.176,'K')), NASAPolynomial(coeffs=[-0.23913,0.0324218,-1.89341e-05,3.79472e-09,-2.64715e-13,-60021.2,36.1025], Tmin=(848.176,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-514.261,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(O2s-CsH) + group(CsFOOO) + radical(OCOJ)"""),
)

species(
    label = '[O]C(F)(F)C(O)(F)OF(1945)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u0 p2 c0 {8,S} {10,S}
7  O u1 p2 c0 {9,S}
8  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
9  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-896.611,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,375,361,584,565,722,1474,351,323,533,609,664,892,1120,1201,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.617923,'amu*angstrom^2'), symmetry=1, barrier=(14.2073,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.61919,'amu*angstrom^2'), symmetry=1, barrier=(14.2364,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.61863,'amu*angstrom^2'), symmetry=1, barrier=(14.2235,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (149.021,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0568127,0.0939705,-0.000154858,1.23044e-07,-3.74282e-11,-107702,27.2072], Tmin=(100,'K'), Tmax=(704.006,'K')), NASAPolynomial(coeffs=[14.802,0.0181273,-1.01694e-05,2.04034e-09,-1.44061e-13,-109975,-40.1512], Tmin=(704.006,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-896.611,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(O2s-CsH) + group(CsCFOO) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + radical(O2sj(Cs-CsF1sF1s))"""),
)

species(
    label = '[O]C(O)(OF)C(F)(F)F(1946)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {5,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u0 p2 c0 {8,S} {10,S}
7  O u1 p2 c0 {8,S}
8  C u0 p0 c0 {5,S} {6,S} {7,S} {9,S}
9  C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-922.335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (149.021,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.451049,0.0905202,-0.000171551,1.6094e-07,-5.6823e-11,-110815,28.131], Tmin=(100,'K'), Tmax=(859.016,'K')), NASAPolynomial(coeffs=[8.12461,0.0278289,-1.50046e-05,2.91239e-09,-1.99314e-13,-111139,-1.93368], Tmin=(859.016,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-922.335,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsOs) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-CsOs) + radical(CCOJ)"""),
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
    label = 'O[C](F)C(F)(F)OF(510)',
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
        HinderedRotor(inertia=(0.488313,'amu*angstrom^2'), symmetry=1, barrier=(11.2273,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.48922,'amu*angstrom^2'), symmetry=1, barrier=(11.2481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.487627,'amu*angstrom^2'), symmetry=1, barrier=(11.2115,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (133.022,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.120369,0.110758,-0.000240512,2.31766e-07,-8.07929e-11,-89510.9,24.1495], Tmin=(100,'K'), Tmax=(898.754,'K')), NASAPolynomial(coeffs=[8.5321,0.0233376,-1.29778e-05,2.43594e-09,-1.58548e-13,-89090.7,-5.68019], Tmin=(898.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-745.312,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sF1sO2s)(F1s)(O2s-H))"""),
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
    label = 'O=C(O)C(F)(F)OF(330)',
    structure = adjacencyList("""1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {4,S}
4 O u0 p2 c0 {3,S} {7,S}
5 O u0 p2 c0 {8,S} {9,S}
6 O u0 p2 c0 {8,D}
7 C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
8 C u0 p0 c0 {5,S} {6,D} {7,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-905.812,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,251,367,519,700,855,1175,1303,230.718,230.721,230.726,230.727,230.73,1022.46],'cm^-1')),
        HinderedRotor(inertia=(0.0100482,'amu*angstrom^2'), symmetry=1, barrier=(7.45435,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0666329,'amu*angstrom^2'), symmetry=1, barrier=(49.4327,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30858,'amu*angstrom^2'), symmetry=1, barrier=(49.4327,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.023,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3913.64,'J/mol'), sigma=(5.83574,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=611.30 K, Pc=44.68 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31261,0.0639673,-8.60146e-05,5.94352e-08,-1.65004e-11,-108851,21.9003], Tmin=(100,'K'), Tmax=(875.403,'K')), NASAPolynomial(coeffs=[10.9088,0.0201199,-1.08833e-05,2.21948e-09,-1.60792e-13,-110532,-23.1188], Tmin=(875.403,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-905.812,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-O2d)H) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cds-OdCsOs)"""),
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
    label = 'O=C(F)C(F)(F)OF(516)',
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
    label = 'FO[C](F)F(240)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66913,0.0316873,-4.80465e-05,3.64256e-08,-1.08807e-11,-38142.7,14.3743], Tmin=(100,'K'), Tmax=(821.842,'K')), NASAPolynomial(coeffs=[7.60338,0.00767066,-4.20997e-06,8.64394e-10,-6.26021e-14,-38953.7,-8.46221], Tmin=(821.842,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-317.517,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsF1sF1sO2s)"""),
)

species(
    label = 'O=C(O)F(134)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (-625.647,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,482,664,788,1296,1923],'cm^-1')),
        HinderedRotor(inertia=(1.70004,'amu*angstrom^2'), symmetry=1, barrier=(39.0874,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (64.0158,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3607.92,'J/mol'), sigma=(5.19854,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=563.55 K, Pc=58.27 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.96655,0.00185137,4.47909e-05,-7.96564e-08,4.23622e-11,-75247.4,7.80392], Tmin=(10,'K'), Tmax=(614.343,'K')), NASAPolynomial(coeffs=[2.83475,0.0173246,-1.27762e-05,4.2862e-09,-5.35304e-13,-75261.2,11.4681], Tmin=(614.343,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-625.647,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""ODC(O)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O][C](O)C(F)(F)OF(1956)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {4,S}
4 O u0 p2 c0 {3,S} {7,S}
5 O u0 p2 c0 {8,S} {9,S}
6 O u1 p2 c0 {8,S}
7 C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
8 C u1 p0 c0 {5,S} {6,S} {7,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-511.586,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,253,525,597,667,842,1178,1324,360,370,350,180,180,2452.42],'cm^-1')),
        HinderedRotor(inertia=(0.270005,'amu*angstrom^2'), symmetry=1, barrier=(6.20794,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.269746,'amu*angstrom^2'), symmetry=1, barrier=(6.202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.37189,'amu*angstrom^2'), symmetry=1, barrier=(31.5424,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.023,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.711237,0.082688,-0.000155626,1.45125e-07,-5.15821e-11,-61421.2,27.5793], Tmin=(100,'K'), Tmax=(832.406,'K')), NASAPolynomial(coeffs=[8.81151,0.0238542,-1.37307e-05,2.74949e-09,-1.92301e-13,-62080,-5.87088], Tmin=(832.406,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-511.586,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2sCF) + group(Cs-CsOsOsH) + group(CsCFFO) + radical(CCOJ) + radical(Cs_P)"""),
)

species(
    label = '[O]C(O)(F)[C](F)OF(1957)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {8,S}
3 F u0 p3 c0 {5,S}
4 O u0 p2 c0 {7,S} {9,S}
5 O u0 p2 c0 {3,S} {8,S}
6 O u1 p2 c0 {7,S}
7 C u0 p0 c0 {1,S} {4,S} {6,S} {8,S}
8 C u1 p0 c0 {2,S} {5,S} {7,S}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (-473.863,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,373,336,502,477,708,1220,1477,395,473,707,1436,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.429027,'amu*angstrom^2'), symmetry=1, barrier=(9.86417,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19403,'amu*angstrom^2'), symmetry=1, barrier=(27.4532,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.429281,'amu*angstrom^2'), symmetry=1, barrier=(9.87002,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.023,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.440422,0.0885115,-0.000166234,1.50791e-07,-5.18146e-11,-56874.1,26.5415], Tmin=(100,'K'), Tmax=(852.322,'K')), NASAPolynomial(coeffs=[10.5987,0.0207726,-1.17073e-05,2.30186e-09,-1.58434e-13,-57876.9,-16.5673], Tmin=(852.322,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-473.863,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2sCF) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(O2sj(Cs-F1sO2sCs)) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]C(O)(F)C([O])(F)F(1958)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {8,S}
3 F u0 p3 c0 {8,S}
4 O u0 p2 c0 {7,S} {9,S}
5 O u1 p2 c0 {7,S}
6 O u1 p2 c0 {8,S}
7 C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
8 C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (-772.445,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,373,336,502,477,708,1220,1477,351,323,533,609,664,892,1120,1201,180],'cm^-1')),
        HinderedRotor(inertia=(0.399192,'amu*angstrom^2'), symmetry=1, barrier=(9.17822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.399458,'amu*angstrom^2'), symmetry=1, barrier=(9.18433,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (130.023,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.762758,0.0772105,-0.000130873,1.10807e-07,-3.65169e-11,-92792.8,25.3063], Tmin=(100,'K'), Tmax=(812.786,'K')), NASAPolynomial(coeffs=[11.6211,0.0175941,-9.4486e-06,1.8588e-09,-1.29415e-13,-94353.9,-23.5728], Tmin=(812.786,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-772.445,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(CsCFOO) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + radical(O2sj(Cs-F1sO2sCs)) + radical(O2sj(Cs-CsF1sF1s))"""),
)

species(
    label = '[O][C](F)C(F)(F)OF(282)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {8,S}
4 F u0 p3 c0 {5,S}
5 O u0 p2 c0 {4,S} {7,S}
6 O u1 p2 c0 {8,S}
7 C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8 C u1 p0 c0 {3,S} {6,S} {7,S}
"""),
    E0 = (-518.741,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,253,525,597,667,842,1178,1324,395,473,707,1436,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.358612,'amu*angstrom^2'), symmetry=1, barrier=(8.24519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.362581,'amu*angstrom^2'), symmetry=1, barrier=(8.33646,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.014,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.124802,0.106886,-0.000241496,2.36733e-07,-8.31204e-11,-62270.9,24.1681], Tmin=(100,'K'), Tmax=(903.345,'K')), NASAPolynomial(coeffs=[7.16295,0.0223503,-1.25027e-05,2.33328e-09,-1.50288e-13,-61364.9,2.98162], Tmin=(903.345,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-518.741,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sH)) + radical(Csj(Cs-F1sF1sO2s)(F1s)(O2s-H))"""),
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
    label = '[O]C([O])(F)C(F)(F)OF(1959)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {8,S}
2 F u0 p3 c0 {8,S}
3 F u0 p3 c0 {9,S}
4 F u0 p3 c0 {5,S}
5 O u0 p2 c0 {4,S} {8,S}
6 O u1 p2 c0 {9,S}
7 O u1 p2 c0 {9,S}
8 C u0 p0 c0 {1,S} {2,S} {5,S} {9,S}
9 C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
"""),
    E0 = (-639.232,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,223,363,546,575,694,1179,1410,373,336,502,477,708,1220,1477,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0714638,'amu*angstrom^2'), symmetry=1, barrier=(1.64309,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05186,'amu*angstrom^2'), symmetry=1, barrier=(24.1844,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (148.013,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.36779,0.0881993,-0.000157128,1.38538e-07,-4.72513e-11,-76759,26.4281], Tmin=(100,'K'), Tmax=(816.098,'K')), NASAPolynomial(coeffs=[11.7098,0.0209142,-1.19644e-05,2.39688e-09,-1.68112e-13,-78220.8,-23.6], Tmin=(816.098,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-639.232,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2sCF) + group(CsCFOO) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + radical(O2sj(Cs-F1sO2sCs)) + radical(O2sj(Cs-F1sO2sCs))"""),
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
    label = '[O]C(O)(F)[C](F)F(1960)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {7,S}
4 O u0 p2 c0 {6,S} {8,S}
5 O u1 p2 c0 {6,S}
6 C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
7 C u1 p0 c0 {2,S} {3,S} {6,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-632.104,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,373,336,502,477,708,1220,1477,190,488,555,1236,1407,180],'cm^-1')),
        HinderedRotor(inertia=(0.686752,'amu*angstrom^2'), symmetry=1, barrier=(15.7898,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.670113,'amu*angstrom^2'), symmetry=1, barrier=(15.4072,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.023,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06239,0.073835,-0.000140038,1.29452e-07,-4.52337e-11,-75927.7,21.567], Tmin=(100,'K'), Tmax=(848.314,'K')), NASAPolynomial(coeffs=[8.75717,0.0190976,-1.0619e-05,2.10023e-09,-1.45403e-13,-76569.2,-10.3763], Tmin=(848.314,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-632.104,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-F1sO2sCs)) + radical(Csj(Cs-F1sO2sO2s)(F1s)(F1s))"""),
)

species(
    label = '[O][C](O)F(137)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u1 p2 c0 {4,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (-194.363,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,482,586,761,1411,442.957],'cm^-1')),
        HinderedRotor(inertia=(0.000859644,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (64.0158,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.30412,0.00866747,1.59293e-05,-2.96311e-08,1.26046e-11,-23345.4,14.3614], Tmin=(100,'K'), Tmax=(953.516,'K')), NASAPolynomial(coeffs=[8.84221,0.00217781,-1.99934e-07,6.07776e-11,-8.37542e-15,-25162.6,-16.0842], Tmin=(953.516,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-194.363,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-F1sO2sH)) + radical(CsF1sO2sO2s)"""),
)

species(
    label = '[O]C(=O)C(F)(F)OF(456)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {4,S}
4 O u0 p2 c0 {3,S} {7,S}
5 O u1 p2 c0 {8,S}
6 O u0 p2 c0 {8,D}
7 C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
8 C u0 p0 c0 {5,S} {6,D} {7,S}
"""),
    E0 = (-680.107,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,251,367,519,700,855,1175,1303,180,180,180,820.564,820.672,820.822,820.866],'cm^-1')),
        HinderedRotor(inertia=(0.0087674,'amu*angstrom^2'), symmetry=1, barrier=(4.18797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.33997,'amu*angstrom^2'), symmetry=1, barrier=(53.8006,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (129.015,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39599,0.0648293,-0.000112617,1.05948e-07,-3.92746e-11,-81711.5,21.7481], Tmin=(100,'K'), Tmax=(786.454,'K')), NASAPolynomial(coeffs=[6.66161,0.0257411,-1.45921e-05,2.95545e-09,-2.10167e-13,-82159.2,0.0291418], Tmin=(786.454,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-680.107,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-O2d)H) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cds-OdCsOs) + radical(CCOJ)"""),
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
    label = '[O]C(O)(F)C(=O)F(1961)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {6,S} {8,S}
4 O u1 p2 c0 {6,S}
5 O u0 p2 c0 {7,D}
6 C u0 p0 c0 {1,S} {3,S} {4,S} {7,S}
7 C u0 p0 c0 {2,S} {5,D} {6,S}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (-770.523,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,486,617,768,1157,1926,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.854338,'amu*angstrom^2'), symmetry=1, barrier=(19.6429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.854027,'amu*angstrom^2'), symmetry=1, barrier=(19.6358,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.024,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32753,0.0652469,-0.000115715,1.01916e-07,-3.43502e-11,-92582.5,22.1915], Tmin=(100,'K'), Tmax=(851.808,'K')), NASAPolynomial(coeffs=[9.22452,0.0162738,-8.53784e-06,1.6471e-09,-1.12647e-13,-93496.5,-12.1084], Tmin=(851.808,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-770.523,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(COCsFO) + radical(C=OCOJ)"""),
)

species(
    label = 'O=C(O)[C](F)OF(452)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {3,S}
3 O u0 p2 c0 {2,S} {6,S}
4 O u0 p2 c0 {7,S} {8,S}
5 O u0 p2 c0 {7,D}
6 C u1 p0 c0 {1,S} {3,S} {7,S}
7 C u0 p0 c0 {4,S} {5,D} {6,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-543.762,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,280,501,1494,1531,180,180,180,580.479,924.219,925.434],'cm^-1')),
        HinderedRotor(inertia=(0.0913279,'amu*angstrom^2'), symmetry=1, barrier=(55.4254,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.41065,'amu*angstrom^2'), symmetry=1, barrier=(55.4256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0913443,'amu*angstrom^2'), symmetry=1, barrier=(55.4053,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.024,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06014,0.0471009,-5.77729e-05,4.01711e-08,-1.17946e-11,-65333.5,18.8253], Tmin=(100,'K'), Tmax=(813.352,'K')), NASAPolynomial(coeffs=[7.1478,0.0220803,-1.16296e-05,2.34976e-09,-1.69497e-13,-66161.1,-4.66866], Tmin=(813.352,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-543.762,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-O2d)H) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsOs) + radical(CsCOF1sO2s)"""),
)

species(
    label = 'O[C](OF)C(F)(F)OF(1963)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {5,S}
4  F u0 p3 c0 {7,S}
5  O u0 p2 c0 {3,S} {8,S}
6  O u0 p2 c0 {9,S} {10,S}
7  O u0 p2 c0 {4,S} {9,S}
8  C u0 p0 c0 {1,S} {2,S} {5,S} {9,S}
9  C u1 p0 c0 {6,S} {7,S} {8,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-602.783,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([504,610,1067,1155,3615,1277.5,1000,253,525,597,667,842,1178,1324,360,370,350,252.695,260.095,261.376],'cm^-1')),
        HinderedRotor(inertia=(0.15499,'amu*angstrom^2'), symmetry=1, barrier=(6.92148,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.310833,'amu*angstrom^2'), symmetry=1, barrier=(15.5958,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143099,'amu*angstrom^2'), symmetry=1, barrier=(6.90504,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.598248,'amu*angstrom^2'), symmetry=1, barrier=(26.0705,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (149.021,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0883928,0.0947814,-0.000156774,1.22508e-07,-3.54073e-11,-72365.5,29.153], Tmin=(100,'K'), Tmax=(654.814,'K')), NASAPolynomial(coeffs=[14.3998,0.0193191,-1.1308e-05,2.30229e-09,-1.63782e-13,-74496.1,-35.7895], Tmin=(654.814,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-602.783,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(O2sCF) + group(Cs-CsOsOsH) + group(CsCFFO) + radical(Cs_P)"""),
)

species(
    label = 'OC(F)(OF)[C](F)OF(1964)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {7,S}
5  O u0 p2 c0 {8,S} {10,S}
6  O u0 p2 c0 {3,S} {8,S}
7  O u0 p2 c0 {4,S} {9,S}
8  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
9  C u1 p0 c0 {2,S} {7,S} {8,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-598.028,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,504,610,1067,1155,364,474,579,619,655,1345,395,473,707,1436,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.621824,'amu*angstrom^2'), symmetry=1, barrier=(14.297,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.621078,'amu*angstrom^2'), symmetry=1, barrier=(14.2798,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.67635,'amu*angstrom^2'), symmetry=1, barrier=(38.5425,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.621551,'amu*angstrom^2'), symmetry=1, barrier=(14.2907,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (149.021,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.3813,0.10685,-0.000196921,1.73854e-07,-5.85689e-11,-71778.5,28.8451], Tmin=(100,'K'), Tmax=(840.525,'K')), NASAPolynomial(coeffs=[13.9666,0.0209607,-1.22174e-05,2.43155e-09,-1.68648e-13,-73568.4,-34.1825], Tmin=(840.525,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-598.028,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(O2sCF) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + radical(CsCsF1sO2s)"""),
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
    E0 = (-493.855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-94.3249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-209.485,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-210.779,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-114.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-399.109,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-400.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-493.208,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-50.5688,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-12.8455,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-311.428,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-102.221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-39.3012,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-141.293,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-123.754,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-370.613,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-301.005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-135.864,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-132.071,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-104.855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(O)(F)C(F)(F)OF(1911)'],
    products = ['HF(38)', 'CF2O(48)', '[O]C(=O)F(136)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(0.00285451,'s^-1'), n=4.62568, Ea=(15.9241,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CF2(42)', '[O]C(O)(F)OF(1944)'],
    products = ['[O]C(O)(F)C(F)(F)OF(1911)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.40908e-07,'m^3/(mol*s)'), n=3.45227, Ea=(235.523,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s',), comment="""Estimated from node CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]C(F)(F)C(O)(F)OF(1945)'],
    products = ['[O]C(O)(F)C(F)(F)OF(1911)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]C(O)(F)C(F)(F)OF(1911)'],
    products = ['[O]C(O)(OF)C(F)(F)F(1946)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(6)', 'O[C](F)C(F)(F)OF(510)'],
    products = ['[O]C(O)(F)C(F)(F)OF(1911)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_ter_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['F(37)', 'O=C(O)C(F)(F)OF(330)'],
    products = ['[O]C(O)(F)C(F)(F)OF(1911)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4.08261e+20,'m^3/(mol*s)'), n=-5.07836, Ea=(45.6856,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O"""),
)

reaction(
    label = 'reaction7',
    reactants = ['OH(4)', 'O=C(F)C(F)(F)OF(516)'],
    products = ['[O]C(O)(F)C(F)(F)OF(1911)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.08261e+20,'m^3/(mol*s)'), n=-5.07836, Ea=(44.5291,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O"""),
)

reaction(
    label = 'reaction8',
    reactants = ['FO[C](F)F(240)', 'O=C(O)F(134)'],
    products = ['[O]C(O)(F)C(F)(F)OF(1911)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.4291e-06,'m^3/(mol*s)'), n=3.20779, Ea=(61.8307,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing_Ext-1R!H-R_N-2R!H->C_N-5R!H-u1',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing_Ext-1R!H-R_N-2R!H->C_N-5R!H-u1"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F(37)', '[O][C](O)C(F)(F)OF(1956)'],
    products = ['[O]C(O)(F)C(F)(F)OF(1911)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.55564e+07,'m^3/(mol*s)'), n=-5.63145e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.7121440562946592, var=2.9508506800589083, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F(37)', '[O]C(O)(F)[C](F)OF(1957)'],
    products = ['[O]C(O)(F)C(F)(F)OF(1911)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', '[O]C(O)(F)C([O])(F)F(1958)'],
    products = ['[O]C(O)(F)C(F)(F)OF(1911)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.81e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_N-2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_N-2R->C"""),
)

reaction(
    label = 'reaction12',
    reactants = ['OH(4)', '[O][C](F)C(F)(F)OF(282)'],
    products = ['[O]C(O)(F)C(F)(F)OF(1911)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C_N-2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C_N-2R->C"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(5)', '[O]C([O])(F)C(F)(F)OF(1959)'],
    products = ['[O]C(O)(F)C(F)(F)OF(1911)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.42074e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]F(228)', '[O]C(O)(F)[C](F)F(1960)'],
    products = ['[O]C(O)(F)C(F)(F)OF(1911)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['FO[C](F)F(240)', '[O][C](O)F(137)'],
    products = ['[O]C(O)(F)C(F)(F)OF(1911)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['HF(38)', '[O]C(=O)C(F)(F)OF(456)'],
    products = ['[O]C(O)(F)C(F)(F)OF(1911)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.218312,'m^3/(mol*s)'), n=1.86531, Ea=(202.481,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd',), comment="""Estimated from node HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F2(77)', '[O]C(O)(F)C(=O)F(1961)'],
    products = ['[O]C(O)(F)C(F)(F)OF(1911)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(90.197,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction18',
    reactants = ['F2(77)', 'O=C(O)[C](F)OF(452)'],
    products = ['[O]C(O)(F)C(F)(F)OF(1911)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(28.577,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction19',
    reactants = ['O[C](OF)C(F)(F)OF(1963)'],
    products = ['[O]C(O)(F)C(F)(F)OF(1911)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(82.5868,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction20',
    reactants = ['OC(F)(OF)[C](F)OF(1964)'],
    products = ['[O]C(O)(F)C(F)(F)OF(1911)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(105.048,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

network(
    label = 'PDepNetwork #838',
    isomers = [
        '[O]C(O)(F)C(F)(F)OF(1911)',
    ],
    reactants = [
        ('HF(38)', 'CF2O(48)', '[O]C(=O)F(136)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #838',
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

