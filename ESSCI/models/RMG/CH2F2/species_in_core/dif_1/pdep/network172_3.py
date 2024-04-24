species(
    label = '[O]O[C](F)F(147)',
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
    label = 'FC1(F)OO1(932)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {4,S} {5,S}
4 O u0 p2 c0 {3,S} {5,S}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {4,S}
"""),
    E0 = (-474.414,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([180,180,761.363,1045.57,1045.77,1045.82,1045.97,1046.15,1046.16],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0062,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2976.61,'J/mol'), sigma=(5.20996,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=464.94 K, Pc=47.76 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.95302,0.00256699,5.22953e-05,-9.50593e-08,5.05491e-11,-57057.6,8.08501], Tmin=(10,'K'), Tmax=(633.13,'K')), NASAPolynomial(coeffs=[3.1996,0.0191329,-1.49226e-05,5.17103e-09,-6.57667e-13,-57198.8,9.5068], Tmin=(633.13,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-474.414,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""FC1(F)OO1""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]C([O])(F)F(152)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u1 p2 c0 {5,S}
4 O u1 p2 c0 {5,S}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {4,S}
"""),
    E0 = (-372.098,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,486,586,623,710,938,1161,1294,464.782],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.0062,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3646.55,'J/mol'), sigma=(6.01309,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=569.58 K, Pc=38.06 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.92149,0.00463651,6.40183e-05,-1.39109e-07,8.57731e-11,-44748,10.2194], Tmin=(10,'K'), Tmax=(572.951,'K')), NASAPolynomial(coeffs=[4.84692,0.0160279,-1.25421e-05,4.35666e-09,-5.55356e-13,-45147.1,3.71309], Tmin=(572.951,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-372.098,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""[O]C([O])(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'CF2O(48)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
"""),
    E0 = (-650.702,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17614,0.021087,-2.39216e-05,1.39949e-08,-3.37605e-12,-77951.3,19.1134], Tmin=(298,'K'), Tmax=(1150,'K')), NASAPolynomial(coeffs=[6.79131,0.00338305,-1.43012e-06,2.71226e-10,-1.90524e-14,-79454,-9.48312], Tmin=(1150,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-650.702,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(83.1447,'J/mol/K'), label="""CF2O""", comment="""Thermo library: Fluorine"""),
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
    label = 'CF2(42)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u2 p0 c0 {1,S} {2,S}
"""),
    E0 = (33.7272,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([700.388,993.445,1307.31],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (50.0074,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2405.34,'J/mol'), sigma=(4.22615,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=375.71 K, Pc=72.31 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28591,0.0107608,-1.05382e-05,4.89881e-09,-8.86384e-13,4216.69,13.1348], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[5.33121,0.00197748,-9.60248e-07,2.10704e-10,-1.5954e-14,3366.49,-2.56367], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(33.7272,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CF2(T)""", comment="""Thermo library: halogens"""),
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
    label = '[O]C(=O)F(143)',
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
    label = '[O][C](F)F(178)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u1 p2 c0 {4,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
"""),
    E0 = (-255.477,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([493,600,700,1144,1293,180],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.0068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.105,0.0167699,-1.53787e-05,4.83895e-09,3.91619e-15,-30692.1,11.7073], Tmin=(100,'K'), Tmax=(1048.23,'K')), NASAPolynomial(coeffs=[8.58998,0.00108971,-4.53613e-07,1.24874e-10,-1.13715e-14,-32130.4,-16.3887], Tmin=(1048.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-255.477,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsOJ) + radical(CsF1sF1sO2s)"""),
)

species(
    label = 'CF2(89)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p1 c0 {1,S} {2,S}
"""),
    E0 = (-196.899,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([192,594,627],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (50.0074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28591,0.0107608,-1.05382e-05,4.89881e-09,-8.86384e-13,-23521.2,13.1348], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[5.33121,0.00197748,-9.60248e-07,2.10704e-10,-1.5954e-14,-24371.4,-2.56367], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-196.899,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CF2""", comment="""Thermo library: Fluorine"""),
)

species(
    label = 'O(172)',
    structure = adjacencyList("""1 O u0 p3 c0
"""),
    E0 = (432.331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (15.9994,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.05029e-14,-1.20495e-16,5.23432e-20,-6.99814e-24,51997.4,2.99252], Tmin=(100,'K'), Tmax=(4230.24,'K')), NASAPolynomial(coeffs=[2.5,3.58763e-09,-1.22015e-12,1.84118e-16,-1.04e-20,51997.4,2.99254], Tmin=(4230.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(432.331,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'F[C]1OO1(174)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {4,S}
3 O u0 p2 c0 {2,S} {4,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
"""),
    E0 = (-18.6882,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([572.4,573.055,574.712,1022.08,1278.6,1281.28],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (63.0078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.06033,-0.00588283,7.56384e-05,-1.39706e-07,8.13069e-11,-2246.24,8.76778], Tmin=(10,'K'), Tmax=(571.945,'K')), NASAPolynomial(coeffs=[3.59104,0.0118581,-8.80966e-06,2.92804e-09,-3.59493e-13,-2429.05,8.70221], Tmin=(571.945,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-18.6882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""F[C]1OO1""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O][C](F)OF(200)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {3,S}
3 O u0 p2 c0 {2,S} {5,S}
4 O u1 p2 c0 {5,S}
5 C u1 p0 c0 {1,S} {3,S} {4,S}
"""),
    E0 = (-52.9566,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,482,586,761,1411,355.197,355.229],'cm^-1')),
        HinderedRotor(inertia=(0.00133632,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.0062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.85551,0.0190494,-4.6901e-06,-1.21077e-08,6.97028e-12,-6322.51,16.4643], Tmin=(100,'K'), Tmax=(979.257,'K')), NASAPolynomial(coeffs=[10.6225,0.00181399,-4.85472e-07,1.4067e-10,-1.44207e-14,-8538.46,-24.3915], Tmin=(979.257,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-52.9566,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-F1sO2sH)) + radical(CsF1sO2sO2s)"""),
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
    E0 = (-85.8359,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (97.3618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (134.906,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-82.9076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-31.5506,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-70.7263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-259.364,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (164.009,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (97.3618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-152.934,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-205.973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (132.847,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction4',
    reactants = ['O(6)', 'CF2O(48)'],
    products = ['[O]O[C](F)F(147)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(2.18075e+08,'m^3/(mol*s)'), n=-0.990129, Ea=(212.026,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05768407177011316, var=24.87296733003861, Tref=1000.0, N=33, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C
Ea raised from 210.2 to 212.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction1',
    reactants = ['O(6)', '[O][C](F)F(178)'],
    products = ['[O]O[C](F)F(147)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction2',
    reactants = ['O2(2)', 'CF2(42)'],
    products = ['[O]O[C](F)F(147)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4.0899e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]O[C](F)F(147)'],
    products = ['FC1(F)OO1(932)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Cpri_rad_out_noH]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O2(2)', 'CF2(89)'],
    products = ['[O]O[C](F)F(147)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2560.47,'cm^3/(mol*s)'), n=2.72845, Ea=(64.1698,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), Tmax=(2500,'K'), comment="""Fitted to 50 data points; dA = *|/ 1.0708, dn = +|- 0.00898775, dEa = +|- 0.048911 kJ/molMatched reaction 10 O2 + CF2 <=> CF2O2 in halocarbene_recombination/training
This reaction matched rate rule [Root_N-3R->H_Ext-3BrCClFINOPSSi-R_N-3BrCClFINOPSSi->C_Ext-1C2s-R]
family: halocarbene_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O(172)', 'CF2O(48)'],
    products = ['FC1(F)OO1(932)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(48.1575,'cm^3/(mol*s)'), n=2.76922, Ea=(37.8393,'kJ/mol'), T0=(1,'K'), Tmin=(298,'K'), Tmax=(2500,'K'), comment="""Estimated using template [elec_def;mb_carbonyl_2H] for rate rule [o_atom_singlet;mb_carbonyl_2H]
Euclidian distance = 1.0
family: 1+2_Cycloaddition"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C([O])(F)F(152)'],
    products = ['FC1(F)OO1(932)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F(37)', 'F[C]1OO1(174)'],
    products = ['FC1(F)OO1(932)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.71905e+10,'m^3/(mol*s)'), n=-0.966851, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.07464897564796032, var=0.299509236916578, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing"""),
)

reaction(
    label = 'reaction9',
    reactants = ['O(6)', '[O][C](F)F(178)'],
    products = ['[O]C([O])(F)F(152)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction10',
    reactants = ['F(37)', '[O]C(=O)F(143)'],
    products = ['[O]C([O])(F)F(152)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.08261e+20,'m^3/(mol*s)'), n=-5.07836, Ea=(45.7999,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O"""),
)

reaction(
    label = 'reaction11',
    reactants = ['O(6)', 'CF2O(48)'],
    products = ['[O]C([O])(F)F(152)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.08261e+20,'m^3/(mol*s)'), n=-5.07836, Ea=(91.8889,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O][C](F)OF(200)'],
    products = ['[O]C([O])(F)F(152)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(75.998,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

network(
    label = 'PDepNetwork #172',
    isomers = [
        '[O]O[C](F)F(147)',
        'FC1(F)OO1(932)',
        '[O]C([O])(F)F(152)',
    ],
    reactants = [
        ('O(6)', 'CF2O(48)'),
        ('O2(2)', 'CF2(42)'),
        ('F(37)', '[O]C(=O)F(143)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #172',
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

