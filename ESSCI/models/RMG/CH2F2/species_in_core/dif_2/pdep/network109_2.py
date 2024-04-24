species(
    label = '[O]OC(F)(F)[C](F)F(269)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {8,S}
4 F u0 p3 c0 {8,S}
5 O u0 p2 c0 {6,S} {7,S}
6 O u1 p2 c0 {5,S}
7 C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8 C u1 p0 c0 {3,S} {4,S} {7,S}
"""),
    E0 = (-656.445,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,253,525,597,667,842,1178,1324,190,488,555,1236,1407,180],'cm^-1')),
        HinderedRotor(inertia=(0.423859,'amu*angstrom^2'), symmetry=1, barrier=(9.74534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.427728,'amu*angstrom^2'), symmetry=1, barrier=(9.8343,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.014,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3380.85,'J/mol'), sigma=(5.63196,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=528.08 K, Pc=42.94 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.752133,0.0844089,-0.000171854,1.63425e-07,-5.74727e-11,-78847.6,24.213], Tmin=(100,'K'), Tmax=(870.769,'K')), NASAPolynomial(coeffs=[8.29005,0.020595,-1.16481e-05,2.27713e-09,-1.5498e-13,-79053.8,-4.75611], Tmin=(870.769,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-656.445,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(ROOJ) + radical(Csj(Cs-F1sF1sO2s)(F1s)(F1s))"""),
)

species(
    label = 'FC1(F)OOC1(F)F(248)',
    structure = adjacencyList("""1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {8,S}
4 F u0 p3 c0 {8,S}
5 O u0 p2 c0 {6,S} {7,S}
6 O u0 p2 c0 {5,S} {8,S}
7 C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8 C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
"""),
    E0 = (-900.264,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.014,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3181.58,'J/mol'), sigma=(5.70556,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=496.96 K, Pc=38.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.82235,0.0116335,0.000134205,-3.35264e-07,2.37701e-10,-108274,11.3155], Tmin=(10,'K'), Tmax=(496.751,'K')), NASAPolynomial(coeffs=[5.33001,0.0317908,-2.41881e-05,8.1937e-09,-1.02179e-12,-108822,1.08586], Tmin=(496.751,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-900.264,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), label="""FC1(F)OOC1(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'CF2CF2(60)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {6,S}
4 F u0 p3 c0 {6,S}
5 C u0 p0 c0 {1,S} {2,S} {6,D}
6 C u0 p0 c0 {3,S} {4,S} {5,D}
"""),
    E0 = (-688.535,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(99.9936,'amu')),
        NonlinearRotor(inertia=([91.6969,155.94,247.638],'amu*angstrom^2'), symmetry=4),
        HarmonicOscillator(frequencies=([198.559,203.927,398.069,428.925,524.86,551.48,555.021,803.576,1208.21,1359.53,1361.83,1918.98],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.015,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2113.54,'J/mol'), sigma=(4.647,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.86796,0.00896674,9.20781e-05,-2.50451e-07,1.90637e-10,-82806.8,9.05323], Tmin=(10,'K'), Tmax=(470.327,'K')), NASAPolynomial(coeffs=[5.38595,0.0191802,-1.42425e-05,4.78608e-09,-5.97116e-13,-83205.3,0.155967], Tmin=(470.327,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-688.535,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""FC(F)DC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=O(177)',
    structure = adjacencyList("""1 O u0 p2 c0 {2,D}
2 O u0 p2 c0 {1,D}
"""),
    E0 = (85.6848,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1487.4],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1857.18,'J/mol'), sigma=(4.34667,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=290.09 K, Pc=51.31 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53732,-0.00121572,5.3162e-06,-4.89446e-09,1.45846e-12,10304.5,4.68368], Tmin=(100,'K'), Tmax=(1074.55,'K')), NASAPolynomial(coeffs=[3.15382,0.00167804,-7.69974e-07,1.51275e-10,-1.08782e-14,10302.3,6.16756], Tmin=(1074.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(85.6848,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2(S)""", comment="""Thermo library: primaryThermoLibrary"""),
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
    label = 'FC1(F)OO1(924)',
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
    label = '[O]C(F)(F)[C](F)F(303)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {7,S}
4 F u0 p3 c0 {7,S}
5 O u1 p2 c0 {6,S}
6 C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
7 C u1 p0 c0 {3,S} {4,S} {6,S}
"""),
    E0 = (-648.922,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,525,597,667,842,1178,1324,190,488,555,1236,1407,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.240999,'amu*angstrom^2'), symmetry=1, barrier=(5.54105,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (116.014,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51581,0.0614317,-0.000112934,1.02531e-07,-3.55816e-11,-77964.4,20.7066], Tmin=(100,'K'), Tmax=(832.95,'K')), NASAPolynomial(coeffs=[8.56232,0.0155595,-8.65617e-06,1.72641e-09,-1.20589e-13,-78720.8,-9.49499], Tmin=(832.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-648.922,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sF1s)) + radical(Csj(Cs-F1sF1sO2s)(F1s)(F1s))"""),
)

species(
    label = 'FC1(F)OC1(F)F(776)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {7,S}
4 F u0 p3 c0 {7,S}
5 O u0 p2 c0 {6,S} {7,S}
6 C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
7 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
"""),
    E0 = (-911.165,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.014,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.86803,0.00818545,0.000108116,-2.50383e-07,1.65379e-10,-109579,10.2581], Tmin=(10,'K'), Tmax=(532.564,'K')), NASAPolynomial(coeffs=[5.14834,0.0260802,-1.97726e-05,6.70734e-09,-8.39604e-13,-110105,1.22508], Tmin=(532.564,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-911.165,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""FC1(F)OC1(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]OC(F)=C(F)F(381)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {7,S}
4 O u0 p2 c0 {5,S} {6,S}
5 O u1 p2 c0 {4,S}
6 C u0 p0 c0 {1,S} {4,S} {7,D}
7 C u0 p0 c0 {2,S} {3,S} {6,D}
"""),
    E0 = (-425.108,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,326,540,652,719,1357,182,240,577,636,1210,1413],'cm^-1')),
        HinderedRotor(inertia=(0.397918,'amu*angstrom^2'), symmetry=1, barrier=(9.14891,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (113.015,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.70698,0.0237193,0.000112193,-4.70654e-07,4.8274e-10,-51130.1,12.1246], Tmin=(10,'K'), Tmax=(374.538,'K')), NASAPolynomial(coeffs=[7.23257,0.020879,-1.58533e-05,5.42931e-09,-6.88263e-13,-51638.3,-4.68163], Tmin=(374.538,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-425.108,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""[O]OC(F)DC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[C](F)[C](F)F(230)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {6,S}
4 F u0 p3 c0 {6,S}
5 C u1 p0 c0 {1,S} {2,S} {6,S}
6 C u1 p0 c0 {3,S} {4,S} {5,S}
"""),
    E0 = (-493.512,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([146,234,414,562,504,606,1176,1296,1354,1460,180],'cm^-1')),
        HinderedRotor(inertia=(0.200092,'amu*angstrom^2'), symmetry=1, barrier=(4.6005,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.015,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.04125,0.0516239,-0.000104942,1.02887e-07,-3.74503e-11,-59293.5,18.0917], Tmin=(100,'K'), Tmax=(852.08,'K')), NASAPolynomial(coeffs=[5.53325,0.0162398,-9.21972e-06,1.83691e-09,-1.27509e-13,-59199.1,5.8493], Tmin=(852.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-493.512,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCsF1sF1s) + radical(CsCsF1sF1s)"""),
)

species(
    label = 'FOO[C](F)[C](F)F(777)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {8,S}
3 F u0 p3 c0 {8,S}
4 F u0 p3 c0 {6,S}
5 O u0 p2 c0 {6,S} {7,S}
6 O u0 p2 c0 {4,S} {5,S}
7 C u1 p0 c0 {1,S} {5,S} {8,S}
8 C u1 p0 c0 {2,S} {3,S} {7,S}
"""),
    E0 = (-290.572,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,277,555,632,395,473,707,1436,190,488,555,1236,1407],'cm^-1')),
        HinderedRotor(inertia=(1.30235,'amu*angstrom^2'), symmetry=1, barrier=(29.9436,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.014,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.799659,0.0835489,-0.000169671,1.65945e-07,-6.07233e-11,-34845.2,25.7755], Tmin=(100,'K'), Tmax=(838.183,'K')), NASAPolynomial(coeffs=[7.10837,0.0253387,-1.52048e-05,3.08532e-09,-2.16705e-13,-34915.6,2.34226], Tmin=(838.183,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-290.572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCsF1sO2s) + radical(Csj(Cs-F1sO2sH)(F1s)(F1s))"""),
)

species(
    label = '[O]O[C](F)C(F)(F)F(778)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {7,S}
4 F u0 p3 c0 {8,S}
5 O u0 p2 c0 {6,S} {8,S}
6 O u1 p2 c0 {5,S}
7 C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
8 C u1 p0 c0 {4,S} {5,S} {7,S}
"""),
    E0 = (-683.094,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.014,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.8648,0.0828448,-0.000171391,1.65948e-07,-5.914e-11,-82057.6,23.0672], Tmin=(100,'K'), Tmax=(871.151,'K')), NASAPolynomial(coeffs=[7.04607,0.0225711,-1.26951e-05,2.47979e-09,-1.68832e-13,-81924.5,1.0444], Tmin=(871.151,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-683.094,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(ROOJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'FC(F)(F)C1(F)OO1(268)',
    structure = adjacencyList("""1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {8,S}
3 F u0 p3 c0 {8,S}
4 F u0 p3 c0 {8,S}
5 O u0 p2 c0 {6,S} {7,S}
6 O u0 p2 c0 {5,S} {7,S}
7 C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
8 C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
"""),
    E0 = (-907.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.014,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.74448,0.0213141,7.8072e-05,-2.13069e-07,1.46089e-10,-109162,12.1787], Tmin=(10,'K'), Tmax=(540.568,'K')), NASAPolynomial(coeffs=[6.61394,0.0286041,-2.13036e-05,7.09797e-09,-8.74448e-13,-109889,-3.75431], Tmin=(540.568,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-907.66,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), label="""FC(F)(F)C1(F)OO1""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[C](F)OO[C](F)F(257)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {8,S}
4 F u0 p3 c0 {8,S}
5 O u0 p2 c0 {6,S} {7,S}
6 O u0 p2 c0 {5,S} {8,S}
7 C u1 p0 c0 {1,S} {2,S} {5,S}
8 C u1 p0 c0 {3,S} {4,S} {6,S}
"""),
    E0 = (-589.461,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,444,542,521,679,613,787,1100,1188,1236,1350,180],'cm^-1')),
        HinderedRotor(inertia=(2.30566,'amu*angstrom^2'), symmetry=1, barrier=(53.0118,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.12,'amu*angstrom^2'), symmetry=1, barrier=(48.7431,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.549417,'amu*angstrom^2'), symmetry=1, barrier=(12.6322,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.014,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2940.86,'J/mol'), sigma=(4.81151,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=459.36 K, Pc=59.91 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.804083,0.0772919,-0.000137217,1.19706e-07,-4.0206e-11,-70787.4,22.8649], Tmin=(100,'K'), Tmax=(827.608,'K')), NASAPolynomial(coeffs=[11.2105,0.0168514,-9.28544e-06,1.84218e-09,-1.28378e-13,-72162.4,-23.272], Tmin=(827.608,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-589.461,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(F1s)(F1s)(O2s-O2s)) + radical(Csj(F1s)(F1s)(O2s-O2s))"""),
)

species(
    label = '[O]C(F)(F)C([O])(F)F(255)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {8,S}
4 F u0 p3 c0 {8,S}
5 O u1 p2 c0 {7,S}
6 O u1 p2 c0 {8,S}
7 C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8 C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
"""),
    E0 = (-812.408,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([261,441,281,365,497,569,525,693,557,771,857,927,1026,1214,1153,1249,3050.6],'cm^-1')),
        HinderedRotor(inertia=(0.481645,'amu*angstrom^2'), symmetry=1, barrier=(11.074,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.014,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3805.25,'J/mol'), sigma=(6.44294,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=594.37 K, Pc=32.28 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40291,0.0622121,-9.92085e-05,8.15334e-08,-2.6593e-11,-97621.4,20.5684], Tmin=(100,'K'), Tmax=(752.636,'K')), NASAPolynomial(coeffs=[9.891,0.0170985,-9.29283e-06,1.88426e-09,-1.34904e-13,-98899,-17.9692], Tmin=(752.636,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-812.408,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sF1s)) + radical(O2sj(Cs-CsF1sF1s))"""),
)

species(
    label = 'F[C]1OOC1(F)F(270)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {7,S}
4 O u0 p2 c0 {5,S} {6,S}
5 O u0 p2 c0 {4,S} {7,S}
6 C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
7 C u1 p0 c0 {3,S} {5,S} {6,S}
"""),
    E0 = (-468.911,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (113.015,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.86334,0.00871615,0.000110259,-2.64044e-07,1.80541e-10,-56392.5,12.5663], Tmin=(10,'K'), Tmax=(512.679,'K')), NASAPolynomial(coeffs=[4.9487,0.0268092,-2.03911e-05,6.90197e-09,-8.5981e-13,-56852.9,4.65089], Tmin=(512.679,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-468.911,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""F[C]1OOC1(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'FC1=C(F)OO1(271)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {4,S} {5,S}
4 O u0 p2 c0 {3,S} {6,S}
5 C u0 p0 c0 {1,S} {3,S} {6,D}
6 C u0 p0 c0 {2,S} {4,S} {5,D}
"""),
    E0 = (-207.71,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([523.035,523.036,523.043,523.043,523.043,523.046,523.049,523.071,1497.32,1497.34,1497.34,1497.36],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.0169,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.90579,0.00578655,8.18583e-05,-1.84343e-07,1.1946e-10,-24975.3,9.3143], Tmin=(10,'K'), Tmax=(538.444,'K')), NASAPolynomial(coeffs=[4.6125,0.0205973,-1.52868e-05,5.12972e-09,-6.38784e-13,-25342.2,3.64195], Tmin=(538.444,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-207.71,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""FC1DC(F)OO1""", comment="""Thermo library: CHOF_G4"""),
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
    E0 = (-285.521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-239.757,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-42.7952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (201.179,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-285.068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-200.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (44.4131,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-139.046,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-36.2611,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (158.457,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-97.409,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-97.3682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-238.171,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-537.171,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-218.083,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-441.031,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-32.9267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (146.578,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction7',
    reactants = ['O2(2)', 'CF2CF2(60)'],
    products = ['[O]OC(F)(F)[C](F)F(269)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(21668.1,'m^3/(mol*s)'), n=0.762578, Ea=(48.547,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_N-3R->C_N-2R!H->O_N-4R!H-u0_Ext-1R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_N-3R->C_N-2R!H->O_N-4R!H-u0_Ext-1R!H-R
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction1',
    reactants = ['[O]OC(F)(F)[C](F)F(269)'],
    products = ['O=O(177)', 'CF2CF2(60)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(53.5956,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission
Ea raised from 0.0 to 53.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction2',
    reactants = ['O(6)', '[O]C(F)(F)[C](F)F(303)'],
    products = ['[O]OC(F)(F)[C](F)F(269)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['CF2(42)', '[O]O[C](F)F(171)'],
    products = ['[O]OC(F)(F)[C](F)F(269)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]OC(F)(F)[C](F)F(269)'],
    products = ['FC1(F)OOC1(F)F(248)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;O_rad;Cpri_rad_out_noH]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]OC(F)(F)[C](F)F(269)'],
    products = ['O(6)', 'FC1(F)OC1(F)F(776)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(8.43325e+10,'s^-1'), n=0.54, Ea=(93.1926,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2OO_S;Cs_rad_intra;OOJ] + [R2OO_S;C_ter_rad_intra;OO] for rate rule [R2OO_S;C_ter_rad_intra;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['F(37)', '[O]OC(F)=C(F)F(381)'],
    products = ['[O]OC(F)(F)[C](F)F(269)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(33.5367,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction8',
    reactants = ['O2(2)', 'F[C](F)[C](F)F(230)'],
    products = ['[O]OC(F)(F)[C](F)F(269)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.53688e+08,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CF2(87)', '[O]O[C](F)F(171)'],
    products = ['[O]OC(F)(F)[C](F)F(269)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(786.723,'m^3/(mol*s)'), n=1.25031, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl"""),
)

reaction(
    label = 'reaction10',
    reactants = ['FOO[C](F)[C](F)F(777)'],
    products = ['[O]OC(F)(F)[C](F)F(269)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(85.9356,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]OC(F)(F)[C](F)F(269)'],
    products = ['[O]O[C](F)C(F)(F)F(778)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(8.5608e+11,'s^-1'), n=0.253035, Ea=(195.943,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_N-4R!H->C_Ext-2R!H-R_N-5R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_N-4R!H->C_Ext-2R!H-R_N-5R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CF2(87)', 'FC1(F)OO1(924)'],
    products = ['FC1(F)OOC1(F)F(248)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.81816e-07,'m^3/(mol*s)'), n=3.45227, Ea=(217.666,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s',), comment="""Estimated from node CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction3',
    reactants = ['FC1(F)OOC1(F)F(248)'],
    products = ['FC(F)(F)C1(F)OO1(268)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction1',
    reactants = ['CF2O(48)', 'CF2O(48)'],
    products = ['FC1(F)OOC1(F)F(248)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.36636e-32,'m^3/(mol*s)'), n=10.9581, Ea=(336.956,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1COCSCdCdd->Cd',), comment="""Estimated from node Root_N-1COCSCdCdd->Cd
Ea raised from 333.8 to 337.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['F[C](F)OO[C](F)F(257)'],
    products = ['FC1(F)OOC1(F)F(248)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_noH;Cpri_rad_out_noH]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]C(F)(F)C([O])(F)F(255)'],
    products = ['FC1(F)OOC1(F)F(248)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F(37)', 'F[C]1OOC1(F)F(270)'],
    products = ['FC1(F)OOC1(F)F(248)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.71905e+10,'m^3/(mol*s)'), n=-0.966851, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.07464897564796032, var=0.299509236916578, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_N-2R-inRing"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F2(106)', 'FC1=C(F)OO1(271)'],
    products = ['FC1(F)OOC1(F)F(248)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

network(
    label = 'PDepNetwork #109',
    isomers = [
        '[O]OC(F)(F)[C](F)F(269)',
        'FC1(F)OOC1(F)F(248)',
    ],
    reactants = [
        ('O2(2)', 'CF2CF2(60)'),
        ('O=O(177)', 'CF2CF2(60)'),
        ('CF2(42)', '[O]O[C](F)F(171)'),
        ('CF2(87)', '[O]O[C](F)F(171)'),
        ('CF2(87)', 'FC1(F)OO1(924)'),
        ('CF2O(48)', 'CF2O(48)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #109',
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

