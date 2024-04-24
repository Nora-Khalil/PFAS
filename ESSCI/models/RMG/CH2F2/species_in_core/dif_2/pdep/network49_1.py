species(
    label = '[O]OOC(=O)F(205)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {3,S} {6,S}
3 O u0 p2 c0 {2,S} {5,S}
4 O u0 p2 c0 {6,D}
5 O u1 p2 c0 {3,S}
6 C u0 p0 c0 {1,S} {2,S} {4,D}
"""),
    E0 = (-411.953,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([482,664,788,1296,1923,449.889,451.278,451.819,452.782,452.884],'cm^-1')),
        HinderedRotor(inertia=(0.287391,'amu*angstrom^2'), symmetry=1, barrier=(41.4839,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000839319,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.0066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.4296,0.0198105,2.4457e-05,-6.03731e-08,2.84566e-11,-49476.3,20.4864], Tmin=(100,'K'), Tmax=(920.993,'K')), NASAPolynomial(coeffs=[17.2536,-0.00628798,4.61052e-06,-8.72873e-10,5.39385e-14,-53830.5,-58.6256], Tmin=(920.993,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-411.953,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-OsOs) + group(O2s-OsH) + group(COFOO) + radical(ROOJ)"""),
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
    label = '[O]OOF(210)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {3,S}
3 O u0 p2 c0 {2,S} {4,S}
4 O u1 p2 c0 {3,S}
"""),
    E0 = (135.84,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([277,555,632,654.187,3258.91],'cm^-1')),
        HinderedRotor(inertia=(0.00635282,'amu*angstrom^2'), symmetry=1, barrier=(1.92928,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (66.9966,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.97115,0.0165468,-1.94483e-05,1.02773e-08,-1.93124e-12,16380.1,14.443], Tmin=(100,'K'), Tmax=(1607.92,'K')), NASAPolynomial(coeffs=[7.78642,-0.00055688,1.28819e-06,-3.02601e-10,2.19178e-14,15494,-9.01478], Tmin=(1607.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(135.84,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(78.9875,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsOs) + group(O2sFO) + group(O2s-OsH) + radical(ROOJ)"""),
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
    label = '[O]OC(=O)F(211)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u0 p2 c0 {5,D}
4 O u1 p2 c0 {2,S}
5 C u0 p0 c0 {1,S} {2,S} {3,D}
"""),
    E0 = (-339.818,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,482,664,788,1296,1923],'cm^-1')),
        HinderedRotor(inertia=(1.35484,'amu*angstrom^2'), symmetry=1, barrier=(31.1504,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.0072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.92964,0.00448623,6.20575e-05,-1.45434e-07,9.86982e-11,-40867.8,10.2015], Tmin=(10,'K'), Tmax=(509.958,'K')), NASAPolynomial(coeffs=[4.25304,0.0157357,-1.15832e-05,3.84885e-09,-4.74388e-13,-41080,7.10146], Tmin=(509.958,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-339.818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""[O]OC(DO)F""", comment="""Thermo library: 2-BTP_G4"""),
)

species(
    label = 'F[C]1OOOO1(212)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {4,S} {6,S}
3 O u0 p2 c0 {5,S} {6,S}
4 O u0 p2 c0 {2,S} {5,S}
5 O u0 p2 c0 {3,S} {4,S}
6 C u1 p0 c0 {1,S} {2,S} {3,S}
"""),
    E0 = (-50.9951,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.0066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.45041,-0.0144644,0.000127147,-1.66872e-07,6.58791e-11,-6088.94,18.2172], Tmin=(100,'K'), Tmax=(929.623,'K')), NASAPolynomial(coeffs=[17.7224,-0.00764369,6.04724e-06,-1.07433e-09,5.94158e-14,-11690.7,-65.4528], Tmin=(929.623,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-50.9951,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsOs) + group(CsFHOO) + ring(Cyclopentane) + radical(CsF1sO2sO2s)"""),
)

species(
    label = '[O]C1(F)OOO1(213)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {4,S} {6,S}
3 O u0 p2 c0 {4,S} {6,S}
4 O u0 p2 c0 {2,S} {3,S}
5 O u1 p2 c0 {6,S}
6 C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
"""),
    E0 = (-213.346,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.0066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.15503,0.0168419,4.76625e-06,-1.21943e-08,3.69849e-12,-25628,13.8359], Tmin=(100,'K'), Tmax=(1401.71,'K')), NASAPolynomial(coeffs=[8.54219,0.0158965,-9.66135e-06,2.01069e-09,-1.44689e-13,-28555.6,-19.0291], Tmin=(1401.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-213.346,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsOs) + group(CsFOOO) + ring(Cyclobutane) + radical(OCOJ)"""),
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
    label = '[O]OO[C]=O(214)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {5,S}
2 O u0 p2 c0 {1,S} {3,S}
3 O u1 p2 c0 {2,S}
4 O u0 p2 c0 {5,D}
5 C u1 p0 c0 {1,S} {4,D}
"""),
    E0 = (7.09814,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,528.404,528.405,528.406,528.407],'cm^-1')),
        HinderedRotor(inertia=(0.22529,'amu*angstrom^2'), symmetry=1, barrier=(44.6382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.225294,'amu*angstrom^2'), symmetry=1, barrier=(44.6383,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (76.0082,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.79698,0.0134039,2.85349e-05,-5.97612e-08,2.76175e-11,909.249,17.108], Tmin=(100,'K'), Tmax=(909.965,'K')), NASAPolynomial(coeffs=[15.1763,-0.00664661,4.9365e-06,-9.69033e-10,6.26932e-14,-2766.53,-49.2653], Tmin=(909.965,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(7.09814,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(99.7737,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-OsOs) + group(O2s-OsH) + group(Cds-OdOsH) + radical(ROOJ) + radical((O)CJOC)"""),
)

species(
    label = '[O]O[O](157)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {3,S}
2 O u1 p2 c0 {1,S}
3 O u1 p2 c0 {1,S}
"""),
    E0 = (192.544,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([767.866,3898.78,3898.78],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (47.9982,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.4747,0.00468953,-6.33899e-06,3.99278e-09,-7.67859e-13,23182.7,11.3223], Tmin=(100,'K'), Tmax=(1873.27,'K')), NASAPolynomial(coeffs=[2.72483,0.000671287,1.37825e-06,-3.55017e-10,2.60932e-14,24449.6,18.0423], Tmin=(1873.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.544,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsOs) + group(O2s-OsH) + group(O2s-OsH) + radical(ROOJ) + radical(ROOJ)"""),
)

species(
    label = 'CFO(50)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u1 p0 c0 {1,S} {2,D}
"""),
    E0 = (-190.359,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(46.9933,'amu')),
        NonlinearRotor(inertia=([2.65864,43.9175,46.5761],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([636.046,1085.22,1918.24],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (47.0084,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(7150.45,'J/mol'), sigma=(4,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02994,-0.00208576,2.20825e-05,-3.30628e-08,1.5841e-11,-22895,6.85532], Tmin=(10,'K'), Tmax=(668.186,'K')), NASAPolynomial(coeffs=[3.39014,0.00554951,-3.6e-06,1.08417e-09,-1.23809e-13,-22894.5,9.04838], Tmin=(668.186,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-190.359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""OD[C]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=[C]OOOF(215)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {6,S}
3 O u0 p2 c0 {2,S} {4,S}
4 O u0 p2 c0 {1,S} {3,S}
5 O u0 p2 c0 {6,D}
6 C u1 p0 c0 {2,S} {5,D}
"""),
    E0 = (-49.6059,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([277,555,632,1855,455,950,476.627,476.633,476.639],'cm^-1')),
        HinderedRotor(inertia=(0.266571,'amu*angstrom^2'), symmetry=1, barrier=(42.9836,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.266609,'amu*angstrom^2'), symmetry=1, barrier=(42.9867,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.266607,'amu*angstrom^2'), symmetry=1, barrier=(42.9894,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.0066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.14419,0.0270087,9.26842e-06,-4.52448e-08,2.27648e-11,-5886.65,20.0731], Tmin=(100,'K'), Tmax=(934.601,'K')), NASAPolynomial(coeffs=[17.6914,-0.00433731,3.09191e-06,-5.46765e-10,3.04158e-14,-10329.8,-62.1047], Tmin=(934.601,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-49.6059,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(120.56,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-OsOs) + group(O2sFO) + group(Cds-OdOsH) + radical((O)CJOC)"""),
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
    E0 = (-286.34,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (185.903,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-8.83736,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (36.951,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-125.4,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (167.936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (90.1313,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (102.697,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O2(2)', '[O]C(=O)F(136)'],
    products = ['[O]OOC(=O)F(205)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(9.96604e+08,'m^3/(mol*s)'), n=0.0472428, Ea=(15.7709,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]OOF(210)', 'CO(12)'],
    products = ['[O]OOC(=O)F(205)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.000742267,'m^3/(mol*s)'), n=2.83796, Ea=(80.8577,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.2632766423807829, var=145.9379134136138, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-1COCbCdCsCtHNOSSidSis->Cs',), comment="""Estimated from node Root_N-1COCbCdCsCtHNOSSidSis->Cs"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O(6)', '[O]OC(=O)F(211)'],
    products = ['[O]OOC(=O)F(205)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]OOC(=O)F(205)'],
    products = ['F[C]1OOOO1(212)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.12332e+09,'s^-1'), n=0.5388, Ea=(360.958,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra] for rate rule [R5_SS_CO;carbonyl_intra;radadd_intra_O]
Euclidian distance = 1.7320508075688772
family: Intra_R_Add_Endocyclic
Ea raised from 356.6 to 361.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]OOC(=O)F(205)'],
    products = ['[O]C1(F)OOO1(213)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(198.607,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_O] for rate rule [R5_SS_CO;carbonylbond_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic
Ea raised from 198.2 to 198.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction6',
    reactants = ['F(37)', '[O]OO[C]=O(214)'],
    products = ['[O]OOC(=O)F(205)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.505e+06,'m^3/(mol*s)'), n=-1.3397e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Sp-4R!H-2C_Ext-2C-R_N-5R!H->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Sp-4R!H-2C_Ext-2C-R_N-5R!H->C"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]O[O](157)', 'CFO(50)'],
    products = ['[O]OOC(=O)F(205)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.01e+06,'m^3/(mol*s)'), n=-1.3397e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Sp-4R!H-2C_Ext-2C-R_N-5R!H->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Sp-4R!H-2C_Ext-2C-R_N-5R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction8',
    reactants = ['O=[C]OOOF(215)'],
    products = ['[O]OOC(=O)F(205)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(64.3569,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

network(
    label = 'PDepNetwork #49',
    isomers = [
        '[O]OOC(=O)F(205)',
    ],
    reactants = [
        ('O2(2)', '[O]C(=O)F(136)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #49',
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

