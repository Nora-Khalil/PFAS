species(
    label = '[O]OOC([O])(F)F(947)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {4,S} {7,S}
4 O u0 p2 c0 {3,S} {6,S}
5 O u1 p2 c0 {7,S}
6 O u1 p2 c0 {4,S}
7 C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
"""),
    E0 = (-352.939,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,486,586,623,710,938,1161,1294,180,180,241.911,1040.92,1041.92],'cm^-1')),
        HinderedRotor(inertia=(0.0027375,'amu*angstrom^2'), symmetry=1, barrier=(2.10692,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.579821,'amu*angstrom^2'), symmetry=1, barrier=(13.3312,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74127,0.0449978,-5.13138e-05,2.60632e-08,-4.66931e-12,-42363.3,23.7156], Tmin=(100,'K'), Tmax=(1018.18,'K')), NASAPolynomial(coeffs=[13.9748,0.00404505,-1.45275e-06,2.72227e-10,-2.01704e-14,-45222.9,-37.3339], Tmin=(1018.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-352.939,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsOs) + group(O2s-OsH) + group(CsFFOO) + radical(O2sj(Cs-F1sF1sO2s)) + radical(ROOJ)"""),
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
    label = 'FC1(F)OO1(934)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {4,S} {5,S}
4 O u0 p2 c0 {3,S} {5,S}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {4,S}
"""),
    E0 = (-474.414,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([180,180,761.365,1045.89,1045.89,1045.9,1045.91,1045.91,1045.93],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0062,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2976.61,'J/mol'), sigma=(5.20996,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=464.94 K, Pc=47.76 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.95302,0.00256699,5.22953e-05,-9.50593e-08,5.05491e-11,-57057.6,8.08501], Tmin=(10,'K'), Tmax=(633.13,'K')), NASAPolynomial(coeffs=[3.1996,0.0191329,-1.49226e-05,5.17103e-09,-6.57667e-13,-57198.8,9.5068], Tmin=(633.13,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-474.414,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""FC1(F)OO1""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]C([O])(F)F(144)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,1.49298e-14,-2.05327e-17,9.23927e-21,-1.28064e-24,29230.2,5.12616], Tmin=(100,'K'), Tmax=(3959.16,'K')), NASAPolynomial(coeffs=[2.5,1.19009e-10,-4.23987e-14,6.68965e-18,-3.94352e-22,29230.2,5.12617], Tmin=(3959.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.034,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[O]OC([O])(F)F(242)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {5,S} {6,S}
4 O u1 p2 c0 {6,S}
5 O u1 p2 c0 {3,S}
6 C u0 p0 c0 {1,S} {2,S} {3,S} {4,S}
"""),
    E0 = (-388.916,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,375,486,586,623,710,938,1161,1294],'cm^-1')),
        HinderedRotor(inertia=(0.563774,'amu*angstrom^2'), symmetry=1, barrier=(12.9623,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.30344,0.038995,-5.35335e-05,3.63417e-08,-9.66956e-12,-46716.1,17.5892], Tmin=(100,'K'), Tmax=(922.714,'K')), NASAPolynomial(coeffs=[9.26777,0.00880528,-4.45747e-06,8.85014e-10,-6.32285e-14,-48001.3,-15.4497], Tmin=(922.714,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-388.916,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(CsFFOO) + radical(O2sj(Cs-F1sF1sO2s)) + radical(ROOJ)"""),
)

species(
    label = '[O]OO[C](F)F(241)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {4,S} {6,S}
4 O u0 p2 c0 {3,S} {5,S}
5 O u1 p2 c0 {4,S}
6 C u1 p0 c0 {1,S} {2,S} {3,S}
"""),
    E0 = (-188.302,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([493,600,700,1144,1293,180,180,887.477,887.489,887.523],'cm^-1')),
        HinderedRotor(inertia=(0.117559,'amu*angstrom^2'), symmetry=1, barrier=(2.70292,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.117522,'amu*angstrom^2'), symmetry=1, barrier=(2.70205,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.29738,0.0349472,-4.16785e-05,2.29342e-08,-4.55973e-12,-22584,20.7893], Tmin=(100,'K'), Tmax=(964.794,'K')), NASAPolynomial(coeffs=[10.8918,0.00414046,-1.28441e-06,2.0607e-10,-1.36549e-14,-24466.9,-21.5296], Tmin=(964.794,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-188.302,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsH) + group(CsFFHO) + radical(ROOJ) + radical(Csj(F1s)(F1s)(O2s-O2s))"""),
)

species(
    label = 'FC1(F)OOOO1(948)',
    structure = adjacencyList("""1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {5,S} {7,S}
4 O u0 p2 c0 {6,S} {7,S}
5 O u0 p2 c0 {3,S} {6,S}
6 O u0 p2 c0 {4,S} {5,S}
7 C u0 p0 c0 {1,S} {2,S} {3,S} {4,S}
"""),
    E0 = (-489.304,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.82883,-0.000299021,0.000102878,-1.46982e-07,5.95113e-11,-58783.5,18.9233], Tmin=(100,'K'), Tmax=(935.246,'K')), NASAPolynomial(coeffs=[19.7736,-0.00600519,4.94725e-06,-8.43427e-10,4.28809e-14,-64873,-77.3018], Tmin=(935.246,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-489.304,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsOs) + group(CsFFOO) + ring(Cyclopentane)"""),
)

species(
    label = 'FC1(F)OOO1(949)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {5,S} {6,S}
4 O u0 p2 c0 {5,S} {6,S}
5 O u0 p2 c0 {3,S} {4,S}
6 C u0 p0 c0 {1,S} {2,S} {3,S} {4,S}
"""),
    E0 = (-443.514,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.86059,0.00715214,5.99531e-05,-9.0459e-08,3.64106e-11,-53285.2,14.5026], Tmin=(100,'K'), Tmax=(959.033,'K')), NASAPolynomial(coeffs=[15.4873,-0.000215728,6.29759e-07,2.83426e-11,-1.57331e-14,-57790.1,-56.7462], Tmin=(959.033,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-443.514,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(CsFFOO) + ring(Cyclobutane)"""),
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
    label = '[O]OOC(=O)F(146)',
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
        HarmonicOscillator(frequencies=([482,664,788,1296,1923,589.236,589.237,589.241,589.241,589.242],'cm^-1')),
        HinderedRotor(inertia=(0.116089,'amu*angstrom^2'), symmetry=1, barrier=(28.6019,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116086,'amu*angstrom^2'), symmetry=1, barrier=(28.6018,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.0066,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3633.72,'J/mol'), sigma=(5.52188,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=567.58 K, Pc=48.97 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.42965,0.0198099,2.44594e-05,-6.03764e-08,2.84581e-11,-49476.3,20.4862], Tmin=(100,'K'), Tmax=(920.985,'K')), NASAPolynomial(coeffs=[17.2534,-0.00628772,4.61037e-06,-8.72837e-10,5.39354e-14,-53830.4,-58.6248], Tmin=(920.985,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-411.953,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-OsOs) + group(O2s-OsH) + group(COFOO) + radical(ROOJ)"""),
)

species(
    label = '[O]O[O](152)',
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3633.48,'J/mol'), sigma=(5.77645,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=567.54 K, Pc=42.77 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.47471,0.00468947,-6.33886e-06,3.99267e-09,-7.67832e-13,23182.7,11.3223], Tmin=(100,'K'), Tmax=(1873.31,'K')), NASAPolynomial(coeffs=[2.72432,0.000671973,1.37793e-06,-3.54951e-10,2.60883e-14,24449.9,18.0454], Tmin=(1873.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.544,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsOs) + group(O2s-OsH) + group(O2s-OsH) + radical(ROOJ) + radical(ROOJ)"""),
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
    label = '[O][C](F)F(257)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.105,0.0167699,-1.53788e-05,4.83907e-09,3.8682e-15,-30692.1,11.7073], Tmin=(100,'K'), Tmax=(1048.23,'K')), NASAPolynomial(coeffs=[8.58999,0.00108969,-4.53602e-07,1.24872e-10,-1.13713e-14,-32130.4,-16.3888], Tmin=(1048.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-255.477,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsOJ) + radical(CsF1sF1sO2s)"""),
)

species(
    label = '[O]OO[C](F)OF(950)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {5,S} {7,S}
4 O u0 p2 c0 {2,S} {7,S}
5 O u0 p2 c0 {3,S} {6,S}
6 O u1 p2 c0 {5,S}
7 C u1 p0 c0 {1,S} {3,S} {4,S}
"""),
    E0 = (-46.6344,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,482,586,761,1411,180,180,180,447.348,840.717,840.749],'cm^-1')),
        HinderedRotor(inertia=(0.439883,'amu*angstrom^2'), symmetry=1, barrier=(10.1138,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00748901,'amu*angstrom^2'), symmetry=1, barrier=(3.75499,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13973,'amu*angstrom^2'), symmetry=1, barrier=(70.0722,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71741,0.0508606,-6.90479e-05,4.46856e-08,-1.11665e-11,-5527.11,24.6035], Tmin=(100,'K'), Tmax=(987.507,'K')), NASAPolynomial(coeffs=[12.2826,0.00806626,-4.04626e-06,8.04267e-10,-5.77587e-14,-7613.82,-26.2352], Tmin=(987.507,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-46.6344,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(145.503,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2sCF) + group(O2s-OsOs) + group(O2s-OsH) + group(CsFHOO) + radical(ROOJ) + radical(CsF1sO2sO2s)"""),
)

species(
    label = '[O][C](F)OOOF(951)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {4,S} {7,S}
4 O u0 p2 c0 {3,S} {5,S}
5 O u0 p2 c0 {2,S} {4,S}
6 O u1 p2 c0 {7,S}
7 C u1 p0 c0 {1,S} {3,S} {6,S}
"""),
    E0 = (17.6848,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([277,555,632,482,586,761,1411,463.023,463.095,463.103,463.111,463.147],'cm^-1')),
        HinderedRotor(inertia=(0.000786031,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000786192,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000785958,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76666,0.0424752,-4.45126e-05,2.16357e-08,-4.01694e-12,2213.01,26.4149], Tmin=(100,'K'), Tmax=(1329.95,'K')), NASAPolynomial(coeffs=[14.2369,0.00496991,-2.21242e-06,4.32125e-10,-3.12044e-14,-1103.98,-37.3025], Tmin=(1329.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(17.6848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(145.503,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsOs) + group(O2sFO) + group(CsFHOO) + radical(O2sj(Cs-F1sO2sH)) + radical(CsF1sO2sO2s)"""),
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
    E0 = (-113.563,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-4.19647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (196.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-204.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-57.6223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-130.433,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-164.073,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-198.629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (78.7517,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (174.269,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (222.084,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]OOC([O])(F)F(947)'],
    products = ['O2(2)', 'FC1(F)OO1(934)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(4.11041e+14,'s^-1'), n=-0.372113, Ea=(97.6903,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO_S;Y_rad_intra;OO_intra]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O(6)', '[O]OC([O])(F)F(242)'],
    products = ['[O]OOC([O])(F)F(947)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['O(6)', '[O]OO[C](F)F(241)'],
    products = ['[O]OOC([O])(F)F(947)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_ter_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]OOC([O])(F)F(947)'],
    products = ['FC1(F)OOOO1(948)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SSSS;Y_rad_out;Ypri_rad_out] for rate rule [R5_SSSS;O_rad;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]OOC([O])(F)F(947)'],
    products = ['O(6)', 'FC1(F)OOO1(949)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.67697e+10,'s^-1'), n=0.54, Ea=(153.631,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnOO;Y_rad_intra;OOJ] + [R3OO;Y_rad_intra;OO] for rate rule [R3OO;Y_rad_intra;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['F(37)', '[O]OOC(=O)F(146)'],
    products = ['[O]OOC([O])(F)F(947)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4.08261e+20,'m^3/(mol*s)'), n=-5.07836, Ea=(66.943,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]O[O](152)', 'CF2O(48)'],
    products = ['[O]OOC([O])(F)F(947)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(8.49e-08,'m^3/(mol*s)'), n=3.486, Ea=(120.308,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_N-3R->C_2R!H->O_Ext-1R!H-R_Ext-1R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_N-3R->C_2R!H->O_Ext-1R!H-R_Ext-1R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction8',
    reactants = ['O2(2)', '[O]C([O])(F)F(144)'],
    products = ['[O]OOC([O])(F)F(947)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(9.96604e+08,'m^3/(mol*s)'), n=0.0472428, Ea=(40.4108,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C
Multiplied by reaction path degeneracy 4.0"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]O[O](152)', '[O][C](F)F(257)'],
    products = ['[O]OOC([O])(F)F(947)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.01e+06,'m^3/(mol*s)'), n=-1.3397e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Sp-4R!H-2C_Ext-2C-R_N-5R!H->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Sp-4R!H-2C_Ext-2C-R_N-5R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]OO[C](F)OF(950)'],
    products = ['[O]OOC([O])(F)F(947)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(79.218,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O][C](F)OOOF(951)'],
    products = ['[O]OOC([O])(F)F(947)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(62.714,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

network(
    label = 'PDepNetwork #336',
    isomers = [
        '[O]OOC([O])(F)F(947)',
    ],
    reactants = [
        ('O2(2)', 'FC1(F)OO1(934)'),
        ('O2(2)', '[O]C([O])(F)F(144)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #336',
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

