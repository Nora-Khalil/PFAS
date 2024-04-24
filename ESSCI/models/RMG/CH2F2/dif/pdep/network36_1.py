species(
    label = '[O]OOOC(F)F(177)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {4,S} {7,S}
4 O u0 p2 c0 {3,S} {5,S}
5 O u0 p2 c0 {4,S} {6,S}
6 O u1 p2 c0 {5,S}
7 C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
8 H u0 p0 c0 {7,S}
"""),
    E0 = (-370.339,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([519,593,830,1115,1166,1389,1413,3114,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.013,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51475,0.0450478,-4.82382e-05,2.42478e-08,-4.56694e-12,-44443.8,27.2206], Tmin=(100,'K'), Tmax=(1452.69,'K')), NASAPolynomial(coeffs=[14.7734,0.00366626,-4.7656e-07,1.95607e-11,5.57074e-17,-47781.7,-39.9254], Tmin=(1452.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-370.339,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsOs) + group(O2s-OsH) + group(CsFFHO) + radical(ROOJ)"""),
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
    label = '[O]OC(F)F(168)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {4,S} {5,S}
4 O u1 p2 c0 {3,S}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
6 H u0 p0 c0 {5,S}
"""),
    E0 = (-427.954,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,519,593,830,1115,1166,1389,1413,3114],'cm^-1')),
        HinderedRotor(inertia=(0.586675,'amu*angstrom^2'), symmetry=1, barrier=(13.4888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0142,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3217.18,'J/mol'), sigma=(5.19785,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=502.52 K, Pc=51.98 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.92725,0.00470123,7.24333e-05,-1.66658e-07,1.13367e-10,-51470.2,11.559], Tmin=(10,'K'), Tmax=(497.316,'K')), NASAPolynomial(coeffs=[3.77619,0.0196241,-1.39226e-05,4.52949e-09,-5.50671e-13,-51624.7,10.478], Tmin=(497.316,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-427.954,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""[O]OC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]OOO(169)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,S} {3,S}
2 O u0 p2 c0 {1,S} {5,S}
3 O u0 p2 c0 {1,S} {4,S}
4 O u1 p2 c0 {3,S}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (76.5163,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([598.567,599.352,600.754,600.834,601.406,2834.58,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00046881,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00046714,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (65.0056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.25054,0.00709423,2.63212e-05,-5.01068e-08,2.3332e-11,9238.71,16.175], Tmin=(100,'K'), Tmax=(881.876,'K')), NASAPolynomial(coeffs=[11.8239,-0.00611983,5.13017e-06,-1.07632e-09,7.44321e-14,6728.3,-29.7686], Tmin=(881.876,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(76.5163,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(99.7737,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsOs) + group(O2s-OsOs) + group(O2s-OsH) + group(O2s-OsH) + radical(ROOJ)"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28591,0.0107608,-1.05382e-05,4.89881e-09,-8.86384e-13,-24340.7,13.1348], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[5.33121,0.00197748,-9.60248e-07,2.10704e-10,-1.5954e-14,-25190.9,-2.56367], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-203.712,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CF2""", comment="""Thermo library: halogens"""),
)

species(
    label = '[O]OOOF(221)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {3,S} {4,S}
3 O u0 p2 c0 {1,S} {2,S}
4 O u0 p2 c0 {2,S} {5,S}
5 O u1 p2 c0 {4,S}
"""),
    E0 = (171.817,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([277,555,632,744.041,748.244,750.242,3896.98],'cm^-1')),
        HinderedRotor(inertia=(0.0065871,'amu*angstrom^2'), symmetry=1, barrier=(2.66102,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0367119,'amu*angstrom^2'), symmetry=1, barrier=(15.0512,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (82.996,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.88455,0.0168411,3.25546e-06,-2.74102e-08,1.51779e-11,20712.3,18.8695], Tmin=(100,'K'), Tmax=(888.361,'K')), NASAPolynomial(coeffs=[12.5058,-0.00504901,4.03081e-06,-8.36267e-10,5.74009e-14,18157.2,-31.1679], Tmin=(888.361,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(171.817,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(99.7737,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsOs) + group(O2s-OsOs) + group(O2sFO) + group(O2s-OsH) + radical(ROOJ)"""),
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
    label = 'O(6)',
    structure = adjacencyList("""multiplicity 3
1 O u2 p2 c0
"""),
    E0 = (243.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (15.9994,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,1.49298e-14,-2.05327e-17,9.23927e-21,-1.28064e-24,29230.2,5.12616], Tmin=(100,'K'), Tmax=(3959.16,'K')), NASAPolynomial(coeffs=[2.5,1.19009e-10,-4.23987e-14,6.68965e-18,-3.94352e-22,29230.2,5.12617], Tmin=(3959.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.034,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[O]OOC(F)F(222)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {4,S} {6,S}
4 O u0 p2 c0 {3,S} {5,S}
5 O u1 p2 c0 {4,S}
6 C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
7 H u0 p0 c0 {6,S}
"""),
    E0 = (-406.316,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([519,593,830,1115,1166,1389,1413,3114,369.567,370.141,371.394,371.669,1852.3],'cm^-1')),
        HinderedRotor(inertia=(0.00121271,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0887881,'amu*angstrom^2'), symmetry=1, barrier=(8.66489,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.0136,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.30892,0.0364485,-4.19925e-05,2.43451e-08,-5.53863e-12,-48806.9,20.2522], Tmin=(100,'K'), Tmax=(1076.17,'K')), NASAPolynomial(coeffs=[9.6708,0.0090853,-3.85275e-06,7.18309e-10,-5.00064e-14,-50391.5,-15.805], Tmin=(1076.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-406.316,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsH) + group(CsFFHO) + radical(ROOJ)"""),
)

species(
    label = '[O-][O+]=O(223)',
    structure = adjacencyList("""1 O u0 p2 c0 {2,S} {3,S}
2 O u0 p2 c0 {1,S} {3,S}
3 O u0 p2 c0 {1,S} {2,S}
"""),
    E0 = (221.344,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (47.9982,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.68624,-0.0442492,0.000174348,-2.12652e-07,8.37384e-11,26624.3,8.96412], Tmin=(100,'K'), Tmax=(903.391,'K')), NASAPolynomial(coeffs=[16.9701,-0.0234422,1.49415e-05,-2.87692e-09,1.87776e-13,21336.4,-66.0334], Tmin=(903.391,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(221.344,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsOs) + group(O2s-OsOs) + group(O2s-OsOs) + ring(Cyclopropane)"""),
)

species(
    label = '[O]C(F)F(170)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u1 p2 c0 {4,S}
4 C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
5 H u0 p0 c0 {4,S}
"""),
    E0 = (-426.241,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(66.9995,'amu')),
        NonlinearRotor(inertia=([47.2408,48.1768,88.2639],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([480.226,502.831,615.399,942.972,1118.44,1159.86,1274.53,1324.03,2842.12],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.0148,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.07235,-0.00674902,8.67846e-05,-1.53143e-07,8.6196e-11,-51263.6,9.19392], Tmin=(10,'K'), Tmax=(577.889,'K')), NASAPolynomial(coeffs=[2.80054,0.016572,-1.1432e-05,3.6345e-09,-4.33761e-13,-51359,12.5348], Tmin=(577.889,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-426.241,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""[O]C(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]OOO[CH]F(139)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {3,S} {6,S}
3 O u0 p2 c0 {2,S} {4,S}
4 O u0 p2 c0 {3,S} {5,S}
5 O u1 p2 c0 {4,S}
6 C u1 p0 c0 {1,S} {2,S} {7,S}
7 H u0 p0 c0 {6,S}
"""),
    E0 = (61.1616,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([580,1155,1237,1373,3147,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0146,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90205,0.0381398,-4.23237e-05,2.19275e-08,-4.22941e-12,7438.42,25.1123], Tmin=(100,'K'), Tmax=(1436.66,'K')), NASAPolynomial(coeffs=[13.034,0.00233038,9.20934e-08,-8.80929e-11,7.61256e-15,4736.79,-30.897], Tmin=(1436.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(61.1616,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(145.503,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsOs) + group(O2s-OsH) + group(CsFHHO) + radical(ROOJ) + radical(CsF1sHO2s)"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.47471,0.00468947,-6.33886e-06,3.99267e-09,-7.67832e-13,23182.7,11.3223], Tmin=(100,'K'), Tmax=(1873.31,'K')), NASAPolynomial(coeffs=[2.72432,0.000671973,1.37793e-06,-3.54951e-10,2.60883e-14,24449.9,18.0454], Tmin=(1873.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.544,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsOs) + group(O2s-OsH) + group(O2s-OsH) + radical(ROOJ) + radical(ROOJ)"""),
)

species(
    label = '[O]OO[O](224)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {3,S}
2 O u0 p2 c0 {1,S} {4,S}
3 O u1 p2 c0 {1,S}
4 O u1 p2 c0 {2,S}
"""),
    E0 = (228.521,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([584.658,585.157,593.212,3787.6,3789.38],'cm^-1')),
        HinderedRotor(inertia=(0.0036014,'amu*angstrom^2'), symmetry=1, barrier=(0.878668,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (63.9976,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.48985,0.0151896,-1.78846e-05,9.27154e-09,-1.58491e-12,27555.2,19.6951], Tmin=(100,'K'), Tmax=(1976.17,'K')), NASAPolynomial(coeffs=[3.91599,0.000510603,2.20836e-06,-5.2646e-10,3.66138e-14,29294.2,17.6693], Tmin=(1976.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(228.521,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(78.9875,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsOs) + group(O2s-OsOs) + group(O2s-OsH) + group(O2s-OsH) + radical(ROOJ) + radical(ROOJ)"""),
)

species(
    label = 'CHF2(81)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u1 p0 c0 {1,S} {2,S} {4,S}
4 H u0 p0 c0 {3,S}
"""),
    E0 = (-256.71,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(51.0046,'amu')),
        NonlinearRotor(inertia=([7.43413,45.9439,52.5803],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([549.125,1005.77,1195.1,1212.61,1359.42,3085.19],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (51.0154,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2178.39,'J/mol'), sigma=(4.123,'angstroms'), dipoleMoment=(1.8,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05476,-0.0040567,3.90133e-05,-5.51349e-08,2.50461e-11,-30875.2,7.58714], Tmin=(10,'K'), Tmax=(697.139,'K')), NASAPolynomial(coeffs=[2.58942,0.0108145,-6.89144e-06,2.06262e-09,-2.34597e-13,-30827.9,13.0014], Tmin=(697.139,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-256.71,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""F[CH]F""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]OOO[C](F)F(225)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {4,S} {7,S}
4 O u0 p2 c0 {3,S} {5,S}
5 O u0 p2 c0 {4,S} {6,S}
6 O u1 p2 c0 {5,S}
7 C u1 p0 c0 {1,S} {2,S} {3,S}
"""),
    E0 = (-152.325,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([493,600,700,1144,1293,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44467,0.0442207,-5.02317e-05,2.57938e-08,-4.84236e-12,-18218.2,27.9691], Tmin=(100,'K'), Tmax=(1525.91,'K')), NASAPolynomial(coeffs=[15.1759,-7.14584e-05,1.4646e-06,-3.55603e-10,2.56916e-14,-21442.8,-40.9138], Tmin=(1525.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-152.325,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(145.503,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsOs) + group(O2s-OsH) + group(CsFFHO) + radical(ROOJ) + radical(Csj(F1s)(F1s)(O2s-O2s))"""),
)

species(
    label = 'OOOO[C](F)F(226)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {4,S} {7,S}
4 O u0 p2 c0 {3,S} {5,S}
5 O u0 p2 c0 {4,S} {6,S}
6 O u0 p2 c0 {5,S} {8,S}
7 C u1 p0 c0 {1,S} {2,S} {3,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-304.33,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([493,600,700,1144,1293,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.013,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63942,0.0424687,-2.70941e-05,-7.13849e-09,8.90263e-12,-36508.8,26.5071], Tmin=(100,'K'), Tmax=(941.739,'K')), NASAPolynomial(coeffs=[16.7313,0.00182747,2.70754e-07,-5.68561e-11,1.72327e-16,-40391.7,-50.9201], Tmin=(941.739,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-304.33,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(166.289,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsOs) + group(O2s-OsH) + group(CsFFHO) + radical(Csj(F1s)(F1s)(O2s-O2s))"""),
)

species(
    label = 'F[CH]OOOOF(227)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {4,S} {7,S}
4 O u0 p2 c0 {3,S} {5,S}
5 O u0 p2 c0 {4,S} {6,S}
6 O u0 p2 c0 {2,S} {5,S}
7 C u1 p0 c0 {1,S} {3,S} {8,S}
8 H u0 p0 c0 {7,S}
"""),
    E0 = (4.45752,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([277,555,632,580,1155,1237,1373,3147,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.013,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50501,0.0487395,-5.10874e-05,2.27409e-08,-3.15801e-12,631.327,27.1588], Tmin=(100,'K'), Tmax=(1027.27,'K')), NASAPolynomial(coeffs=[15.2654,0.00516455,-2.06989e-06,4.11363e-10,-3.11902e-14,-2723.75,-42.1672], Tmin=(1027.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(4.45752,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(166.289,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsOs) + group(O2sFO) + group(CsFHHO) + radical(CsF1sHO2s)"""),
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
    E0 = (-306.264,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (5.41093,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (385.495,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-100.407,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-142.022,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (196.928,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-170.822,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (34.6856,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (122.354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-212.946,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (128.258,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O2(2)', '[O]OC(F)F(168)'],
    products = ['[O]OOOC(F)F(177)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(4.98302e+08,'m^3/(mol*s)'), n=0.0472428, Ea=(67.442,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]OOO(169)', 'CF2(42)'],
    products = ['[O]OOOC(F)F(177)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.76395e-12,'m^3/(mol*s)'), n=5.02686, Ea=(69.7321,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='OH_N-2Br1sCl1sF1sHI1s->H_2Br1sCl1sF1sI1s->F1s',), comment="""Estimated from node OH_N-2Br1sCl1sF1sHI1s->H_2Br1sCl1sF1sI1s->F1s"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]OOOF(221)', 'CHF(40)'],
    products = ['[O]OOOC(F)F(177)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.13916e+09,'m^3/(mol*s)'), n=-1.00842, Ea=(12.0472,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17406617270594374, var=16.167420070870673, Tref=1000.0, N=4, data_mean=0.0, correlation='OHOY',), comment="""Estimated from node OHOY"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(6)', '[O]OOC(F)F(222)'],
    products = ['[O]OOOC(F)F(177)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]OOOC(F)F(177)'],
    products = ['[O-][O+]=O(223)', '[O]C(F)F(170)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.88184e+13,'s^-1'), n=-0.182188, Ea=(165.442,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO;Y_rad_intra;OOR]
Euclidian distance = 0
family: Cyclic_Ether_Formation
Ea raised from 164.1 to 165.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction6',
    reactants = ['F(37)', '[O]OOO[CH]F(139)'],
    products = ['[O]OOOC(F)F(177)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.8422e+07,'m^3/(mol*s)'), n=-0.361029, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H_Ext-2C-R"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]O[O](152)', '[O]C(F)F(170)'],
    products = ['[O]OOOC(F)F(177)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.98302e+08,'m^3/(mol*s)'), n=0.0472428, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]OO[O](224)', 'CHF2(81)'],
    products = ['[O]OOOC(F)F(177)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.01e+06,'m^3/(mol*s)'), n=-1.3397e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Sp-4R!H-2C_Ext-2C-R_N-5R!H->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Sp-4R!H-2C_Ext-2C-R_N-5R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(5)', '[O]OOO[C](F)F(225)'],
    products = ['[O]OOOC(F)F(177)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(43950,'m^3/(mol*s)'), n=1, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_N-3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction10',
    reactants = ['OOOO[C](F)F(226)'],
    products = ['[O]OOOC(F)F(177)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2044.71,'s^-1'), n=2.0375, Ea=(28.5087,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;C_rad_out_single;XH_out] for rate rule [R5H_SSSS;C_rad_out_noH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F[CH]OOOOF(227)'],
    products = ['[O]OOOC(F)F(177)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(60.9255,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

network(
    label = 'PDepNetwork #36',
    isomers = [
        '[O]OOOC(F)F(177)',
    ],
    reactants = [
        ('O2(2)', '[O]OC(F)F(168)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #36',
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

