species(
    label = 'FC(F)OOOOC(F)F(346)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {8,S} {10,S}
7  O u0 p2 c0 {5,S} {8,S}
8  O u0 p2 c0 {6,S} {7,S}
9  C u0 p0 c0 {1,S} {2,S} {5,S} {11,S}
10 C u0 p0 c0 {3,S} {4,S} {6,S} {12,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {10,S}
"""),
    E0 = (-969.199,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([480,558,523,663,745,915,1086,1144,1147,1185,1338,1440,1394,1432,3056,3172,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (166.028,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.380488,0.0764445,-8.27139e-05,4.32521e-08,-8.84679e-12,-116435,35.3444], Tmin=(100,'K'), Tmax=(1191.23,'K')), NASAPolynomial(coeffs=[17.9501,0.0174468,-8.42219e-06,1.67417e-09,-1.20768e-13,-120620,-52.4928], Tmin=(1191.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-969.199,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsOs) + group(CsFFHO) + group(CsFFHO)"""),
)

species(
    label = '[O]OC(F)F(128)',
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
    label = 'OOOOC(F)F(327)',
    structure = adjacencyList("""1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {4,S} {7,S}
4 O u0 p2 c0 {3,S} {5,S}
5 O u0 p2 c0 {4,S} {6,S}
6 O u0 p2 c0 {5,S} {9,S}
7 C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
8 H u0 p0 c0 {7,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-522.344,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([519,593,830,1115,1166,1389,1413,3114,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (116.021,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3517.11,'J/mol'), sigma=(5.59388,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=549.36 K, Pc=45.59 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55601,0.0451149,-3.15334e-05,-1.91076e-10,5.47313e-12,-62727.7,26.3089], Tmin=(100,'K'), Tmax=(971.768,'K')), NASAPolynomial(coeffs=[15.739,0.00638143,-2.07153e-06,4.01892e-10,-3.1738e-14,-66411.8,-46.4827], Tmin=(971.768,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-522.344,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(191.233,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsOs) + group(O2s-OsH) + group(CsFFHO)"""),
)

species(
    label = 'CHF(88)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 C u0 p1 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (153.289,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1169.21,1416.41,2978.06],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (32.017,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2178.39,'J/mol'), sigma=(4.123,'angstroms'), dipoleMoment=(1.8,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.34484,0.00235461,1.93983e-06,-2.65251e-09,7.91169e-13,18514,7.05286], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[4.48366,0.00174964,-5.0479e-07,1.08953e-10,-9.87898e-15,17958.1,0.289222], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(153.289,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CHF""", comment="""Thermo library: Fluorine"""),
)

species(
    label = 'FOOOOC(F)F(412)',
    structure = adjacencyList("""1 F u0 p3 c0 {8,S}
2 F u0 p3 c0 {8,S}
3 F u0 p3 c0 {7,S}
4 O u0 p2 c0 {5,S} {8,S}
5 O u0 p2 c0 {4,S} {6,S}
6 O u0 p2 c0 {5,S} {7,S}
7 O u0 p2 c0 {3,S} {6,S}
8 C u0 p0 c0 {1,S} {2,S} {4,S} {9,S}
9 H u0 p0 c0 {8,S}
"""),
    E0 = (-427.043,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([277,555,632,519,593,830,1115,1166,1389,1413,3114,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (134.011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15584,0.055234,-5.57187e-05,2.36167e-08,-2.96841e-12,-51252.6,29.1278], Tmin=(100,'K'), Tmax=(1037.97,'K')), NASAPolynomial(coeffs=[16.7402,0.00691965,-2.86795e-06,5.71106e-10,-4.29381e-14,-55120.4,-49.6856], Tmin=(1037.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-427.043,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(191.233,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsOs) + group(O2sFO) + group(CsFFHO)"""),
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
    label = 'F[CH]OOOOC(F)F(413)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {6,S} {8,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {4,S} {7,S}
7  O u0 p2 c0 {5,S} {6,S}
8  C u0 p0 c0 {1,S} {2,S} {4,S} {10,S}
9  C u1 p0 c0 {3,S} {5,S} {11,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-537.698,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([519,593,830,1115,1166,1389,1413,3114,580,1155,1237,1373,3147,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (147.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.738245,0.0698417,-7.76872e-05,4.18689e-08,-8.83217e-12,-64551.1,33.3453], Tmin=(100,'K'), Tmax=(1155.88,'K')), NASAPolynomial(coeffs=[16.3026,0.015979,-7.7874e-06,1.55259e-09,-1.12162e-13,-68149.1,-43.9982], Tmin=(1155.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-537.698,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(236.962,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsOs) + group(CsFFHO) + group(CsFHHO) + radical(CsF1sHO2s)"""),
)

species(
    label = '[O]C(F)F(146)',
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3063.56,'J/mol'), sigma=(4.99758,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=478.52 K, Pc=55.69 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.07235,-0.00674902,8.67846e-05,-1.53143e-07,8.6196e-11,-51263.6,9.19392], Tmin=(10,'K'), Tmax=(577.889,'K')), NASAPolynomial(coeffs=[2.80054,0.016572,-1.1432e-05,3.6345e-09,-4.33761e-13,-51359,12.5348], Tmin=(577.889,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-426.241,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""[O]C(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]OOC(F)F(226)',
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
        HarmonicOscillator(frequencies=([519,593,830,1115,1166,1389,1413,3114,371.044,371.074,371.083,371.09,1852.3],'cm^-1')),
        HinderedRotor(inertia=(0.00122402,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0886911,'amu*angstrom^2'), symmetry=1, barrier=(8.66488,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.0136,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.30892,0.0364485,-4.19925e-05,2.43451e-08,-5.53864e-12,-48806.9,20.2522], Tmin=(100,'K'), Tmax=(1076.17,'K')), NASAPolynomial(coeffs=[9.6708,0.0090853,-3.85275e-06,7.18309e-10,-5.00065e-14,-50391.5,-15.805], Tmin=(1076.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-406.316,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsH) + group(CsFFHO) + radical(ROOJ)"""),
)

species(
    label = 'CHF2(81)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u1 p0 c0 {1,S} {2,S} {4,S}
4 H u0 p0 c0 {3,S}
"""),
    E0 = (-258.033,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83276,0.0134164,-9.54328e-06,2.9593e-09,-2.792e-13,-30857.1,16.7589], Tmin=(298,'K'), Tmax=(1150,'K')), NASAPolynomial(coeffs=[5.44828,0.00441165,-1.74189e-06,3.092e-10,-2.01762e-14,-31961,-2.29455], Tmin=(1150,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-258.032,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(83.1447,'J/mol/K'), label="""CHF2""", comment="""Thermo library: Fluorine"""),
)

species(
    label = '[O]OOOC(F)F(206)',
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3517.11,'J/mol'), sigma=(5.59388,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=549.36 K, Pc=45.59 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51475,0.0450478,-4.82382e-05,2.42478e-08,-4.56695e-12,-44443.8,27.2206], Tmin=(100,'K'), Tmax=(1452.69,'K')), NASAPolynomial(coeffs=[14.7734,0.00366623,-4.76548e-07,1.95581e-11,5.59098e-17,-47781.7,-39.9255], Tmin=(1452.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-370.339,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsOs) + group(O2s-OsH) + group(CsFFHO) + radical(ROOJ)"""),
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
    label = 'F[C](F)OOOOC(F)F(414)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {10,S}
4  F u0 p3 c0 {10,S}
5  O u0 p2 c0 {7,S} {9,S}
6  O u0 p2 c0 {8,S} {10,S}
7  O u0 p2 c0 {5,S} {8,S}
8  O u0 p2 c0 {6,S} {7,S}
9  C u0 p0 c0 {1,S} {2,S} {5,S} {11,S}
10 C u1 p0 c0 {3,S} {4,S} {6,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-751.185,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([519,593,830,1115,1166,1389,1413,3114,493,600,700,1144,1293,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (165.02,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.299367,0.075743,-8.51053e-05,4.5229e-08,-9.26533e-12,-90208.7,36.1326], Tmin=(100,'K'), Tmax=(1200.74,'K')), NASAPolynomial(coeffs=[19.3915,0.0121414,-5.65195e-06,1.11529e-09,-8.06136e-14,-94793.6,-59.4688], Tmin=(1200.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-751.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(236.962,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-OsOs) + group(CsFFHO) + group(CsFFHO) + radical(Csj(F1s)(F1s)(O2s-O2s))"""),
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
    E0 = (-524.969,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-320.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (68.1156,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-133.867,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-501.618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-297.432,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-208.441,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]OC(F)F(128)', '[O]OC(F)F(128)'],
    products = ['FC(F)OOOOC(F)F(346)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.24576e+08,'m^3/(mol*s)'), n=0.0472428, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CF2(89)', 'OOOOC(F)F(327)'],
    products = ['FC(F)OOOOC(F)F(346)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.76395e-12,'m^3/(mol*s)'), n=5.02686, Ea=(67.7731,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='OH_N-2Br1sCl1sF1sHI1s->H_2Br1sCl1sF1sI1s->F1s',), comment="""Estimated from node OH_N-2Br1sCl1sF1sHI1s->H_2Br1sCl1sF1sI1s->F1s"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CHF(88)', 'FOOOOC(F)F(412)'],
    products = ['FC(F)OOOOC(F)F(346)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.13916e+09,'m^3/(mol*s)'), n=-1.00842, Ea=(10.9307,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17406617270594374, var=16.167420070870673, Tref=1000.0, N=4, data_mean=0.0, correlation='OHOY',), comment="""Estimated from node OHOY"""),
)

reaction(
    label = 'reaction4',
    reactants = ['F(37)', 'F[CH]OOOOC(F)F(413)'],
    products = ['FC(F)OOOOC(F)F(346)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]C(F)F(146)', '[O]OOC(F)F(226)'],
    products = ['FC(F)OOOOC(F)F(346)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.49151e+08,'m^3/(mol*s)'), n=0.0472428, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_N-2R->C"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CHF2(81)', '[O]OOOC(F)F(206)'],
    products = ['FC(F)OOOOC(F)F(346)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.505e+06,'m^3/(mol*s)'), n=-1.3397e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Sp-4R!H-2C_Ext-2C-R_N-5R!H->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Sp-4R!H-2C_Ext-2C-R_N-5R!H->C"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(5)', 'F[C](F)OOOOC(F)F(414)'],
    products = ['FC(F)OOOOC(F)F(346)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(4.1766,'m^3/(mol*s)'), n=1.94174, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_N-4R!H->C_Ext-2C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_N-4R!H->C_Ext-2C-R
Ea raised from -2.7 to 0.0 kJ/mol."""),
)

network(
    label = 'PDepNetwork #151',
    isomers = [
        'FC(F)OOOOC(F)F(346)',
    ],
    reactants = [
        ('[O]OC(F)F(128)', '[O]OC(F)F(128)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #151',
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

