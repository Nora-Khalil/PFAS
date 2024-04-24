species(
    label = 'OC[C](F)F(1002)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {4,S} {8,S}
4 C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
5 C u1 p0 c0 {1,S} {2,S} {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (-455.869,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,190,488,555,1236,1407,180],'cm^-1')),
        HinderedRotor(inertia=(0.173691,'amu*angstrom^2'), symmetry=1, barrier=(3.99349,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.173679,'amu*angstrom^2'), symmetry=1, barrier=(3.99321,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0414,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.85138,0.0182191,0.000196015,-1.45993e-06,3.09104e-09,-54828.1,10.7199], Tmin=(10,'K'), Tmax=(169.558,'K')), NASAPolynomial(coeffs=[4.37639,0.0240011,-1.58538e-05,5.00357e-09,-6.02456e-13,-54872.1,8.34856], Tmin=(169.558,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-455.869,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""OC[C](F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'CH2CF2(56)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {4,D} {5,S} {6,S}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-361.616,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(64.0125,'amu')),
        NonlinearRotor(inertia=([45.7027,48.2614,93.9642],'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([437.293,557.015,653.832,726.079,816.319,956,966.438,1345.56,1413.22,1792.31,3202.97,3303.55],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2091.09,'J/mol'), sigma=(4.442,'angstroms'), dipoleMoment=(1.4,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.10281,-0.0101072,0.000121983,-2.28108e-07,1.37933e-10,-43490.6,7.77929], Tmin=(10,'K'), Tmax=(534.293,'K')), NASAPolynomial(coeffs=[2.52167,0.0198841,-1.31824e-05,4.13929e-09,-4.93215e-13,-43580.8,11.9914], Tmin=(534.293,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-361.616,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""CDC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'OC[C]F(1878)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {7,S}
3 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
4 C u2 p0 c0 {1,S} {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {2,S}
"""),
    E0 = (-65.6582,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,513.312,513.312,513.338,513.349,513.354,4000],'cm^-1')),
        HinderedRotor(inertia=(0.000639776,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.107457,'amu*angstrom^2'), symmetry=1, barrier=(20.0962,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (62.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.45168,0.0281351,-1.74296e-05,-1.76808e-09,3.60003e-12,-7835.89,12.9138], Tmin=(100,'K'), Tmax=(1007.81,'K')), NASAPolynomial(coeffs=[11.4873,0.00495459,-1.80386e-06,3.81704e-10,-3.06105e-14,-10301.2,-33.9435], Tmin=(1007.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-65.6582,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCCl_triplet)"""),
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
    label = 'CH2OH(81)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,S} {5,S}
2 C u1 p0 c0 {1,S} {3,S} {4,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {1,S}
"""),
    E0 = (-28.7184,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3003.59,4000],'cm^-1')),
        HinderedRotor(inertia=(0.057913,'amu*angstrom^2'), symmetry=1, barrier=(25.9304,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (31.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.47834,-0.0013507,2.78485e-05,-3.64869e-08,1.47907e-11,-3500.73,3.30913], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.09314,0.00594761,-2.06497e-06,3.23008e-10,-1.88126e-14,-4034.1,-1.84691], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-28.7184,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""CH2OH""", comment="""Thermo library: FFCM1(-)"""),
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
    label = 'OC=C(F)F(1217)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {4,S} {7,S}
4 C u0 p0 c0 {3,S} {5,D} {6,S}
5 C u0 p0 c0 {1,S} {2,S} {4,D}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (-502.96,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413],'cm^-1')),
        HinderedRotor(inertia=(0.552091,'amu*angstrom^2'), symmetry=1, barrier=(12.6937,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.85967,0.00956102,9.95844e-05,-2.79091e-07,2.18422e-10,-60486.4,9.26794], Tmin=(10,'K'), Tmax=(463.499,'K')), NASAPolynomial(coeffs=[6.00669,0.0174593,-1.1501e-05,3.70008e-09,-4.58213e-13,-60969.3,-2.50148], Tmin=(463.499,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-502.96,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""OCDC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[CH2][C](F)F(185)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 C u1 p0 c0 {4,S} {5,S} {6,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-109.048,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,190,488,555,1236,1407],'cm^-1')),
        HinderedRotor(inertia=(0.00258864,'amu*angstrom^2'), symmetry=1, barrier=(7.63529,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.25819,0.0186259,-1.56413e-05,8.223e-09,-2.04577e-12,-13090.7,13.5552], Tmin=(100,'K'), Tmax=(887.666,'K')), NASAPolynomial(coeffs=[4.4701,0.0131647,-6.4128e-06,1.29206e-09,-9.37476e-14,-13305.9,7.85287], Tmin=(887.666,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-109.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Cs_P) + radical(CsCsF1sF1s)"""),
)

species(
    label = '[O]C[C](F)F(728)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u1 p2 c0 {4,S}
4 C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
5 C u1 p0 c0 {1,S} {2,S} {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-233.591,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,190,488,555,1236,1407,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0983613,'amu*angstrom^2'), symmetry=1, barrier=(2.26152,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.3035,0.0417665,-6.88725e-05,6.30488e-08,-2.27983e-11,-28037.7,16.0779], Tmin=(100,'K'), Tmax=(796.61,'K')), NASAPolynomial(coeffs=[5.91256,0.0167111,-8.63847e-06,1.7145e-09,-1.20926e-13,-28392.7,0.867685], Tmin=(796.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-233.591,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCOJ) + radical(Csj(Cs-O2sHH)(F1s)(F1s))"""),
)

species(
    label = 'O[CH][C](F)F(1445)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {4,S} {7,S}
4 C u1 p0 c0 {3,S} {5,S} {6,S}
5 C u1 p0 c0 {1,S} {2,S} {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (-278.999,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,190,488,555,1236,1407,180],'cm^-1')),
        HinderedRotor(inertia=(0.636525,'amu*angstrom^2'), symmetry=1, barrier=(14.635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.636282,'amu*angstrom^2'), symmetry=1, barrier=(14.6294,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72308,0.0525948,-8.51294e-05,6.61096e-08,-1.96969e-11,-33476.2,16.8524], Tmin=(100,'K'), Tmax=(832.861,'K')), NASAPolynomial(coeffs=[11.1819,0.00716935,-3.32205e-06,6.30394e-10,-4.3141e-14,-35051.8,-27.0516], Tmin=(832.861,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-278.999,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCsJOH) + radical(Csj(Cs-O2sHH)(F1s)(F1s))"""),
)

species(
    label = 'OC[C]F-2(1879)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {7,S}
3 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
4 C u0 p1 c0 {1,S} {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {2,S}
"""),
    E0 = (-78.9184,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,617,898,1187],'cm^-1')),
        HinderedRotor(inertia=(0.32109,'amu*angstrom^2'), symmetry=1, barrier=(7.38249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.321048,'amu*angstrom^2'), symmetry=1, barrier=(7.38152,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (62.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.88071,0.0254385,-2.53338e-05,1.4826e-08,-3.65566e-12,-9452.08,14.087], Tmin=(100,'K'), Tmax=(964.949,'K')), NASAPolynomial(coeffs=[5.99226,0.0125404,-5.28411e-06,9.74209e-10,-6.69774e-14,-10052.6,-0.813474], Tmin=(964.949,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-78.9184,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsOsHH) + group(CJ2_singlet-FCs)"""),
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
    label = 'O[CH]C(F)F(1880)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {5,S} {8,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u1 p0 c0 {3,S} {4,S} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (-481.391,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (81.0414,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.89233,0.00668522,0.000104616,-2.31512e-07,1.50155e-10,-57892.5,11.6044], Tmin=(10,'K'), Tmax=(527.066,'K')), NASAPolynomial(coeffs=[3.95879,0.0284142,-2.04983e-05,6.77504e-09,-8.32891e-13,-58208.3,8.39674], Tmin=(527.066,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-481.391,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""O[CH]C(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]CC(F)F(1554)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u1 p2 c0 {5,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-436.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([235,523,627,1123,1142,1372,1406,3097,2750,2850,1437.5,1250,1305,750,350,295.454,295.928],'cm^-1')),
        HinderedRotor(inertia=(0.00197665,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0414,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.86274,0.0109615,6.89758e-05,-1.67737e-07,1.18262e-10,-52460.1,11.1968], Tmin=(10,'K'), Tmax=(469.189,'K')), NASAPolynomial(coeffs=[3.25175,0.0271511,-1.78876e-05,5.56681e-09,-6.58859e-13,-52523.6,12.3941], Tmin=(469.189,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-436.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""[O]CC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    E0 = (-227.718,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (112.737,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (110.512,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-185.652,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (24.8499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (83.7171,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (40.9085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (99.4766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-124.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-338.257,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-330.685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['OH(4)', 'CH2CF2(56)'],
    products = ['OC[C](F)F(1002)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(565.393,'m^3/(mol*s)'), n=1.28408, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.06277471317878429, var=1.0052381843076323, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_N-3BrFOS->F_Ext-2CS-R',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_N-3BrFOS->F_Ext-2CS-R"""),
)

reaction(
    label = 'reaction2',
    reactants = ['F(37)', 'OC[C]F(1878)'],
    products = ['OC[C](F)F(1002)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [H/Val7_rad;Birad] for rate rule [Val7_rad;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CF2(42)', 'CH2OH(81)'],
    products = ['OC[C](F)F(1002)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/O;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(5)', 'OC=C(F)F(1217)'],
    products = ['OC[C](F)F(1002)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7.09572,'m^3/(mol*s)'), n=2.03221, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.053734944950677085, var=3.1017726857793777, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_N-4R!H->C',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_N-4R!H->C"""),
)

reaction(
    label = 'reaction5',
    reactants = ['OH(4)', '[CH2][C](F)F(185)'],
    products = ['OC[C](F)F(1002)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7.7e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R
Ea raised from -1.7 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(5)', '[O]C[C](F)F(728)'],
    products = ['OC[C](F)F(1002)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(5)', 'O[CH][C](F)F(1445)'],
    products = ['OC[C](F)F(1002)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(8.43711e+23,'m^3/(mol*s)'), n=-5.99271, Ea=(2.59921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25884854535692575, var=10.966526674850428, Tref=1000.0, N=9, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F(37)', 'OC[C]F-2(1879)'],
    products = ['OC[C](F)F(1002)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(518506,'m^3/(mol*s)'), n=0.4717, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->H_N-2Br1sCl1sF1s->Cl1s_3BrCClFINOPSSi->F',), comment="""Estimated from node Root_N-3R->H_N-2Br1sCl1sF1s->Cl1s_3BrCClFINOPSSi->F"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CF2(87)', 'CH2OH(81)'],
    products = ['OC[C](F)F(1002)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(2.66112,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction10',
    reactants = ['OC[C](F)F(1002)'],
    products = ['O[CH]C(F)F(1880)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.245e+06,'s^-1'), n=1.223, Ea=(12.1085,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 442 used for R2H_S;C_rad_out_noH;Cs_H_out_H/NonDeO
Exact match found for rate rule [R2H_S;C_rad_out_noH;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['OC[C](F)F(1002)'],
    products = ['[O]CC(F)F(1554)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5.248e+06,'s^-1'), n=0.992, Ea=(19.6812,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 441 used for R3H_SS_Cs;C_rad_out_noH;O_H_out
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_noH;O_H_out]
Euclidian distance = 0
family: intra_H_migration
Ea raised from 17.8 to 19.7 kJ/mol to match endothermicity of reaction."""),
)

network(
    label = 'PDepNetwork #420',
    isomers = [
        'OC[C](F)F(1002)',
    ],
    reactants = [
        ('OH(4)', 'CH2CF2(56)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #420',
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

