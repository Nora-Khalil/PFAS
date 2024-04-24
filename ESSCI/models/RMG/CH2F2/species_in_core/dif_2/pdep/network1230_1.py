species(
    label = '[O]C(=O)C[C]=O(3274)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {5,S}
2 O u0 p2 c0 {5,D}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
5 C u0 p0 c0 {1,S} {2,D} {4,S}
6 C u1 p0 c0 {3,D} {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-152.809,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,180,180,180,180,1436.45,1436.69],'cm^-1')),
        HinderedRotor(inertia=(0.115274,'amu*angstrom^2'), symmetry=1, barrier=(2.65037,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.03809,'amu*angstrom^2'), symmetry=1, barrier=(55.916,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0461,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22207,0.0472258,-8.40982e-05,8.54789e-08,-3.29542e-11,-18322.6,18.797], Tmin=(100,'K'), Tmax=(836.379,'K')), NASAPolynomial(coeffs=[2.538,0.0269226,-1.39828e-05,2.72688e-09,-1.89182e-13,-17718.2,21.2586], Tmin=(836.379,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-152.809,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)(Cds-O2d)HH) + group(Cds-OdCsOs) + group(Cds-OdCsH) + radical(CCOJ) + radical(C=OCCJ=O)"""),
)

species(
    label = 'CO2(13)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,D}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,D} {2,D}
"""),
    E0 = (-403.138,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([664.558,664.681,1368.23,2264.78],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0094,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(2028.74,'J/mol'), sigma=(3.763,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(2.65,'angstroms^3'), rotrelaxcollnum=2.1, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.35681,0.00898413,-7.12206e-06,2.4573e-09,-1.42885e-13,-48372,9.9009], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.63651,0.00274146,-9.95898e-07,1.60387e-10,-9.16199e-15,-49024.9,-1.9349], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-403.138,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), label="""CO2""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'CH2CO(75)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,D}
2 C u0 p0 c0 {3,D} {4,S} {5,S}
3 C u0 p0 c0 {1,D} {2,D}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (-60.8183,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2120,512.5,787.5],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0366,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13241,0.0181319,-1.74093e-05,9.35336e-09,-2.01725e-12,-7148.09,13.3808], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.75871,0.00635124,-2.25955e-06,3.62322e-10,-2.15856e-14,-8085.33,-4.9649], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-60.8183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""CH2CO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[CH2]C(=O)C([O])=O(3272)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {4,D}
2 O u1 p2 c0 {6,S}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {1,D} {5,S} {6,S}
5 C u1 p0 c0 {4,S} {7,S} {8,S}
6 C u0 p0 c0 {2,S} {3,D} {4,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-71.7875,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,3000,3100,440,815,1455,1000,180,180,180,912.956,913.782,914.459],'cm^-1')),
        HinderedRotor(inertia=(0.000216936,'amu*angstrom^2'), symmetry=1, barrier=(2.4631,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.26952,'amu*angstrom^2'), symmetry=1, barrier=(29.1887,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0461,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.01751,0.0479637,-6.98526e-05,5.60541e-08,-1.84278e-11,-8566.69,20.8696], Tmin=(100,'K'), Tmax=(738.778,'K')), NASAPolynomial(coeffs=[7.46474,0.0184721,-9.9768e-06,2.02576e-09,-1.45833e-13,-9371.59,-3.76124], Tmin=(738.778,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-71.7875,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)O2s) + radical(C=OC=OOJ) + radical(CJCC=O)"""),
)

species(
    label = '[C]=O(222)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,D}
2 C u2 p0 c0 {1,D}
"""),
    E0 = (439.086,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3054.48],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.01,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.08916,0.0020042,-1.61663e-05,2.55061e-08,-1.16425e-11,52802.7,4.52506], Tmin=(100,'K'), Tmax=(856.108,'K')), NASAPolynomial(coeffs=[0.961633,0.00569044,-3.48043e-06,7.192e-10,-5.0804e-14,53738.7,21.4663], Tmin=(856.108,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(439.086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH2]C([O])=O(882)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {4,S}
2 O u0 p2 c0 {4,D}
3 C u1 p0 c0 {4,S} {5,S} {6,S}
4 C u0 p0 c0 {1,S} {2,D} {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-8.20105,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,840.297,843.048,844.532,847.195,847.24],'cm^-1')),
        HinderedRotor(inertia=(0.00482699,'amu*angstrom^2'), symmetry=1, barrier=(2.45134,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.60873,0.0103835,4.16408e-06,-6.72691e-09,1.73623e-12,-973.928,12.8424], Tmin=(100,'K'), Tmax=(1564.28,'K')), NASAPolynomial(coeffs=[5.70751,0.0133436,-6.65904e-06,1.28862e-09,-8.86359e-14,-2649.33,-1.47843], Tmin=(1564.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.20105,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCOJ) + radical(CJCO)"""),
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
    label = 'O=[C]C[C]=O(3399)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {4,D}
2 O u0 p2 c0 {5,D}
3 C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
4 C u1 p0 c0 {1,D} {3,S}
5 C u1 p0 c0 {2,D} {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (47.1985,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1850,1860,440,470,900,1000],'cm^-1')),
        HinderedRotor(inertia=(0.556553,'amu*angstrom^2'), symmetry=1, barrier=(12.7962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.555979,'amu*angstrom^2'), symmetry=1, barrier=(12.7831,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0467,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96236,0.0530439,-0.000104387,9.76059e-08,-3.34749e-11,5742.3,16.165], Tmin=(100,'K'), Tmax=(905.3,'K')), NASAPolynomial(coeffs=[6.16616,0.0143356,-6.89086e-06,1.24225e-09,-7.98356e-14,5806.23,0.859136], Tmin=(905.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(47.1985,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-O2d)HH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(C=OCCJ=O) + radical(C=OCCJ=O)"""),
)

species(
    label = 'O=C1CC(=O)O1(3278)',
    structure = adjacencyList("""1 O u0 p2 c0 {5,S} {6,S}
2 O u0 p2 c0 {5,D}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
5 C u0 p0 c0 {1,S} {2,D} {4,S}
6 C u0 p0 c0 {1,S} {3,D} {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-404.475,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0461,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.07139,0.0106757,3.47794e-05,-5.07218e-08,1.92855e-11,-48605,17.498], Tmin=(100,'K'), Tmax=(989.816,'K')), NASAPolynomial(coeffs=[8.67131,0.0124785,-4.97923e-06,9.95048e-10,-7.54418e-14,-50910.5,-15.5071], Tmin=(989.816,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-404.475,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-O2d)) + group(Cs-(Cds-O2d)(Cds-O2d)HH) + group(Cds-OdCsOs) + group(Cds-OdCsOs) + ring(Cyclobutane)"""),
)

species(
    label = 'O=C=CC(=O)O(3271)',
    structure = adjacencyList("""1 O u0 p2 c0 {5,S} {8,S}
2 O u0 p2 c0 {5,D}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {5,S} {6,D} {7,S}
5 C u0 p0 c0 {1,S} {2,D} {4,S}
6 C u0 p0 c0 {3,D} {4,D}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {1,S}
"""),
    E0 = (-335.056,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0461,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2527,0.0440707,-7.40444e-05,7.03395e-08,-2.60598e-11,-40240.5,19.6058], Tmin=(100,'K'), Tmax=(819.479,'K')), NASAPolynomial(coeffs=[4.7508,0.0206648,-1.06779e-05,2.0929e-09,-1.46082e-13,-40273.4,10.3484], Tmin=(819.479,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-335.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + missing(O2d-Cdd) + group(Cds-(Cdd-O2d)CsH) + group(Cds-O2d(Cds-Cds)O2s) + missing(Cdd-CdO2d)"""),
)

species(
    label = 'O=[C]C[C]1OO1(3400)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {5,S}
2 O u0 p2 c0 {1,S} {5,S}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
5 C u1 p0 c0 {1,S} {2,S} {4,S}
6 C u1 p0 c0 {3,D} {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (196.036,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0461,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64681,0.0415318,-4.32253e-05,2.19828e-08,-4.15992e-12,23671.2,22.757], Tmin=(100,'K'), Tmax=(1511.92,'K')), NASAPolynomial(coeffs=[12.5856,0.00541113,-2.65617e-07,-1.00979e-10,1.11027e-14,21184.2,-31.8244], Tmin=(1511.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(196.036,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + ring(dioxirane) + radical(Cs_P) + radical(CCCJ=O)"""),
)

species(
    label = '[O][C]1CC(=O)O1(3396)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {5,S} {6,S}
2 O u1 p2 c0 {5,S}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
5 C u1 p0 c0 {1,S} {2,S} {4,S}
6 C u0 p0 c0 {1,S} {3,D} {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-65.7444,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0461,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.1553,0.011433,3.10553e-05,-4.50419e-08,1.74519e-11,-7870.31,18.5519], Tmin=(100,'K'), Tmax=(944.245,'K')), NASAPolynomial(coeffs=[6.20988,0.0164763,-5.52369e-06,9.53486e-10,-6.59918e-14,-9248.85,-0.254458], Tmin=(944.245,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-65.7444,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsOs) + ring(Beta-Propiolactone) + radical(CCOJ) + radical(Cs_P)"""),
)

species(
    label = '[O]C1([O])CC1=O(3401)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {5,S}
2 O u1 p2 c0 {5,S}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
5 C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6 C u0 p0 c0 {3,D} {4,S} {5,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (106.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0461,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81154,0.0304372,2.18034e-05,-7.24256e-08,3.68788e-11,12894.5,19.0218], Tmin=(100,'K'), Tmax=(893.126,'K')), NASAPolynomial(coeffs=[21.6097,-0.00987417,8.28992e-06,-1.71552e-09,1.16627e-13,7429.41,-85.0523], Tmin=(893.126,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(106.417,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + ring(cyclopropanone) + radical(C=OCOJ) + radical(C=OCOJ)"""),
)

species(
    label = 'O=[C][CH]C(=O)O(3402)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {5,S} {8,S}
2 O u0 p2 c0 {5,D}
3 O u0 p2 c0 {6,D}
4 C u1 p0 c0 {5,S} {6,S} {7,S}
5 C u0 p0 c0 {1,S} {2,D} {4,S}
6 C u1 p0 c0 {3,D} {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {1,S}
"""),
    E0 = (-178.613,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,1855,455,950,368.299,368.329,368.358,368.426,1561.34],'cm^-1')),
        HinderedRotor(inertia=(0.00124251,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00518374,'amu*angstrom^2'), symmetry=1, barrier=(45.0269,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.467354,'amu*angstrom^2'), symmetry=1, barrier=(45.0216,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0461,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.16067,0.0433599,-5.75286e-05,4.23187e-08,-1.26901e-11,-21418.5,21.3558], Tmin=(100,'K'), Tmax=(810.241,'K')), NASAPolynomial(coeffs=[7.5822,0.0165937,-7.97419e-06,1.54361e-09,-1.08388e-13,-22297,-3.65895], Tmin=(810.241,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-178.613,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)(Cds-O2d)HH) + group(Cds-OdCsOs) + group(Cds-OdCsH) + radical(CCJCO) + radical(C=OCCJ=O)"""),
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
    label = '[O][C]=O(188)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u1 p0 c0 {1,S} {2,D}
"""),
    E0 = (31.5354,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0094,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.75939,0.00186756,1.032e-05,-1.52371e-08,5.8053e-12,3804.5,8.40408], Tmin=(100,'K'), Tmax=(1021.27,'K')), NASAPolynomial(coeffs=[6.36179,0.000422714,-4.06588e-07,1.52497e-10,-1.51971e-14,2816.75,-6.43929], Tmin=(1021.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.5354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(OJC=O) + radical((O)CJOH)"""),
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
    label = '[O]C(=O)C=C=O(2913)',
    structure = adjacencyList("""multiplicity 2
1 O u1 p2 c0 {5,S}
2 O u0 p2 c0 {5,D}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {5,S} {6,D} {7,S}
5 C u0 p0 c0 {1,S} {2,D} {4,S}
6 C u0 p0 c0 {3,D} {4,D}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-109.351,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2120,512.5,787.5,180,180,180,2038.16,2039.17,2039.64],'cm^-1')),
        HinderedRotor(inertia=(0.159867,'amu*angstrom^2'), symmetry=1, barrier=(3.67565,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0381,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.58336,0.0417198,-8.78371e-05,9.75405e-08,-3.91465e-11,-13111.1,18.5848], Tmin=(100,'K'), Tmax=(844.628,'K')), NASAPolynomial(coeffs=[0.136194,0.0269466,-1.47828e-05,2.92506e-09,-2.03602e-13,-11757.4,35.5443], Tmin=(844.628,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-109.351,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + missing(O2d-Cdd) + group(Cds-(Cdd-O2d)CsH) + group(Cds-O2d(Cds-Cds)O2s) + missing(Cdd-CdO2d) + radical(CCOJ)"""),
)

species(
    label = '[CH2][C]=O(2901)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,D}
2 C u1 p0 c0 {3,S} {4,S} {5,S}
3 C u1 p0 c0 {1,D} {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
"""),
    E0 = (160.185,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,539.612,539.669],'cm^-1')),
        HinderedRotor(inertia=(0.000578908,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0366,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2980.68,'J/mol'), sigma=(5.03063,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=465.58 K, Pc=53.12 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39564,0.0101365,2.30772e-06,-8.97606e-09,3.68258e-12,19290.3,10.0702], Tmin=(100,'K'), Tmax=(1068.89,'K')), NASAPolynomial(coeffs=[6.35051,0.00638958,-2.69372e-06,5.4222e-10,-4.02484e-14,18240.9,-6.33577], Tmin=(1068.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CJC=O) + radical(CsCJ=O)"""),
)

species(
    label = '[O]C([O])=CC=O(3403)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {6,S}
2 O u1 p2 c0 {6,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {5,S} {6,D} {7,S}
5 C u0 p0 c0 {3,D} {4,S} {8,S}
6 C u0 p0 c0 {1,S} {2,S} {4,D}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-164.343,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0461,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.30698,0.0417881,-4.47187e-05,2.47147e-08,-5.71105e-12,-19708.7,16.9935], Tmin=(100,'K'), Tmax=(1015.3,'K')), NASAPolynomial(coeffs=[8.20541,0.0185496,-1.03857e-05,2.17064e-09,-1.59869e-13,-20906.4,-11.5524], Tmin=(1015.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-164.343,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cd-Cd(CO)H) + group(Cds-CdsCsCs) + group(Cds-O2d(Cds-Cds)H) + radical(C=COJ) + radical(C=COJ)"""),
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
    E0 = (-119.014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (102.042,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (464.68,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (324.028,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-110.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-55.6141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (229.832,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (123.156,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (140.212,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (89.7948,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-59.3921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (27.9941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (136.249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-97.8642,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (225.516,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-1.33345,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(=O)C[C]=O(3274)'],
    products = ['CO2(13)', 'CH2CO(75)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(=O)C([O])=O(3272)'],
    products = ['[O]C(=O)C[C]=O(3274)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.08731e+10,'s^-1'), n=0.796, Ea=(140.034,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCCJ;CsJ;CO] + [cCCJ;CsJ-HH;C] for rate rule [cCCJ;CsJ-HH;CO]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[C]=O(222)', '[CH2]C([O])=O(882)'],
    products = ['[O]C(=O)C[C]=O(3274)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/CO;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(6)', 'O=[C]C[C]=O(3399)'],
    products = ['[O]C(=O)C[C]=O(3274)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3335.46,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [CO_rad/NonDe;O_birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]C(=O)C[C]=O(3274)'],
    products = ['O=C1CC(=O)O1(3278)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]C(=O)C[C]=O(3274)'],
    products = ['O=C=CC(=O)O(3271)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.9748e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C(=O)C[C]=O(3274)'],
    products = ['O=[C]C[C]1OO1(3400)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.55936e+11,'s^-1'), n=0.551275, Ea=(348.846,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra] for rate rule [R3_CO;carbonyl_intra_Nd;radadd_intra_O]
Euclidian distance = 2.449489742783178
family: Intra_R_Add_Endocyclic
Ea raised from 347.7 to 348.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C(=O)C[C]=O(3274)'],
    products = ['[O][C]1CC(=O)O1(3396)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.006e+11,'s^-1'), n=0.221, Ea=(242.17,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_CO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C(=O)C[C]=O(3274)'],
    products = ['[O]C1([O])CC1=O(3401)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.1935e+10,'s^-1'), n=0.672833, Ea=(259.226,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra] for rate rule [R4_S_CO;carbonylbond_intra;radadd_intra_CO]
Euclidian distance = 1.7320508075688772
family: Intra_R_Add_Exocyclic
Ea raised from 257.8 to 259.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction10',
    reactants = ['O=[C][CH]C(=O)O(3402)'],
    products = ['[O]C(=O)C[C]=O(3274)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(234.612,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_Ext-1R!H-R',), comment="""Estimated from node Root_3R!H->C_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CO(12)', '[CH2]C([O])=O(882)'],
    products = ['[O]C(=O)C[C]=O(3274)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(65100,'m^3/(mol*s)'), n=0.45, Ea=(33.7546,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->O_N-3BrCClFHINPSSi-inRing_3BrCClFHINPSSi->C_Ext-3C-R_Sp-4R!H-3C_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R->O_N-3BrCClFHINPSSi-inRing_3BrCClFHINPSSi->C_Ext-3C-R_Sp-4R!H-3C_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O][C]=O(188)', 'CH2CO(75)'],
    products = ['[O]C(=O)C[C]=O(3274)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(153.031,'m^3/(mol*s)'), n=1.16366, Ea=(23.4817,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.022037706214473284, var=2.3701416358838845, Tref=1000.0, N=230, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(5)', '[O]C(=O)C=C=O(2913)'],
    products = ['[O]C(=O)C[C]=O(3274)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1306.26,'m^3/(mol*s)'), n=1.34286, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.2935929756704599, var=0.2435564102543282, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_Sp-6R!H=4R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_N-Sp-5R!H-2CS_Sp-4R!H-1COS_Ext-4R!H-R_Sp-6R!H=4R!H"""),
)

reaction(
    label = 'reaction14',
    reactants = ['CO2(13)', '[CH2][C]=O(2901)'],
    products = ['[O]C(=O)C[C]=O(3274)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(111.293,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O][C]=O(188)', '[CH2][C]=O(2901)'],
    products = ['[O]C(=O)C[C]=O(3274)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]C(=O)C[C]=O(3274)'],
    products = ['[O]C([O])=CC=O(3403)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.27137e+08,'s^-1'), n=1.53496, Ea=(117.681,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/OneDe] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/OneDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #1230',
    isomers = [
        '[O]C(=O)C[C]=O(3274)',
    ],
    reactants = [
        ('CO2(13)', 'CH2CO(75)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1230',
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

