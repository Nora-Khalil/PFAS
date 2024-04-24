species(
    label = '[CH]=CC([O])=O(2797)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {4,S}
2 O u0 p2 c0 {4,D}
3 C u0 p0 c0 {4,S} {5,D} {6,S}
4 C u0 p0 c0 {1,S} {2,D} {3,S}
5 C u1 p0 c0 {3,D} {7,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (185.212,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3120,650,792.5,1650,4000,4000,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.000717289,'amu*angstrom^2'), symmetry=1, barrier=(8.14414,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0467,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36229,0.0216391,-5.09949e-05,6.30506e-08,-2.69907e-11,22291.4,13.0249], Tmin=(100,'K'), Tmax=(834.916,'K')), NASAPolynomial(coeffs=[-0.32402,0.0192644,-1.07333e-05,2.16077e-09,-1.52239e-13,23605.2,34.3261], Tmin=(834.916,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(185.212,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + missing(Cd-COCdH) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_P)"""),
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
    label = 'C2H2(70)',
    structure = adjacencyList("""1 C u0 p0 c0 {2,T} {3,S}
2 C u0 p0 c0 {1,T} {4,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (217.784,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,584.389,584.389,2772.01],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (26.0372,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.80868,0.0233616,-3.55172e-05,2.80153e-08,-8.50075e-12,26429,13.9397], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.65878,0.00488397,-1.60829e-06,2.46975e-10,-1.38606e-14,25759.4,-3.99838], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(217.784,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(87.302,'J/(mol*K)'), label="""C2H2""", comment="""Thermo library: FFCM1(-)"""),
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
    label = '[CH]C=C=O(2817)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {4,D}
2 C u0 p0 c0 {3,S} {4,D} {5,S}
3 C u2 p0 c0 {2,S} {6,S}
4 C u0 p0 c0 {1,D} {2,D}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (292.962,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2120,512.5,787.5,246.972,246.978,246.987],'cm^-1')),
        HinderedRotor(inertia=(1.17177,'amu*angstrom^2'), symmetry=1, barrier=(50.7182,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0473,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.00496,0.0183623,-5.86396e-07,-8.19936e-09,3.38635e-12,35273.9,13.1075], Tmin=(100,'K'), Tmax=(1106.37,'K')), NASAPolynomial(coeffs=[5.96291,0.0150316,-6.05408e-06,1.11097e-09,-7.67693e-14,34168.7,-3.49853], Tmin=(1106.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(292.962,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(AllylJ2_triplet)"""),
)

species(
    label = 'O=C1C=CO1(2809)',
    structure = adjacencyList("""1 O u0 p2 c0 {4,S} {5,S}
2 O u0 p2 c0 {5,D}
3 C u0 p0 c0 {4,D} {5,S} {6,S}
4 C u0 p0 c0 {1,S} {3,D} {7,S}
5 C u0 p0 c0 {1,S} {2,D} {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-101.797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (70.0467,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.54788,0.00321963,2.3881e-05,-3.4328e-08,1.3354e-11,-12221,7.94787], Tmin=(100,'K'), Tmax=(978.736,'K')), NASAPolynomial(coeffs=[8.00577,0.00209315,-5.88325e-07,1.82588e-10,-1.88174e-14,-13912.3,-17.6453], Tmin=(978.736,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-101.797,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + missing(Cd-COCdH) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)O2s) + ring(Oxetene)"""),
)

species(
    label = '[CH]=C[C]1OO1(2818)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {3,S}
2 O u0 p2 c0 {1,S} {3,S}
3 C u1 p0 c0 {1,S} {2,S} {4,S}
4 C u0 p0 c0 {3,S} {5,D} {6,S}
5 C u1 p0 c0 {4,D} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (492.619,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (70.0467,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.59041,0.021209,1.26231e-05,-4.13922e-08,2.0639e-11,59308.3,16.7293], Tmin=(100,'K'), Tmax=(896.74,'K')), NASAPolynomial(coeffs=[13.07,0.00100843,2.01114e-06,-4.92979e-10,3.40151e-14,56361.6,-38.6374], Tmin=(896.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(492.619,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(dioxirane) + radical(Cs_P) + radical(Cds_P)"""),
)

species(
    label = 'O=C1[CH][CH]O1(2819)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {4,S} {5,S}
2 O u0 p2 c0 {5,D}
3 C u1 p0 c0 {4,S} {5,S} {6,S}
4 C u1 p0 c0 {1,S} {3,S} {7,S}
5 C u0 p0 c0 {1,S} {2,D} {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (97.3646,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (70.0467,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.24962,0.00275567,6.26876e-05,-9.17663e-08,3.87219e-11,11750.3,16.1714], Tmin=(100,'K'), Tmax=(887.237,'K')), NASAPolynomial(coeffs=[11.5876,0.00156839,3.1498e-06,-7.85164e-10,5.5276e-14,8837.93,-31.1314], Tmin=(887.237,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(97.3646,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + ring(Beta-Propiolactone) + radical(CCJCO) + radical(CCsJOC(O))"""),
)

species(
    label = '[O]C1([O])C=C1(2820)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {3,S}
2 O u1 p2 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
4 C u0 p0 c0 {3,S} {5,D} {6,S}
5 C u0 p0 c0 {3,S} {4,D} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (350.455,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (70.0467,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.86399,0.0190419,-8.3405e-06,9.99786e-10,1.87671e-14,42138.8,9.14415], Tmin=(100,'K'), Tmax=(2568.47,'K')), NASAPolynomial(coeffs=[20.7137,-6.02635e-05,-1.35376e-06,2.68425e-10,-1.53722e-14,31128.6,-92.6241], Tmin=(2568.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(350.455,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopropene) + radical(CC(C)(O)OJ) + radical(CC(C)(O)OJ)"""),
)

species(
    label = '[CH]=C=C([O])O(2821)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,S} {6,S}
2 O u1 p2 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {4,D}
4 C u0 p0 c0 {3,D} {5,D}
5 C u1 p0 c0 {4,D} {7,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (112.738,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,540,610,2055,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(1.40902,'amu*angstrom^2'), symmetry=1, barrier=(32.396,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0467,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54222,0.0471399,-6.12035e-05,3.64213e-08,-7.92208e-12,13653.8,16.7163], Tmin=(100,'K'), Tmax=(1313.34,'K')), NASAPolynomial(coeffs=[14.1536,0.000294517,1.9339e-06,-5.18304e-10,4.04407e-14,11068.6,-44.7943], Tmin=(1313.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(112.738,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=C=CJ)"""),
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
    label = 'C#CC([O])=O(1114)',
    structure = adjacencyList("""multiplicity 2
1 O u1 p2 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 C u0 p0 c0 {3,S} {5,T}
5 C u0 p0 c0 {4,T} {6,S}
6 H u0 p0 c0 {5,S}
"""),
    E0 = (87.8995,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,574.466,574.702,575.363,575.372,575.44],'cm^-1')),
        HinderedRotor(inertia=(0.000509499,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0387,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.42264,0.0164477,1.21129e-05,-6.8773e-08,6.81838e-11,10588.7,10.6553], Tmin=(100,'K'), Tmax=(386.387,'K')), NASAPolynomial(coeffs=[3.92522,0.0164858,-8.38132e-06,1.69263e-09,-1.22498e-13,10510.7,8.20225], Tmin=(386.387,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(87.8995,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo library: 2-BTP_G4 + radical(CCOJ)"""),
)

species(
    label = '[CH]=[CH](2414)',
    structure = adjacencyList("""multiplicity 3
1 C u1 p0 c0 {2,D} {3,S}
2 C u1 p0 c0 {1,D} {4,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (590.758,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1227.31,2941.72],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (26.0372,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88059,-0.000843518,2.04097e-05,-2.51948e-08,9.47595e-12,71059.3,4.72113], Tmin=(100,'K'), Tmax=(935.363,'K')), NASAPolynomial(coeffs=[4.92139,0.00358275,-9.24454e-07,1.57276e-10,-1.19543e-14,70476.2,-2.30644], Tmin=(935.363,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(590.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""C2H2(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=C=C([O])[O](2822)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {4,S}
2 O u1 p2 c0 {4,S}
3 C u0 p0 c0 {5,D} {6,S} {7,S}
4 C u0 p0 c0 {1,S} {2,S} {5,D}
5 C u0 p0 c0 {3,D} {4,D}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (99.7239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (70.0467,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.23483,0.0355648,-3.65588e-05,1.84632e-08,-3.63376e-12,12060.3,15.0232], Tmin=(100,'K'), Tmax=(1243.07,'K')), NASAPolynomial(coeffs=[10.8099,0.00797113,-3.26137e-06,6.05288e-10,-4.22209e-14,9928.46,-28.2123], Tmin=(1243.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(99.7239,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=COJ)"""),
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
    E0 = (16.7441,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (367.528,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (25.0284,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (325.589,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (154.799,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (184.266,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (201.209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (122.087,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (149.053,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (113.383,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (453.826,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (184.702,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=CC([O])=O(2797)'],
    products = ['CO2(13)', 'C2H2(70)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O(6)', '[CH]C=C=O(2817)'],
    products = ['[CH]=CC([O])=O(2797)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [CO_rad/OneDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=CC([O])=O(2797)'],
    products = ['O=C1C=CO1(2809)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4;CdsingleH_rad_out;Opri_rad]
Euclidian distance = 2.23606797749979
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]=CC([O])=O(2797)'],
    products = ['[CH]=C[C]1OO1(2818)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.55936e+11,'s^-1'), n=0.551275, Ea=(308.845,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra] for rate rule [R3_CO;carbonyl_intra_De;radadd_intra_O]
Euclidian distance = 2.449489742783178
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=CC([O])=O(2797)'],
    products = ['O=C1[CH][CH]O1(2819)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.906e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R4_D_CO;carbonyl_intra;radadd_intra_cdsingleH]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=CC([O])=O(2797)'],
    products = ['[O]C1([O])C=C1(2820)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(9.36651e+10,'s^-1'), n=0.5685, Ea=(167.522,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R4_D_CO;carbonylbond_intra;radadd_intra_cdsingleH]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=C=C([O])O(2821)'],
    products = ['[CH]=CC([O])=O(2797)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(256.939,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_Ext-1R!H-R',), comment="""Estimated from node Root_3R!H->C_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O][C]=O(188)', 'C2H2(70)'],
    products = ['[CH]=CC([O])=O(2797)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(8.38778,'m^3/(mol*s)'), n=1.76036, Ea=(41.2352,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.10151511122778119, var=2.139043046412222, Tref=1000.0, N=264, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(5)', 'C#CC([O])=O(1114)'],
    products = ['[CH]=CC([O])=O(2797)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1348.96,'m^3/(mol*s)'), n=1.432, Ea=(17.8168,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_Sp-5R!H-4CClOS_Ext-4CClOS-R_N-Sp-6R!H-4CClOS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_Sp-5R!H-4CClOS_Ext-4CClOS-R_N-Sp-6R!H-4CClOS"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CO2(13)', '[CH]=[CH](2414)'],
    products = ['[CH]=CC([O])=O(2797)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.04e+06,'m^3/(mol*s)'), n=-1.07934e-11, Ea=(94.2299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O][C]=O(188)', '[CH]=[CH](2414)'],
    products = ['[CH]=CC([O])=O(2797)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.47663e+07,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=CC([O])=O(2797)'],
    products = ['C=C=C([O])[O](2822)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(481900,'s^-1'), n=2.375, Ea=(167.958,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 121 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleDe
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #1063',
    isomers = [
        '[CH]=CC([O])=O(2797)',
    ],
    reactants = [
        ('CO2(13)', 'C2H2(70)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1063',
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

