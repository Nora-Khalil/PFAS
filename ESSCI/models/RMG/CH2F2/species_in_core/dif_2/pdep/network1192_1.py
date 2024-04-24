species(
    label = '[O]C(=O)C([O])C(=O)F(3131)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 O u1 p2 c0 {6,S}
3 O u1 p2 c0 {8,S}
4 O u0 p2 c0 {7,D}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {2,S} {7,S} {8,S} {9,S}
7 C u0 p0 c0 {1,S} {4,D} {6,S}
8 C u0 p0 c0 {3,S} {5,D} {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-505.253,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,486,617,768,1157,1926,180,972.249,972.269,972.273,972.281,972.412,972.478,972.496],'cm^-1')),
        HinderedRotor(inertia=(0.187624,'amu*angstrom^2'), symmetry=1, barrier=(4.31386,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.187539,'amu*angstrom^2'), symmetry=1, barrier=(4.31189,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9211,0.0522175,-8.11518e-05,7.88291e-08,-3.0912e-11,-60699.4,26.1769], Tmin=(100,'K'), Tmax=(778.94,'K')), NASAPolynomial(coeffs=[3.76547,0.0306125,-1.61812e-05,3.22493e-09,-2.28492e-13,-60618.6,20.1025], Tmin=(778.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-505.253,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-CsCsOsH) + group(Cds-OdCsOs) + group(COCsFO) + radical(C=OCOJ) + radical(CCOJ)"""),
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
    label = 'O=CC(=O)F(557)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,D}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {2,D} {5,S} {6,S}
5 C u0 p0 c0 {1,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-483.116,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,286,619,818,1246,1924],'cm^-1')),
        HinderedRotor(inertia=(0.753127,'amu*angstrom^2'), symmetry=1, barrier=(17.3159,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (76.0265,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3426.83,'J/mol'), sigma=(5.15159,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=535.26 K, Pc=56.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88796,0.00794544,8.17071e-05,-2.34537e-07,1.90378e-10,-58102.5,9.57665], Tmin=(10,'K'), Tmax=(438.71,'K')), NASAPolynomial(coeffs=[4.97623,0.0166174,-1.15194e-05,3.74172e-09,-4.59031e-13,-58377,3.18363], Tmin=(438.71,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-483.116,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""ODCC(DO)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]C(=O)C([O])F(574)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u1 p2 c0 {5,S}
3 O u1 p2 c0 {6,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6 C u0 p0 c0 {3,S} {4,D} {5,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-365.212,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([241,322,589,815,1096,1220,1302,2892,180,180,1399.52,1402.3,1406.64,1407.54],'cm^-1')),
        HinderedRotor(inertia=(0.207982,'amu*angstrom^2'), symmetry=1, barrier=(4.7819,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0259,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4306.78,'J/mol'), sigma=(6.22203,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=672.71 K, Pc=40.57 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.43792,0.0452694,-9.27235e-05,1.00355e-07,-3.93593e-11,-43879.1,17.7902], Tmin=(100,'K'), Tmax=(857.523,'K')), NASAPolynomial(coeffs=[0.391153,0.0277345,-1.4678e-05,2.85099e-09,-1.95991e-13,-42532.3,33.156], Tmin=(857.523,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-365.212,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsOs) + radical(O2sj(Cs-F1sCOH)) + radical(CCOJ)"""),
)

species(
    label = '[O]C(=O)O[C](F)C=O(3134)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {6,S} {8,S}
3 O u0 p2 c0 {7,D}
4 O u1 p2 c0 {8,S}
5 O u0 p2 c0 {8,D}
6 C u1 p0 c0 {1,S} {2,S} {7,S}
7 C u0 p0 c0 {3,D} {6,S} {9,S}
8 C u0 p0 c0 {2,S} {4,S} {5,D}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-557.265,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([280,501,1494,1531,2782.5,750,1395,475,1775,1000,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4260.67,'J/mol'), sigma=(6.15227,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=665.51 K, Pc=41.52 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.12206,0.0442759,-3.84894e-05,1.69369e-08,-3.11929e-12,-66958.6,22.2179], Tmin=(100,'K'), Tmax=(1243.51,'K')), NASAPolynomial(coeffs=[9.20066,0.0215063,-1.10235e-05,2.21205e-09,-1.58981e-13,-68719,-13.475], Tmin=(1243.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-557.265,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-(Cds-O2d)H) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(Cds-OdOsOs) + radical(OC=OOJ) + radical(CsCOF1sO2s)"""),
)

species(
    label = '[O]C=C([O])OC(=O)F(3177)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {8,S}
2 O u0 p2 c0 {6,S} {8,S}
3 O u1 p2 c0 {6,S}
4 O u1 p2 c0 {7,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {2,S} {3,S} {7,D}
7 C u0 p0 c0 {4,S} {6,D} {9,S}
8 C u0 p0 c0 {1,S} {2,S} {5,D}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-599.806,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.547266,0.0695696,-8.56496e-05,4.8098e-08,-1.02194e-11,-72010.2,24.6753], Tmin=(100,'K'), Tmax=(1170.55,'K')), NASAPolynomial(coeffs=[19.678,0.00419595,-1.87657e-06,3.86505e-10,-2.94532e-14,-76488.9,-70.6323], Tmin=(1170.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-599.806,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(COFOO) + radical(C=COJ) + radical(C=COJ)"""),
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
    label = '[O]C([O])=CC(=O)F(2666)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u1 p2 c0 {7,S}
3 O u1 p2 c0 {7,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {6,S} {7,D} {8,S}
6 C u0 p0 c0 {1,S} {4,D} {5,S}
7 C u0 p0 c0 {2,S} {3,S} {5,D}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-403.058,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,255,533,799,832,1228,350,440,435,1725,270.033,270.035,270.042],'cm^-1')),
        HinderedRotor(inertia=(0.00231182,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79706,0.0514572,-6.67732e-05,4.40719e-08,-1.15975e-11,-48399.9,20.1972], Tmin=(100,'K'), Tmax=(925.29,'K')), NASAPolynomial(coeffs=[10.215,0.0150662,-7.7787e-06,1.56624e-09,-1.12949e-13,-49957.7,-19.7609], Tmin=(925.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-403.058,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cd-Cd(CO)H) + group(Cds-CdsCsCs) + group(COCFO) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[O]C([C]=O)C(=O)F(3224)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u1 p2 c0 {5,S}
3 O u0 p2 c0 {6,D}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
6 C u0 p0 c0 {1,S} {3,D} {5,S}
7 C u1 p0 c0 {4,D} {5,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-305.245,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,486,617,768,1157,1926,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.112014,'amu*angstrom^2'), symmetry=1, barrier=(2.57543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113488,'amu*angstrom^2'), symmetry=1, barrier=(2.6093,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63798,0.0582828,-0.000102082,9.1236e-08,-3.11703e-11,-36633.5,24.3232], Tmin=(100,'K'), Tmax=(859.865,'K')), NASAPolynomial(coeffs=[7.64252,0.0175891,-8.83191e-06,1.67854e-09,-1.1396e-13,-37194.4,-0.995509], Tmin=(859.865,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-305.245,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(COCsFO) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(C=OCCJ=O)"""),
)

species(
    label = 'O=C(F)C1OOC1=O(3138)',
    structure = adjacencyList("""1 F u0 p3 c0 {8,S}
2 O u0 p2 c0 {3,S} {6,S}
3 O u0 p2 c0 {2,S} {7,S}
4 O u0 p2 c0 {7,D}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {2,S} {7,S} {8,S} {9,S}
7 C u0 p0 c0 {3,S} {4,D} {6,S}
8 C u0 p0 c0 {1,S} {5,D} {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-636.359,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22982,0.0163599,7.13581e-05,-1.10621e-07,4.42201e-11,-76452.2,22.8179], Tmin=(100,'K'), Tmax=(969.127,'K')), NASAPolynomial(coeffs=[18.1006,0.00754001,-2.72739e-06,6.96834e-10,-6.53583e-14,-82190.3,-66.9854], Tmin=(969.127,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-636.359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(Cs-CsCsOsH) + group(Cds-OdCsOs) + group(COCsFO) + ring(Cyclobutane)"""),
)

species(
    label = 'O=C(O)C(=O)C(=O)F(3225)',
    structure = adjacencyList("""1 F u0 p3 c0 {8,S}
2 O u0 p2 c0 {7,S} {9,S}
3 O u0 p2 c0 {6,D}
4 O u0 p2 c0 {7,D}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {3,D} {7,S} {8,S}
7 C u0 p0 c0 {2,S} {4,D} {6,S}
8 C u0 p0 c0 {1,S} {5,D} {6,S}
9 H u0 p0 c0 {2,S}
"""),
    E0 = (-839.802,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11357,0.0696813,-0.000117111,1.02123e-07,-3.49881e-11,-100907,23.9927], Tmin=(100,'K'), Tmax=(800.123,'K')), NASAPolynomial(coeffs=[9.47272,0.0204958,-1.1037e-05,2.18801e-09,-1.53526e-13,-102008,-12.9919], Tmin=(800.123,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-839.802,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-O2d(Cds-O2d)(Cds-O2d)) + group(Cds-O2d(Cds-O2d)O2s) + group(COCFO)"""),
)

species(
    label = '[O]C([C]1OO1)C(=O)F(3226)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {8,S}
2 O u0 p2 c0 {3,S} {7,S}
3 O u0 p2 c0 {2,S} {7,S}
4 O u1 p2 c0 {6,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
7 C u1 p0 c0 {2,S} {3,S} {6,S}
8 C u0 p0 c0 {1,S} {5,D} {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-144.066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47865,0.0448387,-2.17942e-05,-1.67309e-08,1.33657e-11,-17226.5,29.0598], Tmin=(100,'K'), Tmax=(921.671,'K')), NASAPolynomial(coeffs=[17.4154,0.0034232,4.47918e-07,-1.53288e-10,8.59827e-15,-21342.7,-52.9196], Tmin=(921.671,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-144.066,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsOsH) + group(COCsFO) + ring(dioxirane) + radical(C=OCOJ) + radical(Cs_P)"""),
)

species(
    label = '[O]C(=O)C1OO[C]1F(3227)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {3,S} {6,S}
3 O u0 p2 c0 {2,S} {7,S}
4 O u1 p2 c0 {8,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {2,S} {7,S} {8,S} {9,S}
7 C u1 p0 c0 {1,S} {3,S} {6,S}
8 C u0 p0 c0 {4,S} {5,D} {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-200.566,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.57016,0.0355012,-2.1541e-05,5.4915e-09,-5.33736e-13,-24075.7,23.8755], Tmin=(100,'K'), Tmax=(2345.58,'K')), NASAPolynomial(coeffs=[16.1012,0.012426,-6.78426e-06,1.29727e-09,-8.66967e-14,-30423.2,-52.9394], Tmin=(2345.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-200.566,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsOsH) + group(CsCFHO) + group(Cds-OdCsOs) + ring(Cs-Cs(F)-O2s-O2s) + radical(CCOJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O][C](F)C1OOC1=O(3162)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {8,S}
2 O u0 p2 c0 {3,S} {6,S}
3 O u0 p2 c0 {2,S} {7,S}
4 O u0 p2 c0 {7,D}
5 O u1 p2 c0 {8,S}
6 C u0 p0 c0 {2,S} {7,S} {8,S} {9,S}
7 C u0 p0 c0 {3,S} {4,D} {6,S}
8 C u1 p0 c0 {1,S} {5,S} {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-264.967,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.86863,0.0267462,4.30571e-05,-8.38939e-08,3.5577e-11,-31773.3,26.4346], Tmin=(100,'K'), Tmax=(967.167,'K')), NASAPolynomial(coeffs=[19.037,0.00572363,-1.85678e-06,4.9839e-10,-4.90962e-14,-37431.9,-67.9049], Tmin=(967.167,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-264.967,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCFHO) + group(Cds-OdCsOs) + ring(Cyclobutane) + radical(O2sj(Cs-CsF1sH)) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]C1[C](F)OOC1=O(3144)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {3,S} {7,S}
3 O u0 p2 c0 {2,S} {8,S}
4 O u1 p2 c0 {6,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
7 C u1 p0 c0 {1,S} {2,S} {6,S}
8 C u0 p0 c0 {3,S} {5,D} {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-329.575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.18631,0.0168182,7.14683e-05,-1.14084e-07,4.67007e-11,-39552.4,25.3637], Tmin=(100,'K'), Tmax=(953.94,'K')), NASAPolynomial(coeffs=[19.0616,0.00431382,-4.72971e-07,2.09883e-10,-2.9658e-14,-45422.6,-69.1475], Tmin=(953.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-329.575,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCFHO) + group(Cds-OdCsOs) + ring(Cyclopentane) + radical(C=OCOJ) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]C(=O)C1OC1([O])F(3228)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {6,S} {7,S}
3 O u1 p2 c0 {7,S}
4 O u1 p2 c0 {8,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {2,S} {7,S} {8,S} {9,S}
7 C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
8 C u0 p0 c0 {4,S} {5,D} {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-381.426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13111,0.0360207,-1.83345e-05,1.58324e-09,6.96572e-13,-45803.6,23.8019], Tmin=(100,'K'), Tmax=(1433.84,'K')), NASAPolynomial(coeffs=[11.9885,0.0178149,-9.01121e-06,1.76885e-09,-1.23975e-13,-49585.7,-30.6379], Tmin=(1433.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-381.426,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsOsH) + group(CsCFOO) + group(Cds-OdCsOs) + ring(Cs(F)(O2)-O2s-Cs) + radical(O2sj(Cs-F1sO2sCs)) + radical(CCOJ)"""),
)

species(
    label = '[O]C1([O])OC1C(=O)F(3206)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {8,S}
2 O u0 p2 c0 {6,S} {7,S}
3 O u1 p2 c0 {7,S}
4 O u1 p2 c0 {7,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {2,S} {7,S} {8,S} {9,S}
7 C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
8 C u0 p0 c0 {1,S} {5,D} {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-366.833,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90837,0.0492429,-6.67174e-05,5.3052e-08,-1.64446e-11,-44047.5,25.7227], Tmin=(100,'K'), Tmax=(954.061,'K')), NASAPolynomial(coeffs=[6.18639,0.0215786,-7.92798e-06,1.28419e-09,-7.89335e-14,-44421,7.60526], Tmin=(954.061,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-366.833,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsOsOs) + group(COCsFO) + ring(Ethylene_oxide) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[O]C1C(=O)OC1([O])F(3171)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {7,S} {8,S}
3 O u1 p2 c0 {6,S}
4 O u1 p2 c0 {7,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
7 C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
8 C u0 p0 c0 {2,S} {5,D} {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-398.103,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.14853,0.0264898,2.81361e-05,-6.4637e-08,2.95935e-11,-47800.8,24.0903], Tmin=(100,'K'), Tmax=(917.259,'K')), NASAPolynomial(coeffs=[15.1681,0.00674413,-1.29442e-07,-8.13914e-11,3.16951e-15,-51747.1,-46.0887], Tmin=(917.259,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-398.103,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCFOO) + group(Cds-OdCsOs) + ring(Beta-Propiolactone) + radical(C=OCOJ) + radical(O2sj(Cs-F1sO2sCs))"""),
)

species(
    label = '[O]C(=O)C(=O)[C](O)F(3229)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {7,S} {9,S}
3 O u0 p2 c0 {6,D}
4 O u1 p2 c0 {8,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {3,D} {7,S} {8,S}
7 C u1 p0 c0 {1,S} {2,S} {6,S}
8 C u0 p0 c0 {4,S} {5,D} {6,S}
9 H u0 p0 c0 {2,S}
"""),
    E0 = (-518.054,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,375,552.5,462.5,1710,280,501,1494,1531,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09107,0.0720908,-0.000126792,1.1678e-07,-4.18738e-11,-62210.6,25.1635], Tmin=(100,'K'), Tmax=(809.667,'K')), NASAPolynomial(coeffs=[7.87979,0.0248137,-1.37537e-05,2.74896e-09,-1.93429e-13,-62859.6,-3.37391], Tmin=(809.667,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-518.054,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)O2s) + radical(C=OC=OOJ) + radical(CsCOF1sO2s)"""),
)

species(
    label = '[O][C](C(=O)O)C(=O)F(3230)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {8,S}
2 O u0 p2 c0 {7,S} {9,S}
3 O u1 p2 c0 {6,S}
4 O u0 p2 c0 {7,D}
5 O u0 p2 c0 {8,D}
6 C u1 p0 c0 {3,S} {7,S} {8,S}
7 C u0 p0 c0 {2,S} {4,D} {6,S}
8 C u0 p0 c0 {1,S} {5,D} {6,S}
9 H u0 p0 c0 {2,S}
"""),
    E0 = (-554.331,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,360,370,350,611,648,830,1210,1753,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40147,0.0636932,-0.000105596,9.49066e-08,-3.35255e-11,-66583.3,27.3012], Tmin=(100,'K'), Tmax=(814.654,'K')), NASAPolynomial(coeffs=[7.32508,0.023845,-1.24078e-05,2.42883e-09,-1.69306e-13,-67191.3,2.12947], Tmin=(814.654,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-554.331,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-CsCsOsH) + group(Cds-OdCsOs) + group(COCsFO) + radical(C=OCOJ) + radical(C2CsJOH)"""),
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
    label = '[O]C(=O)C=O(656)',
    structure = adjacencyList("""multiplicity 2
1 O u1 p2 c0 {4,S}
2 O u0 p2 c0 {4,D}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {2,D} {5,S}
5 C u0 p0 c0 {3,D} {4,S} {6,S}
6 H u0 p0 c0 {5,S}
"""),
    E0 = (-220.22,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,180,554.278,676.839,684.263,3881.32],'cm^-1')),
        HinderedRotor(inertia=(1.37032,'amu*angstrom^2'), symmetry=1, barrier=(31.5064,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (73.0274,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.94681,0.0265976,-2.97423e-05,1.90555e-08,-5.3726e-12,-26451.4,14.3847], Tmin=(100,'K'), Tmax=(825.894,'K')), NASAPolynomial(coeffs=[5.38716,0.0147781,-8.27521e-06,1.72677e-09,-1.27055e-13,-26854.5,3.07825], Tmin=(825.894,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-220.22,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-O2d(Cds-O2d)O2s) + group(Cds-O2d(Cds-O2d)H) + radical(C=OC=OOJ)"""),
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
    label = '[O]C(=O)C(=O)C(=O)F(3231)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {6,D}
3 O u1 p2 c0 {8,S}
4 O u0 p2 c0 {7,D}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {2,D} {7,S} {8,S}
7 C u0 p0 c0 {1,S} {4,D} {6,S}
8 C u0 p0 c0 {3,S} {5,D} {6,S}
"""),
    E0 = (-578.213,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,286,619,818,1246,1924,245.988,245.992,245.993,245.993,246.003,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.553334,'amu*angstrom^2'), symmetry=1, barrier=(23.7638,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00278614,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (119.028,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07757,0.0729091,-0.000138151,1.26283e-07,-4.37284e-11,-69446,25.2055], Tmin=(100,'K'), Tmax=(846.979,'K')), NASAPolynomial(coeffs=[9.31222,0.0172552,-9.89863e-06,1.96452e-09,-1.36047e-13,-70239.6,-9.60458], Tmin=(846.979,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-578.213,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-O2d(Cds-O2d)(Cds-O2d)) + group(Cds-O2d(Cds-O2d)O2s) + group(COCFO) + radical(C=OC=OOJ)"""),
)

species(
    label = '[O][CH]C(=O)F(438)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u1 p2 c0 {4,S}
3 O u0 p2 c0 {5,D}
4 C u1 p0 c0 {2,S} {5,S} {6,S}
5 C u0 p0 c0 {1,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-214.477,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,611,648,830,1210,1753,380.326,381.637],'cm^-1')),
        HinderedRotor(inertia=(0.483418,'amu*angstrom^2'), symmetry=1, barrier=(49.7921,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (76.0265,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.80089,0.0290438,-3.0283e-05,1.66613e-08,-3.81625e-12,-25754.7,13.8243], Tmin=(100,'K'), Tmax=(1028.25,'K')), NASAPolynomial(coeffs=[6.95398,0.0128879,-6.71484e-06,1.38085e-09,-1.01081e-13,-26608.7,-6.32768], Tmin=(1028.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-214.477,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(C=OCOJ) + radical(OCJC=O)"""),
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
    label = '[O]C(=O)C(=O)[C]=O(3220)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {5,D}
2 O u1 p2 c0 {6,S}
3 O u0 p2 c0 {6,D}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {1,D} {6,S} {7,S}
6 C u0 p0 c0 {2,S} {3,D} {5,S}
7 C u1 p0 c0 {4,D} {5,S}
"""),
    E0 = (-201.565,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,1855,455,950,180,180,180,180,180,971.52],'cm^-1')),
        HinderedRotor(inertia=(0.00303399,'amu*angstrom^2'), symmetry=1, barrier=(2.03453,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0680634,'amu*angstrom^2'), symmetry=1, barrier=(45.5076,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69401,0.056308,-9.79267e-05,8.78166e-08,-3.08432e-11,-24165,20.5409], Tmin=(100,'K'), Tmax=(796.505,'K')), NASAPolynomial(coeffs=[8.01747,0.0169742,-9.58141e-06,1.92806e-09,-1.36259e-13,-24932,-7.01846], Tmin=(796.505,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-201.565,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-O2d(Cds-O2d)(Cds-O2d)) + group(Cds-O2d(Cds-O2d)O2s) + group(Cds-O2d(Cds-O2d)H) + radical(C=OC=OOJ) + radical(C=OC=OCJ=O)"""),
)

species(
    label = '[O]C([O])=C(O)C(=O)F(3232)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {6,S} {9,S}
3 O u1 p2 c0 {8,S}
4 O u1 p2 c0 {8,S}
5 O u0 p2 c0 {7,D}
6 C u0 p0 c0 {2,S} {7,S} {8,D}
7 C u0 p0 c0 {1,S} {5,D} {6,S}
8 C u0 p0 c0 {3,S} {4,S} {6,D}
9 H u0 p0 c0 {2,S}
"""),
    E0 = (-609.615,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.348221,0.0758484,-0.000103619,6.44962e-08,-1.51076e-11,-73184.4,24.0532], Tmin=(100,'K'), Tmax=(1066.41,'K')), NASAPolynomial(coeffs=[19.8795,0.00258736,-5.69864e-07,7.37518e-11,-4.75846e-15,-77350,-71.4296], Tmin=(1066.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-609.615,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsCsCs) + group(COCFO) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[O]C(=O)C([C]=O)OF(3233)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {6,S}
3 O u1 p2 c0 {7,S}
4 O u0 p2 c0 {7,D}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {2,S} {7,S} {8,S} {9,S}
7 C u0 p0 c0 {3,S} {4,D} {6,S}
8 C u1 p0 c0 {5,D} {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-200.459,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,1855,455,950,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23359,0.0700164,-0.00012458,1.18589e-07,-4.36029e-11,-24018.9,25.7135], Tmin=(100,'K'), Tmax=(820.682,'K')), NASAPolynomial(coeffs=[6.03545,0.0284812,-1.55256e-05,3.08004e-09,-2.15701e-13,-24196.5,7.21621], Tmin=(820.682,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-200.459,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-O2d)H) + group(Cs-CsCsOsH) + group(Cds-OdCsOs) + group(Cds-OdCsH) + radical(CCOJ) + radical(C=OCCJ=O)"""),
)

species(
    label = '[O]C([C]=O)C(=O)OF(3234)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {7,S}
3 O u1 p2 c0 {6,S}
4 O u0 p2 c0 {7,D}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
7 C u0 p0 c0 {2,S} {4,D} {6,S}
8 C u1 p0 c0 {5,D} {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-182.431,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,1380,1390,370,380,2900,435,1855,455,950,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11028,0.0698267,-0.000115849,1.00969e-07,-3.4569e-11,-21843.4,26.7713], Tmin=(100,'K'), Tmax=(811.48,'K')), NASAPolynomial(coeffs=[9.09658,0.0219132,-1.14834e-05,2.24916e-09,-1.56711e-13,-22858.1,-8.35576], Tmin=(811.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-182.431,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(Cs-CsCsOsH) + group(Cds-OdCsOs) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(C=OCCJ=O)"""),
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
    E0 = (-190.882,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (204.092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-54.8482,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-72.0535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (154.347,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (252.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-182.598,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-127.482,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (170.305,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (113.805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (51.2877,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-15.2036,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-67.0551,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-52.4616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-68.8392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (22.7642,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (5.93527,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-40.4889,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-65.9153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-47.8746,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-186.016,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (131.429,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (83.4196,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-65.4221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (217.174,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (205.725,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(=O)C([O])C(=O)F(3131)'],
    products = ['CO2(13)', 'O=CC(=O)F(557)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(12)', '[O]C(=O)C([O])F(574)'],
    products = ['[O]C(=O)C([O])C(=O)F(3131)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.0026956,'m^3/(mol*s)'), n=2.93313, Ea=(373.673,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_2Br1sCl1sF1sH->F1s',), comment="""Estimated from node Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_2Br1sCl1sF1sH->F1s"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O]C(=O)C([O])C(=O)F(3131)'],
    products = ['[O]C(=O)O[C](F)C=O(3134)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(136.034,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]C(=O)C([O])C(=O)F(3131)'],
    products = ['[O]C=C([O])OC(=O)F(3177)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.87873e+12,'s^-1'), n=0.324012, Ea=(118.829,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(6)', '[O]C([O])=CC(=O)F(2666)'],
    products = ['[O]C(=O)C([O])C(=O)F(3131)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/TwoDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['O(6)', '[O]C([C]=O)C(=O)F(3224)'],
    products = ['[O]C(=O)C([O])C(=O)F(3131)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [CO_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C(=O)C([O])C(=O)F(3131)'],
    products = ['O=C(F)C1OOC1=O(3138)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C(=O)C([O])C(=O)F(3131)'],
    products = ['O=C(O)C(=O)C(=O)F(3225)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C(=O)C([O])C(=O)F(3131)'],
    products = ['[O]C([C]1OO1)C(=O)F(3226)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.55936e+11,'s^-1'), n=0.551275, Ea=(361.188,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra] for rate rule [R3_CO;carbonyl_intra_Nd;radadd_intra_O]
Euclidian distance = 2.449489742783178
family: Intra_R_Add_Endocyclic
Ea raised from 360.6 to 361.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C(=O)C([O])C(=O)F(3131)'],
    products = ['[O]C(=O)C1OO[C]1F(3227)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(304.687,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 303.2 to 304.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]C(=O)C([O])C(=O)F(3131)'],
    products = ['[O][C](F)C1OOC1=O(3162)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.006e+11,'s^-1'), n=0.221, Ea=(242.17,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_O]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C(=O)C([O])C(=O)F(3131)'],
    products = ['[O]C1[C](F)OOC1=O(3144)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(6.24664e+09,'s^-1'), n=0.5388, Ea=(175.679,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra] for rate rule [R5_SS_CO;carbonyl_intra;radadd_intra_O]
Euclidian distance = 1.7320508075688772
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 171.7 to 175.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]C(=O)C([O])C(=O)F(3131)'],
    products = ['[O]C(=O)C1OC1([O])F(3228)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.785e+11,'s^-1'), n=0.342, Ea=(123.827,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonylbond_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic
Ea raised from 121.9 to 123.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]C(=O)C([O])C(=O)F(3131)'],
    products = ['[O]C1([O])OC1C(=O)F(3206)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.785e+11,'s^-1'), n=0.342, Ea=(138.421,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonylbond_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic
Ea raised from 138.0 to 138.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O]C(=O)C([O])C(=O)F(3131)'],
    products = ['[O]C1C(=O)OC1([O])F(3171)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(5.448e+10,'s^-1'), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_O] for rate rule [R5_SS_CO;carbonylbond_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]C(=O)C(=O)[C](O)F(3229)'],
    products = ['[O]C(=O)C([O])C(=O)F(3131)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(226.447,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_Ext-1R!H-R',), comment="""Estimated from node Root_3R!H->C_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O][C](C(=O)O)C(=O)F(3230)'],
    products = ['[O]C(=O)C([O])C(=O)F(3131)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(245.895,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_Ext-1R!H-R',), comment="""Estimated from node Root_3R!H->C_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CFO(50)', '[O]C(=O)C=O(656)'],
    products = ['[O]C(=O)C([O])C(=O)F(3131)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=-1.07934e-11, Ea=(55.7191,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O][C]=O(188)', 'O=CC(=O)F(557)'],
    products = ['[O]C(=O)C([O])C(=O)F(3131)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=-1.07934e-11, Ea=(71.2947,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(5)', '[O]C(=O)C(=O)C(=O)F(3231)'],
    products = ['[O]C(=O)C([O])C(=O)F(3131)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(39.7,'m^3/(mol*s)'), n=1.88, Ea=(4.16296,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CO2(13)', '[O][CH]C(=O)F(438)'],
    products = ['[O]C(=O)C([O])C(=O)F(3131)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(117.228,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O][C]=O(188)', '[O][CH]C(=O)F(438)'],
    products = ['[O]C(=O)C([O])C(=O)F(3131)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction23',
    reactants = ['HF(38)', '[O]C(=O)C(=O)[C]=O(3220)'],
    products = ['[O]C(=O)C([O])C(=O)F(3131)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(251.727,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]C(=O)C([O])C(=O)F(3131)'],
    products = ['[O]C([O])=C(O)C(=O)F(3232)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.66008e+10,'s^-1'), n=0.776642, Ea=(125.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;Cs_H_out_noH] + [R2H_S;O_rad_out;Cs_H_out] for rate rule [R2H_S;O_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O]C(=O)C([C]=O)OF(3233)'],
    products = ['[O]C(=O)C([O])C(=O)F(3131)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(103.262,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O]C([C]=O)C(=O)OF(3234)'],
    products = ['[O]C(=O)C([O])C(=O)F(3131)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(73.7853,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

network(
    label = 'PDepNetwork #1192',
    isomers = [
        '[O]C(=O)C([O])C(=O)F(3131)',
    ],
    reactants = [
        ('CO2(13)', 'O=CC(=O)F(557)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1192',
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

