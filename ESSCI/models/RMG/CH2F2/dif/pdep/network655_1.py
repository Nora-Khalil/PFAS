species(
    label = '[CH2]OC[O](1582)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,S} {4,S}
2 O u1 p2 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
4 C u1 p0 c0 {1,S} {7,S} {8,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (25.0348,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,421.349,421.522,421.584],'cm^-1')),
        HinderedRotor(inertia=(0.000948531,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000947498,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.74078,0.0211544,3.72298e-06,-1.92868e-08,8.49804e-12,3062.05,15.0221], Tmin=(100,'K'), Tmax=(1016.31,'K')), NASAPolynomial(coeffs=[8.92519,0.0118917,-4.86002e-06,9.41349e-10,-6.87754e-14,1026.31,-18.7452], Tmin=(1016.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(25.0348,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-OsOsHH) + group(Cs-OsHHH) + radical(OCOJ) + radical(CsJOCH3)"""),
)

species(
    label = 'CH2O(19)',
    structure = adjacencyList("""1 O u0 p2 c0 {2,D}
2 C u0 p0 c0 {1,D} {3,S} {4,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (-119.055,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.79372,-0.00990833,3.7322e-05,-3.79285e-08,1.31773e-11,-14379.2,0.602798], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.16953,0.00619321,-2.25056e-06,3.65976e-10,-2.20149e-14,-14548.7,6.04208], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-119.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH2O""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'CH2(T)(17)',
    structure = adjacencyList("""multiplicity 3
1 C u2 p0 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (381.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1066.91,2790.99,3622.37],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.01192,-0.000154979,3.26298e-06,-2.40422e-09,5.69497e-13,45867.7,0.5332], Tmin=(100,'K'), Tmax=(1104.57,'K')), NASAPolynomial(coeffs=[3.14983,0.00296674,-9.76055e-07,1.54115e-10,-9.50337e-15,46058.1,4.77807], Tmin=(1104.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[O]C[O](1545)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {3,S}
2 O u1 p2 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (53.4121,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1457.88,1457.91],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (46.0254,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[5.67714,-0.00771635,2.7117e-05,-1.79271e-08,3.39818e-12,6337.74,0.967298], Tmin=(100,'K'), Tmax=(1894.84,'K')), NASAPolynomial(coeffs=[6.52055,0.015002,-1.02611e-05,2.04713e-09,-1.37409e-13,1620.04,-15.2461], Tmin=(1894.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(53.4121,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo library: thermo_DFT_CCSDTF12_BAC + radical(OCOJ) + radical(OCOJ)"""),
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
    label = '[CH2]O[CH2](1706)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {3,S}
2 C u1 p0 c0 {1,S} {4,S} {5,S}
3 C u1 p0 c0 {1,S} {6,S} {7,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (177.918,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,393.551],'cm^-1')),
        HinderedRotor(inertia=(0.00108843,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00108842,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.05159,0.0165784,-1.00609e-07,-1.04835e-08,4.7882e-12,21436.5,11.5123], Tmin=(100,'K'), Tmax=(1037.89,'K')), NASAPolynomial(coeffs=[7.39316,0.00943956,-3.64811e-06,7.00835e-10,-5.09569e-14,20018.6,-12.084], Tmin=(1037.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(177.918,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(CsJOCH3) + radical(CsJOCH3)"""),
)

species(
    label = 'C1OCO1(1765)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,S} {4,S}
2 O u0 p2 c0 {3,S} {4,S}
3 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
4 C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-233.427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (60.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.95173,-0.0183985,0.000183532,-2.55042e-07,1.05038e-10,-27998.3,9.76036], Tmin=(100,'K'), Tmax=(910.358,'K')), NASAPolynomial(coeffs=[29.0324,-0.0252591,1.73212e-05,-3.3283e-09,2.13468e-13,-37211.1,-138.134], Tmin=(910.358,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-233.427,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(Cs-OsOsHH) + group(Cs-OsOsHH) + ring(Cyclobutane)"""),
)

species(
    label = 'COC=O(1722)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,S} {4,S}
2 O u0 p2 c0 {4,D}
3 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
4 C u0 p0 c0 {1,S} {2,D} {8,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-371.674,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.0105058,'amu*angstrom^2'), symmetry=1, barrier=(5.31906,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0943221,'amu*angstrom^2'), symmetry=1, barrier=(47.8214,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (60.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.28125,0.0113178,1.84485e-05,-2.41954e-08,7.75171e-12,-44672.6,12.0363], Tmin=(100,'K'), Tmax=(1141.76,'K')), NASAPolynomial(coeffs=[4.90974,0.0187532,-8.58342e-06,1.6685e-09,-1.18615e-13,-45900.9,0.213053], Tmin=(1141.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-371.674,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""methyl_formate""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2][O](1485)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {2,S}
2 C u1 p0 c0 {1,S} {3,S} {4,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (192.903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88413,-0.00363934,3.28562e-05,-4.13637e-08,1.59643e-11,23210.8,7.47969], Tmin=(100,'K'), Tmax=(933.045,'K')), NASAPolynomial(coeffs=[6.69323,0.000290203,8.6129e-07,-1.56321e-10,7.33527e-15,21991.3,-9.60361], Tmin=(933.045,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: FFCM1(-) + radical(H3COJ) + radical(CsJOH)"""),
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
    label = '[CH2]OC=O(1638)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {3,S} {4,S}
2 O u0 p2 c0 {4,D}
3 C u1 p0 c0 {1,S} {5,S} {6,S}
4 C u0 p0 c0 {1,S} {2,D} {7,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-170.733,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,1010.94],'cm^-1')),
        HinderedRotor(inertia=(0.395837,'amu*angstrom^2'), symmetry=1, barrier=(33.6604,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.46387,'amu*angstrom^2'), symmetry=1, barrier=(33.6572,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (59.044,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.79262,0.0222941,-9.08445e-06,-2.403e-09,1.82246e-12,-20487.5,13.1527], Tmin=(100,'K'), Tmax=(1184.75,'K')), NASAPolynomial(coeffs=[8.2588,0.0116085,-5.39242e-06,1.05479e-09,-7.5224e-14,-22328,-16.4464], Tmin=(1184.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-170.733,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), label="""CH2OCHO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]O[CH]O(1766)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,S} {4,S}
2 O u0 p2 c0 {3,S} {8,S}
3 C u1 p0 c0 {1,S} {2,S} {5,S}
4 C u1 p0 c0 {1,S} {6,S} {7,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {2,S}
"""),
    E0 = (-13.4029,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (60.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.198,0.0396779,-5.03099e-05,3.31155e-08,-8.53808e-12,-1547.15,16.2847], Tmin=(100,'K'), Tmax=(955.28,'K')), NASAPolynomial(coeffs=[9.27221,0.0100563,-3.79755e-06,6.55633e-10,-4.32053e-14,-2898.72,-17.5206], Tmin=(955.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-13.4029,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-OsOsHH) + group(Cs-OsHHH) + radical(OCJO) + radical(CsJOCH3)"""),
)

species(
    label = 'CO[CH][O](1767)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {3,S} {4,S}
2 O u1 p2 c0 {4,S}
3 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
4 C u1 p0 c0 {1,S} {2,S} {8,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (26.1241,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,180,1584.75,1586.78],'cm^-1')),
        HinderedRotor(inertia=(0.178718,'amu*angstrom^2'), symmetry=1, barrier=(4.10908,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.178409,'amu*angstrom^2'), symmetry=1, barrier=(4.10197,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.854,0.0354422,-7.02733e-05,8.23427e-08,-3.43378e-11,3173.34,13.6109], Tmin=(100,'K'), Tmax=(849.231,'K')), NASAPolynomial(coeffs=[-1.87232,0.0312402,-1.61085e-05,3.12781e-09,-2.16053e-13,4930.36,41.2588], Tmin=(849.231,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(26.1241,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-OsOsHH) + group(Cs-OsHHH) + radical(OCOJ) + radical(OCJO)"""),
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
    E0 = (-38.5193,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (371.228,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (357.398,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-30.235,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (24.8809,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (53.5263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-17.7096,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (322.253,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (65.0095,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (139.428,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]OC[O](1582)'],
    products = ['CH2O(19)', 'CH2O(19)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(T)(17)', '[O]C[O](1545)'],
    products = ['[CH2]OC[O](1582)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(87544.3,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['O(6)', '[CH2]O[CH2](1706)'],
    products = ['[CH2]OC[O](1582)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3335.46,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H2/O;O_birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]OC[O](1582)'],
    products = ['C1OCO1(1765)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]OC[O](1582)'],
    products = ['COC=O(1722)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2O(19)', '[CH2][O](1485)'],
    products = ['[CH2]OC[O](1582)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(8.90647e+06,'m^3/(mol*s)'), n=-0.198276, Ea=(43.2315,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3273107057168037, var=7.214547064829838, Tref=1000.0, N=427, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(5)', '[CH2]OC=O(1638)'],
    products = ['[CH2]OC[O](1582)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.42118e-08,'m^3/(mol*s)'), n=4.56427, Ea=(4.77331,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.27792374648988344, var=0.37052311097009455, Tref=1000.0, N=7, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_4R!H-u0',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_4R!H-u0"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2][O](1485)', '[CH2][O](1485)'],
    products = ['[CH2]OC[O](1582)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]OC[O](1582)'],
    products = ['[CH2]O[CH]O(1766)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.41377e+07,'s^-1'), n=1.57383, Ea=(103.529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;Cs_H_out_H/NonDeO] + [R2H_S;O_rad_out;Cs_H_out_1H] for rate rule [R2H_S;O_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CO[CH][O](1767)'],
    products = ['[CH2]OC[O](1582)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.71e+08,'s^-1'), n=1.45, Ea=(176.858,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_O;Y_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #655',
    isomers = [
        '[CH2]OC[O](1582)',
    ],
    reactants = [
        ('CH2O(19)', 'CH2O(19)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #655',
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

