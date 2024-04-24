species(
    label = '[O]C(F)(C=O)O[C]=O(3136)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {6,S} {8,S}
3 O u1 p2 c0 {6,S}
4 O u0 p2 c0 {7,D}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
7 C u0 p0 c0 {4,D} {6,S} {9,S}
8 C u1 p0 c0 {2,S} {5,D}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-471.129,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.65271,'amu*angstrom^2'), symmetry=1, barrier=(37.999,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.65258,'amu*angstrom^2'), symmetry=1, barrier=(37.996,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.65282,'amu*angstrom^2'), symmetry=1, barrier=(38.0017,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05475,0.0710608,-0.000117336,1.02145e-07,-3.51246e-11,-56563.7,24.5811], Tmin=(100,'K'), Tmax=(796.233,'K')), NASAPolynomial(coeffs=[9.26934,0.0224632,-1.19746e-05,2.36675e-09,-1.65992e-13,-57639.4,-11.7186], Tmin=(796.233,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-471.129,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(Cds-OdOsH) + radical(C=OCOJ) + radical((O)CJOC)"""),
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
    label = '[O]C(F)O[C]=O(576)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {5,S} {6,S}
3 O u1 p2 c0 {5,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
6 C u1 p0 c0 {2,S} {4,D}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-345.883,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([451,553,637,1069,1180,1265,1301,3056,1855,455,950,372.978,996.618],'cm^-1')),
        HinderedRotor(inertia=(2.2228,'amu*angstrom^2'), symmetry=1, barrier=(51.1066,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.22354,'amu*angstrom^2'), symmetry=1, barrier=(51.1235,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0259,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3728.91,'J/mol'), sigma=(5.84942,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=582.45 K, Pc=42.28 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.38443,0.0331248,-2.9166e-05,1.25996e-08,-2.17133e-12,-41540.1,18.0887], Tmin=(100,'K'), Tmax=(1381.77,'K')), NASAPolynomial(coeffs=[10.0482,0.0109393,-5.08175e-06,9.79413e-10,-6.88921e-14,-43657.9,-21.3624], Tmin=(1381.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-345.883,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(CsFHOO) + group(Cds-OdOsH) + radical(O2sj(Cs-F1sO2sH)) + radical((O)CJOC)"""),
)

species(
    label = 'O=[C]OC(=O)[CH]OF(3166)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {6,S} {8,S}
3 O u0 p2 c0 {1,S} {7,S}
4 O u0 p2 c0 {6,D}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {2,S} {4,D} {7,S}
7 C u1 p0 c0 {3,S} {6,S} {9,S}
8 C u1 p0 c0 {2,S} {5,D}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-216.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3025,407.5,1350,352.5,1855,455,950,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.789271,0.0779211,-0.000136618,1.19471e-07,-4.02333e-11,-25914.7,27.4857], Tmin=(100,'K'), Tmax=(839.588,'K')), NASAPolynomial(coeffs=[10.5344,0.0193709,-1.03551e-05,2.01685e-09,-1.39004e-13,-27123.8,-15.2807], Tmin=(839.588,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-216.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(191.233,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-O2d)) + group(O2sCF) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsOs) + group(Cds-OdOsH) + radical(OCJC=O) + radical((O)CJOC)"""),
)

species(
    label = 'O=[C]OO[CH]C(=O)F(2966)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {3,S} {6,S}
3 O u0 p2 c0 {2,S} {8,S}
4 O u0 p2 c0 {7,D}
5 O u0 p2 c0 {8,D}
6 C u1 p0 c0 {2,S} {7,S} {9,S}
7 C u0 p0 c0 {1,S} {4,D} {6,S}
8 C u1 p0 c0 {3,S} {5,D}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-400.604,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3025,407.5,1350,352.5,611,648,830,1210,1753,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(1.92943,'amu*angstrom^2'), symmetry=1, barrier=(44.3615,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.93265,'amu*angstrom^2'), symmetry=1, barrier=(44.4355,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.92835,'amu*angstrom^2'), symmetry=1, barrier=(44.3365,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.93285,'amu*angstrom^2'), symmetry=1, barrier=(44.44,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3856.9,'J/mol'), sigma=(5.85905,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=602.44 K, Pc=43.51 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.81019,0.0621594,-6.34708e-05,3.01558e-08,-5.51998e-12,-48060,26.3492], Tmin=(100,'K'), Tmax=(1338.64,'K')), NASAPolynomial(coeffs=[18.3076,0.00987488,-4.8835e-06,9.78125e-10,-7.08111e-14,-52744.5,-63.1688], Tmin=(1338.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-400.604,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(191.233,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + group(Cds-OdOsH) + radical(OCJC=O) + radical((O)CJOC)"""),
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
    label = '[O]C([O])(F)C=O(3041)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u1 p2 c0 {5,S}
3 O u1 p2 c0 {5,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
6 C u0 p0 c0 {4,D} {5,S} {7,S}
7 H u0 p0 c0 {6,S}
"""),
    E0 = (-272.797,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.620008,'amu*angstrom^2'), symmetry=1, barrier=(14.2552,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0259,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9391,0.0520954,-9.61723e-05,8.85863e-08,-3.07158e-11,-32742.1,19.1588], Tmin=(100,'K'), Tmax=(871.375,'K')), NASAPolynomial(coeffs=[6.59297,0.0158185,-8.05219e-06,1.52657e-09,-1.02949e-13,-32987,0.596161], Tmin=(871.375,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-272.797,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(C=OCOJ)"""),
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
    label = 'O=[C]O[C](F)C=O(2659)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {5,S} {7,S}
3 O u0 p2 c0 {6,D}
4 O u0 p2 c0 {7,D}
5 C u1 p0 c0 {1,S} {2,S} {6,S}
6 C u0 p0 c0 {3,D} {5,S} {8,S}
7 C u1 p0 c0 {2,S} {4,D}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-367.182,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([280,501,1494,1531,2782.5,750,1395,475,1775,1000,1855,455,950,430.83,430.986],'cm^-1')),
        HinderedRotor(inertia=(0.400039,'amu*angstrom^2'), symmetry=1, barrier=(52.8047,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.400454,'amu*angstrom^2'), symmetry=1, barrier=(52.8011,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.400272,'amu*angstrom^2'), symmetry=1, barrier=(52.8029,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.30428,0.0414898,-4.59937e-05,3.10574e-08,-9.28818e-12,-44104.4,18.8058], Tmin=(100,'K'), Tmax=(783.828,'K')), NASAPolynomial(coeffs=[5.74318,0.0239417,-1.24143e-05,2.49899e-09,-1.80131e-13,-44643.6,3.05257], Tmin=(783.828,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-367.182,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(Cds-OdOsH) + radical(CsCOF1sO2s) + radical((O)CJOC)"""),
)

species(
    label = 'O=CC1(F)OC(=O)O1(3167)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {6,S} {8,S}
3 O u0 p2 c0 {6,S} {8,S}
4 O u0 p2 c0 {7,D}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
7 C u0 p0 c0 {4,D} {6,S} {9,S}
8 C u0 p0 c0 {2,S} {3,S} {5,D}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-814.503,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.04605,0.0245677,4.08116e-05,-7.51649e-08,3.06216e-11,-97875.3,20.1356], Tmin=(100,'K'), Tmax=(993.627,'K')), NASAPolynomial(coeffs=[17.0707,0.0093329,-4.49904e-06,1.06762e-09,-9.02432e-14,-103095,-63.494], Tmin=(993.627,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-814.503,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-Cs(Cds-O2d)) + group(CsCFOO) + group(Cds-OdCsH) + group(Cds-OdOsOs) + ring(Cyclobutane)"""),
)

species(
    label = 'O=[C]OC1(F)[CH]OO1(3168)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {4,S} {6,S}
3 O u0 p2 c0 {6,S} {8,S}
4 O u0 p2 c0 {2,S} {7,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
7 C u1 p0 c0 {4,S} {6,S} {9,S}
8 C u1 p0 c0 {3,S} {5,D}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-214.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72815,0.0548652,-6.60169e-05,4.35781e-08,-1.1995e-11,-25673.9,21.563], Tmin=(100,'K'), Tmax=(868.639,'K')), NASAPolynomial(coeffs=[8.45266,0.0238994,-1.25438e-05,2.53826e-09,-1.83448e-13,-26842.1,-9.93189], Tmin=(868.639,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-214.11,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-OsCs) + group(O2s-OsCs) + group(CsCFOO) + group(Cs-CsOsHH) + group(Cds-OdOsH) + ring(Cs-Cs(F)-O2s-O2s) + radical(CCsJOO) + radical((O)CJOC)"""),
)

species(
    label = '[O]C1(F)[CH]OC(=O)O1(3169)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {6,S} {8,S}
3 O u0 p2 c0 {7,S} {8,S}
4 O u1 p2 c0 {6,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
7 C u1 p0 c0 {3,S} {6,S} {9,S}
8 C u0 p0 c0 {2,S} {3,S} {5,D}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-522.487,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73449,0.0327237,2.09362e-05,-5.97165e-08,2.69414e-11,-62743.7,21.5276], Tmin=(100,'K'), Tmax=(970.084,'K')), NASAPolynomial(coeffs=[18.4842,0.00556111,-1.85516e-06,4.72874e-10,-4.49473e-14,-67965.1,-68.9335], Tmin=(970.084,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-522.487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(CsCFOO) + group(Cs-CsOsHH) + group(Cds-OdOsOs) + ring(Cyclopentane) + radical(O2sj(Cs-F1sO2sCs)) + radical(CCsJOC(O))"""),
)

species(
    label = '[O]C1OC1(F)O[C]=O(3170)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {6,S} {7,S}
3 O u0 p2 c0 {6,S} {8,S}
4 O u1 p2 c0 {7,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
7 C u0 p0 c0 {2,S} {4,S} {6,S} {9,S}
8 C u1 p0 c0 {3,S} {5,D}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-412.609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.15738,0.0465333,-4.56434e-05,2.54255e-08,-6.34649e-12,-49564.1,21.3799], Tmin=(100,'K'), Tmax=(911.673,'K')), NASAPolynomial(coeffs=[6.38472,0.0279857,-1.51265e-05,3.10987e-09,-2.27057e-13,-50334.9,1.37638], Tmin=(911.673,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-412.609,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(CsCFOO) + group(Cs-CsOsOsH) + group(Cds-OdOsH) + ring(Cs(F)(O2)-O2s-Cs) + radical(CCOJ) + radical((O)CJOC)"""),
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
    label = 'O=[C]OC(=O)C=O(3172)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {5,S} {7,S}
2 O u0 p2 c0 {5,D}
3 O u0 p2 c0 {6,D}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {1,S} {2,D} {6,S}
6 C u0 p0 c0 {3,D} {5,S} {8,S}
7 C u1 p0 c0 {1,S} {4,D}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-382.325,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,1855,455,950,259.074,259.116,259.173,259.193,259.195,1954.19],'cm^-1')),
        HinderedRotor(inertia=(0.00250927,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.489059,'amu*angstrom^2'), symmetry=1, barrier=(23.306,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00251,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.037,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59042,0.0592769,-0.000103303,9.37539e-08,-3.31762e-11,-45902.4,22.8515], Tmin=(100,'K'), Tmax=(809.262,'K')), NASAPolynomial(coeffs=[7.66167,0.019266,-1.06015e-05,2.11454e-09,-1.48635e-13,-46557.5,-3.13035], Tmin=(809.262,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-382.325,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-O2d)) + group(Cds-O2d(Cds-O2d)O2s) + group(Cds-O2d(Cds-O2d)H) + group(Cds-OdOsH) + radical((O)CJOC)"""),
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
    label = 'HCO(14)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,D}
2 C u1 p0 c0 {1,D} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (32.4782,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1131.19,1955.83,1955.83],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.23755,-0.00332075,1.4003e-05,-1.3424e-08,4.37416e-12,3872.41,3.30835], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.92002,0.00252279,-6.71004e-07,1.05616e-10,-7.43798e-15,3653.43,3.58077], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(32.4782,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HCO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'O=[C]OC(=O)F(284)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {5,S} {6,S}
3 O u0 p2 c0 {5,D}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {2,S} {3,D}
6 C u1 p0 c0 {2,S} {4,D}
"""),
    E0 = (-510.842,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([482,664,788,1296,1923,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.284615,'amu*angstrom^2'), symmetry=1, barrier=(6.54387,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.28454,'amu*angstrom^2'), symmetry=1, barrier=(6.54214,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (91.0179,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.16703,0.043103,-7.09532e-05,5.65753e-08,-1.75004e-11,-61376.7,19.026], Tmin=(100,'K'), Tmax=(799.084,'K')), NASAPolynomial(coeffs=[9.27835,0.00750518,-4.12988e-06,8.24649e-10,-5.80989e-14,-62513.2,-13.6869], Tmin=(799.084,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-510.842,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-O2d)) + group(COFOO) + group(Cds-OdOsH) + radical((O)CJOC)"""),
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
    label = 'O=[C]OC(=O)[C]=O(3173)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {5,S} {7,S}
2 O u0 p2 c0 {5,D}
3 O u0 p2 c0 {6,D}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {1,S} {2,D} {6,S}
6 C u1 p0 c0 {3,D} {5,S}
7 C u1 p0 c0 {1,S} {4,D}
"""),
    E0 = (-223.798,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1850,1860,440,470,900,1000,310.927,310.928,310.929,310.929,2053.66,4000],'cm^-1')),
        HinderedRotor(inertia=(0.189947,'amu*angstrom^2'), symmetry=1, barrier=(13.0312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.366287,'amu*angstrom^2'), symmetry=1, barrier=(25.1284,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.366283,'amu*angstrom^2'), symmetry=1, barrier=(25.1284,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56141,0.0592819,-0.000107157,9.41765e-08,-3.17272e-11,-26834.3,23.7331], Tmin=(100,'K'), Tmax=(836.047,'K')), NASAPolynomial(coeffs=[9.40999,0.0126356,-7.14705e-06,1.41604e-09,-9.82691e-14,-27828.8,-10.8253], Tmin=(836.047,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-223.798,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(145.503,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-O2d)) + group(Cds-O2d(Cds-O2d)O2s) + group(Cds-O2d(Cds-O2d)H) + group(Cds-OdOsH) + radical(OC=OCJ=O) + radical((O)CJOC)"""),
)

species(
    label = 'O=[C]OC(O)(F)[C]=O(3174)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {6,S} {8,S}
3 O u0 p2 c0 {6,S} {9,S}
4 O u0 p2 c0 {7,D}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
7 C u1 p0 c0 {4,D} {6,S}
8 C u1 p0 c0 {2,S} {5,D}
9 H u0 p0 c0 {3,S}
"""),
    E0 = (-560.191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.916695,0.0737063,-0.00012069,1.02247e-07,-3.43401e-11,-67270,24.9767], Tmin=(100,'K'), Tmax=(766.773,'K')), NASAPolynomial(coeffs=[10.553,0.0205946,-1.123e-05,2.24274e-09,-1.5843e-13,-68664.2,-18.409], Tmin=(766.773,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-560.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(191.233,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(Cds-OdOsH) + radical(CsCJ=O) + radical((O)CJOC)"""),
)

species(
    label = '[O]C(F)([C]=O)OC=O(3175)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {6,S} {7,S}
3 O u1 p2 c0 {6,S}
4 O u0 p2 c0 {7,D}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
7 C u0 p0 c0 {2,S} {4,D} {9,S}
8 C u1 p0 c0 {5,D} {6,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-512.906,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42428,0.0622731,-9.34844e-05,7.80071e-08,-2.67037e-11,-61600.9,25.3127], Tmin=(100,'K'), Tmax=(710.352,'K')), NASAPolynomial(coeffs=[8.1737,0.0242656,-1.32241e-05,2.68002e-09,-1.92349e-13,-62559.7,-4.94089], Tmin=(710.352,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-512.906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(CsCFOO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(Cds-OdOsH) + radical(C=OCOJ) + radical(CsCJ=O)"""),
)

species(
    label = '[O]C=C(OF)O[C]=O(3176)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {6,S} {8,S}
3 O u0 p2 c0 {1,S} {6,S}
4 O u1 p2 c0 {7,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {2,S} {3,S} {7,D}
7 C u0 p0 c0 {4,S} {6,D} {9,S}
8 C u1 p0 c0 {2,S} {5,D}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-164.412,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,350,440,435,1725,3010,987.5,1337.5,450,1655,1855,455,950,210.45,210.451,210.451,210.454],'cm^-1')),
        HinderedRotor(inertia=(0.884025,'amu*angstrom^2'), symmetry=1, barrier=(27.784,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.88402,'amu*angstrom^2'), symmetry=1, barrier=(27.784,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.884042,'amu*angstrom^2'), symmetry=1, barrier=(27.784,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.57501,0.0794504,-0.000121575,8.96833e-08,-2.56495e-11,-19654.5,23.58], Tmin=(100,'K'), Tmax=(863.398,'K')), NASAPolynomial(coeffs=[14.7697,0.0136869,-7.31999e-06,1.46034e-09,-1.03687e-13,-22105.5,-42.816], Tmin=(863.398,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-164.412,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2sCF) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-OdOsH) + radical(C=COJ) + radical((O)CJOC)"""),
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
    E0 = (-199.723,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (112.379,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (251.123,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (67.8231,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (437.694,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (147.258,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-191.439,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (57.2951,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-154.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-97.4326,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-119.246,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-111.572,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-187.855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-20.5373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-132.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-107.058,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (88.4638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (47.6179,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-124.442,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-69.5986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (185.943,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-31.2407,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(F)(C=O)O[C]=O(3136)'],
    products = ['CO2(13)', 'O=CC(=O)F(557)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(12)', '[O]C(F)O[C]=O(576)'],
    products = ['[O]C(F)(C=O)O[C]=O(3136)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.81927e-16,'m^3/(mol*s)'), n=6.80628, Ea=(305.597,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_N-2Br1sCl1sF1sH->F1s_N-5R!H->C',), comment="""Estimated from node Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_N-2Br1sCl1sF1sH->F1s_N-5R!H->C"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O=[C]OC(=O)[CH]OF(3166)'],
    products = ['[O]C(F)(C=O)O[C]=O(3136)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(196.087,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O=[C]OO[CH]C(=O)F(2966)'],
    products = ['[O]C(F)(C=O)O[C]=O(3136)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(197.022,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[C]=O(222)', '[O]C([O])(F)C=O(3041)'],
    products = ['[O]C(F)(C=O)O[C]=O(3136)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(87544.3,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['O(6)', 'O=[C]O[C](F)C=O(2659)'],
    products = ['[O]C(F)(C=O)O[C]=O(3136)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_ter_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C(F)(C=O)O[C]=O(3136)'],
    products = ['O=CC1(F)OC(=O)O1(3167)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C(F)(C=O)O[C]=O(3136)'],
    products = ['O=[C]OC1(F)[CH]OO1(3168)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(257.018,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 255.4 to 257.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C(F)(C=O)O[C]=O(3136)'],
    products = ['[O]C1(F)[CH]OC(=O)O1(3169)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.12332e+09,'s^-1'), n=0.5388, Ea=(45.5334,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra] for rate rule [R5_SS_CO;carbonyl_intra_H;radadd_intra_CO]
Euclidian distance = 2.449489742783178
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C(F)(C=O)O[C]=O(3136)'],
    products = ['[O]C1OC1(F)O[C]=O(3170)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.785e+11,'s^-1'), n=0.342, Ea=(102.29,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]C(F)(C=O)O[C]=O(3136)'],
    products = ['[O]C1C(=O)OC1([O])F(3171)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.98674e+07,'s^-1'), n=1.31443, Ea=(80.4773,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra] for rate rule [R5_SS_CO;carbonylbond_intra_H;radadd_intra_CO]
Euclidian distance = 2.449489742783178
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CO(12)', '[O]C([O])(F)C=O(3041)'],
    products = ['[O]C(F)(C=O)O[C]=O(3136)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(68.2,'m^3/(mol*s)'), n=8.73864e-09, Ea=(8.56048,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->O',), comment="""Estimated from node Root_3R->O
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CO2(13)', '[O][CH]C(=O)F(438)'],
    products = ['[O]C(F)(C=O)O[C]=O(3136)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.81516e-05,'m^3/(mol*s)'), n=3.04336, Ea=(158.354,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3757377757886876, var=2.242054186761003, Tref=1000.0, N=1042, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F(37)', 'O=[C]OC(=O)C=O(3172)'],
    products = ['[O]C(F)(C=O)O[C]=O(3136)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.08261e+20,'m^3/(mol*s)'), n=-5.07836, Ea=(17.4907,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O][C]=O(188)', 'O=CC(=O)F(557)'],
    products = ['[O]C(F)(C=O)O[C]=O(3136)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.82043e-07,'m^3/(mol*s)'), n=3.71185, Ea=(48.0996,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction16',
    reactants = ['HCO(14)', 'O=[C]OC(=O)F(284)'],
    products = ['[O]C(F)(C=O)O[C]=O(3136)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=-1.07934e-11, Ea=(99.9,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O][C]=O(188)', '[O][CH]C(=O)F(438)'],
    products = ['[O]C(F)(C=O)O[C]=O(3136)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.808e+07,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction18',
    reactants = ['HF(38)', 'O=[C]OC(=O)[C]=O(3173)'],
    products = ['[O]C(F)(C=O)O[C]=O(3136)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(281.124,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O]C(F)(C=O)O[C]=O(3136)'],
    products = ['O=[C]OC(O)(F)[C]=O(3174)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;O_rad_out;XH_out] for rate rule [R3H_SS_Cs;O_rad_out;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O]C(F)(C=O)O[C]=O(3136)'],
    products = ['[O]C(F)([C]=O)OC=O(3175)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(9.82783e+06,'s^-1'), n=1.57155, Ea=(130.124,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS_OCs;Y_rad_out;XH_out] for rate rule [R4H_SSS_OCs;CO_rad_out;CO_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]C=C(OF)O[C]=O(3176)'],
    products = ['[O]C(F)(C=O)O[C]=O(3136)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(78.9489,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O]C(F)(C=O)O[C]=O(3136)'],
    products = ['[O]C=C([O])OC(=O)F(3177)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(168.482,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

network(
    label = 'PDepNetwork #1197',
    isomers = [
        '[O]C(F)(C=O)O[C]=O(3136)',
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
    label = 'PDepNetwork #1197',
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

