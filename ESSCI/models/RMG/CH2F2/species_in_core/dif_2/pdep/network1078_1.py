species(
    label = '[CH]=CC(=[CH])F(2851)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 C u0 p0 c0 {3,S} {4,D} {6,S}
3 C u0 p0 c0 {1,S} {2,S} {5,D}
4 C u1 p0 c0 {2,D} {7,S}
5 C u1 p0 c0 {3,D} {8,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (404.209,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([250,446,589,854,899,3115,3125,620,680,785,800,1600,1700,180,1045.76,1045.81,1045.87],'cm^-1')),
        HinderedRotor(inertia=(0.161872,'amu*angstrom^2'), symmetry=1, barrier=(3.72175,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0649,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.15481,0.0341741,-2.22905e-05,2.10318e-09,1.94559e-12,48687,16.2324], Tmin=(100,'K'), Tmax=(1080.55,'K')), NASAPolynomial(coeffs=[11.9658,0.00946869,-4.11576e-06,8.36068e-10,-6.2431e-14,45888.8,-34.9968], Tmin=(1080.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(404.209,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Cds_P) + radical(Cdj(Cd-CdF1s)(H))"""),
)

species(
    label = 'C2HF(57)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 C u0 p0 c0 {3,T} {4,S}
3 C u0 p0 c0 {1,S} {2,T}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (95.331,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(44.0062,'amu')),
        LinearRotor(inertia=(51.6236,'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([429.793,429.793,596.357,596.357,1107.96,2365.05,3506.88],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0277,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1870.76,'J/mol'), sigma=(4.25,'angstroms'), dipoleMoment=(1,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.4498,0.0030263,3.99146e-05,-8.9615e-08,5.74336e-11,11468.6,5.90915], Tmin=(10,'K'), Tmax=(555.749,'K')), NASAPolynomial(coeffs=[4.23833,0.0086714,-5.87678e-06,1.96876e-09,-2.53031e-13,11206.2,0.995297], Tmin=(555.749,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(95.331,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(87.302,'J/(mol*K)'), label="""C#CF""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'FC1=CC=C1(2865)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 C u0 p0 c0 {3,S} {4,D} {6,S}
3 C u0 p0 c0 {2,S} {5,D} {7,S}
4 C u0 p0 c0 {2,D} {5,S} {8,S}
5 C u0 p0 c0 {1,S} {3,D} {4,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (219.465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (70.0649,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.94301,0.00320423,8.21834e-05,-1.46383e-07,7.9123e-11,26396.6,8.41104], Tmin=(10,'K'), Tmax=(596.917,'K')), NASAPolynomial(coeffs=[1.77313,0.0311244,-2.15992e-05,7.07667e-09,-8.75758e-13,26417.3,15.7632], Tmin=(596.917,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(219.465,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), label="""FC1DCCDC1""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[CH]=[C]F(530)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 C u1 p0 c0 {3,D} {4,S}
3 C u1 p0 c0 {1,S} {2,D}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (350.447,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([167,640,1190,1142.58,1502.03,3807.5],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0277,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.55591,0.0076799,-5.30098e-07,-4.54651e-09,2.16658e-12,42166.8,9.46363], Tmin=(100,'K'), Tmax=(1043.42,'K')), NASAPolynomial(coeffs=[5.86818,0.00351556,-1.29996e-06,2.62251e-10,-1.98915e-14,41428.4,-3.016], Tmin=(1043.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(350.447,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Cds_P) + radical(CdCdF1s)"""),
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
    label = '[CH]=C(F)C#C(1851)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {4,D}
3 C u0 p0 c0 {2,S} {5,T}
4 C u1 p0 c0 {2,D} {6,S}
5 C u0 p0 c0 {3,T} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (347.649,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([322,421,608,607,2175,525,3120,650,792.5,1650,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(1.60429,'amu*angstrom^2'), symmetry=1, barrier=(36.8858,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88681,0.00696908,9.16324e-05,-2.11268e-07,1.39251e-10,41819.5,9.68815], Tmin=(10,'K'), Tmax=(537.835,'K')), NASAPolynomial(coeffs=[5.3922,0.0203287,-1.41115e-05,4.69461e-09,-5.92628e-13,41302.4,0.0574202], Tmin=(537.835,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(347.649,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""[CH]DC(F)C#C""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[CH]=CC#C(2871)',
    structure = adjacencyList("""multiplicity 2
1 C u0 p0 c0 {2,S} {3,D} {5,S}
2 C u0 p0 c0 {1,S} {4,T}
3 C u1 p0 c0 {1,D} {6,S}
4 C u0 p0 c0 {2,T} {7,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (529.735,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,3120,650,792.5,1650,750,770,3400,2100,180,661.203,662.753,665.261],'cm^-1')),
        HinderedRotor(inertia=(0.0106145,'amu*angstrom^2'), symmetry=1, barrier=(3.32525,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (51.0665,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.75147,0.0228456,-8.30515e-06,-7.9138e-09,5.42953e-12,63761.5,11.8343], Tmin=(100,'K'), Tmax=(950.061,'K')), NASAPolynomial(coeffs=[8.93984,0.00798653,-2.52108e-06,4.30933e-10,-3.01739e-14,62080.4,-20.3631], Tmin=(950.061,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(529.735,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""HCCCHCH""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH]=C=C=[CH](1852)',
    structure = adjacencyList("""multiplicity 3
1 C u0 p0 c0 {2,D} {3,D}
2 C u0 p0 c0 {1,D} {4,D}
3 C u1 p0 c0 {1,D} {5,S}
4 C u1 p0 c0 {2,D} {6,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (614.901,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,563.333,586.667,610,1970,2140,3120,650,792.5,1650,180,1333.04],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (50.0585,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.53265,0.0269957,-3.0487e-05,1.76604e-08,-3.72328e-12,74012.9,13.5099], Tmin=(100,'K'), Tmax=(1422.35,'K')), NASAPolynomial(coeffs=[8.18718,0.00427875,6.57284e-07,-3.05855e-10,2.66578e-14,73093.7,-13.3387], Tmin=(1422.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(614.901,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=C=CJ) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]C(F)=C=C(2872)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {4,D} {5,S}
3 C u0 p0 c0 {4,D} {6,S} {7,S}
4 C u0 p0 c0 {2,D} {3,D}
5 C u2 p0 c0 {2,S} {8,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (330.845,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (70.0649,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.49574,0.0331483,-2.03113e-05,6.20312e-09,-7.87555e-13,39845.1,15.2364], Tmin=(100,'K'), Tmax=(1712.17,'K')), NASAPolynomial(coeffs=[8.39363,0.0193697,-8.24022e-06,1.50307e-09,-1.01293e-13,37825.4,-16.3893], Tmin=(1712.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(330.845,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C=C([CH2])F(1518)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {4,D}
3 C u1 p0 c0 {2,S} {6,S} {7,S}
4 C u0 p0 c0 {2,D} {5,D}
5 C u1 p0 c0 {4,D} {8,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (266.137,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,540,610,2055,3120,650,792.5,1650,180,310.756,1088.14,1088.25],'cm^-1')),
        HinderedRotor(inertia=(2.71384,'amu*angstrom^2'), symmetry=1, barrier=(62.3966,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0649,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.18603,0.0354707,-3.30721e-05,1.62519e-08,-3.1442e-12,32077.9,16.1523], Tmin=(100,'K'), Tmax=(1264.08,'K')), NASAPolynomial(coeffs=[10.1277,0.0103402,-3.25125e-06,5.24542e-10,-3.37397e-14,30070.1,-24.0229], Tmin=(1264.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(266.137,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Allyl_P) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=C=C[CH]F(2870)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 C u0 p0 c0 {3,S} {4,D} {6,S}
3 C u1 p0 c0 {1,S} {2,S} {7,S}
4 C u0 p0 c0 {2,D} {5,D}
5 C u1 p0 c0 {4,D} {8,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (259.991,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (70.0649,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28673,0.0334736,-2.96111e-05,1.40662e-08,-2.65621e-12,31334.9,16.8008], Tmin=(100,'K'), Tmax=(1287.46,'K')), NASAPolynomial(coeffs=[9.48251,0.0111172,-3.56419e-06,5.78752e-10,-3.7236e-14,29482,-19.7329], Tmin=(1287.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(259.991,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCdF1sH) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]C=C=CF(2868)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 C u0 p0 c0 {4,D} {5,S} {6,S}
3 C u0 p0 c0 {1,S} {4,D} {7,S}
4 C u0 p0 c0 {2,D} {3,D}
5 C u2 p0 c0 {2,S} {8,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (344.912,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (70.0649,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.43928,0.0318872,-1.81679e-05,5.00046e-09,-5.54693e-13,41541.4,16.1445], Tmin=(100,'K'), Tmax=(2008.56,'K')), NASAPolynomial(coeffs=[10.2725,0.016287,-6.51742e-06,1.13343e-09,-7.33628e-14,38394.7,-27.1093], Tmin=(2008.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(344.912,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(AllylJ2_triplet)"""),
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
    E0 = (89.8215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (98.1058,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (263.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (255.174,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (347.465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (374.793,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (626.818,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (319.559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (257.78,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (281.95,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (251.729,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (243.874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=CC(=[CH])F(2851)'],
    products = ['C2HF(57)', 'C2H2(70)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=CC(=[CH])F(2851)'],
    products = ['FC1=CC=C1(2865)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_DSD;CdsingleH_rad_out;CdsinglepriH_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=[C]F(530)', 'C2H2(70)'],
    products = ['[CH]=CC(=[CH])F(2851)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(0.0486,'m^3/(mol*s)'), n=2.41, Ea=(9.26598,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_N-Sp-2R!H=1R!H_Ext-3R-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_N-Sp-2R!H=1R!H_Ext-3R-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(5)', '[CH]=C(F)C#C(1851)'],
    products = ['[CH]=CC(=[CH])F(2851)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1348.96,'m^3/(mol*s)'), n=1.432, Ea=(10.1066,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_Sp-5R!H-4CClOS_Ext-4CClOS-R_N-Sp-6R!H-4CClOS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_N-4R!H->N_Ext-4CClOS-R_Sp-5R!H-4CClOS_Ext-4CClOS-R_N-Sp-6R!H-4CClOS"""),
)

reaction(
    label = 'reaction5',
    reactants = ['F(37)', '[CH]=CC#C(2871)'],
    products = ['[CH]=CC(=[CH])F(2851)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(59.2256,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C2HF(57)', '[CH]=[CH](2414)'],
    products = ['[CH]=CC(=[CH])F(2851)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.201622,'m^3/(mol*s)'), n=2.45643, Ea=(3.09123,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.17813889912811182, var=0.789831941087029, Tref=1000.0, N=34, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_N-Sp-2R!H=1R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_N-Sp-2R!H=1R!H
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=[C]F(530)', '[CH]=[CH](2414)'],
    products = ['[CH]=CC(=[CH])F(2851)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(408.814,'m^3/(mol*s)'), n=1.30994, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.006495679117668562, var=10.327386067777251, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_N-4R!H->O_N-4BrCClF->F',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_N-4R!H->O_N-4BrCClF->F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction8',
    reactants = ['HF(38)', '[CH]=C=C=[CH](1852)'],
    products = ['[CH]=CC(=[CH])F(2851)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(5353.26,'m^3/(mol*s)'), n=0.732206, Ea=(300.158,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=CC(=[CH])F(2851)'],
    products = ['[CH]C(F)=C=C(2872)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(481900,'s^-1'), n=2.375, Ea=(167.958,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 121 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleDe
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=CC(=[CH])F(2851)'],
    products = ['[CH]=C=C([CH2])F(1518)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=CC(=[CH])F(2851)'],
    products = ['[CH]=C=C[CH]F(2870)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(161.907,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=CC(=[CH])F(2851)'],
    products = ['[CH]C=C=CF(2868)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(154.053,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

network(
    label = 'PDepNetwork #1078',
    isomers = [
        '[CH]=CC(=[CH])F(2851)',
    ],
    reactants = [
        ('C2HF(57)', 'C2H2(70)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1078',
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

