species(
    label = 'C=C(OF)OC=O(3338)',
    structure = adjacencyList("""1  F u0 p3 c0 {3,S}
2  O u0 p2 c0 {5,S} {7,S}
3  O u0 p2 c0 {1,S} {5,S}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {2,S} {3,S} {6,D}
6  C u0 p0 c0 {5,D} {8,S} {9,S}
7  C u0 p0 c0 {2,S} {4,D} {10,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-293.53,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,350,440,435,1725,2950,3100,1380,975,1025,1650,2782.5,750,1395,475,1775,1000,274.163,274.179,1951.7],'cm^-1')),
        HinderedRotor(inertia=(0.54554,'amu*angstrom^2'), symmetry=1, barrier=(29.111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.189,'amu*angstrom^2'), symmetry=1, barrier=(10.0877,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.545804,'amu*angstrom^2'), symmetry=1, barrier=(29.113,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70957,0.0582504,-6.88155e-05,2.90522e-08,7.09059e-12,-35228.7,19.6709], Tmin=(100,'K'), Tmax=(535.393,'K')), NASAPolynomial(coeffs=[7.52885,0.0267623,-1.41842e-05,2.84962e-09,-2.03665e-13,-36023.7,-6.37281], Tmin=(535.393,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-293.53,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2sCF) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cds-OdOsH)"""),
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
    label = 'C=C(O)OF(3492)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {4,S} {8,S}
3 O u0 p2 c0 {1,S} {4,S}
4 C u0 p0 c0 {2,S} {3,S} {5,D}
5 C u0 p0 c0 {4,D} {6,S} {7,S}
6 H u0 p0 c0 {5,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {2,S}
"""),
    E0 = (-200.378,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,231,791,350,440,435,1725,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(1.31372,'amu*angstrom^2'), symmetry=1, barrier=(30.205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.31042,'amu*angstrom^2'), symmetry=1, barrier=(30.129,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.0424,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.87494,0.00722414,0.000106356,-2.21274e-07,1.3152e-10,-24092.6,10.0683], Tmin=(10,'K'), Tmax=(591.099,'K')), NASAPolynomial(coeffs=[5.13645,0.0281847,-2.16884e-05,7.56312e-09,-9.70859e-13,-24757.1,0.286506], Tmin=(591.099,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-200.378,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""CDC(O)OF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'C=COF(3283)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {3,S}
3 C u0 p0 c0 {2,S} {4,D} {5,S}
4 C u0 p0 c0 {3,D} {6,S} {7,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-4.99841,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,536.063],'cm^-1')),
        HinderedRotor(inertia=(0.610082,'amu*angstrom^2'), symmetry=1, barrier=(14.027,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (62.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.84473,0.0161967,-8.81466e-08,-9.83866e-09,4.52656e-12,-601.523,8.34737], Tmin=(10,'K'), Tmax=(953.713,'K')), NASAPolynomial(coeffs=[5.07492,0.0164695,-9.06135e-06,2.40633e-09,-2.48859e-13,-1083.23,1.17548], Tmin=(953.713,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-4.99841,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CDCOF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=CCC(=O)OF(3343)',
    structure = adjacencyList("""1  F u0 p3 c0 {2,S}
2  O u0 p2 c0 {1,S} {6,S}
3  O u0 p2 c0 {6,D}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {6,S} {7,S} {8,S} {9,S}
6  C u0 p0 c0 {2,S} {3,D} {5,S}
7  C u0 p0 c0 {4,D} {5,S} {10,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-404.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.12613,0.044749,-3.68511e-05,1.57935e-08,-2.89113e-12,-48529.7,20.5468], Tmin=(100,'K'), Tmax=(1230.91,'K')), NASAPolynomial(coeffs=[8.36915,0.0244616,-1.21287e-05,2.40385e-09,-1.71668e-13,-50066.6,-10.8692], Tmin=(1230.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-404.034,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-O2d)(Cds-O2d)HH) + group(Cds-OdCsOs) + group(Cds-OdCsH)"""),
)

species(
    label = 'O=COC(=O)CF(3502)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  O u0 p2 c0 {6,S} {7,S}
3  O u0 p2 c0 {6,D}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
6  C u0 p0 c0 {2,S} {3,D} {5,S}
7  C u0 p0 c0 {2,S} {4,D} {10,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-704.424,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.09605,0.0463704,-4.95504e-05,3.08262e-08,-8.31546e-12,-84657.9,23.156], Tmin=(100,'K'), Tmax=(869.547,'K')), NASAPolynomial(coeffs=[6.74547,0.0249828,-1.26562e-05,2.54018e-09,-1.83112e-13,-85466.5,1.37513], Tmin=(869.547,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-704.424,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-O2d)) + group(CsCFHH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsOs) + group(Cds-OdOsH)"""),
)

species(
    label = 'FO[C]1CO[CH]O1(3503)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {4,S}
2  O u0 p2 c0 {5,S} {7,S}
3  O u0 p2 c0 {6,S} {7,S}
4  O u0 p2 c0 {1,S} {6,S}
5  C u0 p0 c0 {2,S} {6,S} {8,S} {9,S}
6  C u1 p0 c0 {3,S} {4,S} {5,S}
7  C u1 p0 c0 {2,S} {3,S} {10,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (2.48387,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2950,3150,900,1000,1100,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.2617,0.0615877,-7.73456e-05,4.45331e-08,-8.96576e-12,451.21,23.445], Tmin=(100,'K'), Tmax=(1542.46,'K')), NASAPolynomial(coeffs=[15.3627,-0.00113794,6.56953e-06,-1.64044e-09,1.23295e-13,-1404.05,-46.8659], Tmin=(1542.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.48387,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2sCF) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-OsOsHH) + ring(1,3-Dioxolane) + radical(Cs_P) + radical(OCJO)"""),
)

species(
    label = '[CH2]C(OF)O[C]=O(3504)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  O u0 p2 c0 {5,S} {7,S}
3  O u0 p2 c0 {1,S} {5,S}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {2,S} {3,S} {6,S} {8,S}
6  C u1 p0 c0 {5,S} {9,S} {10,S}
7  C u1 p0 c0 {2,S} {4,D}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-65.3036,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,1855,455,950,358.443,358.493,358.507],'cm^-1')),
        HinderedRotor(inertia=(0.182535,'amu*angstrom^2'), symmetry=1, barrier=(16.6436,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.459067,'amu*angstrom^2'), symmetry=1, barrier=(41.8494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00754448,'amu*angstrom^2'), symmetry=1, barrier=(24.6864,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.458955,'amu*angstrom^2'), symmetry=1, barrier=(41.8493,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38724,0.0615131,-7.52968e-05,4.82186e-08,-1.25168e-11,-7763.52,23.5104], Tmin=(100,'K'), Tmax=(930.48,'K')), NASAPolynomial(coeffs=[10.6491,0.0216973,-1.11104e-05,2.23025e-09,-1.60639e-13,-9487.11,-20.5053], Tmin=(930.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-65.3036,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2sCF) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-OdOsH) + radical(CJCO) + radical((O)CJOC)"""),
)

species(
    label = 'C[C](OF)O[C]=O(3505)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  O u0 p2 c0 {6,S} {7,S}
3  O u0 p2 c0 {1,S} {6,S}
4  O u0 p2 c0 {7,D}
5  C u0 p0 c0 {6,S} {8,S} {9,S} {10,S}
6  C u1 p0 c0 {2,S} {3,S} {5,S}
7  C u1 p0 c0 {2,S} {4,D}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-71.6463,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,1855,455,950,297.125,297.141,297.186],'cm^-1')),
        HinderedRotor(inertia=(0.00190892,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.689153,'amu*angstrom^2'), symmetry=1, barrier=(43.1787,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0428473,'amu*angstrom^2'), symmetry=1, barrier=(43.1786,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.273887,'amu*angstrom^2'), symmetry=1, barrier=(17.1566,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62134,0.0565264,-6.41137e-05,3.88162e-08,-9.70452e-12,-8534.91,22.7247], Tmin=(100,'K'), Tmax=(955.083,'K')), NASAPolynomial(coeffs=[9.55216,0.0233115,-1.19487e-05,2.40443e-09,-1.7355e-13,-10049.8,-15.1725], Tmin=(955.083,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-71.6463,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2sCF) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-OdOsH) + radical(Cs_P) + radical((O)CJOC)"""),
)

species(
    label = '[CH]=C(OF)OC[O](3506)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  O u0 p2 c0 {5,S} {6,S}
3  O u0 p2 c0 {1,S} {6,S}
4  O u1 p2 c0 {5,S}
5  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
6  C u0 p0 c0 {2,S} {3,S} {7,D}
7  C u1 p0 c0 {6,D} {10,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (135.248,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3120,650,792.5,1650,320.801,320.805,320.811,320.816],'cm^-1')),
        HinderedRotor(inertia=(0.00163802,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.367328,'amu*angstrom^2'), symmetry=1, barrier=(26.8268,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.367335,'amu*angstrom^2'), symmetry=1, barrier=(26.8267,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.930584,0.0657781,-7.77982e-05,4.4355e-08,-9.81697e-12,16378.6,22.9164], Tmin=(100,'K'), Tmax=(1109.21,'K')), NASAPolynomial(coeffs=[15.663,0.0126504,-5.9528e-06,1.17383e-09,-8.45227e-14,13110.4,-49.6859], Tmin=(1109.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(135.248,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2sCF) + group(Cs-OsOsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(OCOJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(OF)O[CH]O(3507)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {3,S}
2  O u0 p2 c0 {5,S} {6,S}
3  O u0 p2 c0 {1,S} {5,S}
4  O u0 p2 c0 {6,S} {9,S}
5  C u0 p0 c0 {2,S} {3,S} {7,D}
6  C u1 p0 c0 {2,S} {4,S} {8,S}
7  C u1 p0 c0 {5,D} {10,S}
8  H u0 p0 c0 {6,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (96.8099,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,3615,1277.5,1000,350,440,435,1725,3025,407.5,1350,352.5,3120,650,792.5,1650,373.567,373.567,373.567],'cm^-1')),
        HinderedRotor(inertia=(0.144149,'amu*angstrom^2'), symmetry=1, barrier=(14.275,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144149,'amu*angstrom^2'), symmetry=1, barrier=(14.275,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144149,'amu*angstrom^2'), symmetry=1, barrier=(14.275,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144149,'amu*angstrom^2'), symmetry=1, barrier=(14.275,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (106.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.307116,0.0851313,-0.000134075,9.85142e-08,-2.69599e-11,11773,24.4758], Tmin=(100,'K'), Tmax=(758.266,'K')), NASAPolynomial(coeffs=[16.2042,0.0105068,-4.72239e-06,8.50236e-10,-5.59269e-14,9096.66,-49.57], Tmin=(758.266,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(96.8099,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2sCF) + group(Cs-OsOsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(OCJO) + radical(Cds_P)"""),
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
    label = '[CH2]C(=O)OC=O(3508)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {4,S} {6,S}
2 O u0 p2 c0 {4,D}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {1,S} {2,D} {5,S}
5 C u1 p0 c0 {4,S} {7,S} {8,S}
6 C u0 p0 c0 {1,S} {3,D} {9,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-328.968,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,211.901,1099,1099.01,1099.01,1099.01,4000],'cm^-1')),
        HinderedRotor(inertia=(0.195085,'amu*angstrom^2'), symmetry=1, barrier=(6.2161,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.598819,'amu*angstrom^2'), symmetry=1, barrier=(19.0806,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.195081,'amu*angstrom^2'), symmetry=1, barrier=(6.2161,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (87.054,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.4321,0.0370197,-3.48417e-05,1.80201e-08,-3.96144e-12,-39511.4,22.4016], Tmin=(100,'K'), Tmax=(1058.62,'K')), NASAPolynomial(coeffs=[7.24702,0.0188264,-9.06277e-06,1.78571e-09,-1.27567e-13,-40530.8,-1.10185], Tmin=(1058.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-328.968,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-O2d)) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsOs) + group(Cds-OdOsH) + radical(CJCO)"""),
)

species(
    label = '[O]C=O(135)',
    structure = adjacencyList("""multiplicity 2
1 O u1 p2 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 H u0 p0 c0 {3,S}
"""),
    E0 = (-138.762,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (45.0174,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.61,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51545,0.0112656,-9.58082e-06,4.36317e-09,-8.44578e-13,-16672.2,7.37034], Tmin=(100,'K'), Tmax=(1184.08,'K')), NASAPolynomial(coeffs=[5.09941,0.00591471,-2.80224e-06,5.46639e-10,-3.87736e-14,-17047.4,-0.538996], Tmin=(1184.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-138.762,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""formyloxy""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=[C]OF(3262)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {4,S}
3 C u0 p0 c0 {4,D} {5,S} {6,S}
4 C u1 p0 c0 {2,S} {3,D}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (235.262,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,2950,3100,1380,975,1025,1650,1685,370,1858.57],'cm^-1')),
        HinderedRotor(inertia=(0.190472,'amu*angstrom^2'), symmetry=1, barrier=(8.45471,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (61.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.27376,0.0120477,-4.21846e-06,6.03974e-11,9.31488e-14,28267.1,10.4556], Tmin=(100,'K'), Tmax=(2667.03,'K')), NASAPolynomial(coeffs=[19.9694,-0.00396407,5.52786e-07,-7.38524e-11,6.5216e-15,17217.5,-85.6832], Tmin=(2667.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(235.262,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(C=CJO)"""),
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
    label = '[CH2]C(=O)OF(1338)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {4,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {2,S} {3,D} {5,S}
5 C u1 p0 c0 {4,S} {6,S} {7,S}
6 H u0 p0 c0 {5,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-97.2878,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,3000,3100,440,815,1455,1000,180,180,738.025,738.027,2257.93],'cm^-1')),
        HinderedRotor(inertia=(1.09055,'amu*angstrom^2'), symmetry=1, barrier=(25.0738,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109118,'amu*angstrom^2'), symmetry=1, barrier=(42.1763,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (77.0344,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.81691,0.0179647,1.8901e-05,-5.41702e-08,3.05982e-11,-11699.4,9.60013], Tmin=(10,'K'), Tmax=(683.94,'K')), NASAPolynomial(coeffs=[5.55193,0.020299,-1.35926e-05,4.18552e-09,-4.85861e-13,-12228.7,-0.245396], Tmin=(683.94,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-97.2878,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), label="""[CH2]C(DO)OF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]F(154)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 O u1 p2 c0 {1,S}
"""),
    E0 = (102.686,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(34.9933,'amu')),
        LinearRotor(inertia=(15.6231,'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([1151.69],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (34.9978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51779,-0.00115859,7.7708e-06,-9.55983e-09,3.74809e-12,12350.1,5.42233], Tmin=(10,'K'), Tmax=(823.913,'K')), NASAPolynomial(coeffs=[3.16214,0.00226288,-1.54389e-06,4.73852e-10,-5.4015e-14,12351.2,6.72012], Tmin=(823.913,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(102.686,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""[O]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'C=[C]OC=O(3509)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {3,S} {5,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {6,S}
4 C u0 p0 c0 {5,D} {7,S} {8,S}
5 C u1 p0 c0 {1,S} {4,D}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-54.4552,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,1685,370,180,877.067],'cm^-1')),
        HinderedRotor(inertia=(1.02066,'amu*angstrom^2'), symmetry=1, barrier=(23.467,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.344321,'amu*angstrom^2'), symmetry=1, barrier=(23.4666,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0546,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.43676,0.0295366,-1.86126e-05,3.82258e-09,2.59809e-13,-6489.17,17.8736], Tmin=(100,'K'), Tmax=(1220.15,'K')), NASAPolynomial(coeffs=[9.68764,0.0123157,-5.49363e-06,1.05395e-09,-7.43111e-14,-8746.14,-20.5483], Tmin=(1220.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-54.4552,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-OdOsH) + radical(C=CJO)"""),
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
    label = '[CH]=C(OF)OC=O(3510)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {5,S} {6,S}
3 O u0 p2 c0 {1,S} {5,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {2,S} {3,S} {7,D}
6 C u0 p0 c0 {2,S} {4,D} {8,S}
7 C u1 p0 c0 {5,D} {9,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-46.4338,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,350,440,435,1725,2782.5,750,1395,475,1775,1000,3120,650,792.5,1650,180,1838.11],'cm^-1')),
        HinderedRotor(inertia=(1.24925,'amu*angstrom^2'), symmetry=1, barrier=(28.7226,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.507822,'amu*angstrom^2'), symmetry=1, barrier=(11.6758,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24962,'amu*angstrom^2'), symmetry=1, barrier=(28.7312,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.044,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22161,0.0679686,-0.000115528,1.04078e-07,-3.69101e-11,-5491.28,21.7534], Tmin=(100,'K'), Tmax=(795.679,'K')), NASAPolynomial(coeffs=[8.20716,0.0231271,-1.26622e-05,2.53098e-09,-1.7856e-13,-6295.11,-8.41696], Tmin=(795.679,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-46.4338,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2sCF) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cds-OdOsH) + radical(Cds_P)"""),
)

species(
    label = 'C=C(OF)O[C]=O(3511)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {5,S} {7,S}
3 O u0 p2 c0 {1,S} {5,S}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {2,S} {3,S} {6,D}
6 C u0 p0 c0 {5,D} {8,S} {9,S}
7 C u1 p0 c0 {2,S} {4,D}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-97.0815,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,350,440,435,1725,2950,3100,1380,975,1025,1650,1855,455,950,239.87,244.236,256.018],'cm^-1')),
        HinderedRotor(inertia=(0.611695,'amu*angstrom^2'), symmetry=1, barrier=(25.7081,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115494,'amu*angstrom^2'), symmetry=1, barrier=(4.74644,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.934169,'amu*angstrom^2'), symmetry=1, barrier=(42.5489,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.044,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06456,0.0724427,-0.000127636,1.15455e-07,-4.03026e-11,-11578.1,20.5421], Tmin=(100,'K'), Tmax=(831.978,'K')), NASAPolynomial(coeffs=[8.49494,0.0225349,-1.20829e-05,2.37087e-09,-1.64565e-13,-12323.6,-10.9882], Tmin=(831.978,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-97.0815,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2sCF) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cds-OdOsH) + radical((O)CJOC)"""),
)

species(
    label = 'O=CO(196)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,S} {5,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {1,S}
"""),
    E0 = (-389.483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (46.0254,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.7662,-0.000867303,3.32874e-05,-3.95806e-08,1.39937e-11,-46830.2,7.569], Tmin=(100,'K'), Tmax=(997.633,'K')), NASAPolynomial(coeffs=[5.52477,0.00779102,-3.35072e-06,6.86787e-10,-5.23885e-14,-47962.8,-4.82894], Tmin=(997.633,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-389.483,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""formic_acid""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C#COF(2902)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {3,S}
3 C u0 p0 c0 {2,S} {4,T}
4 C u0 p0 c0 {3,T} {5,S}
5 H u0 p0 c0 {4,S}
"""),
    E0 = (221.918,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (60.0271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.9094,0.00583843,6.25691e-05,-1.62148e-07,1.16725e-10,26695.5,7.9678], Tmin=(10,'K'), Tmax=(505.498,'K')), NASAPolynomial(coeffs=[5.5535,0.011202,-7.86683e-06,2.64783e-09,-3.37635e-13,26294.6,-1.16432], Tmin=(505.498,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(221.918,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""C#COF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[CH]C(OF)OC=O(3512)',
    structure = adjacencyList("""1  F u0 p3 c0 {3,S}
2  O u0 p2 c0 {5,S} {6,S}
3  O u0 p2 c0 {1,S} {5,S}
4  O u0 p2 c0 {6,D}
5  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {4,D} {9,S}
7  C u0 p1 c0 {5,S} {10,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-22.2414,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (106.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.82981,0.053767,-5.97122e-05,3.69996e-08,-9.85519e-12,-2601.99,22.8105], Tmin=(100,'K'), Tmax=(880.959,'K')), NASAPolynomial(coeffs=[7.63206,0.0274227,-1.48575e-05,3.05689e-09,-2.23201e-13,-3624.33,-4.44683], Tmin=(880.959,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-22.2414,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2sCF) + group(Cs-CsOsOsH) + group(Cds-OdOsH) + group(CsJ2_singlet-CsH)"""),
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
    E0 = (-203.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (19.4817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-9.32886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-101.268,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-69.8715,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (140.832,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (79.9657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (44.7408,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (242.09,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (203.652,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-155.555,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (178.37,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (19.0734,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (130.1,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (247.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (196.592,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-40.9405,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (82.0713,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C(OF)OC=O(3338)'],
    products = ['HF(38)', 'CO2(13)', 'CH2CO(75)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(0.00285451,'s^-1'), n=4.62568, Ea=(7.90742,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(12)', 'C=C(O)OF(3492)'],
    products = ['C=C(OF)OC=O(3338)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.27e-07,'m^3/(mol*s)'), n=3.7, Ea=(256.731,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1COCbCdCsCtHNOSSidSis->Cs_2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->H_1COCbCdCtHNOSSidSis->O',), comment="""Estimated from node Root_N-1COCbCdCsCtHNOSSidSis->Cs_2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->H_1COCbCdCtHNOSSidSis->O"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CO2(13)', 'C=COF(3283)'],
    products = ['C=C(OF)OC=O(3338)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(10.2406,'m^3/(mol*s)'), n=1.86833, Ea=(316.938,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO2;R_H] for rate rule [CO2_Od;Cd/H/NonDeO]
Euclidian distance = 3.1622776601683795
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_CO2"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=C(OF)OC=O(3338)'],
    products = ['O=CCC(=O)OF(3343)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(110.393,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=C(OF)OC=O(3338)'],
    products = ['O=COC(=O)CF(3502)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(141.789,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction6',
    reactants = ['FO[C]1CO[CH]O1(3503)'],
    products = ['C=C(OF)OC=O(3338)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(6.43612e+09,'s^-1'), n=0.0758668, Ea=(56.4791,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for R5JJ
Exact match found for rate rule [R5JJ]
Euclidian distance = 0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C(OF)O[C]=O(3504)'],
    products = ['C=C(OF)OC=O(3338)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C[C](OF)O[C]=O(3505)'],
    products = ['C=C(OF)OC=O(3338)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.328e+09,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=C(OF)OC[O](3506)'],
    products = ['C=C(OF)OC=O(3338)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=C(OF)O[CH]O(3507)'],
    products = ['C=C(OF)OC=O(3338)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', '[CH2]C(=O)OC=O(3508)'],
    products = ['C=C(OF)OC=O(3338)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.55564e+07,'m^3/(mol*s)'), n=-5.63145e-08, Ea=(18.6529,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.7121440562946592, var=2.9508506800589083, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C=O(135)', 'C=[C]OF(3262)'],
    products = ['C=C(OF)OC=O(3338)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.47663e+07,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction13',
    reactants = ['HCO(14)', '[CH2]C(=O)OF(1338)'],
    products = ['C=C(OF)OC=O(3338)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(2.0138,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]F(154)', 'C=[C]OC=O(3509)'],
    products = ['C=C(OF)OC=O(3338)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(5)', '[CH]=C(OF)OC=O(3510)'],
    products = ['C=C(OF)OC=O(3338)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(6.6015e+19,'m^3/(mol*s)'), n=-4.65728, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_N-4R!H->C_N-Sp-3C-2C_Ext-3C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_N-3C-inRing_Ext-3C-R_N-4R!H->C_N-Sp-3C-2C_Ext-3C-R
Ea raised from -10.3 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(5)', 'C=C(OF)O[C]=O(3511)'],
    products = ['C=C(OF)OC=O(3338)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.31504e+12,'m^3/(mol*s)'), n=-2.20453, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_N-3BrCClINOPSSi->C_N-3ClOS->S_3ClO->O',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_N-3BrCClINOPSSi->C_N-3ClOS->S_3ClO->O"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=C(OF)OC=O(3338)'],
    products = ['O=CO(196)', 'C#COF(2902)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.32702e+12,'s^-1'), n=0, Ea=(170.72,'kJ/mol'), T0=(1,'K'), Tmin=(500,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R!H->C_N-5R!H->N_Ext-2R!H-R_1BrClFINOPSSi->O_7R!H->O',), comment="""Estimated from node Root_N-1R!H->C_N-5R!H->N_Ext-2R!H-R_1BrClFINOPSSi->O_7R!H->O
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]C(OF)OC=O(3512)'],
    products = ['C=C(OF)OC=O(3338)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.20443e+17,'s^-1'), n=-1.43042, Ea=(22.4435,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.514991623114097, var=15.884634944160005, Tref=1000.0, N=7, data_mean=0.0, correlation='CCH',), comment="""Estimated from node CCH"""),
)

network(
    label = 'PDepNetwork #1299',
    isomers = [
        'C=C(OF)OC=O(3338)',
    ],
    reactants = [
        ('HF(38)', 'CO2(13)', 'CH2CO(75)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1299',
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

