species(
    label = 'O=C(CO)OF(1588)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {5,S} {9,S}
3 O u0 p2 c0 {1,S} {6,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
6 C u0 p0 c0 {3,S} {4,D} {5,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {2,S}
"""),
    E0 = (-454.339,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,990,1113,2750,2850,1437.5,1250,1305,750,350,366.594,366.594,366.594,366.594,366.594,366.594],'cm^-1')),
        HinderedRotor(inertia=(0.31383,'amu*angstrom^2'), symmetry=1, barrier=(29.9291,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0370302,'amu*angstrom^2'), symmetry=1, barrier=(42.3454,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00630388,'amu*angstrom^2'), symmetry=1, barrier=(29.9291,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.0418,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.99043,0.0417459,-3.58036e-05,1.51651e-08,-2.57813e-12,-54570.1,18.7214], Tmin=(100,'K'), Tmax=(1391.45,'K')), NASAPolynomial(coeffs=[11.3021,0.0149783,-6.94819e-06,1.34024e-09,-9.42742e-14,-57161.5,-29.278], Tmin=(1391.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-454.339,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsOs)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.56838,-0.000852123,2.48917e-06,-1.5633e-09,3.13593e-13,-14284.3,3.57912], Tmin=(100,'K'), Tmax=(1571.64,'K')), NASAPolynomial(coeffs=[2.91307,0.00164657,-6.88609e-07,1.21036e-10,-7.84009e-15,-14180.9,6.71041], Tmin=(1571.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-118.741,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""CO""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'OCOF(1675)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {4,S} {7,S}
3 O u0 p2 c0 {1,S} {4,S}
4 C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
5 H u0 p0 c0 {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {2,S}
"""),
    E0 = (-285.2,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,2750,2850,1437.5,1250,1305,750,350,412.755],'cm^-1')),
        HinderedRotor(inertia=(0.059005,'amu*angstrom^2'), symmetry=1, barrier=(7.13393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.229339,'amu*angstrom^2'), symmetry=1, barrier=(27.7274,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (66.0318,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.92989,0.00447927,7.97445e-05,-1.7605e-07,1.16782e-10,-34301.2,9.66319], Tmin=(10,'K'), Tmax=(500.76,'K')), NASAPolynomial(coeffs=[3.33668,0.0227736,-1.5661e-05,5.0238e-09,-6.05838e-13,-34411.7,10.4177], Tmin=(500.76,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-285.2,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), label="""OCOF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OF(173)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (-95.2653,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(36.0011,'amu')),
        NonlinearRotor(inertia=([0.860315,18.4105,19.2708],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([1005.07,1417.2,3730.42],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (36.0058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02848,-0.00192233,1.33972e-05,-1.61993e-08,6.45714e-12,-11458,4.36077], Tmin=(10,'K'), Tmax=(768.674,'K')), NASAPolynomial(coeffs=[3.2293,0.00415811,-2.21832e-06,5.96331e-10,-6.32051e-14,-11391.9,7.63678], Tmin=(768.674,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-95.2653,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""OF""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C=CO(475)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,S} {6,S}
2 O u0 p2 c0 {4,D}
3 C u0 p0 c0 {1,S} {4,D} {5,S}
4 C u0 p0 c0 {2,D} {3,D}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {1,S}
"""),
    E0 = (-165.709,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,2120,512.5,787.5],'cm^-1')),
        HinderedRotor(inertia=(0.80813,'amu*angstrom^2'), symmetry=1, barrier=(18.5805,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (58.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.86771,0.0228707,-2.2339e-05,1.13145e-08,-2.26468e-12,-19887.6,11.2463], Tmin=(100,'K'), Tmax=(1217.46,'K')), NASAPolynomial(coeffs=[7.78054,0.00672919,-2.45109e-06,4.23933e-10,-2.82988e-14,-21083.9,-13.4219], Tmin=(1217.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-165.709,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""hydroxyketene""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=C(OO)OF(1687)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {4,S} {5,S}
3 O u0 p2 c0 {1,S} {5,S}
4 O u0 p2 c0 {2,S} {9,S}
5 C u0 p0 c0 {2,S} {3,S} {6,D}
6 C u0 p0 c0 {5,D} {7,S} {8,S}
7 H u0 p0 c0 {6,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (-41.4791,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,231,791,350,440,435,1725,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.158837,'amu*angstrom^2'), symmetry=1, barrier=(3.65198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158174,'amu*angstrom^2'), symmetry=1, barrier=(3.63674,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.48089,'amu*angstrom^2'), symmetry=1, barrier=(34.0486,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.0418,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56951,0.0618152,-0.00011051,1.04497e-07,-3.77322e-11,-4909.32,20.2372], Tmin=(100,'K'), Tmax=(844.095,'K')), NASAPolynomial(coeffs=[5.72528,0.024455,-1.27241e-05,2.47033e-09,-1.70403e-13,-4981.51,4.62051], Tmin=(844.095,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-41.4791,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2sCF) + group(O2s-OsH) + group(Cds-CdsCsCs) + group(Cds-CdsHH)"""),
)

species(
    label = '[O]C([CH]O)OF(1688)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {5,S}
3 O u0 p2 c0 {6,S} {9,S}
4 O u1 p2 c0 {5,S}
5 C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
6 C u1 p0 c0 {3,S} {5,S} {8,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {3,S}
"""),
    E0 = (-76.1187,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,180,180,910.562],'cm^-1')),
        HinderedRotor(inertia=(0.182158,'amu*angstrom^2'), symmetry=1, barrier=(4.18817,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.186178,'amu*angstrom^2'), symmetry=1, barrier=(4.28059,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.190072,'amu*angstrom^2'), symmetry=1, barrier=(4.37013,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.0418,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34052,0.0687408,-0.000130274,1.24347e-07,-4.45271e-11,-9069.03,23.8062], Tmin=(100,'K'), Tmax=(860.848,'K')), NASAPolynomial(coeffs=[6.01681,0.0241457,-1.27246e-05,2.45687e-09,-1.67859e-13,-9026.88,6.86761], Tmin=(860.848,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-76.1187,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[O]CC([O])OF(1689)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {5,S}
3 O u1 p2 c0 {5,S}
4 O u1 p2 c0 {6,S}
5 C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
6 C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-30.7115,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,180,180,1644.36,1645.56],'cm^-1')),
        HinderedRotor(inertia=(0.171673,'amu*angstrom^2'), symmetry=1, barrier=(3.9471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171957,'amu*angstrom^2'), symmetry=1, barrier=(3.95363,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.0418,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92793,0.0578175,-0.000113584,1.20513e-07,-4.71679e-11,-3630.98,23.0069], Tmin=(100,'K'), Tmax=(844.536,'K')), NASAPolynomial(coeffs=[0.779297,0.0336335,-1.80139e-05,3.53681e-09,-2.45389e-13,-2380.5,34.609], Tmin=(844.536,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-30.7115,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[O]C[C](O)OF(1690)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {6,S}
3 O u0 p2 c0 {6,S} {9,S}
4 O u1 p2 c0 {5,S}
5 C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
6 C u1 p0 c0 {2,S} {3,S} {5,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {3,S}
"""),
    E0 = (-51.1702,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,360,370,350,214.947,2200.28,2200.29],'cm^-1')),
        HinderedRotor(inertia=(0.225722,'amu*angstrom^2'), symmetry=1, barrier=(7.40483,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.941608,'amu*angstrom^2'), symmetry=1, barrier=(30.8878,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.225786,'amu*angstrom^2'), symmetry=1, barrier=(7.40507,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.0418,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62101,0.0610574,-0.000110746,1.06102e-07,-3.86626e-11,-6077.12,23.9375], Tmin=(100,'K'), Tmax=(844.345,'K')), NASAPolynomial(coeffs=[5.27668,0.0249013,-1.30484e-05,2.53992e-09,-1.75402e-13,-6022.96,10.896], Tmin=(844.345,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-51.1702,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + radical(CCOJ) + radical(Cs_P)"""),
)

species(
    label = 'OC=C(O)OF(1691)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {5,S}
3 O u0 p2 c0 {5,S} {8,S}
4 O u0 p2 c0 {6,S} {9,S}
5 C u0 p0 c0 {2,S} {3,S} {6,D}
6 C u0 p0 c0 {4,S} {5,D} {7,S}
7 H u0 p0 c0 {6,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (-374.767,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,3580,3650,1210,1345,900,1100,350,440,435,1725,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(0.981115,'amu*angstrom^2'), symmetry=1, barrier=(22.5578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.979488,'amu*angstrom^2'), symmetry=1, barrier=(22.5204,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.983847,'amu*angstrom^2'), symmetry=1, barrier=(22.6206,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.0418,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.276708,0.068768,-8.88609e-05,5.0848e-08,-1.05578e-11,-44928.6,20.6877], Tmin=(100,'K'), Tmax=(1378.37,'K')), NASAPolynomial(coeffs=[20.8191,-0.00376884,4.14057e-06,-9.35551e-10,6.80913e-14,-49363.9,-80.556], Tmin=(1378.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-374.767,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsOsH)"""),
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
    label = '[O]C(=O)CO(1692)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {4,S} {8,S}
2 O u1 p2 c0 {5,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5 C u0 p0 c0 {2,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {1,S}
"""),
    E0 = (-363.141,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,180,180,1070.72,1073.93,1074.03,1075.02],'cm^-1')),
        HinderedRotor(inertia=(0.00595127,'amu*angstrom^2'), symmetry=1, barrier=(4.89671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.42185,'amu*angstrom^2'), symmetry=1, barrier=(55.6831,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (75.0434,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.0615,0.0238333,-1.12922e-05,1.84024e-09,-4.47566e-14,-43645,14.8778], Tmin=(100,'K'), Tmax=(2108.69,'K')), NASAPolynomial(coeffs=[12.6007,0.0103872,-5.03438e-06,9.07315e-10,-5.81038e-14,-48701.6,-40.711], Tmin=(2108.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-363.141,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsOs) + radical(CCOJ)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51457,2.92734e-05,-5.32151e-07,1.01948e-09,-3.85939e-13,3414.25,2.10435], Tmin=(100,'K'), Tmax=(1145.76,'K')), NASAPolynomial(coeffs=[3.07194,0.00060402,-1.39806e-08,-2.13441e-11,2.48061e-15,3579.39,4.57801], Tmin=(1145.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.3945,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH(D)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH2]C(=O)OF(1366)',
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
        HarmonicOscillator(frequencies=([990,1113,3000,3100,440,815,1455,1000,180,180,180,738.181,738.499],'cm^-1')),
        HinderedRotor(inertia=(0.064981,'amu*angstrom^2'), symmetry=1, barrier=(25.0738,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0116579,'amu*angstrom^2'), symmetry=1, barrier=(42.1763,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (77.0344,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.81691,0.0179647,1.8901e-05,-5.41702e-08,3.05982e-11,-11699.4,9.60013], Tmin=(10,'K'), Tmax=(683.94,'K')), NASAPolynomial(coeffs=[5.55193,0.020299,-1.35926e-05,4.18552e-09,-4.85861e-13,-12228.7,-0.245396], Tmin=(683.94,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-97.2878,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), label="""[CH2]C(DO)OF""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]CC(=O)OF(1693)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {6,S}
3 O u1 p2 c0 {5,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
6 C u0 p0 c0 {2,S} {4,D} {5,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-210.605,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,2750,2850,1437.5,1250,1305,750,350,294.074,294.458,295.723,297.197,1050.01,1050.47,2598.41],'cm^-1')),
        HinderedRotor(inertia=(0.124157,'amu*angstrom^2'), symmetry=1, barrier=(33.778,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.845507,'amu*angstrom^2'), symmetry=1, barrier=(52.0307,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.0338,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.38315,0.0370496,-3.35537e-05,1.57506e-08,-3.07185e-12,-25273.1,18.0922], Tmin=(100,'K'), Tmax=(1194.53,'K')), NASAPolynomial(coeffs=[8.4009,0.0168986,-8.24959e-06,1.6284e-09,-1.16251e-13,-26710.8,-12.0097], Tmin=(1194.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-210.605,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsOs) + radical(C=OCOJ)"""),
)

species(
    label = '[O]F(228)',
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
    label = 'O=[C]CO(1694)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {3,S} {7,S}
2 O u0 p2 c0 {4,D}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4 C u1 p0 c0 {2,D} {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {1,S}
"""),
    E0 = (-177.837,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.20959,'amu*angstrom^2'), symmetry=1, barrier=(27.8108,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21065,'amu*angstrom^2'), symmetry=1, barrier=(27.8352,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (59.044,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.8007,0.0149395,2.57533e-05,-4.93228e-08,2.10358e-11,-21335.3,13.9997], Tmin=(100,'K'), Tmax=(958.142,'K')), NASAPolynomial(coeffs=[12.515,0.00304265,-4.8669e-07,1.51278e-10,-1.80956e-14,-24512.3,-39.3154], Tmin=(958.142,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-177.837,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CsCJ=O)"""),
)

species(
    label = 'O=[C]OF(142)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {2,S}
2 O u0 p2 c0 {1,S} {4,S}
3 O u0 p2 c0 {4,D}
4 C u1 p0 c0 {2,S} {3,D}
"""),
    E0 = (-35.6837,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.77737,'amu*angstrom^2'), symmetry=1, barrier=(40.8652,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (63.0078,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.98175,0.0254561,-4.70439e-05,4.25357e-08,-1.44417e-11,-4258,9.5318], Tmin=(100,'K'), Tmax=(878.493,'K')), NASAPolynomial(coeffs=[5.67742,0.00648683,-3.22263e-06,6.05559e-10,-4.04983e-14,-4473.27,-1.65362], Tmin=(878.493,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-35.6837,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical((O)CJOC)"""),
)

species(
    label = 'CH2OH(33)',
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
    label = 'O=C([CH]O)OF(482)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {5,S} {8,S}
3 O u0 p2 c0 {1,S} {6,S}
4 O u0 p2 c0 {6,D}
5 C u1 p0 c0 {2,S} {6,S} {7,S}
6 C u0 p0 c0 {3,S} {4,D} {5,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {2,S}
"""),
    E0 = (-315.854,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,990,1113,3025,407.5,1350,352.5,335.949,336.192,336.733,989.126,1250.66,1250.94],'cm^-1')),
        HinderedRotor(inertia=(0.617611,'amu*angstrom^2'), symmetry=1, barrier=(49.7597,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.618396,'amu*angstrom^2'), symmetry=1, barrier=(49.7548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.618613,'amu*angstrom^2'), symmetry=1, barrier=(49.7619,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (93.0338,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.17558,0.0422309,-4.34208e-05,2.3242e-08,-5.11561e-12,-37924.5,17.8181], Tmin=(100,'K'), Tmax=(1077.81,'K')), NASAPolynomial(coeffs=[8.91242,0.0172292,-8.62594e-06,1.72017e-09,-1.23628e-13,-39376.7,-15.1881], Tmin=(1077.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-315.854,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsOs) + radical(OCJC=O)"""),
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
    E0 = (-317.265,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-73.596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (25.3233,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (149.812,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (60.9353,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (147.442,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (98.1009,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-77.6366,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-175.497,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (45.8599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (115.952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (39.6022,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (50.351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (14.7273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=C(CO)OF(1588)'],
    products = ['HF(38)', 'CO2(13)', 'CH2O(19)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(0.00285451,'s^-1'), n=4.62568, Ea=(22.3203,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(12)', 'OCOF(1675)'],
    products = ['O=C(CO)OF(1588)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.000742267,'m^3/(mol*s)'), n=2.83796, Ea=(215.592,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.2632766423807829, var=145.9379134136138, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-1COCbCdCsCtHNOSSidSis->Cs',), comment="""Estimated from node Root_N-1COCbCdCsCtHNOSSidSis->Cs"""),
)

reaction(
    label = 'reaction3',
    reactants = ['OF(173)', 'O=C=CO(475)'],
    products = ['O=C(CO)OF(1588)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(51.5,'cm^3/(mol*s)'), n=3.05, Ea=(171.544,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 5 used for cco_HNd;H_OH
Exact match found for rate rule [cco_HNd;H_OH]
Euclidian distance = 0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=C(OO)OF(1687)'],
    products = ['O=C(CO)OF(1588)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(76.5385,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]C([CH]O)OF(1688)'],
    products = ['O=C(CO)OF(1588)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]CC([O])OF(1689)'],
    products = ['O=C(CO)OF(1588)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C[C](O)OF(1690)'],
    products = ['O=C(CO)OF(1588)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7.76e+08,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['OC=C(O)OF(1691)'],
    products = ['O=C(CO)OF(1588)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(182.377,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_Ext-1R!H-R',), comment="""Estimated from node Root_3R!H->C_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F(37)', '[O]C(=O)CO(1692)'],
    products = ['O=C(CO)OF(1588)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.11128e+07,'m^3/(mol*s)'), n=-5.63145e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.7121440562946592, var=2.9508506800589083, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_3BrCClFINPSSi->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction10',
    reactants = ['OH(4)', '[CH2]C(=O)OF(1366)'],
    products = ['O=C(CO)OF(1588)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.7e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R
Ea raised from -2.6 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(5)', '[O]CC(=O)OF(1693)'],
    products = ['O=C(CO)OF(1588)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]F(228)', 'O=[C]CO(1694)'],
    products = ['O=C(CO)OF(1588)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction13',
    reactants = ['O=[C]OF(142)', 'CH2OH(33)'],
    products = ['O=C(CO)OF(1588)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.81979e+07,'m^3/(mol*s)'), n=-0.126529, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_2CF->C"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(5)', 'O=C([CH]O)OF(482)'],
    products = ['O=C(CO)OF(1588)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.88633e+07,'m^3/(mol*s)'), n=0.213913, Ea=(4.02328,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=8.197906922440378e-05, var=0.01210657880115639, Tref=1000.0, N=9, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_Ext-5R!H-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_N-3R!H->F_3BrCClINOPSSi->C_Ext-2C-R_N-3C-inRing_Ext-3C-R_Ext-5R!H-R"""),
)

network(
    label = 'PDepNetwork #661',
    isomers = [
        'O=C(CO)OF(1588)',
    ],
    reactants = [
        ('HF(38)', 'CO2(13)', 'CH2O(19)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #661',
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

