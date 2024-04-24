species(
    label = 'O=C(O)C(F)C(F)(F)F(1188)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {9,S} {11,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
9  C u0 p0 c0 {5,S} {6,D} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-1277.6,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,206,363,518,1175,1317,1481,3059,193,295,551,588,656,1146,1192,1350,444.331,444.35,444.429,444.457,444.483,444.509],'cm^-1')),
        HinderedRotor(inertia=(0.0357326,'amu*angstrom^2'), symmetry=1, barrier=(5.00687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0135376,'amu*angstrom^2'), symmetry=1, barrier=(41.7924,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.298187,'amu*angstrom^2'), symmetry=1, barrier=(41.7927,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.971163,0.0713453,-8.62449e-05,5.40962e-08,-1.37493e-11,-153555,23.8045], Tmin=(100,'K'), Tmax=(949.066,'K')), NASAPolynomial(coeffs=[11.9721,0.0249798,-1.29636e-05,2.61982e-09,-1.89482e-13,-155643,-28.6934], Tmin=(949.066,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1277.6,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)3-R-CO) + group(Cds-OdCsOs)"""),
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
    label = 'CHFCF2(54)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 C u0 p0 c0 {1,S} {5,D} {6,S}
5 C u0 p0 c0 {2,S} {3,S} {4,D}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-505.825,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(82.003,'amu')),
        NonlinearRotor(inertia=([47.283,130.848,178.131],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([226.38,309.801,488.171,585.738,628.102,776.692,950.625,1195.27,1295.73,1386.74,1841.68,3239.84],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0245,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2091.09,'J/mol'), sigma=(4.442,'angstroms'), dipoleMoment=(1.4,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35521,0.0323113,-3.57904e-05,2.09758e-08,-5.1422e-12,-60645.8,19.2938], Tmin=(298,'K'), Tmax=(1100,'K')), NASAPolynomial(coeffs=[9.29782,0.00669066,-2.71841e-06,4.9759e-10,-3.35867e-14,-62705.1,-20.9391], Tmin=(1100,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-505.825,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(133.032,'J/mol/K'), label="""CHFCF2""", comment="""Thermo library: Fluorine"""),
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
    label = 'OC(F)C(F)(F)F(1271)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {7,S}
4 F u0 p3 c0 {7,S}
5 O u0 p2 c0 {6,S} {9,S}
6 C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
7 C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-1118.65,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,261,493,600,1152,1365,1422,3097,193,295,551,588,656,1146,1192,1350,240.25],'cm^-1')),
        HinderedRotor(inertia=(0.423807,'amu*angstrom^2'), symmetry=1, barrier=(17.4357,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.423252,'amu*angstrom^2'), symmetry=1, barrier=(17.4435,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (118.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.79637,0.0126887,0.000144956,-3.5267e-07,2.39248e-10,-134530,11.5974], Tmin=(10,'K'), Tmax=(531.953,'K')), NASAPolynomial(coeffs=[7.25506,0.0293464,-2.23224e-05,7.74543e-09,-9.93381e-13,-135502,-8.5796], Tmin=(531.953,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1118.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""OC(F)C(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CHF3(41)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 F u0 p3 c0 {4,S}
4 C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
5 H u0 p0 c0 {4,S}
"""),
    E0 = (-708.855,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(70.003,'amu')),
        NonlinearRotor(inertia=([48.8559,48.8561,89.1073],'amu*angstrom^2'), symmetry=3),
        HarmonicOscillator(frequencies=([506.193,506.206,703.466,1151.07,1179.84,1179.87,1414.02,1414.02,3086.64],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (70.0138,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2178.39,'J/mol'), sigma=(4.123,'angstroms'), dipoleMoment=(1.8,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.680758,0.0239859,-2.18336e-05,9.98026e-09,-1.85179e-12,-84989.4,21.0617], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[7.52417,0.00523005,-2.036e-06,3.57296e-10,-2.30997e-14,-87022.5,-14.6112], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-708.855,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(108.088,'J/mol/K'), label="""CHF3""", comment="""Thermo library: Fluorine"""),
)

species(
    label = 'O=C(O)[C]F(826)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,S} {6,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {2,S} {3,D} {5,S}
5 C u0 p1 c0 {1,S} {4,S}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (-238.142,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,262,1290,469.037,469.042,469.058,469.058,469.059],'cm^-1')),
        HinderedRotor(inertia=(0.000766252,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000766216,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (76.0265,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.93204,0.0255218,-2.44778e-05,1.25723e-08,-2.73821e-12,-28605.1,14.1558], Tmin=(100,'K'), Tmax=(1066.47,'K')), NASAPolynomial(coeffs=[6.35429,0.0126859,-6.42393e-06,1.28648e-09,-9.26034e-14,-29335,-2.57485], Tmin=(1066.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-238.142,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-OdCsOs) + group(CJ2_singlet-FCO)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(897.963,'J/mol'), sigma=(3.977,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28591,0.0107608,-1.05382e-05,4.89881e-09,-8.86384e-13,-23521.2,13.1348], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[5.33121,0.00197748,-9.60248e-07,2.10704e-10,-1.5954e-14,-24371.4,-2.56367], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-196.899,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CF2""", comment="""Thermo library: Fluorine"""),
)

species(
    label = 'O=C(O)C(F)F(222)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {6,S} {8,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6 C u0 p0 c0 {3,S} {4,D} {5,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (-811.366,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,195,270,1147,1130,1359,1388,1409,3075,180,180,180,948.084,948.124],'cm^-1')),
        HinderedRotor(inertia=(0.00779228,'amu*angstrom^2'), symmetry=1, barrier=(4.97101,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.066269,'amu*angstrom^2'), symmetry=1, barrier=(42.2671,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0328,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.00261,0.0220897,-2.96049e-06,-9.04426e-09,3.91266e-12,-97592.5,10.9495], Tmin=(10,'K'), Tmax=(1184.8,'K')), NASAPolynomial(coeffs=[10.5537,0.0121968,-5.91178e-06,1.32461e-09,-1.12746e-13,-100003,-25.3872], Tmin=(1184.8,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-811.366,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""ODC(O)C(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FCC(F)(F)F(255)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {6,S}
4 F u0 p3 c0 {6,S}
5 C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
6 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-912.682,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([528,1116,1182,1331,1402,1494,3075,3110,193,295,551,588,656,1146,1192,1350,263.317],'cm^-1')),
        HinderedRotor(inertia=(0.465748,'amu*angstrom^2'), symmetry=1, barrier=(22.9164,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (102.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.219524,0.048987,-5.35013e-05,2.97849e-08,-6.73935e-12,-109444,26.8259], Tmin=(298,'K'), Tmax=(1150,'K')), NASAPolynomial(coeffs=[11.789,0.0101279,-4.32428e-06,8.25437e-10,-5.76051e-14,-112514,-33.8751], Tmin=(1150,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-912.68,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(178.761,'J/mol/K'), label="""CH2F-CF3""", comment="""Thermo library: Fluorine"""),
)

species(
    label = 'OF(155)',
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
    label = 'O=C=CC(F)(F)F(1290)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
6 C u0 p0 c0 {5,S} {7,D} {8,S}
7 C u0 p0 c0 {4,D} {6,D}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-730.861,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([219,296,586,564,718,793,1177,1228,3010,987.5,1337.5,450,1655,2120,512.5,787.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.324383,'amu*angstrom^2'), symmetry=1, barrier=(7.45821,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.77415,0.0175058,9.32954e-05,-2.77798e-07,2.20547e-10,-87900.8,11.5723], Tmin=(10,'K'), Tmax=(455.873,'K')), NASAPolynomial(coeffs=[5.52564,0.0269493,-1.94176e-05,6.4241e-09,-7.93356e-13,-88318.3,1.67048], Tmin=(455.873,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-730.861,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""ODCDCC(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'H2O(3)',
    structure = adjacencyList("""1 O u0 p2 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (-251.755,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1601.58,3620.23,4000],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (18.0153,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4759.22,'J/mol'), sigma=(2.605,'angstroms'), dipoleMoment=(1.844,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=4.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.05764,-0.000787936,2.90877e-06,-1.47519e-09,2.12843e-13,-30281.6,-0.311364], Tmin=(100,'K'), Tmax=(1130.24,'K')), NASAPolynomial(coeffs=[2.84325,0.00275109,-7.81031e-07,1.07244e-10,-5.79392e-15,-29958.6,5.91042], Tmin=(1130.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-251.755,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""H2O""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'CF3CFCO(118)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {6,S}
4 F u0 p3 c0 {7,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
7 C u0 p0 c0 {4,S} {6,S} {8,D}
8 C u0 p0 c0 {5,D} {7,D}
"""),
    E0 = (-842.994,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([219,296,586,564,718,793,1177,1228,145,326,398,834,1303,2120,512.5,787.5,1180.72],'cm^-1')),
        HinderedRotor(inertia=(0.404546,'amu*angstrom^2'), symmetry=1, barrier=(9.3013,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (128.025,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.6191,0.0319905,4.01861e-05,-1.82953e-07,1.52089e-10,-101386,12.436], Tmin=(10,'K'), Tmax=(486.049,'K')), NASAPolynomial(coeffs=[7.46197,0.0258289,-1.93826e-05,6.53779e-09,-8.15302e-13,-102060,-6.42535], Tmin=(486.049,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-842.994,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""ODCDC(F)C(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OC(F)(F)F(254)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 O u0 p2 c0 {5,S} {6,S}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {4,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-924.469,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,431,527,613,668,835,1199,1245,1322],'cm^-1')),
        HinderedRotor(inertia=(0.246742,'amu*angstrom^2'), symmetry=1, barrier=(8.25205,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (86.0132,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3135.61,'J/mol'), sigma=(5.37414,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=489.78 K, Pc=45.84 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93145,0.00398703,6.6773e-05,-1.35516e-07,8.01823e-11,-111184,9.48547], Tmin=(10,'K'), Tmax=(580.755,'K')), NASAPolynomial(coeffs=[4.01947,0.0193555,-1.41814e-05,4.77794e-09,-5.99414e-13,-111464,6.78935], Tmin=(580.755,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-924.469,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""OC(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C=CF(789)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {4,D}
3 C u0 p0 c0 {1,S} {4,D} {5,S}
4 C u0 p0 c0 {2,D} {3,D}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (-160.19,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(60.0011,'amu')),
        NonlinearRotor(inertia=([9.01649,110.348,119.365],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([238.894,460.044,539.472,686.431,1048.64,1208.61,1432.92,2235.18,3235.83],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (60.0271,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2914.22,'J/mol'), sigma=(4.906,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96237,0.0225075,-2.56389e-05,1.67917e-08,-4.73482e-12,-19104,15.6743], Tmin=(298,'K'), Tmax=(1050,'K')), NASAPolynomial(coeffs=[6.77764,0.00623359,-2.55958e-06,4.75314e-10,-3.25951e-14,-20336.9,-8.59113], Tmin=(1050,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-160.19,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(108.088,'J/mol/K'), label="""CHFCO""", comment="""Thermo library: Fluorine"""),
)

species(
    label = 'OC(=CC(F)(F)F)OF(1291)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {6,S}
5  O u0 p2 c0 {9,S} {11,S}
6  O u0 p2 c0 {4,S} {9,S}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
8  C u0 p0 c0 {7,S} {9,D} {10,S}
9  C u0 p0 c0 {5,S} {6,S} {8,D}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-844.651,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,231,791,219,296,586,564,718,793,1177,1228,3010,987.5,1337.5,450,1655,350,440,435,1725,307.726,308.254],'cm^-1')),
        HinderedRotor(inertia=(0.231664,'amu*angstrom^2'), symmetry=1, barrier=(15.4775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127734,'amu*angstrom^2'), symmetry=1, barrier=(8.54632,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128286,'amu*angstrom^2'), symmetry=1, barrier=(8.54772,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.221541,0.0882637,-0.000135182,1.03628e-07,-3.10924e-11,-101457,25.2171], Tmin=(100,'K'), Tmax=(821.207,'K')), NASAPolynomial(coeffs=[14.3069,0.0196631,-9.89139e-06,1.92562e-09,-1.34623e-13,-103770,-39.9636], Tmin=(821.207,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-844.651,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2sCF) + group(CsCdFFF) + group(Cds-CdsCsH) + group(Cds-CdsCsCs)"""),
)

species(
    label = 'OC(=CF)OC(F)(F)F(1292)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {8,S}
6  O u0 p2 c0 {8,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
8  C u0 p0 c0 {5,S} {6,S} {9,D}
9  C u0 p0 c0 {4,S} {8,D} {10,S}
10 H u0 p0 c0 {9,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-1150.22,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,431,527,613,668,835,1199,1245,1322,350,440,435,1725,194,682,905,1196,1383,3221,265.336,265.336,265.336],'cm^-1')),
        HinderedRotor(inertia=(0.291434,'amu*angstrom^2'), symmetry=1, barrier=(14.56,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.291433,'amu*angstrom^2'), symmetry=1, barrier=(14.56,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.291434,'amu*angstrom^2'), symmetry=1, barrier=(14.56,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.301824,0.0814292,-0.000105831,6.61391e-08,-1.59865e-11,-138207,24.5728], Tmin=(100,'K'), Tmax=(1019.84,'K')), NASAPolynomial(coeffs=[17.4929,0.0140017,-6.6554e-06,1.30749e-09,-9.36584e-14,-141713,-58.7015], Tmin=(1019.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1150.22,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(CsFFFO) + group(Cds-CdsCsCs) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs)"""),
)

species(
    label = '[O]C([O])C(F)C(F)(F)F(1293)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u1 p2 c0 {9,S}
6  O u1 p2 c0 {9,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
9  C u0 p0 c0 {5,S} {6,S} {7,S} {11,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-870.927,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([250,417,511,1155,1315,1456,3119,193,295,551,588,656,1146,1192,1350,1380,1390,370,380,2900,435,180,180,1736.15,1736.61],'cm^-1')),
        HinderedRotor(inertia=(0.172558,'amu*angstrom^2'), symmetry=1, barrier=(3.96744,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172277,'amu*angstrom^2'), symmetry=1, barrier=(3.96098,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (146.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.680024,0.0865318,-0.000161403,1.59278e-07,-5.96491e-11,-104642,28.3505], Tmin=(100,'K'), Tmax=(832.579,'K')), NASAPolynomial(coeffs=[4.66398,0.0379163,-2.07122e-05,4.10242e-09,-2.86375e-13,-104284,15.9953], Tmin=(832.579,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-870.927,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCsCsFH) + group(Cs-CsOsOsH) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[O]C(O)[C](F)C(F)(F)F(1294)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {11,S}
6  O u1 p2 c0 {7,S}
7  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
8  C u0 p0 c0 {1,S} {2,S} {3,S} {9,S}
9  C u1 p0 c0 {4,S} {7,S} {8,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-906.158,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,253,522,585,610,849,1160,1215,1402,212,367,445,1450,180,180,2417.56],'cm^-1')),
        HinderedRotor(inertia=(0.211239,'amu*angstrom^2'), symmetry=1, barrier=(4.85679,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.2111,'amu*angstrom^2'), symmetry=1, barrier=(4.8536,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.211148,'amu*angstrom^2'), symmetry=1, barrier=(4.85471,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (146.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.393498,0.0929348,-0.000175913,1.68652e-07,-6.10096e-11,-108869,31.2231], Tmin=(100,'K'), Tmax=(847.865,'K')), NASAPolynomial(coeffs=[6.82752,0.0331323,-1.80142e-05,3.53367e-09,-2.44216e-13,-108901,7.48655], Tmin=(847.865,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-906.158,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCsCsFH) + group(Cs-CsOsOsH) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)) + radical(CCOJ) + radical(CsCsCsF1s)"""),
)

species(
    label = 'OC(O)=C(F)C(F)(F)F(1295)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {9,S} {11,S}
6  O u0 p2 c0 {9,S} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
8  C u0 p0 c0 {4,S} {7,S} {9,D}
9  C u0 p0 c0 {5,S} {6,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-1157.82,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,219,296,586,564,718,793,1177,1228,323,467,575,827,1418,350,440,435,1725,180],'cm^-1')),
        HinderedRotor(inertia=(0.599805,'amu*angstrom^2'), symmetry=1, barrier=(13.7907,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.600683,'amu*angstrom^2'), symmetry=1, barrier=(13.8109,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.600003,'amu*angstrom^2'), symmetry=1, barrier=(13.7953,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.04,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0308202,0.0902056,-0.000131785,9.27975e-08,-2.52312e-11,-139113,24.4582], Tmin=(100,'K'), Tmax=(909.425,'K')), NASAPolynomial(coeffs=[17.2204,0.0145956,-7.0681e-06,1.36753e-09,-9.59391e-14,-142240,-56.8386], Tmin=(909.425,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1157.82,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(CsCdFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cds(F)) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsCsCs)"""),
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
    label = 'O=C(O)[CH]C(F)(F)F(1296)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  O u0 p2 c0 {8,S} {10,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
7  C u1 p0 c0 {6,S} {8,S} {9,S}
8  C u0 p0 c0 {4,S} {5,D} {7,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (-919.799,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,253,522,585,610,849,1160,1215,1402,3025,407.5,1350,352.5,436.068,436.072,436.073,436.073,436.073,436.074],'cm^-1')),
        HinderedRotor(inertia=(0.000886516,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000886553,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000886489,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55396,0.0581248,-6.84028e-05,4.32856e-08,-1.12801e-11,-110542,23.7721], Tmin=(100,'K'), Tmax=(920.312,'K')), NASAPolynomial(coeffs=[9.52159,0.023494,-1.19577e-05,2.39647e-09,-1.72492e-13,-112008,-14.0053], Tmin=(920.312,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-919.799,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-R-CO) + group(Cds-OdCsOs) + radical(CCJCO)"""),
)

species(
    label = 'O=C(O)C(F)[C](F)F(1297)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {7,S} {10,S}
5  O u0 p2 c0 {7,D}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {4,S} {5,D} {6,S}
8  C u1 p0 c0 {2,S} {3,S} {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (-847.02,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,180,233,1109,1254,1325,1339,3290,190,488,555,1236,1407,408.271,408.272,408.277,408.278,408.283,1201.55],'cm^-1')),
        HinderedRotor(inertia=(0.398688,'amu*angstrom^2'), symmetry=1, barrier=(47.1601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.398683,'amu*angstrom^2'), symmetry=1, barrier=(47.1601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0646655,'amu*angstrom^2'), symmetry=1, barrier=(7.64907,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.042,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14462,0.0684097,-9.45845e-05,6.99716e-08,-2.10195e-11,-101775,23.904], Tmin=(100,'K'), Tmax=(808.526,'K')), NASAPolynomial(coeffs=[10.0495,0.0243572,-1.28611e-05,2.59042e-09,-1.86002e-13,-103215,-17.1646], Tmin=(808.526,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-847.02,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-OdCsOs) + radical(Csj(Cs-F1sCOH)(F1s)(F1s))"""),
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
    label = 'O=[C]C(F)C(F)(F)F(1298)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {7,S}
4 F u0 p3 c0 {7,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7 C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
8 C u1 p0 c0 {5,D} {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-862.493,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([255,365,576,1163,1267,1448,3060,193,295,551,588,656,1146,1192,1350,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.238312,'amu*angstrom^2'), symmetry=1, barrier=(5.47926,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.886609,'amu*angstrom^2'), symmetry=1, barrier=(20.3849,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (129.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.61376,0.0321622,5.72137e-05,-2.2209e-07,1.80281e-10,-103731,13.6163], Tmin=(10,'K'), Tmax=(482.045,'K')), NASAPolynomial(coeffs=[7.2644,0.0297382,-2.19644e-05,7.34826e-09,-9.12111e-13,-104406,-4.69047], Tmin=(482.045,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-862.493,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""OD[C]C(F)C(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]C(=O)C(F)C(F)(F)F(1299)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u1 p2 c0 {9,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
9  C u0 p0 c0 {5,S} {6,D} {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-1051.9,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([206,363,518,1175,1317,1481,3059,193,295,551,588,656,1146,1192,1350,505.278,505.307,505.311,505.311,505.315,505.329,505.332],'cm^-1')),
        HinderedRotor(inertia=(0.000660231,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000660233,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (145.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6697,0.0619935,-5.67119e-05,-2.23985e-08,5.63238e-11,-126441,21.6009], Tmin=(100,'K'), Tmax=(475.105,'K')), NASAPolynomial(coeffs=[7.54864,0.030938,-1.68835e-05,3.40869e-09,-2.43444e-13,-127208,-4.57667], Tmin=(475.105,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1051.9,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)3-R-CO) + group(Cds-OdCsOs) + radical(CCOJ)"""),
)

species(
    label = 'CF3(44)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 F u0 p3 c0 {4,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
"""),
    E0 = (-479.453,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(68.9952,'amu')),
        NonlinearRotor(inertia=([46.5949,46.5949,90.0934],'amu*angstrom^2'), symmetry=3),
        HarmonicOscillator(frequencies=([504.804,504.824,700.362,1092.66,1279.94,1279.95],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0058,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1006.05,'J/mol'), sigma=(4.32,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35729,0.0215758,-2.3902e-05,1.29611e-08,-2.80173e-12,-57391.3,18.6474], Tmin=(298,'K'), Tmax=(1250,'K')), NASAPolynomial(coeffs=[7.73167,0.00227115,-8.90496e-07,1.52938e-10,-9.50838e-15,-59145.7,-14.0202], Tmin=(1250,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-479.453,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(83.1447,'J/mol/K'), label="""CF3""", comment="""Thermo library: Fluorine"""),
)

species(
    label = 'O=C(O)[CH]F(1300)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,S} {7,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {2,S} {3,D} {5,S}
5 C u1 p0 c0 {1,S} {4,S} {6,S}
6 H u0 p0 c0 {5,S}
7 H u0 p0 c0 {2,S}
"""),
    E0 = (-437.6,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,235,1215,1347,1486,3221,452.319,452.415,453.012,453.411,453.532],'cm^-1')),
        HinderedRotor(inertia=(0.324557,'amu*angstrom^2'), symmetry=1, barrier=(47.3219,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000822581,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (77.0344,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.91283,0.00519456,8.72792e-05,-1.82058e-07,1.113e-10,-52625.5,10.0943], Tmin=(10,'K'), Tmax=(558.971,'K')), NASAPolynomial(coeffs=[3.89109,0.024878,-1.79443e-05,5.9387e-09,-7.34256e-13,-52928.2,7.45768], Tmin=(558.971,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-437.6,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""ODC(O)[CH]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=[C]O(191)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {3,S} {4,S}
2 O u0 p2 c0 {3,D}
3 C u1 p0 c0 {1,S} {2,D}
4 H u0 p0 c0 {1,S}
"""),
    E0 = (-192.093,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,1244.76,1684.94],'cm^-1')),
        HinderedRotor(inertia=(0.00103939,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (45.0174,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.92208,0.00762454,3.29884e-06,-1.07135e-08,5.11587e-12,-23028.2,11.2926], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.39206,0.00411221,-1.48195e-06,2.39875e-10,-1.43903e-14,-23860.7,-2.23529], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-192.093,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(78.9875,'J/(mol*K)'), label="""HOCO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = 'CF3-CHF(94)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 F u0 p3 c0 {6,S}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
6 C u1 p0 c0 {4,S} {5,S} {7,S}
7 H u0 p0 c0 {6,S}
"""),
    E0 = (-697.945,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,522,585,610,849,1160,1215,1402,334,575,1197,1424,3202,516.186],'cm^-1')),
        HinderedRotor(inertia=(0.380754,'amu*angstrom^2'), symmetry=1, barrier=(8.75429,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.023,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.17776,0.0376992,-4.16073e-05,2.31106e-08,-5.13915e-12,-83881.2,17.2676], Tmin=(298,'K'), Tmax=(1250,'K')), NASAPolynomial(coeffs=[11.7531,0.00713118,-3.12251e-06,6.15301e-10,-4.44146e-14,-86403.9,-31.3333], Tmin=(1250,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-697.945,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(153.818,'J/mol/K'), label="""CF3-CHF""", comment="""Thermo library: Fluorine"""),
)

species(
    label = 'O=C(O)[C](F)C(F)(F)F(1301)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {9,S} {10,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
8  C u1 p0 c0 {4,S} {7,S} {9,S}
9  C u0 p0 c0 {5,S} {6,D} {8,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-1121.23,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,253,522,585,610,849,1160,1215,1402,234,347,1316,1464,391.601,391.602,391.602,391.602,391.602,3074.35],'cm^-1')),
        HinderedRotor(inertia=(0.484184,'amu*angstrom^2'), symmetry=1, barrier=(52.6897,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.484184,'amu*angstrom^2'), symmetry=1, barrier=(52.6897,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.289624,'amu*angstrom^2'), symmetry=1, barrier=(31.5174,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (145.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03656,0.0686013,-8.47026e-05,5.29355e-08,-1.32121e-11,-134749,24.065], Tmin=(100,'K'), Tmax=(972.668,'K')), NASAPolynomial(coeffs=[12.7249,0.0205342,-1.05763e-05,2.12945e-09,-1.53766e-13,-137023,-32.0008], Tmin=(972.668,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1121.23,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)3-R-CO) + group(Cds-OdCsOs) + radical(CsCOCsF1s)"""),
)

species(
    label = 'F2(108)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 F u0 p3 c0 {1,S}
"""),
    E0 = (-8.80492,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(37.9968,'amu')),
        LinearRotor(inertia=(18.2987,'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([1076.6],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (37.9968,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.66955,0.00529418,-6.47091e-06,3.91833e-09,-9.38451e-13,-980.801,7.82659], Tmin=(298,'K'), Tmax=(1150,'K')), NASAPolynomial(coeffs=[4.05774,0.000600034,-2.19218e-07,4.31508e-11,-3.12588e-15,-1324.39,0.863214], Tmin=(1150,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-8.80492,'kJ/mol'), Cp0=(29.1007,'J/mol/K'), CpInf=(37.4151,'J/mol/K'), label="""F2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'O=C(O)C=C(F)F(947)',
    structure = adjacencyList("""1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {6,S} {9,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {6,S} {7,D} {8,S}
6 C u0 p0 c0 {3,S} {4,D} {5,S}
7 C u0 p0 c0 {1,S} {2,S} {5,D}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {3,S}
"""),
    E0 = (-680.412,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,4000,4000,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.000791977,'amu*angstrom^2'), symmetry=1, barrier=(8.99215,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.39274,'amu*angstrom^2'), symmetry=1, barrier=(9.02986,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.35334,0.0410013,-6.93338e-05,6.5514e-08,-2.44534e-11,-81780,17.9381], Tmin=(100,'K'), Tmax=(783.735,'K')), NASAPolynomial(coeffs=[5.45567,0.0174191,-9.3693e-06,1.89146e-09,-1.34737e-13,-82028.3,5.24533], Tmin=(783.735,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-680.412,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + missing(Cd-COCdH) + group(Cds-O2d(Cds-Cds)O2s) + group(CdCFF)"""),
)

species(
    label = 'O=C(O)C(F)=C(F)F(445)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {8,S}
3 F u0 p3 c0 {8,S}
4 O u0 p2 c0 {7,S} {9,S}
5 O u0 p2 c0 {7,D}
6 C u0 p0 c0 {1,S} {7,S} {8,D}
7 C u0 p0 c0 {4,S} {5,D} {6,S}
8 C u0 p0 c0 {2,S} {3,S} {6,D}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (-804.661,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,288,410,724,839,1320,182,240,577,636,1210,1413,180,180,1712.27,1714.57,1715.17],'cm^-1')),
        HinderedRotor(inertia=(0.213533,'amu*angstrom^2'), symmetry=1, barrier=(4.90954,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.212628,'amu*angstrom^2'), symmetry=1, barrier=(4.88874,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4364,0.063091,-0.000106068,9.60819e-08,-3.40779e-11,-96692.5,23.2772], Tmin=(100,'K'), Tmax=(816.633,'K')), NASAPolynomial(coeffs=[7.20996,0.0234135,-1.22525e-05,2.40411e-09,-1.67712e-13,-97255.5,-1.08048], Tmin=(816.633,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-804.661,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CdCCF) + longDistanceInteraction_noncyclic(Cd(F)-CO) + group(Cds-O2d(Cds-Cds)O2s) + group(CdCFF) + longDistanceInteraction_noncyclic(Cds(F)2=Cds(F))"""),
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
    E0 = (-466.325,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-457.262,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-332.926,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-199.323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-434.079,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-95.0087,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-366.142,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-275.928,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-149.913,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-486.521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-288.492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-324.284,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-432.896,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-283.056,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-214.555,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-274.524,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-280.517,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-357.479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-330.464,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-344.625,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-129.643,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-408.42,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=C(O)C(F)C(F)(F)F(1188)'],
    products = ['HF(38)', 'CO2(13)', 'CHFCF2(54)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(3.92484e+07,'s^-1'), n=1.83354, Ea=(251.702,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.49335582114639515, var=20.078174816455903, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_1R!H->C_N-5Br1sCl1sF1sH->H_5Br1sCl1sF1s->F1s_Ext-2C-R_7R!H->O',), comment="""Estimated from node Root_1R!H->C_N-5Br1sCl1sF1sH->H_5Br1sCl1sF1s->F1s_Ext-2C-R_7R!H->O
Multiplied by reaction path degeneracy 3.0"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(12)', 'OC(F)C(F)(F)F(1271)'],
    products = ['O=C(O)C(F)C(F)(F)F(1188)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.000742267,'m^3/(mol*s)'), n=2.83796, Ea=(220.557,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.2632766423807829, var=145.9379134136138, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-1COCbCdCsCtHNOSSidSis->Cs',), comment="""Estimated from node Root_N-1COCbCdCsCtHNOSSidSis->Cs"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CHF3(41)', 'O=C(O)[C]F(826)'],
    products = ['O=C(O)C(F)C(F)(F)F(1188)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.33872e-07,'m^3/(mol*s)'), n=3.38172, Ea=(54.4967,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_N-3Br1sCCl1sF1sHI1s->F1s_N-3Br1sCCl1sH->Cl1s_2Br1sCl1sF1sHI1s->F1s',), comment="""Estimated from node CH_N-3Br1sCCl1sF1sHI1s->F1s_N-3Br1sCCl1sH->Cl1s_2Br1sCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CF2(89)', 'O=C(O)C(F)F(222)'],
    products = ['O=C(O)C(F)C(F)(F)F(1188)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.67164e-06,'m^3/(mol*s)'), n=3.3552, Ea=(249.367,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CY_N-2Br1sCl1sF1sHI1s->H_Ext-4Cs-R',), comment="""Estimated from node CY_N-2Br1sCl1sF1sHI1s->H_Ext-4Cs-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CO2(13)', 'FCC(F)(F)F(255)'],
    products = ['O=C(O)C(F)C(F)(F)F(1188)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(424000,'cm^3/(mol*s)'), n=2.13, Ea=(322.168,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [CO2_Cdd;C_sec]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: 1,3_Insertion_CO2"""),
)

reaction(
    label = 'reaction6',
    reactants = ['OF(155)', 'O=C=CC(F)(F)F(1290)'],
    products = ['O=C(O)C(F)C(F)(F)F(1188)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(51.5,'cm^3/(mol*s)'), n=3.05, Ea=(171.544,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 5 used for cco_HNd;H_OH
Exact match found for rate rule [cco_HNd;H_OH]
Euclidian distance = 0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H2O(3)', 'CF3CFCO(118)'],
    products = ['O=C(O)C(F)C(F)(F)F(1188)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.000454348,'m^3/(mol*s)'), n=2.96333, Ea=(169.034,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_Cdd;H_OH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction8',
    reactants = ['OC(F)(F)F(254)', 'O=C=CF(789)'],
    products = ['O=C(O)C(F)C(F)(F)F(1188)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(6.37684e-08,'m^3/(mol*s)'), n=3.46667, Ea=(249.157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [doublebond;R_OH] + [Cd_Cdd;R_OR] for rate rule [Cd_Cdd;R_OH]
Euclidian distance = 1.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction9',
    reactants = ['OC(=CC(F)(F)F)OF(1291)'],
    products = ['O=C(O)C(F)C(F)(F)F(1188)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(135.165,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction10',
    reactants = ['OC(=CF)OC(F)(F)F(1292)'],
    products = ['O=C(O)C(F)C(F)(F)F(1188)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(104.13,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]C([O])C(F)C(F)(F)F(1293)'],
    products = ['O=C(O)C(F)C(F)(F)F(1188)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C(O)[C](F)C(F)(F)F(1294)'],
    products = ['O=C(O)C(F)C(F)(F)F(1188)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['OC(O)=C(F)C(F)(F)F(1295)'],
    products = ['O=C(O)C(F)C(F)(F)F(1188)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(410000,'s^-1'), n=2.37, Ea=(165.351,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_Ext-1R!H-R',), comment="""Estimated from node Root_3R!H->C_Ext-1R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F(37)', 'O=C(O)[CH]C(F)(F)F(1296)'],
    products = ['O=C(O)C(F)C(F)(F)F(1188)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(5e+07,'m^3/(mol*s)'), n=0, Ea=(4.27748,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_N-2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_N-2CF->C"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F(37)', 'O=C(O)C(F)[C](F)F(1297)'],
    products = ['O=C(O)C(F)C(F)(F)F(1188)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C"""),
)

reaction(
    label = 'reaction16',
    reactants = ['OH(4)', 'O=[C]C(F)C(F)(F)F(1298)'],
    products = ['O=C(O)C(F)C(F)(F)F(1188)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.7e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(5)', '[O]C(=O)C(F)C(F)(F)F(1299)'],
    products = ['O=C(O)C(F)C(F)(F)F(1188)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.42074e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CF3(44)', 'O=C(O)[CH]F(1300)'],
    products = ['O=C(O)C(F)C(F)(F)F(1188)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -5.4 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['O=[C]O(191)', 'CF3-CHF(94)'],
    products = ['O=C(O)C(F)C(F)(F)F(1188)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_4R!H->O',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_4R!H->O"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(5)', 'O=C(O)[C](F)C(F)(F)F(1301)'],
    products = ['O=C(O)C(F)C(F)(F)F(1188)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.21692e+34,'m^3/(mol*s)'), n=-8.80473, Ea=(5.22808,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_N-5CF->C_Sp-4C-2C_Ext-4C-R_Ext-2C-R_Ext-4C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_N-5CF->C_Sp-4C-2C_Ext-4C-R_Ext-2C-R_Ext-4C-R"""),
)

reaction(
    label = 'reaction21',
    reactants = ['F2(108)', 'O=C(O)C=C(F)F(947)'],
    products = ['O=C(O)C(F)C(F)(F)F(1188)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction22',
    reactants = ['HF(38)', 'O=C(O)C(F)=C(F)F(445)'],
    products = ['O=C(O)C(F)C(F)(F)F(1188)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(26.4943,'m^3/(mol*s)'), n=1.22463, Ea=(117.781,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-4COCdCddCtO2d-R_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-4COCdCddCtO2d-R_Ext-3COCdCddCtO2d-R"""),
)

network(
    label = 'PDepNetwork #450',
    isomers = [
        'O=C(O)C(F)C(F)(F)F(1188)',
    ],
    reactants = [
        ('HF(38)', 'CO2(13)', 'CHFCF2(54)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #450',
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

