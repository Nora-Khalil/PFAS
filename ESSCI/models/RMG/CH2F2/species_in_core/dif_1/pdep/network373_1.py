species(
    label = 'OC(F)(F)C(F)C(F)(F)F(883)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  O u0 p2 c0 {9,S} {12,S}
8  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
9  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
10 C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-1567.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,250,417,511,1155,1315,1456,3119,223,363,546,575,694,1179,1410,193,295,551,588,656,1146,1192,1350,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.193525,'amu*angstrom^2'), symmetry=1, barrier=(4.44951,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193498,'amu*angstrom^2'), symmetry=1, barrier=(4.4489,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04397,'amu*angstrom^2'), symmetry=1, barrier=(24.003,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (168.038,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.28937,0.0619178,-2.50251e-05,-7.29557e-08,7.09904e-11,-188498,14.5709], Tmin=(10,'K'), Tmax=(558.885,'K')), NASAPolynomial(coeffs=[10.105,0.0379128,-2.70916e-05,8.82668e-09,-1.07264e-12,-189647,-17.8065], Tmin=(558.885,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1567.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), label="""OC(F)(F)C(F)C(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'CF2O(48)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
"""),
    E0 = (-650.702,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(65.9917,'amu')),
        NonlinearRotor(inertia=([42.7382,43.0674,85.8056],'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([584.127,619.74,793.429,989.231,1281.66,1989.03],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (66.0068,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2914.22,'J/mol'), sigma=(4.906,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17614,0.021087,-2.39216e-05,1.39949e-08,-3.37605e-12,-77951.3,19.1134], Tmin=(298,'K'), Tmax=(1150,'K')), NASAPolynomial(coeffs=[6.79131,0.00338305,-1.43012e-06,2.71226e-10,-1.90524e-14,-79454,-9.48312], Tmin=(1150,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-650.702,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(83.1447,'J/mol/K'), label="""CF2O""", comment="""Thermo library: Fluorine"""),
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
    label = 'OC(F)(F)[C]F(576)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {6,S}
4 O u0 p2 c0 {5,S} {7,S}
5 C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6 C u0 p1 c0 {3,S} {5,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-539.334,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,617,898,1187],'cm^-1')),
        HinderedRotor(inertia=(1.38402,'amu*angstrom^2'), symmetry=1, barrier=(31.8212,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.39593,'amu*angstrom^2'), symmetry=1, barrier=(32.0951,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0239,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03035,0.046278,-6.6989e-05,4.89416e-08,-1.41201e-11,-64798.6,17.5213], Tmin=(100,'K'), Tmax=(850.047,'K')), NASAPolynomial(coeffs=[9.35662,0.011803,-6.15344e-06,1.22947e-09,-8.77629e-14,-66044.1,-16.6334], Tmin=(850.047,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-539.334,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCFFO) + group(CJ2_singlet-FCs)"""),
)

species(
    label = 'OC(F)F(164)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {4,S} {6,S}
4 C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
5 H u0 p0 c0 {4,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-687.011,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,519,593,830,1115,1166,1389,1413,3114],'cm^-1')),
        HinderedRotor(inertia=(0.761553,'amu*angstrom^2'), symmetry=1, barrier=(17.5096,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (68.0228,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3063.56,'J/mol'), sigma=(4.99758,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=478.52 K, Pc=55.69 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93609,0.00366784,6.01076e-05,-1.20272e-07,7.01069e-11,-82625.1,8.12343], Tmin=(10,'K'), Tmax=(593.851,'K')), NASAPolynomial(coeffs=[4.25321,0.0168559,-1.19109e-05,4.02989e-09,-5.14929e-13,-82933,4.48375], Tmin=(593.851,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-687.011,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""OC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[C]C(F)(F)F(307)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 F u0 p3 c0 {6,S}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
6 C u0 p1 c0 {4,S} {5,S}
"""),
    E0 = (-531.432,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([261,346,529,1231,1306,617,898,1187,310.428,1028.96,1029.05],'cm^-1')),
        HinderedRotor(inertia=(0.085726,'amu*angstrom^2'), symmetry=1, barrier=(5.82341,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (100.015,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88935,0.00712406,8.594e-05,-2.10572e-07,1.46542e-10,-63912.6,11.0368], Tmin=(10,'K'), Tmax=(506.569,'K')), NASAPolynomial(coeffs=[4.91752,0.0201688,-1.53534e-05,5.20649e-09,-6.50075e-13,-64288.3,4.0956], Tmin=(506.569,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-531.432,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""F[C]C(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CHF(88)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 C u0 p1 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (153.289,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1169.21,1416.41,2978.06],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (32.017,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2178.39,'J/mol'), sigma=(4.123,'angstroms'), dipoleMoment=(1.8,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.34484,0.00235461,1.93983e-06,-2.65251e-09,7.91169e-13,18514,7.05286], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[4.48366,0.00174964,-5.0479e-07,1.08953e-10,-9.87898e-15,17958.1,0.289222], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(153.289,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CHF""", comment="""Thermo library: Fluorine"""),
)

species(
    label = 'OC(F)(F)C(F)(F)F(1019)',
    structure = adjacencyList("""1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {8,S}
4 F u0 p3 c0 {8,S}
5 F u0 p3 c0 {8,S}
6 O u0 p2 c0 {7,S} {9,S}
7 C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
8 C u0 p0 c0 {3,S} {4,S} {5,S} {7,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-1350.04,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,223,363,546,575,694,1179,1410,193,295,551,588,656,1146,1192,1350,180],'cm^-1')),
        HinderedRotor(inertia=(0.138227,'amu*angstrom^2'), symmetry=1, barrier=(3.17811,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.407754,'amu*angstrom^2'), symmetry=1, barrier=(9.37507,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (136.021,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.57082,0.0364322,3.32484e-05,-1.72022e-07,1.4359e-10,-162369,12.2295], Tmin=(10,'K'), Tmax=(492.269,'K')), NASAPolynomial(coeffs=[7.65639,0.0287673,-2.11978e-05,7.0786e-09,-8.77171e-13,-163080,-7.72772], Tmin=(492.269,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1350.04,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""OC(F)(F)C(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'F[C]C(F)C(F)(F)F(1272)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {7,S}
4 F u0 p3 c0 {7,S}
5 F u0 p3 c0 {8,S}
6 C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7 C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
8 C u0 p1 c0 {5,S} {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-771.338,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([353,444,1253,3145,193,295,551,588,656,1146,1192,1350,617,898,1187,180,180,967.474,968.416],'cm^-1')),
        HinderedRotor(inertia=(0.170131,'amu*angstrom^2'), symmetry=1, barrier=(3.91165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169269,'amu*angstrom^2'), symmetry=1, barrier=(3.89182,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38836,0.0622948,-8.70898e-05,6.34576e-08,-1.85826e-11,-92680.8,21.6378], Tmin=(100,'K'), Tmax=(831.815,'K')), NASAPolynomial(coeffs=[10.2118,0.0198659,-1.05802e-05,2.13957e-09,-1.53996e-13,-94148.7,-19.3056], Tmin=(831.815,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-771.338,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCCFH) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-Cs(F)) + group(CJ2_singlet-FCs)"""),
)

species(
    label = 'OC(F)(F)C(F)F(435)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {7,S}
4 F u0 p3 c0 {7,S}
5 O u0 p2 c0 {6,S} {9,S}
6 C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
7 C u0 p0 c0 {3,S} {4,S} {6,S} {8,S}
8 H u0 p0 c0 {7,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-1115.54,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,223,363,546,575,694,1179,1410,235,523,627,1123,1142,1372,1406,3097,180],'cm^-1')),
        HinderedRotor(inertia=(0.143577,'amu*angstrom^2'), symmetry=1, barrier=(3.30111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.610859,'amu*angstrom^2'), symmetry=1, barrier=(14.0448,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (118.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.63191,0.0316143,2.37687e-05,-1.18291e-07,9.27336e-11,-134166,12.5198], Tmin=(10,'K'), Tmax=(516.15,'K')), NASAPolynomial(coeffs=[6.39527,0.0284113,-1.98498e-05,6.40832e-09,-7.7625e-13,-134694,-1.33427], Tmin=(516.15,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1115.54,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""OC(F)(F)C(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OC(C(F)(F)F)C(F)(F)F(1273)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  O u0 p2 c0 {8,S} {12,S}
8  C u0 p0 c0 {7,S} {9,S} {10,S} {11,S}
9  C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
10 C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-1573.75,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (168.038,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.5596,0.0320358,0.00022874,-7.46347e-07,6.43169e-10,-189273,14.1177], Tmin=(10,'K'), Tmax=(432.114,'K')), NASAPolynomial(coeffs=[9.91238,0.039932,-3.02158e-05,1.0401e-08,-1.32427e-12,-190445,-18.4062], Tmin=(432.114,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1573.75,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), label="""OC(C(F)(F)F)C(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'FC(F)=C(F)C(F)(F)F(879)',
    structure = adjacencyList("""1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {7,S}
4 F u0 p3 c0 {8,S}
5 F u0 p3 c0 {9,S}
6 F u0 p3 c0 {9,S}
7 C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
8 C u0 p0 c0 {4,S} {7,S} {9,D}
9 C u0 p0 c0 {5,S} {6,S} {8,D}
"""),
    E0 = (-1171.65,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([219,296,586,564,718,793,1177,1228,323,467,575,827,1418,182,240,577,636,1210,1413,2018.52],'cm^-1')),
        HinderedRotor(inertia=(0.370476,'amu*angstrom^2'), symmetry=1, barrier=(8.51796,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (150.022,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.48756,0.0445189,-1.69994e-05,-4.62355e-08,4.1721e-11,-140912,13.4084], Tmin=(10,'K'), Tmax=(600.481,'K')), NASAPolynomial(coeffs=[8.64988,0.0279311,-2.00276e-05,6.49183e-09,-7.83507e-13,-141852,-11.5361], Tmin=(600.481,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1171.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""FC(F)DC(F)C(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC(F)=CC(F)(F)F(857)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {6,S}
4 F u0 p3 c0 {8,S}
5 F u0 p3 c0 {8,S}
6 C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
7 C u0 p0 c0 {6,S} {8,D} {9,S}
8 C u0 p0 c0 {4,S} {5,S} {7,D}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-1028.4,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([219,296,586,564,718,793,1177,1228,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,1134.66],'cm^-1')),
        HinderedRotor(inertia=(0.133162,'amu*angstrom^2'), symmetry=1, barrier=(3.06166,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (132.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.71093,0.0287928,1.91777e-05,-6.84733e-08,3.95583e-11,-123686,13.1392], Tmin=(10,'K'), Tmax=(697.143,'K')), NASAPolynomial(coeffs=[7.03966,0.0281026,-1.8947e-05,5.86295e-09,-6.82491e-13,-124598,-4.92755], Tmin=(697.143,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1028.4,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""FC(F)DCC(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'OC(F)(F)[CH]C(F)(F)F(1274)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {8,S}
6  O u0 p2 c0 {7,S} {11,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
8  C u0 p0 c0 {3,S} {4,S} {5,S} {9,S}
9  C u1 p0 c0 {7,S} {8,S} {10,S}
10 H u0 p0 c0 {9,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-1198.84,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,253,525,597,667,842,1178,1324,253,522,585,610,849,1160,1215,1402,3025,407.5,1350,352.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.784803,'amu*angstrom^2'), symmetry=1, barrier=(18.0442,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.235227,'amu*angstrom^2'), symmetry=1, barrier=(5.40833,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.235293,'amu*angstrom^2'), symmetry=1, barrier=(5.40984,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (149.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.3783,0.054081,-1.52178e-05,-8.14086e-08,7.56819e-11,-144182,17.033], Tmin=(10,'K'), Tmax=(546.341,'K')), NASAPolynomial(coeffs=[9.2864,0.0341403,-2.44821e-05,8.00658e-09,-9.76148e-13,-145176,-11.0832], Tmin=(546.341,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1198.84,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), label="""OC(F)(F)[CH]C(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O[C](F)C(F)C(F)(F)F(1275)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {9,S}
6  O u0 p2 c0 {9,S} {11,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
9  C u1 p0 c0 {5,S} {6,S} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-1147.68,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,259,529,569,1128,1321,1390,3140,193,295,551,588,656,1146,1192,1350,395,473,707,1436,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.581505,'amu*angstrom^2'), symmetry=1, barrier=(13.3699,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.581985,'amu*angstrom^2'), symmetry=1, barrier=(13.381,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.580777,'amu*angstrom^2'), symmetry=1, barrier=(13.3532,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (149.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.52172,0.0368576,0.000181779,-6.75324e-07,6.24722e-10,-138034,14.2779], Tmin=(10,'K'), Tmax=(411.977,'K')), NASAPolynomial(coeffs=[9.73647,0.0353683,-2.7076e-05,9.39485e-09,-1.20235e-12,-139045,-16.2549], Tmin=(411.977,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1147.68,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), label="""O[C](F)C(F)C(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OC(F)(F)C(F)[C](F)F(1276)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  O u0 p2 c0 {8,S} {11,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
9  C u1 p0 c0 {4,S} {5,S} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-1135.36,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,259,529,569,1128,1321,1390,3140,223,363,546,575,694,1179,1410,190,488,555,1236,1407,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.273615,'amu*angstrom^2'), symmetry=1, barrier=(6.29094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.841097,'amu*angstrom^2'), symmetry=1, barrier=(19.3385,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.273468,'amu*angstrom^2'), symmetry=1, barrier=(6.28757,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (149.039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.22333,0.0696089,-9.75326e-05,7.23587e-08,-2.18061e-11,-136548,14.9479], Tmin=(10,'K'), Tmax=(785.896,'K')), NASAPolynomial(coeffs=[11.2707,0.02865,-1.93563e-05,6.04254e-09,-7.10357e-13,-137813,-21.9371], Tmin=(785.896,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1135.36,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), label="""OC(F)(F)C(F)[C](F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'F[C](F)C(F)C(F)(F)F(1277)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {9,S}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
9  C u1 p0 c0 {5,S} {6,S} {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-1151.57,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,193,295,551,588,656,1146,1192,1350,190,488,555,1236,1407,180,2422.61],'cm^-1')),
        HinderedRotor(inertia=(0.536938,'amu*angstrom^2'), symmetry=1, barrier=(12.3453,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.26266,'amu*angstrom^2'), symmetry=1, barrier=(29.0311,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (151.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.4642,0.0569969,-6.78565e-05,4.2609e-08,-1.10759e-11,-138501,14.3457], Tmin=(10,'K'), Tmax=(878.417,'K')), NASAPolynomial(coeffs=[9.61596,0.0289839,-2.00211e-05,6.30476e-09,-7.43629e-13,-139582,-14.5355], Tmin=(878.417,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1151.57,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""F[C](F)C(F)C(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]C(F)(F)C(F)C(F)(F)F(1278)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {10,S}
5  F u0 p3 c0 {10,S}
6  F u0 p3 c0 {10,S}
7  O u1 p2 c0 {9,S}
8  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
9  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
10 C u0 p0 c0 {4,S} {5,S} {6,S} {8,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-1310.09,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([250,417,511,1155,1315,1456,3119,351,323,533,609,664,892,1120,1201,193,295,551,588,656,1146,1192,1350,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0971642,'amu*angstrom^2'), symmetry=1, barrier=(2.234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.709412,'amu*angstrom^2'), symmetry=1, barrier=(16.3108,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (167.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.29587,0.0634404,-6.01909e-05,1.69587e-08,4.03472e-12,-157562,14.7031], Tmin=(10,'K'), Tmax=(748.672,'K')), NASAPolynomial(coeffs=[12.2704,0.0294788,-2.01725e-05,6.27903e-09,-7.32297e-13,-159298,-28.6143], Tmin=(748.672,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1310.09,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), label="""[O]C(F)(F)C(F)C(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O[C](F)F(220)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {4,S} {5,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (-470.293,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,493,600,700,1144,1293],'cm^-1')),
        HinderedRotor(inertia=(0.307072,'amu*angstrom^2'), symmetry=1, barrier=(7.06019,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.0148,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3063.56,'J/mol'), sigma=(4.99758,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=478.52 K, Pc=55.69 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.94331,0.003562,5.16413e-05,-1.18327e-07,7.89485e-11,-56559.5,10.5484], Tmin=(10,'K'), Tmax=(519.044,'K')), NASAPolynomial(coeffs=[4.2538,0.012858,-9.00307e-06,2.95251e-09,-3.63989e-13,-56749.2,7.73732], Tmin=(519.044,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-470.293,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""O[C](F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'OC(F)(F)[CH]F(637)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {6,S}
4 O u0 p2 c0 {5,S} {8,S}
5 C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6 C u1 p0 c0 {3,S} {5,S} {7,S}
7 H u0 p0 c0 {6,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-703.673,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,253,525,597,667,842,1178,1324,334,575,1197,1424,3202,180],'cm^-1')),
        HinderedRotor(inertia=(0.516577,'amu*angstrom^2'), symmetry=1, barrier=(11.8771,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.37678,'amu*angstrom^2'), symmetry=1, barrier=(31.6549,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (99.0318,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.79577,0.0139836,0.000131184,-3.74612e-07,2.93964e-10,-84624.9,11.6058], Tmin=(10,'K'), Tmax=(464.833,'K')), NASAPolynomial(coeffs=[6.85848,0.023816,-1.73218e-05,5.85757e-09,-7.39626e-13,-85300.6,-5.02903], Tmin=(464.833,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-703.673,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""OC(F)(F)[CH]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OC(F)(F)[C](F)C(F)(F)F(1279)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  F u0 p3 c0 {10,S}
7  O u0 p2 c0 {8,S} {11,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {10,S}
9  C u0 p0 c0 {3,S} {4,S} {5,S} {10,S}
10 C u1 p0 c0 {6,S} {8,S} {9,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-1367.81,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,253,525,597,667,842,1178,1324,253,522,585,610,849,1160,1215,1402,212,367,445,1450,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.209374,'amu*angstrom^2'), symmetry=1, barrier=(4.81393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.209341,'amu*angstrom^2'), symmetry=1, barrier=(4.81317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.209356,'amu*angstrom^2'), symmetry=1, barrier=(4.81351,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (167.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.41756,0.0717297,0.00020306,-2.74937e-06,6.58302e-09,-164509,15.6032], Tmin=(10,'K'), Tmax=(175.741,'K')), NASAPolynomial(coeffs=[7.27623,0.0415986,-3.22097e-05,1.11931e-08,-1.43327e-12,-164734,1.16173], Tmin=(175.741,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1367.81,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), label="""OC(F)(F)[C](F)C(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C(F)C(F)C(F)(F)F(513)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {9,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
9  C u0 p0 c0 {5,S} {6,D} {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-1271.85,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([206,363,518,1175,1317,1481,3059,193,295,551,588,656,1146,1192,1350,486,617,768,1157,1926,180,624.731],'cm^-1')),
        HinderedRotor(inertia=(0.229656,'amu*angstrom^2'), symmetry=1, barrier=(5.28024,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.588251,'amu*angstrom^2'), symmetry=1, barrier=(13.5251,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (148.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.45138,0.0554271,-6.36074e-05,3.79523e-08,-9.28041e-12,-152967,13.2558], Tmin=(10,'K'), Tmax=(938.495,'K')), NASAPolynomial(coeffs=[10.1572,0.0268458,-1.79258e-05,5.50192e-09,-6.36127e-13,-154226,-18.6702], Tmin=(938.495,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1271.85,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""ODC(F)C(F)C(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'OC(F)=CC(F)(F)F(1280)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {8,S} {10,S}
6  C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
7  C u0 p0 c0 {6,S} {8,D} {9,S}
8  C u0 p0 c0 {4,S} {5,S} {7,D}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-1032.68,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,219,296,586,564,718,793,1177,1228,3010,987.5,1337.5,450,1655,326,540,652,719,1357,180],'cm^-1')),
        HinderedRotor(inertia=(0.0978146,'amu*angstrom^2'), symmetry=1, barrier=(16.5142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.59857,'amu*angstrom^2'), symmetry=1, barrier=(36.7542,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.78812,0.0136763,0.000164941,-4.02256e-07,2.79349e-10,-124200,12.7636], Tmin=(10,'K'), Tmax=(505.974,'K')), NASAPolynomial(coeffs=[5.53615,0.0396626,-3.01038e-05,1.02169e-08,-1.27579e-12,-124886,0.486096], Tmin=(505.974,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1032.68,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""OC(F)DCC(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OC(F)=C(F)C(F)(F)F(1281)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {8,S}
5  F u0 p3 c0 {9,S}
6  O u0 p2 c0 {9,S} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
8  C u0 p0 c0 {4,S} {7,S} {9,D}
9  C u0 p0 c0 {5,S} {6,S} {8,D}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-1175.06,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,219,296,586,564,718,793,1177,1228,323,467,575,827,1418,326,540,652,719,1357,180],'cm^-1')),
        HinderedRotor(inertia=(0.371381,'amu*angstrom^2'), symmetry=1, barrier=(8.53878,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.370822,'amu*angstrom^2'), symmetry=1, barrier=(8.52593,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (148.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.62611,0.0272005,0.000179856,-5.83927e-07,4.96534e-10,-141325,13.4436], Tmin=(10,'K'), Tmax=(438.608,'K')), NASAPolynomial(coeffs=[8.85209,0.0334496,-2.5878e-05,9.00461e-09,-1.15195e-12,-142302,-13.373], Tmin=(438.608,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1175.06,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""OC(F)DC(F)C(F)(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OC(F)(F)C=C(F)F(912)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {6,S} {10,S}
6  C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
7  C u0 p0 c0 {6,S} {8,D} {9,S}
8  C u0 p0 c0 {3,S} {4,S} {7,D}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {5,S}
"""),
    E0 = (-1007.4,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,275,321,533,585,746,850,1103,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,180],'cm^-1')),
        HinderedRotor(inertia=(0.0170999,'amu*angstrom^2'), symmetry=1, barrier=(5.81748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.506554,'amu*angstrom^2'), symmetry=1, barrier=(11.6467,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (130.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.7482,0.0172572,0.00016515,-4.54346e-07,3.48046e-10,-121159,14.6152], Tmin=(10,'K'), Tmax=(468.635,'K')), NASAPolynomial(coeffs=[6.60543,0.034901,-2.58576e-05,8.75963e-09,-1.10002e-12,-121888,-1.92801], Tmin=(468.635,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1007.4,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""OC(F)(F)CDC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'OC(F)(F)C(F)=C(F)F(382)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  O u0 p2 c0 {7,S} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
8  C u0 p0 c0 {3,S} {7,S} {9,D}
9  C u0 p0 c0 {4,S} {5,S} {8,D}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-1154.51,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,275,321,533,585,746,850,1103,323,467,575,827,1418,182,240,577,636,1210,1413,180],'cm^-1')),
        HinderedRotor(inertia=(0.102433,'amu*angstrom^2'), symmetry=1, barrier=(2.35514,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0987329,'amu*angstrom^2'), symmetry=1, barrier=(2.27006,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (148.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.56187,0.035635,0.000138624,-5.6967e-07,5.63035e-10,-138858,13.4425], Tmin=(10,'K'), Tmax=(387.981,'K')), NASAPolynomial(coeffs=[8.14435,0.0331468,-2.4791e-05,8.45085e-09,-1.06863e-12,-139551,-8.66771], Tmin=(387.981,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1154.51,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""OC(F)(F)C(F)DC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    E0 = (-670.4,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-548.84,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-526.053,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-336.958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-467.063,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-214.092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-414.987,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-626.49,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-549.788,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-250.053,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-459.192,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-480.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-432.967,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-420.65,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-481.351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-456.461,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-526.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-541.305,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-514.185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-749.504,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-392.874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-638.887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-369.259,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (-652.957,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['OC(F)(F)C(F)C(F)(F)F(883)'],
    products = ['HF(38)', 'CF2O(48)', 'CHFCF2(54)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5.95626e+10,'s^-1'), n=0.978818, Ea=(255.089,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R!H->C_N-5Br1sCl1sF1sH->H_5Br1sCl1sF1s->F1s_Ext-2C-R_N-7R!H->O_Sp-7BrCClFINPSSi-2C_Ext-2C-R',), comment="""Estimated from node Root_1R!H->C_N-5Br1sCl1sF1sH->H_5Br1sCl1sF1s->F1s_Ext-2C-R_N-7R!H->O_Sp-7BrCClFINPSSi-2C_Ext-2C-R
Multiplied by reaction path degeneracy 3.0"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CHF3(41)', 'OC(F)(F)[C]F(576)'],
    products = ['OC(F)(F)C(F)C(F)(F)F(883)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.33872e-07,'m^3/(mol*s)'), n=3.38172, Ea=(57.5279,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_N-3Br1sCCl1sF1sHI1s->F1s_N-3Br1sCCl1sH->Cl1s_2Br1sCl1sF1sHI1s->F1s',), comment="""Estimated from node CH_N-3Br1sCCl1sF1sHI1s->F1s_N-3Br1sCCl1sH->Cl1s_2Br1sCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction3',
    reactants = ['OC(F)F(164)', 'F[C]C(F)(F)F(307)'],
    products = ['OC(F)(F)C(F)C(F)(F)F(883)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.33872e-07,'m^3/(mol*s)'), n=3.38172, Ea=(50.5689,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CH_N-3Br1sCCl1sF1sHI1s->F1s_N-3Br1sCCl1sH->Cl1s_2Br1sCl1sF1sHI1s->F1s',), comment="""Estimated from node CH_N-3Br1sCCl1sF1sHI1s->F1s_N-3Br1sCCl1sH->Cl1s_2Br1sCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CHF(88)', 'OC(F)(F)C(F)(F)F(1019)'],
    products = ['OC(F)(F)C(F)C(F)(F)F(883)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.7066e-06,'m^3/(mol*s)'), n=3.19155, Ea=(217.971,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CC_2Br1sCl1sF1sHI1s->H_3Br1sCCl1sF1sHI1s->F1s',), comment="""Estimated from node CC_2Br1sCl1sF1sHI1s->H_3Br1sCCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CF2(89)', 'OC(F)C(F)(F)F(1271)'],
    products = ['OC(F)(F)C(F)C(F)(F)F(883)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.40908e-07,'m^3/(mol*s)'), n=3.45227, Ea=(206.667,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s',), comment="""Estimated from node CO_2Br1sCl1sF1sHI1s->F1s_3Br1sCCl1sF1sHI1s->F1s"""),
)

reaction(
    label = 'reaction6',
    reactants = ['OF(155)', 'F[C]C(F)C(F)(F)F(1272)'],
    products = ['OC(F)(F)C(F)C(F)(F)F(883)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.13916e+09,'m^3/(mol*s)'), n=-1.00842, Ea=(10.6903,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17406617270594374, var=16.167420070870673, Tref=1000.0, N=4, data_mean=0.0, correlation='OHOY',), comment="""Estimated from node OHOY"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CF2(89)', 'OC(F)(F)C(F)F(435)'],
    products = ['OC(F)(F)C(F)C(F)(F)F(883)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.67164e-06,'m^3/(mol*s)'), n=3.3552, Ea=(255.635,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='CY_N-2Br1sCl1sF1sHI1s->H_Ext-4Cs-R',), comment="""Estimated from node CY_N-2Br1sCl1sF1sHI1s->H_Ext-4Cs-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction8',
    reactants = ['OC(F)(F)C(F)C(F)(F)F(883)'],
    products = ['OC(C(F)(F)F)C(F)(F)F(1273)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H2O(3)', 'FC(F)=C(F)C(F)(F)F(879)'],
    products = ['OC(F)(F)C(F)C(F)(F)F(883)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.00057209,'m^3/(mol*s)'), n=2.81783, Ea=(231.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OH] for rate rule [Cd/disub_Cd/disub;H_OH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction10',
    reactants = ['OF(155)', 'FC(F)=CC(F)(F)F(857)'],
    products = ['OC(F)(F)C(F)C(F)(F)F(883)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(0.000286045,'m^3/(mol*s)'), n=2.81783, Ea=(231.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;H_OH] for rate rule [Cd/monosub_Cd/disub;H_OH]
Euclidian distance = 1.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction11',
    reactants = ['OC(F)(F)F(254)', 'CHFCF2(54)'],
    products = ['OC(F)(F)C(F)C(F)(F)F(883)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.79e-11,'m^3/(mol*s)'), n=3.97, Ea=(329.281,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_Cd;R_OH] for rate rule [Cd/monosub_Cd/disub;R_OH]
Euclidian distance = 1.0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', 'OC(F)(F)[CH]C(F)(F)F(1274)'],
    products = ['OC(F)(F)C(F)C(F)(F)F(883)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(5e+07,'m^3/(mol*s)'), n=0, Ea=(3.21324,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_N-2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_N-3R!H->F_N-3BrCClINOPSSi->Cl_Ext-1C-R_Ext-3BrCO-R_N-2CF->C"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', 'O[C](F)C(F)C(F)(F)F(1275)'],
    products = ['OC(F)(F)C(F)C(F)(F)F(883)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C
Ea raised from -0.8 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction14',
    reactants = ['F(37)', 'OC(F)(F)C(F)[C](F)F(1276)'],
    products = ['OC(F)(F)C(F)C(F)(F)F(883)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_N-2CF->C"""),
)

reaction(
    label = 'reaction15',
    reactants = ['OH(4)', 'F[C](F)C(F)C(F)(F)F(1277)'],
    products = ['OC(F)(F)C(F)C(F)(F)F(883)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.7e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(5)', '[O]C(F)(F)C(F)C(F)(F)F(1278)'],
    products = ['OC(F)(F)C(F)C(F)(F)F(883)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O[C](F)F(220)', 'CF3-CHF(94)'],
    products = ['OC(F)(F)C(F)C(F)(F)F(883)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -16.1 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['CF3(44)', 'OC(F)(F)[CH]F(637)'],
    products = ['OC(F)(F)C(F)C(F)(F)F(883)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -12.2 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(5)', 'OC(F)(F)[C](F)C(F)(F)F(1279)'],
    products = ['OC(F)(F)C(F)C(F)(F)F(883)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.21692e+34,'m^3/(mol*s)'), n=-8.80473, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_N-5CF->C_Sp-4C-2C_Ext-4C-R_Ext-2C-R_Ext-4C-R',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_N-5CF->C_Sp-4C-2C_Ext-4C-R_Ext-2C-R_Ext-4C-R
Ea raised from -0.9 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['HF(38)', 'O=C(F)C(F)C(F)(F)F(513)'],
    products = ['OC(F)(F)C(F)C(F)(F)F(883)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.109156,'m^3/(mol*s)'), n=1.86531, Ea=(161.64,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd',), comment="""Estimated from node HF_N-3COCdCddCtO2d->Ct_N-3CdO2d->Cd_N-4COCdCddCtO2d->Cdd"""),
)

reaction(
    label = 'reaction21',
    reactants = ['F2(108)', 'OC(F)=CC(F)(F)F(1280)'],
    products = ['OC(F)(F)C(F)C(F)(F)F(883)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(6.79458,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction22',
    reactants = ['HF(38)', 'OC(F)=C(F)C(F)(F)F(1281)'],
    products = ['OC(F)(F)C(F)C(F)(F)F(883)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.14111,'m^3/(mol*s)'), n=1.29695, Ea=(175.47,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction23',
    reactants = ['F2(108)', 'OC(F)(F)C=C(F)F(912)'],
    products = ['OC(F)(F)C(F)C(F)(F)F(883)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(5.1235,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction24',
    reactants = ['HF(38)', 'OC(F)(F)C(F)=C(F)F(382)'],
    products = ['OC(F)(F)C(F)C(F)(F)F(883)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(26.4943,'m^3/(mol*s)'), n=1.22463, Ea=(140.849,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-4COCdCddCtO2d-R_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-4COCdCddCtO2d-R_Ext-3COCdCddCtO2d-R"""),
)

network(
    label = 'PDepNetwork #373',
    isomers = [
        'OC(F)(F)C(F)C(F)(F)F(883)',
    ],
    reactants = [
        ('HF(38)', 'CF2O(48)', 'CHFCF2(54)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #373',
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

