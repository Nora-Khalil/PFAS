species(
    label = '[O]C(F)C([O])(F)F(429)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {7,S}
4 O u1 p2 c0 {6,S}
5 O u1 p2 c0 {7,S}
6 C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
7 C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-613.996,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([391,562,707,872,1109,1210,1289,3137,351,323,533,609,664,892,1120,1201,180],'cm^-1')),
        HinderedRotor(inertia=(1.5766,'amu*angstrom^2'), symmetry=1, barrier=(36.249,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.023,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41151,0.0470069,-3.44331e-05,1.94681e-09,4.4207e-12,-73744.7,20.6885], Tmin=(100,'K'), Tmax=(1020.07,'K')), NASAPolynomial(coeffs=[17.415,0.00468127,-2.23431e-06,5.36432e-10,-4.5404e-14,-78072.5,-62.0468], Tmin=(1020.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-613.996,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sH)) + radical(O2sj(Cs-CsF1sF1s))"""),
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
    label = 'CHFO(46)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 H u0 p0 c0 {3,S}
"""),
    E0 = (-386.87,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(48.0011,'amu')),
        NonlinearRotor(inertia=([5.46358,42.889,48.3526],'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([674.139,1050.3,1121.98,1395.32,1906.92,3070.52],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (48.0164,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2914.22,'J/mol'), sigma=(4.906,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.1185,0.010774,-5.51582e-06,5.81519e-10,2.17042e-13,-46356.7,14.6524], Tmin=(298,'K'), Tmax=(1150,'K')), NASAPolynomial(coeffs=[4.81925,0.00505418,-2.01421e-06,3.63591e-10,-2.44724e-14,-47263.1,0.0972251], Tmin=(1150,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-386.87,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(83.1447,'J/mol/K'), label="""CHFO""", comment="""Thermo library: Fluorine"""),
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
    label = '[O]C(F)(F)[CH]F(474)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {6,S}
4 O u1 p2 c0 {5,S}
5 C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6 C u1 p0 c0 {3,S} {5,S} {7,S}
7 H u0 p0 c0 {6,S}
"""),
    E0 = (-445.565,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,525,597,667,842,1178,1324,334,575,1197,1424,3202,180,395.38],'cm^-1')),
        HinderedRotor(inertia=(0.422093,'amu*angstrom^2'), symmetry=1, barrier=(9.70474,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0239,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.99253,0.0473334,-6.8592e-05,5.03378e-08,-1.45679e-11,-53519.7,18.0995], Tmin=(100,'K'), Tmax=(847.727,'K')), NASAPolynomial(coeffs=[9.47106,0.0120471,-6.15724e-06,1.23976e-09,-8.91188e-14,-54787.7,-16.745], Tmin=(847.727,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-445.565,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sF1s)) + radical(Csj(Cs-F1sF1sO2s)(F1s)(H))"""),
)

species(
    label = '[O]C(F)[C](F)F(475)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {6,S}
4 O u1 p2 c0 {5,S}
5 C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
6 C u1 p0 c0 {2,S} {3,S} {5,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-445.931,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([391,562,707,872,1109,1210,1289,3137,190,488,555,1236,1407,180],'cm^-1')),
        HinderedRotor(inertia=(0.813996,'amu*angstrom^2'), symmetry=1, barrier=(18.7154,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0239,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.14839,0.0397111,-4.32916e-05,2.26821e-08,-4.63769e-12,-53565.5,18.7797], Tmin=(100,'K'), Tmax=(1191.22,'K')), NASAPolynomial(coeffs=[11.353,0.00880293,-4.37161e-06,9.00449e-10,-6.63946e-14,-55758.5,-27.2377], Tmin=(1191.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-445.931,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sH)) + radical(Csj(Cs-F1sO2sH)(F1s)(F1s))"""),
)

species(
    label = 'FC1OOC1(F)F(433)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {7,S}
4 O u0 p2 c0 {5,S} {6,S}
5 O u0 p2 c0 {4,S} {7,S}
6 C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
7 C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-675.311,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.023,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88387,0.00730685,0.000115597,-2.58497e-07,1.70074e-10,-81217.5,12.094], Tmin=(10,'K'), Tmax=(516.337,'K')), NASAPolynomial(coeffs=[3.72286,0.0317812,-2.29787e-05,7.54622e-09,-9.21745e-13,-81510.4,9.76617], Tmin=(516.337,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-675.311,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), label="""FC1OOC1(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C(F)C(O)(F)F(385)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {7,S}
4 O u0 p2 c0 {6,S} {8,S}
5 O u0 p2 c0 {7,D}
6 C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
7 C u0 p0 c0 {3,S} {5,D} {6,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-1029.36,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,251,367,519,700,855,1175,1303,486,617,768,1157,1926,180],'cm^-1')),
        HinderedRotor(inertia=(1.01043,'amu*angstrom^2'), symmetry=1, barrier=(23.2317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01274,'amu*angstrom^2'), symmetry=1, barrier=(23.2848,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.023,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.77229,0.0142182,0.000141818,-3.61218e-07,2.49514e-10,-123789,11.5845], Tmin=(10,'K'), Tmax=(533.644,'K')), NASAPolynomial(coeffs=[8.70117,0.0245747,-2.02509e-05,7.35176e-09,-9.67885e-13,-124989,-15.4095], Tmin=(533.644,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1029.36,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""ODC(F)C(O)(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]C(F)(F)C=O(476)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u1 p2 c0 {5,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
6 C u0 p0 c0 {4,D} {5,S} {7,S}
7 H u0 p0 c0 {6,S}
"""),
    E0 = (-518.002,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([377,639,914,1335,2782.5,750,1395,475,1775,1000,385.069,385.701,385.833,385.973],'cm^-1')),
        HinderedRotor(inertia=(0.0011382,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.0249,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.51133,0.0280011,-1.4965e-05,9.19161e-10,8.40297e-13,-62243.8,16.8477], Tmin=(100,'K'), Tmax=(1344.01,'K')), NASAPolynomial(coeffs=[10.7213,0.0119437,-6.3927e-06,1.30435e-09,-9.39321e-14,-65207.2,-28.0024], Tmin=(1344.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-518.002,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-F1sF1sCO))"""),
)

species(
    label = '[O][C](F)F(268)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u1 p2 c0 {4,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
"""),
    E0 = (-255.477,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([493,600,700,1144,1293,180],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.0068,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.105,0.0167699,-1.53787e-05,4.83895e-09,3.91619e-15,-30692.1,11.7073], Tmin=(100,'K'), Tmax=(1048.23,'K')), NASAPolynomial(coeffs=[8.58998,0.00108971,-4.53613e-07,1.24874e-10,-1.13715e-14,-32130.4,-16.3887], Tmin=(1048.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-255.477,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsOJ) + radical(CsF1sF1sO2s)"""),
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
    label = '[O]C(F)(F)C(=O)F(290)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {7,S}
4 O u1 p2 c0 {6,S}
5 O u0 p2 c0 {7,D}
6 C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
7 C u0 p0 c0 {3,S} {5,D} {6,S}
"""),
    E0 = (-763.591,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([377,639,914,1335,486,617,768,1157,1926,318.885,318.886,318.887,318.889,318.89],'cm^-1')),
        HinderedRotor(inertia=(0.00165779,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (113.015,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.74064,0.0210189,6.45452e-05,-2.04623e-07,1.54562e-10,-91835.6,12.9926], Tmin=(10,'K'), Tmax=(503.021,'K')), NASAPolynomial(coeffs=[6.77249,0.0227039,-1.73971e-05,5.91859e-09,-7.40694e-13,-92466.9,-2.79481], Tmin=(503.021,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-763.591,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""[O]C(F)(F)C(DO)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O]C(F)C(=O)F(477)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u1 p2 c0 {5,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
6 C u0 p0 c0 {2,S} {4,D} {5,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-568.509,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([241,322,589,815,1096,1220,1302,2892,486,617,768,1157,1926,180],'cm^-1')),
        HinderedRotor(inertia=(0.398498,'amu*angstrom^2'), symmetry=1, barrier=(9.16225,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.0249,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.78101,0.0174468,7.58508e-05,-2.6043e-07,2.29389e-10,-68375.7,11.5523], Tmin=(10,'K'), Tmax=(418.42,'K')), NASAPolynomial(coeffs=[5.46396,0.022233,-1.61424e-05,5.37722e-09,-6.6812e-13,-68699.3,2.71568], Tmin=(418.42,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-568.509,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""[O]C(F)C(DO)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O][CH]F(135)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 O u1 p2 c0 {3,S}
3 C u1 p0 c0 {1,S} {2,S} {4,S}
4 H u0 p0 c0 {3,S}
"""),
    E0 = (-19.6796,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([580,1155,1237,1373,3147,180],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (48.0164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53592,0.00799825,-5.22271e-07,-4.56612e-09,2.08752e-12,-2348.33,9.31392], Tmin=(100,'K'), Tmax=(1079.18,'K')), NASAPolynomial(coeffs=[5.97187,0.0038822,-1.62979e-06,3.36437e-10,-2.54188e-14,-3160.18,-3.94923], Tmin=(1079.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-19.6796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsOJ) + radical(CsF1sHO2s)"""),
)

species(
    label = 'F2(77)',
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
    label = '[O][CH]C(=O)F(478)',
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
    label = '[O][C](F)C(=O)F(291)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {5,S}
3 O u1 p2 c0 {5,S}
4 O u0 p2 c0 {6,D}
5 C u1 p0 c0 {2,S} {3,S} {6,S}
6 C u0 p0 c0 {1,S} {4,D} {5,S}
"""),
    E0 = (-436.241,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([280,501,1494,1531,611,648,830,1210,1753,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.353696,'amu*angstrom^2'), symmetry=1, barrier=(8.13216,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.0169,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.23918,0.0462587,-9.06344e-05,8.75347e-08,-3.13139e-11,-52411.4,15.9354], Tmin=(100,'K'), Tmax=(871.842,'K')), NASAPolynomial(coeffs=[5.13072,0.0156267,-8.05455e-06,1.54248e-09,-1.04518e-13,-52255.6,6.16692], Tmin=(871.842,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-436.241,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-F1sCOH)) + radical(CsCOF1sO2s)"""),
)

species(
    label = '[O]C(F)(F)[C](O)F(479)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {7,S}
4 O u0 p2 c0 {7,S} {8,S}
5 O u1 p2 c0 {6,S}
6 C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
7 C u1 p0 c0 {3,S} {4,S} {6,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-644.913,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.023,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.807835,0.0793785,-0.000151236,1.36197e-07,-4.60249e-11,-77458.9,21.0273], Tmin=(100,'K'), Tmax=(869.465,'K')), NASAPolynomial(coeffs=[10.5251,0.0156721,-8.54762e-06,1.65485e-09,-1.12161e-13,-78430.4,-20.3633], Tmin=(869.465,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-644.913,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sF1s)) + radical(Csj(Cs-F1sF1sO2s)(F1s)(O2s-H))"""),
)

species(
    label = '[O][C](F)C(O)(F)F(480)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {7,S}
4 O u0 p2 c0 {6,S} {8,S}
5 O u1 p2 c0 {7,S}
6 C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
7 C u1 p0 c0 {3,S} {5,S} {6,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-678.307,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (114.023,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.824466,0.079627,-0.000153815,1.39251e-07,-4.69627e-11,-81476.4,21.2815], Tmin=(100,'K'), Tmax=(882.453,'K')), NASAPolynomial(coeffs=[10.1935,0.0155755,-8.25199e-06,1.56697e-09,-1.04433e-13,-82289.6,-17.9854], Tmin=(882.453,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-678.307,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sH)) + radical(Csj(Cs-F1sF1sO2s)(F1s)(O2s-H))"""),
)

species(
    label = '[O]C(F)(F)[CH]OF(481)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {4,S}
4 O u0 p2 c0 {3,S} {7,S}
5 O u1 p2 c0 {6,S}
6 C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
7 C u1 p0 c0 {4,S} {6,S} {8,S}
8 H u0 p0 c0 {7,S}
"""),
    E0 = (-311.514,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,253,525,597,667,842,1178,1324,3025,407.5,1350,352.5,180,180,681.494],'cm^-1')),
        HinderedRotor(inertia=(0.409356,'amu*angstrom^2'), symmetry=1, barrier=(9.41191,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.771866,'amu*angstrom^2'), symmetry=1, barrier=(17.7467,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.023,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17354,0.0672757,-0.000107886,8.50348e-08,-2.62117e-11,-37369.4,22.1676], Tmin=(100,'K'), Tmax=(799.119,'K')), NASAPolynomial(coeffs=[11.8133,0.0140166,-7.91135e-06,1.62832e-09,-1.17535e-13,-39069.8,-26.7767], Tmin=(799.119,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-311.514,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sF1s)) + radical(CCsJO)"""),
)

species(
    label = '[O][CH]C(F)(F)OF(482)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {4,S}
4 O u0 p2 c0 {3,S} {6,S}
5 O u1 p2 c0 {7,S}
6 C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
7 C u1 p0 c0 {5,S} {6,S} {8,S}
8 H u0 p0 c0 {7,S}
"""),
    E0 = (-327.033,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,253,525,597,667,842,1178,1324,3025,407.5,1350,352.5,180,180,1884.04],'cm^-1')),
        HinderedRotor(inertia=(0.42436,'amu*angstrom^2'), symmetry=1, barrier=(9.75687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.423509,'amu*angstrom^2'), symmetry=1, barrier=(9.7373,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.023,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.953032,0.0886844,-0.000205361,2.1527e-07,-8.06298e-11,-39244,21.5344], Tmin=(100,'K'), Tmax=(873.814,'K')), NASAPolynomial(coeffs=[1.50426,0.033166,-1.90853e-05,3.74671e-09,-2.55453e-13,-37317.1,30.5264], Tmin=(873.814,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-327.033,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCOJ) + radical(Csj(Cs-F1sF1sO2s)(O2s-H)(H))"""),
)

species(
    label = '[O][C](F)C(F)OF(483)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {4,S}
4 O u0 p2 c0 {3,S} {6,S}
5 O u1 p2 c0 {7,S}
6 C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
7 C u1 p0 c0 {2,S} {5,S} {6,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-298.509,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,487,638,688,1119,1325,1387,3149,395,473,707,1436,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.01269,'amu*angstrom^2'), symmetry=1, barrier=(23.2837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01163,'amu*angstrom^2'), symmetry=1, barrier=(23.2593,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.023,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14979,0.0650998,-9.28953e-05,6.27312e-08,-1.63412e-11,-35801.7,20.8471], Tmin=(100,'K'), Tmax=(947.517,'K')), NASAPolynomial(coeffs=[14.2553,0.0097737,-5.3089e-06,1.10571e-09,-8.13777e-14,-38285.3,-41.673], Tmin=(947.517,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-298.509,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sH)) + radical(Csj(Cs-F1sO2sH)(F1s)(O2s-H))"""),
)

species(
    label = '[O]C(F)[C](F)OF(484)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {4,S}
4 O u0 p2 c0 {3,S} {7,S}
5 O u1 p2 c0 {6,S}
6 C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
7 C u1 p0 c0 {2,S} {4,S} {6,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-301.03,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,391,562,707,872,1109,1210,1289,3137,395,473,707,1436,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.48813,'amu*angstrom^2'), symmetry=1, barrier=(11.2231,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.44975,'amu*angstrom^2'), symmetry=1, barrier=(33.3326,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.023,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25311,0.0660976,-0.000112216,9.73066e-08,-3.3139e-11,-36112,22.1542], Tmin=(100,'K'), Tmax=(781.852,'K')), NASAPolynomial(coeffs=[9.84582,0.0173536,-9.52259e-06,1.91753e-09,-1.3591e-13,-37309.5,-16.2511], Tmin=(781.852,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-301.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sH)) + radical(CsCsF1sO2s)"""),
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
    E0 = (-296.021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (115.444,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (115.079,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-287.736,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-232.62,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-110.893,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-250.841,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-207.871,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-150.053,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-284.512,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (42.8185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (115.324,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-119.138,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-170.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-220.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (86.2938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (99.1251,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (119.366,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (93.9497,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(F)C([O])(F)F(429)'],
    products = ['CF2O(48)', 'CHFO(46)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O(6)', '[O]C(F)(F)[CH]F(474)'],
    products = ['[O]C(F)C([O])(F)F(429)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_sec_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['O(6)', '[O]C(F)[C](F)F(475)'],
    products = ['[O]C(F)C([O])(F)F(429)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_ter_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]C(F)C([O])(F)F(429)'],
    products = ['FC1OOC1(F)F(433)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]C(F)C([O])(F)F(429)'],
    products = ['O=C(F)C(O)(F)F(385)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['F(37)', '[O]C(F)(F)C=O(476)'],
    products = ['[O]C(F)C([O])(F)F(429)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4.08261e+20,'m^3/(mol*s)'), n=-5.07836, Ea=(16.2418,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O][C](F)F(268)', 'CHFO(46)'],
    products = ['[O]C(F)C([O])(F)F(429)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.310487,'m^3/(mol*s)'), n=1.79149, Ea=(73.5311,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.27712170538623165, var=2.20453978455454, Tref=1000.0, N=288, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(5)', '[O]C(F)(F)C(=O)F(290)'],
    products = ['[O]C(F)C([O])(F)F(429)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(15.9,'m^3/(mol*s)'), n=1.84, Ea=(25.9398,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_N-Sp-6R!H=5R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_2COS->O_Ext-1COS-R_Ext-1COS-R_Ext-5R!H-R_N-Sp-6R!H=5R!H"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F(37)', '[O]C(F)C(=O)F(477)'],
    products = ['[O]C(F)C([O])(F)F(429)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.08261e+20,'m^3/(mol*s)'), n=-5.07836, Ea=(27.5882,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CF2O(48)', '[O][CH]F(135)'],
    products = ['[O]C(F)C([O])(F)F(429)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.4291e-06,'m^3/(mol*s)'), n=3.20779, Ea=(67.894,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing_Ext-1R!H-R_N-2R!H->C_N-5R!H-u1',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing_Ext-1R!H-R_N-2R!H->C_N-5R!H-u1"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O][C](F)F(268)', '[O][CH]F(135)'],
    products = ['[O]C(F)C([O])(F)F(429)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F2(77)', '[O][CH]C(=O)F(478)'],
    products = ['[O]C(F)C([O])(F)F(429)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(20.6301,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction13',
    reactants = ['HF(38)', '[O][C](F)C(=O)F(291)'],
    products = ['[O]C(F)C([O])(F)F(429)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(8.28222,'m^3/(mol*s)'), n=1.29695, Ea=(280.241,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]C(F)C([O])(F)F(429)'],
    products = ['[O]C(F)(F)[C](O)F(479)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.66008e+10,'s^-1'), n=0.776642, Ea=(125.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;Cs_H_out_noH] + [R2H_S;O_rad_out;Cs_H_out] for rate rule [R2H_S;O_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O]C(F)C([O])(F)F(429)'],
    products = ['[O][C](F)C(O)(F)F(480)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]C(F)(F)[CH]OF(481)'],
    products = ['[O]C(F)C([O])(F)F(429)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(79.8321,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O][CH]C(F)(F)OF(482)'],
    products = ['[O]C(F)C([O])(F)F(429)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(108.183,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O][C](F)C(F)OF(483)'],
    products = ['[O]C(F)C([O])(F)F(429)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(99.9001,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O]C(F)[C](F)OF(484)'],
    products = ['[O]C(F)C([O])(F)F(429)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(77.0037,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

network(
    label = 'PDepNetwork #158',
    isomers = [
        '[O]C(F)C([O])(F)F(429)',
    ],
    reactants = [
        ('CF2O(48)', 'CHFO(46)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #158',
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

