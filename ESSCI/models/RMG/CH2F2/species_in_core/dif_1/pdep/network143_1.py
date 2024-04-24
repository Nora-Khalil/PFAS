species(
    label = 'F[C](F)OC(F)[C](F)F(340)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  O u0 p2 c0 {7,S} {9,S}
7  C u0 p0 c0 {1,S} {6,S} {8,S} {10,S}
8  C u1 p0 c0 {2,S} {3,S} {7,S}
9  C u1 p0 c0 {4,S} {5,S} {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-889.846,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([487,638,688,1119,1325,1387,3149,190,488,555,1236,1407,493,600,700,1144,1293,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.212653,'amu*angstrom^2'), symmetry=1, barrier=(4.88931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.211676,'amu*angstrom^2'), symmetry=1, barrier=(4.86685,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.807636,'amu*angstrom^2'), symmetry=1, barrier=(18.5691,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (148.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.292316,0.0917377,-0.000166909,1.51854e-07,-5.29813e-11,-106900,29.2184], Tmin=(100,'K'), Tmax=(830.101,'K')), NASAPolynomial(coeffs=[10.3791,0.0250888,-1.38691e-05,2.75968e-09,-1.92769e-13,-107953,-13.822], Tmin=(830.101,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-889.846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sO2sH)(F1s)(F1s)) + radical(Csj(F1s)(F1s)(O2s-Cs))"""),
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
    label = 'CF2(42)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u2 p0 c0 {1,S} {2,S}
"""),
    E0 = (33.7272,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([700.388,993.445,1307.31],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (50.0074,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2405.34,'J/mol'), sigma=(4.22615,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=375.71 K, Pc=72.31 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28591,0.0107608,-1.05382e-05,4.89881e-09,-8.86384e-13,4216.69,13.1348], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[5.33121,0.00197748,-9.60248e-07,2.10704e-10,-1.5954e-14,3366.49,-2.56367], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(33.7272,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CF2(T)""", comment="""Thermo library: halogens"""),
)

species(
    label = '[O]C(F)[C](F)F(371)',
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
    label = 'F[CH]O[C](F)F(372)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {6,S}
4 O u0 p2 c0 {5,S} {6,S}
5 C u1 p0 c0 {1,S} {4,S} {7,S}
6 C u1 p0 c0 {2,S} {3,S} {4,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-460.67,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([580,1155,1237,1373,3147,493,600,700,1144,1293,180,180,2051.26],'cm^-1')),
        HinderedRotor(inertia=(0.385137,'amu*angstrom^2'), symmetry=1, barrier=(8.85506,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.385421,'amu*angstrom^2'), symmetry=1, barrier=(8.86158,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0239,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.02944,0.0516258,-0.000100303,9.81351e-08,-3.59049e-11,-55342.8,19.4656], Tmin=(100,'K'), Tmax=(849.823,'K')), NASAPolynomial(coeffs=[5.06211,0.0189691,-1.0215e-05,2.01039e-09,-1.39201e-13,-55194.5,9.23362], Tmin=(849.823,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-460.67,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(F1s)(O2s-Cs)(H)) + radical(Csj(F1s)(F1s)(O2s-Cs))"""),
)

species(
    label = 'FC1OC(F)(F)C1(F)F(342)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  O u0 p2 c0 {8,S} {9,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u0 p0 c0 {3,S} {6,S} {7,S} {10,S}
9  C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-1180.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (148.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.71256,0.0257385,7.45682e-05,-1.86832e-07,1.16264e-10,-141941,13.1956], Tmin=(10,'K'), Tmax=(589.954,'K')), NASAPolynomial(coeffs=[6.25413,0.0368921,-2.59639e-05,8.33009e-09,-9.99094e-13,-142735,-1.91138], Tmin=(589.954,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1180.2,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), label="""FC1OC(F)(F)C1(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'FC(F)=C(F)OC(F)F(373)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  O u0 p2 c0 {7,S} {8,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
8  C u0 p0 c0 {3,S} {6,S} {9,D}
9  C u0 p0 c0 {4,S} {5,S} {8,D}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-1110.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (148.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36493,0.0623433,-8.83008e-05,6.87972e-08,-2.22766e-11,-133606,13.7167], Tmin=(10,'K'), Tmax=(722.323,'K')), NASAPolynomial(coeffs=[9.20657,0.0299941,-2.11236e-05,6.79619e-09,-8.17738e-13,-134450,-12.5657], Tmin=(722.323,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1110.88,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), label="""FC(F)DC(F)OC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[CH][C](F)F(185)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 C u1 p0 c0 {1,S} {5,S} {6,S}
5 C u1 p0 c0 {2,S} {3,S} {4,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-282.469,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([334,575,1197,1424,3202,190,488,555,1236,1407,1399.55],'cm^-1')),
        HinderedRotor(inertia=(0.462103,'amu*angstrom^2'), symmetry=1, barrier=(10.6247,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.0245,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.7193,0.0310345,-4.23424e-05,3.17512e-08,-9.73807e-12,-33929.6,16.1924], Tmin=(100,'K'), Tmax=(789.997,'K')), NASAPolynomial(coeffs=[6.47706,0.0120087,-6.21892e-06,1.26852e-09,-9.20647e-14,-34523.4,-1.05094], Tmin=(789.997,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-282.469,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: Fluorine + radical(CsCsF1sH) + radical(CsCsF1sF1s)"""),
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
    label = 'F[C](F)OC=C(F)F(374)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {8,S}
4 F u0 p3 c0 {8,S}
5 O u0 p2 c0 {6,S} {8,S}
6 C u0 p0 c0 {5,S} {7,D} {9,S}
7 C u0 p0 c0 {1,S} {2,S} {6,D}
8 C u1 p0 c0 {3,S} {4,S} {5,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-726.812,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,493,600,700,1144,1293,180,180,404.655],'cm^-1')),
        HinderedRotor(inertia=(0.141077,'amu*angstrom^2'), symmetry=1, barrier=(3.24364,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.140425,'amu*angstrom^2'), symmetry=1, barrier=(3.22865,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (129.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.64848,0.0449749,0.000258617,-2.82936e-06,7.04004e-09,-87415.3,12.568], Tmin=(10,'K'), Tmax=(159.375,'K')), NASAPolynomial(coeffs=[5.78799,0.0332054,-2.52296e-05,8.6608e-09,-1.09978e-12,-87536.8,4.50474], Tmin=(159.375,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-726.812,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""F[C](F)OCDC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[O][C](F)F(178)',
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
    label = 'F[C](F)OC(F)=C(F)F(375)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {8,S}
3 F u0 p3 c0 {8,S}
4 F u0 p3 c0 {9,S}
5 F u0 p3 c0 {9,S}
6 O u0 p2 c0 {7,S} {9,S}
7 C u0 p0 c0 {1,S} {6,S} {8,D}
8 C u0 p0 c0 {2,S} {3,S} {7,D}
9 C u1 p0 c0 {4,S} {5,S} {6,S}
"""),
    E0 = (-898.019,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([326,540,652,719,1357,182,240,577,636,1210,1413,493,600,700,1144,1293,180,180,1431.87],'cm^-1')),
        HinderedRotor(inertia=(0.278373,'amu*angstrom^2'), symmetry=1, barrier=(6.40035,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.27701,'amu*angstrom^2'), symmetry=1, barrier=(6.369,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (147.023,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.35146,0.0692691,-0.000148581,1.91607e-07,-1.02765e-10,-108010,14.6241], Tmin=(10,'K'), Tmax=(449.722,'K')), NASAPolynomial(coeffs=[7.50946,0.0322861,-2.52272e-05,8.74762e-09,-1.11356e-12,-108384,-2.11306], Tmin=(449.722,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-898.019,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""F[C](F)OC(F)DC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28591,0.0107608,-1.05382e-05,4.89881e-09,-8.86384e-13,-23521.2,13.1348], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[5.33121,0.00197748,-9.60248e-07,2.10704e-10,-1.5954e-14,-24371.4,-2.56367], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-196.899,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CF2""", comment="""Thermo library: Fluorine"""),
)

species(
    label = 'F[C](F)O[C](F)C(F)F(376)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  O u0 p2 c0 {8,S} {9,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {10,S}
8  C u1 p0 c0 {3,S} {6,S} {7,S}
9  C u1 p0 c0 {4,S} {5,S} {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-894.278,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (148.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.44169,0.0902753,-0.000170139,1.60155e-07,-5.71344e-11,-107440,29.5696], Tmin=(100,'K'), Tmax=(841.662,'K')), NASAPolynomial(coeffs=[8.41917,0.0279017,-1.53836e-05,3.04701e-09,-2.11825e-13,-107917,-2.39504], Tmin=(841.662,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-894.278,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sF1sH)(F1s)(O2s-Cs)) + radical(Csj(F1s)(F1s)(O2s-Cs))"""),
)

species(
    label = 'F[C](F)[C](F)OC(F)F(377)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  O u0 p2 c0 {7,S} {8,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {10,S}
8  C u1 p0 c0 {3,S} {6,S} {9,S}
9  C u1 p0 c0 {4,S} {5,S} {8,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-898.214,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (148.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.441467,0.0895925,-0.000166491,1.5563e-07,-5.54624e-11,-107913,29.863], Tmin=(100,'K'), Tmax=(833.86,'K')), NASAPolynomial(coeffs=[8.77698,0.0275696,-1.5277e-05,3.04066e-09,-2.12344e-13,-108537,-4.24248], Tmin=(833.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-898.214,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sF1sH)(F1s)(O2s-Cs)) + radical(Csj(Cs-F1sO2sH)(F1s)(F1s))"""),
)

species(
    label = 'F[C](F)O[CH]C(F)(F)F(378)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  O u0 p2 c0 {8,S} {9,S}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
8  C u1 p0 c0 {6,S} {7,S} {10,S}
9  C u1 p0 c0 {4,S} {5,S} {6,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-948.881,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (148.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.731995,0.0771861,-0.000111984,8.14898e-08,-2.3383e-11,-114011,28.3216], Tmin=(100,'K'), Tmax=(854.547,'K')), NASAPolynomial(coeffs=[13.1239,0.0191842,-1.01774e-05,2.06979e-09,-1.49515e-13,-116129,-29.515], Tmin=(854.547,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-948.881,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sF1sF1s)(O2s-Cs)(H)) + radical(Csj(F1s)(F1s)(O2s-Cs))"""),
)

species(
    label = 'F[C](F)[CH]OC(F)(F)F(379)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  O u0 p2 c0 {7,S} {8,S}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
8  C u1 p0 c0 {6,S} {9,S} {10,S}
9  C u1 p0 c0 {4,S} {5,S} {8,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-952.469,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (148.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.706458,0.0791698,-0.000129723,1.11735e-07,-3.82827e-11,-114443,29.3091], Tmin=(100,'K'), Tmax=(756.959,'K')), NASAPolynomial(coeffs=[10.5225,0.0237422,-1.28391e-05,2.5859e-09,-1.83986e-13,-115828,-14.6411], Tmin=(756.959,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-952.469,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sF1sH)(O2s-Cs)(H)) + radical(Csj(Cs-O2sHH)(F1s)(F1s))"""),
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
    E0 = (-395.915,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (81.7277,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (66.9881,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-387.631,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-332.515,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-352.252,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-111.164,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-253.672,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-189.556,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-44.0152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-163.638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-142.813,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-262.416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-254.891,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-236.403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-199.085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['F[C](F)OC(F)[C](F)F(340)'],
    products = ['CF2O(48)', 'CHFCF2(54)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CF2(42)', '[O]C(F)[C](F)F(371)'],
    products = ['F[C](F)OC(F)[C](F)F(340)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['CF2(42)', 'F[CH]O[C](F)F(372)'],
    products = ['F[C](F)OC(F)[C](F)F(340)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_sec_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['F[C](F)OC(F)[C](F)F(340)'],
    products = ['FC1OC(F)(F)C1(F)F(342)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_noH;Cpri_rad_out_noH]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['F[C](F)OC(F)[C](F)F(340)'],
    products = ['FC(F)=C(F)OC(F)F(373)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CF2O(48)', 'F[CH][C](F)F(185)'],
    products = ['F[C](F)OC(F)[C](F)F(340)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(9.07578e-06,'m^3/(mol*s)'), n=3.04336, Ea=(86.9875,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3757377757886876, var=2.242054186761003, Tref=1000.0, N=1042, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F(37)', 'F[C](F)OC=C(F)F(374)'],
    products = ['F[C](F)OC(F)[C](F)F(340)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(48.8257,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O][C](F)F(178)', 'CHFCF2(54)'],
    products = ['F[C](F)OC(F)[C](F)F(340)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(13.699,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(5)', 'F[C](F)OC(F)=C(F)F(375)'],
    products = ['F[C](F)OC(F)[C](F)F(340)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(70.704,'m^3/(mol*s)'), n=1.71182, Ea=(2.72724,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_N-4R!H->C_Ext-1COS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_N-4R!H->C_Ext-1COS-R"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O][C](F)F(178)', 'F[CH][C](F)F(185)'],
    products = ['F[C](F)OC(F)[C](F)F(340)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CF2(89)', 'F[CH]O[C](F)F(372)'],
    products = ['F[C](F)OC(F)[C](F)F(340)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(786.723,'m^3/(mol*s)'), n=1.25031, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CF2(89)', '[O]C(F)[C](F)F(371)'],
    products = ['F[C](F)OC(F)[C](F)F(340)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(128827,'m^3/(mol*s)'), n=0.469398, Ea=(6.08571,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.009898522871762284, var=0.3372166703721302, Tref=1000.0, N=6, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F[C](F)OC(F)[C](F)F(340)'],
    products = ['F[C](F)O[C](F)C(F)F(376)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.70914e+10,'s^-1'), n=0.735063, Ea=(133.5,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_noH;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['F[C](F)OC(F)[C](F)F(340)'],
    products = ['F[C](F)[C](F)OC(F)F(377)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.95012e+07,'s^-1'), n=1.37671, Ea=(141.025,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_noH;XH_out] for rate rule [R3H_SS_O;C_rad_out_noH;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F[C](F)OC(F)[C](F)F(340)'],
    products = ['F[C](F)O[CH]C(F)(F)F(378)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(159.512,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F[C](F)OC(F)[C](F)F(340)'],
    products = ['F[C](F)[CH]OC(F)(F)F(379)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(196.83,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

network(
    label = 'PDepNetwork #143',
    isomers = [
        'F[C](F)OC(F)[C](F)F(340)',
    ],
    reactants = [
        ('CF2O(48)', 'CHFCF2(54)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #143',
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

