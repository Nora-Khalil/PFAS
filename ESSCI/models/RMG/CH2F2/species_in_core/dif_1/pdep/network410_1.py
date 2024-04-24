species(
    label = '[CH2]C(F)(F)C(F)[C](F)F(1035)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {3,S} {6,S} {9,S} {10,S}
8  C u1 p0 c0 {6,S} {11,S} {12,S}
9  C u1 p0 c0 {4,S} {5,S} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-757.743,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,259,529,569,1128,1321,1390,3140,3000,3100,440,815,1455,1000,190,488,555,1236,1407,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.203576,'amu*angstrom^2'), symmetry=1, barrier=(4.68061,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.34677,'amu*angstrom^2'), symmetry=1, barrier=(7.97293,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203825,'amu*angstrom^2'), symmetry=1, barrier=(4.68634,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (146.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0278302,0.101325,-0.000183911,1.71829e-07,-6.14353e-11,-91002.6,29.436], Tmin=(100,'K'), Tmax=(834.073,'K')), NASAPolynomial(coeffs=[8.73366,0.0343111,-1.84398e-05,3.63908e-09,-2.534e-13,-91594.7,-6.03132], Tmin=(834.073,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-757.743,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsF1sF1s)(H)(H)) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
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
    label = 'CH2CF2(56)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 C u0 p0 c0 {4,D} {5,S} {6,S}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-349.675,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(64.0125,'amu')),
        NonlinearRotor(inertia=([45.7027,48.2614,93.9642],'amu*angstrom^2'), symmetry=2),
        HarmonicOscillator(frequencies=([437.293,557.015,653.832,726.079,816.319,956,966.438,1345.56,1413.22,1792.31,3202.97,3303.55],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2091.09,'J/mol'), sigma=(4.442,'angstroms'), dipoleMoment=(1.4,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.544641,0.0361013,-4.24522e-05,2.68678e-08,-7.0049e-12,-41578.9,25.9139], Tmin=(298,'K'), Tmax=(1100,'K')), NASAPolynomial(coeffs=[7.52692,0.00816727,-3.24648e-06,5.85466e-10,-3.90185e-14,-43575.5,-14.4929], Tmin=(1100,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-349.675,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(133.032,'J/mol/K'), label="""CH2CF2""", comment="""Thermo library: Fluorine"""),
)

species(
    label = 'F[C](F)CC(F)[C](F)F(1036)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {7,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {1,S} {6,S} {9,S} {12,S}
8  C u1 p0 c0 {2,S} {3,S} {6,S}
9  C u1 p0 c0 {4,S} {5,S} {7,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-764.719,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,259,529,569,1128,1321,1390,3140,146,234,414,562,504,606,1176,1296,1354,1460,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.244406,'amu*angstrom^2'), symmetry=1, barrier=(5.61938,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.245165,'amu*angstrom^2'), symmetry=1, barrier=(5.63682,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244788,'amu*angstrom^2'), symmetry=1, barrier=(5.62816,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (146.058,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2957.56,'J/mol'), sigma=(4.97116,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=461.96 K, Pc=54.63 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0473045,0.0987328,-0.000176616,1.61912e-07,-5.66573e-11,-91843.5,30.1928], Tmin=(100,'K'), Tmax=(848.639,'K')), NASAPolynomial(coeffs=[9.25008,0.0317055,-1.63395e-05,3.16361e-09,-2.1754e-13,-92553.8,-7.67715], Tmin=(848.639,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-764.719,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsHH)(F1s)(F1s)) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
)

species(
    label = '[CH2]C(F)(F)C(F)(F)[CH]F(1033)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {7,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {3,S} {4,S} {6,S} {9,S}
8  C u1 p0 c0 {6,S} {10,S} {11,S}
9  C u1 p0 c0 {5,S} {7,S} {12,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-766.86,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.127825,0.101871,-0.000179645,1.62811e-07,-5.69258e-11,-92094,29.2006], Tmin=(100,'K'), Tmax=(828.867,'K')), NASAPolynomial(coeffs=[10.4315,0.0313668,-1.66802e-05,3.28546e-09,-2.28818e-13,-93173,-15.7096], Tmin=(828.867,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-766.86,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsF1sF1s)(H)(H)) + radical(Csj(Cs-CsF1sF1s)(F1s)(H))"""),
)

species(
    label = 'CH2(T)(68)',
    structure = adjacencyList("""multiplicity 3
1 C u2 p0 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
    E0 = (381.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1066.91,2790.97,3622.38],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.01192,-0.000154978,3.26297e-06,-2.40422e-09,5.69496e-13,45867.7,0.533201], Tmin=(100,'K'), Tmax=(1104.65,'K')), NASAPolynomial(coeffs=[3.14983,0.00296674,-9.76056e-07,1.54115e-10,-9.50339e-15,46058.1,4.77808], Tmin=(1104.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'F[C](F)C(F)[C](F)F(272)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {7,S}
4 F u0 p3 c0 {8,S}
5 F u0 p3 c0 {8,S}
6 C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
7 C u1 p0 c0 {2,S} {3,S} {6,S}
8 C u1 p0 c0 {4,S} {5,S} {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-711.656,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,146,234,414,562,504,606,1176,1296,1354,1460,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.130631,'amu*angstrom^2'), symmetry=1, barrier=(3.00347,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130993,'amu*angstrom^2'), symmetry=1, barrier=(3.01179,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (132.032,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03372,0.0738643,-0.000133097,1.21982e-07,-4.27798e-11,-85493.9,26.1783], Tmin=(100,'K'), Tmax=(838.788,'K')), NASAPolynomial(coeffs=[8.34888,0.0225509,-1.19533e-05,2.34634e-09,-1.62743e-13,-86143.1,-4.38197], Tmin=(838.788,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-711.656,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsF1sH)(F1s)(F1s)) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
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
    label = '[CH2]C(F)(F)[CH]F(1055)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 F u0 p3 c0 {6,S}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u1 p0 c0 {4,S} {7,S} {8,S}
6 C u1 p0 c0 {3,S} {4,S} {9,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-330.121,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,3000,3100,440,815,1455,1000,334,575,1197,1424,3202,1541.11],'cm^-1')),
        HinderedRotor(inertia=(0.283803,'amu*angstrom^2'), symmetry=1, barrier=(6.52519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.283994,'amu*angstrom^2'), symmetry=1, barrier=(6.52958,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69302,0.0567507,-9.31345e-05,8.46039e-08,-3.02362e-11,-39627.2,21.0625], Tmin=(100,'K'), Tmax=(810.485,'K')), NASAPolynomial(coeffs=[6.60538,0.0225535,-1.14235e-05,2.24022e-09,-1.56745e-13,-40096.6,0.412073], Tmin=(810.485,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-330.121,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-CsF1sF1s)(H)(H)) + radical(Csj(Cs-CsF1sF1s)(F1s)(H))"""),
)

species(
    label = 'FC1C(F)(F)CC1(F)F(1078)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {3,S} {6,S} {9,S} {10,S}
8  C u0 p0 c0 {6,S} {9,S} {11,S} {12,S}
9  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-1025.09,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.7046,0.0255133,0.000100803,-2.37145e-07,1.48195e-10,-123285,12.5947], Tmin=(10,'K'), Tmax=(566.91,'K')), NASAPolynomial(coeffs=[4.98583,0.0459474,-3.12505e-05,9.85458e-09,-1.17108e-12,-123904,2.96333], Tmin=(566.91,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-1025.09,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), label="""FC1C(F)(F)CC1(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CC(F)(F)C(F)=C(F)F(1079)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {6,S} {10,S} {11,S} {12,S}
8  C u0 p0 c0 {3,S} {6,S} {9,D}
9  C u0 p0 c0 {4,S} {5,S} {8,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-993.428,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36003,0.0555655,-3.52303e-05,-2.09724e-08,2.87479e-11,-119477,14.0697], Tmin=(10,'K'), Tmax=(582.476,'K')), NASAPolynomial(coeffs=[7.92099,0.0379367,-2.50932e-05,7.78265e-09,-9.15279e-13,-120241,-7.4631], Tmin=(582.476,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-993.428,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), label="""CC(F)(F)C(F)DC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'C=C(F)C(F)[C](F)F(1080)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {7,S}
5  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {2,S} {5,S} {8,D}
7  C u1 p0 c0 {3,S} {4,S} {5,S}
8  C u0 p0 c0 {6,D} {10,S} {11,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-605.066,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,323,467,575,827,1418,190,488,555,1236,1407,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.310086,'amu*angstrom^2'), symmetry=1, barrier=(7.1295,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.506856,'amu*angstrom^2'), symmetry=1, barrier=(11.6536,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.41677,0.0579465,-7.76503e-05,6.50706e-08,-2.40911e-11,-72772.4,13.8092], Tmin=(10,'K'), Tmax=(618.545,'K')), NASAPolynomial(coeffs=[6.80705,0.0360223,-2.44828e-05,7.76666e-09,-9.30228e-13,-73191.8,-0.918315], Tmin=(618.545,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-605.066,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), label="""CDC(F)C(F)[C](F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[CH2]C(F)(F)C=C(F)F(1081)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {5,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {8,S}
5  C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,D} {9,S}
7  C u1 p0 c0 {5,S} {10,S} {11,S}
8  C u0 p0 c0 {3,S} {4,S} {6,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-630.415,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,182,240,577,636,1210,1413],'cm^-1')),
        HinderedRotor(inertia=(0.0753768,'amu*angstrom^2'), symmetry=1, barrier=(1.73306,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0749441,'amu*angstrom^2'), symmetry=1, barrier=(1.72311,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (127.06,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.56047,0.0419245,-1.63376e-05,-1.62633e-08,1.22541e-11,-75819.4,16.3194], Tmin=(10,'K'), Tmax=(803.079,'K')), NASAPolynomial(coeffs=[8.30406,0.0308751,-1.9192e-05,5.60834e-09,-6.25668e-13,-76986.9,-8.0506], Tmin=(803.079,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-630.415,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), label="""[CH2]C(F)(F)CDC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[CH2][C](F)F(187)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 C u1 p0 c0 {4,S} {5,S} {6,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-107.299,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,190,488,555,1236,1407],'cm^-1')),
        HinderedRotor(inertia=(0.00459569,'amu*angstrom^2'), symmetry=1, barrier=(7.47371,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.09702,0.01871,-1.30857e-05,4.49229e-09,-6.17295e-13,-12871.8,14.3472], Tmin=(100,'K'), Tmax=(1681.06,'K')), NASAPolynomial(coeffs=[7.67292,0.00782184,-3.37022e-06,6.39393e-10,-4.43074e-14,-14410.3,-10.1057], Tmin=(1681.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-107.299,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: Fluorine + radical(Cs_P) + radical(CsCsF1sF1s)"""),
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
    label = '[CH2]C(F)(F)C(F)=C(F)F(1082)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {3,S} {6,S} {9,D}
8  C u1 p0 c0 {6,S} {10,S} {11,S}
9  C u0 p0 c0 {4,S} {5,S} {7,D}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-783.14,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,323,467,575,827,1418,3000,3100,440,815,1455,1000,182,240,577,636,1210,1413],'cm^-1')),
        HinderedRotor(inertia=(0.21578,'amu*angstrom^2'), symmetry=1, barrier=(4.9612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216309,'amu*angstrom^2'), symmetry=1, barrier=(4.97337,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (145.051,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.23982,0.0725785,-0.000122692,1.20807e-07,-4.94594e-11,-94188.9,14.0031], Tmin=(10,'K'), Tmax=(582.831,'K')), NASAPolynomial(coeffs=[8.83967,0.0341462,-2.37803e-05,7.6675e-09,-9.28836e-13,-94841.6,-9.98979], Tmin=(582.831,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-783.14,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), label="""[CH2]C(F)(F)C(F)DC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[CH2]C(F)=C[C](F)F(1083)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  C u0 p0 c0 {5,D} {6,S} {8,S}
5  C u0 p0 c0 {1,S} {4,D} {7,S}
6  C u1 p0 c0 {2,S} {3,S} {4,S}
7  C u1 p0 c0 {5,S} {9,S} {10,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-337.068,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,271,519,563,612,1379,161,297,490,584,780,1358,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.589615,'amu*angstrom^2'), symmetry=1, barrier=(13.5564,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152294,'amu*angstrom^2'), symmetry=1, barrier=(56.4778,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.062,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52101,0.0580989,-6.749e-05,4.15435e-08,-1.03976e-11,-40453.7,20.877], Tmin=(100,'K'), Tmax=(962.156,'K')), NASAPolynomial(coeffs=[10.3009,0.0215977,-1.0584e-05,2.11362e-09,-1.52272e-13,-42143.2,-21.1419], Tmin=(962.156,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-337.068,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cd-F1sCd)(H)(H)) + radical(Csj(Cd-CdH)(F1s)(F1s))"""),
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
    label = '[CH2]C(F)=C(F)[C](F)F(1084)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {7,S}
5  C u0 p0 c0 {1,S} {6,D} {7,S}
6  C u0 p0 c0 {2,S} {5,D} {8,S}
7  C u1 p0 c0 {3,S} {4,S} {5,S}
8  C u1 p0 c0 {6,S} {9,S} {10,S}
9  H u0 p0 c0 {8,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-485.921,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([206,336,431,607,515,611,528,696,1312,1446,161,297,490,584,780,1358,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.0969427,'amu*angstrom^2'), symmetry=1, barrier=(2.2289,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.56412,'amu*angstrom^2'), symmetry=1, barrier=(58.9542,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.677246,0.0800675,-0.000132147,1.14743e-07,-3.92523e-11,-58329.9,24.8944], Tmin=(100,'K'), Tmax=(797.466,'K')), NASAPolynomial(coeffs=[10.1795,0.0245038,-1.2772e-05,2.52235e-09,-1.77015e-13,-59594.2,-17.2226], Tmin=(797.466,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-485.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cd-F1sCd)(H)(H)) + radical(Csj(Cd-F1sCd)(F1s)(F1s))"""),
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
    label = '[CH2]C(F)(F)[C](F)C(F)F(1085)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {7,S}
5  F u0 p3 c0 {8,S}
6  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
7  C u0 p0 c0 {3,S} {4,S} {8,S} {10,S}
8  C u1 p0 c0 {5,S} {6,S} {7,S}
9  C u1 p0 c0 {6,S} {11,S} {12,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-766.094,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0299843,0.102984,-0.00019228,1.83068e-07,-6.6076e-11,-92008.5,30.2048], Tmin=(100,'K'), Tmax=(841.843,'K')), NASAPolynomial(coeffs=[7.74486,0.036014,-1.94498e-05,3.83183e-09,-2.66004e-13,-92253.6,0.353831], Tmin=(841.843,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-766.094,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(F1s)) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
)

species(
    label = 'CC(F)(F)[C](F)[C](F)F(1086)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {6,S} {10,S} {11,S} {12,S}
8  C u1 p0 c0 {3,S} {6,S} {9,S}
9  C u1 p0 c0 {4,S} {5,S} {8,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-776.785,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0968556,0.10424,-0.000207341,2.06731e-07,-7.65165e-11,-93302.8,29.4171], Tmin=(100,'K'), Tmax=(853.682,'K')), NASAPolynomial(coeffs=[4.54218,0.0412039,-2.24169e-05,4.40244e-09,-3.03994e-13,-92523.8,17.6821], Tmin=(853.682,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-776.785,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(F1s)) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
)

species(
    label = 'FC[C](F)C(F)[C](F)F(1087)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
7  C u0 p0 c0 {2,S} {8,S} {11,S} {12,S}
8  C u1 p0 c0 {3,S} {6,S} {7,S}
9  C u1 p0 c0 {4,S} {5,S} {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
"""),
    E0 = (-729.613,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,551,1088,1226,1380,1420,1481,3057,3119,212,367,445,1450,190,488,555,1236,1407,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.37334,'amu*angstrom^2'), symmetry=1, barrier=(31.5758,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.220582,'amu*angstrom^2'), symmetry=1, barrier=(5.07161,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.219673,'amu*angstrom^2'), symmetry=1, barrier=(5.05072,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (146.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.195389,0.0938847,-0.000161241,1.47058e-07,-5.22723e-11,-87625,30.5808], Tmin=(100,'K'), Tmax=(814.586,'K')), NASAPolynomial(coeffs=[9.09058,0.0330016,-1.74504e-05,3.45158e-09,-2.41839e-13,-88503.4,-7.0057], Tmin=(814.586,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-729.613,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sCsH)(Cs-F1sHH)(F1s)) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
)

species(
    label = '[CH2][C](F)C(F)C(F)(F)F(1088)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {7,S}
5  F u0 p3 c0 {8,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
8  C u1 p0 c0 {5,S} {6,S} {9,S}
9  C u1 p0 c0 {8,S} {11,S} {12,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-787.815,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.332635,0.0886042,-0.000143631,1.25395e-07,-4.32853e-11,-94627.8,29.5664], Tmin=(100,'K'), Tmax=(801.009,'K')), NASAPolynomial(coeffs=[10.0314,0.0301421,-1.53709e-05,3.01503e-09,-2.11007e-13,-95859.8,-13.0642], Tmin=(801.009,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-787.815,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sCsH)(Cs-HHH)(F1s)) + radical(Csj(Cs-F1sCsH)(H)(H))"""),
)

species(
    label = 'FCC(F)(F)[CH][C](F)F(1089)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  F u0 p3 c0 {9,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {3,S} {6,S} {10,S} {11,S}
8  C u1 p0 c0 {6,S} {9,S} {12,S}
9  C u1 p0 c0 {4,S} {5,S} {8,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {8,S}
"""),
    E0 = (-763.182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.107879,0.0977903,-0.000175158,1.64159e-07,-5.91987e-11,-91661.2,29.793], Tmin=(100,'K'), Tmax=(825.835,'K')), NASAPolynomial(coeffs=[8.21765,0.0351355,-1.88988e-05,3.74408e-09,-2.61818e-13,-92203.6,-2.95421], Tmin=(825.835,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-763.182,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sH)(H)) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = '[CH2]C(F)(F)[CH]C(F)(F)F(1090)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {7,S}
5  F u0 p3 c0 {7,S}
6  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
7  C u0 p0 c0 {3,S} {4,S} {5,S} {8,S}
8  C u1 p0 c0 {6,S} {7,S} {10,S}
9  C u1 p0 c0 {6,S} {11,S} {12,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {9,S}
12 H u0 p0 c0 {9,S}
"""),
    E0 = (-822.133,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (146.058,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.725272,0.0791518,-0.000110802,8.60111e-08,-2.74062e-11,-98768.3,29.4229], Tmin=(100,'K'), Tmax=(760.961,'K')), NASAPolynomial(coeffs=[9.83248,0.0312768,-1.64255e-05,3.32459e-09,-2.39497e-13,-100154,-12.0257], Tmin=(760.961,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-822.133,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sF1sCs)(Cs-F1sF1sF1s)(H)) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
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
    E0 = (-327.723,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-170.404,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-167.788,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (99.7342,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (133.626,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-319.438,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-264.323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-50.4416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-195.209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-68.0397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-179.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-140.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (40.2521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (102.213,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-105.834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-97.0003,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-194.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-233.718,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-159.195,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-115.61,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-103.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-177.077,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(F)(F)C(F)[C](F)F(1035)'],
    products = ['CHFCF2(54)', 'CH2CF2(56)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(F)(F)C(F)[C](F)F(1035)'],
    products = ['F[C](F)CC(F)[C](F)F(1036)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(F)(F)C(F)[C](F)F(1035)'],
    products = ['[CH2]C(F)(F)C(F)(F)[CH]F(1033)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH2(T)(68)', 'F[C](F)C(F)[C](F)F(272)'],
    products = ['[CH2]C(F)(F)C(F)[C](F)F(1035)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4.0899e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['CF2(42)', '[CH2]C(F)(F)[CH]F(1055)'],
    products = ['[CH2]C(F)(F)C(F)[C](F)F(1035)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_sec_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C(F)(F)C(F)[C](F)F(1035)'],
    products = ['FC1C(F)(F)CC1(F)F(1078)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_noH]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C(F)(F)C(F)[C](F)F(1035)'],
    products = ['CC(F)(F)C(F)=C(F)F(1079)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F(37)', 'C=C(F)C(F)[C](F)F(1080)'],
    products = ['[CH2]C(F)(F)C(F)[C](F)F(1035)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(51.7125,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction9',
    reactants = ['F[CH][C](F)F(185)', 'CH2CF2(56)'],
    products = ['[CH2]C(F)(F)C(F)[C](F)F(1035)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(6.91486,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F(37)', '[CH2]C(F)(F)C=C(F)F(1081)'],
    products = ['[CH2]C(F)(F)C(F)[C](F)F(1035)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(59.4636,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction11',
    reactants = ['CHFCF2(54)', '[CH2][C](F)F(187)'],
    products = ['[CH2]C(F)(F)C(F)[C](F)F(1035)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(0.000143481,'m^3/(mol*s)'), n=2.73966, Ea=(4.09045,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24334233861756036, var=1.6451213851449848, Tref=1000.0, N=114, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing_4R!H->C_Sp-2R!H=1R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing_4R!H->C_Sp-2R!H=1R!H"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(5)', '[CH2]C(F)(F)C(F)=C(F)F(1082)'],
    products = ['[CH2]C(F)(F)C(F)[C](F)F(1035)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(56.8734,'m^3/(mol*s)'), n=1.75834, Ea=(0.884802,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.2648472359903826, var=0.02782886759889551, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_4R!H->C_Sp-4C-1COS_Ext-1COS-R_N-6R!H-inRing_N-7R!H-inRing_Ext-4C-R_N-Sp-8R!H#4C',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_4R!H->C_Sp-4C-1COS_Ext-1COS-R_N-6R!H-inRing_N-7R!H-inRing_Ext-4C-R_N-Sp-8R!H#4C"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F[CH][C](F)F(185)', '[CH2][C](F)F(187)'],
    products = ['[CH2]C(F)(F)C(F)[C](F)F(1035)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -7.3 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction14',
    reactants = ['F2(108)', '[CH2]C(F)=C[C](F)F(1083)'],
    products = ['[CH2]C(F)(F)C(F)[C](F)F(1035)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(18.0659,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction15',
    reactants = ['HF(38)', '[CH2]C(F)=C(F)[C](F)F(1084)'],
    products = ['[CH2]C(F)(F)C(F)[C](F)F(1035)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.14111,'m^3/(mol*s)'), n=1.29695, Ea=(231.181,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction16',
    reactants = ['CF2(89)', '[CH2]C(F)(F)[CH]F(1055)'],
    products = ['[CH2]C(F)(F)C(F)[C](F)F(1035)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(786.723,'m^3/(mol*s)'), n=1.25031, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(F)(F)C(F)[C](F)F(1035)'],
    products = ['[CH2]C(F)(F)[C](F)C(F)F(1085)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.70914e+10,'s^-1'), n=0.735063, Ea=(133.5,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_noH;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(F)(F)C(F)[C](F)F(1035)'],
    products = ['CC(F)(F)[C](F)[C](F)F(1086)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.4313e+06,'s^-1'), n=1.72764, Ea=(94.0044,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['FC[C](F)C(F)[C](F)F(1087)'],
    products = ['[CH2]C(F)(F)C(F)[C](F)F(1035)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.36622e+11,'s^-1'), n=0.48217, Ea=(140.398,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-2R!H-R',), comment="""Estimated from node R2F_Ext-2R!H-R"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(F)(F)C(F)[C](F)F(1035)'],
    products = ['[CH2][C](F)C(F)C(F)(F)F(1088)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(212.113,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(F)(F)C(F)[C](F)F(1035)'],
    products = ['FCC(F)(F)[CH][C](F)F(1089)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(224.688,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(F)(F)C(F)[C](F)F(1035)'],
    products = ['[CH2]C(F)(F)[CH]C(F)(F)F(1090)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(150.645,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

network(
    label = 'PDepNetwork #410',
    isomers = [
        '[CH2]C(F)(F)C(F)[C](F)F(1035)',
    ],
    reactants = [
        ('CHFCF2(54)', 'CH2CF2(56)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #410',
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

