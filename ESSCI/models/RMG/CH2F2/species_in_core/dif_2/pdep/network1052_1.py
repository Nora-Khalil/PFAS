species(
    label = 'F[C](F)CC1(F)[CH]O1(2714)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {7,S}
5  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
7  C u1 p0 c0 {4,S} {5,S} {11,S}
8  C u1 p0 c0 {2,S} {3,S} {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-340.951,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,1000,190,488,555,1236,1407,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43565,0.0614648,-7.49608e-05,5.2603e-08,-1.55024e-11,-40919.1,24.7533], Tmin=(100,'K'), Tmax=(813.47,'K')), NASAPolynomial(coeffs=[8.13613,0.0285171,-1.4207e-05,2.81309e-09,-2.00721e-13,-42009.2,-6.18931], Tmin=(813.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-340.951,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCFO) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(CsCsFFH) + ring(Cs(F)(C)-Cs-O2s) + radical(Csj(Cs-F1sO2sCs)(O2s-Cs)(H)_ring) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = 'FC1=CO1(390)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {4,S}
3 C u0 p0 c0 {2,S} {4,D} {5,S}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (74.4547,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,180,180,180,787.093,1742.58,1743.13,1743.63],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (60.0271,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3010.76,'J/mol'), sigma=(4.86933,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=470.27 K, Pc=59.17 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.102,0.0258723,-5.29145e-05,5.76369e-08,-2.28588e-11,8981.24,10.2596], Tmin=(100,'K'), Tmax=(846.703,'K')), NASAPolynomial(coeffs=[1.99159,0.016003,-8.6525e-06,1.7026e-09,-1.18158e-13,9711.08,18.6315], Tmin=(846.703,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(74.4547,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-CdsOsH) + group(CdCFO) + ring(oxirene)"""),
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
    E0 = (-361.616,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.10281,-0.0101072,0.000121983,-2.28108e-07,1.37933e-10,-43490.6,7.77929], Tmin=(10,'K'), Tmax=(534.293,'K')), NASAPolynomial(coeffs=[2.52167,0.0198841,-1.31824e-05,4.13929e-09,-4.93215e-13,-43580.8,11.9914], Tmin=(534.293,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-361.616,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""CDC(F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'F[C](F)CC1O[C]1F(2715)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {7,S}
5  C u0 p0 c0 {4,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {5,S} {8,S} {10,S} {11,S}
7  C u1 p0 c0 {1,S} {4,S} {5,S}
8  C u1 p0 c0 {2,S} {3,S} {6,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-325.984,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,2750,2850,1437.5,1250,1305,750,350,190,488,555,1236,1407,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3202.91,'J/mol'), sigma=(5.34861,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=500.29 K, Pc=47.5 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68179,0.0556779,-5.71586e-05,3.30253e-08,-8.15742e-12,-39127.2,24.6862], Tmin=(100,'K'), Tmax=(949.435,'K')), NASAPolynomial(coeffs=[8.14629,0.0284428,-1.41303e-05,2.81207e-09,-2.01841e-13,-40354.7,-6.16586], Tmin=(949.435,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-325.984,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCFHO) + group(CsCsFFH) + ring(O2s-Cs(F)-Cs(C)) + radical(Csj(Cs-O2sCsH)(F1s)(O2s-Cs)_ring) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = '[CH2]C(F)(F)C1(F)[CH]O1(2712)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  O u0 p2 c0 {5,S} {7,S}
5  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
6  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
7  C u1 p0 c0 {4,S} {5,S} {9,S}
8  C u1 p0 c0 {6,S} {10,S} {11,S}
9  H u0 p0 c0 {7,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-343.439,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,2950,1000,3000,3100,440,815,1455,1000,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1538,0.0646076,-7.21001e-05,4.16278e-08,-9.68112e-12,-41205.2,24.4281], Tmin=(100,'K'), Tmax=(1037.29,'K')), NASAPolynomial(coeffs=[12.1912,0.022046,-1.05537e-05,2.07254e-09,-1.47944e-13,-43495.1,-29.225], Tmin=(1037.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-343.439,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCFO) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Cs(F)(C)-Cs-O2s) + radical(Csj(Cs-F1sO2sCs)(O2s-Cs)(H)_ring) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
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
    label = '[CH2]C1(F)[CH]O1(2727)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {3,S} {4,S}
3 C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
4 C u1 p0 c0 {2,S} {3,S} {6,S}
5 C u1 p0 c0 {3,S} {7,S} {8,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (81.1941,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,3000,3100,440,815,1455,1000,180,180,759.86,762.322,762.438,762.64,763.149,764.966,765.227],'cm^-1')),
        HinderedRotor(inertia=(0.0106634,'amu*angstrom^2'), symmetry=1, barrier=(4.40251,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (74.0536,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28424,0.0319132,-1.91397e-05,7.08133e-10,2.14765e-12,9832.15,16.033], Tmin=(100,'K'), Tmax=(1072.76,'K')), NASAPolynomial(coeffs=[10.8874,0.0105825,-4.34238e-06,8.51854e-10,-6.23669e-14,7367.88,-28.9591], Tmin=(1072.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(81.1941,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sO2sCs)(O2s-Cs)(H)_ring) + radical(Csj(Cs-F1sO2sCs)(H)(H)_1914_ring)"""),
)

species(
    label = 'FC1(F)CC2(F)OC12(2728)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {6,S}
5  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
6  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
7  C u0 p0 c0 {5,S} {8,S} {10,S} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-657.354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05351,0.0551824,-4.2884e-05,1.50853e-08,-2.05568e-12,-78947.2,16.2846], Tmin=(100,'K'), Tmax=(1751.57,'K')), NASAPolynomial(coeffs=[19.5205,0.0130099,-6.7685e-06,1.3394e-09,-9.37354e-14,-85416.4,-83.159], Tmin=(1751.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-657.354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCFO) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCsCsFF) + polycyclic(s2_3_4_ane)"""),
)

species(
    label = 'FC(F)=CC1(F)CO1(2729)',
    structure = adjacencyList("""1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {6,S}
5  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
6  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
7  C u0 p0 c0 {5,S} {8,D} {11,S}
8  C u0 p0 c0 {2,S} {3,S} {7,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-629.179,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60991,0.0536878,-4.86461e-05,2.27488e-08,-4.37866e-12,-75587.8,21.5696], Tmin=(100,'K'), Tmax=(1218.24,'K')), NASAPolynomial(coeffs=[10.9151,0.0231347,-1.10264e-05,2.1617e-09,-1.53871e-13,-77854.9,-25.1593], Tmin=(1218.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-629.179,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCFO) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(CdCFF) + ring(Cs-Cs(F)-O2s)"""),
)

species(
    label = '[O]C=C(F)C[C](F)F(2730)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  O u1 p2 c0 {8,S}
5  C u0 p0 c0 {6,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {5,S} {8,D}
7  C u1 p0 c0 {2,S} {3,S} {5,S}
8  C u0 p0 c0 {4,S} {6,D} {11,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-474.627,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,323,467,575,827,1418,190,488,555,1236,1407,3010,987.5,1337.5,450,1655,180,536.333,536.37],'cm^-1')),
        HinderedRotor(inertia=(0.0894315,'amu*angstrom^2'), symmetry=1, barrier=(18.2676,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.794527,'amu*angstrom^2'), symmetry=1, barrier=(18.2677,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03591,0.0630067,-6.38134e-05,3.1763e-08,-6.25634e-12,-56975.9,26.352], Tmin=(100,'K'), Tmax=(1225.75,'K')), NASAPolynomial(coeffs=[14.8874,0.0178052,-8.49857e-06,1.67812e-09,-1.20344e-13,-60371.6,-43.2929], Tmin=(1225.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-474.627,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFFH) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Csj(Cs-CdHH)(F1s)(F1s))"""),
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
    label = 'F[C](F)CC1=CO1(2731)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {5,S} {6,S}
4  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
5  C u0 p0 c0 {3,S} {4,S} {6,D}
6  C u0 p0 c0 {3,S} {5,D} {10,S}
7  C u1 p0 c0 {1,S} {2,S} {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-13.8466,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,1000,190,488,555,1236,1407,363.003,363.003,363.003,363.004,363.004,1488.71,1488.71,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0851986,'amu*angstrom^2'), symmetry=1, barrier=(7.96679,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.247249,'amu*angstrom^2'), symmetry=1, barrier=(23.1195,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (105.063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45195,0.0610551,-8.77441e-05,7.07674e-08,-2.33835e-11,-1578.29,23.357], Tmin=(100,'K'), Tmax=(735.917,'K')), NASAPolynomial(coeffs=[8.25855,0.0240574,-1.23305e-05,2.44834e-09,-1.73958e-13,-2580.07,-7.39355], Tmin=(735.917,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-13.8466,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(CsCsFFH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(oxirene) + radical(Csj(Cs-CdHH)(F1s)(F1s))"""),
)

species(
    label = '[CH2][C](F)F(185)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 C u1 p0 c0 {4,S} {5,S} {6,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-109.048,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,190,488,555,1236,1407],'cm^-1')),
        HinderedRotor(inertia=(0.00258864,'amu*angstrom^2'), symmetry=1, barrier=(7.63529,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.25819,0.0186259,-1.56413e-05,8.223e-09,-2.04577e-12,-13090.7,13.5552], Tmin=(100,'K'), Tmax=(887.666,'K')), NASAPolynomial(coeffs=[4.4701,0.0131647,-6.4128e-06,1.29206e-09,-9.37476e-14,-13305.9,7.85287], Tmin=(887.666,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-109.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Cs_P) + radical(CsCsF1sF1s)"""),
)

species(
    label = 'F[C]1[CH]O1(2636)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {4,S}
3 C u1 p0 c0 {2,S} {4,S} {5,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
5 H u0 p0 c0 {3,S}
"""),
    E0 = (97.7552,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,180,778.222,778.25,778.304,778.371,778.444,3391.61],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.0271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.18921,0.013453,-5.29115e-07,-1.06379e-08,5.3496e-12,11790.3,13.5601], Tmin=(100,'K'), Tmax=(993.338,'K')), NASAPolynomial(coeffs=[8.19759,0.00359119,-1.19995e-06,2.57134e-10,-2.11238e-14,10286.9,-13.1283], Tmin=(993.338,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(97.7552,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCsJO) + radical(CsCsF1sO2s)"""),
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
    label = 'FC(F)=CC1(F)[CH]O1(2732)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {7,S}
5  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
6  C u0 p0 c0 {5,S} {8,D} {9,S}
7  C u1 p0 c0 {4,S} {5,S} {10,S}
8  C u0 p0 c0 {2,S} {3,S} {6,D}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-445.457,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,1000,182,240,577,636,1210,1413,496.943,496.943,496.948,496.948,496.948,496.948,496.949,496.949,496.968,2205.98],'cm^-1')),
        HinderedRotor(inertia=(0.000682619,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30906,0.05842,-6.19793e-05,3.25065e-08,-6.75353e-12,-53478.5,22.853], Tmin=(100,'K'), Tmax=(1164.42,'K')), NASAPolynomial(coeffs=[13.5152,0.0164893,-7.96392e-06,1.58079e-09,-1.13745e-13,-56321.1,-37.8925], Tmin=(1164.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-445.457,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCFO) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(CdCFF) + ring(Cs-Cs(F)-O2s) + radical(CCsJO)"""),
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
    label = 'F[C](F)C=C1[CH]O1(2733)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {4,S} {6,S}
4 C u0 p0 c0 {3,S} {5,D} {6,S}
5 C u0 p0 c0 {4,D} {7,S} {8,S}
6 C u1 p0 c0 {3,S} {4,S} {9,S}
7 C u1 p0 c0 {1,S} {2,S} {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-129.682,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,1000,161,297,490,584,780,1358,533.01,533.627,533.664,533.725,534.026,534.04,534.106],'cm^-1')),
        HinderedRotor(inertia=(0.000591904,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.055,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6809,0.0321344,2.58233e-05,-7.17406e-08,3.36552e-11,-15496.3,19.1094], Tmin=(100,'K'), Tmax=(936.705,'K')), NASAPolynomial(coeffs=[20.768,-0.00135318,2.55202e-06,-4.49503e-10,2.13687e-14,-21178.8,-82.972], Tmin=(936.705,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-129.682,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(CsCFFH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(methyleneoxirane) + radical(C=CCJO) + radical(Csj(Cd-CdH)(F1s)(F1s))"""),
)

species(
    label = 'CF2(87)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 F u0 p3 c0 {3,S}
3 C u0 p1 c0 {1,S} {2,S}
"""),
    E0 = (-203.712,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([192,594,627],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (50.0074,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(897.963,'J/mol'), sigma=(3.977,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28591,0.0107608,-1.05382e-05,4.89881e-09,-8.86384e-13,-24340.7,13.1348], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[5.33121,0.00197748,-9.60248e-07,2.10704e-10,-1.5954e-14,-25190.9,-2.56367], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-203.712,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CF2""", comment="""Thermo library: halogens"""),
)

species(
    label = 'FC(F)[CH]C1(F)[CH]O1(2734)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  O u0 p2 c0 {5,S} {8,S}
5  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
6  C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
7  C u1 p0 c0 {5,S} {6,S} {10,S}
8  C u1 p0 c0 {4,S} {5,S} {11,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-341.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44907,0.0552202,-5.08538e-05,2.38228e-08,-4.51938e-12,-40991.3,25.8248], Tmin=(100,'K'), Tmax=(1250.85,'K')), NASAPolynomial(coeffs=[12.1881,0.0208789,-9.67229e-06,1.87429e-09,-1.3266e-13,-43677.9,-28.3882], Tmin=(1250.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-341.59,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCFO) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(CsCsFFH) + ring(Cs(F)(C)-Cs-O2s) + radical(CCJCO) + radical(Csj(Cs-F1sO2sCs)(O2s-Cs)(H)_ring)"""),
)

species(
    label = 'F[C](F)[CH]C1(F)CO1(2735)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {5,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {6,S}
5  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
6  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
7  C u1 p0 c0 {5,S} {8,S} {11,S}
8  C u1 p0 c0 {2,S} {3,S} {7,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-355.68,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76294,0.0524582,-4.92517e-05,2.51908e-08,-5.45294e-12,-42700.6,25.7883], Tmin=(100,'K'), Tmax=(1076.46,'K')), NASAPolynomial(coeffs=[8.84871,0.0261286,-1.25631e-05,2.46926e-09,-1.76095e-13,-44226.1,-8.91863], Tmin=(1076.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-355.68,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCFO) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(CsCsFFH) + ring(Cs(F)(C)-Cs-O2s) + radical(CCJCO) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = 'F[C](F)C[C]1OC1F(2723)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {6,S} {7,S}
5  C u0 p0 c0 {7,S} {8,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {4,S} {7,S} {11,S}
7  C u1 p0 c0 {4,S} {5,S} {6,S}
8  C u1 p0 c0 {2,S} {3,S} {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-368.974,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41381,0.0633602,-9.41538e-05,8.55235e-08,-3.16798e-11,-44290.5,25.9582], Tmin=(100,'K'), Tmax=(782.893,'K')), NASAPolynomial(coeffs=[5.44493,0.0328213,-1.6592e-05,3.25415e-09,-2.28702e-13,-44617,9.44305], Tmin=(782.893,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-368.974,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCFHO) + group(CsCsFFH) + ring(O2s-Cs(F)-Cs(C)) + radical(C2CsJO) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = 'FC(F)(F)C[C]1[CH]O1(2736)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {6,S}
4  O u0 p2 c0 {7,S} {8,S}
5  C u0 p0 c0 {6,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
7  C u1 p0 c0 {4,S} {5,S} {8,S}
8  C u1 p0 c0 {4,S} {7,S} {11,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-420.681,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.061,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.615155,0.0723317,-9.47893e-05,6.34468e-08,-1.61273e-11,-50472.4,24.7987], Tmin=(100,'K'), Tmax=(1092.56,'K')), NASAPolynomial(coeffs=[14.0354,0.0144636,-3.34867e-06,3.33335e-10,-1.12292e-14,-52883.6,-38.7486], Tmin=(1092.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-420.681,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(CsCsFFF) + ring(Ethylene_oxide) + radical(C2CsJO) + radical(CCsJO)"""),
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
    E0 = (-104.165,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (16.9482,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-3.12354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (297.918,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-149.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-94.554,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-109.837,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (259.426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (148.403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-39.3866,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-50.6556,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (171.704,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (71.7679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (63.4084,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (22.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-15.2802,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (10.4511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (32.265,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['F[C](F)CC1(F)[CH]O1(2714)'],
    products = ['FC1=CO1(390)', 'CH2CF2(56)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(53.7893,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission
Ea raised from 0.0 to 53.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction2',
    reactants = ['F[C](F)CC1O[C]1F(2715)'],
    products = ['F[C](F)CC1(F)[CH]O1(2714)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(F)(F)C1(F)[CH]O1(2712)'],
    products = ['F[C](F)CC1(F)[CH]O1(2714)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CF2(42)', '[CH2]C1(F)[CH]O1(2727)'],
    products = ['F[C](F)CC1(F)[CH]O1(2714)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['F[C](F)CC1(F)[CH]O1(2714)'],
    products = ['FC1(F)CC2(F)OC12(2728)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_H/NonDeO;Cpri_rad_out_noH]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['F[C](F)CC1(F)[CH]O1(2714)'],
    products = ['FC(F)=CC1(F)CO1(2729)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C=C(F)C[C](F)F(2730)'],
    products = ['F[C](F)CC1(F)[CH]O1(2714)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(181.793,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F(37)', 'F[C](F)CC1=CO1(2731)'],
    products = ['F[C](F)CC1(F)[CH]O1(2714)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(17.3848,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction9',
    reactants = ['FC1=CO1(390)', '[CH2][C](F)F(185)'],
    products = ['F[C](F)CC1(F)[CH]O1(2714)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction10',
    reactants = ['F[C]1[CH]O1(2636)', 'CH2CF2(56)'],
    products = ['F[C](F)CC1(F)[CH]O1(2714)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(0.001023,'m^3/(mol*s)'), n=2.607, Ea=(41.4777,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_4R!H->O',), comment="""Estimated from node Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_4R!H->O"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(5)', 'FC(F)=CC1(F)[CH]O1(2732)'],
    products = ['F[C](F)CC1(F)[CH]O1(2714)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(8.42983,'m^3/(mol*s)'), n=1.92039, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.021512250841106063, var=1.4639131859069563, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_Ext-4R!H-R_Sp-5R!H-4R!H_Sp-2CS=1CCOSS',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_Ext-4R!H-R_Sp-5R!H-4R!H_Sp-2CS=1CCOSS"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F[C]1[CH]O1(2636)', '[CH2][C](F)F(185)'],
    products = ['F[C](F)CC1(F)[CH]O1(2714)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(5e+07,'m^3/(mol*s)'), n=2.59562e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing"""),
)

reaction(
    label = 'reaction13',
    reactants = ['HF(38)', 'F[C](F)C=C1[CH]O1(2733)'],
    products = ['F[C](F)CC1(F)[CH]O1(2714)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(299.567,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['CF2(87)', '[CH2]C1(F)[CH]O1(2727)'],
    products = ['F[C](F)CC1(F)[CH]O1(2714)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(2.93008,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F[C](F)CC1(F)[CH]O1(2714)'],
    products = ['FC(F)[CH]C1(F)[CH]O1(2734)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.76836e+09,'s^-1'), n=1.1815, Ea=(180.499,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_noH;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F[C](F)CC1(F)[CH]O1(2714)'],
    products = ['F[C](F)[CH]C1(F)CO1(2735)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4e-15,'s^-1'), n=8.23, Ea=(142.674,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_H/NonDeO;XH_out] for rate rule [R3H_SS_12cy3;C_rad_out_H/NonDeO;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F[C](F)CC1(F)[CH]O1(2714)'],
    products = ['F[C](F)C[C]1OC1F(2723)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.00763e+12,'s^-1'), n=0.18834, Ea=(168.405,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C"""),
)

reaction(
    label = 'reaction18',
    reactants = ['F[C](F)CC1(F)[CH]O1(2714)'],
    products = ['FC(F)(F)C[C]1[CH]O1(2736)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(190.219,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

network(
    label = 'PDepNetwork #1052',
    isomers = [
        'F[C](F)CC1(F)[CH]O1(2714)',
    ],
    reactants = [
        ('FC1=CO1(390)', 'CH2CF2(56)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1052',
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

