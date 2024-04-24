species(
    label = '[CH2]C(F)(F)O[O](317)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {4,S} {5,S}
4 O u1 p2 c0 {3,S}
5 C u0 p0 c0 {1,S} {2,S} {3,S} {6,S}
6 C u1 p0 c0 {5,S} {7,S} {8,S}
7 H u0 p0 c0 {6,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-283.94,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,253,525,597,667,842,1178,1324,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.377869,'amu*angstrom^2'), symmetry=1, barrier=(8.68795,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.378453,'amu*angstrom^2'), symmetry=1, barrier=(8.70138,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0328,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3366.69,'J/mol'), sigma=(5.8238,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=525.87 K, Pc=38.67 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62595,0.0580252,-0.000101027,8.97264e-08,-3.07049e-11,-34070.2,20.3244], Tmin=(100,'K'), Tmax=(835.988,'K')), NASAPolynomial(coeffs=[8.26983,0.0163097,-8.36704e-06,1.63108e-09,-1.12917e-13,-34834.2,-8.46377], Tmin=(835.988,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-283.94,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(ROOJ) + radical(CJCOOH)"""),
)

species(
    label = 'O2(2)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {2,S}
2 O u1 p2 c0 {1,S}
"""),
    E0 = (-8.62683,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1487.4],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.6,'angstroms^3'), rotrelaxcollnum=3.8, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53732,-0.00121572,5.3162e-06,-4.89446e-09,1.45846e-12,-1038.59,4.68368], Tmin=(100,'K'), Tmax=(1074.55,'K')), NASAPolynomial(coeffs=[3.15382,0.00167804,-7.69974e-07,1.51275e-10,-1.08782e-14,-1040.82,6.16756], Tmin=(1074.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.62683,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2""", comment="""Thermo library: primaryThermoLibrary"""),
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
    label = 'O=O(177)',
    structure = adjacencyList("""1 O u0 p2 c0 {2,D}
2 O u0 p2 c0 {1,D}
"""),
    E0 = (85.6848,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1487.4],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1857.18,'J/mol'), sigma=(4.34667,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=290.09 K, Pc=51.31 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53732,-0.00121572,5.3162e-06,-4.89446e-09,1.45846e-12,10304.5,4.68368], Tmin=(100,'K'), Tmax=(1074.55,'K')), NASAPolynomial(coeffs=[3.15382,0.00167804,-7.69974e-07,1.51275e-10,-1.08782e-14,10302.3,6.16756], Tmin=(1074.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(85.6848,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2(S)""", comment="""Thermo library: primaryThermoLibrary"""),
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
    label = 'FC1(F)CO1(1188)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {4,S} {5,S}
4 C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
5 C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
6 H u0 p0 c0 {5,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-503.243,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,588.2,588.296,588.378,588.4,588.421,589.068,800.385,1507.87,1507.92,1507.95,1508.05],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2932.76,'J/mol'), sigma=(5.28172,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=458.09 K, Pc=45.16 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.09124,-0.00918197,0.00013802,-2.55242e-07,1.50557e-10,-60523.8,9.31183], Tmin=(10,'K'), Tmax=(551.22,'K')), NASAPolynomial(coeffs=[2.29346,0.0259294,-1.7572e-05,5.55773e-09,-6.63159e-13,-60660.8,13.8735], Tmin=(551.22,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-503.243,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""FC1(F)CO1""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[CH2]C([O])(F)F(732)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u1 p2 c0 {4,S}
4 C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
5 C u1 p0 c0 {4,S} {6,S} {7,S}
6 H u0 p0 c0 {5,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-268.405,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([253,525,597,667,842,1178,1324,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.181121,'amu*angstrom^2'), symmetry=1, barrier=(4.16432,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.24484,0.0356831,-3.90797e-05,2.09126e-08,-4.3132e-12,-32216,16.3892], Tmin=(100,'K'), Tmax=(1195.24,'K')), NASAPolynomial(coeffs=[10.9836,0.00643772,-2.37702e-06,4.40781e-10,-3.12138e-14,-34304.9,-27.3284], Tmin=(1195.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-268.405,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sF1s)) + radical(Csj(Cs-F1sF1sO2s)(H)(H))"""),
)

species(
    label = '[O]O[C](F)F(171)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {4,S} {5,S}
4 O u1 p2 c0 {3,S}
5 C u1 p0 c0 {1,S} {2,S} {3,S}
"""),
    E0 = (-195.642,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,493,600,700,1144,1293],'cm^-1')),
        HinderedRotor(inertia=(0.711545,'amu*angstrom^2'), symmetry=1, barrier=(16.3598,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.0062,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3217.18,'J/mol'), sigma=(5.19785,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=502.52 K, Pc=51.98 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.8245,0.013883,6.3822e-05,-2.51617e-07,2.43597e-10,-23530.5,10.5618], Tmin=(10,'K'), Tmax=(394.515,'K')), NASAPolynomial(coeffs=[5.8912,0.0131528,-1.02966e-05,3.57092e-09,-4.54892e-13,-23851,0.518718], Tmin=(394.515,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-195.642,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""[O]O[C](F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'CH2(T)(66)',
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
    label = 'FC1(F)COO1(895)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {4,S} {5,S}
4 O u0 p2 c0 {3,S} {6,S}
5 C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
6 C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-471.012,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0328,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.92793,0.00417742,9.37145e-05,-1.78011e-07,1.01769e-10,-56646.3,10.4485], Tmin=(10,'K'), Tmax=(573.559,'K')), NASAPolynomial(coeffs=[2.28708,0.0321309,-2.2569e-05,7.33691e-09,-8.95255e-13,-56729.7,15.0851], Tmin=(573.559,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-471.012,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), label="""FC1(F)COO1""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'C=C(F)O[O](829)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 O u0 p2 c0 {3,S} {4,S}
3 O u1 p2 c0 {2,S}
4 C u0 p0 c0 {1,S} {2,S} {5,D}
5 C u0 p0 c0 {4,D} {6,S} {7,S}
6 H u0 p0 c0 {5,S}
7 H u0 p0 c0 {5,S}
"""),
    E0 = (-96.7758,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,326,540,652,719,1357,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.725193,'amu*angstrom^2'), symmetry=1, barrier=(16.6736,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (77.0344,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.91162,0.00529573,8.44313e-05,-1.79248e-07,1.11122e-10,-11633.8,10.2654], Tmin=(10,'K'), Tmax=(555.797,'K')), NASAPolynomial(coeffs=[4.21889,0.0228301,-1.61816e-05,5.35569e-09,-6.6581e-13,-11972.9,6.21972], Tmin=(555.797,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-96.7758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CDC(F)O[O]""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[CH2][C](F)OOF(898)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {4,S} {5,S}
4 O u0 p2 c0 {2,S} {3,S}
5 C u1 p0 c0 {1,S} {3,S} {6,S}
6 C u1 p0 c0 {5,S} {7,S} {8,S}
7 H u0 p0 c0 {6,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (90.8727,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([277,555,632,395,473,707,1436,3000,3100,440,815,1455,1000,180,1986.48],'cm^-1')),
        HinderedRotor(inertia=(0.352669,'amu*angstrom^2'), symmetry=1, barrier=(8.10856,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.353261,'amu*angstrom^2'), symmetry=1, barrier=(8.12216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.63175,'amu*angstrom^2'), symmetry=1, barrier=(37.5171,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0328,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28882,0.0697418,-0.000136157,1.29417e-07,-4.60463e-11,11017.4,22.1512], Tmin=(100,'K'), Tmax=(855.36,'K')), NASAPolynomial(coeffs=[7.20094,0.0202847,-1.118e-05,2.20098e-09,-1.51778e-13,10803.8,-0.78386], Tmin=(855.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(90.8727,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCsF1sO2s) + radical(CJCOOH)"""),
)

species(
    label = '[O]O[C](F)CF(899)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {4,S} {6,S}
4 O u1 p2 c0 {3,S}
5 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
6 C u1 p0 c0 {2,S} {3,S} {5,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-233.159,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,551,1088,1226,1380,1420,1481,3057,3119,395,473,707,1436,180],'cm^-1')),
        HinderedRotor(inertia=(0.39753,'amu*angstrom^2'), symmetry=1, barrier=(9.14,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.396475,'amu*angstrom^2'), symmetry=1, barrier=(9.11574,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0328,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77128,0.0580337,-0.000110663,1.05771e-07,-3.75987e-11,-27970.8,20.3106], Tmin=(100,'K'), Tmax=(874.09,'K')), NASAPolynomial(coeffs=[5.5126,0.0201695,-1.00883e-05,1.91315e-09,-1.29079e-13,-27832.5,7.29722], Tmin=(874.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-233.159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(ROOJ) + radical(CsCsF1sO2s)"""),
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
    E0 = (-189.042,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-181.072,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (69.4887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (280.588,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-180.796,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-55.6104,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (114.516,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-9.76211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (268.755,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-8.73557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction7',
    reactants = ['O2(2)', 'CH2CF2(56)'],
    products = ['[CH2]C(F)(F)O[O](317)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(10834.1,'m^3/(mol*s)'), n=0.762578, Ea=(86.3416,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_N-3R->C_N-2R!H->O_N-4R!H-u0_Ext-1R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_N-3R->C_N-2R!H->O_N-4R!H-u0_Ext-1R!H-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(F)(F)O[O](317)'],
    products = ['O=O(177)', 'CH2CF2(56)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(8.00887,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission
Ea raised from 0.0 to 8.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction2',
    reactants = ['O(6)', '[CH2]C([O])(F)F(732)'],
    products = ['[CH2]C(F)(F)O[O](317)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(43772.1,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]O[C](F)F(171)', 'CH2(T)(66)'],
    products = ['[CH2]C(F)(F)O[O](317)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(F)(F)O[O](317)'],
    products = ['FC1(F)COO1(895)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R4_SSS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]C(F)(F)O[O](317)'],
    products = ['O(6)', 'FC1(F)CO1(1188)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+09,'s^-1'), n=1.1, Ea=(133.47,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 35 used for R2OO_S;C_pri_rad_intra;OOJ
Exact match found for rate rule [R2OO_S;C_pri_rad_intra;OOJ]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['F(37)', 'C=C(F)O[O](829)'],
    products = ['[CH2]C(F)(F)O[O](317)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(43.54,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction8',
    reactants = ['O2(2)', '[CH2][C](F)F(185)'],
    products = ['[CH2]C(F)(F)O[O](317)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.62995e+14,'m^3/(mol*s)'), n=-2.27223, Ea=(13.053,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.03642541003782715, var=3.51219038874387, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_3R!H->O_2R->C_N-2C-inRing_Ext-2C-R_Ext-4R!H-R_Sp-5R!H-4R!H
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2][C](F)OOF(898)'],
    products = ['[CH2]C(F)(F)O[O](317)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(83.0222,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]O[C](F)CF(899)'],
    products = ['[CH2]C(F)(F)O[O](317)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.36622e+11,'s^-1'), n=0.48217, Ea=(129.563,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-2R!H-R',), comment="""Estimated from node R2F_Ext-2R!H-R"""),
)

network(
    label = 'PDepNetwork #108',
    isomers = [
        '[CH2]C(F)(F)O[O](317)',
    ],
    reactants = [
        ('O2(2)', 'CH2CF2(56)'),
        ('O=O(177)', 'CH2CF2(56)'),
        ('O(6)', 'FC1(F)CO1(1188)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #108',
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

