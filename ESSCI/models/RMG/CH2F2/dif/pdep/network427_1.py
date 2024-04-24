species(
    label = '[O]C(=O)C[C](F)F(1165)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 O u1 p2 c0 {7,S}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {6,S} {7,S} {8,S} {9,S}
6 C u1 p0 c0 {1,S} {2,S} {5,S}
7 C u0 p0 c0 {3,S} {4,D} {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-452.327,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,190,488,555,1236,1407,180,180,180,1377.35,1378.04,1378.16,1378.89],'cm^-1')),
        HinderedRotor(inertia=(0.221778,'amu*angstrom^2'), symmetry=1, barrier=(5.0991,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.42182,'amu*angstrom^2'), symmetry=1, barrier=(55.6825,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7044,0.061873,-0.00011628,1.1942e-07,-4.59804e-11,-54330.7,20.7415], Tmin=(100,'K'), Tmax=(840.042,'K')), NASAPolynomial(coeffs=[2.1008,0.033344,-1.77664e-05,3.48496e-09,-2.41975e-13,-53457.3,24.4933], Tmin=(840.042,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-452.327,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFFH) + group(Cds-OdCsOs) + radical(CCOJ) + radical(Csj(Cs-COHH)(F1s)(F1s))"""),
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
    label = '[CH2]C(F)(F)C([O])=O(1164)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u1 p2 c0 {7,S}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6 C u1 p0 c0 {5,S} {8,S} {9,S}
7 C u0 p0 c0 {3,S} {4,D} {5,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-428.505,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([249,341,381,651,928,1217,1251,3000,3100,440,815,1455,1000,180,679.632,679.713,679.852,679.938,680.189],'cm^-1')),
        HinderedRotor(inertia=(0.00932817,'amu*angstrom^2'), symmetry=1, barrier=(3.06296,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133474,'amu*angstrom^2'), symmetry=1, barrier=(3.06882,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67946,0.0592279,-0.00010078,9.9451e-08,-3.85058e-11,-51461.7,21.3936], Tmin=(100,'K'), Tmax=(796.923,'K')), NASAPolynomial(coeffs=[3.99557,0.0311419,-1.69327e-05,3.38955e-09,-2.398e-13,-51308.1,14.0249], Tmin=(796.923,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-428.505,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(CCOJ) + radical(Csj(Cs-F1sF1sCO)(H)(H))"""),
)

species(
    label = 'F[C]F(166)',
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28591,0.0107608,-1.05382e-05,4.89881e-09,-8.86384e-13,4216.69,13.1348], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[5.33121,0.00197748,-9.60248e-07,2.10704e-10,-1.5954e-14,3366.49,-2.56367], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(33.7272,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CF2(T)""", comment="""Thermo library: halogens"""),
)

species(
    label = '[CH2]C([O])=O(1203)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {4,S}
2 O u0 p2 c0 {4,D}
3 C u1 p0 c0 {4,S} {5,S} {6,S}
4 C u0 p0 c0 {1,S} {2,D} {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-8.20105,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,844.726,844.726,844.734,844.77,844.77],'cm^-1')),
        HinderedRotor(inertia=(0.00484689,'amu*angstrom^2'), symmetry=1, barrier=(2.45426,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.60873,0.0103835,4.1641e-06,-6.72692e-09,1.73624e-12,-973.928,12.8424], Tmin=(100,'K'), Tmax=(1564.28,'K')), NASAPolynomial(coeffs=[5.7075,0.0133436,-6.65905e-06,1.28862e-09,-8.86361e-14,-2649.32,-1.47836], Tmin=(1564.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.20105,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCOJ) + radical(CJCO)"""),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,1.49298e-14,-2.05327e-17,9.23927e-21,-1.28064e-24,29230.2,5.12616], Tmin=(100,'K'), Tmax=(3959.16,'K')), NASAPolynomial(coeffs=[2.5,1.19009e-10,-4.23987e-14,6.68965e-18,-3.94352e-22,29230.2,5.12617], Tmin=(3959.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.034,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'O=[C]C[C](F)F(1204)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
5 C u1 p0 c0 {1,S} {2,S} {4,S}
6 C u1 p0 c0 {3,D} {4,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-243.768,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,190,488,555,1236,1407,1855,455,950,1192.2],'cm^-1')),
        HinderedRotor(inertia=(0.776288,'amu*angstrom^2'), symmetry=1, barrier=(17.8484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.62152,'amu*angstrom^2'), symmetry=1, barrier=(60.274,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0441,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.83288,0.0293857,-2.14755e-05,7.65592e-09,-1.07841e-12,-29321.8,12.1615], Tmin=(10,'K'), Tmax=(1615.06,'K')), NASAPolynomial(coeffs=[9.89848,0.0143632,-7.52343e-06,1.89685e-09,-1.86956e-13,-31281.1,-20.0092], Tmin=(1615.06,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-243.768,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), label="""OD[C]C[C](F)F""", comment="""Thermo library: 2-BTP_G4"""),
)

species(
    label = 'O=C1CC(F)(F)O1(1168)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {6,S} {7,S}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {6,S} {7,S} {8,S} {9,S}
6 C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
7 C u0 p0 c0 {3,S} {4,D} {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-756.881,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.49755,0.0204008,3.52742e-05,-6.47275e-08,2.7821e-11,-90965.9,17.2757], Tmin=(100,'K'), Tmax=(924.907,'K')), NASAPolynomial(coeffs=[12.1459,0.0111385,-2.35506e-06,3.45719e-10,-2.6022e-14,-94139.3,-36.0257], Tmin=(924.907,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-756.881,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)CsHH) + group(CsCFFO) + group(Cds-OdCsOs) + ring(Beta-Propiolactone)"""),
)

species(
    label = 'O=C(O)C=C(F)F(1087)',
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
        HinderedRotor(inertia=(0.000791703,'amu*angstrom^2'), symmetry=1, barrier=(8.98905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.390386,'amu*angstrom^2'), symmetry=1, barrier=(8.97574,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.35327,0.0410023,-6.93378e-05,6.55202e-08,-2.44567e-11,-81780,17.9383], Tmin=(100,'K'), Tmax=(783.697,'K')), NASAPolynomial(coeffs=[5.45578,0.0174189,-9.36918e-06,1.89143e-09,-1.34735e-13,-82028.3,5.24474], Tmin=(783.697,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-680.412,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + missing(Cd-COCdH) + group(Cds-O2d(Cds-Cds)O2s) + group(CdCFF)"""),
)

species(
    label = 'F[C](F)C[C]1OO1(1205)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {4,S} {6,S}
4 O u0 p2 c0 {3,S} {6,S}
5 C u0 p0 c0 {6,S} {7,S} {8,S} {9,S}
6 C u1 p0 c0 {3,S} {4,S} {5,S}
7 C u1 p0 c0 {1,S} {2,S} {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-90.9247,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31899,0.051773,-5.545e-05,2.92029e-08,-5.89767e-12,-10833.1,24.1692], Tmin=(100,'K'), Tmax=(1278.04,'K')), NASAPolynomial(coeffs=[14.5042,0.00855655,-2.44032e-06,3.57949e-10,-2.17951e-14,-14044.1,-42.0539], Tmin=(1278.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-90.9247,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(CsCsFFH) + ring(dioxirane) + radical(Cs_P) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = '[O][C]1CC(F)(F)O1(1206)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {6,S} {7,S}
4 O u1 p2 c0 {7,S}
5 C u0 p0 c0 {6,S} {7,S} {8,S} {9,S}
6 C u0 p0 c0 {1,S} {2,S} {3,S} {5,S}
7 C u1 p0 c0 {3,S} {4,S} {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-348.017,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.60456,0.0327412,-1.92513e-05,4.80938e-09,-4.56854e-13,-41808.9,20.592], Tmin=(100,'K'), Tmax=(2437.59,'K')), NASAPolynomial(coeffs=[16.4707,0.00998701,-5.24901e-06,9.79783e-10,-6.4083e-14,-48568.7,-58.6586], Tmin=(2437.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-348.017,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(CsCFFO) + ring(Cs-Cs(F)(F)-O2s-Cs) + radical(CCOJ) + radical(Cs_P)"""),
)

species(
    label = '[O]C1([O])CC1(F)F(1207)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 O u1 p2 c0 {7,S}
4 O u1 p2 c0 {7,S}
5 C u0 p0 c0 {6,S} {7,S} {8,S} {9,S}
6 C u0 p0 c0 {1,S} {2,S} {5,S} {7,S}
7 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-243.145,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.38289,0.0316712,-1.70962e-05,3.63448e-09,-2.83076e-13,-29241.1,15.0543], Tmin=(100,'K'), Tmax=(3111.25,'K')), NASAPolynomial(coeffs=[26.5534,0.00188103,-2.7333e-06,5.56754e-10,-3.57623e-14,-43658.5,-123.028], Tmin=(3111.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-243.145,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(CsCsCsFF) + ring(Cs-Cs(F)(F)-Cs(O2)) + radical(CC(C)(O)OJ) + radical(CC(C)(O)OJ)"""),
)

species(
    label = 'O=C(O)[CH][C](F)F(1208)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 O u0 p2 c0 {6,S} {9,S}
4 O u0 p2 c0 {6,D}
5 C u1 p0 c0 {6,S} {7,S} {8,S}
6 C u0 p0 c0 {3,S} {4,D} {5,S}
7 C u1 p0 c0 {1,S} {2,S} {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {3,S}
"""),
    E0 = (-478.13,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,190,488,555,1236,1407,180,180,180,1148.8,1148.85,1149.19],'cm^-1')),
        HinderedRotor(inertia=(0.00380231,'amu*angstrom^2'), symmetry=1, barrier=(3.56396,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.43947,'amu*angstrom^2'), symmetry=1, barrier=(56.0881,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155009,'amu*angstrom^2'), symmetry=1, barrier=(3.56397,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5662,0.0590188,-9.38038e-05,8.24975e-08,-2.88648e-11,-57423.4,24.2624], Tmin=(100,'K'), Tmax=(803.427,'K')), NASAPolynomial(coeffs=[7.29364,0.0227448,-1.15944e-05,2.26178e-09,-1.57786e-13,-58093.3,-0.557207], Tmin=(803.427,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-478.13,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFFH) + group(Cds-OdCsOs) + radical(CCJCO) + radical(Csj(Cs-COHH)(F1s)(F1s))"""),
)

species(
    label = '[O][C]=O(141)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u1 p0 c0 {1,S} {2,D}
"""),
    E0 = (31.5354,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0094,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.75939,0.00186757,1.032e-05,-1.52371e-08,5.80528e-12,3804.5,8.40408], Tmin=(100,'K'), Tmax=(1021.27,'K')), NASAPolynomial(coeffs=[6.3618,0.000422708,-4.06585e-07,1.52496e-10,-1.51971e-14,2816.74,-6.43931], Tmin=(1021.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.5354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(OJC=O) + radical((O)CJOH)"""),
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
    label = '[O]C(=O)C=C(F)F(752)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 O u1 p2 c0 {7,S}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {6,D} {7,S} {8,S}
6 C u0 p0 c0 {1,S} {2,S} {5,D}
7 C u0 p0 c0 {3,S} {4,D} {5,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-454.707,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,4000,4000,4000,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00126454,'amu*angstrom^2'), symmetry=1, barrier=(14.3577,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.98213,0.014039,-5.67082e-06,8.30239e-11,1.27737e-13,-54761,8.46224], Tmin=(100,'K'), Tmax=(2781.16,'K')), NASAPolynomial(coeffs=[36.5915,-0.0182618,4.65194e-06,-6.89871e-10,4.426e-14,-77433.1,-185.518], Tmin=(2781.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-454.707,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + missing(Cd-COCdH) + group(Cds-O2d(Cds-Cds)O2s) + group(CdCFF) + radical(CCOJ)"""),
)

species(
    label = '[CH2][C](F)F(1124)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.25818,0.018626,-1.56415e-05,8.2233e-09,-2.04589e-12,-13090.7,13.5553], Tmin=(100,'K'), Tmax=(887.631,'K')), NASAPolynomial(coeffs=[4.47008,0.0131648,-6.41282e-06,1.29207e-09,-9.3748e-14,-13305.9,7.85297], Tmin=(887.631,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-109.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Cs_P) + radical(CsCsF1sF1s)"""),
)

species(
    label = 'CF2(42)',
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
    label = '[O]C([O])=CC(F)F(1209)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u1 p2 c0 {7,S}
4 O u1 p2 c0 {7,S}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {8,S}
6 C u0 p0 c0 {5,S} {7,D} {9,S}
7 C u0 p0 c0 {3,S} {4,S} {6,D}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-470.846,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77611,0.0521873,-6.08365e-05,3.68585e-08,-9.06917e-12,-56552.3,20.3893], Tmin=(100,'K'), Tmax=(977.613,'K')), NASAPolynomial(coeffs=[9.93123,0.0188195,-9.63807e-06,1.94429e-09,-1.40634e-13,-58146.8,-18.7697], Tmin=(977.613,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-470.846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(CsCFFH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + radical(C=COJ) + radical(C=COJ)"""),
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
    E0 = (-213.025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-31.8852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (264.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (238.567,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-204.741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-149.625,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (148.377,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (29.1446,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-3.84352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-4.21621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-66.3507,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (12.3165,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-181.467,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (161.789,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (27.3879,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-89.202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(=O)C[C](F)F(1165)'],
    products = ['CO2(13)', 'CH2CF2(56)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(F)(F)C([O])=O(1164)'],
    products = ['[O]C(=O)C[C](F)F(1165)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ-HH;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;CO]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['F[C]F(166)', '[CH2]C([O])=O(1203)'],
    products = ['[O]C(=O)C[C](F)F(1165)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/CO;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -1.7 to -1.7 kJ/mol.
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(6)', 'O=[C]C[C](F)F(1204)'],
    products = ['[O]C(=O)C[C](F)F(1165)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [CO_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]C(=O)C[C](F)F(1165)'],
    products = ['O=C1CC(F)(F)O1(1168)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_noH;Opri_rad]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]C(=O)C[C](F)F(1165)'],
    products = ['O=C(O)C=C(F)F(1087)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.9748e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C(=O)C[C](F)F(1165)'],
    products = ['F[C](F)C[C]1OO1(1205)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.55936e+11,'s^-1'), n=0.551275, Ea=(361.402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra] for rate rule [R3_CO;carbonyl_intra_Nd;radadd_intra_O]
Euclidian distance = 2.449489742783178
family: Intra_R_Add_Endocyclic
Ea raised from 360.1 to 361.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C(=O)C[C](F)F(1165)'],
    products = ['[O][C]1CC(F)(F)O1(1206)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.006e+11,'s^-1'), n=0.221, Ea=(242.17,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_cs]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C(=O)C[C](F)F(1165)'],
    products = ['[O]C1([O])CC1(F)F(1207)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.63856e+09,'s^-1'), n=0.755479, Ea=(209.182,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs] for rate rule [R4_S_CO;carbonylbond_intra;radadd_intra_cs]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic
Ea raised from 207.2 to 209.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction10',
    reactants = ['O=C(O)[CH][C](F)F(1208)'],
    products = ['[O]C(=O)C[C](F)F(1165)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(234.612,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_Ext-1R!H-R',), comment="""Estimated from node Root_3R!H->C_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O][C]=O(141)', 'CH2CF2(56)'],
    products = ['[O]C(=O)C[C](F)F(1165)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(153.031,'m^3/(mol*s)'), n=1.16366, Ea=(24.4287,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.022037706214473284, var=2.3701416358838845, Tref=1000.0, N=230, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(5)', '[O]C(=O)C=C(F)F(752)'],
    products = ['[O]C(=O)C[C](F)F(1165)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(306600,'m^3/(mol*s)'), n=0.481, Ea=(15.9169,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_4R!H->C_Sp-4C-1COS_N-6R!H-inRing_Ext-4C-R_Ext-4C-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_4R!H->C_Sp-4C-1COS_N-6R!H-inRing_Ext-4C-R_Ext-4C-R"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CO2(13)', '[CH2][C](F)F(1124)'],
    products = ['[O]C(=O)C[C](F)F(1165)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(91.4175,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O][C]=O(141)', '[CH2][C](F)F(1124)'],
    products = ['[O]C(=O)C[C](F)F(1165)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CF2(42)', '[CH2]C([O])=O(1203)'],
    products = ['[O]C(=O)C[C](F)F(1165)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]C(=O)C[C](F)F(1165)'],
    products = ['[O]C([O])=CC(F)F(1209)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(13570.7,'s^-1'), n=2.68291, Ea=(123.823,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_noH;Cs_H_out_H/OneDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #427',
    isomers = [
        '[O]C(=O)C[C](F)F(1165)',
    ],
    reactants = [
        ('CO2(13)', 'CH2CF2(56)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #427',
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

