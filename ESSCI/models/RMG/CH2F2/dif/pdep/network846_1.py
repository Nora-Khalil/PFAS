species(
    label = 'O=C(OF)C(O)[C](F)F(1919)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {5,S}
4  O u0 p2 c0 {7,S} {11,S}
5  O u0 p2 c0 {3,S} {8,S}
6  O u0 p2 c0 {8,D}
7  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {5,S} {6,D} {7,S}
9  C u1 p0 c0 {1,S} {2,S} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-712.375,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,990,1113,1380,1390,370,380,2900,435,190,488,555,1236,1407,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (143.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.862644,0.0729983,-8.66629e-05,5.16643e-08,-1.23537e-11,-85569.2,29.619], Tmin=(100,'K'), Tmax=(1009.94,'K')), NASAPolynomial(coeffs=[13.5269,0.0228407,-1.21685e-05,2.49121e-09,-1.81637e-13,-88127.3,-31.6044], Tmin=(1009.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-712.375,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(Cs-(Cds-O2d)CsOsH) + group(CsCsFFH) + group(Cds-OdCsOs) + radical(CsCsF1sF1s)"""),
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
    label = 'O=C[C](F)F(1089)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {5,D}
4 C u1 p0 c0 {1,S} {2,S} {5,S}
5 C u0 p0 c0 {3,D} {4,S} {6,S}
6 H u0 p0 c0 {5,S}
"""),
    E0 = (-391.269,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([179,346,818,1406,1524,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.55681,'amu*angstrom^2'), symmetry=1, barrier=(31.8674,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.0255,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2992.84,'J/mol'), sigma=(4.8347,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=467.48 K, Pc=60.09 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.92143,0.00544306,7.55479e-05,-1.88818e-07,1.40422e-10,-47057,10.1435], Tmin=(10,'K'), Tmax=(452.346,'K')), NASAPolynomial(coeffs=[3.74847,0.0196353,-1.35046e-05,4.31321e-09,-5.19454e-13,-47170.9,9.4087], Tmin=(452.346,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-391.269,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""ODC[C](F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'OC(OF)[C](F)F(2000)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {5,S}
4 O u0 p2 c0 {6,S} {9,S}
5 O u0 p2 c0 {3,S} {6,S}
6 C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
7 C u1 p0 c0 {1,S} {2,S} {6,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {4,S}
"""),
    E0 = (-511.149,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,557,1111,1380,1390,370,380,2900,435,190,488,555,1236,1407,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.399974,'amu*angstrom^2'), symmetry=1, barrier=(9.19619,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.396184,'amu*angstrom^2'), symmetry=1, barrier=(9.10906,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.28944,'amu*angstrom^2'), symmetry=1, barrier=(29.6468,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (115.031,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.48571,0.0449182,3.35714e-05,-2.91541e-07,3.19242e-10,-61478.9,12.3546], Tmin=(10,'K'), Tmax=(395.558,'K')), NASAPolynomial(coeffs=[8.52708,0.0264957,-2.00293e-05,6.87674e-09,-8.73401e-13,-62132.5,-10.511], Tmin=(395.558,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-511.149,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), label="""OC(OF)[C](F)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C(OF)C(F)(F)[CH]O(2032)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {9,S}
5  O u0 p2 c0 {8,S} {11,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
8  C u1 p0 c0 {5,S} {7,S} {10,S}
9  C u0 p0 c0 {4,S} {6,D} {7,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-715.325,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,3615,1277.5,1000,249,341,381,651,928,1217,1251,3025,407.5,1350,352.5,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (143.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0710342,0.0948883,-0.000160983,1.40339e-07,-4.7868e-11,-85900.4,28.1944], Tmin=(100,'K'), Tmax=(808.18,'K')), NASAPolynomial(coeffs=[11.6341,0.0268034,-1.44692e-05,2.86107e-09,-2.00154e-13,-87414.9,-22.935], Tmin=(808.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-715.325,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + radical(CCsJOH)"""),
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
    label = 'O=[C]C(O)=C(F)F(2096)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {5,S} {8,S}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {3,S} {6,D} {7,S}
6 C u0 p0 c0 {1,S} {2,S} {5,D}
7 C u1 p0 c0 {4,D} {5,S}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (-484.671,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,182,240,577,636,1210,1413,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(0.994546,'amu*angstrom^2'), symmetry=1, barrier=(22.8666,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.98686,'amu*angstrom^2'), symmetry=1, barrier=(22.6898,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.931299,0.0731613,-0.000124841,1.03489e-07,-3.34526e-11,-58187.3,20.1853], Tmin=(100,'K'), Tmax=(763.583,'K')), NASAPolynomial(coeffs=[12.2562,0.013828,-8.26943e-06,1.69922e-09,-1.2153e-13,-59916.6,-31.3946], Tmin=(763.583,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-484.671,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-O2d)O2s) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJ=O)"""),
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
    label = 'O=C(C=C(F)F)OF(848)',
    structure = adjacencyList("""1 F u0 p3 c0 {8,S}
2 F u0 p3 c0 {8,S}
3 F u0 p3 c0 {4,S}
4 O u0 p2 c0 {3,S} {7,S}
5 O u0 p2 c0 {7,D}
6 C u0 p0 c0 {7,S} {8,D} {9,S}
7 C u0 p0 c0 {4,S} {5,D} {6,S}
8 C u0 p0 c0 {1,S} {2,S} {6,D}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-545.904,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,180,4000,4000,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.252137,'amu*angstrom^2'), symmetry=1, barrier=(5.79713,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.253356,'amu*angstrom^2'), symmetry=1, barrier=(5.82515,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (126.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85389,0.0536526,-9.7139e-05,9.19911e-08,-3.39189e-11,-65586.1,19.9008], Tmin=(100,'K'), Tmax=(790.008,'K')), NASAPolynomial(coeffs=[6.67042,0.0187161,-1.07748e-05,2.20796e-09,-1.57861e-13,-66017.9,-0.117383], Tmin=(790.008,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-545.904,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + missing(Cd-COCdH) + group(Cds-O2d(Cds-Cds)O2s) + group(CdCFF)"""),
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
    label = 'O=C(OF)C(O)[C]F(2097)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {6,S} {10,S}
4  O u0 p2 c0 {2,S} {7,S}
5  O u0 p2 c0 {7,D}
6  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {4,S} {5,D} {6,S}
8  C u2 p0 c0 {1,S} {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-322.645,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,990,1113,1380,1390,370,380,2900,435,90,150,250,247,323,377,431,433,1065,1274],'cm^-1')),
        HinderedRotor(inertia=(0.000700407,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000700668,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.253572,'amu*angstrom^2'), symmetry=1, barrier=(43.3045,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30099,0.060382,-6.41203e-05,3.33318e-08,-6.9142e-12,-38708.9,25.9313], Tmin=(100,'K'), Tmax=(1157.93,'K')), NASAPolynomial(coeffs=[13.4723,0.0183369,-9.6544e-06,1.97362e-09,-1.43882e-13,-41527.6,-34.573], Tmin=(1157.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-322.645,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(Cs-(Cds-O2d)CsOsH) + group(CsCsFHH) + group(Cds-OdCsOs) + radical(CsCCl_triplet)"""),
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
    label = 'OC1[C](OF)OC1(F)F(2098)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {6,S}
4  O u0 p2 c0 {8,S} {9,S}
5  O u0 p2 c0 {7,S} {11,S}
6  O u0 p2 c0 {3,S} {9,S}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
9  C u1 p0 c0 {4,S} {6,S} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {5,S}
"""),
    E0 = (-611.685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (143.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00754,0.0658245,-6.99619e-05,3.65935e-08,-7.60229e-12,-73460.9,27.2995], Tmin=(100,'K'), Tmax=(1161.23,'K')), NASAPolynomial(coeffs=[14.5766,0.0190844,-9.58616e-06,1.93161e-09,-1.39975e-13,-76612.2,-40.1915], Tmin=(1161.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-611.685,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2sCF) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(CsCFFO) + ring(Cs-Cs(F)(F)-O2s-Cs) + radical(Cs_P)"""),
)

species(
    label = '[O]C1(OF)C(O)C1(F)F(2099)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {5,S}
4  O u0 p2 c0 {7,S} {11,S}
5  O u0 p2 c0 {3,S} {8,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {5,S} {6,S} {7,S} {9,S}
9  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-513.098,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (143.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.419117,0.0763808,-8.83109e-05,4.91347e-08,-1.0628e-11,-61580.4,26.3743], Tmin=(100,'K'), Tmax=(1133.02,'K')), NASAPolynomial(coeffs=[17.7617,0.0151542,-7.25236e-06,1.43955e-09,-1.0397e-13,-65510.2,-59.4591], Tmin=(1133.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-513.098,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsOsH) + group(CsCsCsFF) + ring(Cs-Cs(F)(F)-Cs(O2)) + radical(CC(C)(O)OJ)"""),
)

species(
    label = 'OC(OF)=C(O)[C](F)F(2100)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {5,S}
4  O u0 p2 c0 {7,S} {10,S}
5  O u0 p2 c0 {3,S} {8,S}
6  O u0 p2 c0 {8,S} {11,S}
7  C u0 p0 c0 {4,S} {8,D} {9,S}
8  C u0 p0 c0 {5,S} {6,S} {7,D}
9  C u1 p0 c0 {1,S} {2,S} {7,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {6,S}
"""),
    E0 = (-608.917,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3580,3650,1210,1345,900,1100,231,791,325,375,415,465,420,450,1700,1750,161,297,490,584,780,1358,180],'cm^-1')),
        HinderedRotor(inertia=(0.633148,'amu*angstrom^2'), symmetry=1, barrier=(14.5573,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.633294,'amu*angstrom^2'), symmetry=1, barrier=(14.5607,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.63442,'amu*angstrom^2'), symmetry=1, barrier=(14.5866,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.44921,'amu*angstrom^2'), symmetry=1, barrier=(33.3202,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (143.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.605089,0.107602,-0.000182628,1.47692e-07,-4.56379e-11,-73076,28.765], Tmin=(100,'K'), Tmax=(867.372,'K')), NASAPolynomial(coeffs=[17.5179,0.0151999,-7.56702e-06,1.40758e-09,-9.35697e-14,-75887.8,-54.1753], Tmin=(867.372,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-608.917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2sCF) + group(O2s-(Cds-Cd)H) + group(CsCFFH) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsCs) + radical(Csj(Cd-O2sCd)(F1s)(F1s))"""),
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
    label = 'OC=C(F)F(535)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {4,S} {7,S}
4 C u0 p0 c0 {3,S} {5,D} {6,S}
5 C u0 p0 c0 {1,S} {2,S} {4,D}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (-502.96,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413],'cm^-1')),
        HinderedRotor(inertia=(0.552091,'amu*angstrom^2'), symmetry=1, barrier=(12.6937,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3292.1,'J/mol'), sigma=(5.16209,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=514.22 K, Pc=54.31 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.85967,0.00956102,9.95844e-05,-2.79091e-07,2.18422e-10,-60486.4,9.26794], Tmin=(10,'K'), Tmax=(463.499,'K')), NASAPolynomial(coeffs=[6.00669,0.0174593,-1.1501e-05,3.70008e-09,-4.58213e-13,-60969.3,-2.50148], Tmin=(463.499,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-502.96,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""OCDC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'O=C(OF)C(O)=C(F)F(2101)',
    structure = adjacencyList("""1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {5,S}
4  O u0 p2 c0 {7,S} {10,S}
5  O u0 p2 c0 {3,S} {8,S}
6  O u0 p2 c0 {8,D}
7  C u0 p0 c0 {4,S} {8,S} {9,D}
8  C u0 p0 c0 {5,S} {6,D} {7,S}
9  C u0 p0 c0 {1,S} {2,S} {7,D}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (-683.291,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,990,1113,350,440,435,1725,182,240,577,636,1210,1413,348.04,348.04,348.04,348.04,348.041,2204.22],'cm^-1')),
        HinderedRotor(inertia=(0.091116,'amu*angstrom^2'), symmetry=1, barrier=(7.83219,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0911169,'amu*angstrom^2'), symmetry=1, barrier=(7.83229,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.271628,'amu*angstrom^2'), symmetry=1, barrier=(23.3485,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.27855,0.0904548,-0.000159796,1.41042e-07,-4.81887e-11,-82055.1,26.7048], Tmin=(100,'K'), Tmax=(819.027,'K')), NASAPolynomial(coeffs=[11.5664,0.0227305,-1.26937e-05,2.52735e-09,-1.76805e-13,-83481.7,-22.9201], Tmin=(819.027,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-683.291,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2sCF) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-Cds)O2s) + group(CdCFF) + longDistanceInteraction_noncyclic(Cd(F)2=CdOs)"""),
)

species(
    label = '[O]C(=O)C(O)[C](F)F(2102)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  O u0 p2 c0 {6,S} {10,S}
4  O u1 p2 c0 {8,S}
5  O u0 p2 c0 {8,D}
6  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
7  C u1 p0 c0 {1,S} {2,S} {6,S}
8  C u0 p0 c0 {4,S} {5,D} {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-621.178,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,190,488,555,1236,1407,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56576,0.059443,-7.75638e-05,5.88488e-08,-1.8935e-11,-74628.2,27.0946], Tmin=(100,'K'), Tmax=(744.741,'K')), NASAPolynomial(coeffs=[7.32042,0.0285336,-1.53061e-05,3.1158e-09,-2.25421e-13,-75485.3,1.02792], Tmin=(744.741,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-621.178,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsOsH) + group(CsCsFFH) + group(Cds-OdCsOs) + radical(CCOJ) + radical(CsCsF1sF1s)"""),
)

species(
    label = 'O=C([CH][C](F)F)OF(746)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {8,S}
2 F u0 p3 c0 {8,S}
3 F u0 p3 c0 {4,S}
4 O u0 p2 c0 {3,S} {7,S}
5 O u0 p2 c0 {7,D}
6 C u1 p0 c0 {7,S} {8,S} {9,S}
7 C u0 p0 c0 {4,S} {5,D} {6,S}
8 C u1 p0 c0 {1,S} {2,S} {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-343.622,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,3025,407.5,1350,352.5,190,488,555,1236,1407,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06783,0.0716563,-0.000121557,1.08902e-07,-3.82985e-11,-41229.5,26.2213], Tmin=(100,'K'), Tmax=(805.217,'K')), NASAPolynomial(coeffs=[8.50455,0.0240488,-1.30041e-05,2.57931e-09,-1.80997e-13,-42081.4,-5.89856], Tmin=(805.217,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-343.622,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFFH) + group(Cds-OdCsOs) + radical(CCJCO) + radical(Csj(Cs-COHH)(F1s)(F1s))"""),
)

species(
    label = '[O]C([C](F)F)C(=O)OF(2080)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {8,D}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {4,S} {6,D} {7,S}
9  C u1 p0 c0 {1,S} {2,S} {7,S}
10 H u0 p0 c0 {7,S}
"""),
    E0 = (-468.642,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,1380,1390,370,380,2900,435,190,488,555,1236,1407,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02905,0.0708816,-9.29905e-05,6.27456e-08,-1.7069e-11,-56262.4,29.807], Tmin=(100,'K'), Tmax=(891.104,'K')), NASAPolynomial(coeffs=[11.6716,0.0231067,-1.25669e-05,2.57486e-09,-1.87199e-13,-58159,-20.3098], Tmin=(891.104,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-468.642,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(Cs-(Cds-O2d)CsOsH) + group(CsCsFFH) + group(Cds-OdCsOs) + radical(C=OCOJ) + radical(CsCsF1sF1s)"""),
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
    label = 'O=[C]C(O)[C](F)F(2103)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {5,S} {9,S}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
6 C u1 p0 c0 {1,S} {2,S} {5,S}
7 C u1 p0 c0 {4,D} {5,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {3,S}
"""),
    E0 = (-421.236,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,190,488,555,1236,1407,1855,455,950,280.216],'cm^-1')),
        HinderedRotor(inertia=(0.00216244,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.205073,'amu*angstrom^2'), symmetry=1, barrier=(11.3331,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204255,'amu*angstrom^2'), symmetry=1, barrier=(11.3532,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5017,0.0565062,-7.33032e-05,4.71115e-08,-1.18634e-11,-50574.3,26.6575], Tmin=(100,'K'), Tmax=(973.965,'K')), NASAPolynomial(coeffs=[12.0973,0.0129907,-6.28493e-06,1.23815e-09,-8.8477e-14,-52638.2,-24.1805], Tmin=(973.965,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-421.236,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(CsCsFFH) + group(Cds-OdCsH) + radical(CsCsF1sF1s) + radical(CCCJ=O)"""),
)

species(
    label = 'O[CH][C](F)F(1121)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {4,S} {7,S}
4 C u1 p0 c0 {3,S} {5,S} {6,S}
5 C u1 p0 c0 {1,S} {2,S} {4,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {3,S}
"""),
    E0 = (-278.999,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,190,488,555,1236,1407,180],'cm^-1')),
        HinderedRotor(inertia=(0.636341,'amu*angstrom^2'), symmetry=1, barrier=(14.6307,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.636327,'amu*angstrom^2'), symmetry=1, barrier=(14.6304,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0334,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72308,0.0525948,-8.51295e-05,6.61097e-08,-1.96969e-11,-33476.2,16.8524], Tmin=(100,'K'), Tmax=(832.734,'K')), NASAPolynomial(coeffs=[11.1819,0.00716936,-3.32206e-06,6.30396e-10,-4.31411e-14,-35051.8,-27.0516], Tmin=(832.734,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-278.999,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCsJOH) + radical(Csj(Cs-O2sHH)(F1s)(F1s))"""),
)

species(
    label = 'O=C(OF)[C](O)[C](F)F(2104)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {9,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {5,S}
4  O u0 p2 c0 {7,S} {10,S}
5  O u0 p2 c0 {3,S} {8,S}
6  O u0 p2 c0 {8,D}
7  C u1 p0 c0 {4,S} {8,S} {9,S}
8  C u0 p0 c0 {5,S} {6,D} {7,S}
9  C u1 p0 c0 {1,S} {2,S} {7,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (-535.748,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,990,1113,360,370,350,190,488,555,1236,1407,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.457257,0.0847864,-0.000134936,1.09535e-07,-3.52851e-11,-64314.4,31.1858], Tmin=(100,'K'), Tmax=(761.777,'K')), NASAPolynomial(coeffs=[12.2756,0.0227275,-1.27323e-05,2.58537e-09,-1.85001e-13,-66114.9,-22.6146], Tmin=(761.777,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-535.748,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(Cs-(Cds-O2d)CsOsH) + group(CsCsFFH) + group(Cds-OdCsOs) + radical(C2CsJOH) + radical(CsCsF1sF1s)"""),
)

species(
    label = 'O=C(OF)C(O)[C]F-2(2105)',
    structure = adjacencyList("""1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {4,S}
3  O u0 p2 c0 {6,S} {10,S}
4  O u0 p2 c0 {2,S} {7,S}
5  O u0 p2 c0 {7,D}
6  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
7  C u0 p0 c0 {4,S} {5,D} {6,S}
8  C u0 p1 c0 {1,S} {6,S}
9  H u0 p0 c0 {6,S}
10 H u0 p0 c0 {3,S}
"""),
    E0 = (-327.529,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,990,1113,1380,1390,370,380,2900,435,617,898,1187,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (124.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36755,0.0620295,-7.50884e-05,4.77069e-08,-1.23201e-11,-39301.3,25.8189], Tmin=(100,'K'), Tmax=(933.879,'K')), NASAPolynomial(coeffs=[10.6115,0.0224352,-1.14909e-05,2.30599e-09,-1.66072e-13,-41027.8,-18.1452], Tmin=(933.879,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-327.529,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(Cs-CsCsOsH) + group(Cds-OdCsOs) + group(CJ2_singlet-FCs)"""),
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
    label = 'O=C(OF)[C](O)C(F)F(2082)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {5,S}
4  O u0 p2 c0 {8,S} {11,S}
5  O u0 p2 c0 {3,S} {9,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {10,S}
8  C u1 p0 c0 {4,S} {7,S} {9,S}
9  C u0 p0 c0 {5,S} {6,D} {8,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-738.097,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (143.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.721208,0.0776365,-0.000101105,6.72696e-08,-1.79664e-11,-88659.1,28.6216], Tmin=(100,'K'), Tmax=(908.859,'K')), NASAPolynomial(coeffs=[12.8467,0.0242689,-1.30232e-05,2.65742e-09,-1.92881e-13,-90863.1,-28.7177], Tmin=(908.859,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-738.097,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(Cs-(Cds-O2d)CsOsH) + group(CsCsFFH) + group(Cds-OdCsOs) + radical(C2CsJOH)"""),
)

species(
    label = '[O]C(C(=O)OF)C(F)F(1923)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {9,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {1,S} {2,S} {7,S} {11,S}
9  C u0 p0 c0 {4,S} {6,D} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-670.991,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,1380,1390,370,380,2900,435,235,523,627,1123,1142,1372,1406,3097,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (143.041,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3809.89,'J/mol'), sigma=(5.79401,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=595.10 K, Pc=44.44 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07864,0.066586,-7.06698e-05,3.7614e-08,-8.07527e-12,-80598.3,28.6829], Tmin=(100,'K'), Tmax=(1114.47,'K')), NASAPolynomial(coeffs=[13.27,0.0228304,-1.17789e-05,2.38665e-09,-1.73148e-13,-83315.8,-31.4547], Tmin=(1114.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-670.991,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(Cs-(Cds-O2d)CsOsH) + group(CsCsFFH) + group(Cds-OdCsOs) + radical(C=OCOJ)"""),
)

species(
    label = '[O]C(=O)C(O)C(F)(F)F(2106)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {8,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {7,S} {11,S}
5  O u1 p2 c0 {9,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
8  C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
9  C u0 p0 c0 {5,S} {6,D} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {4,S}
"""),
    E0 = (-1051.95,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (143.041,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36573,0.0631797,-6.64736e-05,3.65469e-08,-8.35697e-12,-126430,26.1242], Tmin=(100,'K'), Tmax=(1031.21,'K')), NASAPolynomial(coeffs=[10.5706,0.0274734,-1.45337e-05,2.96731e-09,-2.15906e-13,-128328,-18.5668], Tmin=(1031.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-1051.95,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsOsH) + group(CsCsFFF) + longDistanceInteraction_noncyclic(Cs(F)3-CsOs) + longDistanceInteraction_noncyclic(Cs(F)3-R-CO) + group(Cds-OdCsOs) + radical(CCOJ)"""),
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
    E0 = (-375.598,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-113.387,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-275.479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-93.807,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-53.2289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (67.3423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (34.9691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-153.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-195.831,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-120.308,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-208.992,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-154.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-215.211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (1.86825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (60.2582,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-1.45443,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (2.41325,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-5.14797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (62.4581,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-197.434,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-263.902,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (-353.896,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-325.609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=C(OF)C(O)[C](F)F(1919)'],
    products = ['HF(38)', 'CO2(13)', 'O=C[C](F)F(1089)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(0.00285451,'s^-1'), n=4.62568, Ea=(19.6811,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3917696014098424, var=52.52930589142705, Tref=1000.0, N=14, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(12)', 'OC(OF)[C](F)F(2000)'],
    products = ['O=C(OF)C(O)[C](F)F(1919)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(0.000742267,'m^3/(mol*s)'), n=2.83796, Ea=(199.407,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.2632766423807829, var=145.9379134136138, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-1COCbCdCsCtHNOSSidSis->Cs',), comment="""Estimated from node Root_N-1COCbCdCsCtHNOSSidSis->Cs"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O=C(OF)C(F)(F)[CH]O(2032)'],
    products = ['O=C(OF)C(O)[C](F)F(1919)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.889e+11,'s^-1'), n=0.232, Ea=(122.75,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCsCJ;CsJ;CO_O] for rate rule [cCs(-R!HR!H)CJ;CsJ;CO_O]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['OF(173)', 'O=[C]C(O)=C(F)F(2096)'],
    products = ['O=C(OF)C(O)[C](F)F(1919)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(0.000227174,'m^3/(mol*s)'), n=2.96333, Ea=(169.034,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_Cdd;H_OH]
Euclidian distance = 0
family: 1,3_Insertion_ROR"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O=C(OF)C(O)[C](F)F(1919)'],
    products = ['OH(4)', 'O=C(C=C(F)F)OF(848)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.62709e+20,'s^-1'), n=-1.9758, Ea=(342.051,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.1791595854394145, var=100.97114496904264, Tref=1000.0, N=30, data_mean=0.0, correlation='Root',), comment="""Estimated from node Root"""),
)

reaction(
    label = 'reaction6',
    reactants = ['F(37)', 'O=C(OF)C(O)[C]F(2097)'],
    products = ['O=C(OF)C(O)[C](F)F(1919)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [H/Val7_rad;Birad] for rate rule [Val7_rad;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F[C]F(166)', 'O=C([CH]O)OF(482)'],
    products = ['O=C(OF)C(O)[C](F)F(1919)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination
Ea raised from -1.7 to -1.7 kJ/mol.
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction8',
    reactants = ['O=C(OF)C(O)[C](F)F(1919)'],
    products = ['OC1[C](OF)OC1(F)F(2098)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(242.17,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra_Nd;radadd_intra_cs]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['O=C(OF)C(O)[C](F)F(1919)'],
    products = ['[O]C1(OF)C(O)C1(F)F(2099)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.63856e+09,'s^-1'), n=0.755479, Ea=(199.449,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs] for rate rule [R4_S_CO;carbonylbond_intra_Nd;radadd_intra_cs]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['OC(OF)=C(O)[C](F)F(2100)'],
    products = ['O=C(OF)C(O)[C](F)F(1919)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(205000,'s^-1'), n=2.37, Ea=(171.514,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R!H->C_Ext-1R!H-R',), comment="""Estimated from node Root_3R!H->C_Ext-1R!H-R"""),
)

reaction(
    label = 'reaction11',
    reactants = ['O=[C]OF(142)', 'OC=C(F)F(535)'],
    products = ['O=C(OF)C(O)[C](F)F(1919)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(153.031,'m^3/(mol*s)'), n=1.16366, Ea=(12.5559,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.022037706214473284, var=2.3701416358838845, Tref=1000.0, N=230, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(5)', 'O=C(OF)C(O)=C(F)F(2101)'],
    products = ['O=C(OF)C(O)[C](F)F(1919)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(56.8734,'m^3/(mol*s)'), n=1.75834, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.2648472359903826, var=0.02782886759889551, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_4R!H->C_Sp-4C-1COS_Ext-1COS-R_N-6R!H-inRing_N-7R!H-inRing_Ext-4C-R_N-Sp-8R!H#4C',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R_4R!H->C_Sp-4C-1COS_Ext-1COS-R_N-6R!H-inRing_N-7R!H-inRing_Ext-4C-R_N-Sp-8R!H#4C"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F(37)', '[O]C(=O)C(O)[C](F)F(2102)'],
    products = ['O=C(OF)C(O)[C](F)F(1919)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.808e+07,'m^3/(mol*s)'), n=2.17087e-08, Ea=(15.98,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction14',
    reactants = ['OH(4)', 'O=C([CH][C](F)F)OF(746)'],
    products = ['O=C(OF)C(O)[C](F)F(1919)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.7e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_2R->C_Ext-2C-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(5)', '[O]C([C](F)F)C(=O)OF(2080)'],
    products = ['O=C(OF)C(O)[C](F)F(1919)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.21037e+06,'m^3/(mol*s)'), n=0.349925, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_2BrCClFHNO->O_Ext-2O-R_N-3R!H->Cl_3BrCFINOPSSi->C"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]F(228)', 'O=[C]C(O)[C](F)F(2103)'],
    products = ['O=C(OF)C(O)[C](F)F(1919)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O=[C]OF(142)', 'O[CH][C](F)F(1121)'],
    products = ['O=C(OF)C(O)[C](F)F(1919)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.63131e-11,'m^3/(mol*s)'), n=4.71246, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_N-1BrCFS-inRing_1BrCFS->C_N-2R->S_N-2BrCF->Br_Ext-1C-R_3R!H->F_Ext-2CF-R_Ext-4R!H-R
Ea raised from -17.6 to 0.0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(5)', 'O=C(OF)[C](O)[C](F)F(2104)'],
    products = ['O=C(OF)C(O)[C](F)F(1919)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(1.69912,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2BrCClFHNO-inRing_N-2BrCClFHNO->O_N-2CHN->N_2CH->C_Ext-2C-R_3R!H->F_Ext-2C-R_4R!H->C_Ext-4C-R_N-5R!H->Cl_N-5BrCFINOPSSi->Br_5CF->C"""),
)

reaction(
    label = 'reaction19',
    reactants = ['F(37)', 'O=C(OF)C(O)[C]F-2(2105)'],
    products = ['O=C(OF)C(O)[C](F)F(1919)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(518506,'m^3/(mol*s)'), n=0.4717, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->H_N-2Br1sCl1sF1s->Cl1s_3BrCClFINOPSSi->F',), comment="""Estimated from node Root_N-3R->H_N-2Br1sCl1sF1s->Cl1s_3BrCClFINOPSSi->F"""),
)

reaction(
    label = 'reaction20',
    reactants = ['CF2(42)', 'O=C([CH]O)OF(482)'],
    products = ['O=C(OF)C(O)[C](F)F(1919)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.00976185,'m^3/(mol*s)'), n=2.64543, Ea=(5.03668,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24524072153041726, var=3.368891203868511, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s"""),
)

reaction(
    label = 'reaction21',
    reactants = ['O=C(OF)C(O)[C](F)F(1919)'],
    products = ['O=C(OF)[C](O)C(F)F(2082)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(8.2e+09,'s^-1'), n=0.65, Ea=(131.378,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_noH;Cs_H_out_OneDe] for rate rule [R2H_S;C_rad_out_noH;Cs_H_out_CO]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['O=C(OF)C(O)[C](F)F(1919)'],
    products = ['[O]C(C(=O)OF)C(F)F(1923)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(5.248e+06,'s^-1'), n=0.992, Ea=(41.3841,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 441 used for R3H_SS_Cs;C_rad_out_noH;O_H_out
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_noH;O_H_out]
Euclidian distance = 0
family: intra_H_migration
Ea raised from 40.5 to 41.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction23',
    reactants = ['O=C(OF)C(O)[C](F)F(1919)'],
    products = ['[O]C(=O)C(O)C(F)(F)F(2106)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(69.6704,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

network(
    label = 'PDepNetwork #846',
    isomers = [
        'O=C(OF)C(O)[C](F)F(1919)',
    ],
    reactants = [
        ('HF(38)', 'CO2(13)', 'O=C[C](F)F(1089)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #846',
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

