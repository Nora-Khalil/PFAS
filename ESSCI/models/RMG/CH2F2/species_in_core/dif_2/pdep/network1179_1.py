species(
    label = 'O=C=C(F)[O+]=[C-]O(3004)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 O u0 p1 c+1 {5,S} {6,D}
3 O u0 p2 c0 {6,S} {8,S}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {1,S} {2,S} {7,D}
6 C u0 p1 c-1 {2,D} {3,S}
7 C u0 p0 c0 {4,D} {5,D}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (-55.2119,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,2120,512.5,787.5,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.3178,'amu*angstrom^2'), symmetry=1, barrier=(53.2909,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.35964,'amu*angstrom^2'), symmetry=1, barrier=(54.2528,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45096,0.0728374,-0.000157629,1.63844e-07,-6.13749e-11,-6564.67,33.2499], Tmin=(100,'K'), Tmax=(875.735,'K')), NASAPolynomial(coeffs=[1.36019,0.0320795,-1.72948e-05,3.32685e-09,-2.25226e-13,-4969.99,42.6898], Tmin=(875.735,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-55.2119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + missing(O2d-Cdd) + group(Cd(Cdd-Od)FO) + group(CsJ2_singlet-CsH) + missing(Cdd-CdO2d)"""),
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
    label = 'O=CC(=O)F(557)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 O u0 p2 c0 {4,D}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {2,D} {5,S} {6,S}
5 C u0 p0 c0 {1,S} {3,D} {4,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-483.116,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,286,619,818,1246,1924],'cm^-1')),
        HinderedRotor(inertia=(0.753127,'amu*angstrom^2'), symmetry=1, barrier=(17.3159,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (76.0265,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3426.83,'J/mol'), sigma=(5.15159,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=535.26 K, Pc=56.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88796,0.00794544,8.17071e-05,-2.34537e-07,1.90378e-10,-58102.5,9.57665], Tmin=(10,'K'), Tmax=(438.71,'K')), NASAPolynomial(coeffs=[4.97623,0.0166174,-1.15194e-05,3.74172e-09,-4.59031e-13,-58377,3.18363], Tmin=(438.71,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-483.116,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""ODCC(DO)F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = 'O=C=C=[O+][C-](O)F(3120)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 O u0 p1 c+1 {5,S} {6,D}
3 O u0 p2 c0 {5,S} {8,S}
4 O u0 p2 c0 {7,D}
5 C u0 p1 c-1 {1,S} {2,S} {3,S}
6 C u0 p0 c0 {2,D} {7,D}
7 C u0 p0 c0 {4,D} {6,D}
8 H u0 p0 c0 {3,S}
"""),
    E0 = (-182.277,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.63078,0.0350496,-6.19829e-05,6.01705e-08,-2.24408e-11,-21878.3,13.9843], Tmin=(100,'K'), Tmax=(823.69,'K')), NASAPolynomial(coeffs=[4.51974,0.0156921,-8.18497e-06,1.61743e-09,-1.13283e-13,-21844,7.33475], Tmin=(823.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-182.277,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + missing(O2d-Cdd) + group(CJ2_singlet-FO) + group(Cdd-CdsCds) + missing(Cdd-CddO2d)"""),
)

species(
    label = 'O=C=C(F)[OH+][C-]=O(3121)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 O u0 p1 c+1 {5,S} {6,S} {8,S}
3 O u0 p2 c0 {6,D}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {1,S} {2,S} {7,D}
6 C u0 p1 c-1 {2,S} {3,D}
7 C u0 p0 c0 {4,D} {5,D}
8 H u0 p0 c0 {2,S}
"""),
    E0 = (113.127,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.10794,0.053146,-0.000111821,1.17107e-07,-4.47443e-11,13663.1,18.6481], Tmin=(100,'K'), Tmax=(854.722,'K')), NASAPolynomial(coeffs=[2.18972,0.0250377,-1.3835e-05,2.72749e-09,-1.88615e-13,14661.9,24.1908], Tmin=(854.722,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(113.127,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + missing(O2d-C2dc) + missing(O2d-Cdd) + group(Cd(Cdd-Od)FO) + group(CsJ2_singlet-CsH) + missing(Cdd-CdO2d)"""),
)

species(
    label = 'O=C(F)C(=O)[C-]=[OH+](3122)',
    structure = adjacencyList("""1 F u0 p3 c0 {6,S}
2 O u0 p1 c+1 {7,D} {8,S}
3 O u0 p2 c0 {5,D}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {3,D} {6,S} {7,S}
6 C u0 p0 c0 {1,S} {4,D} {5,S}
7 C u0 p1 c-1 {2,D} {5,S}
8 H u0 p0 c0 {2,S}
"""),
    E0 = (-209.092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (104.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78429,0.0590807,-0.000115448,1.16595e-07,-4.41114e-11,-25078.2,22.5923], Tmin=(100,'K'), Tmax=(837.732,'K')), NASAPolynomial(coeffs=[3.74754,0.0262355,-1.46108e-05,2.90483e-09,-2.02754e-13,-24583.5,18.384], Tmin=(837.732,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-209.092,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cds-OdCsCs) + group(COCFO) + group(CsJ2_singlet-CsH)"""),
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
    E0 = (66.3969,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (306.609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (258.349,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (154.319,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=C=C(F)[O+]=[C-]O(3004)'],
    products = ['CO(12)', 'O=CC(=O)F(557)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(4.02942e+10,'s^-1'), n=0.375329, Ea=(5.99296,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.06684469846581474, var=24.965000021544366, Tref=1000.0, N=36, data_mean=0.0, correlation='Root_N-1R!H->C_N-5R!H->N',), comment="""Estimated from node Root_N-1R!H->C_N-5R!H->N"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O=C=C(F)[O+]=[C-]O(3004)'],
    products = ['O=C=C=[O+][C-](O)F(3120)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(246.205,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O=C=C(F)[O+]=[C-]O(3004)'],
    products = ['O=C=C(F)[OH+][C-]=O(3121)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.87617e-12,'s^-1'), n=7.25159, Ea=(197.945,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.05095623138621764, var=6.9856901089733165, Tref=1000.0, N=14, data_mean=0.0, correlation='Root_N-3R!H->C',), comment="""Estimated from node Root_N-3R!H->C"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O=C=C(F)[O+]=[C-]O(3004)'],
    products = ['O=C(F)C(=O)[C-]=[OH+](3122)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(93.9152,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

network(
    label = 'PDepNetwork #1179',
    isomers = [
        'O=C=C(F)[O+]=[C-]O(3004)',
    ],
    reactants = [
        ('CO(12)', 'O=CC(=O)F(557)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1179',
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

