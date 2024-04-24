species(
    label = '[O]C(=O)O[CH]C(=O)F(3133)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {6,S} {8,S}
3 O u0 p2 c0 {7,D}
4 O u1 p2 c0 {8,S}
5 O u0 p2 c0 {8,D}
6 C u1 p0 c0 {2,S} {7,S} {9,S}
7 C u0 p0 c0 {1,S} {3,D} {6,S}
8 C u0 p0 c0 {2,S} {4,S} {5,D}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-527.907,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,611,648,830,1210,1753,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2419,0.0591955,-6.82679e-05,3.7347e-08,-7.94174e-12,-63392,25.4372], Tmin=(100,'K'), Tmax=(1150.23,'K')), NASAPolynomial(coeffs=[14.987,0.0113965,-5.93461e-06,1.21937e-09,-8.95692e-14,-66554,-42.799], Tmin=(1150.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-527.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + group(Cds-OdOsOs) + radical(OC=OOJ) + radical(CCsJOC(O))"""),
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
    label = '[O]C([O])=O(236)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {4,S}
2 O u1 p2 c0 {4,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
"""),
    E0 = (-60.3782,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([180,620.874,889.154,889.303,890.758,1871.83],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.0088,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36989,0.0105617,-1.27493e-06,-7.37835e-09,3.82327e-12,-7236.21,10.2409], Tmin=(100,'K'), Tmax=(995.953,'K')), NASAPolynomial(coeffs=[7.17796,0.00288735,-1.19278e-06,2.48515e-10,-1.94623e-14,-8372.65,-10.0126], Tmin=(995.953,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-60.3782,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CO3t1""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH]C(=O)F-2(2781)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 C u2 p0 c0 {3,S} {5,S}
5 H u0 p0 c0 {4,S}
"""),
    E0 = (-5.0725,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([486,617,768,1157,1926,180,1655.21,1655.36],'cm^-1')),
        HinderedRotor(inertia=(0.0191737,'amu*angstrom^2'), symmetry=1, barrier=(5.30701,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.0271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.32766,0.0152939,-1.10759e-05,3.91583e-09,-5.61423e-13,-586.537,12.101], Tmin=(100,'K'), Tmax=(1580.4,'K')), NASAPolynomial(coeffs=[6.5355,0.00717482,-3.3698e-06,6.6515e-10,-4.72046e-14,-1600.47,-4.84311], Tmin=(1580.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-5.0725,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCJ2_triplet)"""),
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
    label = 'O=[C]O[CH]C(=O)F(2655)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {5,S} {7,S}
3 O u0 p2 c0 {6,D}
4 O u0 p2 c0 {7,D}
5 C u1 p0 c0 {2,S} {6,S} {8,S}
6 C u0 p0 c0 {1,S} {3,D} {5,S}
7 C u1 p0 c0 {2,S} {4,D}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-335.839,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,611,648,830,1210,1753,1855,455,950,413.205,413.681,413.963],'cm^-1')),
        HinderedRotor(inertia=(0.12461,'amu*angstrom^2'), symmetry=1, barrier=(15.0817,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124353,'amu*angstrom^2'), symmetry=1, barrier=(15.0821,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124586,'amu*angstrom^2'), symmetry=1, barrier=(15.0851,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28611,0.0559816,-7.0658e-05,4.12365e-08,-9.13802e-12,-40291.2,23.5862], Tmin=(100,'K'), Tmax=(1121.26,'K')), NASAPolynomial(coeffs=[15.686,0.00461086,-1.93477e-06,3.75418e-10,-2.74252e-14,-43520.4,-47.5328], Tmin=(1121.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-335.839,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + group(Cds-OdOsH) + radical(CCsJOC(O)H) + radical((O)CJOCC)"""),
)

species(
    label = 'O=C1OC(C(=O)F)O1(3178)',
    structure = adjacencyList("""1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {6,S} {8,S}
3 O u0 p2 c0 {6,S} {8,S}
4 O u0 p2 c0 {7,D}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
7 C u0 p0 c0 {1,S} {4,D} {6,S}
8 C u0 p0 c0 {2,S} {3,S} {5,D}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-834.564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.33837,0.0167539,5.95803e-05,-9.30068e-08,3.66573e-11,-100297,22.2555], Tmin=(100,'K'), Tmax=(985.9,'K')), NASAPolynomial(coeffs=[16.3872,0.0092166,-4.20578e-06,1.01204e-09,-8.72125e-14,-105471,-57.5136], Tmin=(985.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-834.564,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-Cs(Cds-O2d)) + group(Cs-CsOsOsH) + group(COCsFO) + group(Cds-OdOsOs) + ring(Cyclobutane)"""),
)

species(
    label = '[O]C(=O)OC1O[C]1F(3203)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {6,S} {7,S}
3 O u0 p2 c0 {6,S} {8,S}
4 O u1 p2 c0 {8,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
7 C u1 p0 c0 {1,S} {2,S} {6,S}
8 C u0 p0 c0 {3,S} {4,S} {5,D}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-391.292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62428,0.0438776,-2.86949e-05,3.03041e-09,2.16227e-12,-46969,23.5335], Tmin=(100,'K'), Tmax=(1108.28,'K')), NASAPolynomial(coeffs=[14.2386,0.0128026,-6.19708e-06,1.26351e-09,-9.3334e-14,-50652.6,-42.6242], Tmin=(1108.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-391.292,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-O2d)) + group(O2s-(Cds-O2d)H) + group(Cs-CsOsOsH) + group(CsCFHO) + group(Cds-OdOsOs) + ring(Cs(O2)-O2s-Cs(F)) + radical(OC=OOJ) + radical(Csj(Cs-O2sO2sH)(F1s)(O2s-Cs)_ring)"""),
)

species(
    label = 'O=C(F)[CH]O[C]1OO1(3204)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {8,S}
2 O u0 p2 c0 {6,S} {7,S}
3 O u0 p2 c0 {4,S} {6,S}
4 O u0 p2 c0 {3,S} {6,S}
5 O u0 p2 c0 {8,D}
6 C u1 p0 c0 {2,S} {3,S} {4,S}
7 C u1 p0 c0 {2,S} {8,S} {9,S}
8 C u0 p0 c0 {1,S} {5,D} {7,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-158.174,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19024,0.0619337,-8.01886e-05,5.09241e-08,-1.25571e-11,-18922.8,25.7662], Tmin=(100,'K'), Tmax=(999.258,'K')), NASAPolynomial(coeffs=[13.6384,0.0121033,-5.38667e-06,1.01842e-09,-7.12174e-14,-21410.6,-34.2796], Tmin=(999.258,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-158.174,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-OsOsOsH) + group(Cs-(Cds-O2d)OsHH) + group(COCsFO) + ring(dioxirane) + radical(Cs_P) + radical(CCsJOCs)"""),
)

species(
    label = '[O][C]1OC(C(=O)F)O1(3205)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {6,S} {8,S}
3 O u0 p2 c0 {6,S} {8,S}
4 O u0 p2 c0 {7,D}
5 O u1 p2 c0 {8,S}
6 C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
7 C u0 p0 c0 {1,S} {4,D} {6,S}
8 C u1 p0 c0 {2,S} {3,S} {5,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-371.331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.27885,0.0223105,-1.23842e-06,-5.91045e-09,1.44327e-12,-44696.8,17.5212], Tmin=(100,'K'), Tmax=(2066.52,'K')), NASAPolynomial(coeffs=[20.9579,0.0119817,-9.67783e-06,1.95337e-09,-1.30038e-13,-56278.3,-86.3954], Tmin=(2066.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-371.331,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-OsOsOsH) + group(COCsFO) + ring(Cyclobutane) + radical(OCOJ) + radical(Cs_P)"""),
)

species(
    label = 'O=C1O[CH][C](F)OO1(3143)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {6,S} {8,S}
3 O u0 p2 c0 {4,S} {7,S}
4 O u0 p2 c0 {3,S} {8,S}
5 O u0 p2 c0 {8,D}
6 C u1 p0 c0 {2,S} {7,S} {9,S}
7 C u1 p0 c0 {1,S} {3,S} {6,S}
8 C u0 p0 c0 {2,S} {4,S} {5,D}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-406.92,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25187,0.0338959,4.73121e-05,-9.90326e-08,4.20486e-11,-48818.5,20.6314], Tmin=(100,'K'), Tmax=(987.857,'K')), NASAPolynomial(coeffs=[25.4816,0.000917443,-1.51067e-06,6.58897e-10,-7.14131e-14,-56783.6,-112.052], Tmin=(987.857,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-406.92,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(Cs-CsOsHH) + group(CsCFHO) + group(Cds-OdOsOs) + ring(Cyclohexanone) + radical(CCsJOC(O)) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]C1([O])OC1C(=O)F(3206)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {8,S}
2 O u0 p2 c0 {6,S} {7,S}
3 O u1 p2 c0 {7,S}
4 O u1 p2 c0 {7,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {2,S} {7,S} {8,S} {9,S}
7 C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
8 C u0 p0 c0 {1,S} {5,D} {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-366.833,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90837,0.0492429,-6.67174e-05,5.3052e-08,-1.64446e-11,-44047.5,25.7227], Tmin=(100,'K'), Tmax=(954.061,'K')), NASAPolynomial(coeffs=[6.18639,0.0215786,-7.92798e-06,1.28419e-09,-7.89335e-14,-44421,7.60526], Tmin=(954.061,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-366.833,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsOsOs) + group(COCsFO) + ring(Ethylene_oxide) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[O]C1(F)[CH]OC(=O)O1(3169)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {6,S} {8,S}
3 O u0 p2 c0 {7,S} {8,S}
4 O u1 p2 c0 {6,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
7 C u1 p0 c0 {3,S} {6,S} {9,S}
8 C u0 p0 c0 {2,S} {3,S} {5,D}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-522.487,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73449,0.0327237,2.09362e-05,-5.97165e-08,2.69414e-11,-62743.7,21.5276], Tmin=(100,'K'), Tmax=(970.084,'K')), NASAPolynomial(coeffs=[18.4842,0.00556111,-1.85516e-06,4.72874e-10,-4.49473e-14,-67965.1,-68.9335], Tmin=(970.084,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-522.487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(CsCFOO) + group(Cs-CsOsHH) + group(Cds-OdOsOs) + ring(Cyclopentane) + radical(O2sj(Cs-F1sO2sCs)) + radical(CCsJOC(O))"""),
)

species(
    label = '[O][C]=O(188)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.75939,0.00186756,1.032e-05,-1.52371e-08,5.8053e-12,3804.5,8.40408], Tmin=(100,'K'), Tmax=(1021.27,'K')), NASAPolynomial(coeffs=[6.36179,0.000422714,-4.06588e-07,1.52497e-10,-1.51971e-14,2816.75,-6.43929], Tmin=(1021.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.5354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(OJC=O) + radical((O)CJOH)"""),
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
    label = '[O]C(=O)OC=C=O(3207)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {5,S} {6,S}
2 O u1 p2 c0 {6,S}
3 O u0 p2 c0 {6,D}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {1,S} {7,D} {8,S}
6 C u0 p0 c0 {1,S} {2,S} {3,D}
7 C u0 p0 c0 {4,D} {5,D}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-315.151,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2120,512.5,787.5,424.8,424.862,425.179,425.452,425.567,425.922,425.981,427.051],'cm^-1')),
        HinderedRotor(inertia=(0.000930309,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167636,'amu*angstrom^2'), symmetry=1, barrier=(21.506,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (101.037,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34553,0.0525329,-5.92928e-05,2.95345e-08,-5.1467e-12,-37803.2,20.7182], Tmin=(100,'K'), Tmax=(1029.79,'K')), NASAPolynomial(coeffs=[15.9676,0.00434197,-1.63198e-06,3.20485e-10,-2.44253e-14,-41271,-52.4693], Tmin=(1029.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-315.151,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-O2d)H) + missing(O2d-Cdd) + group(Cds-(Cdd-O2d)OsH) + group(Cds-OdOsOs) + missing(Cdd-CdO2d) + radical(OC=OOJ)"""),
)

species(
    label = '[O][CH]C(=O)F(438)',
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
    label = 'O=[C][CH]OC(=O)OF(3208)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {6,S} {7,S}
3 O u0 p2 c0 {1,S} {7,S}
4 O u0 p2 c0 {7,D}
5 O u0 p2 c0 {8,D}
6 C u1 p0 c0 {2,S} {8,S} {9,S}
7 C u0 p0 c0 {2,S} {3,S} {4,D}
8 C u1 p0 c0 {5,D} {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-228.484,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,3025,407.5,1350,352.5,1855,455,950,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.741187,0.0656786,-7.41497e-05,3.80549e-08,-7.46357e-12,-27357.9,26.6847], Tmin=(100,'K'), Tmax=(1255.99,'K')), NASAPolynomial(coeffs=[19.1445,0.00706879,-4.15339e-06,9.0149e-10,-6.83282e-14,-31980.7,-66.2951], Tmin=(1255.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-228.484,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(191.233,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2sCF) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-OdOsOs) + radical(CCsJOC(O)) + radical(CsCJ=O)"""),
)

species(
    label = '[O]C(=O)OC(F)[C]=O(3195)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {6,S} {7,S}
3 O u1 p2 c0 {7,S}
4 O u0 p2 c0 {7,D}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
7 C u0 p0 c0 {2,S} {3,S} {4,D}
8 C u1 p0 c0 {5,D} {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-537.95,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([355,410,600,1181,1341,1420,3056,1855,455,950,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41654,0.0605324,-8.06068e-05,5.47807e-08,-1.48391e-11,-64610.6,24.221], Tmin=(100,'K'), Tmax=(900.074,'K')), NASAPolynomial(coeffs=[11.0688,0.0176362,-9.11801e-06,1.82965e-09,-1.31496e-13,-66348.1,-21.3295], Tmin=(900.074,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-537.95,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-(Cds-O2d)H) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(Cds-OdOsOs) + radical(OC=OOJ) + radical(COj(Cs-F1sO2sH)(O2d))"""),
)

species(
    label = '[O]C(=O)C([O])(F)C=O(3132)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u1 p2 c0 {6,S}
3 O u1 p2 c0 {8,S}
4 O u0 p2 c0 {7,D}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7 C u0 p0 c0 {4,D} {6,S} {9,S}
8 C u0 p0 c0 {3,S} {5,D} {6,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-472.149,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,180,180,180,180,934.393,934.491,934.523],'cm^-1')),
        HinderedRotor(inertia=(2.28346,'amu*angstrom^2'), symmetry=1, barrier=(52.5012,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0847007,'amu*angstrom^2'), symmetry=1, barrier=(52.502,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41538,0.0657411,-0.000114533,1.109e-07,-4.17586e-11,-56701.9,23.9666], Tmin=(100,'K'), Tmax=(810.739,'K')), NASAPolynomial(coeffs=[4.92978,0.0306075,-1.66081e-05,3.30474e-09,-2.32417e-13,-56687,11.3557], Tmin=(810.739,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-472.149,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-O2d)H) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsOs) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCOJ)"""),
)

species(
    label = '[O]C(=O)OC=[C]F(3209)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {5,S} {6,S}
3 O u1 p2 c0 {6,S}
4 O u0 p2 c0 {6,D}
5 C u0 p0 c0 {2,S} {7,D} {8,S}
6 C u0 p0 c0 {2,S} {3,S} {4,D}
7 C u1 p0 c0 {1,S} {5,D}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-193.869,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,167,640,1190,180,180,180,465.849,735.549,735.558,735.57,735.612],'cm^-1')),
        HinderedRotor(inertia=(0.0101776,'amu*angstrom^2'), symmetry=1, barrier=(3.90885,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.727819,'amu*angstrom^2'), symmetry=1, barrier=(16.734,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61749,0.0554552,-7.43608e-05,4.92111e-08,-1.28457e-11,-23233.8,20.9139], Tmin=(100,'K'), Tmax=(936.412,'K')), NASAPolynomial(coeffs=[11.4104,0.0136236,-7.35267e-06,1.50571e-09,-1.09563e-13,-25067.9,-25.6879], Tmin=(936.412,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-193.869,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-O2d)H) + group(Cds-CdsOsH) + group(Cds-OdOsOs) + group(CdCFH) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + radical(OC=OOJ) + radical(Cdj(Cd-O2sH)(F1s))"""),
)

species(
    label = 'O=C1OC=C(F)OO1(3149)',
    structure = adjacencyList("""1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {6,S} {8,S}
3 O u0 p2 c0 {4,S} {7,S}
4 O u0 p2 c0 {3,S} {8,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {2,S} {7,D} {9,S}
7 C u0 p0 c0 {1,S} {3,S} {6,D}
8 C u0 p0 c0 {2,S} {4,S} {5,D}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-604.491,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9693,0.0202651,7.09834e-05,-1.10532e-07,4.29957e-11,-72608.8,16.7804], Tmin=(100,'K'), Tmax=(1001.58,'K')), NASAPolynomial(coeffs=[19.5798,0.0108466,-6.13595e-06,1.52014e-09,-1.29605e-13,-79191.7,-83.4602], Tmin=(1001.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-604.491,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-O2d)) + group(Cds-CdsOsH) + group(CdCFO) + group(Cds-OdOsOs) + longDistanceInteraction_cyclic(Cd(F)=CdOs) + ring(Cyclohexane)"""),
)

species(
    label = '[O]C1([O])OC=C(F)O1(3198)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {6,S} {7,S}
3 O u0 p2 c0 {6,S} {8,S}
4 O u1 p2 c0 {6,S}
5 O u1 p2 c0 {6,S}
6 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
7 C u0 p0 c0 {1,S} {2,S} {8,D}
8 C u0 p0 c0 {3,S} {7,D} {9,S}
9 H u0 p0 c0 {8,S}
"""),
    E0 = (-341.284,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.93694,0.0292058,-9.14873e-06,-1.7991e-09,8.36809e-13,-41013.9,15.8684], Tmin=(100,'K'), Tmax=(1843.44,'K')), NASAPolynomial(coeffs=[13.2563,0.018439,-9.84674e-06,1.87409e-09,-1.25243e-13,-46793.7,-45.5856], Tmin=(1843.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-341.284,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-OsOsOsOs) + group(CdCFO) + group(Cds-CdsOsH) + ring(Cyclopentane) + radical(OCOJ) + radical(OCOJ) + longDistanceInteraction_cyclic(Cd(F)=CdOs)"""),
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
    label = '[O]C(=O)O[C]=C=O(3199)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {5,S} {6,S}
2 O u1 p2 c0 {5,S}
3 O u0 p2 c0 {5,D}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {1,S} {2,S} {3,D}
6 C u1 p0 c0 {1,S} {7,D}
7 C u0 p0 c0 {4,D} {6,D}
"""),
    E0 = (-75.407,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2120,512.5,787.5,433.208,433.574,433.925,434.699,434.738,435.176,435.337,1940.57],'cm^-1')),
        HinderedRotor(inertia=(0.125375,'amu*angstrom^2'), symmetry=1, barrier=(16.7956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000890528,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.03,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78539,0.0514693,-7.583e-05,5.42942e-08,-1.51472e-11,-8992.02,22.2761], Tmin=(100,'K'), Tmax=(882.756,'K')), NASAPolynomial(coeffs=[10.9347,0.0100097,-5.37836e-06,1.08649e-09,-7.80376e-14,-10607.3,-20.7228], Tmin=(882.756,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-75.407,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-O2d)H) + missing(O2d-Cdd) + group(Cds-OdOsOs) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(OC=OOJ) + radical(C=CJO)"""),
)

species(
    label = '[O]C(=O)O[C]=C(O)F(3210)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 O u0 p2 c0 {7,S} {8,S}
3 O u0 p2 c0 {6,S} {9,S}
4 O u1 p2 c0 {7,S}
5 O u0 p2 c0 {7,D}
6 C u0 p0 c0 {1,S} {3,S} {8,D}
7 C u0 p0 c0 {2,S} {4,S} {5,D}
8 C u1 p0 c0 {2,S} {6,D}
9 H u0 p0 c0 {3,S}
"""),
    E0 = (-415.295,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,293,496,537,1218,1685,370,180,180,180,314.288,827.071,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.149755,'amu*angstrom^2'), symmetry=1, barrier=(3.44316,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149755,'amu*angstrom^2'), symmetry=1, barrier=(3.44316,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149755,'amu*angstrom^2'), symmetry=1, barrier=(3.44316,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11221,0.0680301,-9.76308e-05,6.92402e-08,-1.93425e-11,-49848.3,25.9104], Tmin=(100,'K'), Tmax=(877.303,'K')), NASAPolynomial(coeffs=[12.4908,0.0161494,-8.9243e-06,1.83054e-09,-1.32805e-13,-51844.8,-27.4949], Tmin=(877.303,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-415.295,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-O2d)H) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-OdOsOs) + radical(OC=OOJ) + radical(C=CJO)"""),
)

species(
    label = '[O]C(F)=[C]OC(=O)O(3211)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 O u0 p2 c0 {6,S} {8,S}
3 O u0 p2 c0 {6,S} {9,S}
4 O u0 p2 c0 {6,D}
5 O u1 p2 c0 {7,S}
6 C u0 p0 c0 {2,S} {3,S} {4,D}
7 C u0 p0 c0 {1,S} {5,S} {8,D}
8 C u1 p0 c0 {2,S} {7,D}
9 H u0 p0 c0 {3,S}
"""),
    E0 = (-517.58,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,293,496,537,1218,1685,370,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43345,0.0608845,-7.28467e-05,4.2872e-08,-1.01345e-11,-62161.7,25.6921], Tmin=(100,'K'), Tmax=(1017.99,'K')), NASAPolynomial(coeffs=[12.1345,0.0188368,-1.08896e-05,2.29711e-09,-1.70002e-13,-64340.4,-26.1248], Tmin=(1017.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-517.58,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-O2d)H) + group(Cds-CdsOsH) + group(CdCFO) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-OdOsOs) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = '[O]C(=O)OC=[C]OF(3212)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 O u0 p2 c0 {6,S} {7,S}
3 O u0 p2 c0 {1,S} {8,S}
4 O u1 p2 c0 {7,S}
5 O u0 p2 c0 {7,D}
6 C u0 p0 c0 {2,S} {8,D} {9,S}
7 C u0 p0 c0 {2,S} {4,S} {5,D}
8 C u1 p0 c0 {3,S} {6,D}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-99.0766,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,3010,987.5,1337.5,450,1655,1685,370,180,180,180,297.832,834.067,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.148364,'amu*angstrom^2'), symmetry=1, barrier=(3.41119,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148364,'amu*angstrom^2'), symmetry=1, barrier=(3.41119,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148364,'amu*angstrom^2'), symmetry=1, barrier=(3.41119,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (120.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.894145,0.0673158,-8.77858e-05,5.39029e-08,-1.27061e-11,-11803.4,27.1331], Tmin=(100,'K'), Tmax=(1048.63,'K')), NASAPolynomial(coeffs=[16.1818,0.00900125,-4.37057e-06,8.71731e-10,-6.32045e-14,-15009.6,-47.3469], Tmin=(1048.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-99.0766,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2sCF) + group(O2s-(Cds-O2d)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdOsOs) + radical(OC=OOJ) + radical(C=CJO)"""),
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
    E0 = (-237.882,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (224.574,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (197.22,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-229.974,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-54.427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (132.068,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-81.306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-116.895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-76.8075,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-133.364,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-124.335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (84.8059,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-235.898,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (107.083,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (142.295,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-89.1284,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-47.848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (339.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-230.351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-51.2589,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (218.3,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (25.0756,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-126.378,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (240.215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(=O)O[CH]C(=O)F(3133)'],
    products = ['CO2(13)', 'O=CC(=O)F(557)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]C([O])=O(236)', '[CH]C(=O)F-2(2781)'],
    products = ['[O]C(=O)O[CH]C(=O)F(3133)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(131316,'m^3/(mol*s)'), n=0.920148, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_sec_rad;Birad] for rate rule [O_rad/OneDe;Birad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Birad_R_Recombination
Ea raised from -3.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['O(6)', 'O=[C]O[CH]C(=O)F(2655)'],
    products = ['[O]C(=O)O[CH]C(=O)F(3133)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [CO_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]C(=O)O[CH]C(=O)F(3133)'],
    products = ['O=C1OC(C(=O)F)O1(3178)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.6e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 2.23606797749979
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]C(=O)O[CH]C(=O)F(3133)'],
    products = ['[O]C(=O)OC1O[C]1F(3203)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(9.27213e+10,'s^-1'), n=0.543712, Ea=(183.455,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra_cs] for rate rule [R3_CO;carbonyl_intra;radadd_intra_csHO]
Euclidian distance = 2.449489742783178
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]C(=O)O[CH]C(=O)F(3133)'],
    products = ['O=C(F)[CH]O[C]1OO1(3204)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.55936e+11,'s^-1'), n=0.551275, Ea=(369.951,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra] for rate rule [R3_CO;carbonyl_intra_Nd;radadd_intra_O]
Euclidian distance = 2.449489742783178
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C(=O)O[CH]C(=O)F(3133)'],
    products = ['[O][C]1OC(C(=O)F)O1(3205)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.37972e+08,'s^-1'), n=1.13751, Ea=(156.576,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_csHCO]
Euclidian distance = 2.449489742783178
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 153.6 to 156.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C(=O)O[CH]C(=O)F(3133)'],
    products = ['O=C1O[CH][C](F)OO1(3143)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(6.65022e+10,'s^-1'), n=0.274726, Ea=(120.987,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra] for rate rule [R6_linear;carbonyl_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 118.3 to 121.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C(=O)O[CH]C(=O)F(3133)'],
    products = ['[O]C1([O])OC1C(=O)F(3206)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.90568e+10,'s^-1'), n=0.237, Ea=(161.075,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_csHDe] for rate rule [R4_S_CO;carbonylbond_intra;radadd_intra_csHDe]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic
Ea raised from 159.2 to 161.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C(=O)O[CH]C(=O)F(3133)'],
    products = ['[O]C1(F)[CH]OC(=O)O1(3169)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6.33774e+10,'s^-1'), n=0.492453, Ea=(104.519,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;multiplebond_intra;radadd_intra_O] + [R6;multiplebond_intra;radadd_intra] for rate rule [R6;carbonylbond_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O][C]=O(188)', 'O=CC(=O)F(557)'],
    products = ['[O]C(=O)O[CH]C(=O)F(3133)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(153.031,'m^3/(mol*s)'), n=1.16366, Ea=(37.2212,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.022037706214473284, var=2.3701416358838845, Tref=1000.0, N=230, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', '[O]C(=O)OC=C=O(3207)'],
    products = ['[O]C(=O)O[CH]C(=O)F(3133)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(37.0404,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CO2(13)', '[O][CH]C(=O)F(438)'],
    products = ['[O]C(=O)O[CH]C(=O)F(3133)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.30416e-16,'m^3/(mol*s)'), n=6.08142, Ea=(91.6915,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.42494241855262277, var=2.081425725945606, Tref=1000.0, N=380, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_N-Sp-2R!H#1R!H_Sp-4R!H-3R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O][C]=O(188)', '[O][CH]C(=O)F(438)'],
    products = ['[O]C(=O)O[CH]C(=O)F(3133)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O=[C][CH]OC(=O)OF(3208)'],
    products = ['[O]C(=O)O[CH]C(=O)F(3133)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.64853e+07,'s^-1'), n=1.15307, Ea=(80.7542,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R5nF',), comment="""Estimated from node R5nF"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]C(=O)O[CH]C(=O)F(3133)'],
    products = ['[O]C(=O)OC(F)[C]=O(3195)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.36622e+11,'s^-1'), n=0.48217, Ea=(148.754,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-2R!H-R',), comment="""Estimated from node R2F_Ext-2R!H-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]C(=O)C([O])(F)C=O(3132)'],
    products = ['[O]C(=O)O[CH]C(=O)F(3133)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(134.276,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['O(6)', '[O]C(=O)OC=[C]F(3209)'],
    products = ['[O]C(=O)O[CH]C(=O)F(3133)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_sec_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O]C(=O)O[CH]C(=O)F(3133)'],
    products = ['O=C1OC=C(F)OO1(3149)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Ypri_rad_out] for rate rule [R6_SSSDS;O_rad;Opri_rad]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O]C(=O)O[CH]C(=O)F(3133)'],
    products = ['[O]C1([O])OC=C(F)O1(3198)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.49749e+08,'s^-1'), n=0.656505, Ea=(186.623,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SMS;multiplebond_intra;radadd_intra] for rate rule [R6_SMS_CO;carbonylbond_intra;radadd_intra_O]
Euclidian distance = 1.7320508075688772
family: Intra_R_Add_Exocyclic
Ea raised from 182.9 to 186.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction21',
    reactants = ['HF(38)', '[O]C(=O)O[C]=C=O(3199)'],
    products = ['[O]C(=O)O[CH]C(=O)F(3133)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.64483e+40,'m^3/(mol*s)'), n=-10.004, Ea=(284.795,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_N-3COCdCddCtO2d->Cd"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O]C(=O)O[C]=C(O)F(3210)'],
    products = ['[O]C(=O)O[CH]C(=O)F(3133)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleNd;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[O]C(F)=[C]OC(=O)O(3211)'],
    products = ['[O]C(=O)O[CH]C(=O)F(3133)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]C(=O)OC=[C]OF(3212)'],
    products = ['[O]C(=O)O[CH]C(=O)F(3133)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(49.267,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

network(
    label = 'PDepNetwork #1194',
    isomers = [
        '[O]C(=O)O[CH]C(=O)F(3133)',
    ],
    reactants = [
        ('CO2(13)', 'O=CC(=O)F(557)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1194',
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

