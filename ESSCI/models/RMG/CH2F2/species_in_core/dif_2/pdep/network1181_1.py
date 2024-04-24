species(
    label = '[O]C(F)(F)C([O])(F)C=O(3006)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u1 p2 c0 {7,S}
5  O u1 p2 c0 {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
8  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
9  C u0 p0 c0 {6,D} {7,S} {10,S}
10 H u0 p0 c0 {9,S}
"""),
    E0 = (-703.023,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,351,323,533,609,664,892,1120,1201,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.927134,'amu*angstrom^2'), symmetry=1, barrier=(21.3166,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.92739,'amu*angstrom^2'), symmetry=1, barrier=(21.3225,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.615694,0.0817469,-0.000126536,9.62049e-08,-2.71992e-11,-84439,27.3639], Tmin=(100,'K'), Tmax=(655.573,'K')), NASAPolynomial(coeffs=[12.0905,0.0216189,-1.15791e-05,2.30506e-09,-1.62853e-13,-86156,-24.7705], Tmin=(655.573,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-703.023,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(O2sj(Cs-CsF1sF1s))"""),
)

species(
    label = 'CF2O(48)',
    structure = adjacencyList("""1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {1,S} {2,S} {3,D}
"""),
    E0 = (-618.61,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.06775,-0.00614966,6.73615e-05,-1.17623e-07,6.56735e-11,-74400.6,7.68563], Tmin=(10,'K'), Tmax=(584.627,'K')), NASAPolynomial(coeffs=[3.15981,0.0116963,-8.27581e-06,2.66621e-09,-3.20585e-13,-74493.3,9.87819], Tmin=(584.627,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-618.61,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""ODC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = '[O]C(F)C([O])(F)F(348)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {7,S}
4 O u1 p2 c0 {6,S}
5 O u1 p2 c0 {7,S}
6 C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
7 C u0 p0 c0 {3,S} {5,S} {6,S} {8,S}
8 H u0 p0 c0 {7,S}
"""),
    E0 = (-613.996,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([351,323,533,609,664,892,1120,1201,391,562,707,872,1109,1210,1289,3137,180],'cm^-1')),
        HinderedRotor(inertia=(1.5766,'amu*angstrom^2'), symmetry=1, barrier=(36.249,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.023,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3807.04,'J/mol'), sigma=(6.24717,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=594.65 K, Pc=35.43 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41151,0.0470069,-3.44331e-05,1.94681e-09,4.4207e-12,-73744.7,20.6885], Tmin=(100,'K'), Tmax=(1020.07,'K')), NASAPolynomial(coeffs=[17.415,0.00468127,-2.23431e-06,5.36432e-10,-4.5404e-14,-78072.5,-62.0468], Tmin=(1020.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-613.996,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-CsF1sF1s)) + radical(O2sj(Cs-CsF1sH))"""),
)

species(
    label = '[O]C(=COF)C([O])(F)F(3092)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {9,S}
5  O u1 p2 c0 {7,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8  C u0 p0 c0 {6,S} {7,S} {9,D}
9  C u0 p0 c0 {4,S} {8,D} {10,S}
10 H u0 p0 c0 {9,S}
"""),
    E0 = (-472.943,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,526,555,698,907,1200,1145,1227,350,440,435,1725,3010,987.5,1337.5,450,1655,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.975029,'amu*angstrom^2'), symmetry=1, barrier=(22.4178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.97269,'amu*angstrom^2'), symmetry=1, barrier=(22.3641,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.315376,0.0827449,-0.000116667,7.91567e-08,-2.07276e-11,-56750.6,27.8681], Tmin=(100,'K'), Tmax=(943.985,'K')), NASAPolynomial(coeffs=[16.7049,0.0132945,-6.30695e-06,1.215e-09,-8.53789e-14,-59844.8,-50.2568], Tmin=(943.985,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-472.943,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2sCF) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(O2sj(Cs-F1sF1sCd)) + radical(C=C(C)OJ)"""),
)

species(
    label = '[O]C(F)(F)O[CH]C(=O)F(3007)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
8  C u1 p0 c0 {4,S} {9,S} {10,S}
9  C u0 p0 c0 {3,S} {6,D} {8,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-769.353,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,486,586,623,710,938,1161,1294,3025,407.5,1350,352.5,611,648,830,1210,1753,180,180,234.34,741.902],'cm^-1')),
        HinderedRotor(inertia=(0.0882243,'amu*angstrom^2'), symmetry=1, barrier=(2.02845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.794674,'amu*angstrom^2'), symmetry=1, barrier=(18.2711,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.794545,'amu*angstrom^2'), symmetry=1, barrier=(18.2682,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3928.52,'J/mol'), sigma=(6.25769,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=613.63 K, Pc=36.38 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.576176,0.0770894,-0.000102379,6.53e-08,-1.6183e-11,-92409.9,28.0187], Tmin=(100,'K'), Tmax=(992.174,'K')), NASAPolynomial(coeffs=[16.1523,0.0142934,-7.44255e-06,1.50962e-09,-1.09643e-13,-95500.8,-47.0049], Tmin=(992.174,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-769.353,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)OsHH) + group(CsFFOO) + group(COCsFO) + radical(O2sj(Cs-F1sF1sO2s)) + radical(CCsJOCs)"""),
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
    label = '[O]C=C(F)C([O])(F)F(2673)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {7,S}
4 O u1 p2 c0 {6,S}
5 O u1 p2 c0 {8,S}
6 C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
7 C u0 p0 c0 {3,S} {6,S} {8,D}
8 C u0 p0 c0 {5,S} {7,D} {9,S}
9 H u0 p0 c0 {8,S}
"""),
    E0 = (-598.505,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([526,555,698,907,1200,1145,1227,323,467,575,827,1418,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.49432,'amu*angstrom^2'), symmetry=1, barrier=(34.3575,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.969481,0.0605051,-6.66811e-05,3.46144e-08,-6.89829e-12,-71869.2,24.6419], Tmin=(100,'K'), Tmax=(1237.02,'K')), NASAPolynomial(coeffs=[16.9965,0.00868079,-3.83962e-06,7.47408e-10,-5.38492e-14,-75834.3,-56.0882], Tmin=(1237.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-598.505,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cds(F)) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cd(F)=CdOs) + group(Cds-CdsOsH) + radical(O2sj(Cs-F1sF1sCd)) + radical(C=COJ)"""),
)

species(
    label = '[O]C(F)(C=O)[C](F)F(3093)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {7,S}
4 O u1 p2 c0 {6,S}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
7 C u1 p0 c0 {2,S} {3,S} {6,S}
8 C u0 p0 c0 {5,D} {6,S} {9,S}
9 H u0 p0 c0 {8,S}
"""),
    E0 = (-550.739,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,190,488,555,1236,1407,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.1643,'amu*angstrom^2'), symmetry=1, barrier=(26.7696,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16288,'amu*angstrom^2'), symmetry=1, barrier=(26.7368,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (126.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.847262,0.0763394,-0.000131578,1.15083e-07,-3.91043e-11,-66131.8,25.4981], Tmin=(100,'K'), Tmax=(821.506,'K')), NASAPolynomial(coeffs=[10.2223,0.0205074,-1.10383e-05,2.17222e-09,-1.51193e-13,-67328.5,-15.796], Tmin=(821.506,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-550.739,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CsCsF1sF1s)"""),
)

species(
    label = 'O=CC1(F)OOC1(F)F(3013)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
8  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
9  C u0 p0 c0 {6,D} {7,S} {10,S}
10 H u0 p0 c0 {9,S}
"""),
    E0 = (-811.021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32027,0.0564643,-5.24432e-05,2.33494e-08,-4.12478e-12,-97445,23.0906], Tmin=(100,'K'), Tmax=(1350.52,'K')), NASAPolynomial(coeffs=[14.623,0.0170639,-8.68173e-06,1.74705e-09,-1.25879e-13,-101038,-45.0845], Tmin=(1350.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-811.021,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCCFO) + group(CsCFFO) + group(Cds-OdCsH) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + ring(O2s-Cs(F)(F)-Cs-O2s)"""),
)

species(
    label = '[O]C(F)(F)C1(F)[CH]OO1(3094)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {7,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
8  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
9  C u1 p0 c0 {5,S} {7,S} {10,S}
10 H u0 p0 c0 {9,S}
"""),
    E0 = (-446.004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12765,0.0678858,-8.58445e-05,5.60956e-08,-1.47834e-11,-53542.5,24.9001], Tmin=(100,'K'), Tmax=(918.87,'K')), NASAPolynomial(coeffs=[11.5418,0.0225503,-1.18355e-05,2.39882e-09,-1.7361e-13,-55456.3,-24.4607], Tmin=(918.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-446.004,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsOsHH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + ring(Cs-Cs(F)-O2s-O2s) + radical(O2sj(Cs-CsF1sF1s)) + radical(CCsJOO)"""),
)

species(
    label = '[O]C1(F)[CH]OOC1(F)F(3031)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {8,S}
5  O u0 p2 c0 {4,S} {9,S}
6  O u1 p2 c0 {7,S}
7  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
8  C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
9  C u1 p0 c0 {5,S} {7,S} {10,S}
10 H u0 p0 c0 {9,S}
"""),
    E0 = (-557.158,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.376114,0.0693917,-7.98967e-05,4.49545e-08,-9.47733e-12,-66871.3,24.8422], Tmin=(100,'K'), Tmax=(1324.03,'K')), NASAPolynomial(coeffs=[17.4435,0.00868746,-7.67217e-07,-1.03396e-10,1.51008e-14,-70589.5,-59.2624], Tmin=(1324.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-557.158,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(CsCCFO) + group(CsCFFO) + group(Cs-CsOsHH) + ring(12dioxolane) + radical(O2sj(Cs-F1sCsCs)) + radical(CCsJOOC) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F))"""),
)

species(
    label = '[O]C1OC1(F)C([O])(F)F(3095)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {7,S} {8,S}
5  O u1 p2 c0 {8,S}
6  O u1 p2 c0 {9,S}
7  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
8  C u0 p0 c0 {4,S} {5,S} {7,S} {10,S}
9  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-641.057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31563,0.0656018,-8.56919e-05,6.23084e-08,-1.89352e-11,-77010.5,24.163], Tmin=(100,'K'), Tmax=(789.995,'K')), NASAPolynomial(coeffs=[8.59993,0.0287191,-1.56612e-05,3.21053e-09,-2.33264e-13,-78161.5,-9.26248], Tmin=(789.995,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-641.057,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsOsOsH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + ring(Cs(O2)-O2s-Cs(F)) + radical(CCOJ) + radical(O2sj(Cs-CsF1sF1s))"""),
)

species(
    label = '[O]C1OC(F)(F)C1([O])F(3059)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  O u0 p2 c0 {8,S} {9,S}
5  O u1 p2 c0 {7,S}
6  O u1 p2 c0 {8,S}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
8  C u0 p0 c0 {4,S} {6,S} {7,S} {10,S}
9  C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-712.026,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56519,0.0595324,-7.49928e-05,5.46227e-08,-1.69057e-11,-85554.6,25.1751], Tmin=(100,'K'), Tmax=(771.81,'K')), NASAPolynomial(coeffs=[7.48323,0.0288609,-1.53823e-05,3.13214e-09,-2.26899e-13,-86468.1,-1.84309], Tmin=(771.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-712.026,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(CsCCFO) + group(Cs-CsOsOsH) + group(CsCFFO) + ring(Cs-Cs(F)(F)-O2s-Cs) + radical(O2sj(Cs-F1sCsCs)) + radical(CCOJ) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F))"""),
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
    label = '[O]C(F)(F)C(=O)C=O(3096)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 O u1 p2 c0 {6,S}
4 O u0 p2 c0 {7,D}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
7 C u0 p0 c0 {4,D} {6,S} {8,S}
8 C u0 p0 c0 {5,D} {7,S} {9,S}
9 H u0 p0 c0 {8,S}
"""),
    E0 = (-616.428,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([377,639,914,1335,375,552.5,462.5,1710,2782.5,750,1395,475,1775,1000,180,180,994.119,995.197,996.096],'cm^-1')),
        HinderedRotor(inertia=(0.241639,'amu*angstrom^2'), symmetry=1, barrier=(5.55575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.241549,'amu*angstrom^2'), symmetry=1, barrier=(5.55368,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45267,0.064515,-0.000112719,1.07737e-07,-3.99257e-11,-74055.7,24.4368], Tmin=(100,'K'), Tmax=(817.993,'K')), NASAPolynomial(coeffs=[5.43764,0.0281624,-1.51283e-05,2.99453e-09,-2.09763e-13,-74143.4,9.46129], Tmin=(817.993,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-616.428,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(O2sj(Cs-F1sF1sCO))"""),
)

species(
    label = '[O][C](F)F(176)',
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
    label = 'HCO(14)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,D}
2 C u1 p0 c0 {1,D} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (32.4782,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1131.19,1955.83,1955.83],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.23755,-0.00332075,1.4003e-05,-1.3424e-08,4.37416e-12,3872.41,3.30835], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[3.92002,0.00252279,-6.71004e-07,1.05616e-10,-7.43798e-15,3653.43,3.58077], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(32.4782,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HCO""", comment="""Thermo library: FFCM1(-)"""),
)

species(
    label = '[O]C(F)(F)C(=O)F(304)',
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
    label = '[O]C(F)(C=O)C(=O)F(3097)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {8,S}
3 O u1 p2 c0 {6,S}
4 O u0 p2 c0 {7,D}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
7 C u0 p0 c0 {4,D} {6,S} {9,S}
8 C u0 p0 c0 {2,S} {5,D} {6,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-686.16,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,486,617,768,1157,1926,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.16676,'amu*angstrom^2'), symmetry=1, barrier=(26.8262,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17597,'amu*angstrom^2'), symmetry=1, barrier=(27.0379,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14244,0.0688121,-0.000113218,9.80766e-08,-3.35513e-11,-82428.8,24.399], Tmin=(100,'K'), Tmax=(796.52,'K')), NASAPolynomial(coeffs=[9.25056,0.0213979,-1.13173e-05,2.23316e-09,-1.56519e-13,-83508.1,-11.5399], Tmin=(796.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-686.16,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cds-OdCsH) + group(COCsFO) + radical(C=OCOJ)"""),
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
    label = 'F2(106)',
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
    label = '[O]C=C([O])C(=O)F(3098)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 O u1 p2 c0 {5,S}
3 O u1 p2 c0 {6,S}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {2,S} {6,D} {7,S}
6 C u0 p0 c0 {3,S} {5,D} {8,S}
7 C u0 p0 c0 {1,S} {4,D} {5,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-350.208,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,255,533,799,832,1228,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.98945,'amu*angstrom^2'), symmetry=1, barrier=(22.7494,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (104.036,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.376829,0.062269,-7.56558e-05,4.02788e-08,-7.73099e-12,-41974.9,23.536], Tmin=(100,'K'), Tmax=(1514.91,'K')), NASAPolynomial(coeffs=[20.3496,-0.00437702,4.1067e-06,-8.8305e-10,6.19861e-14,-46430.2,-75.8492], Tmin=(1514.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-350.208,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsOsH) + group(COCFO) + radical(C=OCOJ) + radical(C=COJ)"""),
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
    label = '[O]C(F)(F)C(=O)[C]=O(3099)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 O u1 p2 c0 {6,S}
4 O u0 p2 c0 {7,D}
5 O u0 p2 c0 {8,D}
6 C u0 p0 c0 {1,S} {2,S} {3,S} {7,S}
7 C u0 p0 c0 {4,D} {6,S} {8,S}
8 C u1 p0 c0 {5,D} {7,S}
"""),
    E0 = (-456.467,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([377,639,914,1335,375,552.5,462.5,1710,1855,455,950,180,180,180,1270.81,1271.28],'cm^-1')),
        HinderedRotor(inertia=(0.193326,'amu*angstrom^2'), symmetry=1, barrier=(4.44495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193738,'amu*angstrom^2'), symmetry=1, barrier=(4.45442,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.027,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36504,0.0676227,-0.000128037,1.21855e-07,-4.37658e-11,-54814.8,25.5714], Tmin=(100,'K'), Tmax=(849.84,'K')), NASAPolynomial(coeffs=[6.4264,0.0229974,-1.25541e-05,2.45942e-09,-1.69678e-13,-54923.8,6.39659], Tmin=(849.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-456.467,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(O2sj(Cs-F1sF1sCO)) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C(F)(F)C(O)(F)[C]=O(3100)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {7,S} {10,S}
5  O u1 p2 c0 {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
8  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
9  C u1 p0 c0 {6,D} {7,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (-786.795,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.199616,0.0891494,-0.000144861,1.15041e-07,-3.55097e-11,-94497.9,28.8779], Tmin=(100,'K'), Tmax=(800.605,'K')), NASAPolynomial(coeffs=[14.7248,0.0165613,-8.82917e-06,1.74025e-09,-1.21765e-13,-96823.2,-37.9641], Tmin=(800.605,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-786.795,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-OdCsH) + radical(O2sj(Cs-CsF1sF1s)) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C(F)([C]=O)C(O)(F)F(3101)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {8,S} {10,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
8  C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
9  C u1 p0 c0 {6,D} {7,S}
10 H u0 p0 c0 {4,S}
"""),
    E0 = (-803.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.121845,0.0915203,-0.000154833,1.27903e-07,-4.08547e-11,-96448.4,29.3435], Tmin=(100,'K'), Tmax=(826.127,'K')), NASAPolynomial(coeffs=[14.3342,0.0169425,-8.95761e-06,1.7398e-09,-1.19854e-13,-98600,-35.3177], Tmin=(826.127,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-803.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C=C(OF)C([O])(F)F(3102)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u1 p2 c0 {7,S}
6  O u1 p2 c0 {9,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {8,S}
8  C u0 p0 c0 {4,S} {7,S} {9,D}
9  C u0 p0 c0 {6,S} {8,D} {10,S}
10 H u0 p0 c0 {9,S}
"""),
    E0 = (-469.285,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,526,555,698,907,1200,1145,1227,350,440,435,1725,3010,987.5,1337.5,450,1655,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.28087,'amu*angstrom^2'), symmetry=1, barrier=(29.4497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.28385,'amu*angstrom^2'), symmetry=1, barrier=(29.5183,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.448935,0.0786983,-0.000102007,6.32767e-08,-1.52117e-11,-56314.4,27.2982], Tmin=(100,'K'), Tmax=(1023.62,'K')), NASAPolynomial(coeffs=[17.0378,0.0138744,-7.01572e-06,1.41076e-09,-1.02217e-13,-59710.6,-53.121], Tmin=(1023.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-469.285,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(O2s-(Cds-Cd)H) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(O2sj(Cs-F1sF1sCd)) + radical(C=COJ)"""),
)

species(
    label = '[O]C=C([O])C(F)(F)OF(3103)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {7,S}
5  O u1 p2 c0 {8,S}
6  O u1 p2 c0 {9,S}
7  C u0 p0 c0 {1,S} {2,S} {4,S} {8,S}
8  C u0 p0 c0 {5,S} {7,S} {9,D}
9  C u0 p0 c0 {6,S} {8,D} {10,S}
10 H u0 p0 c0 {9,S}
"""),
    E0 = (-602.792,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,275,321,533,585,746,850,1103,350,440,435,1725,3010,987.5,1337.5,450,1655,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.975033,'amu*angstrom^2'), symmetry=1, barrier=(22.4179,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.977107,'amu*angstrom^2'), symmetry=1, barrier=(22.4656,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.048164,0.0891519,-0.000129822,8.94979e-08,-2.37143e-11,-72358.7,26.4335], Tmin=(100,'K'), Tmax=(934.255,'K')), NASAPolynomial(coeffs=[18.0467,0.0120881,-6.08669e-06,1.19909e-09,-8.51705e-14,-75721.6,-59.1745], Tmin=(934.255,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-602.792,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(CsCFFO) + longDistanceInteraction_noncyclic(Cs(F)2-CdOs) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = '[O][C](F)C(F)(C=O)OF(3104)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {7,S}
5  O u1 p2 c0 {8,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
8  C u1 p0 c0 {2,S} {5,S} {7,S}
9  C u0 p0 c0 {6,D} {7,S} {10,S}
10 H u0 p0 c0 {9,S}
"""),
    E0 = (-421.6,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,395,473,707,1436,2782.5,750,1395,475,1775,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.21413,'amu*angstrom^2'), symmetry=1, barrier=(27.9153,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00124,'amu*angstrom^2'), symmetry=1, barrier=(23.0205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.995654,'amu*angstrom^2'), symmetry=1, barrier=(22.892,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.146992,0.0942379,-0.000168766,1.49903e-07,-5.11392e-11,-50577.1,29.2238], Tmin=(100,'K'), Tmax=(836.836,'K')), NASAPolynomial(coeffs=[11.5212,0.0235964,-1.29736e-05,2.55263e-09,-1.76829e-13,-51911,-20.2197], Tmin=(836.836,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-421.6,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(O2s-CsH) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-OdCsH) + radical(O2sj(Cs-CsF1sH)) + radical(CsCsF1sO2s)"""),
)

species(
    label = '[O]C(F)(C=O)[C](F)OF(3105)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {4,S}
4  O u0 p2 c0 {3,S} {8,S}
5  O u1 p2 c0 {7,S}
6  O u0 p2 c0 {9,D}
7  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
8  C u1 p0 c0 {2,S} {4,S} {7,S}
9  C u0 p0 c0 {6,D} {7,S} {10,S}
10 H u0 p0 c0 {9,S}
"""),
    E0 = (-404.44,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([557,1111,1380,1390,370,380,2900,435,395,473,707,1436,2782.5,750,1395,475,1775,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.07823,'amu*angstrom^2'), symmetry=1, barrier=(24.7906,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08198,'amu*angstrom^2'), symmetry=1, barrier=(24.8769,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.09008,'amu*angstrom^2'), symmetry=1, barrier=(25.063,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.033,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.111412,0.0955988,-0.000173158,1.55292e-07,-5.34056e-11,-48512.6,29.2285], Tmin=(100,'K'), Tmax=(835.881,'K')), NASAPolynomial(coeffs=[11.2976,0.0243697,-1.35748e-05,2.68309e-09,-1.86291e-13,-49764.3,-19.0343], Tmin=(835.881,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-404.44,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2sCF) + group(CsCCFO) + longDistanceInteraction_noncyclic(Cs(F)-CO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFHO) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CsCsF1sO2s)"""),
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
    E0 = (-280.798,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (37.1066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (155.114,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-151.249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (66.7541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (114.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-272.513,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-23.7795,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-134.933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-178.507,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-158.755,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-103.129,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-261.707,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-197.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (-154.468,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-264.767,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-41.7918,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (91.6407,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-33.4474,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-205.517,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-153.404,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (54.8008,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-0.415286,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (111.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (99.0938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(F)(F)C([O])(F)C=O(3006)'],
    products = ['CF2O(48)', 'O=CC(=O)F(557)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CO(12)', '[O]C(F)C([O])(F)F(348)'],
    products = ['[O]C(F)(F)C([O])(F)C=O(3006)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(88.9,'m^3/(mol*s)'), n=1.51, Ea=(347.619,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_N-2Br1sCl1sF1sH->F1s_5R!H->C_Ext-1Cs-R_Ext-1Cs-R',), comment="""Estimated from node Root_1COCbCdCsCtHNOSSidSis->Cs_N-2Br1sCbCdCl1sCsCtF1sHI1sNSSidSis->Cs_Ext-1Cs-R_N-2Br1sCl1sF1sH->F1s_5R!H->C_Ext-1Cs-R_Ext-1Cs-R"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]C(=COF)C([O])(F)F(3092)'],
    products = ['[O]C(F)(F)C([O])(F)C=O(3006)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(205.832,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]C(F)(F)C([O])(F)C=O(3006)'],
    products = ['[O]C(F)(F)O[CH]C(=O)F(3007)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(9.39365e+11,'s^-1'), n=0.324012, Ea=(129.549,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.014478493324023197, var=15.997960675483611, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_N-1R!H-inRing_Ext-4R!H-R',), comment="""Estimated from node Root_N-1R!H-inRing_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(6)', '[O]C=C(F)C([O])(F)F(2673)'],
    products = ['[O]C(F)(F)C([O])(F)C=O(3006)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_ter_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['O(6)', '[O]C(F)(C=O)[C](F)F(3093)'],
    products = ['[O]C(F)(F)C([O])(F)C=O(3006)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_ter_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C(F)(F)C([O])(F)C=O(3006)'],
    products = ['O=CC1(F)OOC1(F)F(3013)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C(F)(F)C([O])(F)C=O(3006)'],
    products = ['[O]C(F)(F)C1(F)[CH]OO1(3094)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(257.018,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 255.4 to 257.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C(F)(F)C([O])(F)C=O(3006)'],
    products = ['[O]C1(F)[CH]OOC1(F)F(3031)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.12332e+09,'s^-1'), n=0.5388, Ea=(145.864,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra] for rate rule [R5_SS_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 2.449489742783178
family: Intra_R_Add_Endocyclic
Ea raised from 143.6 to 145.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C(F)(F)C([O])(F)C=O(3006)'],
    products = ['[O]C1OC1(F)C([O])(F)F(3095)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.785e+11,'s^-1'), n=0.342, Ea=(102.29,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]C(F)(F)C([O])(F)C=O(3006)'],
    products = ['[O]C1OC(F)(F)C1([O])F(3059)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_O] for rate rule [R5_SS_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', '[O]C(F)(F)C(=O)C=O(3096)'],
    products = ['[O]C(F)(F)C([O])(F)C=O(3006)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.08261e+20,'m^3/(mol*s)'), n=-5.07836, Ea=(18.1825,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O][C](F)F(176)', 'O=CC(=O)F(557)'],
    products = ['[O]C(F)(F)C([O])(F)C=O(3006)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.4291e-06,'m^3/(mol*s)'), n=3.20779, Ea=(54.6622,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing_Ext-1R!H-R_N-2R!H->C_N-5R!H-u1',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing_Ext-1R!H-R_N-2R!H->C_N-5R!H-u1"""),
)

reaction(
    label = 'reaction14',
    reactants = ['HCO(14)', '[O]C(F)(F)C(=O)F(304)'],
    products = ['[O]C(F)(F)C([O])(F)C=O(3006)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=-1.07934e-11, Ea=(111.087,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H_N-2R!H->C"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F(37)', '[O]C(F)(C=O)C(=O)F(3097)'],
    products = ['[O]C(F)(F)C([O])(F)C=O(3006)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.08261e+20,'m^3/(mol*s)'), n=-5.07836, Ea=(36.5757,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O"""),
)

reaction(
    label = 'reaction16',
    reactants = ['CF2O(48)', '[O][CH]C(=O)F(438)'],
    products = ['[O]C(F)(F)C([O])(F)C=O(3006)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(9.07578e-06,'m^3/(mol*s)'), n=3.04336, Ea=(146.095,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.3757377757886876, var=2.242054186761003, Tref=1000.0, N=1042, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O][C](F)F(176)', '[O][CH]C(=O)F(438)'],
    products = ['[O]C(F)(F)C([O])(F)C=O(3006)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(5.93777,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction18',
    reactants = ['F2(106)', '[O]C=C([O])C(=O)F(3098)'],
    products = ['[O]C(F)(F)C([O])(F)C=O(3006)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.000118654,'m^3/(mol*s)'), n=2.63647, Ea=(28.4291,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='YY',), comment="""Estimated from node YY
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction19',
    reactants = ['HF(38)', '[O]C(F)(F)C(=O)[C]=O(3099)'],
    products = ['[O]C(F)(F)C([O])(F)C=O(3006)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(281.909,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O]C(F)(F)C([O])(F)C=O(3006)'],
    products = ['[O]C(F)(F)C(O)(F)[C]=O(3100)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;O_rad_out;XH_out] for rate rule [R3H_SS_Cs;O_rad_out;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]C(F)(F)C([O])(F)C=O(3006)'],
    products = ['[O]C(F)([C]=O)C(O)(F)F(3101)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.75172e+06,'s^-1'), n=1.80068, Ea=(127.394,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;O_rad_out;XH_out] for rate rule [R4H_SSS;O_rad_out;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O]C=C(OF)C([O])(F)F(3102)'],
    products = ['[O]C(F)(F)C([O])(F)C=O(3006)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(101.861,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[O]C=C([O])C(F)(F)OF(3103)'],
    products = ['[O]C(F)(F)C([O])(F)C=O(3006)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(180.152,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O][C](F)C(F)(C=O)OF(3104)'],
    products = ['[O]C(F)(F)C([O])(F)C=O(3006)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(111.115,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O]C(F)(C=O)[C](F)OF(3105)'],
    products = ['[O]C(F)(F)C([O])(F)C=O(3006)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(81.3093,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

network(
    label = 'PDepNetwork #1181',
    isomers = [
        '[O]C(F)(F)C([O])(F)C=O(3006)',
    ],
    reactants = [
        ('CF2O(48)', 'O=CC(=O)F(557)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1181',
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

