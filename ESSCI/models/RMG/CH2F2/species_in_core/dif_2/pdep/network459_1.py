species(
    label = '[O]C(F)(F)C=[C]F(1052)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {7,S}
4 O u1 p2 c0 {5,S}
5 C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6 C u0 p0 c0 {5,S} {7,D} {8,S}
7 C u1 p0 c0 {3,S} {6,D}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-314.043,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([526,555,698,907,1200,1145,1227,3010,987.5,1337.5,450,1655,167,640,1190,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.17659,'amu*angstrom^2'), symmetry=1, barrier=(27.0522,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87986,0.0461386,-4.88669e-05,2.53519e-08,-5.18932e-12,-37693.9,20.3488], Tmin=(100,'K'), Tmax=(1182.34,'K')), NASAPolynomial(coeffs=[11.8453,0.0124248,-6.09565e-06,1.23551e-09,-9.00874e-14,-40050.4,-29.3979], Tmin=(1182.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-314.043,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-F1sF1sCd)) + radical(Cdj(Cd-CsH)(F1s))"""),
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
    label = 'C2HF(57)',
    structure = adjacencyList("""1 F u0 p3 c0 {3,S}
2 C u0 p0 c0 {3,T} {4,S}
3 C u0 p0 c0 {1,S} {2,T}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (95.331,'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(44.0062,'amu')),
        LinearRotor(inertia=(51.6236,'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([429.793,429.793,596.357,596.357,1107.96,2365.05,3506.88],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0277,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1870.76,'J/mol'), sigma=(4.25,'angstroms'), dipoleMoment=(1,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.4498,0.0030263,3.99146e-05,-8.9615e-08,5.74336e-11,11468.6,5.90915], Tmin=(10,'K'), Tmax=(555.749,'K')), NASAPolynomial(coeffs=[4.23833,0.0086714,-5.87678e-06,1.96876e-09,-2.53031e-13,11206.2,0.995297], Tmin=(555.749,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(95.331,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(87.302,'J/(mol*K)'), label="""C#CF""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'F[C]C=C(F)F(1096)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {6,S}
4 C u0 p0 c0 {5,D} {6,S} {7,S}
5 C u0 p0 c0 {1,S} {2,S} {4,D}
6 C u2 p0 c0 {3,S} {4,S}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-199.193,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,290.002,290.011,1499.69],'cm^-1')),
        HinderedRotor(inertia=(0.00200412,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.0351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.52272,0.0350047,-3.94019e-05,2.35736e-08,-5.77123e-12,-23906.3,16.3068], Tmin=(100,'K'), Tmax=(978.242,'K')), NASAPolynomial(coeffs=[7.71577,0.0137706,-6.84235e-06,1.3845e-09,-1.00585e-13,-24922.3,-8.6324], Tmin=(978.242,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-199.193,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCCl_triplet)"""),
)

species(
    label = 'FC1=CC(F)(F)O1(1093)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {7,S}
4 O u0 p2 c0 {5,S} {7,S}
5 C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6 C u0 p0 c0 {5,S} {7,D} {8,S}
7 C u0 p0 c0 {3,S} {4,S} {6,D}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-607.666,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.88726,0.00660468,0.000109124,-2.24183e-07,1.34132e-10,-73083,10.7267], Tmin=(10,'K'), Tmax=(574.596,'K')), NASAPolynomial(coeffs=[4.06958,0.0312823,-2.30324e-05,7.73795e-09,-9.66989e-13,-73532.3,6.22091], Tmin=(574.596,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-607.666,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), label="""FC1DCC(F)(F)O1""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'O=C(F)C=[C]F(1097)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {5,D}
4 C u0 p0 c0 {5,S} {6,D} {7,S}
5 C u0 p0 c0 {1,S} {3,D} {4,S}
6 C u1 p0 c0 {2,S} {4,D}
7 H u0 p0 c0 {4,S}
"""),
    E0 = (-278.413,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,255,533,799,832,1228,167,640,1190,264.238],'cm^-1')),
        HinderedRotor(inertia=(0.00241086,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (91.0361,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.77182,0.0239976,-9.59278e-06,-6.41263e-09,4.53523e-12,-33485.8,11.0457], Tmin=(10,'K'), Tmax=(939.953,'K')), NASAPolynomial(coeffs=[7.42042,0.0163751,-1.00426e-05,2.85286e-09,-3.08611e-13,-34520.8,-8.18818], Tmin=(939.953,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-278.413,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""ODC(F)CD[C]F""", comment="""Thermo library: CHOF_G4"""),
)

species(
    label = '[CH]=[C]F(530)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {3,S}
2 C u1 p0 c0 {3,D} {4,S}
3 C u1 p0 c0 {1,S} {2,D}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (350.447,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([167,640,1190,1142.58,1502.03,3807.5],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0277,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.55591,0.0076799,-5.30098e-07,-4.54651e-09,2.16658e-12,42166.8,9.46363], Tmin=(100,'K'), Tmax=(1043.42,'K')), NASAPolynomial(coeffs=[5.86818,0.00351556,-1.29996e-06,2.62251e-10,-1.98915e-14,41428.4,-3.016], Tmin=(1043.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(350.447,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Cds_P) + radical(CdCdF1s)"""),
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
    label = '[O]C(F)(F)C#CF(1098)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {7,S}
4 O u1 p2 c0 {5,S}
5 C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6 C u0 p0 c0 {5,S} {7,T}
7 C u0 p0 c0 {3,S} {6,T}
"""),
    E0 = (-285.758,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([170,411,408,596,588,668,1120,1103,2175,525,239,401,1367,280.403],'cm^-1')),
        HinderedRotor(inertia=(0.00214484,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.027,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.75826,0.0184656,9.41181e-05,-3.02626e-07,2.50749e-10,-34367.4,12.3782], Tmin=(10,'K'), Tmax=(448.533,'K')), NASAPolynomial(coeffs=[6.44741,0.0231744,-1.75763e-05,5.99733e-09,-7.55539e-13,-34897.2,-1.65647], Tmin=(448.533,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-285.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""[O]C(F)(F)C#CF""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'O=C(F)[C]=[C]F(1099)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {4,D}
4 C u0 p0 c0 {1,S} {3,D} {5,S}
5 C u1 p0 c0 {4,S} {6,D}
6 C u1 p0 c0 {2,S} {5,D}
"""),
    E0 = (-40.9128,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([102,183,524,698,842,1211,1685,370,167,640,1190],'cm^-1')),
        HinderedRotor(inertia=(0.0426076,'amu*angstrom^2'), symmetry=1, barrier=(13.2258,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (90.0281,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.31363,0.0408089,-5.82627e-05,4.14337e-08,-1.17122e-11,-4863.15,17.557], Tmin=(100,'K'), Tmax=(862.732,'K')), NASAPolynomial(coeffs=[8.74598,0.0109858,-6.41039e-06,1.36549e-09,-1.01374e-13,-5973.03,-12.5257], Tmin=(862.732,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-40.9128,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(C=CJC=O) + radical(Cdj(Cd-COH)(F1s))"""),
)

species(
    label = '[O]C(F)(F)[C]=CF(1100)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {6,S}
4 O u1 p2 c0 {5,S}
5 C u0 p0 c0 {1,S} {2,S} {4,S} {7,S}
6 C u0 p0 c0 {3,S} {7,D} {8,S}
7 C u1 p0 c0 {5,S} {6,D}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-300.868,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([526,555,698,907,1200,1145,1227,615,860,1140,1343,3152,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.813742,'amu*angstrom^2'), symmetry=1, barrier=(18.7095,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.97056,0.0474072,-5.58214e-05,3.38039e-08,-8.23917e-12,-36115.3,20.5538], Tmin=(100,'K'), Tmax=(990.627,'K')), NASAPolynomial(coeffs=[9.79057,0.0158305,-8.0071e-06,1.62531e-09,-1.18216e-13,-37664.6,-17.0994], Tmin=(990.627,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-300.868,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(O2sj(Cs-F1sF1sCd)) + radical(Cdj(Cs-F1sF1sO2s)(Cd-F1sH))"""),
)

species(
    label = 'OC(F)(F)[C]=[C]F(1101)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {7,S}
4 O u0 p2 c0 {5,S} {8,S}
5 C u0 p0 c0 {1,S} {2,S} {4,S} {6,S}
6 C u1 p0 c0 {5,S} {7,D}
7 C u1 p0 c0 {3,S} {6,D}
8 H u0 p0 c0 {4,S}
"""),
    E0 = (-298.765,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,140,424,499,621,667,843,876,1082,1685,370,167,640,1190],'cm^-1')),
        HinderedRotor(inertia=(0.345909,'amu*angstrom^2'), symmetry=1, barrier=(7.95312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.34914,'amu*angstrom^2'), symmetry=1, barrier=(8.02743,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31045,0.0669288,-0.000122917,1.12332e-07,-3.91095e-11,-35843.7,21.8916], Tmin=(100,'K'), Tmax=(841.849,'K')), NASAPolynomial(coeffs=[8.44442,0.018253,-9.85364e-06,1.94256e-09,-1.34748e-13,-36521.1,-8.18704], Tmin=(841.849,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-298.765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Cdj(Cs-F1sF1sO2s)(Cd-F1sH)) + radical(Cdj(Cd-CsH)(F1s))"""),
)

species(
    label = 'F[C]C=C(F)OF(1102)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {4,S}
4 O u0 p2 c0 {3,S} {6,S}
5 C u0 p0 c0 {6,D} {7,S} {8,S}
6 C u0 p0 c0 {1,S} {4,S} {5,D}
7 C u2 p0 c0 {2,S} {5,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-46.5992,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([231,791,3010,987.5,1337.5,450,1655,326,540,652,719,1357,198.295,198.379,198.443,2072.56],'cm^-1')),
        HinderedRotor(inertia=(1.0981,'amu*angstrom^2'), symmetry=1, barrier=(30.6759,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0100638,'amu*angstrom^2'), symmetry=1, barrier=(30.6753,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96866,0.0519184,-6.84054e-05,3.63522e-08,1.0995e-12,-5538.55,19.6608], Tmin=(100,'K'), Tmax=(533.506,'K')), NASAPolynomial(coeffs=[7.42393,0.0214366,-1.19984e-05,2.47367e-09,-1.79227e-13,-6268.92,-4.61995], Tmin=(533.506,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-46.5992,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCCl_triplet)"""),
)

species(
    label = 'O=C(F)[CH][C](F)F(464)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {6,S}
3 F u0 p3 c0 {7,S}
4 O u0 p2 c0 {7,D}
5 C u1 p0 c0 {6,S} {7,S} {8,S}
6 C u1 p0 c0 {1,S} {2,S} {5,S}
7 C u0 p0 c0 {3,S} {4,D} {5,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-515.422,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,190,488,555,1236,1407,611,648,830,1210,1753,329.065,1473.94],'cm^-1')),
        HinderedRotor(inertia=(0.144345,'amu*angstrom^2'), symmetry=1, barrier=(11.4621,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.687173,'amu*angstrom^2'), symmetry=1, barrier=(53.4543,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98659,0.0477077,-5.8135e-05,3.68603e-08,-9.44956e-12,-61921.3,21.5871], Tmin=(100,'K'), Tmax=(941.925,'K')), NASAPolynomial(coeffs=[9.32482,0.0165442,-8.50634e-06,1.73363e-09,-1.26212e-13,-63303.6,-13.3764], Tmin=(941.925,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-515.422,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCJC=O) + radical(Csj(Cs-COHH)(F1s)(F1s))"""),
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
    E0 = (-171.734,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (186.151,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-163.449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-32.6982,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (-88.6599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (2.38877,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (70.8835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (237.279,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (87.7417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (45.5785,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-62.7495,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (186.412,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-32.0077,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C(F)(F)C=[C]F(1052)'],
    products = ['CF2O(48)', 'C2HF(57)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O(6)', 'F[C]C=C(F)F(1096)'],
    products = ['[O]C(F)(F)C=[C]F(1052)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_ter_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]C(F)(F)C=[C]F(1052)'],
    products = ['FC1=CC(F)(F)O1(1093)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;O_rad;Cdsinglepri_rad_out]
Euclidian distance = 1.7320508075688772
family: Birad_recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['F(37)', 'O=C(F)C=[C]F(1097)'],
    products = ['[O]C(F)(F)C=[C]F(1052)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4.08261e+20,'m^3/(mol*s)'), n=-5.07836, Ea=(30.5133,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_2R!H->O"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CF2O(48)', '[CH]=[C]F(530)'],
    products = ['[O]C(F)(F)C=[C]F(1052)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(9.10216e-08,'m^3/(mol*s)'), n=3.71185, Ea=(37.1938,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.339286171633947, var=3.7349511333863634, Tref=1000.0, N=1489, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O][C](F)F(176)', 'C2HF(57)'],
    products = ['[O]C(F)(F)C=[C]F(1052)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.000105193,'m^3/(mol*s)'), n=3.55802, Ea=(20.2256,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-3C-R_2R!H->C_1R!H->C_N-4R!H->C',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-3C-R_2R!H->C_1R!H->C_N-4R!H->C"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(5)', '[O]C(F)(F)C#CF(1098)'],
    products = ['[O]C(F)(F)C=[C]F(1052)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2440,'m^3/(mol*s)'), n=1.64, Ea=(2.52728,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-4R!H-R_N-Sp-6R!H=4R!H',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_N-Sp-2CS=1CCOSS_N-5R!H-inRing_Ext-4R!H-R_N-Sp-6R!H=4R!H"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O][C](F)F(176)', '[CH]=[C]F(530)'],
    products = ['[O]C(F)(F)C=[C]F(1052)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.38316e+06,'m^3/(mol*s)'), n=1.31229e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.016021952005170214, var=0.3543710496450803, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C"""),
)

reaction(
    label = 'reaction9',
    reactants = ['HF(38)', 'O=C(F)[C]=[C]F(1099)'],
    products = ['[O]C(F)(F)C=[C]F(1052)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2676.63,'m^3/(mol*s)'), n=0.732206, Ea=(267.458,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.13214706218215122, var=14.38614219007748, Tref=1000.0, N=10, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C(F)(F)[C]=CF(1100)'],
    products = ['[O]C(F)(F)C=[C]F(1052)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_single]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['OC(F)(F)[C]=[C]F(1101)'],
    products = ['[O]C(F)(F)C=[C]F(1052)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.8074e+07,'s^-1'), n=1.23467, Ea=(93.7055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS_Cs;Y_rad_out;O_H_out] + [R3H_SS_Cs;Cd_rad_out;XH_out] for rate rule [R3H_SS_Cs;Cd_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F[C]C=C(F)OF(1102)'],
    products = ['[O]C(F)(F)C=[C]F(1052)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.69879e+07,'s^-1'), n=1.42748, Ea=(90.7018,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.05505882546177618, var=61.63242994454046, Tref=1000.0, N=11, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]C(F)(F)C=[C]F(1052)'],
    products = ['O=C(F)[CH][C](F)F(464)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(139.726,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

network(
    label = 'PDepNetwork #459',
    isomers = [
        '[O]C(F)(F)C=[C]F(1052)',
    ],
    reactants = [
        ('CF2O(48)', 'C2HF(57)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #459',
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

