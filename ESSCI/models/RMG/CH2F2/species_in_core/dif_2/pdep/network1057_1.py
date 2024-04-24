species(
    label = 'F[C](F)C(F)C1O[C]1F(2762)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {8,S}
6  C u0 p0 c0 {5,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {1,S} {6,S} {9,S} {11,S}
8  C u1 p0 c0 {2,S} {5,S} {6,S}
9  C u1 p0 c0 {3,S} {4,S} {7,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-517.543,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,259,529,569,1128,1321,1390,3140,190,488,555,1236,1407,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07832,0.070088,-9.12071e-05,6.53218e-08,-1.92646e-11,-62146,27.1611], Tmin=(100,'K'), Tmax=(818.676,'K')), NASAPolynomial(coeffs=[9.63898,0.0282623,-1.4575e-05,2.92026e-09,-2.09493e-13,-63547.7,-12.4266], Tmin=(818.676,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-517.543,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCFHO) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + ring(O2s-Cs(F)-Cs(C)) + radical(Csj(Cs-O2sCsH)(F1s)(O2s-Cs)_ring) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
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
    E0 = (-511.455,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93201,0.00413887,7.03233e-05,-1.49827e-07,9.42397e-11,-61509.2,9.50863], Tmin=(10,'K'), Tmax=(540.182,'K')), NASAPolynomial(coeffs=[3.82032,0.0197002,-1.38027e-05,4.49179e-09,-5.49698e-13,-61712.1,7.9889], Tmin=(540.182,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-511.455,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""FCDC(F)F""", comment="""Thermo library: CHOF_G4"""),
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
    label = 'F[C]1OC([C](F)F)C1F(2763)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {8,S}
6  C u0 p0 c0 {5,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {1,S} {6,S} {8,S} {11,S}
8  C u1 p0 c0 {2,S} {5,S} {7,S}
9  C u1 p0 c0 {3,S} {4,S} {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-566.252,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26922,0.0626859,-6.49342e-05,3.4404e-08,-7.42028e-12,-68008.2,25.1933], Tmin=(100,'K'), Tmax=(1104.21,'K')), NASAPolynomial(coeffs=[12.0438,0.0236552,-1.19137e-05,2.39293e-09,-1.72806e-13,-70387.7,-27.8558], Tmin=(1104.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-566.252,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCFHO) + group(CsCsFFH) + ring(O2s-Cs-Cs-Cs(F)) + radical(Csj(Cs-F1sCsH)(F1s)(O2s-Cs)_ring) + radical(Csj(Cs-O2sCsH)(F1s)(F1s)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[C](F)C(F)C1(F)[CH]O1(2761)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {8,S}
6  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
7  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
8  C u1 p0 c0 {5,S} {6,S} {11,S}
9  C u1 p0 c0 {3,S} {4,S} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-525.174,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.897728,0.0740384,-0.000102579,7.75054e-08,-2.38371e-11,-63057.4,27.1737], Tmin=(100,'K'), Tmax=(790.288,'K')), NASAPolynomial(coeffs=[10.1169,0.0273741,-1.40044e-05,2.78238e-09,-1.98123e-13,-64514.5,-15.1332], Tmin=(790.288,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-525.174,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCFO) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cs-CsOsHH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + ring(Cs(F)(C)-Cs-O2s) + radical(Csj(Cs-F1sO2sCs)(O2s-Cs)(H)_ring) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
)

species(
    label = 'F[CH]C(F)(F)C1O[C]1F(2760)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {8,S}
6  C u0 p0 c0 {5,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {6,S} {9,S}
8  C u1 p0 c0 {3,S} {5,S} {6,S}
9  C u1 p0 c0 {4,S} {7,S} {11,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-530.132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.69488,0.0791066,-0.000120582,1.01033e-07,-3.43069e-11,-63647.2,27.3241], Tmin=(100,'K'), Tmax=(738.886,'K')), NASAPolynomial(coeffs=[9.7708,0.0285018,-1.48626e-05,2.95019e-09,-2.08838e-13,-64948.2,-13.4436], Tmin=(738.886,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-530.132,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCFHO) + group(CsCsFHH) + ring(O2s-Cs(F)-Cs(C)) + radical(Csj(Cs-O2sCsH)(F1s)(O2s-Cs)_ring) + radical(Csj(Cs-CsF1sF1s)(F1s)(H))"""),
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
    label = 'F[CH]C1O[C]1F(2764)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u0 p2 c0 {4,S} {5,S}
4 C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
5 C u1 p0 c0 {1,S} {3,S} {4,S}
6 C u1 p0 c0 {2,S} {4,S} {8,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {6,S}
"""),
    E0 = (-93.3689,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,334,575,1197,1424,3202,180,180,180,180,890.941,891.174,891.327,891.928,892.002,896.304],'cm^-1')),
        HinderedRotor(inertia=(0.00386014,'amu*angstrom^2'), symmetry=1, barrier=(2.17249,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0441,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.261,0.0385558,-3.66292e-05,1.77682e-08,-3.48358e-12,-11167.4,18.2126], Tmin=(100,'K'), Tmax=(1212.89,'K')), NASAPolynomial(coeffs=[9.59237,0.0143772,-6.72668e-06,1.33197e-09,-9.57128e-14,-12945.8,-18.572], Tmin=(1212.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-93.3689,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-O2sCsH)(F1s)(O2s-Cs)_ring) + radical(Csj(Cs-O2sCsH)(F1s)(H)_1935_ring)"""),
)

species(
    label = 'FC1C2OC2(F)C1(F)F(2765)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {7,S}
6  C u0 p0 c0 {5,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
8  C u0 p0 c0 {2,S} {6,S} {9,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-822.429,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.897799,0.0629918,-5.69925e-05,2.42559e-08,-4.05594e-12,-98799.5,18.048], Tmin=(100,'K'), Tmax=(1431.98,'K')), NASAPolynomial(coeffs=[17.4241,0.0168281,-8.63558e-06,1.74289e-09,-1.25513e-13,-103533,-67.6156], Tmin=(1431.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-822.429,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(CsCCFO) + group(CsCsCsFH) + group(CsCsCsFF) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + polycyclic(s2_3_4_ane)"""),
)

species(
    label = 'FC1=C(C(F)C(F)F)O1(2766)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {8,S} {9,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {2,S} {3,S} {6,S} {11,S}
8  C u0 p0 c0 {5,S} {6,S} {9,D}
9  C u0 p0 c0 {4,S} {5,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-580.563,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.482897,0.0940522,-0.000186347,1.88114e-07,-7.08991e-11,-69715.1,26.449], Tmin=(100,'K'), Tmax=(840.81,'K')), NASAPolynomial(coeffs=[3.79739,0.0402683,-2.25772e-05,4.48963e-09,-3.13003e-13,-68928.7,19.0242], Tmin=(840.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-580.563,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(CsCCFH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsOs) + group(CdCFO) + ring(oxirene)"""),
)

species(
    label = 'FC(F)=C(F)C1OC1F(2767)',
    structure = adjacencyList("""1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {7,S}
6  C u0 p0 c0 {5,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
8  C u0 p0 c0 {2,S} {6,S} {9,D}
9  C u0 p0 c0 {3,S} {4,S} {8,D}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-779.366,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02501,0.0612101,-5.88561e-05,2.74782e-08,-5.05471e-12,-93625.5,26.0447], Tmin=(100,'K'), Tmax=(1311.45,'K')), NASAPolynomial(coeffs=[15.6523,0.0165964,-7.82857e-06,1.53889e-09,-1.0997e-13,-97462.1,-48.4895], Tmin=(1311.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-779.366,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(CsCFHO) + group(CdCsCdF) + group(CdCFF) + longDistanceInteraction_noncyclic(Cds(F)2=Cds(F)) + ring(Cs-Cs(F)-O2s)"""),
)

species(
    label = 'O=C(F)[CH]C(F)[C](F)F(2768)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {10,S}
7  C u1 p0 c0 {6,S} {9,S} {11,S}
8  C u1 p0 c0 {2,S} {3,S} {6,S}
9  C u0 p0 c0 {4,S} {5,D} {7,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-720.326,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([259,529,569,1128,1321,1390,3140,3025,407.5,1350,352.5,190,488,555,1236,1407,611,648,830,1210,1753,235.941,238.501,2389.59],'cm^-1')),
        HinderedRotor(inertia=(0.944063,'amu*angstrom^2'), symmetry=1, barrier=(37.6654,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.266844,'amu*angstrom^2'), symmetry=1, barrier=(10.517,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.971262,'amu*angstrom^2'), symmetry=1, barrier=(37.6486,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.615297,0.0810901,-0.000123022,1.00906e-07,-3.34542e-11,-86519.6,28.4952], Tmin=(100,'K'), Tmax=(736.757,'K')), NASAPolynomial(coeffs=[10.4008,0.0279499,-1.4805e-05,2.9596e-09,-2.10498e-13,-87961.2,-15.7221], Tmin=(736.757,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-720.326,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFH) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(COCsFO) + radical(CCJC=O) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
)

species(
    label = 'F[CH][C](F)F(183)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {5,S}
3 F u0 p3 c0 {5,S}
4 C u1 p0 c0 {1,S} {5,S} {6,S}
5 C u1 p0 c0 {2,S} {3,S} {4,S}
6 H u0 p0 c0 {4,S}
"""),
    E0 = (-285.254,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([334,575,1197,1424,3202,190,488,555,1236,1407,1765.02],'cm^-1')),
        HinderedRotor(inertia=(0.510429,'amu*angstrom^2'), symmetry=1, barrier=(11.7358,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.0245,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.54626,0.0360262,-6.18601e-05,5.6897e-08,-2.03341e-11,-34259.7,17.1201], Tmin=(100,'K'), Tmax=(816.913,'K')), NASAPolynomial(coeffs=[5.79757,0.0130515,-6.72074e-06,1.3276e-09,-9.30775e-14,-34555.5,3.53266], Tmin=(816.913,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-285.254,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CsCsF1sH) + radical(CsCsF1sF1s)"""),
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
    label = 'F[C](F)C(F)C1=C(F)O1(2769)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {9,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {7,S} {8,S}
6  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
7  C u0 p0 c0 {5,S} {6,S} {8,D}
8  C u0 p0 c0 {4,S} {5,S} {7,D}
9  C u1 p0 c0 {2,S} {3,S} {6,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-377.257,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([174,267,591,721,1107,1278,1348,3273,190,488,555,1236,1407,180,180,180,275.995,957.347,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.149957,'amu*angstrom^2'), symmetry=1, barrier=(3.4478,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149957,'amu*angstrom^2'), symmetry=1, barrier=(3.4478,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.044,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.634231,0.0931848,-0.000196744,2.04179e-07,-7.78118e-11,-45270.9,27.714], Tmin=(100,'K'), Tmax=(847.14,'K')), NASAPolynomial(coeffs=[2.20986,0.0404968,-2.3332e-05,4.65916e-09,-3.24522e-13,-43914.2,29.9567], Tmin=(847.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-377.257,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(CsCCFH) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cds-CdsCsOs) + group(CdCFO) + ring(oxirene) + radical(Csj(Cs-F1sCdH)(F1s)(F1s))"""),
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
    label = 'F[C]1OC1C=C(F)F(2722)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {7,S}
5  C u0 p0 c0 {4,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {5,S} {8,D} {10,S}
7  C u1 p0 c0 {1,S} {4,S} {5,S}
8  C u0 p0 c0 {2,S} {3,S} {6,D}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-415.688,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,180,180,180,180,941.008,941.016,941.019,941.021,941.021,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00434442,'amu*angstrom^2'), symmetry=1, barrier=(2.72988,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63305,0.0497536,-4.57309e-05,2.08991e-08,-3.81786e-12,-49908.7,25.64], Tmin=(100,'K'), Tmax=(1307,'K')), NASAPolynomial(coeffs=[12.4648,0.0166031,-7.68489e-06,1.49259e-09,-1.05793e-13,-52740.1,-29.5171], Tmin=(1307,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-415.688,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(CsCFHO) + group(Cds-CdsCsH) + group(CdCFF) + ring(Cs-Cs(F)-O2s) + radical(CsCsF1sO2s)"""),
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
    label = 'F[C]1OC1C(F)=C(F)F(2770)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {8,S}
6  C u0 p0 c0 {5,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {1,S} {6,S} {9,D}
8  C u1 p0 c0 {2,S} {5,S} {6,S}
9  C u0 p0 c0 {3,S} {4,S} {7,D}
10 H u0 p0 c0 {6,S}
"""),
    E0 = (-584.812,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,323,467,575,827,1418,182,240,577,636,1210,1413,180,180,180,180,845.319,845.321,845.321,845.322,2157.22,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00981905,'amu*angstrom^2'), symmetry=1, barrier=(4.97894,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.044,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20567,0.0658717,-8.51617e-05,5.69025e-08,-1.52731e-11,-70239.7,26.9041], Tmin=(100,'K'), Tmax=(904.794,'K')), NASAPolynomial(coeffs=[11.3344,0.0210928,-1.09245e-05,2.20256e-09,-1.58891e-13,-72072.6,-20.9478], Tmin=(904.794,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-584.812,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(CsCFHO) + group(CdCsCdF) + group(CdCFF) + longDistanceInteraction_noncyclic(Cds(F)2=Cds(F)) + ring(Cs-Cs(F)-O2s) + radical(CsCsF1sO2s)"""),
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
    label = 'F[C]1O[C]1C=C(F)F(2771)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {7,S}
2 F u0 p3 c0 {8,S}
3 F u0 p3 c0 {8,S}
4 O u0 p2 c0 {5,S} {7,S}
5 C u1 p0 c0 {4,S} {6,S} {7,S}
6 C u0 p0 c0 {5,S} {8,D} {9,S}
7 C u1 p0 c0 {1,S} {4,S} {5,S}
8 C u0 p0 c0 {2,S} {3,S} {6,D}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-235.583,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,182,240,577,636,1210,1413,329.526,329.542,329.568,329.615,329.617,858.929,1665.14,1665.18,3853.78],'cm^-1')),
        HinderedRotor(inertia=(0.00155278,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71145,0.0547939,-7.64759e-05,5.9145e-08,-1.87123e-11,-28255.8,26.4804], Tmin=(100,'K'), Tmax=(767.687,'K')), NASAPolynomial(coeffs=[8.1557,0.0212143,-1.08602e-05,2.16021e-09,-1.53882e-13,-29245.2,-2.9053], Tmin=(767.687,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-235.583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(CsCFHO) + group(Cds-CdsCsH) + group(CdCFF) + ring(Cs-Cs(F)-O2s) + radical(C2CsJO) + radical(CsCsF1sO2s)"""),
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
    label = 'F[C](F)C(F)[C]1OC1F(2772)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {8,S}
6  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
7  C u0 p0 c0 {2,S} {5,S} {8,S} {11,S}
8  C u1 p0 c0 {5,S} {6,S} {7,S}
9  C u1 p0 c0 {3,S} {4,S} {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-560.533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.871145,0.0770896,-0.000125955,1.14962e-07,-4.15156e-11,-67311.9,28.2125], Tmin=(100,'K'), Tmax=(804.957,'K')), NASAPolynomial(coeffs=[7.1069,0.032317,-1.68338e-05,3.31151e-09,-2.31949e-13,-67869.2,2.25577], Tmin=(804.957,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-560.533,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCFHO) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + ring(O2s-Cs(F)-Cs(C)) + radical(C2CsJO) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
)

species(
    label = 'F[C]1OC1[C](F)C(F)F(2773)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {9,S}
6  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {8,S} {11,S}
8  C u1 p0 c0 {3,S} {6,S} {7,S}
9  C u1 p0 c0 {4,S} {5,S} {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-528.444,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06635,0.0704837,-9.42022e-05,7.01386e-08,-2.15422e-11,-63456.8,27.9091], Tmin=(100,'K'), Tmax=(787.157,'K')), NASAPolynomial(coeffs=[9.25487,0.028873,-1.49089e-05,2.98251e-09,-2.13513e-13,-64746,-9.63601], Tmin=(787.157,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-528.444,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCFHO) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + ring(O2s-Cs(F)-Cs(C)) + radical(CsCsCsF1s) + radical(Csj(Cs-O2sCsH)(F1s)(O2s-Cs)_ring)"""),
)

species(
    label = 'F[C](F)[C](F)C1OC1F(2774)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {7,S}
6  C u0 p0 c0 {5,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
8  C u1 p0 c0 {2,S} {6,S} {9,S}
9  C u1 p0 c0 {3,S} {4,S} {8,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-550.164,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.86962,0.0590089,-2.84758e-05,-1.01878e-07,1.29672e-10,-66105.6,25.7752], Tmin=(100,'K'), Tmax=(439.691,'K')), NASAPolynomial(coeffs=[7.50012,0.0316341,-1.64421e-05,3.22819e-09,-2.25567e-13,-66831.3,0.616402], Tmin=(439.691,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-550.164,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCFHO) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + ring(O2s-Cs(F)-Cs(C)) + radical(CsCsCsF1s) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
)

species(
    label = 'F[C]1O[C]1C(F)C(F)F(2775)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {8,S} {9,S}
6  C u0 p0 c0 {1,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {2,S} {3,S} {6,S} {11,S}
8  C u1 p0 c0 {5,S} {6,S} {9,S}
9  C u1 p0 c0 {4,S} {5,S} {8,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-538.813,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.916128,0.0738467,-0.000105733,8.43725e-08,-2.75484e-11,-64698.8,27.5449], Tmin=(100,'K'), Tmax=(744.783,'K')), NASAPolynomial(coeffs=[9.32533,0.0286791,-1.4756e-05,2.93007e-09,-2.08178e-13,-65951.3,-10.5459], Tmin=(744.783,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-538.813,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCFHO) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + ring(O2s-Cs(F)-Cs(C)) + radical(C2CsJO) + radical(Csj(Cs-O2sCsH)(F1s)(O2s-Cs)_ring)"""),
)

species(
    label = 'F[C](F)[CH]C1OC1(F)F(2776)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {7,S}
6  C u0 p0 c0 {5,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
8  C u1 p0 c0 {6,S} {9,S} {11,S}
9  C u1 p0 c0 {3,S} {4,S} {8,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-566.078,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35679,0.0599135,-6.14121e-05,3.2801e-08,-7.14268e-12,-67989.8,28.6684], Tmin=(100,'K'), Tmax=(1095.73,'K')), NASAPolynomial(coeffs=[11.429,0.0231444,-1.1077e-05,2.17595e-09,-1.55308e-13,-70197.1,-20.8452], Tmin=(1095.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-566.078,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCFFO) + group(CsCsFFH) + ring(Cs(C)-Cs(F)(F)-O2s) + radical(CCJCO) + radical(Csj(Cs-CsHH)(F1s)(F1s))"""),
)

species(
    label = 'F[C]1OC1[CH]C(F)(F)F(2777)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {9,S}
6  C u0 p0 c0 {5,S} {8,S} {9,S} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {3,S} {8,S}
8  C u1 p0 c0 {6,S} {7,S} {11,S}
9  C u1 p0 c0 {4,S} {5,S} {6,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-573.579,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2425,0.0579282,-5.32815e-05,2.42994e-08,-4.42955e-12,-68884.1,27.9948], Tmin=(100,'K'), Tmax=(1309.48,'K')), NASAPolynomial(coeffs=[13.9026,0.0192558,-8.98226e-06,1.74618e-09,-1.2376e-13,-72199.8,-36.4964], Tmin=(1309.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-573.579,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(CsCFHO) + group(CsCsFFF) + ring(O2s-Cs(F)-Cs(C)) + radical(CCJCO) + radical(Csj(Cs-O2sCsH)(F1s)(O2s-Cs)_ring)"""),
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
    E0 = (-113.089,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (105.368,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-33.6966,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-33.6966,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (264.269,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-185.347,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-130.232,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-130.232,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-193.632,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (113.112,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (158.458,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (48.0367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (-57.1832,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-29.4446,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (136.413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (69.5993,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (26.8296,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-19.4552,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-60.1321,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (-52.6069,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (-71.2748,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (9.81244,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (-32.996,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['F[C](F)C(F)C1O[C]1F(2762)'],
    products = ['CHFCF2(54)', 'FC1=CO1(390)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(80.5425,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission
Ea raised from 0.0 to 80.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction2',
    reactants = ['F[C](F)C(F)C1O[C]1F(2762)'],
    products = ['F[C]1OC([C](F)F)C1F(2763)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction3',
    reactants = ['F[C](F)C(F)C1O[C]1F(2762)'],
    products = ['F[C](F)C(F)C1(F)[CH]O1(2761)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['F[C](F)C(F)C1O[C]1F(2762)'],
    products = ['F[CH]C(F)(F)C1O[C]1F(2760)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CF2(42)', 'F[CH]C1O[C]1F(2764)'],
    products = ['F[C](F)C(F)C1O[C]1F(2762)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_sec_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['F[C](F)C(F)C1O[C]1F(2762)'],
    products = ['FC1C2OC2(F)C1(F)F(2765)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_noH;Cpri_rad_out_noH]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F[C](F)C(F)C1O[C]1F(2762)'],
    products = ['FC1=C(C(F)C(F)F)O1(2766)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['F[C](F)C(F)C1O[C]1F(2762)'],
    products = ['FC(F)=C(F)C1OC1F(2767)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['O=C(F)[CH]C(F)[C](F)F(2768)'],
    products = ['F[C](F)C(F)C1O[C]1F(2762)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5.92717e+11,'s^-1'), n=0.412677, Ea=(202.784,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 201.5 to 202.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction10',
    reactants = ['F[CH][C](F)F(183)', 'FC1=CO1(390)'],
    products = ['F[C](F)C(F)C1O[C]1F(2762)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(5)', 'F[C](F)C(F)C1=C(F)O1(2769)'],
    products = ['F[C](F)C(F)C1O[C]1F(2762)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.97324,'m^3/(mol*s)'), n=2.12322, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0984061311673169, var=1.772613507237397, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R"""),
)

reaction(
    label = 'reaction12',
    reactants = ['F(37)', 'F[C]1OC1C=C(F)F(2722)'],
    products = ['F[C](F)C(F)C1O[C]1F(2762)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(66.9222,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CHFCF2(54)', 'F[C]1[CH]O1(2636)'],
    products = ['F[C](F)C(F)C1O[C]1F(2762)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.001023,'m^3/(mol*s)'), n=2.607, Ea=(32.6056,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_4R!H->O',), comment="""Estimated from node Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_4R!H->O"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(5)', 'F[C]1OC1C(F)=C(F)F(2770)'],
    products = ['F[C](F)C(F)C1O[C]1F(2762)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(50675.7,'m^3/(mol*s)'), n=0.856706, Ea=(19.6514,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.3149223058909152, var=1.98898808922157, Tref=1000.0, N=74, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_N-2CS-inRing_Ext-1COS-R_N-4R!H-inRing_4R!H-u0_Ext-2CS-R_Sp-2CS=1CCOSS_N-5R!H-inRing_Sp-5R!H-2CS_Ext-2CS-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['F[CH][C](F)F(183)', 'F[C]1[CH]O1(2636)'],
    products = ['F[C](F)C(F)C1O[C]1F(2762)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(5e+07,'m^3/(mol*s)'), n=2.59562e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing"""),
)

reaction(
    label = 'reaction16',
    reactants = ['HF(38)', 'F[C]1O[C]1C=C(F)F(2771)'],
    products = ['F[C](F)C(F)C1O[C]1F(2762)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(281.116,'m^3/(mol*s)'), n=1.03051, Ea=(262.385,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17160344105964337, var=17.74244203879562, Tref=1000.0, N=7, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction17',
    reactants = ['CF2(87)', 'F[CH]C1O[C]1F(2764)'],
    products = ['F[C](F)C(F)C1O[C]1F(2762)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(786.723,'m^3/(mol*s)'), n=1.25031, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl"""),
)

reaction(
    label = 'reaction18',
    reactants = ['F[C](F)C(F)C1O[C]1F(2762)'],
    products = ['F[C](F)C(F)[C]1OC1F(2772)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(9.06381e+10,'s^-1'), n=0.647667, Ea=(174.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_noH;Cs_H_out_NonDe] for rate rule [R2H_S_cy3;C_rad_out_noH;Cs_H_out_NDMustO]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['F[C](F)C(F)C1O[C]1F(2762)'],
    products = ['F[C]1OC1[C](F)C(F)F(2773)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.70914e+10,'s^-1'), n=0.735063, Ea=(133.5,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_noH;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['F[C](F)C(F)C1O[C]1F(2762)'],
    products = ['F[C](F)[C](F)C1OC1F(2774)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.95012e+07,'s^-1'), n=1.37671, Ea=(141.025,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_noH;XH_out] for rate rule [R3H_SS_12cy3;C_rad_out_noH;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['F[C](F)C(F)C1O[C]1F(2762)'],
    products = ['F[C]1O[C]1C(F)C(F)F(2775)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.50974e+07,'s^-1'), n=1.33047, Ea=(122.357,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_noH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['F[C](F)C(F)C1O[C]1F(2762)'],
    products = ['F[C](F)[CH]C1OC1(F)F(2776)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.00930803,'s^-1'), n=4.16824, Ea=(203.444,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F"""),
)

reaction(
    label = 'reaction23',
    reactants = ['F[C](F)C(F)C1O[C]1F(2762)'],
    products = ['F[C]1OC1[CH]C(F)(F)F(2777)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(5.38157e+14,'s^-1'), n=-0.447076, Ea=(160.636,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08474952276158404, var=26.08378321593064, Tref=1000.0, N=4, data_mean=0.0, correlation='R2F_Ext-1R!H-R',), comment="""Estimated from node R2F_Ext-1R!H-R"""),
)

network(
    label = 'PDepNetwork #1057',
    isomers = [
        'F[C](F)C(F)C1O[C]1F(2762)',
    ],
    reactants = [
        ('CHFCF2(54)', 'FC1=CO1(390)'),
    ],
    bathGas = {
        'N2': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1057',
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

