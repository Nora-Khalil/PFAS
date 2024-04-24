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
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,215,315,519,588,595,1205,1248,334,575,1197,1424,3202,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.69488,0.0791066,-0.000120582,1.01033e-07,-3.43069e-11,-63647.2,27.3241], Tmin=(100,'K'), Tmax=(738.886,'K')), NASAPolynomial(coeffs=[9.7708,0.0285018,-1.48626e-05,2.95019e-09,-2.08838e-13,-64948.2,-13.4436], Tmin=(738.886,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-530.132,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCFHO) + group(CsCsFHH) + ring(O2s-Cs(F)-Cs(C)) + radical(Csj(Cs-O2sCsH)(F1s)(O2s-Cs)_ring) + radical(Csj(Cs-CsF1sF1s)(F1s)(H))"""),
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
    label = 'F[CH]C1(F)O[C](F)C1F(2893)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {8,S}
6  C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
7  C u0 p0 c0 {2,S} {6,S} {8,S} {10,S}
8  C u1 p0 c0 {3,S} {5,S} {7,S}
9  C u1 p0 c0 {4,S} {6,S} {11,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-558.49,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.852431,0.0756252,-0.000103945,7.73261e-08,-2.34691e-11,-67063.3,24.4578], Tmin=(100,'K'), Tmax=(799.155,'K')), NASAPolynomial(coeffs=[10.3371,0.0281517,-1.48385e-05,2.99168e-09,-2.15008e-13,-68579.2,-19.1737], Tmin=(799.155,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-558.49,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCFO) + group(CsCsCsFH) + group(CsCFHO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + ring(O2s-Cs-Cs-Cs(F)) + radical(Csj(Cs-F1sCsH)(F1s)(O2s-Cs)_ring) + radical(Csj(Cs-F1sO2sCs)(F1s)(H)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F))"""),
)

species(
    label = 'F[CH]C(F)(F)C1(F)[CH]O1(2759)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {8,S}
6  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
7  C u0 p0 c0 {2,S} {3,S} {6,S} {9,S}
8  C u1 p0 c0 {5,S} {6,S} {10,S}
9  C u1 p0 c0 {4,S} {7,S} {11,S}
10 H u0 p0 c0 {8,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-525.208,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,2950,1000,334,575,1197,1424,3202,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3279.08,'J/mol'), sigma=(5.91083,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=512.18 K, Pc=36.03 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.796721,0.07799,-0.000110064,7.75942e-08,-1.90732e-11,-63059.6,26.6411], Tmin=(100,'K'), Tmax=(626.751,'K')), NASAPolynomial(coeffs=[10.492,0.0269341,-1.37689e-05,2.71143e-09,-1.91217e-13,-64487.4,-17.2986], Tmin=(626.751,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-525.208,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCFO) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-CsOsHH) + group(CsCsFHH) + ring(Cs(F)(C)-Cs-O2s) + radical(Csj(Cs-F1sO2sCs)(O2s-Cs)(H)_ring) + radical(Csj(Cs-CsF1sF1s)(F1s)(H))"""),
)

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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3221.19,'J/mol'), sigma=(5.26728,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=503.14 K, Pc=50.02 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07832,0.070088,-9.12071e-05,6.53218e-08,-1.92646e-11,-62146,27.1611], Tmin=(100,'K'), Tmax=(818.676,'K')), NASAPolynomial(coeffs=[9.63898,0.0282623,-1.4575e-05,2.92026e-09,-2.09493e-13,-63547.7,-12.4266], Tmin=(818.676,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-517.543,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCFHO) + group(CsCsFFH) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + ring(O2s-Cs(F)-Cs(C)) + radical(Csj(Cs-O2sCsH)(F1s)(O2s-Cs)_ring) + radical(Csj(Cs-CsF1sH)(F1s)(F1s))"""),
)

species(
    label = 'CHF(40)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {2,S}
2 C u2 p0 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (214.928,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([787.277,1376.72,4000],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (32.017,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2402.79,'J/mol'), sigma=(4.33176,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=375.31 K, Pc=67.08 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.93333,-0.00026336,8.89187e-06,-1.03032e-08,3.50809e-12,25853.7,4.33729], Tmin=(100,'K'), Tmax=(1056.12,'K')), NASAPolynomial(coeffs=[4.72426,0.00164131,-7.73114e-07,1.90987e-10,-1.59925e-14,25413.4,-0.815528], Tmin=(1056.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.928,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CH2_triplet)"""),
)

species(
    label = 'F[C](F)C1O[C]1F(2645)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {7,S}
4 O u0 p2 c0 {5,S} {6,S}
5 C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
6 C u1 p0 c0 {1,S} {4,S} {5,S}
7 C u1 p0 c0 {2,S} {3,S} {5,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-298.903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,190,488,555,1236,1407,180,180,324.861,971.33,971.551,971.584,971.63,971.812,971.919,2089.26],'cm^-1')),
        HinderedRotor(inertia=(1.6237,'amu*angstrom^2'), symmetry=1, barrier=(37.332,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.035,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.82538,0.0424736,-3.98931e-05,1.79937e-08,-3.15273e-12,-35867,21.5195], Tmin=(100,'K'), Tmax=(1389.13,'K')), NASAPolynomial(coeffs=[13.3538,0.00927683,-4.04634e-06,7.90013e-10,-5.65461e-14,-39069.9,-37.8877], Tmin=(1389.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-298.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-O2sCsH)(F1s)(O2s-Cs)_ring) + radical(Csj(Cs-O2sCsH)(F1s)(F1s)_1952_ring)"""),
)

species(
    label = 'FC1C(F)(F)C2OC12F(2784)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {8,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {7,S}
6  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
7  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
8  C u0 p0 c0 {2,S} {6,S} {9,S} {11,S}
9  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-820.986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.95001,0.0625318,-5.60538e-05,2.35902e-08,-3.90908e-12,-98628.5,18.1528], Tmin=(100,'K'), Tmax=(1439.98,'K')), NASAPolynomial(coeffs=[17.1882,0.0174253,-9.0672e-06,1.8369e-09,-1.3242e-13,-103305,-66.1079], Tmin=(1439.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-820.986,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsCCFO) + group(Cs-CsCsOsH) + group(CsCsCsFH) + group(CsCsCsFF) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)-Cs(F)) + longDistanceInteraction_cyclic(Cs(F)2-Cs(F)) + polycyclic(s2_3_4_ane)"""),
)

species(
    label = 'FCC(F)(F)C1=C(F)O1(2894)',
    structure = adjacencyList("""1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {8,S} {9,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {3,S} {6,S} {10,S} {11,S}
8  C u0 p0 c0 {5,S} {6,S} {9,D}
9  C u0 p0 c0 {4,S} {5,S} {8,D}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-569.243,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.860696,0.0788276,-0.000135489,1.28703e-07,-4.78692e-11,-68360.6,25.0312], Tmin=(100,'K'), Tmax=(803.249,'K')), NASAPolynomial(coeffs=[6.18308,0.0342762,-1.85912e-05,3.7109e-09,-2.61723e-13,-68633.5,4.14396], Tmin=(803.249,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-569.243,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + group(Cds-CdsCsOs) + group(CdCFO) + ring(Cd(C-FF)-Cd-O2s)"""),
)

species(
    label = 'O=C(F)[CH]C(F)(F)[CH]F(2895)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {9,D}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u1 p0 c0 {6,S} {9,S} {10,S}
8  C u1 p0 c0 {3,S} {6,S} {11,S}
9  C u0 p0 c0 {4,S} {5,D} {7,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {8,S}
"""),
    E0 = (-732.915,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([215,315,519,588,595,1205,1248,3025,407.5,1350,352.5,334,575,1197,1424,3202,611,648,830,1210,1753,180,180,2549.34],'cm^-1')),
        HinderedRotor(inertia=(0.257082,'amu*angstrom^2'), symmetry=1, barrier=(5.91081,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.6024,'amu*angstrom^2'), symmetry=1, barrier=(36.8423,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.60752,'amu*angstrom^2'), symmetry=1, barrier=(36.9601,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.321971,0.0889305,-0.00014764,1.29281e-07,-4.46975e-11,-88024.7,28.3419], Tmin=(100,'K'), Tmax=(794.764,'K')), NASAPolynomial(coeffs=[10.466,0.0283071,-1.51623e-05,3.00626e-09,-2.11247e-13,-89334.9,-16.3657], Tmin=(794.764,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-732.915,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFHH) + group(COCsFO) + radical(CCJC=O) + radical(Csj(Cs-CsF1sF1s)(F1s)(H))"""),
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
    label = 'F[CH]C(F)(F)C1=C(F)O1(2896)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {9,S}
4  F u0 p3 c0 {8,S}
5  O u0 p2 c0 {7,S} {8,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {9,S}
7  C u0 p0 c0 {5,S} {6,S} {8,D}
8  C u0 p0 c0 {4,S} {5,S} {7,D}
9  C u1 p0 c0 {3,S} {6,S} {10,S}
10 H u0 p0 c0 {9,S}
"""),
    E0 = (-368.137,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([248,333,466,604,684,796,1061,1199,334,575,1197,1424,3202,180,180,180,180,1183.74,1600,1740.63,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.136228,'amu*angstrom^2'), symmetry=1, barrier=(3.13216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136228,'amu*angstrom^2'), symmetry=1, barrier=(3.13216,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (141.044,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.969766,0.0769377,-0.000138375,1.32827e-07,-4.9213e-11,-44177.5,27.672], Tmin=(100,'K'), Tmax=(815.889,'K')), NASAPolynomial(coeffs=[6.13303,0.0313563,-1.73115e-05,3.45846e-09,-2.43254e-13,-44345.5,7.94689], Tmin=(815.889,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-368.137,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCsFHH) + group(Cds-CdsCsOs) + group(CdCFO) + ring(Cd(C-FF)-Cd-O2s) + radical(Csj(Cs-F1sF1sCd)(F1s)(H))"""),
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
    label = 'FC=C(F)C1O[C]1F(2897)',
    structure = adjacencyList("""multiplicity 2
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  O u0 p2 c0 {5,S} {7,S}
5  C u0 p0 c0 {4,S} {6,S} {7,S} {9,S}
6  C u0 p0 c0 {1,S} {5,S} {8,D}
7  C u1 p0 c0 {2,S} {4,S} {5,S}
8  C u0 p0 c0 {3,S} {6,D} {10,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {8,S}
"""),
    E0 = (-398.235,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,323,467,575,827,1418,194,682,905,1196,1383,3221,180,180,180,840.733,840.775,840.933,841.219,841.282,1978.18,4000],'cm^-1')),
        HinderedRotor(inertia=(0.233788,'amu*angstrom^2'), symmetry=1, barrier=(5.37524,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (123.053,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5969,0.0548078,-5.90083e-05,3.26908e-08,-7.32718e-12,-47811.7,24.7823], Tmin=(100,'K'), Tmax=(1070.38,'K')), NASAPolynomial(coeffs=[11.0343,0.0195396,-9.58362e-06,1.90697e-09,-1.37138e-13,-49832,-21.3896], Tmin=(1070.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-398.235,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(CsCFHO) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + ring(Cs-Cs(F)-O2s) + radical(CsCsF1sO2s)"""),
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
    label = 'FC=C(F)[C]1O[C]1F(2898)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {6,S}
2 F u0 p3 c0 {7,S}
3 F u0 p3 c0 {8,S}
4 O u0 p2 c0 {5,S} {7,S}
5 C u1 p0 c0 {4,S} {6,S} {7,S}
6 C u0 p0 c0 {1,S} {5,S} {8,D}
7 C u1 p0 c0 {2,S} {4,S} {5,S}
8 C u0 p0 c0 {3,S} {6,D} {9,S}
9 H u0 p0 c0 {8,S}
"""),
    E0 = (-218.13,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([271,519,563,612,1379,194,682,905,1196,1383,3221,252.785,252.8,252.803,252.804,252.811,1296.1,1815.68,1815.7,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00263814,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (122.045,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3646,0.063734,-0.00010436,9.13202e-08,-3.158e-11,-26145.7,26.722], Tmin=(100,'K'), Tmax=(800.345,'K')), NASAPolynomial(coeffs=[8.37325,0.0213517,-1.11445e-05,2.19304e-09,-1.53559e-13,-27032,-4.05831], Tmin=(800.345,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-218.13,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(CsCFHO) + group(CdCsCdF) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + group(CdCFH) + longDistanceInteraction_noncyclic(Cds(F)=Cds(F)) + ring(Cs-Cs(F)-O2s) + radical(C2CsJO) + radical(CsCsF1sO2s)"""),
)

species(
    label = 'CHF(86)',
    structure = adjacencyList("""1 F u0 p3 c0 {2,S}
2 C u0 p1 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (138.756,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1169.21,1416.41,2978.06],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (32.017,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2178.39,'J/mol'), sigma=(4.123,'angstroms'), dipoleMoment=(1.8,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NIST_Fluorine"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.34484,0.00235461,1.93983e-06,-2.65251e-09,7.91169e-13,16766.1,7.05286], Tmin=(298,'K'), Tmax=(1300,'K')), NASAPolynomial(coeffs=[4.48366,0.00174964,-5.0479e-07,1.08953e-10,-9.87898e-15,16210.2,0.289222], Tmin=(1300,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(138.756,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(58.2013,'J/mol/K'), label="""CHF""", comment="""Thermo library: halogens"""),
)

species(
    label = 'F[CH]C(F)(F)[C]1OC1F(2881)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {7,S} {8,S}
6  C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
7  C u0 p0 c0 {3,S} {5,S} {8,S} {10,S}
8  C u1 p0 c0 {5,S} {6,S} {7,S}
9  C u1 p0 c0 {4,S} {6,S} {11,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-573.122,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.659224,0.0838439,-0.000146119,1.36485e-07,-4.92843e-11,-68820.4,27.7749], Tmin=(100,'K'), Tmax=(828.987,'K')), NASAPolynomial(coeffs=[6.98065,0.0330256,-1.74049e-05,3.41066e-09,-2.3718e-13,-69170.4,2.67375], Tmin=(828.987,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-573.122,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCFHO) + group(CsCsFHH) + ring(O2s-Cs(F)-Cs(C)) + radical(C2CsJO) + radical(Csj(Cs-CsF1sF1s)(F1s)(H))"""),
)

species(
    label = 'FCC(F)(F)[C]1O[C]1F(2899)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {6,S}
2  F u0 p3 c0 {6,S}
3  F u0 p3 c0 {7,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {8,S} {9,S}
6  C u0 p0 c0 {1,S} {2,S} {7,S} {8,S}
7  C u0 p0 c0 {3,S} {6,S} {10,S} {11,S}
8  C u1 p0 c0 {5,S} {6,S} {9,S}
9  C u1 p0 c0 {4,S} {5,S} {8,S}
10 H u0 p0 c0 {7,S}
11 H u0 p0 c0 {7,S}
"""),
    E0 = (-549.268,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.708593,0.0790471,-0.000122642,1.05059e-07,-3.63405e-11,-65949.5,27.2514], Tmin=(100,'K'), Tmax=(760.759,'K')), NASAPolynomial(coeffs=[9.33603,0.0291185,-1.51934e-05,3.01019e-09,-2.12569e-13,-67130.1,-11.1434], Tmin=(760.759,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-549.268,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(CsCsCsFF) + longDistanceInteraction_noncyclic(Cs(F)2-Cs(F)) + group(CsCFHO) + group(CsCsFHH) + ring(O2s-Cs(F)-Cs(C)) + radical(C2CsJO) + radical(Csj(Cs-O2sCsH)(F1s)(O2s-Cs)_ring)"""),
)

species(
    label = 'F[CH][C](F)C1OC1(F)F(2900)',
    structure = adjacencyList("""multiplicity 3
1  F u0 p3 c0 {7,S}
2  F u0 p3 c0 {7,S}
3  F u0 p3 c0 {8,S}
4  F u0 p3 c0 {9,S}
5  O u0 p2 c0 {6,S} {7,S}
6  C u0 p0 c0 {5,S} {7,S} {8,S} {10,S}
7  C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
8  C u1 p0 c0 {3,S} {6,S} {9,S}
9  C u1 p0 c0 {4,S} {8,S} {11,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {9,S}
"""),
    E0 = (-559.897,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (142.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.952145,0.0741149,-0.000102632,7.30024e-08,-1.87341e-11,-67236.9,28.304], Tmin=(100,'K'), Tmax=(629.923,'K')), NASAPolynomial(coeffs=[9.75032,0.02758,-1.40456e-05,2.77052e-09,-1.95888e-13,-68530.5,-11.5458], Tmin=(629.923,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-559.897,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(CsCsCsFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(CsCFFO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + ring(Cs(C)-Cs(F)(F)-O2s) + radical(CsCsCsF1s) + radical(Csj(Cs-F1sCsH)(F1s)(H))"""),
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
    E0 = (-133.553,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (72.3156,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (-64.4422,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (-54.1603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (219.472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-218.4,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (-163.284,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-226.684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (92.6483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (147.115,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (35.8983,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-81.02,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (115.949,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (15.6411,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (143.3,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (-52.5079,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (-70.8914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-13.8289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-17.1343,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['F[CH]C(F)(F)C1O[C]1F(2760)'],
    products = ['CHFCF2(54)', 'FC1=CO1(390)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(93.1315,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission
Ea raised from 0.0 to 93.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction2',
    reactants = ['F[CH]C(F)(F)C1O[C]1F(2760)'],
    products = ['F[CH]C1(F)O[C](F)C1F(2893)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2e+13,'s^-1'), n=0, Ea=(299,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [OF]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_XY_interchange"""),
)

reaction(
    label = 'reaction3',
    reactants = ['F[CH]C(F)(F)C1(F)[CH]O1(2759)'],
    products = ['F[CH]C(F)(F)C1O[C]1F(2760)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ;C]
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
    reactants = ['CHF(40)', 'F[C](F)C1O[C]1F(2645)'],
    products = ['F[CH]C(F)(F)C1O[C]1F(2760)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['F[CH]C(F)(F)C1O[C]1F(2760)'],
    products = ['FC1C(F)(F)C2OC12F(2784)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_single] for rate rule [R4_SSS;C_rad_out_noH;Cpri_rad_out_1H]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['F[CH]C(F)(F)C1O[C]1F(2760)'],
    products = ['FCC(F)(F)C1=C(F)O1(2894)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['O=C(F)[CH]C(F)(F)[CH]F(2895)'],
    products = ['F[CH]C(F)(F)C1O[C]1F(2760)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(5.92717e+11,'s^-1'), n=0.412677, Ea=(202.784,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 201.5 to 202.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['F[CH][C](F)F(183)', 'FC1=CO1(390)'],
    products = ['F[CH]C(F)(F)C1O[C]1F(2760)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.04353e-06,'m^3/(mol*s)'), n=3.0961, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.25403387200380967, var=0.12219132316784118, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Ext-4R!H-R_Ext-3R-R_N-Sp-5R!H=4R!H"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(5)', 'F[CH]C(F)(F)C1=C(F)O1(2896)'],
    products = ['F[CH]C(F)(F)C1O[C]1F(2760)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.97324,'m^3/(mol*s)'), n=2.12322, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0984061311673169, var=1.772613507237397, Tref=1000.0, N=8, data_mean=0.0, correlation='Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R',), comment="""Estimated from node Root_3R->H_N-2R!H->N_N-1R!H->N_N-2COS->O_2CS-inRing_N-Sp-2CS-=1CCOSS_1COS-inRing_Ext-2CS-R_Ext-4R!H-R_Ext-1COS-R"""),
)

reaction(
    label = 'reaction11',
    reactants = ['F(37)', 'FC=C(F)C1O[C]1F(2897)'],
    products = ['F[CH]C(F)(F)C1O[C]1F(2760)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(57.7946,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CHFCF2(54)', 'F[C]1[CH]O1(2636)'],
    products = ['F[CH]C(F)(F)C1O[C]1F(2760)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(0.001023,'m^3/(mol*s)'), n=2.607, Ea=(29.2324,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_4R!H->O',), comment="""Estimated from node Root_3R-inRing_Ext-3R-R_Sp-4R!H-3R_4R!H->O"""),
)

reaction(
    label = 'reaction13',
    reactants = ['F[CH][C](F)F(183)', 'F[C]1[CH]O1(2636)'],
    products = ['F[CH]C(F)(F)C1O[C]1F(2760)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(5e+07,'m^3/(mol*s)'), n=2.59562e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_N-1BrCFOS->O_1BrCFS-inRing_Ext-1BrCFS-R_Sp-3R!H-1BrCFS_2R-inRing"""),
)

reaction(
    label = 'reaction14',
    reactants = ['HF(38)', 'FC=C(F)[C]1O[C]1F(2898)'],
    products = ['F[CH]C(F)(F)C1O[C]1F(2760)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.14111,'m^3/(mol*s)'), n=1.29695, Ea=(211.438,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R',), comment="""Estimated from node HF_Ext-3COCdCddCtO2d-R_Ext-4COCdCddCtO2d-R_6R!H->F_Ext-3COCdCddCtO2d-R"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CHF(86)', 'F[C](F)C1O[C]1F(2645)'],
    products = ['F[CH]C(F)(F)C1O[C]1F(2760)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(786.723,'m^3/(mol*s)'), n=1.25031, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl',), comment="""Estimated from node Root_N-3R->H_Ext-3BrCClFINOPSSi-R_3BrCClFINOPSSi->C_2Br1sCl1sF1s->F1s_N-4R!H->Br_N-4ClF->Cl"""),
)

reaction(
    label = 'reaction16',
    reactants = ['F[CH]C(F)(F)C1O[C]1F(2760)'],
    products = ['F[CH]C(F)(F)[C]1OC1F(2881)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(9.06381e+10,'s^-1'), n=0.647667, Ea=(174.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_noH;Cs_H_out_NonDe] for rate rule [R2H_S_cy3;C_rad_out_noH;Cs_H_out_NDMustO]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['F[CH]C(F)(F)C1O[C]1F(2760)'],
    products = ['FCC(F)(F)[C]1O[C]1F(2899)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(91.367,'s^-1'), n=3.04268, Ea=(155.793,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_1H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['F[CH]C(F)(F)C1O[C]1F(2760)'],
    products = ['F[CH][C](F)C1OC1(F)F(2900)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.0186161,'s^-1'), n=4.16824, Ea=(212.856,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R3F',), comment="""Estimated from node R3F
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction19',
    reactants = ['F[CH]C(F)(F)C1O[C]1F(2760)'],
    products = ['F[C]1OC1[C](F)C(F)F(2773)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(8.5608e+11,'s^-1'), n=0.253035, Ea=(209.55,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_N-4R!H->C_Ext-2R!H-R_N-5R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_N-4R!H->C_Ext-2R!H-R_N-5R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

network(
    label = 'PDepNetwork #1055',
    isomers = [
        'F[CH]C(F)(F)C1O[C]1F(2760)',
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
    label = 'PDepNetwork #1055',
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

