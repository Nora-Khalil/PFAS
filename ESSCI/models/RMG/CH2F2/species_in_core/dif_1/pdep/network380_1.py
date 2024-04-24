species(
    label = '[CH2]C(F)(F)C([O])=O(948)',
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
        HarmonicOscillator(frequencies=([249,341,381,651,928,1217,1251,3000,3100,440,815,1455,1000,180,180,180,1157.28,1158.09,1158.78],'cm^-1')),
        HinderedRotor(inertia=(0.00332538,'amu*angstrom^2'), symmetry=1, barrier=(3.16153,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.33341,'amu*angstrom^2'), symmetry=1, barrier=(53.6498,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67959,0.0592261,-0.000100773,9.94399e-08,-3.85e-11,-51461.7,21.3931], Tmin=(100,'K'), Tmax=(796.957,'K')), NASAPolynomial(coeffs=[3.99536,0.0311423,-1.69329e-05,3.38961e-09,-2.39805e-13,-51308.1,14.0261], Tmin=(796.957,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-428.505,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCCFF) + longDistanceInteraction_noncyclic(Cs(F)2-CO) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(CCOJ) + radical(Csj(Cs-F1sF1sCO)(H)(H))"""),
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
    E0 = (-349.675,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.544641,0.0361013,-4.24522e-05,2.68678e-08,-7.0049e-12,-41578.9,25.9139], Tmin=(298,'K'), Tmax=(1100,'K')), NASAPolynomial(coeffs=[7.52692,0.00816727,-3.24648e-06,5.85466e-10,-3.90185e-14,-43575.5,-14.4929], Tmin=(1100,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(-349.675,'kJ/mol'), Cp0=(33.2579,'J/mol/K'), CpInf=(133.032,'J/mol/K'), label="""CH2CF2""", comment="""Thermo library: Fluorine"""),
)

species(
    label = '[O]C(=O)C[C](F)F(949)',
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
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,190,488,555,1236,1407,180,180,180,999.422,999.424,999.502,999.741],'cm^-1')),
        HinderedRotor(inertia=(0.118687,'amu*angstrom^2'), symmetry=1, barrier=(2.72884,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118139,'amu*angstrom^2'), symmetry=1, barrier=(2.71625,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.043,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3898.35,'J/mol'), sigma=(5.74734,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=608.91 K, Pc=46.59 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70445,0.0618723,-0.000116278,1.19416e-07,-4.59784e-11,-54330.7,20.7413], Tmin=(100,'K'), Tmax=(840.051,'K')), NASAPolynomial(coeffs=[2.10069,0.0333442,-1.77665e-05,3.48499e-09,-2.41977e-13,-53457.3,24.4938], Tmin=(840.051,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-452.327,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsHH) + group(CsCsFFH) + group(Cds-OdCsOs) + radical(CCOJ) + radical(Csj(Cs-COHH)(F1s)(F1s))"""),
)

species(
    label = '[CH2][C](F)C(=O)OF(998)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {3,S}
3 O u0 p2 c0 {2,S} {6,S}
4 O u0 p2 c0 {6,D}
5 C u1 p0 c0 {1,S} {6,S} {7,S}
6 C u0 p0 c0 {3,S} {4,D} {5,S}
7 C u1 p0 c0 {5,S} {8,S} {9,S}
8 H u0 p0 c0 {7,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-154.977,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([990,1113,234,347,1316,1464,3000,3100,440,815,1455,1000,391.611,391.902,391.967,392.192,392.408,1674.28],'cm^-1')),
        HinderedRotor(inertia=(0.403845,'amu*angstrom^2'), symmetry=1, barrier=(43.749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.399776,'amu*angstrom^2'), symmetry=1, barrier=(43.7338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.40205,'amu*angstrom^2'), symmetry=1, barrier=(43.8137,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59493,0.0578276,-7.41678e-05,5.11708e-08,-1.44973e-11,-18557.2,21.6382], Tmin=(100,'K'), Tmax=(850.963,'K')), NASAPolynomial(coeffs=[9.10557,0.0225237,-1.19378e-05,2.41852e-09,-1.74681e-13,-19835.4,-13.3843], Tmin=(850.963,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-154.977,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2sCF) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(CsCOCsF1s) + radical(Csj(Cs-F1sCOH)(H)(H))"""),
)

species(
    label = 'CH2(T)(68)',
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
    label = '[O]C(=O)[C](F)F(301)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u1 p2 c0 {6,S}
4 O u0 p2 c0 {6,D}
5 C u1 p0 c0 {1,S} {2,S} {6,S}
6 C u0 p0 c0 {3,S} {4,D} {5,S}
"""),
    E0 = (-403.431,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([179,346,818,1406,1524,552.36,552.367,552.375,552.376,552.437,552.438],'cm^-1')),
        HinderedRotor(inertia=(0.000552515,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.0169,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.50118,0.0377573,-6.33905e-05,6.35209e-08,-2.57029e-11,-48472.2,17.3926], Tmin=(100,'K'), Tmax=(738.699,'K')), NASAPolynomial(coeffs=[4.33901,0.0200027,-1.14933e-05,2.38477e-09,-1.73033e-13,-48530.9,10.5238], Tmin=(738.699,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-403.431,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(CCOJ) + radical(Csj(CO-O2sO2d)(F1s)(F1s))"""),
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
    label = '[CH2]C(F)(F)[C]=O(999)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {1,S} {2,S} {5,S} {6,S}
5 C u1 p0 c0 {4,S} {7,S} {8,S}
6 C u1 p0 c0 {3,D} {4,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-243.455,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([169,350,401,508,725,1243,1223,3000,3100,440,815,1455,1000,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(0.281761,'amu*angstrom^2'), symmetry=1, barrier=(6.47823,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.282515,'amu*angstrom^2'), symmetry=1, barrier=(6.49557,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0441,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77667,0.0557146,-9.88107e-05,9.22915e-08,-3.32112e-11,-29207.4,20.2577], Tmin=(100,'K'), Tmax=(825.037,'K')), NASAPolynomial(coeffs=[6.4684,0.0200859,-1.06135e-05,2.09934e-09,-1.46847e-13,-29543.2,1.18213], Tmin=(825.037,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-243.455,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: CHOF_G4 + radical(Csj(Cs-F1sF1sCO)(H)(H)) + radical(COj(Cs-F1sF1sCs)(O2d))"""),
)

species(
    label = 'O=C1OCC1(F)F(953)',
    structure = adjacencyList("""1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {6,S} {7,S}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6 C u0 p0 c0 {3,S} {5,S} {8,S} {9,S}
7 C u0 p0 c0 {3,S} {4,D} {5,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-690.493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.45145,0.023581,2.3141e-05,-5.13415e-08,2.31834e-11,-82981.6,18.5747], Tmin=(100,'K'), Tmax=(911.73,'K')), NASAPolynomial(coeffs=[11.4033,0.0117117,-2.41857e-06,3.16177e-10,-2.13996e-14,-85752.9,-30.0321], Tmin=(911.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-690.493,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(CsCCFF) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + longDistanceInteraction_cyclic(Cs(F)2-CO) + ring(Beta-Propiolactone)"""),
)

species(
    label = '[CH2]C(F)(F)[C]1OO1(1000)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {4,S} {6,S}
4 O u0 p2 c0 {3,S} {6,S}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6 C u1 p0 c0 {3,S} {4,S} {5,S}
7 C u1 p0 c0 {5,S} {8,S} {9,S}
8 H u0 p0 c0 {7,S}
9 H u0 p0 c0 {7,S}
"""),
    E0 = (-83.2056,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77855,0.0469602,-4.6292e-05,2.27348e-08,-4.43816e-12,-9925.77,22.9292], Tmin=(100,'K'), Tmax=(1234.11,'K')), NASAPolynomial(coeffs=[11.8642,0.0142707,-6.55967e-06,1.27144e-09,-9.02274e-14,-12415.1,-27.8497], Tmin=(1234.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-83.2056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(CsCsCsFF) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + ring(O2s-O2s-Cs(C-FF)) + radical(Cs_P) + radical(Csj(Cs-CsF1sF1s)(H)(H))"""),
)

species(
    label = '[O][C]1OCC1(F)F(1001)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {5,S}
3 O u0 p2 c0 {6,S} {7,S}
4 O u1 p2 c0 {7,S}
5 C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
6 C u0 p0 c0 {3,S} {5,S} {8,S} {9,S}
7 C u1 p0 c0 {3,S} {4,S} {5,S}
8 H u0 p0 c0 {6,S}
9 H u0 p0 c0 {6,S}
"""),
    E0 = (-290.047,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.51421,0.0334035,-1.96658e-05,4.82242e-09,-4.40809e-13,-34831.8,20.2975], Tmin=(100,'K'), Tmax=(2265.8,'K')), NASAPolynomial(coeffs=[17.979,0.00847589,-4.73477e-06,8.91618e-10,-5.81154e-14,-42449.2,-68.3051], Tmin=(2265.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-290.047,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(CsCsCsFF) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + ring(Cs-O2s-Cs-Cs(F)(F)) + radical(CCOJ) + radical(Cs_P)"""),
)

species(
    label = '[O]C1([O])CC1(F)F(995)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.38289,0.0316712,-1.70962e-05,3.63448e-09,-2.83076e-13,-29241.1,15.0543], Tmin=(100,'K'), Tmax=(3111.18,'K')), NASAPolynomial(coeffs=[26.5529,0.00188163,-2.73355e-06,5.568e-10,-3.57653e-14,-43658.1,-123.025], Tmin=(3111.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-243.145,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(CsCsCsFF) + ring(Cs-Cs(F)(F)-Cs(O2)) + radical(CC(C)(O)OJ) + radical(CC(C)(O)OJ)"""),
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
    label = 'C=C(F)C([O])=O(1002)',
    structure = adjacencyList("""multiplicity 2
1 F u0 p3 c0 {4,S}
2 O u1 p2 c0 {6,S}
3 O u0 p2 c0 {6,D}
4 C u0 p0 c0 {1,S} {5,D} {6,S}
5 C u0 p0 c0 {4,D} {7,S} {8,S}
6 C u0 p0 c0 {2,S} {3,D} {4,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {5,S}
"""),
    E0 = (-225.509,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([288,410,724,839,1320,2950,3100,1380,975,1025,1650,180,180,1637.67,1639.75,1642.38,1643.3],'cm^-1')),
        HinderedRotor(inertia=(0.200618,'amu*angstrom^2'), symmetry=1, barrier=(4.61261,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (89.0451,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.6618,0.0356335,-5.66103e-05,6.05029e-08,-2.51047e-11,-27080.3,17.5644], Tmin=(100,'K'), Tmax=(806.761,'K')), NASAPolynomial(coeffs=[1.30485,0.0276006,-1.42303e-05,2.80341e-09,-1.96967e-13,-26381,26.7966], Tmin=(806.761,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-225.509,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CdCCF) + longDistanceInteraction_noncyclic(Cd(F)-CO) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsHH) + radical(CCOJ)"""),
)

species(
    label = '[O][C]=O(190)',
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
    label = '[CH2][C](F)F(187)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {4,S}
2 F u0 p3 c0 {4,S}
3 C u1 p0 c0 {4,S} {5,S} {6,S}
4 C u1 p0 c0 {1,S} {2,S} {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
"""),
    E0 = (-107.299,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,190,488,555,1236,1407],'cm^-1')),
        HinderedRotor(inertia=(0.00459569,'amu*angstrom^2'), symmetry=1, barrier=(7.47371,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (64.034,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.09702,0.01871,-1.30857e-05,4.49229e-09,-6.17295e-13,-12871.8,14.3472], Tmin=(100,'K'), Tmax=(1681.06,'K')), NASAPolynomial(coeffs=[7.67292,0.00782184,-3.37022e-06,6.39393e-10,-4.43074e-14,-14410.3,-10.1057], Tmin=(1681.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-107.299,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: Fluorine + radical(Cs_P) + radical(CsCsF1sF1s)"""),
)

species(
    label = '[O]C(=O)[C](F)CF(1003)',
    structure = adjacencyList("""multiplicity 3
1 F u0 p3 c0 {5,S}
2 F u0 p3 c0 {6,S}
3 O u1 p2 c0 {7,S}
4 O u0 p2 c0 {7,D}
5 C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
6 C u1 p0 c0 {2,S} {5,S} {7,S}
7 C u0 p0 c0 {3,S} {4,D} {6,S}
8 H u0 p0 c0 {5,S}
9 H u0 p0 c0 {5,S}
"""),
    E0 = (-446.824,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.043,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.75231,0.0346645,4.75762e-06,-9.19369e-08,9.52855e-11,-53703.2,18.4956], Tmin=(100,'K'), Tmax=(411.575,'K')), NASAPolynomial(coeffs=[4.34516,0.029713,-1.55705e-05,3.14808e-09,-2.27236e-13,-53923.5,11.1418], Tmin=(411.575,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-446.824,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(CsCCFH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + longDistanceInteraction_noncyclic(Cs(F)-CO) + group(CsCsFHH) + longDistanceInteraction_noncyclic(Cs(F)-Cs(F)) + group(Cds-OdCsOs) + radical(CCOJ) + radical(CsCOCsF1s)"""),
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
    E0 = (-217.654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-60.3359,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (244.984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (188.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (210.429,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-209.37,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (127.645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-79.1962,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-32.2943,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (97.916,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-80.0083,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-203.692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (135.087,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (-45.2151,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(F)(F)C([O])=O(948)'],
    products = ['CO2(13)', 'CH2CF2(56)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(F)(F)C([O])=O(948)'],
    products = ['[O]C(=O)C[C](F)F(949)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-R!HR!H)CJ;CsJ-HH;C] for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;CO]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2][C](F)C(=O)OF(998)'],
    products = ['[CH2]C(F)(F)C([O])=O(948)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.88952e+10,'s^-1'), n=0.725184, Ea=(189.11,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.005830439029566195, var=4.831276094152293, Tref=1000.0, N=2, data_mean=0.0, correlation='F',), comment="""Estimated from node F"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH2(T)(68)', '[O]C(=O)[C](F)F(301)'],
    products = ['[CH2]C(F)(F)C([O])=O(948)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.04495e+06,'m^3/(mol*s)'), n=0.382229, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_ter_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -1.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(6)', '[CH2]C(F)(F)[C]=O(999)'],
    products = ['[CH2]C(F)(F)C([O])=O(948)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [CO_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C(F)(F)C([O])=O(948)'],
    products = ['O=C1OCC1(F)F(953)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C(F)(F)C([O])=O(948)'],
    products = ['[CH2]C(F)(F)[C]1OO1(1000)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.55936e+11,'s^-1'), n=0.551275, Ea=(345.299,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra] for rate rule [R3_CO;carbonyl_intra_Nd;radadd_intra_O]
Euclidian distance = 2.449489742783178
family: Intra_R_Add_Endocyclic
Ea raised from 343.9 to 345.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C(F)(F)C([O])=O(948)'],
    products = ['[O][C]1OCC1(F)F(1001)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(6.54148e+08,'s^-1'), n=0.924088, Ea=(138.458,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_cs2H]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 135.3 to 138.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(F)(F)C([O])=O(948)'],
    products = ['[O]C1([O])CC1(F)F(995)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.34239e+09,'s^-1'), n=0.889391, Ea=(185.36,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonylbond_intra;radadd_intra_cs2H]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic
Ea raised from 183.5 to 185.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction10',
    reactants = ['F(37)', 'C=C(F)C([O])=O(1002)'],
    products = ['[CH2]C(F)(F)C([O])=O(948)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.575e+07,'m^3/(mol*s)'), n=3.11585e-09, Ea=(39.6824,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F',), comment="""Estimated from node Root_N-3R-inRing_N-3R->C_N-1R!H->N_N-2R!H->O_N-3BrClFNOS->Cl_N-2CNS->N_N-3BrFNOS->N_3BrFOS->F"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O][C]=O(190)', 'CH2CF2(56)'],
    products = ['[CH2]C(F)(F)C([O])=O(948)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(153.031,'m^3/(mol*s)'), n=1.16366, Ea=(27.281,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.022037706214473284, var=2.3701416358838845, Tref=1000.0, N=230, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_Sp-4R!H=3R_Sp-2R!H=1R!H"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CO2(13)', '[CH2][C](F)F(187)'],
    products = ['[CH2]C(F)(F)C([O])=O(948)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(0.000143481,'m^3/(mol*s)'), n=2.73966, Ea=(95.8938,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.24334233861756036, var=1.6451213851449848, Tref=1000.0, N=114, data_mean=0.0, correlation='Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing_4R!H->C_Sp-2R!H=1R!H',), comment="""Estimated from node Root_N-3R-inRing_Ext-3R-R_N-Sp-4R!H=3R_3R->C_Ext-1R!H-R_N-5R!H-inRing_4R!H->C_Sp-2R!H=1R!H"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O][C]=O(190)', '[CH2][C](F)F(187)'],
    products = ['[CH2]C(F)(F)C([O])=O(948)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(9.04e+06,'m^3/(mol*s)'), n=2.17087e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R',), comment="""Estimated from node Root_N-1R->H_N-1BrCClFINOPSSi->N_N-1BrCClFOS->Cl_1BrCFOS->O_Ext-1O-R_N-3R!H->O_Ext-2R-R_2R->C_Ext-2C-R_Ext-2C-R"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(F)(F)C([O])=O(948)'],
    products = ['[O]C(=O)[C](F)CF(1003)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.01526e+12,'s^-1'), n=0.18834, Ea=(172.439,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='R2F_Ext-1R!H-R_4R!H->C',), comment="""Estimated from node R2F_Ext-1R!H-R_4R!H->C
Multiplied by reaction path degeneracy 2.0"""),
)

network(
    label = 'PDepNetwork #380',
    isomers = [
        '[CH2]C(F)(F)C([O])=O(948)',
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
    label = 'PDepNetwork #380',
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

