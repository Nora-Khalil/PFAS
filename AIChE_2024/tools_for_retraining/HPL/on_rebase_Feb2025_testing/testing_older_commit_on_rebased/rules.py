#!/usr/bin/env python
# encoding: utf-8

name = "1+2_Cycloaddition/rules"
shortDesc = ""
longDesc = """

"""
entry(
    index = 1,
    label = "Root",
    kinetics = ArrheniusBM(A=(3.65062e+43,'m^3/(mol*s)'), n=-10.9254, w0=(533238,'J/mol'), E0=(211624,'J/mol'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.3583955782505395, var=162.7615503188367, Tref=1000.0, N=42, data_mean=0.0, correlation='Root',), comment="""BM rule fitted to 42 training reactions at node Root
    Total Standard Deviation in ln(k): 28.989070417350188"""),
    rank = 11,
    shortDesc = """BM rule fitted to 42 training reactions at node Root
Total Standard Deviation in ln(k): 28.989070417350188""",
    longDesc = 
"""
BM rule fitted to 42 training reactions at node Root
Total Standard Deviation in ln(k): 28.989070417350188
""",
)

entry(
    index = 2,
    label = "Root_Ext-3R-R",
    kinetics = ArrheniusBM(A=(5.9786e+44,'m^3/(mol*s)'), n=-11.3045, w0=(555517,'J/mol'), E0=(213149,'J/mol'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.4085912481084455, var=180.0845526459928, Tref=1000.0, N=29, data_mean=0.0, correlation='Root_Ext-3R-R',), comment="""BM rule fitted to 29 training reactions at node Root_Ext-3R-R
    Total Standard Deviation in ln(k): 30.441833810939148"""),
    rank = 11,
    shortDesc = """BM rule fitted to 29 training reactions at node Root_Ext-3R-R
Total Standard Deviation in ln(k): 30.441833810939148""",
    longDesc = 
"""
BM rule fitted to 29 training reactions at node Root_Ext-3R-R
Total Standard Deviation in ln(k): 30.441833810939148
""",
)

entry(
    index = 3,
    label = "Root_1R->C",
    kinetics = ArrheniusBM(A=(2.94176e+08,'m^3/(mol*s)'), n=-0.463393, w0=(476125,'J/mol'), E0=(47612.5,'J/mol'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.05349156708281672, var=2.5349342802848747, Tref=1000.0, N=12, data_mean=0.0, correlation='Root_1R->C',), comment="""BM rule fitted to 12 training reactions at node Root_1R->C
    Total Standard Deviation in ln(k): 3.326235253360758"""),
    rank = 11,
    shortDesc = """BM rule fitted to 12 training reactions at node Root_1R->C
Total Standard Deviation in ln(k): 3.326235253360758""",
    longDesc = 
"""
BM rule fitted to 12 training reactions at node Root_1R->C
Total Standard Deviation in ln(k): 3.326235253360758
""",
)

entry(
    index = 4,
    label = "Root_N-1R->C",
    kinetics = ArrheniusBM(A=(0.53862,'m^3/(mol*s)'), n=1.86213, w0=(572500,'J/mol'), E0=(98533.3,'J/mol'), Tmin=(303.03,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_N-1R->C',), comment="""BM rule fitted to 1 training reactions at node Root_N-1R->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_N-1R->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_N-1R->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 5,
    label = "Root_Ext-3R-R_Sp-4R!H-3R",
    kinetics = ArrheniusBM(A=(2.53975e+62,'m^3/(mol*s)'), n=-16.2979, w0=(543529,'J/mol'), E0=(213324,'J/mol'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.8837658861059445, var=335.0753577820745, Tref=1000.0, N=17, data_mean=0.0, correlation='Root_Ext-3R-R_Sp-4R!H-3R',), comment="""BM rule fitted to 17 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R
    Total Standard Deviation in ln(k): 41.429883448650315"""),
    rank = 11,
    shortDesc = """BM rule fitted to 17 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R
Total Standard Deviation in ln(k): 41.429883448650315""",
    longDesc = 
"""
BM rule fitted to 17 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R
Total Standard Deviation in ln(k): 41.429883448650315
""",
)

entry(
    index = 6,
    label = "Root_Ext-3R-R_N-Sp-4R!H-3R",
    kinetics = ArrheniusBM(A=(1.06287e-07,'m^3/(mol*s)'), n=3.48388, w0=(572500,'J/mol'), E0=(124514,'J/mol'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.019764276166514334, var=9.541692004373376, Tref=1000.0, N=12, data_mean=0.0, correlation='Root_Ext-3R-R_N-Sp-4R!H-3R',), comment="""BM rule fitted to 12 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R
    Total Standard Deviation in ln(k): 6.242211330847662"""),
    rank = 11,
    shortDesc = """BM rule fitted to 12 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R
Total Standard Deviation in ln(k): 6.242211330847662""",
    longDesc = 
"""
BM rule fitted to 12 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R
Total Standard Deviation in ln(k): 6.242211330847662
""",
)

entry(
    index = 7,
    label = "Root_1R->C_Ext-1C-R",
    kinetics = ArrheniusBM(A=(2.5716e+09,'m^3/(mol*s)'), n=-0.704858, w0=(474200,'J/mol'), E0=(47420,'J/mol'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.003178602591253677, var=0.5008544564983712, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_1R->C_Ext-1C-R',), comment="""BM rule fitted to 5 training reactions at node Root_1R->C_Ext-1C-R
    Total Standard Deviation in ln(k): 1.4267589342073508"""),
    rank = 11,
    shortDesc = """BM rule fitted to 5 training reactions at node Root_1R->C_Ext-1C-R
Total Standard Deviation in ln(k): 1.4267589342073508""",
    longDesc = 
"""
BM rule fitted to 5 training reactions at node Root_1R->C_Ext-1C-R
Total Standard Deviation in ln(k): 1.4267589342073508
""",
)

entry(
    index = 8,
    label = "Root_1R->C_Ext-2R-R",
    kinetics = ArrheniusBM(A=(3.33384e+09,'m^3/(mol*s)'), n=-0.775738, w0=(480000,'J/mol'), E0=(48000,'J/mol'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.03181634044488337, var=0.1433731634963656, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_1R->C_Ext-2R-R',), comment="""BM rule fitted to 4 training reactions at node Root_1R->C_Ext-2R-R
    Total Standard Deviation in ln(k): 0.8390264519529201"""),
    rank = 11,
    shortDesc = """BM rule fitted to 4 training reactions at node Root_1R->C_Ext-2R-R
Total Standard Deviation in ln(k): 0.8390264519529201""",
    longDesc = 
"""
BM rule fitted to 4 training reactions at node Root_1R->C_Ext-2R-R
Total Standard Deviation in ln(k): 0.8390264519529201
""",
)

entry(
    index = 9,
    label = "Root_1R->C_Sp-2R=1C",
    kinetics = ArrheniusBM(A=(0.00650771,'m^3/(mol*s)'), n=2.42295, w0=(480000,'J/mol'), E0=(48000,'J/mol'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.587997258796586, var=0.16964704494037178, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_1R->C_Sp-2R=1C',), comment="""BM rule fitted to 2 training reactions at node Root_1R->C_Sp-2R=1C
    Total Standard Deviation in ln(k): 4.815657794528946"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_1R->C_Sp-2R=1C
Total Standard Deviation in ln(k): 4.815657794528946""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_1R->C_Sp-2R=1C
Total Standard Deviation in ln(k): 4.815657794528946
""",
)

entry(
    index = 10,
    label = "Root_1R->C_N-Sp-2R=1C",
    kinetics = ArrheniusBM(A=(1.77e+09,'m^3/(mol*s)'), n=-0.662, w0=(462500,'J/mol'), E0=(46250,'J/mol'), Tmin=(200,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->C_N-Sp-2R=1C',), comment="""BM rule fitted to 1 training reactions at node Root_1R->C_N-Sp-2R=1C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_1R->C_N-Sp-2R=1C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_1R->C_N-Sp-2R=1C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 11,
    label = "Root_Ext-3R-R_Sp-4R!H-3R_4R!H->O",
    kinetics = ArrheniusBM(A=(2.39334e-46,'m^3/(mol*s)'), n=14.7713, w0=(572500,'J/mol'), E0=(-110286,'J/mol'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-17.24787716633618, var=46.96062662985535, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_Ext-3R-R_Sp-4R!H-3R_4R!H->O',), comment="""BM rule fitted to 5 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_4R!H->O
    Total Standard Deviation in ln(k): 57.07438804600418"""),
    rank = 11,
    shortDesc = """BM rule fitted to 5 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_4R!H->O
Total Standard Deviation in ln(k): 57.07438804600418""",
    longDesc = 
"""
BM rule fitted to 5 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_4R!H->O
Total Standard Deviation in ln(k): 57.07438804600418
""",
)

entry(
    index = 12,
    label = "Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O",
    kinetics = ArrheniusBM(A=(1.52958e-35,'m^3/(mol*s)'), n=11.6327, w0=(531458,'J/mol'), E0=(31554.8,'J/mol'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-1.1734088086985555, var=75.48094170605887, Tref=1000.0, N=12, data_mean=0.0, correlation='Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O',), comment="""BM rule fitted to 12 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O
    Total Standard Deviation in ln(k): 20.36535527659792"""),
    rank = 11,
    shortDesc = """BM rule fitted to 12 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O
Total Standard Deviation in ln(k): 20.36535527659792""",
    longDesc = 
"""
BM rule fitted to 12 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O
Total Standard Deviation in ln(k): 20.36535527659792
""",
)

entry(
    index = 13,
    label = "Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R",
    kinetics = ArrheniusBM(A=(1.74208e-07,'m^3/(mol*s)'), n=3.40879, w0=(572500,'J/mol'), E0=(124664,'J/mol'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.037589497325199166, var=10.587001202588679, Tref=1000.0, N=11, data_mean=0.0, correlation='Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R',), comment="""BM rule fitted to 11 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R
    Total Standard Deviation in ln(k): 6.617387277631007"""),
    rank = 11,
    shortDesc = """BM rule fitted to 11 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R
Total Standard Deviation in ln(k): 6.617387277631007""",
    longDesc = 
"""
BM rule fitted to 11 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R
Total Standard Deviation in ln(k): 6.617387277631007
""",
)

entry(
    index = 14,
    label = "Root_1R->C_Ext-1C-R_Ext-1C-R",
    kinetics = ArrheniusBM(A=(3.18e+07,'m^3/(mol*s)'), n=0, w0=(486000,'J/mol'), E0=(48600,'J/mol'), Tmin=(298,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->C_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 1 training reactions at node Root_1R->C_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_1R->C_Ext-1C-R_Ext-1C-R
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_1R->C_Ext-1C-R_Ext-1C-R
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 15,
    label = "Root_1R->C_Ext-1C-R_Sp-2R=1C",
    kinetics = ArrheniusBM(A=(9.41381e+08,'m^3/(mol*s)'), n=-0.607357, w0=(480000,'J/mol'), E0=(48000,'J/mol'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08355106414541179, var=0.031111968581513858, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_1R->C_Ext-1C-R_Sp-2R=1C',), comment="""BM rule fitted to 2 training reactions at node Root_1R->C_Ext-1C-R_Sp-2R=1C
    Total Standard Deviation in ln(k): 0.5635342003474927"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_1R->C_Ext-1C-R_Sp-2R=1C
Total Standard Deviation in ln(k): 0.5635342003474927""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_1R->C_Ext-1C-R_Sp-2R=1C
Total Standard Deviation in ln(k): 0.5635342003474927
""",
)

entry(
    index = 16,
    label = "Root_1R->C_Ext-1C-R_N-Sp-2R=1C",
    kinetics = ArrheniusBM(A=(4.63231e+09,'m^3/(mol*s)'), n=-0.7664, w0=(462500,'J/mol'), E0=(46250,'J/mol'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.0009309671137264789, var=1.1655949953965052, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_1R->C_Ext-1C-R_N-Sp-2R=1C',), comment="""BM rule fitted to 2 training reactions at node Root_1R->C_Ext-1C-R_N-Sp-2R=1C
    Total Standard Deviation in ln(k): 2.1667057286443416"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_1R->C_Ext-1C-R_N-Sp-2R=1C
Total Standard Deviation in ln(k): 2.1667057286443416""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_1R->C_Ext-1C-R_N-Sp-2R=1C
Total Standard Deviation in ln(k): 2.1667057286443416
""",
)

entry(
    index = 17,
    label = "Root_1R->C_Ext-2R-R_3R->C",
    kinetics = ArrheniusBM(A=(4.13595e+09,'m^3/(mol*s)'), n=-0.801251, w0=(474000,'J/mol'), E0=(47400,'J/mol'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.0395530227295251, var=0.0335984334446083, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_1R->C_Ext-2R-R_3R->C',), comment="""BM rule fitted to 2 training reactions at node Root_1R->C_Ext-2R-R_3R->C
    Total Standard Deviation in ln(k): 0.46684489712018756"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_1R->C_Ext-2R-R_3R->C
Total Standard Deviation in ln(k): 0.46684489712018756""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_1R->C_Ext-2R-R_3R->C
Total Standard Deviation in ln(k): 0.46684489712018756
""",
)

entry(
    index = 18,
    label = "Root_1R->C_Ext-2R-R_N-3R->C",
    kinetics = ArrheniusBM(A=(3.86082e+06,'m^3/(mol*s)'), n=0.0243327, w0=(486000,'J/mol'), E0=(48600,'J/mol'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.02516094828390788, var=1.7607258136569388, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_1R->C_Ext-2R-R_N-3R->C',), comment="""BM rule fitted to 2 training reactions at node Root_1R->C_Ext-2R-R_N-3R->C
    Total Standard Deviation in ln(k): 2.7233484267253916"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_1R->C_Ext-2R-R_N-3R->C
Total Standard Deviation in ln(k): 2.7233484267253916""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_1R->C_Ext-2R-R_N-3R->C
Total Standard Deviation in ln(k): 2.7233484267253916
""",
)

entry(
    index = 19,
    label = "Root_1R->C_Sp-2R=1C_3R->C",
    kinetics = ArrheniusBM(A=(1.98e+06,'m^3/(mol*s)'), n=0, w0=(474000,'J/mol'), E0=(47400,'J/mol'), Tmin=(296,'K'), Tmax=(728,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->C_Sp-2R=1C_3R->C',), comment="""BM rule fitted to 1 training reactions at node Root_1R->C_Sp-2R=1C_3R->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_1R->C_Sp-2R=1C_3R->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_1R->C_Sp-2R=1C_3R->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 20,
    label = "Root_1R->C_Sp-2R=1C_N-3R->C",
    kinetics = ArrheniusBM(A=(700000,'m^3/(mol*s)'), n=0, w0=(486000,'J/mol'), E0=(48600,'J/mol'), Tmin=(300,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->C_Sp-2R=1C_N-3R->C',), comment="""BM rule fitted to 1 training reactions at node Root_1R->C_Sp-2R=1C_N-3R->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_1R->C_Sp-2R=1C_N-3R->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_1R->C_Sp-2R=1C_N-3R->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 21,
    label = "Root_Ext-3R-R_Sp-4R!H-3R_4R!H->O_Ext-4O-R_Ext-3R-R_Ext-5R!H-R",
    kinetics = ArrheniusBM(A=(6.43209e+53,'m^3/(mol*s)'), n=-15.5793, w0=(572500,'J/mol'), E0=(50764.7,'J/mol'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.831024435251295, var=0.34690475935680176, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_Ext-3R-R_Sp-4R!H-3R_4R!H->O_Ext-4O-R_Ext-3R-R_Ext-5R!H-R',), comment="""BM rule fitted to 3 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_4R!H->O_Ext-4O-R_Ext-3R-R_Ext-5R!H-R
    Total Standard Deviation in ln(k): 5.781325229399296"""),
    rank = 11,
    shortDesc = """BM rule fitted to 3 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_4R!H->O_Ext-4O-R_Ext-3R-R_Ext-5R!H-R
Total Standard Deviation in ln(k): 5.781325229399296""",
    longDesc = 
"""
BM rule fitted to 3 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_4R!H->O_Ext-4O-R_Ext-3R-R_Ext-5R!H-R
Total Standard Deviation in ln(k): 5.781325229399296
""",
)

entry(
    index = 22,
    label = "Root_Ext-3R-R_Sp-4R!H-3R_4R!H->O_Ext-4O-R_Ext-5R!H-R",
    kinetics = ArrheniusBM(A=(5.34467e-08,'m^3/(mol*s)'), n=4.75559, w0=(572500,'J/mol'), E0=(57250,'J/mol'), Tmin=(303.03,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-3R-R_Sp-4R!H-3R_4R!H->O_Ext-4O-R_Ext-5R!H-R',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_4R!H->O_Ext-4O-R_Ext-5R!H-R
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_4R!H->O_Ext-4O-R_Ext-5R!H-R
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_4R!H->O_Ext-4O-R_Ext-5R!H-R
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 23,
    label = "Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C",
    kinetics = ArrheniusBM(A=(1.09626e+10,'m^3/(mol*s)'), n=-1.18034, w0=(490417,'J/mol'), E0=(144079,'J/mol'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.329607791278391, var=5.390846928468328, Tref=1000.0, N=6, data_mean=0.0, correlation='Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C',), comment="""BM rule fitted to 6 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C
    Total Standard Deviation in ln(k): 5.482793765866567"""),
    rank = 11,
    shortDesc = """BM rule fitted to 6 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C
Total Standard Deviation in ln(k): 5.482793765866567""",
    longDesc = 
"""
BM rule fitted to 6 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C
Total Standard Deviation in ln(k): 5.482793765866567
""",
)

entry(
    index = 24,
    label = "Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C",
    kinetics = ArrheniusBM(A=(2.44892e-39,'m^3/(mol*s)'), n=12.717, w0=(572500,'J/mol'), E0=(1343.97,'J/mol'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-1.5181274212804807, var=82.16873793629212, Tref=1000.0, N=6, data_mean=0.0, correlation='Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C',), comment="""BM rule fitted to 6 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C
    Total Standard Deviation in ln(k): 21.98670723698378"""),
    rank = 11,
    shortDesc = """BM rule fitted to 6 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C
Total Standard Deviation in ln(k): 21.98670723698378""",
    longDesc = 
"""
BM rule fitted to 6 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C
Total Standard Deviation in ln(k): 21.98670723698378
""",
)

entry(
    index = 25,
    label = "Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_5R!H->F",
    kinetics = ArrheniusBM(A=(0.0160101,'m^3/(mol*s)'), n=2.09578, w0=(572500,'J/mol'), E0=(116559,'J/mol'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.060407312939327225, var=3.4113773501814233, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_5R!H->F',), comment="""BM rule fitted to 4 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_5R!H->F
    Total Standard Deviation in ln(k): 3.8545056797289043"""),
    rank = 11,
    shortDesc = """BM rule fitted to 4 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_5R!H->F
Total Standard Deviation in ln(k): 3.8545056797289043""",
    longDesc = 
"""
BM rule fitted to 4 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_5R!H->F
Total Standard Deviation in ln(k): 3.8545056797289043
""",
)

entry(
    index = 26,
    label = "Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F",
    kinetics = ArrheniusBM(A=(1.35394e-17,'m^3/(mol*s)'), n=6.24951, w0=(572500,'J/mol'), E0=(110567,'J/mol'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.16823401916699962, var=3.254104310757667, Tref=1000.0, N=7, data_mean=0.0, correlation='Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F',), comment="""BM rule fitted to 7 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F
    Total Standard Deviation in ln(k): 4.039067430295244"""),
    rank = 11,
    shortDesc = """BM rule fitted to 7 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F
Total Standard Deviation in ln(k): 4.039067430295244""",
    longDesc = 
"""
BM rule fitted to 7 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F
Total Standard Deviation in ln(k): 4.039067430295244
""",
)

entry(
    index = 27,
    label = "Root_1R->C_Ext-1C-R_Sp-2R=1C_Ext-2R-R",
    kinetics = ArrheniusBM(A=(1.54e+07,'m^3/(mol*s)'), n=0, w0=(486000,'J/mol'), E0=(48600,'J/mol'), Tmin=(298,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->C_Ext-1C-R_Sp-2R=1C_Ext-2R-R',), comment="""BM rule fitted to 1 training reactions at node Root_1R->C_Ext-1C-R_Sp-2R=1C_Ext-2R-R
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_1R->C_Ext-1C-R_Sp-2R=1C_Ext-2R-R
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_1R->C_Ext-1C-R_Sp-2R=1C_Ext-2R-R
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 28,
    label = "Root_1R->C_Ext-1C-R_N-Sp-2R=1C_Ext-2R-R",
    kinetics = ArrheniusBM(A=(4.7e+09,'m^3/(mol*s)'), n=-0.823, w0=(462500,'J/mol'), E0=(46250,'J/mol'), Tmin=(200,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->C_Ext-1C-R_N-Sp-2R=1C_Ext-2R-R',), comment="""BM rule fitted to 1 training reactions at node Root_1R->C_Ext-1C-R_N-Sp-2R=1C_Ext-2R-R
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_1R->C_Ext-1C-R_N-Sp-2R=1C_Ext-2R-R
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_1R->C_Ext-1C-R_N-Sp-2R=1C_Ext-2R-R
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 29,
    label = "Root_1R->C_Ext-2R-R_3R->C_Ext-4R!H-R",
    kinetics = ArrheniusBM(A=(1.85e+09,'m^3/(mol*s)'), n=-0.7, w0=(474000,'J/mol'), E0=(47400,'J/mol'), Tmin=(200,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->C_Ext-2R-R_3R->C_Ext-4R!H-R',), comment="""BM rule fitted to 1 training reactions at node Root_1R->C_Ext-2R-R_3R->C_Ext-4R!H-R
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_1R->C_Ext-2R-R_3R->C_Ext-4R!H-R
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_1R->C_Ext-2R-R_3R->C_Ext-4R!H-R
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 30,
    label = "Root_1R->C_Ext-2R-R_N-3R->C_Ext-2R-R",
    kinetics = ArrheniusBM(A=(7.6e+06,'m^3/(mol*s)'), n=0, w0=(486000,'J/mol'), E0=(48600,'J/mol'), Tmin=(298,'K'), Tmax=(410,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->C_Ext-2R-R_N-3R->C_Ext-2R-R',), comment="""BM rule fitted to 1 training reactions at node Root_1R->C_Ext-2R-R_N-3R->C_Ext-2R-R
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_1R->C_Ext-2R-R_N-3R->C_Ext-2R-R
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_1R->C_Ext-2R-R_N-3R->C_Ext-2R-R
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 31,
    label = "Root_Ext-3R-R_Sp-4R!H-3R_4R!H->O_Ext-4O-R_Ext-3R-R_Ext-5R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-8R!H-R_8R!H->C",
    kinetics = ArrheniusBM(A=(8.37235e+53,'m^3/(mol*s)'), n=-15.5905, w0=(572500,'J/mol'), E0=(54687.9,'J/mol'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-19.812270417415302, var=0.19227652450927105, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_Ext-3R-R_Sp-4R!H-3R_4R!H->O_Ext-4O-R_Ext-3R-R_Ext-5R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-8R!H-R_8R!H->C',), comment="""BM rule fitted to 2 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_4R!H->O_Ext-4O-R_Ext-3R-R_Ext-5R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-8R!H-R_8R!H->C
    Total Standard Deviation in ln(k): 50.658637156539335"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_4R!H->O_Ext-4O-R_Ext-3R-R_Ext-5R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-8R!H-R_8R!H->C
Total Standard Deviation in ln(k): 50.658637156539335""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_4R!H->O_Ext-4O-R_Ext-3R-R_Ext-5R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-8R!H-R_8R!H->C
Total Standard Deviation in ln(k): 50.658637156539335
""",
)

entry(
    index = 32,
    label = "Root_Ext-3R-R_Sp-4R!H-3R_4R!H->O_Ext-4O-R_Ext-3R-R_Ext-5R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-8R!H-R_N-8R!H->C",
    kinetics = ArrheniusBM(A=(9.85919e-16,'m^3/(mol*s)'), n=4.09382, w0=(572500,'J/mol'), E0=(57250,'J/mol'), Tmin=(303.03,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-3R-R_Sp-4R!H-3R_4R!H->O_Ext-4O-R_Ext-3R-R_Ext-5R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-8R!H-R_N-8R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_4R!H->O_Ext-4O-R_Ext-3R-R_Ext-5R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-8R!H-R_N-8R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_4R!H->O_Ext-4O-R_Ext-3R-R_Ext-5R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-8R!H-R_N-8R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_4R!H->O_Ext-4O-R_Ext-3R-R_Ext-5R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-8R!H-R_N-8R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 33,
    label = "Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_4BrCClF->F",
    kinetics = ArrheniusBM(A=(2.80502e-20,'m^3/(mol*s)'), n=7.14725, w0=(498625,'J/mol'), E0=(38962.2,'J/mol'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.4121423006753344, var=3.114469141446709, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_4BrCClF->F',), comment="""BM rule fitted to 4 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_4BrCClF->F
    Total Standard Deviation in ln(k): 4.573461541974407"""),
    rank = 11,
    shortDesc = """BM rule fitted to 4 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_4BrCClF->F
Total Standard Deviation in ln(k): 4.573461541974407""",
    longDesc = 
"""
BM rule fitted to 4 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_4BrCClF->F
Total Standard Deviation in ln(k): 4.573461541974407
""",
)

entry(
    index = 34,
    label = "Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_N-4BrCClF->F",
    kinetics = ArrheniusBM(A=(2.95403e-06,'m^3/(mol*s)'), n=3.41638, w0=(474000,'J/mol'), E0=(47400,'J/mol'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.4682509003578339, var=1.625277061065213, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_N-4BrCClF->F',), comment="""BM rule fitted to 2 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_N-4BrCClF->F
    Total Standard Deviation in ln(k): 3.73227346956574"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_N-4BrCClF->F
Total Standard Deviation in ln(k): 3.73227346956574""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_N-4BrCClF->F
Total Standard Deviation in ln(k): 3.73227346956574
""",
)

entry(
    index = 35,
    label = "Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-3R-R",
    kinetics = ArrheniusBM(A=(0.00524266,'m^3/(mol*s)'), n=2.4366, w0=(572500,'J/mol'), E0=(122370,'J/mol'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.8163698570402267, var=195.64042736584582, Tref=1000.0, N=4, data_mean=0.0, correlation='Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-3R-R',), comment="""BM rule fitted to 4 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-3R-R
    Total Standard Deviation in ln(k): 30.09171524367362"""),
    rank = 11,
    shortDesc = """BM rule fitted to 4 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-3R-R
Total Standard Deviation in ln(k): 30.09171524367362""",
    longDesc = 
"""
BM rule fitted to 4 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-3R-R
Total Standard Deviation in ln(k): 30.09171524367362
""",
)

entry(
    index = 36,
    label = "Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-4BrCClF-R",
    kinetics = ArrheniusBM(A=(1.8498e-07,'m^3/(mol*s)'), n=3.40152, w0=(572500,'J/mol'), E0=(144447,'J/mol'), Tmin=(303.03,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-4BrCClF-R',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-4BrCClF-R
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-4BrCClF-R
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-4BrCClF-R
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 37,
    label = "Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_5R!H->F_Ext-2R-R_6R!H->C",
    kinetics = ArrheniusBM(A=(2.03833e-05,'m^3/(mol*s)'), n=2.91521, w0=(572500,'J/mol'), E0=(114127,'J/mol'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.006177771675028784, var=0.2821070671969892, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_5R!H->F_Ext-2R-R_6R!H->C',), comment="""BM rule fitted to 3 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_5R!H->F_Ext-2R-R_6R!H->C
    Total Standard Deviation in ln(k): 1.080312060508416"""),
    rank = 11,
    shortDesc = """BM rule fitted to 3 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_5R!H->F_Ext-2R-R_6R!H->C
Total Standard Deviation in ln(k): 1.080312060508416""",
    longDesc = 
"""
BM rule fitted to 3 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_5R!H->F_Ext-2R-R_6R!H->C
Total Standard Deviation in ln(k): 1.080312060508416
""",
)

entry(
    index = 38,
    label = "Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_5R!H->F_Ext-2R-R_N-6R!H->C",
    kinetics = ArrheniusBM(A=(0.000352549,'m^3/(mol*s)'), n=2.48352, w0=(572500,'J/mol'), E0=(87444.7,'J/mol'), Tmin=(303.03,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_5R!H->F_Ext-2R-R_N-6R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_5R!H->F_Ext-2R-R_N-6R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_5R!H->F_Ext-2R-R_N-6R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_5R!H->F_Ext-2R-R_N-6R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 39,
    label = "Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R",
    kinetics = ArrheniusBM(A=(7.76196e-43,'m^3/(mol*s)'), n=13.5623, w0=(572500,'J/mol'), E0=(26099.1,'J/mol'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.8786799159153253, var=0.3634958394327543, Tref=1000.0, N=5, data_mean=0.0, correlation='Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R',), comment="""BM rule fitted to 5 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R
    Total Standard Deviation in ln(k): 3.4164056122994793"""),
    rank = 11,
    shortDesc = """BM rule fitted to 5 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R
Total Standard Deviation in ln(k): 3.4164056122994793""",
    longDesc = 
"""
BM rule fitted to 5 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R
Total Standard Deviation in ln(k): 3.4164056122994793
""",
)

entry(
    index = 40,
    label = "Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-5BrCClINOPSSi-R",
    kinetics = ArrheniusBM(A=(4.4056e-09,'m^3/(mol*s)'), n=3.70948, w0=(572500,'J/mol'), E0=(135210,'J/mol'), Tmin=(303.03,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-5BrCClINOPSSi-R',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-5BrCClINOPSSi-R
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-5BrCClINOPSSi-R
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-5BrCClINOPSSi-R
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 41,
    label = "Root_Ext-3R-R_Sp-4R!H-3R_4R!H->O_Ext-4O-R_Ext-3R-R_Ext-5R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-8R!H-R_8R!H->C_Ext-8C-R_9R!H->C",
    kinetics = ArrheniusBM(A=(2.09967e-13,'m^3/(mol*s)'), n=3.46188, w0=(572500,'J/mol'), E0=(57250,'J/mol'), Tmin=(303.03,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-3R-R_Sp-4R!H-3R_4R!H->O_Ext-4O-R_Ext-3R-R_Ext-5R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-8R!H-R_8R!H->C_Ext-8C-R_9R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_4R!H->O_Ext-4O-R_Ext-3R-R_Ext-5R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-8R!H-R_8R!H->C_Ext-8C-R_9R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_4R!H->O_Ext-4O-R_Ext-3R-R_Ext-5R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-8R!H-R_8R!H->C_Ext-8C-R_9R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_4R!H->O_Ext-4O-R_Ext-3R-R_Ext-5R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-8R!H-R_8R!H->C_Ext-8C-R_9R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 42,
    label = "Root_Ext-3R-R_Sp-4R!H-3R_4R!H->O_Ext-4O-R_Ext-3R-R_Ext-5R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-8R!H-R_8R!H->C_Ext-8C-R_N-9R!H->C",
    kinetics = ArrheniusBM(A=(5.99786e-14,'m^3/(mol*s)'), n=3.67717, w0=(572500,'J/mol'), E0=(57250,'J/mol'), Tmin=(303.03,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-3R-R_Sp-4R!H-3R_4R!H->O_Ext-4O-R_Ext-3R-R_Ext-5R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-8R!H-R_8R!H->C_Ext-8C-R_N-9R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_4R!H->O_Ext-4O-R_Ext-3R-R_Ext-5R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-8R!H-R_8R!H->C_Ext-8C-R_N-9R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_4R!H->O_Ext-4O-R_Ext-3R-R_Ext-5R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-8R!H-R_8R!H->C_Ext-8C-R_N-9R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_4R!H->O_Ext-4O-R_Ext-3R-R_Ext-5R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-7R!H-R_Ext-5R!H-R_Ext-5R!H-R_Ext-8R!H-R_8R!H->C_Ext-8C-R_N-9R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 43,
    label = "Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_4BrCClF->F_Ext-1C-R",
    kinetics = ArrheniusBM(A=(4.10298e-06,'m^3/(mol*s)'), n=3.07477, w0=(474000,'J/mol'), E0=(76277.8,'J/mol'), Tmin=(298,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_4BrCClF->F_Ext-1C-R',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_4BrCClF->F_Ext-1C-R
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_4BrCClF->F_Ext-1C-R
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_4BrCClF->F_Ext-1C-R
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 44,
    label = "Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_4BrCClF->F_Ext-2R-R",
    kinetics = ArrheniusBM(A=(1.01799e-05,'m^3/(mol*s)'), n=2.87159, w0=(474000,'J/mol'), E0=(98608.6,'J/mol'), Tmin=(298,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_4BrCClF->F_Ext-2R-R',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_4BrCClF->F_Ext-2R-R
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_4BrCClF->F_Ext-2R-R
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_4BrCClF->F_Ext-2R-R
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 45,
    label = "Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_4BrCClF->F_2R->C",
    kinetics = ArrheniusBM(A=(2.07011e-05,'m^3/(mol*s)'), n=2.98446, w0=(474000,'J/mol'), E0=(47400,'J/mol'), Tmin=(298,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_4BrCClF->F_2R->C',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_4BrCClF->F_2R->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_4BrCClF->F_2R->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_4BrCClF->F_2R->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 46,
    label = "Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_4BrCClF->F_N-2R->C",
    kinetics = ArrheniusBM(A=(4.81575e-05,'m^3/(mol*s)'), n=2.76922, w0=(572500,'J/mol'), E0=(107328,'J/mol'), Tmin=(298,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_4BrCClF->F_N-2R->C',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_4BrCClF->F_N-2R->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_4BrCClF->F_N-2R->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_4BrCClF->F_N-2R->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 47,
    label = "Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_N-4BrCClF->F_4BrCCl->Br",
    kinetics = ArrheniusBM(A=(0.000145611,'m^3/(mol*s)'), n=2.95653, w0=(474000,'J/mol'), E0=(47400,'J/mol'), Tmin=(298,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_N-4BrCClF->F_4BrCCl->Br',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_N-4BrCClF->F_4BrCCl->Br
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_N-4BrCClF->F_4BrCCl->Br
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_N-4BrCClF->F_4BrCCl->Br
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 48,
    label = "Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_N-4BrCClF->F_N-4BrCCl->Br",
    kinetics = ArrheniusBM(A=(8.68219e-05,'m^3/(mol*s)'), n=2.97056, w0=(474000,'J/mol'), E0=(47400,'J/mol'), Tmin=(298,'K'), Tmax=(2500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_N-4BrCClF->F_N-4BrCCl->Br',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_N-4BrCClF->F_N-4BrCCl->Br
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_N-4BrCClF->F_N-4BrCCl->Br
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_1R->C_N-4BrCClF->F_N-4BrCCl->Br
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 49,
    label = "Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-3R-R_4BrCClF->C",
    kinetics = ArrheniusBM(A=(2644.15,'m^3/(mol*s)'), n=0.797762, w0=(572500,'J/mol'), E0=(149630,'J/mol'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.17830343602060747, var=294.759867097508, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-3R-R_4BrCClF->C',), comment="""BM rule fitted to 3 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-3R-R_4BrCClF->C
    Total Standard Deviation in ln(k): 34.866437250665626"""),
    rank = 11,
    shortDesc = """BM rule fitted to 3 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-3R-R_4BrCClF->C
Total Standard Deviation in ln(k): 34.866437250665626""",
    longDesc = 
"""
BM rule fitted to 3 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-3R-R_4BrCClF->C
Total Standard Deviation in ln(k): 34.866437250665626
""",
)

entry(
    index = 50,
    label = "Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-3R-R_N-4BrCClF->C",
    kinetics = ArrheniusBM(A=(0.00491008,'m^3/(mol*s)'), n=2.38401, w0=(572500,'J/mol'), E0=(80383.6,'J/mol'), Tmin=(303.03,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-3R-R_N-4BrCClF->C',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-3R-R_N-4BrCClF->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-3R-R_N-4BrCClF->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-3R-R_N-4BrCClF->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 51,
    label = "Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_5R!H->F_Ext-2R-R_6R!H->C_Ext-6C-R_7R!H->C",
    kinetics = ArrheniusBM(A=(2.6436e-05,'m^3/(mol*s)'), n=2.88986, w0=(572500,'J/mol'), E0=(114713,'J/mol'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.0002285734005305904, var=0.984252878074767, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_5R!H->F_Ext-2R-R_6R!H->C_Ext-6C-R_7R!H->C',), comment="""BM rule fitted to 2 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_5R!H->F_Ext-2R-R_6R!H->C_Ext-6C-R_7R!H->C
    Total Standard Deviation in ln(k): 1.9894623447118933"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_5R!H->F_Ext-2R-R_6R!H->C_Ext-6C-R_7R!H->C
Total Standard Deviation in ln(k): 1.9894623447118933""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_5R!H->F_Ext-2R-R_6R!H->C_Ext-6C-R_7R!H->C
Total Standard Deviation in ln(k): 1.9894623447118933
""",
)

entry(
    index = 52,
    label = "Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_5R!H->F_Ext-2R-R_6R!H->C_Ext-6C-R_N-7R!H->C",
    kinetics = ArrheniusBM(A=(2.81961e-05,'m^3/(mol*s)'), n=2.73363, w0=(572500,'J/mol'), E0=(106031,'J/mol'), Tmin=(303.03,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_5R!H->F_Ext-2R-R_6R!H->C_Ext-6C-R_N-7R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_5R!H->F_Ext-2R-R_6R!H->C_Ext-6C-R_N-7R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_5R!H->F_Ext-2R-R_6R!H->C_Ext-6C-R_N-7R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_5R!H->F_Ext-2R-R_6R!H->C_Ext-6C-R_N-7R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 53,
    label = "Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R_Ext-6R!H-R_5BrCClINOPSSi->O",
    kinetics = ArrheniusBM(A=(8.17503e-12,'m^3/(mol*s)'), n=4.475, w0=(572500,'J/mol'), E0=(104343,'J/mol'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-0.034757139199668025, var=0.09256908687007695, Tref=1000.0, N=3, data_mean=0.0, correlation='Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R_Ext-6R!H-R_5BrCClINOPSSi->O',), comment="""BM rule fitted to 3 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R_Ext-6R!H-R_5BrCClINOPSSi->O
    Total Standard Deviation in ln(k): 0.6972735176051794"""),
    rank = 11,
    shortDesc = """BM rule fitted to 3 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R_Ext-6R!H-R_5BrCClINOPSSi->O
Total Standard Deviation in ln(k): 0.6972735176051794""",
    longDesc = 
"""
BM rule fitted to 3 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R_Ext-6R!H-R_5BrCClINOPSSi->O
Total Standard Deviation in ln(k): 0.6972735176051794
""",
)

entry(
    index = 54,
    label = "Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R_Ext-6R!H-R_N-5BrCClINOPSSi->O",
    kinetics = ArrheniusBM(A=(2.27425e-12,'m^3/(mol*s)'), n=4.81028, w0=(572500,'J/mol'), E0=(141758,'J/mol'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=1.405864762151905e-17, var=0.11495646616874614, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R_Ext-6R!H-R_N-5BrCClINOPSSi->O',), comment="""BM rule fitted to 2 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R_Ext-6R!H-R_N-5BrCClINOPSSi->O
    Total Standard Deviation in ln(k): 0.6797100508055892"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R_Ext-6R!H-R_N-5BrCClINOPSSi->O
Total Standard Deviation in ln(k): 0.6797100508055892""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R_Ext-6R!H-R_N-5BrCClINOPSSi->O
Total Standard Deviation in ln(k): 0.6797100508055892
""",
)

entry(
    index = 55,
    label = "Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-3R-R_4BrCClF->C_Ext-4C-R_6R!H->C",
    kinetics = ArrheniusBM(A=(3.2922e-05,'m^3/(mol*s)'), n=3.02248, w0=(572500,'J/mol'), E0=(159345,'J/mol'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.003984533462247416, var=701.6587757381756, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-3R-R_4BrCClF->C_Ext-4C-R_6R!H->C',), comment="""BM rule fitted to 2 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-3R-R_4BrCClF->C_Ext-4C-R_6R!H->C
    Total Standard Deviation in ln(k): 53.11312306256407"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-3R-R_4BrCClF->C_Ext-4C-R_6R!H->C
Total Standard Deviation in ln(k): 53.11312306256407""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-3R-R_4BrCClF->C_Ext-4C-R_6R!H->C
Total Standard Deviation in ln(k): 53.11312306256407
""",
)

entry(
    index = 56,
    label = "Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-3R-R_4BrCClF->C_Ext-4C-R_N-6R!H->C",
    kinetics = ArrheniusBM(A=(3.87031e-05,'m^3/(mol*s)'), n=3.04873, w0=(572500,'J/mol'), E0=(62012.5,'J/mol'), Tmin=(303.03,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-3R-R_4BrCClF->C_Ext-4C-R_N-6R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-3R-R_4BrCClF->C_Ext-4C-R_N-6R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-3R-R_4BrCClF->C_Ext-4C-R_N-6R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-3R-R_4BrCClF->C_Ext-4C-R_N-6R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 57,
    label = "Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_5R!H->F_Ext-2R-R_6R!H->C_Ext-6C-R_7R!H->C_Ext-7C-R_8R!H->C",
    kinetics = ArrheniusBM(A=(6.58773e-05,'m^3/(mol*s)'), n=2.61294, w0=(572500,'J/mol'), E0=(108780,'J/mol'), Tmin=(303.03,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_5R!H->F_Ext-2R-R_6R!H->C_Ext-6C-R_7R!H->C_Ext-7C-R_8R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_5R!H->F_Ext-2R-R_6R!H->C_Ext-6C-R_7R!H->C_Ext-7C-R_8R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_5R!H->F_Ext-2R-R_6R!H->C_Ext-6C-R_7R!H->C_Ext-7C-R_8R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_5R!H->F_Ext-2R-R_6R!H->C_Ext-6C-R_7R!H->C_Ext-7C-R_8R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 58,
    label = "Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_5R!H->F_Ext-2R-R_6R!H->C_Ext-6C-R_7R!H->C_Ext-7C-R_N-8R!H->C",
    kinetics = ArrheniusBM(A=(4.67866e-05,'m^3/(mol*s)'), n=2.73161, w0=(572500,'J/mol'), E0=(106740,'J/mol'), Tmin=(303.03,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_5R!H->F_Ext-2R-R_6R!H->C_Ext-6C-R_7R!H->C_Ext-7C-R_N-8R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_5R!H->F_Ext-2R-R_6R!H->C_Ext-6C-R_7R!H->C_Ext-7C-R_N-8R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_5R!H->F_Ext-2R-R_6R!H->C_Ext-6C-R_7R!H->C_Ext-7C-R_N-8R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_5R!H->F_Ext-2R-R_6R!H->C_Ext-6C-R_7R!H->C_Ext-7C-R_N-8R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 59,
    label = "Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R_Ext-6R!H-R_5BrCClINOPSSi->O_Ext-6R!H-R_Ext-6R!H-R_Ext-8R!H-R_9R!H->C",
    kinetics = ArrheniusBM(A=(7.68887e-12,'m^3/(mol*s)'), n=4.4869, w0=(572500,'J/mol'), E0=(105480,'J/mol'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=-2.81172952430381e-17, var=0.03505081054756626, Tref=1000.0, N=2, data_mean=0.0, correlation='Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R_Ext-6R!H-R_5BrCClINOPSSi->O_Ext-6R!H-R_Ext-6R!H-R_Ext-8R!H-R_9R!H->C',), comment="""BM rule fitted to 2 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R_Ext-6R!H-R_5BrCClINOPSSi->O_Ext-6R!H-R_Ext-6R!H-R_Ext-8R!H-R_9R!H->C
    Total Standard Deviation in ln(k): 0.3753237286206311"""),
    rank = 11,
    shortDesc = """BM rule fitted to 2 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R_Ext-6R!H-R_5BrCClINOPSSi->O_Ext-6R!H-R_Ext-6R!H-R_Ext-8R!H-R_9R!H->C
Total Standard Deviation in ln(k): 0.3753237286206311""",
    longDesc = 
"""
BM rule fitted to 2 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R_Ext-6R!H-R_5BrCClINOPSSi->O_Ext-6R!H-R_Ext-6R!H-R_Ext-8R!H-R_9R!H->C
Total Standard Deviation in ln(k): 0.3753237286206311
""",
)

entry(
    index = 60,
    label = "Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R_Ext-6R!H-R_5BrCClINOPSSi->O_Ext-6R!H-R_Ext-6R!H-R_Ext-8R!H-R_N-9R!H->C",
    kinetics = ArrheniusBM(A=(2.22077e-12,'m^3/(mol*s)'), n=4.39186, w0=(572500,'J/mol'), E0=(85496.4,'J/mol'), Tmin=(303.03,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R_Ext-6R!H-R_5BrCClINOPSSi->O_Ext-6R!H-R_Ext-6R!H-R_Ext-8R!H-R_N-9R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R_Ext-6R!H-R_5BrCClINOPSSi->O_Ext-6R!H-R_Ext-6R!H-R_Ext-8R!H-R_N-9R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R_Ext-6R!H-R_5BrCClINOPSSi->O_Ext-6R!H-R_Ext-6R!H-R_Ext-8R!H-R_N-9R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R_Ext-6R!H-R_5BrCClINOPSSi->O_Ext-6R!H-R_Ext-6R!H-R_Ext-8R!H-R_N-9R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 61,
    label = "Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R_Ext-6R!H-R_N-5BrCClINOPSSi->O_Ext-7R!H-R",
    kinetics = ArrheniusBM(A=(4.89931e-15,'m^3/(mol*s)'), n=5.3217, w0=(572500,'J/mol'), E0=(118750,'J/mol'), Tmin=(303.03,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R_Ext-6R!H-R_N-5BrCClINOPSSi->O_Ext-7R!H-R',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R_Ext-6R!H-R_N-5BrCClINOPSSi->O_Ext-7R!H-R
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R_Ext-6R!H-R_N-5BrCClINOPSSi->O_Ext-7R!H-R
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R_Ext-6R!H-R_N-5BrCClINOPSSi->O_Ext-7R!H-R
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 62,
    label = "Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-3R-R_4BrCClF->C_Ext-4C-R_6R!H->C_Ext-6C-R_7R!H->C",
    kinetics = ArrheniusBM(A=(0.00828286,'m^3/(mol*s)'), n=2.34779, w0=(572500,'J/mol'), E0=(87031.3,'J/mol'), Tmin=(303.03,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-3R-R_4BrCClF->C_Ext-4C-R_6R!H->C_Ext-6C-R_7R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-3R-R_4BrCClF->C_Ext-4C-R_6R!H->C_Ext-6C-R_7R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-3R-R_4BrCClF->C_Ext-4C-R_6R!H->C_Ext-6C-R_7R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-3R-R_4BrCClF->C_Ext-4C-R_6R!H->C_Ext-6C-R_7R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 63,
    label = "Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-3R-R_4BrCClF->C_Ext-4C-R_6R!H->C_Ext-6C-R_N-7R!H->C",
    kinetics = ArrheniusBM(A=(2.67496e-08,'m^3/(mol*s)'), n=3.75287, w0=(572500,'J/mol'), E0=(220685,'J/mol'), Tmin=(303.03,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-3R-R_4BrCClF->C_Ext-4C-R_6R!H->C_Ext-6C-R_N-7R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-3R-R_4BrCClF->C_Ext-4C-R_6R!H->C_Ext-6C-R_N-7R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-3R-R_4BrCClF->C_Ext-4C-R_6R!H->C_Ext-6C-R_N-7R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-3R-R_Sp-4R!H-3R_N-4R!H->O_N-1R->C_Ext-3R-R_4BrCClF->C_Ext-4C-R_6R!H->C_Ext-6C-R_N-7R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 64,
    label = "Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R_Ext-6R!H-R_5BrCClINOPSSi->O_Ext-6R!H-R_Ext-6R!H-R_Ext-8R!H-R_9R!H->C_Ext-9C-R_10R!H->C",
    kinetics = ArrheniusBM(A=(1.69502e-13,'m^3/(mol*s)'), n=4.7063, w0=(572500,'J/mol'), E0=(85246.3,'J/mol'), Tmin=(303.03,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R_Ext-6R!H-R_5BrCClINOPSSi->O_Ext-6R!H-R_Ext-6R!H-R_Ext-8R!H-R_9R!H->C_Ext-9C-R_10R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R_Ext-6R!H-R_5BrCClINOPSSi->O_Ext-6R!H-R_Ext-6R!H-R_Ext-8R!H-R_9R!H->C_Ext-9C-R_10R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R_Ext-6R!H-R_5BrCClINOPSSi->O_Ext-6R!H-R_Ext-6R!H-R_Ext-8R!H-R_9R!H->C_Ext-9C-R_10R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R_Ext-6R!H-R_5BrCClINOPSSi->O_Ext-6R!H-R_Ext-6R!H-R_Ext-8R!H-R_9R!H->C_Ext-9C-R_10R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

entry(
    index = 65,
    label = "Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R_Ext-6R!H-R_5BrCClINOPSSi->O_Ext-6R!H-R_Ext-6R!H-R_Ext-8R!H-R_9R!H->C_Ext-9C-R_N-10R!H->C",
    kinetics = ArrheniusBM(A=(2.27066e-13,'m^3/(mol*s)'), n=4.68033, w0=(572500,'J/mol'), E0=(85070.3,'J/mol'), Tmin=(303.03,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R_Ext-6R!H-R_5BrCClINOPSSi->O_Ext-6R!H-R_Ext-6R!H-R_Ext-8R!H-R_9R!H->C_Ext-9C-R_N-10R!H->C',), comment="""BM rule fitted to 1 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R_Ext-6R!H-R_5BrCClINOPSSi->O_Ext-6R!H-R_Ext-6R!H-R_Ext-8R!H-R_9R!H->C_Ext-9C-R_N-10R!H->C
    Total Standard Deviation in ln(k): 11.540182761524994"""),
    rank = 11,
    shortDesc = """BM rule fitted to 1 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R_Ext-6R!H-R_5BrCClINOPSSi->O_Ext-6R!H-R_Ext-6R!H-R_Ext-8R!H-R_9R!H->C_Ext-9C-R_N-10R!H->C
Total Standard Deviation in ln(k): 11.540182761524994""",
    longDesc = 
"""
BM rule fitted to 1 training reactions at node Root_Ext-3R-R_N-Sp-4R!H-3R_Ext-2R-R_N-5R!H->F_Ext-2R-R_Ext-6R!H-R_5BrCClINOPSSi->O_Ext-6R!H-R_Ext-6R!H-R_Ext-8R!H-R_9R!H->C_Ext-9C-R_N-10R!H->C
Total Standard Deviation in ln(k): 11.540182761524994
""",
)

