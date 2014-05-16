from math import log10

def get_x_section_from_slha_file(slhafile):
    with open(slhafile,'r') as f:
        is_in_block_mass=False
        for line in f:
            words=line.split()
            if len(words)>=2:
                if words[0].upper()=='BLOCK' and words[1].upper()=='MASS':
                    is_in_block_mass=True
                elif words[0].upper=='BLOCK' and not words[1].upper=='MASS':
                    is_in_block_mass=False
                if is_in_block_mass:
                    try:
                        if int(words[0])==1000006:
                            return get_x_section_from_mass(float(words[1]))
                    except ValueError:
                        pass

def get_x_section_from_mass(mstop):
    for i, mstop_xsec in enumerate(data):
        m, xsec = mstop_xsec
        if m > mstop:
            m1=data[i-1][0]
            m2=m
            log_xsec_1=log10(data[i-1][1])
            log_xsec_2=log10(xsec)
            x=(mstop-m1)/(m2-m1)
            log_xsec=x*log_xsec_2 + (1-x)*log_xsec_1
            return (10**log_xsec)*1e-12

#data from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SUSYCrossSections8TeVstopsbottom
data=[
[100.0  ,    559.757 ],
[105.0  ,    448.456 ],
[110.0  ,    361.917 ],
[115.0  ,    293.281 ],
[120.0  ,    240.077 ],
[125.0  ,    197.122 ],
[130.0  ,    163.376 ],
[135.0  ,    135.791 ],
[140.0  ,    113.319 ],
[145.0  ,    95.0292 ],
[150.0  ,    80.268 ],
[155.0  ,    68.0456 ],
[160.0  ,    58.01 ],
[165.0  ,    49.6639 ],
[170.0  ,    42.6441 ],
[175.0  ,    36.7994 ],
[180.0  ,    31.8695 ],
[185.0  ,    27.7028 ],
[190.0  ,    24.1585 ],
[195.0  ,    21.1597 ],
[200.0  ,    18.5245 ],
[205.0  ,    16.2439 ],
[210.0  ,    14.3201 ],
[215.0  ,    12.6497 ],
[220.0  ,    11.1808 ],
[225.0  ,    9.90959 ],
[230.0  ,    8.78125 ],
[235.0  ,    7.81646 ],
[240.0  ,    6.96892 ],
[245.0  ,    6.22701 ],
[250.0  ,    5.57596 ],
[255.0  ,    5.00108 ],
[260.0  ,    4.48773 ],
[265.0  ,    4.03416 ],
[270.0  ,    3.63085 ],
[275.0  ,    3.2781 ],
[280.0  ,    2.95613 ],
[285.0  ,    2.67442 ],
[290.0  ,    2.42299 ],
[295.0  ,    2.19684 ],
[300.0  ,    1.99608 ],
[305.0  ,    1.81486 ],
[310.0  ,    1.64956 ],
[315.0  ,    1.50385 ],
[320.0  ,    1.3733 ],
[325.0  ,    1.25277 ],
[330.0  ,    1.14277 ],
[335.0  ,    1.04713 ],
[340.0  ,    0.959617 ],
[345.0  ,    0.879793 ],
[350.0  ,    0.807323 ],
[355.0  ,    0.74141 ],
[360.0  ,    0.681346 ],
[365.0  ,    0.626913 ],
[370.0  ,    0.576882 ],
[375.0  ,    0.531443 ],
[380.0  ,    0.489973 ],
[385.0  ,    0.452072 ],
[390.0  ,    0.4176 ],
[395.0  ,    0.385775 ],
[400.0  ,    0.35683 ],
[405.0  ,    0.329881 ],
[410.0  ,    0.305512 ],
[415.0  ,    0.283519 ],
[420.0  ,    0.262683 ],
[425.0  ,    0.243755 ],
[430.0  ,    0.226367 ],
[435.0  ,    0.209966 ],
[440.0  ,    0.195812 ],
[445.0  ,    0.181783 ],
[450.0  ,    0.169668 ],
[455.0  ,    0.158567 ],
[460.0  ,    0.147492 ],
[465.0  ,    0.137392 ],
[470.0  ,    0.128326 ],
[475.0  ,    0.119275 ],
[480.0  ,    0.112241 ],
[485.0  ,    0.104155 ],
[490.0  ,    0.0977878 ],
[495.0  ,    0.091451 ],
[500.0  ,    0.0855847 ],
[505.0  ,    0.0801322 ],
[510.0  ,    0.0751004 ],
[515.0  ,    0.0703432 ],
[520.0  ,    0.0660189 ],
[525.0  ,    0.0618641 ],
[530.0  ,    0.0580348 ],
[535.0  ,    0.0545113 ],
[540.0  ,    0.0511747 ],
[545.0  ,    0.0481537 ],
[550.0  ,    0.0452067 ],
[555.0  ,    0.0424781 ],
[560.0  ,    0.0399591 ],
[565.0  ,    0.0376398 ],
[570.0  ,    0.0354242 ],
[575.0  ,    0.0333988 ],
[580.0  ,    0.0313654 ],
[585.0  ,    0.0295471 ],
[590.0  ,    0.0279395 ],
[595.0  ,    0.0263263 ],
[600.0  ,    0.0248009 ],
[605.0  ,    0.0233806 ],
[610.0  ,    0.0220672 ],
[615.0  ,    0.0208461 ],
[620.0  ,    0.0196331 ],
[625.0  ,    0.0185257 ],
[630.0  ,    0.0175075 ],
[635.0  ,    0.0164955 ],
[640.0  ,    0.0155809 ],
[645.0  ,    0.0147721 ],
[650.0  ,    0.0139566 ],
[655.0  ,    0.0132456 ],
[660.0  ,    0.0125393 ],
[665.0  ,    0.0118287 ],
[670.0  ,    0.0112223 ],
[675.0  ,    0.0106123 ],
[680.0  ,    0.0100516 ],
[685.0  ,    0.0095256 ],
[690.0  ,    0.0090306 ],
[695.0  ,    0.00856339 ],
[700.0  ,    0.0081141 ],
[705.0  ,    0.00769525 ],
[710.0  ,    0.00730084 ],
[715.0  ,    0.00692243 ],
[720.0  ,    0.00656729 ],
[725.0  ,    0.00623244 ],
[730.0  ,    0.00591771 ],
[735.0  ,    0.00561049 ],
[740.0  ,    0.00532605 ],
[745.0  ,    0.00506044 ],
[750.0  ,    0.00480639 ],
[755.0  ,    0.00455979 ],
[760.0  ,    0.00433688 ],
[765.0  ,    0.00412174 ],
[770.0  ,    0.00391839 ],
[775.0  ,    0.00372717 ],
[780.0  ,    0.00354211 ],
[785.0  ,    0.00336904 ],
[790.0  ,    0.00320476 ],
[795.0  ,    0.00304935 ],
[800.0  ,    0.00289588 ],
[805.0  ,    0.00275424 ],
[810.0  ,    0.0026184 ],
[815.0  ,    0.00249291 ],
[820.0  ,    0.00237168 ],
[825.0  ,    0.00226163 ],
[830.0  ,    0.00214607 ],
[835.0  ,    0.00204589 ],
[840.0  ,    0.00195172 ],
[845.0  ,    0.0018573 ],
[850.0  ,    0.00176742 ],
[855.0  ,    0.00168383 ],
[860.0  ,    0.00160403 ],
[865.0  ,    0.00153063 ],
[870.0  ,    0.00145772 ],
[875.0  ,    0.0013878 ],
[880.0  ,    0.00132077 ],
[885.0  ,    0.00126234 ],
[890.0  ,    0.00120568 ],
[895.0  ,    0.00114627 ],
[900.0  ,    0.00109501 ],
[905.0  ,    0.001044 ],
[910.0  ,    0.000996193 ],
[915.0  ,    0.00095071 ],
[920.0  ,    0.000907494 ],
[925.0  ,    0.000866391 ],
[930.0  ,    0.000826533 ],
[935.0  ,    0.000789573 ],
[940.0  ,    0.000753768 ],
[945.0  ,    0.000719675 ],
[950.0  ,    0.000687022 ],
[955.0  ,    0.000656279 ],
[960.0  ,    0.000626876 ],
[965.0  ,    0.000598955 ],
[970.0  ,    0.000571551 ],
[975.0  ,    0.000546728 ],
[980.0  ,    0.000522495 ],
[985.0  ,    0.000499017 ],
[990.0  ,    0.000476255 ],
[995.0  ,    0.000455959 ],
[1000.0 ,    0.000435488 ],
[1005.0 ,    0.000416116 ],
[1010.0 ,    0.00039791 ],
[1015.0 ,    0.000379994 ],
[1020.0 ,    0.000363934 ],
[1025.0 ,    0.000347646 ],
[1030.0 ,    0.00033204 ],
[1035.0 ,    0.000318049 ],
[1040.0 ,    0.000303756 ],
[1045.0 ,    0.000290392 ],
[1050.0 ,    0.000277943 ],
[1055.0 ,    0.000265929 ],
[1060.0 ,    0.000254659 ],
[1065.0 ,    0.000243251 ],
[1070.0 ,    0.00023289 ],
[1075.0 ,    0.000222651 ],
[1080.0 ,    0.000213396 ],
[1085.0 ,    0.000204211 ],
[1090.0 ,    0.000196038 ],
[1095.0 ,    0.000187913 ],
[1100.0 ,    0.000179699 ],
[1105.0 ,    0.000172125 ],
[1110.0 ,    0.000165045 ],
[1115.0 ,    0.000157905 ],
[1120.0 ,    0.000151236 ],
[1125.0 ,    0.000144737 ],
[1130.0 ,    0.000138657 ],
[1135.0 ,    0.000133343 ],
[1140.0 ,    0.000127478 ],
[1145.0 ,    0.00012234 ],
[1150.0 ,    0.000117215 ],
[1155.0 ,    0.000112199 ],
[1160.0 ,    0.000107256 ],
[1165.0 ,    0.000103046 ],
[1170.0 ,    9.86633e-05 ],
[1175.0 ,    9.44977e-05 ],
[1180.0 ,    9.05131e-05 ],
[1185.0 ,    8.67972e-05 ],
[1190.0 ,    8.31669e-05 ],
[1195.0 ,    7.96568e-05 ],
[1200.0 ,    7.63052e-05 ],
[1205.0 ,    7.31846e-05 ],
[1210.0 ,    7.00987e-05 ],
[1215.0 ,    6.72666e-05 ],
[1220.0 ,    6.44711e-05 ],
[1225.0 ,    6.18495e-05 ],
[1230.0 ,    5.93775e-05 ],
[1235.0 ,    5.69258e-05 ],
[1240.0 ,    5.46007e-05 ],
[1245.0 ,    5.23902e-05 ],
[1250.0 ,    5.02385e-05 ],
[1255.0 ,    4.81984e-05 ],
[1260.0 ,    4.6268e-05 ],
[1265.0 ,    4.43332e-05 ],
[1270.0 ,    4.25134e-05 ],
[1275.0 ,    4.07881e-05 ],
[1280.0 ,    3.91417e-05 ],
[1285.0 ,    3.75567e-05 ],
[1290.0 ,    3.60147e-05 ],
[1295.0 ,    3.4617e-05 ],
[1300.0 ,    3.31856e-05 ],
[1305.0 ,    3.18643e-05 ],
[1310.0 ,    3.05669e-05 ],
[1315.0 ,    2.93562e-05 ],
[1320.0 ,    2.81652e-05 ],
[1325.0 ,    2.70203e-05 ],
[1330.0 ,    2.59296e-05 ],
[1335.0 ,    2.49015e-05 ],
[1340.0 ,    2.39221e-05 ],
[1345.0 ,    2.29308e-05 ],
[1350.0 ,    2.20166e-05 ],
[1355.0 ,    2.11331e-05 ],
[1360.0 ,    2.02909e-05 ],
[1365.0 ,    1.94752e-05 ],
[1370.0 ,    1.8714e-05 ],
[1375.0 ,    1.79681e-05 ],
[1380.0 ,    1.72697e-05 ],
[1385.0 ,    1.65664e-05 ],
[1390.0 ,    1.58669e-05 ],
[1395.0 ,    1.52759e-05 ],
[1400.0 ,    1.46673e-05 ],
[1405.0 ,    1.40698e-05 ],
[1410.0 ,    1.35484e-05 ],
[1415.0 ,    1.30386e-05 ],
[1420.0 ,    1.253e-05 ],
[1425.0 ,    1.20335e-05 ],
[1430.0 ,    1.15351e-05 ],
[1435.0 ,    1.11127e-05 ],
[1440.0 ,    1.06301e-05 ],
[1445.0 ,    1.02226e-05 ],
[1450.0 ,    9.84287e-06 ],
[1455.0 ,    9.45889e-06 ],
[1460.0 ,    9.09124e-06 ],
[1465.0 ,    8.73208e-06 ],
[1470.0 ,    8.38911e-06 ],
[1475.0 ,    8.05993e-06 ],
[1480.0 ,    7.74413e-06 ],
[1485.0 ,    7.44538e-06 ],
[1490.0 ,    7.14898e-06 ],
[1495.0 ,    6.86798e-06 ],
[1500.0 ,    6.59698e-06 ],
[1505.0 ,    6.33638e-06 ],
[1510.0 ,    6.08509e-06 ],
[1515.0 ,    5.84265e-06 ],
[1520.0 ,    5.6093e-06 ],
[1525.0 ,    5.38651e-06 ],
[1530.0 ,    5.1742e-06 ],
[1535.0 ,    4.97273e-06 ],
[1540.0 ,    4.77117e-06 ],
[1545.0 ,    4.58178e-06 ],
[1550.0 ,    4.39875e-06 ],
[1555.0 ,    4.22724e-06 ],
[1560.0 ,    4.05582e-06 ],
[1565.0 ,    3.89113e-06 ],
[1570.0 ,    3.72963e-06 ],
[1575.0 ,    3.58102e-06 ],
[1580.0 ,    3.43895e-06 ],
[1585.0 ,    3.29773e-06 ],
[1590.0 ,    3.16675e-06 ],
[1595.0 ,    3.03589e-06 ],
[1600.0 ,    2.91221e-06 ],
[1605.0 ,    2.79219e-06 ],
[1610.0 ,    2.66816e-06 ],
[1615.0 ,    2.55964e-06 ],
[1620.0 ,    2.44758e-06 ],
[1625.0 ,    2.34644e-06 ],
[1630.0 ,    2.24557e-06 ],
[1635.0 ,    2.15574e-06 ],
[1640.0 ,    2.0669e-06 ],
[1645.0 ,    1.98515e-06 ],
[1650.0 ,    1.90249e-06 ],
[1655.0 ,    1.83115e-06 ],
[1660.0 ,    1.76963e-06 ],
[1665.0 ,    1.70036e-06 ],
[1670.0 ,    1.63964e-06 ],
[1675.0 ,    1.58714e-06 ],
[1680.0 ,    1.52765e-06 ],
[1685.0 ,    1.47645e-06 ],
[1690.0 ,    1.42533e-06 ],
[1695.0 ,    1.37529e-06 ],
[1700.0 ,    1.32467e-06 ],
[1705.0 ,    1.27389e-06 ],
[1710.0 ,    1.22407e-06 ],
[1715.0 ,    1.18167e-06 ],
[1720.0 ,    1.13218e-06 ],
[1725.0 ,    1.09086e-06 ],
[1730.0 ,    1.04182e-06 ],
[1735.0 ,    1.00428e-06 ],
[1740.0 ,    9.64964e-07 ],
[1745.0 ,    9.26859e-07 ],
[1750.0 ,    8.90605e-07 ],
[1755.0 ,    8.56319e-07 ],
[1760.0 ,    8.2246e-07 ],
[1765.0 ,    7.90884e-07 ],
[1770.0 ,    7.59995e-07 ],
[1775.0 ,    7.3087e-07 ],
[1780.0 ,    7.02001e-07 ],
[1785.0 ,    6.75402e-07 ],
[1790.0 ,    6.49168e-07 ],
[1795.0 ,    6.24339e-07 ],
[1800.0 ,    6.00052e-07 ],
[1805.0 ,    5.76531e-07 ],
[1810.0 ,    5.5474e-07 ],
[1815.0 ,    5.33316e-07 ],
[1820.0 ,    5.129e-07 ],
[1825.0 ,    4.93338e-07 ],
[1830.0 ,    4.74104e-07 ],
[1835.0 ,    4.55409e-07 ],
[1840.0 ,    4.37889e-07 ],
[1845.0 ,    4.21436e-07 ],
[1850.0 ,    4.04973e-07 ],
[1855.0 ,    3.89723e-07 ],
[1860.0 ,    3.7399e-07 ],
[1865.0 ,    3.59977e-07 ],
[1870.0 ,    3.45452e-07 ],
[1875.0 ,    3.32155e-07 ],
[1880.0 ,    3.19667e-07 ],
[1885.0 ,    3.06274e-07 ],
[1890.0 ,    2.95073e-07 ],
[1895.0 ,    2.83743e-07 ],
[1900.0 ,    2.724e-07 ],
[1905.0 ,    2.62039e-07 ],
[1910.0 ,    2.51734e-07 ],
[1915.0 ,    2.4149e-07 ],
[1920.0 ,    2.32237e-07 ],
[1925.0 ,    2.23897e-07 ],
[1930.0 ,    2.14605e-07 ],
[1935.0 ,    2.06292e-07 ],
[1940.0 ,    1.99019e-07 ],
[1945.0 ,    1.90756e-07 ],
[1950.0 ,    1.83446e-07 ],
[1955.0 ,    1.7615e-07 ],
[1960.0 ,    1.69813e-07 ],
[1965.0 ,    1.62682e-07 ],
[1970.0 ,    1.56389e-07 ],
[1975.0 ,    1.51029e-07 ],
[1980.0 ,    1.44807e-07 ],
[1985.0 ,    1.39563e-07 ],
[1990.0 ,    1.34259e-07 ],
[1995.0 ,    1.29026e-07 ],
[2000.0 ,    1.23774e-07 ],
]

