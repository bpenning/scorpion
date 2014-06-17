from math import log10

def get_gluino_x_section_from_slha_file(slhafile):
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
                        if int(words[0])==1000021:
                            return get_x_section_from_mass(float(words[1]))
                    except ValueError:
                        pass

def get_x_section_from_mass(mg):
    for i, mg_xsec in enumerate(data):
        m, xsec = mg_xsec
        if m > mg:
            m1=data[i-1][0]
            m2=m
            log_xsec_1=log10(data[i-1][1])
            log_xsec_2=log10(xsec)
            x=(mg-m1)/(m2-m1)
            log_xsec=x*log_xsec_2 + (1-x)*log_xsec_1
            return (10**log_xsec)*1e-12


data=[
[200 ,     1010.1 ],
[205 ,     888.038 ],
[210 ,     782.82 ],
[215 ,     690.67 ],
[220 ,     610.177 ],
[225 ,     540.437 ],
[230 ,     479.688 ],
[235 ,     426.694 ],
[240 ,     379.598 ],
[245 ,     338.228 ],
[250 ,     302.3 ],
[255 ,     270.569 ],
[260 ,     242.843 ],
[265 ,     217.628 ],
[270 ,     195.926 ],
[275 ,     176.213 ],
[280 ,     158.621 ],
[285 ,     143.51 ],
[290 ,     129.421 ],
[295 ,     117.331 ],
[300 ,     106.269 ],
[305 ,     96.4393 ],
[310 ,     87.6653 ],
[315 ,     79.7744 ],
[320 ,     72.687 ],
[325 ,     66.2653 ],
[330 ,     60.5336 ],
[335 ,     55.369 ],
[340 ,     50.6463 ],
[345 ,     46.4307 ],
[350 ,     42.5663 ],
[355 ,     39.0611 ],
[360 ,     35.9491 ],
[365 ,     33.0444 ],
[370 ,     30.435 ],
[375 ,     28.0354 ],
[380 ,     25.88 ],
[385 ,     23.8782 ],
[390 ,     22.073 ],
[395 ,     20.4283 ],
[400 ,     18.8721 ],
[405 ,     17.4227 ],
[410 ,     16.1707 ],
[415 ,     14.9804 ],
[420 ,     13.8804 ],
[425 ,     12.8817 ],
[430 ,     11.9789 ],
[435 ,     11.0855 ],
[440 ,     10.3321 ],
[445 ,     9.59091 ],
[450 ,     8.93257 ],
[455 ,     8.32507 ],
[460 ,     7.75814 ],
[465 ,     7.24166 ],
[470 ,     6.76031 ],
[475 ,     6.31553 ],
[480 ,     5.89981 ],
[485 ,     5.51473 ],
[490 ,     5.16376 ],
[495 ,     4.82913 ],
[500 ,     4.52485 ],
[505 ,     4.24081 ],
[510 ,     3.96713 ],
[515 ,     3.71672 ],
[520 ,     3.49333 ],
[525 ,     3.27008 ],
[530 ,     3.0669 ],
[535 ,     2.88611 ],
[540 ,     2.70371 ],
[545 ,     2.54097 ],
[550 ,     2.38869 ],
[555 ,     2.24732 ],
[560 ,     2.11612 ],
[565 ,     1.99418 ],
[570 ,     1.87218 ],
[575 ,     1.76075 ],
[580 ,     1.65934 ],
[585 ,     1.56783 ],
[590 ,     1.47661 ],
[595 ,     1.39575 ],
[600 ,     1.31429 ],
[605 ,     1.24324 ],
[610 ,     1.17213 ],
[615 ,     1.10066 ],
[620 ,     1.04006 ],
[625 ,     0.985407 ],
[630 ,     0.931497 ],
[635 ,     0.880476 ],
[640 ,     0.832847 ],
[645 ,     0.786926 ],
[650 ,     0.744472 ],
[655 ,     0.704916 ],
[660 ,     0.667058 ],
[665 ,     0.631306 ],
[670 ,     0.597843 ],
[675 ,     0.56624 ],
[680 ,     0.536593 ],
[685 ,     0.50942 ],
[690 ,     0.482843 ],
[695 ,     0.457432 ],
[700 ,     0.433971 ],
[705 ,     0.411416 ],
[710 ,     0.390086 ],
[715 ,     0.370645 ],
[720 ,     0.352514 ],
[725 ,     0.33408 ],
[730 ,     0.317737 ],
[735 ,     0.301462 ],
[740 ,     0.286089 ],
[745 ,     0.271852 ],
[750 ,     0.258598 ],
[755 ,     0.24633 ],
[760 ,     0.234109 ],
[765 ,     0.222826 ],
[770 ,     0.211617 ],
[775 ,     0.201354 ],
[780 ,     0.191111 ],
[785 ,     0.18202 ],
[790 ,     0.173801 ],
[795 ,     0.165647 ],
[800 ,     0.157399 ],
[805 ,     0.150172 ],
[810 ,     0.142838 ],
[815 ,     0.135881 ],
[820 ,     0.129292 ],
[825 ,     0.123451 ],
[830 ,     0.117155 ],
[835 ,     0.11198 ],
[840 ,     0.106724 ],
[845 ,     0.101468 ],
[850 ,     0.0966803 ],
[855 ,     0.0922606 ],
[860 ,     0.0878851 ],
[865 ,     0.0837802 ],
[870 ,     0.0799422 ],
[875 ,     0.0762312 ],
[880 ,     0.072704 ],
[885 ,     0.0693676 ],
[890 ,     0.0661977 ],
[895 ,     0.0631751 ],
[900 ,     0.060276 ],
[905 ,     0.0575238 ],
[910 ,     0.054922 ],
[915 ,     0.052471 ],
[920 ,     0.0500914 ],
[925 ,     0.0478794 ],
[930 ,     0.0457437 ],
[935 ,     0.0436802 ],
[940 ,     0.0417248 ],
[945 ,     0.0398759 ],
[950 ,     0.0381246 ],
[955 ,     0.0364316 ],
[960 ,     0.0348221 ],
[965 ,     0.033282 ],
[970 ,     0.0318672 ],
[975 ,     0.0304177 ],
[980 ,     0.0290889 ],
[985 ,     0.0278154 ],
[990 ,     0.0265993 ],
[995 ,     0.0254762 ],
[1000 ,    0.0243547 ],
[1005 ,    0.0233412 ],
[1010 ,    0.0222816 ],
[1015 ,    0.0213092 ],
[1020 ,    0.0204008 ],
[1025 ,    0.0195439 ],
[1030 ,    0.0186399 ],
[1035 ,    0.0178964 ],
[1040 ,    0.0170845 ],
[1045 ,    0.0163849 ],
[1050 ,    0.0156931 ],
[1055 ,    0.0149871 ],
[1060 ,    0.0143428 ],
[1065 ,    0.0137535 ],
[1070 ,    0.0131578 ],
[1075 ,    0.0126249 ],
[1080 ,    0.0120916 ],
[1085 ,    0.0115882 ],
[1090 ,    0.0110906 ],
[1095 ,    0.0105972 ],
[1100 ,    0.0101744 ],
[1105 ,    0.00975322 ],
[1110 ,    0.00934364 ],
[1115 ,    0.00895424 ],
[1120 ,    0.00858597 ],
[1125 ,    0.00823026 ],
[1130 ,    0.00788948 ],
[1135 ,    0.00756671 ],
[1140 ,    0.00724917 ],
[1145 ,    0.00695169 ],
[1150 ,    0.00666673 ],
[1155 ,    0.00639182 ],
[1160 ,    0.00613149 ],
[1165 ,    0.00588262 ],
[1170 ,    0.0056425 ],
[1175 ,    0.0054104 ],
[1180 ,    0.00519563 ],
[1185 ,    0.00498468 ],
[1190 ,    0.00477732 ],
[1195 ,    0.00458673 ],
[1200 ,    0.00440078 ],
[1205 ,    0.00422593 ],
[1210 ,    0.0040514 ],
[1215 ,    0.00389138 ],
[1220 ,    0.00372785 ],
[1225 ,    0.00357858 ],
[1230 ,    0.00343744 ],
[1235 ,    0.00330131 ],
[1240 ,    0.0031691 ],
[1245 ,    0.00303793 ],
[1250 ,    0.00291565 ],
[1255 ,    0.00280625 ],
[1260 ,    0.00269022 ],
[1265 ,    0.00258061 ],
[1270 ,    0.00248068 ],
[1275 ,    0.00238091 ],
[1280 ,    0.00228827 ],
[1285 ,    0.00219864 ],
[1290 ,    0.00210573 ],
[1295 ,    0.00202657 ],
[1300 ,    0.00194443 ],
[1305 ,    0.00186692 ],
[1310 ,    0.00179418 ],
[1315 ,    0.00172271 ],
[1320 ,    0.00165064 ],
[1325 ,    0.00158945 ],
[1330 ,    0.00152847 ],
[1335 ,    0.0014675 ],
[1340 ,    0.0014106 ],
[1345 ,    0.0013492 ],
[1350 ,    0.00129951 ],
[1355 ,    0.00124801 ],
[1360 ,    0.0011984 ],
[1365 ,    0.00114871 ],
[1370 ,    0.00110711 ],
[1375 ,    0.00106538 ],
[1380 ,    0.00102485 ],
[1385 ,    0.00098222 ],
[1390 ,    0.000943754 ],
[1395 ,    0.000906445 ],
[1400 ,    0.000871201 ],
[1405 ,    0.00083687 ],
[1410 ,    0.000804842 ],
[1415 ,    0.000773344 ],
[1420 ,    0.000743009 ],
[1425 ,    0.000714083 ],
[1430 ,    0.000686632 ],
[1435 ,    0.000659782 ],
[1440 ,    0.000634212 ],
[1445 ,    0.000609136 ],
[1450 ,    0.00058599 ],
[1455 ,    0.00056309 ],
[1460 ,    0.000541725 ],
[1465 ,    0.000520546 ],
[1470 ,    0.000500264 ],
[1475 ,    0.000481062 ],
[1480 ,    0.000462071 ],
[1485 ,    0.000443875 ],
[1490 ,    0.000426659 ],
[1495 ,    0.000410737 ],
[1500 ,    0.000394612 ],
[1505 ,    0.000379584 ],
[1510 ,    0.000364688 ],
[1515 ,    0.000350346 ],
[1520 ,    0.000337197 ],
[1525 ,    0.000323984 ],
[1530 ,    0.000311214 ],
[1535 ,    0.000299031 ],
[1540 ,    0.000287604 ],
[1545 ,    0.00027663 ],
[1550 ,    0.000265533 ],
[1555 ,    0.000255837 ],
[1560 ,    0.000245873 ],
[1565 ,    0.000236606 ],
[1570 ,    0.000227404 ],
[1575 ,    0.000218539 ],
[1580 ,    0.000209696 ],
[1585 ,    0.00020156 ],
[1590 ,    0.000194171 ],
[1595 ,    0.000186477 ],
[1600 ,    0.000179423 ],
[1605 ,    0.000172361 ],
[1610 ,    0.000165332 ],
[1615 ,    0.000159471 ],
[1620 ,    0.000153392 ],
[1625 ,    0.000147446 ],
[1630 ,    0.000141507 ],
[1635 ,    0.000136462 ],
[1640 ,    0.000131372 ],
[1645 ,    0.000126312 ],
[1650 ,    0.000121361 ],
[1655 ,    0.000116488 ],
[1660 ,    0.000112313 ],
[1665 ,    0.000107505 ],
[1670 ,    0.000103502 ],
[1675 ,    9.95545e-05 ],
[1680 ,    9.56648e-05 ],
[1685 ,    9.19785e-05 ],
[1690 ,    8.84717e-05 ],
[1695 ,    8.5105e-05 ],
[1700 ,    8.17761e-05 ],
[1705 ,    7.86327e-05 ],
[1710 ,    7.55929e-05 ],
[1715 ,    7.26295e-05 ],
[1720 ,    6.98621e-05 ],
[1725 ,    6.71419e-05 ],
[1730 ,    6.4556e-05 ],
[1735 ,    6.20851e-05 ],
[1740 ,    5.97211e-05 ],
[1745 ,    5.74255e-05 ],
[1750 ,    5.51497e-05 ],
[1755 ,    5.30377e-05 ],
[1760 ,    5.10074e-05 ],
[1765 ,    4.90166e-05 ],
[1770 ,    4.71201e-05 ],
[1775 ,    4.52786e-05 ],
[1780 ,    4.35033e-05 ],
[1785 ,    4.18723e-05 ],
[1790 ,    4.01971e-05 ],
[1795 ,    3.86643e-05 ],
[1800 ,    3.71771e-05 ],
[1805 ,    3.57622e-05 ],
[1810 ,    3.43583e-05 ],
[1815 ,    3.29627e-05 ],
[1820 ,    3.1738e-05 ],
[1825 ,    3.04796e-05 ],
[1830 ,    2.93288e-05 ],
[1835 ,    2.81348e-05 ],
[1840 ,    2.70614e-05 ],
[1845 ,    2.60099e-05 ],
[1850 ,    2.50085e-05 ],
[1855 ,    2.40021e-05 ],
[1860 ,    2.31231e-05 ],
[1865 ,    2.21893e-05 ],
[1870 ,    2.12977e-05 ],
[1875 ,    2.04823e-05 ],
[1880 ,    1.97133e-05 ],
[1885 ,    1.89005e-05 ],
[1890 ,    1.81887e-05 ],
[1895 ,    1.74813e-05 ],
[1900 ,    1.67835e-05 ],
[1905 ,    1.60881e-05 ],
[1910 ,    1.54787e-05 ],
[1915 ,    1.4876e-05 ],
[1920 ,    1.42691e-05 ],
[1925 ,    1.37589e-05 ],
[1930 ,    1.32437e-05 ],
[1935 ,    1.2654e-05 ],
[1940 ,    1.21604e-05 ],
[1945 ,    1.17397e-05 ],
[1950 ,    1.12426e-05 ],
[1955 ,    1.08304e-05 ],
[1960 ,    1.04217e-05 ],
[1965 ,    9.97664e-06 ],
[1970 ,    9.58718e-06 ],
[1975 ,    9.20837e-06 ],
[1980 ,    8.84663e-06 ],
[1985 ,    8.49858e-06 ],
[1990 ,    8.15342e-06 ],
[1995 ,    7.84163e-06 ],
[2000 ,    7.52864e-06 ],
]