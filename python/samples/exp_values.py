import libjad_DelphesAnalysis as j
from python.core.types import *

bg_at = DoubleVector([787.0, 310.0, 202.0, 60.4, 20.3, 7.7, 3.2, 2.8])
bgunc_at = DoubleVector([0.05, 0.05, 0.05, 0.07, 0.09, 0.10, 0.13, 0.15])
data_at = IntVector([782, 321, 196, 62, 21, 6, 3, 1])

bg_atlas5 = DoubleVector([7.0, 64.8, 115.0, 5.39, 5.68, 38.6, 214.0, 6.84, 12.1, 34.0, 119.0])
bgunc_atlas5 = DoubleVector([0.35, 0.19, 0.19, 0.42, 0.32, 0.13, 0.06, 0.34, 0.26, 0.17, 0.10])
data_atlas5 = IntVector([1, 59, 85, 1, 14, 36, 210, 9, 13, 25, 148])

bg_at7 = DoubleVector([3735.5, 1478.9, 1049.6, 370.1, 116.4, 52.6, 14.0, 8.6])
bgunc_at7 = DoubleVector([0.017, 0.03, 0.03, 0.04, 0.062, 0.09, 0.21, 0.22])
#data_at7 = IntVector([3735, 1478, 1049, 370, 154, 52, 14, 8])
data_at7 = IntVector([3703, 1536, 1043, 346, 122, 44, 14, 6])

bg_at7b = DoubleVector([2933.0, 1139.0, 783.0, 261.0, 81.5, 34.2, 10.4, 5.3,
                        630.0, 271.0, 202.0, 78.0, 24.2, 10.6, 2.9, 2.2,
                        162.0, 61.8, 58.8, 28.0, 9.0, 7.1, 0.6, 0.9,
                        10.5, 7.1, 5.8, 3.1, 1.7, 0.7, 0.1, 0.2])
bgunc_at7b = DoubleVector([0.02, 0.035, 0.034, 0.054, 0.080, 0.117, 0.269, 0.321,
                           0.041, 0.059, 0.050, 0.088, 0.120, 0.160, 0.310, 0.318,
                           0.080, 0.102, 0.082, 0.125, 0.156, 0.197, 0.500, 0.444,
                           0.333, 0.310, 0.241, 0.322, 0.294, 0.714, 1.00, 0.500])
data_at7b = IntVector([2919, 1166, 769, 255, 91, 31, 10, 4,
                       614, 294, 214, 71, 20, 6, 4, 0,
                       160, 68, 52, 19, 11, 7, 0, 2,
                       10, 8, 8, 1, 0, 0, 0, 0])

bg_at8 = DoubleVector([4136.5, 1765.0, 1310.5, 431.4, 154.8, 60.2, 25.2, 21.2])
bgunc_at8 = DoubleVector([0.018, 0.02, 0.02, 0.04, 0.06, 0.09, 0.16, 0.18])
#data_at8 = IntVector([4136, 1765, 1310, 431, 154, 60, 25, 21])
data_at8 = IntVector([4093, 1803, 1269, 462, 158, 64, 14, 31])

bg_at8b = DoubleVector([6235.0, 2900.0, 1955.0, 558.0, 186.0, 51.3, 21.2, 16.1,
                        1010.0, 447.0, 390.0, 250.0, 111.0, 53.3, 18.5, 19.4,
                        1162.0, 481.0, 341.0, 86.7, 24.8, 7.2, 3.3, 2.1,
                        521.0, 232.0, 188.0, 106.0, 42.1, 17.9, 9.8, 6.8,
                        224.0, 98.2, 59.0, 12.8, 3.0, 0.5, 0.1, 0.1,
                        208.0, 103.0, 85.9, 51.7, 19.9, 6.8, 1.7, 1.3,
                        25.3, 11.7, 6.7, 3.9, 2.3, 1.2, 0.3, 0.1,
                        0.9, 0.3, 0.6])
bgunc_at8b = DoubleVector([0.016, 0.021, 0.020, 0.027, 0.059, 0.074, 0.108, 0.106,
                           0.034, 0.043, 0.049, 0.048, 0.081, 0.081, 0.130, 0.139,
                           0.032, 0.040, 0.047, 0.065, 0.113, 0.153, 0.212, 0.238,
                           0.048, 0.065, 0.064, 0.057, 0.105, 0.123, 0.153, 0.176,
                           0.067, 0.086, 0.102, 0.125, 0.300, 0.400, 1.000, 1.000,
                           0.082, 0.087, 0.084, 0.091, 0.171, 0.191, 0.412, 0.308,
                           0.198, 0.154, 0.209, 0.205, 0.261, 0.333, 0.667, 1.000,
                           0.778, 0.667, 0.500])
data_at8b = IntVector([6232, 2904, 1965, 552, 177, 58, 16, 25,
                       1009, 452, 375, 274, 113, 56, 16, 27,
                       1164, 473, 329, 95, 23, 8, 4, 1,
                       515, 236, 204, 92, 51, 13, 13, 6,
                       222, 107, 58, 12, 5, 1, 0, 0,
                       204, 107, 84, 59, 24, 5, 1, 2,
                       25, 13, 4, 2, 2, 3, 0, 0,
                       1, 0, 2])

bg_alphat20b = DoubleVector([12391.0, 5566.0, 3424.0, 2473.0, 684.0, 231.0, 69.8, 27.9, 11.0, 5.4, 3.7,
                             106.0, 520.0, 448.0, 369.0, 224.0, 109.0, 46.8, 19.2, 9.2, 3.7, 3.2,
			     1626.0, 862.0, 548.0, 413.0, 103.0, 26.6, 9.5, 3.1, 2.6, 0.4, 0.2,
			     39.8, 223.0, 211.0, 171.0, 98.7, 35.5, 14.8, 6.8, 3.2, 1.1, 1.2,
			     177.0, 113.0, 97.3, 66.0, 14.6, 2.9, 1.1, 0.2, 0.1,
			     13.0, 75.5, 96.1, 69.5, 43.9, 16.1, 5.2, 1.7, 1.9,
			     1.0, 8.1, 11.7, 8.0, 5.3, 2.3, 0.8, 0.2, 0.2,
			     0.0, 0.1, 0.6, 0.7])
			     
bgunc_alphat20b = DoubleVector([0.030, 0.058, 0.060, 0.053, 0.056, 0.082, 0.095, 0.140, 0.191, 0.185, 0.243,
                                0.094, 0.065, 0.121, 0.084, 0.076, 0.128, 0.130, 0.146, 0.185, 0.297, 0.219,
				0.037, 0.056, 0.064, 0.068, 0.078, 0.120, 0.168, 0.290, 0.346, 0.500, 0.500,
				0.088, 0.063, 0.118, 0.105, 0.102, 0.169, 0.149, 0.221, 0.281, 0.455, 0.333,
				0.051, 0.062, 0.068, 0.095, 0.116, 0.172, 0.273, 0.500, 1.000,
                                0.077, 0.064, 0.119, 0.119, 0.130, 0.230, 0.231, 0.235, 0.263,
				0.200, 0.099, 0.162, 0.175, 0.226, 0.261, 0.250, 0.500, 0.500,
				0.000, 1.000, 0.500, 0.429])

bg_alphat20b_valid = DoubleVector([12391.0, 5566.0, 3424.0, 2473.0, 684.0, 231.0, 69.8, 27.9, 11.0, 5.4, 3.7,
                             106.0, 520.0, 448.0, 369.0, 224.0, 109.0, 46.8, 19.2, 9.2, 3.7, 3.2,
			     1626.0, 862.0, 548.0, 413.0, 103.0, 26.6, 9.5, 3.1, 2.6, 0.4, 0.2,
			     39.8, 223.0, 211.0, 171.0, 98.7, 35.5, 14.8, 6.8, 3.2, 1.1, 1.2
			     ]) 
bgunc_alphat20b_valid = DoubleVector([0.030, 0.058, 0.060, 0.053, 0.056, 0.082, 0.095, 0.140, 0.191, 0.185, 0.243,
                                0.094, 0.065, 0.121, 0.084, 0.076, 0.128, 0.130, 0.146, 0.185, 0.297, 0.219,
				0.037, 0.056, 0.064, 0.068, 0.078, 0.120, 0.168, 0.290, 0.346, 0.500, 0.500,
				0.088, 0.063, 0.118, 0.105, 0.102, 0.169, 0.149, 0.221, 0.281, 0.455, 0.333
				])
data_alphat20b_valid = IntVector([13310, 5404, 3580, 2475, 698, 224, 81, 28, 12, 9, 3,
                            103, 590, 455, 391, 274, 149, 56, 18, 12, 7, 8,
			    1802, 835, 597, 416, 97, 29, 7, 4, 1, 0, 0,
			    38, 207, 250, 181, 108, 41, 10, 14, 6, 1, 1
			    ])
# For expected limits:
#data_alphat20b = IntVector([12391, 5566, 3424, 2473, 684, 231, 70, 28, 11, 5, 4,
#                             106, 520, 448, 369, 224, 109, 47, 19, 9, 4, 3,
#			     1626, 862, 548, 413, 103, 27, 10, 3, 3, 0, 0,
#			     40, 223, 211, 171, 99, 35, 15, 7, 3, 1, 1,
#			     177, 113, 97, 66, 15, 3, 1, 0, 0,
#			     13, 75, 96, 69, 44, 16, 5, 2, 2,
#			     1, 8, 12, 8, 5, 2, 1, 0, 0,
#			     0, 0, 1, 1])

# For observed limits:
data_alphat20b = IntVector([13310, 5404, 3580, 2475, 698, 224, 81, 28, 12, 9, 3,
                            103, 590, 455, 391, 274, 149, 56, 18, 12, 7, 8,
			    1802, 835, 597, 416, 97, 29, 7, 4, 1, 0, 0,
			    38, 207, 250, 181, 108, 41, 10, 14, 6, 1, 1,
			    174, 119, 117, 70, 18, 7, 1, 0, 0,
			    17, 82, 95, 86, 48, 19, 10, 2, 1,
			    0, 9, 8, 8, 9, 1, 3, 0, 0,
			    0, 0, 0, 3])
	
#bg_alphat20b = DoubleVector([13235.0, 5417.0, 3562.0, 2482.0, 689.0, 231.0, 72.2, 28.1, 11.2, 6.0, 3.7,
                        #     104.0, 567.0, 454.0, 397.0, 250.0, 134.0, 55.5, 19.2, 9.8, 4.6, 4.2,
			#     1743.0, 843.0, 580.0, 413.0, 102.0, 27.0, 9.1, 3.3, 2.3, 0.3, 0.2,
			#     39.4, 216.0, 238.0, 179.0, 104.0, 36.4, 14.2, 8.6, 3.9, 1.1, 1.2,
			#     176.0, 115.0, 104.0, 68.1, 15.2, 3.3, 1.1, 0.2, 0.1,
			#     13.3, 77.1, 95.5, 77.4, 48.2, 18.4, 6.2, 1.7, 1.8,
			#     1.0, 8.2, 10.9, 8.2, 5.8, 2.2, 0.9, 0.2, 0.2,
			#     0.0, 0.1, 0.5, 0.9])
			     
#bgunc_alphat20b = DoubleVector([0.009, 0.013, 0.019, 0.016, 0.022, 0.056, 0.072, 0.114, 0.143, 0.183, 0.243,
                        #        0.077, 0.037, 0.044, 0.035, 0.044, 0.075, 0.079, 0.130, 0.163, 0.217, 0.238,
			#	0.018, 0.031, 0.040, 0.039, 0.049, 0.107, 0.132, 0.212, 0.304, 0.667, 0.500,
			#	0.076, 0.042, 0.059, 0.056, 0.058, 0.096, 0.127, 0.174, 0.231, 0.364, 0.333,
			#	0.045, 0.052, 0.058, 0.069, 0.105, 0.182, 0.273, 0.500, 0.000,
                        #       0.075, 0.057, 0.081, 0.076, 0.077, 0.168, 0.161, 0.235, 0.222,
			#	0.200, 0.098, 0.147, 0.134, 0.172, 0.273, 0.333, 0.500, 0.500,
			#	0.000, 1.000, 0.600, 0.333])


		    

bg_lp5 = DoubleVector([84.1, 21.3, 8.2, 33.3, 7.5, 2.9, 15.7, 5.9, 4.6, 65.3, 11.2, 4.8, 19.4, 6.3, 2.6, 9.2, 6.8, 2.6])
bgunc_lp5 = DoubleVector([0.10, 0.15, 0.22, 0.12, 0.24, 0.43, 0.19, 0.27, 0.34, 0.11, 0.24, 0.39, 0.19, 0.30, 0.58, 0.28, 0.33, 0.61])
data_lp5 = IntVector([78, 22, 8, 23, 8, 1, 16, 7, 2, 71, 13, 8, 29, 5, 1, 11, 5, 1])

bg_monojet5 = DoubleVector([7842.0, 2757.0, 1225.0, 573.0])
bgunc_monojet5 = DoubleVector([0.05, 0.06, 0.08, 0.11])
data_monojet5 = IntVector([7584, 2774, 1142, 522])

#bg_multi8 = DoubleVector([6023.1, 2280.0, 419.5, 57.6])
#bgunc_multi8 = DoubleVector([0.05, 0.06, 0.08, 0.11])
#data_multi8 = IntVector([7584, 2774, 1142, 522])
bg_monojet20 = DoubleVector([49154.0, 18506.0, 7875.0, 3663.0,1931.0,949.0,501.0])
bgunc_monojet20 = DoubleVector([0.034, 0.037,0.043,0.054,0.068,0.087,0.12])
data_monojet20 = IntVector([50419, 19108, 8056, 3677,1772,894,508])


bg_DMbSR1 = DoubleVector([385])
bgunc_DMbSR1 = DoubleVector([35])
data_DMbSR1 = IntVector([440])

bg_monojet20ss = DoubleVector([35862.0, 17409.0, 8064.0, 3907.0,2098.0,1096.0,563.0])
bgunc_monojet20ss = DoubleVector([0.0411, 0.046,0.054,0.064,0.076,0.097,0.126])
data_monojet20ss = IntVector([36582, 17646, 8119, 3896, 1898, 1003, 565])

bg_zerolepmt2_8_20_gluino_test = DoubleVector([
      #CMS PAS SUS-13-019 Table 2
      #low HT; >=6 jets; 1 b-jet
      32.0, 14.7, 4.8, 
      #low HT; >=6 jets; 2 b-jet
      12.0, 4.6, 2.8, 
      #low HT; >=6 jets; >=3 b-jet
      16.1, 4.6,
      #medium HT; >=6 jets; 1 b-jet
      38.7, 21.1, 10.5, 3.0,
      #medium HT; >=6 jets; 2 b-jet
      41.0,19.4, 10.4, 4.3, 
      #medium HT; >=6 jets; >=3 b-jet
      31.9, 16.1, 6.1,
      #high HT; >=6 jets; 1 b-jet
      7.3, 5.1, 2.3, 
      #high HT; >=6 jets; 2 b-jet
      10.6, 4.7, 
      #high HT; >=6 jets; >=3 b-jet
      4.5
      ])

data_zerolepmt2_8_20_gluino_test = IntVector([
      #low HT; >=6 jets; 1 b-jet
      31, 23, 11,  
      #low HT; >=6 jets; 2 b-jet
      15, 13,6, 
      #low HT; >=6 jets; >=3 b-jet
      16,  7, 
      #medium HT; >=6 jets; 1 b-jet
      38, 21, 13, 4,
      #medium HT; >=6 jets; 2 b-jet
      54,28,8,6,
      #medium HT; >=6 jets; >=3 b-jet
      17,13, 1,
      #high HT; >=6 jets; 1 b-jet
      6,5,1,
      #high HT; >=6 jets; 2 b-jet
      10,2,
      #high HT; >=6 jets; >=3 b-jet
      3
      ])

bgunc_zerolepmt2_8_20_gluino_test = DoubleVector([
#low HT; >=6 jets; 1 b-jet
# 6.7,3.1,1.5,
#low HT; >=6 jets; 2 b-jet
# 4.3,1.6,1.0,
#low HT; >=6 jets; >=3 b-jet
# 6.2,1.7,
#medium HT; >=6 jets; 1 b-jet
# 5.9, 3.5, 1.9, 0.8, 
#medium HT; >=6 jets; 2 b-jet
# 7.0, 3.8, 2.1, 0.8, 
#medium HT; >=6 jets; >=3 b-jet
# 11.4, 6.3, 2.4,
#high HT; >=6 jets; 1 b-jet
# 3.2,2.4,1.1,
#high HT; >=6 jets; 2 b-jet
# 6.0,2.9,
#high HT; >=6 jets; >=3 b-jet                                        #
# 2.1
  0.209375   , 0.21088435 , 0.3125     , 0.35833333 , 0.34782609 , 0.35714286  ,
  0.38509317 , 0.36956522 , 0.15245478 , 0.16587678 , 0.18095238 , 0.26666667  ,
  0.17073171 , 0.19587629 , 0.20192308 , 0.18604651 , 0.35736677 , 0.39130435  ,
  0.39344262 , 0.43835616 , 0.47058824 , 0.47826087 , 0.56603774 , 0.61702128  ,
  0.46666667
      ])
bg_zerolepmt2_8_20_stop_test = DoubleVector([
      #CMS PAS SUS-13-019 Table 2
      #low HT; 3-5 jets; 1 b-jet
      305, 167, 103, 43.6, 17.9, 4.0, 
      #low HT; 3-5 jets; 2 b-jet
      91.1, 52.7, 18.6, 4.5, 
      #low HT; >=6 jets; 1 b-jet
      32.0, 14.7, 4.8, 
      #low HT; >=6 jets; 2 b-jet
      12.0, 4.6, 2.8, 
      #low HT; >=6 jets; >=3 b-jet
      16.1, 4.6,
      #medium HT; 3-5 jets; 1 b-jet
      93.4, 69.5, 52.8, 38.6, 15.9, 3.6, 
      #medium HT; 3-5 jets; 2 b-jet
      42.4, 26.5, 15.4, 5.5, 2.9, 
      #medium HT; >=6 jets; 1 b-jet
      38.7, 21.1, 10.5, 3.0,
      #medium HT; >=6 jets; 2 b-jet
      41.0,19.4, 10.4, 4.3, 
      #medium HT; >=6 jets; >=3 b-jet
      31.9, 16.1, 6.1,
      #high HT; 3-5 jets; 1 b-jet
      13.5, 8.7, 6.2, 3.5, 
      #high HT; 3-5 jets; 2 b-jet
      6.8, 2.9, 
      #high HT; >=6 jets; 1 b-jet
      7.3, 5.1, 2.3, 
      #high HT; >=6 jets; 2 b-jet
      10.6, 4.7, 
      #high HT; >=6 jets; >=3 b-jet
      4.5
      ])

data_zerolepmt2_8_20_stop_test = IntVector([
      #low HT; 3-5 jets; 1 b-jet
      300, 172, 98, 47, 19, 4, 
      #low HT; 3-5 jets; 2 b-jet
      97, 39, 16, 11, 
      #low HT; >=6 jets; 1 b-jet
      31, 23, 11,  
      #low HT; >=6 jets; 2 b-jet
      15, 13,6, 
      #low HT; >=6 jets; >=3 b-jet
      16,  7, 
      #medium HT; 3-5 jets; 1 b-jet
      87, 71, 63, 47, 19, 4, 
      #medium HT; 3-5 jets; 2 b-jet
      53, 29, 19,  11, 5, 
      #medium HT; >=6 jets; 1 b-jet
      38, 21, 13, 4,
      #medium HT; >=6 jets; 2 b-jet
      54,28,8,6,
      #medium HT; >=6 jets; >=3 b-jet
      17,13, 1,
      #high HT; 3-5 jets; 1 b-jet
      28,7,9,3,
      #high HT; 3-5 jets; 2 b-jet
      9,6,
      #high HT; >=6 jets; 1 b-jet
      6,5,1,
      #high HT; >=6 jets; 2 b-jet
      10,2,
      #high HT; >=6 jets; >=3 b-jet
      3
      ])

bgunc_zerolepmt2_8_20_stop_test = DoubleVector([
#low HT; 3-5 jets; 1 b-jet
# 34,21,16,8.7,4.1,1.1,
#low HT; 3-5 jets; 2 b-jet
# 22.0,13.7,5.8,1.9,
#low HT; >=6 jets; 1 b-jet
# 6.7,3.1,1.5,
#low HT; >=6 jets; 2 b-jet
# 4.3,1.6,1.0,
#low HT; >=6 jets; >=3 b-jet
# 6.2,1.7,
#medium HT; 3-5 jets; 1 b-jet
# 10.7, 8.7, 6.8, 5.1, 3.2, 0.9, 
#medium HT; 3-5 jets; 2 b-jet
# 7.5, 5.5, 3.7, 1.7, 1.1, 
#medium HT; >=6 jets; 1 b-jet
# 5.9, 3.5, 1.9, 0.8, 
#medium HT; >=6 jets; 2 b-jet
# 7.0, 3.8, 2.1, 0.8, 
#medium HT; >=6 jets; >=3 b-jet
# 11.4, 6.3, 2.4,
#high HT; 3-5 jets; 1 b-jet
# 3.1,2.2,1.6,1.0,
#high HT; 3-5 jets; 2 b-jet
# 2.3,1.1,
#high HT; >=6 jets; 1 b-jet
# 3.2,2.4,1.1,
#high HT; >=6 jets; 2 b-jet
# 6.0,2.9,
#high HT; >=6 jets; >=3 b-jet                                        #
# 2.1
  0.11147541 , 0.1257485   ,0.15533981  ,0.19954128  ,0.22905028  ,0.275      ,
  0.24149286 , 0.25996205  ,0.31182796  ,0.42222222  ,0.209375    ,0.21088435  ,
  0.3125     , 0.35833333  ,0.34782609  ,0.35714286  ,0.38509317  ,0.36956522  ,
  0.11456103 , 0.12517986  ,0.12878788  ,0.13212435  ,0.20125786  ,0.25        ,
  0.17688679 , 0.20754717  ,0.24025974  ,0.30909091  ,0.37931034  ,0.15245478  ,
  0.16587678 , 0.18095238  ,0.26666667  ,0.17073171  ,0.19587629  ,0.20192308  ,
  0.18604651 , 0.35736677  ,0.39130435  ,0.39344262  ,0.22962963  ,0.25287356  ,
  0.25806452 , 0.28571429  ,0.33823529  ,0.37931034  ,0.43835616  ,0.47058824  ,
  0.47826087 , 0.56603774  ,0.61702128  ,0.46666667              
      ])

bg_zerolepmt2_8_20_t2qq = DoubleVector([
      #CMS PAS SUS-13-019 Table 2
      #low HT; 2 jets; 0 b-jet
      553, 395, 288, 236, 165, 68.9, 17.3, 4.1, 
      #low HT; 2 jets; >= 1 b-jet
      56.4, 34.2, 25.9, 19.9, 12.6, 2.6, 
      #low HT; 3-5 jets; 0 b-jet
      979, 711, 492, 280, 138, 60.0, 13.8, 3.6,
      #low HT; 3-5 jets; 1 b-jet
      305, 167, 103, 43.6, 17.9, 4.0, 
      #low HT; >=6 jets; 0 b-jet
      50.8, 14.7, 7.3, 
      #medium HT; 2 jets; 0 b-jet
      167, 128, 85.8, 70.0, 38.1, 43.4, 21.3, 20.8, 3.5, 
      #medium HT; 2 jets; >= 1 b-jet
      27.4, 21.1, 13.4, 7.3, 3.4, 
      #medium HT; 3-5 jets; 0 b-jet
      243, 180, 134, 112, 89.0, 67.0, 35.0, 10.0, 3.4, 
      #medium HT; 3-5 jets; 1 b-jet
      93.4, 69.5, 52.8, 38.6, 15.9, 3.6, 
      #medium HT; >=6 jets; 0 b-jet
      38.5, 19.3, 14.1, 5.8, 2.3, 
      #high HT; 2 jets; 0 b-jet
      21.9, 19.4, 14.5, 6.3, 4.3, 3.0, 
      #high HT; 2 jets; >= 1 b-jet
      11.4, 4.4, 
      #high HT; 3-5 jets; 0 b-jet
      34.9, 31.1, 25.5, 19.3, 9.1, 5.0, 4.4, 
      #high HT; 3-5 jets; 1 b-jet
      13.5, 8.7, 6.2, 3.5, 
      #high HT; >=6 jets; 0 b-jet
      12.1, 10.1, 4.5, 
      ])

data_zerolepmt2_8_20_t2qq = IntVector([
      #low HT; 2 jets; 0 b-jet
      588, 451, 318, 232, 162, 61, 19, 1, 
      #low HT; 2 jets; >= 1 b-jet
      56, 44, 29, 13, 15, 3, 
      #low HT; 3-5 jets; 0 b-jet
      1041, 827, 522, 333, 145, 66, 21, 2, 
      #low HT; 3-5 jets; 1 b-jet
      300, 172, 98, 47, 19, 4, 
      #low HT; >=6 jets; 0 b-jet
      56, 16, 8,  
      #medium HT; 2 jets; 0 b-jet
      171, 104, 91, 78,  48, 45, 29, 10, 2, 
      #medium HT; 2 jets; >= 1 b-jet
      30, 19, 15, 7, 9,  
      #medium HT; 3-5 jets; 0 b-jet
      234, 203, 152, 119, 91,  75, 40, 16, 4, 
      #medium HT; 3-5 jets; 1 b-jet
      87, 71, 63, 47, 19, 4, 
      #medium HT; >=6 jets; 0 b-jet
      44, 34,  23, 9, 4, 
      #high HT; 2 jets; 0 b-jet
      18, 18, 10, 9, 8, 6, 
      #high HT; 2 jets; >= 1 b-jet
      2, 2, 
      #high HT; 3-5 jets; 0 b-jet
      39,32,25,19,6,5,5,
      #high HT; 3-5 jets; 1 b-jet
      28,7,9,3,
      #high HT; >=6 jets; 0 b-jet
      12,7,2,
      ])

bgunc_zerolepmt2_8_20_t2qq = DoubleVector(
#low HT; 2 jets; 0 b-jet
# 70,53,40,52,36,15.5,4.3,1.6,
#low HT; 2 jets; >= 1 b-jet
# 12.8,8.1,7.4,5.8,3.8,0.8,
#low HT; 3-5 jets; 0 b-jet
# 108,86,65,57,29,13.6,3.9,1.5,
#low HT; 3-5 jets; 1 b-jet
# 34,21,16,8.7,4.1,1.1,
#low HT; >=6 jets; 0 b-jet
# 8.9,3.1,2.3,
#medium HT; 2 jets; 0 b-jet
# 21, 17, 11.3, 10.3, 5.8, 10.1, 4.7, 5.6, 1.4, 
#medium HT; 2 jets; >= 1 b-jet
# 9.6, 7.5, 5.4, 3.5, 1.7, 
#medium HT; 3-5 jets; 0 b-jet
# 23, 19, 16, 14, 12.2, 14.2, 8.0, 2.7, 1.5, 
#medium HT; 3-5 jets; 1 b-jet
# 10.7, 8.7, 6.8, 5.1, 3.2, 0.9, 
#medium HT; >=6 jets; 0 b-jet
# 6.2, 3.6, 2.8, 1.9, 0.8, 
#high HT; 2 jets; 0 b-jet
# 4.9,4.3,3.4,1.8,1.6,1.4,
#high HT; 2 jets; >= 1 b-jet
# 8.1,2.6,
#high HT; 3-5 jets; 0 b-jet
# 4.7,4.7,4.3,3.5,2.5,1.6,1.6,
#high HT; 3-5 jets; 1 b-jet
# 3.1,2.2,1.6,1.0,
#high HT; >=6 jets; 0 b-jet
# 2.9,3.2,1.7,
[ 0.12658228 , 0.13417722 , 0.13888889  ,0.22033898  ,0.21818182 , 0.22496372,
  0.24855491 , 0.3902439  , 0.22695035  ,0.23684211  ,0.28571429 , 0.29145729,
  0.3015873  , 0.30769231 , 0.11031665  ,0.1209564   ,0.13211382 , 0.20357143,
  0.21014493 , 0.22666667 , 0.2826087   ,0.41666667  ,0.11147541 , 0.1257485 ,
  0.15533981 , 0.19954128 , 0.22905028  ,0.275       ,0.17519685 , 0.21088435,
  0.31506849 , 0.1257485  , 0.1328125   ,0.13170163  ,0.14714286 , 0.15223097,
  0.23271889 , 0.22065728 , 0.26923077  ,0.4         ,0.35036496 , 0.35545024,
  0.40298507 , 0.47945205 , 0.5         ,0.09465021  ,0.10555556 , 0.11940299,
  0.125      , 0.13707865 , 0.2119403   ,0.22857143  ,0.27       , 0.44117647,
  0.11456103 , 0.12517986 , 0.12878788  ,0.13212435  ,0.20125786 , 0.25      ,
  0.16103896 , 0.1865285  , 0.19858156  ,0.32758621  ,0.34782609 , 0.22374429,
  0.22164948 , 0.23448276 , 0.28571429  ,0.37209302  ,0.46666667 , 0.71052632,
  0.59090909 , 0.13467049 , 0.1511254   ,0.16862745  ,0.18134715 , 0.27472527,
  0.32       , 0.36363636 , 0.22962963  ,0.25287356  ,0.25806452 , 0.28571429,
  0.23966942 , 0.31683168 , 0.37777778]                        
      )

bg_zerolepmt2_8_20 = DoubleVector([
      #CMS PAS SUS-13-019 Table 2
      #low HT; 2 jets; 0 b-jet
      553, 395, 288, 236, 165, 68.9, 17.3, 4.1, 
      #low HT; 2 jets; >= 1 b-jet
      56.4, 34.2, 25.9, 19.9, 12.6, 2.6, 
      #low HT; 3-5 jets; 0 b-jet
      979, 711, 492, 280, 138, 60.0, 13.8, 3.6,
      #low HT; 3-5 jets; 1 b-jet
      305, 167, 103, 43.6, 17.9, 4.0, 
      #low HT; 3-5 jets; 2 b-jet
      91.1, 52.7, 18.6, 4.5, 
      #low HT; >=6 jets; 0 b-jet
      50.8, 14.7, 7.3, 
      #low HT; >=6 jets; 1 b-jet
      32.0, 14.7, 4.8, 
      #low HT; >=6 jets; 2 b-jet
      12.0, 4.6, 2.8, 
      #low HT; >=6 jets; >=3 b-jet
      16.1, 4.6,
      #medium HT; 2 jets; 0 b-jet
      167, 128, 85.8, 70.0, 38.1, 43.4, 21.3, 20.8, 3.5, 
      #medium HT; 2 jets; >= 1 b-jet
      27.4, 21.1, 13.4, 7.3, 3.4, 
      #medium HT; 3-5 jets; 0 b-jet
      243, 180, 134, 112, 89.0, 67.0, 35.0, 10.0, 3.4, 
      #medium HT; 3-5 jets; 1 b-jet
      93.4, 69.5, 52.8, 38.6, 15.9, 3.6, 
      #medium HT; 3-5 jets; 2 b-jet
      42.4, 26.5, 15.4, 5.5, 2.9, 
      #medium HT; >=6 jets; 0 b-jet
      38.5, 19.3, 14.1, 5.8, 2.3, 
      #medium HT; >=6 jets; 1 b-jet
      38.7, 21.1, 10.5, 3.0,
      #medium HT; >=6 jets; 2 b-jet
      41.0,19.4, 10.4, 4.3, 
      #medium HT; >=6 jets; >=3 b-jet
      31.9, 16.1, 6.1,
      #high HT; 2 jets; 0 b-jet
      21.9, 19.4, 14.5, 6.3, 4.3, 3.0, 
      #high HT; 2 jets; >= 1 b-jet
      11.4, 4.4, 
      #high HT; 3-5 jets; 0 b-jet
      34.9, 31.1, 25.5, 19.3, 9.1, 5.0, 4.4, 
      #high HT; 3-5 jets; 1 b-jet
      13.5, 8.7, 6.2, 3.5, 
      #high HT; 3-5 jets; 2 b-jet
      6.8, 2.9, 
      #high HT; >=6 jets; 0 b-jet
      12.1, 10.1, 4.5, 
      #high HT; >=6 jets; 1 b-jet
      7.3, 5.1, 2.3, 
      #high HT; >=6 jets; 2 b-jet
      10.6, 4.7, 
      #high HT; >=6 jets; >=3 b-jet
      4.5
      ])

data_zerolepmt2_8_20 = IntVector([
      #low HT; 2 jets; 0 b-jet
      588, 451, 318, 232, 162, 61, 19, 1, 
      #low HT; 2 jets; >= 1 b-jet
      56, 44, 29, 13, 15, 3, 
      #low HT; 3-5 jets; 0 b-jet
      1041, 827, 522, 333, 145, 66, 21, 2, 
      #low HT; 3-5 jets; 1 b-jet
      300, 172, 98, 47, 19, 4, 
      #low HT; 3-5 jets; 2 b-jet
      97, 39, 16, 11, 
      #low HT; >=6 jets; 0 b-jet
      56, 16, 8,  
      #low HT; >=6 jets; 1 b-jet
      31, 23, 11,  
      #low HT; >=6 jets; 2 b-jet
      15, 13,6, 
      #low HT; >=6 jets; >=3 b-jet
      16,  7, 
      #medium HT; 2 jets; 0 b-jet
      171, 104, 91, 78,  48, 45, 29, 10, 2, 
      #medium HT; 2 jets; >= 1 b-jet
      30, 19, 15, 7, 9,  
      #medium HT; 3-5 jets; 0 b-jet
      234, 203, 152, 119, 91,  75, 40, 16, 4, 
      #medium HT; 3-5 jets; 1 b-jet
      87, 71, 63, 47, 19, 4, 
      #medium HT; 3-5 jets; 2 b-jet
      53, 29, 19,  11, 5, 
      #medium HT; >=6 jets; 0 b-jet
      44, 34,  23, 9, 4, 
      #medium HT; >=6 jets; 1 b-jet
      38, 21, 13, 4,
      #medium HT; >=6 jets; 2 b-jet
      54,28,8,6,
      #medium HT; >=6 jets; >=3 b-jet
      17,13, 1,
      #high HT; 2 jets; 0 b-jet
      18, 18, 10, 9, 8, 6, 
      #high HT; 2 jets; >= 1 b-jet
      2, 2, 
      #high HT; 3-5 jets; 0 b-jet
      39,32,25,19,6,5,5,
      #high HT; 3-5 jets; 1 b-jet
      28,7,9,3,
      #high HT; 3-5 jets; 2 b-jet
      9,6,
      #high HT; >=6 jets; 0 b-jet
      12,7,2,
      #high HT; >=6 jets; 1 b-jet
      6,5,1,
      #high HT; >=6 jets; 2 b-jet
      10,2,
      #high HT; >=6 jets; >=3 b-jet
      3
      ])

bgunc_zerolepmt2_8_20 = DoubleVector([
#low HT; 2 jets; 0 b-jet
# 70,53,40,52,36,15.5,4.3,1.6,
#low HT; 2 jets; >= 1 b-jet
# 12.8,8.1,7.4,5.8,3.8,0.8,
#low HT; 3-5 jets; 0 b-jet
# 108,86,65,57,29,13.6,3.9,1.5,
#low HT; 3-5 jets; 1 b-jet
# 34,21,16,8.7,4.1,1.1,
#low HT; 3-5 jets; 2 b-jet
# 22.0,13.7,5.8,1.9,
#low HT; >=6 jets; 0 b-jet
# 8.9,3.1,2.3,
#low HT; >=6 jets; 1 b-jet
# 6.7,3.1,1.5,
#low HT; >=6 jets; 2 b-jet
# 4.3,1.6,1.0,
#low HT; >=6 jets; >=3 b-jet
# 6.2,1.7,
#medium HT; 2 jets; 0 b-jet
# 21, 17, 11.3, 10.3, 5.8, 10.1, 4.7, 5.6, 1.4, 
#medium HT; 2 jets; >= 1 b-jet
# 9.6, 7.5, 5.4, 3.5, 1.7, 
#medium HT; 3-5 jets; 0 b-jet
# 23, 19, 16, 14, 12.2, 14.2, 8.0, 2.7, 1.5, 
#medium HT; 3-5 jets; 1 b-jet
# 10.7, 8.7, 6.8, 5.1, 3.2, 0.9, 
#medium HT; 3-5 jets; 2 b-jet
# 7.5, 5.5, 3.7, 1.7, 1.1, 
#medium HT; >=6 jets; 0 b-jet
# 6.2, 3.6, 2.8, 1.9, 0.8, 
#medium HT; >=6 jets; 1 b-jet
# 5.9, 3.5, 1.9, 0.8, 
#medium HT; >=6 jets; 2 b-jet
# 7.0, 3.8, 2.1, 0.8, 
#medium HT; >=6 jets; >=3 b-jet
# 11.4, 6.3, 2.4,
#high HT; 2 jets; 0 b-jet
# 4.9,4.3,3.4,1.8,1.6,1.4,
#high HT; 2 jets; >= 1 b-jet
# 8.1,2.6,
#high HT; 3-5 jets; 0 b-jet
# 4.7,4.7,4.3,3.5,2.5,1.6,1.6,
#high HT; 3-5 jets; 1 b-jet
# 3.1,2.2,1.6,1.0,
#high HT; 3-5 jets; 2 b-jet
# 2.3,1.1,
#high HT; >=6 jets; 0 b-jet
# 2.9,3.2,1.7,
#high HT; >=6 jets; 1 b-jet
# 3.2,2.4,1.1,
#high HT; >=6 jets; 2 b-jet
# 6.0,2.9,
#high HT; >=6 jets; >=3 b-jet                                        #
# 2.1
        0.12658228,  0.13417722,  0.13888889,  0.22033898,  0.21818182,
        0.22496372,  0.24855491,  0.3902439 ,  0.22695035,  0.23684211,
        0.28571429,  0.29145729,  0.3015873 ,  0.30769231,  0.11031665,
        0.1209564 ,  0.13211382,  0.20357143,  0.21014493,  0.22666667,
        0.2826087 ,  0.41666667,  0.11147541,  0.1257485 ,  0.15533981,
        0.19954128,  0.22905028,  0.275     ,  0.24149286,  0.25996205,
        0.31182796,  0.42222222,  0.17519685,  0.21088435,  0.31506849,
        0.209375  ,  0.21088435,  0.3125    ,  0.35833333,  0.34782609,
        0.35714286,  0.38509317,  0.36956522,  0.1257485 ,  0.1328125 ,
        0.13170163,  0.14714286,  0.15223097,  0.23271889,  0.22065728,
        0.26923077,  0.4       ,  0.35036496,  0.35545024,  0.40298507,
        0.47945205,  0.5       ,  0.09465021,  0.10555556,  0.11940299,
        0.125     ,  0.13707865,  0.2119403 ,  0.22857143,  0.27      ,
        0.44117647,  0.11456103,  0.12517986,  0.12878788,  0.13212435,
        0.20125786,  0.25      ,  0.17688679,  0.20754717,  0.24025974,
        0.30909091,  0.37931034,  0.16103896,  0.1865285 ,  0.19858156,
        0.32758621,  0.34782609,  0.15245478,  0.16587678,  0.18095238,
        0.26666667,  0.17073171,  0.19587629,  0.20192308,  0.18604651,
        0.35736677,  0.39130435,  0.39344262,  0.22374429,  0.22164948,
        0.23448276,  0.28571429,  0.37209302,  0.46666667,  0.71052632,
        0.59090909,  0.13467049,  0.1511254 ,  0.16862745,  0.18134715,
        0.27472527,  0.32      ,  0.36363636,  0.22962963,  0.25287356,
        0.25806452,  0.28571429,  0.33823529,  0.37931034,  0.23966942,
        0.31683168,  0.37777778,  0.43835616,  0.47058824,  0.47826087,
        0.56603774,  0.61702128,  0.46666667
      ])



#tables 4 and 6
bg_lp8_20b_all=DoubleVector([251., 83., 31., 11.5, 29., 17., 9.5, 4.7, 1662., 537., 180., 66., 79., 38., 19., 9.9])
#absolute uncerainties are=[50., 21., 8., 3.6, 7., 5., 2.8, 1.4, 203., 75., 28., 13., 12., 7., 5., 2.7 ]
bgunc_lp8_20b_all=DoubleVector([0.199, 0.253, 0.258, 0.313, 0.241, 0.294, 0.295, 0.298, 0.122, 0.140, 0.156, 0.197, 0.152, 0.184, 0.263, 0.273])
data_lp8_20b_all=IntVector([227, 69, 21, 9, 23, 11, 3, 2, 1624, 487, 151, 52, 90, 39, 18, 5])
bg_lp8_20b_t2tt=DoubleVector([251., 83., 31., 11.5, 29., 17., 9.5, 4.7])
#absolute uncerainties are=[50., 21., 8., 3.6, 7., 5., 2.8, 1.4, 203., 75., 28., 13., 12., 7., 5., 2.7 ]
bgunc_lp8_20b_t2tt=DoubleVector([0.199, 0.253, 0.258, 0.313, 0.241, 0.294, 0.295, 0.298])
data_lp8_20b_t2tt=IntVector([227, 69, 21, 9, 23, 11, 3, 2])
bg_lp8_20b_t2bbww=DoubleVector([ 1662., 537., 180., 66., 79., 38., 19., 9.9])
##absolute uncerainties are=[50., 21., 8., 3.6, 7., 5., 2.8, 1.4, 203., 75., 28., 13., 12., 7., 5., 2.7 ]
bgunc_lp8_20b_t2bbww=DoubleVector([ 0.122, 0.140, 0.156, 0.197, 0.152, 0.184, 0.263, 0.273])
data_lp8_20b_t2bbww=IntVector([ 1624, 487, 151, 52, 90, 39, 18, 5])

bg_os5 = DoubleVector([5.7, 5.2, 5.6, 5.7, 5.2, 5.6])
bgunc_os5 = DoubleVector([1.02, 0.87, 0.71, 1.02, 0.87, 0.71])
data_os5 = IntVector([9, 6, 5, 10, 5, 13])

bg_ss5 = DoubleVector([1.1, 1.2, 2.6])
bgunc_ss5 = DoubleVector([1.0, 1.0, 0.54])
data_ss5 = IntVector([1, 0, 3])

#use this bin for the NS studies as it is optimal
bg_ssb8 = DoubleVector([2.2])
bgunc_ssb8 = DoubleVector([0.45])
data_ssb8 = IntVector([1])

bg_ssb8full = DoubleVector([40.0, 32.0, 2.2, 8.1, 5.7, 1.7, 1.2, 8.1])
bgunc_ssb8full = DoubleVector([0.35, 0.34, 0.45, 0.42, 0.42, 0.41, 0.50, 0.41])
data_ssb8full = IntVector([43, 38, 1, 10, 7, 1, 1, 9])

bg_ss8LowPt = DoubleVector([44.0, 12.0, 12.0, 9.1, 21.0, 13.0, 3.5, 5.8,
                             32.0, 6.0, 17.0, 10.0, 13.0, 5.5, 4.2, 6.8,
                             7.6, 1.5, 7.1, 4.4, 2.8, 1.3, 1.8, 3.4])
bgunc_ss8LowPt = DoubleVector([16.0, 4.0, 5.0, 3.4, 8.0, 5.0, 1.4, 2.1,
                                13.0, 2.2, 7.0, 4.0, 5.0, 2.0, 1.6, 2.5,
                                2.8, 0.7, 2.7, 1.7, 1.1, 0.6, 0.8, 1.3])
data_ss8LowPt = IntVector([50, 17, 13, 4, 22, 18, 2, 4,
                            40, 5, 15, 6, 9, 5, 3, 11,
                            10, 1, 6, 11, 1, 2, 0, 3])

bg_ss8HighPt = DoubleVector([51.0, 9.0, 8.0, 5.6, 20.0, 9.0, 2.4, 3.6,
                             36.0, 3.8, 10.0, 5.9, 11.0, 3.9, 2.8, 4.0,
                             7.1, 1.0, 3.8, 2.8, 2.9, 0.8, 1.2, 2.2])
#bgunc_ss8HighPt = DoubleVector([18.0, 3.5, 3.1, 2.1, 7.0, 4.0, 1.0, 1.5,
#                                14.0, 1.4, 4.0, 2.2, 4.0, 1.5, 1.1, 1.5,
#                                2.5, 0.5, 1.4, 1.2, 1.1, 0.5, 0.6, 1.0])
bgunc_ss8HighPt = DoubleVector([ 
    0.35294118,  0.38888889,  0.3875    ,  0.375     ,  0.35      ,
    0.44444444,  0.41666667,  0.41666667,  0.38888889,  0.36842105,
    0.4       ,  0.37288136,  0.36363636,  0.38461538,  0.39285714,
    0.375     ,  0.35211268,  0.5       ,  0.36842105,  0.42857143,
    0.37931034,  0.625     ,  0.5       ,  0.45454545])
data_ss8HighPt = IntVector([48, 11, 5, 2, 12, 11, 1, 3,
                            29, 5, 6, 2, 11, 2, 3, 7,
                            12, 1, 3, 7, 4, 1, 0, 2])
#only the 2 b-tag signal regions
bg_ss8HighPtBtag2 = DoubleVector([ 7.1, 1.0, 3.8, 2.8, 2.9, 0.8, 1.2, 2.2])
#bgunc_ss8HighPtBtag2 = DoubleVector([18.0, 3.5, 3.1, 2.1, 7.0, 4.0, 1.0, 1.5,
#                                14.0, 1.4, 4.0, 2.2, 4.0, 1.5, 1.1, 1.5,
#                                2.5, 0.5, 1.4, 1.2, 1.1, 0.5, 0.6, 1.0])
bgunc_ss8HighPtBtag2 = DoubleVector([ 0.35211268,  0.5       ,  0.36842105,  0.42857143,
    0.37931034,  0.625     ,  0.5       ,  0.45454545])
data_ss8HighPtBtag2 = IntVector([ 12, 1, 3, 7, 4, 1, 0, 2])

bg_cms3l8=DoubleVector(
        [1.00000000e-02,   1.00000000e-02,   2.00000000e-02,
         1.10000000e-01,   4.00000000e-03,   1.00000000e-02,
         4.00000000e-03,   1.20000000e-01,   4.00000000e-03,
         7.00000000e-02,   4.00000000e-03,   2.00000000e-02,
         1.00000000e-02,   2.50000000e-01,   1.30000000e-01,
         1.20000000e-01,   1.00000000e-01,   5.00000000e-01,
         4.20000000e-01,   4.20000000e-01,   7.00000000e-02,
         2.90000000e-01,   4.00000000e-02,   2.30000000e-01,
         2.30000000e-01,   7.00000000e-01,   2.30000000e-01,
         3.40000000e-01,   2.00000000e-02,   2.70000000e-01,
         3.00000000e-02,   3.10000000e-01,   2.00000000e-01,
         1.30000000e+00,   6.00000000e-02,   4.90000000e-01,
         1.00000000e-02,   1.00000000e-02,   1.50000000e-01,
         3.40000000e-01,   3.00000000e-02,   1.30000000e-01,
         8.00000000e-01,   3.60000000e-01,   2.70000000e-01,
         8.00000000e-02,   7.40000000e+00,   8.00000000e-01,
         1.10000000e-01,   1.70000000e-01,   3.00000000e-02,
         4.00000000e-02,   1.00000000e-02,   7.00000000e-01,
         4.00000000e-03,   2.80000000e-01,   1.00000000e-02,
         7.00000000e-01,   4.00000000e-03,   1.30000000e-01,
         6.00000000e-02,   6.00000000e-01,   2.00000000e-02,
         3.20000000e-01,   5.00000000e-01,   2.50000000e+00,
         3.80000000e-01,   2.10000000e-01,   1.80000000e-01,
         2.10000000e+00,   1.60000000e-01,   4.50000000e-01,
         1.20000000e+00,   9.60000000e+00,   4.20000000e-01,
         5.00000000e-01,   4.60000000e-01,   7.50000000e+00,
         9.00000000e-02,   7.00000000e-01,   3.00000000e+00,
         4.00000000e+01,   3.10000000e-01,   1.50000000e+00,
         4.00000000e-02,   5.00000000e-02,   3.40000000e-01,
         4.60000000e-01,   1.80000000e-01,   2.00000000e-02,
         3.90000000e+00,   5.00000000e-01,   8.90000000e+00,
         2.30000000e-01,   1.60000000e+02,   2.90000000e+00,
         3.70000000e+00,   3.30000000e+01,   5.50000000e+00,
         6.10000000e+01,   3.50000000e+00,   3.60000000e+01,
         7.70000000e+00,   9.10000000e+01,   2.10000000e+00,
         2.50000000e+01,   3.60000000e+00,   5.90000000e+01,
         3.60000000e+00,   1.00000000e+01,   4.70000000e+00,
         2.20000000e+01,   9.70000000e+00,   1.40000000e+01,
         9.10000000e+00,   2.30000000e+01,   6.10000000e+01,
         1.50000000e+01,   1.40000000e+01,   1.20000000e+01,
         5.00000000e+00,   1.10000000e+01,   6.80000000e+00,
         3.00000000e+01,   1.10000000e+01,   1.90000000e+01,
         9.90000000e+00,   3.20000000e+01,   8.00000000e+01,
         5.00000000e+01,   2.20000000e+01,   2.40000000e+01,
         7.30000000e+00,   3.30000000e+01,   5.30000000e+00,
         2.30000000e+01,   2.50000000e+01,   8.60000000e+01,
         1.00000000e+01,   2.60000000e+01,   1.30000000e+02,
         5.40000000e+02,   3.20000000e+01,   7.50000000e+01,
         1.10000000e+01,   1.11000000e+02,   1.00000000e+01,
         1.19000000e+02,   3.80000000e+01,   4.02000000e+02,
         2.60000000e+01,   2.98000000e+02,   5.10000000e+01,
         1.03500000e+03,   2.30000000e+01,   2.40000000e+02,
         1.30000000e+01,   3.80000000e+01,   6.50000000e+00,
         3.50000000e+01,   2.40000000e+01,   5.00000000e+01,
         2.00000000e+01,   5.40000000e+01,   1.50000000e+02,
         4.80000000e+01,   1.40000000e+01,   2.30000000e+01,
         4.60000000e+01,   1.40000000e+02,   1.80000000e+01,
         9.30000000e+01,   1.30000000e+02,   3.60000000e+02,
         4.80000000e+01,   1.33000000e+02,   7.80000000e+02,
         1.20000000e+03,   4.70000000e+01,   7.50000000e+01,
         2.00000000e+02,   1.90000000e+03,   1.80000000e+01,
         9.40000000e+01,   5.60000000e+02,   9.00000000e+03,
         4.20000000e+01,   2.28000000e+02,   4.10000000e+03,
         5.00000000e+04,   1.56000000e+02,   9.25000000e+02]
        )

bgunc_cms3l8=DoubleVector(  
      [ 3.        ,  6.        ,  2.        ,  0.72727273,  5.        ,
        6.        ,  7.5       ,  0.58333333,  5.        ,  1.42857143,
        5.        ,  1.        ,  2.        ,  0.44      ,  0.61538462,
        1.        ,  0.6       ,  0.54      ,  0.52380952,  0.45238095,
        0.85714286,  0.44827586,  1.        ,  0.56521739,  0.47826087,
        0.44285714,  0.56521739,  0.47058824,  1.5       ,  0.44444444,
        1.33333333,  0.48387097,  0.4       ,  0.38461538,  0.66666667,
        0.3877551 ,  2.        ,  6.        ,  1.06666667,  0.52941176,
        0.66666667,  0.69230769,  0.5       ,  0.52777778,  0.48148148,
        0.625     ,  0.47297297,  0.5       ,  0.72727273,  0.58823529,
        1.33333333,  1.        ,  3.        ,  0.47142857,  5.        ,
        0.57142857,  2.        ,  0.42857143,  5.        ,  0.61538462,
        0.66666667,  0.4       ,  2.        ,  0.625     ,  0.36      ,
        0.2       ,  0.52631579,  0.47619048,  0.33333333,  0.23809524,
        0.5       ,  0.53333333,  0.25      ,  0.16666667,  0.54761905,
        0.32      ,  0.39130435,  0.26666667,  0.66666667,  0.44285714,
        0.26666667,  0.25      ,  0.48387097,  0.31333333,  0.75      ,
        0.8       ,  0.44117647,  0.54347826,  0.72222222,  1.5       ,
        0.64102564,  0.42      ,  0.26966292,  0.39130435,  0.2125    ,
        0.27586207,  0.43243243,  0.42424242,  0.4       ,  0.49180328,
        0.4       ,  0.44444444,  0.35064935,  0.50549451,  0.38095238,
        0.4       ,  0.41666667,  0.49152542,  0.33333333,  0.48      ,
        0.34042553,  0.5       ,  0.34020619,  0.45714286,  0.37362637,
        0.47826087,  0.37704918,  0.32666667,  0.31428571,  0.48333333,
        0.32      ,  0.47272727,  0.35294118,  0.5       ,  0.34545455,
        0.33684211,  0.37373737,  0.5       ,  0.4       ,  0.22      ,
        0.28636364,  0.40833333,  0.2739726 ,  0.26363636,  0.28301887,
        0.47826087,  0.272     ,  0.26744186,  0.25      ,  0.42307692,
        0.31538462,  0.2962963 ,  0.203125  ,  0.25333333,  0.44545455,
        0.48648649,  0.53      ,  0.51260504,  0.39473684,  0.37810945,
        0.5       ,  0.50671141,  0.21568627,  0.24637681,  0.43478261,
        0.47083333,  0.26923077,  0.47368421,  0.44615385,  0.51428571,
        0.375     ,  0.5       ,  0.5       ,  0.51851852,  0.17333333,
        0.27083333,  0.34285714,  0.47826087,  0.21086957,  0.34285714,
        0.44444444,  0.50537634,  0.20769231,  0.25555556,  0.47916667,
        0.5112782 ,  0.15384615,  0.25833333,  0.27659574,  0.42666667,
        0.175     ,  0.28421053,  0.37222222,  0.44680851,  0.15535714,
        0.3       ,  0.26190476,  0.27631579,  0.16341463,  0.3       ,
        0.15384615,  0.28432432]
        )

data_cms3l8=IntVector([
     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     1,     0,     0,     1,     0,
     0,     0,     0,     1,     0,     0,     0,     1,     0,
     1,     0,     0,     0,     0,     0,     0,     0,     1,
     0,     0,     1,     0,     0,     0,     0,     0,     1,
     0,     5,     2,     0,     0,     0,     0,     0,     2,
     0,     0,     0,     1,     0,     0,     0,     3,     0,
     0,     1,     2,     1,     0,     0,     4,     0,     1,
     2,     9,     2,     0,     2,    15,     0,     0,     4,
    41,     1,     2,     0,     0,     0,     0,     2,     0,
     4,     0,     7,     1,   156,     4,     5,    35,     1,
    47,     3,    34,     8,    82,     4,    25,     1,    52,
     5,     2,     3,    19,     7,    18,     8,    21,    39,
    17,     9,    10,     4,    14,     6,    32,    10,    24,
    10,    25,    78,    70,    22,    36,     3,    41,     4,
    15,    26,   110,     5,    24,   135,   542,    31,    86,
     7,   101,    13,    87,    35,   406,    29,   269,    53,
   910,    29,   237,    18,    25,    10,    24,    21,    41,
    14,    42,   150,    39,    15,    19,    50,   169,    20,
    85,   142,   353,    48,   140,   773,  1276,    56,    81,
   178,  1676,    17,   115,   510,  9939,    34,   226,  3869,
 50188,   148,   906])
bg_zerolep8 = DoubleVector([
                           ])

bgunc_zerolep8 = DoubleVector([
                           ])

data_zerolep8 = IntVector([
                           ])
