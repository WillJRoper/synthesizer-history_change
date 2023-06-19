

import numpy as np
from dataclasses import dataclass
from .units import Quantity
from .utils import fnu_to_flam
from . import exceptions


interesting_lines = [
'H 1 1215.67A', # LyA
'H 1 6564.62A', # Ha
'H 1 4862.69A', # Hb
'H 1 4341.68A', # Hg
'H 1 4102.89A', # Hd
'H 1 3971.19A',
'H 1 3890.15A',
'H 1 3836.47A',
'H 1 1.87561m', # Pa
'H 1 1.28215m', # Pb
'H 1 1.09410m', # Pg
'H 1 1.00521m', # Pd
'H 1 9548.54A',
'H 1 9231.50A',
'H 1 4.05224m', # Bra
'H 1 2.62585m', # Brb
'H 1 2.16611m', # Brg
'H 1 1.94507m', # Brd
'H 1 1.81790m',
He 2 1640.43A
He 2 1215.13A
He 2 1084.94A
He 2 1025.27A
He 2 933.444A
He 2 930.337A
He 2 927.846A
He 2 4686.95A
He 2 3203.97A
He 2 1.01261m
He 1 1.08332m
He 1 3889.74A
He 1 3188.66A
He 1 2945.96A
He 1 2829.91A
He 1 2761.25A
He 1 2721.81A
He 1 2695.41A
He 1 2676.84A
He 1 2663.26A
He 1 2653.03A
He 1 2645.12A
He 1 2638.88A
He 1 2633.86A
He 1 2629.77A
He 1 2626.39A
He 1 2623.57A
He 1 2621.18A
He 1 2619.15A
He 1 2617.40A
He 1 2615.88A
He 1 2614.56A
He 1 2613.41A
He 1 2612.39A
He 1 2611.49A
He 1 2.05869m
He 1 5017.08A
He 1 3965.85A
He 1 3614.67A
He 1 3448.58A
He 1 3356.38A
He 1 3298.28A
He 1 3259.60A
He 1 3232.48A
He 1 3212.70A
He 1 3197.82A
He 1 3186.34A
He 1 3177.28A
He 1 3170.02A
He 1 3137.66A
He 1 7067.16A
He 1 5877.26A
He 1 4714.35A
He 1 4472.74A
He 1 4121.99A
He 1 4027.34A
He 1 3868.58A
He 1 3820.70A
He 1 3705.82A
He 1 3635.12A
He 1 3588.19A
He 1 3555.35A
He 1 3531.44A
He 1 3513.47A
He 1 3499.61A
He 1 3488.69A
He 1 3479.93A
He 1 3472.80A
He 1 3466.91A
He 1 3461.98A
He 1 3457.83A
He 1 3454.29A
He 1 3451.25A
He 1 3448.62A
He 1 3446.33A
He 1 3444.32A
He 1 3442.55A
He 1 3440.98A
He 1 7283.36A
He 1 6680.00A
He 1 5049.05A
He 1 4923.31A
He 1 4438.80A
He 1 4389.16A
He 1 4170.15A
He 1 4144.93A
He 1 4010.18A
He 1 3927.52A
He 1 3872.79A
He 1 3834.57A
He 1 3806.77A
He 1 3785.90A
He 1 3769.81A
He 1 3757.14A
He 1 3746.99A
He 1 3738.71A
He 1 3731.88A
He 1 3726.18A
He 1 3721.37A
He 1 3717.27A
He 1 3713.75A
He 1 3710.70A
He 1 3708.05A
He 1 3705.72A
He 1 3703.67A
He 1 3701.86A
He 1 4.29576m
He 1 1.25309m
He 1 9466.17A
He 1 8364.03A
He 1 7791.41A
He 1 7485.33A
He 1 7289.01A
He 1 1.50878m
He 1 1.10161m
He 1 9606.08A
He 1 2.11238m
He 1 1.70073m
He 1 1.28496m
He 1 1.19725m
He 1 1.06707m
He 1 1.03142m
He 1 9517.67A
He 1 9064.88A
He 1 8778.55A
He 1 8584.59A
He 1 8446.51A
He 1 8344.43A
He 1 8266.68A
He 1 8206.01A
He 1 8157.71A
He 1 8118.60A
He 1 8086.47A
He 1 8059.74A
He 1 8037.25A
He 1 8018.16A
He 1 7946.79A
He 1 1.95484m
He 1 1.86905m
He 1 1.29884m
He 1 1.27884m
He 1 1.09997m
He 1 1.09160m
He 1 1.00302m
He 1 9528.62A
He 1 9212.76A
He 1 8999.37A
He 1 8847.75A
He 1 8735.80A
He 1 8650.62A
He 1 8584.21A
He 1 8531.37A
He 1 8488.60A
He 1 8453.48A
He 1 8424.28A
He 1 8399.71A
He 1 8378.86A
He 1 8300.96A
He 1 1.87023m
He 1 1.85606m
He 1 1.27940m
He 1 1.09201m
He 1 1.00336m
He 1 9531.71A
He 1 9215.65A
He 1 9002.13A
He 1 8850.41A
He 1 8738.40A
He 1 8653.17A
He 1 8303.31A
He 1 2.11361m
He 1 1.90946m
He 1 1.34154m
He 1 1.29720m
He 1 1.10480m
He 1 1.01399m
He 1 9627.52A
He 1 9305.17A
He 1 9087.54A
He 1 8932.95A
He 1 8818.85A
He 1 8375.92A
He 1 2.04328m
He 1 3.70361m
He 1 2.47342m
He 1 2.06001m
He 1 1.85902m
He 1 1.74247m
He 1 1.66768m
He 1 1.61635m
He 1 1.57937m
He 1 1.44274m
He 1 4.24405m
He 1 4.03774m
He 1 2.66790m
He 1 2.61921m
He 1 2.16125m
He 1 1.94108m
He 1 1.81436m
He 1 1.73342m
He 1 1.67803m
He 1 1.63821m
He 1 1.60851m
He 1 1.49168m
He 1 4.04093m
He 1 2.62056m
He 1 2.16217m
He 1 1.94182m
He 1 1.81501m
He 1 4.04902m
He 1 2.62409m
He 1 2.16471m
He 1 1.94387m
He 1 1.81680m
He 1 1.73564m
He 1 1.68011m
He 1 1.64020m
He 1 1.49332m
He 1 4.04903m
He 1 2.62410m
He 1 2.16472m
He 1 1.94388m
He 1 4.12273m
He 1 2.18401m
He 1 4.64191m
He 1 3.73259m
He 1 3.29064m
He 1 4.65048m
He 1 3.73813m
He 1 3.29495m
He 1 3.03737m
He 1 4.65155m
He 1 3.73882m
He 1 3.29548m
He 1 4.65156m
He 1 4.81482m
Al 3 1854.72A
Al 3 1862.79A
Ar 1 1048.22A
Ar 1 1066.66A
Ar 3 3110.08A
Ar 3 5193.27A
Ar 3 7137.76A
Ar 3 7753.24A
Ar 4 2854.50A
Ar 4 2869.06A
Ar 4 4712.58A
Ar 4 4741.45A
Ar 4 7172.68A
C 1 1657.01A
C 1 9852.96A
C 2 1036.34A
C 2 1037.02A
C 2 1334.53A
C 2 1335.66A
C 2 1335.71A
C 2 1760.82A
C 2 2324.21A
C 2 2325.40A
C 2 2326.11A
C 2 2327.65A
C 2 2328.84A
C 2 7233.33A
C 3 1908.73A
C 3 1906.68A
C 3 977.020A
C 3 1175.26A
C 3 1175.59A
C 3 1176.37A
C 3 2297.58A
C 3 4651.55A
Ca 4 3.20698m
Co 2 2.98463m
Co 2 1.01907m
Co 3 6578.12A
Co 3 5890.11A
Cr 2 8002.28A
Cr 2 8127.53A
Cr 2 8231.93A
Fe 2 1.25702m
Fe 2 4417.51A
Fe 2 4288.60A
Fe 2 2600.17A
Fe 2 2586.65A
Fe 2 2382.77A
Fe 2 2344.21A
Fe 2 1.32092m
Fe 2 1.24888m
Fe 2 7639.64A
Fe 2 4890.98A
Fe 2 4493.89A
Fe 2 4459.20A
Fe 2 4360.56A
Fe 2 2626.45A
Fe 2 2612.65A
Fe 2 2599.15A
Fe 2 2396.36A
Fe 2 2389.36A
Fe 2 2365.55A
Fe 2 1.37219m
Fe 2 1.29462m
Fe 2 7689.05A
Fe 2 4729.39A
Fe 2 4415.02A
Fe 2 2632.11A
Fe 2 2618.40A
Fe 2 2607.87A
Fe 2 2405.62A
Fe 2 2381.49A
Fe 2 1.32814m
Fe 2 1.27913m
Fe 2 4453.35A
Fe 2 2631.83A
Fe 2 2614.60A
Fe 2 1.29813m
Fe 2 1.27069m
Fe 2 4476.16A
Fe 2 1.64400m
Fe 2 1.53389m
Fe 2 8619.32A
Fe 2 7157.13A
Fe 2 5274.81A
Fe 2 5160.21A
Fe 2 5113.05A
Fe 2 4815.88A
Fe 2 4776.05A
Fe 2 4245.16A
Fe 2 3377.17A
Fe 2 2494.05A
Fe 2 2360.72A
Fe 2 2348.83A
Fe 2 1.80989m
Fe 2 1.67733m
Fe 2 1.59991m
Fe 2 9054.43A
Fe 2 8894.35A
Fe 2 7454.59A
Fe 2 5434.64A
Fe 2 5263.09A
Fe 2 5221.51A
Fe 2 5159.44A
Fe 2 4948.75A
Fe 2 4906.71A
Fe 2 4875.85A
Fe 2 4348.07A
Fe 2 4278.03A
Fe 2 1.80051m
Fe 2 1.71159m
Fe 2 1.66422m
Fe 2 9401.62A
Fe 2 9229.15A
Fe 2 9035.98A
Fe 2 5335.13A
Fe 2 4354.00A
Fe 2 4320.83A
Fe 2 1.79759m
Fe 2 1.74541m
Fe 2 9473.53A
Fe 2 9270.11A
Fe 2 1.81189m
Fe 2 3278.29A
Fe 2 2881.60A
Fe 2 2756.55A
Fe 2 2740.36A
Fe 2 2750.13A
Fe 3 4926.04A
Fe 3 4882.48A
Fe 3 4659.31A
Fe 3 4608.40A
Fe 3 5271.87A
Fe 3 4988.73A
Fe 3 4755.97A
Fe 3 4702.94A
Fe 3 4668.25A
Fe 3 5012.56A
Fe 3 4770.86A
Fe 3 4735.16A
Fe 3 5413.48A
Fe 3 4931.82A
Fe 3 4778.95A
Fe 3 5086.08A
K 3 4.61808m
Mg 1 2026.48A
Mg 1 2852.96A
Mg 1 4563.88A
Mg 1 4572.38A
Mg 1 4576.57A
Mn 2 2576.88A
Mn 2 2594.50A
Mn 2 2606.46A
N 1 952.303A
N 1 952.415A
N 1 953.415A
N 1 953.655A
N 1 953.970A
N 1 954.104A
N 1 963.990A
N 1 964.626A
N 1 965.041A
N 1 1134.17A
N 1 1134.41A
N 1 1134.98A
N 1 1168.33A
N 1 1168.54A
N 1 1199.55A
N 1 1200.22A
N 1 1200.71A
N 1 5199.35A
N 1 5201.71A
N 1 7470.37A
N 1 8190.26A
N 1 8218.59A
N 1 8225.39A
N 1 8244.66A
N 1 8682.67A
N 1 1.04006m
N 5 1238.82A
N 5 1242.80A
Ne 3 3869.86A
Ne 3 3968.59A
Ne 3 1814.56A
Ne 3 3343.14A
O 1 1039.23A
O 1 1040.94A
O 1 1041.69A
O 1 1302.17A
O 1 1304.86A
O 1 1306.03A
O 1 5578.89A
O 1 6302.05A
O 1 6365.54A
O 1 8448.68A
O 2 2470.97A
O 2 2471.09A
O 2 7320.94A
O 2 7322.01A
O 2 7331.69A
O 2 7332.75A
O 2 3729.88A
O 2 3727.09A
O 3 4932.60A
O 3 4960.29A
O 3 5008.24A
O 3 2321.66A
O 3 4364.44A
O 3 1660.81A
O 3 1666.15A
P 2 1152.82A
P 2 1.14713m
P 2 1156.97A
P 2 1155.01A
P 2 1.18861m
P 2 1159.09A
P 2 7878.16A
P 3 1334.81A
P 3 1344.85A
P 4 950.657A
P 4 1028.09A
P 4 1030.51A
P 4 1035.52A
S 2 6732.67A
S 2 6718.29A
S 2 4077.50A
S 2 4069.75A
S 2 1.03392m
S 2 1.02895m
S 2 1.03733m
S 2 1.03233m
S 2 1259.52A
S 2 1253.81A
S 2 1250.58A
S 2 912.736A
S 3 8831.82A
S 3 9071.11A
S 3 9533.23A
S 3 3722.69A
S 3 3798.25A
S 3 6313.81A
S 3 1713.11A
S 3 1728.94A
S 3 1190.20A
S 3 1194.45A
S 3 1202.12A
S 3 1194.06A
S 3 1201.73A
S 3 1200.97A
S 3 1015.50A
S 3 1012.50A
S 3 1015.57A
S 3 1021.11A
S 3 1015.78A
S 3 1021.32A
Si 2 2335.12A
Si 2 2329.23A
Si 2 1808.01A
Si 2 1526.71A
Si 2 1304.37A
Si 2 1260.42A
Si 2 1193.29A
Si 2 1190.42A
Si 2 2350.89A
Si 2 2344.92A
Si 2 2335.32A
Si 2 1817.45A
Si 2 1816.93A
Si 2 1533.43A
Si 2 1309.28A
Si 2 1265.00A
Si 2 1264.74A
Si 2 1197.39A
Si 2 1194.50A
Si 3 1892.03A
Si 3 1882.71A
Si 3 1206.50A
Si 4 1122.48A
Si 4 1128.34A
Si 4 1393.75A
Si 4 1402.77A
Si 4 1722.53A
Si 4 1727.38A
Al 2 2669.95A
Al 2 2661.15A
Al 2 1670.79A
C 4 1550.78A
C 4 1548.19A
Ca 2 7325.90A
Ca 2 7293.48A
Ca 2 3969.59A
Ca 2 8664.52A
Ca 2 3934.78A
Ca 2 8544.44A
Cl 2 8581.05A
Cl 2 9126.10A
Cl 3 5539.41A
Cl 3 5519.24A
Cl 3 3354.14A
Cl 3 8502.35A
Cl 3 3343.76A
Cl 3 8435.98A
Cl 3 8483.19A
Cl 4 7532.62A
Cl 4 8047.84A
Fe 4 3095.86A
Fe 4 3095.43A
Fe 4 2836.57A
Fe 4 2830.19A
Fe 4 2824.33A
Fe 4 2578.69A
Fe 4 2.83640m
Fe 4 2.80631m
Fe 4 2.86525m
Fe 4 2568.38A
Fe 4 2.71643m
Fe 4 2568.17A
Fe 4 2.71415m
Fe 4 2.77400m
Fe 4 6741.68A
Fe 4 6763.25A
Fe 4 6736.27A
Fe 4 6999.03A
Fe 4 6794.36A
Fe 4 5235.22A
Fe 4 5034.99A
Fe 4 4907.93A
Fe 4 4919.35A
Fe 4 7224.75A
Fe 4 4901.34A
Fe 4 4904.44A
Fe 4 7185.97A
Fe 4 4869.52A
Fe 4 7193.21A
Fe 4 4869.31A
Fe 4 4207.78A
Fe 4 4210.06A
Fe 4 4145.37A
Mg 2 2803.53A
Mg 2 2796.35A
Mg 2 2937.37A
Mg 2 2929.49A
Mg 2 2798.82A
Mg 2 2791.60A
N 2 6549.86A
N 2 6585.27A
N 2 3063.79A
N 2 5756.21A
N 2 1083.99A
N 2 915.612A
N 2 2139.68A
N 2 1084.58A
N 2 1084.56A
N 2 916.012A
N 2 916.020A
N 2 915.962A
N 2 2143.46A
N 2 1085.70A
N 2 1085.55A
N 2 1085.53A
N 2 916.701A
N 2 916.710A
N 3 1748.65A
N 3 1753.99A
N 3 1746.82A
N 3 1752.16A
N 3 1749.67A
N 3 989.799A
N 3 991.511A
N 3 991.577A
N 4 1486.50A
N 4 1483.32A
Ni 2 1.19134m
Ni 2 1.07181m
Ni 2 7379.86A
Ni 2 8303.27A
Ni 2 1.93930m
Ni 2 6668.64A
Ni 2 7413.65A
Ni 2 4327.45A
Ni 2 4629.34A
Ni 2 7257.82A
Ni 2 4202.35A
Ni 2 3994.19A
Ni 2 3439.86A
P 5 1128.01A
P 5 1117.98A
S 4 1416.89A
S 4 1406.02A
S 4 1062.66A
S 4 1073.52A
S 4 1072.97A
S 6 944.524A
S 6 933.380A
C 1 943.939A
C 1 1247.50A
C 1 1314.77A
C 1 1165.98A
N 1 942.258A
O 1 976.791A
O 1 969.049A
O 1 946.387A
O 1 934.501A
O 1 973.577A
O 1 927.472A
O 1 948.506A
O 1 922.969A
O 1 919.912A
O 1 935.672A
O 1 917.726A
O 1 918.704A
O 1 928.186A
O 1 917.473A
O 1 923.438A
O 1 920.243A
O 1 917.970A
Mg 1 1817.31A
Mg 1 1739.09A
Mg 2 959.306A
Mg 2 919.202A
Al 2 929.014A
Si 2 1179.59A
Si 2 993.354A
Si 2 1024.03A
Ar 1 918.493A
P 2 967.202A
P 2 963.828A
P 3 918.147A
Cl 2 1071.32A
Cl 3 1011.29A
Cl 4 981.258A
Ti 2 3365.71A
Ti 2 3238.08A
Ti 3 1288.59A
Cr 2 2060.44A
Cr 3 1033.06A
Cr 3 1035.74A
Cr 3 924.009A
Mn 2 1198.85A
Co 2 1470.18A
Co 3 942.738A
Ni 2 1751.83A
Ni 2 1744.24A
Ni 2 1484.97A
Ni 2 1472.74A
Ni 2 1400.37A
Ni 2 1375.73A
Ni 2 1324.11A
C 2R 1335.00A
C 3R 977.000A
C 3R 1909.00A
C 3R 1175.00A
N 2R 6584.00A
N 2R 1085.00A
N 3R 990.000A
O 2R 3729.00A
O 2R 4651.00A
O 2R 4342.22A
O 2R 3737.06A


]



def get_line_id(id):
    """
    A class representing a spectral line or set of lines (e.g. a doublet)

    Arguments
    ----------
    id : str, list, tuple
        a str, list, or tuple containing the id(s) of the lines

    Returns
    -------
    string
        string representation of the id

    """'H 1 6564.62A'

    if isinstance(id, list):
        return ','.join(id)
    else:
        return id


def get_roman_numeral(number):
    """
    Function to convert an integer into a roman numeral str.

    Used for renaming emission lines from the cloudy defaults.

    Returns
    ---------
    str
        string reprensentation of the roman numeral
    """

    num = [1, 4, 5, 9, 10, 40, 50, 90,
           100, 400, 500, 900, 1000]
    sym = ["I", "IV", "V", "IX", "X", "XL",
           "L", "XC", "C", "CD", "D", "CM", "M"]
    i = 12

    roman = ''
    while number:
        div = number // num[i]
        number %= num[i]

        while div:
            roman += sym[i]
            div -= 1
        i -= 1
    return roman


def get_fancy_line_id(id):
    """
    Function to get a nicer line id, e.g. "O 3 5008.24A" -> "OIII5008"

    Returns
    ---------
    str
        the fancy line id
    """

    element, ion, wavelength = id.split(' ')

    wavelength = float(wavelength[:-1])

    return f'{element}{get_roman_numeral(int(ion))}{wavelength: .4g}'


class LineRatios:

    """
    A dataclass holding useful line ratios (e.g. R23) and diagrams (pairs of ratios), e.g. BPT.
    """
    # short-hand


    def __init__(self):

        O3 = ['O 3 4960.29A', 'O 3 5008.24A']
        O2 = ['O 2 3727.09A', 'O 2 3729.88A']
        Hb = 'H 1 4862.69A'
        Ha = 'H 1 6564.62A'

        self.ratios = {}

        # Balmer decrement, should be ~2.86 for dust free
        self.ratios['BalmerDecrement'] = [[Ha], [Hb]]
        self.ratios['N2'] = [['N 2 6585.27A'], [Ha]]  #  add reference
        self.ratios['S2'] = [['S 2 6732.67A', 'S 2 6718.29A'], [Ha]]  #  add reference
        self.ratios['O1'] = [['O 1 6302.05A'], [Ha]]  #  add reference
        self.ratios['R2'] = [['O 2 3727.09A'], [Hb]]  #  add reference
        self.ratios['R3'] = R3 = [['O 3 5008.24A'], [Hb]]  #  add reference
        self.ratios['R23'] = [O3+O2, [Hb]]  #  add reference
        self.ratios['O32'] = [['O 3 5008.24A'], ['O 2 3727.09A']]  #  add reference
        self.ratios['Ne3O2'] = [['Ne 3 3968.59A'], ['O 2 3727.09A']]  #  add reference

        self.available_ratios = tuple(self.ratios.keys())

        self.diagrams = {}
        self.diagrams['OHNO'] = [R3, [['Ne 3 3869.86A'], O2]]  #  add reference
        self.diagrams['BPT-NII'] = [[['N 2 6585.27A'], [Ha]], R3]  #  add reference
        # diagrams['VO78'] = [[], []]
        # diagrams['unVO78'] = [[], []]

        self.available_diagrams = tuple(self.diagrams.keys())

    def get_ratio_label_(self, ab, fancy=False):
        """
        Get a line ratio label

        Arguments
        -------
        ab
            a list of lists of lines, e.g. [[l1,l2], [l3]]
        fancy
            flag to return fancy label instead

        Returns
        -------
        str
            a label
        """

        a, b = ab

        if fancy:
            a = map(get_fancy_line_id, a)
            b = map(get_fancy_line_id, b)

        return f"({'+'.join(a)})/({'+'.join(b)})"

    def get_ratio_label(self, ratio_id, fancy=False):
        """
        Get a line ratio label

        Arguments
        -------
        ratio_id
            a ratio_id where the ratio lines are defined in LineRatios

        Returns
        -------
        str
            a label
        """

        ab = self.ratios[ratio_id]

        return f'{ratio_id}={self.get_ratio_label_(ab, fancy = fancy)}'


class LineCollection:

    """
    A class holding a collection of emission lines

    Attributes
    ----------
    lines : dictionary of Line objects

    Methods
    -------

    """

    def __init__(self, lines):

        self.lines = lines
        self.line_ids = list(self.lines.keys())

        # these should be filtered to only show ones that are available for the availalbe line_ids

        self.lineratios = LineRatios()

        self.available_ratios = self.lineratios.available_ratios
        self.available_diagrams = self.lineratios.available_diagrams

    def __getitem__(self, line_id):

        return self.lines[line_id]

    def __str__(self):
        """Function to print a basic summary of the LineCollection object.

        Returns a string containing the id, wavelength, luminosity, equivalent width, and flux if generated.

        Returns
        -------
        str
            Summary string containing the total mass formed and lists of the available SEDs, lines, and images.
        """

        # Set up string for printing
        pstr = ""

        # Add the content of the summary to the string to be printed
        pstr += "-"*10 + "\n"
        pstr += f"LINE COLLECTION\n"
        pstr += f"lines: {self.line_ids}\n"
        pstr += f"available ratios: {self.available_ratios}\n"
        pstr += f"available diagrams: {self.available_diagrams}\n"
        pstr += "-"*10

        return pstr

    def get_ratio_(self, ab):
        """
        Measure (and return) a line ratio

        Arguments
        -------
        ab
            a list of lists of lines, e.g. [[l1,l2], [l3]]

        Returns
        -------
        float
            a line ratio
        """

        a, b = ab

        return np.sum([self.lines[l].luminosity for l in a]) / \
            np.sum([self.lines[l].luminosity for l in b])

    def get_ratio(self, ratio_id):
        """
        Measure (and return) a line ratio

        Arguments
        -------
        ratio_id
            a ratio_id where the ratio lines are defined in LineRatios

        Returns
        -------
        float
            a line ratio
        """

        ab = self.lineratios.ratios[ratio_id]

        return self.get_ratio_(ab)

    def get_ratio_label_(self, ab, fancy=False):
        """
        Get a line ratio label

        Arguments
        -------
        ab
            a list of lists of lines, e.g. [[l1,l2], [l3]]
        fancy
            flag to return fancy label instead

        Returns
        -------
        str
            a label
        """

        a, b = ab

        if fancy:
            a = map(get_fancy_line_id, a)
            b = map(get_fancy_line_id, b)

        return f"({','.join(a)})/({','.join(b)})"

    def get_ratio_label(self, ratio_id, fancy=False):
        """
        Get a line ratio label

        Arguments
        -------
        ratio_id
            a ratio_id where the ratio lines are defined in LineRatios

        Returns
        -------
        str
            a label
        """
        ab = self.lineratios.ratios[ratio_id]

        return f'{ratio_id}={self.get_ratio_label_(ab, fancy = fancy)}'

    def get_diagram(self, diagram_id):
        """
        Return a pair of line ratios for a given diagram_id (E.g. BPT)

        Arguments
        -------
        ratdiagram_idio_id
            a diagram_id where the pairs of ratio lines are defined in LineRatios

        Returns
        -------
        tuple (float)
            a pair of line ratios
        """
        ab, cd = self.lineratios.diagrams[diagram_id]

        return self.get_ratio_(ab), self.get_ratio_(cd)

    def get_diagram_label(self, diagram_id, fancy=False):
        """
        Get a line ratio label

        Arguments
        -------
        ab
            a list of lists of lines, e.g. [[l1,l2], [l3]]

        Returns
        -------
        str
            a label
        """
        ab, cd = self.lineratios.diagrams[diagram_id]

        return self.get_ratio_label_(ab, fancy=fancy), self.get_ratio_label_(cd, fancy=fancy)


class Line:

    """
    A class representing a spectral line or set of lines (e.g. a doublet)

    Attributes
    ----------
    lam : wavelength of the line

    Methods
    -------

    """

    wavelength = Quantity()
    continuum = Quantity()
    luminosity = Quantity()
    flux = Quantity()
    ew = Quantity()

    def __init__(self, id_, wavelength_, luminosity_, continuum_):

        self.id_ = id_

        # --- these are maintained because we may want to hold on to the individual lines of a doublet
        self.wavelength_ = wavelength_
        self.luminosity_ = luminosity_
        self.continuum_ = continuum_

        self.id = get_line_id(id_)
        self.continuum = np.mean(continuum_)  #  mean continuum value in units of erg/s/Hz
        self.wavelength = np.mean(wavelength_)  # mean wavelength of the line in units of AA
        self.luminosity = np.sum(luminosity_)  # total luminosity of the line in units of erg/s/Hz
        self.flux = None  # line flux in erg/s/cm2, generated by method

        # continuum at line wavelength, erg/s/AA
        self._continuum_lam = fnu_to_flam(self._wavelength, self._continuum)
        self.ew = self._luminosity / self._continuum_lam  # AA

    def __str__(self):
        """Function to print a basic summary of the Line object.

        Returns a string containing the id, wavelength, luminosity, equivalent width, and flux if generated.

        Returns
        -------
        str
            Summary string containing the total mass formed and lists of the available SEDs, lines, and images.
        """

        # Set up string for printing
        pstr = ""

        # Add the content of the summary to the string to be printed
        pstr += "-"*10 + "\n"
        pstr += f"SUMMARY OF {self.id}" + "\n"
        pstr += f"wavelength: {self.wavelength:.1f}" + "\n"
        pstr += f"log10(luminosity/{self.luminosity.units}): {np.log10(self.luminosity):.2f}" + "\n"
        pstr += f"equivalent width: {self.ew:.0f}" + "\n"
        if self._flux:
            pstr += f"log10(flux/{self.flux.units}): {np.log10(self.flux):.2f}"
        pstr += "-"*10

        return pstr

    def __add__(self, second_line):
        """
        Function allowing adding of two Line objects together. This should NOT be used to add different lines together.

        Returns
        -------
        obj (Line)
            New instance of Line
        """

        if second_line.id == self.id:

            return Line(self.id, self._wavelength, self._luminosity + second_line._luminosity, self._continuum + second_line._continuum)

        else:

            exceptions.InconsistentAddition('Wavelength grids must be identical')

    def get_flux(self, cosmo, z):
        """Calculate the line flux in units of erg/s/cm2

        Returns the line flux and (optionally) updates the line object.

        Parameters
        -------
        cosmo: obj
            Astropy cosmology object

        z: float
            Redshift

        Returns
        -------
        flux: float
            Flux of the line in units of erg/s/cm2
            """

        luminosity_distance = cosmo.luminosity_distance(
            z).to('cm').value  # the luminosity distance in cm

        self.flux = self._luminosity / (4 * np.pi * luminosity_distance**2)

        return self.flux
