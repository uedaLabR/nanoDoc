import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

path = "C:\\Users\\hirok\\Desktop\\covid19\\detection200730\\*.txt"
l = glob.glob(path)

def isNumeric(s):
    s = s.replace('.', '').replace('0', '').replace('-', '')
    if len(s) == 0:
        return True
    return s.isnumeric()

def getCoord(p):

    if "16S" in p:
        return [516, 527,966,967,1207,1402,1407,1498,1516,1518,1519]
    if "23S" in p:
        #return [745,746,747,955,1618,1835,1911,1915,1917,1939,1962,2030,2069,2251,2445,2449,2457,2498,2501,2503,2504,2552,2580,2604,2605]
        return [745, 746, 747, 955, 1618+1, 1835+1, 1911+1, 1915+1, 1917+1, 1939+1, 1962+1, 2030+1, 2069+1, 2251+1, 2445+1, 2449+1, 2457+1, 2498+1,
                2501+1, 2503+1, 2504+1, 2552+1, 2580+1, 2604+1, 2605+1]
    if "18S" in p:
        return [28,100,106,120,211,302,414,420,436,466,541,562,578,619,632,759,766,796,974,999,1007,1126,1181,1187,1191,1269,1271,1280,1290,1415,1428,1572,1575,1639,1773,1781,1782]
    if "25S" in p:
        return [645,649,650,663,776,805,807,817,867,876,898,908,956,960,966,986,990,1004,1042,1052,1056,1110,1124,1133,1437,1449,1450,1888,2129,2133,2142,2191,2197,2220,2256,2258,2260,2264,2266,2278,2280,2281,2288,2314,2337,2340,2347,2349,2351,2416,2417,2421,2619,2634,2640,2724,2729,2735,2791,2793,2815,2826,2843,2865,2870,2880,2921,2922,2923,2944,2946,2948,2959,2975]

    return []

def dataKSStat(p):

    path = ""
    r = []
    r.append(0)
    r.append(0)
    r.append(0)
    r.append(0)
    r.append(0)
    r.append(0)
    if "16S" in p:
        path = "C:\\Users\\hirok\\Desktop\\covid19\\ksstats\\16S.csv"
    if "23S" in p:
        path = "C:\\Users\\hirok\\Desktop\\covid19\\ksstats\\23S.csv"
    if "18S" in p:
        path = "C:\\Users\\hirok\\Desktop\\covid19\\ksstats\\18S.csv"
    if "25S" in p:
        path = "C:\\Users\\hirok\\Desktop\\covid19\\ksstats\\25S.csv"

    f = open(path,encoding="utf-8")
    lines = f.readlines()
    for line in lines:
        line = line.split(",")
        r.append(float(line[1]))
    f.close()
    return r

def wavelet(data,i):

    # maxv= 0
    # for sigma in (0.1,0.2,0.3,0.4,0.5):
    #     innerproduct = 0
    #     normfac = 0
    #
    #     if innerproduct > maxv:
    #         maxv = innerproduct
    innerproduct = 0
    normfac = 0
    for n in range(-7, 8):
        dist = norm.pdf(n/2.5, 0, 1.0)
        normfac = normfac + dist
        innerproduct = innerproduct + (dist * data[i + n])

    innerproduct = innerproduct / normfac

    return innerproduct

data =[]
rnum = 0
import math

for p in l:

    print(p)
    f = open(p)
    lines = f.readlines()
    data = []
    cddata = []
    depths = []
    data.append(0)
    data.append(0)
    data.append(0)

    cddata.append(0)
    cddata.append(0)
    cddata.append(0)

    for line in lines:

       line = line.split("\t")
       line[9] = line[9].replace('\n', '')
       pos = int(line[0])
       depth = int(line[2])
       score = float(line[11])

       currentdiff = abs(float(line[4])- float(line[6]))
       cddata.append(currentdiff)

       data.append(score)
       depths.append(depth)
       # print(pos,depth,score)

    plt.rcParams['figure.figsize'] = (50, 10)

    xcoords = getCoord(p)
    ksdata = dataKSStat(p)

    fig = plt.figure()
    ax1 = fig.add_subplot(2, 1, 1)

    for xc in xcoords:
       ax1.axvline(xc,color='red', alpha= .25)

    #normalize
    maxv = 500
    data = np.clip(data, 0, maxv)
    data = (data / maxv)
    ax1.plot(data,lw=1)
    ax1.set_ylim(-0.1,1)

    ax2 = fig.add_subplot(2, 1, 2)

    for xc in xcoords:
       ax2.axvline(xc,color='red', alpha= .25)
    print(ksdata)
    ax2.plot(ksdata, lw=1)
    ax2.set_ylim(-0.1, 1)

    maxv = np.max(cddata)
    cddata = (cddata / maxv)*0.8

    # ax3 = fig.add_subplot(3, 1, 3)
    # print(len(cddata),cddata)
    # for xc in xcoords:
    #    ax3.axvline(xc,color='red', alpha= .25)
    #
    # ax3.plot(cddata, lw=1)
    # ax3.set_ylim(-0.1, 1)
    #fig.tight_layout()
    title=""
    if "16S" in p:
        title="16S"
    if "23S" in p:
        title="23S"
    if "18S" in p:
        title="18S"
    if "25S" in p:
        title = "25S"

    plt.show()


    # dataSmooth = np.zeros(len(data))
    # waveunit = 7
    # for i in range(len(data)):
    #
    #     if i < waveunit or i > len(data)-waveunit-1:
    #         dataSmooth[i] = data[i]
    #     else:
    #         dataSmooth[i] = wavelet(data,i)
    #
    # for xc in xcoords:
    #    plt.axvline(xc,color='red', alpha= .25)
    # plt.plot(dataSmooth,lw=1)
    # plt.show()


