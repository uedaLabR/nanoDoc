import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from sklearn.metrics import roc_curve
from sklearn.metrics import auc

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
        return [745,746,747,955,1618,1835,1911,1915,1917,1939,1962,2030,2069,2251,2445,2449,2457,2498,2501,2503,2504,2552,2580,2605]
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
    cnt = 0
    for line in lines:
        line = line.split(",")
        r.append(float(line[1]))
        cnt = cnt +1
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

def getAnswerLabel(dl,xcoords):
    r =[]
    for i in range(dl):
        if   (i-1 in xcoords) or (i in xcoords) or (i+1 in xcoords):
            r.append(1)
        else:
            r.append(0)
    return r

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

    cnt = 0
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
       cnt = cnt +1


    plt.rcParams['figure.figsize'] = (50, 10)

    #normalize
    maxv = 500
    data = np.clip(data, 0, maxv)
    data = (data / maxv)

    xcoords = getCoord(p)
    ksdata = dataKSStat(p)

    minlen = min(len(data), len(ksdata))
    y_true = getAnswerLabel(minlen, xcoords)

    rangestart = 0
    rangesend = 0
    title = ""
    if "16S" in p:
        rangestart = 500
        rangesend = 1538
        title = "E.coli 16S"
    if "23S" in p:
        rangestart = 700
        rangesend = 2700
        title = "E.coli 23S"
    if "18S" in p:
        rangestart = 0
        rangesend = 1796
        title = "Yeast 18S"
    if "25S" in p:
        rangestart = 600
        rangesend = 3000
        title = "Yeast 25S"
    print(minlen)
    ksdata = ksdata[rangestart :rangesend]
    data = data[rangestart : rangesend]
    cddata = cddata[rangestart : rangesend]
    y_true = y_true[rangestart:rangesend]

    y_score = data

    fpr, tpr, thresholds = roc_curve(y_true, y_score)
    auc_c = auc(fpr, tpr)
    plt.title(title)
    plt.plot(fpr, tpr,marker='o',markersize=4)
    plt.xlabel('FPR: False positive rate')
    plt.ylabel('TPR: True positive rate')
    plt.grid()
    print(auc_c)

    maxv = np.max(cddata)
    cddata = (cddata / maxv)*0.8

    y_score = cddata
    fpr, tpr, thresholds = roc_curve(y_true, y_score)
    auc_c = auc(fpr, tpr)
    plt.plot(fpr, tpr,marker='.',markersize=4)
    plt.xlabel('FPR: False positive rate')
    plt.ylabel('TPR: True positive rate')
    plt.grid()
    print(auc_c)



    y_score = ksdata
    fpr, tpr, thresholds = roc_curve(y_true, y_score)
    auc_c = auc(fpr, tpr)
    plt.plot(fpr, tpr,marker='x',markersize=4)
    plt.xlabel('FPR: False positive rate')
    plt.ylabel('TPR: True positive rate')
    plt.grid()
    print(auc_c)
    g = plt.subplot()
    g.set_aspect('equal')
    plt.show()
