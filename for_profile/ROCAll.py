import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from sklearn.metrics import roc_curve
from sklearn.metrics import auc
from sklearn.metrics import matthews_corrcoef

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


def getAnswerLabel(dl,xcoords):
    r =[]
    for i in range(dl):
        if  (i-1 in xcoords) or (i in xcoords) or (i+1 in xcoords) :
            r.append(1)
        else:
            r.append(0)
    return r


import math

dataall = []
cddataall = []
depths = []
answerlabel = []
kslist =[]

for p in l:

    data = []
    cddata = []
    print(p)
    f = open(p)
    lines = f.readlines()

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

    xcoords = getCoord(p)
    ksdata = dataKSStat(p)

    rangestart = 0
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

    minlen = min(len(data), len(ksdata))
    print(minlen)
    y_true = getAnswerLabel(minlen, xcoords)

    ksdata = ksdata[rangestart :rangesend]
    data = data[rangestart : rangesend]
    cddata = cddata[rangestart : rangesend]
    y_true = y_true[rangestart:rangesend]

    plt.plot(data)
    plt.show()

    kslist.extend(ksdata)
    dataall.extend(data)
    cddataall.extend(cddata)
    answerlabel.extend(y_true)

g = plt.subplot()
g.set_aspect('equal')
plt.xlabel('FPR: False positive rate')
plt.ylabel('TPR: True positive rate')
plt.grid()


fpr, tpr, thresholds = roc_curve(answerlabel, dataall)
auc_c = auc(fpr, tpr)
plt.plot(fpr, tpr,marker='o',markersize=3)
print(auc_c)

for n in range(len(fpr)):
    print(fpr[n],thresholds[n])


fpr, tpr, thresholds = roc_curve(answerlabel, cddataall)
auc_c = auc(fpr, tpr)
plt.plot(fpr, tpr,marker='.',markersize=3)
print(auc_c)

fpr, tpr, thresholds = roc_curve(answerlabel, kslist)
auc_c = auc(fpr, tpr)
plt.plot(fpr, tpr,marker='x',markersize=3)
print(auc_c)

plt.show()

maxd = max(dataall)
plt.plot(dataall)
plt.show()
dataall = np.array(dataall)
dataall = dataall/maxd
plt.plot(dataall)
plt.show()

