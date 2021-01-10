import glob
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd

def exp_func(x, a, b):
    return b*(x**a)

def exp_fit(val1_quan, val2_quan):
    l_popt, l_pcov = curve_fit(exp_func, val1_quan, val2_quan, maxfev=10000, check_finite=False)
    return exp_func(val1_quan, *l_popt),l_popt


path = "C:\\Users\\hirok\\Desktop\\covid19\\depthstat\\*"
l = glob.glob(path)

dic={}
for p in l:
    # print(p)
    f = open(p)
    lines = f.readlines()

    xv = []
    yv_1 = []
    yv_2 = []
    nuc =""
    for line in lines:
        data = line.split("\t")
        nuc = data[0]
        key = data[1]

        xv.append(float(data[1]))
        xv.append(float(data[1]))

        thres1_1 = int(data[3])
        thres1_2 = int(data[4])
        yv_1.append(thres1_1)
        yv_1.append(thres1_2)

        thres2_1 = int(data[5])
        thres2_2 = int(data[6])

        yv_2.append(thres2_1)
        yv_2.append(thres2_2)

        #


        if key in dic:
            ls = dic[key]
        else:
            dic[key] = []
        dic[key].append(float(data[2]))
    # print(xv)
    # print(yv_1)
    y_fit, l_popt1 = exp_fit(xv, yv_1)
    y_fit, l_popt2 = exp_fit(xv, yv_2)

    print(nuc,l_popt1[0],l_popt1[1],l_popt2[0],l_popt2[1])

    f.close()



sds = []
x=[]
for key,val in dic.items():

    # plt.plot(val)
    # plt.show()
    mean = np.mean(val)
    sd = np.std(val, ddof=1)
    # print(key,mean,sd)
    x.append(int(key))
    sds.append(sd)

y_fit,l_popt=exp_fit(x,sds)

ax=plt.subplot(1,1,1)
ax.plot(x,sds,label='obs')
ax.plot(x,y_fit,label='model')
plt.show()
print('a : {},   b : {}'.format(l_popt[0],l_popt[1])) #求めたパラメータa,bを確認