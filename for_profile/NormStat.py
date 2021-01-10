import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import math

for time in (0,5,10,15,20,25,30):

    x =[]
    y =[]
    max = 0
    for v in range(-50,50):

        vv = (v/50) *4
        yadd = time - ((time/4) * np.abs(vv))
        yy = norm.pdf(vv, 0, 1) * yadd
        x.append(vv)
        y.append(yy)
        if yy > max:
            max = yy

    y = np.array(y)
    y = y /max

    # cnt = 0
    # for yy in y:
    #     cnt = cnt +1
    #     if yy > 0.05:
    #         print(cnt)

    plt.plot(x,y)
    plt.show()