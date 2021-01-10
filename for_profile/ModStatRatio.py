import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
from matplotlib_venn import venn3
from matplotlib_venn import venn2
from Bio import SeqIO
import seaborn as sns
import statistics
import math

modf = "C:\\Users\\hirok\\Desktop\\covid19\\modwithr.txt"
dfmod = pd.read_csv(modf, sep='\t')

print(dfmod)

#print(dfmod.query('chrom == "ec16s"'))

path = "C:\\Users\\hirok\\Desktop\\covid19\\detection200730\\*.txt"
l = glob.glob(path)

def thinouted(lst, newsize):
    """間引いたリストを返す"""
    result = []
    cnt = 0
    for x in lst:
        cnt -= newsize
        if cnt < 0:
            cnt += len(lst)
            result.append(x)
    return result

holder = []
holdertmp = []
nomodscore = []
yl = []
Nm = []
m6Al = []
x = []
y = []
for p in l:

    print(p)
    chrom = ""
    if "16S" in p:
        chrom = "ec16s"
        continue
    if "23S" in p:
        chrom = "ec23s"
        continue
    if "18S" in p:
        chrom = "sc18s"
    if "25S" in p:
        chrom = "sc25s"

    f = open(p)
    lines = f.readlines()
    posin = 2
    for line in lines:
        line = line.split("\t")
        posin = posin +1
        line[9] = line[9].replace('\n', '')
        # pos = int(line[0])
        pos = posin
        score = float(line[11])
        if score > 500:
            score = 500
        score = score /500

        df = dfmod.query('chrom == "' +chrom+'" & pos == '+str(pos))
        if df.empty:

           if posin %30 == 0:
               df0 = dfmod.query('chrom == "' + chrom + '" & pos == ' + str(pos - 3))
               df1 = dfmod.query('chrom == "' + chrom + '" & pos == ' + str(pos - 2))
               df2 = dfmod.query('chrom == "' +chrom+'" & pos == '+str(pos-1))
               df3 = dfmod.query('chrom == "' + chrom + '" & pos == ' + str(pos +1))
               df4 = dfmod.query('chrom == "' + chrom + '" & pos == ' + str(pos +2))
               df5 = dfmod.query('chrom == "' + chrom + '" & pos == ' + str(pos + 3))
               if (df0.empty and df1.empty and df2.empty and df3.empty and df4.empty and df5.empty):
                    holdertmp.append(("No Mod",score))
                    nomodscore.append(score)

        else:

            mod = df["mod"].values[0]
            ratio = int(df["ratio"].values[0])/100
            print(pos,mod,score,ratio)
            x.append(ratio)
            y.append(score)
            holder.append((mod, score))
            #
            if mod == "Y":
                yl.append(score)
            if mod == "Am" or mod == "Cm" or mod == "Gm" or mod == "Um" :
                Nm.append(score)
            if mod == "m6A":
                m6Al.append(score)

    f.close()

plt.scatter(x,y)
plt.xlabel('Ratio')
plt.ylabel('Score')
plt.grid()


plt.show()