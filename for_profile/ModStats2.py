import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
from matplotlib_venn import venn3
from matplotlib_venn import venn2
from Bio import SeqIO
import seaborn as sns
import statistics
import math

modf = "C:\\Users\\hirok\\Desktop\\covid19\\mod.txt"
dfmod = pd.read_csv(modf, sep='\t')

print(dfmod)

path = "C:\\Users\\hirok\\Desktop\\covid19\\detection200730\\*.txt"
l = glob.glob(path)
ndict = {}

nuclist = {"CGCCC","CGTAA","GGAAC","GAACC","TTAAA","ACTAT","ACTCT"}

for p in l:

    print(p)
    chrom = ""
    if "16S" in p:
        chrom = "ec16s"
    if "23S" in p:
        chrom = "ec23s"
    if "18S" in p:
        chrom = "sc18s"
    if "25S" in p:
        chrom = "sc25s"

    f = open(p)
    lines = f.readlines()
    posin = 3
    if (chrom == "sc25s" or chrom == "sc18s"):
        posin = 2
    for line in lines:
        line = line.split("\t")
        posin = posin +1
        fmer = line[1]
        line[9] = line[9].replace('\n', '')
        # pos = int(line[0])
        pos = posin
        if(chrom == "ec23s" and pos > 1000):
            pos = posin -1
        score = float(line[11])
        if score > 500:
            score = 500
        score = score /500

        df = dfmod.query('chrom == "' + chrom + '" & pos == ' + str(pos))
        if not df.empty:

            mod = df["mod"].values[0]
            #print(pos,mod,score,fmer)
            if fmer in ndict:
                ndict[fmer].add(mod +" "+chrom +" "+str(pos))
            else:
                se = {mod+" "+chrom +" "+str(pos)}
                ndict[fmer] = se

            if fmer in nuclist:
                print(mod,line)

for kv in ndict:

    if len(ndict[kv]) > 1 :
        print(kv,ndict[kv])