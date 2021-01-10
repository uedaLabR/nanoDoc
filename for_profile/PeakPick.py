import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
from matplotlib_venn import venn3
from matplotlib_venn import venn2
from Bio import SeqIO

def peakPick(df,snp):

   data = []
   for n in range(len(df)):

       pos = df.iat[n, 0]
       value = df.iat[n, 1]
       if value > 500:
           value = 500
       value = value / 500
       if value > 0.95 and pos >= 29751:
           continue

       lastidx = len(data) - 1
       if value > 0.05:

           if pos in snp or pos-1 in snp or pos+1 in snp:
               print(pos,value)
               continue
           if lastidx<0:
               data.append((pos, value))
           else:
               last = data[lastidx]
               lastpos = last[0]
               if(pos-lastpos <= 2):
                   if (last[1]<value):
                        data[lastidx] = (pos,value)
               else :
                   data.append((pos,value))
   s = set()
   for tp in data:
       s.add(tp[0])
   return pd.DataFrame(data),s

def todict(df):
    dict = {}
    for n in range(len(df)):

        pos = df.iat[n, 0]
        value = df.iat[n, 1]
        if value > 500:
            value = 500
        value = value / 500
        dict[pos] = value
    return dict



wigKorea = "C:\\Users\\hirok\\Desktop\\covid19\\wig\\window\\Korea_window.wig"
wigEng = "C:\\Users\\hirok\\Desktop\\covid19\\wig\\window\\England_window.wig"
wigAus = "C:\\Users\\hirok\\Desktop\\covid19\\wig\\window\\Australia_window.wig"

wigKoreao = "C:\\Users\\hirok\\Desktop\\covid19\\wig\\KoreaPeak.tsv"
wigEngo = "C:\\Users\\hirok\\Desktop\\covid19\\wig\\EnglandPeak.tsv"
wigAuso = "C:\\Users\\hirok\\Desktop\\covid19\\wig\\Austraria.tsv"

dfKorea = pd.read_csv(wigKorea, sep='\t')
dfEng = pd.read_csv(wigEng, sep='\t')
dfAus = pd.read_csv(wigAus, sep='\t')

dfKoreaDict = todict(dfKorea)
dfEngDict = todict(dfEng)
dfAusDict = todict(dfAus)

koreasnp = []
engsnp = [4402, 5062,18488,23605,29596]
aussnp = [4402, 5062,8782,19065,22303,26144,28144,29750]

dfKorea,dfKoreaSet = peakPick(dfKorea,koreasnp)
dfEng,dfEngSet = peakPick(dfEng,engsnp)
dfAus,dfAusSet = peakPick(dfAus,aussnp)



print("korea =" + str(len(dfKorea)))
print("eng =" + str(len(dfEng)))
print("aus =" + str(len(dfAus)))

ke = dfKoreaSet.intersection(dfEngSet)
ka = dfKoreaSet.intersection(dfAusSet)
ea = dfEngSet.intersection(dfAusSet)

all = dfKoreaSet.intersection(dfEngSet).intersection(dfAusSet)

allunion = dfKoreaSet.union(dfEngSet).union(dfAusSet)

print(len(ke))
print(len(ka))
print(len(ea))

print(len(all))

kcount = len(dfKorea)
ecount = len(dfEng)
acount = len(dfAus)

kecount = len(ke)
kacount = len(ka)
eacount = len(ea)
allintercect = len(all)

from matplotlib_venn import venn3
venn3(subsets = (kcount-(kecount+kacount-allintercect),ecount-(kecount+eacount-allintercect),(kecount-allintercect),acount-(kacount+eacount-allintercect), (kacount-allintercect), (eacount-allintercect), allintercect), set_labels = ('Korea', 'England', 'Australia'))
plt.show()



ic = 0
kimpos = [6891,12788,13750,15948,17585,25357,25460,26058,26297,26731,27106,27165,27268,27305,27486,27948,28591,28612,28653,28669,28701,28789,28805,28860,28879,28898,28931,28949,28958,29016,29041,29088,29127,29155,29170,29298,29313,29378,29408,29444,29776]
for n in dfKoreaSet:
    if n in kimpos or n-1 in kimpos or  n+1 in kimpos or  n+2 in kimpos or n-2 in kimpos:
        ic = ic +1

print(len(kimpos))
print("ic="+str(ic))

venn2(subsets=(kcount-ic,len(kimpos)-ic, ic),set_labels = ('nanoDOC', 'Kim et.al.'))
plt.show()

print(allunion)

records = SeqIO.parse("C:\\Users\\hirok\\Desktop\\covid19\\data\\Korea2.txt", 'fasta')
kseq = ""
for record in records:
  print(record.name)
  kseq = record.seq

counter={}
for pos in allunion:
    idx = pos -1
    nuc = kseq[idx]
    if nuc in counter:
        c = counter[nuc]
        c = c+1
        counter[nuc] = c
    else:
        counter[nuc] = 1

print(counter)

for pos in allunion:
    subseq = kseq[pos-5:pos+4]
    print(subseq)

poslist = []
for pos in allunion:

    poslist.append(pos)
poslist = sorted(poslist)


for pos in poslist:

    korea = "-"
    england = "-"
    austraria = "-"

    if pos in dfKoreaDict:
        korea = '{:.3f}'.format(dfKoreaDict[pos])
    if pos in dfEngDict:
        england = '{:.3f}'.format(dfEngDict[pos])
    if pos in dfAusDict:
        austraria = '{:.3f}'.format(dfAusDict[pos])

    print(pos,korea,england,austraria)

path = "C:\\Users\\hirok\\Desktop\\covid19\\wig\\abb5813.wig"
f = open(path)
lines = f.readlines()
pos2 = set()
for l in lines:
    if l.startswith("v"):
        continue
    print(l)
    pos2.add(int(l.split("\t")[0]))
f.close()

ic = 0
for n in pos2:
    if (n in poslist) or (n-1 in poslist) or  (n+1 in poslist) or  (n+2 in poslist) or (n-2 in poslist):
        ic = ic +1
        print(n)

venn2(subsets=(len(poslist)-ic,len(pos2), ic),set_labels = ('nanoDOC','mismatch','ic'))
plt.show()

print(len(poslist))

for pos in poslist:
    nuc = kseq[pos]

    korea = "-"
    england = "-"
    austraria = "-"

    if pos in dfKoreaDict:
        korea = '{:.3f}'.format(dfKoreaDict[pos])
    if pos in dfEngDict:
        england = '{:.3f}'.format(dfEngDict[pos])
    if pos in dfAusDict:
        austraria = '{:.3f}'.format(dfAusDict[pos])

    print(str(pos) +"\t" + nuc +"\t" + korea +"\t" + england +"\t" + austraria)

# dfKorea.to_csv(wigKoreao, sep='\t',index=False)
# dfEng.to_csv(wigEngo, sep='\t',index=False)
# dfAus.to_csv(wigAuso, sep='\t',index=False)