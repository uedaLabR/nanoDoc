# from Bio import SeqIO
#
# # seq = "GGTTAAGCGACTAAGCGTACACGGTGGATGCCCTGGCAGTCAGAGGCGATGAAGGACGTGCTAATCTGCGATAAGCGTCGGTAAGGTGATATGAACCGTTATAACCGGCGATTTCCGAATGGGGAAACCCAGTGTGTTTCGACACACTATCATTAACTGAATCCATAGGTTAATGAGGCGAACCGGGGGAACTGAAACATCTAAGTACCCCGAGGAAAAGAAATCAACCGAGATTCCCCCAGTAGCGGCGAGCGAACGGGGAGCAGCCCAGAGCCTGAATCAGTGTGTGTGTTAGTGGAAGCGTCTGGAAAGGCGCGCGATACAGGGTGACAGCCCCGTACACAAAAATGCACATGCTGTGAGCTCGATGAGTAGGGCGGGACACGTGGTATCCTGTCTGAATATGGGGGGACCATCCTCCAAGGCTAAATACTCCTGACTGACCGATAGTGAACCAGTACCGTGAGGGAAAGGCGAAAAGAACCCCGGCGAGGGGAGTGAAAAAGAACCTGAAACCGTGTACGTACAAGCAGTGGGAGCACGCTTAGGCGTGTGACTGCGTACCTTTTGTATAATGGGTCAGCGACTTATATTCTGTAGCAAGGTTAACCGAATAGGGGAGCCGAAGGGAAACCGAGTCTTAACTGGGCGTTAAGTTGCAGGGTATAGACCCGAAACCCGGTGATCTAGCCATGGGCAGGTTGAAGGTTGGGTAACACTAACTGGAGGACCGAACCGACTAATGTTGAAAAATTAGCGGATGACTTGTGGCTGGGGGTGAAAGGCCAATCAAACCGGGAGATAGCTGGTTCTCCCCGAAAGCTATTTAGGTAGCGCCTCGTGAATTCATCTCCGGGGGTAGAGCACTGTTTCGGCAAGGGGGTCATCCCGACTTACCAACCCGATGCAAACTGCGAATACCGGAGAATGTTATCACGGGAGACACACGGCGGGTGCTAACGTCCGTCGTGAAGAGGGAAACAACCCAGACCGCCAGCTAAGGTCCCAAAGTCATGGTTAAGTGGGAAACGATGTGGGAAGGCCCAGACAGCCAGGATGTTGGCTTAGAAGCAGCCATCATTTAAAGAAAGCGTAATAGCTCACTGGTCGAGTCGGCCTGCGCGGAAGATGTAACGGGGCTAAACCATGCACCGAAGCTGCGGCAGCGACGCTTATGCGTTGTTGGGTAGGGGAGCGTTCTGTAAGCCTGCGAAGGTGTGCTGTGAGGCATGCTGGAGGTATCAGAAGTGCGAATGCTGACATAAGTAACGATAAAGCGGGTGAAAAGCCCGCTCGCCGGAAGACCAAGGGTTCCTGTCCAACGTTAATCGGGGCAGGGTGAGTCGACCCCTAAGGCGAGGCCGAAAGGCGTAGTCGATGGGAAACAGGTTAATATTCCTGTACTTGGTGTTACTGCGAAGGGGGGACGGAGAAGGCTATGTTGGCCGGGCGACGGTTGTCCCGGTTTAAGCGTGTAGGCTGGTTTTCCAGGCAAATCCGGAAAATCAAGGCTGAGGCGTGATGACGAGGCACTACGGTGCTGAAGCAACAAATGCCCTGCTTCCAGGAAAAGCCTCTAAGCATCAGGTAACATCAAATCGTACCCCAAACCGACACAGGTGGTCAGGTAGAGAATACCAAGGCGCTTGAGAGAACTCGGGTGAAGGAACTAGGCAAAATGGTGCCGTAACTTCGGGAGAAGGCACGCTGATATGTAGGTGAGGTCCCTCGCGGATGGAGCTGAAATCAGTCGAAGATACCAGCTGGCTGCAACTGTTTATTAAAAACACAGCACTGTGCAAACACGAAAGTGGACGTATACGGTGTGACGCCTGCCCGGTGCCGGAAGGTTAATTGATGGGGTTAGCGCAAGCGAAGCTCTTGATCGAAGCCCCGGTAAACGGCGGCCGTAACTATAACGGTCCTAAGGTAGCGAAATTCCTTGTCGGGTAAGTTCCGACCTGCACGAATGGCGTAATGATGGCCAGGCTGTCTCCACCCGAGACTCAGTGAAATTGAACTCGCTGTGAAGATGCAGTGTACCCGCGGCAAGACGGAAAGACCCCGTGAACCTTTACTATAGCTTGACACTGAACATTGAGCCTTGATGTGTAGGATAGGTGGGAGGCTTTGAAGTGTGGACGCCAGTCTGCATGGAGCCGACCTTGAAATACCACCCTTTAATGTTTGATGTTCTAACGTTGACCCGTAATCCGGGTTGCGGACAGTGTCTGGTGGGTAGTTTGACTGGGGCGGTCTCCTCCTAAAGAGTAACGGAGGAGCACGAAGGTTGGCTAATCCTGGTCGGACATCAGGAGGTTAGTGCAATGGCATAAGCCAGCTTGACTGCGAGCGTGACGGCGCGAGCAGGTGCGAAAGCAGGTCATAGTGATCCGGTGGTTCTGAATGGAAGGGCCATCGCTCAACGGATAAAAGGTACTCCGGGGATAACAGGCTGATACCGCCCAAGAGTTCATATCGACGGCGGTGTTTGGCACCTCGATGTCGGCTCATCACATCCTGGGGCTGAAGTAGGTCCCAAGGGTATGGCTGTTCGCCATTTAAAGTGGTACGCGAGCTGGGTTTAGAACGTCGTGAGACAGTTCGGTCCCTATCTGCCGTGGGCGCTGGAGAACTGAGGGGGGCTGCTCCTAGTACGAGAGGACCGGAGTGGACGCATCACTGGTGTTCGGGTTGTCATGCCAATGGCACTGCCCGGTAGCTAAATGCGGAAGAGATAAGTGCTGAAAGCATCTAAGCACGAAACTTGCCCCGAGATGAGTTCTCCCTGACCCTTTAAGGGTCCTGAAGGAACGTTGAAGACGACGACGTTGATAGGCCGGGTGTGTAAGCGCAGCGATGCGTTGAGCTAACCGGTACTAATGAACCGTGAGGCTTAACCTT"
# #
# # print(seq[1618-3:1618+2])
# # print(seq[2030-3:2030+2])
# #
# # f = open("C:\\Users\\hirok\\Desktop\\covid19\\detection200730\\23S_4037519_4040423.txt")
# # l = f.readlines()
# # data = []
# # for line in l:
# #     nuc = line.split("\t")[1][2]
# #     data.append(nuc)
# # f.close()
# # s = "GGT" + "".join(data)
# # print(seq)
# # print(s)
# #
# # def readFasta(path):
# #     f = open(path)
# #     l = f.readlines()
# #     L = []
# #     for line in l:
# #         if not line.startswith(">"):
# #             L.append(line.replace('\n',''))
# #     f.close()
# #     s = ''.join(L)
# #     return s
# #
# # def subseq(seq,start,end):
# #
# #     seqr = ""
# #     cnt = 0
# #     for n in seq:
# #         if cnt >= start and cnt <= end:
# #             seqr = seqr+n
# #         cnt = cnt +1
# #
# #     return seqr
#
# records = SeqIO.parse("C:\\Users\\hirok\\Desktop\\covid19\\data\\Curlcake.fa", 'fasta')
# rec = ""
#
# for record in records:
#     print(record.name)
#     rec = record.seq[4168641:4171544]
#     rec2 = record.seq[4037519:4040423]
#
# print(rec)
# print(rec2)
#
# # rfile = readFasta("C:\\Users\\hirok\\Desktop\\covid19\\data\\NC000913.fa")
# # rec = rfile
# # #print(rec)
# # print(len(rec))
# seq23s = rec
# # print(4040423-4037519)
# # print(len(seq23s))
# # print(seq23s)
# # #ec23s	1618	m6A
# # #ec23s	2030	m6A
# print(seq23s[1616-5:1617+5])
# print(seq23s[2028:2029])
#
# seq23s = rec2
# # print(4040423-4037519)
# # print(len(seq23s))
# # print(seq23s)
# # #ec23s	1618	m6A
# # #ec23s	2030	m6A
# print(seq23s[1616-5:1617+5])
# print(seq23s[2029-5:2030+5])
#
# idx = 0
# for n in rec:
#     if n != rec2[idx]:
#         print(idx)
#     idx = idx +1
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
from sklearn.metrics import roc_curve
from sklearn.metrics import auc


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

def isDarch(s):
    if s[0] == 'C':
        return False
    if s[1] == 'C' or s[1] == 'U':
        return False
    if s[3] != 'C':
        return False
    if s[4] == 'G':
        return False
    return True

path = "C:\\Users\\hirok\\Desktop\\covid19\\detection200730\\m6A\\*.txt"
l = glob.glob(path)

holdertmp = []

nomodscore = []
holder = []
highest = {}

holderD = []
highestD = {}

for p in l:

    print(p)

    f = open(p)
    lines = f.readlines()
    posin = 3
    holder.append(("No Mod", 0))
    fmerbf = ""
    fmerafter = ""
    for cn  in range(1,len(lines)-1):

        line = lines[cn]

        line = line.split("\t")
        line[9] = line[9].replace('\n', '')
        # pos = int(line[0])
        pos = posin
        score = float(line[11])
        if score > 500:
            score = 500
        score = score /500
        if score < 0:
            score = 0
        fivemer = line[1].replace('T','U')

        fmerbf = lines[cn-1].split("\t")[1].replace('T','U')
        fmerafter = lines[cn+1].split("\t")[1].replace('T','U')
        if not ("A" in fivemer) and not ("A" in  fmerbf) and not ("A" in fmerafter) :

           print("test",fivemer,fmerbf)
           holder.append(("No Mod",score))
           nomodscore.append(score)

        elif fivemer[2:3]=="A":
            if(fivemer.count('A')==1):

                threemer = fivemer[1:4]
                holder.append((threemer,score))



    f.close()

c_array = np.percentile(nomodscore, q=[25,50,75])

disporder =['No Mod' ,'GAC','CAG', 'GAU'  ,'UAG', 'GAG', 'UAU', 'UAC','CAU','CAC']


#print(holder)
df = pd.DataFrame(holder,columns=['threemer', 'score'])
# plot_order = df.sort_values(by='score', ascending=False).fivemer.values
# print(plot_order)
u = df['threemer'].unique()
print(u)

plt.figure(figsize=(24, 12))
#ax = sns.violinplot(x=df["modification"], y =df["score"], inner=None, color="0.8", linewidth=1)
ax = sns.boxplot(x=df["threemer"], y =df["score"],color="0.95",linewidth=0.5 ,order =disporder )
ax = sns.swarmplot(x=df["threemer"], y =df["score"],palette="deep" ,order =disporder )
ax.axhline(c_array[0], ls='--',linewidth=0.5)
ax.axhline(c_array[1], ls='--',linewidth=0.5)
ax.axhline(c_array[2], ls='--',linewidth=0.5)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45,fontsize='15')

plt.show()
