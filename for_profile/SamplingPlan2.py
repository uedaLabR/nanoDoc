import numpy as np
from Bio import SeqIO

class Counter:

    cnt = 1
    s = set()
    l = []

    def __init__(self, v,name,n,depth,idxs):
       self.s = set()
       self.l = []
       self.s.add(v)
       self.l.append((name,n,depth,idxs))

    def inc(self,v,name,n,depth,idxs):
       self.cnt = self.cnt + 1
       self.s.add(v)
       self.l.append((name, n, depth,idxs))

    def getS(self):
        return self.s
    def getCnt(self):
        return self.cnt
    def getList(self):
        return self.l

def depthG(readrange,n):

    a = readrange[readrange[:,0] <= n]
    b = a[a[:,1] >= n+5]
    return b

def loadDepthDB(name):

    if name == "hCoV-19/South":
        indexf = "C:\\Users\\hirok\\Desktop\\covid19\\coverageCalc\\indexKoraIVT_2.txt"
        readrange = np.loadtxt(indexf, delimiter=',',dtype='int16', usecols=[2,3,4])
        return readrange
    else:

        xlist=[]
#        indexf = "C:\\Users\\hirok\\Desktop\\covid19\\coverageCalc\\indexM6a_2.txt"
        indexf = "C:\\Users\\hirok\\Desktop\\covid19\\coverageCalc\\indexM6Amod_2.txt"
        f = open(indexf)
        line = True
        while line:
            line = f.readline()
            data = line.split(",")
            if(len(data)<3):
                break
            if data[1] == name:
                xlist.append((data[2],data[3],data[4]))

        f.close
        return np.array(xlist,dtype='int16')


fivemerDict={}
recordL=[]
# records = SeqIO.parse("C:\\Users\\hirok\\Desktop\\covid19\\data\\Korea2.txt", 'fasta')
# for record in records:
#     recordL.append(record)
records = SeqIO.parse("C:\\Users\\hirok\\Desktop\\covid19\\data\\Curlcake.fa", 'fasta')
for record in records:
    recordL.append(record)

for record in recordL:
   print(record.name)
   print(record.seq)
   readrange = loadDepthDB(record.name)
   seq = str(record.seq)

   for n in range(10, len(seq)-16):

      fmer = seq[n:n+5]
      b4 = seq[n-1]
      after = seq[n+6]
      v = b4+after
      m1 = depthG(readrange,n)
      depth = len(m1)

      idxs = set(m1[:,2].tolist())
#      print(idxs)
      print(record.name,n,depth)
      if fivemerDict.get(fmer) == None:
          fivemerDict[fmer] = Counter(v,record.name,n,depth,idxs)
      else:
          fivemerDict[fmer].inc(v,record.name,n,depth,idxs)


posList=[]
samplenum = 10000
ret = sorted(fivemerDict.items(), key=lambda x: x[0])
for r in ret:
    ls = r[1].getList()
    llen = len(ls)
    divn = samplenum//llen
    modn = samplenum%llen
    #
    dlist=[]
    cnt = 0
    for d in ls:
        if cnt ==0:
            dlist.append((d[0],d[1],d[2],divn+modn,r[0],d[3]))
        else:
            dlist.append((d[0], d[1], d[2], divn,r[0],d[3]))
        cnt = cnt+1

    posList.extend(dlist)
    print(r[0],r[1].getCnt(),len(r[1].getS()),dlist)

from operator import itemgetter, attrgetter
posList = sorted(posList, key=itemgetter(0,1))
#path_w="C:\\Users\\hirok\\Desktop\\covid19\\coverageCalc\\readPlan\\g_"
path_w="C:\\Users\\hirok\\Desktop\\covid19\\coverageCalc\\readPlanm6A\\g_"
sidx = 0
sb4 = None

f = None
for d in posList:
    s = d[5]
    if s != sb4:
       sidx = sidx+1
       if f!=None:
           f.close()
           f = open(path_w+str(sidx)+".txt", mode='w')
       else:
           f = open(path_w + str(sidx) + ".txt", mode='w')
    f.write(",".join(map(str, d))+"\n")
    sb4 = s
f.close()