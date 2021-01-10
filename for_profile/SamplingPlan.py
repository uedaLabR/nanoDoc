import numpy as np
from Bio import SeqIO

def depth(readrange,n):

    a = readrange[readrange[:,0] <= n]
    b = a[a[:,1] >= n+5]
    # print(n)
    # print(b)
    print(n,len(b))
    return len(b)

fivemerDict = {}

indexf = "C:\\Users\\hirok\\Desktop\\covid19\\coverageCalc\\indexKoraIVT.txt"
readrange = np.loadtxt(indexf, delimiter=',',dtype='int16', usecols=[2,3])
print(readrange)

records = SeqIO.parse("C:\\Users\\hirok\\Desktop\\covid19\\data\\Korea2.txt", 'fasta')
for record in records:
  print(record.name)
  print(record.seq)
  seq = str(record.seq)
  for n in range(5, len(seq)):
      fmer = seq[n - 5:n]

      if fivemerDict.get(fmer) == None:
          fivemerDict[fmer] = depth(readrange,n)
      else:
          fivemerDict[fmer] = fivemerDict[fmer] + depth(readrange,n)


# records = SeqIO.parse("C:\\Users\\hirok\\Desktop\\covid19\\data\\Curlcake.fa", 'fasta')
# for record in records:
#   print(record.name)
#   print(record.seq)
#   seq = str(record.seq)
#   for n in range(5,len(seq)):
#       fmer = seq[n-5:n]
#
#       if fivemerDict.get(fmer) == None:
#           fivemerDict[fmer] = 1
#       else:
#           fivemerDict[fmer] = fivemerDict[fmer]+1



  ret = sorted(fivemerDict.items(), key=lambda x: x[0])
  for r in ret:
      print(r)