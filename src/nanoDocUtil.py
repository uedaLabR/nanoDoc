import glob
from sklearn.model_selection import train_test_split
import pyarrow.parquet as pq
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import random
import tensorflow as tf # add
import numpy as np
import cnn_network
import itertools
import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import LocalOutlierFactor
from sklearn import metrics
from sklearn.preprocessing import MinMaxScaler
from Bio import SeqIO
import faiss
from numpy import mean, absolute
import math
import gc
from scipy.optimize import curve_fit

DATA_LENGTH_UNIT = 60
DATA_LENGTH = DATA_LENGTH_UNIT*5+20

class PqReader:
    
    def __init__(self,indexfile,minreadlen):
        
        tb = []
        self.indexfile = indexfile
        self.minreadlen = minreadlen
        indexfile = indexfile+"/index.txt"
        #print(indexfile)
        f = open(indexfile)
        idxno =0
        rowcnt = 0
        chrb4 = ""
        for s_line in f:
            s_line = s_line.split(",")
            
            
            if chrb4 != s_line[1]:
                idxno = idxno+1
                rowcnt = 0
                
            chrb4 =  s_line[1]
            
            rowcnt = rowcnt+1
            if(rowcnt%4000==0):
                idxno = idxno+1                
                
            tb.append((idxno,s_line[1],int(s_line[2]),int(s_line[3])))    
            
        f.close()
        self.df = pd.DataFrame(tb,columns=('idx','chr','start','end'))
        #print(self.df)
        self.currentPqs = set()
    
    def loadpq(self,set0,chr,pos,strand):
        
        ss = sorted(set0)
        cnt = 0
        if len(self.currentPqs) > 0 and self.currentPqs <= set0:
            ss = set0-self.currentPqs
            cnt = len(self.currentPqs)
            totaldf = self.pq
        
 
        for s in ss:
           
            s = self.indexfile+"/algined"+str(s)+".pq"
            print(s)
            
            table = pq.read_table(s, columns=['nucb4After','mapped_chrom','mapped_strand','position','mapped_start','mapped_end','signal','originalsize'])
            df = table.to_pandas()  
            if strand == "+":
                df = df.query('(mapped_end - 40) >= '+str(pos)+  ' & mapped_chrom=="'+chr+'" & mapped_strand =="'+strand +'" & (mapped_end - mapped_start) > ' +str(self.minreadlen) )
            else:
                df = df.query('(mapped_start +40) <= '+str(pos)+  ' & mapped_chrom=="'+chr+'" & mapped_strand =="'+strand +'" & (mapped_end - mapped_start) > ' +str(self.minreadlen) )
            gc.collect()
            if cnt == 0:
                totaldf = df
            else:
                totaldf = pd.concat([totaldf, df], axis=0)     
                
                
            cnt = cnt +1
        
        return totaldf
    
    def freeUnsued(self,chr,pos,strand):
        
            if strand == "+":
                df = df.query('(mapped_end - 40) >= '+str(pos)+  ' & mapped_chrom=="'+chr+'" & mapped_strand =="'+strand +'" & (mapped_end - mapped_start) > ' +str(self.minreadlen) )
            else:
                df = df.query('(mapped_start +40) <= '+str(pos)+  ' & mapped_chrom=="'+chr+'" & mapped_strand =="'+strand +'" & (mapped_end - mapped_start) > ' +str(self.minreadlen) )
   
    def checkUpdate(self,chr,pos,strand):
        
        df = self.df.query('start <= '+str(pos) +' & end >= '+str(pos+5)+  ' & chr == "'+chr+'"')
#         print(df)
        df = df['idx'].drop_duplicates()
        s = set(df)
        if(s!= self.currentPqs or self.pq.empty):
            self.pq =self.loadpq(s,chr,pos,strand)        
            self.currentPqs =s
        
        #if pos%10 ==0:    
        #    self.freeUnsued(chr,pos,strand)
        
#             print(pos,s)
            
    def getData(self,chr,strand,pos,maxtake):
        
        unitwidth = DATA_LENGTH_UNIT
        self.checkUpdate(chr,pos,strand)
        pqlen = len(self.pq)    
        _pq = self.pq
        #print("pqlen=",len(_pq))    
        data = []    
        takecnt = 0
        for index, row in _pq.iterrows():

            mapped_chrom = row['mapped_chrom']
            mapped_start = row['mapped_start']
            mapped_end = row['mapped_end']

            if strand == "-":

                relpos = mapped_end - pos     
                if relpos < 0:
                    continue
                if pos-6 < mapped_start:
                    continue            


                signal = row['signal'][relpos*unitwidth:(relpos+5)*unitwidth]
                originalsize = row['originalsize'][relpos:relpos+5]        
            

                if(len(signal) == unitwidth*5):
                    data.append((signal,originalsize))
                    takecnt = takecnt+1

            else:
		
                relpos = pos-mapped_start     
                if relpos < 0:
                    continue
                if pos+5 > mapped_end:
                    continue            


                signal = row['signal'][relpos*unitwidth:(relpos+5)*unitwidth]
                originalsize = row['originalsize'][relpos:relpos+5]        
            

                if(len(signal) == unitwidth*5):
                    data.append((signal,originalsize))
                    takecnt = takecnt+1


            
            if takecnt == maxtake:
                break
        
        return data,len(data)

from Bio.Seq import Seq
def getSeq(ref,chrtgt,start,end,strand):

    records = SeqIO.parse(ref, 'fasta')
    for r in records:
        if r.id == chrtgt:
            record = r
    seq = record.seq[start:end]
    if strand =="-":
        seq2 = Seq(str(seq))
        seq = seq2.reverse_complement()
    return seq

def getFirstChrom(ref):
    
    records = SeqIO.parse(ref, 'fasta')
    for r in records:
        return  r.id
    return ""

def readParam(paramf):
    
    f = open(paramf)
    for line in f:
        if line.startswith("a="):
            a = float(line.replace("a=",""))
        if line.startswith("b="):
            b = float(line.replace("b=",""))            
        if line.startswith("uplimit="):
            uplimit =int(line.replace("uplimit=",""))
        if line.startswith("takeparcentile="):
            takeparcentile = float(line.replace("takeparcentile=",""))            
       
    f.close()
    return a,b,uplimit,takeparcentile

def normalize(v, axis=-1, order=2):
    l2 = np.linalg.norm(v, ord = order, axis=axis, keepdims=True)
    l2[l2==0] = 1
    return v/l2

def histG(dist):
    plt.xlim(0, 3000)
    plt.hist(dist, range=(10,3000),bins=2900)
    plt.show()
    
def reducesize(data,size):
    
    if len(data) <= size:
        return data,len(data)
    
    a_train, a_test = train_test_split(data, test_size=size)
    return a_test,len(a_test)



def mad(a, axis=None):
    """
    Compute *Median Absolute Deviation* of an array along given axis.
    """

    # Median along given axis, but *keeping* the reduced axis so that
    # result can still broadcast against a.
    med = np.median(a, axis=axis, keepdims=True)
    mad = np.median(np.absolute(a - med), axis=axis)  # MAD along given axis

    return mad
    
def getFormat(lt):
    
    flattensignal = []        
    train_x = []
    cnt = 0
    for row in lt:

        cnt = cnt+1
        signal = np.array(list(row[0]))
        
#        for signal mean 
        flattensignal.extend(signal[DATA_LENGTH_UNIT*2:DATA_LENGTH_UNIT*3])
        
        signal = zeropadding10(signal)
        signal = np.array(signal)
        signal = signal.astype('float32')/255.  

        originalsize = extendAry(row[1])
        originalsize = zeropadding10(originalsize)
        
        

        if cnt<=0:
            
 #           print(len(signal))
            plt.plot(signal)
            plt.show()
            print(len(originalsize))
            plt.plot(originalsize)
 #           plt.show()
        
        train_x.extend(signal)
        train_x.extend(originalsize)
     
    train_x = np.reshape(train_x, (-1, DATA_LENGTH,2))    
    med = np.median(flattensignal)
    _mad = mad(flattensignal)
    
    return train_x,med,_mad



def getDist(a,b):
    
    index = faiss.IndexFlatL2(a.shape[1])   # build the index
    index.add(a)                  # add vectors to the index
    dists, result = index.search(b, k=10)     # actual search
    
    index = faiss.IndexFlatL2(b.shape[1])   # build the index
    index.add(b)                  # add vectors to the index
    dists2, result2 = index.search(a, k=10)     # actual search
    
    dists = np.mean(dists, axis=1)
    dists2 = np.mean(dists2, axis=1)
    
    dists = np.clip(dists, 0, 3000)
    dists2 = np.clip(dists2, 0, 3000)
    
    return dists,dists2

def acc(dist):
    
    hist,idx = np.histogram(dist, range=(10,3000),bins=2990)    
       
    acc = np.add.accumulate(hist)
    maxv = np.max(acc)
    acc = acc/maxv
    totalcnt = np.sum(hist)
      
    return acc,totalcnt


def extendAry(ar):
    
#     print(ar)
    ar_ret = []
    for i in ar:
        
        if i >= 500:
            i = 500
        i = i/500
        ar = [i]*DATA_LENGTH_UNIT
        ar_ret.extend(ar)        
    
#     print(ar_ret)
    return addnoise(ar_ret)    


def zeropadding10(sig):

    b4 =  addnoise([0]*10)
    after =  addnoise([0]*10)
    sig2 = np.append(b4,sig)
    return np.append(sig2,after)

def addnoise(npa):
    
    noise_unit= 0.005
    noise = np.random.normal(-noise_unit,noise_unit,len(npa))
    return npa + noise

def addnoise_old(ls):
    ls = list(map(lambda x: x + random.uniform(-1.0, +1.0) *0.005, ls))
    return ls