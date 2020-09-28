import pyarrow.parquet as pq
import pandas as pd
from Bio import SeqIO
from sklearn.model_selection import train_test_split
import os
import shutil
import nanoDocUtil

def loadParquet(p,chr, idxs):


    cnt = 0
    for idx in idxs:
        idx = idx.replace('}', '')
        idx = idx.replace('{', '')
        idx = idx.replace(' ', '')
        pp = p +"/algined" +idx + ".pq"
        table = pq.read_table(pp)
        df = table.to_pandas()
        if cnt == 0:
            totaldf = df
        else:
            totaldf = pd.concat([totaldf, df], axis=0)
        cnt = cnt + 1
#    print(totaldf)
    return totaldf


#     table2 = pq.read_table('example.parquet')

def loadRef(ref,chr):
    record = None
    records = SeqIO.parse(ref, 'fasta')
    for r in records:
       record = r
       if r.id == chr:
          break
    return record


def getData(df1, reference, position, samplenum):

    df1 = df1[(df1.mapped_start < position) & (df1.mapped_end > position + 5)]
    train, test = train_test_split(df1, test_size=samplenum)
    print(test)
    cnt = 0
    unitwidth = 60

    #         df = pd.DataFrame(data, columns=['nucb4After','mapped_chrom','position','mapped_start','mapped_end','signal','originalsize'])
    data = []
    for index, row in test.iterrows():
        mapped_chrom = row['mapped_chrom']
        mapped_start = row['mapped_start']
        mapped_end = row['mapped_end']
        nucb4After = reference.seq[position - 1] + reference.seq[position + 5]

        relpos = position - mapped_start
        signal = row['signal'][relpos * unitwidth:(relpos + 5) * unitwidth]

        #         if cnt <= 1:
        #             isignal = list(signal)
        #             plt.figure(figsize=(20, 5))
        #             plt.plot(isignal)

        originalsize = row['originalsize'][relpos:relpos + 5]
        cnt = cnt + 1
        data.append((nucb4After, mapped_chrom, position, mapped_start, mapped_end, signal, originalsize))

    return data


def mkpq(planf,ref,p,outdir):

    f = open(planf)
    line = True
    parquetLoaded = False
    matrix = None
    reference = None
    lastreference = None
    cnt = 0
    minreadlen = 100
    refpr = nanoDocUtil.PqReader(p, minreadlen)
    uplimit = 100000

    while line:

        line = f.readline().rstrip('\n')
        data = line.split(",")
        if len(data) < 2:
            break

        chr = data[0]
        position = int(data[1])
        samplenum = int(data[3])
        fivemer = data[4]

        if chr != lastreference :
            reference = loadRef(ref,chr)
	
        outpath = outdir + "/" + fivemer + "/" + chr + "_" + data[1] + ".pq"
        if os.path.exists(outpath):
            continue

        strand = '+'
        rawdatas, cnt = refpr.getData(chr, strand, position, uplimit)
        
        if not os.path.exists(outdir  + "/" + fivemer):
            os.makedirs(outdir + "/" + fivemer)



        df = pd.DataFrame(rawdatas,
                          columns=['signal','originalsize'])
        df.to_parquet(outpath)
        cnt = cnt + 1
        lastreference = chr

    f.close
    
    #merge files
    files = []    
    for x in os.listdir(outdir):  
        if os.path.isdir(outdir + "/" + x):
            files.append(x)
    
    totaldf =None
    for dir in files: 
        for each in os.listdir(outdir+ "/" +dir):  
            s = outdir+ "/" +dir+ "/"+each
            table = pq.read_table(s, columns=['signal','originalsize'])
            df = table.to_pandas()  
            if cnt == 0:
                totaldf = df
            else:
                totaldf = pd.concat([totaldf, df], axis=0)  
        outpath = outdir+ "/" +dir +".pq"
        print(outpath)
        totaldf.to_parquet(outpath)
        shutil.rmtree(outdir+ "/" +dir)

import sys

if __name__ == "__main__":
    outdir = "/groups2/gac50430/nanopore/dataset4DL/fivemer"
    args = sys.argv
    path = args[1]
    main(path, outdir)
    print(path)

