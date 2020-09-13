import pyarrow.parquet as pq
import pandas as pd
from Bio import SeqIO
from sklearn.model_selection import train_test_split
import os


def loadParquet(chr, idxs):

    p = "/groups2/gac50430/nanopore/equalbinnedpq/koreaIVT/algined"
    if chr.startswith("cc6m"):
        p = "/groups2/gac50430/nanopore/equalbinnedpq/m6aIVT/algined"

    cnt = 0
    for idx in idxs:
        idx = idx.replace('}', '')
        idx = idx.replace('{', '')
        idx = idx.replace(' ', '')
        pp = p + idx + ".pq"
        table = pq.read_table(pp)
        df = table.to_pandas()
        if cnt == 0:
            totaldf = df
        else:
            totaldf = pd.concat([totaldf, df], axis=0)
        cnt = cnt + 1

    return totaldf


#     table2 = pq.read_table('example.parquet')

def loadRef(chr):
    record = None
    if chr.startswith("cc6m"):

        records = SeqIO.parse("/groups2/gac50430/nanopore/reference/Curlcake.fa", 'fasta')
        for r in records:
            print(r.name, chr, chr == r.name)
            if (chr == r.name):
                record = r
                break

    else:

        records = SeqIO.parse("/groups2/gac50430/nanopore/reference/Cov2_Korea.fa", 'fasta')
        for r in records:
            record = r
            break
    return record


def getData(df1, reference, position, samplenum):
    df1 = df1[(df1.mapped_start < position) & (df1.mapped_end > position + 5)]
    train, test = train_test_split(df1, test_size=samplenum)
    # print(test)
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


def mkpq(planf,ref,pq, outdir):

    f = open(planf)
    line = True
    parquetLoaded = False
    matrix = None
    reference = None
    cnt = 0
    while line:

        line = f.readline().rstrip('\n')
        data = line.split(",")
        if len(data) < 2:
            break
        if parquetLoaded == False:
            matrix = loadParquet(data[0], data[5:])
            reference = loadRef(data[0])
            parquetLoaded = True
        #         print(data)
        #        position = int(data[1])+1
        position = int(data[1])
        samplenum = int(data[3])
        fivemer = data[4]
        if not os.path.exists(outdir + fivemer):
            os.makedirs(outdir + fivemer)
        chr = data[0]
        if chr.startswith("hCoV"):
            chr = "hCoV-19_South"
        outpath = outdir + fivemer + "/" + chr + "_" + data[1] + ".pq"
        if os.path.exists(outpath):
            os.remove(outpath)
        data = getData(matrix, reference, position, samplenum)
        df = pd.DataFrame(data,
                          columns=['nucb4After', 'mapped_chrom', 'position', 'mapped_start', 'mapped_end', 'signal',
                                   'originalsize'])
        df.to_parquet(outpath)
        cnt = cnt + 1

    f.close


import sys

if __name__ == "__main__":
    outdir = "/groups2/gac50430/nanopore/dataset4DL/fivemer"
    args = sys.argv
    path = args[1]
    main(path, outdir)
    print(path)

