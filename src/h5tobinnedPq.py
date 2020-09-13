import h5py
from statsmodels import robust
import numpy as np
import random
from scipy import interpolate
import pyarrow as pa
import pyarrow.parquet
import pandas as pd
# import mappy as mp
import tqdm
from functools import cmp_to_key
from multiprocessing import Pool
from tqdm import tqdm

def equalize(bindata):
    binsize = 60
    l = len(bindata)

    if l == 0:
        return addnoise([0] * binsize)
    elif l == 1:
        return addnoise([int(bindata[0])] * binsize)
    elif l == 2:
        datal1 = [int(bindata[0])] * int(binsize / 2)
        datal2 = [int(bindata[1])] * int(binsize / 2)
        datal1.extend(datal2)
        return addnoise(datal1)
    elif l == binsize:
        return bindata
    elif l < binsize:
        return addnoise(unsampling(bindata, binsize))
    else:
        return downsampling(bindata, binsize)


def addnoise(npa):
    
    npa = np.array(npa)
    noise_unit= 0.005
    noise = np.random.normal(-noise_unit,noise_unit,len(npa))
    return npa + noise


def unsampling(bindata, binsize):
    lsize = len(bindata)

    x_observed = []
    y_observed = []

    unit = binsize / lsize
    for i in range(0, lsize):
        x_observed.append(i * unit)
        y_observed.append(bindata[i])

    method = lambda x, y: interpolate.interp1d(x, y, kind="quadratic")
    try:
        fitted_curve = method(x_observed, y_observed)
    except ValueError:
        print(x_observed, y_observed)
    x_latent = np.linspace(min(x_observed), max(x_observed), binsize)
    ls = fitted_curve(x_latent)
    #    print("ups",bindata,ls)
    ls = np.array(ls)
    ls = ls.astype(np.int8)
    return ls


def downsampling(bindata, binsize):
    lsize = len(bindata)
    indexs = random.sample(range(0, lsize), binsize)
    ls = []
    for i in range(0, lsize):
        if i in indexs:
            ls.append(bindata[i])
    return ls


def toByteRange(x):
    y = int(x) + 128
    if y > 255:
        #        print(y)
        y = 255
    if y < 0:
        #        print(y)
        y = 0
    return y


def indexFiles(path):

    unit = 4000
    cnt = 0
    idxcnt = 0
    ref = ""
    lastref = ""
    blist = []
    indexlist = []
    with open(path) as f:
        for s_line in f:
            ref = s_line.split(",")[1]
            if ref != lastref or cnt % 4000 == 0:
                if len(indexlist) > 0:
                    blist.append(indexlist)
                indexlist = []
                cnt = 0
            indexlist.append(s_line)
            cnt = cnt + 1
            lastref = ref
        if len(indexlist) > 0:
            blist.append(indexlist)

    return blist


def makeParquet(path, pathout, refg, thread):

    indexlist = indexFiles(path)
    blist = []
    idx = 1
    for d in indexlist:
        pathoutpq = pathout +"/"+ str(idx) + ".pq"
        tp = (d,pathoutpq,refg)
        blist.append(tp)
        idx = idx + 1

    size = len(blist)
    print(size)
    
#    for cont in blist:
#        _makeParquet(cont)
    with Pool(thread) as p:
        p.map(_makeParquet, blist)
#        tqdm(imap,total=size)


def _makeParquet(tp):

    (indexcontent, pathout, refg) = tp
    # aligner = mp.Aligner(refg, preset="map-ont")

    wholedata = []

    cnt_1 = 0
    for s_line in indexcontent:
        fdata = s_line.split(",")
        onerecord = getRecord(fdata)
        if onerecord != None:
           wholedata.append(onerecord)
           cnt_1 = cnt_1+1
           print(cnt_1)
        # (mapped_chrom, mapped_strand, mapped_start, mapped_end, clipped_bases_start, clipped_bases_end, fastq, bytedata,binsizelist)
    df = pd.DataFrame(wholedata,
                   columns=['mapped_chrom', 'mapped_strand', 'mapped_start', 'mapped_end', 'clipped_bases_start',
                                'clipped_bases_end', 'fastq', 'ctg', 'st', 'ed', 'cigar', 'md', 'signal',
                                 'originalsize'])
    df.to_parquet(pathout)


def normalise(signal):
    med = np.median(signal)  # take median
    mad = robust.mad(signal)  # take mad
    signal = ((signal - med) / mad)  # normalize
    signal = (((signal - 1) / 4) * 128)  # fit into byte
    signal = signal.astype(np.int8)  # type conversion
    return signal


def jointByteArray(binnedlist):
    elist = []
    for ll in binnedlist:
        ll = list(map(toByteRange, ll))
        elist.extend(ll)
    return bytes(elist)


def getFastQ(f, key):
    if key in f:
        group = f[key]
        datalist = list(group.values())
        fastq = datalist[0]
        fastq = fastq[()]
        return fastq.decode('utf-8')
    return None

def getRecord(fdata):
    binsize = 60

    with h5py.File(fdata[0], 'r') as f:

        group = f['/Analyses/RawGenomeCorrected_000/BaseCalled_template']
        datalist = list(group.values())
        alignment = datalist[0]
#        attr = list(alignment.attrs.items())
#        print(attr)
        event = datalist[1]
        allevent = event[()]
        raw_start = event.attrs.get('read_start_rel_to_raw')
        dlen = len(allevent)-1
        end = allevent[dlen][2]

        group = f['/Raw/Reads/']
        readid = list(group.values())[0].name
        signalkey = readid + "/Signal"
        signal = f[signalkey][()]

        signal = signal[::-1] #reverse for rna
        signal = signal[0:end] #upto tombo segmantation
        signal = normalise(signal)

        # seq = []
        # for n in range(len(event.value)):
        #     seq.append(event.value[n][4])
        # seq = np.array(seq, dtype=np.unicode)
        # a_str = ''.join(str(x) for x in seq)

        mapped_chrom = alignment.attrs.get('mapped_chrom')
        mapped_strand = alignment.attrs.get('mapped_strand')
        mapped_start = alignment.attrs.get('mapped_start')
        mapped_end = alignment.attrs.get('mapped_end')

        clipped_bases_start = alignment.attrs.get('clipped_bases_start')
        clipped_bases_end = alignment.attrs.get('clipped_bases_end')

        binnedlist = []
        binsizelist = []
        for n in range(len(allevent)):
            #           print(event.value[n][0], event.value[n][1], event.value[n][2], event.value[n][3], event.value[n][4])
            start = raw_start + allevent[n][2]
            size = allevent[n][3]
            binnedlist.append(signal[start:start + size])
            binsizelist.append(size)

        binnedlist = list(map(equalize, binnedlist))
        bytedata = jointByteArray(binnedlist)

        fastq = getFastQ(f,'/Analyses/Basecall_1D_001/BaseCalled_template')
        if fastq == None:
            fastq = getFastQ(f,'/Analyses/Basecall_1D_000/BaseCalled_template')

        ctg,st,ed,cigar,md = "",0,0,"",""


        tp = (mapped_chrom,mapped_strand,mapped_start,mapped_end,clipped_bases_start,clipped_bases_end,fastq,ctg,st,ed,cigar,md,bytedata,binsizelist)
        return tp

    return None


def _getRecord(fdata):
    binsize = 60

    try:

        f = h5py.File(fdata[0], 'r') 
        group = f['/Analyses/RawGenomeCorrected_000/BaseCalled_template']
        fastq = getFastQ(f,'/Analyses/Basecall_1D_001/BaseCalled_template')
        if fastq == None:
            fastq = getFastQ(f,'/Analyses/Basecall_1D_000/BaseCalled_template')

        datalist = list(group.values())
        alignment = datalist[0]

        mapped_chrom = alignment.attrs.get('mapped_chrom')
        mapped_strand = alignment.attrs.get('mapped_strand')
        mapped_start = alignment.attrs.get('mapped_start')
        mapped_end = alignment.attrs.get('mapped_end')

        clipped_bases_start = alignment.attrs.get('clipped_bases_start')
        clipped_bases_end = alignment.attrs.get('clipped_bases_end')


        event = datalist[1]
        allevent = event[()]
        raw_start = event.attrs.get('read_start_rel_to_raw')
        dlen = len(allevent)-1
        end = allevent[dlen][2]

        group = f['/Raw/Reads/']
        readid = list(group.values())[0].name
        signalkey = readid + "/Signal"
        signal = f[signalkey][()]

        signal = signal[::-1] #reverse for rna
        signal = signal[0:end] #upto tombo segmantation
        signal = normalise(signal)

        size = len(allevent)
        binnedlist = []
        binsizelist = []
        for n in range(size):

#            print(event.value[n][0], event.value[n][1], event.value[n][2], event.value[n][3], event.value[n][4])
            start = raw_start +allevent[n][2]
            size = allevent[n][3]
#            print(n,size,len(signal),start,size)
            binnedlist.append(signal[start:start + size])
            binsizelist.append(size)

        f.close()

        binnedlist = list(map(equalize, binnedlist))
        bytedata = jointByteArray(binnedlist)
        ctg, st, ed, cigar, md = "", 0, 0, "", ""
        bytedata =[]
        tp = (mapped_chrom,mapped_strand,mapped_start,mapped_end,clipped_bases_start,clipped_bases_end,fastq,ctg,st,ed,cigar,md,bytedata,binsizelist)
        print(mapped_chrom,mapped_strand,mapped_start,mapped_end)

        return tp



    except OSError:
        pass
    except KeyError:
        pass
    print("None")
    return None




import sys

# if __name__ == "__main__":
    #    path = 'C:\\Users\\hirok\\Desktop\\covid19\\kimivt\\index.txt'
    #    pathout = 'C:\\Users\\hirok\\Desktop\\covid19\\kimivt\\kinivt'
    #    index = str(1)
    # args = sys.argv
    # main(args[1], args[2], args[3], args[4])