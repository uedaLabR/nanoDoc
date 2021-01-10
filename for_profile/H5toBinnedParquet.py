import h5py
from statsmodels import robust
import numpy as np
import random
from scipy import interpolate
import pyarrow as pa
import pyarrow.parquet
import pandas as pd
import mappy as mp

def equalize(bindata):
    binsize = 60
    l = len(bindata)

    if l == 0:
        return addnoise([0] * binsize)
    elif l == 1:
        return addnoise([int(bindata[0])] * binsize)
    elif l == 2:
        datal1 = [int(bindata[0])] * int(binsize / 2 )
        datal2 = [int(bindata[1])] *  int(binsize / 2)
        datal1.extend(datal2)
        return addnoise(datal1)
    elif l == binsize:
        return bindata
    elif l < binsize:
        return addnoise(unsampling(bindata, binsize))
    else:
        return downsampling(bindata, binsize)


def addnoise(ls):
    ls = list(map(lambda x: x + random.uniform(-1.0, +1.0) * 5, ls))
    return ls

def unsampling(bindata, binsize):
    lsize = len(bindata)

    x_observed = []
    y_observed = []

    unit = binsize / lsize
    for i in range(0, lsize):
        x_observed.append(i * unit)
        y_observed.append(bindata[i])

    # "2次元"
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

    y = int(x) +128
    if y > 255:
#        print(y)
        y = 255
    if y < 0:
#        print(y)
        y = 0
    return y

def indexFiles(path,index):

    unit = 4000
    cnt = 0
    idxcnt = 0
    ref = ""
    lastref = ""
    indexlist = []
    with open(path) as f:
        for s_line in f:
            ref = s_line.split(",")[1]
            if ref != lastref or cnt % 4000 == 0:
                idxcnt = idxcnt + 1
                # print(idxcnt)

            if idxcnt == index:
                indexlist.append(s_line) #calculate only indexed files for distributed processing
            if idxcnt > index:
                break

            cnt = cnt + 1
            lastref = ref
    return indexlist

def main(path,pathout,refg,index):

    indexlist = indexFiles(path,int(index))
    aligner = mp.Aligner(refg,preset="map-ont",MD=True)

    wholedata = []
    for s_line in indexlist:
        fdata = s_line.split(",")
        onerecord = getRecord(fdata,aligner)
        if onerecord != None:
            wholedata.append(onerecord)

    #(mapped_chrom, mapped_strand, mapped_start, mapped_end, clipped_bases_start, clipped_bases_end, fastq, bytedata,binsizelist)
    df = pd.DataFrame(wholedata, columns=['mapped_chrom', 'mapped_strand' ,'mapped_start','mapped_end','clipped_bases_start', 'clipped_bases_end','fastq','ctg','st','ed','cigar','md','signal','originalsize' ])
    pathout = pathout+str(index)+".pq"
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
        ll = list(map((lambda x: int(x)+128), ll))
        elist.extend(ll)
    return bytes(elist)

def getFastQ(f,key):

    if key in f:
        group = f[key]
        datalist = list(group.values())
        fastq = datalist[0]
        fastq =  fastq.value
        return fastq
    return None

def getRecord(fdata,aligner):
    binsize = 60
    # print(fl)
    with h5py.File(fdata[0], 'r') as f:

        group = f['/Analyses/RawGenomeCorrected_000/BaseCalled_template']
        datalist = list(group.values())
        alignment = datalist[0]
#        attr = list(alignment.attrs.items())
#        print(attr)
        event = datalist[1]
        raw_start = event.attrs.get('read_start_rel_to_raw')
        dlen = len(event[()])-1
        end = event.value[dlen][2]

        group = f['/Raw/Reads/']
        readid = list(group.values())[0].name
        signalkey = readid + "/Signal"
        signal = f[signalkey].value

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
        for n in range(len(event.value)):
            #           print(event.value[n][0], event.value[n][1], event.value[n][2], event.value[n][3], event.value[n][4])
            start = raw_start + event.value[n][2]
            size = event.value[n][3]
            binnedlist.append(signal[start:start + size])
            binsizelist.append(size)

        binnedlist = list(map(equalize, binnedlist))
        bytedata = jointByteArray(binnedlist)

        fastq = getFastQ(f,'/Analyses/Basecall_1D_001/BaseCalled_template')
        if fastq == None:
            fastq = getFastQ(f,'/Analyses/Basecall_1D_000/BaseCalled_template')

        ctg,st,ed,cigar,md = "",0,0,"",""
        try:
            for hit in aligner.map(fastq.split("\n")[1]):
                ctg,st,ed,cigar,md = hit.ctg, hit.r_st, hit.r_en, hit.cigar_str,hit.MD
                break #take first hit
        except:
            pass

        tp = (mapped_chrom,mapped_strand,mapped_start,mapped_end,clipped_bases_start,clipped_bases_end,fastq,ctg,st,ed,cigar,md,bytedata,binsizelist)
        return tp

    return None


if __name__ == "__main__":
    path = 'C:\\Users\\hirok\\Desktop\\covid19\\kimivt\\index.txt'
    pathout = 'C:\\Users\\hirok\\Desktop\\covid19\\kimivt\\kinivt'
    refg = 'C:\\Users\\hirok\\Desktop\\covid19\\data\\Korea2.txt'
    index = str(1)
    main(path,pathout,refg,index)