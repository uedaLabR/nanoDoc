import glob
from sklearn.model_selection import train_test_split

import tensorflow as tf  # add
import numpy as np
import cnn_network
import nanoDocUtil
from tensorflow.keras.layers import GlobalAveragePooling1D
import numpy as np
from tensorflow.keras import Model


DATA_LENGTH_UNIT = 60
DATA_LENGTH = DATA_LENGTH_UNIT * 5 + 20


def getModel():
    shape1 = (None, DATA_LENGTH, 2)
    num_classes_org = 1024
    with tf.device("/cpu:0"):
        model = cnn_network.build_network(shape=shape1, num_classes=num_classes_org)
        model._layers.pop() #remove last layer
        model._layers.pop() #remove last layer
        model._layers.pop() #remove last layer
        flat = GlobalAveragePooling1D()(model.layers[-1].output)
        model_t = Model(inputs=model.input, outputs=flat)
        model_t.summary()
    return model_t


def getSD(cnt, coeffA, coeffB):
    # y=b*(x**a) to approximate std value for the score
    sd = coeffB * (cnt ** coeffA)
    minv = 0.0001
    if sd < minv:
        sd = minv
    return sd


def callMinusStrand(wfile,coeffA, coeffB, uplimit, takeparcentile, seq, refpr, targetpr, model_t, fw, chrom, chromtgt, start,
                    end,refpq, targetpq, minreadlen):
    n = end
    idx = 0
    strand = "-"
    while n > start:
        subs = seq[idx:idx + 5]
        cnt,cntref = eachProcess(wfile,n, subs, strand, coeffA, coeffB, uplimit, takeparcentile, seq, refpr, targetpr, model_t, fw, chrom,
                    chromtgt)


        n = n - 1
        idx = idx + 1


def callPlusStrand(wfile,coeffA, coeffB, uplimit, takeparcentile, seq, refpr, targetpr, model_t, fw, chrom, chromtgt, start,
                   end,refpq, targetpq,minreadlen):
    strand = "+"
    for n in range(start, end):
        subs = seq[(n - start):(n - start) + 5]
        cnt,cntref = eachProcess(wfile,n, subs, strand, coeffA, coeffB, uplimit, takeparcentile, seq, refpr, targetpr, model_t, fw, chrom,
                    chromtgt)


def eachProcess(wfile,n, subs, strand, coeffA, coeffB, uplimit, takeparcentile, seq, refpr, targetpr, model_t, fw, chrom,
                chromtgt):
    weight_path = wfile + str(subs) + "/model_t_ep_1.h5"
    model_t.load_weights(weight_path)
    cnt = 0
    cntref = 0
    # target signal
    tgd = targetpr.getData(chromtgt, strand, n, uplimit)
    if tgd is not None:    
        rawdatas, cnt = tgd

    # reference signal
    refd = refpr.getData(chrom, strand, n, cnt * 2)
    if refd is not None:
        refdatas, cntref = refd
    

    if (cnt < 10 or cntref < 10 or (tgd is None) or (refd is None)):
        infos = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}".format(n, str(subs), cnt, cnt, 0, 0, 0,
                                                                                      0, 0, "0", "0", "0")
        print(infos)
        fw.writelines(infos + "\n")
        fw.flush()
        return (cnt,cntref)

    print(cnt,cntref,chrom,chromtgt)
    rawdatas, cnt = nanoDocUtil.reducesize(rawdatas, cntref // 2)  # raw data have at least half size as reference data

    refdata1, refdata2 = train_test_split(refdatas, test_size=(cntref // 2))
    refdata1, medr, madr = nanoDocUtil.getFormat(refdata1)
    refdata2, medr, madr = nanoDocUtil.getFormat(refdata2)

    xref = model_t.predict(refdata1)
    xref2 = model_t.predict(refdata2)
    dist1, dist2 = nanoDocUtil.getDist(xref, xref2)
    #
    acc1r, totolcnt1r = nanoDocUtil.acc(dist1)
    acc2r, totolcnt2r = nanoDocUtil.acc(dist2)
    accr = (acc1r + acc2r) / 2

    thres1 = np.abs(acc1r - takeparcentile).argmin()  # percentile one direction
    thres2 = np.abs(acc2r - takeparcentile).argmin()  # percentile reverse direction
    thres = int((thres1 + thres2) / 2)
    totalcnt = (totolcnt1r + totolcnt2r) / 2

    #
    rowdata, med, mad = nanoDocUtil.getFormat(rawdatas)
    medratio = med / medr
    medratios = '{:.5f}'.format(medratio)

    xrow = model_t.predict(rowdata)

    dist1, dist2 = nanoDocUtil.getDist(xref, xrow)
    acc1, totalcnt1 = nanoDocUtil.acc(dist1)
    acc2, totalcnt2 = nanoDocUtil.acc(dist2)

    diff1 = accr - acc1
    overthrs1 = sum(diff1[thres:]) / totalcnt1

    diff2 = accr - acc2
    overthrs2 = sum(diff2[thres:]) / totalcnt2

    if cnt < 50:
        overthres1 = 0
        overthres2 = 0

    scoreDisplay1 = '{:.7f}'.format(overthrs1)
    scoreDisplay2 = '{:.7f}'.format(overthrs2)

    sd = getSD(cnt, coeffA, coeffB)
    if overthrs1 < 0:
        overthrs1 = 0
    if overthrs2 < 0:
        overthrs2 = 0

    score = (overthrs1 + overthrs2) / sd

    if score > 1:
        score = 1
    if score < 0:
        score = 0


#    if score > 500:
#        score = 500
#    score = score/500     #normalize
    
    scoreDisplay3 = '{:.7f}'.format(score)

    infos = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}".format(n, str(subs), cnt, cntref, med, mad,
                                                                                  medr, madr, medratios, scoreDisplay1,
                                                                                  scoreDisplay2, scoreDisplay3)
    print(infos)
    fw.writelines(infos + "\n")
    fw.flush()
    return (cnt,cntref)


def modCall(wfile, paramf, ref, refpq, targetpq, out, chrom, chromtgt, start, end, strand, minreadlen):
    coeffA, coeffB, uplimit, takeparcentile = nanoDocUtil.readParam(paramf)
    

    if chrom == "" :
        chrom = nanoDocUtil.getFirstChrom(ref)
        chromtgt = chrom
        print("modcallinit",chrom)      
    seq = nanoDocUtil.getSeq(ref, chrom, start, end, strand)                         
    if end < 0:
        end = len(seq)

    coeffA,coeffB,uplimit,takeparcentile = nanoDocUtil.readParam(paramf)


    refpr = nanoDocUtil.PqReader(refpq, minreadlen)
    targetpr = nanoDocUtil.PqReader(targetpq, minreadlen)
    model_t = getModel()

    fw = open(out, mode='w')

    if strand == "-":
        callMinusStrand(wfile,coeffA, coeffB, uplimit, takeparcentile, seq, refpr, targetpr, model_t, fw, chrom, chromtgt,
                        start, end,refpq, targetpq,minreadlen)
    else:
        callPlusStrand(wfile,coeffA, coeffB, uplimit, takeparcentile, seq, refpr, targetpr, model_t, fw, chrom, chromtgt,
                       start, end,refpq, targetpq,minreadlen)

    fw.close()


import sys

if __name__ == '__main__':

    #    wfile = "/groups2/gac50430/nanopore/dataset4DL/weight5merm6A/"
    #    paramf = "/groups2/gac50430/nanopore/shell/modcall/param.txt"
    #    ref ="/groups2/gac50430/nanopore/reference/NC000913.fa"
    #    refpq = "/groups2/gac50430/nanopore/equalbinnedpq/ecrRnaIvt"
    #    targetpq = "/groups2/gac50430/nanopore/equalbinnedpq/ecrRnaNative"
    #    out = "/groups2/gac50430/nanopore/detection/ecoli/23S\m6aplus.txt"
    #    chrom = "NC_000913.3"
    #    start = 4037519+1600
    #    end =  4037519+1730

    #    modCall(wfile,paramf,ref,refpq,targetpq,out,chrom,start,end)
    wfile = sys.argv[1]
    paramf = sys.argv[2]
    ref = sys.argv[3]
    refpq = sys.argv[4]
    targetpq = sys.argv[5]
    out = sys.argv[6]
    chrom = sys.argv[7]
    start = int(sys.argv[8])
    end = int(sys.argv[9])
    strand = sys.argv[10]
    #    minreadlen = 700
    minreadlen = 200
    #    if len(sys.argv) > 11 :
    #        minreadlen = int(sys.argv[11])
    chromtgt = chrom
    # for covid analysis
    if "england" in out:
        chromtgt = "hCoV-19/England/02/2020|EPI_ISL_407073"
    if "austraria" in out:
        chromtgt = "MT007544.1"
    modCall(wfile, paramf, ref, refpq, targetpq, out, chrom, chromtgt, start, end, strand, minreadlen)


