import os

def getIndex(path):

    unit = 4000
    cnt = 0
    idxcnt = 0
    ref = ""
    lastref = ""
    indexlist = []
    with open(path) as f:
        for s_line in f:
            ref = s_line.split(",")[1]
            if ref != lastref or cnt%4000 ==0:
                idxcnt = idxcnt +1
                # print(idxcnt)
                indexlist.append(idxcnt)
            cnt = cnt +1
            lastref = ref
    return indexlist

import sys
if __name__ == "__main__":

    args = sys.argv
    indexfile = args[1]
    outfile = args[2]
    indexes = getIndex(indexfile)
    for i in indexes:
        os.system("qsub -g gac50430 -l rt_C.small=1 -l h_rt=12:00:00 /groups2/gac50430/nanopore/shell/toequalbinparquet/parquet.sh " + indexfile + " " +outfile+" "+ str(i))