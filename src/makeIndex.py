import glob
import h5py
from multiprocessing import Pool
import tqdm
from functools import cmp_to_key
from multiprocessing import Pool
import os, fnmatch
from tqdm import tqdm

def find_files(directory, pattern):
	
    retl =[]
    for root, dirs, files in os.walk(directory):
        for basename in files:
            if fnmatch.fnmatch(basename, pattern):
                filename = os.path.join(root, basename)
                retl.append(filename)
    return retl

def cmp(a, b):

    # print(type(a))
    # print(a)
    if a[1] == b[1]:# chromsome
        if a[2] == b[2]:#mapped start
            return -1 if a[3] < b[3] else 1 #mapped end
        else:
            return -1 if a[2] < b[2] else 1
    else:
        return -1 if a[1] < b[1] else 1

def extractMapandFileInfo(file):

    # cnt = 0
    try:

        with h5py.File(file, 'r') as f:
            e = '/Analyses/RawGenomeCorrected_000/BaseCalled_template/Alignment' in f
            if e:

                group = f['/Analyses/RawGenomeCorrected_000/BaseCalled_template']
                datalist = list(group.values())
                alignment = datalist[0]
                mapped_chrom = alignment.attrs.get('mapped_chrom')
                mapped_start = alignment.attrs.get('mapped_start')
                mapped_end = alignment.attrs.get('mapped_end')
                tp = (file, mapped_chrom, mapped_start, mapped_end)
                return tp

                #
                # if cnt % 10000 == 0:
                #     print("read " + str(cnt) + " files")
                # cnt = cnt + 1
    except OSError:
        pass
    except KeyError:
        pass
    return None


def mkIdx(basedir,path_w,thread):

    l = find_files(basedir, '*.fast5')
    size = len(l)
    print(str(size)+" files to process")
    with Pool(thread) as p:

        imap = p.imap(extractMapandFileInfo, l)
        filelist = list(tqdm(imap,total=size))

    filelist = list(filter(None, filelist))
    print("start sorting")
    filelist = sorted(filelist, key=cmp_to_key(cmp))
    print("end sorting")
    print("writing file")
    with open(path_w, mode='w') as f:
        for tp in filelist:
            f.writelines(','.join(map(str, tp))+"\n")
    print("end writing file")


# if __name__ == "__main__":
#     basedir = "C:\\Users\\hirok\\Desktop\\covid19\\kimivt\\*\\*.fast5"
#     path_w = 'C:\\Users\\hirok\\Desktop\\covid19\\kimivt\\index.txt'
#     import sys
#     args = sys.argv
#     main(basedir,path_w)