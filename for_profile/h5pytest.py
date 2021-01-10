import h5py
import numpy as np
import os
import glob


def valid(fs):

    try:
        f = h5py.File(fs, 'r')
        dset = f['Raw/Reads/']
        for key, val in dset.items():
            d2 = f['Raw/Reads/' + key]
            for key, val in d2.items():
                n1 = np.array(val[:])
    except:
        f.close()
        return False

    f.close()
    return True


path = 'C:\\Users\\hirok\\Desktop\\covid19'

for s in os.listdir(path):
    if os.path.isdir(os.path.join(path, s)):
        pp = path + "/" + s
        print(pp)
        files = glob.glob(pp + "/*.fast5")
        cnt, cnta = (0, 0)
        for file in files:
            if valid(file):
                print(file)
                cnt = cnt + 1
            cnta = cnta + 1
        print(cnt, cnta)