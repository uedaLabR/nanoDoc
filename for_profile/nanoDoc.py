import glob
from sklearn.model_selection import train_test_split
import pyarrow.parquet as pq
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import random
import tensorflow as tf  # add
import numpy as np
import cnn_network
import ndoc_utils

DATA_LENGTH_UNIT = 60
DATA_LENGTH = DATA_LENGTH_UNIT * 5 + 20


def loadData(s, count):
    table = pq.read_table(s)
    df = table.to_pandas()
    df1 = df[['nucb4After', 'signal', 'originalsize']]

    data_x = []
    cnt = 0
    for index, row in df1.iterrows():

        signal = np.array(list(row[1]))

        signal = zeropadding10(signal)
        data_x.append(signal)

        #         originalsize = np.array(extendAry(row[2]))
        #         originalsize = zeropadding10(originalsize)
        #         data_x.append(originalsize)
        cnt = cnt + 1
        if cnt % 1000 == 0:
            print(cnt)
        if cnt == count:
            break

    data_x = np.array(data_x)
    data_x = data_x.astype('float32') / 255.
    data_x = np.reshape(data_x, (count, DATA_LENGTH, 1))
    return data_x


def getNearExculudeSelf(nuc):
    first = nuc[0]
    last = nuc[4]
    nucs = ('A', 'T', 'C', 'G')
    gen = []
    for n1, n2, n3 in product(nucs, nucs, nucs):

        nc = first + n1 + n2 + n3 + last
        if nc != nuc:
            gen.append(nc)

    return gen


def prepData(s_data):
    train_x = []
    test_x = []
    train_y = []
    test_y = []
    p_flg = 0
    flg = 0
    totalcnt = 0
    df = pq.read_table(s_data).to_pandas()
    df = df[['signal', 'originalsize']]

    samplecnt = 0
    for path in paths:

        table = pq.read_table(path)
        df = table.to_pandas()
        #         print(df)
        df = df[['signal', 'originalsize']]

        for idx, row in df.iterrows():

            flg = samplecnt

            signal = np.array(list(row[0]))
            signal = zeropadding10(signal)
            signal = np.array(signal)
            signal = signal.astype('float32') / 255.

            originalsize = np.array(extendAry(row[1]))
            originalsize = zeropadding10(originalsize)

            if cnt == 0:
                print(path.replace(s_data + "/", ""))
                fmer = path.replace(s_data + "/", "").replace(".pq", "")
                plt.title(fmer)
                plt.plot(signal)
                #                 plt.show()
                #                 fig = plt.figure(fmer)
                plt.savefig("/groups2/gac50430/nanopore/dataset4DL/figs/" + fmer + ".png")
                plt.clf()
            #                 plt.plot(originalsize)
            #                 plt.show()

            testidx = (idx % 12 >= 10)
            if testidx:
                test_x.append(signal)
                test_x.append(originalsize)
                test_y.append(flg)
            else:
                train_x.append(signal)
                train_x.append(originalsize)
                train_y.append(flg)

            cnt = cnt + 1
            totalcnt = totalcnt + 1

            if cnt % 12000 == 0:
                print(samplecnt, totalcnt, path, totalcnt, idx, row)
            if cnt == 36000:
                break

        samplecnt = samplecnt + 1

    print("totalcnt", totalcnt)

    train_x = np.array(train_x)
    test_x = np.array(test_x)
    train_y = np.array(train_y)
    test_y = np.array(test_y)
    num_classes = np.unique(train_y).size

    print("train_x.shape", train_x.shape)
    print("test_x.shape", train_x.shape)
    print("train_y.shape", train_x.shape)
    print("test_y.shape", train_x.shape)

    print(num_classes, 'classes')

    print('y_train shape:', train_y.shape)
    print('y_test shape:', test_y.shape)

    train_x = np.reshape(train_x, (-1, DATA_LENGTH, 2))
    test_x = np.reshape(test_x, (-1, DATA_LENGTH, 2))
    train_y = np.reshape(train_y, (-1, 1,))
    test_y = np.reshape(test_y, (-1, 1,))

    test_y = test_y - 1
    train_y = train_y - 1
    train_y = keras.utils.to_categorical(train_y, num_classes)
    test_y = keras.utils.to_categorical(test_y, num_classes)

    print('train_x:', train_x.shape)
    print('train_y:', train_y.shape)
    print('test_x shape:', test_x.shape)
    print('test_y shape:', test_y.shape)

    return train_x, test_x, train_y, test_y, num_classes


def main(s_data, s_out):

    train_x, test_x, train_y, test_y, num_classes = prepData(s_data)


if __name__ == '__main__':
    s_data = "/fs2/groups2/gac50430/nanopore/dataset4DL/1200each.pq"
    s_out = "/groups2/gac50430/nanopore/dataset4DL/1200weight.h5"
    main(s_data, s_out)