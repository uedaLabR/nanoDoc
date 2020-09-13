import pyarrow.parquet as pq

import numpy as np
import cnn_network
import tensorflow
from tensorflow.keras.layers import GlobalAveragePooling1D, Dense
from tensorflow.keras import Model
#from tensorflow.keras import Network
from tf_agents.networks.network import Network
from tensorflow.keras.optimizers import SGD
from tensorflow.keras import backend as K
import itertools
import nanoDocUtil


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
        signal = nanoDocUtil.zeropadding10(signal)
        data_x.append(signal)

        originalsize = np.array(nanoDocUtil.extendAry(row[2]))
        originalsize = nanoDocUtil.zeropadding10(originalsize)
        data_x.append(originalsize)

        cnt = cnt + 1
        if cnt % 1000 == 0:
            print(cnt)
        if cnt == count:
            break

    data_x = np.array(data_x)
    data_x = np.reshape(data_x, (count, DATA_LENGTH, 2))
    return data_x


def getNearExculudeSelf(nuc):
    first = nuc[0]
    last = nuc[4]
    nucs = ('A', 'T', 'C', 'G')
    gen = []
    for n1, n2, n3 in itertools.product(nucs, nucs, nucs):

        nc = first + n1 + n2 + n3 + last
        if nc != nuc:
            gen.append(nc)

    return gen


def prepDataNear(s_data, nuc):
    train_x = []
    test_x = []
    train_y = []
    test_y = []
    p_flg = 0
    flg = 0
    totalcnt = 0

    nucs = getNearExculudeSelf(nuc)
    paths = list(map(lambda x: s_data + "/" + x + ".pq", nucs))

    samplecnt = 0
    for path in paths:

        table = pq.read_table(path)
        df = table.to_pandas()
        #         print(df)
        df = df[['signal', 'originalsize']]
        #         print(df)
        cnt = 0
        for idx, row in df.iterrows():

            flg = samplecnt

            signal = np.array(list(row[0]))
            signal = nanoDocUtil.zeropadding10(signal)
            signal = np.array(signal)
            signal = signal.astype('float32') / 255.

            originalsize = np.array(nanoDocUtil.extendAry(row[1]))
            originalsize = nanoDocUtil.zeropadding10(originalsize)

            #             if cnt ==0:
            #                 print(path.replace(s_data+"/",""))
            #                 fmer = path.replace(s_data+"/","").replace(".pq","")
            #                 plt.title(fmer)
            #                 plt.plot(signal)
            #                 plt.show()
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
            if cnt == 12000:
                break

        samplecnt = samplecnt + 1

    #     print("totalcnt",totalcnt)

    train_x = np.array(train_x)
    test_x = np.array(test_x)
    train_y = np.array(train_y)
    test_y = np.array(test_y)
    num_classes = np.unique(train_y).size

    #     print("train_x.shape",train_x.shape)
    #     print("test_x.shape",train_x.shape)
    #     print("train_y.shape",train_x.shape)
    #     print("test_y.shape",train_x.shape)

    #     print(num_classes, 'classes')

    #     print('y_train shape:',train_y.shape)
    #     print('y_test shape:', test_y.shape)

    train_x = np.reshape(train_x, (-1, DATA_LENGTH, 2))
    test_x = np.reshape(test_x, (-1, DATA_LENGTH, 2))
    train_y = np.reshape(train_y, (-1, 1,))
    test_y = np.reshape(test_y, (-1, 1,))

    test_y = test_y - 1
    train_y = train_y - 1
    train_y = tensorflow.keras.utils.to_categorical(train_y, num_classes)
    test_y = tensorflow.keras.utils.to_categorical(test_y, num_classes)

    #     print('train_x:', train_x.shape)
    #     print('train_y:', train_y.shape)
    #     print('test_x shape:',test_x.shape)
    #     print('test_y shape:', test_y.shape)

    return train_x, test_x, train_y, test_y, num_classes


def prepData(s_data, nuc):
    train_x = []
    test_x = []
    train_y = []
    test_y = []
    p_flg = 0
    flg = 0
    totalcnt = 0
    path = s_data + "/" + nuc + ".pq"

    table = pq.read_table(path)
    df = table.to_pandas()
    df = df[['signal', 'originalsize']]
    cnt = 0
    for idx, row in df.iterrows():

        signal = np.array(list(row[0]))
        signal = nanoDocUtil.zeropadding10(signal)
        signal = np.array(signal)
        signal = signal.astype('float32') / 255.

        originalsize = np.array(nanoDocUtil.extendAry(row[1]))
        originalsize = nanoDocUtil.zeropadding10(originalsize)

        #         if cnt ==0:
        #             print(path.replace(s_data+"/",""))
        #             fmer = path.replace(s_data+"/","").replace(".pq","")
        #             plt.title(fmer)
        #             plt.plot(signal)
        #             plt.show()
        #             plt.plot(originalsize)
        #             plt.show()

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
            print(totalcnt, path, totalcnt, idx, row)
        if cnt == 12000:
            break

    #     print("totalcnt",totalcnt)

    train_x = np.array(train_x)
    test_x = np.array(test_x)
    train_y = np.array(train_y)
    test_y = np.array(test_y)
    num_classes = np.unique(train_y).size


    train_x = np.reshape(train_x, (-1, DATA_LENGTH, 2))
    test_x = np.reshape(test_x, (-1, DATA_LENGTH, 2))
    train_y = np.reshape(train_y, (-1, 1,))
    test_y = np.reshape(test_y, (-1, 1,))

    test_y = test_y - 1
    train_y = train_y - 1
    train_y = tensorflow.keras.utils.to_categorical(train_y, num_classes)
    test_y = tensorflow.keras.utils.to_categorical(test_y, num_classes)

    #     print('train_x:', train_x.shape)
    #     print('train_y:', train_y.shape)
    #     print('test_x shape:',test_x.shape)
    #     print('test_y shape:', test_y.shape)

    return train_x, test_x, train_y, test_y, num_classes


batchsize = 128
feature_out = 1280  # secondary network out for MobileNet
alpha = 0.5  # for MobileNetV2
lambda_ = 0.1  # for compact loss

#bst = "/fs2/groups2/gac50430/nanopore/dataset4DL/weight2/340.hdf"
#save_path = "/fs2/groups2/gac50430/nanopore/dataset4DL/bestweight.hdf"



def train(s_data, s_out, nuc,bestwight, epoch_num):
    batch_size = 1024
    num_classes_org = 1024
    num_classes = 63
    shape1 = (None, DATA_LENGTH, 2)

    model = cnn_network.build_network(shape=shape1, num_classes=num_classes_org)
    model.load_weights(bestwight)
    model.layers.pop()  # remove last layer
    model.layers.pop()  # remove last layer
    model.layers.pop()  # remove last layer

    for layer in model.layers:
        if layer.name == "conv1d_20":
            break
        else:
            layer.trainable = False

    flat = GlobalAveragePooling1D()(model.layers[-1].output)
    model_t = Model(inputs=model.input, outputs=flat)
    model_r = Network(inputs=model_t.input,
                      outputs=model_t.output,
                      name="shared_layer")

    prediction = Dense(num_classes, activation='softmax')(model_t.output)
    model_r = Model(inputs=model_r.input, outputs=prediction)

    optimizer = SGD(lr=5e-5, decay=0.00005)
    model_r.compile(optimizer=optimizer, loss="categorical_crossentropy")
    model_t.compile(optimizer=optimizer, loss=original_loss)

    x_ref, test_x_r, y_ref, test_y_r, num_classes_r = prepDataNear(s_data, nuc)
    x_target, test_x, train_y, test_y, num_classes = prepData(s_data, nuc)

    ref_samples = np.arange(x_ref.shape[0])

    loss, loss_c = [], []
    epochs = []
    print("training...")

    for epochnumber in range(epoch_num):
        x_r, y_r, lc, ld = [], [], [], []

        np.random.shuffle(x_target)
        np.random.shuffle(ref_samples)

        for i in range(len(x_target)):
            x_r.append(x_ref[ref_samples[i]])
            y_r.append(y_ref[ref_samples[i]])
        x_r = np.array(x_r)
        y_r = np.array(y_r)

        for i in range(int(len(x_target) / batchsize)):
            batch_target = x_target[i * batchsize:i * batchsize + batchsize]
            batch_ref = x_r[i * batchsize:i * batchsize + batchsize]
            batch_y = y_r[i * batchsize:i * batchsize + batchsize]
            # target data
            lc.append(model_t.train_on_batch(batch_target, np.zeros((batchsize, feature_out))))

            # reference data
            ld.append(model_r.train_on_batch(batch_ref, batch_y))

        loss.append(np.mean(ld))
        loss_c.append(np.mean(lc))
        epochs.append(epochnumber)

        print("epoch : {} ,Descriptive loss : {}, Compact loss : {}".format(epochnumber + 1, loss[-1], loss_c[-1]))
        if epochnumber % 1 == 0:
            model_t.save_weights(s_out +'/'+ nuc + '/model_t_ep_{}.h5'.format(epochnumber))
            model_r.save_weights(s_out +'/'+ nuc + '/model_r_ep_{}.h5'.format(epochnumber))


batchsize = 64
# feature_out = 512 #secondary network out for VGG16
# feature_out = 1280 #secondary network out for MobileNet
feature_out = 1536  # secondary network out for Inception Resnetv2
alpha = 0.5  # for MobileNetV2
lambda_ = 0.1  # for compact loss
classes = 63


def original_loss(y_true, y_pred):
    lc = 1 / (classes * batchsize) * batchsize ** 2 * K.sum((y_pred - K.mean(y_pred, axis=0)) ** 2, axis=[1]) / (
                (batchsize - 1) ** 2)
    return lc


import sys

# if __name__ == '__main__':
#     train(sys.argv[2], sys.argv[3], sys.argv[1], 10)
