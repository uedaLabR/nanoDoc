import pyarrow.parquet as pq
import matplotlib.pyplot as plt
import numpy as np
# from tensorflow.keras.utils.training_utils import multi_gpu_model
import tensorflow
import cnn_network
import nanoDocUtil

DATA_LENGTH_UNIT = 60
DATA_LENGTH = DATA_LENGTH_UNIT * 5 + 20


def prepData(s_data):
    train_x = []
    test_x = []
    train_y = []
    test_y = []
    p_flg = 0
    flg = 0
    totalcnt = 0
    df = pq.read_table(s_data).to_pandas()
    df1 = df[['labelidx', 'signal', 'originalsize']]

    for idx, row in df1.iterrows():

        flg = row[0]
        signal = np.array(list(row[1]))
        signal = nanoDocUtil.zeropadding10(signal)
        signal = np.array(signal)
        signal = signal.astype('float32') / 255.

        originalsize = np.array(nanoDocUtil.extendAry(row[2]))
        originalsize = nanoDocUtil.zeropadding10(originalsize)

        testidx = (idx % 12 >= 10)
        if testidx:
            test_x.extend(signal)
            test_x.extend(originalsize)
            test_y.append(flg)
        else:
            train_x.extend(signal)
            train_x.extend(originalsize)
            train_y.append(flg)

        if idx == 10:
            plt.plot(signal)
            plt.show()
            plt.plot(originalsize)
            plt.show()

        totalcnt = totalcnt + 1

        if totalcnt % 12000 == 0:
            print(totalcnt, idx, row)

    print("totalcnt", totalcnt)
    train_x = np.array(train_x)
    test_x = np.array(test_x)
    train_y = np.array(train_y)
    test_y = np.array(test_y)
    num_classes = np.unique(train_y).size
    print(type(train_x))

    print("train_x.shape", train_x.shape)
    print("test_x.shape", test_x.shape)

    print(num_classes, 'classes')

    print('y_train shape:', train_y.shape)
    print('y_test shape:', test_y.shape)

    train_x = np.reshape(train_x, (-1, DATA_LENGTH, 2))
    test_x = np.reshape(test_x, (-1, DATA_LENGTH, 2))
    train_y = np.reshape(train_y, (-1, 1,))
    test_y = np.reshape(test_y, (-1, 1,))

    train_y = tensorflow.keras.utils.to_categorical(train_y, num_classes)
    test_y = tensorflow.keras.utils.to_categorical(test_y, num_classes)

    print('train_x:', train_x.shape)
    print('train_y:', train_y.shape)
    print('test_x shape:', test_x.shape)
    print('test_y shape:', test_y.shape)

    return train_x, test_x, train_y, test_y, num_classes


def main(s_data, s_out,epochs):

    batch_size = 1024
    num_classes = 1024
    gpu_count = 4
    shape1 = (None, DATA_LENGTH, 2)

    # with tf.device("/cpu:0"):
    model = cnn_network.build_network(shape=shape1, num_classes=num_classes)
    model.summary()
    # model = multi_gpu_model(model, gpus=gpu_count)  # add

    train_x, test_x, train_y, test_y, num_classes = prepData(s_data)
    model.compile(loss='categorical_crossentropy',
                  optimizer=tensorflow.keras.optimizers.Adam(lr=0.0005, beta_1=0.9, beta_2=0.999, epsilon=None, decay=0.0,
                                                  amsgrad=False),
                  # optimizer=opt,
                  metrics=['accuracy'])

    for ep in range(epochs):
        model.fit(train_x, train_y, epochs=1, batch_size=batch_size, verbose=1,
                  shuffle=True, validation_data=(test_x, test_y))
        outpath = s_out + str(ep) + ".hdf"
        model.save_weights(outpath)
    #        score = model.evaluate(test_x, test_y, verbose=0)
    #        print('epoch',ep,'Test loss:', score[0],'Test accuracy:', score[1])

    #     checkpoint = ModelCheckpoint(s_out+"/best_all.hdf5", monitor='val_acc', verbose=1,
    #     save_best_only=True, mode='auto', period=1)

    #     fit = model.fit(train_x, train_y, epochs=epochs, batch_size=batch_size,verbose=1,
    #                    shuffle=True, validation_data=(test_x, test_y),callbacks=[checkpoint])

    #     score = model.evaluate(test_x, test_y, verbose=0)

    # ----------------------------------------------
    # Some plots
    # ----------------------------------------------
    fig, (axL, axR) = plt.subplots(ncols=2, figsize=(10, 4))

    # loss
    def plot_history_loss(fit):
        # Plot the loss in the history
        axL.plot(fit.history['loss'], label="loss for training")
        axL.plot(fit.history['val_loss'], label="loss for validation")
        axL.set_title('model loss')
        axL.set_xlabel('epoch')
        axL.set_ylabel('loss')
        axL.legend(loc='upper right')

    # acc
    def plot_history_acc(fit):
        # Plot the loss in the history
        axR.plot(fit.history['acc'], label="loss for training")
        axR.plot(fit.history['val_acc'], label="loss for validation")
        axR.set_title('model accuracy')
        axR.set_xlabel('epoch')
        axR.set_ylabel('accuracy')
        axR.legend(loc='upper right')

    # plot_history_loss(fit)
    # plot_history_acc(fit)
    # fig.savefig('/groups2/gac50430/nanopore/dataset4DL/learnigcurve.png')
    # plt.close()


if __name__ == '__main__':
    s_data = "/fs2/groups2/gac50430/nanopore/dataset4DL/2400each.pq"
    #     s_data ="/fs2/groups2/gac50430/nanopore/dataset4DL/6000each.pq"
    s_out = "/fs2/groups2/gac50430/nanopore/dataset4DL/weight/"
    main(s_data, s_out)