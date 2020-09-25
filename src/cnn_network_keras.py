from funcy import concat, identity, juxt, partial, rcompose, repeat, take
#from tensorflow.keras.layers import Activation, Add, BatchNormalization, Concatenate, Conv1D, Dense, Dropout, concatenate,  \
#    GlobalAveragePooling1D, Input, MaxPooling1D, GaussianNoise,AveragePooling1D
#from tensorflow.keras.models import Model, save_model
#from tensorflow.keras.regularizers import l2
import keras
from keras.models import Sequential,Model
from keras.layers import Dense, Activation, Flatten, Conv1D, Input, Dropout,BatchNormalization
from keras.layers import Conv1D, MaxPooling1D, AveragePooling1D, Dropout, concatenate,\
    BatchNormalization, GaussianNoise, GlobalAveragePooling1D, Softmax
import random
from keras.layers import Reshape
from keras.layers import Multiply, multiply
import tensorflow as tf # add
import numpy as np
from funcy import concat, identity, juxt, partial, rcompose, repeat, take
from keras.callbacks import LearningRateScheduler
from keras.layers import Activation, Add, BatchNormalization, Concatenate, Conv1D, Dense, Dropout, GlobalAveragePooling1D, Input, MaxPooling1D
from keras.models import Model, save_model
from keras.optimizers import SGD
from keras.preprocessing.image import ImageDataGenerator
from keras.regularizers import l2
from keras.utils import plot_model
from operator import getitem
from keras.layers import Conv1D, MaxPooling1D, AveragePooling1D, Dropout, concatenate,BatchNormalization, GaussianNoise, Softmax
from keras import backend as K
from itertools import product
from keras.callbacks import ModelCheckpoint



def build_network(shape, num_classes):
    # Utility functions.

    def ljuxt(*fs):
        return rcompose(juxt(*fs), list)

    def conv1D(filters, kernel_size):
        return Conv1D(filters, kernel_size, padding='same', kernel_initializer='he_normal',
                      kernel_regularizer=l2(0.0001))

    def conv1D_halve(filters, kernel_size):
        return Conv1D(filters, kernel_size, padding='same', strides=2, kernel_initializer='he_normal',
                      kernel_regularizer=l2(0.0001))

    def dense(units, activation):
        return Dense(units, activation=activation, kernel_regularizer=l2(0.0001))

    # Define SqueezeNet.
    def fire_module(filters_squeeze, filters_expand):
        return rcompose(BatchNormalization(),
                        Activation('relu'),
                        conv1D(filters_squeeze, 1),
                        BatchNormalization(),
                        Activation('relu'),
                        ljuxt(conv1D(filters_expand // 2, 1),
                              conv1D(filters_expand // 2, 3)),
                        Concatenate())

    def fire_module_with_shortcut(filters_squeeze, filters_expand):
        return rcompose(ljuxt(fire_module(filters_squeeze, filters_expand),
                              identity),
                        Add())

    # Inception　Module.
    def inception():
        u1 = rcompose(AveragePooling1D(pool_size=3, strides=1, padding='same'),
                      conv1D(48, 1))
        u2 = conv1D(48, 1)
        u3 = rcompose(conv1D(16, 1),
                      conv1D(48, 3))
        u4 = rcompose(conv1D(16, 1),
                      conv1D(48, 3),
                      conv1D(48, 3))

        return rcompose(ljuxt(u1, u2, u3, u4),
                        Concatenate(axis=2))

    def convBlock(f1, k1, f2, k2, f3, k3, do_r):
        return rcompose(conv1D(f1, k1),
                        conv1D(f2, k2),
                        conv1D(f3, k3),
                        MaxPooling1D(pool_size=2),
                        BatchNormalization(),
                        Dropout(do_r))

    do_r = 0.30
    # x = inputs
    input = Input(batch_shape=shape)
    nnBlock = rcompose(GaussianNoise(stddev=0.005),
                       conv1D_halve(32, 3),
                       BatchNormalization(),
                       Dropout(do_r),
                       convBlock(48, 3, 48, 3, 48, 3, do_r),
                       convBlock(16, 1, 32, 3, 32, 3, do_r),
                       fire_module(8, 32),
                       fire_module_with_shortcut(8, 32),
                       fire_module(16, 64),
                       MaxPooling1D(),
                       fire_module_with_shortcut(16, 64),
                       BatchNormalization(),
                       Dropout(do_r),
                       Activation('relu'),
                       conv1D(16, 3),
                       conv1D_halve(16, 3),
                       convBlock(16, 3, 16, 1, 16, 3, do_r),
                       conv1D(num_classes, 1),
                       GlobalAveragePooling1D(),
                       Activation('softmax'))

    return Model(inputs=input, outputs=nnBlock(input))