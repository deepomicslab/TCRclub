#!/usr/bin/env python
# coding: utf-8

import warnings
warnings.filterwarnings('ignore')
from numpy import *
import sys
from collections import Counter
import pandas as pd
import random
from sklearn import *
from scipy import stats
set_printoptions(threshold=sys.maxsize, linewidth = 100)
from sklearn.metrics import classification_report
import tensorflow as tf
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Layer, Dense, Activation, Conv1D, Flatten, Dropout, Input, BatchNormalization, \
                                    MaxPooling1D, UpSampling1D, Lambda, Reshape, GlobalAveragePooling1D, \
                                    Concatenate, GlobalMaxPooling1D, Concatenate, GlobalMaxPooling2D
from tensorflow.keras import backend as K
from tensorflow.keras.callbacks import TensorBoard
from tensorflow.keras.preprocessing.image import ImageDataGenerator
from tensorflow import keras
import logging
logging.basicConfig()
import struct
print(keras.__version__, tf.__version__)
# use keras backend (K) to force channels-last ordering
K.set_image_data_format('channels_last')
tf.compat.v1.disable_eager_execution()
import pickle, re


class Sampling1(Layer):
    """Uses (z_mean, z_log_var) to sample z, the vector encoding a digit."""
    def call(self, inputs):
        z_mean, z_log_var = inputs
        batch = tf.shape(z_mean)[0]
        dim = tf.shape(z_mean)[1]
        epsilon = tf.keras.backend.random_normal(shape=(batch, dim), mean=0., stddev=1e-12)
        return z_mean + tf.exp(0.5 * z_log_var) * epsilon


class CenterLossLayer(Layer):
    def __init__(self, mat_shape=(1032,32), alpha=0.5, **kwargs):
        super().__init__(**kwargs)
        self.alpha = alpha
        self.mat_shape = mat_shape

    def build(self, input_shape):
        self.centers = self.add_weight(name='centers', shape=self.mat_shape, initializer='uniform', trainable=False)
        super().build(input_shape)

    def call(self, inputs, mask=None):
        latent = inputs[0] # x[0] is Nx30, x[1] is Nx838 onehot, self.centers is 838x30
        labels = inputs[1]
        delta_centers = K.dot(K.transpose(labels), (K.dot(labels, self.centers) - latent))  # 838x30
        center_counts = K.sum(K.transpose(labels), axis=1, keepdims=True) + 1  # 838x1
        delta_centers /= center_counts
        new_centers = self.centers - self.alpha * delta_centers
        self.add_update((self.centers, new_centers), inputs)
        self.result = latent - K.dot(labels, self.centers)
        self.result = K.sum(self.result ** 2, axis=1, keepdims=True) #/ K.dot(x[1], center_counts)
        return self.result # Nx1

    def compute_output_shape(self, input_shape):
        return K.int_shape(self.result)
    
    def get_config(self):
        config = super().get_config().copy()
        config.update({'alpha': self.alpha})
        return config



def zero_loss(y_true, y_pred):
    return 0.5 * K.sum(y_pred, axis=0)


if __name__ == '__main__':

    K.clear_session()
    random.seed(12345);tf.random.set_random_seed(12345)

    latent_dim1 = 32
    lambda_centerloss = 0.1
    input_img1 = Input(shape=x_dim, name='input1')
    aux_input = Input(shape=(e_dim,), name = 'input2')
    x = Conv1D(128, 3, activation='relu', padding='same', name='enc_conv1')(input_img1)
    x = Conv1D(128, 3, activation='relu', padding='same', name='enc_conv2')(x)
    x = Conv1D(128, 3, activation='relu', padding='same', name='enc_conv3')(x)
    x = Conv1D(256, 5, activation='relu', padding='same', name='enc_conv4')(x)
    x = Conv1D(256, 5, activation='relu', padding='same', name='enc_conv5')(x)
    x = Conv1D(256, 5, activation='relu', padding='same', name='enc_conv6')(x)
    x = Flatten(name='enc_hidden')(x)
    z_mean = Dense(latent_dim1, name='mu')(x)
    z_log_sigma = Dense(latent_dim1, name='log_sigma')(x)
    z_sampled = Sampling1()([z_mean, z_log_sigma])
    #x_sampled = Lambda(sampling, name='x_sampled')([x_mean, x_log_sigma])
    encoder1 = Model(input_img1, z_sampled, name='encoder1')
    side = CenterLossLayer(mat_shape=(e_dim,latent_dim1),alpha=0.5, name='center_side')([z_sampled, aux_input])
    center_loss = Model([input_img1,aux_input],side, name='center')

    def vae_loss(y_true, y_pred, beta=0.1):
        # cross-entropy loss (recontruction loss)
        #reconstruction_loss = keras.losses.categorical_crossentropy(y_true, y_pred)
        reconstruction_loss = tf.reduce_mean(keras.losses.categorical_crossentropy(y_true, y_pred))*x_dim[0]
        #reconstruction_loss = K.sum(K.binary_crossentropy(y_pred, y_true), axis=1)
        # KL divergence loss (regularizer)
        kl_loss = 0.5 * K.sum(K.exp(z_log_sigma) + K.square(z_mean) - 1. - z_log_sigma, axis=1)
        return reconstruction_loss + beta*kl_loss

    encoded_input1 = Input(shape=(latent_dim1,), name='encoded_input')

    d1 = Dense(900, name='dec_hidden1')(encoded_input1)
    d1 = BatchNormalization()(d1)
    d1 = Activation("relu")(d1)
    decoded_len = Dense(l_dim, activation='softmax', name='dec_len')(d1)

    decoded_x = Reshape((30,30), name='dec_reshape')(d1)
    decoded_x = Conv1D(256, 5, activation='relu', padding='same', name='enc_conv7')(decoded_x)
    decoded_x = Conv1D(256, 5, activation='relu', padding='same', name='enc_conv8')(decoded_x)
    decoded_x = Conv1D(256, 5, activation='relu', padding='same', name='enc_conv9')(decoded_x)
    decoded_x = Conv1D(128, 3, activation='relu', padding='same', name='enc_conv10')(decoded_x)
    decoded_x = Conv1D(128, 3, activation='relu', padding='same', name='enc_conv11')(decoded_x)
    decoded_x = Conv1D(128, 3, activation='relu', padding='same', name='enc_conv12')(decoded_x)
    decoded_seq = Conv1D(21, 1, activation='softmax', padding='same', name='d1_output')(decoded_x)

    decoder_x = Model(encoded_input1,decoded_seq,name='seq_recon')
    decoder_y = Model(encoded_input1,decoded_len,name='length')

    AE = Model([input_img1,aux_input], [center_loss([input_img1,aux_input]), decoder_x(encoder1(input_img1)), 
                                        decoder_y(encoder1(input_img1))])
    AE.compile(optimizer=keras.optimizers.Adam(),loss=[zero_loss,vae_loss,'categorical_crossentropy'],
            loss_weights=[lambda_centerloss,1,1],metrics={'seq_recon':'accuracy','length':'accuracy'})
    AE.summary()

    dummy1 = zeros((trainX.shape[0], 1))
    dummy2 = zeros((validX.shape[0], 1))
    history = AE.fit([trainX, trainE],[dummy1,trainX,trainL],shuffle=True,epochs=6,batch_size=128, 
                    validation_data=([validX,validE],[dummy2,validX,validL]))

    AE.save('cVAE_len_centerloss_ept_6epoch.h5')
