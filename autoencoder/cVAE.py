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
from tensorflow.keras.layers import Layer
from tensorflow.keras import backend as K
from tensorflow import keras
import logging
logging.basicConfig()
import struct
print(keras.__version__, tf.__version__)
# use keras backend (K) to force channels-last ordering
K.set_image_data_format('channels_last')
tf.compat.v1.disable_eager_execution()
import pickle, re

amino_acids=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','-']
def amino_onehot_encoding(seqs, aa=amino_acids, max_length=40):
    n=len(seqs)
    matrix=zeros((n,max_length,21),dtype=int8)
    for x in range(n):
        seq=seqs[x]+(max_length-len(seqs[x]))*'-'
        for i in range(max_length):
        #for i in range(len(seq)):
            matrix[x][i][aa.index(seq[i])] = 1
    return matrix

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


