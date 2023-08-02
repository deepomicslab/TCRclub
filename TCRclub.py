#!/usr/bin/env python
# coding: utf-8

import argparse
import tensorflow as tf
import os
import numpy as np
import random
from numba import cuda
import pandas as pd
from tensorflow import keras
from tensorflow.keras import backend as K
from autoencoder.cVAE import Sampling1, CenterLossLayer, amino_onehot_encoding
from collections import Counter, defaultdict
import argparse
import torch
from pyseat.SEAT import SEAT
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

def setup_seed(seed):
     torch.manual_seed(seed)
     os.environ['PYTHONHASHSEED'] = str(seed)
     torch.manual_seed(seed)
     torch.cuda.manual_seed(seed)
     torch.cuda.manual_seed_all(seed)
     np.random.seed(seed)
     random.seed(seed)
     torch.backends.cudnn.deterministic = True
     tf.random.set_seed(seed)


'''
    k k neighbours
    T is the matrix of the TCR embeddings n x f, f is the dimensions of embedding
    R is the matrix of sc rnq seq expression, n x g, g means gene
    W: l x f, diagonal matrix as vector
    A: 1 x n, diagonal matrix as vector
    C: n x n, adjecent matrix
''' 

def consensus_seat(consensus_matrix, cutoff=0.002):

    seat = SEAT(affinity="precomputed",
                sparsification="affinity",
                objective="SE",
                split_se_cutoff = cutoff, 
                strategy="top_down")
    seat.fit_predict(consensus_matrix) 

    label = seat.clubs
    groups = []
    for i in set(label):
        groups.append(np.argwhere(np.array(label) == i).reshape(-1))
    
    return groups


class TCRclub():
    def __init__(self, T, RR, k=5, alpha = 0, beta=0, fixed_ini=False):
        self.TCR = T
        #self.R = R,
        self.RR = RR #n x n 
        n, f = self.TCR.shape
        self.k = k

        self.W = torch.ones((1, f))
        self.A = torch.ones((n, 1))
        self.W_grads = torch.zeros(self.W.shape)
  
        if fixed_ini:
            self.C = self.initialize_C()
        else:
            self.C = self.random_initialize() #default
            
        #for optim_W
        self.Hpq = torch.zeros((f, f))
        self.Tip_Tjp = self.TCR.T.unsqueeze(-1) * self.TCR.T.unsqueeze(1) # f x n x n
        self.alpha = torch.ones(n)*alpha
        self.beta = torch.diag(torch.ones(f)*beta)
    
    def initialize_C(self):

        C = torch.zeros(self.RR.shape)
        TWT = torch.matmul(torch.mul(self.TCR, self.W), self.TCR.T)  # n x n
        similarity = (self.RR - torch.mul(self.A, TWT)).pow(2)  # n x n
        values, indices = similarity.topk(self.k, dim=1, largest=False)

        for i in np.arange(self.RR.shape[0]):
             C[i, indices[i]] = 1

        temp_C = C + C.T
        temp_C[temp_C>1] = 1
        
        C.zero_()
        for row, temp_s in enumerate(temp_C):
            ref_indices = torch.nonzero(temp_s).reshape(-1).cpu().numpy()
            for i in (set(ref_indices)&set(indices[row].cpu().numpy())):
                C[row][i] = 1
        
        #C.fill_diagonal_(0) #revise here

        return C
        
    def random_initialize(self):

        C = torch.zeros(self.RR.shape)
        TWT = torch.matmul(torch.mul(self.TCR, self.W), self.TCR.T)  # n x n
        similarity = (self.RR - torch.mul(self.A, TWT)).pow(2)  # n x n
        values, indices = similarity.topk(self.k, dim=1, largest=False)

        for i in np.arange(C.shape[0]):
            C[i, indices[i]] = 1
        
        temp_C = C + C.T 
        temp_C[temp_C>1] = 1
        
        top_values, indices = similarity.topk(self.k, dim=1, largest=False)
        probility = temp_C * similarity

        for i in np.arange(top_values.shape[0]):
            Q1, Q2, Q3, Q4 = torch.quantile(top_values[i], q=torch.tensor([0,0.25,0.75,1]))
            upperbound = (Q3-Q2)*1.5 + Q3
            probility[i][probility[i] > upperbound] = 0

        probility[probility!=0] = 1 / probility[probility!=0]
        probility = probility / torch.sum(probility, dim=1, keepdim=True)
        
        for i, row in enumerate(probility):
            sample_num = self.k if self.k <= torch.count_nonzero(row).item() else torch.count_nonzero(row>0).item()
            row_indice = torch.multinomial(row, sample_num)
            C[i][row_indice] = 1
        
        return C

    
    def to(self):
        for name, param in self.__dict__.items():
            if name not in ['k']:
                self.__dict__[name] = param.cuda().detach()
        
    def loss(self):
        
        loss_value = 0
        TWT = torch.matmul(torch.mul(self.TCR, self.W), self.TCR.T)  # n x n
        similarity = (self.RR - torch.mul(self.A, TWT)).pow(2)  # n x n
        loss_value = torch.sum(torch.mul(self.C, similarity)) + self.alpha[0]*torch.sum(self.A.pow(2)) + self.beta[0][0]*torch.sum(self.W.pow(2))

        return loss_value.item()

    
    def updateC(self):

        TWT = torch.matmul(torch.mul(self.TCR, self.W), self.TCR.T)  # n x n
        similarity = (self.RR - torch.mul(self.A, TWT)).pow(2)  # n x n
        self.C.zero_()
        values, indices = similarity.topk(self.k, dim=1, largest=False)
        for i in np.arange(similarity.shape[0]):
            self.C[i, indices[i]] = 1

        self.undirected_C(similarity)
    
    def undirected_C(self, similarity):

        temp_C = self.C + self.C.T
        temp_C[temp_C>1] = 1

        values, indices = similarity.topk(self.k, dim=1, largest=False)
        self.C.zero_()

        for row, temp_s in enumerate(temp_C):
            ref_indices = torch.nonzero(temp_s).reshape(-1).cpu().numpy()
            for i in (set(ref_indices)&set(indices[row].cpu().numpy())):
                self.C[row][i] = 1
        #self.C.fill_diagonal_(0)


    def optim_W(self):
        '''
        Hpq: Coefficient matrix of W
        Bp: HpqW = Bp
        '''
        C = self.C.detach()
        A = self.A.detach()
        Bp = []
        for i in range(self.W.shape[-1]): # f at row 
            Hpqk = self.Tip_Tjp[i].unsqueeze(0) * self.Tip_Tjp * C * (self.A.pow(2))   # f x n x n
            self.Hpq[i,:] = torch.sum(Hpqk, dim = [1,2]) # f x 1
        Bp = torch.sum(self.Tip_Tjp * self.RR * C * self.A, dim = [1,2]).reshape(-1, 1)
        self.Hpq = self.Hpq + self.beta
        try:
            self.W = torch.linalg.solve(self.Hpq, Bp).reshape(1,-1)
        except RuntimeError:
            print("W is singular!")
            #self.W = torch.matmul(torch.linalg.pinv(self.Hpq), Bp).reshape(1,-1)
            self.W = torch.linalg.lstsq(self.Hpq.cpu(), Bp.cpu(), driver='gelsd').solution.reshape(1,-1).cuda()

        return self.W
    
    
    def optim_A(self):
        W = self.W.detach()
        TWT = torch.matmul(torch.mul(self.TCR, W), self.TCR.T)
        left = torch.sum(self.C * TWT.pow(2), dim=1) + self.alpha # n x 1
        right = torch.sum(self.C * TWT * self.RR, dim=1) # n x 1
        if torch.sum(left==0) != 0:
            return "break"

        self.A = (right / left).reshape(-1, 1)

        return self.A

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--tcr_file", type=str)
    parser.add_argument("--rna_file", type=str)
    parser.add_argument("--k", type=int, default=5)
    parser.add_argument("--repeat_times", type=int, default=50)
    parser.add_argument("--alpha", type=float, default=1e-7)
    parser.add_argument("--beta", type=float, default=1e-7)
    parser.add_argument("--single_cutoff", type=float, default=0.0001)
    parser.add_argument("--con_cutoff", type=float, default=0.002)
    parser.add_argument("--con_topk", type=int, default=25, help="the topk results with smallest loss.")
    parser.add_argument("--out", type=str, default="outputs")
    parser.add_argument("--multiple_sample", action='store_true')
    parser.add_argument("--fixed_initialization", action='store_true')
    args = parser.parse_args()

    setup_seed(0)
    tcr_path = args.tcr_file
    rna_path = args.rna_file
    tcr_file = pd.read_csv(tcr_path, sep = ',', index_col='barcode')  #the first column has to be "barcode"
    if args.multiple_sample:
        print("please note that multiple_sample is selected.")
        tcr_file['cdr3'] = tcr_file['cdr3'].str.cat(tcr_file['sample'], sep=':')
        unique_tcr_file = tcr_file.drop_duplicates(subset='cdr3', keep='first')
        onehot = amino_onehot_encoding(unique_tcr_file['cdr3'].map(lambda x: str(x).split(':')[0]).tolist(), max_length=30)
    else:
        unique_tcr_file = tcr_file.drop_duplicates(subset='cdr3', keep='first')
        onehot = amino_onehot_encoding(unique_tcr_file.cdr3.tolist(), max_length=30)
    rna_file = pd.read_csv(rna_path, sep = ',', index_col='barcode') #the first column has to be "barcode"

    if not os.path.exists(args.out):
        os.mkdir(args.out)

    #load the autoencoder
    K.clear_session()
    model = keras.models.load_model("autoencoder/new_classifier_k_20_lambda_1_softmax_adam_blosum_20epoch.h5",
            custom_objects={'Sampling1':Sampling1,'CenterLossLayer':CenterLossLayer}, compile=False).layers[2]
    T = model.predict(onehot)
    cuda.select_device(0)
    cuda.close()

    #obtain the pair-wise expression distances between TCR clones:RR
    R = np.float32(rna_file.values)
    RR = np.matmul(R, R.T)
    rna = pd.DataFrame(RR) 
    rna['cdr3'] = tcr_file.reset_index()['cdr3']
    RR_txn = rna.groupby('cdr3').agg('mean')
    rna = pd.DataFrame(RR_txn.values.T)
    rna['cdr3'] = tcr_file.reset_index()['cdr3']
    RR_txt = rna.groupby('cdr3').agg('mean')
    u, sigma, vt = np.linalg.svd(RR_txt.values)
    if len(sigma) < len(RR_txt):
        S = np.diag(np.hstack(sigma, np.zeros(len(RR_txt)-len(sigma))))
    else:
        S = np.diag(sigma)
    filter_RR = np.dot(np.dot(u,S),vt)

    if filter_RR.shape[0] != T.shape[0]:
        raise Exception("T cells in scRNA file are not matched with scTCR file. Please check your input files.")
    
    TCR = torch.from_numpy(T/np.std(T,axis=0))
    RNA = torch.from_numpy(filter_RR)
    
    inputs = {}
    inputs['TCR'] = TCR.numpy()
    inputs['RR'] = RNA.numpy()
    #np.save(os.path.join(args.out, "inputs.npy"), inputs)

    results = defaultdict(dict)
    for repeat_time in np.arange(args.repeat_times):
        epochs = 1000
        clubproducer = TCRclub(TCR, RNA, k=args.k, alpha=args.alpha, beta=args.beta, fixed_ini=args.fixed_initialization)
        clubproducer.to()
        loss = clubproducer.loss()
        print("The start loss is {}".format(loss))
        minloss = 1e+25
        prevloss = 0
        with torch.no_grad():
            for epoch in np.arange(epochs):
                W = clubproducer.optim_W()
                loss = clubproducer.loss()
                #print('W', loss)
                A = clubproducer.optim_A()
                if A == 'break':
                    break
                loss = clubproducer.loss()
                #print('A', loss)
                clubproducer.updateC()
                loss = clubproducer.loss()
                #print('C', loss),
                if loss < minloss:
                    minloss = loss
                    if minloss < 1e-7:
                        break
                if np.abs(loss-prevloss) < 1e-3:
                    break
                prevloss = loss
                #torch.cuda.empty_cache()
        print("The final loss is {}".format(minloss))
        
        results[repeat_time]['loss'] = minloss
        
        TWT = torch.matmul(torch.mul(clubproducer.TCR, clubproducer.W), clubproducer.TCR.T)  # n x n
        distance = (clubproducer.RR - torch.mul(clubproducer.A, TWT)).pow(2)
        results[repeat_time]['A'] = clubproducer.A.cpu().numpy()
        results[repeat_time]['W'] = clubproducer.W.cpu().numpy()
        
        X = distance + distance.T
        #X = X / X.max()
        
        seat = SEAT(affinity="precomputed",
                    sparsification="affinity",
                    objective="SE",
                    split_se_cutoff = args.single_cutoff, 
                    strategy="top_down")
        seat.fit_predict(X.cpu().numpy())
        label = seat.clubs
        groups = []
        for i in set(label):
            groups.append(np.argwhere(np.array(label) == i).reshape(-1))
        results[repeat_time]['clustering_result'] = groups
        

    #np.save(os.path.join(args.out, "k{}_maxk{}_beta{}.npy".format(args.k, args.k, args.beta)), results)
    
    tcr_file = tcr_file.reset_index()
    out_file = tcr_file.copy()
    clustered_idx = set(sum([list(row) for row in results[0]["clustering_result"]],[]))
    clustered_idx = list(clustered_idx)
    clustered_idx = np.array(clustered_idx)

    top_loss = {}
    puritys = defaultdict(list)
    consensus_purity = {}
    consensus_matrice = {}
    consensus_judge = {}
    combined_clusters = []

    results = dict(sorted(results.items(), key=lambda x: x[1]["loss"]))
    top_results = {k:v for i, (k, v) in enumerate(results.items()) if i in range(0,args.con_topk)}
    
    consensus_matrix = np.zeros((len(clustered_idx),len(clustered_idx)))
    for key, values in top_results.items():
        #print("cluster purity:", len(sum([list(cluster) for cluster in values["clustering_result"] if len(cluster)>1], []))/len(unique_tcr_file))
        for cluster in values["clustering_result"]:
            for a in np.arange(len(cluster)):
                consensus_matrix[cluster[a]][cluster[a]] += 1
                for b in np.arange(a+1, len(cluster)):
                    consensus_matrix[cluster[a]][cluster[b]] += 1
                    consensus_matrix[cluster[b]][cluster[a]] += 1
    
    groups = consensus_seat(consensus_matrix, cutoff=args.con_cutoff)


    out_file["club"] = np.nan
    for cluster_ID, cluster in enumerate(groups):
        tcrs = unique_tcr_file.iloc[cluster].cdr3.tolist()  #tcr: cdr3:sample
        expand_tcrs_idxs = tcr_file[tcr_file.cdr3.isin(tcrs)].index.to_list() 
        out_file.iloc[expand_tcrs_idxs, out_file.columns.to_list().index("club")] = cluster_ID
    if out_file['club'].isnull().any():
        raise ValueError("Some cells are not assigned clusterID!")

    if args.multiple_sample:
        out_file['cdr3'] = out_file['cdr3'].str.replace(':.*', '', regex=True)
    out_file.to_csv(os.path.join(args.out, "consensus_result.csv"), index=False)

