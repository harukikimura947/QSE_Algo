import time
import numpy as np
import pandas as pd
import numba
import scipy.optimize as optimize
import scipy.sparse as spsp
import matplotlib.pyplot as plt
import csv
import os

from scipy.sparse import csr_matrix
from scipy.spatial import distance
from scipy.optimize import linprog
from scipy.optimize import minimize
from collections import defaultdict
from mpl_toolkits.mplot3d import Axes3D

class Long:
    '''長期均衡問題クラス

    Attributes
    ----------
    prm: Parameter
       パラメータクラス 
    '''

    def __init__(self, prm, algprm, short):
        '''

        Parameters
        ----------
        prm: Parameter
            パラメータクラス
        long: Long
            長期問題クラス
        '''
        self.prm = prm
        self.algprm = algprm
        self.short = short
        
    def bond(self, m, n):
        '''企業・家計分布を一つの結合リストにする関数

        Parameters
        ----------
        m: numpy.ndarray (K, )
            企業分布
        n: numpy.ndarray (K, K)
            家計分布

        Returns
        ----------
        mn: numpy.ndarray (K + K * K, )
            企業・家計分布の結合リスト
        '''
        
        # 行列を1次元のリストに変換
        n_list = n.ravel()
        
        mn = np.concatenate((m, n_list))

        return mn
    
    def breakingdown(self, mn):
        '''結合リストを企業・家計分布に分解する関数

        Parameters
        ----------
        mn: numpy.ndarray (K + K * K, )
            企業・家計分布の結合リスト
        
        Returns
        ----------
        m: numpy.ndarray (K, )
            企業・家計分布の結合リスト
        n: numpy.ndarray (K, K)
            家計分布
        '''
        
        m = np.array(mn[:self.prm.K])
        
        n_list = mn[self.prm.K:]
        n = np.array([n_list[i:i+self.prm.K] for i in range(0, len(n_list), self.prm.K)])

        return m, n
    
    def Z_LP(self, mn, RW):
        '''目的関数を計算する関数

        Parameters
        ----------
        mn: numpy.ndarray (K + K * K, )
            企業・家計分布の結合リスト

        Returns
        ----------
        F: float
            目的関数値
        '''

        m, n = self.breakingdown(mn)
        
        Short = self.short.Z_SD(RW, m, n)
        
        F_value = - Short \
        - (1 / 2) * m@self.prm.D@m\
        + (1 / self.prm.theta_firm) *  (m @ np.log(m / self.prm.M)) \
        + (1 / self.prm.theta_house) * np.sum(n * np.log(n / self.prm.N))

        # print("len(F_value)", len(F_value))

        return F_value
    
    def long_df(self, mn, RW):
        '''目的関数の勾配を計算する関数

        Parameters
        ----------
        mn: numpy.ndarray (K + K * K, )
            企業・家計分布の結合リスト
        RW: numpy.ndarray (2 * K, )
            価格変数の結合リスト

        Returns
        ----------
        dF: numpy.ndarray (K + K * K, )
            目的関数の勾配
        '''

        m, n = self.breakingdown(mn)

        dF_m = - self.short.pi(R, W, m) + (1/self.prm.theta_firm) * (np.log(m / self.prm.M) + 1) # numpy.ndarray(K, )
        dF_n = - self.short.v(R, W) + (1/self.prm.theta_house) * (np.log(n / self.prm.N) + 1) # numpy.ndarray(K, K)
        
        dF = self.bond(dF_m, dF_n)

        return dF

    def solve(self, m0, n0, RW_ini, err_short, err_long, obj_corr, long_itr):
        '''長期均衡を解く関数

        Parameters
        ----------
        m0: numpy.ndarray (K, )
            企業の初期分布
        n0: numpy.ndarray (K, K)
            家計の初期分布
        max_itr: int
            最大反復回数
        err_long: float
            収束判定の閾値
            
        Returns
        -------
        mn_k: numpy.ndarray (K + K * K, )
            企業・家計分布の結合リスト
        '''

        max_value = 0

        obj_before = 0
        long_iteration = []
        obj_list = []
        obj_rel_list = []
        
        beta = 10
        
        #Step.0 初期実行解を決める
        m_before = m0
        n_before = n0

        long_iteration_append = long_iteration.append
        obj_list_append = obj_list.append
        obj_rel_list_append = obj_rel_list.append

        for k in range(1, long_itr):
            # print("="*40)
            # print("klong:", k)
            
            #Step.1 短期均衡を解く
            R, W, iteration, short_obj_rel = \
            self.short.short_solve(RW_ini, m_before, n_before, err_short, short_itr=100000, rel=0)
            
            RW = np.concatenate((R, W))
            
            #Step.2 探索方向を決める
            m_d, n_d = self.logit(R, W, m_before)
            
            #Step.3 解を更新する
            long_iteration_append(k)

            #Step.4 収束判定
#             max_value = k
#             if obj_rel < 1e-4:
#                 break
            
#             if k > 1:
#                 if obj > obj_before:
#                     beta += 2
            
            # alpha = 1 / beta
            # print("alpha:", alpha)

            alpha = 0.05
            m = (1 - alpha) * m_before + alpha * m_d
            n = (1 - alpha) * n_before + alpha * n_d
            
            if np.max(abs((m - m_before) / m_before)) < err_long and\
               np.max(abs((n - n_before) / n_before)) < err_long:
                break
            
            m_before = m
            n_before = n
        
        print("long_max_value:", max_value)
        print("long_max_value:", k)
        
        m_true, n_true = m, n
        
        print("Lyapunov:", self.Lyapunov(m_true, n_true, RW))

        return m_true, n_true, RW, long_iteration, obj_list, obj_rel_list
    
    def logit(self, R, W, m):
        '''長期均衡の探索方向を計算する関数

        Parameters
        ----------
        mn: numpy.ndarray(K + K * K, )
            家計分布と企業分布の結合リスト
        RW: numpy.ndarray (K, )
            地代と賃金の結合リスト

        Returns
        -------
        mn_d: numpy.ndarray (K + K * K, )
            企業分布の探索方向と家計分布の探索方向の結合リスト
        '''
        
        pi = self.short.pi(R, W, m)
        v = self.short.v(R, W)
        
#         theta_v = self.prm.theta_house * v
#         max_v = np.max(theta_v)
#         diff_v = theta_v - max_v
        
#         log_logit = theta_v - np.log(np.sum(np.exp(diff_v))) - max_v
#         n_logit = np.exp(log_logit)

#         n_d = self.prm.N * n_logit

        m_d = self.prm.M \
        * (np.exp(self.prm.theta_firm * pi) \
        / np.sum(np.exp(self.prm.theta_firm * pi)))

        n_d = self.prm.N \
        * np.exp(self.prm.theta_house * v) \
        / np.sum(np.exp(self.prm.theta_house * v))

        return m_d, n_d

    def long_armijo(self, Z_LP, dZ, mn, mn_d, RW, c_1=0.5, beta=0.6):
        
        t = 1.0

        while True:
            print("Armijo")
            if Z_LP(mn + t * mn_d, RW) < Z_LP(mn, RW) + c_1 * t * dZ(mn, RW)@mn_d:
                break
            t *= beta
            
        return t
    
    def Lyapunov(self, m, n, RW):
        
        R = RW[:self.prm.K]
        W = RW[self.prm.K:]
        
        pi = self.short.pi(R, W, m)
        v = self.short.v(R, W)
        
        G = np.sum(pi * m) + np.sum(v * n)\
        - (self.prm.M * np.log(np.sum(np.exp(self.prm.theta_firm * pi))) / self.prm.theta_firm) - (self.prm.N * np.log(np.sum(np.exp(self.prm.theta_house * v))) / self.prm.theta_house) \
        - (np.sum(m * np.log(m / self.prm.M)) / self.prm.theta_firm) - (np.sum(n * np.log(n / self.prm.N)) / self.prm.theta_house)
        
        return G