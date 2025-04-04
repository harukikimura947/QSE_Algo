import time
import numpy as np
import pandas as pd
import numba
# import cProfile
import scipy.optimize as optimize
import scipy.sparse as spsp
import matplotlib.pyplot as plt
import csv
import os
from dir_exp_sci import exp_Short_sci

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

    def __init__(self, prm):
        '''

        Parameters
        ----------
        prm: Parameter
            パラメータクラス
        long: Long
            長期問題クラス
        '''
        self.prm = prm
        
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
        n_list = n.flatten()
        
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
        
        Z_short = exp_Short_sci.Short(self.prm, m, n)
        
        Short = Z_short.Z_SD(RW)
        
        # print("Z_SD:", Short)
        
        F_value = - Short \
        - (1 / 2) * np.dot(m, np.dot(self.prm.D, m))\
        + (1 / self.prm.theta_firm) * np.dot(m, np.log(m / self.prm.M)) \
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
        
        df_short = exp_Short_sci.Short(self.prm, m, n)

        dF_m = - df_short.pi(R, W) + (1/self.prm.theta_firm) * (np.log(m / self.prm.M) + 1) # numpy.ndarray(K, )
        dF_n = - df_short.v(R, W) + (1/self.prm.theta_house) * (np.log(n / self.prm.N) + 1) # numpy.ndarray(K, K)
        
        dF = self.bond(dF_m, dF_n)

        return dF

    def solve(self, m0, n0, RW_ini, long_itr, err_long):
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

        mn0 = self.bond(m0, n0)
        RW_before = RW_ini * np.ones(2 * self.prm.K)
        obj_before = 0
        
        beta = 10
        
        #Step.0 初期実行解を決める
        mn_before = mn0

        for k in range(1, long_itr):
            print("="*40)
            print("klong:", k)
            
            m_before, n_before = self.breakingdown(mn_before)
            
            #Step.1 短期均衡を解く
            short = exp_Short_sci.Short(self.prm, m_before, n_before)
            
            R, W = short.hanyou(sci_alg="L-BFGS-B", ftol=1e-8)
            
            RW = np.concatenate((R, W))
            
            print("RW_true:", np.mean(RW))
            
            #Step.3 探索方向を決める
            mn_d = self.logit(R, W, short)
            
            #Step.4 解を更新する
            obj = self.Z_LP(mn_before, RW)
            print("obj_before:", obj_before)
            print("obj:", obj)
            
            if k > 1:
                if obj > obj_before:
                    beta += 2
            
            alpha = 1 / beta
            mn = (1 - alpha) * mn_before + alpha * mn_d
            
            print("beta:", beta)
            
            #Step.2 収束判定
            if np.max(abs((mn - mn_before) / mn_before)) < err_long:
                break
            
            max_value = k
                
            mn_before = mn
            obj_before = obj
        
        print("long_max_value:", max_value)
        
        m_true, n_true = self.breakingdown(mn)
        
        print("Lyapunov:", self.Lyapunov(m_true, n_true, RW))

        return m_true, n_true, RW
    
    def logit(self, R, W, short):
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
        
        pi = short.pi(R, W)
        v = short.v(R, W)

        m_d = self.prm.M \
        * (np.exp(self.prm.theta_firm * pi) \
        / np.sum(np.exp(self.prm.theta_firm * pi)))

        n_d = self.prm.N \
        * (np.exp(self.prm.theta_house * v) \
        / np.sum(np.exp(self.prm.theta_house * v)))
        
        mn_d = self.bond(m_d, n_d)

        return mn_d

    def long_armijo(self, Z_LP, dZ, mn, mn_d, RW, c_1=0.5, beta=0.6):
        
        t = 1.0

        while True:
            print("Armijo")
            if Z_LP(mn + t * mn_d, RW) < Z_LP(mn, RW) + c_1 * t * np.dot(dZ(mn, RW), mn_d):
                break
            t *= beta
            
        return t
    
    def Lyapunov(self, m, n, RW):
        
        R = RW[:self.prm.K]
        W = RW[self.prm.K:]
        
        Lyap_short = exp_Short_sci.Short(self.prm, m, n)
        
        pi = Lyap_short.pi(R, W)
        v = Lyap_short.v(R, W)
        
        G = np.sum(pi * m) + np.sum(v * n)\
        - (self.prm.M * np.log(np.sum(np.exp(self.prm.theta_firm * pi))) / self.prm.theta_firm) - (self.prm.N * np.log(np.sum(np.exp(self.prm.theta_house * v))) / self.prm.theta_house) \
        - (np.sum(m * np.log(m / self.prm.M)) / self.prm.theta_firm) - (np.sum(n * np.log(n / self.prm.N)) / self.prm.theta_house)
        
        return G