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

from scipy.special import logsumexp
from scipy.spatial import distance
from scipy.optimize import linprog
from scipy.optimize import minimize
from collections import defaultdict
from mpl_toolkits.mplot3d import Axes3D

class Parameter:
    '''モデルのパラメータを保持するクラス

    Attributes
    ----------
    K: int
        立地点数
    tau: float
        企業間交流費用パラメータ
    t: float 
        通勤費用パラメータ
    S_total: float 
        都市の総床面積
    theta_firm: float
        企業のランダム項パラメータ 
    theta_house: 
        家計のランダム項パラメータ
    T: numpy.ndarray (K, K)
        通勤費用行列
    Ker: numpy.ndarray (K, K)
        Tのカーネル行列
    D: numpy.ndarray (K, K)
        空間割引行列
    S: numpy.ndarray (K, )
        各立地点の床面積
    M: float
        総企業数
    N: float
        総家計数
    '''

    def __init__(self, Col, distance_matrix, t, tau, Scaling, S_total, S_bar,
                 theta_firm, theta_house, E, alpha_1, alpha_2, beta_1, beta_2, M, N):
        '''

        Parameters
        ----------
        Num_cols: int
            立地点の列の数
        K: int
            立地点数
        distance_matrix: numpy.ndarray (K, K)
            距離行列
        tau: float
            企業間交流費用パラメータ
        t: float 
            通勤費用パラメータ
        S_total: float 
            都市の総床面積
        theta_firm: float
            企業のランダム項パラメータ 
        theta_house: float
            家計のランダム項パラメータ
        alpha_1: float
            間接効用関数のパラメータ1
        alpha_2: float
            間接効用関数のパラメータ2
        beta_1: float
            利潤関数のパラメータ1
        beta_1: float
            利潤関数のパラメータ2
        '''
        self.Col = Col
        self.K = Col * Col
        self.tau = tau
        self.t = t 
        self.Scaling = Scaling
        self.S_total = S_total
        self.S_bar = S_bar
        self.theta_firm = theta_firm
        self.theta_house = theta_house
        self.E = E
        self.alpha_1 = alpha_1
        self.alpha_2 = alpha_2
        self.beta_1 = beta_1
        self.beta_2 = beta_2
        self.M = M
        self.N = N
        self.distance_matrix = distance_matrix
        
        self.T = t * distance_matrix
        self.Ker = np.exp(-self.theta_house*self.T)
        self.D = np.exp(-tau*distance_matrix)

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
        # R, W = self.short.solve(m, n, max_itr=self.itr_max_short, err_short=self.err_short)
        # RW = np.concatenate((R, W))
        Short = self.short.Z_SD(RW)
        
        F_value = - Short \
        - (1 / 2) * np.dot(m, np.dot(self.prm.D, m))\
        + (1/self.prm.theta_firm)*np.dot(m, np.log(m/self.prm.M)) \
        + (1/self.prm.theta_house)*np.sum(n * np.log(n/self.prm.N))

        # print("len(F_value)", len(F_value))

        return F_value
    
    def long_primal_df(self, mn, RW):
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

        dF_m = - self.short.pi(R, W) + (1/self.prm.theta_firm) * (np.log(m / self.prm.M) + 1) # numpy.ndarray(K, )
        dF_n = - self.short.v(R, W) + (1/self.prm.theta_house) * (np.log(n / self.prm.N) + 1) # numpy.ndarray(K, K)
        
        dF = self.bond(dF_m, dF_n)

        return dF

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
    
    
    def logit(self, pi, v):
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
        
        max_pi = np.max(self.prm.theta_firm * pi)
        max_v = np.max(self.prm.theta_house * v)

        m_d = self.prm.M \
        * (np.exp(self.prm.theta_firm * pi - max_pi) \
        / np.sum(np.exp(self.prm.theta_firm * pi - max_pi)))

        n_d = self.prm.N \
        * np.exp(self.prm.theta_house * v - max_v) \
        / np.sum(np.exp(self.prm.theta_house * v - max_v))
        
        mn_d = self.bond(m_d, n_d)

        return mn_d


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

    def solve(self, m0, n0, long_itr, err_long):
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
        
        #Step.1 短期均衡の解を求める
        RW_hist = {}
        
        mn_hist = {}
        mn_d_hist = {}
        #Step.0 初期実行解を決める
        mn_before = mn0

        #Step.2 反復
        for k in range(1, long_itr):
            # print("="*40)
            # print("klong:", k)
            
            m_before, n_before = self.breakingdown(mn_before)
            
            #Step.1 短期均衡を解く
            R_true = ((self.prm.alpha_1 * np.sum(n_before, axis=1) + self.prm.beta_1 * m_before) / self.prm.S_bar)
            W_true = (1 / self.prm.E) * (self.prm.beta_2 * (m_before / np.sum(n_before, axis=0)) + self.prm.alpha_2)
            v_true = np.maximum(0.1, self.prm.E * W_true * np.ones((self.prm.K, self.prm.K))\
                    - (self.prm.alpha_1 * np.log(R_true) * np.ones((self.prm.K, self.prm.K))).T\
                    - self.prm.alpha_2 * np.log(W_true) * np.ones((self.prm.K, self.prm.K))\
                    + self.prm.alpha_1 * np.log(self.prm.alpha_1) + self.prm.alpha_2 * np.log(self.prm.alpha_2)\
                    - self.prm.alpha_1 - self.prm.alpha_2 - self.prm.T)
            pi_true = np.maximum(0.1, np.dot(self.prm.D, m_before)\
                    - self.prm.beta_1 * np.log(R_true) - self.prm.beta_2 * np.log(W_true)\
                    + self.prm.beta_1 * np.log(self.prm.beta_1) + self.prm.beta_2 * np.log(self.prm.beta_2)\
                    - self.prm.beta_1 - self.prm.beta_2)
            
            RW_hist[k] = np.concatenate((R_true, W_true))
            
            #Step.2 探索方向を決める
            mn_d = self.logit(pi_true, v_true)
            
            #Step.3 ステップサイズを決める
            alpha = 0.1
            
            #Step.4 解を更新する
            mn = (1 - alpha) * mn_before + alpha * mn_d
            
            max_value = k
            
            #Step.5 収束判定
            if np.max(abs((mn - mn_before) / mn_before)) < err_long:
                break
            
            mn_before = mn
        
        print("long_max_value:", max_value)
        
        m_true, n_true = self.breakingdown(mn)
        
        max_pi = np.max(self.prm.theta_firm * pi_true)
        max_v = np.max(self.prm.theta_house * v_true)
        
        G = np.sum(pi_true * m_true) + np.sum(v_true * n_true) \
          - (self.prm.M * logsumexp(self.prm.theta_firm * pi_true) / self.prm.theta_firm) \
          - (self.prm.N * logsumexp(self.prm.theta_house * v_true) / self.prm.theta_house) \
          - (np.sum(m_true * np.log(m_true / self.prm.M)) / self.prm.theta_firm) \
        - (np.sum(n_true * np.log(n_true / self.prm.N)) / self.prm.theta_house)
        
        print("Lyapunov:", G)

        return m_true, n_true, RW_hist[max_value]