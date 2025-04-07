"""Module summary.

This is Short-Algo.py.
"""

import time
import numpy as np
import pandas as pd
import scipy.optimize as optimize
import scipy.sparse as spsp
import matplotlib.pyplot as plt
import csv
import os

from scipy.spatial import distance
from scipy.optimize import linprog
from scipy.optimize import minimize
from collections import defaultdict
from mpl_toolkits.mplot3d import Axes3D

class Short:
    '''短期均衡問題クラス

    Attributes
    ----------
    prm: Parameter
       パラメータクラス 
    '''

    def __init__(self, prm, algprm):
        '''

        Parameters
        ----------
        prm: Parameter
            モデルのパラメータクラス
        algprm: Algo-Parameter
            アルゴリズムのパラメータクラス
        '''
        self.prm = prm
        self.algprm = algprm
    
    def demand(self, R, W):
        '''需要関数を設定する関数

        Parameters
        ----------
        R: numpy.ndarray(matrix) (K, )
            地代リスト
        W: numpy.ndarray(matrix) (K, )
            賃金リスト
        '''
        
        S_H = ((1 / self.prm.T) ** (1 / (1 - self.prm.alpha_1 - self.prm.alpha_2))) \
                * ((self.prm.alpha_1 / R) ** ((1 - self.prm.alpha_2) / (1 - self.prm.alpha_1 - self.prm.alpha_2)))[:, np.newaxis]\
                * ((self.prm.alpha_2 / W) ** (self.prm.alpha_2 / (1 - self.prm.alpha_1 - self.prm.alpha_2)))
        
        L_H = ((1 / self.prm.T) ** (1 / (1 - self.prm.alpha_1 - self.prm.alpha_2))) \
                * ((self.prm.alpha_1 / R) ** (self.prm.alpha_1 / (1 - self.prm.alpha_1 - self.prm.alpha_2)))[:, np.newaxis]\
                * ((self.prm.alpha_2 / W) ** ((1 - self.prm.alpha_1) / (1 - self.prm.alpha_1 - self.prm.alpha_2)))
        
        S_F = (self.prm.beta_1 / R) ** ((1 - self.prm.beta_2) / (1 - self.prm.beta_1 - self.prm.beta_2)) \
            * (self.prm.beta_2 / W) ** (self.prm.beta_2 / (1 - self.prm.beta_1 - self.prm.beta_2))
        
        L_F = (self.prm.beta_1 / R) ** (self.prm.beta_1 / (1 - self.prm.beta_1 - self.prm.beta_2)) \
            * (self.prm.beta_2 / W) ** ((1 - self.prm.beta_1) / (1 - self.prm.beta_1 - self.prm.beta_2))
        
        return S_H, L_H, S_F, L_F

    def pi(self, R, W, m):
        '''利潤関数を計算する関数

        Parameters
        ----------
        RW: numpy.ndarray (K, )
            地代と賃金の結合リスト

        Returns
        -------
        pi : numpy.ndarray (K, )
            利潤関数値
        '''

        pi = np.dot(self.prm.D, m)\
           + (1 - self.prm.beta_1 - self.prm.beta_2)\
           * ((self.prm.beta_1 / R) ** (self.prm.beta_1 / (1 - self.prm.beta_1 - self.prm.beta_2)))\
           * ((self.prm.beta_2 / W) ** (self.prm.beta_2 / (1 - self.prm.beta_1 - self.prm.beta_2)))
        
        return pi
    
    def v(self, R, W):
        '''間接効用関数を計算する関数

        Parameters
        ----------
        RW: numpy.ndarray (K, )
            地代と賃金の結合リスト

        Returns
        -------
        v : numpy.ndarray (K, K)
            間接効用関数値
        '''
        
        v = self.prm.E * W \
          + (1 - self.prm.alpha_1 - self.prm.alpha_2)\
          * (1 / self.prm.T) ** (1 / (1 - self.prm.alpha_1 - self.prm.alpha_2))\
          * ((self.prm.alpha_1 / R) ** (self.prm.alpha_1 / (1 - self.prm.alpha_1 - self.prm.alpha_2)))[:, np.newaxis]\
          * ((self.prm.alpha_2 / W) ** (self.prm.alpha_2 / (1 - self.prm.alpha_1 - self.prm.alpha_2)))

        return v
    
    def v_flat(self, R, W):
        '''間接効用関数を計算する関数

        Parameters
        ----------
        RW: numpy.ndarray (K, )
            地代と賃金の結合リスト

        Returns
        -------
        v : numpy.ndarray (K, K)
            間接効用関数値
        '''
        
        v = self.prm.E * W \
          + (1 - self.prm.alpha_1 - self.prm.alpha_2)\
          * (1 / self.prm.T) ** (1 / (1 - self.prm.alpha_1 - self.prm.alpha_2))\
          * ((self.prm.alpha_1 / R) ** (self.prm.alpha_1 / (1 - self.prm.alpha_1 - self.prm.alpha_2)))[:, np.newaxis]\
          * ((self.prm.alpha_2 / W) ** (self.prm.alpha_2 / (1 - self.prm.alpha_1 - self.prm.alpha_2)))

        return v.flatten()

    def Z_SD(self, RW, m, n):
        '''目的関数を計算する関数

        Parameters
        ----------
        RW: numpy.ndarray (K, )
            地代と賃金の結合リスト
        m: numpy.ndarray (K, )
            企業分布
        n: numpy.ndarray (K, K)
            家計分布

        Returns
        -------
        Z_SD : float
            短期均衡の目的関数値
        '''

        R = RW[:self.prm.K]
        W = RW[self.prm.K:]
        
        # pi_noex = \
        #         + (1 - self.prm.beta_1 - self.prm.beta_2)\
        #         * ((self.prm.beta_1 / R) ** (self.prm.beta_1 / (1 - self.prm.beta_1 - self.prm.beta_2)))\
        #         * ((self.prm.beta_2 / W) ** (self.prm.beta_2 / (1 - self.prm.beta_1 - self.prm.beta_2)))
        
        # F_value = np.sum(self.v(R, W) * n)\
        #         + np.sum((self.pi(R, W, m) - np.dot(self.prm.D, m)) * m)\
        #         + self.prm.S_bar * np.sum(R)
        
        #flatten()遅い！！
        F_value = self.v_flat(R, W) @ n.flatten()\
                + (self.pi(R, W, m) - (self.prm.D @ m)) @ m\
                + self.prm.S_bar * np.sum(R)
        
        # print("v.flatten():", self.v(R, W).flatten())
        
        # F_value = np.sum(self.v(R, W) * n)\
        #         + pi_noex @ m\
        #         + self.prm.S_bar * np.sum(R)
        
        # F_value = np.dot(self.v(R, W).flatten(), n.flatten())\
        #         + np.sum((self.pi(R, W, m) - np.dot(self.prm.D, m)) * m)\
        #         + self.prm.S_bar * np.sum(R)
        
        return F_value
    
    def short_dual_df(self, RW, m, n):
        '''目的関数の勾配を計算する関数

        Parameters
        ----------
        R: numpy.ndarray (K, )
            地代
        W: numpy.ndarray (K, )
            賃金
        m: numpy.ndarray (K, )
            企業分布
        n: numpy.ndarray (K, K)
            家計分布

        Returns
        -------
        dRW : numpy.ndarray (2 * K, )
            短期均衡問題の勾配
        '''
        
        R = RW[:self.prm.K]
        W = RW[self.prm.K:]
        
        S_H, L_H, S_F, L_F = self.demand(R, W)
        
        dR = self.prm.S_bar - np.sum(S_H * n, axis=1) - S_F * m  # forall i
        dW = self.prm.E * np.sum(n, axis=0) - np.sum(L_H * n, axis=0) - L_F * m # forall j

        #ベクトルくっつけるの遅い可能性・・・
        return np.concatenate((dR, dW))

    def short_solve(self, RW_ini, m_fixed, n_fixed, err_short, short_itr, rel):
        '''短期均衡を解く関数

        Parameters
        ----------
        RW0: numpy.ndarray(2 * K, )
            価格変数の初期値
        m: numpy.ndarray(K, )
            長期均衡から与えられる企業分布
        n: numpy.ndarray(K, K)
            長期均衡から与えられる家計分布
        max_itr: int
            最大反復回数
        err_short: float
            収束判定の閾値
        Returns
        -------
        R: numpy.ndarray (K, )
            地代
        W: numpy.ndarray (K, )
            賃金
        '''
        #Step.0 初期値を決める

        RW_before = RW_ini * np.ones(2 * self.prm.K, dtype=np.float64)
        
        p_bar_before = np.ones(2 * self.prm.K, dtype=np.float64)
        
        L_before = self.algprm.L
        
        t_before = 1
        
        max_value = 0
        obj_hist = {}
        obj_hist[0] = 1000
        
        iteration = []
        obj_rel_list = []
        
        R_ini = RW_before[:self.prm.K]
        W_ini = RW_before[self.prm.K:]
        
        if np.min(self.v(R_ini, W_ini)) < 0:
            raise ValueError("v must be non-negative")

        for k in range(1, short_itr):
            iteration.append(k)
            
            #Step.1 Backtracking
            i = self.backtracking(p_bar_before, L_before, m_fixed, n_fixed)
            
            #Lは 1/Lで扱った方が計算回数は少ない
            L = (self.algprm.eta ** i) * L_before
            
            dZ_SD = self.short_dual_df(p_bar_before, m_fixed, n_fixed)
            
            #step.2 解の更新
            RW = np.maximum(self.prm.RW_proj, p_bar_before - (dZ_SD / L), dtype=np.float64)
            
            #Step.3 収束判定
            if np.max(abs((RW - RW_before) / RW_before)) < err_short:
                break
            # obj = self.Z_SD(RW, m_fixed, n_fixed)
            # obj_rel = abs(obj - obj_corr) / abs(obj_corr)
            # obj_rel_list.append(obj_rel)
            
            # if obj_rel < 1e-4:
            #     break
            
            # if k > 1:
            #     max_obj = max(obj_hist.values())
            #     if ((obj_hist[k-1] - obj_hist[k]) / np.maximum(max_obj, 1)) < err_short:
            #         break
            
            #Step.4 Adaptive restart
            if np.dot(self.short_dual_df(RW_before, m_fixed, n_fixed), (RW - RW_before)) > 0:
                t_before = 1
            
            #step.5 momentum項の計算
            t = (1 + np.sqrt(1 + 4 * (t_before ** 2))) * 0.5
            
            p_bar = np.maximum(self.prm.RW_proj, RW + ((t_before - 1) / t) * (RW - RW_before), dtype=np.float64)
                
            max_value = k
            
            RW_before = RW
            p_bar_before = p_bar
            L_before = L
            t_before = t
            
        print("short_max_value:", max_value)
        
        R_hist = RW[:self.prm.K]
        W_hist = RW[self.prm.K:]
        
        RW_hist = np.concatenate((R_hist, W_hist))
        
        equ_R, equ_W = self.equilibrium(R_hist, W_hist, m_fixed, n_fixed)
        
        return R_hist, W_hist, iteration, obj_rel_list
    
    def Q(self, p, p_bar, L_bar, m_fixed, n_fixed):
        
        #np.linalg.normではなく，内積の方が速いかも．
        Q = self.Z_SD(p_bar, m_fixed, n_fixed) + np.dot(p - p_bar, self.short_dual_df(p_bar, m_fixed, n_fixed)) + (L_bar / 2) * np.linalg.norm(p - p_bar) ** 2
        
        # Q_2 = self.Z_SD(p_bar) - (1 / (2 * L_bar)) * (np.linalg.norm(self.short_dual_df(p_bar)) ** 2)
        
        return Q

    def backtracking(self, p_bar, L, m_fixed, n_fixed):
        
        i = 0
        
        while True:
            L_bar = (self.algprm.eta ** i) * L
            p = np.maximum(self.algprm.p_proj, p_bar - (self.short_dual_df(p_bar, m_fixed, n_fixed) / L_bar))
            if self.Z_SD(p, m_fixed, n_fixed) <= self.Q(p, p_bar, L_bar, m_fixed, n_fixed):
                break
            
            i = i + 1
            
        return i
    
    def equilibrium(self, R, W, m_fixed, n_fixed):
        
        RW = np.concatenate((R, W))
        
        dRW = self.short_dual_df(RW, m_fixed, n_fixed)
        
        dR = dRW[:self.prm.K]
        dW = dRW[self.prm.K:]
        
        equ_R = np.max(R * dR)
        equ_W = np.max(W * dW)
        
        return equ_R, equ_W