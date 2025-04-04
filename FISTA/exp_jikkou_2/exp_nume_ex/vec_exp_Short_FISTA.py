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
# import btra2
# import gradient
# import demand

from scipy.sparse import csr_matrix
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
        
        self.power_S_H_R = ((1 - self.prm.alpha_2) / (1 - self.prm.alpha_1 - self.prm.alpha_2))
        self.power_S_H_W = (self.prm.alpha_2 / (1 - self.prm.alpha_1 - self.prm.alpha_2))
        self.power_L_H_R = (self.prm.alpha_1 / (1 - self.prm.alpha_1 - self.prm.alpha_2))
        self.power_L_H_W = ((1 - self.prm.alpha_1) / (1 - self.prm.alpha_1 - self.prm.alpha_2))
        
        self.power_S_F_R = ((1 - self.prm.beta_2) / (1 - self.prm.beta_1 - self.prm.beta_2))
        self.power_S_F_W = (self.prm.beta_2 / (1 - self.prm.beta_1 - self.prm.beta_2))
        self.power_L_F_R = (self.prm.beta_1 / (1 - self.prm.beta_1 - self.prm.beta_2))
        self.power_L_F_W = ((1 - self.prm.beta_1) / (1 - self.prm.beta_1 - self.prm.beta_2))
        
        self.T_power = (1 / self.prm.T) ** (1 / (1 - self.prm.alpha_1 - self.prm.alpha_2))
        self.alpha_S_R = self.prm.alpha_1 ** self.power_S_H_R
        self.alpha_S_W = self.prm.alpha_2 ** self.power_S_H_W
        self.alpha_L_R = self.prm.alpha_1 ** self.power_L_H_R
        self.alpha_L_W = self.prm.alpha_2 ** self.power_L_H_W
        
        self.beta_S_R = self.prm.beta_1 ** self.power_S_F_R
        self.beta_S_W = self.prm.beta_2 ** self.power_S_F_W
        self.beta_L_R = self.prm.beta_1 ** self.power_L_F_R
        self.beta_L_W = self.prm.beta_2 ** self.power_L_F_W
        
        Z_I = np.zeros((self.prm.K, self.prm.K * self.prm.K))
        Z_J = np.zeros((self.prm.K, self.prm.K * self.prm.K))
        
#         #Rのためのsparse行列
#         for i in range(self.prm.K):
#             for j in range(i * self.prm.K, i * self.prm.K + self.prm.K):
#                 Z_I[i][j] = 1.0
                
#         #Wのためのsparse行列
#         for i in range(self.prm.K):
#             for j in range(self.prm.K):
#                 Z_J[i][j * self.prm.K + i] = 1.0
        
        rows_I = np.arange(self.prm.K)
        cols_I = np.arange(self.prm.K) * self.prm.K
        Z_I[rows_I[:, None], cols_I[:, None] + np.arange(self.prm.K)] = 1.0
        
        rows_J = np.arange(self.prm.K)[:, None] # 行インデックス
        cols_J = np.arange(self.prm.K) * self.prm.K + np.arange(self.prm.K)[:, None]  # 列インデックス
        Z_J[rows_J, cols_J] = 1.0
        
        self.sparse_I = csr_matrix(Z_I)
        self.sparse_J = csr_matrix(Z_J)
    
    def demand_sparse(self, R, W):
        '''需要関数を設定する関数

        Parameters
        ----------
        R: numpy.ndarray(matrix) (K, )
            地代リスト
        W: numpy.ndarray(matrix) (K, )
            賃金リスト
        '''
        #同じような形が何回も現れているので，それは一番最初に求めておく．
        
        S_H = self.T_power\
            * ((self.alpha_S_R * (1.0 / R) ** self.power_S_H_R) @ self.sparse_I)\
            * ((self.alpha_S_W * (1.0 / W) ** self.power_S_H_W) @ self.sparse_J)
        
        L_H = self.T_power\
            * ((self.alpha_L_R * (1.0 / R) ** self.power_L_H_R) @ self.sparse_I)\
            * ((self.alpha_L_W * (1.0 / W) ** self.power_L_H_W) @ self.sparse_J)
        
        S_F = self.beta_S_R * (1.0 / R) ** self.power_S_F_R \
            * self.beta_S_W * (1.0 / W) ** self.power_S_F_W
        
        L_F = self.beta_L_R * (1.0 / R) ** self.power_L_F_R \
            * self.beta_L_W * (1.0 / W) ** self.power_L_F_W
        
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
           * (self.beta_L_R * (1.0 / R) ** self.power_L_F_R)\
           * (self.beta_S_W * (1.0 / W) ** self.power_S_F_W)
        
        # コード確認用
        # pi = np.dot(self.prm.D, m)\
        #    + (1 - self.prm.beta_1 - self.prm.beta_2)\
        #    * (self.beta_L_R * (1.0 / R) ** self.power_L_F_R)\ # (self.prm.beta_1 / (1 - self.prm.beta_1 - self.prm.beta_2))
        #    * (self.beta_S_W * (1.0 / W) ** self.power_S_F_W)  # (self.prm.beta_2 / (1 - self.prm.beta_1 - self.prm.beta_2))
        
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
        
        v = self.prm.E * (W @ self.sparse_J) \
          + (1 - self.prm.alpha_1 - self.prm.alpha_2)\
          * self.T_power\
          * ((self.alpha_L_R * (1.0 / R) ** self.power_L_H_R) @ self.sparse_I)\
          * ((self.alpha_S_W * (1.0 / W) ** self.power_S_H_W) @ self.sparse_J)
        
        # コード確認用
        # v = self.prm.E * W \
        #   + (1 - self.prm.alpha_1 - self.prm.alpha_2)\
        #   * self.T_power\
        #   * (self.alpha_L_R * (1.0 / R) ** self.power_L_H_R)[:, np.newaxis]\ #  (self.prm.alpha_1 / (1 - self.prm.alpha_1 - self.prm.alpha_2))
        #   * (self.alpha_S_W * (1.0 / W) ** self.power_S_H_W)                 #  (self.prm.alpha_2 / (1 - self.prm.alpha_1 - self.prm.alpha_2))

        return v

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
        F_value = self.v(R, W) @ n\
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
        
        S_H, L_H, S_F, L_F = self.demand_sparse(R, W)
        
        dR = self.prm.S_bar - (self.sparse_I @ (S_H * n)) - S_F * m  # forall i
        dW = self.prm.E * (self.sparse_J @ n) - (self.sparse_J @ (L_H * n)) - L_F * m # forall j

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
        R_err = []
        R_list = []
        W_err = []
        W_list = []
        v_err = []
        v_list = []
        pi_err = []
        pi_list = []
        obj_err = []
        Z_list = []
        grad_mean = []
        
        obj_list = []
        obj_rel_list = []
        
        R_ini = RW_before[:self.prm.K]
        W_ini = RW_before[self.prm.K:]
        
        if np.min(self.v(R_ini, W_ini)) < 0:
            raise ValueError("v must be non-negative")

        for k in range(1, short_itr):
            # print("="*40)
            # print("kshort:", max_value)
            iteration.append(k)
            
            #Step.1 Backtracking
            i = self.backtracking(p_bar_before, L_before, m_fixed, n_fixed)
            # i = btra2.backtracking(L_before, p_bar_before, m_fixed, n_fixed)
            
            #Lは 1/Lで扱った方が計算回数は少ない
            L = (self.algprm.eta ** i) * L_before
            
            dZ_SD = self.short_dual_df(p_bar_before, m_fixed, n_fixed)
            
            #step.2 解の更新
            RW = np.maximum(self.prm.RW_proj, p_bar_before - (dZ_SD / L), dtype=np.float64)
            
            #Step.3 収束判定
            # obj = self.Z_SD(RW, m_fixed, n_fixed)
            # obj_rel = abs(obj - obj_corr) / abs(obj_corr)
            # obj_rel_list.append(obj_rel)
            
            # if obj_rel < 1e-4:
            #     break
            
            # if k > 1:
            #     max_obj = max(obj_hist.values())
            #     if ((obj_hist[k-1] - obj_hist[k]) / np.maximum(max_obj, 1)) < err_short:
            #         break
                
            if np.max(abs((RW - RW_before) / RW_before)) < err_short:
                break
            
            #Step.4 Adaptive restart
            if np.dot(self.short_dual_df(RW_before, m_fixed, n_fixed), (RW - RW_before)) > 0:
                t_before = 1
            
            #step.5 momentum項の計算
            t = (1 + np.sqrt(1 + 4 * (t_before ** 2))) * 0.5
            
            p_bar = np.maximum(self.prm.RW_proj, RW + ((t_before - 1) / t) * (RW - RW_before), dtype=np.float64)
            
#             if rel == 1:
#                 R = RW[:self.prm.K]
#                 W = RW[self.prm.K:]
                
#                 v = self.v(R, W)
#                 pi = self.pi(R, W, m_fixed)
                
#                 R_list.append(np.mean(R))
#                 W_list.append(np.mean(W))
#                 pi_list.append(np.mean(pi))
#                 v_list.append(np.mean(v))
#                 Z_list.append(self.Z_SD(RW, m_fixed, n_fixed))
                
#                 iteration.append(k)
#                 R_err.append(np.mean(abs((R - R_corr) / R_corr)))
#                 W_err.append(np.mean(abs((W - W_corr) / W_corr)))
#                 v_err.append(np.mean(abs((v - v_corr) / v_corr)))
#                 pi_err.append(np.mean(abs((pi - pi_corr) / pi_corr)))
#                 obj_err.append(np.mean(abs((self.Z_SD(RW, m_fixed, n_fixed) - obj_corr) / obj_corr)))
#                 grad_mean.append(np.mean(self.short_dual_df(RW, m_fixed, n_fixed)))
                
            max_value = k
            
            RW_before = RW
            p_bar_before = p_bar
            L_before = L
            t_before = t
            
        # print("short_max_value:", max_value)
        
        R_hist = RW[:self.prm.K]
        W_hist = RW[self.prm.K:]
        
        RW_hist = np.concatenate((R_hist, W_hist))
        
        equ_R, equ_W = self.equilibrium(R_hist, W_hist, m_fixed, n_fixed)
        
        # print("jac:", np.mean(self.short_dual_df(RW_hist, m_fixed, n_fixed)))
        # print("equ_R:", equ_R)
        # print("equ_W:", equ_W)
        
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
    
#     def backtracking(self, p_bar, L, m_fixed, n_fixed):
        
#         i = 0
#         # L_bar = L
        
#         while True:
#             # L_bar = (self.algprm.eta ** i) * L
#             p = np.maximum(self.algprm.p_proj, p_bar - (self.short_dual_df(p_bar, m_fixed, n_fixed) / L_bar))
#             if self.Z_SD(p, m_fixed, n_fixed) <= self.Q(p, p_bar, L_bar, m_fixed, n_fixed):
#                 break
            
#             L_bar = self.algprm.eta * L_bar
            
#        return i 後で直す
    
    def equilibrium(self, R, W, m_fixed, n_fixed):
        
        RW = np.concatenate((R, W))
        
        dRW = self.short_dual_df(RW, m_fixed, n_fixed)
        
        dR = dRW[:self.prm.K]
        dW = dRW[self.prm.K:]
        
        equ_R = np.max(R * dR)
        equ_W = np.max(W * dW)
        
        # G = np.sum(R * dR) + np.sum(W * dW)
        
        return equ_R, equ_W