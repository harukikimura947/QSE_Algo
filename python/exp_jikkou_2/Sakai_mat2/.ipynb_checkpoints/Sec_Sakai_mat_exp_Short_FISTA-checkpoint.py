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

from scipy.spatial import distance
from scipy.optimize import linprog
from scipy.optimize import minimize
from collections import defaultdict
from mpl_toolkits.mplot3d import Axes3D
from numba import njit

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
        
        self.T_power_v = (1 - self.prm.alpha_1 - self.prm.alpha_2) * self.T_power
        self.pi_beta = 1 - self.prm.beta_1 - self.prm.beta_2

        #Rのためのsparse行列
        # Z_I = np.zeros((self.prm.K, self.prm.K * self.prm.K))

        #Wのためのsparse行列
        # Z_J = np.zeros((self.prm.K, self.prm.K * self.prm.K))
        
#         for i in range(self.prm.K):
#             for j in range(i * self.prm.K, i * self.prm.K + self.prm.K):
#                 Z_I[i][j] = 1.0
                
#         for i in range(self.prm.K):
#             for j in range(self.prm.K):
#                 Z_J[i][j * self.prm.K + i] = 1.0
        
#         rows_I = np.arange(self.prm.K)
#         cols_I = np.arange(self.prm.K) * self.prm.K
#         Z_I[rows_I[:, None], cols_I[:, None] + np.arange(self.prm.K)] = 1.0
        
#         rows_J = np.arange(self.prm.K)[:, None] # 行インデックス
#         cols_J = np.arange(self.prm.K) * self.prm.K + np.arange(self.prm.K)[:, None]  # 列インデックス
#         Z_J[rows_J, cols_J] = 1.0

        self.one_K = np.ones(self.prm.K, dtype=np.float64)
        self.one_2K = np.ones(2 * self.prm.K, dtype=np.float64)

    def demand_sparse(self, R, W):
        '''需要関数を設定する関数

        Parameters
        ----------
        R: numpy.ndarray(matrix) (K, )
            地代リスト
        W: numpy.ndarray(matrix) (K, )
            賃金リスト
        '''


        # 外積案 np.outer
        S_H = self.T_power\
            * np.outer((self.alpha_S_R * np.power(R, -self.power_S_H_R)), (self.alpha_S_W * np.power(W, -self.power_S_H_W)))
        
        L_H = self.T_power\
            * np.outer((self.alpha_L_R * np.power(R, -self.power_L_H_R)), (self.alpha_L_W * np.power(W, -self.power_L_H_W)))

        # S_F = self.beta_S_R * ((R_inv) ** self.power_S_F_R) \
        #     * self.beta_S_W * ((W_inv) ** self.power_S_F_W)
        
        # L_F = self.beta_L_R * ((R_inv) ** self.power_L_F_R) \
        #     * self.beta_L_W * ((W_inv) ** self.power_L_F_W)
        
        S_F = self.beta_S_R * (np.power(R, - self.power_S_F_R)) \
            * self.beta_S_W * (np.power(W, - self.power_S_F_W))
        L_F = self.beta_L_R * (np.power(R, - self.power_L_F_R)) \
            * self.beta_L_W * (np.power(W, - self.power_L_F_W))

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
        
        pi = self.prm.D@m \
           + (1 - self.prm.beta_1 - self.prm.beta_2)\
           * (self.beta_L_R * np.power(R, -self.power_L_F_R))\
           * (self.beta_S_W * np.power(W, -self.power_S_F_W))
        
        # コード確認用
        # pi = np.dot(self.prm.D, m)\
        #    + (1 - self.prm.beta_1 - self.prm.beta_2)\
        #    * (self.beta_L_R * (1.0 / R) ** self.power_L_F_R)\ # (self.prm.beta_1 / (1 - self.prm.beta_1 - self.prm.beta_2))
        #    * (self.beta_S_W * (1.0 / W) ** self.power_S_F_W)  # (self.prm.beta_2 / (1 - self.prm.beta_1 - self.prm.beta_2))
        
        return pi
    
    def pi_noex(self, R, W, m):
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
        
        pi = self.pi_beta\
           * (self.beta_L_R * np.power(R, -self.power_L_F_R))\
           * (self.beta_S_W * np.power(W, -self.power_S_F_W))
        
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

        v = self.prm.E * W \
          + self.T_power_v\
          * np.outer((self.alpha_L_R * np.power(R, -self.power_L_H_R)), (self.alpha_S_W * np.power(W, -self.power_S_H_W)))
        
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
        
        
        #flatten()はコピーを生成するので，ravel()よりメモリを食う．
        F_value = self.v(R, W).ravel() @ n.ravel()\
                + self.pi_noex(R, W, m) @ m\
                + self.prm.S_bar * np.sum(R)
        
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
        n: numpy.ndarray (K, K) ???
            家計分布

        Returns
        -------
        dRW : numpy.ndarray (2 * K, )
            短期均衡問題の勾配
        '''

        R = RW[:self.prm.K]
        W = RW[self.prm.K:]
        
        S_H, L_H, S_F, L_F = self.demand_sparse(R, W)
        
        dR = self.prm.S_bar - self.sum_n1(S_H * n, axis=1) - S_F * m  # forall i
        dW = self.sum_n0((self.prm.E-L_H) * n, axis=0) - L_F * m # forall j
        # dW = self.prm.E * self.sum_n0(n, axis=0) - self.sum_n0(L_H * n, axis=0) - L_F * m # forall j

        #ベクトルくっつけるの遅い可能性・・・
        return np.concatenate((dR, dW))
    
    # https://www.procrasist.com/entry/einsum#google_vignette
    # 書き換えやすいようにaxis引数は残した
    def sum_n0(self, n, axis):
        return np.einsum('ij->j', n)
    
    def sum_n1(self, n, axis):
        return np.einsum('ij->i', n)

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

        RW_before = RW_ini * self.one_2K
        
        p_bar_before = self.one_2K
        
        L_before = self.algprm.L
        
        t_before = 1
        
        max_value = 0
        obj_hist = {}
        obj_hist[0] = 1000
        
        iteration = []
#         R_err = []
#         R_list = []
#         W_err = []
#         W_list = []
#         v_err = []
#         v_list = []
#         pi_err = []
#         pi_list = []
#         obj_err = []
#         Z_list = []
#         grad_mean = []
        
#         obj_list = []
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
            L = self.backtracking(p_bar_before, L_before, m_fixed, n_fixed)
            
            dZ_SD = self.short_dual_df(p_bar_before, m_fixed, n_fixed)
            
            #step.2 解の更新
            # RW = np.maximum(self.prm.RW_proj, p_bar_before - (dZ_SD / L), dtype=np.float64)
            RW = np.maximum(self.prm.RW_proj, p_bar_before - (dZ_SD / L))
            
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
            if self.short_dual_df(RW_before, m_fixed, n_fixed)@(RW - RW_before) > 0:
                t_before = 1
            
            #step.5 momentum項の計算
            t = (1.0 + np.sqrt(1.0 + 4.0 * (np.power(t_before, 2)))) * 0.5
            
            # p_bar = np.maximum(self.prm.RW_proj, RW + ((t_before - 1.0) / t) * (RW - RW_before), dtype=np.float64)
            p_bar = np.maximum(self.prm.RW_proj, RW + ((t_before - 1.0) / t) * (RW - RW_before))
            
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
        
        Q = self.Z_SD(p_bar, m_fixed, n_fixed) + (p - p_bar)@self.short_dual_df(p_bar, m_fixed, n_fixed) + (0.5*L_bar) * (p - p_bar)@(p - p_bar)
                
        return Q
    
    
    def backtracking(self, p_bar, L, m_fixed, n_fixed):
        
        i = 0
        L_bar = L
        dRdW = self.short_dual_df(p_bar, m_fixed, n_fixed)
        Z_SD_p_bar = self.Z_SD(p_bar, m_fixed, n_fixed)
        m_flat = m_fixed.ravel()
        n_flat = n_fixed.ravel()
        
        # while True:
        for k in range(10**6):
            p = np.maximum(self.algprm.p_proj, p_bar - dRdW / L_bar)

            # 重複したZ_SDとQの内部で重複した計算が多い気がする．
            # Qの計算の中身では，このwhile文直前のdRdWを使えるはず
            # そもそもif文で評価するZ_SD-Q自体が解析的にもっと整理できそう（self.Qの中でself.Z_SDで呼んでいるので，例えば(self.prm.D @ m)) @ mとかは，少なくともこのif文の評価では計算しなくてよい）
            # その上でpとL_barに依存しない部分は，while文の外で計算しておくのが良さそう

            # if self.Z_SD(p, m_fixed, n_fixed) -self.Q(p, p_bar, L_bar, m_fixed, n_fixed) <= 0.0:
            #     break
            
            p_R = p[:self.prm.K]
            p_W = p[self.prm.K:]
            
            v = self.prm.E * p_W \
              + self.T_power_v\
              * np.outer((self.alpha_L_R * np.power(p_R, -self.power_L_H_R)), (self.alpha_S_W * np.power(p_W, -self.power_S_H_W)))
            
            pi = self.pi_beta\
               * (self.beta_L_R * np.power(p_R, -self.power_L_F_R))\
               * (self.beta_S_W * np.power(p_W, -self.power_S_F_W))
            
            Z = v.ravel() @ n_flat\
              + pi.ravel() @ m_flat\
              + self.prm.S_bar * np.sum(p_R)
            
            # ちょっとだけ修正：
            # if self.Z_SD(p, m_fixed, n_fixed) - (Z_SD_p_bar +  (p - p_bar)@dRdW + (0.5*L_bar) * (p - p_bar)@(p - p_bar) )<= 0.0:
            if Z - (Z_SD_p_bar +  (p - p_bar)@dRdW + (0.5*L_bar) * (p - p_bar)@(p - p_bar) )<= 0.0:
                break
            
            L_bar = L_bar*self.algprm.eta
            i += 1

        return L_bar
    
    def equilibrium(self, R, W, m_fixed, n_fixed):
        
        RW = np.concatenate((R, W))
        
        dRW = self.short_dual_df(RW, m_fixed, n_fixed)
        
        dR = dRW[:self.prm.K]
        dW = dRW[self.prm.K:]
        
        equ_R = np.max(R * dR)
        equ_W = np.max(W * dW)
        
        # G = np.sum(R * dR) + np.sum(W * dW)
        
        return equ_R, equ_W