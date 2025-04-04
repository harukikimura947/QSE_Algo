"""Module summary.

This is Short-Algo.py.
"""

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

from scipy.spatial import distance
from scipy.optimize import linprog
from scipy.optimize import minimize
from collections import defaultdict
from mpl_toolkits.mplot3d import Axes3D
from numba import njit
# from typing import Tuple
# from line_profiler import LineProfiler
# profile = LineProfiler()

short_solve = []
short_solve_lst = []
short_iteration = []
Lipchitz = []

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
        
        self.T = t*distance_matrix
        self.Ker = np.exp(-self.theta_house*self.T)
        self.D = np.exp(-tau*distance_matrix)
        
class Algo_Parameter:
    '''アルゴリズム中のパラメータを保持するクラス

    Attributes
    ----------
    L: float
            FISTAのパラメータ
    eta: float
        FISTAのパラメータ
    p_proj:
        パラメータpの射影後の値
    '''

    def __init__(self, L, eta, p_proj):
        '''

        Parameters
        ----------
        L: float
            FISTAのパラメータ
        eta: float
            FISTAのパラメータ
        p_proj:
            パラメータpの射影後の値
        '''
        self.L = L
        self.eta = eta
        self.p_proj = p_proj

class Short:
    '''短期均衡問題クラス

    Attributes
    ----------
    prm: Parameter
       パラメータクラス 
    '''

    def __init__(self, prm, algprm, m, n):
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
        self.m = m
        self.n = n
    
    def demand(self, R, W):
        '''需要関数を設定する関数

        Parameters
        ----------
        R: numpy.ndarray(matrix) (K, )
            地代リスト
        W: numpy.ndarray(matrix) (K, )
            賃金リスト
        '''
        
        #最小距離の半分の長さ
        
        S_H = np.maximum(0, self.prm.alpha_1 / R)
        L_H = np.maximum(0, self.prm.alpha_2 / W)
        S_F = np.maximum(0, self.prm.beta_1 / R)
        L_F = np.maximum(0, self.prm.beta_2 / W)
        
        return S_H, L_H, S_F, L_F

    def pi(self, R, W):
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

        pi = np.maximum(0.1, np.dot(self.prm.D, self.m)\
        - self.prm.beta_1 * np.log(R) - self.prm.beta_2 * np.log(W)\
        + self.prm.beta_1 * np.log(self.prm.beta_1) + self.prm.beta_2 * np.log(self.prm.beta_2)\
        - self.prm.beta_1 - self.prm.beta_2)
        
#         pi = np.dot(self.prm.D, m)\
#         - self.prm.beta_1 * np.log(R) - self.prm.beta_2 * np.log(W)\
#         + self.prm.beta_1 * np.log(self.prm.beta_1) + self.prm.beta_2 * np.log(self.prm.beta_2)\
#         - self.prm.beta_1 - self.prm.beta_2
        
#         if np.min(pi) < 0:
#             print("end")
#             break
        
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

        v = np.maximum(0.1, self.prm.E * W * np.ones((self.prm.K, self.prm.K))\
        - (self.prm.alpha_1 * np.log(R) * np.ones((self.prm.K, self.prm.K))).T\
        - self.prm.alpha_2 * np.log(W) * np.ones((self.prm.K, self.prm.K))\
        + self.prm.alpha_1 * np.log(self.prm.alpha_1) + self.prm.alpha_2 * np.log(self.prm.alpha_2)\
        - self.prm.alpha_1 - self.prm.alpha_2 - self.prm.T)
        
#         v = self.prm.E * W * np.ones((self.prm.K, self.prm.K))\
#         - (self.prm.alpha_1 * np.log(R) * np.ones((self.prm.K, self.prm.K))).T\
#         - self.prm.alpha_2 * np.log(W) * np.ones((self.prm.K, self.prm.K))\
#         + self.prm.alpha_1 * np.log(self.prm.alpha_1) + self.prm.alpha_2 * np.log(self.prm.alpha_2)\
#         - self.prm.alpha_1 - self.prm.alpha_2 - self.prm.T
        
#         if np.min(v) < 0:
#             print("end")
#             break

        return v

    def Z_SD(self, RW):
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
        
        pi_noex = np.maximum(0.1, 
        - self.prm.beta_1 * np.log(R) - self.prm.beta_2 * np.log(W)\
        + self.prm.beta_1 * np.log(self.prm.beta_1) + self.prm.beta_2 * np.log(self.prm.beta_2)\
        - self.prm.beta_1 - self.prm.beta_2)
        
        F_value = np.sum(self.v(R, W) * self.n)\
        + np.sum(pi_noex * self.m)\
        + self.prm.S_bar * np.sum(R)
        
        return F_value
    
    def short_dual_df(self, RW):
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
        
        # dR = np.maximum(1e-3, self.prm.S_bar - S_H * np.sum(n, axis=1) - S_F * m)
        # dW = np.maximum(1e-3, self.prm.E * np.sum(n, axis=0) - L_H * np.sum(n, axis=0) - L_F * m)
        
        dR = self.prm.S_bar - S_H * np.sum(self.n, axis=1) - S_F * self.m
        dW = self.prm.E * np.sum(self.n, axis=0) - L_H * np.sum(self.n, axis=0) - L_F * self.m

        return np.concatenate((dR, dW))

    def short_solve(self, RW_ini, R_corr, W_corr, pi_corr, v_corr, obj_corr,
                    short_itr, err_short, rel):
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
        RW_before = RW_ini * np.ones(2 * self.prm.K)

        B_before = np.eye(len(RW_before))
        
        max_value = 0
        
        iteration = []
        R_err = []
        W_err = []
        v_err = []
        pi_err = []
        obj_err = []
        grad_mean = []

        for k in range(1, short_itr):
            
            #勾配計算
            dRW_before = self.short_dual_df(RW_before)

            #Step.1 探索方向を決める
            B_inv_before = np.linalg.inv(B_before)
            RW_d = -np.dot(B_inv_before, dRW_before)

            #Step.2 ステップサイズを決める
            # alpha = 1 / k
            alpha = self.short_armijo(RW_before, RW_d)

            #Step.3 解を更新する
            RW = np.maximum(1e-3, RW_before + alpha * RW_d)
            
            if rel == 1:
                R = RW[:self.prm.K]
                W = RW[self.prm.K:]
                
                v = self.v(R, W)
                pi = self.pi(R, W)
                
                iteration.append(k)
                R_err.append(np.mean(abs((R - R_corr) / R_corr)))
                W_err.append(np.mean(abs((W - W_corr) / W_corr)))
                v_err.append(np.mean(abs((v - v_corr) / v_corr)))
                pi_err.append(np.mean(abs((pi - pi_corr) / pi_corr)))
                obj_err.append(np.mean(abs((self.Z_SD(RW) - obj_corr) / obj_corr)))
                grad_mean.append(np.mean(self.short_dual_df(RW)))
            
            #勾配計算
            dRW = self.short_dual_df(RW)

            # BFGS更新式
            s = RW - RW_before
            y = dRW - dRW_before
            B = self.BFGS(B_before, s, y)
            
            #Step.4 収束判定
            if np.max(abs((RW - RW_before) / RW_before)) < err_short:
                break
    
            max_value = k
        
            RW_before = RW
            B_before = B
            
        print("short_max_value:", max_value)
        
        R_hist = RW[:self.prm.K]
        W_hist = RW[self.prm.K:]
        
        S_H, L_H, S_F, L_F = self.demand(R_hist, W_hist)
        
        print("Lyapunov:", np.sum(R_hist * (self.prm.S_bar - S_H * np.sum(self.n, axis=1) - S_F * self.m))\
              + np.sum(W_hist * (self.prm.E * np.sum(self.n, axis=0) - L_H * np.sum(self.n, axis=0) - L_F * self.m)))
        
        return R_hist, W_hist, iteration, R_err, W_err, v_err, pi_err, obj_err, grad_mean

    def short_armijo(self, RW, RW_d, alpha=0.8, beta=0.8):
        """
        Armijoのルール

        Parameters
        ----------
        grad: function
            勾配ベクトルを計算する関数
        RW: numpy.ndarray (2 * K, )
            地代と賃金の結合リスト
        RW_d: numpy.ndarray (2 * K, )
            探索方向
        alpha: float
            Armijoの定数 (0 < alpha < 1)
        beta: float
            ステップサイズを減少させる定数 (0 < beta < 1)

        Returns
        ----------
        - step_size: 適切なステップサイズ
        """
        t = 1.0
        
        while self.Z_SD(np.maximum(1e-3, RW + (t * RW_d))) > self.Z_SD(RW) + alpha * t * np.dot(self.short_dual_df(RW), RW_d):
            t *= beta

        return t

    def BFGS(self, B_k, s, y):
        """

        Parameters
        ----------
        B_k: numpy.ndarray
            現在の近似行列 (2 * K, 2 * K)
        s: numpy.ndarray (2 * K, )
            価格変数リストの変化量
        y: numpy.ndarray (2 * K, )
            勾配ベクトルの変化量

        Returns
        ----------
        B_k1: numpy.ndarray
            更新後の近似行列 (2 * K, 2 * K)
        """
        B_k1 = B_k - (np.outer(np.dot(B_k, s), np.dot(B_k, s)) / np.dot(np.dot(B_k, s), s)) + np.outer(y, y) / np.dot(y, s)

        return B_k1