"""Module summary.

This is Short-Algo.py.
"""

import time
import numpy as np
import pandas as pd
import numba
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.tri import Triangulation
# import cProfile
import scipy.optimize as optimize
import scipy.sparse as spsp
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
        
        #最小距離の半分の長さ
        
        S_H = np.maximum(0, self.prm.alpha_1 / R)
        L_H = np.maximum(0, self.prm.alpha_2 / W)
        S_F = np.maximum(0, self.prm.beta_1 / R)
        L_F = np.maximum(0, self.prm.beta_2 / W)
        
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

        pi = np.maximum(0.1, np.dot(self.prm.D, m)\
        - self.prm.beta_1 * np.log(R) - self.prm.beta_2 * np.log(W)\
        + self.prm.beta_1 * np.log(self.prm.beta_1) + self.prm.beta_2 * np.log(self.prm.beta_2)\
        - self.prm.beta_1 - self.prm.beta_2)
        
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
        
        pi_noex = np.maximum(0.1, 
        - self.prm.beta_1 * np.log(R) - self.prm.beta_2 * np.log(W)\
        + self.prm.beta_1 * np.log(self.prm.beta_1) + self.prm.beta_2 * np.log(self.prm.beta_2)\
        - self.prm.beta_1 - self.prm.beta_2)
        
        # print("pi_noex:", pi_noex)
        
        R_sum = np.sum(R)
        F_value = np.sum(self.v(R, W) * n)\
        + np.sum(pi_noex * m)\
        + self.prm.S_bar * R_sum
        
        return F_value
    
    def short_dual_df(self, RW: np.ndarray, m: np.ndarray, n: np.ndarray):
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
        
        dR = self.prm.S_bar - S_H * np.sum(n, axis=1) - S_F * m
        dW = self.prm.E * np.sum(n, axis=0) - L_H * np.sum(n, axis=0) - L_F * m

        return np.concatenate((dR, dW))

    def short_solve(self, m, n, short_itr, err_short):
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

        RW_before = np.ones(2 * self.prm.K)
        
        p_bar_before = np.ones(2 * self.prm.K)
        
        L_before = self.algprm.L
        
        t_before = 1
        
        max_value = 0
        
        iteration = []
        R_err = []
        W_err = []
        v_err = []
        pi_err = []
        obj_err = []

        for k in range(1, short_itr):
            # print("="*40)
            # print("kshort:", max_value)
            
            #Step.2 勾配を計算する
            i = self.backtracking(p_bar_before, m, n, L_before)
            
            L = (self.algprm.eta ** i) * L_before
            
            dZ_SD = self.short_dual_df(p_bar_before, m, n)
            
            RW = np.maximum(1e-3, p_bar_before - (dZ_SD / L))
            
            t = (1 + np.sqrt(1 + 4 * (t_before ** 2))) / 2
            
            p_bar = np.maximum(1e-3, RW + ((t_before - 1) / t) * (RW - RW_before))
            
            if np.dot(self.short_dual_df(RW, m, n), (RW - RW_before)) > 0:
                t = 1
            
            if k == 1 or k % 30 == 0:
                
                R = RW[:self.prm.K]
                W = RW[self.prm.K:]
                
                pi = self.pi(R, W, m)

                # 3次元行列の生成
                matrix = np.reshape(pi, (self.prm.Col, self.prm.Col))

                # 行列の座標を作成
                x, y = np.meshgrid(range(self.prm.Col), range(self.prm.Col))

                # 3Dグラフの作成
                fig = plt.figure()
                ax = fig.add_subplot(111, projection='3d')

                # plot_surfaceを使用して3Dグラフを描画
                # ax.plot_surface(x, y, matrix, cmap='viridis')
                ax.plot_wireframe(x, y, matrix, cmap='viridis')

                # グラフの軸ラベル設定
                # ax.set_xlabel('Col')
                # ax.set_ylabel('Row')
                ax.set_zlabel('pi')
                
                ax.text2D(0.5, -0.1, f'iteration={k}', transform=ax.transAxes, ha='center')

                plt.gca().invert_xaxis()

                # グラフ表示
                plt.show()
            
            max_value = k
            
            #Step.5 収束判定
            if np.max(abs((RW - RW_before) / RW_before)) < err_short:
                break
            
            RW_before = RW
            p_bar_before = p_bar
            L_before = L
            t_before = t
            
        print("short_max_value:", max_value)
        
        R_hist = RW[:self.prm.K]
        W_hist = RW[self.prm.K:]
        
        S_H, L_H, S_F, L_F = self.demand(R_hist, W_hist)
        
        print("Lyapunov:", np.sum(R_hist * (self.prm.S_bar - S_H * np.sum(n, axis=1) - S_F * m))\
              + np.sum(W_hist * (self.prm.E * np.sum(n, axis=0) - L_H * np.sum(n, axis=0) - L_F * m)))
        
        return None
    
    def Q(self, p, p_bar, m, n, L_bar):
        
        Q = self.Z_SD(p_bar, m, n) + np.dot(p - p_bar, self.short_dual_df(p_bar, m, n)) + (L_bar / 2) * np.linalg.norm(p - p_bar)
        
        # Q = self.Z_SD(y, m, n) - (1 / (2 * L_bar)) * np.linalg.norm(self.short_dual_df(y, m, n))
        
        return Q
    
    def backtracking(self, p_bar, m, n, L):
        
        i = 0
        
        while True:
            L_bar = (self.algprm.eta ** i) * L
            p = np.maximum(self.algprm.p_proj, p_bar - (self.short_dual_df(p_bar, m, n) / L_bar))
            if self.Z_SD(p, m, n) <= self.Q(p, p_bar, m, n, L_bar):
                break
            
            i = i + 1
            
        return i