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
        
class Short:
    '''短期均衡問題クラス

    Attributes
    ----------
    prm: Parameter
       パラメータクラス 
    '''

    def __init__(self, prm, m, n):
        '''

        Parameters
        ----------
        prm: Parameter
            モデルのパラメータクラス
        algprm: Algo-Parameter
            アルゴリズムのパラメータクラス
        '''
        self.prm = prm
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
        
        dR = self.prm.S_bar - S_H * np.sum(self.n, axis=1) - S_F * self.m
        dW = self.prm.E * np.sum(self.n, axis=0) - L_H * np.sum(self.n, axis=0) - L_F * self.m

        return np.concatenate((dR, dW))