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

    def __init__(self,
                 Col, distance_matrix, t, tau, Scaling, S_total, S_bar,
                 theta_firm, theta_house, E, RW_proj, alter_T_num,
                 alpha_1, alpha_2, beta_1, beta_2, M, N):
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
        self.RW_proj = RW_proj
        self.alter_T = Scaling * t * alter_T_num
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
        
        S_H = ((1 / np.maximum(self.prm.alter_T, self.prm.T)) ** (1 / (1 - self.prm.alpha_1 - self.prm.alpha_2))) \
            * (((self.prm.alpha_1 / R) * np.ones((self.prm.K, self.prm.K))) ** ((1 - self.prm.alpha_2) / (1 - self.prm.alpha_1 - self.prm.alpha_2))).T \
            * (((self.prm.alpha_2 / W) * np.ones((self.prm.K, self.prm.K))) ** (self.prm.alpha_2 / (1 - self.prm.alpha_1 - self.prm.alpha_2)))
        
        L_H = ((1 / np.maximum(self.prm.alter_T, self.prm.T)) ** (1 / (1 - self.prm.alpha_1 - self.prm.alpha_2))) \
            * (((self.prm.alpha_1 / R) * np.ones((self.prm.K, self.prm.K))) ** (self.prm.alpha_1 / (1 - self.prm.alpha_1 - self.prm.alpha_2))).T \
            * (((self.prm.alpha_2 / W) * np.ones((self.prm.K, self.prm.K))) ** ((1 - self.prm.alpha_1) / (1 - self.prm.alpha_1 - self.prm.alpha_2)))
        
        S_F = (self.prm.beta_1 / R) ** ((1 - self.prm.beta_2) / (1 - self.prm.beta_1 - self.prm.beta_2)) *\
              (self.prm.beta_2 / W) ** (self.prm.beta_2 / (1 - self.prm.beta_1 - self.prm.beta_2)) #(K, )
        
        L_F = (self.prm.beta_1 / R) ** (self.prm.beta_1 / (1 - self.prm.beta_1 - self.prm.beta_2)) *\
              (self.prm.beta_2 / W) ** ((1 - self.prm.beta_1) / (1 - self.prm.beta_1 - self.prm.beta_2)) #(K, )
        
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

        pi = np.dot(self.prm.D, self.m)\
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
        
        v = self.prm.E * W * np.ones((self.prm.K, self.prm.K))\
          + (1 - self.prm.alpha_1 - self.prm.alpha_2)\
          * (1 / np.maximum(self.prm.alter_T, self.prm.T)) ** (1 / (1 - self.prm.alpha_1 - self.prm.alpha_2))\
          * (((self.prm.alpha_1 / R) * np.ones((self.prm.K, self.prm.K))) ** (self.prm.alpha_1 / (1 - self.prm.alpha_1 - self.prm.alpha_2))).T\
          * (((self.prm.alpha_2 / W) * np.ones((self.prm.K, self.prm.K))) ** (self.prm.alpha_2 / (1 - self.prm.alpha_1 - self.prm.alpha_2)))

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
        
        F_value = np.sum(self.v(R, W) * self.n)\
                + np.sum((self.pi(R, W) - np.dot(self.prm.D, self.m)) * self.m)\
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
        
        dR = self.prm.S_bar - np.sum(S_H * self.n, axis=1) - S_F * self.m  # forall i
        dW = self.prm.E * np.sum(self.n, axis=0) - np.sum(L_H * self.n, axis=0) - L_F * self.m # forall j

        return np.concatenate((dR, dW))
    
    def callback(self, RW):
        
        R = RW[:self.prm.K]
        W = RW[self.prm.K:]
        
        # print("grad:", self.short_dual_df(RW))
        # print("v:", self.v(R, W))
        # print("R:", R)
        # print("W:", W)
        # print("obj:", self.Z_SD(RW))

    def hanyou(self, sci_alg, ftol):

        R_ini = np.ones(self.prm.K)
        W_ini = np.ones(self.prm.K)

        ini = np.concatenate((R_ini, W_ini))
        
        if np.min(self.v(R_ini, W_ini)) < 0:
            raise ValueError("v must be non-negative")

        bounds = [(self.prm.RW_proj, None)] * len(ini)
        
        options = {"ftol":ftol, "maxfun":10000, "maxls":40000, "maxiter":40000, "iprint":100}

        result = minimize(self.Z_SD, ini, bounds = bounds, method = f'{sci_alg}', jac = self.short_dual_df, constraints = None, callback = self.callback, options = options)

        RW = result.x
        
        R_corr = RW[:self.prm.K]
        W_corr = RW[self.prm.K:]
        
        equ_R, equ_W = self.equilibrium(R_corr, W_corr)
        
        # print("success:", result.success)
        # print("iteration:", result.nit)
        # print("status:", result.status)
        # print("message:", result.message)
        # print("jac:", np.mean(result.jac))
        # print("equ_R:", equ_R)
        # print("equ_W:", equ_W)

        return R_corr, W_corr
    
    def equilibrium(self, R, W):
        
        RW = np.concatenate((R, W))
        
        dRW = self.short_dual_df(RW)
        
        dR = dRW[:self.prm.K]
        dW = dRW[self.prm.K:]
        
        equ_R = np.max(R * dR)
        equ_W = np.max(W * dW)
        
        # G = np.sum(R * dR) + np.sum(W * dW)
        
        return equ_R, equ_W