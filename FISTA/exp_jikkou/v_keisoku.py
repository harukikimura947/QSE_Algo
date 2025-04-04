from exp_nume_ex import lattice
from exp_nume_ex import exp_kakunou
from exp_nume_ex import exp_qseplot
from exp_nume_ex import exp_Short_FISTA
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
from matplotlib.ticker import ScalarFormatter
from matplotlib.font_manager import FontProperties


E = 5

Col = 20

K = Col * Col

M = 0.1 * K
N = 0.1 * K

#         dF_long = pd.read_excel(fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/log/excelfile/' \
#             fr'shortlong/kaisekikai/theta=5.0/mesh/E={E}/v_proj=0.1/K={Col}^2/truevalue/jikken_F.xlsx')
#         dH_long = pd.read_excel(fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/log/excelfile/' \
#             fr'shortlong/kaisekikai/theta=5.0/mesh/E={E}/v_proj=0.1/K={Col}^2/truevalue/jikken_H.xlsx')

#         m_fixed = np.array(dF_long['m'].tolist())
#         n_fixed = np.array(dH_long['n'].tolist()).reshape(K, K)

m_per = M / K
m_fixed = np.array([m_per] * K)
n_fixed = np.full((K, K), N / (K * K))

# Scaling * tが格子の最小距離となる。

Scaling = 10

alter_T_num = 0.5

city_network = lattice.make_lattice(Col)

# 通勤費用パラメータ
t = 0.1

# 距離抵抗パラメータ
tau = 0.01

# 総土地供給量
S_total = K
S_bar = S_total / K

Coordinate_Data = \
np.array([(city_network['node_dic'][str(i)]['x_pos']*Scaling,
           city_network['node_dic'][str(i)]['y_pos']*Scaling) for i in range(K)])
distance_matrix = distance.squareform(distance.pdist(Coordinate_Data))

theta_firm = 5.0
theta_house = 5.0

alpha_1 = 0.4
alpha_2 = 0.4
beta_1 = 0.4
beta_2 = 0.4

L = 0.2
eta = 1.2
p_proj = 1e-3

RW_ini = 1.0
RW_proj = 1e-3
err_short = 1e-5

method = "FISTA"
long = "ichiyou"
dic = "Scaling=30divCol"

#パラメータ設定
prm = exp_Short_FISTA.Parameter(Col, distance_matrix, t, tau, Scaling, S_total, S_bar,
                theta_firm, theta_house, E, RW_proj, alter_T_num,
                alpha_1, alpha_2, beta_1, beta_2, M, N)

algprm = exp_Short_FISTA.Algo_Parameter(L, eta, p_proj)
short = exp_Short_FISTA.Short(prm, algprm)

alter_T = alter_T_num * Scaling * t

def v():
    
    R = np.ones(prm.K)
    W = np.ones(prm.K)
    
    v = prm.E * W * np.ones((prm.K, prm.K))\
      + (1 - prm.alpha_1 - prm.alpha_2)\
      * (1 / np.maximum(prm.alter_T, prm.T)) ** (1 / (1 - prm.alpha_1 - prm.alpha_2))\
      * ((prm.alpha_1 / R) ** (prm.alpha_1  / ((1 - prm.alpha_1 - prm.alpha_2))))[:, np.newaxis]\
      * ((prm.alpha_2 / W) ** (prm.alpha_2  / ((1 - prm.alpha_1 - prm.alpha_2))))

    return v

    # v = prm.E * W * np.ones((prm.K, prm.K))\
    #   + (1 - prm.alpha_1 - prm.alpha_2)\
    #   * (1 / np.maximum(prm.alter_T, prm.T)) ** (1 / (1 - prm.alpha_1 - prm.alpha_2))\
    #   * (((prm.alpha_1 / R) * np.ones((prm.K, prm.K))) ** (prm.alpha_1 / (1 - prm.alpha_1 - prm.alpha_2))).T\
    #   * (((prm.alpha_2 / W) * np.ones((prm.K, prm.K))) ** (prm.alpha_2 / (1 - prm.alpha_1 - prm.alpha_2)))

# print("T:", np.maximum(alter_T, prm.T))

# #step.1-1
# start_time = time.process_time()
# R_true, W_true, \
# iteration, R_err, W_err, v_err, pi_err, obj_err, grad_mean, \
# R_list, W_list, pi_list, v_list, Z_list = \
# short.short_solve(RW_ini, m_fixed, n_fixed, err_short,
#                   R_corr=None, W_corr=None, pi_corr=None, v_corr=None, obj_corr=None,
#                   short_itr=10001, rel=0)
# end_time = time.process_time()
# elapsed_time = end_time - start_time

# pi_true = short.pi(R_true, W_true, m_fixed)
# v_true = short.v(R_true, W_true)

# RW_true = np.concatenate((R_true, W_true))

# print("CPU_time:", elapsed_time, "秒")