from exp_nume_ex import exp_Short_FISTA
from exp_nume_ex import lattice
from exp_nume_ex import exp_kakunou
from exp_nume_ex import exp_qseplot
from exp_nume_ex import exp_Long_FISTA
from exp_nume_ex import exp_parameter
from dir_exp_sci import exp_lattice_sci
from dir_exp_sci import exp_kakunou_sci
from dir_exp_sci import exp_Short_sci
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

E = 5
Col = 20
switch = 1

K = Col * Col

M = 0.1 * K
N = 0.1 * K

# Scaling * tが格子の最小距離となる。
Scaling = 30
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

theta_firm = 1.0
theta_house = 1.0

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
err_long = 1e-3

m_per = M / K
m0 = np.array([m_per] * K)
n0 = np.full((K, K), N / (K * K))

method = "FISTA"
dic = "Scaling=300divCol"

#パラメータ設定
prm = exp_parameter.Parameter(
            Col, distance_matrix, t, tau, Scaling, S_total, S_bar,
            theta_firm, theta_house, E, RW_proj, alter_T_num,
            alpha_1, alpha_2, beta_1, beta_2, M, N)

algprm = exp_parameter.Algo_Parameter(L, eta, p_proj)
short = exp_Short_FISTA.Short(prm, algprm)
long = exp_Long_FISTA.Long(prm, algprm, short)

start_time = time.process_time()
m, n, RW, long_iteration, obj_list = long.solve(m0, n0, RW_ini, err_short, err_long, long_itr=100000)
end_time = time.process_time()
elapsed_time = end_time - start_time

print("CPU_time:", elapsed_time, "秒")

R = RW[:prm.K]
W = RW[prm.K:]

pi = short.pi(R, W, m)
v = short.v(R, W)

# exp_kakunou.long_kakunou_true(m, n, R, W, pi, v, prm, algprm, method, err_short, err_long, dic)
# exp_kakunou.long_kakunou_nume_num(long_iteration, prm, algprm, method, err_short, err_long, dic)