"""Module summary.

This is Short-Algo.py.
"""
# from exp_nume_ex import exp_Short_FISTA
# from exp_nume_ex import shusei_exp_Short_FISTA
# from exp_nume_ex import exp_Long_FISTA
# from exp_nume_ex import lattice
# from exp_nume_ex import exp_kakunou
# from exp_nume_ex import exp_qseplot
# from exp_nume_ex import exp_parameter
# from exp_nume_ex import vec_exp_parameter
# from exp_nume_ex import vec_exp_Short_FISTA
# from exp_nume_ex import vec_exp_Long_FISTA
from Sakai_mat import Sakai_mat_exp_parameter
from Sakai_mat import Sakai_mat_exp_Short_FISTA
from Sakai_mat import Sakai_mat_exp_Long_FISTA
from Sakai_mat import lattice
import time
import numpy as np
import pandas as pd
import scipy.optimize as optimize
import scipy.sparse as spsp
import matplotlib.pyplot as plt
import csv
import os
import cProfile
import pstats

from scipy.spatial import distance
from scipy.optimize import linprog
from scipy.optimize import minimize
from collections import defaultdict
from mpl_toolkits.mplot3d import Axes3D

def main():

    E = 5
    Col = 10

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
    obj_corr = 1.0

    m_per = M / K
    m0 = np.array([m_per] * K)
    n0 = np.full((K, K), N / (K * K)).flatten()

    method = "FISTA"
    dic = "Scaling=30"

    #パラメータ設定
#     prm = exp_parameter.Parameter(
#                 Col, distance_matrix, t, tau, Scaling, S_total, S_bar,
#                 theta_firm, theta_house, E, RW_proj, alter_T_num,
#                 alpha_1, alpha_2, beta_1, beta_2, M, N)

#     algprm = exp_parameter.Algo_Parameter(L, eta, p_proj)
#     short = shusei_exp_Short_FISTA.Short(prm, algprm)
#     long = exp_Long_FISTA.Long(prm, algprm, short)

    prm = Sakai_mat_exp_parameter.Parameter(
                Col, distance_matrix, t, tau, Scaling, S_total, S_bar,
                theta_firm, theta_house, E, RW_proj, alter_T_num,
                alpha_1, alpha_2, beta_1, beta_2, M, N)

    algprm = Sakai_mat_exp_parameter.Algo_Parameter(L, eta, p_proj)
    short = Sakai_mat_exp_Short_FISTA.Short(prm, algprm)
    long = Sakai_mat_exp_Long_FISTA.Long(prm, algprm, short)

    m, n, RW, long_iteration, obj_list, obj_rel_list = long.solve(m0, n0, RW_ini, err_short, err_long, obj_corr, long_itr = 100000)

if __name__ == "__main__":
    cProfile.run("main()", filename="main_mat.prof")
    # sts = pstats.Stats("main_vec.prof")
    
#     sv = snakeviz.stats.SnakevizStats(sts)
#     html = sv.html()
    
#     with open(output_file, "w") as f:
#         f.write(html)
    
#     save_snakeviz_html(stats, "profile_results.html")
    
    # sts.strip_dirs().sort_stats("cumulative").print_stats()