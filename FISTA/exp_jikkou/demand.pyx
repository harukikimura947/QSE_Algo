# cython: language_level=3, boundscheck=False, wraparound=False

from exp_nume_ex import lattice
import numpy as np
import scipy.sparse as spsp
from scipy.spatial import distance
from scipy.optimize import linprog
cimport numpy as cnp
ctypedef cnp.float64_t DTYPE_T

Col = 20

K = Col * Col

E = 5

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

alpha_1 = 0.4
alpha_2 = 0.4
beta_1 = 0.4
beta_2 = 0.4

L = 0.2
eta = 1.2
p_proj = 1e-3

T = np.maximum(Scaling * t * alter_T_num , t * distance_matrix)
D = np.exp(-tau * distance_matrix)
 
cpdef cnp.ndarray[DTYPE_T, ndim=2] S_H_def(cnp.ndarray[DTYPE_T, ndim=1] R, cnp.ndarray[DTYPE_T, ndim=1] W):
    '''需要関数を設定する関数

    Parameters
    ----------
    R: numpy.ndarray(matrix) (K, )
        地代リスト
    W: numpy.ndarray(matrix) (K, )
        賃金リスト
    '''
    cdef cnp.ndarray[DTYPE_T, ndim=2] S_H

    S_H = ((1 / T) ** (1 / (1 - alpha_1 - alpha_2))) \
        * (((alpha_1 / R) * np.ones((K, K))) ** ((1 - alpha_2) / (1 - alpha_1 - alpha_2))).T \
        * (((alpha_2 / W) * np.ones((K, K))) ** (alpha_2 / (1 - alpha_1 - alpha_2)))

    return S_H

cpdef cnp.ndarray[DTYPE_T, ndim=2] L_H_def(cnp.ndarray[DTYPE_T, ndim=1] R, cnp.ndarray[DTYPE_T, ndim=1] W):
    '''需要関数を設定する関数

    Parameters
    ----------
    R: numpy.ndarray(matrix) (K, )
        地代リスト
    W: numpy.ndarray(matrix) (K, )
        賃金リスト
    '''
    cdef cnp.ndarray[DTYPE_T, ndim=2] L_H

    L_H = ((1 / T) ** (1 / (1 - alpha_1 - alpha_2))) \
        * (((alpha_1 / R) * np.ones((K, K))) ** (alpha_1 / (1 - alpha_1 - alpha_2))).T \
        * (((alpha_2 / W) * np.ones((K, K))) ** ((1 - alpha_1) / (1 - alpha_1 - alpha_2)))

    return L_H

cpdef cnp.ndarray[DTYPE_T, ndim=1] S_F_def(cnp.ndarray[DTYPE_T, ndim=1] R, cnp.ndarray[DTYPE_T, ndim=1] W):
    '''需要関数を設定する関数

    Parameters
    ----------
    R: numpy.ndarray(matrix) (K, )
        地代リスト
    W: numpy.ndarray(matrix) (K, )
        賃金リスト
    '''
    cdef cnp.ndarray[DTYPE_T, ndim=1] S_F

    S_F = (beta_1 / R) ** ((1 - beta_2) / (1 - beta_1 - beta_2)) *\
          (beta_2 / W) ** (beta_2 / (1 - beta_1 - beta_2)) #(K, )

    return S_F

cpdef cnp.ndarray[DTYPE_T, ndim=1] L_F_def(cnp.ndarray[DTYPE_T, ndim=1] R, cnp.ndarray[DTYPE_T, ndim=1] W):
    '''需要関数を設定する関数

    Parameters
    ----------
    R: numpy.ndarray(matrix) (K, )
        地代リスト
    W: numpy.ndarray(matrix) (K, )
        賃金リスト
    '''
    cdef cnp.ndarray[DTYPE_T, ndim=1] L_F

    L_F = (beta_1 / R) ** (beta_1 / (1 - beta_1 - beta_2)) *\
          (beta_2 / W) ** ((1 - beta_1) / (1 - beta_1 - beta_2)) #(K, )

    return L_F