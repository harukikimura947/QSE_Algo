# cython: language_level=3, boundscheck=False, wraparound=False

from exp_nume_ex import lattice
import numpy as np
import scipy.sparse as spsp
from scipy.spatial import distance
from scipy.optimize import linprog
cimport numpy as cnp
ctypedef cnp.float64_t DTYPE_T #numpyの詳細な型指定

# cdef int Col
# cdef int K
# cdef double Scaling
# cdef double alter_T_num
# cdef double t
# cdef double tau
# cdef double S_total
# cdef double S_bar
# cdef double alpha_1
# cdef double alpha_2
# cdef double beta_1
# cdef double beta_2
# cdef double L
# cdef double eta
# cdef double p_proj
# cdef cnp.ndarray[double, ndim=2] T
# cdef cnp.ndarray[double, ndim=2] D

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
 
cdef cnp.ndarray[DTYPE_T, ndim=2] S_H_def(cnp.ndarray[DTYPE_T, ndim=1] R, cnp.ndarray[DTYPE_T, ndim=1] W):
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

cdef cnp.ndarray[DTYPE_T, ndim=2] L_H_def(cnp.ndarray[DTYPE_T, ndim=1] R, cnp.ndarray[DTYPE_T, ndim=1] W):
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

cdef cnp.ndarray[DTYPE_T, ndim=1] S_F_def(cnp.ndarray[DTYPE_T, ndim=1] R, cnp.ndarray[DTYPE_T, ndim=1] W):
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

cdef cnp.ndarray[DTYPE_T, ndim=1] L_F_def(cnp.ndarray[DTYPE_T, ndim=1] R, cnp.ndarray[DTYPE_T, ndim=1] W):
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

cdef cnp.ndarray[DTYPE_T, ndim=2] v(cnp.ndarray[DTYPE_T, ndim=1] R, cnp.ndarray[DTYPE_T, ndim=1] W):
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
    cdef cnp.ndarray[DTYPE_T, ndim=2] v

    v = E * W * np.ones((K, K))\
      + (1 - alpha_1 - alpha_2)\
      * (1 / T) ** (1 / (1 - alpha_1 - alpha_2))\
      * (((alpha_1 / R) * np.ones((K, K))) ** (alpha_1 / (1 - alpha_1 - alpha_2))).T\
      * (((alpha_2 / W) * np.ones((K, K))) ** (alpha_2 / (1 - alpha_1 - alpha_2)))

    return v

cdef cnp.ndarray[DTYPE_T, ndim=1] pi(cnp.ndarray[DTYPE_T, ndim=1] R, cnp.ndarray[DTYPE_T, ndim=1] W, 
                                    cnp.ndarray[DTYPE_T, ndim=1] m):
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
    cdef cnp.ndarray[DTYPE_T, ndim=1] pi

    pi = np.dot(D, m)\
       + (1 - beta_1 - beta_2)\
       * ((beta_1 / R) ** (beta_1 / (1 - beta_1 - beta_2)))\
       * ((beta_2 / W) ** (beta_2 / (1 - beta_1 - beta_2)))

    return pi

cdef cnp.ndarray[DTYPE_T, ndim=1] short_dual_df(cnp.ndarray[DTYPE_T, ndim=1] RW, cnp.ndarray[DTYPE_T, ndim=1] m, cnp.ndarray[DTYPE_T, ndim=2] n):
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

    cdef cnp.ndarray[DTYPE_T, ndim=1] R
    cdef cnp.ndarray[DTYPE_T, ndim=1] W
    cdef cnp.ndarray[DTYPE_T, ndim=2] S_H
    cdef cnp.ndarray[DTYPE_T, ndim=2] L_H
    cdef cnp.ndarray[DTYPE_T, ndim=1] S_F
    cdef cnp.ndarray[DTYPE_T, ndim=1] L_F
    cdef cnp.ndarray[DTYPE_T, ndim=1] dR
    cdef cnp.ndarray[DTYPE_T, ndim=1] dW
    cdef cnp.ndarray[DTYPE_T, ndim=1] dRW

    R = RW[:K]
    W = RW[K:]

    S_H = S_H_def(R, W)
    L_H = L_H_def(R, W)
    S_F = S_F_def(R, W)
    L_F = L_F_def(R, W)

    dR = S_bar - np.sum(S_H * n, axis=1) - S_F * m  # forall i
    dW = E * np.sum(n, axis=0) - np.sum(L_H * n, axis=0) - L_F * m # forall j

    dRW = np.concatenate((dR, dW))

    return dRW

cdef double Z_SD(cnp.ndarray[DTYPE_T, ndim=1] RW, cnp.ndarray[DTYPE_T, ndim=1] m, cnp.ndarray[DTYPE_T, ndim=2] n):
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

    cdef cnp.ndarray[DTYPE_T, ndim=1] R
    cdef cnp.ndarray[DTYPE_T, ndim=1] W
    cdef DTYPE_T F_value

    R = RW[:K]
    W = RW[K:]

    F_value = np.sum(v(R, W) * n)\
            + np.sum((pi(R, W, m) - np.dot(D, m)) * m)\
            + S_bar * np.sum(R)

    return F_value

cdef double Q(cnp.ndarray[DTYPE_T, ndim=1] p, cnp.ndarray[DTYPE_T, ndim=1] p_bar,
              DTYPE_T L_bar, cnp.ndarray[DTYPE_T, ndim=1] m_fixed, cnp.ndarray[DTYPE_T, ndim=2] n_fixed):

    cdef DTYPE_T Q

    Q = Z_SD(p_bar, m_fixed, n_fixed) + np.dot(p - p_bar, short_dual_df(p_bar, m_fixed, n_fixed)) \
      + (L_bar / 2) * np.linalg.norm(p - p_bar) ** 2

    # Q_2 = self.Z_SD(p_bar) - (1 / (2 * L_bar)) * (np.linalg.norm(self.short_dual_df(p_bar)) ** 2)

    return Q

cpdef double backtracking(DTYPE_T L,
                          cnp.ndarray[DTYPE_T, ndim=1] p_bar, 
                          cnp.ndarray[DTYPE_T, ndim=1] m_fixed, 
                          cnp.ndarray[DTYPE_T, ndim=2] n_fixed):

    cdef DTYPE_T ita
    cdef DTYPE_T L_bar
    cdef cnp.ndarray[DTYPE_T, ndim=1] p
    cdef cnp.ndarray[DTYPE_T, ndim=1] p_bar_R
    cdef cnp.ndarray[DTYPE_T, ndim=1] p_bar_W
    cdef cnp.ndarray[DTYPE_T, ndim=2] v_Q
    cdef cnp.ndarray[DTYPE_T, ndim=1] pi_Q
    cdef cnp.ndarray[DTYPE_T, ndim=1] p_R
    cdef cnp.ndarray[DTYPE_T, ndim=1] p_W
    cdef cnp.ndarray[DTYPE_T, ndim=2] v
    cdef cnp.ndarray[DTYPE_T, ndim=1] pi
    cdef cnp.ndarray[DTYPE_T, ndim=2] S_H
    cdef cnp.ndarray[DTYPE_T, ndim=2] L_H
    cdef cnp.ndarray[DTYPE_T, ndim=1] S_F
    cdef cnp.ndarray[DTYPE_T, ndim=1] L_F
    cdef cnp.ndarray[DTYPE_T, ndim=1] dR
    cdef cnp.ndarray[DTYPE_T, ndim=1] dW
    cdef cnp.ndarray[DTYPE_T, ndim=1] dRW
    cdef DTYPE_T Z
    cdef DTYPE_T Z_Q
    cdef DTYPE_T Q
    cdef cnp.ndarray[cnp.int_t, ndim=1] index
    
    index = np.arange(0, 101, 1)
    
    p_bar_R = p_bar[:K]
    p_bar_W = p_bar[K:]
    
    v_Q = E * p_bar_W * np.ones((K, K))\
            + (1 - alpha_1 - alpha_2)\
            * (1 / T) ** (1 / (1 - alpha_1 - alpha_2))\
            * (((alpha_1 / p_bar_R) * np.ones((K, K))) ** (alpha_1 / (1 - alpha_1 - alpha_2))).T\
            * (((alpha_2 / p_bar_W) * np.ones((K, K))) ** (alpha_2 / (1 - alpha_1 - alpha_2)))
        
    pi_Q = np.dot(D, m_fixed)\
         + (1 - beta_1 - beta_2)\
         * ((beta_1 / p_bar_R) ** (beta_1 / (1 - beta_1 - beta_2)))\
         * ((beta_2 / p_bar_W) ** (beta_2 / (1 - beta_1 - beta_2)))

    Z_Q = np.sum(v_Q * n_fixed) + np.sum((pi_Q - np.dot(D, m_fixed)) * m_fixed) + S_bar * np.sum(p_bar_R)

    S_H = ((1 / T) ** (1 / (1 - alpha_1 - alpha_2))) \
        * (((alpha_1 / p_bar_R) * np.ones((K, K))) ** ((1 - alpha_2) / (1 - alpha_1 - alpha_2))).T \
        * (((alpha_2 / p_bar_W) * np.ones((K, K))) ** (alpha_2 / (1 - alpha_1 - alpha_2)))
    
    L_H = ((1 / T) ** (1 / (1 - alpha_1 - alpha_2))) \
        * (((alpha_1 / p_bar_R) * np.ones((K, K))) ** (alpha_1 / (1 - alpha_1 - alpha_2))).T \
        * (((alpha_2 / p_bar_W) * np.ones((K, K))) ** ((1 - alpha_1) / (1 - alpha_1 - alpha_2)))
    
    S_F = (beta_1 / p_bar_R) ** ((1 - beta_2) / (1 - beta_1 - beta_2)) *\
          (beta_2 / p_bar_W) ** (beta_2 / (1 - beta_1 - beta_2)) #(K, )
    
    L_F = (beta_1 / p_bar_R) ** (beta_1 / (1 - beta_1 - beta_2)) *\
          (beta_2 / p_bar_W) ** ((1 - beta_1) / (1 - beta_1 - beta_2)) #(K, )
    
    dR = S_bar - np.sum(S_H * n_fixed, axis=1) - S_F * m_fixed  # forall i
    dW = E * np.sum(n_fixed, axis=0) - np.sum(L_H * n_fixed, axis=0) - L_F * m_fixed # forall j
    
    dRW = np.concatenate((dR, dW))
    
    

    # while True:
    for i in range(index.shape[0]):
        
        L_bar = (eta ** index[i]) * L

        p = np.maximum(p_proj, p_bar - (short_dual_df(p_bar, m_fixed, n_fixed) / L_bar))
        
        p_R = p_bar[:K]
        p_W = p_bar[K:]
        
        v = E * p_W * np.ones((K, K))\
            + (1 - alpha_1 - alpha_2)\
            * (1 / T) ** (1 / (1 - alpha_1 - alpha_2))\
            * (((alpha_1 / p_R) * np.ones((K, K))) ** (alpha_1 / (1 - alpha_1 - alpha_2))).T\
            * (((alpha_2 / p_W) * np.ones((K, K))) ** (alpha_2 / (1 - alpha_1 - alpha_2)))
        
        pi = np.dot(D, m_fixed)\
             + (1 - beta_1 - beta_2)\
             * ((beta_1 / p_R) ** (beta_1 / (1 - beta_1 - beta_2)))\
             * ((beta_2 / p_W) ** (beta_2 / (1 - beta_1 - beta_2)))

        Z = np.sum(v * n_fixed) + np.sum((pi - np.dot(D, m_fixed)) * m_fixed) + S_bar * np.sum(p_R)
        
        Q = Z_Q + np.dot(p - p_bar, dRW) + (L_bar / 2) * np.linalg.norm(p - p_bar) ** 2

        if Z <= Q:
            ita = i
            break

    return ita