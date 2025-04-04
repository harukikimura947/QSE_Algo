# cython: language_level=3, boundscheck=False, wraparound=False

from exp_nume_ex import exp_parameter
import numpy as np
cimport numpy as cnp

cdef class btra:
    
    cdef public:
        
        exp_parameter.Parameter prm
        exp_parameter.Algo_Parameter algprm

        def __init__(self, exp_parameter.Parameter prm, exp_parameter.Algo_Parameter algprm):
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

        cdef cnp.ndarray[double, ndim=2] S_H(self, cnp.ndarray[double, ndim=1] R, cnp.ndarray[double, ndim=1] W):
            '''需要関数を設定する関数

            Parameters
            ----------
            R: numpy.ndarray(matrix) (K, )
                地代リスト
            W: numpy.ndarray(matrix) (K, )
                賃金リスト
            '''
            cdef cnp.ndarray[double, ndim=2] S_H

            S_H = ((1 / self.prm.T) ** (1 / (1 - self.prm.alpha_1 - self.prm.alpha_2))) \
                * (((self.prm.alpha_1 / R) * np.ones((self.prm.K, self.prm.K))) ** ((1 - self.prm.alpha_2) / (1 - self.prm.alpha_1 - self.prm.alpha_2))).T \
                * (((self.prm.alpha_2 / W) * np.ones((self.prm.K, self.prm.K))) ** (self.prm.alpha_2 / (1 - self.prm.alpha_1 - self.prm.alpha_2)))

            return S_H

        cdef cnp.ndarray[double, ndim=2] L_H(self, cnp.ndarray[double, ndim=1] R, cnp.ndarray[double, ndim=1] W):
            '''需要関数を設定する関数

            Parameters
            ----------
            R: numpy.ndarray(matrix) (K, )
                地代リスト
            W: numpy.ndarray(matrix) (K, )
                賃金リスト
            '''
            cdef cnp.ndarray[double, ndim=2] L_H

            L_H = ((1 / self.prm.T) ** (1 / (1 - self.prm.alpha_1 - self.prm.alpha_2))) \
                * (((self.prm.alpha_1 / R) * np.ones((self.prm.K, self.prm.K))) ** (self.prm.alpha_1 / (1 - self.prm.alpha_1 - self.prm.alpha_2))).T \
                * (((self.prm.alpha_2 / W) * np.ones((self.prm.K, self.prm.K))) ** ((1 - self.prm.alpha_1) / (1 - self.prm.alpha_1 - self.prm.alpha_2)))

            return L_H

        cdef cnp.ndarray[double, ndim=1] S_F(self, cnp.ndarray[double, ndim=1] R, cnp.ndarray[double, ndim=1] W):
            '''需要関数を設定する関数

            Parameters
            ----------
            R: numpy.ndarray(matrix) (K, )
                地代リスト
            W: numpy.ndarray(matrix) (K, )
                賃金リスト
            '''
            cdef cnp.ndarray[double, ndim=1] S_F

            S_F = (self.prm.beta_1 / R) ** ((1 - self.prm.beta_2) / (1 - self.prm.beta_1 - self.prm.beta_2)) *\
                  (self.prm.beta_2 / W) ** (self.prm.beta_2 / (1 - self.prm.beta_1 - self.prm.beta_2)) #(K, )

            return S_F

        cdef cnp.ndarray[double, ndim=1] L_F(self, cnp.ndarray[double, ndim=1] R, cnp.ndarray[double, ndim=1] W):
            '''需要関数を設定する関数

            Parameters
            ----------
            R: numpy.ndarray(matrix) (K, )
                地代リスト
            W: numpy.ndarray(matrix) (K, )
                賃金リスト
            '''
            cdef cnp.ndarray[double, ndim=1] L_F

            L_F = (self.prm.beta_1 / R) ** (self.prm.beta_1 / (1 - self.prm.beta_1 - self.prm.beta_2)) *\
                  (self.prm.beta_2 / W) ** ((1 - self.prm.beta_1) / (1 -self.prm.beta_1 - self.prm.beta_2)) #(K, )

            return L_F

        cdef cnp.ndarray[double, ndim=2] v(self, cnp.ndarray[double, ndim=1] R, cnp.ndarray[double, ndim=1] W):
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
            cdef cnp.ndarray[double, ndim=2] v

            v = self.prm.E * W * np.ones((self.prm.K, self.prm.K))\
              + (1 - self.prm.alpha_1 - self.prm.alpha_2)\
              * (1 / self.prm.T) ** (1 / (1 - self.prm.alpha_1 - self.prm.alpha_2))\
              * (((self.prm.alpha_1 / R) * np.ones((self.prm.K, self.prm.K))) ** (self.prm.alpha_1 / (1 - self.prm.alpha_1 - self.prm.alpha_2))).T\
              * (((self.prm.alpha_2 / W) * np.ones((self.prm.K, self.prm.K))) ** (self.prm.alpha_2 / (1 - self.prm.alpha_1 - self.prm.alpha_2)))

            return v

        cdef cnp.ndarray[double, ndim=1] pi(self, 
                                            cnp.ndarray[double, ndim=1] R, cnp.ndarray[double, ndim=1] W, 
                                            cnp.ndarray[double, ndim=1] m):
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
            cdef cnp.ndarray[double, ndim=1] pi

            pi = np.dot(self.prm.D, m)\
               + (1 - self.prm.beta_1 - self.prm.beta_2)\
               * ((self.prm.beta_1 / R) ** (self.prm.beta_1 / (1 - self.prm.beta_1 - self.prm.beta_2)))\
               * ((self.prm.beta_2 / W) ** (self.prm.beta_2 / (1 - self.prm.beta_1 - self.prm.beta_2)))

            return pi

        cdef cnp.ndarray[double, ndim=1] short_dual_df(self, 
                                                       cnp.ndarray[double, ndim=1] RW, cnp.ndarray[double, ndim=1] m, cnp.ndarray[double, ndim=2] n):
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

            cdef cnp.ndarray[double, ndim=1] R
            cdef cnp.ndarray[double, ndim=1] W
            cdef cnp.ndarray[double, ndim=2] S_H
            cdef cnp.ndarray[double, ndim=2] L_H
            cdef cnp.ndarray[double, ndim=1] S_F
            cdef cnp.ndarray[double, ndim=1] L_F
            cdef cnp.ndarray[double, ndim=1] dR
            cdef cnp.ndarray[double, ndim=1] dW
            cdef cnp.ndarray[double, ndim=1] dRW

            R = RW[:self.prm.K]
            W = RW[self.prm.K:]

            S_H = self.S_H(R, W)
            L_H = self.L_H(R, W)
            S_F = self.S_F(R, W)
            L_F = self.L_F(R, W)

            dR = self.prm.S_bar - np.sum(S_H * n, axis=1) - S_F * m  # forall i
            dW = self.prm.E * np.sum(n, axis=0) - np.sum(L_H * n, axis=0) - L_F * m # forall j

            dRW = np.concatenate((dR, dW))

            return dRW

        cdef double Z_SD(self, 
                         cnp.ndarray[double, ndim=1] RW, cnp.ndarray[double, ndim=1] m, cnp.ndarray[double, ndim=2] n):
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

            cdef cnp.ndarray[double, ndim=1] R
            cdef cnp.ndarray[double, ndim=1] W
            cdef double F_value

            R = RW[:self.prm.K]
            W = RW[self.prm.K:]

            F_value = np.sum(self.v(R, W) * n)\
                    + np.sum((self.pi(R, W, m) - np.dot(self.prm.D, m)) * m)\
                    + self.prm.S_bar * np.sum(R)

            return F_value

        cdef double Q(self, 
                      cnp.ndarray[double, ndim=1] p, cnp.ndarray[double, ndim=1] p_bar,
                      double L_bar, cnp.ndarray[double, ndim=1] m_fixed, cnp.ndarray[double, ndim=2] n_fixed):

            cdef double Q

            Q = self.Z_SD(p_bar, m_fixed, n_fixed) + np.dot(p - p_bar, self.short_dual_df(p_bar, m_fixed, n_fixed)) \
              + (L_bar / 2) * np.linalg.norm(p - p_bar) ** 2

            # Q_2 = self.Z_SD(p_bar) - (1 / (2 * L_bar)) * (np.linalg.norm(self.short_dual_df(p_bar)) ** 2)

            return Q

        cpdef double backtracking(self, double L,
                                  cnp.ndarray[double, ndim=1] p_bar, 
                                  cnp.ndarray[double, ndim=1] m_fixed, 
                                  cnp.ndarray[double, ndim=2] n_fixed):

            cdef double i = 0
            cdef double L_bar
            cdef cnp.ndarray[double, ndim=1] p

            while True:
                L_bar = (self.algprm.eta ** i) * L

                p = np.maximum(self.algprm.p_proj, p_bar - (self.short_dual_df(p_bar, m_fixed, n_fixed) / L_bar))

                if self.Z_SD(p, m_fixed, n_fixed) <= self.Q(p, p_bar, L_bar, m_fixed, n_fixed):
                    break

                i = i + 1

            return i