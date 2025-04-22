import time
import numpy as np
import pandas as pd
import numba
import scipy.optimize as optimize
import scipy.sparse as spsp
import matplotlib.pyplot as plt
import csv
import os

from scipy.sparse import csr_matrix
from scipy.spatial import distance
from scipy.optimize import linprog
from scipy.optimize import minimize
from collections import defaultdict
from mpl_toolkits.mplot3d import Axes3D

class Long:
    '''長期均衡問題クラス

    Attributes
    ----------
    prm: Parameter
       パラメータクラス 
    '''

    def __init__(self, prm, algprm, short):
        '''

        Parameters
        ----------
        prm: Parameter
            パラメータクラス
        long: Long
            長期問題クラス
        '''
        self.prm = prm
        self.algprm = algprm
        self.short = short
        
    def bond(self, m, n):
        '''企業・家計分布を一つの結合リストにする関数

        Parameters
        ----------
        m: numpy.ndarray (K, )
            企業分布
        n: numpy.ndarray (K, K)
            家計分布

        Returns
        ----------
        mn: numpy.ndarray (K + K * K, )
            企業・家計分布の結合リスト
        '''
        
        # 行列を1次元のリストに変換
        n_list = n.ravel()
        
        mn = np.concatenate((m, n_list))

        return mn
    
    def breakingdown(self, mn):
        '''結合リストを企業・家計分布に分解する関数

        Parameters
        ----------
        mn: numpy.ndarray (K + K * K, )
            企業・家計分布の結合リスト
        
        Returns
        ----------
        m: numpy.ndarray (K, )
            企業・家計分布の結合リスト
        n: numpy.ndarray (K, K)
            家計分布
        '''
        
        m = np.array(mn[:self.prm.K])
        
        n_list = mn[self.prm.K:]
        n = np.array([n_list[i:i+self.prm.K] for i in range(0, len(n_list), self.prm.K)])

        return m, n
    
    def Z_LP(self, mn, RW):
        '''目的関数を計算する関数

        Parameters
        ----------
        mn: numpy.ndarray (K + K * K, )
            企業・家計分布の結合リスト

        Returns
        ----------
        F: float
            目的関数値
        '''

        m, n = self.breakingdown(mn)
        
        Short = self.short.Z_SD(RW, m, n)
        
        F_value = - Short \
        - (1 / 2) * m@self.prm.D@m\
        + (1 / self.prm.theta_firm) *  (m @ np.log(m / self.prm.M)) \
        + (1 / self.prm.theta_house) * np.sum(n * np.log(n / self.prm.N))

        # print("len(F_value)", len(F_value))

        return F_value
    
    def long_df(self, mn, RW):
        '''目的関数の勾配を計算する関数

        Parameters
        ----------
        mn: numpy.ndarray (K + K * K, )
            企業・家計分布の結合リスト
        RW: numpy.ndarray (2 * K, )
            価格変数の結合リスト

        Returns
        ----------
        dF: numpy.ndarray (K + K * K, )
            目的関数の勾配
        '''

        m, n = self.breakingdown(mn)

        dF_m = - self.short.pi(R, W, m) + (1/self.prm.theta_firm) * (np.log(m / self.prm.M) + 1) # numpy.ndarray(K, )
        dF_n = - self.short.v(R, W) + (1/self.prm.theta_house) * (np.log(n / self.prm.N) + 1) # numpy.ndarray(K, K)
        
        dF = self.bond(dF_m, dF_n)

        return dF

    def solve(self, m0, n0, RW_ini, err_short, err_long, obj_corr, long_itr):
        '''長期均衡を解く関数

        Parameters
        ----------
        m0: numpy.ndarray (K, )
            企業の初期分布
        n0: numpy.ndarray (K, K)
            家計の初期分布
        max_itr: int
            最大反復回数
        err_long: float
            収束判定の閾値
            
        Returns
        -------
        mn_k: numpy.ndarray (K + K * K, )
            企業・家計分布の結合リスト
        '''

        max_value = 0

        obj_before = 0
        long_iteration = []
        obj_list = []
        obj_rel_list = []
        
        beta = 10
        
        #Step.0 初期実行解を決める
        m_before = m0
        n_before = n0

        long_iteration_append = long_iteration.append
        obj_list_append = obj_list.append
        obj_rel_list_append = obj_rel_list.append

        for k in range(1, long_itr):
            # print("="*40)
            # print("klong:", k)
            
            #Step.1 短期均衡を解く
            R, W, iteration, short_obj_rel = \
            self.short.short_solve(RW_ini, m_before, n_before, err_short, short_itr=100000, rel=0)
            
            RW = np.concatenate((R, W))
            
            #Step.2 探索方向を決める
            m_d, n_d = self.logit(R, W, m_before)
            
            #Step.3 解を更新する
            long_iteration_append(k)

            #Step.4 収束判定
#             max_value = k
#             if obj_rel < 1e-4:
#                 break
            
#             if k > 1:
#                 if obj > obj_before:
#                     beta += 2
            
            # alpha = 1 / beta
            # print("alpha:", alpha)

            alpha = 0.05
            m = (1 - alpha) * m_before + alpha * m_d
            n = (1 - alpha) * n_before + alpha * n_d
            
            if np.max(abs((m - m_before) / m_before)) < err_long and\
               np.max(abs((n - n_before) / n_before)) < err_long:
                break
            
            m_before = m
            n_before = n
        
        print("long_max_value:", max_value)
        print("long_max_value:", k)
        
        m_true, n_true = m, n
        
        print("Lyapunov:", self.Lyapunov(m_true, n_true, RW))

        return m_true, n_true, RW, long_iteration, obj_list, obj_rel_list
    
    def logit(self, R, W, m):
        '''長期均衡の探索方向を計算する関数

        Parameters
        ----------
        mn: numpy.ndarray(K + K * K, )
            家計分布と企業分布の結合リスト
        RW: numpy.ndarray (K, )
            地代と賃金の結合リスト

        Returns
        -------
        mn_d: numpy.ndarray (K + K * K, )
            企業分布の探索方向と家計分布の探索方向の結合リスト
        '''
        
        pi = self.short.pi(R, W, m)
        v = self.short.v(R, W)
        
#         theta_v = self.prm.theta_house * v
#         max_v = np.max(theta_v)
#         diff_v = theta_v - max_v
        
#         log_logit = theta_v - np.log(np.sum(np.exp(diff_v))) - max_v
#         n_logit = np.exp(log_logit)

#         n_d = self.prm.N * n_logit

        m_d = self.prm.M \
        * (np.exp(self.prm.theta_firm * pi) \
        / np.sum(np.exp(self.prm.theta_firm * pi)))

        n_d = self.prm.N \
        * np.exp(self.prm.theta_house * v) \
        / np.sum(np.exp(self.prm.theta_house * v))

        return m_d, n_d

    def long_armijo(self, Z_LP, dZ, mn, mn_d, RW, c_1=0.5, beta=0.6):
        
        t = 1.0

        while True:
            print("Armijo")
            if Z_LP(mn + t * mn_d, RW) < Z_LP(mn, RW) + c_1 * t * dZ(mn, RW)@mn_d:
                break
            t *= beta
            
        return t
    
    def Lyapunov(self, m, n, RW):
        
        R = RW[:self.prm.K]
        W = RW[self.prm.K:]
        
        pi = self.short.pi(R, W, m)
        v = self.short.v(R, W)
        
        G = np.sum(pi * m) + np.sum(v * n)\
        - (self.prm.M * np.log(np.sum(np.exp(self.prm.theta_firm * pi))) / self.prm.theta_firm) - (self.prm.N * np.log(np.sum(np.exp(self.prm.theta_house * v))) / self.prm.theta_house) \
        - (np.sum(m * np.log(m / self.prm.M)) / self.prm.theta_firm) - (np.sum(n * np.log(n / self.prm.N)) / self.prm.theta_house)
        
        return G


# void short_solve(
#     size_t K, double S_bar,
#     double alpha_P1, double alpha_P2, double beta_P1, double beta_P2,
#     double coef_R_alpha, double coef_W_alpha, double coef_R_beta, double coef_W_beta,
#     double coef_pi, double coef_v,
#     double nE[], double m[], double **T_n,
#     double RW_proj, double p_proj, double par_L, double eta,
#     double R_hist[], double W_hist[], int *g_out,
#     double err_short, int short_itr)
# {
#     double *RW_before = (double *)malloc(2 * K * sizeof(double));
#     double *RW = (double *)malloc(2 * K * sizeof(double));
#     double *data_R = (double *)malloc(K * sizeof(double));
#     double *data_W = (double *)malloc(K * sizeof(double));
#     double *R = (double *)malloc(K * sizeof(double));
#     double *W = (double *)malloc(K * sizeof(double));
#     double *data_RW = NULL;
#     double data_Z = 0;
#     double *p_bar_before = (double *)malloc(2 * K * sizeof(double));
#     double *p_bar = (double *)malloc(2 * K * sizeof(double));
#     double L_before = par_L;
#     double t_before = 1.0;

#     for (size_t i = 0; i < 2 * K; i++)
#     {
#         RW_before[i] = 1.0; // 初期値を 1.0 に設定
#         RW[i] = 0.0;
#         p_bar_before[i] = 1.0;
#     }

#     int g = 1;

#     double L = 0;
#     double Z_SD_p_bar = 0;

#     double *dR = (double *)malloc(K * sizeof(double));
#     double *dW = (double *)malloc(K * sizeof(double));
#     for (size_t i = 0; i < K; i++)
#     {
#         dR[i] = 0.0;
#         dW[i] = 0.0;
#     }
#     double *dRdW = (double *)malloc(2 * K * sizeof(double));

#     double *dR_dot = (double *)malloc(K * sizeof(double));
#     double *dW_dot = (double *)malloc(K * sizeof(double));
#     double *dRdW_dot = (double *)malloc(2 * K * sizeof(double));

#     double *RW_diff = (double *)malloc(2 * K * sizeof(double));

#     // 要素数取得&データ読み込み
#     int valid_elements = read_RW(K, par_L, eta, &data_RW, &data_Z);

#     for (size_t i = 0; i < K; ++i)
#     {
#         data_R[i] = data_RW[i];
#         data_W[i] = data_RW[K + i];
#     }

#     double shiken = 0.0;
#     for (size_t i = 0; i < 2 * K; i++)
#     {
#         if (data_RW[i] > shiken)
#         {
#             shiken = data_RW[i];
#         }
#     }

#     printf("shiken: %.10f\n", shiken);
#     printf("data_Z: %.10f\n", data_Z);

#     // ディレクトリ作成
#     // char dirpath[256];
#     // if (create_parameter_directory("C:\\Users", K, par_L, eta, dirpath) != 0)
#     // {
#     //     fprintf(stderr, "Error creating directory\n");
#     //     exit(EXIT_FAILURE);
#     // }

#     // // CSVファイルパス
#     // char csv_filepath[512];
#     // snprintf(csv_filepath, sizeof(csv_filepath), "%s\\iteration_data.csv", dirpath);

#     // // CSVファイルパス
#     // // char csv_filepath[512];
#     // // snprintf(csv_filepath, sizeof(csv_filepath), "%s\\10000iteration.csv", dirpath);

#     // // CSVファイルを開く
#     // FILE *fp = fopen(csv_filepath, "w"); // "a" は追記モード, "w"は上書きモード
#     // if (fp == NULL)
#     // {
#     //     perror("Error opening file");
#     //     exit(EXIT_FAILURE);
#     // }

#     // ★★★ ヘッダーの書き込みは、ファイルが新規作成された場合のみ ★★★
#     // ファイルポインタの位置が0 (ファイルの先頭) なら、ファイルは空
#     // if (ftell(fp) == 0)
#     // { 
#     //     fprintf(fp, "RW,Z\n");
#     // }

#     // if (ftell(fp) == 0)
#     // {
#     //     fprintf(fp,"iteration,max_diff,R_max_diff,W_max_diff,Z_diff\n");
#     // }

#     for (int k = 0; k < short_itr; k++)
#     {
#         for (size_t i = 0; i < K; i++)
#         {
#             dR[i] = 0.0;
#             dW[i] = 0.0;
#         }

#         // Step 1: Backtracking
#         Z_SD_p_bar = effi_grad_Z(K, S_bar,
#                                  alpha_P1, alpha_P2, beta_P1, beta_P2,
#                                  coef_R_alpha, coef_W_alpha, coef_R_beta, coef_W_beta,
#                                  coef_pi, coef_v,
#                                  p_bar_before, nE, m, T_n, dR, dW);

#         L = backtracking(K, S_bar,
#                          alpha_P1, alpha_P2, beta_P1, beta_P2,
#                          coef_R_alpha, coef_W_alpha, coef_R_beta, coef_W_beta,
#                          coef_pi, coef_v,
#                          nE, m, T_n,
#                          Z_SD_p_bar, p_proj, eta, L_before, dR, dW, p_bar_before);

#         for (size_t i = 0; i < K; i++)
#         {
#             dRdW[i] = dR[i];
#             dRdW[i + K] = dW[i];
#         }

#         // Step 2: 解の更新
#         for (size_t i = 0; i < 2 * K; i++)
#         {
#             RW[i] = fmax(RW_proj, p_bar_before[i] - (dRdW[i] / L));
#         }

#         // Step 3: 収束判定
#         double max_diff = 0.0;
#         double R_max_diff = 0.0;
#         double W_max_diff = 0.0;

#         double Z_now = Z_SD(K, S_bar, coef_pi, coef_v,
#                         alpha_P1, alpha_P2, beta_P1, beta_P2,
#                         RW, nE, m, T_n, dR, dW);

#         for (size_t i = 0; i < K; ++i)
#         {
#             R[i] = RW[i];
#             W[i] = RW[K + i];
#         }

#         for (size_t i = 0; i < 2 * K; i++)
#         {
#             double diff = fabs((data_RW[i] - RW[i]) / RW[i]);
#             if (diff > max_diff)
#             {
#                 max_diff = diff;
#             }
#         }

#         for (size_t i = 0; i < K; i++)
#         {
#             double R_diff = fabs((data_R[i] - R[i]) / R[i]);
#             if (R_diff > R_max_diff)
#             {
#                 R_max_diff = R_diff;
#             }
#         }

#         for (size_t i = 0; i < K; i++)
#         {
#             double W_diff = fabs((data_W[i] - W[i]) / W[i]);
#             if (W_diff > W_max_diff)
#             {
#                 W_max_diff = W_diff;
#             }
#         }

#         double Z_diff = fabs((data_Z - Z_now) / Z_now);

#         // fprintf(fp, "%d,%.16f,%.16f,%.16f,%.16f\n",
#         //     k + 1, max_diff, R_max_diff, W_max_diff, Z_diff);

#         // Step 4: Adaptive restart
#         short_dual_df(
#             K, S_bar,
#             alpha_P1, alpha_P2, beta_P1, beta_P2,
#             coef_R_alpha, coef_W_alpha, coef_R_beta, coef_W_beta,
#             RW_before, nE, m, T_n, dR_dot, dW_dot);

#         for (size_t i = 0; i < K; i++)
#         {
#             dRdW_dot[i] = dR_dot[i];
#             dRdW_dot[i + K] = dW_dot[i];
#         }

#         for (size_t i = 0; i < 2 * K; i++)
#         {
#             RW_diff[i] = RW[i] - RW_before[i];
#         }

#         double RW_dot = dot_product(2 * K, dRdW_dot, RW_diff);

#         if (RW_dot > 0)
#         {
#             t_before = 1.0;
#         }

#         // Step 5: momentum項の計算
#         double t = (1.0 + sqrt(1.0 + 4.0 * t_before * t_before)) / 2.0;

#         for (size_t i = 0; i < 2 * K; i++)
#         {
#             p_bar[i] = fmax(RW_proj, RW[i] + ((t_before - 1.0) / t) * (RW[i] - RW_before[i]));
#         }

#         // RW_before, t_before の更新
#         for (size_t i = 0; i < 2 * K; i++)
#         {
#             RW_before[i] = RW[i];
#         }
#         for (size_t i = 0; i < 2 * K; i++)
#         {
#             p_bar_before[i] = p_bar[i];
#         }
#         t_before = t;
#         L_before = L;

#         g++;

#         if (g == short_itr + 1)
#         {
#             for (size_t i = 0; i < 2 * K; i++)
#             {
#                 if (i < K)
#                 {
#                     R_hist[i] = RW[i];
#                 }
#                 else
#                 {
#                     W_hist[i - K] = RW[i];
#                 }
#             }
#         }
#     }

#     // double Z_result = Z_SD(K, S_bar, coef_pi, coef_v,
#     //                        alpha_P1, alpha_P2, beta_P1, beta_P2,
#     //                        RW, nE, m, T_n, dR, dW);

#     // // RWの値を縦に記録
#     // for (size_t i = 0; i < 2 * K; i++)
#     // {
#     //     fprintf(fp, "%.16f", RW[i]); // RWの値を記録

#     //     // 最初の行のみ Z を記録
#     //     if (i == 0)
#     //     {
#     //         fprintf(fp, ",%.16f", Z_result); // Z を記録
#     //     }

#     //     fprintf(fp, "\n"); // 改行
#     // }

#     // fclose(fp);

#     *g_out = g;

#     equilibrium(
#         K, S_bar,
#         alpha_P1, alpha_P2, beta_P1, beta_P2,
#         coef_R_alpha, coef_W_alpha, coef_R_beta, coef_W_beta,
#         RW,
#         nE, m,
#         T_n);

#     printf("g: %d\n", g);

#     // メモリの解放
#     free(RW_before);
#     free(RW);
#     free(data_RW);
#     free(p_bar_before);
#     free(p_bar);
#     free(dR);
#     free(dW);
#     free(dRdW);
#     free(dR_dot);
#     free(dW_dot);
#     free(dRdW_dot);
#     free(RW_diff);
# }

# int main()
# {
#     // 定数パラメータの定義
#     const int E = 5;
#     const double M = 1.0;
#     const double N = 1.0;
#     const double alter_T_num = 0.5;
#     const double S_total = 100;
#     const double t = 0.1;
#     const double alpha_1 = 0.4;
#     const double alpha_2 = 0.4;
#     const double beta_1 = 0.4;
#     const double beta_2 = 0.4;
#     const double p_proj = 1e-5;
#     const double RW_proj = 1e-5;
#     const double err_short = 1e-5;
#     const int short_itr = 1000;

#     //--- 数値実験のパラメータと結果の配列 (例) ---
#     int K_list[] = {100};
#     double par_L[] = {0.200};
#     double eta[] = {1.5000};

#     int num_K = sizeof(K_list) / sizeof(K_list[0]);
#     int num_L = sizeof(par_L) / sizeof(par_L[0]);
#     int num_eta = sizeof(eta) / sizeof(eta[0]);

#     clock_t start, end;
#     double cpu_time_used;

#     for (int pk = 0; pk < num_K; pk++)
#     {
#         const size_t Col = sqrt(K_list[pk]);
#         const size_t K = Col * Col;
#         const int int_Col = Col;
#         const int int_K = int_Col * int_Col;
#         const double Scaling = 10.0 / int_Col;
#         const double S_bar = S_total / int_K;

#         // 動的配列の確保
#         double *R_hist = (double *)malloc(K * sizeof(double));
#         double *W_hist = (double *)malloc(K * sizeof(double));
#         double *RW_hist = (double *)malloc(2 * K * sizeof(double));
#         double *data_RW = (double *)malloc(2 * K * sizeof(double));
#         double *m0 = (double *)malloc(K * sizeof(double));

#         // メモリ確保のエラーチェック
#         if (R_hist == NULL || W_hist == NULL || RW_hist == NULL || data_RW == NULL || m0 == NULL)
#         {
#             fprintf(stderr, "メモリの確保に失敗しました\n");
#             exit(EXIT_FAILURE);
#         }

#         // 動的領域に確保
#         double **Coordinate_Data = malloc(K * sizeof(double *));
#         double **distance_matrix = malloc(K * sizeof(double *));
#         double **T = malloc(K * sizeof(double *));
#         double **T_n = malloc(K * sizeof(double *));
#         double **n0 = malloc(K * sizeof(double *));

#         for (size_t i = 0; i < K; i++)
#         {
#             Coordinate_Data[i] = malloc(2 * sizeof(double));
#             distance_matrix[i] = malloc(K * sizeof(double));
#             T[i] = malloc(K * sizeof(double));
#             T_n[i] = malloc(K * sizeof(double));
#             n0[i] = malloc(K * sizeof(double));
#         }

#         // 座標データの生成
#         for (int i = 0; i < int_K; i++)
#         {
#             Coordinate_Data[i][0] = (i % int_Col) * Scaling;
#             Coordinate_Data[i][1] = (i / int_Col) * Scaling;
#         }

#         // 距離行列の作成
#         for (size_t i = 0; i < K; i++)
#         {
#             for (size_t j = 0; j < K; j++)
#             {
#                 double dx = Coordinate_Data[i][0] - Coordinate_Data[j][0];
#                 double dy = Coordinate_Data[i][1] - Coordinate_Data[j][1];
#                 distance_matrix[i][j] = sqrt(dx * dx + dy * dy);
#             }
#         }

#         // Compute alpha and beta related parameters
#         double alpha_12 = 1.0 - alpha_1 - alpha_2;
#         double alpha_inv_12 = 1.0 / alpha_12;
#         double alpha_P1 = alpha_1 * alpha_inv_12;
#         double alpha_P2 = alpha_2 * alpha_inv_12;

#         double alpha_1_P1 = pow(alpha_1, alpha_P1);
#         double alpha_2_P2 = pow(alpha_2, alpha_P2);
#         double alpha_1_PP1 = pow(alpha_1, alpha_P1 + 1);
#         double alpha_2_PP2 = pow(alpha_2, alpha_P2 + 1);

#         double coef_R_alpha = alpha_2_P2 * alpha_1_PP1;
#         double coef_W_alpha = alpha_1_P1 * alpha_2_PP2;

#         double beta_12 = 1.0 - beta_1 - beta_2;
#         double beta_inv_12 = 1.0 / beta_12;
#         double beta_P1 = beta_1 * beta_inv_12;
#         double beta_P2 = beta_2 * beta_inv_12;

#         double beta_1_P1 = pow(beta_1, beta_P1);
#         double beta_2_P2 = pow(beta_2, beta_P2);
#         double beta_1_PP1 = pow(beta_1, beta_P1 + 1);
#         double beta_2_PP2 = pow(beta_2, beta_P2 + 1);

#         double coef_R_beta = beta_2_P2 * beta_1_PP1;
#         double coef_W_beta = beta_1_P1 * beta_2_PP2;

#         double coef_pi = beta_12 * beta_1_P1 * beta_2_P2;
#         double coef_v = alpha_12 * alpha_1_P1 * alpha_2_P2;

#         double *nE = (double *)malloc(K * sizeof(double));

#         // メモリ確保のエラーチェック
#         if (nE == NULL)
#         {
#             fprintf(stderr, "メモリの確保に失敗しました\n");
#             exit(EXIT_FAILURE);
#         }

#         for (size_t i = 0; i < K; i++)
#         {
#             for (size_t j = 0; j < K; j++)
#             {
#                 n0[i][j] = N / (int_K * int_K);
#             }
#         }

#         // T行列とD行列の計算
#         for (size_t i = 0; i < K; i++)
#         {
#             for (size_t j = 0; j < K; j++)
#             {
#                 T[i][j] = fmax(Scaling * t * alter_T_num, t * distance_matrix[i][j]);
#             }
#         }

#         // Compute T_n
#         for (size_t i = 0; i < K; i++)
#         {
#             for (size_t j = 0; j < K; j++)
#             {
#                 T_n[i][j] = pow(1.0 / exp(T[i][j]), alpha_inv_12) * n0[i][j];
#             }
#         }

#         // Compute nE
#         for (size_t j = 0; j < K; j++)
#         {
#             nE[j] = 0.0;
#             for (size_t i = 0; i < K; i++)
#             {
#                 nE[j] += E * n0[i][j];
#             }
#         }

#         // 初期化
#         double m_per = M / int_K;
#         for (size_t i = 0; i < K; i++)
#         {
#             m0[i] = m_per;
#         }

#         for (size_t i = 0; i < K; i++)
#         {
#             R_hist[i] = 0.0;
#             W_hist[i] = 0.0;
#         }

#         int g_out = 0;
#         int *g_out_ptr = &g_out;

#         for (int p1 = 0; p1 < num_L; p1++)
#         {
#             for (int p2 = 0; p2 < num_eta; p2++)
#             {
#                 printf("K: %d\n", K_list[pk]);
#                 printf("par_L: %f\n", par_L[p1]);
#                 printf("eta: %f\n", eta[p2]);

#                 start = clock();

#                 short_solve(
#                     K, S_bar,
#                     alpha_P1, alpha_P2, beta_P1, beta_P2,
#                     coef_R_alpha, coef_W_alpha, coef_R_beta, coef_W_beta,
#                     coef_pi, coef_v,
#                     nE, m0, T_n,
#                     RW_proj, p_proj, par_L[p1], eta[p2],
#                     R_hist, W_hist, g_out_ptr,
#                     err_short, short_itr);

#                 end = clock();
#                 cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
#                 printf("CPU time used: %f seconds\n", cpu_time_used);
#             }
#         }

#         // メモリの解放
#         for (size_t i = 0; i < K; i++)
#         {
#             free(Coordinate_Data[i]);
#             free(distance_matrix[i]);
#             free(T[i]);
#             free(T_n[i]);
#             free(n0[i]);
#         }

#         free(Coordinate_Data);
#         free(distance_matrix);
#         free(T);
#         free(T_n);
#         free(n0);

#         free(R_hist);
#         free(W_hist);
#         free(RW_hist);
#         free(data_RW);
#         free(m0);
#         free(nE);
#     }

#     return 0;
# }