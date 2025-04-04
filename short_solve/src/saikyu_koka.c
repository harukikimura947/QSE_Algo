#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <direct.h>
#include <errno.h>
#include <stdbool.h>

#define MAX_ELEMENTS 1800
#define LINE_LENGTH 65536 // 1行分の最大長（大きめにする）
#define FILENAME_TEMPLATE "C://Users//kimura//Downloads//short_solve//data//K=900//par_L=%.4f,eta=%.4f//20000iteration.csv"

// ベクトル要素ごとの累乗
void elementwise_pow(size_t size, double result[], const double vec[], double power)
{
    for (size_t i = 0; i < size; ++i)
    {
        result[i] = pow(vec[i], power);
    }
}

// ベクトル要素ごとの積
void elementwise_mul_vec(size_t size, double result[], const double vec1[], const double vec2[])
{
    for (size_t i = 0; i < size; ++i)
    {
        result[i] = vec1[i] * vec2[i];
    }
}

// 行列要素ごとの積 (2次元配列の引数)ただし，同じ大きさの行列同士に限る．
void elementwise_mul_mat(size_t size, double **result, double **mat1, double **mat2)
{
    for (size_t i = 0; i < size; ++i)
    {
        for (size_t j = 0; j < size; ++j)
        {
            result[i][j] = mat1[i][j] * mat2[i][j];
        }
    }
}

// 行列とベクトルの行方向の要素積と総和
double sum_row_product(size_t size, double **mat, const double vec[])
{
    double sum = 0.0;
    for (size_t j = 0; j < size; ++j)
    {
        sum += mat[size][j] * vec[j];
    }
    return sum;
}

// 行列とベクトルの列方向の要素積と総和
double sum_col_product(size_t size, double **mat, const double vec[])
{
    double sum = 0.0;
    for (size_t i = 0; i < size; ++i)
    {
        sum += mat[i][size] * vec[i];
    }
    return sum;
}

// ベクトルの総和
double sum_vector(size_t size, const double vec[])
{
    double sum = 0.0;
    for (size_t i = 0; i < size; ++i)
    {
        sum += vec[i];
    }
    return sum;
}

// ベクトルの内積
double dot_product(size_t size, const double a[], const double b[])
{
    double result = 0.0;
    for (size_t i = 0; i < size; i++)
    {
        result += a[i] * b[i];
    }
    return result;
}

// ベクトルの外積（行ベクトル×列ベクトル）
void outer_product(size_t size, const double vec1[], const double vec2[], double **result)
{
    for (size_t i = 0; i < size; ++i)
    {
        for (size_t j = 0; j < size; ++j)
        {
            result[i][j] = vec1[i] * vec2[j];
        }
    }
}

void short_dual_df(
    size_t K, double S_bar,
    double alpha_P1, double alpha_P2, double beta_P1, double beta_P2,
    double coef_R_alpha, double coef_W_alpha, double coef_R_beta, double coef_W_beta,
    double RW[],
    double nE[], double m[],
    double **T_n,
    double dR[], double dW[])
{
    double R[K];
    double W[K];
    double R_inv[K];
    double W_inv[K];
    double R_inv_alpha[K];
    double W_inv_alpha[K];
    double R_inv_beta[K];
    double W_inv_beta[K];
    double R_T[K];
    double W_T[K];

    for (size_t i = 0; i < K; ++i)
    {
        R[i] = RW[i];
        W[i] = RW[K + i];
    }

    // Compute R_inv, W_inv
    for (size_t i = 0; i < K; i++)
    {
        R_inv[i] = 1.0 / R[i];
    }
    for (size_t j = 0; j < K; j++)
    {
        W_inv[j] = 1.0 / W[j];
    }

    // Compute R_inv_alpha, W_inv_alpha, R_inv_beta, W_inv_beta
    for (size_t i = 0; i < K; i++)
    {
        R_inv_alpha[i] = pow(R_inv[i], alpha_P1);
        R_inv_beta[i] = pow(R_inv[i], beta_P1);
    }
    for (size_t j = 0; j < K; j++)
    {
        W_inv_alpha[j] = pow(W_inv[j], alpha_P2);
        W_inv_beta[j] = pow(W_inv[j], beta_P2);
    }

    // Compute R_T and W_T
    for (size_t j = 0; j < K; j++)
    {
        R_T[j] = 0.0;
        for (size_t i = 0; i < K; i++)
        {
            R_T[j] += R_inv_alpha[i] * T_n[i][j];
        }
    }

    for (size_t i = 0; i < K; i++)
    {
        W_T[i] = 0.0;
        for (size_t j = 0; j < K; j++)
        {
            W_T[i] += W_inv_alpha[j] * T_n[i][j];
        }
    }

    // Compute gradient dZ/dR
    for (size_t i = 0; i < K; i++)
    {
        dR[i] = S_bar
        - (coef_R_alpha * R_inv[i] * R_inv_alpha[i] * W_T[i])
        - (coef_R_beta * R_inv[i] * R_inv_beta[i] * W_inv_beta[i] * m[i]);
    }

    // Compute gradient dZ/dW
    for (size_t j = 0; j < K; j++)
    {
        dW[j] = nE[j]
        - (coef_W_alpha * W_inv[j] * W_inv_alpha[j] * R_T[j])
        - (coef_W_beta * W_inv[j] * W_inv_beta[j] * R_inv_beta[j] * m[j]);
    }
}

void equilibrium(
    size_t K, double S_bar,
    double alpha_P1, double alpha_P2, double beta_P1, double beta_P2,
    double coef_R_alpha, double coef_W_alpha, double coef_R_beta, double coef_W_beta,
    double RW[],
    double nE[], double m[],
    double **T_n)
{
    double dR_equ[K];
    double dW_equ[K];

    short_dual_df(
        K, S_bar,
        alpha_P1, alpha_P2, beta_P1, beta_P2,
        coef_R_alpha, coef_W_alpha, coef_R_beta, coef_W_beta,
        RW,
        nE, m,
        T_n,
        dR_equ, dW_equ);

    int dR_size = sizeof(dR_equ) / sizeof(dR_equ[0]);
    int dW_size = sizeof(dW_equ) / sizeof(dW_equ[0]);

    // 最大値を求める
    double dR_max = fabs(dR_equ[0]); // 最初の要素を仮の最大値とする
    for (int i = 1; i < dR_size; i++)
    {
        if (fabs(dR_equ[i]) > dR_max)
        {
            dR_max = fabs(dR_equ[i]);
        }
    }
    printf("dRの最大値: %.16f\n", dR_max);

    // 平均値を求める
    double dR_sum = 0.0;
    for (int i = 0; i < dR_size; i++)
    {
        dR_sum += dR_equ[i];
    }
    double dR_average = dR_sum / dR_size;
    printf("dRの平均値: %.16f\n", dR_average); // 小数点以下2桁まで表示

    // 最大値を求める
    double dW_max = fabs(dW_equ[0]); // 最初の要素を仮の最大値とする
    for (int i = 1; i < dW_size; i++)
    {
        if (fabs(dW_equ[i]) > dW_max)
        {
            dW_max = fabs(dW_equ[i]);
        }
    }
    printf("dWの最大値: %.16f\n", dW_max);

    // 平均値を求める
    double dW_sum = 0.0;
    for (int i = 0; i < dW_size; i++)
    {
        dW_sum += dW_equ[i];
    }
    double dW_average = dW_sum / dW_size;
    printf("dWの平均値: %.16f\n", dW_average); // 小数点以下2桁まで表示
}

double effi_grad_Z(
    size_t K, double S_bar,
    double alpha_P1, double alpha_P2, double beta_P1, double beta_P2,
    double coef_R_alpha, double coef_W_alpha, double coef_R_beta, double coef_W_beta,
    double coef_pi, double coef_v,
    double RW[],
    double nE[], double m[],
    double **T_n,
    double dR[], double dW[])
{

    double R[K];
    double W[K];
    double R_inv[K];
    double W_inv[K];
    double R_inv_alpha[K];
    double W_inv_alpha[K];
    double R_inv_beta[K];
    double W_inv_beta[K];
    double R_T[K];
    double W_T[K];

    for (size_t i = 0; i < K; ++i)
    {
        R[i] = RW[i];
        W[i] = RW[K + i];
    }

    // Compute R_inv, W_inv
    for (size_t i = 0; i < K; i++)
    {
        R_inv[i] = 1.0 / R[i];
    }
    for (size_t j = 0; j < K; j++)
    {
        W_inv[j] = 1.0 / W[j];
    }

    // Compute R_inv_alpha, W_inv_alpha, R_inv_beta, W_inv_beta
    for (size_t i = 0; i < K; i++)
    {
        R_inv_alpha[i] = pow(R_inv[i], alpha_P1);
        R_inv_beta[i] = pow(R_inv[i], beta_P1);
    }

    for (size_t j = 0; j < K; j++)
    {
        W_inv_alpha[j] = pow(W_inv[j], alpha_P2);
        W_inv_beta[j] = pow(W_inv[j], beta_P2);
    }

    // Compute R_T and W_T
    for (size_t j = 0; j < K; j++)
    {
        R_T[j] = 0.0;
        for (size_t i = 0; i < K; i++)
        {
            R_T[j] += R_inv_alpha[i] * T_n[i][j];
        }
    }

    for (size_t i = 0; i < K; i++)
    {
        W_T[i] = 0.0;
        for (size_t j = 0; j < K; j++)
        {
            W_T[i] += W_inv_alpha[j] * T_n[i][j];
        }
    }

    for (size_t i = 0; i < K; i++)
    {
        dR[i] = S_bar
        - (coef_R_alpha * R_inv[i] * R_inv_alpha[i] * W_T[i])
        - (coef_R_beta * R_inv[i] * R_inv_beta[i] * W_inv_beta[i] * m[i]);
    }

    // Compute gradient dZ/dW
    for (size_t j = 0; j < K; j++)
    {
        dW[j] = nE[j]
        - (coef_W_alpha * W_inv[j] * W_inv_alpha[j] * R_T[j])
        - (coef_W_beta * W_inv[j] * W_inv_beta[j] * R_inv_beta[j] * m[j]);
    }

    // Compute Z_SD
    double Z_SD = 0.0;
    double Z_SD_1 = 0.0;
    double Z_SD_2 = 0.0;
    double Z_SD_3 = 0.0;
    double Z_SD_4 = 0.0;

    // First summation: beta-related term
    for (size_t i = 0; i < K; i++)
    {
        Z_SD_1 += R_inv_beta[i] * W_inv_beta[i] * m[i];
    }

    Z_SD_1 *= coef_pi;

    // Second summation: S_bar * R
    for (size_t i = 0; i < K; i++)
    {
        Z_SD_2 += R[i];
    }

    Z_SD_2 *= S_bar;

    // Third summation: W * nE
    for (size_t j = 0; j < K; j++)
    {
        Z_SD_3 += W[j] * nE[j];
    }

    // Fourth summation: alpha-related term
    for (size_t i = 0; i < K; i++)
    {
        Z_SD_4 += R_inv_alpha[i] * W_T[i];
    }

    Z_SD_4 *= coef_v;

    Z_SD = Z_SD_1 + Z_SD_2 + Z_SD_3 + Z_SD_4;

    return Z_SD;
}

double Z_SD(
    size_t K, double S_bar,
    double coef_pi, double coef_v,
    double alpha_P1, double alpha_P2, double beta_P1, double beta_P2,
    double RW[],
    double nE[], double m[],
    double **T_n)
{

    double R[K];
    double W[K];
    double R_inv[K];
    double W_inv[K];
    double R_inv_alpha[K];
    double W_inv_alpha[K];
    double R_inv_beta[K];
    double W_inv_beta[K];
    double R_T[K];
    double W_T[K];

    for (size_t i = 0; i < K; ++i)
    {
        R[i] = RW[i];
        W[i] = RW[K + i];
    }

    // Compute R_inv, W_inv
    for (size_t i = 0; i < K; i++)
    {
        R_inv[i] = 1.0 / R[i];
    }
    for (size_t j = 0; j < K; j++)
    {
        W_inv[j] = 1.0 / W[j];
    }

    // Compute R_inv_alpha, W_inv_alpha, R_inv_beta, W_inv_beta
    for (size_t i = 0; i < K; i++)
    {
        R_inv_alpha[i] = pow(R_inv[i], alpha_P1);
        R_inv_beta[i] = pow(R_inv[i], beta_P1);
    }
    for (size_t j = 0; j < K; j++)
    {
        W_inv_alpha[j] = pow(W_inv[j], alpha_P2);
        W_inv_beta[j] = pow(W_inv[j], beta_P2);
    }

    for (size_t i = 0; i < K; i++)
    {
        W_T[i] = 0.0;
        for (size_t j = 0; j < K; j++)
        {
            W_T[i] += W_inv_alpha[j] * T_n[i][j];
        }
    }

    // Compute Z_SD
    double Z_SD = 0.0;
    double Z_SD_1 = 0.0;
    double Z_SD_2 = 0.0;
    double Z_SD_3 = 0.0;
    double Z_SD_4 = 0.0;

    // First summation: beta-related term
    for (size_t i = 0; i < K; i++)
    {
        Z_SD_1 += R_inv_beta[i] * W_inv_beta[i] * m[i];
    }

    Z_SD_1 *= coef_pi;

    // Second summation: S_bar * R
    for (size_t i = 0; i < K; i++)
    {
        Z_SD_2 += R[i];
    }

    Z_SD_2 *= S_bar;

    // Third summation: W * nE
    for (size_t j = 0; j < K; j++)
    {
        Z_SD_3 += W[j] * nE[j];
    }

    // Fourth summation: alpha-related term
    for (size_t i = 0; i < K; i++)
    {
        Z_SD_4 += R_inv_alpha[i] * W_T[i];
    }

    Z_SD_4 *= coef_v;

    Z_SD = Z_SD_1 + Z_SD_2 + Z_SD_3 + Z_SD_4;

    return Z_SD;
}

double armijo(
    size_t K, double S_bar,
    double alpha_P1, double alpha_P2, double beta_P1, double beta_P2,
    double coef_pi, double coef_v,
    double nE[], double m[], double **T_n,
    double RW[], double dRdW[], double Z_before)
{

    double step = 1.0;
    double c_1 = 0.5;
    double beta = 0.8;
    double RW_next[2 * K];
    double d[2 * K];
    double Z_next = 0.0;

    for (size_t i = 0; i < 2 * K; ++i)
    {
        d[i] = - dRdW[i];
    }

    double d_dot = dot_product(2 * K, dRdW, d);

    while(true)
    {
    for (size_t i = 0; i < 2 * K; ++i)
    {
        RW_next[i] = RW[i] + step * d[i];
    }

    Z_next = Z_SD(
        K, S_bar,
        coef_pi, coef_v,
        alpha_P1, alpha_P2, beta_P1, beta_P2,
        RW_next,
        nE, m,
        T_n);

    if (Z_next < Z_before + c_1 * step * d_dot)
    {
        break;
    }
    step *= beta;
    }

    return step;
}

void short_solve(
    size_t K, double S_bar,
    double alpha_P1, double alpha_P2, double beta_P1, double beta_P2,
    double coef_R_alpha, double coef_W_alpha, double coef_R_beta, double coef_W_beta,
    double coef_pi, double coef_v,
    double nE[], double m[], double **T_n,
    double RW_proj, double p_proj, double par_L, double eta,
    double R_hist[], double W_hist[], int *g_out,
    double err_short, int short_itr)
{
    double RW_before[2 * K];
    double RW[2 * K];

    double Z_before = 0.0;

    for (size_t i = 0; i < 2 * K; i++)
    {
        RW_before[i] = 1.0; // 初期値を 1.0 に設定
        RW[i] = 0.0;
    }

    int g = 1;

    // v_cal の結果を確認
    // if (v_cal_minimum(R_ini, W_ini, K) < 0)
    // {
    //     fprintf(stderr, "Error: v must be non-negative\n");
    //     exit(EXIT_FAILURE);
    // }

    double dR[K];
    double dW[K];
    for (size_t i = 0; i < K; i++)
    {
        dR[i] = 0.0;
        dW[i] = 0.0;
    }
    double dRdW[2 * K];

    double RW_diff[2 * K];

    for (int k = 0; k < short_itr; k++)
    {
        for (size_t i = 0; i < K; i++)
        {
            dR[i] = 0.0;
            dW[i] = 0.0;
        }

        Z_before = effi_grad_Z(K, S_bar,
                            alpha_P1, alpha_P2, beta_P1, beta_P2,
                            coef_R_alpha, coef_W_alpha, coef_R_beta, coef_W_beta,
                            coef_pi, coef_v,
                            RW_before, nE, m, T_n, dR, dW);

        for (size_t i = 0; i < K; i++)
        {
            dRdW[i] = dR[i];
            dRdW[i + K] = dW[i];
        }

        // double step; // 事前に宣言しておく

        double step = armijo(
            K, S_bar,
            alpha_P1, alpha_P2, beta_P1, beta_P2,
            coef_pi, coef_v,
            nE, m, T_n,
            RW_before, dRdW, Z_before);

        // if (g <= 200)
        // {
        //     step = armijo(
        //         K, S_bar,
        //         alpha_P1, alpha_P2, beta_P1, beta_P2,
        //         coef_pi, coef_v,
        //         nE, m, T_n,
        //         RW_before, dRdW, Z_before);
        // }
        // else
        // { // else を使うと、すべてのケースをカバーできる
        //     step = 0.01;
        // }

        // Step 2: 解の更新
        for (size_t i = 0; i < 2 * K; i++)
        {
            RW[i] = fmax(RW_proj, RW_before[i] - step * dRdW[i]);
        }

        // Step 3: 収束判定
        double max_diff = 0.0;
        for (size_t i = 0; i < 2 * K; i++)
        {
            double diff = fabs((RW[i] - RW_before[i]) / RW_before[i]);

            if (diff > max_diff)
            {
                max_diff = diff;
            }
        }

        if (max_diff < err_short)
        {
            // RW_box に RW を追加 (ここでは RW_hist を仮使用)
            for (size_t i = 0; i < 2 * K; i++)
            {
                if (i < K)
                {
                    R_hist[i] = RW[i];
                }
                else
                {
                    W_hist[i - K] = RW[i];
                }
            }
            break;
        }

        // RW_beforeの更新
        for (size_t i = 0; i < 2 * K; i++)
        {
            RW_before[i] = RW[i];
        }

        g++;

        if (g == short_itr + 1)
        {
            // RW_box に RW を追加 (ここでは RW_hist を仮使用)
            for (size_t i = 0; i < 2 * K; i++)
            {
                if (i < K)
                {
                    R_hist[i] = RW[i];
                }
                else
                {
                    W_hist[i - K] = RW[i];
                }
            }
        }
    }

    *g_out = g;

    equilibrium(
        K, S_bar,
        alpha_P1, alpha_P2, beta_P1, beta_P2,
        coef_R_alpha, coef_W_alpha, coef_R_beta, coef_W_beta,
        RW,
        nE, m,
        T_n);

    printf("g: %d\n", g);
}

int main()
{
    // 定数パラメータの定義
    const int E = 5;
    const size_t Col = 4;
    const size_t K = Col * Col;
    const int int_Col = 4;
    const int int_K = int_Col * int_Col;
    const double M = 1.0;
    const double N = 1.0;
    const double Scaling = 10.0 / int_Col;
    const double alter_T_num = 0.5;
    const double S_total = 100;
    const double S_bar = S_total / int_K;
    const double t = 0.1;
    const double alpha_1 = 0.4;
    const double alpha_2 = 0.4;
    const double beta_1 = 0.4;
    const double beta_2 = 0.4;
    // const double par_L = 0.2;
    // const double eta = 1.8;
    const double p_proj = 1e-5;
    const double RW_proj = 1e-5;
    const double err_short = 1e-6;
    const int short_itr = 1000;

    double R_hist[K];
    double W_hist[K];
    double RW_hist[2 * K];
    double data_RW[2 * K];
    double m0[K];

    // 動的領域に確保
    double **Coordinate_Data = malloc(K * sizeof(double *));
    double **distance_matrix = malloc(K * sizeof(double *));
    double **T = malloc(K * sizeof(double *));
    double **T_n = malloc(K * sizeof(double *));
    double **n0 = malloc(K * sizeof(double *));

    for (size_t i = 0; i < K; i++)
    {
        Coordinate_Data[i] = malloc(2 * sizeof(double));
        distance_matrix[i] = malloc(K * sizeof(double));
        T[i] = malloc(K * sizeof(double));
        T_n[i] = malloc(K * sizeof(double));
        n0[i] = malloc(K * sizeof(double));
    }

    // 座標データの生成
    for (int i = 0; i < int_K; i++)
    {
        Coordinate_Data[i][0] = (i % int_Col) * Scaling;
        Coordinate_Data[i][1] = (i / int_Col) * Scaling;
    }

    // 距離行列の作成
    for (size_t i = 0; i < K; i++)
    {
        for (size_t j = 0; j < K; j++)
        {
            double dx = Coordinate_Data[i][0] - Coordinate_Data[j][0];
            double dy = Coordinate_Data[i][1] - Coordinate_Data[j][1];
            distance_matrix[i][j] = sqrt(dx * dx + dy * dy);
        }
    }

    // Compute alpha and beta related parameters
    double alpha_12 = 1.0 - alpha_1 - alpha_2;
    double alpha_inv_12 = 1.0 / alpha_12;
    double alpha_P1 = alpha_1 * alpha_inv_12;
    double alpha_P2 = alpha_2 * alpha_inv_12;
    
    double alpha_1_P1 = pow(alpha_1, alpha_P1);
    double alpha_2_P2 = pow(alpha_2, alpha_P2);
    double alpha_1_PP1 = pow(alpha_1, alpha_P1 + 1);
    double alpha_2_PP2 = pow(alpha_2, alpha_P2 + 1);

    double coef_R_alpha = alpha_2_P2 * alpha_1_PP1;
    double coef_W_alpha = alpha_1_P1 * alpha_2_PP2;

    double beta_12 = 1.0 - beta_1 - beta_2;
    double beta_inv_12 = 1.0 / beta_12;
    double beta_P1 = beta_1 * beta_inv_12;
    double beta_P2 = beta_2 * beta_inv_12;

    double beta_1_P1 = pow(beta_1, beta_P1);
    double beta_2_P2 = pow(beta_2, beta_P2);
    double beta_1_PP1 = pow(beta_1, beta_P1 + 1);
    double beta_2_PP2 = pow(beta_2, beta_P2 + 1);

    double coef_R_beta = beta_2_P2 * beta_1_PP1;
    double coef_W_beta = beta_1_P1 * beta_2_PP2;

    double coef_pi = beta_12 * beta_1_P1 * beta_2_P2;
    double coef_v = alpha_12 * alpha_1_P1 * alpha_2_P2;

    double nE[K];

    for (size_t i = 0; i < K; i++)
    {
        for (size_t j = 0; j < K; j++)
        {
            n0[i][j] = N / (int_K * int_K);
        }
    }

    // T行列とD行列の計算
    for (size_t i = 0; i < K; i++)
    {
        for (size_t j = 0; j < K; j++)
        {
            T[i][j] = fmax(Scaling * t * alter_T_num, t * distance_matrix[i][j]);
        }
    }

    // Compute T_n
    for (size_t i = 0; i < K; i++)
    {
        for (size_t j = 0; j < K; j++)
        {
            T_n[i][j] = pow(1.0 / exp(T[i][j]), alpha_inv_12) * n0[i][j];
        }
    }

    // Compute nE
    for (size_t j = 0; j < K; j++)
    {
        nE[j] = 0.0;
        for (size_t i = 0; i < K; i++)
        {
            nE[j] += E * n0[i][j];
        }
    }

    // 初期化
    double m_per = M / int_K;
    for (size_t i = 0; i < K; i++)
    {
        m0[i] = m_per;
    }

    for (size_t i = 0; i < K; i++)
    {
        R_hist[i] = 0.0;
        W_hist[i] = 0.0;
    }

    int g_out = 0;
    int *g_out_ptr = &g_out;

    //--- 数値実験のパラメータと結果の配列 (例) ---
    // double par_L[] = {0.0050, 0.0080, 0.0100, 0.0200, 0.0300, 0.0500, 0.0800};
    // double eta[] = {1.5000, 1.8000, 2.0000};

    double par_L[] = {0.050};
    double eta[] = {1.8000};

    int num_L = sizeof(par_L) / sizeof(par_L[0]);
    int num_eta = sizeof(eta) / sizeof(eta[0]);

    clock_t start, end;
    double cpu_time_used;

    start = clock();

    for (int p1 = 0; p1 < num_L; p1++)
    {
        for (int p2 = 0; p2 < num_eta; p2++)
        {

            short_solve(
                K, S_bar,
                alpha_P1, alpha_P2, beta_P1, beta_P2,
                coef_R_alpha, coef_W_alpha, coef_R_beta, coef_W_beta,
                coef_pi, coef_v,
                nE, m0, T_n,
                RW_proj, p_proj, par_L[p1], eta[p2],
                R_hist, W_hist, g_out_ptr,
                err_short, short_itr);

            for (size_t i = 0; i < K; i++)
            {
                printf("R: %.10f\n", R_hist[i]);
            }

            for (size_t i = 0; i < K; i++)
            {
                printf("W: %.10f\n", W_hist[i]);
            }

            for (size_t i = 0; i < K; i++)
            {
                RW_hist[i] = R_hist[i];
                RW_hist[i + K] = W_hist[i];
            }
        }
    }

    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("CPU time used: %f seconds\n", cpu_time_used);

    // メモリの解放
    for (size_t i = 0; i < K; i++)
    {
        free(Coordinate_Data[i]);
        free(distance_matrix[i]);
        free(T[i]);
        free(T_n[i]);
        free(n0[i]);
    }

    free(Coordinate_Data);
    free(distance_matrix);
    free(T);
    free(T_n);
    free(n0);

    return 0;
}