#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <direct.h>
#include <errno.h>

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

int main()
{
    // 定数パラメータの定義
    const int E = 5;
    const double M = 1.0;
    const double N = 1.0;
    const double alter_T_num = 0.5;
    const double S_total = 100;
    const double t = 0.1;
    const double alpha_1 = 0.4;
    const double alpha_2 = 0.4;
    const double beta_1 = 0.4;
    const double beta_2 = 0.4;
    const double p_proj = 1e-5;
    const double RW_proj = 1e-5;
    const double err_short = 1e-5;
    const int short_itr = 1000;

    //--- 数値実験のパラメータと結果の配列 (例) ---
    int K_list[] = {100};
    double par_L[] = {0.20};
    double eta[] = {1.5000};

    int num_K = sizeof(K_list) / sizeof(K_list[0]);
    int num_L = sizeof(par_L) / sizeof(par_L[0]);
    int num_eta = sizeof(eta) / sizeof(eta[0]);

    clock_t start, end;
    double cpu_time_used;

    for (int pk = 0; pk < num_K; pk++)
    {
        const size_t Col = sqrt(K_list[pk]);
        const size_t K = Col * Col;
        printf("K: %d\n", K);
        printf("2 * K: %d\n", 2 * K);
        const int int_Col = Col;
        const int int_K = int_Col * int_Col;
        const double Scaling = 10.0 / int_Col;
        const double S_bar = S_total / int_K;

        // 動的配列の確保
        double *R_hist = (double *)malloc(K * sizeof(double));
        double *W_hist = (double *)malloc(K * sizeof(double));
        double *RW_hist = (double *)malloc(2 * K * sizeof(double));
        double *m0 = (double *)malloc(K * sizeof(double));

        // メモリ確保のエラーチェック
        if (R_hist == NULL || W_hist == NULL || RW_hist == NULL || m0 == NULL)
        {
            fprintf(stderr, "メモリの確保に失敗しました\n");
            exit(EXIT_FAILURE);
        }

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

        double *nE = (double *)malloc(K * sizeof(double));

        // メモリ確保のエラーチェック
        if (nE == NULL)
        {
            fprintf(stderr, "メモリの確保に失敗しました\n");
            exit(EXIT_FAILURE);
        }

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

        for (int p1 = 0; p1 < num_L; p1++)
        {
            for (int p2 = 0; p2 < num_eta; p2++)
            {
                printf("K: %d\n", K_list[pk]);
                printf("par_L: %f\n", par_L[p1]);
                printf("eta: %f\n", eta[p2]);

                start = clock();

                double *RW_before = (double *)malloc(2 * K * sizeof(double));
                double *RW = (double *)malloc(2 * K * sizeof(double));
                double *R = (double *)malloc(K * sizeof(double));
                double *W = (double *)malloc(K * sizeof(double));
                double *p_bar_before = (double *)malloc(2 * K * sizeof(double));
                double *p_bar = (double *)malloc(2 * K * sizeof(double));
                double L_before = par_L[p1];
                double t_before = 1.0;

                for (size_t i = 0; i < 2 * K; i++)
                {
                    RW_before[i] = 1.0; // 初期値を 1.0 に設定
                    RW[i] = 0.0;
                    p_bar_before[i] = 1.0;
                }

                int g = 1;

                double L = 0;
                double Z_SD_p_bar = 0;

                double *dR = (double *)malloc(K * sizeof(double));
                double *dW = (double *)malloc(K * sizeof(double));
                for (size_t i = 0; i < K; i++)
                {
                    dR[i] = 0.0;
                    dW[i] = 0.0;
                }
                double *dRdW = (double *)malloc(2 * K * sizeof(double));

                double *dR_dot = (double *)malloc(K * sizeof(double));
                double *dW_dot = (double *)malloc(K * sizeof(double));
                double *dRdW_dot = (double *)malloc(2 * K * sizeof(double));

                double *RW_diff = (double *)malloc(2 * K * sizeof(double));

                for (int k = 0; k < short_itr; k++)
                {

                    // effi_grad_Z関数====================================================================
                    // 動的配列の確保

                    // double *R_func = (double *)malloc(K * sizeof(double));
                    // double *W_func = (double *)malloc(K * sizeof(double));
                    // double *R_inv = (double *)malloc(K * sizeof(double));
                    // double *W_inv = (double *)malloc(K * sizeof(double));
                    // double *R_inv_alpha = (double *)malloc(K * sizeof(double));
                    // double *R_inv_beta = (double *)malloc(K * sizeof(double));
                    // double *W_inv_alpha = (double *)malloc(K * sizeof(double));
                    // double *W_inv_beta = (double *)malloc(K * sizeof(double));
                    // double *R_T = (double *)malloc(K * sizeof(double));
                    // double *W_T = (double *)malloc(K * sizeof(double));

                    // // メモリ確保のエラーチェック
                    // if (R_func == NULL || W_func == NULL || R_inv == NULL || W_inv == NULL ||
                    //     R_inv_alpha == NULL || R_inv_beta == NULL || W_inv_alpha == NULL ||
                    //     W_inv_beta == NULL || R_T == NULL || W_T == NULL)
                    // {
                    //     fprintf(stderr, "メモリの確保に失敗しました\n");
                    //     exit(EXIT_FAILURE);
                    // }

                    double R_func[K];
                    double W_func[K];
                    double R_inv[K];
                    double W_inv[K];
                    double R_inv_alpha[K];
                    double R_inv_beta[K];
                    double W_inv_alpha[K];
                    double W_inv_beta[K];
                    double R_T[K];
                    double W_T[K];

                    // R と W を RW から初期化
                    for (size_t i = 0; i < K; ++i)
                    {
                        R_func[i] = p_bar_before[i];
                        W_func[i] = p_bar_before[K + i];
                    }

                    // Compute R_inv, W_inv
                    for (size_t i = 0; i < K; i++)
                    {
                        R_inv[i] = 1.0 / R_func[i];
                    }
                    for (size_t j = 0; j < K; j++)
                    {
                        W_inv[j] = 1.0 / W_func[j];
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
                        dR[i] = S_bar - (coef_R_alpha * R_inv[i] * R_inv_alpha[i] * W_T[i]) - (coef_R_beta * R_inv[i] * R_inv_beta[i] * W_inv_beta[i] * m0[i]);
                    }

                    // Compute gradient dZ/dW
                    for (size_t j = 0; j < K; j++)
                    {
                        dW[j] = nE[j] - (coef_W_alpha * W_inv[j] * W_inv_alpha[j] * R_T[j]) - (coef_W_beta * W_inv[j] * W_inv_beta[j] * R_inv_beta[j] * m0[j]);
                    }

                    // Compute Z_SD
                    double Z_SD_p_bar = 0.0;
                    double Z_SD_1 = 0.0;
                    double Z_SD_2 = 0.0;
                    double Z_SD_3 = 0.0;
                    double Z_SD_4 = 0.0;

                    // First summation: beta-related term
                    for (size_t i = 0; i < K; i++)
                    {
                        Z_SD_1 += R_inv_beta[i] * W_inv_beta[i] * m0[i];
                    }

                    Z_SD_1 *= coef_pi;

                    // Second summation: S_bar * R
                    for (size_t i = 0; i < K; i++)
                    {
                        Z_SD_2 += R_func[i];
                    }

                    Z_SD_2 *= S_bar;

                    // Third summation: W * nE
                    for (size_t j = 0; j < K; j++)
                    {
                        Z_SD_3 += W_func[j] * nE[j];
                    }

                    // Fourth summation: alpha-related term
                    for (size_t i = 0; i < K; i++)
                    {
                        Z_SD_4 += R_inv_alpha[i] * W_T[i];
                    }

                    Z_SD_4 *= coef_v;

                    Z_SD_p_bar = Z_SD_1 + Z_SD_2 + Z_SD_3 + Z_SD_4;

                    // free(R_func);
                    // free(W_func);
                    // free(R_inv);
                    // free(W_inv);
                    // free(R_inv_alpha);
                    // free(W_inv_alpha);
                    // free(R_inv_beta);
                    // free(W_inv_beta);
                    // free(R_T);
                    // free(W_T);

                    // effi_grad_Z関数終了====================================================================

                    // backtracking関数======================================================================
                    double L_bar = L_before;

                    for (size_t i = 0; i < 2 * K; i++)
                    {
                        p_bar[i] = p_bar_before[i];
                    }

                    // 動的配列の確保
                    double *dRdW_back = (double *)malloc(2 * K * sizeof(double));
                    double *p = (double *)malloc(2 * K * sizeof(double));
                    double *diff_p = (double *)malloc(2 * K * sizeof(double));

                    if (dRdW_back == NULL || p == NULL || diff_p == NULL)
                    {
                        fprintf(stderr, "メモリの確保に失敗しました\n");
                        exit(EXIT_FAILURE);
                    }

                    int itr = 0;

                    for (size_t i = 0; i < K; i++)
                    {
                        dRdW_back[i] = dR[i];
                        dRdW_back[K + i] = dW[i];
                    }

                    for (size_t k = 0; k < 1000000; k++)
                    {
                        for (size_t j = 0; j < 2 * K; j++)
                        {
                            p[j] = fmax(p_proj, p_bar[j] - dRdW_back[j] / L_bar);
                        }

                        // 判定条件
                        // Z_SD関数==================================================================
                        // double *R_func2 = (double *)malloc(K * sizeof(double));
                        // double *W_func2 = (double *)malloc(K * sizeof(double));
                        // double *R_inv2 = (double *)malloc(K * sizeof(double));
                        // double *W_inv2 = (double *)malloc(K * sizeof(double));
                        // double *R_inv_alpha2 = (double *)malloc(K * sizeof(double));
                        // double *R_inv_beta2 = (double *)malloc(K * sizeof(double));
                        // double *W_inv_alpha2 = (double *)malloc(K * sizeof(double));
                        // double *W_inv_beta2 = (double *)malloc(K * sizeof(double));
                        // double *R_T2 = (double *)malloc(K * sizeof(double));
                        // double *W_T2 = (double *)malloc(K * sizeof(double));

                        // // メモリ確保のエラーチェック
                        // if (R_func2 == NULL || W_func2 == NULL || R_inv2 == NULL || W_inv2 == NULL ||
                        //     R_inv_alpha2 == NULL || R_inv_beta2 == NULL || W_inv_alpha2 == NULL ||
                        //     W_inv_beta2 == NULL || R_T2 == NULL || W_T2 == NULL)
                        // {
                        //     fprintf(stderr, "メモリの確保に失敗しました\n");
                        //     exit(EXIT_FAILURE);
                        // }

                        double R_func2[K];
                        double W_func2[K];
                        double R_inv2[K];
                        double W_inv2[K];
                        double R_inv_alpha2[K];
                        double R_inv_beta2[K];
                        double W_inv_alpha2[K];
                        double W_inv_beta2[K];
                        double R_T2[K];
                        double W_T2[K];

                        // R と W を RW から初期化
                        for (size_t i = 0; i < K; ++i)
                        {
                            R_func2[i] = p[i];
                            W_func2[i] = p[K + i];
                        }

                        // Compute R_inv, W_inv
                        for (size_t i = 0; i < K; i++)
                        {
                            R_inv2[i] = 1.0 / R_func2[i];
                        }
                        for (size_t j = 0; j < K; j++)
                        {
                            W_inv2[j] = 1.0 / W_func2[j];
                        }

                        // Compute R_inv_alpha, W_inv_alpha, R_inv_beta, W_inv_beta
                        for (size_t i = 0; i < K; i++)
                        {
                            R_inv_alpha2[i] = pow(R_inv2[i], alpha_P1);
                            R_inv_beta2[i] = pow(R_inv2[i], beta_P1);
                        }
                        for (size_t j = 0; j < K; j++)
                        {
                            W_inv_alpha2[j] = pow(W_inv2[j], alpha_P2);
                            W_inv_beta2[j] = pow(W_inv2[j], beta_P2);
                        }

                        // Compute W_T
                        for (size_t i = 0; i < K; i++)
                        {
                            W_T2[i] = 0.0;
                            for (size_t j = 0; j < K; j++)
                            {
                                W_T2[i] += W_inv_alpha2[j] * T_n[i][j];
                            }
                        }

                        // Compute Z_SD
                        double Z_SD_p = 0.0;
                        double Z_SD_1 = 0.0;
                        double Z_SD_2 = 0.0;
                        double Z_SD_3 = 0.0;
                        double Z_SD_4 = 0.0;

                        // First summation: beta-related term
                        for (size_t i = 0; i < K; i++)
                        {
                            Z_SD_1 += R_inv_beta2[i] * W_inv_beta2[i] * m0[i];
                        }

                        Z_SD_1 *= coef_pi;

                        // Second summation: S_bar * R
                        for (size_t i = 0; i < K; i++)
                        {
                            Z_SD_2 += R_func2[i];
                        }

                        Z_SD_2 *= S_bar;

                        // Third summation: W * nE
                        for (size_t j = 0; j < K; j++)
                        {
                            Z_SD_3 += W_func2[j] * nE[j];
                        }

                        // Fourth summation: alpha-related term
                        for (size_t i = 0; i < K; i++)
                        {
                            Z_SD_4 += R_inv_alpha2[i] * W_T2[i];
                        }

                        Z_SD_4 *= coef_v;

                        Z_SD_p = Z_SD_1 + Z_SD_2 + Z_SD_3 + Z_SD_4;

                        // free(R_func2);
                        // free(W_func2);
                        // free(R_inv2);
                        // free(W_inv2);
                        // free(R_inv_alpha2);
                        // free(W_inv_alpha2);
                        // free(R_inv_beta2);
                        // free(W_inv_beta2);
                        // free(R_T2);
                        // free(W_T2);
                        // Z_SD関数終了===============================================================

                        double norm_squared = 0.0;
                        for (size_t j = 0; j < 2 * K; j++)
                        {
                            double diff = p[j] - p_bar[j];
                            norm_squared += diff * diff;
                            diff_p[j] = diff;
                        }

                        double dot = dot_product(2 * K, diff_p, dRdW_back);

                        if (Z_SD_p - (Z_SD_p_bar + dot + 0.5 * L_bar * norm_squared) <= 0.0)
                        {
                            break;
                        }

                        // 更新
                        L_bar *= eta[p2]; // algprm.eta として最後の要素を使用
                        itr++;
                    }

                    // 動的配列の解放
                    free(dRdW_back);
                    free(p);
                    free(diff_p);

                    L = L_bar;

                    // backtracking関数終了==================================================================

                    for (size_t i = 0; i < K; i++)
                    {
                        dRdW[i] = dR[i];
                        dRdW[i + K] = dW[i];
                    }

                    printf("L: %f\n", L);

                    // Step 2: 解の更新
                    for (size_t i = 0; i < 2 * K; i++)
                    {
                        RW[i] = fmax(RW_proj, p_bar_before[i] - (dRdW[i] / L));
                    }

                    // Step 3: 収束判定
                    // 数値実験の結果として、変数の相対誤差を測定する
                    double max_diff = 0.0;

                    for (size_t i = 0; i < 2 * K; i++)
                    {
                        double diff = fabs((RW[i] - RW_before[i]) / RW_before[i]);
                        if (diff > max_diff)
                        {
                            max_diff = diff;
                        }
                    }

                    if (max_diff < 1e-4)
                    {
                        break;
                    }

                    // Step 4: Adaptive restart
                    // short_dual_df関数===================================================================

                    // double *R_func3 = (double *)malloc(K * sizeof(double));
                    // double *W_func3 = (double *)malloc(K * sizeof(double));
                    // double *R_inv3 = (double *)malloc(K * sizeof(double));
                    // double *W_inv3 = (double *)malloc(K * sizeof(double));
                    // double *R_inv_alpha3 = (double *)malloc(K * sizeof(double));
                    // double *R_inv_beta3 = (double *)malloc(K * sizeof(double));
                    // double *W_inv_alpha3 = (double *)malloc(K * sizeof(double));
                    // double *W_inv_beta3 = (double *)malloc(K * sizeof(double));
                    // double *R_T3 = (double *)malloc(K * sizeof(double));
                    // double *W_T3 = (double *)malloc(K * sizeof(double));

                    // // メモリ確保のエラーチェック
                    // if (R_func3 == NULL || W_func3 == NULL || R_inv3 == NULL || W_inv3 == NULL ||
                    //     R_inv_alpha3 == NULL || R_inv_beta3 == NULL || W_inv_alpha3 == NULL ||
                    //     W_inv_beta3 == NULL || R_T3 == NULL || W_T3 == NULL)
                    // {
                    //     fprintf(stderr, "メモリの確保に失敗しました\n");
                    //     exit(EXIT_FAILURE);
                    // }

                    double R_func3[K];
                    double W_func3[K];
                    double R_inv3[K];
                    double W_inv3[K];
                    double R_inv_alpha3[K];
                    double R_inv_beta3[K];
                    double W_inv_alpha3[K];
                    double W_inv_beta3[K];
                    double R_T3[K];
                    double W_T3[K];

                    // R と W を RW から初期化
                    for (size_t i = 0; i < K; ++i)
                    {
                        R_func3[i] = RW_before[i];
                        W_func3[i] = RW_before[K + i];
                    }

                    // Compute R_inv, W_inv
                    for (size_t i = 0; i < K; i++)
                    {
                        R_inv3[i] = 1.0 / R_func3[i];
                    }
                    for (size_t j = 0; j < K; j++)
                    {
                        W_inv3[j] = 1.0 / W_func3[j];
                    }

                    // Compute R_inv_alpha, W_inv_alpha, R_inv_beta, W_inv_beta
                    for (size_t i = 0; i < K; i++)
                    {
                        R_inv_alpha3[i] = pow(R_inv3[i], alpha_P1);
                        R_inv_beta3[i] = pow(R_inv3[i], beta_P1);
                    }
                    for (size_t j = 0; j < K; j++)
                    {
                        W_inv_alpha3[j] = pow(W_inv3[j], alpha_P2);
                        W_inv_beta3[j] = pow(W_inv3[j], beta_P2);
                    }

                    // Compute R_T and W_T
                    for (size_t j = 0; j < K; j++)
                    {
                        R_T3[j] = 0.0;
                        for (size_t i = 0; i < K; i++)
                        {
                            R_T3[j] += R_inv_alpha3[i] * T_n[i][j];
                        }
                    }

                    for (size_t i = 0; i < K; i++)
                    {
                        W_T3[i] = 0.0;
                        for (size_t j = 0; j < K; j++)
                        {
                            W_T3[i] += W_inv_alpha3[j] * T_n[i][j];
                        }
                    }

                    // Compute gradient dZ/dR
                    for (size_t i = 0; i < K; i++)
                    {
                        dR[i] = S_bar - (coef_R_alpha * R_inv3[i] * R_inv_alpha3[i] * W_T3[i]) - (coef_R_beta * R_inv3[i] * R_inv_beta3[i] * W_inv_beta3[i] * m0[i]);
                    }

                    // Compute gradient dZ/dW
                    for (size_t j = 0; j < K; j++)
                    {
                        dW[j] = nE[j] - (coef_W_alpha * W_inv3[j] * W_inv_alpha3[j] * R_T3[j]) - (coef_W_beta * W_inv3[j] * W_inv_beta3[j] * R_inv_beta3[j] * m0[j]);
                    }

                    // 動的配列の解放
                    // free(R_func3);
                    // free(W_func3);
                    // free(R_inv3);
                    // free(W_inv3);
                    // free(R_inv_alpha3);
                    // free(W_inv_alpha3);
                    // free(R_inv_beta3);
                    // free(W_inv_beta3);
                    // free(R_T3);
                    // free(W_T3);

                    // short_dual_df関数終了============================================================

                    for (size_t i = 0; i < K; i++)
                    {
                        dRdW_dot[i] = dR_dot[i];
                        dRdW_dot[i + K] = dW_dot[i];
                    }

                    for (size_t i = 0; i < 2 * K; i++)
                    {
                        RW_diff[i] = RW[i] - RW_before[i];
                    }

                    double RW_dot = dot_product(2 * K, dRdW_dot, RW_diff);

                    if (RW_dot > 0)
                    {
                        t_before = 1.0;
                        printf("restart on \n");
                    }

                    // Step 5: momentum項の計算
                    double t = (1.0 + sqrt(1.0 + 4.0 * t_before * t_before)) / 2.0;

                    for (size_t i = 0; i < 2 * K; i++)
                    {
                        p_bar[i] = fmax(RW_proj, RW[i] + ((t_before - 1.0) / t) * (RW[i] - RW_before[i]));
                    }

                    // RW_before, t_before の更新
                    for (size_t i = 0; i < 2 * K; i++)
                    {
                        RW_before[i] = RW[i];
                    }
                    for (size_t i = 0; i < 2 * K; i++)
                    {
                        p_bar_before[i] = p_bar[i];
                    }
                    t_before = t;
                    L_before = L;

                    g++;

                    printf("g: %d\n", g);

                    if (g == short_itr + 1)
                    {
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

                // *g_out = g;

                printf("g: %d\n", g);

                // メモリの解放
                free(RW_before);
                free(R);
                free(W);
                free(RW);
                free(p_bar_before);
                free(p_bar);
                free(dR);
                free(dW);
                free(dRdW);
                free(dR_dot);
                free(dW_dot);
                free(dRdW_dot);
                free(RW_diff);

                end = clock();
                cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
                printf("CPU time used: %f seconds\n", cpu_time_used);
            }
        }

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

        free(R_hist);
        free(W_hist);
        free(RW_hist);
        free(m0);
        free(nE);
    }

    return 0;
}