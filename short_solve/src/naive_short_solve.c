#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

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

// `short_dual_df`関数
void short_dual_df(
    size_t K, double S_bar, double E,
    double alpha1, double alpha2, double beta1, double beta2,
    const double RW[], double m[], double **n,
    double **T,
    double dR[], double dW[])
{
    double R[K], W[K];
    double S_F[K], L_F[K];

    // RWをRとWに分割
    for (size_t i = 0; i < K; ++i)
    {
        R[i] = RW[i];
        W[i] = RW[K + i];
    }

    for (int i = 0; i < K; i++)
    {
        dR[i] = S_bar;
        double sum1 = 0.0, sum2 = 0.0;

        for (int j = 0; j < K; j++)
        {
            double T_factor = pow(1.0 / exp(T[i][j]), 1.0 / (1.0 - alpha1 - alpha2));
            double R_factor = pow(alpha1 / R[i], (1.0 - alpha2) / (1.0 - alpha1 - alpha2));
            double W_factor = pow(alpha2 / W[j], alpha2 / (1.0 - alpha1 - alpha2));

            sum1 += T_factor * R_factor * W_factor * n[i][j];
        }
        S_F[i] =
        pow(beta1 / R[i], (1 - beta2) / (1 - beta1 - beta2))
        * pow(beta2 / W[i], beta2 / (1 - beta1 - beta2));
        dR[i] -= sum1 + S_F[i] * m[i];
    }

    for (int j = 0; j < K; j++)
    {
        dW[j] = 0.0;
        double sum3 = 0.0, sum4 = 0.0;

        for (int i = 0; i < K; i++)
        {
            double T_factor = pow(1.0 / exp(T[i][j]), 1.0 / (1.0 - alpha1 - alpha2));
            double R_factor = pow(alpha1 / R[i], alpha1 / (1.0 - alpha1 - alpha2));
            double W_factor = pow(alpha2 / W[j], (1.0 - alpha1) / (1.0 - alpha1 - alpha2));

            sum3 += T_factor * R_factor * W_factor * n[i][j];
            sum4 += n[i][j];
        }
        L_F[j] =
        pow(beta1 / R[j], beta1 / (1 - beta1 - beta2))
        * pow(beta2 / W[j], (1 - beta2) / (1 - beta1 - beta2));

        dW[j] = E * sum4 - sum3 - L_F[j] * m[j];
    }
}

void pi_noex(size_t K, double pi[], double R[], double W[],
             double beta_L_R, double beta_S_W, double power_L_F_R, double power_S_F_W, double pi_beta)
{
    for (size_t i = 0; i < K; i++)
    {
        double R_pi = beta_L_R * pow(R[i], -power_L_F_R);
        double W_pi = beta_S_W * pow(W[i], -power_S_F_W);
        pi[i] = pi_beta * R_pi * W_pi;
    }
}

double Z_SD(size_t K, double E, double S_bar, double RW[], double m[], double **n, double **T,
            double alpha1, double alpha2, double beta1, double beta2)
{
    double Z_SD = 0.0;
    double R[K], W[K];

    // RWをRとWに分割
    for (size_t i = 0; i < K; ++i)
    {
        R[i] = RW[i];
        W[i] = RW[K + i];
    }

    // 第一の和
    for (int i = 0; i < K; i++)
    {
        double term1 = (1 - beta1 - beta2) * pow(beta1 / R[i], beta1 / (1 - beta1 - beta2)) * pow(beta2 / W[i], beta2 / (1 - beta1 - beta2));
        Z_SD += term1 * m[i] + S_bar * R[i];
    }

    // 第二の二重和
    for (int i = 0; i < K; i++)
    {
        for (int j = 0; j < K; j++)
        {
            double term2 = (1 - alpha1 - alpha2) * pow(1 / exp(T[i][j]), 1 / (1 - alpha1 - alpha2)) * pow(alpha1 / R[i], alpha1 / (1 - alpha1 - alpha2)) * pow(alpha2 / W[j], alpha2 / (1 - alpha1 - alpha2));
            Z_SD += (E * W[j] + term2) * n[i][j];
        }
    }

    return Z_SD;
}

double backtracking(
    size_t K, double S_bar, double E,
    double alpha1, double alpha2, double beta1, double beta2,
    double **T, double Z_SD_p_bar,
    double p_proj, double eta,
    double L, double dR[], double dW[], double p_bar[], double m_fixed[], double **n_fixed)
{
    double L_bar = L;
    double dRdW[2 * K];
    double p[2 * K];
    double diff_p[2 * K];
    int itr = 0;

    for (size_t i = 0; i < K; i++)
    {
        dRdW[i] = dR[i];
        dRdW[K + i] = dW[i];
    }

    for (size_t k = 0; k < 1000000; k++)
    {
        for (size_t j = 0; j < 2 * K; j++)
        {
            p[j] = fmax(p_proj, p_bar[j] - dRdW[j] / L_bar);
        }

        // 判定条件
        double Z_SD_p = Z_SD(K, E, S_bar, p, m_fixed, n_fixed, T,
                             alpha1, alpha2, beta1, beta2);

        double norm_squared = 0.0;
        for (size_t j = 0; j < 2 * K; j++)
        {
            double diff = p[j] - p_bar[j];
            norm_squared += diff * diff;
            diff_p[j] = diff;
        }

        double dot = dot_product(2 * K, diff_p, dRdW);

        if (Z_SD_p - (Z_SD_p_bar + dot + 0.5 * L_bar * norm_squared) <= 0.0)
        {
            break;
        }

        // 更新
        L_bar *= eta; // algprm.eta として最後の要素を使用
        itr++;
    }

    return L_bar;
}

void short_solve(
    size_t K, double E, double S_bar, double RW_proj, double p_proj, double par_L, double eta,
    double R_hist[], double W_hist[], int *g_out,
    double m_fixed[], double **n_fixed,
    double alpha1, double alpha2, double beta1, double beta2,
    double **T,
    double err_short, int short_itr)
{
    double RW_before[2 * K];
    double RW[2 * K];
    double p_bar_before[2 * K];
    double p_bar[2 * K];
    double L_before = par_L;

    double t_before = 1.0;

    for (size_t i = 0; i < 2 * K; i++)
    {
        RW_before[i] = 1.0; // 初期値を 1.0 に設定
        RW[i] = 0.0;
        p_bar_before[i] = 1.0;
    }

    int g = 1;

    // v_cal の結果を確認
    // if (v_cal_minimum(R_ini, W_ini, K) < 0)
    // {
    //     fprintf(stderr, "Error: v must be non-negative\n");
    //     exit(EXIT_FAILURE);
    // }
    double L = 0;
    double Z_SD_p_bar = 0;

    double dR[K];
    double dW[K];
    for (size_t i = 0; i < K; i++)
    {
        dR[i] = 0.0;
        dW[i] = 0.0;
    }
    double dRdW[2 * K];

    double dR_dot[K];
    double dW_dot[K];
    double dRdW_dot[2 * K];

    double RW_diff[2 * K];

    for (int k = 0; k < short_itr; k++)
    {
        for (size_t i = 0; i < K; i++)
        {
            dR[i] = 0.0;
            dW[i] = 0.0;
        }

        // short_dual_df の結果を取得
        short_dual_df(K, S_bar, E,
                      alpha1, alpha2, beta1, beta2,
                      p_bar_before, m_fixed, n_fixed, T, dR, dW);

        // short_dual_df の結果を取得
        Z_SD_p_bar = Z_SD(K, E, S_bar, p_bar_before, m_fixed, n_fixed, T,
                        alpha1, alpha2, beta1, beta2);

        // Step 1: Backtracking
        L = backtracking(
            K, S_bar, E,
            alpha1, alpha2, beta1, beta2,
            T, Z_SD_p_bar, p_proj, eta, L_before,
            dR, dW, p_bar_before, m_fixed, n_fixed);

        for (size_t i = 0; i < K; i++)
        {
            dRdW[i] = dR[i];
            dRdW[i + K] = dW[i];
        }

        // Step 2: 解の更新
        for (size_t i = 0; i < 2 * K; i++)
        {
            RW[i] = fmax(RW_proj, p_bar_before[i] - (dRdW[i] / L));
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

        // Step 4: Adaptive restart
        short_dual_df(K, S_bar, E,
                    alpha1, alpha2, beta1, beta2,
                    RW_before, m_fixed, n_fixed, T, dR_dot, dW_dot);

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

    printf("g: %d\n", g);
}

int main()
{
    // 定数パラメータの定義
    const int E = 5;
    const size_t Col = 30;
    const size_t K = Col * Col;
    const int int_Col = 30;
    const int int_K = int_Col * int_Col;
    const double M = 1.0;
    const double N = 1.0;
    const double Scaling = 10.0 / int_Col;
    const double alter_T_num = 0.5;
    const double S_total = 100;
    const double S_bar = S_total / int_K;
    const double t = 0.1;
    const double alpha1 = 0.4;
    const double alpha2 = 0.4;
    const double beta1 = 0.4;
    const double beta2 = 0.4;
    const double par_L = 0.2;
    const double eta = 1.2;
    const double p_proj = 1e-3;
    const double RW_proj = 1e-3;
    const double err_short = 1e-5;
    const int short_itr = 50;

    double R_hist[K];
    double W_hist[K];
    double m0[K];

    // 動的領域に確保
    double **Coordinate_Data = malloc(K * sizeof(double *));
    double **distance_matrix = malloc(K * sizeof(double *));
    double **T = malloc(K * sizeof(double *));
    double **n0 = malloc(K * sizeof(double *));

    for (size_t i = 0; i < K; i++)
    {
        Coordinate_Data[i] = malloc(2 * sizeof(double));
        distance_matrix[i] = malloc(K * sizeof(double));
        T[i] = malloc(K * sizeof(double));
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

    // T行列とD行列の計算
    for (size_t i = 0; i < K; i++)
    {
        for (size_t j = 0; j < K; j++)
        {
            T[i][j] = fmax(Scaling * t * alter_T_num, t * distance_matrix[i][j]);
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
        for (size_t j = 0; j < K; j++)
        {
            n0[i][j] = N / (int_K * int_K);
        }
    }

    for (size_t i = 0; i < K; i++)
    {
        R_hist[i] = 0.0;
        W_hist[i] = 0.0;
    }

    int g_out = 0;
    int *g_out_ptr = &g_out;

    clock_t start, end;
    double cpu_time_used;

    start = clock();

    short_solve(
        K, E, S_bar, RW_proj, p_proj, par_L, eta,
        R_hist, W_hist, g_out_ptr,
        m0, n0,
        alpha1, alpha2, beta1, beta2,
        T,
        err_short, short_itr);

    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("CPU time used: %f seconds\n", cpu_time_used);

    // メモリの解放
    for (size_t i = 0; i < K; i++)
    {
        free(Coordinate_Data[i]);
        free(distance_matrix[i]);
        free(T[i]);
        free(n0[i]);
    }

    free(Coordinate_Data);
    free(distance_matrix);
    free(T);
    free(n0);

    return 0;
}