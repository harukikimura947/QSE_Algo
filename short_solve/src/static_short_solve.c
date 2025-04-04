#include <stdio.h>
#include <math.h>
#include <time.h>

void compute_T_power(size_t K, double T_power[K][K], double T[K][K], double alpha_1, double alpha_2)
{
    double exponent = 1.0 / (1.0 - alpha_1 - alpha_2);

    for (size_t i = 0; i < K; i++)
    {
        for (size_t j = 0; j < K; j++)
        {
            T_power[i][j] = pow(1.0 / exp(T[i][j]), exponent);
        }
    }
}

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
void elementwise_mul_mat(size_t size, double result[size][size], const double mat1[size][size], const double mat2[size][size])
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
double sum_row_product(size_t size, const double mat[size][size], const double vec[])
{
    double sum = 0.0;
    for (size_t j = 0; j < size; ++j)
    {
        sum += mat[size][j] * vec[j];
    }
    return sum;
}

// 行列とベクトルの列方向の要素積と総和
double sum_col_product(size_t size, const double mat[size][size], const double vec[])
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

// // ベクトルの要素ごとの加算（結合）
// void concatenate_vectors(const double vec1[], size_t size1, const double vec2[], size_t size2, double result[], size_t *result_size)
// {
//     for (size_t i = 0; i < size1; ++i)
//     {
//         result[i] = vec1[i];
//     }
//     for (size_t i = 0; i < size2; ++i)
//     {
//         result[size1 + i] = vec2[i];
//     }
//     *result_size = size1 + size2;
// }

double dot_product_kahan(size_t size, const double a[], const double b[])
{
    double sum = 0.0;
    double c = 0.0; // 誤差補正用変数
    for (size_t i = 0; i < size; i++)
    {
        double prod = a[i] * b[i];
        double y = prod - c; // 誤差を補正した値
        double t = sum + y;  // 正しい加算結果
        c = (t - sum) - y;   // 新しい誤差
        sum = t;             // 更新
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
void outer_product(size_t size, const double vec1[], const double vec2[], double result[size][size])
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
    double alpha_S_R, double alpha_S_W, double alpha_L_R, double alpha_L_W,
    double power_S_H_R, double power_S_H_W, double power_L_H_R, double power_L_H_W,
    double beta_S_R, double beta_S_W, double beta_L_R, double beta_L_W,
    double power_S_F_R, double power_S_F_W, double power_L_F_R, double power_L_F_W,
    const double RW[], double m[], double n[K][K],
    double T_power[K][K],
    double dR[], double dW[])
{

    double R[K], W[K];
    double S_H_R[K], S_H_W[K], L_H_R[K], L_H_W[K];
    double S_F[K], L_F[K];

    // RWをRとWに分割
    for (size_t i = 0; i < K; ++i)
    {
        R[i] = RW[i];
        W[i] = RW[K + i];
    }

    // S_H_Rの計算
    elementwise_pow(K, S_H_R, R, -power_S_H_R);
    for (size_t i = 0; i < K; ++i)
    {
        S_H_R[i] *= alpha_S_R;
    }

    // S_H_Wの計算
    elementwise_pow(K, S_H_W, W, -power_S_H_W);

    // L_H_Rの計算
    elementwise_pow(K, L_H_R, R, -power_L_H_R);

    // L_H_Wの計算
    elementwise_pow(K, L_H_W, W, -power_L_H_W);
    for (size_t i = 0; i < K; ++i)
    {
        L_H_W[i] *= alpha_L_W;
    }

    // S_Fの計算
    double temp1[K], temp2[K];
    elementwise_pow(K, temp1, R, -power_S_F_R);
    elementwise_pow(K, temp2, W, -power_S_F_W);
    elementwise_mul_vec(K, S_F, temp1, temp2);
    for (size_t i = 0; i < K; ++i)
    {
        S_F[i] *= beta_S_R * beta_S_W;
    }

    // L_Fの計算
    elementwise_pow(K, temp1, R, -power_L_F_R);
    elementwise_pow(K, temp2, W, -power_L_F_W);
    elementwise_mul_vec(K, L_F, temp1, temp2);
    for (size_t i = 0; i < K; ++i)
    {
        L_F[i] *= beta_L_R * beta_L_W;
    }

    // dRの計算
    for (size_t i = 0; i < K; i++)
    {
        double Sn = 0.0;
        for (size_t j = 0; j < K; j++)
        {
            Sn += T_power[i][j] * S_H_W[j] * n[i][j];
        }
        dR[i] = S_bar - S_H_R[i] * alpha_S_W * Sn - S_F[i] * m[i];
    }

    // dWの計算
    for (size_t j = 0; j < K; j++)
    {
        double Ln = 0.0;
        double sum_n = 0.0;
        for (size_t i = 0; i < K; i++)
        {
            Ln += T_power[i][j] * L_H_R[i] * n[i][j];
            sum_n += n[i][j];
        }
        dW[j] = E * sum_n - L_H_W[j] * alpha_L_R * Ln - L_F[j] * m[j];
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

double Z_SD(size_t K, double E, double S_bar, double RW[], double m[], double n[K][K], double T_power_v[K][K],
            double beta_L_R, double beta_S_W, double power_L_F_R, double power_S_F_W, double pi_beta,
            double alpha_L_R, double alpha_S_W, double power_L_H_R, double power_S_H_W)
{
    double R[K], W[K], R_v[K], W_v[K], vn = 0;

    // RWからRとWに分離
    for (size_t i = 0; i < K; i++)
    {
        R[i] = RW[i];
        W[i] = RW[K + i];
    }

    // R_vとW_vの計算
    for (size_t i = 0; i < K; i++)
    {
        R_v[i] = alpha_L_R * pow(R[i], -power_L_H_R);
        W_v[i] = alpha_S_W * pow(W[i], -power_S_H_W);
    }

    // vnの計算
    for (size_t i = 0; i < K; i++)
    {
        for (size_t j = 0; j < K; j++)
        {
            vn += (E * W[j] + T_power_v[i][j] * R_v[i] * W_v[j]) * n[i][j];
        }
    }

    // F_valueの計算
    double F_value = vn;
    double pi_noex_result[K]; // pi_noex関数の戻り値を格納する配列

    pi_noex(K, pi_noex_result, R, W,
            beta_L_R, beta_S_W, power_L_F_R, power_S_F_W, pi_beta);

    // 総和の計算
    for (size_t i = 0; i < K; i++)
    {
        F_value += pi_noex_result[i] * m[i];
    }

    // 総和の計算
    for (size_t i = 0; i < K; i++)
    {
        F_value += S_bar * R[i];
    }

    return F_value;
}

double backtracking(
    size_t K, double S_bar, double E,
    double alpha_S_R, double alpha_S_W, double alpha_L_R, double alpha_L_W,
    double power_S_H_R, double power_S_H_W, double power_L_H_R, double power_L_H_W,
    double beta_S_R, double beta_S_W, double beta_L_R, double beta_L_W,
    double power_S_F_R, double power_S_F_W, double power_L_F_R, double power_L_F_W,
    double pi_beta,
    double T_power[K][K], double T_power_v[K][K],
    double p_proj, double eta,
    double L, double dR[], double dW[], double p_bar[], double m_fixed[], double n_fixed[K][K])
{
    double L_bar = L;
    double Z_SD_p_bar;
    double dRdW[2 * K];
    double p[2 * K];
    double diff_p[2 * K];
    int itr = 0;

    // dRdW 計算
    short_dual_df(K, S_bar, E,
                  alpha_S_R, alpha_S_W, alpha_L_R, alpha_L_W,
                  power_S_H_R, power_S_H_W, power_L_H_R, power_L_H_W,
                  beta_S_R, beta_S_W, beta_L_R, beta_L_W,
                  power_S_F_R, power_S_F_W, power_L_F_R, power_L_F_W,
                  p_bar, m_fixed, n_fixed, T_power, dR, dW);

    for (size_t i = 0; i < K; i++)
    {
        dRdW[i] = dR[i];
        dRdW[K + i] = dW[i];
    }

    // for (size_t i = 0; i < 2 * K; i++)
    // {
    //     printf("dRdW[i]: %f\n", dRdW[i]);
    // }

    // 初期 Z_SD(p_bar)
    Z_SD_p_bar = Z_SD(K, E, S_bar, p_bar, m_fixed, n_fixed, T_power_v,
                      beta_L_R, beta_S_W, power_L_F_R, power_S_F_W, pi_beta,
                      alpha_L_R, alpha_S_W, power_L_H_R, power_S_H_W);

    // printf("Z_SD_p_bar: %f\n", Z_SD_p_bar);

    for (size_t k = 0; k < 1000000; k++)
    {

        for (size_t j = 0; j < 2 * K; j++)
        {
            p[j] = fmax(p_proj, p_bar[j] - dRdW[j] / L_bar);
        }

        // 判定条件
        double Z_SD_p = Z_SD(K, E, S_bar, p, m_fixed, n_fixed, T_power_v,
                             beta_L_R, beta_S_W, power_L_F_R, power_S_F_W, pi_beta,
                             alpha_L_R, alpha_S_W, power_L_H_R, power_S_H_W);
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
            // printf("btra_iteration: %d\n", itr);
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
    double m_fixed[], double n_fixed[K][K],
    double alpha_S_R, double alpha_S_W, double alpha_L_R, double alpha_L_W,
    double power_S_H_R, double power_S_H_W, double power_L_H_R, double power_L_H_W,
    double beta_S_R, double beta_S_W, double beta_L_R, double beta_L_W,
    double power_S_F_R, double power_S_F_W, double power_L_F_R, double power_L_F_W,
    double pi_beta,
    double T_power[K][K], double T_power_v[K][K],
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
        // printf("========================================\n");
        // printf("iteration: %d\n", k);

        // printf("L_before: %f\n", L_before);

        // for (size_t i = 0; i < 2 * K; i++)
        // {
        //     printf("p_bar_before: %f\n", p_bar_before[i]);
        // }

        for (size_t i = 0; i < K; i++)
        {
            dR[i] = 0.0;
            dW[i] = 0.0;
        }

        // Step 1: Backtracking
        L = backtracking(K, S_bar, E,
                         alpha_S_R, alpha_S_W, alpha_L_R, alpha_L_W,
                         power_S_H_R, power_S_H_W, power_L_H_R, power_L_H_W,
                         beta_S_R, beta_S_W, beta_L_R, beta_L_W,
                         power_S_F_R, power_S_F_W, power_L_F_R, power_L_F_W, pi_beta, T_power, T_power_v,
                         p_proj, eta, L_before, dR, dW, p_bar_before, m_fixed, n_fixed);

        // printf("L: %f\n", L);

        // short_dual_df の結果を取得
        short_dual_df(K, S_bar, E,
                      alpha_S_R, alpha_S_W, alpha_L_R, alpha_L_W,
                      power_S_H_R, power_S_H_W, power_L_H_R, power_L_H_W,
                      beta_S_R, beta_S_W, beta_L_R, beta_L_W,
                      power_S_F_R, power_S_F_W, power_L_F_R, power_L_F_W,
                      p_bar_before, m_fixed, n_fixed, T_power, dR, dW);

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

        // for (size_t i = 0; i < 2 * K; i++)
        // {
        //     printf("dRdW: %.17f\n", dRdW[i]);
        // }

        // if (k == 0)
        // {
        // printf("配列の要素: ");
        // for (size_t i = 0; i < 2 * K; i++)
        // {
        //     printf("RW_1iteration:%f ", RW[i]); // 各要素を出力
        // }
        // for (size_t i = 0; i < 2 * K; i++)
        // {
        //     printf("dRdW_1iteration:%f ", dRdW[i]); // 各要素を出力
        // }
        // printf("\n");
        // }

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
                      alpha_S_R, alpha_S_W, alpha_L_R, alpha_L_W,
                      power_S_H_R, power_S_H_W, power_L_H_R, power_L_H_W,
                      beta_S_R, beta_S_W, beta_L_R, beta_L_W,
                      power_S_F_R, power_S_F_W, power_L_F_R, power_L_F_W,
                      RW_before, m_fixed, n_fixed, T_power, dR_dot, dW_dot);

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
            // printf("Adaptive_on\n");
            t_before = 1.0;
        }

        // Step 5: momentum項の計算
        double t = (1.0 + sqrt(1.0 + 4.0 * t_before * t_before)) / 2.0;

        // printf("t: %f\n", t);

        for (size_t i = 0; i < 2 * K; i++)
        {
            p_bar[i] = fmax(RW_proj, RW[i] + ((t_before - 1.0) / t) * (RW[i] - RW_before[i]));
        }

        // RW_before, t_before の更新
        for (size_t i = 0; i < 2 * K; i++)
        {
            RW_before[i] = RW[i];
            // printf("RW_before: %f\n", RW_before[i]);
        }
        for (size_t i = 0; i < 2 * K; i++)
        {
            p_bar_before[i] = p_bar[i];
        }
        t_before = t;
        L_before = L;

        // printf("t_before: %f\n", t_before);

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

#define Size 900

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
    const double alpha_1 = 0.4;
    const double alpha_2 = 0.4;
    const double beta_1 = 0.4;
    const double beta_2 = 0.4;
    const double par_L = 0.2;
    const double eta = 1.2;
    const double p_proj = 1e-3;
    const double RW_proj = 1e-3;
    const double err_short = 1e-5;
    const int short_itr = 10000;
    double R_hist[K];
    double W_hist[K];
    double m0[K];

    // 静的領域に確保
    static double Coordinate_Data[Size][2];
    static double distance_matrix[Size][Size];
    static double T[Size][Size];
    static double T_power[Size][Size];
    static double T_power_v[Size][Size];
    static double n0[Size][Size];

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

    // モデルパラメータ計算
    double power_S_H_R = (1 - alpha_2) / (1 - alpha_1 - alpha_2);
    double power_S_H_W = alpha_2 / (1 - alpha_1 - alpha_2);
    double power_L_H_R = alpha_1 / (1 - alpha_1 - alpha_2);
    double power_L_H_W = (1 - alpha_1) / (1 - alpha_1 - alpha_2);

    double power_S_F_R = (1 - beta_2) / (1 - beta_1 - beta_2);
    double power_S_F_W = beta_2 / (1 - beta_1 - beta_2);
    double power_L_F_R = beta_1 / (1 - beta_1 - beta_2);
    double power_L_F_W = (1 - beta_1) / (1 - beta_1 - beta_2);

    double alpha_S_R = pow(alpha_1, power_S_H_R);
    double alpha_S_W = pow(alpha_2, power_S_H_W);
    double alpha_L_R = pow(alpha_1, power_L_H_R);
    double alpha_L_W = pow(alpha_2, power_L_H_W);

    double beta_S_R = pow(beta_1, power_S_F_R);
    double beta_S_W = pow(beta_2, power_S_F_W);
    double beta_L_R = pow(beta_1, power_L_F_R);
    double beta_L_W = pow(beta_2, power_L_F_W);

    double pi_beta = 1 - beta_1 - beta_2;

    // T行列とD行列の計算
    for (size_t i = 0; i < K; i++)
    {
        for (size_t j = 0; j < K; j++)
        {
            T[i][j] = fmax(Scaling * t * alter_T_num, t * distance_matrix[i][j]);
        }
    }

    compute_T_power(K, T_power, T, alpha_1, alpha_2);

    for (size_t i = 0; i < K; i++)
    {
        for (size_t j = 0; j < K; j++)
        {
            T_power_v[i][j] = (1 - alpha_1 - alpha_2) * T_power[i][j];
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
        alpha_S_R, alpha_S_W, alpha_L_R, alpha_L_W,
        power_S_H_R, power_S_H_W, power_L_H_R, power_L_H_W,
        beta_S_R, beta_S_W, beta_L_R, beta_L_W,
        power_S_F_R, power_S_F_W, power_L_F_R, power_L_F_W,
        pi_beta,
        T_power, T_power_v,
        err_short, short_itr);

    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("CPU time used: %f seconds\n", cpu_time_used);

    return 0;
}