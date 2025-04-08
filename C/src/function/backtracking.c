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

    for (size_t i = 0; i < 2 * K; i++)
    {
        printf("RW:%f ", RW[i]);
    }

    // RWからRとWに分離
    for (size_t i = 0; i < K; i++)
    {
        R[i] = RW[i];
        W[i] = RW[K + i];
        printf("W:%f ", W[i]);
    }

    // R_vとW_vの計算
    for (size_t i = 0; i < K; i++)
    {
        R_v[i] = alpha_L_R * pow(R[i], -power_L_H_R);
        W_v[i] = alpha_S_W * pow(W[i], -power_S_H_W);
        printf("R_v:%f ", R_v[i]);
        printf("W_v:%f ", W_v[i]);
    }

    // vnの計算
    for (size_t i = 0; i < K; i++)
    {
        for (size_t j = 0; j < K; j++)
        {
            vn += (E * W[j] + T_power_v[i][j] * R_v[i] * W_v[j]) * n[i][j];
        }
    }

    printf("vn:%f ", vn);

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

    // 初期 Z_SD(p_bar)
    Z_SD_p_bar = Z_SD(K, E, S_bar, p_bar, m_fixed, n_fixed, T_power_v,
                      beta_L_R, beta_S_W, power_L_F_R, power_S_F_W, pi_beta,
                      alpha_L_R, alpha_S_W, power_L_H_R, power_S_H_W);

    for (size_t k = 0; k < 10000; k++)
    {
        for (size_t j = 0; j < 2 * K; j++)
        {
            p[j] = fmax(p_proj, p_bar[j] - dRdW[j] / L_bar);
        }

        // 判定条件
        double Z_SD_p = Z_SD(K, E, S_bar, p, m_fixed, n_fixed, T_power_v,
                             beta_L_R, beta_S_W, power_L_F_R, power_S_F_W, pi_beta,
                             alpha_L_R, alpha_S_W, power_L_H_R, power_S_H_W);
        double dot = dot_product(2 * K, p, dRdW);
        double norm_squared = 0.0;
        for (size_t j = 0; j < 2 * K; j++)
        {
            double diff = p[j] - p_bar[j];
            norm_squared += diff * diff;
        }

        if (Z_SD_p - (Z_SD_p_bar + dot + 0.5 * L_bar * norm_squared) <= 0.0)
        {
            break;
        }

        // 更新
        L_bar *= eta; // algprm.eta として最後の要素を使用
    }

    return L_bar;
}

// メイン関数
int main()
{
    // 定数パラメータの定義
    const int E = 5;
    const size_t Col = 2;
    const size_t K = Col * Col;
    const int int_Col = 2;
    const int int_K = int_Col * int_Col;
    const double M = 1.0;
    const double N = 1.0;
    const double Scaling = 10.0 / int_Col;
    const double alter_T_num = 0.5;
    const double S_total = 100;
    const double S_bar = S_total / int_K;
    const double t = 0.1;
    // const double tau = 0.5;
    // const double theta_firm = 1.0;
    // const double theta_house = 1.0;
    const double alpha_1 = 0.4;
    const double alpha_2 = 0.4;
    const double beta_1 = 0.4;
    const double beta_2 = 0.4;
    const double par_L = 0.2;
    const double eta = 1.2;
    const double p_proj = 1e-3;
    const double RW_proj = 1e-3;
    const double err_short = 1e-5;
    // const double err_long = 1e-3;
    // const double obj_corr = 1.0;
    const int short_itr = 1000;
    // const int long_itr = 1000;

    // 座標データの生成
    double Coordinate_Data[K][2];
    for (int i = 0; i < int_K; i++)
    {
        Coordinate_Data[i][0] = (i % int_Col) * Scaling;
        Coordinate_Data[i][1] = (i / int_Col) * Scaling;
    }

    // 距離行列の作成
    double distance_matrix[K][K];
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
    double T[K][K];
    // double D[K][K];
    for (size_t i = 0; i < K; i++)
    {
        for (size_t j = 0; j < K; j++)
        {
            T[i][j] = fmax(Scaling * t * alter_T_num, t * distance_matrix[i][j]);
            // D[i][j] = exp(-tau * distance_matrix[i][j]);
        }
    }

    double T_power[K][K];
    compute_T_power(K, T_power, T, alpha_1, alpha_2);

    double T_power_v[K][K];
    for (size_t i = 0; i < K; i++)
    {
        for (size_t j = 0; j < K; j++)
        {
            T_power_v[i][j] = (1 - alpha_1 - alpha_2) * T_power[i][j];
        }
    }

    // 初期化
    double m_per = M / int_K;
    double m0[K];
    for (size_t i = 0; i < K; i++)
    {
        m0[i] = m_per;
    }

    double n0[K][K];
    for (size_t i = 0; i < K; i++)
    {
        for (size_t j = 0; j < K; j++)
        {
            n0[i][j] = N / (int_K * int_K);
        }
    }

    double RW[2 * K];
    for (size_t i = 0; i < 2 * K; i++)
    {
        RW[i] = 1.0;
    }

    double dR[K];
    double dW[K];
    for (size_t i = 0; i < K; i++)
    {
        dR[i] = 0;
        dW[i] = 0;
    }

    double L_bar = backtracking(K, S_bar, E,
                 alpha_S_R, alpha_S_W, alpha_L_R, alpha_L_W,
                 power_S_H_R, power_S_H_W, power_L_H_R, power_L_H_W,
                 beta_S_R, beta_S_W, beta_L_R, beta_L_W,
                 power_S_F_R, power_S_F_W, power_L_F_R, power_L_F_W, pi_beta, T_power, T_power_v,
                 p_proj, eta, par_L, dR, dW, RW, m0, n0);

    printf("L_bar:%f ", L_bar);

    return 0;
}