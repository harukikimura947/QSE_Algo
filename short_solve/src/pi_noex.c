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

    double pi_noex_result[K];
    for (size_t i = 0; i < K; i++)
    {
        pi_noex_result[i] = 0.0;
    }

    double R[K];
    double W[K];
    for (size_t i = 0; i < K; i++)
    {
        R[i] = 1.0;
        W[i] = 1.0;
    }

    pi_noex(K, pi_noex_result, R, W,
            beta_L_R, beta_S_W, power_L_F_R, power_S_F_W, pi_beta);

    for (size_t i = 0; i < K; ++i)
    {
        printf("pi_noex:%f ", pi_noex_result[i]);
    }

    return 0;
}