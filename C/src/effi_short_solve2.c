#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <direct.h>
#include <errno.h>

#define LINE_LENGTH 1000000 // 1行分の最大長（大きめにする）
#define FILENAME_TEMPLATE "C://Users//kimura//Downloads//short_solve//data//K=%d//par_L=%.4f,eta=%.4f//10000iteration.csv"

// RW配列とZを読み取る関数（実際に読み取った要素数を返す）
int read_RW(size_t K, double A, double B, double **RW, double *Z)
{
    // MAX_ELEMENTS を 2 * K として動的に決定
    size_t MAX_ELEMENTS = 2 * K;

    // RW配列を動的に確保
    *RW = (double *)malloc(MAX_ELEMENTS * sizeof(double));
    if (*RW == NULL)
    {
        printf("メモリの確保に失敗しました。\n");
        return 0;
    }

    char filename[256];
    snprintf(filename, sizeof(filename), FILENAME_TEMPLATE, K, A, B);

    FILE *file = fopen(filename, "r");
    if (!file)
    {
        printf("ファイル %s を開けませんでした。\n", filename);
        free(*RW); // メモリを解放
        return 0;
    }

    char *line = (char *)malloc(LINE_LENGTH * sizeof(char));
    if (line == NULL)
    {
        printf("メモリの確保に失敗しました。\n");
        fclose(file);
        free(*RW); // メモリを解放
        return 0;
    }

    // 1行目（ヘッダー）を読み飛ばす
    if (fgets(line, LINE_LENGTH, file) == NULL)
    {
        printf("ファイル %s のヘッダーを読み取れませんでした。\n", filename);
        fclose(file);
        free(*RW);  // メモリを解放
        free(line); // lineを解放
        return 0;
    }

    // データ行を読み取る
    int i = 0;
    while (fgets(line, LINE_LENGTH, file) != NULL && i < MAX_ELEMENTS)
    {
        double rw_value, z_value;
        if (i == 0)
        {
            // 最初の行では RW と Z の両方を読み取る
            if (sscanf(line, "%lf,%lf", &rw_value, &z_value) == 2)
            {
                (*RW)[i] = rw_value;
                *Z = z_value; // Zの値を記録
                i++;
            }
        }
        else
        {
            // 2行目以降では RW のみを読み取る
            if (sscanf(line, "%lf", &rw_value) == 1)
            {
                (*RW)[i] = rw_value;
                i++;
            }
        }
    }

    fclose(file);
    free(line); // lineを解放
    return i;   // 実際に読み取った要素数を返す
}

int create_parameter_directory(const char *base_dir, int K, double par_L, double eta, char *dirpath_out)
{
    char *username = getenv("USERNAME");
    if (username == NULL)
    {
        fprintf(stderr, "ユーザー名を取得できませんでした。\n");
        return -1;
    }

    // 親ディレクトリ (Downloads) のパスを生成 (Windows用)
    char parent_dir[256];
    snprintf(parent_dir, sizeof(parent_dir), "C:\\Users\\%s\\Downloads", username); // ここを修正

    // 親ディレクトリが存在しない場合は作成
    if (_mkdir(parent_dir) == -1)
    {
        if (errno != EEXIST)
        {
            perror("_mkdir (parent)");
            return -1;
        }
    }

    // short_solve ディレクトリのパスを生成 (Windows用)
    char short_dir[256];
    snprintf(short_dir, sizeof(short_dir), "%s\\short_solve", parent_dir);

    // data ディレクトリのパスを生成 (Windows用)
    char data_dir[256];
    snprintf(data_dir, sizeof(data_dir), "%s\\data", short_dir);

    // Kディレクトリのパスを生成 (Windows用)
    char k_dir[256];
    snprintf(k_dir, sizeof(k_dir), "%s\\K=%d", data_dir, K);

    // Kディレクトリが存在しない場合は作成
    if (_mkdir(k_dir) == -1)
    {
        if (errno != EEXIST)
        {
            perror("_mkdir (K dir)");
            return -1;
        }
    }

    // パラメータ付きディレクトリのパスを生成 (Windows用, par_L と eta)
    char dirpath[256];
    snprintf(dirpath, sizeof(dirpath), "%s\\par_L=%.4f,eta=%.4f", k_dir, par_L, eta);

    // パラメータ付きディレクトリを作成
    if (_mkdir(dirpath) == -1)
    {
        if (errno != EEXIST)
        {
            perror("_mkdir");
            return -1;
        }
    }

    // 作成したディレクトリのパスを返す
    if (dirpath_out != NULL)
    {
        strncpy(dirpath_out, dirpath, 256);
        dirpath_out[255] = '\0';
    }

    return 0;
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
    // 動的配列の確保
    double *R = (double *)malloc(K * sizeof(double));
    double *W = (double *)malloc(K * sizeof(double));
    double *R_inv = (double *)malloc(K * sizeof(double));
    double *W_inv = (double *)malloc(K * sizeof(double));
    double *R_inv_alpha = (double *)malloc(K * sizeof(double));
    double *W_inv_alpha = (double *)malloc(K * sizeof(double));
    double *R_inv_beta = (double *)malloc(K * sizeof(double));
    double *W_inv_beta = (double *)malloc(K * sizeof(double));
    double *R_T = (double *)malloc(K * sizeof(double));
    double *W_T = (double *)malloc(K * sizeof(double));

    // メモリ確保のエラーチェック
    if (R == NULL || W == NULL || R_inv == NULL || W_inv == NULL ||
        R_inv_alpha == NULL || W_inv_alpha == NULL || R_inv_beta == NULL ||
        W_inv_beta == NULL || R_T == NULL || W_T == NULL)
    {
        fprintf(stderr, "メモリの確保に失敗しました\n");
        exit(EXIT_FAILURE);
    }

    // R と W を RW から初期化
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
        dR[i] = S_bar - (coef_R_alpha * R_inv[i] * R_inv_alpha[i] * W_T[i]) - (coef_R_beta * R_inv[i] * R_inv_beta[i] * W_inv_beta[i] * m[i]);
    }

    // Compute gradient dZ/dW
    for (size_t j = 0; j < K; j++)
    {
        dW[j] = nE[j] - (coef_W_alpha * W_inv[j] * W_inv_alpha[j] * R_T[j]) - (coef_W_beta * W_inv[j] * W_inv_beta[j] * R_inv_beta[j] * m[j]);
    }

    // 動的配列の解放
    free(R);
    free(W);
    free(R_inv);
    free(W_inv);
    free(R_inv_alpha);
    free(W_inv_alpha);
    free(R_inv_beta);
    free(W_inv_beta);
    free(R_T);
    free(W_T);
}

void equilibrium(
    size_t K, double S_bar,
    double alpha_P1, double alpha_P2, double beta_P1, double beta_P2,
    double coef_R_alpha, double coef_W_alpha, double coef_R_beta, double coef_W_beta,
    double RW[],
    double nE[], double m[],
    double **T_n)
{
    // 動的配列の確保
    double *dR_equ = (double *)malloc(K * sizeof(double));
    double *dW_equ = (double *)malloc(K * sizeof(double));

    // メモリ確保のエラーチェック
    if (dR_equ == NULL || dW_equ == NULL)
    {
        fprintf(stderr, "メモリの確保に失敗しました\n");
        exit(EXIT_FAILURE);
    }

    short_dual_df(
        K, S_bar,
        alpha_P1, alpha_P2, beta_P1, beta_P2,
        coef_R_alpha, coef_W_alpha, coef_R_beta, coef_W_beta,
        RW,
        nE, m,
        T_n,
        dR_equ, dW_equ);

    // 最大値を求める
    double dR_max = fabs(dR_equ[0]); // 最初の要素を仮の最大値とする
    for (size_t i = 1; i < K; i++)
    {
        if (fabs(dR_equ[i]) > dR_max)
        {
            dR_max = fabs(dR_equ[i]);
        }
    }
    printf("dRの最大値: %.16f\n", dR_max);

    // 平均値を求める
    double dR_sum = 0.0;
    for (size_t i = 0; i < K; i++)
    {
        dR_sum += dR_equ[i];
    }
    double dR_average = dR_sum / K;
    printf("dRの平均値: %.16f\n", dR_average);

    // 最大値を求める
    double dW_max = fabs(dW_equ[0]); // 最初の要素を仮の最大値とする
    for (size_t i = 1; i < K; i++)
    {
        if (fabs(dW_equ[i]) > dW_max)
        {
            dW_max = fabs(dW_equ[i]);
        }
    }
    printf("dWの最大値: %.16f\n", dW_max);

    // 平均値を求める
    double dW_sum = 0.0;
    for (size_t i = 0; i < K; i++)
    {
        dW_sum += dW_equ[i];
    }
    double dW_average = dW_sum / K;
    printf("dWの平均値: %.16f\n", dW_average);

    // 動的配列の解放
    free(dR_equ);
    free(dW_equ);
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
    // 動的配列の確保
    double *R = (double *)malloc(K * sizeof(double));
    double *W = (double *)malloc(K * sizeof(double));
    double *R_inv = (double *)malloc(K * sizeof(double));
    double *W_inv = (double *)malloc(K * sizeof(double));
    double *R_inv_alpha = (double *)malloc(K * sizeof(double));
    double *R_inv_beta = (double *)malloc(K * sizeof(double));
    double *W_inv_alpha = (double *)malloc(K * sizeof(double));
    double *W_inv_beta = (double *)malloc(K * sizeof(double));
    double *R_T = (double *)malloc(K * sizeof(double));
    double *W_T = (double *)malloc(K * sizeof(double));

    // メモリ確保のエラーチェック
    if (R == NULL || W == NULL || R_inv == NULL || W_inv == NULL ||
        R_inv_alpha == NULL || R_inv_beta == NULL || W_inv_alpha == NULL ||
        W_inv_beta == NULL || R_T == NULL || W_T == NULL)
    {
        fprintf(stderr, "メモリの確保に失敗しました\n");
        exit(EXIT_FAILURE);
    }

    // R と W を RW から初期化
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
        dR[i] = S_bar - (coef_R_alpha * R_inv[i] * R_inv_alpha[i] * W_T[i]) - (coef_R_beta * R_inv[i] * R_inv_beta[i] * W_inv_beta[i] * m[i]);
    }

    // Compute gradient dZ/dW
    for (size_t j = 0; j < K; j++)
    {
        dW[j] = nE[j] - (coef_W_alpha * W_inv[j] * W_inv_alpha[j] * R_T[j]) - (coef_W_beta * W_inv[j] * W_inv_beta[j] * R_inv_beta[j] * m[j]);
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

    // 動的配列の解放
    free(R);
    free(W);
    free(R_inv);
    free(W_inv);
    free(R_inv_alpha);
    free(R_inv_beta);
    free(W_inv_alpha);
    free(W_inv_beta);
    free(R_T);
    free(W_T);

    return Z_SD;
}

double Z_SD(
    size_t K, double S_bar,
    double coef_pi, double coef_v,
    double alpha_P1, double alpha_P2, double beta_P1, double beta_P2,
    double RW[],
    double nE[], double m[],
    double **T_n,
    double dR[], double dW[])
{
    // 動的配列の確保
    double *R = (double *)malloc(K * sizeof(double));
    double *W = (double *)malloc(K * sizeof(double));
    double *R_inv = (double *)malloc(K * sizeof(double));
    double *W_inv = (double *)malloc(K * sizeof(double));
    double *R_inv_alpha = (double *)malloc(K * sizeof(double));
    double *W_inv_alpha = (double *)malloc(K * sizeof(double));
    double *R_inv_beta = (double *)malloc(K * sizeof(double));
    double *W_inv_beta = (double *)malloc(K * sizeof(double));
    double *W_T = (double *)malloc(K * sizeof(double));

    // メモリ確保のエラーチェック
    if (R == NULL || W == NULL || R_inv == NULL || W_inv == NULL ||
        R_inv_alpha == NULL || W_inv_alpha == NULL || R_inv_beta == NULL ||
        W_inv_beta == NULL || W_T == NULL)
    {
        fprintf(stderr, "メモリの確保に失敗しました\n");
        exit(EXIT_FAILURE);
    }

    // R と W を RW から初期化
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

    // Compute W_T
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

    // 動的配列の解放
    free(R);
    free(W);
    free(R_inv);
    free(W_inv);
    free(R_inv_alpha);
    free(W_inv_alpha);
    free(R_inv_beta);
    free(W_inv_beta);
    free(W_T);

    return Z_SD;
}

double backtracking(
    size_t K, double S_bar,
    double alpha_P1, double alpha_P2, double beta_P1, double beta_P2,
    double coef_R_alpha, double coef_W_alpha, double coef_R_beta, double coef_W_beta,
    double coef_pi, double coef_v,
    double nE[], double m[],
    double **T_n,
    double Z_SD_p_bar, double p_proj, double eta,
    double L, double dR[], double dW[], double p_bar[])
{
    double L_bar = L;

    // 動的配列の確保
    double *dRdW = (double *)malloc(2 * K * sizeof(double));
    double *p = (double *)malloc(2 * K * sizeof(double));
    double *diff_p = (double *)malloc(2 * K * sizeof(double));

    if (dRdW == NULL || p == NULL || diff_p == NULL)
    {
        fprintf(stderr, "メモリの確保に失敗しました\n");
        exit(EXIT_FAILURE);
    }

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
        double Z_SD_p = Z_SD(K, S_bar, coef_pi, coef_v,
                             alpha_P1, alpha_P2, beta_P1, beta_P2,
                             p, nE, m, T_n, dR, dW);

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

    // 動的配列の解放
    free(dRdW);
    free(p);
    free(diff_p);

    return L_bar;
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
    double *RW_before = (double *)malloc(2 * K * sizeof(double));
    double *RW = (double *)malloc(2 * K * sizeof(double));
    double *data_R = (double *)malloc(K * sizeof(double));
    double *data_W = (double *)malloc(K * sizeof(double));
    double *R = (double *)malloc(K * sizeof(double));
    double *W = (double *)malloc(K * sizeof(double));
    double *data_RW = NULL;
    double data_Z = 0;
    double *p_bar_before = (double *)malloc(2 * K * sizeof(double));
    double *p_bar = (double *)malloc(2 * K * sizeof(double));
    double L_before = par_L;
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

    // 要素数取得&データ読み込み
    int valid_elements = read_RW(K, par_L, eta, &data_RW, &data_Z);

    for (size_t i = 0; i < K; ++i)
    {
        data_R[i] = data_RW[i];
        data_W[i] = data_RW[K + i];
    }

    double shiken = 0.0;
    for (size_t i = 0; i < 2 * K; i++)
    {
        if (data_RW[i] > shiken)
        {
            shiken = data_RW[i];
        }
    }

    printf("shiken: %.10f\n", shiken);
    printf("data_Z: %.10f\n", data_Z);

    // ディレクトリ作成
    char dirpath[256];
    if (create_parameter_directory("C:\\Users", K, par_L, eta, dirpath) != 0)
    {
        fprintf(stderr, "Error creating directory\n");
        exit(EXIT_FAILURE);
    }

    // CSVファイルパス
    char csv_filepath[512];
    snprintf(csv_filepath, sizeof(csv_filepath), "%s\\iteration_data.csv", dirpath);

    // CSVファイルパス
    // char csv_filepath[512];
    // snprintf(csv_filepath, sizeof(csv_filepath), "%s\\10000iteration.csv", dirpath);

    // CSVファイルを開く
    FILE *fp = fopen(csv_filepath, "w"); // "a" は追記モード, "w"は上書きモード
    if (fp == NULL)
    {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    // ★★★ ヘッダーの書き込みは、ファイルが新規作成された場合のみ ★★★
    // ファイルポインタの位置が0 (ファイルの先頭) なら、ファイルは空
    // if (ftell(fp) == 0)
    // { 
    //     fprintf(fp, "RW,Z\n");
    // }

    if (ftell(fp) == 0)
    {
        fprintf(fp,"iteration,max_diff,R_max_diff,W_max_diff,Z_diff\n");
    }

    for (int k = 0; k < short_itr; k++)
    {
        for (size_t i = 0; i < K; i++)
        {
            dR[i] = 0.0;
            dW[i] = 0.0;
        }

        // Step 1: Backtracking
        Z_SD_p_bar = effi_grad_Z(K, S_bar,
                                 alpha_P1, alpha_P2, beta_P1, beta_P2,
                                 coef_R_alpha, coef_W_alpha, coef_R_beta, coef_W_beta,
                                 coef_pi, coef_v,
                                 p_bar_before, nE, m, T_n, dR, dW);

        L = backtracking(K, S_bar,
                         alpha_P1, alpha_P2, beta_P1, beta_P2,
                         coef_R_alpha, coef_W_alpha, coef_R_beta, coef_W_beta,
                         coef_pi, coef_v,
                         nE, m, T_n,
                         Z_SD_p_bar, p_proj, eta, L_before, dR, dW, p_bar_before);

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
        double R_max_diff = 0.0;
        double W_max_diff = 0.0;

        double Z_now = Z_SD(K, S_bar, coef_pi, coef_v,
                        alpha_P1, alpha_P2, beta_P1, beta_P2,
                        RW, nE, m, T_n, dR, dW);

        for (size_t i = 0; i < K; ++i)
        {
            R[i] = RW[i];
            W[i] = RW[K + i];
        }

        for (size_t i = 0; i < 2 * K; i++)
        {
            double diff = fabs((data_RW[i] - RW[i]) / RW[i]);
            if (diff > max_diff)
            {
                max_diff = diff;
            }
        }

        for (size_t i = 0; i < K; i++)
        {
            double R_diff = fabs((data_R[i] - R[i]) / R[i]);
            if (R_diff > R_max_diff)
            {
                R_max_diff = R_diff;
            }
        }

        for (size_t i = 0; i < K; i++)
        {
            double W_diff = fabs((data_W[i] - W[i]) / W[i]);
            if (W_diff > W_max_diff)
            {
                W_max_diff = W_diff;
            }
        }

        double Z_diff = fabs((data_Z - Z_now) / Z_now);

        fprintf(fp, "%d,%.16f,%.16f,%.16f,%.16f\n",
            k + 1, max_diff, R_max_diff, W_max_diff, Z_diff);

        // Step 4: Adaptive restart
        short_dual_df(
            K, S_bar,
            alpha_P1, alpha_P2, beta_P1, beta_P2,
            coef_R_alpha, coef_W_alpha, coef_R_beta, coef_W_beta,
            RW_before, nE, m, T_n, dR_dot, dW_dot);

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

    // double Z_result = Z_SD(K, S_bar, coef_pi, coef_v,
    //                        alpha_P1, alpha_P2, beta_P1, beta_P2,
    //                        RW, nE, m, T_n, dR, dW);

    // // RWの値を縦に記録
    // for (size_t i = 0; i < 2 * K; i++)
    // {
    //     fprintf(fp, "%.16f", RW[i]); // RWの値を記録

    //     // 最初の行のみ Z を記録
    //     if (i == 0)
    //     {
    //         fprintf(fp, ",%.16f", Z_result); // Z を記録
    //     }

    //     fprintf(fp, "\n"); // 改行
    // }

    fclose(fp);

    *g_out = g;

    equilibrium(
        K, S_bar,
        alpha_P1, alpha_P2, beta_P1, beta_P2,
        coef_R_alpha, coef_W_alpha, coef_R_beta, coef_W_beta,
        RW,
        nE, m,
        T_n);

    printf("g: %d\n", g);

    // メモリの解放
    free(RW_before);
    free(RW);
    free(data_RW);
    free(p_bar_before);
    free(p_bar);
    free(dR);
    free(dW);
    free(dRdW);
    free(dR_dot);
    free(dW_dot);
    free(dRdW_dot);
    free(RW_diff);
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
    int K_list[] = {2500};
    double par_L[] = {0.0050};
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
        const int int_Col = Col;
        const int int_K = int_Col * int_Col;
        const double Scaling = 10.0 / int_Col;
        const double S_bar = S_total / int_K;

        // 動的配列の確保
        double *R_hist = (double *)malloc(K * sizeof(double));
        double *W_hist = (double *)malloc(K * sizeof(double));
        double *RW_hist = (double *)malloc(2 * K * sizeof(double));
        double *data_RW = (double *)malloc(2 * K * sizeof(double));
        double *m0 = (double *)malloc(K * sizeof(double));

        // メモリ確保のエラーチェック
        if (R_hist == NULL || W_hist == NULL || RW_hist == NULL || data_RW == NULL || m0 == NULL)
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

                short_solve(
                    K, S_bar,
                    alpha_P1, alpha_P2, beta_P1, beta_P2,
                    coef_R_alpha, coef_W_alpha, coef_R_beta, coef_W_beta,
                    coef_pi, coef_v,
                    nE, m0, T_n,
                    RW_proj, p_proj, par_L[p1], eta[p2],
                    R_hist, W_hist, g_out_ptr,
                    err_short, short_itr);

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
        free(data_RW);
        free(m0);
        free(nE);
    }

    return 0;
}