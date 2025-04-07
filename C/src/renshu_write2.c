#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <direct.h>
#include <errno.h>

// 1次元配列をCSVファイルに書き込む関数 (double型対応)
int write_array1d_to_csv(const char *dirpath, const char *filename, const double *arr, int size)
{
    char filepath[256];
    FILE *fp;

    // ファイルパスを生成
    snprintf(filepath, sizeof(filepath), "%s\\%s", dirpath, filename);

    fp = fopen(filepath, "w");
    if (fp == NULL)
    {
        perror("fopen");
        return -1; // エラー
    }

    // 1次元配列を書き込む
    for (int i = 0; i < size; i++)
    {
        fprintf(fp, "%.1f", arr[i]); // double型を出力 (精度を指定)
        if (i < size - 1)
        {
            fprintf(fp, ",");
        }
    }
    fprintf(fp, "\n");

    fclose(fp);
    return 0; // 成功
}

// ディレクトリを作成する関数 (Windows専用, パラメータ付き, 親ディレクトリも作成, double対応)
int create_parameter_directory(const char *base_dir, double param1, double param2)
{
    char dirpath[256];
    char parent_dir[256];
    char *username = getenv("USERNAME");
    if (username == NULL)
    {
        fprintf(stderr, "ユーザー名を取得できませんでした。\n");
        return -1;
    }

    // 親ディレクトリ (results) のパスを生成
    snprintf(parent_dir, sizeof(parent_dir), "%s\\%s\\Downloads\\results", base_dir, username);

    // 親ディレクトリが存在しない場合は作成
    if (_mkdir(parent_dir) == -1)
    {
        if (errno != EEXIST)
        {
            perror("_mkdir (parent)");
            return -1;
        }
    }

    // パラメータ付きディレクトリのパスを生成 (double型に対応)
    snprintf(dirpath, sizeof(dirpath), "%s\\para1=%.1f,para2=%.1f", parent_dir, param1, param2);

    // パラメータ付きディレクトリを作成
    if (_mkdir(dirpath) == -1)
    {
        if (errno != EEXIST)
        {
            perror("_mkdir");
            return -1;
        }
    }
    return 0;
}

int main()
{
    //--- 数値実験のパラメータと結果の配列 (例) ---
    double param1_values[] = {10.0, 20.5};    // パラメータ1の値
    double param2_values[] = {1.2, 2.3, 3.4}; // パラメータ2の値
    int num_param1 = sizeof(param1_values) / sizeof(param1_values[0]);
    int num_param2 = sizeof(param2_values) / sizeof(param2_values[0]);

    //--- 書き込み先の基本ディレクトリ ---
    const char *base_directory = "C:\\Users";

    for (int p1 = 0; p1 < num_param1; p1++)
    {
        for (int p2 = 0; p2 < num_param2; p2++)
        {
            double current_param1 = param1_values[p1];
            double current_param2 = param2_values[p2];

            //--- パラメータディレクトリを作成 ---
            if (create_parameter_directory(base_directory, current_param1, current_param2) != 0)
            {
                fprintf(stderr, "ディレクトリ作成に失敗しました。\n");
                return 1;
            }

            //--- ディレクトリパスを生成 ---
            char *username = getenv("USERNAME");
            if (username == NULL)
            {
                fprintf(stderr, "ユーザー名を取得できませんでした。\n");
                return 1;
            }
            char dirpath[256];
            snprintf(dirpath, sizeof(dirpath), "C:\\Users\\%s\\Downloads\\results\\para1=%.1f,para2=%.1f", username, current_param1, current_param2);

            //--- ここで数値実験を行い、結果を配列に格納する ---
            // (例: 1次元配列)
            int size1d = 5;
            double result_arr1d[size1d]; // double型に変更, サイズは固定
            for (int i = 0; i < size1d; i++)
            {
                result_arr1d[i] = current_param1 * (i + 1) + current_param2; // パラメータを使った計算例
            }

            //--- 配列をCSVファイルに書き込む ---
            // ファイル名は固定
            if (write_array1d_to_csv(dirpath, "result1d.csv", result_arr1d, size1d) != 0)
            {
                fprintf(stderr, "1次元配列の書き込みに失敗しました.\n");
                return 1;
            }
            printf("パラメータ (para1=%.1f, para2=%.1f) の結果を保存しました。\n", current_param1, current_param2);
        }
    }

    printf("すべてのパラメータの結果を保存しました。\n");
    return 0;
}