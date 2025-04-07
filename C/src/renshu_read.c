// #include <stdio.h>
// #include <stdlib.h>
// #include <string.h>
// #include <direct.h>

// // CSVファイルから整数を読み込む関数 (エラー処理付き)
// int read_int_from_csv(FILE *fp, int *value)
// {
//     if (fscanf(fp, "%d", value) != 1)
//     {
//         if (feof(fp))
//         {
//             fprintf(stderr, "Error: End of file reached unexpectedly.\n");
//         }
//         else
//         {
//             fprintf(stderr, "Error: Invalid integer format in CSV.\n");
//             // CSVファイル中に数字でないものがある場合など
//         }
//         // ファイルポインタを先頭に戻す．
//         fseek(fp, 0, SEEK_SET);
//         return 0; // 読み込み失敗
//     }
//     return 1; // 読み込み成功
// }

// int main()
// {
//     // 読み込む配列のサイズ (書き込み時と同じにする必要がある)
//     // サイズが不明な場合、動的に確保する必要がある（後述）
//     int arr1d_size = 5;
//     int arr1d[arr1d_size];

//     int rows = 3;
//     int cols = 3;
//     int arr2d[rows][cols]; // VLA (可変長配列) - C99以降

//     FILE *fp;
//     char filepath[256];

//     snprintf(filepath, sizeof(filepath), "C:\\Users\\%s\\Downloads\\aaa\\Array.csv", getenv("USERNAME"));

//     fp = fopen(filepath, "r"); // "r" モード (読み込み)
//     if (fp == NULL)
//     {
//         fprintf(stderr, "ファイルを開けませんでした: %s\n", filepath);
//         perror("fopen");
//         return 1;
//     }

//     // 1次元配列の読み込み
//     char line[1024]; // 各行を読み込むバッファ
//                      // "1D Array:" の行を読み飛ばす
//     if (!fgets(line, sizeof(line), fp))
//     {
//         perror("fgets");
//         return 1;
//     }
//     // 数字が始まる行まで読み飛ばす.
//     while (fgets(line, sizeof(line), fp) != NULL)
//     {
//         if (isdigit(line[0]) || (line[0] == '-' && isdigit(line[1])))
//         { // 数字か，-数字ではじまる場合
//             break;
//         }
//     }

//     // 1次元配列の読み込み
//     char *token = strtok(line, ","); // 最初のトークンを取得
//     for (int i = 0; i < arr1d_size; i++)
//     {
//         if (token == NULL)
//         {
//             fprintf(stderr, "Error: Not enough elements in 1D array data.\n");
//             fclose(fp);
//             return 1;
//         }
//         arr1d[i] = atoi(token);    // 文字列を整数に変換
//         token = strtok(NULL, ","); // 次のトークンを取得
//     }

//     // "2D Array:" の行を読み飛ばす
//     while (fgets(line, sizeof(line), fp) != NULL)
//     {
//         if (strstr(line, "2D Array:") != NULL) // 行に"2D Array:"が含まれているか確認
//         {
//             break;
//         }
//     }

//     // 2次元配列の読み込み
//     for (int i = 0; i < rows; i++)
//     {
//         if (fgets(line, sizeof(line), fp) == NULL)
//         {
//             fprintf(stderr, "Error reading 2D array row %d\n", i);
//             fclose(fp);
//             return 1;
//         }

//         char *token = strtok(line, ",");
//         for (int j = 0; j < cols; j++)
//         {
//             if (token == NULL)
//             {
//                 fprintf(stderr, "Error: Not enough elements in row %d\n", i);
//                 fclose(fp);
//                 return 1;
//             }
//             arr2d[i][j] = atoi(token); // 文字列を整数に変換
//             token = strtok(NULL, ","); // 次のトークンを取得. strtokは内部で状態を保持するので、NULLを指定する
//         }
//     }

//     fclose(fp);

//     // 読み込んだデータの確認 (表示)
//     printf("1D Array:\n");
//     for (int i = 0; i < arr1d_size; i++)
//     {
//         printf("%d ", arr1d[i]);
//     }
//     printf("\n");

//     printf("2D Array:\n");
//     for (int i = 0; i < rows; i++)
//     {
//         for (int j = 0; j < cols; j++)
//         {
//             printf("%d ", arr2d[i][j]);
//         }
//         printf("\n");
//     }

//     return 0;
// }

// int main()
// {
//     //--- 数値実験のパラメータと結果の配列 (例) ---
//     int param_values[] = {10, 20, 30}; // パラメータの値
//     int num_params = sizeof(param_values) / sizeof(param_values[0]);

//     //--- 書き込み先の基本ディレクトリ ---
//     const char *base_directory = "C:\\Users\\%s\\Downloads\\results";

//     for (int p = 0; p < num_params; p++)
//     {
//         int current_param = param_values[p];

//         //--- ディレクトリを作成 ---
//         char *username = getenv("USERNAME");
//         if (username == NULL)
//         {
//             fprintf(stderr, "ユーザー名を取得できませんでした。\n");
//             return 1;
//         }

//         // ディレクトリを作成 (例: results_10, results_20, results_30)
//         if (create_directory(base_directory, current_param) != 0)
//         {
//             fprintf(stderr, "ディレクトリ作成に失敗しました。\n");
//             return 1;
//         }

//         //--- ここで数値実験を行い、結果を配列に格納する ---
//         // (例: 1次元配列)
//         int result_arr1d[5];
//         for (int i = 0; i < 5; i++)
//         {
//             result_arr1d[i] = current_param * (i + 1); // パラメータを使った簡単な計算
//         }

//         // (例: 2次元配列)
//         int result_arr2d[2][3];
//         for (int i = 0; i < 2; i++)
//         {
//             for (int j = 0; j < 3; j++)
//             {
//                 result_arr2d[i][j] = current_param * (i + 1) * (j + 1);
//             }
//         }

//         //--- 配列をCSVファイルに書き込む ---
//         char base_filename[256];
//         snprintf(base_filename, sizeof(base_filename), "C:\\Users\\%s\\Downloads\\results\\%s_%d\\result", getenv("USERNAME"), getenv("USERNAME"), current_param); //  ファイルパス+ファイル名を指定

//         if (write_array1d_to_csv(base_filename, current_param, result_arr1d, 5) != 0)
//         {
//             fprintf(stderr, "1次元配列の書き込みに失敗しました。\n");
//             return 1;
//         }
//         snprintf(base_filename, sizeof(base_filename), "C:\\Users\\%s\\Downloads\\results\\%s_%d\\result2d", getenv("USERNAME"), getenv("USERNAME"), current_param);

//         if (write_array2d_to_csv(base_filename, current_param, result_arr2d, 2, 3) != 0)
//         {
//             fprintf(stderr, "2次元配列の書き込みに失敗しました。\n");
//             return 1;
//         }
//         printf("パラメータ %d の結果を保存しました。\n", current_param);
//     }

//     printf("すべてのパラメータの結果を保存しました。\n");

//     return 0;
// }

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_ELEMENTS 1800 // 期待する配列の最大要素数
#define LINE_LENGTH 65536 // 1行分の最大長（大きめにする）
#define FILENAME_TEMPLATE "C://Users//kimura//Downloads//short_solve//data//K=900//par_L=%.4f,eta=%.4f//20000iteration.csv"

// RW配列を読み取る関数（実際に読み取った要素数を返す）
int read_RW(double A, double B, double RW[])
{
    char filename[256];
    snprintf(filename, sizeof(filename), FILENAME_TEMPLATE, A, B);

    FILE *file = fopen(filename, "r");
    if (!file)
    {
        printf("ファイル %s を開けませんでした。\n", filename);
        return 0;
    }

    char line[LINE_LENGTH];

    // 1行目（ヘッダー）を読み飛ばす
    if (fgets(line, sizeof(line), file) == NULL)
    {
        printf("ファイル %s のヘッダーを読み取れませんでした。\n", filename);
        fclose(file);
        return 0;
    }

    // 2行目（データ）を読み取る
    if (fgets(line, sizeof(line), file) == NULL)
    {
        printf("ファイル %s のデータ行を読み取れませんでした。\n", filename);
        fclose(file);
        return 0;
    }

    fclose(file); // ファイルを閉じる（メモリ節約）

    // **修正: `sscanf()` を使って安全にデータを解析**
    int i = 0;
    char *ptr = line;
    while (*ptr && i < MAX_ELEMENTS)
    {
        if (sscanf(ptr, "%lf", &RW[i]) == 1)
        { // 1つのdouble値を読み取る
            i++;
        }
        while (*ptr && *ptr != ',')
            ptr++; // カンマまで進める
        if (*ptr == ',')
            ptr++; // 次の値へ
    }

    return i; // 実際に読み取った要素数を返す
}

int main()
{
    double RW[MAX_ELEMENTS]; // RW配列

    double par_L[] = {0.0050};
    double eta[] = {1.8000};

    int num_L = sizeof(par_L) / sizeof(par_L[0]);
    int num_eta = sizeof(eta) / sizeof(eta[0]);

    for (int p1 = 0; p1 < num_L; p1++)
    {
        for (int p2 = 0; p2 < num_eta; p2++)
        {
            printf("Reading file for par_L=%.4f, eta=%.4f...\n", par_L[p1], eta[p2]);

            int valid_elements = read_RW(par_L[p1], eta[p2], RW); // 実際の要素数を取得

            // 読み取った要素数を表示
            printf("読み取った要素数: %d\n", valid_elements);

            // RW 配列のデータを表示（valid_elements のみ）
            // for (int i = 0; i < valid_elements; i++)
            // {
            //     printf("%lf ", RW[i]);
            // }
            // printf("\n");
        }
    }

    return 0;
}
