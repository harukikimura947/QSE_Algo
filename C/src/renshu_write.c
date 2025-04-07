#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <direct.h> // _mkdir()用

int main()
{
    // 1次元配列
    int arr1d[] = {1, 2, 3, 4, 5};
    int arr1d_size = sizeof(arr1d) / sizeof(arr1d[0]);

    // 2次元配列
    int arr2d[][3] = {
        {1, 2, 3},
        {4, 5, 6},
        {7, 8, 9}};
    int rows = sizeof(arr2d) / sizeof(arr2d[0]);
    int cols = sizeof(arr2d[0]) / sizeof(arr2d[0][0]);

    FILE *fp;
    char filepath[256];

    snprintf(filepath, sizeof(filepath), "C:\\Users\\%s\\Downloads\\aaa\\Array.csv", getenv("USERNAME"));

    char dirpath[256];
    snprintf(dirpath, sizeof(dirpath), "C:\\Users\\%s\\Downloads\\aaa", getenv("USERNAME"));
    if (_mkdir(dirpath) == -1)
    {
        if (errno != EEXIST)
        {
            fprintf(stderr, "ディレクトリを作成できませんでした。: %s\n", dirpath);
            // return 1; // エラーが起きても続行する場合
        }
    }

    fp = fopen(filepath, "w"); // "w" モード (上書き)
    if (fp == NULL)
    {
        fprintf(stderr, "ファイルを開けませんでした: %s\n", filepath);
        perror("fopen");
        return 1;
    }

    // 1次元配列を書き込む
    fprintf(fp, "1D Array:\n");
    for (int i = 0; i < arr1d_size; i++)
    {
        fprintf(fp, "%d", arr1d[i]);
        if (i < arr1d_size - 1)
        {
            fprintf(fp, ",");
        }
    }
    fprintf(fp, "\n");

    // 2次元配列を書き込む
    fprintf(fp, "2D Array:\n");
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            fprintf(fp, "%d", arr2d[i][j]);
            if (j < cols - 1)
            {
                fprintf(fp, ",");
            }
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
    printf("%s にデータが保存されました。\n", filepath);
    return 0;
}