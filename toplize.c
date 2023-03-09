#include <stdio.h>
#include <math.h>

#define ROWS 6
#define COLS 6

float toeplitz[ROWS][COLS];

int main() {
    // 第一行元素
    float input[] = {0, 0.05, 0.0866025, 0.1, 0.0866025, 0.2};

    // 生成 Toeplitz 矩阵
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLS; j++) {
            int k = abs(i - j);
            toeplitz[i][j] = input[k];
        }
    }

    // 打印矩阵
    printf("Toeplitz matrix:\n");
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLS; j++) {
            printf("%.6f ", toeplitz[i][j]);
        }
        printf("\n");
    }

    return 0;
}
