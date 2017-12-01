//simplify input matrix to three-diagonal matrix

#include "task_08_04.h"

/**
 * Функция, проверяющая входную матрицу на симметричность
 * @param  n размерность матрицы
 * @param  A указатель на матрицу
 * @return   0 - матрица не симметрична
 *           1 - матрица симметрична
 */
int is_simmetric(int n, double* A) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (A[i * n + j] != A[j * n + i]) return 0;
        }
    }
    return 1;
}

int sim_08_04(int n, double* A, double* tmp, double precision) {
    if (!is_simmetric(n, A))
        return -1;

    for (int i = 1; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            double temp = sqrt(A[i*n+i-1] * A[i*n+i-1] + A[j*n+i-1] * A[j*n+i-1]);
            
            if (temp <= precision) continue;

            double alpha = A[i*n+i-1] / temp;
            double betta = -A[j*n+i-1] / temp;

            for (int k = 0; k < n; k++) {
                double a_ik = alpha * A[i*n+k] - betta * A[j*n+k];
                double a_jk = betta * A[i*n+k] + alpha * A[j*n+k];
                A[i*n+k] = a_ik;
                A[j*n+k] = a_jk;
            }

            for (int k = 0; k < n; k++) {
                double a_ki = alpha * A[k*n+i] - betta * A[k*n+j];
                double a_kj = betta * A[k*n+i] + alpha * A[k*n+j];
                A[k*n+i] = a_ki;
                A[k*n+j] = a_kj;
            }
        }
    }
    return 0;
}

int sim_memsize_08_04(int n) {
    return 0;
}
