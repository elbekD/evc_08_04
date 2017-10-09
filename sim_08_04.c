//simplify input matrix to three-diagonal matrix

#include "task_08_04.h"

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

    double alpha, betta;
    for (int i = 1; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            double tmp = sqrt(A[i*n+i-1] * A[i*n+i-1] + A[j*n+i-1] * A[j*n+i-1]);
            if (tmp < precision) {
                alpha = betta = .0;
            } else {
                alpha = A[i*n+i-1] / tmp;
                betta = -A[j*n+i-1] / tmp;
            }

            for (int k = i - 1; k < n; k++) {
                double a_ik = alpha * A[i*n+k] - betta * A[j*n+k];
                double a_jk = betta * A[i*n+k] + alpha * A[j*n+k];
                A[i*n+k] = a_ik;
                A[j*n+k] = a_jk;
            }

            if (debug) {
                printf("a=%1.9lf, b=%1.9lf\n", alpha, betta);
                printf("Before. i=%d j=%d: \n", i, j);
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++)
                        printf("%1.3lf ", A[i*n+j]);
                    printf("\n");
                }
            }

            tmp = sqrt(A[(i-1)*n+i] * A[(i-1)*n+i] + A[(i-1)*n+j] * A[(i-1)*n+j]);
            if (tmp < precision) {
                alpha = betta = .0;
            } else {
                alpha = A[(i-1)*n+i] / tmp;
                betta = -A[(i-1)*n+j] / tmp;
            }

            for (int k = i - 1; k < n; k++) {
                double a_ki = alpha * A[k*n+i] - betta * A[k*n+j];
                double a_kj = betta * A[k*n+i] + alpha * A[k*n+j];
                A[k*n+i] = a_ki;
                A[k*n+j] = a_kj;
            }

            if (debug) {
                printf("After. i=%d j=%d: \n", i, j);
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++)
                        printf("%1.3lf ", A[i*n+j]);
                    printf("\n");
                }
            }
        }
    }

    return 0;
}

int sim_memsize_08_04(int n) {
    return 0;
}
