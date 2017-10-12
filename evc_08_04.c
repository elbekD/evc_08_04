//find eigenvalue of matrix using LR method
#include "task_08_04.h"

void LR_decomposition(int n, double* A, double precision) {
    for (int i = 1; i < n; i++) {
        double l = fabs(A[(i-1)*n+i-1]) < precision ? .0 : A[i*n+i-1] / A[(i-1)*n+i-1];
        for (int k = i; k <= i + 1; k++) {
            double r = A[i*n+k] - l * A[(i-1)*n+k];
            A[i*n+k] = r;
        }
        A[i*n+i-1] = l;
    }
}

void compute_next_A(int n, double* A) {
    A[0*n+0] = A[0*n+0] + A[0*n+1] * A[1*n+0];
    A[0*n+1] = A[0*n+1];
    for (int i = 1; i < n - 1; i++) {
        A[i*n+i-1] = A[i*n+i] * A[i*n+i-1];
        A[i*n+i] = A[i*n+i] + A[i*n+i+1] * A[(i+1)*n+i];
        A[i*n+i+1] = A[i*n+i+1];
    }
    A[(n-1)*n+n-2] = A[(n-1)*n+n-2] * A[(n-1)*n+n-1];
    A[(n-1)*n+n-1] = A[(n-1)*n+n-1];
}

int is_epsilon_reached(int n, double epsilon, double* A, double* E) {
    double tmp = .0;
    for (int i = 0; i < n; i++)
        tmp += (A[i*n+i] - E[i]) * (A[i*n+i] - E[i]);
    return sqrt(tmp) < epsilon;
}

int evc_08_04(int n, int max_iterations, double epsilon, double* A, double* E, double* tmp, double precision) {
    if (max_iterations <= 0) {
        while (!is_epsilon_reached(n, epsilon, A, E)) {
            for (int i = 0; i < n; i++)
                E[i] = A[i*n+i];
            LR_decomposition(n, A, precision);
            compute_next_A(n, A);
        }
        return 0;
    } else {
        for (int i = 0; i < max_iterations; i++) {
            for (int i = 0; i < n; i++)
                E[i] = A[i*n+i];
            LR_decomposition(n, A, precision);
            compute_next_A(n, A);
            if (is_epsilon_reached(n, epsilon, A, E))
                return 0;
        }
        return 1;
    }
    return -1;
}

int evc_memsize_08_04(int n) {
    return 0;
}
