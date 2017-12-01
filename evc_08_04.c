//find eigenvalue of matrix using LR method
#include "task_08_04.h"

/**
 * Функция LR разложения трехдиагональной матрицы.
 * Матрицы L и R двухдиагональны, при этом у матрицы L на главной диагонали единицы.
 * Матрица L хранится под главной диагональю исходной матрицы.
 * Матрица R хранится на главной диагонали и над главной диагональю исходной матрицы.
 * @param n         размерность матрицы
 * @param A         указатель на матрицу
 * @param precision "точность локализации нуля"
 */
void LR_decomposition(int n, double* A, double precision) {
    for (int i = 1; i < n; i++) {
        double l = fabs(A[(i-1)*n+i-1]) <= precision || fabs(A[i*n+i-1]) <= precision ? .0 : A[i*n+i-1] / A[(i-1)*n+i-1];
        for (int k = i; k <= i + 1 && k < n; k++) {
            double r = A[i*n+k] - l * A[(i-1)*n+k];
            A[i*n+k] = r;
        }
        A[i*n+i-1] = l;
    } 
}

/**
 * Функция, вычисляющая произведение матрицы R на L, получившихся в результате LR разложения.
 * Диагональные элементы, получившиеся в результате данного произведения, будут приближениями
 * к собственным значениям матрицы
 * @param n размерность матрицы
 * @param A указатель на матрицу
 */
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

/**
 * Функция, определяющая, достигла ли точность вычисления собственных значений указанного epsilon.
 * Под точностью подразумевается сходимость матрицы L к единичной матрице.
 * @param  n       размерность матрицы
 * @param  epsilon точность вычисления
 * @param  A       указатель на матрицу
 */
int is_epsilon_reached(int n, double epsilon, double* A) {
    for (int i = 1; i < n; i++) {
        if (fabs(A[i*n+i-1]) > epsilon) {
            return 0;
        }
    }
    return 1;
}

int evc_08_04(int n, int max_iterations, double epsilon, double* A, double* E, double* tmp, double precision) {
    if (max_iterations <= 0) {
        while (!is_epsilon_reached(n, epsilon, A)) {
            LR_decomposition(n, A, precision);
            compute_next_A(n, A);
        }
        for (int i = 0; i < n; i++)
            E[i] = A[i*n+i];
        return 0;
    } else {
        for (int i = 0; i < max_iterations; i++) {
            LR_decomposition(n, A, precision);
            compute_next_A(n, A);
            if (is_epsilon_reached(n, epsilon, A)) {
                for (int j = 0; j < n; j++)
                    E[j] = A[j*n+j];
                return 0;
            }
        }
        return 1;
    }
    return -1;
}

int evc_memsize_08_04(int n) {
    return 0;
}
