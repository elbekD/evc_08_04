#include "task_08_04.h"

#define DEFAULT_IN_FILENAME "08_04_in.txt"
#define DEFAULT_OUT_FILENAME "08_04_out.txt"
#define PRECISION_KEY "-prec="
#define EPSILON_KEY "-eps="
#define ITERATION_KEY "-max_iter="

void hint() {
    printf("Usage: evc [input_file_name] [output_file_name] [options]\n"
           "  Where options include:\n"
           "  -d    print debug messages [default OFF]\n"
           "  -e    print errors [default OFF]\n"
           "  -p    print matrix [default OFF]\n"
           "  -t    print execution time [default OFF]\n"
           "  -prec=<num>       precision [default - 1e-14]\n"
           "  -eps=<num>        'epsilon' [default - 1e-10]\n"
           "  -max_iter=<num>   limit number of iterations [default - 0, i.e. not limit]\n"
           "  -h, -?     print this and exit\n");
}

int main(int argc, char const *argv[])
{
    debug = 0;
    error = 0;

    char print_matrix = 0;
    char print_time = 0;
    double precision = 1e-14;
    double epsilon = 1e-10;
    int max_iterations = 0;
    char h = 0;

    double* A = NULL;
    double* E = NULL;
    double* sim_tmp = NULL;
    double* evc_tmp = NULL;

    char input_filename[128] = DEFAULT_IN_FILENAME;
    char output_filename[128] = DEFAULT_OUT_FILENAME;

    if (argc > 1) {
        for (int i = 1; i < argc; i++) {
            if (i == 1 && argv[i][0] != '-')
                strcpy(input_filename, argv[i]);
            else if (i == 2 && argv[i][0] != '-')
                strcpy(output_filename, argv[i]);
            else if (strcmp("-d", argv[i]) == 0) debug = 1;
            else if (strcmp("-e", argv[i]) == 0) error = 1;
            else if (strcmp("-p", argv[i]) == 0) print_matrix = 1;
            else if (strcmp("-t", argv[i]) == 0) print_time = 1;
            else if (strcmp("-h", argv[i]) == 0 || strcmp("-?", argv[i]) == 0) { h = 1; break; }
            else if (strncmp(PRECISION_KEY, argv[i], strlen(PRECISION_KEY)) == 0) {
                int n = sscanf(argv[i] + strlen(PRECISION_KEY), "%lf", &precision);
                if (n != 1) return __LINE__;
            } else if (strncmp(EPSILON_KEY, argv[i], strlen(EPSILON_KEY)) == 0) {
                int n = sscanf(argv[i] + strlen(EPSILON_KEY), "%lf", &epsilon);
                if (n != 1) return __LINE__;
            } else if (strncmp(ITERATION_KEY, argv[i], strlen(ITERATION_KEY)) == 0) {
                int n = sscanf(argv[i] + strlen(ITERATION_KEY), "%d", &max_iterations);
                if (n != 1) return __LINE__;
            } else {
                if (error) fprintf(stderr, "Error occured during argument parsing\n");
                return __LINE__;
            }
        }
    }

    if (h) {
        hint();
        return 0;
    }

    int dim;
    FILE* in = fopen(input_filename, "r");
    FILE* out = fopen(output_filename, "w");
    if (in == NULL || out == NULL) {
        if (error) fprintf(stderr, "Input or output file open error\n");
        fclose(in);
        fclose(out);
        return __LINE__;
    }

    if (fscanf(in, "%d", &dim) != 1) {
        if (error) fprintf(stderr, "Can not read matrix dimension\n");
        return __LINE__;
    }

    A = (double*)malloc(dim * dim * sizeof(double));
    E = (double*)malloc(dim * sizeof(double));
    sim_tmp = (double*)malloc(sizeof(double) * sim_memsize_08_04(dim));
    evc_tmp = (double*)malloc(sizeof(double) * evc_memsize_08_04(dim));

    if (A == NULL || E == NULL || sim_tmp == NULL || evc_tmp == NULL) {
        if (error) fprintf(stderr, "Memory allocation error\nA=%p, E=%p, sim_tmp=%p, evc_tmp=%p", A, E, sim_tmp, evc_tmp);
        free(A);
        free(E);
        free(sim_tmp);
        free(evc_tmp);
        return __LINE__;
    }

    for (int i = 0; i < dim * dim; i++) {
        if (fscanf(in, "%lf", &A[i]) != 1) {
            if (error) fprintf(stderr, "Error occured during reading matrix A\n");
            return __LINE__;
        }
    }

    time_t overall_time = 0;
    
    time_t s_time = clock();
    int sim_res = sim_08_04(dim, A, sim_tmp, precision);
    time_t e_time = clock();
    
    overall_time += e_time - s_time;

    if (sim_res == -1) {
        if (debug) fprintf(stderr, "The simplification method is not applicable\n");
        fprintf(out, "0\n");
    } else {
        if (print_matrix) {
            printf("\n");
            if (debug) printf("SIM result\n");
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++)
                    printf("%1.3lf ", A[i*dim+j]);
                printf("\n");
            }
        }

        for (int i = 0; i < dim; i++) E[i] = .0;
        s_time = clock();
        int evc_res = evc_08_04(dim, max_iterations, epsilon, A, E, evc_tmp, precision);
        e_time = clock();
        overall_time += e_time - s_time;

        if (print_matrix) {
            printf("\n");
            if (debug) printf("EVC result\n");
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++)
                    printf("%1.3lf ", A[i*dim+j]);
                printf("\n");
            }
        }

        if (evc_res == -1) {
            if (debug) fprintf(stderr, "EVC method is not applicable\n");
            fprintf(out, "0\n");
        } else if (evc_res == 1) {
            if (debug) fprintf(stderr, "EVC method does not converge. Iterations: %d\n", max_iterations);
            fprintf(out, "1\n");
        } else if (evc_res == 0) {
            fprintf(out, "%d\n", dim);
            for (int i = 0; i < dim; i++)
                fprintf(out, "%1.9lf\n", E[i]);
        }
    }

    if (print_time)
        printf("Overall time: %lf\n", (double)overall_time / CLOCKS_PER_SEC);

    free(A);
    free(E);
    free(sim_tmp);
    free(evc_tmp);
    fclose(in);
    fclose(out);
    return 0;
}