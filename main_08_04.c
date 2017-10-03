#include "task_08_04.h"

#define DEFAULT_IN_FILENAME "08_04_in.txt"
#define DEFAULT_OUT_FILENAME "08_04_out.txt"
#define PRECISION_KEY "-prec="
#define EPSILON_KEY "-eps="
#define ITERATION_KEY "-max_iter="
//interface
//parse args, alloc mem, write/read to/from file
//must print overall solving time

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
    double* tmp = NULL;

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

    int n;
    FILE* in = fopen(input_filename, "r");
    FILE* out = fopen(output_filename, "w");
    if (in == NULL || out == NULL) {
        if (error) fprintf(stderr, "Input or output file open error\n");
        return __LINE__;
    }

    if (fscanf(in, "%d", &n) != 1) {
        if (error) fprintf(stderr, "Matrix dimension '%d' incorrect\n", n);
        return __LINE__;
    }

    A = (double*)malloc(n * n * sizeof(double));
    E = (double*)malloc(n * sizeof(double));

    if (A == NULL || E == NULL) {
        if (error) fprintf(stderr, "Memory allocation error\n");
        return __LINE__;
    }

    for (int i = 0; i < n * n; i++) {
        if (fscanf(in, "%lf", &A[i]) != 1) {
            if (error) fprintf(stderr, "Error occured during reading matrix A\n");
            return __LINE__;
        }
    }

    time_t overall_time = 0;
    time_t s_time = clock();
    int sim_res = sim_08_04(n, A, tmp, precision);
    time_t e_time = clock();
    overall_time += e_time - s_time;
    if (debug) printf("Simplification time: %lf\n", (double)(e_time - s_time) / CLOCKS_PER_SEC);

    if (sim_res == -1) {
        if (debug) fprintf(stderr, "The simplification method is not applicable\n");
        fprintf(out, "0\n");
    } else {
        s_time = clock();
        int evc_res = evc_08_04(n, max_iterations, epsilon, A, E, tmp, precision);
        e_time = clock();
        overall_time += e_time - s_time;
        if (evc_res == -1) {
            if (debug) fprintf(stderr, "EVC method is not applicable\n");
            fprintf(out, "0\n");
        } else if (evc_res == 1) {
            if (debug) fprintf(stderr, "EVC method does not converge. Iterations: %d\n", max_iterations);
            fprintf(out, "1\n");
        } else if (evc_res == 0) {
            fprintf(out, "%d\n", n);
            for (int i = 0; i < n; i++)
                fprintf(out, "%1.9lf\n", E[i]);
        }
    }

    if (print_time)
        printf("Overall time: %lf\n", (double)overall_time / CLOCKS_PER_SEC);

    free(A);
    free(E);
    free(tmp);
    fclose(in);
    fclose(out);
    return 0;
}