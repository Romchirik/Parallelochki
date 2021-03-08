#include <cmath>
#include <ctime>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define MAX_EL_SIZE 100
#define SIZE 4
#define EPS 0.00001

void print_matrix(const double* result, size_t height, size_t width);

void mul_vector_number(const double* matrix, double* result, double mul, size_t size)
{
    for (int i = 0; i < size; i++) {
        result[i] = matrix[i] * mul;
    }
}

void sub_vector_vector(const double* first, const double* second, double* result, size_t size)
{
    for (int i = 0; i < size; i++) {
        result[i] = first[i] - second[i];
    }
}

double mul_scalar_vector_linear(const double* first, const double* second, size_t size)
{
    double result = 0;
    for (int i = 0; i < size; i++) {
        result += first[i] * second[i];
    }
    return result;
}

double mul_scalar_vector(const double* first, const double* second, size_t size, int proc_num, int rank)
{
    //    double curr_p_result = 0, result = 0;

    double result = 0;
    for (int i = 0; i < size; i++) {
        result += first[i] * second[i];
    }
    if (rank == 0) {
        printf("scalar mul result: %lf\n", result);
        if (result < 0.000000001) {
            printf("AAAAAAAAAALARM!!!\n");
            print_matrix(first, 1, SIZE);
            print_matrix(second, 1, SIZE);
        }
    }

    return result;

    // size_t num_lines_to_count = size / proc_num;
    // size_t shift = num_lines_to_count * rank;
    // for (int i = shift; i < shift + num_lines_to_count; i++) {
    //     curr_p_result += first[i] * second[i];
    // }
    // MPI_Allreduce(&curr_p_result, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    // return result;
}

double count_mod_vector(const double* vector, size_t size, int num_processes, int rank)
{
    return sqrt(mul_scalar_vector(vector, vector, size, num_processes, rank));
}

double* create_matrix(size_t width, size_t height)
{
    return (double*)calloc(width * height, sizeof(double));
}

void print_matrix(const double* result, size_t height, size_t width)
{
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            printf("%lf  ", result[i * width + j]);
        }
        printf("\n");
    }
}

void swap(double** a, double** b)
{
    double* tmp = *a;
    *a = *b;
    *b = tmp;
}

void mul_matrix_vector_linear(const double* matrix, const double* vector, size_t size, double* result)
{
    for (int i = 0; i < size; i++) {
        result[i] = mul_scalar_vector_linear(matrix + i, vector, size);
    }
}

void mul_matrix_vector(const double* matrix, const double* vector, size_t size, double* result, int proc_num, int rank)
{
    double* buffer = (double*)malloc(sizeof(double) * size);
    size_t num_lines_to_count = size / proc_num;
    size_t shift = num_lines_to_count * rank;

    for (int i = 0; i < num_lines_to_count; i++) {
        for (int j = 0; j < size; j++) {
            result[i + shift] = matrix[i * size + j] * vector[j];
        }
    }

    MPI_Allgather(result + shift, num_lines_to_count, MPI_DOUBLE, buffer, num_lines_to_count, MPI_DOUBLE, MPI_COMM_WORLD);
    for (int i = 0; i < size; i++) {
        result[i] = buffer[i];
    }
    //memcpy(result, buffer, sizeof(double) * size);
}

double rand_double()
{
    return abs((double)rand() / RAND_MAX * 4.0 - 2.0);
}

void generate_matrix(double* matrix, size_t width, size_t height, size_t proc_num, size_t rank)
{
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            matrix[i * width + j] = rand_double();
        }
    }

    for (int i = 0; i < height; i++) {
        matrix[i * width + height * rank + i] = rand() % MAX_EL_SIZE;
    }
}

void generate_right_side(double* result, size_t size)
{
    for (int i = 0; i < SIZE; i++) {
        result[i] = rand() % MAX_EL_SIZE;
    }
}

int main(int argc, char* argv[])
{
    struct timespec start, end;
    double total_time;

    //mpi stuff
    int proc_num, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num); // get number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //get process idenifier

    //math stuff
    double t_n = 0;
    double* tmp_vector = create_matrix(SIZE, 1);
    double* b = create_matrix(SIZE, 1);
    double* y_n = create_matrix(SIZE, 1);
    double* x_n = create_matrix(SIZE, 1);
    double* main_matrix = create_matrix(SIZE, SIZE / proc_num); //A

    generate_matrix(main_matrix, SIZE, SIZE / proc_num, proc_num, rank);
    if (rank == 0) {
        generate_right_side(b, SIZE);
    }

    // start_processing
    if (rank == 0) {
        clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    }
    printf("kavo1\n");
    fflush(stdout);

    MPI_Bcast(b, SIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double b_mod = count_mod_vector(b, SIZE, proc_num, rank);
    double a;
    do {

        mul_matrix_vector(main_matrix, x_n, SIZE, tmp_vector, proc_num, rank);
        // if (rank == 0) {
        //     printf("tmp_matrix:\n");
        //     print_matrix(tmp_vector, 1, SIZE);
        //     printf("b:\n");
        //     print_matrix(b, 1, SIZE);
        // }

        sub_vector_vector(tmp_vector, b, y_n, SIZE);

        mul_matrix_vector(main_matrix, y_n, SIZE, tmp_vector, proc_num, rank);
        if (rank == 0) {
            printf("tmp_matrix:\n");
            print_matrix(tmp_vector, 1, SIZE);
            printf("y_n:\n");
            print_matrix(y_n, 1, SIZE);
        }
        t_n = mul_scalar_vector(y_n, tmp_vector, SIZE, proc_num, rank) / mul_scalar_vector(tmp_vector, tmp_vector, SIZE, proc_num, rank);
        usleep(1000000);
        mul_vector_number(y_n, tmp_vector, t_n, SIZE);
        sub_vector_vector(x_n, tmp_vector, x_n, SIZE);
        a = mul_scalar_vector(y_n, y_n, SIZE, proc_num, rank) / b_mod * b_mod;

        if (rank == 0) {
            printf("a = %lf\n", a);
            printf("x_n:\n");
            print_matrix(x_n, 1, SIZE);
            printf("y_n:\n");
            print_matrix(y_n, 1, SIZE);
            printf("t_n %lf:\n ", t_n);
            printf("tmp_matrix:\n");
            print_matrix(tmp_vector, 1, SIZE);
            printf("========================================\n");
        }
    } while (EPS * EPS < a);

    //finish processing

    if (rank == 0) {
        clock_gettime(CLOCK_MONOTONIC_RAW, &end);
        total_time = end.tv_sec - start.tv_sec + 0.000000001 * (double)(end.tv_nsec - start.tv_nsec);
        // print_matrix(x_n, 1, SIZE);
        // printf("Time elapsed: %lf\n", total_time);
    }

    fflush(stdout);
    free(tmp_vector);
    free(b);
    free(y_n);
    free(x_n);
    free(main_matrix);

    MPI_Finalize();
    return 0;
}
