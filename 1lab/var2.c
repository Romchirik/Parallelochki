#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAX_U_EL 10
#define SIZE 100
#define EPS 0.00001

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

double mul_scalar_vector(const double* first, const double* second, size_t size)
{
    double tmp = 0;
    for (int i = 0; i < size; i++) {
        tmp += first[i] * second[i];
    }
    return tmp;
}

double count_mod_vector(const double* vector, size_t size)
{
    double tmp = 0;
    for (int i = 0; i < size; i++) {
        tmp = tmp + vector[i] * vector[i];
    }
    return sqrt(tmp);
}

double* create_matrix(size_t width, size_t height)
{
    return (double*)calloc(width * height, sizeof(double));
}

void print_matrix(double* result, size_t height, size_t width)
{
    for (int i = 0; i < height; i++) {
        for (int j = 90; j < width; j++) {
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

void mul_matrix_vector(const double* matrix, const double* vector, size_t size, double* result)
{
    for (int i = 0; i < size; i++) {
        result[i] = mul_scalar_vector(matrix + i, vector, size);
    }
}

void generate_right_side(double* main_matrix, double* result, size_t size)
{
    double* tmp = create_matrix(SIZE, 1);
    for (int i = 0; i < SIZE; i++) {
        tmp[i] = rand() % MAX_U_EL;
    }
    mul_matrix_vector(main_matrix, tmp, SIZE, result);
    //printf("generated result:\n");
    //print_matrix(tmp, 1, SIZE);
    free(tmp);
}

void fill_diag(double* matrix, size_t size)
{
    for (int i = 0; i < SIZE; i++) {
        for (int j = i; j < SIZE; j++) {
            matrix[i * SIZE + j] = 1;
            matrix[j * SIZE + i] = 1;
        }
    }
    for (int i = 0; i < size; i++) {
        matrix[i * size + i] = 2;
    }
}

int main(int argc, char** argv)
{
    //time stuff
    struct timespec start, end;
    double total_time;

    //math stuff
    double t_n = 0;
    double* tmp_vector = create_matrix(SIZE, 1);
    double* b = create_matrix(SIZE, 1);
    double* y_n = create_matrix(SIZE, 1);
    double* x_n = create_matrix(SIZE, 1);
    double* main_matrix = create_matrix(SIZE, SIZE); //A

    //mpi stuff
    int size, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Получение числа процессов
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    printf("my_rank: %d", rank);
    if (rank == 0) {
        srand(time(NULL));
        fill_diag(main_matrix, SIZE);
        generate_right_side(main_matrix, b, SIZE);
    }

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    MPI_Bcast(main_matrix, SIZE * SIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(b, SIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //print_matrix(main_matrix, SIZE, SIZE);
    //     mul_matrix_vector(main_matrix, x_n, SIZE, tmp_vector);
    //     sub_vector_vector(tmp_vector, b, y_n, SIZE);

    //     mul_matrix_vector(main_matrix, y_n, SIZE, tmp_vector);
    //     t_n = mul_scalar_vector(y_n, tmp_vector, SIZE) / mul_scalar_vector(tmp_vector, tmp_vector, SIZE);
    //     mul_vector_number(y_n, tmp_vector, t_n, SIZE);
    //     sub_vector_vector(x_n, tmp_vector, x_n, SIZE);
    // }
    // while (EPS < (count_mod_vector(y_n, SIZE) / b_mod))
    //     ;

    // //finish processing
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    total_time = end.tv_sec - start.tv_sec + 0.000000001 * (double)(end.tv_nsec - start.tv_nsec);

    //print_matrix(x_n, 1, SIZE);
    //printf("Time elapsed: %lf\n", total_time);

    free(tmp_vector);
    free(b);
    free(y_n);
    free(x_n);
    free(main_matrix);

    MPI_Finalize();
    return 0;
}
