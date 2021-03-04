#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define SIZE 100
#define EPS 0.00001


void mul_vector_number(const double *matrix, double *result, double mul, size_t size) {
    for (int i = 0; i < size; i++) {
        result[i] = matrix[i] * mul;
    }
}


void sub_vector_vector(const double *first, const double *second, double *result, size_t size) {
    for (int i = 0; i < size; i++) {
        result[i] = first[i] - second[i];
    }
}

double mul_scalar_vector(const double *first, const double *second, size_t size) {
    double tmp = 0;
    for (int i = 0; i < size; i++) {
        tmp += first[i] * second[i];
    }
    return tmp;
}

double count_mod_vector(const double *vector, size_t size) {
    double tmp = 0;
    for (int i = 0; i < size; i++) {
        tmp = tmp + vector[i] * vector[i];
    }
    return sqrt(tmp);
}

double *create_matrix(size_t width, size_t height) {
    return (double *) calloc(width * height, sizeof(double));
}

void print_matrix(double *result, size_t height, size_t width) {
    for (int i = 0; i < height; i++) {
        for (int j = 90; j < width; j++) {
            printf("%lf  ", result[i * width + j]);
        }
        printf("\n");
    }
}

void swap(double **a, double **b) {
    double *tmp = *a;
    *a = *b;
    *b = tmp;
}

void mul_matrix_vector(const double *matrix, const double *vector, size_t size, double *result) {
    for (int i = 0; i < size; i++) {
        result[i] = mul_scalar_vector(matrix + i, vector, size);
    }
}


int main(int argc, char **argv) {
    struct timespec start, end;
    double total_time;

    //srand(time(NULL));

    double t_n = 0;
    double *u = create_matrix(SIZE, 1);
    double *b = create_matrix(SIZE, 1);
    double *y_n = create_matrix(SIZE, 1);
    double *ay_n = create_matrix(SIZE, 1);
    double *x_n = create_matrix(SIZE, 1);
    double *t_y_n = create_matrix(SIZE, 1);
    double *ax_n = create_matrix(SIZE, 1);
    double *main_matrix = create_matrix(SIZE, SIZE);

    for (int i = 0; i < SIZE; i++) {
        u[i] = rand() % 10;
    }

    double main_diag = 2;
    for (int i = 0; i < SIZE; i++) {
        for (int j = i; j < SIZE; j++) {
            main_matrix[i * SIZE + j] = 1;
            main_matrix[j * SIZE + i] = 1;
            if (i == j) {
                continue;
            }
        }
        main_matrix[i * SIZE + i] = main_diag;
    }


    //print_matrix(main_matrix, SIZE, SIZE);
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    mul_matrix_vector(main_matrix, u, SIZE, b);
    double b_mod = count_mod_vector(b, SIZE);

    int counter = 0;
    do {
        mul_matrix_vector(main_matrix, x_n, SIZE, ax_n);
        sub_vector_vector(ax_n, b, y_n, SIZE);

        mul_matrix_vector(main_matrix, y_n, SIZE, ay_n);
        t_n = mul_scalar_vector(y_n, ay_n, SIZE) / mul_scalar_vector(ay_n, ay_n, SIZE);
        mul_vector_number(y_n, t_y_n, t_n, SIZE);
        sub_vector_vector(x_n, t_y_n, x_n, SIZE);
        printf("%d\n", counter);
        counter++;
    } while (EPS < (count_mod_vector(y_n, SIZE) / b_mod));

    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    total_time = end.tv_sec - start.tv_sec + 0.000000001 * (double) (end.tv_nsec - start.tv_nsec);

    printf("Time elapsed: %lf\n", total_time);
    print_matrix(u, 1, SIZE);
    //print_matrix(b, 1, SIZE);
    print_matrix(x_n, 1, SIZE);

    double t_n = 0;
    double *u = create_matrix(SIZE, 1);
    double *b = create_matrix(SIZE, 1);
    double *y_n = create_matrix(SIZE, 1);
    double *ay_n = create_matrix(SIZE, 1);
    double *x_n = create_matrix(SIZE, 1);
    double *t_y_n = create_matrix(SIZE, 1);
    double *ax_n = create_matrix(SIZE, 1);
    double *main_matrix = create_matrix(SIZE, SIZE);

    free(u);
    free(b);
    free
}
