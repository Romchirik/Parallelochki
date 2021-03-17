#include <cmath>
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define DEBUG false

#define MAIN_DIAG 2000
#define SIZE 9600
#define EPS 0.001
#define MAX_EL_SIZE 1

double rand_double(double min, double max)
{
    double range = (max - min);
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

void mul_vector_number(const double* vector, double* result, double mul, size_t size)
{
    for (int i = 0; i < size; i++) {
        result[i] = vector[i] * mul;
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

double* create_matrix(size_t width, size_t height)
{
    return (double*)calloc(width * height, sizeof(double));
}

void print_matrix(double* result, size_t height, size_t width)
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

void mul_matrix_vector(const double* matrix, const double* vector, size_t size, double* result)
{
    for (int i = 0; i < size; i++) {
        double tmp = 0;
        for (int j = 0; j < size; j++) {
            tmp += matrix[i * size + j] * vector[j];
        }
        result[i] = tmp;
    }
}

void generate_task(double* matrix, double* right_side, size_t size)
{
    double* tmp_vector = (double*)malloc(size * sizeof(double));
    for (size_t i = 0; i < size; i++) {
        for (size_t j = i; j < size; j++) {
            double tmp = rand_double(0, MAX_EL_SIZE);
            matrix[i * size + j] = tmp;
            matrix[j * size + i] = tmp;
        }
    }
    for (int i = 0; i < size; i++) {
        matrix[i * size + i] += MAIN_DIAG;
        right_side[i] = rand_double(0, MAX_EL_SIZE);
    }

    if (DEBUG) {
        printf("generated right side:\n");
        print_matrix(right_side, 1, SIZE);
        printf("generated matrix:\n");
        print_matrix(matrix, SIZE, SIZE);
    }

    free(tmp_vector);
}
bool next_step(double* main_matrix, double* solution, double* right_side, size_t size, double accuracy)
{
    double* tmp = create_matrix(SIZE, 1);
    mul_matrix_vector(main_matrix, solution, size, tmp);
    sub_vector_vector(tmp, right_side, tmp, size);

    double b_mod = sqrt(mul_scalar_vector(right_side, right_side, SIZE));
    double ax_b_mod = sqrt(mul_scalar_vector(tmp, tmp, SIZE));
    printf("a = %lf\n", ax_b_mod / b_mod);
    if (ax_b_mod / b_mod < accuracy) {
        return false;
    } else {
        return true;
    }
    free(tmp);
}

int main(int argc, char** argv)
{

    struct timespec start, end;
    double total_time;

    double t_n = 0;

    double* tmp_vector = create_matrix(SIZE, 1);
    double* b = create_matrix(SIZE, 1);
    double* y_n = create_matrix(SIZE, 1);
    double* x_n = create_matrix(SIZE, 1);
    double* main_matrix = create_matrix(SIZE, SIZE); //A

    generate_task(main_matrix, b, SIZE);

    //start processing
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);

    do {
        mul_matrix_vector(main_matrix, x_n, SIZE, tmp_vector); //count Ax
        if (DEBUG) {
            printf("tmp_vector1:\n");
            print_matrix(tmp_vector, 1, SIZE);
        }

        sub_vector_vector(tmp_vector, b, y_n, SIZE); //count Ax - b = y
        if (DEBUG) {
            printf("tmp_vector:\n");
            print_matrix(tmp_vector, 1, SIZE);
        }

        mul_matrix_vector(main_matrix, y_n, SIZE, tmp_vector); //count Ay
        if (DEBUG) {
            printf("tmp_vector2:\n");
            print_matrix(tmp_vector, 1, SIZE);
        }

        t_n = mul_scalar_vector(y_n, tmp_vector, SIZE) / mul_scalar_vector(tmp_vector, tmp_vector, SIZE);
        if (DEBUG) {
            printf("y * Ay: %lf\n", mul_scalar_vector(y_n, tmp_vector, SIZE));
            printf("yA *Ay: %lf\n", mul_scalar_vector(tmp_vector, tmp_vector, SIZE));
            printf("t_n = %lf\n", t_n);
        }

        mul_vector_number(y_n, tmp_vector, t_n, SIZE);
        if (DEBUG) {
            printf("t_y * t_n:\n");
            print_matrix(tmp_vector, 1, SIZE);
        }

        sub_vector_vector(x_n, tmp_vector, x_n, SIZE);
        if (DEBUG) {
            printf("x_n:\n");
            print_matrix(x_n, 1, SIZE);
            printf("===========================================\n");
            sleep(1);
        }

    } while (next_step(main_matrix, x_n, b, SIZE, EPS));

    //finish processing
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    total_time = end.tv_sec - start.tv_sec + 0.000000001 * (double)(end.tv_nsec - start.tv_nsec);

    printf("Time elapsed: %lf\n", total_time);
    if (DEBUG) {
        printf("solution:\n");
        print_matrix(x_n, 1, SIZE);
    }
    free(tmp_vector);
    free(b);
    free(y_n);
    free(x_n);
    free(main_matrix);

    return 0;
}
