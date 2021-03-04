#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define SIZE 3
#define EPS 0.00001

void mul_matrix_matrix(const double* first, size_t height1, size_t width1, const double* second, size_t height2,
    size_t width2, double* result)
{
    for (int i = 0; i < height1; i++) {
        for (int j = 0; j < width2; j++) {
            double sum = 0;
            for (int k = 0; k < width1; k++) {
                sum += first[i * width1 + k] * second[k * width2 + j];
            }
            result[i * width2 + j] = sum;
        }
    }
}

void mul_matrix_number(const double* matrix, size_t height, size_t width, double number, double* result)
{
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            result[i * width + j] = matrix[i * width + j] * number;
        }
    }
}

void sum_vector(const double* first, const double* second, double* result, size_t size)
{
    for (int i = 0; i < size; i++) {
        result[i] = first[i] + second[i];
    }
}

void sub_matrix_matrix(const double* first, size_t height, size_t width, const double* second, double* result)
{
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            result[i * width + j] = first[i * width + j] - second[i * width + j];
        }
    }
}

double mul_scalar_vector(const double* first, size_t size, const double* second)
{
    double tmp = 0;
    for (int i = 0; i < size; i++) {
        tmp += first[i] * second[i];
    }
    return tmp;
}

double count_mod_vector(double* vector, size_t size)
{
    double tmp = 0;
    for (int i = 0; i < size; i++) {
        tmp += vector[i] * vector[i];
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

int main(int argc, char** argv)
{
    struct timespec start {
    }, end {};
    double total_time;

    srand(time(NULL));
    double b_mod = 0;
    double t_n = 0;
    double* u = create_matrix(SIZE, 1);
    double* b = create_matrix(SIZE, 1);
    double* y_n = create_matrix(SIZE, 1);
    double* ay_n = create_matrix(SIZE, 1);
    double* x_n = create_matrix(SIZE, 1);
    double* t_y_n = create_matrix(SIZE, 1);
    double* ax_n = create_matrix(SIZE, 1);
    double* main_matrix = create_matrix(SIZE, SIZE);

    for (int i = 0; i < SIZE; i++) {
        u[i] = rand() % 20;
    }

    for (int i = 0; i < SIZE; i++) {
        for (int j = i; j < SIZE; j++) {
            double tmp = rand() % 20;
            main_matrix[i * SIZE + j] = tmp;
            main_matrix[j * SIZE + i] = tmp;
            if (i == j) {
                main_matrix[i * SIZE + j] += 200;
            }
        }
    }

    b_mod = count_mod_vector(b, SIZE);
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    mul_matrix_matrix(main_matrix, SIZE, SIZE, u, SIZE, 1, b);

    do {
        mul_matrix_matrix(main_matrix, SIZE, SIZE, x_n, SIZE, 1, ax_n);
        sub_matrix_matrix(ax_n, SIZE, 1, b, y_n);

        mul_matrix_matrix(main_matrix, SIZE, SIZE, y_n, SIZE, 1, ay_n);
        t_n = mul_scalar_vector(y_n, SIZE, ay_n) / mul_scalar_vector(ay_n, SIZE, ay_n);
        mul_matrix_number(y_n, SIZE, 1, t_n, t_y_n);
        sub_matrix_matrix(x_n, SIZE, 1, t_y_n, x_n);

        print_matrix(u, 1, SIZE);
        print_matrix(x_n, 1, SIZE);

    } while (count_mod_vector(y_n, SIZE) / b_mod > EPS);

    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    total_time = end.tv_sec - start.tv_sec + 0.000000001 * (double)(end.tv_nsec - start.tv_nsec);

    printf("Time elapsed: %lf\n", total_time);
    print_matrix(u, 1, SIZE);
    //print_matrix(b, 1, SIZE);
    print_matrix(x_n, 1, SIZE);
}
