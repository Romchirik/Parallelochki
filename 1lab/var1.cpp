#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define SIZE 3
#define EPS 0.00001

void mul_matrix_matrix(double* first, size_t height1, size_t width1, double* second, size_t height2, size_t width2, double* result)
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

void mul_matrix_number(double* matrix, size_t height, size_t width, double number, double* result)
{
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            result[i * width + j] = matrix[i * width + j] * number;
        }
    }
}

void sum_vector(double* first, double* second, double* result, size_t size)
{
    for (int i = 0; i < size; i++) {
        result[i] = first[i] + second[i];
    }
}

void sub_matrix_matrix(double* first, size_t height, size_t width, double* second, double* result)
{
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            result[i * width + j] = first[i * width + j] - second[i * width + j];
        }
    }
}

double mul_scalar_vector(double* first, size_t size, double* second)
{
    double tmp = 0;
    for (int i = 0; i < size; i++) {
        tmp += first[i] * second[i];
    }
    return tmp;
}

double count_mod_vector(double* vector, size_t size)
{
    double sum = 0;
    for (int i = 0; i < size; i++) {
        sum += vector[i] * vector[i];
    }
    return sqrt(sum);
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

int main(int argc, char** argv)
{
    srand(0);
    double t_n = 0;
    double* b = create_matrix(SIZE, 1);
    double* y_n = create_matrix(SIZE, 1);
    double* x_n = create_matrix(SIZE, 1);
    double* tmp_vector = create_matrix(SIZE, 1);
    double* axn_b = create_matrix(SIZE, 1);
    double* main_matrix = create_matrix(SIZE, SIZE);

    for (int i = 0; i < SIZE; i++) {
        b[i] = 1;
    }

    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            if (i == j) {
                main_matrix[i * SIZE + j] = 2;
            } else {
                main_matrix[i * SIZE + j] = 1;
            }
        }
    }

    do {
        mul_matrix_matrix(main_matrix, SIZE, SIZE, x_n, SIZE, 1, y_n);
        sub_matrix_matrix(y_n, SIZE, 1, b, y_n);
        mul_matrix_matrix(main_matrix, SIZE, SIZE, y_n, SIZE, 1, axn_b);
        t_n = mul_scalar_vector(y_n, SIZE, tmp_vector) / mul_scalar_vector(tmp_vector, SIZE, tmp_vector);

        mul_matrix_number(y_n, SIZE, 1, t_n, tmp_vector);
        sub_matrix_matrix(x_n, SIZE, 1, tmp_vector, x_n);

    } while ((count_mod_vector(axn_b, SIZE) / count_mod_vector(b, SIZE)) > EPS);

    print_matrix(b, 1, SIZE);
    print_matrix(x_n, 1, SIZE);
}
