#include <cmath>
#include <ctime>
#include <stdio.h>
#include <stdlib.h>

#define MAX_EL_SIZE 50
#define SIZE 1000
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

void mul_matrix_vector(const double* matrix, const double* vector, size_t size, double* result)
{
    for (int i = 0; i < size; i++) {
        result[i] = mul_scalar_vector(matrix + i, vector, size);
    }
}

void generate_right_side(double* result, size_t size)
{
    for (int i = 0; i < SIZE; i++) {
        result[i] = rand() % MAX_EL_SIZE;
    }
}

void fill_diag(double* matrix, size_t size)
{
    for (int i = 0; i < SIZE; i++) {
        for (int j = i; j < SIZE; j++) {
            matrix[i * SIZE + j] = 1;
            matrix[j * SIZE + i] = 1;
            if (i == j) {
                matrix[j * SIZE + i] = 2;
            }
        }
    }
}
double rand_double()
{
    return (double)rand() / RAND_MAX * 4.0 - 2.0;
}

void generate_matrix(double* matrix, size_t size)
{
    double rand_value = 0;

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < i; j++) {
            matrix[i * size + j] = matrix[j * size + i];
        }
        for (int j = i; j < size; j++) {
            matrix[i * size + j] = rand_double();
            if (i == j) {
                matrix[i * size + j] = fabs(matrix[i * size + j]) + 110.0;
            }
        }
    }
}

int main(int argc, char** argv)
{
    struct timespec start, end;
    double total_time;

    srand(time(NULL));

    double t_n = 0;

    double* tmp_vector = create_matrix(SIZE, 1);
    double* b = create_matrix(SIZE, 1);
    double* y_n = create_matrix(SIZE, 1);
    double* x_n = create_matrix(SIZE, 1);
    double* x_n_1 = create_matrix(SIZE, 1);
    double* main_matrix = create_matrix(SIZE, SIZE); //A

    fill_diag(main_matrix, SIZE);
    generate_right_side(b, SIZE);

    generate_matrix(main_matrix, SIZE);

    //start processing
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);

    double b_mod = count_mod_vector(b, SIZE);

    do {
        mul_matrix_vector(main_matrix, x_n, SIZE, tmp_vector);
        sub_vector_vector(tmp_vector, b, y_n, SIZE);

        mul_matrix_vector(main_matrix, y_n, SIZE, tmp_vector);
        t_n = mul_scalar_vector(y_n, tmp_vector, SIZE) / mul_scalar_vector(tmp_vector, tmp_vector, SIZE);

        mul_vector_number(y_n, tmp_vector, t_n, SIZE);
        sub_vector_vector(x_n, tmp_vector, x_n, SIZE);

        printf("%lf\n", count_mod_vector(y_n, SIZE) / b_mod);
    } while (EPS < (count_mod_vector(y_n, SIZE) / b_mod));

    //finish processing
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    total_time = end.tv_sec - start.tv_sec + 0.000000001 * (double)(end.tv_nsec - start.tv_nsec);

    printf("Time elapsed: %lf\n", total_time);
    print_matrix(x_n, 1, SIZE);
    free(tmp_vector);
    free(b);
    free(y_n);
    free(x_n);
    free(main_matrix);

    return 0;
}
