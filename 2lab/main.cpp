#include <math.h>
#include <memory.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/times.h>
#include <time.h>
#include <unistd.h>
#include <xmmintrin.h>

#define MAX_EL_SIZE 10
#define DEBUG false
#define NUM_THREADS 4
#define N 2048
const int M = 10;

void sum_matrix_matrix(float* first, float* second, float* res, int size)
{
    __m128 sum;
#pragma omp parallel for private(sum)
    for (int i = 0; i < size; ++i) {
        __m128 *AA, *BB;
        AA = (__m128*)(first + i * size);
        BB = (__m128*)(second + i * size);
        for (int j = 0; j < N / 4; ++j) {
            sum = _mm_add_ps(AA[j], BB[j]);
            _mm_store_ps((res + i * size + j * 4), sum);
        }
    }
}

void sub_matrix_matrix(float* first, float* second, float* res, int size)
{
    __m128 sub;
#pragma omp parallel for private(sub)
    for (int i = 0; i < size; ++i) {
        __m128 *AA, *BB;
        AA = (__m128*)(first + i * size);
        BB = (__m128*)(second + i * size);

        for (int j = 0; j < size / 4; ++j) {
            sub = _mm_sub_ps(AA[j], BB[j]);
            _mm_store_ps((res + i * size + j * 4), sub);
        }
    }
}

void mul_matrix_matrix(float* first, float* second, float* res, size_t size)
{
#pragma omp parallel
    {
        auto sum = (float*)calloc(size, sizeof(float));
        auto tmp = (float*)calloc(size, sizeof(float));

#pragma omp for
        for (size_t i = 0; i < size; i++) {
            memset(sum, 0, size * sizeof(float));
            for (int j = 0; j < size; j++) {
                __m128 Aij = _mm_set1_ps(first[i * size + j]);

                for (size_t k = 0; k < size; k += 4) {
                    __m128 tmp_vec = _mm_load_ps(second + (j * size + k));
                    __m128 to_store = _mm_mul_ps(Aij, tmp_vec);
                    _mm_store_ps(tmp + k, to_store);
                }

                for (size_t k = 0; k < size; k += 4) {
                    __m128 a = _mm_load_ps(sum + k);
                    __m128 b = _mm_load_ps(tmp + k);
                    __m128 r = _mm_add_ps(a, b);
                    _mm_store_ps(sum + k, r);
                }
            }
            memcpy(res + i * size, sum, sizeof(float) * size);
        }
        free(sum);
        free(tmp);
    }
}

float max_row_sum(float* matrix, int size)
{
    float max = 0;
    for (int i = 0; i < size; i++) {
        float tmp = 0;
#pragma omp parallel reduction(+ \
                               : tmp)
        {
#pragma omp for
            for (int j = 0; j < size; j++) {
                tmp += fabsf(matrix[i * size + j]);
            }
        }

        if (tmp > max)
            max = tmp;
    }
    return max;
}

float max_column_sum(float* A, int size)
{
    float max = 0;
    for (int i = 0; i < size; ++i) {
        float tmp = 0;
#pragma omp parallel reduction(+ \
                               : tmp)
        {
#pragma omp for
            for (int j = 0; j < size; ++j) {
                tmp += fabsf(A[j * size + i]);
            }
        }

        if (tmp > max)
            max = tmp;
    }
    return max;
}

float* reverse_matrix(float* matrix, int size)
{
    float *b_matrix, *i_matrix, *ba, *r, *inv, *buffer;
    float max_sum_row = max_row_sum(matrix, size);
    float max_sum_column = max_column_sum(matrix, size);

    b_matrix = (float*)calloc(size * size, sizeof(float));
    i_matrix = (float*)calloc(size * size, sizeof(float));
    ba = (float*)calloc(size * size, sizeof(float));
    r = (float*)calloc(size * size, sizeof(float));
    inv = (float*)calloc(size * size, sizeof(float));
    buffer = (float*)calloc(size * size, sizeof(float));

    for (int i = 0; i < size; ++i) {
        i_matrix[i * size + i] = 1;
        buffer[i * size + i] = 1;
    }

    //count B

    for (int i = 0; i < size; ++i) {
#pragma omp parallel for
        for (int j = 0; j < size; ++j) {
            b_matrix[i * size + j] = matrix[j * size + i] / (max_sum_column * max_sum_row);
        }
    }

    //count R
    mul_matrix_matrix(matrix, b_matrix, ba, size);
    sub_matrix_matrix(i_matrix, ba, r, size);

    for (int k = 0; k < M - 1; ++k) {
        printf("progress:%d\n", k + 1);
        mul_matrix_matrix(i_matrix, r, ba, size);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                i_matrix[i * size + j] = ba[i * size + j];
            }
        }
        sum_matrix_matrix(buffer, i_matrix, buffer, size);
    }
    mul_matrix_matrix(buffer, b_matrix, inv, size);

    free(b_matrix);
    free(ba);
    free(i_matrix);
    free(r);
    free(buffer);

    return inv;
}

void print_matrix(float* result, size_t height, size_t width)
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
    omp_set_num_threads(NUM_THREADS);

#pragma omp parallel
    {
        printf("I'm process\n");
    }
    struct timespec start, end;

    float *matrix, *reversed_matrix, *proof;
    proof = (float*)calloc(N * N, sizeof(float));
    matrix = (float*)calloc(N * N, sizeof(float));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i != j)
                matrix[i * N + j] = 2;
            else
                matrix[i * N + j] = 25;
        }
    }

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    reversed_matrix = reverse_matrix(matrix, N);
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);

    printf("Total time taken: %f\n", end.tv_sec - start.tv_sec + 0.000000001 * (end.tv_nsec - start.tv_nsec));

    free(matrix);
    free(reversed_matrix);
    return 0;
}
