#include "mpi.h"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>

#define HEIGHT_1 4
#define WIDTH_HEIGHT 4
#define WIDTH_2 4
#define MAX_MATRIX_EL 100

int mul_matrix_local(const double* first, const double* second, double* result, int height1, int width_height, int height2)
{
    int tmp = 0;
    for (int i = 0; i < height1; i++) {
        for (int k = 0; k < width_height; k++) {
            for (int j = 0; j < height2; j++) {
                result[i * height2 + j] += first[i * width_height + k] * second[k * height2 + j];
                tmp++;
            }
        }
    }
    return tmp;
}

double* create_matrix(const int height, const int width, const bool fill)
{
    double* tmp = (double*)calloc(height * width, sizeof(double));
    if (fill) {
        for (int i = 0; i < height * width; i++) {
            tmp[i] = rand() % MAX_MATRIX_EL;
        }
    }
    return tmp;
}

void print_matrix(const double* matrix, const int height, const int width)
{
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            printf("%lf", matrix[i * width + j]);
        }
        printf("\n");
    }
}

int main(int argc, char* argv[])
{
    struct timespec start, end;
    double total_time;

    int rank, new_rank, size;
    int dims[2], periods[2];

    MPI_Comm col_comm, row_comm, new_comm;
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int proc_grid_width = 4;
    const int proc_grid_height = 4;

    double* first;
    double* second;
    double* result;

    if (rank == 0) {
        first = create_matrix(HEIGHT_1, WIDTH_HEIGHT, true);
        second = create_matrix(WIDTH_HEIGHT, WIDTH_2, true);
        result = create_matrix(HEIGHT_1, WIDTH_2, false);
    }

    int ndims = 2;
    dims[0] = proc_grid_width;
    dims[1] = proc_grid_height;

    periods[0] = 0;
    periods[1] = 0;

    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, 1, &new_comm);
    int coords[2];

    MPI_Cart_coords(new_comm, rank, 2, coords);
    MPI_Comm_split(new_comm, coords[1], coords[0], &col_comm);
    MPI_Comm_split(new_comm, coords[0], coords[1], &row_comm);

    int num_rows_for_proc = HEIGHT_1 / dims[0];
    int num_cols_for_proc = WIDTH_2 / dims[1];

    double* rows_first = create_matrix(WIDTH_HEIGHT, num_rows_for_proc, false);
    double* columns_second = create_matrix(num_cols_for_proc, WIDTH_HEIGHT, false);
    double* part_result = create_matrix(num_rows_for_proc, num_cols_for_proc, false);

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);

    if (coords[1] == 0) {
        printf("row\n");
        MPI_Scatter(first, num_rows_for_proc * WIDTH_HEIGHT, MPI_DOUBLE, rows_first, num_rows_for_proc * WIDTH_HEIGHT, MPI_DOUBLE, 0, col_comm);
    }

    if (coords[0] == 0) {
        printf("col\n");
        MPI_Datatype column_type;
        MPI_Datatype column_resized_type;

        MPI_Type_vector(WIDTH_HEIGHT, num_cols_for_proc, WIDTH_2, MPI_DOUBLE, &column_type);
        MPI_Type_commit(&column_type);
        MPI_Type_create_resized(column_type, 0, num_cols_for_proc * sizeof(double), &column_resized_type);
        MPI_Type_commit(&column_resized_type);
        MPI_Scatter(second, 1, column_resized_type, columns_second, WIDTH_HEIGHT * num_cols_for_proc, MPI_DOUBLE, 0, row_comm);
        MPI_Type_free(&column_resized_type);
        MPI_Type_free(&column_type);
    }

    MPI_Bcast(rows_first, num_rows_for_proc * WIDTH_HEIGHT, MPI_DOUBLE, 0, row_comm);
    MPI_Bcast(columns_second, num_cols_for_proc * WIDTH_HEIGHT, MPI_DOUBLE, 0, col_comm);

    printf("Proc %d iters: %d\n", rank, mul_matrix_local(rows_first, columns_second, part_result, num_rows_for_proc, WIDTH_HEIGHT, num_cols_for_proc));

    if (rank == 0) {
        for (int i = 0; i < dims[0]; i++) {
            for (int j = 0; j < dims[1]; j++) {

                if (i != 0 || j != 0) {
                    int proc_coordinates[2];
                    proc_coordinates[0] = i;
                    proc_coordinates[1] = j;

                    int sending_process_rank;

                    MPI_Cart_rank(new_comm, proc_coordinates, &sending_process_rank);
                    MPI_Recv(part_result, num_rows_for_proc * num_cols_for_proc, MPI_DOUBLE, sending_process_rank, 0, new_comm, MPI_STATUS_IGNORE);
                }

                int start_x = num_rows_for_proc * i;
                int start_y = num_cols_for_proc * j;

                for (int x = 0; x < num_rows_for_proc; x++) {
                    for (int y = 0; y < num_cols_for_proc; y++) {
                        result[(start_x + x) * WIDTH_2 + (start_y + y)] = part_result[x * num_cols_for_proc + y];
                    }
                }
            }
        }

    } else {
        MPI_Send(part_result, num_rows_for_proc * num_cols_for_proc, MPI_DOUBLE, 0, 0, new_comm);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    if (rank == 0) {
        total_time = end.tv_sec - start.tv_sec + 0.000000001 * (double)(end.tv_nsec - start.tv_nsec);
        printf("Time elapsed: %lf\n", total_time);
        free(first);
        free(second);
        free(result);
    }

    free(rows_first);
    free(columns_second);
    free(part_result);

    MPI_Finalize();

    return 0;
}

// void mul_matrix_matrix(float* first, float* second, float* res, size_t size)
// {

//     {
//         auto sum = (float*)calloc(size, sizeof(float));
//         auto tmp = (float*)calloc(size, sizeof(float));

//         for (size_t i = 0; i < size; i++) {
//             memset(sum, 0, size * sizeof(float));
//             for (int j = 0; j < size; j++) {
//                 __m128 Aij = _mm_set1_ps(first[i * size + j]);

//                 for (size_t k = 0; k < size; k += 4) {
//                     __m128 tmp_vec = _mm_load_ps(second + (j * size + k));
//                     __m128 to_store = _mm_mul_ps(Aij, tmp_vec);
//                     _mm_store_ps(tmp + k, to_store);
//                 }

//                 for (size_t k = 0; k < size; k += 4) {
//                     __m128 a = _mm_load_ps(sum + k);
//                     __m128 b = _mm_load_ps(tmp + k);
//                     __m128 r = _mm_add_ps(a, b);
//                     _mm_store_ps(sum + k, r);
//                 }
//             }
//             memcpy(res + i * size, sum, sizeof(float) * size);
//         }
//         free(sum);
//         free(tmp);
//     }
// }
