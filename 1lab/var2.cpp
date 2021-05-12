#include <cmath>
#include <ctime>
#include <mpi/mpi.h>
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

// double mul_scalar_vector(const double* first, const double* second, size_t size, size_t proc_num, size_t rank)
// {
//     double tmp = 0;
//     for (int i = 0; i < size; i++) {
//         tmp += first[i] * second[i];
//     }
//     return tmp;
// }

double mul_scalar_vector(const double* first, const double* second, size_t size, size_t proc_num, size_t rank)
{
    double local_result = 0, result = 0;
    size_t num_lines_to_count = size / proc_num;
    size_t shift = num_lines_to_count * rank;

    for (int i = 0; i < num_lines_to_count; i++) {
        local_result += first[i + shift] * second[i + shift];
    }
    MPI_Allreduce(&local_result, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return result;
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

void mul_matrix_vector(const double* matrix, const double* vector, size_t size, double* result, size_t proc_num, size_t rank)
{
    double* buffer = (double*)calloc(size, sizeof(double));
    size_t num_lines_to_count = size / proc_num;
    size_t shift = num_lines_to_count * rank;

    for (int i = 0; i < num_lines_to_count; i++) {
        for (int j = 0; j < size; j++) {
            result[i + shift] += matrix[i * size + j] * vector[j];
        }
    }

    MPI_Allgather(result + shift, num_lines_to_count, MPI_DOUBLE, buffer, num_lines_to_count, MPI_DOUBLE, MPI_COMM_WORLD);
    memcpy(result, buffer, size * sizeof(double));
    free(buffer);
}

void copy_column(double* matrix, double* dest, size_t size, size_t proc_num, size_t rank)
{
    size_t dest_width = SIZE / proc_num;
    size_t shift = dest_width * rank;

    printf("%d\n", dest_width);
    printf("%d\n", shift);
    for (int i = 0; i < size; i++) {
        //memcpy(dest + (i * dest_width), matrix + (i * size) + shift, dest_width * sizeof(double));
        for (int j = 0; j < dest_width; j++) {
            dest[i * size + j] = matrix[shift + (i * size) + j];
        }
    }
}

void cut_recv_matrix(double* maitix, double* dest, size_t size, size_t proc_num, size_t rank)
{
    size_t dest_width = SIZE / proc_num;
    size_t shift = dest_width * rank;
    double* buffer = (double*)calloc(size, sizeof(double));
    if (rank == 0) {
        copy_column(maitix, dest, size, proc_num, rank);
        printf("Process %d recieved:\n", rank);
        print_matrix(dest, size, dest_width);
        for (int i = 1; i < proc_num; i++) {
            copy_column(maitix, buffer, size, proc_num, i);
            MPI_Send(buffer, size * dest_width, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }

    } else {
        for (int i = 0; i < size; i++) {
            //sleep(rank);
            MPI_Recv(dest, size * dest_width, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("Process %d recieved:\n", rank);
            print_matrix(dest, size, dest_width);
        }
    }
    free(buffer);
    printf("kacnofnvls");
    // MPI_Barrier(MPI_COMM_WORLD);
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

bool next_step(double* main_matrix, double* solution, double* right_side, size_t size, double accuracy, size_t proc_num, size_t rank)
{
    double* tmp = create_matrix(SIZE, 1);
    mul_matrix_vector(main_matrix, solution, size, tmp, proc_num, rank);
    sub_vector_vector(tmp, right_side, tmp, size);

    double b_mod = sqrt(mul_scalar_vector(right_side, right_side, SIZE, proc_num, rank));
    double ax_b_mod = sqrt(mul_scalar_vector(tmp, tmp, SIZE, proc_num, rank));
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

    int proc_num, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num); // get number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //get process identifier

    double* tmp_vector = create_matrix(SIZE, 1);
    double* b = create_matrix(SIZE, 1);
    double* y_n = create_matrix(SIZE, 1);
    double* x_n = create_matrix(SIZE, 1);
    double* main_matrix = create_matrix(SIZE, SIZE / proc_num); //A
    double* main_matrix_full;

    if (rank == 0) {
        main_matrix_full = create_matrix(SIZE, SIZE);
        generate_task(main_matrix_full, b, SIZE);
    }

    size_t lines_to_count = SIZE / proc_num;
    size_t shift = lines_to_count * rank;
    MPI_Scatter(main_matrix_full, lines_to_count * SIZE, MPI_DOUBLE, main_matrix, lines_to_count * SIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   
    //start processing
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);

    do {
        mul_matrix_vector(main_matrix, x_n, SIZE, tmp_vector, proc_num, rank); //count Ax
        if (DEBUG && rank == 0) {
            printf("tmp_vector1:\n");
            print_matrix(tmp_vector, 1, SIZE);
        }

        sub_vector_vector(tmp_vector, b, y_n, SIZE); //count Ax - b = y
        if (DEBUG && rank == 0) {
            printf("tmp_vector:\n");
            print_matrix(tmp_vector, 1, SIZE);
        }

        mul_matrix_vector(main_matrix, y_n, SIZE, tmp_vector, proc_num, rank); //count Ay
        if (DEBUG && rank == 0) {
            printf("tmp_vector2:\n");
            print_matrix(tmp_vector, 1, SIZE);
        }

        t_n = mul_scalar_vector(y_n, tmp_vector, SIZE, proc_num, rank) / mul_scalar_vector(tmp_vector, tmp_vector, SIZE, proc_num, rank);
        if (DEBUG && rank == 0) {
            printf("y * Ay: %lf\n", mul_scalar_vector(y_n, tmp_vector, SIZE, proc_num, rank));
            printf("yA *Ay: %lf\n", mul_scalar_vector(tmp_vector, tmp_vector, SIZE, proc_num, rank));
            printf("t_n = %lf\n", t_n);
        }

        mul_vector_number(y_n, tmp_vector, t_n, SIZE);
        if (DEBUG && rank == 0) {
            printf("t_y * t_n:\n");
            print_matrix(tmp_vector, 1, SIZE);
        }

        sub_vector_vector(x_n, tmp_vector, x_n, SIZE);
        if (DEBUG && rank == 0) {
            printf("x_n:\n");
            print_matrix(x_n, 1, SIZE);
            printf("===========================================\n");
            sleep(1);
        }

    } while (next_step(main_matrix, x_n, b, SIZE, EPS, proc_num, rank));

    //finish processing
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    total_time = end.tv_sec - start.tv_sec + 0.000000001 * (double)(end.tv_nsec - start.tv_nsec);
    if (rank == 0) {
        printf("Time elapsed: %lf\n", total_time);
        if (DEBUG) {
            printf("solution:\n");
            print_matrix(x_n, 1, SIZE);
        }
    }
    free(tmp_vector);
    free(b);
    free(y_n);
    free(x_n);
    free(main_matrix);

    MPI_Finalize();
    return 0;
}
