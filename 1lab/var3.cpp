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

void print_matrix(double* result, size_t height, size_t width);

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

void sub_vector_vector(const double* first, const double* second, double* result, size_t size, size_t proc_num)
{
    size_t num_lines_to_count = size / proc_num;

    double* local_result = (double*)calloc(num_lines_to_count, sizeof(double));

    for (int i = 0; i < num_lines_to_count; i++) {
        local_result[i] = first[i] - second[i];
    }

    // printf("Local result from %d\n");
    // print_matrix(local_result, 1, num_lines_to_count);
    MPI_Allgather(local_result, num_lines_to_count, MPI_DOUBLE, result, num_lines_to_count, MPI_DOUBLE, MPI_COMM_WORLD);
    free(local_result);
}
void sub_vector_vector_partial(const double* first, const double* second, double* result, size_t size, size_t proc_num)
{
    size_t num_lines_to_count = size / proc_num;

    double* local_result = (double*)calloc(num_lines_to_count, sizeof(double));

    for (int i = 0; i < num_lines_to_count; i++) {
        result[i] = first[i] - second[i];
    }

    // printf("Local result from %d\n");
    // print_matrix(local_result, 1, num_lines_to_count);
    free(local_result);
}

double mul_scalar_vector(const double* first, const double* second, size_t size, size_t proc_num, size_t rank)
{
    double local_result = 0, result = 0;
    size_t num_lines_to_count = size / proc_num;
    size_t shift = num_lines_to_count * rank;

    for (int i = 0; i < size; i++) {
        local_result += first[i] * second[i];
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

void extract_my_part(double* matrix, double* my_part, size_t size, size_t proc_num, size_t rank)
{
    size_t my_part_width = SIZE / proc_num;
    size_t shift = my_part_width * rank;

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < my_part_width; j++) {
            my_part[i * my_part_width + j] = matrix[i * size + shift + j];
        }
    }
    if (DEBUG) {
        sleep(rank);
        printf("Process %d recieved:\n", rank);
        print_matrix(my_part, size, my_part_width);
    }
}

void mul_matrix_vector(const double* matrix, const double* vector, size_t size, double* result, size_t proc_num, size_t rank)
{
    size_t num_lines_to_count = size / proc_num;
    size_t shift = num_lines_to_count * rank;

    double* local_result = (double*)calloc(size, sizeof(double));

    for (int i = 0; i < size; i++) {
        double tmp = 0;
        for (int j = 0; j < num_lines_to_count; j++) {
            tmp += matrix[i * num_lines_to_count + j] * vector[j];
        }
        local_result[i] = tmp;
    }

    MPI_Allreduce(local_result, result, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    free(local_result);
}

bool next_step(double* main_matrix, double* solution, double* right_side, size_t size, double accuracy, size_t proc_num, size_t rank)
{
    size_t num_lines_to_count = SIZE / proc_num;
    size_t shift = num_lines_to_count * rank;

    double* tmp = (double*)calloc(size, sizeof(double));
    mul_matrix_vector(main_matrix, solution, size, tmp, proc_num, rank);
    sub_vector_vector(tmp + shift, right_side, tmp, size, proc_num);

    MPI_Barrier(MPI_COMM_WORLD);
    double b_mod = sqrt(mul_scalar_vector(right_side, right_side, num_lines_to_count, proc_num, rank));
    double ax_b_mod = sqrt(mul_scalar_vector(tmp, tmp, SIZE, proc_num, rank));

    printf("a = %lf\n", ax_b_mod / b_mod);
    if (ax_b_mod / b_mod < accuracy) {
        return false;
    } else {
        return true;
    }
    free(tmp);
    return false;
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
    double* b = create_matrix(SIZE / proc_num, 1);
    double* b_full = create_matrix(SIZE, 1);
    double* y_n = create_matrix(SIZE, 1);
    double* x_n = create_matrix(SIZE / proc_num, 1);
    double* main_matrix = create_matrix(SIZE, SIZE / proc_num); //A
    double* main_matrix_full = create_matrix(SIZE, SIZE);

    if (rank == 0) {
        generate_task(main_matrix_full, b_full, SIZE);
    }

    MPI_Bcast(main_matrix_full, SIZE * SIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(b_full, SIZE / proc_num, MPI_DOUBLE, b, SIZE / proc_num, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    extract_my_part(main_matrix_full, main_matrix, SIZE, proc_num, rank);
    //start processing
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);

    size_t num_lines_to_count = SIZE / proc_num;
    size_t shift = num_lines_to_count * rank;

    do {
        mul_matrix_vector(main_matrix, x_n, SIZE, tmp_vector, proc_num, rank); //count Ax
        if (DEBUG && rank == 0) {
            printf("tmp_vector0:\n");
            print_matrix(tmp_vector, 1, SIZE);
        }

        sub_vector_vector(tmp_vector + shift, b, y_n, SIZE, proc_num); //count Ax - b = y
        if (DEBUG && rank == 0) {
            printf("y_n:\n");
            print_matrix(y_n, 1, SIZE);
        }

        mul_matrix_vector(main_matrix, y_n + shift, SIZE, tmp_vector, proc_num, rank); //count Ay
        if (DEBUG && rank == 0) {
            printf("tmp_vector2:\n");
            print_matrix(tmp_vector, 1, SIZE);
        }

        double up = mul_scalar_vector(y_n, tmp_vector, SIZE, proc_num, rank);
        double down = mul_scalar_vector(tmp_vector, tmp_vector, SIZE, proc_num, rank);
        t_n = up / down;
        if (DEBUG && rank == 0) {
            printf("y * Ay: %lf\n", up);
            printf("yA *Ay: %lf\n", down);
            printf("t_n = %lf\n", t_n);
        }

        mul_vector_number(y_n, tmp_vector, t_n, SIZE);
        if (DEBUG && rank == 0) {
            printf("t_y * t_n:\n");
            print_matrix(tmp_vector, 1, SIZE);
        }

        sub_vector_vector_partial(x_n, tmp_vector + shift, x_n, SIZE, proc_num);
        if (DEBUG && rank == 0) {
            printf("x_n:\n");
            print_matrix(x_n, 1, SIZE);
            printf("===========================================\n");
            sleep(1);
        }
        printf("%lf\n", t_n);
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
    free(b_full);
    free(y_n);
    free(x_n);
    free(main_matrix);
    free(main_matrix_full);

    MPI_Finalize();

    return 0;
}
