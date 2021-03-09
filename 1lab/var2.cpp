#include <cmath>
#include <ctime>
#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>

#define MAX_EL_SIZE 50
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

double* create_matrix(size_t width, size_t height)
{
    return (double*)calloc(width * height, sizeof(double));
}

void print_matrix(double* result, size_t height, size_t width)
{
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            printf("%f  ", result[i * width + j]);
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

void copy(const double* source, double* dest, size_t count)
{
    for (int i = 0; i < count; i++) {
        dest[i] = source[i];
    }
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
    copy(buffer, result, size);
    free(buffer);
}

void generate_right_side(double* result, size_t size)
{
    for (int i = 0; i < SIZE; i++) {
        result[i] = rand() % MAX_EL_SIZE;
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
    int proc_num, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num); // get number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //get process identifier

    double t_n = 0;
    double* tmp_vector = create_matrix(SIZE, 1);
    double* b = create_matrix(SIZE, 1);
    double* y_n = create_matrix(SIZE, 1);
    double* x_n = create_matrix(SIZE, 1);
    double* x_n_1 = create_matrix(SIZE, 1);
    double* main_matrix = create_matrix(SIZE, SIZE / proc_num); //A cut
    double* main_matrix_full;

    if (rank == 0) {
        main_matrix_full = create_matrix(SIZE, SIZE); //A
        generate_matrix(main_matrix_full, SIZE);
        generate_right_side(b, SIZE);
    }

    size_t lines_to_count = SIZE / proc_num;
    size_t shift = lines_to_count * rank;

    print_matrix(main_matrix, SIZE, SIZE);
    print_matrix(b, 1, SIZE);

    MPI_Scatter(main_matrix_full, lines_to_count * SIZE, MPI_DOUBLE, main_matrix, lines_to_count * SIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(b, SIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //start processing
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);

    double b_mod = sqrt(mul_scalar_vector(b, b, SIZE));

    do {
        mul_matrix_vector(main_matrix, x_n, SIZE, tmp_vector, proc_num, rank);
        sub_vector_vector(tmp_vector, b, y_n, SIZE);
//        if (rank == 0) {
//            printf("y_n:\n");
//            print_matrix(y_n, 1, SIZE);
//        }

        mul_matrix_vector(main_matrix, y_n, SIZE, tmp_vector, proc_num, rank);
//        if (rank == 0) {
//            printf("tmp_vector:\n");
//            print_matrix(tmp_vector, 1, SIZE);
//        }

        t_n = mul_scalar_vector(y_n, tmp_vector, SIZE) / mul_scalar_vector(tmp_vector, tmp_vector, SIZE);
//        if (rank == 0) {
//            printf("scalar_mul y_n tmp_vector = %lf\n", mul_scalar_vector(y_n, tmp_vector, SIZE));
//            printf("scalar_mul tmp_vector tmp_vector = %lf\n", mul_scalar_vector(tmp_vector, tmp_vector, SIZE));
//            printf("t_n = %lf\n", t_n);
//            print_matrix(tmp_vector, 1, SIZE);
//        }

        mul_vector_number(y_n, tmp_vector, t_n, SIZE);
        sub_vector_vector(x_n, tmp_vector, x_n, SIZE);
        printf("a = %lf\n", sqrt(mul_scalar_vector(y_n, y_n, SIZE)) / b_mod);
//
//        if (rank == 0) {
//            printf("============================================================\n");
//        }
//        sleep(1);
    } while (EPS < (sqrt(mul_scalar_vector(y_n, y_n, SIZE)) / b_mod));

    //finish processing
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    total_time = end.tv_sec - start.tv_sec + 0.000000001 * (double)(end.tv_nsec - start.tv_nsec);

    if(rank == 0){
        printf("Time elapsed: %f\n", total_time);
        print_matrix(x_n, 1, SIZE);
        free(main_matrix_full);

    }


    free(tmp_vector);
    free(b);
    free(y_n);
    free(x_n);
    free(main_matrix);

    MPI_Finalize();
    return 0;
}
