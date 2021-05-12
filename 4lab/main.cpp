#include <iostream>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#define NUM_PARTICLES 1024
#define ACCURACY 0.0000001

double* create_cube(const size_t size_x, const size_t size_y, const size_t size_z)
{
    return (double*)calloc(size_x * size_y * size_z, sizeof(double));
}

double func(const double x, const double y, const double z)
{
    return (x * x + y * y + z * z);
}

double right_side(const double x, const double y, const double z)
{
    return (6.0 - 100000 * func(x, y, z));
}

double recount_val(const int x, const int y, const int z, double* cube_part, const double neigh_dist, const int proc_num, const int rank)
{
    int x_starting_coord = rank * NUM_PARTICLES / proc_num;
    double first_coef = 1.0 / (6 / (neigh_dist * neigh_dist) + 100000);
    double second_coef = ((cube_part[((x - 1) * NUM_PARTICLES * NUM_PARTICLES) + (y * NUM_PARTICLES) + z] + cube_part[((x + 1) * NUM_PARTICLES * NUM_PARTICLES) + (y * NUM_PARTICLES) + z]
                              + cube_part[(x * NUM_PARTICLES * NUM_PARTICLES) + ((y - 1) * NUM_PARTICLES) + z] + cube_part[(x * NUM_PARTICLES * NUM_PARTICLES) + ((y + 1) * NUM_PARTICLES) + z]
                              + cube_part[(x * NUM_PARTICLES * NUM_PARTICLES) + (y * NUM_PARTICLES) + z - 1] + cube_part[(x * NUM_PARTICLES * NUM_PARTICLES) + (y * NUM_PARTICLES) + z + 1])
            / neigh_dist
        + right_side((-1.0 + (x_starting_coord + x)) * neigh_dist, -1.0 + y * neigh_dist, -1.0 + z * neigh_dist));

    return first_coef * second_coef;
}

int main(int argc, char** argv)
{

    MPI_Request receive_request[2];
    MPI_Request send_request[2];

    struct timespec start, end;
    double total_time;

    int local_stop_flag = 0;
    int global_stop_flag = 0;
    int proc_num, rank;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &proc_num); // get number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //get process identifier

    double* area_part = create_cube(NUM_PARTICLES / proc_num + 2, NUM_PARTICLES, NUM_PARTICLES);
    int x_starting_coord = rank * NUM_PARTICLES / proc_num;
    double neigh_dist = 2.0 / NUM_PARTICLES;

    double borders_x[2] = { -1.0, 1.0 };
    double borders_y[2] = { -1.0, 1.0 };
    double borders_z[2] = { -1.0, 1.0 };

    for (int z = 0; z < 2; z++) {
        for (int x = 1; x < NUM_PARTICLES / proc_num + 1; x++) {
            for (int y = 0; y < NUM_PARTICLES; y++) {
                area_part[(x * NUM_PARTICLES * NUM_PARTICLES) + (NUM_PARTICLES * (z % 2)) + (y * NUM_PARTICLES)]
                    = func(-1.0 + neigh_dist * (x_starting_coord + x), -1.0 + neigh_dist * y, borders_z[z]);
            }
        }
    }

    for (int y = 0; y < 2; y++) {
        for (int x = 1; x < NUM_PARTICLES / proc_num + 1; x++) {
            for (int z = 0; z < NUM_PARTICLES; z++) {
                area_part[(x * NUM_PARTICLES * NUM_PARTICLES) + ((NUM_PARTICLES - 1) * (y % 2)) + z]
                    = func(-1.0 + neigh_dist * (x_starting_coord + x), borders_y[y], -1.0 + neigh_dist * z);
            }
        }
    }

    if (rank == 0) {
        for (int z = 0; z < NUM_PARTICLES; z++) {
            for (int y = 0; y < NUM_PARTICLES; y++) {
                area_part[(1 * NUM_PARTICLES * NUM_PARTICLES) + (NUM_PARTICLES * y) + z]
                    = func(-1, -1.0 + neigh_dist * y, -1.0 + neigh_dist * y);
            }
        }
    } else if (rank == proc_num - 1) {
        for (int z = 0; z < NUM_PARTICLES; z++) {
            for (int y = 0; y < NUM_PARTICLES; y++) {
                area_part[((NUM_PARTICLES / proc_num) * NUM_PARTICLES * NUM_PARTICLES) + (NUM_PARTICLES * y) + z]
                    = func(1, -1.0 + neigh_dist * y, -1.0 + neigh_dist * y);
            }
        }
    }

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);

    while (true) {
        double prev_function_estimation = area_part[NUM_PARTICLES / proc_num * NUM_PARTICLES * NUM_PARTICLES + NUM_PARTICLES * (NUM_PARTICLES - 5) + (NUM_PARTICLES - 5)];

        //sending data
        if (rank > 0) {
            MPI_Isend(area_part + (NUM_PARTICLES * NUM_PARTICLES), NUM_PARTICLES * NUM_PARTICLES, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &send_request[0]);
            MPI_Irecv(area_part, NUM_PARTICLES * NUM_PARTICLES, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &receive_request[1]);
        }
        if (rank < proc_num - 1) {
            MPI_Isend(area_part + (proc_num / NUM_PARTICLES * NUM_PARTICLES * NUM_PARTICLES), NUM_PARTICLES * NUM_PARTICLES, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &send_request[1]);
            MPI_Irecv(area_part + (NUM_PARTICLES * NUM_PARTICLES * (proc_num / NUM_PARTICLES + 1)), NUM_PARTICLES * NUM_PARTICLES, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &receive_request[0]);
        }

        //calculating values
        for (int x = (NUM_PARTICLES / proc_num + 2) / 2; x > 1; x--) {
            for (int y = 1; y < NUM_PARTICLES; y++) {
                for (int z = 1; z < NUM_PARTICLES; z++) {
                    area_part[x * NUM_PARTICLES * NUM_PARTICLES + y * NUM_PARTICLES + z] = recount_val(x, y, z, area_part, neigh_dist, proc_num, rank);
                }
            }
        }

        for (int x = (NUM_PARTICLES / proc_num + 2) / 2; x < NUM_PARTICLES / proc_num; x++) {
            for (int y = 1; y < NUM_PARTICLES; y++) {
                for (int z = 1; z < NUM_PARTICLES; z++) {
                    area_part[x * NUM_PARTICLES * NUM_PARTICLES + y * NUM_PARTICLES + z] = recount_val(x, y, z, area_part, neigh_dist, proc_num, rank);
                }
            }
        }

        if (rank != 0) {
            MPI_Wait(&send_request[0], MPI_STATUS_IGNORE);
            MPI_Wait(&receive_request[1], MPI_STATUS_IGNORE);
        }
        if (rank != proc_num - 1) {
            MPI_Wait(&send_request[1], MPI_STATUS_IGNORE);
            MPI_Wait(&receive_request[0], MPI_STATUS_IGNORE);
        }

        if (rank > 0) {
            for (int i = 1; i < NUM_PARTICLES - 1; i++) {
                for (int j = 1; j < NUM_PARTICLES - 1; j++) {
                    area_part[(1 * NUM_PARTICLES * NUM_PARTICLES) + (NUM_PARTICLES * j) + i] = recount_val(1, j, i, area_part, neigh_dist, proc_num, rank);
                }
            }
        }
        if (rank < proc_num - 1) {
            for (int z = 1; z < NUM_PARTICLES - 1; z++) {
                for (int y = 1; y < NUM_PARTICLES - 1; y++) {
                    area_part[((NUM_PARTICLES / proc_num) * NUM_PARTICLES * NUM_PARTICLES) + (NUM_PARTICLES * y) + z] = recount_val(NUM_PARTICLES / proc_num, y, z, area_part, neigh_dist, proc_num, rank);
                }
            }
        }

        double check_value = fabs(area_part[NUM_PARTICLES / proc_num * NUM_PARTICLES * NUM_PARTICLES + NUM_PARTICLES * (NUM_PARTICLES - 5) + (NUM_PARTICLES - 5)] - prev_function_estimation);

        if (check_value < ACCURACY) {
            local_stop_flag = 0;
        } else {
            local_stop_flag = 1;
        }

        MPI_Allreduce(&local_stop_flag, &global_stop_flag, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if (global_stop_flag == 0) {
            break;
        }
    }

    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    if (rank == 0) {
        total_time = end.tv_sec - start.tv_sec + 0.000000001 * (double)(end.tv_nsec - start.tv_nsec);
        printf("Time elapsed: %lf\n", total_time);
    }

    free(area_part);
    MPI_Finalize();
    return 0;
}