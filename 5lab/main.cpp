#include <math.h>
#include <mpi.h>
#include <pthread.h>
#include <random>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define TASK_LIMITER 5000
#define NUM_TASKS_PACKS 5
#define NUM_TASKS_FOR_PAC 100

int tasks[NUM_TASKS_FOR_PAC];

pthread_mutex_t mutex;
pthread_mutex_t rem_task_mutex;

int tasks_remaining;
int tasks_completed;

double task_res = 0;
int proc_num, rank;

void generateTasks(int* tasksArray, int numTasks, int i)
{
}

void* tasks_executor(void* args)
{

    int providedLevel;
    struct timespec start, end;

    for (int i = 0; i < NUM_TASKS_PACKS; ++i) {

        clock_gettime(CLOCK_MONOTONIC_RAW, &start);

        for (int j = 0; j < NUM_TASKS_FOR_PAC; j++) {
            tasks[j] = rand() % TASK_LIMITER;
        }

        tasks_remaining = NUM_TASKS_FOR_PAC;
        tasks_completed = 0;

        pthread_mutex_lock(&rem_task_mutex);
        int loc_rem_tasks = tasks_remaining;
        pthread_mutex_unlock(&rem_task_mutex);

        for (int j = 0; j < loc_rem_tasks; j++) {
            pthread_mutex_lock(&mutex);
            int repeat = tasks[j];
            pthread_mutex_unlock(&mutex);

            for (int k = 0; k < repeat; k++) {
                task_res += sqrt(k);
            }

            tasks_completed++;
            pthread_mutex_lock(&rem_task_mutex);
            loc_rem_tasks = tasks_remaining - tasks_completed;
            pthread_mutex_unlock(&rem_task_mutex);
        }
        pthread_mutex_lock(&rem_task_mutex);
        tasks_remaining = 0;
        pthread_mutex_unlock(&rem_task_mutex);

        for (int k = 0; k < proc_num; ++k) {
            if (k != rank) {
                int numNewTasks;
                MPI_Send(&rank, 1, MPI_INT, k, 0, MPI_COMM_WORLD);
                MPI_Recv(&numNewTasks, 1, MPI_INT, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if (numNewTasks > 0) {
                    MPI_Recv(tasks, numNewTasks, MPI_INT, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    tasks_remaining = numNewTasks;

                    pthread_mutex_lock(&rem_task_mutex);
                    int currentRemainingTasks = tasks_remaining;
                    pthread_mutex_unlock(&rem_task_mutex);

                    for (int j = 0; j < currentRemainingTasks; j++) {
                        pthread_mutex_lock(&mutex);
                        int repeatNum = tasks[j];
                        pthread_mutex_unlock(&mutex);

                        for (int m = 0; m < repeatNum; m++) {
                            task_res += sin(m);
                        }
                        tasks_completed++;
                        pthread_mutex_lock(&rem_task_mutex);
                        currentRemainingTasks = tasks_remaining;
                        pthread_mutex_unlock(&rem_task_mutex);
                    }
                    pthread_mutex_lock(&rem_task_mutex);
                    tasks_remaining = 0;
                    pthread_mutex_unlock(&rem_task_mutex);
                }
            }
        }
        clock_gettime(CLOCK_MONOTONIC_RAW, &end);
        double my_time = end.tv_sec - start.tv_sec + 0.000000001 * (double)(end.tv_nsec - start.tv_nsec);
        double min_time, max_time;
        MPI_Allreduce(&my_time, &min_time, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&my_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        for (int j = 0; j < proc_num; ++j) {
            MPI_Barrier(MPI_COMM_WORLD);
            if (rank == 0 && j == 0) {
                printf("iterarion %d\n", i);
                printf("disbalace in %% = %f\n", ((max_time - min_time) / max_time) * 100);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    int stopJob = -1;
    MPI_Send(&stopJob, 1, MPI_INT, rank, 0, MPI_COMM_WORLD);

    return NULL;
}

void* communicator(void* args)
{
    while (true) {
        int requesting_proc;
        MPI_Status status;
        MPI_Recv(&requesting_proc, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        if (requesting_proc == -1) {
            pthread_exit(NULL);
        }

        pthread_mutex_lock(&mutex);
        if (tasks_remaining > 10) {
            int shared_tasks = tasks_remaining / 2;
            pthread_mutex_lock(&rem_task_mutex);
            tasks_remaining -= shared_tasks;
            pthread_mutex_unlock(&rem_task_mutex);

            MPI_Send(&shared_tasks, 1, MPI_INT, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
            MPI_Send(&tasks[shared_tasks - 1], shared_tasks, MPI_INT, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
        } else {
            int numShareTasks = 0;
            MPI_Send(&numShareTasks, 1, MPI_INT, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
        }
        pthread_mutex_unlock(&mutex);
    }
    return NULL;
}

int main(int argc, char** argv)
{
    double total_time;
    int providedLevel;
    struct timespec start, end;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &providedLevel);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);

    pthread_mutex_init(&mutex, NULL);
    pthread_mutex_init(&rem_task_mutex, NULL);

    pthread_t commiunicator_thread, task_thread;
    pthread_attr_t communicator_attr, tasks_attr;

    pthread_attr_init(&tasks_attr);
    pthread_attr_init(&communicator_attr);

    pthread_attr_setdetachstate(&communicator_attr, PTHREAD_CREATE_JOINABLE);
    pthread_attr_setdetachstate(&tasks_attr, PTHREAD_CREATE_JOINABLE);

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);

    pthread_create(&commiunicator_thread, &communicator_attr, tasks_executor, NULL);
    pthread_create(&task_thread, &tasks_attr, communicator, NULL);

    pthread_join(commiunicator_thread, NULL);
    pthread_join(task_thread, NULL);

    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    if (rank == 0) {
        total_time = end.tv_sec - start.tv_sec + 0.000000001 * (double)(end.tv_nsec - start.tv_nsec);
        printf("Time elapsed: %lf\n", total_time);
    }

    pthread_attr_destroy(&communicator_attr);
    pthread_attr_destroy(&tasks_attr);

    pthread_mutex_destroy(&mutex);
    pthread_mutex_destroy(&rem_task_mutex);
    MPI_Finalize();
    return 0;
}