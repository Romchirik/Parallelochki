#include <mpi.h>
#include <pthread.h>
#include <stdio.h>
#include <unistd.h>

pthread_mutex_t mutex;

void* job_thread(void* args)
{
    return NULL;
}

void* sender_thread(void* args)
{
    uint64_t t = 0;

    return NULL;
}

int main(int argc, char** argv)
{

    int prov, rank, proc_num;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &prov);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);

    // initializing attributes
    pthread_attr_t attr;

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    // 0 - sender
    // 1 - job thread
    pthread_t job_thr;
    pthread_t sender;

    pthread_create(&job_thr, &attr, job_thread, NULL);
    pthread_create(&sender, &attr, sender_thread, NULL);

    if (pthread_join(sender, NULL) != 0 || pthread_join(job_thr, NULL) != 0) {
        perror("Unable to join threads");
        return 0;
    }

    pthread_mutex_destroy(&mutex);
    MPI_Finalize();
    return 0;
}