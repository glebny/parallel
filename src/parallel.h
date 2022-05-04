#ifndef PARALLEL_H
#define PARALLEL_H

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include "brusselator.h"

void check_success(int code)
{
    if (code != MPI_SUCCESS)
    {
        printf("Error.\n");
        MPI_Abort(MPI_COMM_WORLD, code);
    }
}

int rnk, sz; //rank, size
MPI_Status st;

//initializing MPI
int MPI_Part(int argc, char* argv[])
{
    int res = MPI_Init(&argc, &argv);
    check_success(res);

    res = MPI_Comm_size(MPI_COMM_WORLD, &sz);
    check_success(res);

    res = MPI_Comm_rank(MPI_COMM_WORLD, &rnk);
    check_success(res);
}

#define p sz //p = number of processes, as it said in Main.cpp (comments on ~57-58 strings)
#define N_prl N / p //N == 256 from brusselator.h (?)
#define rank rnk
#define MASTER 1


#endif //PARALLEL_H