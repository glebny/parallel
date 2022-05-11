#ifndef REDUCTION_H
#define REDUCTION_H

#include <stdio.h>
#include <math.h>
//#include <C:\Program Files\mingw64\lib\gcc\x86_64-w64-mingw32\8.1.0\include\omp.h>
#include<C:\msys64\mingw64\lib\gcc\x86_64-w64-mingw32\11.3.0\include\omp.h>
#include <stdbool.h>
#include "common.h"
#include "brusselator.h"

int if_power_of_2(int n);

void prl_reduction(matrix_t* subDia, matrix_t* mainDia, matrix_t* supDia, node_t* F, node_t* u1, int size, int sizeInput);

void reduction(matrix_t* Am, matrix_t* Bm, matrix_t* Cm, node_t* Fm, node_t* u1);

void reduction_method(const node_t* u, node_t* u1, int sizeInput, bool flag, int size);

int sizeInput, sizeN, size;

#endif // REDUCTION_H
