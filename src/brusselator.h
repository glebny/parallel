#ifndef BRUSSELATOR_H
#define BRUSSELATOR_H

#include <stdlib.h>

#include "common.h"

/* Параметры брюсселятора. */
#define a 1.8
#define b 1.0

/* Коефициенты диффузии. */
#define d 1
#define D 16

/* Параметры расчетной области. */
#define A 262143.0 // 262143.0
#define N 262143 //262143
#define h (A / N)

/* Число шагов. */
#define S 100000

/* Шаг по времени. */
#define dt 0.1


/* Функции правой части. */
double f(const node_t x)
{
    return b + x.u * (x.u * x.v - 1.0 - a);
}

double g(const node_t x)
{
    return x.u * (a - x.u * x.v);
}

/* Производные. */
double dfdu(const node_t x)
{
    return 2.0 * x.u * x.v - (1.0 + a);
}
double dgdu(const node_t x)
{
    return a - 2.0 * x.u * x.v;
}
double dfdv(const node_t x)
{
    return x.u * x.u;
}
double dgdv(const node_t x)
{
    return - x.u * x.u;
}


/* Начальное условие. */
void init(node_t *u, int sizeInput, int size)
{
    int i;
    if (sizeInput != size){
        //#pragma omp parallel for
        for(i=0;i<sizeInput;i++){
            u[i].u = b;
            u[i].v = a / b;

            /* Возмущения. */
            u[i].u *= (1 + 0.01 * (double)rand() / RAND_MAX);
            u[i].v *= (1 + 0.01 * (double)rand() / RAND_MAX);
        }
    //    printf("tut1\n");
        //#pragma omp parallel for
        for (i=sizeInput; i< size; i++){
            u[i].u = 1;
            u[i].v = 1;
        }
    //    printf("tut1\n");
    }

}

#endif // BRUSSELATOR_H
