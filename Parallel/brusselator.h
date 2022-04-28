#ifndef BRUSSELATOR_H
#define BRUSSELATOR_H

#include <stdlib.h>
#include "common.h"

// in this file:
// parameters and functions for brusselator model
// are defined

/* Параметры брюсселятора. */
#define a 1.8
#define b 1.0

/* Коэффициенты диффузии. */
#define d 1
#define D 16

/* Параметры расчетной области. */
#define A 256.0
#define N 256
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
void init(node_t *u)
{
    int i;
    for (i = 0; i < N; i++) {
        u[i].u = b;
        u[i].v = a / b;

        /* Возмущения. */
        u[i].u *= (1 + 0.01 * (double)rand() / RAND_MAX);
        u[i].v *= (1 + 0.01 * (double)rand() / RAND_MAX);
    }
}

#endif // BRUSSELATOR_H
