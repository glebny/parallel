#ifndef BLOOD_H
#define BLOOD_H

#include <stdlib.h>

#include "common.h"

/* Параметры модели свёртывания крови. */
#define alpha 2.0
#define beta 0.0015
#define gamma 5.0
#define u0 2.95
#define v0 0.0525
#define C 5.0
#define k1 0.05
#define k2 0.35
#define ucr (k1 * u0 / (alpha - k1))

/* Коефициенты диффузии. */
#define d 0.0001
#define D d

/* Параметры расчетной области. */
#define A 1.0
#define N 1024
#define h (A / N)

/* Число шагов. */
#define S 50000

/* Шаг по времени. */
#define dt 0.0015


/* Функции правой части. */
double f(const node_t x)
{
	return alpha * x.u * x.u / (x.u + u0) - k1 * x.u - gamma * x.u * x.v;
}

double g(const node_t x)
{
	return beta * x.u * (1.0 - x.v / C) * (1.0 + x.v * x.v / v0 / v0) - k2 * x.v;
}

/* Производные. */
double dfdu(const node_t x)
{
    return 2 * alpha * x.u / (x.u + u0) - alpha * x.u * x.u / (x.u + u0) / (x.u + u0) - k1 - gamma * x.v;
}
double dgdu(const node_t x)
{
    return beta * (1.0 - x.v / C) * (1.0 + x.v * x.v / v0 / v0);
}
double dfdv(const node_t x)
{
    return -gamma * x.u;
}
double dgdv(const node_t x)
{
    return -beta * x.u * (1.0 + x.v * x.v / v0 / v0) / C + 2 * beta * x.u * x.v * (1.0 - x.v / C) / v0 / v0 - k2;
}

/* Начальное условие. */
void init(node_t *u)
{
	int i;
	for (i = 0; i < N; i++) {
		u[i].v = 0.0;

		if (i > 10 && i < 180) {
			u[i].u = 0.8;
		} else {
			u[i].u = 0.0;
		}
	}
}

#endif // BLOOD_H
