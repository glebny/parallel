#ifndef BRUSSELATOR_H
#define BRUSSELATOR_H

#include <stdlib.h>
#include "common.h"

/* Brusselator parameters */
#define a 1.8
#define b 1.0

/* Diffusion coefficients */
#define d 1
#define D 16

/* Computional domain parameters */
#define A 262143.0
#define N 262143
#define h (A / N)

/* Number of steps */
#define S 100000

/* Timestep */
#define dt 0.1

/* Right side functions */
double f(const node_t x);
double g(const node_t x);
double dfdu(const node_t x);
double dgdu(const node_t x);
double dfdv(const node_t x);
double dgdv(const node_t x);

/* Initial conditions */
void init(node_t* u, int sizeInput, int size);

#endif // BRUSSELATOR_H
