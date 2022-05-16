#ifndef BLOOD_H
#define BLOOD_H

#include <stdlib.h>
#include "common.h"

/* Parameters of the blood coagulation model */
#define alpha 2.0
#define beta 0.0015
#define gamma 5.0
#define u0 2.95
#define v0 0.0525
#define C 5.0
#define k1 0.05
#define k2 0.35
#define ucr (k1 * u0 / (alpha - k1))

/* Diffusion coefficients */
#define d 0.0001
#define D d

/* Computational area parameters */
#define A 1.0
#define N 1024
#define h (A / N)

/* Number of steps */
#define S 50000

/* Timestep */
#define dt 0.0015


/* Right side functions */
double f(const node_t x);
double g(const node_t x);

/* Derivatives */
double dfdu(const node_t x);
double dgdu(const node_t x);
double dfdv(const node_t x);
double dgdv(const node_t x);

/* Initial conditions */
void init(node_t *u);

#endif // BLOOD_H
