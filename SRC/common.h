#ifndef REACTION_DIFFUSION_COMMON_H
#define REACTION_DIFFUSION_COMMON_H

/* BASIC VECTOR & MATRIX OPERATIONS */

/* VECTOR struct */
typedef struct {
    double u;
    double v;
} node_t;

/* MATRIX struct */
typedef struct {
    double a11, a12, a21, a22;
} matrix_t;

/* Determinant */
double mdet(matrix_t m);

/* Scaling matrix by a */
matrix_t mscale(matrix_t m, double a);

/* Scaling vector by a */
node_t vscale(node_t x, double a);

/* Inverse matrix */
matrix_t minv(matrix_t m);

/* matrix x matrix multiplication */
matrix_t mmmul(matrix_t a, matrix_t b);

/* matrix x vector multiplication */
node_t mvmul(matrix_t a, node_t x);

/* matrix + matrix addition */
matrix_t madd(matrix_t a, matrix_t b);

/* vector + vector addition */
node_t vadd(node_t x, double a, node_t y, double b);

#endif // REACTION_DIFFUSION_COMMON_H
