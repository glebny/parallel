#include "common.h"

/* Determinant */
double mdet(matrix_t m)
{
    return m.a11 * m.a22 - m.a12 * m.a21;
}

/* Scaling matrix by a */
matrix_t mscale(matrix_t m, double a)
{
    matrix_t r;
    r.a11 = m.a11 * a;
    r.a12 = m.a12 * a;
    r.a21 = m.a21 * a;
    r.a22 = m.a22 * a;
    return r;
}

/* Scaling vector by a */
node_t vscale(node_t x, double a)
{
    node_t r;
    r.u = x.u * a;
    r.v = x.v * a;
    return r;
}

/* Inverse matrix */
matrix_t minv(matrix_t m)
{
    matrix_t r;
    r.a11 = m.a22;
    r.a12 = -m.a12;
    r.a21 = -m.a21;
    r.a22 = m.a11;
    return mscale(r, 1.0 / mdet(m));
}

/* transpose matrix */

matrix_t mT(matrix_t m)
{
    matrix_t r;
    r.a11 = m.a11;
    r.a12 = m.a21;
    r.a21 = m.a12;
    r.a22 = m.a22;
    return r;
}

/* vector x matrix multiplication */
node_t vTMmul(node_t x, matrix_t a)
{
    node_t r;
    r.u = a.a11 * x.u + a.a21 * x.v;
    r.v = a.a12 * x.u + a.a22 * x.v;
    return r;
}

/* matrix x matrix multiplication */
matrix_t mmmul(matrix_t a, matrix_t b)
{
    matrix_t r;
    r.a11 = a.a11 * b.a11 + a.a12 * b.a21;
    r.a12 = a.a11 * b.a12 + a.a12 * b.a22;
    r.a21 = a.a21 * b.a11 + a.a22 * b.a21;
    r.a22 = a.a21 * b.a12 + a.a22 * b.a22;
    return r;
}

/* matrix x vector multiplication */
node_t mvmul(matrix_t a, node_t x)
{
    node_t r;
    r.u = a.a11 * x.u + a.a12 * x.v;
    r.v = a.a21 * x.u + a.a22 * x.v;
    return r;
}

/* matrix + matrix addition */
matrix_t madd(matrix_t a, matrix_t b)
{
    matrix_t c;
    c.a11 = a.a11 + b.a11;
    c.a12 = a.a12 + b.a12;
    c.a21 = a.a21 + b.a21;
    c.a22 = a.a22 + b.a22;
    return c;
}

/* vector + vector addition */
node_t vadd(node_t x, double a, node_t y, double b)
{
    node_t c;
    c.u = x.u * a + y.u * b;
    c.v = x.v * a + y.v * b;
    return c;
}
