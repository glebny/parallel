#ifndef REACTION_DIFFUSION_COMMON_H
#define REACTION_DIFFUSION_COMMON_H

// in this file:
// matrix, vector and basic operations for them
// are defined

typedef struct //2D vector x
{
	double u;
	double v;
} node_t;


/*
 * a b
 * c d
 */
typedef struct //matrix 2 x 2
{
	double a11, a12, a21, a22;
} matrix_t;

double mdet(matrix_t m) //matrix determinant
{
	return m.a11 * m.a22 - m.a12 * m.a21;
}

matrix_t mscale(matrix_t m, double a) //scaling matrix by a
{
	matrix_t r;
	r.a11 = m.a11 * a;
	r.a12 = m.a12 * a;
	r.a21 = m.a21 * a;
	r.a22 = m.a22 * a;
	return r;
}

matrix_t minv(matrix_t m) //inverse matrix for m
{
	matrix_t r;
	r.a11 = m.a22;
	r.a12 = -m.a12;
	r.a21 = -m.a21;
	r.a22 = m.a11;
	return mscale(r, 1.0 / mdet(m));
}

matrix_t mmmul(matrix_t a, matrix_t b) //matrix * matrix multiplication
{
	matrix_t r;
	r.a11 = a.a11 * b.a11 + a.a12 * b.a21;
	r.a12 = a.a11 * b.a12 + a.a12 * b.a22;
	r.a21 = a.a21 * b.a11 + a.a22 * b.a21;
	r.a22 = a.a21 * b.a12 + a.a22 * b.a22;
	return r;
}

node_t mvmul(matrix_t a, node_t x) //matrix x vector multiplication
{
	node_t r;
	r.u = a.a11 * x.u + a.a12 * x.v;
	r.v = a.a21 * x.u + a.a22 * x.v;
	return r;
}

node_t vscale(node_t x, double a) //scaling vector by a
{
	node_t r;
	r.u = x.u * a;
	r.v = x.v * a;
	return r;
}

matrix_t madd(matrix_t a, matrix_t b) //matrix + matrix addition
{
	matrix_t c;
	c.a11 = a.a11 + b.a11;
	c.a12 = a.a12 + b.a12;
	c.a21 = a.a21 + b.a21;
	c.a22 = a.a22 + b.a22;
	return c;
}

node_t vadd(node_t x, double a, node_t y, double b) //ax + by lineal product of vectors
{
	node_t c;
	c.u = x.u * a + y.u * b;
	c.v = x.v * a + y.v * b;
	return c;
}

#endif // REACTION_DIFFUSION_COMMON_H