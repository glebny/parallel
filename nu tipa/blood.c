#include "blood.h"

/* Right side functions*/
double f(const node_t x)
{
    return alpha * x.u * x.u / (x.u + u0) - k1 * x.u - gamma * x.u * x.v;
}
double g(const node_t x)
{
    return beta * x.u * (1.0 - x.v / C) * (1.0 + x.v * x.v / v0 / v0) - k2 * x.v;
}

/* Derivatives */
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

/* Initial conditions */
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
