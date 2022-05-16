#include "reduction.h"

int if_power_of_2(int n) {
    n = n + 1;
    while ((n & 1) == 0)
        n >>= 1;
    return (n == 1) ? 1 : 0;
}

void prl_reduction(matrix_t* subDia, matrix_t* mainDia, matrix_t* supDia, node_t* F, node_t* u1, int size, int sizeInput)
{
    int i, j;
    int index1, index2, offset;
    matrix_t alpha, gamma;
    int logSize = log2(size + 1) - 1;
    /// Cyclic Reduction Step 1
    for (i = 0; i < logSize; i++) {
        int step = pow(2, i + 1);

#pragma omp parallel shared (subDia, supDia, mainDia,F, size) private(j,index1, index2, alpha, gamma)
        {
#pragma omp for
            for (j = pow(2, i + 1) - 1; j < size; j = j + step) {
                //offset = pow(2, i);
                index1 = j - pow(2, i);
                index2 = j + pow(2, i);

                alpha = mmmul(subDia[j], minv(mainDia[index1]));
                gamma = mmmul(supDia[j], minv(mainDia[index2]));

                subDia[j] = mmmul(mscale(subDia[index1], -1), alpha);

                mainDia[j] = madd(mainDia[j], madd(mscale(mmmul(supDia[index1], alpha), -1)
                    , mscale(mmmul(subDia[index2], gamma), -1)));
                supDia[j] = mmmul(mscale(supDia[index2], -1), gamma);

                F[j] = vadd(F[j], 1.0, vadd(mvmul(alpha, F[index1]), -1.0, mvmul(gamma, F[index2]), -1.0), 1.0);
            }
        }
    }

    int index = (size - 1) / 2;
//    u1[index] = vTMmul(F[index], minv(mainDia[index]));
    u1[index] = mvmul(minv(mainDia[index]), F[index]);
    for (i = log2(size + 1) - 2; i >= 0; i--) {
        int step = pow(2, i + 1);

#pragma omp parallel shared (u1, subDia, supDia, mainDia,F, size) private(j,index1, index2, alpha, gamma)
        {
#pragma omp for
            for (j = pow(2, i + 1) - 1; j < size; j = j + step) {
                offset = pow(2, i);
                index1 = j - offset;
                index2 = j + offset;

                if (index1 - offset < 0) {

//                    u1[index1] = vadd(F[index1], 1.0, vTMmul(mvmul(supDia[index1], u1[index1 + offset]), minv(mainDia[index1])), -1.0);
                     u1[index1] = vadd(F[index1], 1.0, mvmul(minv(mainDia[index1]), mvmul(supDia[index1], u1[index1 + offset])), -1.0);
                }
                else {

//                    u1[index1] = vadd(vadd(F[index1], 1.0, mvmul(subDia[index1], u1[index1 - offset]), -1.0), 1.0,
//                        vTMmul(mvmul(supDia[index1], u1[index1 + offset]), minv(mainDia[index1])), -1.0);
                    u1[index1] = vadd(vadd(F[index1], 1.0, mvmul(subDia[index1], u1[index1 - offset]), -1.0), 1.0,
                        mvmul(minv(mainDia[index1]), mvmul(supDia[index1], u1[index1 + offset])), -1.0);
                }

                if (index2 + offset >= size) {

//                    u1[index2] = vadd(F[index2], 1.0, vTMmul(mvmul(subDia[index2], u1[index2 - offset]), minv(mainDia[index2])), -1.0);
                    u1[index2] = vadd(F[index2], 1.0, mvmul(minv(mainDia[index2]), mvmul(subDia[index2], u1[index2 - offset])), -1.0);
                }
                else {

//                    u1[index2] = vadd(vadd(F[index2], 1.0, mvmul(subDia[index2], u1[index2 - offset]), -1.0), 1.0,
//                        vTMmul(mvmul(supDia[index2], u1[index2 + offset]), minv(mainDia[index2])), -1.0);
                    u1[index2] = vadd(vadd(F[index2], 1.0, mvmul(subDia[index2], u1[index2 - offset]), -1.0), 1.0,
                        mvmul(minv(mainDia[index2]), mvmul(supDia[index2], u1[index2 + offset])), -1.0);
                }
            }
        }
    }
}

void reduction(matrix_t* Am, matrix_t* Bm, matrix_t* Cm, node_t* Fm, node_t* u1)
{
    int s, j;
    for (s = 1; s < N; s *= 2) {
        matrix_t Cx = mmmul(mscale(Cm[0], -1), minv(Bm[s]));
        Bm[0] = madd(Bm[0], mmmul(Cx, Am[s]));
        Fm[0] = vadd(Fm[0], 1.0, mvmul(Cx, Fm[s]), 1.0);
        Cm[0] = mmmul(Cx, Cm[s]);
        matrix_t Ax = mmmul(mscale(Am[N], -1), minv(Bm[N - s]));
        Bm[N] = madd(Bm[N], mmmul(Ax, Cm[N - s]));
        Fm[N] = vadd(Fm[N], 1.0, mvmul(Ax, Fm[N - s]), 1.0);
        Am[N] = mmmul(Ax, Am[N - s]);
        for (j = 2 * s; j < N; j += 2 * s) {
            Ax = mmmul(mscale(Am[j], -1), minv(Bm[j - s]));
            Cx = mmmul(mscale(Cm[j], -1), minv(Bm[j + s]));
            Bm[j] = madd(Bm[j], mmmul(Ax, Cm[j - s]));
            Bm[j] = madd(Bm[j], mmmul(Cx, Cm[j + s]));
            Fm[j] = vadd(Fm[j], 1.0, mvmul(Ax, Fm[j - s]), 1.0);
            Fm[j] = vadd(Fm[j], 1.0, mvmul(Cx, Fm[j + s]), 1.0);
            Am[j] = mmmul(Ax, Am[j - s]);
            Cm[j] = mmmul(Cx, Cm[j + s]);
        }
    }

    node_t r0 = { .u = Fm[0].u, .v = Fm[0].v };
    node_t rN = { .u = Fm[N].u, .v = Fm[N].v };
    matrix_t Ax = mmmul(mscale(Am[N], -1), minv(Bm[0]));
    matrix_t Cx = mmmul(mscale(Cm[0], -1), minv(Bm[N]));
    r0 = vadd(r0, 1.0, mvmul(Cx, Fm[N]), 1.0);
    rN = vadd(rN, 1.0, mvmul(Ax, Fm[0]), -1.0);
    matrix_t S0 = { .a11 = Bm[0].a11, .a12 = Bm[0].a12, .a21 = Bm[0].a21, .a22 = Bm[0].a22 };
    matrix_t SN = { .a11 = Bm[N].a11, .a12 = Bm[N].a12, .a21 = Bm[N].a21, .a22 = Bm[N].a22 };
    S0 = madd(S0, mmmul(Cx, Am[N]));
    SN = madd(SN, mmmul(Ax, Cm[0]));

    u1[0] = mvmul(minv(S0), r0);
    u1[N] = mvmul(minv(SN), rN);

    for (s = N / 2; s > 0; s /= 2) {
        for (j = s; j < N; j += 2 * s) {
            Fm[j] = vadd(Fm[j], 1.0, mvmul(Am[j], u1[j - s]), -1.0);
            Fm[j] = vadd(Fm[j], 1.0, mvmul(Cm[j], u1[j + s]), -1.0);
            u1[j] = mvmul(minv(Bm[j]), Fm[j]);
        }
    }
}

void reduction_method(const node_t* u, node_t* u1, int sizeInput, bool flag, int size)
{
    int i;

    if ((u == NULL) || (u1 == NULL)) {
        fprintf(stderr, "Invalid arguments for reduction_method\n");
        exit(-1);
    }

    node_t* F = (node_t*)malloc(sizeof(node_t) * (size + 1));
    matrix_t* mainDia = (matrix_t*)malloc(sizeof(matrix_t) * (size + 1));
    matrix_t* subDia = (matrix_t*)malloc(sizeof(matrix_t) * (size + 1));
    matrix_t* supDia = (matrix_t*)malloc(sizeof(matrix_t) * (size + 1));

    if ((F == NULL) || (mainDia == NULL) || (subDia == NULL) || (supDia == NULL)) {
        fprintf(stderr, "Cannot alloc enough memory\n");
        exit(-1);
    }

    if (flag == false) {
        for (i = 0; i < sizeInput; i++) {
            subDia[i] = (matrix_t){ .a11 = dt * d / h / h, .a12 = 0.0, .a21 = 0.0, .a22 = dt * D / h / h };
            mainDia[i] = (matrix_t){ .a11 = -1.0 - 2.0 * dt * d / h / h + dt * dfdu(u[i]),
                    .a12 = dt * dfdv(u[i]),
                    .a21 = dt * dgdu(u[i]),
                    .a22 = -1.0 - 2.0 * dt * D / h / h + dt * dgdv(u[i]) };
            supDia[i] = (matrix_t){ .a11 = dt * d / h / h, .a12 = 0.0, .a21 = 0.0, .a22 = dt * D / h / h };
            F[i] = (node_t){ .u = -u[i].u - dt * f(u[i]) +
                    dt * dfdu(u[i]) * u[i].u + dt * dfdv(u[i]) * u[i].v,
                    .v = -u[i].v - dt * g(u[i]) +
                    dt * dgdu(u[i]) * u[i].u + dt * dgdv(u[i]) * u[i].v };
        }

        for (i = sizeInput; i < size; i++) {
            F[i] = (node_t){ .u = 1.0,
                    .v = 1.0 };
            mainDia[i] = (matrix_t){ .a11 = 1.0,
                    .a12 = 0.0,
                    .a21 = 0.0,
                    .a22 = 1.0 };
            subDia[i] = (matrix_t){ .a11 = 0.0, .a21 = 0.0,
                    .a12 = 0.0, .a22 = 0.0 };
            supDia[i] = (matrix_t){ .a11 = 0.0, .a21 = 0.0,
                    .a12 = 0.0, .a22 = 0.0 };
        }

    }
    else {
        for (i = 0; i < size; i++) {
            subDia[i] = (matrix_t){ .a11 = dt * d / h / h, .a12 = 0.0, .a21 = 0.0, .a22 = dt * D / h / h };
            mainDia[i] = (matrix_t){ .a11 = -1.0 - 2.0 * dt * d / h / h + dt * dfdu(u[i]),
                    .a12 = dt * dfdv(u[i]),
                    .a21 = dt * dgdu(u[i]),
                    .a22 = -1.0 - 2.0 * dt * D / h / h + dt * dgdv(u[i]) };
            supDia[i] = (matrix_t){ .a11 = dt * d / h / h, .a12 = 0.0, .a21 = 0.0, .a22 = dt * D / h / h };
            F[i] = (node_t){ .u = -u[i].u - dt * f(u[i]) +
                    dt * dfdu(u[i]) * u[i].u + dt * dfdv(u[i]) * u[i].v,
                    .v = -u[i].v - dt * g(u[i]) +
                    dt * dgdu(u[i]) * u[i].u + dt * dgdv(u[i]) * u[i].v };
        }
    }
    subDia[0] = (matrix_t){ .a11 = 0.0, .a12 = 0.0, .a21 = 0.0, .a22 = 0.0 };
    mainDia[0] = (matrix_t){ .a11 = -1.0, .a12 = 0.0, .a21 = 0.0, .a22 = -1.0 };
    supDia[0] = (matrix_t){ .a11 = 1.0, .a12 = 0.0, .a21 = 0.0, .a22 = 1.0 };
    F[0] = (node_t){ .u = 0, .v = 0 };
    subDia[sizeInput] = (matrix_t){ .a11 = 1.0, .a12 = 0.0, .a21 = 0.0, .a22 = 1.0 };
    mainDia[sizeInput] = (matrix_t){ .a11 = -1.0, .a12 = 0.0, .a21 = 0.0, .a22 = -1.0 };
    supDia[sizeInput] = (matrix_t){ .a11 = 0.0, .a12 = 0.0, .a21 = 0.0, .a22 = 0.0 };
    F[sizeInput] = (node_t){ .u = 0, .v = 0 };

    prl_reduction(subDia, mainDia, supDia, F, u1, size, sizeInput);

    free(F);
    free(mainDia);
    free(subDia);
    free(supDia);
}

void reduction_method1(const node_t *u, node_t *u1)
{
    int i;
    /* Матрицы метода редукции. */
    matrix_t *Am = (matrix_t*)malloc(sizeof(matrix_t) * (N + 1));
    matrix_t *Bm = (matrix_t*)malloc(sizeof(matrix_t) * (N + 1));
    matrix_t *Cm = (matrix_t*)malloc(sizeof(matrix_t) * (N + 1));
    node_t *Fm   = (node_t*)malloc(sizeof(node_t) * (N + 1));

    if ((u == NULL) || (u1 == NULL)){
        fprintf(stderr, "Invalid arguments for reduction_method\n");
        exit(-1);
    }

    if ((Am == NULL) || (Bm == NULL) || (Cm == NULL) || (Fm == NULL)) {
        fprintf(stderr, "Cannot alloc enough memory\n");
        exit(-1);
    }

    /* Постоянные матрицы. */
    for (i = 1; i < N; i++) {
        Am[i] = (matrix_t){ .a11 = dt * d / h / h, .a12 = 0.0, .a21 = 0.0, .a22 = dt * D / h / h };
        Bm[i] = (matrix_t){ .a11 = -1.0 - 2.0 * dt * d / h / h + dt * dfdu(u[i]),
                    .a12 = dt * dfdv(u[i]),
                    .a21 = dt * dgdu(u[i]),
                    .a22 = -1.0 - 2.0 * dt * D / h / h + dt * dgdv(u[i])};
        Cm[i] = (matrix_t){ .a11 = dt * d / h / h, .a12 = 0.0, .a21 = 0.0, .a22 = dt * D / h / h };
        Fm[i] = (node_t){ .u = -u[i].u - dt * f(u[i]) +
                    dt * dfdu(u[i]) * u[i].u + dt * dfdv(u[i]) * u[i].v,
                  .v = -u[i].v - dt * g(u[i]) +
                    dt * dgdu(u[i]) * u[i].u + dt * dgdv(u[i]) * u[i].v };
    }

    /* Граничные условия. */
    Am[0] = (matrix_t){ .a11 = 0.0, .a12 = 0.0, .a21 = 0.0, .a22 = 0.0 };
    Bm[0] = (matrix_t){ .a11 = -1.0, .a12 = 0.0, .a21 = 0.0, .a22 = -1.0 };
    Cm[0] = (matrix_t){ .a11 = 1.0, .a12 = 0.0, .a21 = 0.0, .a22 = 1.0 };
    Fm[0] = (node_t){ .u = 0, .v = 0 };
    Am[N] = (matrix_t){ .a11 = 1.0, .a12 = 0.0, .a21 = 0.0, .a22 = 1.0 };
    Bm[N] = (matrix_t){ .a11 = -1.0, .a12 = 0.0, .a21 = 0.0, .a22 = -1.0 };
    Cm[N] = (matrix_t){ .a11 = 0.0, .a12 = 0.0, .a21 = 0.0, .a22 = 0.0 };
    Fm[N] = (node_t){ .u = 0, .v = 0 };


    reduction(Am, Bm, Cm, Fm, u1);


    free(Am);
    free(Bm);
    free(Cm);
    free(Fm);
}
