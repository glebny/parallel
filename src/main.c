#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <C:\Program Files\mingw64\lib\gcc\x86_64-w64-mingw32\8.1.0\include\omp.h>
#include <time.h>
#include "vtk.h"
#include "common.h"
#include <stdbool.h>
//#include "parallel.h"

/*
#ifdef BLOOD
#include "blood.h"
#endif
#ifdef BRUSSELATOR
#include "brusselator.h"
#endif
*/

#include "brusselator.h" //выберем пока один для определенности
//#include "blood.h"
int if_power_of_2(int n) {
    n = n + 1;
    while ((n & 1) == 0)
        n >>= 1;
    return (n == 1) ? 1 : 0;
}
//void prl_reduction(matrix_t* Am, matrix_t* Bm, matrix_t* Cm, node_t* Fm, node_t* u1) //судя по всему речь о массивах матриц, триндец
//{
//    int s, j;
//    // for each proccesses
//    for (s = 1; s < N_prl; s *= 2) {
//        // define boundary values
//        // for 0 and last 0 and N
//        // for others mpi_send for right boundary
//        matrix_t Cx = mmmul(mscale(Cm[0], -1), minv(Bm[s]));
//        Bm[0] = madd(Bm[0], mmmul(Cx, Am[s]));
//        Fm[0] = vadd(Fm[0], 1.0, mvmul(Cx, Fm[s]), 1.0);
//        Cm[0] = mmmul(Cx, Cm[s]);
//        matrix_t Ax = mmmul(mscale(Am[N], -1), minv(Bm[N - s]));
//        Bm[N] = madd(Bm[N], mmmul(Ax, Cm[N - s]));
//        Fm[N] = vadd(Fm[N], 1.0, mvmul(Ax, Fm[N - s]), 1.0);
//        Am[N] = mmmul(Ax, Am[N - s]);
//        for (j = 2 * s; j < N; j += 2 * s) {
//            Ax = mmmul(mscale(Am[j], -1), minv(Bm[j - s]));
//            Cx = mmmul(mscale(Cm[j], -1), minv(Bm[j + s]));
//            Bm[j] = madd(Bm[j], mmmul(Ax, Cm[j - s]));
//            Bm[j] = madd(Bm[j], mmmul(Cx, Cm[j + s]));
//            Fm[j] = vadd(Fm[j], 1.0, mvmul(Ax, Fm[j - s]), 1.0);
//            Fm[j] = vadd(Fm[j], 1.0, mvmul(Cx, Fm[j + s]), 1.0);
//            Am[j] = mmmul(Ax, Am[j - s]);
//            Cm[j] = mmmul(Cx, Cm[j + s]);
//        }
//    }


//    // Gather on Master all last points
//    // and call reduction (Am_last, Bm_last, ... , u[0]);

///*добавлено мной*/
//    matrix_t Am_last = Am[N];
//    matrix_t Bm_last = Bm[N];
//    matrix_t Cm_last = Cm[N];
//    node_t Fm_last = Fm[N];
///*добавлено мной*/

//    if (rank == MASTER) {
//        //áûëî reduction(Am_last, Bm_last, Cm_last, Fm_last, u1[0]);
//        reduction (&Am_last, &Bm_last, &Cm_last, &Fm_last, &u1[0]);
//    }
//    // and Scatter u[p+1] on all proccesses

//    // rebuild the solution
//    // N_prl = N / p (p - number of proccesses)
//    // for master (rank = 0) N_prl + 1 (+ 0)

//    for (s = N_prl / 2; s > 0; s /= 2) {
//        for (j = s; j < N_prl; j += 2 * s) {
//            Fm[j] = vadd(Fm[j], 1.0, mvmul(Am[j], u1[j - s]), -1.0);
//            Fm[j] = vadd(Fm[j], 1.0, mvmul(Cm[j], u1[j + s]), -1.0);
//            u1[j] = mvmul(minv(Bm[j]), Fm[j]);
//        }
//    }
//}

void prl_reduction(matrix_t* subDia, matrix_t* mainDia, matrix_t* supDia, node_t* F, node_t* u1, int size, int sizeInput)
{
//    clock_t start = clock();
    int i,j;
    int index1,index2,offset;
    matrix_t alpha,gamma;
    int logSize = log2(size+1)-1;
    /// Cyclic Reduction Step 1
    //printf("tut1\n");
    for(i=0;i<logSize;i++){
        int step = pow(2,i+1);
#pragma omp parallel shared (subDia, supDia, mainDia,F, size) private(j,index1, index2, alpha, gamma)
        {
#pragma omp for
            //printf("tut1\n");
            for(j=pow(2,i+1)-1;j<size;j=j+ step){
                offset = pow(2,i);
                index1 = j - offset;
                index2 = j + offset;

                alpha = mmmul(subDia[j],minv(mainDia[index1]));
                gamma = mmmul(supDia[j],minv(mainDia[index2]));
               // printf("tut1\n");
                //#pragma omp atomic capture
                subDia[j] = mmmul(mscale(subDia[index1], -1), alpha);

                mainDia[j] = madd( mainDia[j] , madd( mscale ( mmmul( supDia[index1] , alpha), -1)
                                                      , mscale( mmmul( subDia[index2] , gamma),-1)));
                supDia[j] = mmmul(mscale(supDia[index2], -1), gamma);

                F[j] = vadd(F[j],1.0,vadd(mvmul(alpha,F[index1]),-1.0,mvmul(gamma,F[index2]),-1.0),1.0);
            }

        }

    }

    int index = (size - 1)/2;
    u1[index] = mvmul(minv(mainDia[index]),F[index]);
    //printf("tut1\n");
    for(i=log2(size+1)-2;i>=0;i--){
        int step = pow(2,i+1);
#pragma omp parallel shared(u1,F,subDia, supDia, mainDia, size) private(j,index1, index2, alpha, gamma)

        {
#pragma omp  for
            for(j=pow(2,i+1)-1;j<size;j=j+ step){
                offset = pow(2,i);
                index1 = j - offset;
                index2 = j + offset;

                //printf("Executed by %d \n", omp_get_thread_num());
                if (index1 - offset < 0){

                    u1[index1] =
                            vadd(F[index1],1.0,mvmul(minv(mainDia[index1]),mvmul(supDia[index1],u1[index1+offset])) ,-1.0);
                }
                else{
                    u1[index1] = vadd(vadd(F[index1],1.0,mvmul(subDia[index1],u1[index1-offset]),-1.0),1.0,
                            mvmul(minv(mainDia[index1]),mvmul(supDia[index1],u1[index1+offset])),-1.0);

                }

                if(index2 + offset >= size ){
                    u1[index2] = vadd(F[index2],1.0,mvmul(minv(mainDia[index2]),mvmul(subDia[index2],u1[index2-offset])),-1.0);
                }
                else{
                    u1[index2] = vadd(vadd(F[index2],1.0,mvmul(subDia[index2],u1[index2-offset]),-1.0),1.0,
                            mvmul(minv(mainDia[index2]),mvmul(supDia[index2],u1[index2+offset])),-1.0);

                }
            }

        }
    //printf("tut1\n");

    }



//    for(i=0;i<sizeInput;i++){
//        cout << x[i] << endl;
//    }

    // Stop measuring time and calculate the elapsed time
//    clock_t end = clock();
//    double elapsed = (double)(end - start)/CLOCKS_PER_SEC;

    //    printf("Time measured: %.3f seconds.\n", elapsed);
//    cout << "The Time measured is : " << elapsed << endl;
//    return 0;  // return 0 to the OS.
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

    /*
     * Now the system is reduced to
     * B0 u0 + C0 uN = d0
     * AN u0 + BN uN = dN
     * B0 u0 + C0 uN = d0
     * C0 (BN \ AN) u0 + C0 uN = C0 (BN \ dN)
     * u0 = [B0 - C0 (BN \ AN)] \ [d0 - C0 (BN \ dN)]
     * uN = [BN - AN (B0 \ C0)] \ [dN - AN (B0 \ d0)]
     */
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
    /* Матрицы метода редукции. */
//    matrix_t* Am = (matrix_t*)malloc(sizeof(matrix_t) * (N + 1));
//    matrix_t* Bm = (matrix_t*)malloc(sizeof(matrix_t) * (N + 1));
//    matrix_t* Cm = (matrix_t*)malloc(sizeof(matrix_t) * (N + 1));
//    node_t* Fm = (node_t*)malloc(sizeof(node_t) * (N + 1));

    if ((u == NULL) || (u1 == NULL)) {
        fprintf(stderr, "Invalid arguments for reduction_method\n");
        exit(-1);
    }

//    if ((Am == NULL) || (Bm == NULL) || (Cm == NULL) || (Fm == NULL)) {
//        fprintf(stderr, "Cannot alloc enough memory\n");
//        exit(-1);
//    }

//    /* Постоянные матрицы. */
//    for (i = 1; i < N; i++) {
//        Am[i] = (matrix_t){ .a11 = dt * d / h / h, .a12 = 0.0, .a21 = 0.0, .a22 = dt * D / h / h };
//        Bm[i] = (matrix_t){ .a11 = -1.0 - 2.0 * dt * d / h / h + dt * dfdu(u[i]),
//                .a12 = dt * dfdv(u[i]),
//                .a21 = dt * dgdu(u[i]),
//                .a22 = -1.0 - 2.0 * dt * D / h / h + dt * dgdv(u[i]) };
//        Cm[i] = (matrix_t){ .a11 = dt * d / h / h, .a12 = 0.0, .a21 = 0.0, .a22 = dt * D / h / h };
//        Fm[i] = (node_t){ .u = -u[i].u - dt * f(u[i]) +
//                dt * dfdu(u[i]) * u[i].u + dt * dfdv(u[i]) * u[i].v,
//                .v = -u[i].v - dt * g(u[i]) +
//                dt * dgdu(u[i]) * u[i].u + dt * dgdv(u[i]) * u[i].v };
//    }

//    /* Граничные условия. */
//    Am[0] = (matrix_t){ .a11 = 0.0, .a12 = 0.0, .a21 = 0.0, .a22 = 0.0 };
//    Bm[0] = (matrix_t){ .a11 = -1.0, .a12 = 0.0, .a21 = 0.0, .a22 = -1.0 };
//    Cm[0] = (matrix_t){ .a11 = 1.0, .a12 = 0.0, .a21 = 0.0, .a22 = 1.0 };
//    Fm[0] = (node_t){ .u = 0, .v = 0 };
//    Am[N] = (matrix_t){ .a11 = 1.0, .a12 = 0.0, .a21 = 0.0, .a22 = 1.0 };
//    Bm[N] = (matrix_t){ .a11 = -1.0, .a12 = 0.0, .a21 = 0.0, .a22 = -1.0 };
//    Cm[N] = (matrix_t){ .a11 = 0.0, .a12 = 0.0, .a21 = 0.0, .a22 = 0.0 };
//    Fm[N] = (node_t){ .u = 0, .v = 0 };
    //
  //  printf("tut1\n");
    node_t* F = (node_t*)malloc(sizeof(node_t) * (size+1));
    matrix_t* mainDia = (matrix_t*)malloc(sizeof(matrix_t) * (size+1));
    matrix_t* subDia = (matrix_t*)malloc(sizeof(matrix_t) * (size+1));
    matrix_t* supDia = (matrix_t*)malloc(sizeof(matrix_t) * (size+1));

    if ((F == NULL) || (mainDia == NULL) || (subDia == NULL) || (supDia == NULL)) {
        fprintf(stderr, "Cannot alloc enough memory\n");
        exit(-1);
    }
   // printf("tut2\n, flag: %d\n",flag);
    if (flag == false){
        //#pragma omp parallel for
        for(i=0;i<sizeInput;i++){
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
    //    printf("tut1\n");
        //#pragma omp parallel for
        for (i=sizeInput; i< size; i++){
            F[i] =  (node_t){ .u = 1.0,
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
    //    printf("tut1\n");
    }
    else{
        for(i=0;i<size;i++){
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
   // printf("tut1\n");
    subDia[0] = (matrix_t){ .a11 = 0.0, .a12 = 0.0, .a21 = 0.0, .a22 = 0.0 };
  //  printf("%f, %f, %f %f  A%d\n",subDia[0].a11, subDia[0].a12, subDia[0].a21, subDia[0].a22,i);
    mainDia[0] = (matrix_t){ .a11 = -1.0, .a12 = 0.0, .a21 = 0.0, .a22 = -1.0 };
    supDia[0] = (matrix_t){ .a11 = 1.0, .a12 = 0.0, .a21 = 0.0, .a22 = 1.0 };
    F[0] = (node_t){ .u = 0, .v = 0 };
    subDia[sizeInput-1] = (matrix_t){ .a11 = 1.0, .a12 = 0.0, .a21 = 0.0, .a22 = 1.0 };
    mainDia[sizeInput-1] = (matrix_t){ .a11 = -1.0, .a12 = 0.0, .a21 = 0.0, .a22 = -1.0 };
    supDia[sizeInput-1] = (matrix_t){ .a11 = 0.0, .a12 = 0.0, .a21 = 0.0, .a22 = 0.0 };
    F[sizeInput-1] = (node_t){ .u = 0, .v = 0 };
   // printf("tut1\n");
    //
    //reduction(Am, Bm, Cm, Fm, u1);
//    printf("tut peded reduct\n");
//    prl_reduction(subDia, mainDia, supDia, F, u1, size, sizeInput);
//    for(i=0;i<size;i++){
//            printf("%f, %f,  F%d\n",F[i].u, F[i].v,i);
//    }
//    for(i=0;i<size;i++){
//            printf("%f, %f,  U%d\n",u1[i].u, u1[i].v,i);
//    }
//    for(i=0;i<size;i++){
//            printf("%f, %f, %f %f  A%d\n",subDia[i].a11, subDia[i].a12, subDia[i].a21, subDia[i].a22,i);
//    }
//    for(i=0;i<size;i++){
//            printf("%f, %f, %f %f B%d\n", mainDia[i].a11, mainDia[i].a12, mainDia[i].a21, mainDia[i].a22,i);
//    }
//    for(i=0;i<size;i++){
//            printf("%f, %f, %f %f  C%d\n",supDia[i].a11, supDia[i].a12, supDia[i].a21, supDia[i].a22,i);
//    }
    prl_reduction(subDia, mainDia, supDia, F, u1, size, sizeInput);
//    printf("tut posle reduct\n");
//    for(i=0;i<size;i++){
//            printf("%f, %f,  F%d\n",F[i].u, F[i].v,i);
//    }
//    for(i=0;i<size;i++){
//            printf("%f, %f,  U%d\n",u1[i].u, u1[i].v,i);
//    }
//    for(i=0;i<size;i++){
//            printf("%f, %f, %f %f  A%d\n",subDia[i].a11, subDia[i].a12, subDia[i].a21, subDia[i].a22,i);
//    }
//    for(i=0;i<size;i++){
//            printf("%f, %f, %f %f B%d\n", mainDia[i].a11, mainDia[i].a12, mainDia[i].a21, mainDia[i].a22,i);
//    }
//    for(i=0;i<size;i++){
//            printf("%f, %f, %f %f  C%d\n",supDia[i].a11, supDia[i].a12, supDia[i].a21, supDia[i].a22,i);
//    }

    free(F);
  //  printf("tut posle reduct\n");
    free(mainDia);
  //  printf("tut posle reduct\n");
    free(subDia);
 //   printf("tut posle reduct\n");
    free(supDia);
   // printf("tut posle reduct\n");
    //
}

//
int sizeInput, sizeN, size;
//

int main(int argc, char** argv)
{
    //
    int count =0;
    bool flag;
    sizeInput = N;
    int maxMatsize = 25;
    double * p = (double*)malloc(sizeof(double) * (maxMatsize));

    for (int k = 0; k<maxMatsize; k++ ){

        p[k] = pow(2,k)-1;

        if(p[k] == sizeInput){
            size = p[k];
            flag =  true;

            // cout << sizeN;
            // cout << endl;
        }
        else {
            if (p[k]>= sizeInput ){
                count = count +1;
            }
        }
    }
    sizeN = p[maxMatsize-count];
    if (flag ==  false){
        size = sizeN;
    }

    clock_t start = clock();

    int i,j;
    int index1,index2,offset;
    double alpha,gamma;

    //
    node_t* u, * u1;
    node_t* temp;
    //int i;
    char buf[256];
    const char* save[2] = { "u", "v" };

    u = (node_t*)malloc(sizeof(node_t) * (size + 1));
    u1 = (node_t*)malloc(sizeof(node_t) * (size + 1));
    if ((u == NULL) || (u1 == NULL)) {
        fprintf(stderr, "Cannot alloc enough memory\n");
        exit(-1);
    }

//    if (!if_power_of_2(N)) {
//        fprintf(stderr, "Parameter N must be 2^k-1\n");
//        exit(-1);
//    }

    init(u,sizeInput, size);
    init(u1,sizeInput, size);

    for (i = 0; i < 10; i++) { //10
        /* Сохраняем посчитанные значения. */
        if (i % 2 == 0) {
            sprintf(buf, "C:/Users/Gleb/Desktop/lobanov/S_reaction_diffusion_matrix/src/res1/data_%06d.vtk", i);
            write_to_vtk2((double*)u, buf, save, N, 0.0, h, 2);
        }

        /* Обновляем значение. */
        reduction_method(u, u1,sizeInput,flag, size);
//            for(i=0;i<size;i++){
//                    printf("%f, %f,  F%d\n",u1[i].u, u1[i].v,i);
//            }

     //   printf("tut1\n");
        /* Меняем местами временной слой. */
        temp = u;
        u = u1;
        u1 = temp;
    }

    free(u);
    free(u1);
    clock_t end = clock();
        double elapsed = (double)(end - start)/CLOCKS_PER_SEC;

        printf("Time measured: %.3f seconds.\n", elapsed);
    return 0;
}

//#define _GNU_SOURCE
//#include <stdio.h>
//#include <stdlib.h>

//#include "vtk.h"
//#include "common.h"
////#include "parallel.h"

///*
//#ifdef BLOOD
//#include "blood.h"
//#endif
//#ifdef BRUSSELATOR
//#include "brusselator.h"
//#endif
//*/

//#include "brusselator.h" //выберем пока один для определенности
////#include "blood.h"
//int if_power_of_2(int n) {
//    while ((n & 1) == 0)
//        n >>= 1;
//    return (n == 1) ? 1 : 0;
//}
////void prl_reduction(matrix_t* Am, matrix_t* Bm, matrix_t* Cm, node_t* Fm, node_t* u1) //судя по всему речь о массивах матриц, триндец
////{
////    int s, j;
////    // for each proccesses
////    for (s = 1; s < N_prl; s *= 2) {
////        // define boundary values
////        // for 0 and last 0 and N
////        // for others mpi_send for right boundary
////        matrix_t Cx = mmmul(mscale(Cm[0], -1), minv(Bm[s]));
////        Bm[0] = madd(Bm[0], mmmul(Cx, Am[s]));
////        Fm[0] = vadd(Fm[0], 1.0, mvmul(Cx, Fm[s]), 1.0);
////        Cm[0] = mmmul(Cx, Cm[s]);
////        matrix_t Ax = mmmul(mscale(Am[N], -1), minv(Bm[N - s]));
////        Bm[N] = madd(Bm[N], mmmul(Ax, Cm[N - s]));
////        Fm[N] = vadd(Fm[N], 1.0, mvmul(Ax, Fm[N - s]), 1.0);
////        Am[N] = mmmul(Ax, Am[N - s]);
////        for (j = 2 * s; j < N; j += 2 * s) {
////            Ax = mmmul(mscale(Am[j], -1), minv(Bm[j - s]));
////            Cx = mmmul(mscale(Cm[j], -1), minv(Bm[j + s]));
////            Bm[j] = madd(Bm[j], mmmul(Ax, Cm[j - s]));
////            Bm[j] = madd(Bm[j], mmmul(Cx, Cm[j + s]));
////            Fm[j] = vadd(Fm[j], 1.0, mvmul(Ax, Fm[j - s]), 1.0);
////            Fm[j] = vadd(Fm[j], 1.0, mvmul(Cx, Fm[j + s]), 1.0);
////            Am[j] = mmmul(Ax, Am[j - s]);
////            Cm[j] = mmmul(Cx, Cm[j + s]);
////        }
////    }


////    // Gather on Master all last points
////    // and call reduction (Am_last, Bm_last, ... , u[0]);

/////*добавлено мной*/
////    matrix_t Am_last = Am[N];
////    matrix_t Bm_last = Bm[N];
////    matrix_t Cm_last = Cm[N];
////    node_t Fm_last = Fm[N];
/////*добавлено мной*/

////    if (rank == MASTER) {
////        //áûëî reduction(Am_last, Bm_last, Cm_last, Fm_last, u1[0]);
////        reduction (&Am_last, &Bm_last, &Cm_last, &Fm_last, &u1[0]);
////    }
////    // and Scatter u[p+1] on all proccesses

////    // rebuild the solution
////    // N_prl = N / p (p - number of proccesses)
////    // for master (rank = 0) N_prl + 1 (+ 0)

////    for (s = N_prl / 2; s > 0; s /= 2) {
////        for (j = s; j < N_prl; j += 2 * s) {
////            Fm[j] = vadd(Fm[j], 1.0, mvmul(Am[j], u1[j - s]), -1.0);
////            Fm[j] = vadd(Fm[j], 1.0, mvmul(Cm[j], u1[j + s]), -1.0);
////            u1[j] = mvmul(minv(Bm[j]), Fm[j]);
////        }
////    }
////}

//void reduction(matrix_t* Am, matrix_t* Bm, matrix_t* Cm, node_t* Fm, node_t* u1)
//{
//    int s, j;
//    for (s = 1; s < N; s *= 2) {
//        matrix_t Cx = mmmul(mscale(Cm[0], -1), minv(Bm[s]));
//        Bm[0] = madd(Bm[0], mmmul(Cx, Am[s]));
//        Fm[0] = vadd(Fm[0], 1.0, mvmul(Cx, Fm[s]), 1.0);
//        Cm[0] = mmmul(Cx, Cm[s]);
//        matrix_t Ax = mmmul(mscale(Am[N], -1), minv(Bm[N - s]));
//        Bm[N] = madd(Bm[N], mmmul(Ax, Cm[N - s]));
//        Fm[N] = vadd(Fm[N], 1.0, mvmul(Ax, Fm[N - s]), 1.0);
//        Am[N] = mmmul(Ax, Am[N - s]);
//        for (j = 2 * s; j < N; j += 2 * s) {
//            Ax = mmmul(mscale(Am[j], -1), minv(Bm[j - s]));
//            Cx = mmmul(mscale(Cm[j], -1), minv(Bm[j + s]));
//            Bm[j] = madd(Bm[j], mmmul(Ax, Cm[j - s]));
//            Bm[j] = madd(Bm[j], mmmul(Cx, Cm[j + s]));
//            Fm[j] = vadd(Fm[j], 1.0, mvmul(Ax, Fm[j - s]), 1.0);
//            Fm[j] = vadd(Fm[j], 1.0, mvmul(Cx, Fm[j + s]), 1.0);
//            Am[j] = mmmul(Ax, Am[j - s]);
//            Cm[j] = mmmul(Cx, Cm[j + s]);
//        }
//    }

//    /*
//     * Now the system is reduced to
//     * B0 u0 + C0 uN = d0
//     * AN u0 + BN uN = dN
//     * B0 u0 + C0 uN = d0
//     * C0 (BN \ AN) u0 + C0 uN = C0 (BN \ dN)
//     * u0 = [B0 - C0 (BN \ AN)] \ [d0 - C0 (BN \ dN)]
//     * uN = [BN - AN (B0 \ C0)] \ [dN - AN (B0 \ d0)]
//     */
//    node_t r0 = { .u = Fm[0].u, .v = Fm[0].v };
//    node_t rN = { .u = Fm[N].u, .v = Fm[N].v };
//    matrix_t Ax = mmmul(mscale(Am[N], -1), minv(Bm[0]));
//    matrix_t Cx = mmmul(mscale(Cm[0], -1), minv(Bm[N]));
//    r0 = vadd(r0, 1.0, mvmul(Cx, Fm[N]), 1.0);
//    rN = vadd(rN, 1.0, mvmul(Ax, Fm[0]), -1.0);
//    matrix_t S0 = { .a11 = Bm[0].a11, .a12 = Bm[0].a12, .a21 = Bm[0].a21, .a22 = Bm[0].a22 };
//    matrix_t SN = { .a11 = Bm[N].a11, .a12 = Bm[N].a12, .a21 = Bm[N].a21, .a22 = Bm[N].a22 };
//    S0 = madd(S0, mmmul(Cx, Am[N]));
//    SN = madd(SN, mmmul(Ax, Cm[0]));

//    u1[0] = mvmul(minv(S0), r0);
//    u1[N] = mvmul(minv(SN), rN);

//    for (s = N / 2; s > 0; s /= 2) {
//        for (j = s; j < N; j += 2 * s) {
//            Fm[j] = vadd(Fm[j], 1.0, mvmul(Am[j], u1[j - s]), -1.0);
//            Fm[j] = vadd(Fm[j], 1.0, mvmul(Cm[j], u1[j + s]), -1.0);
//            u1[j] = mvmul(minv(Bm[j]), Fm[j]);
//        }
//    }
//}

//void reduction_method(const node_t* u, node_t* u1)
//{
//    int i;
//    /* Матрицы метода редукции. */
//    matrix_t* Am = (matrix_t*)malloc(sizeof(matrix_t) * (N + 1));
//    matrix_t* Bm = (matrix_t*)malloc(sizeof(matrix_t) * (N + 1));
//    matrix_t* Cm = (matrix_t*)malloc(sizeof(matrix_t) * (N + 1));
//    node_t* Fm = (node_t*)malloc(sizeof(node_t) * (N + 1));

//    if ((u == NULL) || (u1 == NULL)) {
//        fprintf(stderr, "Invalid arguments for reduction_method\n");
//        exit(-1);
//    }

//    if ((Am == NULL) || (Bm == NULL) || (Cm == NULL) || (Fm == NULL)) {
//        fprintf(stderr, "Cannot alloc enough memory\n");
//        exit(-1);
//    }

//    /* Постоянные матрицы. */
//    for (i = 1; i < N; i++) {
//        Am[i] = (matrix_t){ .a11 = dt * d / h / h, .a12 = 0.0, .a21 = 0.0, .a22 = dt * D / h / h };
//        Bm[i] = (matrix_t){ .a11 = -1.0 - 2.0 * dt * d / h / h + dt * dfdu(u[i]),
//                    .a12 = dt * dfdv(u[i]),
//                    .a21 = dt * dgdu(u[i]),
//                    .a22 = -1.0 - 2.0 * dt * D / h / h + dt * dgdv(u[i]) };
//        Cm[i] = (matrix_t){ .a11 = dt * d / h / h, .a12 = 0.0, .a21 = 0.0, .a22 = dt * D / h / h };
//        Fm[i] = (node_t){ .u = -u[i].u - dt * f(u[i]) +
//                    dt * dfdu(u[i]) * u[i].u + dt * dfdv(u[i]) * u[i].v,
//                  .v = -u[i].v - dt * g(u[i]) +
//                    dt * dgdu(u[i]) * u[i].u + dt * dgdv(u[i]) * u[i].v };
//    }

//    /* Граничные условия. */
//    Am[0] = (matrix_t){ .a11 = 0.0, .a12 = 0.0, .a21 = 0.0, .a22 = 0.0 };
//    Bm[0] = (matrix_t){ .a11 = -1.0, .a12 = 0.0, .a21 = 0.0, .a22 = -1.0 };
//    Cm[0] = (matrix_t){ .a11 = 1.0, .a12 = 0.0, .a21 = 0.0, .a22 = 1.0 };
//    Fm[0] = (node_t){ .u = 0, .v = 0 };
//    Am[N] = (matrix_t){ .a11 = 1.0, .a12 = 0.0, .a21 = 0.0, .a22 = 1.0 };
//    Bm[N] = (matrix_t){ .a11 = -1.0, .a12 = 0.0, .a21 = 0.0, .a22 = -1.0 };
//    Cm[N] = (matrix_t){ .a11 = 0.0, .a12 = 0.0, .a21 = 0.0, .a22 = 0.0 };
//    Fm[N] = (node_t){ .u = 0, .v = 0 };

//    reduction(Am, Bm, Cm, Fm, u1);
//    //prl_reduction(Am, Bm, Cm, Fm, u1);
//    int size = N+1;
//    for(i=0;i<size;i++){
//                printf("%f, %f,  F%d\n",Fm[i].u, Fm[i].v,i);
//        }
//        for(i=0;i<N+1;i++){
//                printf("%f, %f,  U%d\n",u1[i].u, u1[i].v,i);
//        }
//        for(i=0;i<size;i++){
//                printf("%f, %f, %f %f  A%d\n",Am[i].a11, Am[i].a12, Am[i].a21, Am[i].a22,i);
//        }
//        for(i=0;i<size;i++){
//                printf("%f, %f, %f %f B%d\n", Bm[i].a11, Bm[i].a12, Bm[i].a21, Bm[i].a22,i);
//        }
//        for(i=0;i<size;i++){
//                printf("%f, %f, %f %f  C%d\n",Cm[i].a11, Cm[i].a12, Cm[i].a21, Cm[i].a22,i);
//        }
//    free(Am);
//    free(Bm);
//    free(Cm);
//    free(Fm);
//}

//int main(int argc, char** argv)
//{
//    node_t* u, * u1;
//    node_t* temp;
//    int i;
//    char buf[256];
//    const char* save[2] = { "u", "v" };

//    u = (node_t*)malloc(sizeof(node_t) * (N + 1));
//    u1 = (node_t*)malloc(sizeof(node_t) * (N + 1));
//    if ((u == NULL) || (u1 == NULL)) {
//        fprintf(stderr, "Cannot alloc enough memory\n");
//        exit(-1);
//    }

//    if (!if_power_of_2(N)) {
//        fprintf(stderr, "Parameter N must be 2^k\n");
//        exit(-1);
//    }

//    init(u);
//    init(u1);

//    for (i = 0; i < 1; i++) {
//        /* Сохраняем посчитанные значения. */
//        if (i % 100 == 0) {
//            sprintf(buf, "C:/Users/Gleb/Desktop/lobanov/S_reaction_diffusion_matrix/src/res/data_%06d.vtk", i);
//            write_to_vtk2((double*)u, buf, save, N, 0.0, h, 2);
//        }

//        /* Обновляем значение. */
//        reduction_method(u, u1);

//        /* Меняем местами временной слой. */
//        temp = u;
//        u = u1;
//        u1 = temp;
//    }

//    free(u);
//    free(u1);
//    return 0;
//}
