#include "reduction.h"
#include <time.h>
#include "vtk.h"
#include "cursedpart.h"
#define THREAD_NUM 8
int main(int argc, char** argv)
{
    int N1 = welcome();
    omp_set_num_threads(THREAD_NUM);

    int sizeInput, sizeN, size;

    int count = 0;
    bool flag;
    sizeInput = N1;
    int maxMatsize = 25;
    double* p = (double*)malloc(sizeof(double) * (maxMatsize));

    for (int k = 0; k < maxMatsize; k++) {

        p[k] = pow(2, k) - 1;

        if (p[k] == sizeInput) {
            size = p[k];
            flag = true;
        }
        else {
            if (p[k] >= sizeInput) {
                count = count + 1;
            }
        }
    }
    sizeN = p[maxMatsize - count];
    if (flag == false) {
        size = sizeN;
    }

    clock_t start = clock();

    int i, j;
    int index1, index2, offset;
    double alpha, gamma;

    node_t* u, * u1;
    node_t* temp;
    char buf[256];
    const char* save[2] = { "u", "v" };

    u = (node_t*)malloc(sizeof(node_t) * (size + 1));
    u1 = (node_t*)malloc(sizeof(node_t) * (size + 1));
    if ((u == NULL) || (u1 == NULL)) {
        fprintf(stderr, "Cannot alloc enough memory\n");
        exit(-1);
    }

    init(u, sizeInput, size);
    init(u1, sizeInput, size);

    for (i = 0; i < 100; i++) { //10

        /* Saving calculated values */
        if (i % 1 == 0) {
           // sprintf(buf, "C:\\Users\\Gleb\\Documents\\parallel-main\\NOv\\res1\\data_%06d.vtk", i);
            sprintf(buf, "C:\\Users\\Gleb\\Documents\\parallel-main1\\NOv\\res1\\data_%06d.vtk", i);
            write_to_vtk2((double*)u, buf, save, N, 0.0, h, 2);
        }

        /* Refreshing the value */
        reduction_method(u, u1, sizeInput, flag, size);
        reduction_method1(u, u1);


        /* Swapping time layer */
        temp = u;
        u = u1;
        u1 = temp;
    }

    free(u);
    free(u1);
    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;

    printWorkTime(elapsed);
    return 0;
}
