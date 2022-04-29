#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>

#include "C:\Users\ilina\source\repos\Parallel\vtk.h"
#include "C:\Users\ilina\source\repos\Parallel\common.h"
#include "C:\Users\ilina\source\repos\Parallel\parallel.h"
#include "C:\Users\ilina\source\repos\Parallel\blood.h" //������� ���� ���� ��� ��������������

/*

#ifdef BLOOD
#include "C:\Users\ilina\source\repos\Parallel\blood.h"
#endif

#ifdef BRUSSELATOR
#include "C:\Users\ilina\source\repos\Parallel\brusselator.h"
#endif

*/


int if_power_of_2(int n) { //�������, �����������, �������� �� �������� �� ��� ����� �������� ������
	while ((n & 1) == 0)
		n >>= 1;
	return (n == 1) ? 1 : 0;
}

void prl_reduction(matrix_t* Am, matrix_t* Bm, matrix_t* Cm, node_t* Fm, node_t* u1)
{
	int s, j;
	// for each proccesses
	for (s = 1; s < N_prl; s *= 2) {
		// define boundary values
		// for 0 and last 0 and N
		// for others mpi_send for right boundary
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

	// Gather on Master all last points
	// and call reduction (Am_last, Bm_last, ... , u[0]);

/*��������� ����*/
	matrix_t Am_last = Am[N];
	matrix_t Bm_last = Bm[N];
	matrix_t Cm_last = Cm[N];
	node_t Fm_last = Fm[N];
/*�� ����, ������*/

	if (rank == MASTER) {
		//���� reduction(Am_last, Bm_last, Cm_last, Fm_last, u1[0]);
		reduction (&Am_last, &Bm_last, &Cm_last, &Fm_last, &u1[0]);
	}
	// and Scatter u[p+1] on all proccesses

	// rebuild the solution
	// N_prl = N / p (p - number of proccesses)
	// for master (rank = 0) N_prl + 1 (+ 0)

	for (s = N_prl / 2; s > 0; s /= 2) {
		for (j = s; j < N_prl; j += 2 * s) {
			Fm[j] = vadd(Fm[j], 1.0, mvmul(Am[j], u1[j - s]), -1.0);
			Fm[j] = vadd(Fm[j], 1.0, mvmul(Cm[j], u1[j + s]), -1.0);
			u1[j] = mvmul(minv(Bm[j]), Fm[j]);
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

void reduction_method(const node_t* u, node_t* u1)
{
	int i;
	/* ������� ������ ��������. */
	matrix_t* Am = (matrix_t*)malloc(sizeof(matrix_t) * (N + 1));
	matrix_t* Bm = (matrix_t*)malloc(sizeof(matrix_t) * (N + 1));
	matrix_t* Cm = (matrix_t*)malloc(sizeof(matrix_t) * (N + 1));
	node_t* Fm = (node_t*)malloc(sizeof(node_t) * (N + 1));

//� ����� �� ������ ��� ��������, ���� �� ������ ��� �� ������?
	if ((u == NULL) || (u1 == NULL)) {
		fprintf(stderr, "Invalid arguments for reduction_method\n");
		exit(-1);
	}
//

	if ((Am == NULL) || (Bm == NULL) || (Cm == NULL) || (Fm == NULL)) {
		fprintf(stderr, "Cannot alloc enough memory\n");
		exit(-1);
	}

	/* ���������� �������. */
	for (i = 1; i < N; i++)
	{
		Am[i] = (matrix_t)
		{ 
			                .a11 = dt * d / h / h,   .a12 = 0.0,

			                .a21 = 0.0,			     .a22 = dt * D / h / h
		};
		//��� �������� ���� ������� ����� �������� ���� �� ��������� � ���, ��� �������� � ���������. ��� ��� �����������? �� ������ �� ������.


		Bm[i] = (matrix_t)
		{
			                .a11 = -1.0 - 2.0 * dt * d / h / h + dt * dfdu(u[i]),    .a12 = dt * dfdv(u[i]),

					        .a21 = dt * dgdu(u[i]),                                  .a22 = -1.0 - 2.0 * dt * D / h / h + dt * dgdv(u[i])
		};
		//��� ���� ����� ��������


		Cm[i] = (matrix_t)
		{ 
			                .a11 = dt * d / h / h,   .a12 = 0.0,

			                .a21 = 0.0,              .a22 = dt * D / h / h
		};
		//� ����� ������� �� �������


		Fm[i] = (node_t)
		{
			        .u = -u[i].u - dt * f(u[i]) + dt * dfdu(u[i]) * u[i].u + dt * dfdv(u[i]) * u[i].v,
				    .v = -u[i].v - dt * g(u[i]) + dt * dgdu(u[i]) * u[i].u + dt * dgdv(u[i]) * u[i].v
		};
		//��� ������ ���-�� ��� �� ������ ���, ��� ����. h^2, �� ������� �� ����� ����� ����� ���������, ���, ��������, ������ ���.
	}

	/* ��������� �������. */
	Am[0] = (matrix_t){ .a11 = 0.0, .a12 = 0.0, .a21 = 0.0, .a22 = 0.0 };
	Bm[0] = (matrix_t){ .a11 = -1.0, .a12 = 0.0, .a21 = 0.0, .a22 = -1.0 };
	Cm[0] = (matrix_t){ .a11 = 1.0, .a12 = 0.0, .a21 = 0.0, .a22 = 1.0 };
	Fm[0] = (node_t){ .u = 0, .v = 0 };
	Am[N] = (matrix_t){ .a11 = 1.0, .a12 = 0.0, .a21 = 0.0, .a22 = 1.0 };
	Bm[N] = (matrix_t){ .a11 = -1.0, .a12 = 0.0, .a21 = 0.0, .a22 = -1.0 };
	Cm[N] = (matrix_t){ .a11 = 0.0, .a12 = 0.0, .a21 = 0.0, .a22 = 0.0 };
	Fm[N] = (node_t){ .u = 0, .v = 0 };

	//reduction(Am, Bm, Cm, Fm, u1);
	prl_reduction(Am, Bm, Cm, Fm, u1);

	free(Am);
	free(Bm);
	free(Cm);
	free(Fm);
}

int main(int argc, char** argv)
{
	node_t* u, * u1; //u - ��� u ��� ��������, �� ���� � ������ ������ �������, u1 - ��� u � ���������, �� ���� � ����� ������ �������
	node_t* temp;
	int i;
	char buf[256];
	const char* save[2] = { "u", "v" };

	u = (node_t*)malloc(sizeof(node_t) * (N + 1)); //�������� ������ ��� ������� �������� � ������ � ����� ������ �������
	u1 = (node_t*)malloc(sizeof(node_t) * (N + 1));

	if ((u == NULL) || (u1 == NULL)) {
		fprintf(stderr, "Cannot alloc enough memory\n");
		exit(-1);
	}

	if (!if_power_of_2(N)) {
		fprintf(stderr, "Parameter N must be 2^k\n"); //N - ���������� ������� �������. h - ������ ����
		exit(-1);
	}

	init(u); //������ ��������� ������� (��. brusselator.h)
	init(u1);

	for (i = 0; i < S; i++) { //����� ����� S = 50 000
		/* ������ ��� ����� ��������� ����������� ��������. */
		if (i % 100 == 0) {
			sprintf(buf, "C:/Users/ilina/source/repos/Parallel/res/data_%06d.vtk", i); //%06d: "6" - ���.����� ������ = 6. "0" - ���� ���.����� ������ 6, ����� ���������� ����� ������. �-� ��� i = 100 ������� "res/data_000100.vtk"
			write_to_vtk2((double*)u, buf, save, N, 0.0, h, 2); //���������� ��� ��� � ����
		}

		/* ��������� ��������. */
		reduction_method(u, u1);

		/* ������ ������� ��������� ����. */
		temp = u;
		u = u1;
		u1 = temp;
	}

	free(u);
	free(u1);
	return 0;
}