#ifndef VTK_H
#define VTK_H

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>

void write_to_vtk(double* u, const char* path, const char* data, const int _N, const double o, const double _h);
void write_to_vtk2(double* u, const char* path, const char** data, const int _N, const double o, const double _h, const int k);
void write_to_vtk2d(double* u, const char* path, const char** data, const int _N[2], const double o[2], const double _h[2], const int gs, const int k);

#endif //VTK_H
