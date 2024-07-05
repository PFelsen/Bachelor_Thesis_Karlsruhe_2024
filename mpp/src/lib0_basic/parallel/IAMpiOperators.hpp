#ifndef MPIOPERATORS_HPP
#define MPIOPERATORS_HPP

#include <mpi.h>

class IAInterval;

class IACInterval;

extern MPI_Op MPI_SUM_IAINTERVAL;
extern MPI_Op MPI_SUM_IACINTERVAL;
extern MPI_Op MPI_HULL_IAINTERVAL;
extern MPI_Op MPI_INTERSECT_IAINTERVAL;

void mpi_sum_iainterval(double *in, double *inout, int *len, MPI_Datatype *dtype);

void mpi_sum_iacinterval(double *in, double *inout, int *len, MPI_Datatype *dtype);

void mpi_hull_iainterval(double *in, double *inout, int *len, MPI_Datatype *dtype);

void mpi_intersect_iainterval(double *in, double *inout, int *len, MPI_Datatype *dtype);

#endif // MPIOPERATORS_HPP
