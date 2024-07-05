#include "IAMpiOperators.hpp"

#include "IACInterval.hpp"

MPI_Op MPI_SUM_IAINTERVAL;
MPI_Op MPI_SUM_IACINTERVAL;
MPI_Op MPI_HULL_IAINTERVAL;
MPI_Op MPI_INTERSECT_IAINTERVAL;

void mpi_sum_iainterval(double *in, double *inout, int *len, MPI_Datatype *dtype) {
  for (int i = 0; i < *len / 2; ++i) {
    IAInterval tmp =
        IAInterval(in[2 * i], in[2 * i + 1]) + IAInterval(inout[2 * i], inout[2 * i + 1]);
    inout[2 * i] = tmp.inf();
    inout[2 * i + 1] = tmp.sup();
  }
}

void mpi_sum_iacinterval(double *in, double *inout, int *len, MPI_Datatype *dtype) {
  for (int i = 0; i < *len / 4; ++i) {
    IACInterval tmp =
        IACInterval(IAInterval(in[4 * i], in[4 * i + 1]), IAInterval(in[4 * i + 2], in[4 * i + 3]))
        + IACInterval(IAInterval(inout[4 * i], inout[4 * i + 1]),
                      IAInterval(inout[4 * i + 2], inout[4 * i + 3]));
    inout[4 * i] = tmp.real().inf();
    inout[4 * i + 1] = tmp.real().sup();
    inout[4 * i + 2] = tmp.imag().inf();
    inout[4 * i + 3] = tmp.imag().sup();
  }
}

void mpi_hull_iainterval(double *in, double *inout, int *len, MPI_Datatype *dtype) {
  for (int i = 0; i < *len / 2; ++i) {
    IAInterval tmp =
        IAInterval(in[2 * i], in[2 * i + 1]) | IAInterval(inout[2 * i], inout[2 * i + 1]);
    inout[2 * i] = tmp.inf();
    inout[2 * i + 1] = tmp.sup();
  }
}

void mpi_intersect_iainterval(double *in, double *inout, int *len, MPI_Datatype *dtype) {
  for (int i = 0; i < *len / 2; ++i) {
    IAInterval tmp =
        IAInterval(in[2 * i], in[2 * i + 1]) & IAInterval(inout[2 * i], inout[2 * i + 1]);
    inout[2 * i] = tmp.inf();
    inout[2 * i + 1] = tmp.sup();
  }
}