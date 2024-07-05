#include "ParallelMatrix.hpp"

#include <cstdlib>


#undef TRUE
#undef FALSE

#ifdef DCOMPLEX
#define Create_Dense_Matrix zCreate_Dense_Matrix
#define Create_CompCol_Matrix zCreate_CompCol_Matrix
#define Print_CompCol_Matrix zPrint_CompCol_Matrix
#define Print_SuperNode_Matrix zPrint_SuperNode_Matrix
#define gstrs zgstrs
#define gstrf zgstrf
#define SLU_DT SLU_Z

typedef complex<double> doublecomplex;

#include "../lib/superlu/include/slu_zdefs.hpp"
#else
#define Create_CompCol_Matrix dCreate_CompCol_Matrix
#define Create_Dense_Matrix dCreate_Dense_Matrix
#define Print_CompCol_Matrix dPrint_CompCol_Matrix
#define Print_SuperNode_Matrix dPrint_SuperNode_Matrix
#define gstrs dgstrs
#define gstrf dgstrf
#define SLU_DT SLU_D

#include "slu_ddefs.h"

#endif

const char trans = 'N';
const Scalar one = 1;
const Scalar minone = -1;

ParallelMatrix::ParallelMatrix(int r, int c, Communicator &PC, int matrixsize, int maxP) :
    C(PC), a(0), IPIV(0), local_col(0), local_row(0), row(r), col(c), dec(0), C_loc(0),
    other_column(0), other_row(0), other_IPIV(0) {
  if (row == 0 || col == 0) return;
  Define_ParallelDistribution(matrixsize, maxP);
}

ParallelMatrix::ParallelMatrix(Communicator &PC) :
    C(PC), a(0), IPIV(0), local_col(0), local_row(0), row(0), col(0), dec(0), C_loc(0),
    other_column(0), other_row(0), other_IPIV(0) {}

ParallelMatrix::ParallelMatrix(int r, const ParallelMatrix &PM, int matrixsize, int maxP) :
    C(PM.C), a(0), IPIV(0), local_col(PM.local_col), local_row(PM.local_row), row(r), col(PM.col),
    dec(0), C_loc(0), other_column(0), other_row(0), other_IPIV(0) {
  if (row == 0 || col == 0) { return; }
  if (PM.row == 0) {
    Define_ParallelDistribution(matrixsize, maxP);
    return;
  }
  dec.resize(PM.dec.size());
  for (unsigned int i = 0; i < dec.size(); ++i) {
    dec[i].size = PM.dec[i].size;
    dec[i].proc = PM.dec[i].proc;
    dec[i].shift = PM.dec[i].shift;
  }
  int color = 0;
  if (local_col > 0) color = 1;

  //     if (dec.size() < C.size())
  C_loc = new Communicator(C, color);
  //     if (!C_loc) C_loc = &C;
  //     else
  if (color == 0) {
    delete C_loc;
    C_loc = 0;
  }

  local_row = row;
  Create();
}

void ParallelMatrix::Create() {
  if (local_col == 0 || local_row == 0) return;
  if (a) delete[] a;
  a = 0;
  a = new Scalar[local_row * local_col];
  Scalar *t = a;
  for (int i = 0; i < local_col; ++i) {
    for (int j = 0; j < local_row; ++j)
      *(t++) = 0;
  }
}

void ParallelMatrix::Create_Distribution(int MINSIZE, int maxP) {
  int P = C.Size();
  if (maxP > 0 && maxP < P) P = maxP;
  int N = col;
  if (N < P) P = N;
  int p = C.Proc();
  if (MINSIZE < 0) MINSIZE = N;
  if (N < MINSIZE) MINSIZE = N;

  int local_size = int(N / P);
  int num_dec = P;
  if (local_size < MINSIZE) {
    local_size = MINSIZE;
    num_dec = int(N / local_size);
    local_size = int(N / num_dec);
  }

  int extra = N - num_dec * local_size;

  dec.resize(num_dec);
  local_col = 0;
  int color = 0;
  int SHIFT = 0;
  for (int i = 0; i < extra; ++i) {
    dec[i].size = local_size + 1;
    dec[i].proc = i;
    if (p == i) {
      local_col = local_size + 1;
      color = 1;
    }
    dec[i].shift = SHIFT;
    SHIFT += dec[i].size;
  }

  for (int i = extra; i < num_dec; ++i) {
    dec[i].size = local_size;
    dec[i].proc = i;
    if (p == i) {
      local_col = local_size;
      color = 1;
    }
    dec[i].shift = SHIFT;
    SHIFT += dec[i].size;
  }

  if ((int)dec.size() < C.Size()) C_loc = new Communicator(C, color);
}

void ParallelMatrix::Define_ParallelDistribution(int MINSIZE, int maxP) {
  Create_Distribution(MINSIZE, maxP);
  local_row = row;
  Create();
}

void ParallelMatrix::makeLU(int i) {
  if (C.Proc() == i) {
    if (IPIV) { THROW("Error in ParallelMatrix::makeLU") }
    int n = dec[i].size;
    IPIV = new int[n];
    GETRF(&n, &n, a + dec[i].shift, &col, IPIV, &info);
    if (info != 0) { THROW("Error in GETRF (LU)") }
  }
}

void ParallelMatrix::Solve(int i, Scalar *rhs, int nrhs, int size_rhs) {
  if (C.Proc() == i) {
    int n = dec[i].size;
    GETRS(&trans, &n, &nrhs, a + dec[i].shift, &col, IPIV, rhs + dec[i].shift, &size_rhs, &info);
  }
}

void ParallelMatrix::Solve_own(int i) {
  if (!other_column) return;
  pardec dec = get_dec(i);
  GETRS(&trans, &dec.size, &local_col, other_column + dec.shift, &row, other_IPIV, a + dec.shift,
        &row, &info);
}

void ParallelMatrix::MMM_own(int i) {
  if (!other_column) return;
  pardec dec = get_dec(i);
  int n = row - dec.shift - dec.size;
  GEMM(&trans, &trans, &n, &local_col, &dec.size, &minone, other_column + dec.shift + dec.size,
       &row, a + dec.shift, &row, &one, a + dec.shift + dec.size, &row);
}

void ParallelMatrix::MMM_right(Scalar *A_other_column, pardec pardec) {
  int n = row - pardec.shift - pardec.size;
  GEMM(&trans, &trans, &n, &local_col, &pardec.size, &minone,
       A_other_column + pardec.shift + pardec.size, &row, a + pardec.shift, &row, &one,
       a + pardec.shift + pardec.size, &row);
}

void ParallelMatrix::MMM_left(Scalar *in, int Arows, pardec pardec) {
  if (!other_column) return;
  GEMM(&trans, &trans, &row, &local_col, &pardec.size, &minone, other_column, &row, in, &Arows,
       &one, a, &row);
  //     LIB_GEMM(&trans, &trans, &leftrows, &local_col, &dec.size, &minone, leftcolumn, &leftrows,
  //              a+dec.shift, &row, &one, left, &leftrows); //GEMM LEFT
}

void ParallelMatrix::MMM_schur(Scalar *left, Scalar *right, int right_row, pardec pardec) {
  GEMM(&trans, &trans, &row, &local_col, &pardec.size, &minone, left, &row, right + pardec.shift,
       &right_row, &one, a, &row); // SCHUR
}

void ParallelMatrix::Solve_right(Scalar *rightref, int local_col_right, pardec pardec) {
  GETRS(&trans, &pardec.size, &local_col_right, other_column + pardec.shift, &row, other_IPIV,
        rightref + pardec.shift, &row, &info); // Solve RIGHT
}

void ParallelMatrix::MV_rhs(int i, Scalar *rhs, int nrhs, int size_rhs) { // call with A->
  pardec dec = get_dec(i);
  int n = row - dec.shift - dec.size;
  GEMM(&trans, &trans, &n, &nrhs, &dec.size, &minone, a + dec.shift + dec.size, &row,
       rhs + dec.shift, &size_rhs, &one, rhs + dec.shift + dec.size, &size_rhs);
}

void ParallelMatrix::MV_rhs_left(int i, int Arow, Scalar *rhs, int nrhs,
                                 int size_rhs) { // call with Left->
  pardec dec = get_dec(i);
  GEMM(&trans, &trans, &row, &nrhs, &dec.size, &minone, a, &row, rhs + dec.shift, &size_rhs, &one,
       rhs + Arow, &size_rhs);
}

void ParallelMatrix::MV_rhs_right(int p, Scalar *rhs, int nrhs, int size_rhs) {
  if (p < get_dec_size()) {
    pardec dec = get_dec(p);
    GEMM(&trans, &trans, &row, &nrhs, &local_col, &minone, a, &row, rhs + row + dec.shift,
         &size_rhs, &one, rhs, &size_rhs);
  }
}

void ParallelMatrix::MV_rhs_Z(int i, Scalar *rhs, int nrhs, int size_rhs) {
  pardec dec = get_dec(i);
  GEMM(&trans, &trans, &dec.shift, &nrhs, &local_col, &minone, a, &row, rhs + dec.shift, &size_rhs,
       &one, rhs, &size_rhs);
}

void ParallelMatrix::Create_other_column(int i, bool withIPIV) {
  if (other_column) delete[] other_column;
  other_column = 0;
  if (other_IPIV) delete[] other_IPIV;
  other_IPIV = 0;
  size_t sz = dec[i].size * row;
  other_column = new Scalar[sz];
  if (withIPIV) other_IPIV = new int[dec[i].size];
}

void ParallelMatrix::Communicate_Column_global(int i, bool withIPIV) {
  Create_other_column(i, withIPIV);
  size_t sz = dec[i].size * row;
  if (C.Proc() == i) {
    memcpy(other_column, a, sz * sizeof(Scalar));
    if (withIPIV) memcpy(other_IPIV, IPIV, dec[i].size * sizeof(int));
  }
  C.Broadcast(other_column, sz, i);
  if (withIPIV) C.Broadcast(other_IPIV, dec[i].size, i);
}

void ParallelMatrix::Communicate_Column_intern(int i, Communicator *C_locloc, bool withIPIV) {
  if (!C_locloc) {
    Communicate_Column_global(i, withIPIV);
    return;
  }
  if (C.Proc() != C_locloc->Proc()) return;
  Create_other_column(i, withIPIV);
  size_t sz = dec[i].size * row;
  if (C_locloc->Proc() == i) {
    memcpy(other_column, a, sz * sizeof(Scalar));
    if (withIPIV) memcpy(other_IPIV, IPIV, dec[i].size * sizeof(int));
  }
  C_locloc->Broadcast(other_column, sz, i);
  if (withIPIV) C_locloc->Broadcast(other_IPIV, dec[i].size, i);
}

void ParallelMatrix::Create_other_row(int num_rows) {
  if (other_row) delete[] other_row;
  other_row = 0;
  size_t sz = cols() * (num_rows);
  other_row = new Scalar[sz];
}

void ParallelMatrix::Communicate_row(int r_begin, int r_end, Communicator *PC) {
  int num_rows = r_end - r_begin;
  Create_other_row(num_rows);
  if (!PC) return;
  if (C.Proc() >= PC->Size()) return;

  if (PC->Size() == 1) {
    for (int i = 0; i < local_col; ++i) {
      memcpy(other_row + i * num_rows, a + r_begin + i * rows(), num_rows * sizeof(Scalar));
    }
    return;
  }

  Scalar *senddata = new Scalar[local_col * num_rows];
  size_t sendsize = 0;
  vector<int> receivesize(PC->Size(), 0);
  vector<int> displs(PC->Size(), 0);
  Scalar *receivedata = other_row;
  for (unsigned int p = 0; p < dec.size(); ++p) {
    receivesize[p] = num_rows * dec[p].size;
    displs[p] = num_rows * dec[p].shift;
  }
  for (int i = 0; i < local_col; ++i) {
    memcpy(senddata + i * num_rows, a + r_begin + i * rows(), num_rows * sizeof(Scalar));
  }
  PC->Allgatherv(senddata, sendsize, receivedata, &receivesize[0], &displs[0]);

  delete[] senddata;
}

void ParallelMatrix::make_LU_multiply_local(Scalar *left, int r_begin, int r_end) {
  if (!a) return;
  int r = row - r_end;
  int k = r_end - r_begin;
  int c = local_col;
  const char dotrans = 'T';
  Scalar *B = a + r_begin;
  Scalar *C = a + r_end;
  GEMM(&dotrans, &trans, &r, &c, &k, &minone, left, &k, B, &row, &one, C, &row);
}

void ParallelMatrix::makeLU() {
  if (!a) return;
  for (int i = 0; i < (int)dec.size(); ++i) {
    makeLU(i);
    if (C_loc->Proc() == C.Proc()) Communicate_Column_intern(i, C_loc);
    if (C.Proc() > i) {
      GETRS(&trans, &dec[i].size, &local_col, other_column + dec[i].shift, &row, other_IPIV,
            a + dec[i].shift, &row, &info);
      int n = row - dec[i].shift - dec[i].size;
      GEMM(&trans, &trans, &n, &local_col, &dec[i].size, &minone,
           other_column + dec[i].shift + dec[i].size, &row, a + dec[i].shift, &row, &one,
           a + dec[i].shift + dec[i].size, &row);
    }
  }
}

void ParallelMatrix::set(int i, int j, Scalar val) {
  if (!a) return;
  int p = C.Proc();
  if (p >= (int)dec.size()) return;
  if ((j < dec[p].shift) || (j >= dec[p].shift + local_col)) return;
  *(a + (j - dec[p].shift) * row + i) = val;
}

void ParallelMatrix::add(int i, int j, Scalar val) {
  if (!a) return;
  int p = C.Proc();
  if (p >= (int)dec.size()) return;
  if ((j < dec[p].shift) || (j >= dec[p].shift + local_col)) return;
  *(a + (j - dec[p].shift) * row + i) += val;
}

Scalar *ParallelMatrix::operator()(int i, int j) {
  if (!a) return NULL;
  int p = C.Proc();
  if (p >= (int)dec.size()) return NULL;
  if ((j < dec[p].shift) || (j >= dec[p].shift + local_col)) return NULL;
  return a + (j - dec[p].shift) * row + i;
}

Scalar *ParallelMatrix::loc(int i, int j) {
  if (!a) return NULL;
  return a + j * row + i;
}

void ParallelMatrix::loc_add(int i, int j, Scalar val) {
  if (!a) return;
  *(a + j * row + i) += val;
}

Scalar *ParallelMatrix::operator()(int i, int j) const {
  if (!a) return NULL;
  int p = C.Proc();
  if (p >= (int)dec.size()) return NULL;
  if ((j < dec[p].shift) || (j >= dec[p].shift + local_col)) return NULL;
  return a + (j - dec[p].shift) * row + i;
}

Scalar *ParallelMatrix::loc(int i, int j) const {
  if (!a) return NULL;
  return a + j * row + i;
}
