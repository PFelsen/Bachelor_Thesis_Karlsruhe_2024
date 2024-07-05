#ifndef PARALLELMATRIX_HPP
#define PARALLELMATRIX_HPP

#include "Communicator.hpp"
#include "Sparse.hpp"

#include "Lapdef.hpp"

struct pardec {
  int size;
  int proc;
  int shift;
};

class ParallelMatrix {
protected:
  Communicator &C; // nur Referenz
  Scalar *a;       // local Matrix
  int *IPIV;       // perhaps std::vector<int*> for different IPIVs
  int local_col;   // sum = col
  int local_row;
  int row;
  int col;
  //     std::vector <double*> local_dec;
  std::vector<pardec> dec; // decomposition

  Communicator *C_loc; // falls weniger kommuniziert werden soll...

  Scalar *other_column; // for Communication
  Scalar *other_row;
  int *other_IPIV;
  int info;

  void Create_Distribution(int MINSIZE, int maxP);

  virtual void Create();
public:
  ParallelMatrix(int r, int c, Communicator &PC, int matrixsize = 0, int maxP = 0);

  ParallelMatrix(int r, const ParallelMatrix &PM, int matrixsize = 0, int maxP = 0);

  ParallelMatrix(Communicator &PC);

  virtual void Define_ParallelDistribution(int MINSIZE = 0, int maxP = 0);

  int rows() const { return row; }

  int cols() const { return col; }

  int cols_loc() const { return local_col; }

  virtual void makeLU(int i);

  virtual void Solve(int i, Scalar *rhs, int nrhs = 1, int size_rhs = 1);

  void Solve_own(int);

  void MMM_own(int); // Matrix-Matrix-Multiplication

  void MMM_left(Scalar *, int, pardec);

  void MMM_right(Scalar *, pardec);

  void MMM_schur(Scalar *, Scalar *, int, pardec);

  void MV_rhs(int, Scalar *, int, int);

  void MV_rhs_left(int, int, Scalar *, int, int);

  void MV_rhs_right(int, Scalar *, int, int);

  void MV_rhs_Z(int, Scalar *, int, int);

  void Solve_right(Scalar *, int, pardec);

  const std::vector<pardec> *get_dec() { return &dec; }

  pardec &get_dec(int i) { return dec[i]; }

  int get_dec_size() const { return dec.size(); }

  int *get_other_IPIV() { return other_IPIV; }

  Scalar *get_other_column() { return other_column; }

  Scalar *get_other_row() { return other_row; }

  Communicator &get_PC() { return C; }

  Communicator &get_PC() const { return C; }

  Communicator *get_PC_loc() const { return C_loc; }

  virtual void Create_other_column(int i, bool withIPIV = true);

  virtual void Communicate_Column_global(int i, bool withIPIV = true);

  virtual void Communicate_Column_intern(int i, Communicator *C_locloc, bool withIPIV = true);

  void Create_other_row(int num_rows);

  void Communicate_row(int r_begin, int r_end, Communicator *PC = NULL);

  void make_LU_multiply_local(Scalar *left, int r_begin, int r_end);

  virtual void makeLU();

  void Copy(const ParallelMatrix &CP) {
    if (CP.a) memcpy(a, CP.a, local_col * row * sizeof(Scalar));
  }

  void Clear() {
    if (a) delete[] a;
    a = 0;
    if (IPIV) delete[] IPIV;
    IPIV = 0;
    if (other_column) delete[] other_column;
    other_column = 0;
    if (other_row) delete[] other_row;
    other_row = 0;
    if (other_IPIV) delete[] other_IPIV;
    other_IPIV = 0;
  }

  virtual ~ParallelMatrix() {
    Clear();
    if (C_loc) delete C_loc;
    C_loc = 0;
  }

  Scalar *ref() { return a; }

  void set(int i, int j, Scalar val);

  void add(int i, int j, Scalar val);

  Scalar *operator()(int i, int j);

  Scalar *loc(int i, int j);

  void loc_add(int i, int j, Scalar val);

  Scalar *operator()(int i, int j) const;

  Scalar *loc(int i, int j) const;
};

#endif // PARALLELMATRIX_HPP
