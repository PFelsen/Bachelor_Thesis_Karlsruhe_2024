#ifndef _PARALLELSOLVERMATRIX_H_
#define _PARALLELSOLVERMATRIX_H_

#include "ParallelMatrix.hpp"
#include "ParallelSolverSteps.hpp"

class ParallelSolverMatrix {
protected:
  const ParallelSolverClusterStep &cluster; // reference to the cluster_step
  const ParallelSolverOneStep &step;        // reference to one step -- needed for SetNextMatrix...
  Communicator &C;
  Communicator *C_loc;

  const VectorProcSet &vps; // reference to vecprocset in cluster_step
  ParallelMatrix *A;
  ParallelMatrix *Right;
  ParallelMatrix *Left;
  ParallelMatrix *Schur;

  Scalar *rhs;
  int size_rhs;
  int nrhs;
  int info;
private:
  void Copy_Schur_Matrix(ParallelSolverMatrix &next);

  void COPY(Scalar *data, size_t size, size_t numrow, Scalar *tmp);

  void SEND(Communicator &Comm, Scalar *data, size_t size, size_t numrow, int dest);

  void RECEIVE(Communicator &Comm, Scalar *tmp, size_t size, size_t numrow, int source);

  void Set_LEFT(ParallelSolverMatrix &next, const ParallelSolverClusterStep *current_cluster,
                const Scalar *tmp, const int numrows, const int numcols, const int shift);

  void set_sr_dec_LEFT(ParallelSolverMatrix &next, const ParallelSolverClusterStep *current_cluster,
                       const std::vector<pardec> *current_dec, int &proc_shift, int &shift_send,
                       int &size, int &shift_receive);

  void SetNext_Matrix_LEFT_CLUSTERED(ParallelSolverMatrix &next,
                                     const ParallelSolverClusterStep *current_cluster,
                                     const std::vector<pardec> *current_dec,
                                     int &current_proc_shift);

  void Set_RIGHT(ParallelSolverMatrix &next, const ParallelSolverClusterStep *current_cluster,
                 const Scalar *tmp, const int numrows, const int numcols, const int shift);

  void set_sr_dec_RIGHT(ParallelSolverMatrix &next,
                        const ParallelSolverClusterStep *current_cluster,
                        const std::vector<pardec> *current_dec, int &proc_shift, int &shift_send,
                        int &size, int &shift_receive);

  void SetNext_Matrix_RIGHT_CLUSTERED(ParallelSolverMatrix &next,
                                      const ParallelSolverClusterStep *current_cluster,
                                      const std::vector<pardec> *current_dec,
                                      int &current_proc_shift);

  void Communicate_other_dec(const ParallelSolverMatrix &next, std::vector<pardec> &other_dec,
                             int shift_own_proc, int shift_other_proc, int &other_cluster_number);
public:
  ParallelSolverMatrix(const ParallelSolverOneStep &one_step, int matrixsize = 0, int maxP = 0);

  void Set(const SparseMatrix &S) { Set(S.nzval(), S.rowind(), S.colptr(), S.size()); }

  virtual void Set(const Scalar *a, const int *d, const int *col, int n);

  void Set_rhs(const Scalar *u, int shift = 0);

  void Create_rhs(int NRHS = 1);

  void Write_rhs(Scalar *u, int shift = 0);

  int min(int a, int b) {
    if (a < b) return a;
    else return b;
  }

  void SetNext_Matrix(ParallelSolverMatrix &next);

  virtual void makeLU();

  void SetNext_rhs_LEFT(ParallelSolverMatrix &next);

  void SetNext_rhs_RIGHT(ParallelSolverMatrix &next);

  virtual void SolveL();

  virtual void SolveU();

  virtual void Destruct() {}

  virtual ~ParallelSolverMatrix();

  ParallelMatrix *getA() const { return A; }

  ParallelMatrix *getRight() const { return Right; }

  ParallelMatrix *getLeft() const { return Left; }

  ParallelMatrix *getSchur() const { return Schur; }

  int solsize() { return A->rows(); }

  int schursize() { return Right->cols(); }

  virtual int get_LU_procs();

  int get_Schur_procs();

  int get_procs() { return C.Size(); }
};

class ParallelSolverMatrix_Sparse : public ParallelSolverMatrix {
  SparseRMatrix *SparseA;
  SparseSolver *S;
  SparseRMatrix *SparseLeft;
public:
  ParallelSolverMatrix_Sparse(const ParallelSolverOneStep &one_step) :
      ParallelSolverMatrix(one_step, -1), SparseA(0), S(0), SparseLeft(0) {}

  void Set(const Scalar *, const int *, const int *, int);

  void makeLU();

  void SolveL();

  void SolveU();

  void Destruct();

  ~ParallelSolverMatrix_Sparse() { Destruct(); }

  virtual int get_LU_procs() { return 1; }
};


#endif // _PARALLELSOLVERMATRIX_H_
