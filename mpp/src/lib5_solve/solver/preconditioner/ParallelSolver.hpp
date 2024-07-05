#ifndef _PARALLELSOLVER_H_
#define _PARALLELSOLVER_H_

#include "ParallelSolverMatrix.hpp"
#include "Preconditioner.hpp"

class ParallelSolver : public Preconditioner {
  ParallelSolverAllSteps *steps;
  vector<ParallelSolverMatrix *> PSM;

  int min_matrix_size;
  int maxP;
  bool PS_cd;
public:
  ParallelSolver();

  void Construct(const Matrix &);

  void Destruct();

  ~ParallelSolver() { Destruct(); }

  void multiply(Vector &, const Vector &) const override;

  void multiply(Vectors &U, const Vectors &B) const override;

  string Name() const override { return "ParallelMatrixSolver"; }

  friend std::ostream &operator<<(std::ostream &s, const ParallelSolver &ds) {
    return s << ds.Name();
  }
};

#endif // of #ifndef _PARALLELSOLVER_H_
