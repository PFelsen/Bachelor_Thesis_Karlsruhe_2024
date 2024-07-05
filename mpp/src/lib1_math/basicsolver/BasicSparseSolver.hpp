#ifndef BASICSPARSESOLVER_HPP
#define BASICSPARSESOLVER_HPP

#include "GlobalDefinitions.hpp"
#include "SparseRMatrix.hpp"

class SparseSolver {
public:
  virtual void Solve(Scalar *, int nrhs = 1){};

  virtual void Solve(Scalar *, const Scalar *){};

  virtual ~SparseSolver() {}

  virtual void test_output(){};
};

SparseSolver *GetSparseSolver(SparseRMatrix &, const string &name = "SuperLU");

#endif // BASICSPARSESOLVER_HPP
