#ifndef RITZGALERKIN_HPP
#define RITZGALERKIN_HPP

#include "IEigenSolver.hpp"

class RitzGalerkin : public IEigenSolver {
  explicit RitzGalerkin(std::unique_ptr<LinearSolver> &&solver,
                        std::unique_ptr<LinearSolver> &&solver2, int verbose, double eps,
                        int maxstep, int printSteps) :
      IEigenSolver("RitzGalerkin", std::move(solver), std::move(solver2), verbose, eps, maxstep,
                   printSteps) {}

  void
  operator()(Eigenfcts &u, Eigenvalues &lambda, Matrix &A, Matrix &B,
             std::function<double(Vector &, Vector &, Eigenvalue, Operator &, Operator &, Vector &)>
                 ritzDefect,
             Vector &tmp, int fev = 0) override {
    compute(u, lambda, A, B, ritzDefect, tmp);
  }

  void
  operator()(Eigenfcts &u, Eigenvalues &lambda, Matrix &A, Matrix &B, Operator &P,
             std::function<double(Vector &, Vector &, Eigenvalue, Operator &, Operator &, Vector &)>
                 ritzDefect,
             Vector &tmp, int fev = 0) override {
    compute(u, lambda, A, B, ritzDefect, tmp);
  }

  void
  compute(Eigenfcts &u, Eigenvalues &lambda, Matrix &A, Matrix &B,
          std::function<double(Vector &, Vector &, Eigenvalue, Operator &, Operator &, Vector &)>
              ritzDefect,
          Vector &tmp);

  friend class EigenSolverCreator;
};

#endif // RITZGALERKIN_HPP
