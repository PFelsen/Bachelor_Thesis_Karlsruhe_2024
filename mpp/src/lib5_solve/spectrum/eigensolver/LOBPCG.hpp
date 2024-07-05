#ifndef LOBPCG_HPP
#define LOBPCG_HPP

#include "IEigenSolver.hpp"

class LOBPCG : public ILOBPCG {
protected:
  explicit LOBPCG(std::unique_ptr<LinearSolver> &&solver, std::unique_ptr<LinearSolver> &&solver2,
                  int verbose, double eps, int maxstep, int printSteps) :
      ILOBPCG("LOBPCG", 1, std::move(solver), std::move(solver2), verbose, eps, maxstep,
              printSteps) {}

  void ritzStepLOBPCG(Eigenfcts &u, std::vector<Vectors> &additional, Eigenvalues &lambda,
                      Operator &A, Operator &B, Vector &tmp1, Vector &tmp2) override;

  void update(Eigenfcts &u, std::vector<Vectors> &additional, Eigenvalues &lambda, Matrix &A,
              Matrix &B, Vector &residual,
              std::function<void(Vector &, Operator &, Operator &, Vector &)> randVector,
              double defect, int r) override;

  friend class EigenSolverCreator;
};

class LOBPCGExtended : public ILOBPCG {
protected:
  explicit LOBPCGExtended(std::unique_ptr<LinearSolver> &&solver,
                          std::unique_ptr<LinearSolver> &&solver2, int verbose, double eps,
                          int maxstep, int printSteps) :
      ILOBPCG("LOBPCGExtended", 2, std::move(solver), std::move(solver2), verbose, eps, maxstep,
              printSteps) {}

  void ritzStepLOBPCG(Eigenfcts &u, std::vector<Vectors> &additional, Eigenvalues &lambda,
                      Operator &A, Operator &B, Vector &tmp1, Vector &tmp2) override;

  void update(Eigenfcts &u, std::vector<Vectors> &additional, Eigenvalues &lambda, Matrix &A,
              Matrix &B, Vector &residual,
              std::function<void(Vector &, Operator &, Operator &, Vector &)> randVector,
              double defect, int r) override;

  friend class EigenSolverCreator;
};

#endif // LOBPCG_HPP
