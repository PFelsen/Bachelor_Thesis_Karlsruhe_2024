#ifndef LOBPCGSELECTIVE_HPP
#define LOBPCGSELECTIVE_HPP

#include "IEigenSolver.hpp"

class LOBPCGSelective : public ILOBPCGSelective {
protected:
  explicit LOBPCGSelective(std::unique_ptr<LinearSolver> &&solver,
                           std::unique_ptr<LinearSolver> &&solver2, int verbose, double eps,
                           int maxstep, int printSteps) :
      ILOBPCGSelective("LOBPCGSelective", 2, std::move(solver), std::move(solver2), verbose, eps,
                       maxstep, printSteps) {}

  void ritzStepLOBPCG(Eigenfcts &u, std::vector<Vectors> &additional, Eigenvalues &lambda,
                      Operator &A, Operator &B, Vector &tmp1, Vector &tmp2) override;

  bool
  update(Eigenfcts &u, std::vector<Vectors> &additional, Eigenvalues &lambda, Matrix &A, Matrix &B,
         Vector &residual, Vector &tmp,
         std::function<void(Vector &, Operator &, Operator &, Vector &)> randVector,
         std::function<double(Vector &, Vector &, Eigenvalue, Operator &, Operator &, Vector &)>
             ritzDefect,
         int fev) override;

  friend class EigenSolverCreator;
};

class LOBPCGSelectiveFast : public ILOBPCGSelective {
protected:
  explicit LOBPCGSelectiveFast(std::unique_ptr<LinearSolver> &&solver,
                               std::unique_ptr<LinearSolver> &&solver2, int verbose, double eps,
                               int maxstep, int printSteps) :
      ILOBPCGSelective("LOBPCGSelectiveFast", 2, std::move(solver), std::move(solver2), verbose,
                       eps, maxstep, printSteps) {}

  void ritzStepLOBPCG(Eigenfcts &u, std::vector<Vectors> &additional, Eigenvalues &lambda,
                      Operator &A, Operator &B, Vector &tmp1, Vector &tmp2) override;

  bool
  update(Eigenfcts &u, std::vector<Vectors> &additional, Eigenvalues &lambda, Matrix &A, Matrix &B,
         Vector &residual, Vector &tmp,
         std::function<void(Vector &, Operator &, Operator &, Vector &)> randVector,
         std::function<double(Vector &, Vector &, Eigenvalue, Operator &, Operator &, Vector &)>
             ritzDefect,
         int fev) override;

  friend class EigenSolverCreator;
};

class LOBPCGSelectiveVeryFast : public ILOBPCGSelective {
protected:
  explicit LOBPCGSelectiveVeryFast(std::unique_ptr<LinearSolver> &&solver,
                                   std::unique_ptr<LinearSolver> &&solver2, int verbose, double eps,
                                   int maxstep, int printSteps) :
      ILOBPCGSelective("LOBPCGSelectiveVeryFast", 2, std::move(solver), std::move(solver2), verbose,
                       eps, maxstep, printSteps) {}

  void setMatricesSmallVeryFast(Eigenfcts &u, Vectors &a1, Operator &A, Operator &B, Vector &tmp1,
                                Vector &tmp2);

  void setMatricesLargeVeryFast(Eigenfcts &u, Vectors &a1, Vectors &a2, Operator &A, Operator &B,
                                Vector &tmp1, Vector &tmp2);

  void ritzStepLOBPCG(Eigenfcts &u, std::vector<Vectors> &additional, Eigenvalues &lambda,
                      Operator &A, Operator &B, Vector &tmp1, Vector &tmp2) override;

  bool
  update(Eigenfcts &u, std::vector<Vectors> &additional, Eigenvalues &lambda, Matrix &A, Matrix &B,
         Vector &residual, Vector &tmp,
         std::function<void(Vector &, Operator &, Operator &, Vector &)> randVector,
         std::function<double(Vector &, Vector &, Eigenvalue, Operator &, Operator &, Vector &)>
             ritzDefect,
         int fev) override;

  friend class EigenSolverCreator;
};


#endif // LOBPCGSELECTIVE_HPP
